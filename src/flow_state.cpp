/*
  FlowState — residual capacities and placement over a fixed graph.

  Maintains per-edge residual capacity and cumulative edge flows. Supports
  two placement strategies when pushing flow along an SPF DAG:
    - Proportional: distribute flow proportionally to residual capacity,
      processing nodes in topological order from source to destination.
    - EqualBalanced: distribute flow equally across available parallel edges,
      *single-pass ECMP admission* on a fixed DAG. We compute a single global
      scale so no edge is oversubscribed, then return. Re-running on updated
      residuals intentionally changes the allowed next-hop set (progressive/TE);
      use place_max_flow() if you want that behavior.
*/
#include "netgraph/core/flow_state.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/constants.hpp"
#include "netgraph/core/profiling.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>
#include <stdexcept>

namespace netgraph::core {
namespace {

// Dinic-style edge for the reversed residual graph (used in Proportional placement).
// Each edge u->v represents a group of parallel forward edges with aggregated capacity.
struct DinicEdge {
  std::int32_t to;        // destination node
  std::int32_t rev;       // index of reverse edge in adj[to]
  double cap;             // current residual capacity
  double init_cap;        // initial capacity (snapshot at augmentation start)
  std::int32_t group;     // index into groups array (-1 if reverse edge)
};

// Workspace for Dinic-like augmentations on the reversed residual graph.
// We reverse the DAG (dst becomes source) to enable topological flow placement.
struct FlowWorkspace {
  std::vector<std::vector<DinicEdge>> adj; // reversed residual graph
  std::vector<std::int32_t> level;         // BFS level for level graph
  std::vector<std::int32_t> it;            // DFS iteration pointer per node

  void reset(std::int32_t n) {
    adj.assign(static_cast<std::size_t>(n), {});
    level.assign(static_cast<std::size_t>(n), -1);
    it.assign(static_cast<std::size_t>(n), 0);
  }
  void add_edge(std::int32_t u, std::int32_t v, double c, std::int32_t group_idx) {
    // Add forward edge u->v and its reverse v->u (with zero initial capacity).
    // Forward edge stores the group index for later proportional distribution.
    DinicEdge a{v, static_cast<std::int32_t>(adj[static_cast<std::size_t>(v)].size()), c, c, group_idx};
    DinicEdge b{u, static_cast<std::int32_t>(adj[static_cast<std::size_t>(u)].size()), 0.0, 0.0, -1};
    adj[static_cast<std::size_t>(u)].push_back(a);
    adj[static_cast<std::size_t>(v)].push_back(b);
  }
  bool bfs(std::int32_t s, std::int32_t t) {
    std::fill(level.begin(), level.end(), -1);
    std::queue<std::int32_t> q;
    level[static_cast<std::size_t>(s)] = 0;
    q.push(s);
    while (!q.empty()) {
      auto u = q.front(); q.pop();
      for (auto const& e : adj[static_cast<std::size_t>(u)]) {
        if (e.cap >= kMinCap && level[static_cast<std::size_t>(e.to)] < 0) {
          level[static_cast<std::size_t>(e.to)] = level[static_cast<std::size_t>(u)] + 1;
          q.push(e.to);
        }
      }
    }
    return level[static_cast<std::size_t>(t)] >= 0;
  }
  double dfs(std::int32_t u, std::int32_t t, double f) {
    if (u == t) return f;
    for (auto& i = it[static_cast<std::size_t>(u)]; i < static_cast<std::int32_t>(adj[static_cast<std::size_t>(u)].size()); ++i) {
      DinicEdge& e = adj[static_cast<std::size_t>(u)][static_cast<std::size_t>(i)];
      if (e.cap < kMinCap) continue;
      if (level[static_cast<std::size_t>(e.to)] != level[static_cast<std::size_t>(u)] + 1) continue;
      double d = dfs(e.to, t, std::min(f, e.cap));
      if (d >= kMinFlow) {
        e.cap -= d;
        DinicEdge& r = adj[static_cast<std::size_t>(e.to)][static_cast<std::size_t>(e.rev)];
        r.cap += d;
        return d;
      }
    }
    return 0.0;
  }
};

// Edges from parent u to child v grouped by (u,v) for proportional/EB logic.
// In the DAG, parent u is a predecessor of child v on a shortest path.
struct EdgeGroup {
  std::int32_t from; // child v (destination of grouped edges)
  std::int32_t to;   // parent u (source of grouped edges)
  std::vector<EdgeId> eids; // underlying forward edges u->v (may be multiple parallel edges)
  Cap sum_cap {0.0};  // sum of residual capacities for Proportional placement
  Cap min_cap {0.0};  // min residual capacity for EqualBalanced placement
};

// Build grouped edges by (parent u, child v) that can reach destination t,
// using the current residual snapshot.
static std::vector<EdgeGroup> build_groups_residual(const StrictMultiDiGraph& g,
                                                    const PredDAG& dag, NodeId t,
                                                    const std::vector<Cap>& residual) {
  std::vector<EdgeGroup> groups;
  const auto& offsets = dag.parent_offsets;
  const auto& parents = dag.parents;
  const auto& via     = dag.via_edges;
  const auto N = g.num_nodes();
  // Compute reachability: BFS backward from destination t to identify nodes on SPF DAG.
  std::vector<char> reach(static_cast<std::size_t>(N), 0);
  if (t >= 0 && t < N) {
    std::queue<std::int32_t> q; q.push(t); reach[static_cast<std::size_t>(t)] = 1;
    while (!q.empty()) {
      auto v = q.front(); q.pop();
      // Iterate over v's predecessors (parents in the DAG).
      std::size_t s = static_cast<std::size_t>(offsets[static_cast<std::size_t>(v)]);
      std::size_t e = static_cast<std::size_t>(offsets[static_cast<std::size_t>(v + 1)]);
      for (std::size_t i = s; i < e; ++i) {
        auto u = parents[i];
        if (!reach[static_cast<std::size_t>(u)]) { reach[static_cast<std::size_t>(u)] = 1; q.push(u); }
      }
    }
  }
  // For each reachable node v, group its incoming DAG edges by parent node u.
  // This creates one group per (u, v) pair, aggregating parallel edges.
  for (std::int32_t v = 0; v < N; ++v) {
    if (!reach[static_cast<std::size_t>(v)]) continue;
    // Small linear grouping by parent (faster than a hash for typical degrees).
    std::vector<std::pair<std::int32_t, std::vector<std::int32_t>>> by_parent;
    const std::size_t s = static_cast<std::size_t>(offsets[static_cast<std::size_t>(v)]);
    const std::size_t e = static_cast<std::size_t>(offsets[static_cast<std::size_t>(v + 1)]);
    for (std::size_t i = s; i < e; ++i) {
      const auto u = parents[i];
      bool found = false;
      for (auto& pr : by_parent) {
        if (pr.first == u) { pr.second.push_back(via[i]); found = true; break; }
      }
      if (!found) by_parent.emplace_back(u, std::vector<std::int32_t>{ via[i] });
    }
    for (auto& kv : by_parent) {
      EdgeGroup gr; gr.from = v; gr.to = kv.first; gr.eids.clear();
      gr.sum_cap = static_cast<Cap>(0.0); gr.min_cap = std::numeric_limits<Cap>::infinity();
      for (auto eid0 : kv.second) {
        const Cap c = residual[static_cast<std::size_t>(eid0)];
        if (c >= kMinCap) {
          gr.eids.push_back(eid0);
          gr.sum_cap += c;
          gr.min_cap = std::min(gr.min_cap, c);
        }
      }
      if (gr.min_cap == std::numeric_limits<Cap>::infinity()) gr.min_cap = static_cast<Cap>(0.0);
      if (!gr.eids.empty()) groups.push_back(std::move(gr));
    }
  }
  return groups;
}

// Construct reversed residual graph for Dinic BFS/DFS using group capacities.
static void build_reversed_residual(FlowWorkspace& ws, std::int32_t N, const std::vector<EdgeGroup>& groups) {
  ws.reset(N);
  for (std::int32_t idx = 0; idx < static_cast<std::int32_t>(groups.size()); ++idx) {
    const auto& gr = groups[static_cast<std::size_t>(idx)];
    if (gr.sum_cap >= kMinCap) {
      ws.add_edge(gr.from, gr.to, gr.sum_cap, idx);
    }
  }
}

// Note: equal-balanced logic is implemented inline in place_on_dag to avoid
// maintaining a separate assignment routine. This keeps behavior explicit and
// consistent with the documented reference behavior.

} // namespace

FlowState::FlowState(const StrictMultiDiGraph& g) : g_(&g) {
  residual_.assign(static_cast<std::size_t>(g.num_edges()), 0.0);
  edge_flow_.assign(static_cast<std::size_t>(g.num_edges()), 0.0);
  // initialize residual to capacity
  auto cap = g.capacity_view();
  for (std::size_t i=0;i<residual_.size();++i) residual_[i] = cap[i];
}

FlowState::FlowState(const StrictMultiDiGraph& g, std::span<const Cap> residual_init) : g_(&g) {
  if (static_cast<std::size_t>(g.num_edges()) != residual_init.size()) {
    throw std::invalid_argument("FlowState: residual_init length must equal num_edges");
  }
  residual_.assign(residual_init.begin(), residual_init.end());
  edge_flow_.assign(static_cast<std::size_t>(g.num_edges()), 0.0);
}

void FlowState::reset() noexcept {
  auto cap = g_->capacity_view();
  for (std::size_t i=0;i<residual_.size();++i) residual_[i] = cap[i];
  std::fill(edge_flow_.begin(), edge_flow_.end(), 0.0);
}

void FlowState::reset(std::span<const Cap> residual_init) {
  if (residual_init.size() != residual_.size()) {
    throw std::invalid_argument("FlowState::reset: residual_init length must equal num_edges");
  }
  std::copy(residual_init.begin(), residual_init.end(), residual_.begin());
  std::fill(edge_flow_.begin(), edge_flow_.end(), 0.0);
}

Flow FlowState::place_on_dag(NodeId src, NodeId dst, const PredDAG& dag,
                             Flow requested_flow, FlowPlacement placement,
                             std::vector<std::pair<EdgeId, Flow>>* trace) {
  NGRAPH_PROFILE_SCOPE("place_on_dag");
  const auto N = g_->num_nodes();
  if (src < 0 || src >= N || dst < 0 || dst >= N || src == dst) return 0.0;

  // Build groups using current residual
  auto groups = build_groups_residual(*g_, dag, dst, residual_);

  Flow placed = static_cast<Flow>(0.0);
  double remaining = static_cast<double>(requested_flow);

  if (placement == FlowPlacement::Proportional) {
    // Proportional placement: use Dinic-like augmentation on reversed DAG.
    // We reverse the DAG (dst -> src) so flow propagates topologically from dst.
    FlowWorkspace ws; build_reversed_residual(ws, N, groups);
    while (remaining > kMinFlow && ws.bfs(dst, src)) {
      std::fill(ws.it.begin(), ws.it.end(), 0);
      Flow pushed_layer = static_cast<Flow>(0.0);
      while (true) {
        double pushed = ws.dfs(dst, src, remaining);
        if (pushed < kMinFlow) break;
        pushed_layer += static_cast<Flow>(pushed);
        remaining -= pushed;
        if (remaining <= kMinFlow) break;
      }
      if (pushed_layer < kMinFlow) break;
      placed += pushed_layer;

      // Distribute the pushed flow onto the underlying parallel edges.
      // For each group, divide the flow proportionally to residual capacity.
      // This ensures fair load balancing across parallel edges.
      for (std::size_t u = 0; u < ws.adj.size(); ++u) {
        for (const auto& e : ws.adj[u]) {
          if (e.group < 0) continue;  // skip reverse edges
          double sent = e.init_cap - e.cap; if (sent < kMinFlow) continue;
          const auto& gr = groups[static_cast<std::size_t>(e.group)];
          // Guard against division by zero when group capacity is numerically zero.
          double denom = gr.sum_cap > kMinCap ? gr.sum_cap : 1.0;
          // Proportional split: each edge gets share = sent * (edge_residual / sum_residual).
          for (auto eid : gr.eids) {
            Cap base = residual_[static_cast<std::size_t>(eid)];
            double share = sent * (static_cast<double>(base) / denom);
            edge_flow_[static_cast<std::size_t>(eid)] += static_cast<Cap>(share);
            residual_[static_cast<std::size_t>(eid)] = std::max(static_cast<Cap>(0.0), static_cast<Cap>(base - share));
            if (trace) { trace->emplace_back(eid, static_cast<Flow>(share)); }
          }
        }
      }
      // Rebuild groups for next tier using updated residual
      groups = build_groups_residual(*g_, dag, dst, residual_);
      build_reversed_residual(ws, N, groups);
    }
  } else {
    // EqualBalanced placement: split flow equally across parallel edges, with
    // topological accumulation to correctly handle reconvergent DAGs.

    // Build forward adjacency from parent u to child v for each group and
    // compute aggregated reverse capacities per group.
    std::vector<std::vector<std::size_t>> succ(static_cast<std::size_t>(N));
    std::vector<double> rev_cap(groups.size(), 0.0);
    for (std::size_t gi = 0; gi < groups.size(); ++gi) {
      const auto& gr = groups[gi];
      if (gr.eids.empty()) continue;
      // EB: group admissible total = min_edge_residual * |edges|
      const double cap_rev = static_cast<double>(gr.min_cap) * static_cast<double>(gr.eids.size());
      if (cap_rev >= kMinCap) {
        succ[static_cast<std::size_t>(gr.to)].push_back(gi); // u -> v (group index)
        rev_cap[gi] = cap_rev;
      }
    }

    // Compute reachability from src on this succ graph (to ignore disconnected parts).
    std::vector<char> reach(static_cast<std::size_t>(N), 0);
    if (src >= 0 && src < N) {
      std::queue<std::int32_t> q; q.push(src); reach[static_cast<std::size_t>(src)] = 1;
      while (!q.empty()) {
        auto u = q.front(); q.pop();
        for (auto gi : succ[static_cast<std::size_t>(u)]) {
          auto v = groups[gi].from;
          if (!reach[static_cast<std::size_t>(v)]) { reach[static_cast<std::size_t>(v)] = 1; q.push(v); }
        }
      }
    }

    // Precompute total outgoing per-edge count per node for equal split.
    std::vector<int> node_split(static_cast<std::size_t>(N), 0);
    for (std::size_t u = 0; u < succ.size(); ++u) {
      if (!reach[u]) continue;
      int s = 0;
      for (auto gi : succ[u]) s += static_cast<int>(groups[gi].eids.size());
      node_split[u] = s;
    }

    // Topological accumulation: indegree counts over reachable subgraph.
    std::vector<int> indeg(static_cast<std::size_t>(N), 0);
    for (std::size_t u = 0; u < succ.size(); ++u) {
      if (!reach[u]) continue;
      for (auto gi : succ[u]) {
        auto v = static_cast<std::size_t>(groups[gi].from);
        if (reach[v]) indeg[v] += 1;
      }
    }

    // Kahn's algorithm over reachable nodes starting from src.
    std::queue<std::int32_t> q;
    std::vector<double> inflow(static_cast<std::size_t>(N), 0.0);
    std::vector<double> assigned(groups.size(), 0.0);
    if (src >= 0 && src < N && reach[static_cast<std::size_t>(src)]) {
      q.push(src);
      inflow[static_cast<std::size_t>(src)] = 1.0; // unit flow
    }

    while (!q.empty()) {
      auto u = q.front(); q.pop();
      double f_in = inflow[static_cast<std::size_t>(u)];
      if (f_in < kEpsilon) continue;
      int split = node_split[static_cast<std::size_t>(u)];
      if (split <= 0) continue;
      for (auto gi : succ[static_cast<std::size_t>(u)]) {
        const auto& gr = groups[gi];
        if (gr.eids.empty()) continue;
        // Group share proportional to number of edges (equal per-edge split).
        double push = f_in * (static_cast<double>(gr.eids.size()) / static_cast<double>(split));
        if (push < kEpsilon) continue;
        assigned[gi] += push;
        auto v = static_cast<std::size_t>(gr.from);
        inflow[v] += push;
        // Decrement indegree and enqueue child when all parents processed.
        if (--indeg[v] == 0) q.push(static_cast<std::int32_t>(v));
      }
    }

    // Early exit if destination is not reachable with capacity.
    // With require_capacity=false, SPF includes zero-capacity edges in the DAG,
    // but those edges are filtered out during group building. If no path with
    // capacity reaches the destination, return 0.
    if (inflow[static_cast<std::size_t>(dst)] < kEpsilon) {
      return static_cast<Flow>(0.0);
    }

    // Single-pass ECMP admission: scale the unit assignment by the smallest
    // per-group headroom so that no edge is oversubscribed under *fixed equal
    // per-edge splits*. Any further injection with the same splits would violate
    // the saturated group(s).
    double ratio = std::numeric_limits<double>::infinity();
    for (std::size_t gi = 0; gi < groups.size(); ++gi) {
      if (assigned[gi] > 0.0) {
        double r = rev_cap[gi] / assigned[gi];
        if (r < ratio) ratio = r;
      }
    }
    if (!std::isfinite(ratio)) ratio = 0.0;

    // Place up to requested amount.
    Flow use = static_cast<Flow>(std::min(ratio, static_cast<double>(remaining)));
    if (use >= kMinFlow) {
      placed += use;
      // Apply scaled group assignments equally over parallel edges.
      for (std::size_t gi = 0; gi < groups.size(); ++gi) {
        const auto& gr = groups[gi]; if (gr.eids.empty()) continue;
        double flow_scaled = assigned[gi] * static_cast<double>(use);
        if (flow_scaled < kMinFlow) continue;
        double per_edge = flow_scaled / static_cast<double>(gr.eids.size());
        for (auto eid : gr.eids) {
          edge_flow_[static_cast<std::size_t>(eid)] += static_cast<Flow>(per_edge);
          double base = static_cast<double>(residual_[static_cast<std::size_t>(eid)]);
          residual_[static_cast<std::size_t>(eid)] = static_cast<Cap>(std::max(0.0, base - per_edge));
          if (trace) trace->emplace_back(eid, static_cast<Flow>(per_edge));
        }
      }
    }
  }
  return placed;
}

Flow FlowState::place_max_flow(NodeId src, NodeId dst, FlowPlacement placement,
                               bool shortest_path, bool require_capacity,
                               std::span<const bool> node_mask,
                               std::span<const bool> edge_mask) {
  // require_capacity controls routing behavior:
  //   - true: Require edges to have capacity, exclude saturated links (SDN/TE, progressive)
  //   - false: Routes based on costs only, ignore capacity (IP/IGP, fixed routing)
  //
  // NOTE: With FlowPlacement::EqualBalanced + require_capacity=true, this behaves as a
  // progressive "fill". For IP ECMP, use require_capacity=false + shortest_path=true.

  // Validate mask lengths
  if (!node_mask.empty() && node_mask.size() != static_cast<std::size_t>(g_->num_nodes())) {
    throw std::invalid_argument("FlowState::place_max_flow: node_mask length mismatch");
  }
  if (!edge_mask.empty() && edge_mask.size() != static_cast<std::size_t>(g_->num_edges())) {
    throw std::invalid_argument("FlowState::place_max_flow: edge_mask length mismatch");
  }

  Flow total = static_cast<Flow>(0.0);
  while (true) {
    EdgeSelection sel;
    sel.multi_edge = true;
    sel.require_capacity = require_capacity;
    sel.tie_break = EdgeTieBreak::Deterministic;
    auto [dist, dag] = shortest_paths(*g_, src, dst, /*multipath=*/true, sel,
                                      require_capacity ? residual_ : std::span<const Cap>{},
                                      node_mask, edge_mask);
    if (static_cast<std::size_t>(dst) >= dag.parent_offsets.size()-1 || dag.parent_offsets[static_cast<std::size_t>(dst)] == dag.parent_offsets[static_cast<std::size_t>(dst)+1]) {
      break;
    }
    Flow placed = place_on_dag(src, dst, dag, std::numeric_limits<double>::infinity(), placement);
    if (placed < kMinFlow) break;
    total += placed;
    if (shortest_path) break;
  }
  return total;
}

MinCut FlowState::compute_min_cut(NodeId src, std::span<const bool> node_mask, std::span<const bool> edge_mask) const {
  MinCut out;
  const auto N = g_->num_nodes();

  // Early return if source is out of range
  if (!(src >= 0 && src < N)) return out;

  // Validate mask lengths
  if (!node_mask.empty() && node_mask.size() != static_cast<std::size_t>(N)) {
    throw std::invalid_argument("FlowState::compute_min_cut: node_mask length mismatch");
  }
  if (!edge_mask.empty() && edge_mask.size() != static_cast<std::size_t>(g_->num_edges())) {
    throw std::invalid_argument("FlowState::compute_min_cut: edge_mask length mismatch");
  }

  const auto row = g_->row_offsets_view();
  const auto col = g_->col_indices_view();
  const auto aei = g_->adj_edge_index_view();
  const auto in_row = g_->in_row_offsets_view();
  const auto in_col = g_->in_col_indices_view();
  const auto in_aei = g_->in_adj_edge_index_view();
  const bool use_node_mask = !node_mask.empty();
  const bool use_edge_mask = !edge_mask.empty();

  // Early return if source is masked out
  if (use_node_mask && !node_mask[static_cast<std::size_t>(src)]) return out;

  std::vector<char> visited(static_cast<std::size_t>(N), 0);
  std::queue<std::int32_t> q;
  visited[static_cast<std::size_t>(src)] = 1;
  q.push(src);
  while (!q.empty()) {
    auto u = q.front(); q.pop();
    if (use_node_mask && !node_mask[static_cast<std::size_t>(u)]) continue;
    // Forward residual arcs
    auto start = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto end   = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
    for (std::size_t j = start; j < end; ++j) {
      auto v = static_cast<std::int32_t>(col[j]);
      auto eid = static_cast<std::size_t>(aei[j]);
      if (use_edge_mask && !edge_mask[eid]) continue;
      if (use_node_mask && !node_mask[static_cast<std::size_t>(v)]) continue;
      if (residual_[eid] > kMinCap && !visited[static_cast<std::size_t>(v)]) {
        visited[static_cast<std::size_t>(v)] = 1;
        q.push(v);
      }
    }
    // Reverse residual arcs
    auto rs = static_cast<std::size_t>(in_row[static_cast<std::size_t>(u)]);
    auto re = static_cast<std::size_t>(in_row[static_cast<std::size_t>(u)+1]);
    for (std::size_t j = rs; j < re; ++j) {
      auto w = static_cast<std::int32_t>(in_col[j]);
      auto eid = static_cast<std::size_t>(in_aei[j]);
      if (use_edge_mask && !edge_mask[eid]) continue;
      if (use_node_mask && !node_mask[static_cast<std::size_t>(w)]) continue;
      double flow_e = g_->capacity_view()[eid] - residual_[eid];
      if (flow_e > kMinFlow && !visited[static_cast<std::size_t>(w)]) {
        visited[static_cast<std::size_t>(w)] = 1;
        q.push(w);
      }
    }
  }
  // Collect cut edges
  for (std::int32_t u = 0; u < N; ++u) {
    if (!visited[static_cast<std::size_t>(u)]) continue;
    if (use_node_mask && !node_mask[static_cast<std::size_t>(u)]) continue;
    auto s3 = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto e3 = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
    for (std::size_t j = s3; j < e3; ++j) {
      auto v = static_cast<std::int32_t>(col[j]);
      if (use_node_mask && !node_mask[static_cast<std::size_t>(v)]) continue;
      if (visited[static_cast<std::size_t>(v)]) continue;
      auto eid = static_cast<std::size_t>(aei[j]);
      if (use_edge_mask && !edge_mask[eid]) continue;
      if (residual_[eid] <= kMinCap) {
        out.edges.push_back(static_cast<EdgeId>(eid));
      }
    }
  }
  return out;
}

void FlowState::apply_deltas(std::span<const std::pair<EdgeId, Flow>> deltas, bool add) noexcept {
  const auto cap = g_->capacity_view();
  for (const auto& pr : deltas) {
    std::size_t eid = static_cast<std::size_t>(pr.first);
    double df = static_cast<double>(pr.second);
    if (df <= 0.0) continue;
    if (eid >= residual_.size()) continue;
    // base_res is not needed; use base_flow and cap directly
    double base_flow = edge_flow_[eid];
    if (add) {
      edge_flow_[eid] = static_cast<Flow>(base_flow + df);
      double new_res = std::max(0.0, static_cast<double>(cap[eid]) - edge_flow_[eid]);
      residual_[eid] = static_cast<Cap>(new_res);
    } else {
      // Clamp flow into [0, cap] and set residual accordingly.
      double unclamped = base_flow - df;
      double clamped = std::min(static_cast<double>(cap[eid]), std::max(0.0, unclamped));
      edge_flow_[eid] = static_cast<Flow>(clamped);
      residual_[eid] = static_cast<Cap>(std::max(0.0, static_cast<double>(cap[eid]) - clamped));
    }
  }
}

} // namespace netgraph::core
