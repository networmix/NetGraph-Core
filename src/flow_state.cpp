/*
  FlowState — residual capacities and placement over a fixed graph.

  Maintains per-edge residual capacity and cumulative edge flows. Supports
  two placement strategies when pushing flow along an SPF DAG:
    - Proportional: distribute flow proportionally to residual capacity,
      processing nodes in topological order from source to destination.
    - EqualBalanced: distribute flow equally across available parallel edges,
      respecting per-edge residual capacity constraints.
*/
#include "netgraph/core/flow_state.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/constants.hpp"

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
  double init_cap;        // initial capacity (snapshot before augmentation)
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
static std::vector<EdgeGroup> build_groups_residual(const StrictMultiDiGraph& g, const PredDAG& dag, NodeId t,
                                                    const std::vector<Cap>& residual) {
  std::vector<EdgeGroup> groups;
  auto offsets = dag.parent_offsets;
  auto parents = dag.parents;
  auto via = dag.via_edges;
  const auto N = g.num_nodes();
  // Compute reachability: BFS backward from destination t to identify nodes on SPF DAG.
  // Python developers: this is like a reverse BFS to find ancestors.
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
    // Group edges by parent node (unordered_map is like Python dict).
    std::unordered_map<std::int32_t, std::vector<std::int32_t>> by_parent;
    std::size_t s = static_cast<std::size_t>(offsets[static_cast<std::size_t>(v)]);
    std::size_t e = static_cast<std::size_t>(offsets[static_cast<std::size_t>(v + 1)]);
    for (std::size_t i = s; i < e; ++i) {
      auto u = parents[i];
      by_parent[u].push_back(via[i]);
    }
    for (auto& kv : by_parent) {
      EdgeGroup gr; gr.from = v; gr.to = kv.first; gr.eids.clear();
      gr.sum_cap = static_cast<Cap>(0.0); gr.min_cap = std::numeric_limits<Cap>::infinity();
      // Include only edges with residual >= kMinCap to match residual filtering.
      for (auto eid0 : kv.second) {
        Cap c = residual[static_cast<std::size_t>(eid0)];
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
                             bool shortest_path,
                             std::vector<std::pair<EdgeId, Flow>>* trace) {
  if (src == dst) return 0.0;
  const auto N = g_->num_nodes();
  if (dst < 0 || dst >= N) return 0.0;

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
      // Always cap per-augment push by remaining requested amount
      auto cap_push = [&](){ return shortest_path ? remaining : remaining; };
      Flow pushed_layer = static_cast<Flow>(0.0);
      while (true) {
        double pushed = ws.dfs(dst, src, cap_push());
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
            if (trace && share >= kMinFlow) {
              trace->emplace_back(eid, static_cast<Flow>(share));
            }
          }
        }
      }
      // Rebuild groups for next tier using updated residual
      groups = build_groups_residual(*g_, dag, dst, residual_);
      build_reversed_residual(ws, N, groups);
    }
  } else {
    // EqualBalanced placement: split flow equally across parallel edges.
    // This ensures uniform load distribution, limited by the min residual capacity.

    // Build reversed adjacency succ: parent u <- child v for each group.
    std::vector<std::vector<std::size_t>> succ(static_cast<std::size_t>(N));
    std::vector<double> rev_cap(groups.size(), 0.0);
    for (std::size_t gi = 0; gi < groups.size(); ++gi) {
      const auto& gr = groups[gi];
      if (gr.eids.empty()) continue;
      // Find minimum residual capacity among parallel edges.
      double min_c = std::numeric_limits<double>::infinity();
      for (auto eid : gr.eids) {
        double c = static_cast<double>(residual_[static_cast<std::size_t>(eid)]);
        if (c < min_c) min_c = c;
      }
      if (!std::isfinite(min_c)) min_c = 0.0;
      // Aggregate capacity for this group is min_cap * num_edges.
      double cap_rev = min_c * static_cast<double>(gr.eids.size());
      if (cap_rev >= kMinCap) {
        succ[static_cast<std::size_t>(gr.to)].push_back(gi);
        rev_cap[gi] = cap_rev;
      }
    }
    // Equal-split BFS from src over reversed graph.
    // We compute a relative assignment for each group based on equal splitting.
    std::vector<double> assigned(groups.size(), 0.0);

    // Precompute total edge count per node (for equal splitting).
    std::vector<int> node_split(static_cast<std::size_t>(N), 0);
    for (std::size_t u = 0; u < succ.size(); ++u) {
      int s = 0; for (auto gi : succ[u]) s += static_cast<int>(groups[gi].eids.size());
      node_split[u] = s;
    }

    // BFS from src with a flow value of 1.0 (unit flow).
    // Each node splits its inflow equally among all outgoing edges.
    std::vector<char> visited(static_cast<std::size_t>(N), 0);
    std::queue<std::pair<std::int32_t,double>> q2; q2.emplace(src, 1.0);
    while (!q2.empty()) {
      auto [u, inflow] = q2.front(); q2.pop();
      if (inflow < kMinFlow) continue;
      visited[static_cast<std::size_t>(u)] = 1;
      int split = node_split[static_cast<std::size_t>(u)];
      if (split <= 0) continue;
      for (auto gi : succ[static_cast<std::size_t>(u)]) {
        const auto& gr = groups[gi]; if (gr.eids.empty()) continue;
        // Each parallel edge in this group gets an equal share.
        double push = inflow * (static_cast<double>(gr.eids.size()) / static_cast<double>(split));
        if (push < kMinFlow) continue;
        assigned[gi] += push;
        if (!visited[static_cast<std::size_t>(gr.from)]) q2.emplace(gr.from, push);
      }
    }
    // Compute the bottleneck scaling ratio: max flow we can push without exceeding capacity.
    // For each group, the ratio is (available_capacity / assigned_relative_flow).
    double ratio = std::numeric_limits<double>::infinity();
    for (std::size_t gi = 0; gi < groups.size(); ++gi) {
      if (assigned[gi] > 0.0) {
        double r = rev_cap[gi] / assigned[gi];
        if (r < ratio) ratio = r;
      }
    }
    if (!std::isfinite(ratio)) ratio = 0.0;

    // Scale the unit flow by the bottleneck ratio to get the actual flow to place.
    // In shortest_path mode, cap by requested amount.
    Flow use = static_cast<Flow>(std::min(ratio, static_cast<double>(remaining)));
    if (use >= kMinFlow) {
      placed += use;
      // Apply scaled flow to each group, distributing equally among parallel edges.
      for (std::size_t gi = 0; gi < groups.size(); ++gi) {
        const auto& gr = groups[gi]; if (gr.eids.empty()) continue;
        double flow_scaled = assigned[gi] * static_cast<double>(use);
        if (flow_scaled < kMinFlow) continue;
        // Equal split: each edge in the group gets the same flow.
        double per_edge = flow_scaled / static_cast<double>(gr.eids.size());
        for (auto eid : gr.eids) {
          edge_flow_[static_cast<std::size_t>(eid)] += static_cast<Flow>(per_edge);
          double base = static_cast<double>(residual_[static_cast<std::size_t>(eid)]);
          residual_[static_cast<std::size_t>(eid)] = static_cast<Cap>(std::max(0.0, base - per_edge));
          if (trace && per_edge >= kMinFlow) trace->emplace_back(eid, static_cast<Flow>(per_edge));
        }
      }
    }
  }
  return placed;
}

Flow FlowState::place_max_flow(NodeId src, NodeId dst, FlowPlacement placement, bool shortest_path) {
  Flow total = static_cast<Flow>(0.0);
  while (true) {
    EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
    auto [dist, dag] = shortest_paths(*g_, src, dst, /*multipath=*/true, sel, residual_);
    if (static_cast<std::size_t>(dst) >= dag.parent_offsets.size()-1 || dag.parent_offsets[static_cast<std::size_t>(dst)] == dag.parent_offsets[static_cast<std::size_t>(dst)+1]) {
      break;
    }
    Flow placed = place_on_dag(src, dst, dag, std::numeric_limits<double>::infinity(), placement, shortest_path);
    if (placed < kMinFlow) break;
    total += placed;
    if (shortest_path) break;
  }
  return total;
}

MinCut FlowState::compute_min_cut(NodeId src, std::span<const bool> node_mask, std::span<const bool> edge_mask) const {
  MinCut out;
  const auto N = g_->num_nodes();
  const auto row = g_->row_offsets_view();
  const auto col = g_->col_indices_view();
  const auto aei = g_->adj_edge_index_view();
  const auto in_row = g_->in_row_offsets_view();
  const auto in_col = g_->in_col_indices_view();
  const auto in_aei = g_->in_adj_edge_index_view();
  std::vector<char> visited(static_cast<std::size_t>(N), 0);
  std::queue<std::int32_t> q;
  if (src >= 0 && src < N) { visited[static_cast<std::size_t>(src)] = 1; q.push(src); }
  const bool use_node_mask = (node_mask.size() == static_cast<std::size_t>(N));
  const bool use_edge_mask = (edge_mask.size() == static_cast<std::size_t>(g_->num_edges()));
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
      double new_flow = std::max(0.0, base_flow - df);
      edge_flow_[eid] = static_cast<Flow>(new_flow);
      double new_res = std::min(static_cast<double>(cap[eid]), static_cast<double>(cap[eid]) - new_flow);
      // guard for numeric drift
      if (new_res < kMinCap && (static_cast<double>(cap[eid]) - new_flow) < kMinCap) new_res = std::max(0.0, static_cast<double>(cap[eid]) - new_flow);
      residual_[eid] = static_cast<Cap>(new_res);
    }
  }
}

} // namespace netgraph::core
