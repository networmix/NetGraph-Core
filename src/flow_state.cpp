/*
  FlowState — residual capacities and placement over a fixed graph.

  This module maintains per-edge residual capacity and cumulative edge flows.
  It supports two placement strategies when pushing flow along an SPF DAG:
    - Proportional: distribute according to current residuals in each (u,v)
      neighbor group, tier-by-tier using a Dinic-like layered approach.
    - EqualBalanced: approximate equal splitting across parallel edges based on
      a reversed residual graph that models per-edge min capacity.

  The goal is deterministic behavior matching the documented semantics while
  keeping internals explicit and efficient in C++.
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

struct DinicEdge {
  std::int32_t to;
  std::int32_t rev;
  double cap;
  double init_cap;
  std::int32_t group; // index into groups (-1 if none)
};

// Workspace for Dinic-like augmentations on the reversed residual graph.
struct FlowWorkspace {
  std::vector<std::vector<DinicEdge>> adj; // reversed residual graph
  std::vector<std::int32_t> level;
  std::vector<std::int32_t> it;

  void reset(std::int32_t n) {
    adj.assign(static_cast<std::size_t>(n), {});
    level.assign(static_cast<std::size_t>(n), -1);
    it.assign(static_cast<std::size_t>(n), 0);
  }
  void add_edge(std::int32_t u, std::int32_t v, double c, std::int32_t group_idx) {
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
struct EdgeGroup {
  std::int32_t from; // child v
  std::int32_t to;   // parent u
  std::vector<EdgeId> eids; // underlying forward edges u->v
  Cap sum_cap {0.0};
  Cap min_cap {0.0};
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
  // reachability to t
  std::vector<char> reach(static_cast<std::size_t>(N), 0);
  if (t >= 0 && t < N) {
    std::queue<std::int32_t> q; q.push(t); reach[static_cast<std::size_t>(t)] = 1;
    while (!q.empty()) {
      auto v = q.front(); q.pop();
      std::size_t s = static_cast<std::size_t>(offsets[static_cast<std::size_t>(v)]);
      std::size_t e = static_cast<std::size_t>(offsets[static_cast<std::size_t>(v + 1)]);
      for (std::size_t i = s; i < e; ++i) {
        auto u = parents[i];
        if (!reach[static_cast<std::size_t>(u)]) { reach[static_cast<std::size_t>(u)] = 1; q.push(u); }
      }
    }
  }
  for (std::int32_t v = 0; v < N; ++v) {
    if (!reach[static_cast<std::size_t>(v)]) continue;
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
      // include only edges with residual >= kMinCap to match residual filtering
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
        if (shortest_path) break;
        if (remaining <= kMinFlow) break;
      }
      if (pushed_layer < kMinFlow) break;
      placed += pushed_layer;
      // Distribute on forward parallels proportional to residual snapshot before update
      for (std::size_t u = 0; u < ws.adj.size(); ++u) {
        for (const auto& e : ws.adj[u]) {
          if (e.group < 0) continue;
          double sent = e.init_cap - e.cap; if (sent < kMinFlow) continue;
          const auto& gr = groups[static_cast<std::size_t>(e.group)];
          double denom = gr.sum_cap > 0.0 ? gr.sum_cap : 1.0;
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
      if (shortest_path) break;
      // Rebuild groups for next tier using updated residual
      groups = build_groups_residual(*g_, dag, dst, residual_);
      build_reversed_residual(ws, N, groups);
    }
  } else { // EqualBalanced (reversed residual)
    // Build reversed adjacency succ: parent u <- child v for each group
    std::vector<std::vector<std::size_t>> succ(static_cast<std::size_t>(N));
    std::vector<double> rev_cap(groups.size(), 0.0);
    for (std::size_t gi = 0; gi < groups.size(); ++gi) {
      const auto& gr = groups[gi];
      if (gr.eids.empty()) continue;
      double min_c = std::numeric_limits<double>::infinity();
      for (auto eid : gr.eids) {
        double c = static_cast<double>(residual_[static_cast<std::size_t>(eid)]);
        if (c < min_c) min_c = c;
      }
      if (!std::isfinite(min_c)) min_c = 0.0;
      double cap_rev = min_c * static_cast<double>(gr.eids.size());
      if (cap_rev >= kMinCap) {
        succ[static_cast<std::size_t>(gr.to)].push_back(gi);
        rev_cap[gi] = cap_rev;
      }
    }
    // Equal-split BFS from src over reversed graph
    std::vector<double> assigned(groups.size(), 0.0);
    std::vector<int> node_split(static_cast<std::size_t>(N), 0);
    for (std::size_t u = 0; u < succ.size(); ++u) {
      int s = 0; for (auto gi : succ[u]) s += static_cast<int>(groups[gi].eids.size());
      node_split[u] = s;
    }
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
        double push = inflow * (static_cast<double>(gr.eids.size()) / static_cast<double>(split));
        if (push < kMinFlow) continue;
        assigned[gi] += push;
        if (!visited[static_cast<std::size_t>(gr.from)]) q2.emplace(gr.from, push);
      }
    }
    // Compute scaling ratio
    double ratio = std::numeric_limits<double>::infinity();
    for (std::size_t gi = 0; gi < groups.size(); ++gi) {
      if (assigned[gi] > 0.0) {
        double r = rev_cap[gi] / assigned[gi];
        if (r < ratio) ratio = r;
      }
    }
    if (!std::isfinite(ratio)) ratio = 0.0;
    // In shortest_path mode, cap the equal-balanced push by the requested amount;
    // otherwise allow multiple tiers per iteration.
    Flow use = static_cast<Flow>(std::min(ratio, static_cast<double>(remaining)));
    if (use >= kMinFlow) {
      placed += use;
      for (std::size_t gi = 0; gi < groups.size(); ++gi) {
        const auto& gr = groups[gi]; if (gr.eids.empty()) continue;
        double flow_scaled = assigned[gi] * static_cast<double>(use);
        if (flow_scaled < kMinFlow) continue;
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

MinCut FlowState::compute_min_cut(NodeId src, const bool* node_mask, const bool* edge_mask) const {
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
  while (!q.empty()) {
    auto u = q.front(); q.pop();
    if (node_mask && !node_mask[static_cast<std::size_t>(u)]) continue;
    // Forward residual arcs
    auto start = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto end   = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
    for (std::size_t j = start; j < end; ++j) {
      auto v = static_cast<std::int32_t>(col[j]);
      auto eid = static_cast<std::size_t>(aei[j]);
      if (edge_mask && !edge_mask[eid]) continue;
      if (node_mask && !node_mask[static_cast<std::size_t>(v)]) continue;
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
      if (edge_mask && !edge_mask[eid]) continue;
      if (node_mask && !node_mask[static_cast<std::size_t>(w)]) continue;
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
    if (node_mask && !node_mask[static_cast<std::size_t>(u)]) continue;
    auto s3 = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto e3 = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
    for (std::size_t j = s3; j < e3; ++j) {
      auto v = static_cast<std::int32_t>(col[j]);
      if (node_mask && !node_mask[static_cast<std::size_t>(v)]) continue;
      if (visited[static_cast<std::size_t>(v)]) continue;
      auto eid = static_cast<std::size_t>(aei[j]);
      if (edge_mask && !edge_mask[eid]) continue;
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
