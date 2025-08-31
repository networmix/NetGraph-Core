#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/shortest_paths.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

namespace netgraph::core {

namespace {
constexpr double kMinCap = 1.0 / 4096.0;
constexpr double kMinFlow = 1.0 / 4096.0;

struct DinicEdge {
  std::int32_t to;
  std::int32_t rev;   // index into adj[to]
  double cap;         // residual capacity
  double init_cap;    // initial capacity (for flow extraction)
  std::int32_t group; // index into groups (-1 if none)
};

// Use public shortest_paths_with_residual from shortest_paths.cpp

struct FlowWorkspace {
  std::vector<std::vector<DinicEdge>> adj; // reversed residual graph for proportional mode
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
        // reduce forward
        e.cap -= d;
        // increase backward
        DinicEdge& r = adj[static_cast<std::size_t>(e.to)][static_cast<std::size_t>(e.rev)];
        r.cap += d;
        return d;
      }
    }
    return 0.0;
  }
};

struct EdgeGroup {
  std::int32_t from; // v (child)
  std::int32_t to;   // u (parent)
  std::vector<std::int32_t> eids; // underlying forward edges u->v
  double sum_cap {0.0};
  double min_cap {0.0};
};

// Build edge groups from PredDAG restricted to nodes that can reach t (forward),
// discovered by a reversed BFS over the DAG starting from t.
static std::vector<EdgeGroup> build_groups(const StrictMultiDiGraph& g, const PredDAG& dag, NodeId t) {
  std::vector<EdgeGroup> groups;
  auto offsets = dag.parent_offsets;
  auto parents = dag.parents;
  auto via = dag.via_edges;
  auto cap = g.capacity_view();
  const auto N = g.num_nodes();

  // Mark nodes reachable to t via DAG (following parents backward)
  std::vector<char> reach(static_cast<std::size_t>(N), 0);
  if (t >= 0 && t < N) {
    std::queue<std::int32_t> q;
    q.push(t);
    reach[static_cast<std::size_t>(t)] = 1;
    while (!q.empty()) {
      auto v = q.front(); q.pop();
      std::size_t s = static_cast<std::size_t>(offsets[static_cast<std::size_t>(v)]);
      std::size_t e = static_cast<std::size_t>(offsets[static_cast<std::size_t>(v + 1)]);
      for (std::size_t i = s; i < e; ++i) {
        auto u = parents[i];
        if (!reach[static_cast<std::size_t>(u)]) {
          reach[static_cast<std::size_t>(u)] = 1;
          q.push(u);
        }
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
      EdgeGroup gr;
      gr.from = v;
      gr.to = kv.first;
      gr.eids = std::move(kv.second);
      gr.sum_cap = 0.0;
      gr.min_cap = std::numeric_limits<double>::infinity();
      for (auto eid : gr.eids) {
        double c = cap[static_cast<std::size_t>(eid)];
        gr.sum_cap += c;
        gr.min_cap = std::min(gr.min_cap, c);
      }
      if (gr.min_cap == std::numeric_limits<double>::infinity()) gr.min_cap = 0.0;
      groups.push_back(std::move(gr));
    }
  }
  return groups;
}

// Build reversed residual network for Dinic from groups (proportional mode)
static void build_reversed_residual(FlowWorkspace& ws, std::int32_t N, const std::vector<EdgeGroup>& groups) {
  ws.reset(N);
  for (std::int32_t idx = 0; idx < static_cast<std::int32_t>(groups.size()); ++idx) {
    const auto& gr = groups[static_cast<std::size_t>(idx)];
    if (gr.sum_cap >= kMinCap) {
      // reversed orientation: edge from v (child) to u (parent)
      ws.add_edge(gr.from, gr.to, gr.sum_cap, idx);
    }
  }
}

struct ForwardAssign { std::int32_t u; std::int32_t v; double flow; std::size_t group_index; };

// Equal-balanced assignment on forward adjacency using groups
static std::pair<double, std::vector<ForwardAssign>> equal_balance_assign(std::int32_t N, std::int32_t src,
                                                                         const std::vector<EdgeGroup>& groups) {
  // Build forward adjacency from groups
  std::vector<std::vector<std::size_t>> adj(static_cast<std::size_t>(N));
  for (std::size_t gi = 0; gi < groups.size(); ++gi) {
    adj[static_cast<std::size_t>(groups[gi].to)].push_back(gi); // u -> (group u->v)
  }
  // Count total parallel edges leaving each node (sum over counts in groups)
  std::vector<int> split(static_cast<std::size_t>(N), 0);
  for (std::size_t u = 0; u < adj.size(); ++u) {
    int s = 0;
    for (auto gi : adj[u]) s += static_cast<int>(groups[gi].eids.size());
    split[u] = s;
  }
  // BFS-like pass from src with nominal 1.0
  std::vector<char> visited(static_cast<std::size_t>(N), 0);
  std::queue<std::pair<std::int32_t,double>> q;
  q.emplace(src, 1.0);
  std::vector<ForwardAssign> assigns;
  while (!q.empty()) {
    auto [u, inflow] = q.front(); q.pop();
    visited[static_cast<std::size_t>(u)] = 1;
    if (inflow < kMinFlow) continue;
    int sc = split[static_cast<std::size_t>(u)];
    if (sc <= 0) continue;
    for (auto gi : adj[static_cast<std::size_t>(u)]) {
      const auto& gr = groups[gi];
      if (gr.eids.empty()) continue;
      double push = inflow * static_cast<double>(gr.eids.size()) / static_cast<double>(sc);
      if (push < kMinFlow) continue;
      assigns.push_back(ForwardAssign{gr.to, gr.from, push, gi});
      if (!visited[static_cast<std::size_t>(gr.from)]) {
        q.emplace(gr.from, push);
      }
    }
  }
  // Compute scaling ratio based on capacities in forward direction
  double min_ratio = std::numeric_limits<double>::infinity();
  for (auto const& asg : assigns) {
    const auto& gr = groups[asg.group_index];
    // Effective capacity for equal-balanced is min(cap) * count
    double cap_eff = gr.min_cap * static_cast<double>(gr.eids.size());
    if (asg.flow > 0.0) {
      double ratio = cap_eff / asg.flow;
      if (ratio < min_ratio) min_ratio = ratio;
    }
  }
  if (min_ratio == std::numeric_limits<double>::infinity()) min_ratio = 0.0;
  return {min_ratio, assigns};
}

} // namespace

std::pair<double, FlowSummary>
calc_max_flow(const StrictMultiDiGraph& g, NodeId s, NodeId t,
              FlowPlacement placement, bool shortest_path,
              double /*eps*/, bool with_edge_flows,
              const bool* node_mask, const bool* edge_mask) {
  FlowSummary summary;
  const auto N = g.num_nodes();
  if (s < 0 || s >= N || t < 0 || t >= N || s == t) {
    return {0.0, std::move(summary)};
  }
  if ((node_mask && !node_mask[static_cast<std::size_t>(s)]) ||
      (node_mask && !node_mask[static_cast<std::size_t>(t)])) {
    return {0.0, std::move(summary)};
  }

  // Residual capacities per edge
  std::vector<double> residual(static_cast<std::size_t>(g.num_edges()));
  for (std::size_t i=0;i<residual.size();++i) residual[i] = g.capacity_view()[i];
  std::vector<double> edge_flows_accum;
  if (with_edge_flows) edge_flows_accum.assign(static_cast<std::size_t>(g.num_edges()), 0.0);
  double total = 0.0;
  std::unordered_map<double,double> cost_dist;

  while (true) {
    auto [dist, dag] = shortest_paths_with_residual(g, s, t, EdgeSelect::AllMinCostWithCapRemaining, /*multipath*/ true, /*eps*/ 1e-12, residual, node_mask, edge_mask);
    if (static_cast<std::size_t>(t) >= dag.parent_offsets.size()-1 || dag.parent_offsets[static_cast<std::size_t>(t)] == dag.parent_offsets[static_cast<std::size_t>(t)+1]) {
      break; // no path
    }
    double path_cost = dist[static_cast<std::size_t>(t)];
    auto groups = build_groups(g, dag, t);

    // Recompute group capacities from residuals
    for (auto& gr : groups) {
      gr.sum_cap = 0.0;
      gr.min_cap = std::numeric_limits<double>::infinity();
      for (auto eid : gr.eids) {
        double c = residual[static_cast<std::size_t>(eid)];
        gr.sum_cap += c;
        gr.min_cap = std::min(gr.min_cap, c);
      }
      if (gr.min_cap == std::numeric_limits<double>::infinity()) gr.min_cap = 0.0;
    }

    if (placement == FlowPlacement::Proportional) {
      FlowWorkspace ws;
      build_reversed_residual(ws, N, groups);
      double pushed_tier = 0.0;
      while (ws.bfs(t, s)) {
        std::fill(ws.it.begin(), ws.it.end(), 0);
        if (shortest_path) {
          double pushed = ws.dfs(t, s, std::numeric_limits<double>::infinity());
          if (pushed < kMinFlow) break;
          pushed_tier += pushed;
          break;
        } else {
          while (true) {
            double pushed = ws.dfs(t, s, std::numeric_limits<double>::infinity());
            if (pushed < kMinFlow) break;
            pushed_tier += pushed;
          }
        }
      }
      if (pushed_tier < kMinFlow) break;
      total += pushed_tier;
      cost_dist[path_cost] += pushed_tier;
      // distribute per-edge and update residual
      for (std::size_t u = 0; u < ws.adj.size(); ++u) {
        for (const auto& e : ws.adj[u]) {
          if (e.group < 0) continue;
          double sent = e.init_cap - e.cap;
          if (sent < kMinFlow) continue;
          const auto& gr = groups[static_cast<std::size_t>(e.group)];
          double denom = gr.sum_cap > 0.0 ? gr.sum_cap : 1.0;
          for (auto eid : gr.eids) {
            // Proportional to residual before update
            double base = residual[static_cast<std::size_t>(eid)];
            double share = sent * (base / denom);
            if (with_edge_flows) edge_flows_accum[static_cast<std::size_t>(eid)] += share;
            residual[static_cast<std::size_t>(eid)] = std::max(0.0, base - share);
          }
        }
      }
    } else { // Equal-balanced
      auto [ratio, assigns] = equal_balance_assign(N, s, groups);
      if (ratio < kMinFlow) break;
      total += ratio;
      cost_dist[path_cost] += ratio;
      // update residual, accumulate flows equally across parallels
      for (auto const& asg : assigns) {
        const auto& gr = groups[asg.group_index];
        if (gr.eids.empty()) continue;
        double flow_scaled = asg.flow * ratio;
        double per_edge = flow_scaled / static_cast<double>(gr.eids.size());
        for (auto eid : gr.eids) {
          if (with_edge_flows) edge_flows_accum[static_cast<std::size_t>(eid)] += per_edge;
          double base = residual[static_cast<std::size_t>(eid)];
          residual[static_cast<std::size_t>(eid)] = std::max(0.0, base - per_edge);
        }
      }
    }
    if (shortest_path) break;
  }

  summary.total_flow = total;
  if (with_edge_flows) summary.edge_flows = edge_flows_accum;
  if (!cost_dist.empty()) {
    std::vector<std::pair<double,double>> pairs(cost_dist.begin(), cost_dist.end());
    std::sort(pairs.begin(), pairs.end(), [](auto const& a, auto const& b){ return a.first < b.first; });
    for (auto const& pr : pairs) { CostBucket b; b.cost = pr.first; b.share = pr.second; summary.cost_distribution.buckets.push_back(b); }
  }
  // Min-cut extraction (proportional/equal-balanced): compute reachability in final residual graph
  // Residual graph has forward residual = residual[e], reverse residual = flow[e] (capacity - residual)
  if (total >= kMinFlow) {
    const auto row = g.row_offsets_view();
    const auto col = g.col_indices_view();
    const auto aei = g.adj_edge_index_view();
    const auto in_row = g.in_row_offsets_view();
    const auto in_col = g.in_col_indices_view();
    const auto in_aei = g.in_adj_edge_index_view();
    std::vector<char> visited(static_cast<std::size_t>(N), 0);
    std::queue<std::int32_t> q;
    q.push(s);
    visited[static_cast<std::size_t>(s)] = 1;
    while (!q.empty()) {
      auto u = q.front(); q.pop();
      if (node_mask && !node_mask[static_cast<std::size_t>(u)]) continue;
      // Forward residual arcs: u -> v if residual > 0
      auto start = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
      auto end   = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
      for (std::size_t j = start; j < end; ++j) {
        auto v = static_cast<std::int32_t>(col[j]);
        auto eid = static_cast<std::size_t>(aei[j]);
        if (edge_mask && !edge_mask[eid]) continue;
        if (node_mask && !node_mask[static_cast<std::size_t>(v)]) continue;
        if (residual[eid] > kMinCap && !visited[static_cast<std::size_t>(v)]) {
          visited[static_cast<std::size_t>(v)] = 1;
          q.push(v);
        }
      }
      // Reverse residual arcs: v -> u if flow(u->v) > 0
      auto rs = static_cast<std::size_t>(in_row[static_cast<std::size_t>(u)]);
      auto re = static_cast<std::size_t>(in_row[static_cast<std::size_t>(u)+1]);
      for (std::size_t j = rs; j < re; ++j) {
        auto w = static_cast<std::int32_t>(in_col[j]);
        auto eid = static_cast<std::size_t>(in_aei[j]);
        if (edge_mask && !edge_mask[eid]) continue;
        if (node_mask && !node_mask[static_cast<std::size_t>(w)]) continue;
        double flow_e = g.capacity_view()[eid] - residual[eid];
        if (flow_e > kMinFlow && !visited[static_cast<std::size_t>(w)]) {
          visited[static_cast<std::size_t>(w)] = 1;
          q.push(w);
        }
      }
    }
    // Collect cut edges: from reachable to non-reachable with no forward residual
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
        if (residual[eid] <= kMinCap) {
          summary.min_cut.edges.push_back(static_cast<std::int32_t>(g.link_id_of(static_cast<std::int32_t>(eid))));
        }
      }
    }
  }
  return {summary.total_flow, std::move(summary)};
}

std::vector<FlowSummary>
batch_max_flow(const StrictMultiDiGraph& g,
               const std::vector<std::pair<NodeId,NodeId>>& pairs,
               FlowPlacement placement, bool shortest_path,
               double eps, bool with_edge_flows,
               int /*threads*/, std::optional<std::uint64_t> /*seed*/,
               const std::vector<const bool*>& node_masks,
               const std::vector<const bool*>& edge_masks) {
  std::vector<FlowSummary> out;
  out.reserve(pairs.size());
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    auto pr = pairs[i];
    const bool* nm = (i < node_masks.size() ? node_masks[i] : nullptr);
    const bool* em = (i < edge_masks.size() ? edge_masks[i] : nullptr);
    auto [val, summary] = calc_max_flow(g, pr.first, pr.second, placement, shortest_path, eps, with_edge_flows, nm, em);
    out.push_back(std::move(summary));
  }
  return out;
}
} // namespace netgraph::core
