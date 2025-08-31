#include "netgraph/core/k_shortest_paths.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>

// Use multipath SPF DAG for spur enumeration
#include "netgraph/core/shortest_paths.hpp"

namespace netgraph::core {

namespace {
struct Path {
  std::vector<std::int32_t> nodes;
  std::vector<std::int32_t> edges;
  double cost{0.0};
};
struct Candidate {
  double cost;
  std::vector<std::int32_t> nodes;
  std::vector<std::int32_t> edges;
  bool operator>(const Candidate& o) const { return cost > o.cost; }
};


static std::optional<Path> dijkstra_single(const StrictMultiDiGraph& g, NodeId s, NodeId t,
                                           const std::vector<unsigned char>* node_mask,
                                           const std::vector<unsigned char>* edge_mask) {
  const auto row = g.row_offsets_view();
  const auto col = g.col_indices_view();
  const auto aei = g.adj_edge_index_view();
  const auto cost = g.cost_view();
  const int N = g.num_nodes();
  if (s < 0 || t < 0 || s >= N || t >= N) return std::nullopt;
  auto ok_node = [&](std::int32_t v){ return (!node_mask) || (*node_mask)[static_cast<std::size_t>(v)] != 0; };
  auto ok_edge = [&](std::size_t e){ return (!edge_mask) || (*edge_mask)[e] != 0; };
  if (!ok_node(s) || !ok_node(t)) return std::nullopt;
  std::vector<double> dist(static_cast<std::size_t>(N), std::numeric_limits<double>::infinity());
  std::vector<std::int32_t> parent(static_cast<std::size_t>(N), -1);
  std::vector<std::int32_t> via(static_cast<std::size_t>(N), -1);
  using QItem = std::pair<double, std::int32_t>;
  auto cmp = [](const QItem& a, const QItem& b){ return a.first > b.first; };
  std::priority_queue<QItem, std::vector<QItem>, decltype(cmp)> pq(cmp);
  dist[static_cast<std::size_t>(s)] = 0.0;
  pq.emplace(0.0, s);
  while (!pq.empty()) {
    auto [d_u, u] = pq.top(); pq.pop();
    if (d_u > dist[static_cast<std::size_t>(u)]) continue;
    if (u == t) break;
    auto start = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto end   = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
    std::size_t i = start;
    while (i < end) {
      std::int32_t v = col[i];
      // Skip neighbor group if node masked
      if (!ok_node(v)) { std::size_t j=i; while (j<end && col[j]==v) ++j; i=j; continue; }
      // Find best edge u->v (min cost; tie-breaker: smallest edge id)
      double min_edge_cost = std::numeric_limits<double>::infinity();
      std::int32_t best_eid = -1;
      std::size_t j = i;
      for (; j < end && col[j] == v; ++j) {
        auto eidx = static_cast<std::size_t>(aei[j]);
        if (!ok_edge(eidx)) continue;
        double ec = cost[eidx];
        if (ec < min_edge_cost || (ec == min_edge_cost && static_cast<std::int32_t>(eidx) < best_eid)) {
          min_edge_cost = ec; best_eid = static_cast<std::int32_t>(eidx);
        }
      }
      if (best_eid >= 0) {
        double nd = d_u + min_edge_cost;
        auto vi = static_cast<std::size_t>(v);
        if (nd < dist[vi]) { dist[vi] = nd; parent[vi] = u; via[vi] = best_eid; pq.emplace(nd, v); }
      }
      i = j;
    }
  }
  if (!std::isfinite(dist[static_cast<std::size_t>(t)])) return std::nullopt;
  // Reconstruct
  std::vector<std::int32_t> nodes_rev;
  std::vector<std::int32_t> edges_rev;
  for (std::int32_t v = t; v != s; v = parent[static_cast<std::size_t>(v)]) {
    if (v < 0) return std::nullopt;
    nodes_rev.push_back(v);
    edges_rev.push_back(via[static_cast<std::size_t>(v)]);
  }
  Path p;
  p.nodes.reserve(nodes_rev.size() + 1);
  p.nodes.push_back(s);
  for (auto it = nodes_rev.rbegin(); it != nodes_rev.rend(); ++it) p.nodes.push_back(*it);
  p.edges.reserve(edges_rev.size());
  for (auto it = edges_rev.rbegin(); it != edges_rev.rend(); ++it) p.edges.push_back(*it);
  p.cost = dist[static_cast<std::size_t>(t)];
  return p;
}

// Enumerate all shortest spur paths from spur->t using PredDAG produced by
// shortest_paths(spur, ..., multipath=true). Returns a list of (nodes, edges)
// sequences in forward order.
static void dfs_spur_paths(NodeId spur, NodeId v,
                           const std::vector<std::int32_t>& off,
                           const std::vector<std::int32_t>& parents,
                           const std::vector<std::int32_t>& via,
                           std::vector<std::int32_t>& nodes_rev,
                           std::vector<std::int32_t>& edges_rev,
                           std::vector<std::vector<std::int32_t>>& out_nodes,
                           std::vector<std::vector<std::int32_t>>& out_edges) {
  if (v == spur) {
    // build forward
    std::vector<std::int32_t> n; n.reserve(nodes_rev.size() + 1);
    n.push_back(spur);
    for (auto it = nodes_rev.rbegin(); it != nodes_rev.rend(); ++it) n.push_back(*it);
    std::vector<std::int32_t> e; e.reserve(edges_rev.size());
    for (auto it = edges_rev.rbegin(); it != edges_rev.rend(); ++it) e.push_back(*it);
    out_nodes.push_back(std::move(n));
    out_edges.push_back(std::move(e));
    // restore vectors for caller
    out_nodes.back(); out_edges.back();
    return;
  }
  auto s = off[static_cast<std::size_t>(v)];
  auto e = off[static_cast<std::size_t>(v) + 1];
  for (std::int32_t i = s; i < e; ++i) {
    auto u = parents[static_cast<std::size_t>(i)];
    auto eid = via[static_cast<std::size_t>(i)];
    edges_rev.push_back(eid);
    nodes_rev.push_back(v);
    dfs_spur_paths(spur, u, off, parents, via, nodes_rev, edges_rev, out_nodes, out_edges);
    nodes_rev.pop_back();
    edges_rev.pop_back();
  }
}

} // namespace

std::vector<std::pair<std::vector<double>, PredDAG>> k_shortest_paths(
    const StrictMultiDiGraph& g, NodeId s, NodeId t,
    int k, std::optional<double> max_cost_factor,
    bool unique, double eps) {
  std::vector<Path> paths;
  if (k <= 0) return {};
  if (s < 0 || t < 0 || s >= g.num_nodes() || t >= g.num_nodes()) return {};

  // Base shortest path
  auto p0 = dijkstra_single(g, s, t, nullptr, nullptr);
  if (!p0) return {};
  double best_cost = p0->cost;
  double max_cost = std::numeric_limits<double>::infinity();
  if (max_cost_factor && *max_cost_factor > 0.0) {
    max_cost = best_cost * (*max_cost_factor);
  }
  if (p0->cost <= max_cost) paths.push_back(*p0);
  if (k == 1) {
    // Convert and return
    std::vector<std::pair<std::vector<double>, PredDAG>> items;
    items.reserve(paths.size());
    auto cost_view = g.cost_view();
    for (auto const& P : paths) {
      std::vector<double> dist(static_cast<std::size_t>(g.num_nodes()), std::numeric_limits<double>::infinity());
      PredDAG dag;
      dag.parent_offsets.assign(static_cast<std::size_t>(g.num_nodes() + 1), 0);
      // Fill distances along path and one-parent predecessors
      if (!P.nodes.empty()) {
        dist[static_cast<std::size_t>(P.nodes.front())] = 0.0;
        for (std::size_t i = 1; i < P.nodes.size(); ++i) {
          auto u = P.nodes[i-1]; auto v = P.nodes[i]; auto e = P.edges[i-1];
          dist[static_cast<std::size_t>(v)] = dist[static_cast<std::size_t>(u)] + cost_view[static_cast<std::size_t>(e)];
          dag.parent_offsets[static_cast<std::size_t>(v+1)] = 1;
        }
        for (std::size_t v = 1; v < dag.parent_offsets.size(); ++v) dag.parent_offsets[v] += dag.parent_offsets[v-1];
        dag.parents.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
        dag.via_edges.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
        std::size_t idx = 0;
        for (std::size_t i = 1; i < P.nodes.size(); ++i) {
          dag.parents[idx] = P.nodes[i-1]; dag.via_edges[idx] = P.edges[i-1]; ++idx;
        }
      }
      items.emplace_back(std::move(dist), std::move(dag));
    }
    return items;
  }

  // Candidate heap across spur deviations
  auto cost_view = g.cost_view();
  std::priority_queue<Candidate, std::vector<Candidate>, std::greater<Candidate>> B;
  auto path_signature = [&](const std::vector<std::int32_t>& edges){ return edges; };
  struct VectorHash { size_t operator()(const std::vector<std::int32_t>& v) const noexcept {
    size_t h = 1469598103934665603ull;
    for (auto x : v) { h ^= static_cast<size_t>(x + 0x9e3779b97f4a7c15ull); h *= 1099511628211ull; }
    return h;
  }};
  std::unordered_set<std::vector<std::int32_t>, VectorHash> visited;
  visited.insert(path_signature(p0->edges));

  for (int i = 1; i < k; ++i) {
    const Path& last = paths.back();
    // Precompute prefix cumulative costs along last path to avoid re-summing
    std::vector<double> prefix_cost;
    prefix_cost.assign(last.nodes.size(), 0.0);
    for (std::size_t idx = 1; idx < last.nodes.size(); ++idx) {
      // edge at idx-1
      auto e = last.edges[idx - 1];
      prefix_cost[idx] = prefix_cost[idx - 1] + cost_view[static_cast<std::size_t>(e)];
    }
    // Spur node positions 0..len-2
    for (std::size_t j = 0; j + 1 < last.nodes.size(); ++j) {
      std::int32_t spur_node = last.nodes[j];
      // Build masks: exclude prefix nodes [0..j-1]
      std::vector<unsigned char> node_mask;
      node_mask.assign(static_cast<std::size_t>(g.num_nodes()), static_cast<unsigned char>(1));
      for (std::size_t r = 0; r < j; ++r) node_mask[static_cast<std::size_t>(last.nodes[r])] = 0;
      // Edge mask: exclude next edges of previous accepted paths that share this prefix
      std::vector<unsigned char> edge_mask;
      edge_mask.assign(static_cast<std::size_t>(g.num_edges()), static_cast<unsigned char>(1));
      for (auto const& P : paths) {
        if (P.nodes.size() > j && std::equal(P.nodes.begin(), P.nodes.begin() + j, last.nodes.begin())) {
          // Exclude edge at position j for this path
          if (P.edges.size() > j) {
            edge_mask[static_cast<std::size_t>(P.edges[j])] = 0;
          }
        }
      }
      // Multipath spur PredDAG from spur_node -> t
      std::unique_ptr<bool[]> nm_ptr;
      std::unique_ptr<bool[]> em_ptr;
      const bool* nm = nullptr; const bool* em = nullptr;
      if (!node_mask.empty()) { nm_ptr = std::unique_ptr<bool[]>(new bool[static_cast<std::size_t>(g.num_nodes())]);
        for (std::size_t idx=0; idx<static_cast<std::size_t>(g.num_nodes()); ++idx) nm_ptr[idx] = node_mask[idx] != 0; nm = nm_ptr.get(); }
      if (!edge_mask.empty()) { em_ptr = std::unique_ptr<bool[]>(new bool[static_cast<std::size_t>(g.num_edges())]);
        for (std::size_t idx=0; idx<static_cast<std::size_t>(g.num_edges()); ++idx) em_ptr[idx] = edge_mask[idx] != 0; em = em_ptr.get(); }
      auto [dist_spur, dag_spur] = shortest_paths(g, spur_node, t, EdgeSelect::AllMinCost, /*multipath*/ true, eps, nm, em);
      // If no spur parents for t, skip
      if (static_cast<std::size_t>(t + 1) >= dag_spur.parent_offsets.size() ||
          dag_spur.parent_offsets[static_cast<std::size_t>(t)] == dag_spur.parent_offsets[static_cast<std::size_t>(t + 1)]) {
        continue;
      }
      std::vector<std::int32_t> nodes_rev;
      std::vector<std::int32_t> edges_rev;
      std::vector<std::vector<std::int32_t>> spur_nodes_list;
      std::vector<std::vector<std::int32_t>> spur_edges_list;
      dfs_spur_paths(spur_node, t, dag_spur.parent_offsets, dag_spur.parents, dag_spur.via_edges,
                     nodes_rev, edges_rev, spur_nodes_list, spur_edges_list);
      for (std::size_t si = 0; si < spur_nodes_list.size(); ++si) {
        const auto& spur_nodes = spur_nodes_list[si];
        const auto& spur_edges = spur_edges_list[si];
        std::vector<std::int32_t> cand_nodes;
        std::vector<std::int32_t> cand_edges;
        cand_nodes.reserve(j + spur_nodes.size());
        cand_edges.reserve(j + spur_edges.size());
        for (std::size_t r = 0; r < j; ++r) cand_nodes.push_back(last.nodes[r]);
        for (auto n : spur_nodes) cand_nodes.push_back(n);
        for (std::size_t r = 0; r < j; ++r) cand_edges.push_back(last.edges[r]);
        for (auto e : spur_edges) cand_edges.push_back(e);
        // Compute candidate cost as prefix_cost[j] + sum(spur_edges)
        double cand_cost = prefix_cost[j];
        for (auto e : spur_edges) cand_cost += cost_view[static_cast<std::size_t>(e)];
        if (cand_cost > max_cost) continue;
        auto sig = path_signature(cand_edges);
        if (unique && visited.find(sig) != visited.end()) continue;
        visited.insert(sig);
        B.push(Candidate{cand_cost, std::move(cand_nodes), std::move(cand_edges)});
      }
    }
    // pick next best candidate not exceeding threshold
    if (B.empty()) break;
    bool accepted = false;
    while (!B.empty()) {
      auto c = B.top(); B.pop();
      if (c.cost <= max_cost) {
        paths.push_back(Path{std::move(c.nodes), std::move(c.edges), c.cost});
        accepted = true;
        break;
      }
    }
    if (!accepted) break;
  }
  // Convert to SPF-compatible outputs
  std::vector<std::pair<std::vector<double>, PredDAG>> items;
  items.reserve(paths.size());
  for (auto const& P : paths) {
    std::vector<double> dist(static_cast<std::size_t>(g.num_nodes()), std::numeric_limits<double>::infinity());
    PredDAG dag;
    dag.parent_offsets.assign(static_cast<std::size_t>(g.num_nodes() + 1), 0);
    if (!P.nodes.empty()) {
      dist[static_cast<std::size_t>(P.nodes.front())] = 0.0;
      for (std::size_t i = 1; i < P.nodes.size(); ++i) {
        auto u = P.nodes[i-1]; auto v = P.nodes[i]; auto e = P.edges[i-1];
        dist[static_cast<std::size_t>(v)] = dist[static_cast<std::size_t>(u)] + cost_view[static_cast<std::size_t>(e)];
        dag.parent_offsets[static_cast<std::size_t>(v+1)] = 1;
      }
      for (std::size_t v = 1; v < dag.parent_offsets.size(); ++v) dag.parent_offsets[v] += dag.parent_offsets[v-1];
      dag.parents.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
      dag.via_edges.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
      std::size_t idx = 0;
      for (std::size_t i = 1; i < P.nodes.size(); ++i) {
        dag.parents[idx] = P.nodes[i-1]; dag.via_edges[idx] = P.edges[i-1]; ++idx;
      }
    }
    items.emplace_back(std::move(dist), std::move(dag));
  }
  return items;
}

} // namespace netgraph::core
