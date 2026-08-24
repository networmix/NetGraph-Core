/*
  k_shortest_paths — Yen-like enumeration with deterministic tie-breaking.

  Returns SPF-compatible outputs for each path: full distance array and a
  PredDAG that encodes a single concrete path as single-parent predecessors.
*/
#include "netgraph/core/k_shortest_paths.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <queue>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

// Use multipath SPF DAG for spur enumeration
#include "netgraph/core/shortest_paths.hpp"

namespace netgraph::core {

namespace {
struct Path {
  std::vector<std::int32_t> nodes;
  std::vector<EdgeId> edges;
  Cost cost{0};
};
struct Candidate {
  Cost cost;
  std::vector<std::int32_t> nodes;
  std::vector<EdgeId> edges;
  bool operator>(const Candidate& o) const { return cost > o.cost; }
};

struct VectorHash { size_t operator()(const std::vector<EdgeId>& v) const noexcept {
  size_t h = 1469598103934665603ull;
  for (auto x : v) { h ^= static_cast<size_t>(x + 0x9e3779b97f4a7c15ull); h *= 1099511628211ull; }
  return h;
}};
using PathSet = std::unordered_set<std::vector<EdgeId>, VectorHash>;


static std::optional<Path> dijkstra_single(const StrictMultiDiGraph& g, NodeId s, NodeId t,
                                           const std::vector<std::uint8_t>* node_mask,
                                           const std::vector<std::uint8_t>* edge_mask) {
  const auto row = g.row_offsets_view();
  const auto col = g.col_indices_view();
  const auto aei = g.adj_edge_index_view();
  const auto cost = g.cost_view();
  const int N = g.num_nodes();
  if (s < 0 || t < 0 || s >= N || t >= N) return std::nullopt;
  auto ok_node = [&](std::int32_t v){ return (!node_mask) || (*node_mask)[static_cast<std::size_t>(v)] != 0; };
  auto ok_edge = [&](std::size_t e){ return (!edge_mask) || (*edge_mask)[e] != 0; };
  if (!ok_node(s) || !ok_node(t)) return std::nullopt;
  std::vector<Cost> dist(static_cast<std::size_t>(N), std::numeric_limits<Cost>::max());
  std::vector<std::int32_t> parent(static_cast<std::size_t>(N), -1);
  std::vector<std::int32_t> via(static_cast<std::size_t>(N), -1);
  using QItem = std::pair<Cost, std::int32_t>;
  auto cmp = [](const QItem& a, const QItem& b){ return a.first > b.first; };
  std::priority_queue<QItem, std::vector<QItem>, decltype(cmp)> pq(cmp);
  dist[static_cast<std::size_t>(s)] = static_cast<Cost>(0);
  pq.emplace(static_cast<Cost>(0), s);
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
      Cost min_edge_cost = std::numeric_limits<Cost>::max();
      std::int32_t best_eid = -1;
      std::size_t j = i;
      for (; j < end && col[j] == v; ++j) {
        auto eidx = static_cast<std::size_t>(aei[j]);
        if (!ok_edge(eidx)) continue;
        Cost ec = static_cast<Cost>(cost[eidx]);
        if (ec < min_edge_cost || (ec == min_edge_cost && static_cast<std::int32_t>(eidx) < best_eid)) {
          min_edge_cost = ec; best_eid = static_cast<std::int32_t>(eidx);
        }
      }
      if (best_eid >= 0) {
        Cost nd = static_cast<Cost>(d_u + min_edge_cost);
        auto vi = static_cast<std::size_t>(v);
        if (nd < dist[vi]) { dist[vi] = nd; parent[vi] = u; via[vi] = best_eid; pq.emplace(nd, v); }
      }
      i = j;
    }
  }
  if (dist[static_cast<std::size_t>(t)] == std::numeric_limits<Cost>::max()) return std::nullopt;
  // Reconstruct
  std::vector<std::int32_t> nodes_rev;
  std::vector<EdgeId> edges_rev;
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

// Enumerate shortest spur paths from spur->t using the PredDAG produced by
// shortest_paths(spur, ..., multipath=true), in deterministic DFS order.
// All spur paths in the DAG share the same (shortest) cost, and at most
// max_count candidates can ever be accepted downstream, so enumeration stops
// after emitting max_count paths. Without the bound this walk materializes
// every equal-cost path — exponential in ECMP fan-out (2^k paths on a k-stage
// ladder). When dedup is supplied (unique mode), paths whose full edge
// sequence (prefix + spur segment) was already visited are skipped without
// consuming the bound, so the cap cannot starve fresh candidates.
static void dfs_spur_paths(NodeId spur, NodeId v,
                           const std::vector<std::int32_t>& off,
                           const std::vector<std::int32_t>& parents,
                           const std::vector<std::int32_t>& via,
                           std::vector<std::int32_t>& nodes_rev,
                           std::vector<EdgeId>& edges_rev,
                           std::vector<std::vector<std::int32_t>>& out_nodes,
                           std::vector<std::vector<EdgeId>>& out_edges,
                           std::size_t max_count,
                           const std::vector<EdgeId>* prefix_edges,
                           const PathSet* dedup) {
  if (out_nodes.size() >= max_count) return;
  if (v == spur) {
    // build forward
    std::vector<std::int32_t> n; n.reserve(nodes_rev.size() + 1);
    n.push_back(spur);
    for (auto it = nodes_rev.rbegin(); it != nodes_rev.rend(); ++it) n.push_back(*it);
    std::vector<EdgeId> e; e.reserve(edges_rev.size());
    for (auto it = edges_rev.rbegin(); it != edges_rev.rend(); ++it) e.push_back(*it);
    if (dedup != nullptr) {
      std::vector<EdgeId> full;
      full.reserve(prefix_edges->size() + e.size());
      full.insert(full.end(), prefix_edges->begin(), prefix_edges->end());
      full.insert(full.end(), e.begin(), e.end());
      if (dedup->find(full) != dedup->end()) return;  // already visited; free the slot
    }
    out_nodes.push_back(std::move(n));
    out_edges.push_back(std::move(e));
    return;
  }
  auto s = off[static_cast<std::size_t>(v)];
  auto e = off[static_cast<std::size_t>(v) + 1];
  for (std::int32_t i = s; i < e; ++i) {
    if (out_nodes.size() >= max_count) return;
    auto u = parents[static_cast<std::size_t>(i)];
    auto eid = via[static_cast<std::size_t>(i)];
    edges_rev.push_back(static_cast<EdgeId>(eid));
    nodes_rev.push_back(v);
    dfs_spur_paths(spur, u, off, parents, via, nodes_rev, edges_rev, out_nodes, out_edges,
                   max_count, prefix_edges, dedup);
    nodes_rev.pop_back();
    edges_rev.pop_back();
  }
}

} // namespace

std::vector<std::pair<std::vector<Cost>, PredDAG>> k_shortest_paths(
    const StrictMultiDiGraph& g, NodeId src, NodeId dst,
    int k, std::optional<double> max_cost_factor,
    bool unique,
    std::span<const bool> node_mask,
    std::span<const bool> edge_mask) {
  std::vector<Path> paths;
  if (k <= 0) return {};
  if (src < 0 || dst < 0 || src >= g.num_nodes() || dst >= g.num_nodes()) return {};

  if (!node_mask.empty() && node_mask.size() != static_cast<std::size_t>(g.num_nodes())) {
    throw std::invalid_argument("k_shortest_paths: node_mask length mismatch");
  }
  if (!edge_mask.empty() && edge_mask.size() != static_cast<std::size_t>(g.num_edges())) {
    throw std::invalid_argument("k_shortest_paths: edge_mask length mismatch");
  }

  // Base shortest path
  std::vector<std::uint8_t> node_mask_vec;
  std::vector<std::uint8_t> edge_mask_vec;
  const std::vector<std::uint8_t>* nm_ptr = nullptr;
  const std::vector<std::uint8_t>* em_ptr = nullptr;
  if (node_mask.size() == static_cast<std::size_t>(g.num_nodes())) {
    node_mask_vec.assign(static_cast<std::size_t>(g.num_nodes()), static_cast<std::uint8_t>(1));
    for (std::size_t i=0;i<node_mask_vec.size();++i) node_mask_vec[i] = node_mask[i] ? 1u : 0u;
    nm_ptr = &node_mask_vec;
  }
  if (edge_mask.size() == static_cast<std::size_t>(g.num_edges())) {
    edge_mask_vec.assign(static_cast<std::size_t>(g.num_edges()), static_cast<std::uint8_t>(1));
    for (std::size_t i=0;i<edge_mask_vec.size();++i) edge_mask_vec[i] = edge_mask[i] ? 1u : 0u;
    em_ptr = &edge_mask_vec;
  }
  auto p0 = dijkstra_single(g, src, dst, nm_ptr, em_ptr);
  if (!p0) return {};
  double best_cost = static_cast<double>(p0->cost);
  double max_cost = std::numeric_limits<double>::infinity();
  if (max_cost_factor && *max_cost_factor > 0.0) {
    max_cost = best_cost * (*max_cost_factor);
  }
  if (p0->cost <= max_cost) paths.push_back(*p0);
  // A max_cost_factor below 1.0 puts the ceiling under the shortest path itself, so
  // nothing was admitted. Return now: the k > 1 loop below opens with paths.back(),
  // which is undefined behaviour on an empty vector. (The k == 1 branch would have
  // masked this, so the guard must sit ahead of it.)
  if (paths.empty()) return {};
  if (k == 1) {
    // Convert and return
    std::vector<std::pair<std::vector<Cost>, PredDAG>> items;
    items.reserve(paths.size());
    for (auto const& P : paths) {
      items.push_back(path_to_pred_dag(g, P.nodes, P.edges));
    }
    return items;
  }

  // Candidate heap across spur deviations
  auto cost_view = g.cost_view();
  std::priority_queue<Candidate, std::vector<Candidate>, std::greater<Candidate>> B;
  auto path_signature = [&](const std::vector<EdgeId>& edges){ return edges; };
  PathSet visited;
  visited.insert(path_signature(p0->edges));

  for (int i = 1; i < k; ++i) {
    const Path& last = paths.back();
    // Precompute prefix cumulative costs along last path to avoid re-summing
    std::vector<Cost> prefix_cost;
    prefix_cost.assign(last.nodes.size(), 0);
    for (std::size_t idx = 1; idx < last.nodes.size(); ++idx) {
      // edge at idx-1
      auto e = last.edges[idx - 1];
      prefix_cost[idx] = prefix_cost[idx - 1] + static_cast<Cost>(cost_view[static_cast<std::size_t>(e)]);
    }
    // Spur node positions 0..len-2
    for (std::size_t j = 0; j + 1 < last.nodes.size(); ++j) {
      std::int32_t spur_node = last.nodes[j];
      // Build masks: exclude prefix nodes [0..j-1]
      std::vector<std::uint8_t> node_mask_local;
      node_mask_local.assign(static_cast<std::size_t>(g.num_nodes()), static_cast<std::uint8_t>(1));
      for (std::size_t r = 0; r < j; ++r) node_mask_local[static_cast<std::size_t>(last.nodes[r])] = 0;
      if (nm_ptr) {
        for (std::size_t idx=0; idx<node_mask_local.size(); ++idx) node_mask_local[idx] = (node_mask_local[idx] && (*nm_ptr)[idx]) ? 1u : 0u;
      }
      // Edge mask: exclude next edges of already-accepted paths that share this prefix
      std::vector<std::uint8_t> edge_mask_local;
      edge_mask_local.assign(static_cast<std::size_t>(g.num_edges()), static_cast<std::uint8_t>(1));
      for (auto const& P : paths) {
        if (P.nodes.size() > j && std::equal(P.nodes.begin(), P.nodes.begin() + static_cast<std::ptrdiff_t>(j + 1), last.nodes.begin())) {
          // Exclude edge at position j for this path
          if (P.edges.size() > j) {
            edge_mask_local[static_cast<std::size_t>(P.edges[j])] = 0;
          }
        }
      }
      if (em_ptr) {
        for (std::size_t idx=0; idx<edge_mask_local.size(); ++idx) edge_mask_local[idx] = (edge_mask_local[idx] && (*em_ptr)[idx]) ? 1u : 0u;
      }
      // Multipath spur PredDAG from spur_node -> t
      std::unique_ptr<bool[]> nm_buf2;
      std::unique_ptr<bool[]> em_buf2;
      std::span<const bool> nm, em;
      if (!node_mask_local.empty()) {
        nm_buf2 = std::unique_ptr<bool[]>(new bool[static_cast<std::size_t>(g.num_nodes())]);
        for (std::size_t idx = 0; idx < static_cast<std::size_t>(g.num_nodes()); ++idx) nm_buf2[idx] = (node_mask_local[idx] != 0);
        nm = std::span<const bool>(nm_buf2.get(), static_cast<std::size_t>(g.num_nodes()));
      }
      if (!edge_mask_local.empty()) {
        em_buf2 = std::unique_ptr<bool[]>(new bool[static_cast<std::size_t>(g.num_edges())]);
        for (std::size_t idx = 0; idx < static_cast<std::size_t>(g.num_edges()); ++idx) em_buf2[idx] = (edge_mask_local[idx] != 0);
        em = std::span<const bool>(em_buf2.get(), static_cast<std::size_t>(g.num_edges()));
      }
      EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
      auto [dist_spur, dag_spur] = shortest_paths(g, spur_node, dst, /*multipath=*/true, sel, std::span<const Cap>(), nm, em);
      // If no spur parents for t, skip
      if (static_cast<std::size_t>(dst + 1) >= dag_spur.parent_offsets.size() ||
          dag_spur.parent_offsets[static_cast<std::size_t>(dst)] == dag_spur.parent_offsets[static_cast<std::size_t>(dst + 1)]) {
        continue;
      }
      // Every spur path in the DAG costs dist_spur[dst]; skip the whole family
      // if the resulting candidates would exceed the cost ceiling.
      if (static_cast<double>(prefix_cost[j]) + static_cast<double>(dist_spur[static_cast<std::size_t>(dst)]) > max_cost) {
        continue;
      }
      // At most (k - accepted) further paths can ever be accepted, so that many
      // fresh candidates from this spur family suffice.
      const std::size_t want = static_cast<std::size_t>(k) - paths.size();
      std::vector<EdgeId> prefix_edges(last.edges.begin(),
                                       last.edges.begin() + static_cast<std::ptrdiff_t>(j));
      std::vector<std::int32_t> nodes_rev;
      std::vector<std::int32_t> edges_rev;
      std::vector<std::vector<std::int32_t>> spur_nodes_list;
      std::vector<std::vector<std::int32_t>> spur_edges_list;
      dfs_spur_paths(spur_node, dst, dag_spur.parent_offsets, dag_spur.parents, dag_spur.via_edges,
                     nodes_rev, edges_rev, spur_nodes_list, spur_edges_list,
                     want, &prefix_edges, unique ? &visited : nullptr);
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
        Cost cand_cost = prefix_cost[j];
        for (auto e : spur_edges) cand_cost += static_cast<Cost>(cost_view[static_cast<std::size_t>(e)]);
        if (cand_cost > max_cost) continue;
        auto sig = path_signature(cand_edges);
        if (unique && visited.find(sig) != visited.end()) continue;
        visited.insert(sig);
        B.push(Candidate{cand_cost, std::move(cand_nodes), std::move(cand_edges)});
      }
    }
    // Pick next best candidate not exceeding threshold
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
  std::vector<std::pair<std::vector<Cost>, PredDAG>> items;
  items.reserve(paths.size());
  for (auto const& P : paths) {
    items.push_back(path_to_pred_dag(g, P.nodes, P.edges));
  }
  return items;
}

} // namespace netgraph::core
