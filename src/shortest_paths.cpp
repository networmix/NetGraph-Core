/*
  shortest_paths — Dijkstra over a StrictMultiDiGraph with flexible selection.

  Features:
    - Optional residual-aware traversal (treat residual as capacity gate).
    - Multipath mode collects all equal-cost predecessor edges per node.
    - Single-path mode chooses one best edge per neighbor group with
      deterministic tie-breaking (or by higher residual if requested).
    - Optional early exit when a specific destination is provided.
*/
#include "netgraph/core/shortest_paths.hpp"
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <queue>
#include <utility>
#include <vector>
#include "netgraph/core/constants.hpp"

namespace netgraph::core {

namespace {
static std::pair<std::vector<Cost>, PredDAG>
shortest_paths_core(const StrictMultiDiGraph& g, NodeId src,
                    std::optional<NodeId> dst,
                    const EdgeSelection& selection,
                    std::span<const Cap> residual,
                    const bool* node_mask,
                    const bool* edge_mask) {
  const auto N = g.num_nodes();
  const auto row = g.row_offsets_view();
  const auto col = g.col_indices_view();
  const auto aei = g.adj_edge_index_view();
  const auto cost = g.cost_view();
  const auto cap  = g.capacity_view();

  std::vector<Cost> dist(static_cast<std::size_t>(N), std::numeric_limits<Cost>::max());
  if (src >= 0 && src < N && (!node_mask || node_mask[static_cast<std::size_t>(src)])) {
    dist[static_cast<std::size_t>(src)] = static_cast<Cost>(0);
  }
  std::vector<std::vector<std::pair<NodeId, std::vector<EdgeId>>>> pred_lists(static_cast<std::size_t>(N));
  if (src >= 0 && src < N) pred_lists[static_cast<std::size_t>(src)] = {};

  using QItem = std::pair<Cost, NodeId>;
  auto cmp = [](const QItem& a, const QItem& b) { return a.first > b.first; };
  std::priority_queue<QItem, std::vector<QItem>, decltype(cmp)> pq(cmp);
  pq.emplace(static_cast<Cost>(0), src);
  Cost best_dst_cost = std::numeric_limits<Cost>::max();
  bool have_best_dst = false;
  const bool early_exit = dst.has_value();
  const NodeId dst_node = dst.value_or(-1);

  const bool has_residual = (residual.size() == static_cast<std::size_t>(g.num_edges()));
  const bool require_cap = selection.require_capacity || has_residual;
  const bool multipath = selection.multipath;

  while (!pq.empty()) {
    auto [d_u, u] = pq.top(); pq.pop();
    if (u < 0 || u >= N) continue;
    if (d_u > dist[static_cast<std::size_t>(u)]) continue;
    if (early_exit && u == dst_node && !have_best_dst) { best_dst_cost = d_u; have_best_dst = true; }
    if (early_exit && u == dst_node) {
      if (pq.empty() || pq.top().first > best_dst_cost) break; else continue;
    }

    auto start = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto end   = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
    std::size_t i = start;
    while (i < end) {
      NodeId v = col[i];
      if (node_mask && !node_mask[static_cast<std::size_t>(v)]) {
        std::size_t j_skip = i; while (j_skip < end && col[j_skip] == v) ++j_skip; i = j_skip; continue;
      }
      Cost min_edge_cost = std::numeric_limits<Cost>::max();
      std::vector<EdgeId> selected_edges;
      double best_rem_for_min_cost = -1.0;
      std::size_t j = i;
      int best_edge_id = -1;
      for (; j < end && col[j] == v; ++j) {
        auto e = static_cast<std::size_t>(aei[j]);
        if (edge_mask && !edge_mask[e]) continue;
        const Cap rem = has_residual ? residual[e] : cap[e];
        if (require_cap && rem < kMinCap) continue;
        const Cost ecost = static_cast<Cost>(cost[e]);
        if (ecost < min_edge_cost) {
          min_edge_cost = ecost;
          selected_edges.clear();
          if (multipath) {
            selected_edges.push_back(static_cast<EdgeId>(aei[j]));
          } else {
            best_edge_id = static_cast<int>(e);
            best_rem_for_min_cost = static_cast<double>(rem);
          }
        } else if (ecost == min_edge_cost) {
          if (multipath) {
            selected_edges.push_back(static_cast<EdgeId>(aei[j]));
          } else {
            // tie-break among equal-cost edges
            if (selection.tie_break == EdgeTieBreak::PreferHigherResidual) {
              if (static_cast<double>(rem) > best_rem_for_min_cost + 1e-18) {
                best_edge_id = static_cast<int>(e);
                best_rem_for_min_cost = static_cast<double>(rem);
              } else if (std::abs(static_cast<double>(rem) - best_rem_for_min_cost) <= 1e-18) {
                // further tie-break deterministically by smaller edge id
                if (best_edge_id < 0 || static_cast<int>(e) < best_edge_id) {
                  best_edge_id = static_cast<int>(e);
                }
              }
            } else {
              // Deterministic: smallest edge id
              if (best_edge_id < 0 || static_cast<int>(e) < best_edge_id) {
                best_edge_id = static_cast<int>(e);
              }
            }
          }
        }
      }
      if (!multipath && best_edge_id >= 0) {
        selected_edges.clear();
        selected_edges.push_back(static_cast<EdgeId>(best_edge_id));
      }
      if (!selected_edges.empty()) {
        Cost new_cost = static_cast<Cost>(d_u + (min_edge_cost==std::numeric_limits<Cost>::max() ? 0 : min_edge_cost));
        auto v_idx = static_cast<std::size_t>(v);
        if (new_cost < dist[v_idx]) { dist[v_idx] = new_cost; pred_lists[v_idx].clear(); pred_lists[v_idx].push_back({u, std::move(selected_edges)}); pq.emplace(new_cost, v); }
        else if (multipath && new_cost == dist[v_idx]) { pred_lists[v_idx].push_back({u, std::move(selected_edges)}); }
      }
      i = j;
    }
    if (have_best_dst) { if (pq.empty() || pq.top().first > best_dst_cost) break; }
  }

  PredDAG dag; dag.parent_offsets.assign(static_cast<std::size_t>(N+1), 0);
  for (std::int32_t v=0; v<N; ++v) { std::size_t c=0; for (auto const& pe : pred_lists[static_cast<std::size_t>(v)]) c += pe.second.size(); dag.parent_offsets[static_cast<std::size_t>(v+1)] = static_cast<std::int32_t>(c); }
  for (std::size_t k=1; k<dag.parent_offsets.size(); ++k) dag.parent_offsets[k] += dag.parent_offsets[k-1];
  dag.parents.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  dag.via_edges.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  for (std::int32_t v=0; v<N; ++v) {
    auto base = static_cast<std::size_t>(dag.parent_offsets[static_cast<std::size_t>(v)]);
    std::size_t k = 0;
    for (auto const& pe : pred_lists[static_cast<std::size_t>(v)]) { NodeId p = pe.first; for (auto e : pe.second) { dag.parents[base+k] = p; dag.via_edges[base+k] = e; ++k; } }
  }
  return {std::move(dist), std::move(dag)};
}
} // namespace

std::pair<std::vector<Cost>, PredDAG>
shortest_paths(const StrictMultiDiGraph& g, NodeId src,
               std::optional<NodeId> dst,
               const EdgeSelection& selection,
               std::span<const Cap> residual,
               const bool* node_mask,
               const bool* edge_mask) {
  return shortest_paths_core(g, src, dst, selection, residual, node_mask, edge_mask);
}

} // namespace netgraph::core
