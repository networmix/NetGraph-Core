#include "netgraph/core/shortest_paths.hpp"
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <queue>
#include <utility>
#include <vector>

namespace netgraph::core {

namespace {
// Threshold for remaining capacity checks (matches NetGraph MIN_CAP default).
constexpr double kMinCap = 1.0 / 4096.0;

inline bool nearly_equal(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}
} // namespace

std::pair<std::vector<double>, PredDAG>
shortest_paths(const StrictMultiDiGraph& g, NodeId src,
               std::optional<NodeId> dst,
               EdgeSelect policy, bool multipath, double eps,
               const bool* node_mask,
               const bool* edge_mask) {
  const auto N = g.num_nodes();
  // Supported policies: ALL_MIN_COST, SINGLE_MIN_COST,
  // ALL_MIN_COST_WITH_CAP_REMAINING, ALL_ANY_COST_WITH_CAP_REMAINING,
  // SINGLE_MIN_COST_WITH_CAP_REMAINING
  // Note: eps controls equality for cost ties. Python bindings default to 1e-10;
  // pass a smaller value (e.g., 1e-12) to mimic NetGraph's isclose in tests.
  if (src < 0 || src >= N) {
    // Return empty structures if invalid source
    return {std::vector<double>(static_cast<std::size_t>(N), std::numeric_limits<double>::infinity()), PredDAG{}};
  }
  // If node mask excludes source, return empty
  if (node_mask && !node_mask[static_cast<std::size_t>(src)]) {
    return {std::vector<double>(static_cast<std::size_t>(N), std::numeric_limits<double>::infinity()), PredDAG{}};
  }

  const auto row = g.row_offsets_view();
  const auto col = g.col_indices_view();
  const auto aei = g.adj_edge_index_view();
  const auto cost = g.cost_view();
  const auto cap  = g.capacity_view();

  std::vector<double> dist(static_cast<std::size_t>(N), std::numeric_limits<double>::infinity());
  dist[static_cast<std::size_t>(src)] = 0.0;

  // predecessor lists per node: vector of pairs (parent node, list of via edges)
  std::vector<std::vector<std::pair<NodeId, std::vector<EdgeId>>>> pred_lists(static_cast<std::size_t>(N));
  pred_lists[static_cast<std::size_t>(src)] = {};

  using QItem = std::pair<double, NodeId>;
  auto cmp = [](const QItem& a, const QItem& b) { return a.first > b.first; };
  std::priority_queue<QItem, std::vector<QItem>, decltype(cmp)> pq(cmp);
  pq.emplace(0.0, src);

  double best_dst_cost = std::numeric_limits<double>::infinity();
  bool have_best_dst = false;
  const bool early_exit = dst.has_value();
  const NodeId dst_node = dst.value_or(-1);

  while (!pq.empty()) {
    auto [d_u, u] = pq.top(); pq.pop();
    if (d_u > dist[static_cast<std::size_t>(u)] + eps) continue;

    // Record best destination distance when popped first time
    if (early_exit && u == dst_node && !have_best_dst) {
      best_dst_cost = d_u;
      have_best_dst = true;
    }

    // Do not expand from destination itself, but allow equal-cost preds to be collected
    if (early_exit && u == dst_node) {
      // check heap head for early termination
      if (pq.empty() || pq.top().first > best_dst_cost + eps) break;
      else continue;
    }

    // Explore neighbors of u
    auto start = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto end   = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);

    // Neighbor selection: compute min edge cost and collect edges per neighbor v.
    // Because edges are compacted/sorted by (src, dst, ...), all parallel edges u->v
    // are consecutive; we group them in a single pass without extra allocations.

    // Iterate consecutive segments by v
    std::size_t i = start;
    while (i < end) {
      NodeId v = col[i];
      // Skip entire neighbor group if node v is masked out
      if (node_mask && !node_mask[static_cast<std::size_t>(v)]) {
        std::size_t j_skip = i;
        while (j_skip < end && col[j_skip] == v) ++j_skip;
        i = j_skip;
        continue;
      }
      double min_edge_cost = std::numeric_limits<double>::infinity();
      std::vector<EdgeId> selected_edges;

      // Scan all parallel edges u->v (consecutive due to compaction)
      std::size_t j = i;
      for (; j < end && col[j] == v; ++j) {
        auto e = static_cast<std::size_t>(aei[j]);
        // Skip masked-out edges
        if (edge_mask && !edge_mask[e]) continue;
        const double ecap = cap[e];
        const double ecost = cost[e];
        // Capacity-screening policies (no flow tracking; treat cap as remaining)
        if (policy == EdgeSelect::AllMinCostWithCapRemaining ||
            policy == EdgeSelect::AllAnyCostWithCapRemaining ||
            policy == EdgeSelect::SingleMinCostWithCapRemaining) {
          if (ecap < kMinCap) continue;
        }

        if (policy == EdgeSelect::AllAnyCostWithCapRemaining) {
          // Collect all edges with cap; min_edge_cost still determines dist
          if (ecost + eps < min_edge_cost) {
            min_edge_cost = ecost;
          }
          selected_edges.push_back(static_cast<EdgeId>(aei[j]));
        } else {
          // Min-cost based policies
          if (ecost + eps < min_edge_cost) {
            min_edge_cost = ecost;
            selected_edges.clear();
            selected_edges.push_back(static_cast<EdgeId>(aei[j]));
          } else if (nearly_equal(ecost, min_edge_cost, eps)) {
            if (policy == EdgeSelect::AllMinCost || policy == EdgeSelect::AllMinCostWithCapRemaining) {
              selected_edges.push_back(static_cast<EdgeId>(aei[j]));
            } else {
              // SINGLE_* keep only one representative (first by compaction order)
              // nothing to do
            }
          }
        }
      }

      if (!selected_edges.empty()) {
        double new_cost = d_u + min_edge_cost;
        auto v_idx = static_cast<std::size_t>(v);
        if (new_cost + eps < dist[v_idx]) {
          dist[v_idx] = new_cost;
          pred_lists[v_idx].clear();
          pred_lists[v_idx].push_back({u, std::move(selected_edges)});
          pq.emplace(new_cost, v);
        } else if (multipath && nearly_equal(new_cost, dist[v_idx], eps)) {
          pred_lists[v_idx].push_back({u, std::move(selected_edges)});
        }
      }

      i = j;
    }

    if (have_best_dst) {
      if (pq.empty() || pq.top().first > best_dst_cost + eps) break;
    }
  }

  // Build PredDAG from pred_lists: flatten each (parent, edges[]) as entries per edge
  PredDAG dag;
  dag.parent_offsets.assign(static_cast<std::size_t>(N + 1), 0);
  // count total entries per node
  for (std::int32_t v = 0; v < N; ++v) {
    std::size_t count = 0;
    for (auto const& pe : pred_lists[static_cast<std::size_t>(v)]) count += pe.second.size();
    dag.parent_offsets[static_cast<std::size_t>(v + 1)] = static_cast<std::int32_t>(count);
  }
  for (std::size_t i2 = 1; i2 < dag.parent_offsets.size(); ++i2) {
    dag.parent_offsets[i2] += dag.parent_offsets[i2 - 1];
  }
  dag.parents.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  dag.via_edges.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  for (std::int32_t v = 0; v < N; ++v) {
    auto base = static_cast<std::size_t>(dag.parent_offsets[static_cast<std::size_t>(v)]);
    std::size_t k = 0;
    for (auto const& pe : pred_lists[static_cast<std::size_t>(v)]) {
      NodeId p = pe.first;
      for (auto e : pe.second) {
        dag.parents[base + k] = p;
        dag.via_edges[base + k] = e;
        ++k;
      }
    }
  }

  return {std::move(dist), std::move(dag)};
}

std::pair<std::vector<double>, PredDAG>
shortest_paths_with_residual(const StrictMultiDiGraph& g, NodeId src,
                             std::optional<NodeId> dst,
                             EdgeSelect policy, bool multipath, double eps,
                             const std::vector<double>& residual,
                             const bool* node_mask,
                             const bool* edge_mask) {
  const auto N = g.num_nodes();
  const auto row = g.row_offsets_view();
  const auto col = g.col_indices_view();
  const auto aei = g.adj_edge_index_view();
  const auto cost = g.cost_view();
  const auto cap  = g.capacity_view();

  std::vector<double> dist(static_cast<std::size_t>(N), std::numeric_limits<double>::infinity());
  if (src >= 0 && src < N && (!node_mask || node_mask[static_cast<std::size_t>(src)])) {
    dist[static_cast<std::size_t>(src)] = 0.0;
  }
  std::vector<std::vector<std::pair<NodeId, std::vector<EdgeId>>>> pred_lists(static_cast<std::size_t>(N));
  if (src >= 0 && src < N) pred_lists[static_cast<std::size_t>(src)] = {};

  using QItem = std::pair<double, NodeId>;
  auto cmp = [](const QItem& a, const QItem& b) { return a.first > b.first; };
  std::priority_queue<QItem, std::vector<QItem>, decltype(cmp)> pq(cmp);
  pq.emplace(0.0, src);
  double best_dst_cost = std::numeric_limits<double>::infinity();
  bool have_best_dst = false;
  const bool early_exit = dst.has_value();
  const NodeId dst_node = dst.value_or(-1);

  auto nearly = [&](double a, double b){ return std::abs(a-b) <= eps; };

  while (!pq.empty()) {
    auto [d_u, u] = pq.top(); pq.pop();
    if (u < 0 || u >= N) continue;
    if (d_u > dist[static_cast<std::size_t>(u)] + eps) continue;
    if (early_exit && u == dst_node && !have_best_dst) { best_dst_cost = d_u; have_best_dst = true; }
    if (early_exit && u == dst_node) {
      if (pq.empty() || pq.top().first > best_dst_cost + eps) break; else continue;
    }

    auto start = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto end   = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
    std::size_t i = start;
    while (i < end) {
      NodeId v = col[i];
      if (node_mask && !node_mask[static_cast<std::size_t>(v)]) {
        std::size_t j_skip = i; while (j_skip < end && col[j_skip] == v) ++j_skip; i = j_skip; continue;
      }
      double min_edge_cost = std::numeric_limits<double>::infinity();
      std::vector<EdgeId> selected_edges;
      std::size_t j = i;
      for (; j < end && col[j] == v; ++j) {
        auto e = static_cast<std::size_t>(aei[j]);
        if (edge_mask && !edge_mask[e]) continue;
        const double rem = residual[e];
        if (rem < kMinCap) continue;
        const double ecost = cost[e];
        if (policy == EdgeSelect::AllAnyCostWithCapRemaining) {
          if (ecost + eps < min_edge_cost) min_edge_cost = ecost;
          selected_edges.push_back(static_cast<EdgeId>(aei[j]));
        } else if (policy == EdgeSelect::SingleMinCostWithCapRemainingLoadFactored) {
          // Compute cost factor: cost*100 + round(flow/capacity*10)
          const double c = cap[e];
          const double flow = std::max(0.0, c - rem);
          const double load_factor = c > 0.0 ? std::round((flow / c) * 10.0) : 999999.0;
          const double eff_cost = ecost * 100.0 + load_factor;
          if (eff_cost + eps < min_edge_cost) {
            min_edge_cost = eff_cost;
            selected_edges.clear();
            selected_edges.push_back(static_cast<EdgeId>(aei[j]));
          }
        } else {
          // Min-cost w/ cap remaining (ALL or SINGLE)
          if (ecost + eps < min_edge_cost) {
            min_edge_cost = ecost;
            selected_edges.clear();
            selected_edges.push_back(static_cast<EdgeId>(aei[j]));
          } else if (nearly(ecost, min_edge_cost)) {
            if (policy == EdgeSelect::AllMinCostWithCapRemaining) {
              selected_edges.push_back(static_cast<EdgeId>(aei[j]));
            }
          }
        }
      }
      if (!selected_edges.empty()) {
        double new_cost;
        if (policy == EdgeSelect::SingleMinCostWithCapRemainingLoadFactored) {
          // Distance accumulates base edge cost, not load factor
          // Recover base min cost among selected edges for distance
          double min_base = std::numeric_limits<double>::infinity();
          for (auto eid : selected_edges) {
            min_base = std::min(min_base, cost[static_cast<std::size_t>(eid)]);
          }
          new_cost = d_u + min_base;
        } else {
          new_cost = d_u + ((min_edge_cost==std::numeric_limits<double>::infinity()) ? 0.0 : min_edge_cost);
        }
        auto v_idx = static_cast<std::size_t>(v);
        if (new_cost + eps < dist[v_idx]) { dist[v_idx] = new_cost; pred_lists[v_idx].clear(); pred_lists[v_idx].push_back({u, std::move(selected_edges)}); pq.emplace(new_cost, v); }
        else if (multipath && nearly(new_cost, dist[v_idx])) { pred_lists[v_idx].push_back({u, std::move(selected_edges)}); }
      }
      i = j;
    }
    if (have_best_dst) { if (pq.empty() || pq.top().first > best_dst_cost + eps) break; }
  }

  PredDAG dag; dag.parent_offsets.assign(static_cast<std::size_t>(N+1), 0);
  for (std::int32_t v=0; v<N; ++v) { std::size_t c=0; for (auto const& pe : pred_lists[static_cast<std::size_t>(v)]) c += pe.second.size(); dag.parent_offsets[static_cast<std::size_t>(v+1)] = static_cast<std::int32_t>(c); }
  for (std::size_t k=1; k<dag.parent_offsets.size(); ++k) dag.parent_offsets[k] += dag.parent_offsets[k-1];
  dag.parents.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  dag.via_edges.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  for (std::int32_t v=0; v<N; ++v) {
    auto base = static_cast<std::size_t>(dag.parent_offsets[static_cast<std::size_t>(v)]);
    std::size_t k = 0;
    for (auto const& pe : pred_lists[static_cast<std::size_t>(v)]) {
      NodeId p = pe.first; for (auto e : pe.second) { dag.parents[base+k] = p; dag.via_edges[base+k] = e; ++k; }
    }
  }

  return {std::move(dist), std::move(dag)};
}

} // namespace netgraph::core
