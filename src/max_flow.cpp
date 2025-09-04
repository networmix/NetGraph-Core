/*
  calc_max_flow — incremental max-flow via repeated SPF tiers over residuals.

  Uses FlowState to track residuals and place flow along SPF predecessor DAGs.
  Accumulates total flow, optional per-edge flows, cost distribution, and
  derives a min-cut from the final residual graph.
*/
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/flow_state.hpp"
#include "netgraph/core/constants.hpp"

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
} // namespace

std::pair<Flow, FlowSummary>
calc_max_flow(const StrictMultiDiGraph& g, NodeId src, NodeId dst,
              FlowPlacement placement, bool shortest_path,
              bool with_edge_flows,
              const bool* node_mask, const bool* edge_mask) {
  FlowSummary summary;
  const auto N = g.num_nodes();
  if (src < 0 || src >= N || dst < 0 || dst >= N || src == dst) {
    return {0.0, std::move(summary)};
  }
  if ((node_mask && !node_mask[static_cast<std::size_t>(src)]) ||
      (node_mask && !node_mask[static_cast<std::size_t>(dst)])) {
    return {0.0, std::move(summary)};
  }
  // Use FlowState for residual tracking and per-edge placement
  FlowState fs(g);
  Flow total = static_cast<Flow>(0.0);
  std::vector<std::pair<Cost,Flow>> cost_dist; // (cost, share)

  // Iterate tiers: SPF over current residual -> place on DAG -> accumulate stats
  while (true) {
    EdgeSelection sel; sel.multipath = true; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
    auto [dist, dag] = shortest_paths(
        g, src, dst,
        sel,
        fs.residual_view(),
        node_mask, edge_mask);

    // No path if t has no parents in DAG
    if (static_cast<std::size_t>(dst) >= dag.parent_offsets.size() - 1 ||
        dag.parent_offsets[static_cast<std::size_t>(dst)] == dag.parent_offsets[static_cast<std::size_t>(dst) + 1]) {
      break;
    }

    Cost path_cost = dist[static_cast<std::size_t>(dst)];
    Flow placed = fs.place_on_dag(src, dst, dag, std::numeric_limits<double>::infinity(), placement, shortest_path);
    if (placed < kMinFlow) break;
    total += placed;
    // Merge by exact cost (integer)
    bool merged = false;
    for (auto& pr : cost_dist) { if (pr.first == path_cost) { pr.second += placed; merged = true; break; } }
    if (!merged) cost_dist.emplace_back(path_cost, placed);
    if (shortest_path) break;
  }

  summary.total_flow = total;
  if (with_edge_flows) {
    auto ef = fs.edge_flow_view();
    summary.edge_flows.assign(ef.begin(), ef.end());
  }
  if (!cost_dist.empty()) {
    std::sort(cost_dist.begin(), cost_dist.end(), [](auto const& a, auto const& b){ return a.first < b.first; });
    for (auto const& pr : cost_dist) { CostBucket b; b.cost = pr.first; b.share = pr.second; summary.cost_distribution.buckets.push_back(b); }
  }
  // Min-cut extraction (proportional/equal-balanced): compute reachability in final residual graph
  // Residual graph has forward residual = residual[e], reverse residual = flow[e] (capacity - residual)
  if (total >= kMinFlow) {
    summary.min_cut = fs.compute_min_cut(src, node_mask, edge_mask);
  }
  return {summary.total_flow, std::move(summary)};
}

std::vector<FlowSummary>
batch_max_flow(const StrictMultiDiGraph& g,
               const std::vector<std::pair<NodeId,NodeId>>& pairs,
               FlowPlacement placement, bool shortest_path,
               bool with_edge_flows,
               const std::vector<const bool*>& node_masks,
               const std::vector<const bool*>& edge_masks) {
  std::vector<FlowSummary> out;
  out.reserve(pairs.size());
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    auto pr = pairs[i];
    const bool* nm = (i < node_masks.size() ? node_masks[i] : nullptr);
    const bool* em = (i < edge_masks.size() ? edge_masks[i] : nullptr);
    auto [val, summary] = calc_max_flow(g, pr.first, pr.second, placement, shortest_path, with_edge_flows, nm, em);
    out.push_back(std::move(summary));
  }
  return out;
}
} // namespace netgraph::core
