/*
  calc_max_flow — incremental max-flow via repeated SPF tiers over residuals.

  Uses FlowState to track residuals and place flow along SPF predecessor DAGs.
  Accumulates total flow, optional per-edge flows, cost distribution, and
  derives a min-cut from the final residual graph.

  Algorithm: Iterative augmentation using Successive Shortest Path on residual graphs,
  pushing flow across full ECMP/WCMP DAGs at each step.
*/
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/flow_state.hpp"
#include "netgraph/core/constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

namespace netgraph::core {

std::pair<Flow, FlowSummary>
calc_max_flow(const StrictMultiDiGraph& g, NodeId src, NodeId dst,
              FlowPlacement placement, bool shortest_path,
              bool require_capacity,
              bool with_edge_flows,
              bool with_reachable,
              bool with_residuals,
              std::span<const bool> node_mask, std::span<const bool> edge_mask) {
  FlowSummary summary;
  const auto N = g.num_nodes();
  if (src < 0 || src >= N || dst < 0 || dst >= N || src == dst) {
    return {0.0, std::move(summary)};
  }
  if (!node_mask.empty() && node_mask.size() != static_cast<std::size_t>(g.num_nodes())) {
    throw std::invalid_argument("calc_max_flow: node_mask length mismatch");
  }
  if (!edge_mask.empty() && edge_mask.size() != static_cast<std::size_t>(g.num_edges())) {
    throw std::invalid_argument("calc_max_flow: edge_mask length mismatch");
  }
  const bool use_node_mask = !node_mask.empty();
  const bool use_edge_mask = !edge_mask.empty();
  if ((use_node_mask && !node_mask[static_cast<std::size_t>(src)]) ||
      (use_node_mask && !node_mask[static_cast<std::size_t>(dst)])) {
    return {0.0, std::move(summary)};
  }
  // Use FlowState for residual tracking and per-edge placement
  FlowState fs(g);
  Flow total = static_cast<Flow>(0.0);
  std::vector<std::pair<Cost,Flow>> cost_dist; // (cost, flow)

  // Iterate tiers: SPF over current residual or costs only
  // require_capacity=true: Require edges to have capacity, exclude saturated links (SDN/TE behavior)
  // require_capacity=false: Routes based on costs only, ignore capacity (IP/IGP behavior)
  while (true) {
    EdgeSelection sel;
    sel.multi_edge = true;
    sel.require_capacity = require_capacity;
    sel.tie_break = EdgeTieBreak::Deterministic;
    auto [dist, dag] = shortest_paths(
        g, src, dst,
        /*multipath=*/true,
        sel,
        require_capacity ? fs.residual_view() : std::span<const Cap>{},
        use_node_mask ? node_mask : std::span<const bool>{},
        use_edge_mask ? edge_mask : std::span<const bool>{});

    // No path if t has no parents in DAG
    if (static_cast<std::size_t>(dst) >= dag.parent_offsets.size() - 1 ||
        dag.parent_offsets[static_cast<std::size_t>(dst)] == dag.parent_offsets[static_cast<std::size_t>(dst) + 1]) {
      break;
    }

    Cost path_cost = dist[static_cast<std::size_t>(dst)];
    Flow placed = fs.place_on_dag(src, dst, dag, std::numeric_limits<double>::infinity(), placement);

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
  if (with_residuals) {
    auto res = fs.residual_view();
    summary.residual_capacity.resize(static_cast<std::size_t>(g.num_edges()));
    for (std::size_t i = 0; i < summary.residual_capacity.size(); ++i) {
      summary.residual_capacity[i] = static_cast<Cap>(res[i]);
    }
  }
  if (!cost_dist.empty()) {
    std::sort(cost_dist.begin(), cost_dist.end(), [](auto const& a, auto const& b){ return a.first < b.first; });
    summary.costs.reserve(cost_dist.size());
    summary.flows.reserve(cost_dist.size());
    for (auto const& pr : cost_dist) { summary.costs.push_back(pr.first); summary.flows.push_back(pr.second); }
  }
  // Min-Cut: Set of edges crossing from the reachable set (S) to the unreachable set (T).
  // Computed by finding all nodes reachable from src in the residual graph.
  if (total >= kMinFlow) {
    auto mc = fs.compute_min_cut(src, node_mask, edge_mask);
    summary.min_cut = mc;
    if (with_reachable) {
      summary.reachable_nodes.assign(static_cast<std::size_t>(g.num_nodes()), 0u);
      auto residual = fs.residual_view();
      auto capv = fs.capacity_view();
      const bool reach_use_node_mask = use_node_mask;
      const bool reach_use_edge_mask = use_edge_mask;
      const auto N = static_cast<std::size_t>(g.num_nodes());
      std::vector<std::int32_t> stack;
      stack.reserve(N);
      stack.push_back(src);
      while (!stack.empty()) {
        auto n = static_cast<std::size_t>(stack.back());
        stack.pop_back();
        if (summary.reachable_nodes[n]) continue;
        if (reach_use_node_mask && !node_mask[n]) continue;
        summary.reachable_nodes[n] = 1u;
        auto ro = g.row_offsets_view();
        auto ci = g.col_indices_view();
        auto ae = g.adj_edge_index_view();
        auto s = static_cast<std::size_t>(ro[n]);
        auto e = static_cast<std::size_t>(ro[n+1]);
        for (std::size_t p = s; p < e; ++p) {
          auto v = static_cast<std::size_t>(ci[p]);
          auto eid = static_cast<std::size_t>(ae[p]);
          if (reach_use_edge_mask && !edge_mask[eid]) continue;
          if (reach_use_node_mask && !node_mask[v]) continue;
          if (residual[eid] > kMinCap && !summary.reachable_nodes[v]) {
            stack.push_back(static_cast<std::int32_t>(v));
          }
        }
        auto iro = g.in_row_offsets_view();
        auto ici = g.in_col_indices_view();
        auto iae = g.in_adj_edge_index_view();
        auto rs = static_cast<std::size_t>(iro[n]);
        auto re = static_cast<std::size_t>(iro[n+1]);
        for (std::size_t p = rs; p < re; ++p) {
          auto u = static_cast<std::size_t>(ici[p]);
          auto eid = static_cast<std::size_t>(iae[p]);
          if (reach_use_edge_mask && !edge_mask[eid]) continue;
          if (reach_use_node_mask && !node_mask[u]) continue;
          auto flow = capv[eid] - residual[eid];
          if (flow > kMinFlow && !summary.reachable_nodes[u]) {
            stack.push_back(static_cast<std::int32_t>(u));
          }
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
               bool require_capacity,
               bool with_edge_flows,
               bool with_reachable,
               bool with_residuals,
               const std::vector<std::span<const bool>>& node_masks,
               const std::vector<std::span<const bool>>& edge_masks) {
  // Validate all masks upfront
  for (const auto& nm : node_masks) {
    if (!nm.empty() && nm.size() != static_cast<std::size_t>(g.num_nodes())) {
      throw std::invalid_argument("batch_max_flow: node_mask length mismatch");
    }
  }
  for (const auto& em : edge_masks) {
    if (!em.empty() && em.size() != static_cast<std::size_t>(g.num_edges())) {
      throw std::invalid_argument("batch_max_flow: edge_mask length mismatch");
    }
  }

  std::vector<FlowSummary> out;
  out.reserve(pairs.size());
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    auto pr = pairs[i];
    std::span<const bool> nm = (i < node_masks.size() ? node_masks[i] : std::span<const bool>{});
    std::span<const bool> em = (i < edge_masks.size() ? edge_masks[i] : std::span<const bool>{});
    auto [val, summary] = calc_max_flow(g, pr.first, pr.second,
                                        placement, shortest_path,
                                        require_capacity,
                                        with_edge_flows,
                                        with_reachable,
                                        with_residuals,
                                        nm, em);
    out.push_back(std::move(summary));
  }
  return out;
}

std::vector<std::pair<EdgeId, Flow>>
sensitivity_analysis(const StrictMultiDiGraph& g, NodeId src, NodeId dst,
                     FlowPlacement placement, bool shortest_path,
                     bool require_capacity,
                     std::span<const bool> node_mask,
                     std::span<const bool> edge_mask) {
  // Step 1: Baseline analysis to identify saturated edges
  // Uses shortest_path mode to match the routing semantics being analyzed.
  auto [baseline_flow, summary] = calc_max_flow(
      g, src, dst, placement,
      shortest_path,
      require_capacity,
      /*with_edge_flows=*/false,
      /*with_reachable=*/false,
      /*with_residuals=*/true,
      node_mask, edge_mask);

  if (baseline_flow < kMinFlow || summary.residual_capacity.empty()) {
    return {};
  }

  std::vector<EdgeId> candidates;
  // Identify saturated edges (residual <= kMinCap)
  // Only consider edges that are not already masked out
  for (EdgeId eid = 0; eid < static_cast<EdgeId>(summary.residual_capacity.size()); ++eid) {
    if (!edge_mask.empty() && !edge_mask[eid]) continue; // Already masked
    if (summary.residual_capacity[eid] <= kMinCap) {
      candidates.push_back(eid);
    }
  }

  if (candidates.empty()) {
    return {};
  }

  // Prepare mutable mask buffer
  const auto N = static_cast<std::size_t>(g.num_edges());
  std::unique_ptr<bool[]> test_mask_buf(new bool[N]);
  if (edge_mask.empty()) {
    std::fill(test_mask_buf.get(), test_mask_buf.get() + N, true);
  } else {
    std::copy(edge_mask.begin(), edge_mask.end(), test_mask_buf.get());
  }
  // View for passing to calc_max_flow
  std::span<const bool> test_mask_span(test_mask_buf.get(), N);

  std::vector<std::pair<EdgeId, Flow>> results;
  results.reserve(candidates.size());

  // Step 2: Iterate candidates, testing flow reduction when each is removed
  for (EdgeId eid : candidates) {
    // Mask out the edge
    test_mask_buf[eid] = false;

    auto [new_flow, _] = calc_max_flow(
        g, src, dst, placement,
        shortest_path,
        require_capacity,
        /*with_edge_flows=*/false,
        /*with_reachable=*/false,
        /*with_residuals=*/false,
        node_mask, test_mask_span);

    double delta = baseline_flow - new_flow;
    if (delta > kMinFlow) {
      results.emplace_back(eid, delta);
    }

    // Restore mask
    test_mask_buf[eid] = true;
  }

  return results;
}

} // namespace netgraph::core
