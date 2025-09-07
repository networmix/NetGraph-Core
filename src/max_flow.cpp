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
              bool with_reachable,
              bool with_residuals,
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
  std::vector<std::pair<Cost,Flow>> cost_dist; // (cost, flow)

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
  // Min-cut extraction (proportional/equal-balanced): compute reachability in final residual graph
  // Residual graph has forward residual = residual[e], reverse residual = flow[e] (capacity - residual)
  if (total >= kMinFlow) {
    auto mc = fs.compute_min_cut(src, node_mask, edge_mask);
    summary.min_cut = mc;
    if (with_reachable) {
      // Reconstruct reachable set as a boolean vector of size N from min-cut computation
      // FlowState::compute_min_cut records only cut edges; recompute reachability quickly
      summary.reachable_nodes.assign(static_cast<std::size_t>(g.num_nodes()), 0u);
      // BFS over residual graph using FlowState views
      auto residual = fs.residual_view();
      auto capv = fs.capacity_view();
      const auto N = static_cast<std::size_t>(g.num_nodes());
      std::vector<std::int32_t> stack;
      stack.reserve(N);
      stack.push_back(src);
      while (!stack.empty()) {
        auto n = static_cast<std::size_t>(stack.back());
        stack.pop_back();
        if (summary.reachable_nodes[n]) continue;
        summary.reachable_nodes[n] = 1u;
        // Forward residual arcs
        auto ro = g.row_offsets_view();
        auto ci = g.col_indices_view();
        auto ae = g.adj_edge_index_view();
        auto s = static_cast<std::size_t>(ro[n]);
        auto e = static_cast<std::size_t>(ro[n+1]);
        for (std::size_t p = s; p < e; ++p) {
          auto v = static_cast<std::size_t>(ci[p]);
          auto eid = static_cast<std::size_t>(ae[p]);
          if (residual[eid] > kMinCap && !summary.reachable_nodes[v]) {
            stack.push_back(static_cast<std::int32_t>(v));
          }
        }
        // Reverse residual arcs (flow > 0)
        auto iro = g.in_row_offsets_view();
        auto ici = g.in_col_indices_view();
        auto iae = g.in_adj_edge_index_view();
        auto rs = static_cast<std::size_t>(iro[n]);
        auto re = static_cast<std::size_t>(iro[n+1]);
        for (std::size_t p = rs; p < re; ++p) {
          auto u = static_cast<std::size_t>(ici[p]);
          auto eid = static_cast<std::size_t>(iae[p]);
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
               bool with_edge_flows,
               bool with_reachable,
               bool with_residuals,
               const std::vector<const bool*>& node_masks,
               const std::vector<const bool*>& edge_masks) {
  std::vector<FlowSummary> out;
  out.reserve(pairs.size());
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    auto pr = pairs[i];
    const bool* nm = (i < node_masks.size() ? node_masks[i] : nullptr);
    const bool* em = (i < edge_masks.size() ? edge_masks[i] : nullptr);
    auto [val, summary] = calc_max_flow(g, pr.first, pr.second,
                                        placement, shortest_path,
                                        with_edge_flows,
                                        with_reachable,
                                        with_residuals,
                                        nm, em);
    out.push_back(std::move(summary));
  }
  return out;
}
} // namespace netgraph::core
