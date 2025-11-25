/* Max-flow utility APIs with summaries and batch evaluation. */
#pragma once

#include <optional>
#include <cstdint>
#include <utility>
#include <vector>

#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

struct MinCut { std::vector<EdgeId> edges; };

struct FlowSummary {
  Flow total_flow {0.0};
  MinCut min_cut {};
  // Parallel arrays: costs[i] has placed flow flows[i] at that total path cost.
  // Costs are ascending and unique; flows are raw volumes (not normalized).
  std::vector<Cost> costs;
  std::vector<Flow> flows;
  std::vector<Flow> edge_flows; // filled if requested
  // Optional large arrays; populated only when requested via API flags
  std::vector<Cap> residual_capacity; // length == g.num_edges()
  std::vector<std::uint8_t> reachable_nodes; // length == g.num_nodes(); 0/1 flags
};

[[nodiscard]] std::pair<Flow, FlowSummary>
calc_max_flow(const StrictMultiDiGraph& g, NodeId src, NodeId dst,
              FlowPlacement placement, bool shortest_path,
              bool require_capacity,
              bool with_edge_flows,
              bool with_reachable,
              bool with_residuals,
              std::span<const bool> node_mask = {},
              std::span<const bool> edge_mask = {});

[[nodiscard]] std::vector<FlowSummary>
batch_max_flow(const StrictMultiDiGraph& g,
               const std::vector<std::pair<NodeId,NodeId>>& pairs,
               FlowPlacement placement, bool shortest_path,
               bool require_capacity,
               bool with_edge_flows,
               bool with_reachable,
               bool with_residuals,
               const std::vector<std::span<const bool>>& node_masks = {},
               const std::vector<std::span<const bool>>& edge_masks = {});

// Performs sensitivity analysis to identify edges that constrain flow.
//
// Computes baseline flow (using full max-flow or shortest-path-only mode),
// then tests removing each saturated edge to measure its criticality.
//
// Arguments:
//   g: The input graph.
//   src, dst: Source and destination nodes.
//   placement: Flow placement strategy.
//   shortest_path: If true, uses single-pass shortest-path flow (IP/IGP mode).
//                  If false, uses full iterative max-flow (SDN/TE mode).
//   require_capacity: If true, excludes saturated edges from routing.
//   node_mask, edge_mask: Optional masks to exclude nodes/edges.
//
// Returns:
//   Pairs of (EdgeId, FlowDelta) for edges whose removal reduces flow.
[[nodiscard]] std::vector<std::pair<EdgeId, Flow>>
sensitivity_analysis(const StrictMultiDiGraph& g, NodeId src, NodeId dst,
                     FlowPlacement placement, bool shortest_path,
                     bool require_capacity,
                     std::span<const bool> node_mask = {},
                     std::span<const bool> edge_mask = {});

} // namespace netgraph::core
