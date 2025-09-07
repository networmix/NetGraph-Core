/*
  FlowState — per-edge residual/flow tracking and placement helpers.
  Implements proportional and equal-balanced placements over SPF DAGs.
*/
#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"
#include "netgraph/core/max_flow.hpp"

namespace netgraph::core {

// FlowState maintains per-edge residual capacity and per-edge placed flow for a
// given immutable StrictMultiDiGraph. It provides efficient placement on a
// provided predecessor DAG (PredDAG) using either proportional or equal-balanced
// strategies, updating internal residuals deterministically.
class FlowState {
public:
  explicit FlowState(const StrictMultiDiGraph& g);
  FlowState(const StrictMultiDiGraph& g, std::span<const Cap> residual_init);
  ~FlowState() noexcept = default;

  // Reset residual to the graph's initial capacities and clear edge_flow.
  void reset() noexcept;
  void reset(std::span<const Cap> residual_init);

  // Views over internal buffers (length == g.num_edges()).
  [[nodiscard]] std::span<const Cap> capacity_view() const noexcept { return g_->capacity_view(); }
  [[nodiscard]] std::span<const Cap> residual_view() const noexcept { return residual_; }
  [[nodiscard]] std::span<const Flow> edge_flow_view() const noexcept { return edge_flow_; }

  // Mutating placement along a given PredDAG tier between src and dst.
  // requested_flow may be +inf. Returns the amount actually placed.
  [[nodiscard]] Flow place_on_dag(NodeId src, NodeId dst,
                    const PredDAG& dag,
                    Flow requested_flow,
                    FlowPlacement placement,
                    bool shortest_path = false,
                    // Optional trace collector to record per-edge deltas applied by this call
                    std::vector<std::pair<EdgeId, Flow>>* trace = nullptr);

  // Convenience: run repeated placements until exhaustion (or single tier when
  // shortest_path=true). Returns total placed flow. Uses internal residual.
  [[nodiscard]] Flow place_max_flow(NodeId src, NodeId dst,
                      FlowPlacement placement,
                      bool shortest_path = false);

  // Compute min-cut with respect to current residual state, starting reachability
  // from source s on the residual graph (forward arcs: residual>MIN; reverse arcs:
  // positive flow). Honors optional masks.
  [[nodiscard]] MinCut compute_min_cut(NodeId src,
                         const bool* node_mask = nullptr,
                         const bool* edge_mask = nullptr) const;

  // Apply or revert a set of edge flow deltas directly.
  // When add==true, treats each (eid, flow) as additional placed flow on the edge.
  // When add==false, removes previously placed flow (reverts), clamping to [0, capacity].
  void apply_deltas(std::span<const std::pair<EdgeId, Flow>> deltas, bool add) noexcept;

private:
  const StrictMultiDiGraph* g_ {nullptr};
  std::vector<Cap> residual_;
  std::vector<Flow> edge_flow_;
};

} // namespace netgraph::core
