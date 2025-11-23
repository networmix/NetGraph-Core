/*
  FlowState — per-edge residual/flow tracking and placement helpers.
  Implements proportional (WCMP-like) and equal-balanced (ECMP) placements
  over shortest-path (Dijkstra) predecessor DAGs.
  ECMP = Equal-Cost Multi-Path; WCMP = Weighted-Cost Multi-Path.
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
// given immutable StrictMultiDiGraph. It places flow on a predecessor DAG (PredDAG)
// using either proportional or equal-balanced strategies, updating internal
// residuals deterministically.
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
  //
  // EqualBalanced semantics here are *single-pass ECMP admission*:
  //  - Use the supplied DAG and fixed equal per-edge splits per (u->v) group.
  //  - Compute one global scale:  min_g cap_rev[g] / assigned[g],
  //      where cap_rev[g] = min_edge_residual(g) * |E_g|
  //        (enforces equal per-edge shares),
  //            assigned[g] = unit-demand load on group g under equal splits.
  //  - Place once and return. We do NOT re-split/recompute after a bottleneck
  //    saturates. Re-invoking this on the updated residuals changes the effective
  //    next-hop set (progressive traffic-engineering behavior) and is outside
  //    "single-pass ECMP admission".
  [[nodiscard]] Flow place_on_dag(NodeId src, NodeId dst,
                    const PredDAG& dag,
                    Flow requested_flow,
                    FlowPlacement placement,
                    // Optional trace collector to record per-edge allocations applied by this call
                    std::vector<std::pair<EdgeId, Flow>>* trace = nullptr);

  // Convenience: run repeated placements until exhaustion (or single tier when
  // shortest_path=true). Returns total placed flow. Uses internal residual.
  //
  // require_capacity: Whether to require edges to have capacity.
  //   - true (default): Routes adapt to residuals (SDN/TE behavior).
  //   - false: Routes based on costs only (IP/IGP behavior).
  //
  // NOTE: With FlowPlacement::EqualBalanced + require_capacity=true, this behaves as a
  // progressive "fill": it recomputes the SPF DAG after each saturation. That is *not*
  // the single-pass ECMP admission used by place_on_dag(). For IP ECMP, use
  // require_capacity=false + shortest_path=true + EqualBalanced.
  [[nodiscard]] Flow place_max_flow(NodeId src, NodeId dst,
                      FlowPlacement placement,
                      bool shortest_path = false,
                      bool require_capacity = true,
                      std::span<const bool> node_mask = {},
                      std::span<const bool> edge_mask = {});

  // Compute min-cut (edges crossing from reachable set S to unreachable set T)
  // based on reachability in the current residual graph.
  // Reachability starts from source s; forward arcs allowed if residual > kMinCap,
  // reverse arcs allowed if flow > kMinFlow.
  [[nodiscard]] MinCut compute_min_cut(NodeId src,
                         std::span<const bool> node_mask = {},
                         std::span<const bool> edge_mask = {}) const;

  // Apply or revert a set of edge flow allocations directly.
  // When add==true, treats each (eid, flow) as additional placed flow on the edge.
  // When add==false, removes placed flow (reverts allocations), clamping to [0, capacity].
  void apply_deltas(std::span<const std::pair<EdgeId, Flow>> deltas, bool add) noexcept;

private:
  const StrictMultiDiGraph* g_ {nullptr};
  std::vector<Cap> residual_;
  std::vector<Flow> edge_flow_;
};

} // namespace netgraph::core
