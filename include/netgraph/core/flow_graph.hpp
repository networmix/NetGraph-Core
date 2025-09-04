/* Shared flow ledger layering over FlowState with per-flow deltas. */
#pragma once

#include <cstdint>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

#include "netgraph/core/flow_state.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// FlowGraph is a shared, authoritative flow ledger over a StrictMultiDiGraph.
// It composes a FlowState for residual and aggregate edge_flow management, and
// maintains per-flow edge deltas to support exact removal/reopt.
class FlowGraph {
public:
  explicit FlowGraph(const StrictMultiDiGraph& g);

  // Views
  std::span<const Cap> capacity_view() const noexcept { return fs_.capacity_view(); }
  std::span<const Cap> residual_view() const noexcept { return fs_.residual_view(); }
  std::span<const Flow> edge_flow_view() const noexcept { return fs_.edge_flow_view(); }

  // Access underlying graph (const)
  const StrictMultiDiGraph& graph() const noexcept { return *g_; }

  // Placement: applies placement and records per-edge deltas for this flow.
  // Returns placed amount.
  Flow place(const FlowIndex& idx, NodeId src, NodeId dst,
             const PredDAG& dag, Flow amount,
             FlowPlacement placement, bool shortest_path = false);

  // Remove a specific flow, reverting its edge deltas from the ledger.
  void remove(const FlowIndex& idx);

  // Remove all flows belonging to a given flowClass.
  void remove_by_class(std::int32_t flowClass);

  // Reset all state to initial capacity and clear ledger.
  void reset();

  // Inspect: return a copy of the flow's edges and amounts.
  std::vector<std::pair<EdgeId, Flow>> get_flow_edges(const FlowIndex& idx) const;

  // Attempt to reconstruct a single path (LSP) for this flow from the ledger.
  // Returns empty vector if the flow does not correspond to a unique simple path
  // (e.g., when placed with multipath/proportional splitting).
  std::vector<EdgeId> get_flow_path(const FlowIndex& idx) const;

private:
  const StrictMultiDiGraph* g_ {nullptr};
  FlowState fs_;
  // Per-flow ledger: only edges with non-zero assigned flow are stored.
  std::unordered_map<FlowIndex, std::vector<std::pair<EdgeId, Flow>>, FlowIndexHash> ledger_;
};

} // namespace netgraph::core
