/* FlowGraph manages per-flow edge allocations over FlowState. */
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

// FlowGraph manages per-flow edge allocations over a StrictMultiDiGraph.
// Composes FlowState for residual/aggregate edge flow management.
class FlowGraph {
public:
  explicit FlowGraph(const StrictMultiDiGraph& g);
  ~FlowGraph() noexcept = default;

  // Views
  [[nodiscard]] std::span<const Cap> capacity_view() const noexcept { return fs_.capacity_view(); }
  [[nodiscard]] std::span<const Cap> residual_view() const noexcept { return fs_.residual_view(); }
  [[nodiscard]] std::span<const Flow> edge_flow_view() const noexcept { return fs_.edge_flow_view(); }

  // Access underlying graph (const)
  [[nodiscard]] const StrictMultiDiGraph& graph() const noexcept { return *g_; }

// Apply placement and record per-edge allocations for this flow. Returns placed amount.
  [[nodiscard]] Flow place(const FlowIndex& idx, NodeId src, NodeId dst,
             const PredDAG& dag, Flow amount,
             FlowPlacement placement);

  // Remove a specific flow, reverting its edge allocations from the ledger.
  void remove(const FlowIndex& idx);

  // Remove all flows belonging to a given flowClass.
  void remove_by_class(FlowClass flowClass);

  // Reset all state to initial capacity and clear ledger.
  void reset() noexcept;

  // Inspect: return a copy of the flow's edges and amounts.
  [[nodiscard]] std::vector<std::pair<EdgeId, Flow>> get_flow_edges(const FlowIndex& idx) const;

// Reconstruct single path for this flow from ledger.
// Returns empty vector if flow uses multipath/proportional splitting.
  [[nodiscard]] std::vector<EdgeId> get_flow_path(const FlowIndex& idx) const;

private:
  const StrictMultiDiGraph* g_ {nullptr};
  FlowState fs_;
  // Per-flow ledger: stores only edges with non-zero flow
  std::unordered_map<FlowIndex, std::vector<std::pair<EdgeId, Flow>>, FlowIndexHash> ledger_;
};

} // namespace netgraph::core
