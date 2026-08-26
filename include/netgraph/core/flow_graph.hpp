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

  // Copies and moves take a FRESH uid (and version 0). uid must be unique per
  // reachable instance: a shared uid would let two objects that diverge after
  // copying present equal StateStamps over different residual content, and
  // FlowPolicy's memo fast path would then serve a DAG cached for the other
  // object's state. Moved-from objects stay valid and mutable, so moves must
  // not share the uid either; the one-time cache miss this costs is nothing.
  FlowGraph(const FlowGraph& other);
  FlowGraph& operator=(const FlowGraph& other);
  FlowGraph(FlowGraph&& other) noexcept;
  FlowGraph& operator=(FlowGraph&& other) noexcept;

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

  // Monotonic stamp identifying the exact residual state of this FlowGraph.
  // uid is process-unique per instance (never reused across FlowGraph
  // lifetimes); version bumps on every mutation that can change residuals.
  // Equal stamps from two observations therefore guarantee identical residual
  // content, which FlowPolicy uses to elide repeated identical SPF runs.
  struct StateStamp {
    std::uint64_t uid;
    std::uint64_t version;
    bool operator==(const StateStamp&) const noexcept = default;
  };
  [[nodiscard]] StateStamp state_stamp() const noexcept { return {uid_, version_}; }

// Reconstruct single path for this flow from ledger.
// Returns empty vector if flow uses multipath/proportional splitting.
  [[nodiscard]] std::vector<EdgeId> get_flow_path(const FlowIndex& idx) const;

private:
  const StrictMultiDiGraph* g_ {nullptr};
  FlowState fs_;
  // Per-flow ledger: stores only edges with non-zero flow
  std::unordered_map<FlowIndex, std::vector<std::pair<EdgeId, Flow>>, FlowIndexHash> ledger_;
  // See state_stamp(). fs_ is private and every mutation path (place, remove,
  // remove_by_class, reset) bumps version_, so the stamp cannot go stale.
  std::uint64_t uid_ {0};
  std::uint64_t version_ {0};
};

} // namespace netgraph::core
