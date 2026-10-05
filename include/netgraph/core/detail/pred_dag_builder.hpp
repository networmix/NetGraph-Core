/* Internal construction of PredDAG values (search output, kernel routes, Path::to_dag). */
#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "netgraph/core/pred_dag.hpp"

namespace netgraph::core::detail {

struct PredDAGBuilder {
  // The DAG of an edge set known to be acyclic and duplicate-free: per node its predecessor
  // entries (edge_src, e) in ascending EdgeId. Not validated beyond what the callers guarantee.
  [[nodiscard]] static PredDAG from_edges(const StrictMultiDiGraph& g, std::span<const EdgeId> edges);
  // Per child its ascending via members, all with the one parent parent_of[child].
  [[nodiscard]] static PredDAG from_children(std::uint64_t graph_uid, const std::vector<std::vector<EdgeId>>& per_child,
                                             const std::vector<NodeId>& parent_of);
  // The DAG's edges in ascending EdgeId, without a copy (kept by the DAG since construction).
  [[nodiscard]] static std::span<const EdgeId> ascending(const PredDAG& dag) noexcept;
};

} // namespace netgraph::core::detail
