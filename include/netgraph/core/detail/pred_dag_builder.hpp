/* Internal construction of PredDAG values from trusted inputs (search output, flow routes,
   Path::to_dag) and read access to their sparse storage. */
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "netgraph/core/pred_dag.hpp"

namespace netgraph::core::detail {

struct PredDAGBuilder {
  // The DAG of an edge set known to be acyclic and duplicate-free: per node its predecessor
  // entries (edge_src, e) in ascending EdgeId. Not validated beyond what the callers guarantee.
  // O(|edges| + N) for ascending input (every caller's), O(|edges| log |edges| + N) otherwise
  // (a sort first); the grouping by child is a stable counting sort over a temporary N + 1
  // array, and nothing retained is proportional to the graph's size.
  [[nodiscard]] static PredDAG of_edges(const Graph& g, std::span<const EdgeId> edges);
  // A Path's DAG: hop h's members (ascending) are the entries of child nodes[h + 1], all with
  // the parent nodes[h]; the nodes are distinct.
  [[nodiscard]] static PredDAG of_path(std::uint64_t graph_uid, std::int32_t num_nodes,
                                       const std::vector<std::vector<EdgeId>>& hops, const std::vector<NodeId>& nodes);
  // The DAG's edges in ascending EdgeId, without a copy (kept by the DAG since construction).
  [[nodiscard]] static std::span<const EdgeId> ascending(const PredDAG& dag) noexcept;
  // The sparse storage: the member nodes (ascending) and their entry offsets (size + 1).
  [[nodiscard]] static std::span<const NodeId> member_nodes(const PredDAG& dag) noexcept;
  [[nodiscard]] static std::span<const std::int32_t> member_offsets(const PredDAG& dag) noexcept;
  // The retained payload in bytes (every stored array; the DagCache's `bytes` statistic).
  [[nodiscard]] static std::size_t payload_bytes(const PredDAG& dag) noexcept;
};

} // namespace netgraph::core::detail
