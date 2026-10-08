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
  // A small set (8 |edges| < N) is grouped by sorting, O(|edges| log |edges|), with nothing
  // proportional to the graph's size; a larger one by a stable counting sort over a temporary
  // N + 1 array, O(|edges| + N) (a sort first when the input is not ascending). Nothing retained
  // is proportional to the graph's size.
  [[nodiscard]] static PredDAG of_edges(const Graph& g, std::span<const EdgeId> edges);
  // The DAG without entries.
  [[nodiscard]] static PredDAG empty(const Graph& g);
  // A DAG from its finished parts (the searches' emission): `nodes` the children with entries
  // (ascending), `offsets` their entry ranges (nodes.size() + 1, from 0), `parents` and `via` the
  // entries (ascending via within a child), `ascending` the via in ascending EdgeId.
  [[nodiscard]] static PredDAG of_parts(const Graph& g, std::vector<NodeId> nodes, std::vector<std::int32_t> offsets,
                                        std::vector<NodeId> parents, std::vector<EdgeId> via,
                                        std::vector<EdgeId> ascending);
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

// Internal construction of Path values from trusted walks (the searches' results): the checks of
// Path's constructor hold by construction, so none is repeated (no per-path O(E) scan).
struct PathBuilder {
  // A walk of single-member hops [[e1], [e2], ...].
  [[nodiscard]] static Path of_walk(const Graph& g, std::span<const EdgeId> walk);
  // Contiguous hops of parallel members (each hop's members share (u, v); no node repeats).
  [[nodiscard]] static Path of_hops(const Graph& g, std::vector<std::vector<EdgeId>> hops);
};

} // namespace netgraph::core::detail
