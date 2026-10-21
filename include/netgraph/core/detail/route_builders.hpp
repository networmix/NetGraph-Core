/* Internal construction of PredecessorDAG values from trusted inputs (search output, flow routes,
   Path::to_dag) and read access to their sparse storage. */
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "netgraph/core/route_structures.hpp"

namespace netgraph::core::detail {

// Input checks shared by native factories and Python's attached argument pass.
void check_csr(const Graph& graph, std::span<const std::int32_t> offsets,
               std::span<const NodeId> parents, std::span<const EdgeId> via);
void check_hops(const Graph& graph, const std::vector<std::vector<EdgeId>>& hops);

struct PredecessorDAGBuilder {
  // The same arrays after check_csr; builds owned storage without repeating validation.
  [[nodiscard]] static PredecessorDAG of_csr(const Graph& graph, std::span<const std::int32_t> offsets,
      std::span<const NodeId> parents, std::span<const EdgeId> via);
  // Build a DAG from acyclic, duplicate-free edges; callers guarantee these preconditions.
  // Group predecessor entries (edge source, e) by node, with ascending EdgeId within each node.
  // For graph node count N and 8 |edges| < N, sort in O(|edges| log |edges|) without graph-sized
  // scratch. Otherwise use stable counting sort with an N + 1 temporary in O(|edges| + N),
  // sorting first if edges are not ascending. Retained storage never scales with graph size.
  [[nodiscard]] static PredecessorDAG of_edges(const Graph& graph, std::span<const EdgeId> edges);
  // The DAG without entries.
  [[nodiscard]] static PredecessorDAG empty(const Graph& graph);
  // A DAG from its finished parts (the searches' emission): `nodes` the children with entries
  // (ascending), `offsets` their entry ranges (nodes.size() + 1, from 0), `parents` and `via` the
  // entries (ascending via within a child), `ascending` the via in ascending EdgeId.
  [[nodiscard]] static PredecessorDAG of_parts(const Graph& graph, std::vector<NodeId> nodes, std::vector<std::int32_t> offsets,
                                        std::vector<NodeId> parents, std::vector<EdgeId> via,
                                        std::vector<EdgeId> ascending);
  // A Path's DAG: hop h's members (ascending) are the entries of child nodes[h + 1], all with
  // the parent nodes[h]; the nodes are distinct.
  [[nodiscard]] static PredecessorDAG of_path(std::uint64_t graph_uid, std::int32_t num_nodes,
                                       const std::vector<std::vector<EdgeId>>& hops, const std::vector<NodeId>& nodes);
  // The DAG's edges in ascending EdgeId, without a copy (kept by the DAG since construction).
  [[nodiscard]] static std::span<const EdgeId> ascending(const PredecessorDAG& dag) noexcept;
  // The sparse storage: the member nodes (ascending) and their entry offsets (size + 1).
  [[nodiscard]] static std::span<const NodeId> member_nodes(const PredecessorDAG& dag) noexcept;
  [[nodiscard]] static std::span<const std::int32_t> member_offsets(const PredecessorDAG& dag) noexcept;
  // The retained payload in bytes (every stored array; the ShortestPathCache's `bytes` statistic).
  [[nodiscard]] static std::size_t payload_bytes(const PredecessorDAG& dag) noexcept;
};

// Internal construction of Path values from trusted walks (the searches' results): the checks of
// Path's constructor hold by construction, so none is repeated (no per-path O(E) scan).
struct PathBuilder {
  // A walk of single-member hops [[e1], [e2], ...].
  [[nodiscard]] static Path of_walk(const Graph& graph, std::span<const EdgeId> walk);
  // Contiguous hops of parallel members (each hop's members share (u, v); no node repeats).
  [[nodiscard]] static Path of_hops(const Graph& graph, std::vector<std::vector<EdgeId>> hops);
};

} // namespace netgraph::core::detail
