/* Shortest paths (Dijkstra) with multipath predecessor DAG support. */
#pragma once

#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// Predecessor DAG: compact representation of all equal-cost predecessors.
// For each node v, entries are stored in [parent_offsets[v], parent_offsets[v+1])
// as pairs (parents[i], via_edges[i]) where via_edges[i] is the compacted EdgeId
// used to reach v from parents[i]. Multiple parallel edges are represented by
// multiple entries with the same parent.
struct PredDAG {
  std::vector<std::int32_t> parent_offsets;
  std::vector<NodeId> parents;
  std::vector<EdgeId> via_edges;
};

// Optional node/edge masks:
// - node_mask[v] == true means node v is allowed; false excludes it from search.
// - edge_mask[e] == true means edge e is allowed; false excludes it from search.
// If masks are nullptr, they are ignored.
std::pair<std::vector<Cost>, PredDAG>
shortest_paths(const StrictMultiDiGraph& g, NodeId src,
               std::optional<NodeId> dst,
               const EdgeSelection& selection,
               std::span<const Cap> residual = {},
               const bool* node_mask = nullptr,
               const bool* edge_mask = nullptr);

} // namespace netgraph::core
