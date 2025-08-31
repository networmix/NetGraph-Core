#pragma once

#include <optional>
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
  std::vector<std::int32_t> parents;
  std::vector<std::int32_t> via_edges;
};

// Optional node/edge masks:
// - node_mask[v] == true means node v is allowed; false excludes it from search.
// - edge_mask[e] == true means edge e is allowed; false excludes it from search.
// If masks are nullptr, they are ignored.
std::pair<std::vector<double>, PredDAG>
shortest_paths(const StrictMultiDiGraph& g, NodeId src,
               std::optional<NodeId> dst,
               EdgeSelect policy, bool multipath, double eps,
               const bool* node_mask = nullptr,
               const bool* edge_mask = nullptr);

// Residual-aware shortest paths: like shortest_paths but considers per-edge
// residual capacities (residual[e] = remaining capacity). Some EdgeSelect
// policies require residual/capacity/load (e.g., load-factored). When residual
// is provided, edges with residual < MIN_CAP are excluded.
std::pair<std::vector<double>, PredDAG>
shortest_paths_with_residual(const StrictMultiDiGraph& g, NodeId src,
                             std::optional<NodeId> dst,
                             EdgeSelect policy, bool multipath, double eps,
                             const std::vector<double>& residual,
                             const bool* node_mask = nullptr,
                             const bool* edge_mask = nullptr);

} // namespace netgraph::core
