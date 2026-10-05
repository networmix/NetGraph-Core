/* The flow core restricted to an edge subset (internal; APPROACH 2.2 and 7). */
#pragma once

#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

struct SubsetFlow {
  Amount placed {0};
  std::vector<Amount> flows {};  // aligned with the input edges
  std::vector<std::pair<Cost, Amount>> cost_distribution {};
};

// The minimum-cost flow of the maximum value under `cap` (nullopt: uncapped) from `sources` to
// `destinations` over the edges `edges` (ascending EdgeId, distinct) with capacities `caps`
// (aligned with `edges`, each in [0, 2^62)). The one flow core of APPROACH 2.2 with the same
// arc order (ascending EdgeId at every node, the virtual arc first, virtual arcs by ascending
// endpoint NodeId), from zero flow; the vertices are compressed to the ones the edges and
// endpoints touch, which changes no flow, curve or augmenting path. The endpoint sets are
// ascending, distinct and disjoint; the caller validated every input (the MaxFlow kernel's
// replacement solve on V_i = W + own_i over the flow's split graph).
[[nodiscard]] SubsetFlow solve_subset(const StrictMultiDiGraph& g, std::span<const EdgeId> edges,
                                      std::span<const Amount> caps, std::span<const NodeId> sources,
                                      std::span<const NodeId> destinations, std::optional<Amount> cap);

} // namespace netgraph::core::detail
