/* Core type aliases and the bounds of the units domain. */
#pragma once

#include <cstdint>
#include <limits>

namespace netgraph::core {

// Node and edge identifiers are signed 32-bit integers. An EdgeId is the rank of the edge in
// the stable order (src, cost, dst, input) of its graph; the input position of an
// edge is Graph::input_index()[e].
using NodeId = std::int32_t;
using EdgeId = std::int32_t;
// Costs: integers, never converted. An edge cost is at least 1, so a path's cost is positive; a
// search's distance is 0 at its roots, and kUnreachable (below) marks no path.
using Cost = std::int64_t;
// An amount: an integer number of units of the graph's resolution. Capacities, views, flows and
// requests are amounts; there is no finer grid.
using Amount = std::int64_t;

// Exclusive bound of the units domain: every capacity, the sum of all capacities and the sum of
// all edge costs stay below 2^62.
inline constexpr std::int64_t kDomainLimit = std::int64_t{1} << 62;

// The path-cost sentinel of an unreachable node or an empty route. Python keeps this value itself
// in the integer distance array SpfResult.costs and shows a single unreachable cost
// (SpfResult.target_cost, FlowRecord.route_cost) as None. Costs are positive and their sum is
// below 2^62, so no path cost reaches it.
inline constexpr Cost kUnreachable = std::numeric_limits<Cost>::max();

} // namespace netgraph::core
