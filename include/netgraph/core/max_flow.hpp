/* The one integer minimum-cost flow core (APPROACH 2.2). */
#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// Potential of a vertex that no Dijkstra run of the call settled. Potentials are otherwise
// non-negative path costs, so the sentinel never collides.
inline constexpr Cost kNeverReached = -1;
// Upper bound of the internal value cap; also the capacity of the virtual source and sink arcs,
// strictly above every attainable value (below 2^62), so a virtual arc never saturates.
inline constexpr Amount kMaxFlowCap = INT64_MAX - 1;

// Operation counters of one call: diagnostics for the performance gate's drift check, identical
// across identical calls, never an oracle. Virtual vertices and arcs count like real ones.
struct FlowCoreCounters {
  std::int64_t dijkstra_runs {0};  // every run, the final failed one included
  std::int64_t dijkstra_pops {0};  // heap pops, stale ones included
  std::int64_t settled {0};        // vertices finalised, summed over runs
  std::int64_t arc_scans {0};      // arcs examined by Dijkstra, Dinic BFS and DFS advance
  std::int64_t dinic_phases {0};   // BFS level builds, the last unsuccessful one of a tier included
  std::int64_t augmentations {0};  // DFS paths pushed
};

// The private payload of _place_global_units.
struct FlowCoreResult {
  Amount placed {0};
  // The marginal cost curve: (cost, units) with costs strictly increasing and units positive.
  std::vector<std::pair<Cost, Amount>> cost_distribution {};
  std::vector<Amount> edge_flows {};  // per edge, id order
  // Edges with positive view from the full residual closure of the real sources to its
  // complement, ascending id; nullopt while the destinations stay connected (a cap stopped the
  // solve) and on an input cap of 0.
  std::optional<std::vector<EdgeId>> min_cut {};
  std::vector<Cost> potentials {};  // per node; kNeverReached when never settled
  FlowCoreCounters counters {};
};

// Minimum-cost flow of the maximum value attainable under `cap` from the source set to the
// destination set over the eligible edges, on the capacities `view` (units, id order). Every
// solve starts from zero flow. `cap` nullopt is uncapped; an input cap of 0 returns the zero
// payload without running a search. Throws std::invalid_argument before any work for a view of
// the wrong length or with an entry outside [0, capacity_units[e]], an empty or out-of-range
// endpoint set, a vertex in both sets, a cap outside [0, 2^63 - 2], or a foreign Eligibility.
// Excluded endpoints are skipped.
[[nodiscard]] FlowCoreResult place_global_units(const StrictMultiDiGraph& g,
                                                std::span<const Amount> view,
                                                std::span<const NodeId> sources,
                                                std::span<const NodeId> destinations,
                                                const Eligibility* elig,
                                                std::optional<Amount> cap);

} // namespace netgraph::core
