/* The public search: spf from a set of roots, forward or reverse,
   with an optional set of targets (anycast); ksp, the K shortest paths in the canonical order;
   and resolve_to_paths, a DAG unfolded into explicit paths. */
#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// The queue holds nodes waiting for their shortest distance to be finalised (C++ only).
// Auto groups them into cyclic distance buckets if max_edge_cost / max(1, edge_cost_gcd) + 1 <= 65,536;
// edge_cost_gcd is the greatest common divisor of edge costs. Otherwise use a binary heap, which
// keeps the smallest key first; Heap always uses it. Both pop in (distance, node) order and
// give identical results. The bucket limit bounds bucket storage.
enum class SpfQueue { Auto, Heap };

// spf requires a direction. Without targets, finalise every reachable node's shortest distance;
// with targets, stop at the nearest targets (anycast). multipath = true returns every equal-cost
// walk as a DAG, false the one canonical path (targets required); parallel_edges = false keeps one
// member of each (u, v) pair by the tie rule; tie_break orders equal-cost choices (RouteTieBreak).
// A Higher* tie rule with multipath and parallel_edges both true has nothing to choose and is
// rejected. A supplied residual gives unused capacity: int64[E] in EdgeId order, each in
// [0, edge_capacity_units], borrowed for this call. Only edges with residual >= min_residual are
// usable, and these values also supply the availability scores. min_residual must be >= 1 and
// applies only when residual is supplied. The Eligibility always applies. Without a residual no
// edge is skipped for its capacity: an edge of capacity 0 is searched like any other. Invalid
// enum values throw std::invalid_argument naming the field. Fields remain assignable, so spf
// validates them again.
struct SpfOptions {
  explicit SpfOptions(Direction direction_);
  Direction direction;
  std::optional<std::vector<NodeId>> targets {};
  RouteTieBreak tie_break {RouteTieBreak::FewerHops};
  bool multipath {true};
  bool parallel_edges {true};
  std::optional<std::span<const Amount>> residual {};
  Amount min_residual {1};
  const Eligibility* eligibility {nullptr};
};

// The result of spf. costs: full searches only (no targets), one entry per node, kUnreachable
// for the nodes the search did not reach. dag: the unique DAG of the member graph's edges that
// are tight in the traversal direction, returned in their original orientation u -> v (d the
// search's distance: d(u) + cost = d(v) for a forward search, d(v) + cost = d(u) for a reverse
// one); with targets, pruned to the walks between the roots and a nearest target (forward: from a
// root to the target; reverse: from the target to a root), or the canonical single path for
// multipath = false. target_cost: searches with targets only, the nearest targets' cost, or
// kUnreachable when no target is reached.
struct SpfResult {
  std::optional<std::vector<Cost>> costs;
  PredecessorDAG dag;
  std::optional<Cost> target_cost;
};

// Validates in this order before any work (std::invalid_argument): the enumerations (direction,
// tie_break, the queue), the root set (non-empty, ids in range; duplicates removed), the target
// set when given (the same), a vertex in both sets, multipath = false without targets, a Higher*
// tie_break with multipath and parallel_edges both true ("spf tie_break requires multipath=False or
// parallel_edges=False"), the residual (length E, then every entry in [0, edge_capacity_units]),
// min_residual >= 1, a foreign Eligibility. A root or target excluded by the Eligibility is skipped; all roots excluded gives
// an empty result (no error).
[[nodiscard]] SpfResult spf(const Graph& graph, std::span<const NodeId> roots, const SpfOptions& options,
                            SpfQueue queue = SpfQueue::Auto);

// ksp requires k >= 1. parallel_edges = true keeps all eligible parallel edges; false keeps one
// per (u, v) pair by the tie rule. tie_break orders paths of equal cost; the paths are always
// single (k sets their number).
// cost_ceiling is inclusive: AbsoluteCost(c), c >= 0, or CostStretch(f), finite f in [1, inf), times
// the static shortest cost under the Eligibility (edges of capacity 0 included, residual
// ignored). f is compared exactly, as for CostCeiling.
// A residual is validated as for spf and borrowed for the call; edges need at least 1 unused
// unit, and it supplies the availability scores. Without a residual no edge is skipped for its capacity.
// Fields remain assignable, so ksp validates them.
struct KspOptions {
  explicit KspOptions(std::int64_t k_);
  std::int64_t k;
  RouteTieBreak tie_break {RouteTieBreak::FewerHops};
  bool parallel_edges {true};
  std::optional<CostCeiling> cost_ceiling {};
  std::optional<std::span<const Amount>> residual {};
  const Eligibility* eligibility {nullptr};
};

// The K shortest source -> destination paths with no repeated node, ordered by (cost, -bottleneck
// under a Higher* tie rule, hops, edge-id sequence), the bottleneck being the minimum edge score
// along a path: the minimal prefix of that order, truncated at k and at the ceiling. Each
// path is a Path of single-member hops with its cost. Validates in this order before any work
// (std::invalid_argument): tie_break, source and destination in range, source == destination, k >= 1,
// the ceiling, the residual, a foreign Eligibility. An excluded endpoint, no static route under a
// CostStretch ceiling, or no path gives an empty list.
[[nodiscard]] std::vector<std::pair<Cost, Path>> ksp(const Graph& graph, NodeId source, NodeId destination, const KspOptions& options);

// The source -> destination walks of a DAG in increasing (cost, hops, edge-id sequence): a hop's
// cost is its cheapest member's, and a grouped hop is represented in the key by its lowest member
// id. split_parallel_edges = false keeps every parallel member of the DAG's (u, v) pair in one hop
// (members ascending); true gives one path per member combination, sorted by that path's own
// (cost, hops, edge-id sequence). max_paths (>= 1) truncates the sorted list; absent, every walk.
// Validates (std::invalid_argument): the DAG's graph, source and destination in range, source !=
// destination, max_paths >= 1. An endpoint outside the DAG, or no walk, gives an empty list.
[[nodiscard]] std::vector<Path> resolve_to_paths(const Graph& graph, const PredecessorDAG& dag, NodeId source, NodeId destination,
                                                 bool split_parallel_edges = false,
                                                 std::optional<std::int64_t> max_paths = std::nullopt);

} // namespace netgraph::core
