/* The public search (APPROACH 5 and 6): spf from a set of roots, forward or reverse, with an
   optional set of targets (anycast); ksp, the K shortest paths in the canonical order; and
   resolve_to_paths, a DAG unfolded into explicit paths. */
#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// The frontier queue of a search (C++ only; the results never depend on it). Auto takes the
// monotone bucket queue when the graph's costs fit max_cost / max(1, cost_gcd) + 1 <= 65,536
// buckets, else the binary heap; Heap always takes the heap. Both pop in (distance, node) order
// and give identical results.
enum class SpfQueue { Auto, Heap };

// The options of spf (APPROACH 6). direction and selection are explicit; targets absent is a
// full search (every reachable node settled), present a search that stops at the nearest
// targets (anycast). A supplied residual (int64[E] in EdgeId order, each entry in
// [0, capacity_units]) is borrowed for the call; it adds residual gating at min_residual units
// and supplies the residual values of the selected priority. min_residual (an integer >= 1) gates
// a supplied residual only: static eligibility (the Eligibility, capacity >= 1 unit) is
// unconditional. The constructor rejects a direction outside its enumeration
// (std::invalid_argument naming the field); the fields stay assignable, so spf checks them again.
struct SpfOptions {
  SpfOptions(Direction direction_, PathSelection selection_);
  Direction direction;
  std::optional<std::vector<NodeId>> targets {};
  PathSelection selection;
  std::optional<std::span<const Amount>> residual {};
  Amount min_residual {1};
  const Eligibility* elig {nullptr};
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
  PredDAG dag;
  std::optional<Cost> target_cost;
};

// Validates in this order before any work (std::invalid_argument): the enumerations (direction,
// the selection, the queue), the root set (non-empty, ids in range; duplicates removed), the
// target set when given (the same), a vertex in both sets, multipath = false without targets,
// the residual (length E, then every entry in [0, capacity_units]), min_residual >= 1, a foreign
// Eligibility. A root or target excluded by the Eligibility is skipped; all roots excluded gives
// an empty result (no error).
[[nodiscard]] SpfResult spf(const Graph& g, std::span<const NodeId> roots, const SpfOptions& opts,
                            SpfQueue queue = SpfQueue::Auto);

// The options of ksp (APPROACH 6, rule 10). k (>= 1) and selection are explicit; the selection's
// multi_edge is the member rule (Expand keeps every eligible parallel member, Collapse one per
// (u, v) pair by select_adjacency) and multipath is ignored. max_path_cost: an inclusive ceiling,
// Absolute(c) (c >= 0) or Stretch(f) (f finite in [1, inf), read as the exact dyadic rational it
// is) times the static shortest cost of the member graph under the Eligibility, residual ignored.
// A supplied residual (validated as for spf, borrowed for the call) gates at 1 unit and supplies
// the priority values. The fields stay assignable; ksp checks them.
struct KspOptions {
  KspOptions(std::int64_t k_, PathSelection selection_);
  std::int64_t k;
  std::optional<CostCeiling> max_path_cost {};
  PathSelection selection;
  std::optional<std::span<const Amount>> residual {};
  const Eligibility* elig {nullptr};
};

// The K shortest simple src -> dst paths of the member graph in increasing (cost, -bottleneck
// under PreferHigherPriority, hops, edge-id sequence), the bottleneck being the minimum priority
// value along a path: the minimal prefix of that order, truncated at k and at the ceiling. Each
// path is a Path of single-member hops with its cost. Validates in this order before any work
// (std::invalid_argument): the selection's enumerations, src and dst in range, src == dst, k >= 1,
// the ceiling, the residual, a foreign Eligibility. An excluded endpoint, no static route under a
// Stretch ceiling, or no path gives an empty list.
[[nodiscard]] std::vector<std::pair<Cost, Path>> ksp(const Graph& g, NodeId src, NodeId dst, const KspOptions& opts);

// The src -> dst walks of a DAG (APPROACH 5) in increasing (cost, hops, edge-id sequence): a hop's
// cost is its cheapest member's, and a grouped hop is represented in the key by its lowest member
// id. split_parallel_edges = false keeps every parallel member of the DAG's (u, v) pair in one hop
// (members ascending); true gives one path per member combination, in the same order. max_paths
// (>= 1) truncates the sorted list; absent, every walk. Validates (std::invalid_argument): a DAG of
// another graph, src and dst in range, src == dst, max_paths >= 1. An endpoint outside the DAG, or
// no walk, gives an empty list.
[[nodiscard]] std::vector<Path> resolve_to_paths(const Graph& g, const PredDAG& dag, NodeId src, NodeId dst,
                                                 bool split_parallel_edges = false,
                                                 std::optional<std::int64_t> max_paths = std::nullopt);

} // namespace netgraph::core
