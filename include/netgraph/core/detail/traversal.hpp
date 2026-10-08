/* The one traversal core of every search (internal; APPROACH 6 and 14): public spf, the
   policies' route searches, the DagCache's supports and KSP's spur searches.

   One Dijkstra loop, parametrised by the direction (Forward scans a node's out-edge id range,
   Reverse its in-edge slots) and the queue (a binary heap, or monotone cyclic buckets when the
   graph's costs fit 65,536 buckets), both popping in (distance, node) order. Positive costs make
   the set of tight edges d(u) + cost(e) = d(v) of a fixed member graph unique (2.4), whatever the
   pop order among equal distances; the loop records them in per-node lists while it runs. An edge
   is a candidate when it passes static eligibility, capacity >= 1 unit and the residual gate;
   Collapse (multi_edge = false) keeps one member per (u, v) pair by select_adjacency when the
   tight edges are emitted (it never changes a distance). With targets the loop stops at the first
   target popped: every node nearer than that distance d* is settled by then, so every target at
   d* has its final distance, and the walks back from the targets at d* are the anycast-pruned DAG.
   The canonical single path is extracted from that DAG in stages (section 6) with the search's
   own distances as the topological order: no second traversal.

   Every call leases the thread's workspace (node arrays reset through touched lists); the lease
   restores it on every exit, exceptions included, and overlapping leases are refused.
*/
#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/search.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

// An exact priority value: num / den with den >= 1 (den == 1 for Residual and Capacity, the
// capacity for Headroom); `infinite` above every finite value (the bottleneck of an empty walk).
struct Priority {
  Amount num {0};
  Amount den {1};
  bool infinite {false};
};
[[nodiscard]] int compare(const Priority& a, const Priority& b) noexcept;  // -1, 0, 1

// The priority value of an edge with `residual` units left of `capacity`: the residual units
// (Residual), the capacity units (Capacity) or the rational residual / capacity (Headroom).
[[nodiscard]] Priority priority_value(EdgePriority priority, Amount residual, Amount capacity) noexcept;

// The per-edge priority values a search reads: `residual` empty (no residual) reads the static
// capacity for every priority (APPROACH 10); else priority_value on the residual.
struct PriorityValues {
  const Graph* g {nullptr};
  EdgePriority priority {EdgePriority::Capacity};
  std::span<const Amount> residual {};
  [[nodiscard]] Priority of(EdgeId e) const noexcept;
};

// Any source of per-edge priority values (a search's PriorityValues, a policy's working view),
// without a virtual call or an allocation.
struct PriorityFn {
  const void* ctx {nullptr};
  Priority (*fn)(const void*, EdgeId) noexcept {nullptr};
  [[nodiscard]] Priority operator()(EdgeId e) const noexcept { return fn(ctx, e); }
  [[nodiscard]] static PriorityFn of(const PriorityValues& values) noexcept;
};

// select_adjacency (APPROACH 6): whether member `e` replaces `cur`, the choice so far among the
// eligible members of one (u, v) pair seen in ascending EdgeId (cur < e): the lower cost, then
// under PreferHigherPriority (`prefer`) the higher priority value; on a tie the lower EdgeId
// (cur) stays. `value(e)` reads the priority value of an edge on the view of the selection.
template <class Value>
[[nodiscard]] bool replaces(const Graph& g, EdgeId e, EdgeId cur, bool prefer, const Value& value) {
  const auto cost = g.cost();
  const Cost ce = cost[static_cast<std::size_t>(e)];
  const Cost cc = cost[static_cast<std::size_t>(cur)];
  if (ce != cc) return ce < cc;
  return prefer && compare(value(e), value(cur)) > 0;
}

// ---------------------------------------------------------------- queue dispatch

// The bucket queue's limit (an engineering cap, not a proven crossover): Auto takes the buckets
// when W = max_cost / max(1, cost_gcd) + 1 <= kMaxBucketWidth.
inline constexpr std::int64_t kMaxBucketWidth = std::int64_t {1} << 16;
// W of g (1 on a graph without edges: max_cost and cost_gcd are 0 there).
[[nodiscard]] std::int64_t bucket_width(const Graph& g) noexcept;
// Whether a search on g with `queue` runs on the buckets (Heap never does).
[[nodiscard]] bool uses_buckets(const Graph& g, SpfQueue queue) noexcept;

// ---------------------------------------------------------------- searches of validated inputs

// The candidates of a search beyond capacity >= 1 unit (validated by the caller): the static
// Eligibility (null: everything allowed) and a residual gate (`residual` empty: none; else E
// entries, an edge needs residual >= min_residual >= 1).
struct Gate {
  const Eligibility* elig {nullptr};
  std::span<const Amount> residual {};
  Amount min_residual {1};
};

// A policy's route: the edges in ascending EdgeId and the cost of its walks.
struct SearchRoute {
  std::vector<EdgeId> edges {};
  Cost cost {kUnreachable};
};

// CapacityAware's search (APPROACH 9.1): forward from `src` toward the destination set (allowed
// ones; disjoint from src) over the candidates of `gate`, Collapse and the staged single path by
// `selection` with the priority values on `gate.residual`; nullopt when no destination is reached.
[[nodiscard]] std::optional<SearchRoute> route_search(const Graph& g, const Gate& gate, const PathSelection& selection,
                                                      NodeId src, std::span<const NodeId> dests);

// CostOnly's projection of a cached full static support (APPROACH 9.1, 10): the support edges
// reachable from `src` (ascending `support`: a reverse multipath DAG rooted at the destination
// set), Collapse by `selection` on the priority values `view` (no residual gating), and for
// multipath = false the staged canonical single path on that part, its distances derived from the
// support itself. An empty route when src reaches no destination of `dests` in the support.
[[nodiscard]] SearchRoute support_route(const Graph& g, std::span<const EdgeId> support, NodeId src,
                                        std::span<const NodeId> dests, const PathSelection& selection, PriorityFn view);

// The full static support of a search from (Forward) or toward (Reverse) the root set (every root
// in range; excluded roots skipped): every tight edge of the statically eligible member graph
// (multipath, multi_edge, no residual) over the nodes it reaches, forward-oriented (the DagCache's
// value, and a CostOnly policy's DAG when it has no cache).
[[nodiscard]] PredDAG static_support(const Graph& g, std::span<const NodeId> roots, Direction direction,
                                     const Eligibility* elig);

// The static cost from each of `sources` to the nearest of `dests` (allowed ones) over the
// statically eligible edges of at least 1 unit, kUnreachable where there is no walk (the
// CapacityAware and Global rosters and the Stretch baselines).
[[nodiscard]] std::vector<Cost> static_costs(const Graph& g, const Eligibility* elig, std::span<const NodeId> dests,
                                             std::span<const NodeId> sources);

// The rule-13 check of a cost ceiling (APPROACH 6; KSP and CapacityAware): Absolute's cost >= 0,
// Stretch's factor finite and in [1, inf); std::invalid_argument naming max_path_cost.
void check_ceiling(const std::optional<CostCeiling>& ceiling);

#if NETGRAPH_CORE_PRIVATE_HOOKS
// Test evidence (C++ tests, hooks builds only): the queue of every search this thread ran.
struct QueueCounts {
  std::int64_t heap {0};
  std::int64_t buckets {0};
};
QueueCounts& thread_queue_counts() noexcept;
#endif

} // namespace netgraph::core::detail
