/* Shared traversal for spf, policy routes, cached shortest-path supports and KSP alternatives.
   A support contains every edge on the eligible shortest paths. A KSP spur search finds an
   alternative suffix from a node on an earlier path.

   Dijkstra finalises the nearest pending node, scanning out-edges for Forward and in-edges
   for Reverse. Both heap and cyclic distance-bucket queues pop in (distance, node) order;
   buckets are used when the graph's costs fit 65,536 buckets. A candidate must pass static
   eligibility and any residual-capacity threshold; its capacity is never read (an edge of
   capacity 0 is searched like any other).
   A tight edge extends a shortest path exactly: d(u) + cost(e) = d(v) in traversal direction.
   Positive costs make these edges unique for a fixed member graph regardless of equal-distance
   pop order. Record them in per-node lists; on output, multi_edge = false (Collapse) selects
   one member per (u, v) by select_adjacency, without changing distances.
   With targets, stop when the first target at distance d* is popped. Every nearer node is
   already final, so positive costs ensure every target at d* has its final distance too.
   Keep only tight edges on walks to these nearest targets (anycast pruning).
   Extract the single path from that DAG by bottleneck under PreferHigherPriority (largest
   minimum edge priority), then fewest hops, then smallest EdgeId sequence. Increasing search
   distance orders every predecessor before its child, so no second traversal is needed.
   Each call exclusively borrows the thread's workspace (a lease). Reset it on every exit,
   including exceptions: the touched nodes only, or the whole arrays when a quarter or more of
   the nodes were touched; refuse overlapping leases.
*/
#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "netgraph/core/detail/hooks_flag.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/search.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

// An exact priority value: num / den with den >= 1 (den == 1 for Residual and Capacity, the
// capacity for Headroom on a positive edge); `infinite` above every finite value (the bottleneck
// of an empty walk).
struct Priority {
  Amount num {0};
  Amount den {1};
  bool infinite {false};
};
[[nodiscard]] int compare(const Priority& a, const Priority& b) noexcept;  // -1, 0, 1

// The priority value of an edge with `residual` units left of `capacity`: the residual units
// (Residual), the capacity units (Capacity) or the rational residual / capacity (Headroom; 0 on
// an edge of capacity 0).
[[nodiscard]] Priority priority_value(EdgePriority priority, Amount residual, Amount capacity) noexcept;

// The per-edge priority values a search reads: priority_value on `residual`, or, when it is
// empty (no residual), on the static capacity as the residual (Residual and Capacity read the
// capacity, Headroom 1 on a positive edge and 0 on an edge of capacity 0).
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

// select_adjacency: whether member `e` replaces `cur`, the choice so far among the
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

// The candidates of a search (validated by the caller): the static Eligibility (null:
// everything allowed) and a residual gate (`residual` empty: none; else E entries, an edge needs
// residual >= min_residual >= 1). Capacity is not a condition.
struct Gate {
  const Eligibility* elig {nullptr};
  std::span<const Amount> residual {};
  Amount min_residual {1};
};

// A policy's route and the cost of its walks: its edges in ascending EdgeId (support_route), or
// a single path's in walk order (route_search).
struct SearchRoute {
  std::vector<EdgeId> edges {};
  Cost cost {kUnreachable};
};

// T1's CSPF: forward from `src` toward the destination set (allowed ones; disjoint from src) over
// the candidates of `gate` (min_residual 0: no bandwidth filter), Collapse and the staged single
// path by `selection` (its multipath is not read) with the priority values on `gate.residual`;
// the path's edges in walk order, nullopt when no destination is reached.
[[nodiscard]] std::optional<SearchRoute> route_search(const Graph& g, const Gate& gate, const PathSelection& selection,
                                                      NodeId src, std::span<const NodeId> dests);

// CostOnly's projection of a cached full static support: the support edges
// reachable from `src` (ascending `support`: a reverse multipath DAG rooted at the destination
// set), Collapse by `selection` on the priority values `view` (no residual gating), and for
// multipath = false the staged canonical single path on that part, its distances derived from the
// support itself. An empty route when src reaches no destination of `dests` in the support.
[[nodiscard]] SearchRoute support_route(const Graph& g, std::span<const EdgeId> support, NodeId src,
                                        std::span<const NodeId> dests, const PathSelection& selection, PriorityFn view);

// The full static support of a search from (Forward) or toward (Reverse) the root set (a
// normalised set: nonempty, in range, ascending, no duplicates; excluded roots skipped): every
// tight edge of the statically eligible member graph (multipath, multi_edge, no residual) over the
// nodes it reaches, forward-oriented (the DagCache's value, and a CostOnly policy's DAG when it has
// no cache). The empty DAG when no root is allowed.
[[nodiscard]] PredDAG static_support(const Graph& g, std::span<const NodeId> roots, Direction direction,
                                     const Eligibility* elig);

// The static cost from each of `sources` to the nearest of `dests` (in range, ascending, no
// duplicates, possibly empty; allowed ones) over the statically eligible edges, whatever their
// capacity, kUnreachable where there is no walk, and everywhere when no destination is allowed
// (the T1 and Global rosters and T1's Stretch baseline).
[[nodiscard]] std::vector<Cost> static_costs(const Graph& g, const Eligibility* elig, std::span<const NodeId> dests,
                                             std::span<const NodeId> sources);

// The validation of a cost ceiling (KSP and T1): Absolute's cost >= 0,
// Stretch's factor finite and in [1, inf); std::invalid_argument naming max_path_cost.
void check_ceiling(const std::optional<CostCeiling>& ceiling);

#if NETGRAPH_CORE_PRIVATE_HOOKS
// Test evidence (C++ tests, hooks builds only): how many of this thread's searches ran on the
// bucket queue; a search that leaves the count unchanged ran on the heap.
struct QueueCounts {
  std::int64_t buckets {0};
};
QueueCounts& thread_queue_counts() noexcept;
#endif

} // namespace netgraph::core::detail
