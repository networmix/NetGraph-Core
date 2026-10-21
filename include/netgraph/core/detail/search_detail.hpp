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
   pop order. Record them in per-node lists; on output, parallel_edges = false (Collapse) selects
   one member per (u, v) by replaces(), without changing distances.
   With targets, stop when the first target at distance d* is popped. Every nearer node is
   already final, so positive costs ensure every target at d* has its final distance too.
   Keep only tight edges on walks to these nearest targets (anycast pruning).
   Extract the single path from that DAG by bottleneck under a Higher* tie rule (largest
   minimum edge score), then fewest hops, then smallest EdgeId sequence. Increasing search
   distance orders every predecessor before its child, so no second traversal is needed.
   Each call exclusively borrows the thread's workspace (a lease). Reset it on every exit,
   including exceptions: the touched nodes only, or the whole arrays when a quarter or more of
   the nodes were touched; refuse overlapping leases.
*/
#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "netgraph/core/detail/hooks_flag.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/search.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

#if NETGRAPH_CORE_PRIVATE_HOOKS
// Native-only full-Yen comparator and observations. The caller owns the probe throughout
// the call; observers run between workspace leases and may throw to test cleanup.
struct KspProbe {
  bool full_yen {false};
  void (*spur)(void*, std::size_t deviation, std::size_t index) {nullptr};
  void (*duplicate)(void*, std::size_t old_deviation, std::size_t new_deviation) {nullptr};
  void* context {nullptr};
};
[[nodiscard]] const KspProbe* current_ksp_probe() noexcept;
void set_ksp_probe(const KspProbe* probe) noexcept;
struct ScopedKspProbe {
  const KspProbe* previous;
  explicit ScopedKspProbe(const KspProbe* probe) noexcept : previous(current_ksp_probe()) { set_ksp_probe(probe); }
  ~ScopedKspProbe() { set_ksp_probe(previous); }
  ScopedKspProbe(const ScopedKspProbe&) = delete;
  ScopedKspProbe& operator=(const ScopedKspProbe&) = delete;
};
#endif

// A route choice as the searches read it: the tie rule, all equal-cost walks (multipath) or the
// canonical single path, and whether every parallel edge of a (u, v) pair is kept (parallel_edges).
// The public options (SpfOptions, KspOptions, ShortestPaths, DistributedTE) each give it its fields:
// KspOptions and DistributedTE always order single paths.
struct RouteChoice {
  RouteTieBreak tie_break {RouteTieBreak::FewerHops};
  bool multipath {true};
  bool parallel_edges {true};
  // A Higher* tie rule: the edge score breaks cost ties.
  [[nodiscard]] bool prefers_score() const noexcept { return tie_break != RouteTieBreak::FewerHops; }
  // The ineffective tie rule (ShortestPaths and spf): a Higher* tie rule over every equal-cost walk
  // and every parallel edge has nothing to choose.
  [[nodiscard]] bool tie_is_ineffective() const noexcept { return prefers_score() && multipath && parallel_edges; }
  [[nodiscard]] bool single_path() const noexcept { return !multipath; }
  friend bool operator==(const RouteChoice&, const RouteChoice&) = default;
};

// An exact edge score: num / den with den >= 1 (den == 1 for the capacity and the availability,
// the capacity for the headroom on a positive edge); `infinite` above every finite value (the
// bottleneck of an empty walk).
struct Priority {
  Amount num {0};
  Amount den {1};
  bool infinite {false};
};
[[nodiscard]] int compare(const Priority& a, const Priority& b) noexcept;  // -1, 0, 1

// The score of an edge with `residual` units available of `capacity` under a Higher* tie rule: the
// capacity units (HigherCapacity), the residual units (HigherAvailableCapacity) or the rational
// residual / capacity (HigherHeadroom; 0 on an edge of capacity 0). FewerHops reads no score (the
// capacity is returned and never compared).
[[nodiscard]] Priority priority_value(RouteTieBreak tie_break, Amount residual, Amount capacity) noexcept;

// The per-edge scores a search reads: priority_value on `residual`, or, when it is empty (no
// residual), on the static capacity as the residual (HigherCapacity and HigherAvailableCapacity read
// the capacity, HigherHeadroom 1 on a positive edge and 0 on an edge of capacity 0).
struct PriorityValues {
  const Graph* graph {nullptr};
  RouteTieBreak tie_break {RouteTieBreak::FewerHops};
  std::span<const Amount> residual {};
  [[nodiscard]] Priority of(EdgeId e) const noexcept;
};

// Any source of per-edge priority values (a search's PriorityValues, a policy's working view),
// without a virtual call or an allocation.
struct PriorityReader {
  const void* context {nullptr};
  Priority (*function)(const void*, EdgeId) noexcept {nullptr};
  [[nodiscard]] Priority operator()(EdgeId e) const noexcept { return function(context, e); }
  [[nodiscard]] static PriorityReader of(const PriorityValues& values) noexcept;
};

// replaces: whether member `e` replaces `cur`, the choice so far among the
// eligible members of one (u, v) pair seen in ascending EdgeId (cur < e): the lower cost, then
// under a Higher* tie rule (`prefer`) the higher score; on a tie the lower EdgeId (cur) stays.
// `value(e)` reads the score of an edge on the view of the route choice.
template <class Value>
[[nodiscard]] bool replaces(const Graph& graph, EdgeId e, EdgeId cur, bool prefer, const Value& value) {
  const auto cost = graph.edge_costs();
  const Cost ce = cost[static_cast<std::size_t>(e)];
  const Cost cc = cost[static_cast<std::size_t>(cur)];
  if (ce != cc) return ce < cc;
  return prefer && compare(value(e), value(cur)) > 0;
}

// ---------------------------------------------------------------- queue dispatch

// The bucket queue's limit (an engineering cap, not a proven crossover): Auto takes the buckets
// when W = max_edge_cost / max(1, edge_cost_gcd) + 1 <= kMaxBucketWidth.
inline constexpr std::int64_t kMaxBucketWidth = std::int64_t {1} << 16;
// W of graph (1 on a graph without edges: max_edge_cost and edge_cost_gcd are 0 there).
[[nodiscard]] std::int64_t bucket_width(const Graph& graph) noexcept;
// Whether a search on graph with `queue` runs on the buckets (Heap never does).
[[nodiscard]] bool uses_buckets(const Graph& graph, SpfQueue queue) noexcept;

// ---------------------------------------------------------------- searches of validated inputs

// The candidates of a search (validated by the caller): the static Eligibility (null:
// everything allowed) and a residual gate (`residual` empty: none; else E entries, an edge needs
// residual >= min_residual >= 1). Capacity is not a condition.
struct SearchFilter {
  const Eligibility* eligibility {nullptr};
  std::span<const Amount> residual {};
  Amount min_residual {1};
};

// A policy's route and the cost of its walks: its edges in ascending EdgeId (support_route), or
// a single path's in walk order (route_search).
struct SearchRoute {
  std::vector<EdgeId> edges {};
  Cost cost {kUnreachable};
};

// T1's CSPF: forward from `source` toward the destination set (allowed ones; disjoint from it) over
// the candidates of `gate` (min_residual 0: no bandwidth filter), Collapse and the staged single
// path by `choice` (its multipath is not read) with the scores on `gate.residual`;
// the path's edges in walk order, nullopt when no destination is reached.
[[nodiscard]] std::optional<SearchRoute> route_search(const Graph& graph, const SearchFilter& gate, const RouteChoice& choice,
                                                      NodeId source, std::span<const NodeId> destinations);

// ShortestPaths's projection of a cached full static support: the support edges
// reachable from `source` (ascending `support`: a reverse multipath DAG rooted at the destination
// set), Collapse by `choice` on the scores `view` (no residual gating), and for
// multipath = false the staged canonical single path on that part, its distances derived from the
// support itself. An empty route when `source` reaches no destination in the support.
[[nodiscard]] SearchRoute support_route(const Graph& graph, std::span<const EdgeId> support, NodeId source,
                                        std::span<const NodeId> destinations, const RouteChoice& choice, PriorityReader view);

// The full static support of a search from (Forward) or toward (Reverse) the root set (a
// normalised set: nonempty, in range, ascending, no duplicates; excluded roots skipped): every
// tight edge of the statically eligible member graph (multipath, parallel edges, no residual) over the
// nodes it reaches, forward-oriented (the ShortestPathCache's value, and a ShortestPaths policy's DAG when it has
// no cache). The empty DAG when no root is allowed.
[[nodiscard]] PredecessorDAG static_support(const Graph& graph, std::span<const NodeId> roots, Direction direction,
                                     const Eligibility* eligibility);

// The static cost from each of `sources` to the nearest of `destinations` (in range, ascending, no
// duplicates, possibly empty; allowed ones) over the statically eligible edges, whatever their
// capacity, kUnreachable where there is no walk, and everywhere when no destination is allowed
// (the T1 and CentralizedTE rosters and T1's CostStretch baseline).
[[nodiscard]] std::vector<Cost> static_costs(const Graph& graph, const Eligibility* eligibility, std::span<const NodeId> destinations,
                                             std::span<const NodeId> sources);

// The validation of a cost ceiling at use (ksp and T1): an AbsoluteCost's cost >= 0, a
// CostStretch's factor >= 1; std::invalid_argument naming cost_ceiling.
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
