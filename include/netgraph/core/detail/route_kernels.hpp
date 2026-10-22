/* Placement kernels assign the amounts of a plan's injections to route edges. Apportionment
   splits an integer total among recipients.

   A plan freezes the route, its split graph and any Weighted weights at creation. It has one or
   more roots: the sources that inject (one for a flow of a single source; a ShortestPaths demand's
   roster for its pooled program). The split graph keeps edges allowed by the kernel that lie on
   a walk from a root to a destination; it removes edges leaving destinations. Its edges are the
   plan's members. A router reached by several roots adds their traffic and splits it once.
   A member-aligned array has one entry per member, in the same order as plan.forwarding_edges. Policies
   combine these allocations; the ledger commits them.
*/
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/detail/search_detail.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

// ---------------------------------------------------------------- apportionment

// Round robin: give floor(total / k) to each of k recipients, plus one to the first
// total mod k in caller order. Edge shares use ascending EdgeId; shares among a source's
// flows use flow order. Used by Equal (Lossless and Lossy) at every node and to split a
// source's amount over its flows.
// Visit only the n positive shares (n = min(total, k)); later entries stay 0 in a zeroed output.
// The callback form computes the split once and calls share(i, x) in index order for i in
// [0, n). Stop if it returns false; return the accepted count: n, or that rejected index i.
// The span form writes out[0 .. n), where k = out.size(), and returns n. round_robin_share is the
// rule itself: the share of recipient i in [0, k) (T1's head designs and tunnel requests use it).
constexpr Amount round_robin_share(Amount total, Amount k, Amount i) noexcept {
  return total / k + (i < total % k ? 1 : 0);
}

template <class Share>
inline std::size_t round_robin(Amount total, std::size_t k, Share&& share) {
  const auto kk = static_cast<Amount>(k);
  const std::size_t n = total / kk == 0 ? static_cast<std::size_t>(total % kk) : k;
  for (std::size_t i = 0; i < n; ++i) {
    if (!share(i, round_robin_share(total, kk, static_cast<Amount>(i)))) return i;
  }
  return n;
}

// The round-robin increments when a cumulative total X (any 128-bit count) grows by v over k =
// out.size() recipients: out[i] = RR(X + v, k, i) - RR(X, k, i). They depend only on r = X mod k:
// floor((r + v) / k) + [i < (r + v) mod k] - [i < r] (ShortestPaths' Lossy offers and a program's
// equal sources).
void round_robin_increments(const Counter128& total, Amount v, std::span<Amount> out) noexcept;

inline std::size_t round_robin(Amount total, std::span<Amount> out) noexcept {
  return round_robin(total, out.size(), [&](std::size_t i, Amount x) noexcept {
    out[i] = x;
    return true;
  });
}

// Highest averages H(total, w) (Jefferson/D'Hondt; Weighted splits and the attribution of
// pooled traffic to sources): give each next unit to the member with largest w / (a + 1), where w
// is its weight and a its allocation so far. Break ties by lower member index. Start each member
// at its lower quota floor(total * w / W), where W is the weight sum; then assign the fewer than k
// remaining units by those quotients, comparing with exact 128-bit cross products. Weights are
// in [0, 2^62) (their sum is exact, any size), total < 2^63: a zero weight gets nothing, and all-zero weights give
// nothing (out all 0). Equal positive weights give round robin. Increasing total never reduces
// any member's allocation (house monotonicity).
void highest_averages(Amount total, std::span<const Amount> weights, std::span<Amount> out);

// The same rule for cumulative offered totals below 2^127. Weight sums stay
// below 2^93 (fewer than 2^31 members). Quotas/priorities use exact 192-bit products.
struct U128;
void highest_averages(U128 total, std::span<const Amount> weights, std::span<U128> out);

// ---------------------------------------------------------------- plans

// Throws std::invalid_argument for a value outside its enumeration and for the unsupported
// combination (Lossy with Optimized: "OPTIMIZED requires admission=LOSSLESS").
void check_kernel(RouteAllocation allocation, Admission admission);
// The combination rule of check_kernel alone, for a caller that has already checked the
// enumerations (check_config).
void check_kernel_combination(RouteAllocation allocation, Admission admission);

// Where a Weighted plan's member weights come from (ignored by the other kernels): Capacity, the
// static capacity; Explicit, `values` per EdgeId (one entry per edge; a positive weight on an edge
// of capacity 0 is 0); ResidualSnapshot, `values` aligned with the route's edges in ascending
// EdgeId (the free capacity read on the creation view). A route edge of weight 0 is no member.
struct WeightSource {
  WeightBasis basis {WeightBasis::Capacity};
  std::span<const Amount> values {};
};

struct PreparedRoute {
  RouteAllocation allocation {RouteAllocation::Optimized};
  std::vector<EdgeId> support_edges {};      // ascending: the plan's route (a single flow's record route)
  std::vector<EdgeId> forwarding_edges {};    // ascending: the split graph
  std::vector<Amount> weights {};    // Weighted: per member (positive), frozen; else empty
  // Destinations ascend. Nodes have a topological order: every member's tail precedes its head
  // (a single root comes first). Position p's out-members are members[outgoing_begin[p] ..
  // outgoing_end[p]); they are consecutive in ascending EdgeId because EdgeIds sort by source first.
  // Each member also stores its head's position. Empty split graph means empty structure. The
  // roots (ascending: the sources that inject) and per root its position, or -1 when it has no
  // member. All views share one block owned by the plan (layout()), so the structure takes one
  // allocation whatever the number of roots.
  std::span<const NodeId> sources {};
  std::span<const std::int32_t> source_positions {};
  std::span<const NodeId> destinations {};
  std::span<const NodeId> topological_order {};
  std::span<const std::int32_t> outgoing_begin {};
  std::span<const std::int32_t> outgoing_end {};
  std::span<const std::int32_t> edge_destination_positions {};
  std::span<const std::uint8_t> is_destination {};

  PreparedRoute() = default;
  PreparedRoute(const PreparedRoute& other);
  PreparedRoute(PreparedRoute&& other) noexcept;
  PreparedRoute& operator=(const PreparedRoute& other);
  PreparedRoute& operator=(PreparedRoute&& other) noexcept;
  ~PreparedRoute() = default;

  // The writable views of a new block: the roots and the destinations copied in, every root's
  // position -1, the structure for `positions` nodes and `num_members` members zeroed. One
  // allocation, made before anything changes (a failure leaves the plan as it was); the previous
  // block is released.
  struct Layout {
    std::span<std::int32_t> source_positions;
    std::span<NodeId> topological_order;
    std::span<std::int32_t> outgoing_begin;
    std::span<std::int32_t> outgoing_end;
    std::span<std::int32_t> edge_destination_positions;
    std::span<std::uint8_t> is_destination;
  };
  Layout layout(std::span<const NodeId> roots, std::span<const NodeId> destinations, std::size_t positions,
                std::size_t num_members);

private:
  // Points the views into store_ for these sizes (the layout's order: sources, source_positions,
  // destinations, topological_order, outgoing_begin, outgoing_end, edge_destination_positions,
  // then is_destination as bytes).
  Layout point(std::size_t num_roots, std::size_t num_dsts, std::size_t positions, std::size_t num_members) noexcept;
  std::vector<std::int32_t> store_ {};
};

// The edges of `dag` that `eligibility` allows (every one when null) and that lie on some walk from
// `source` to the destination set (edges leaving a destination removed: destinations are terminal),
// ascending.
[[nodiscard]] std::vector<EdgeId> prune(const Graph& graph, const PredecessorDAG& dag, NodeId source,
                                        std::span<const NodeId> destinations, const Eligibility* eligibility = nullptr);

// The minimum walk cost from `source` to the destination set over `route` (a record's route:
// its min_route_cost); kUnreachable when the route has no such walk. The split-graph builder's
// cost pass alone: no plan is built.
[[nodiscard]] Cost min_route_cost(const Graph& graph, const PredecessorDAG& route, NodeId source,
                              std::span<const NodeId> destinations);

// The plan of a new flow on `route` (acyclic; the flow's DAG or path, its record's route), its
// Weighted weights from `weights` (a WeightSource; its creation view aligned with the route's edges
// in ascending EdgeId): Optimized and Equal every route edge (an edge of capacity 0 is a
// member that carries nothing: Equal still gives it its round-robin share), Weighted the
// edges of positive weight, with those weights; Lossy keeps its kernel's membership. The members
// are pruned to the walks from `source` to the destination set (none leaving a destination) and
// given their topological structure, in one pass. With `cost`, the route's minimum walk cost
// (min_route_cost) is computed in the same pass when every route edge is a member. `destinations`
// ascending and disjoint from `source`.
[[nodiscard]] PreparedRoute prepare_route(const Graph& graph, RouteAllocation allocation, Admission admission,
                                   const WeightSource& weights, NodeId source, std::span<const NodeId> destinations,
                                   const PredecessorDAG& route, Cost* cost = nullptr);
// The same plan with several roots (a ShortestPaths demand's pooled program): the members are pruned to
// the walks from any root to the destination set; `roots` ascending and disjoint from `destinations`.
[[nodiscard]] PreparedRoute prepare_route(const Graph& graph, RouteAllocation allocation, Admission admission,
                                   const WeightSource& weights, std::span<const NodeId> roots,
                                   std::span<const NodeId> destinations, const PredecessorDAG& route);

// The plan of a supplied split route (from_routes): `dag` pruned to the walks from `source` to the
// destination set (the plan's route, ascending), its members every route edge (Equal) or
// the route edges of positive weight (Weighted: `basis` with `edge_values` one entry per edge, the
// explicit weights or the residual snapshot; a weight on an edge of capacity 0 is 0 unless the
// basis is the snapshot), with `cost` the route's minimum walk cost. One pass over `dag`, and a
// second only when a Weighted route edge has weight 0. Equal to prepare_route on the pruned route.
[[nodiscard]] PreparedRoute prepare_supplied_route(const Graph& graph, RouteAllocation allocation, WeightBasis basis,
                                         std::span<const Amount> edge_values, NodeId source,
                                         std::span<const NodeId> destinations, const PredecessorDAG& dag, Cost* cost);

// A plan's replay under a mask: its members that `eligibility` allows (every one when null),
// pruned to the walks from `roots` (ascending) to `destinations` (ascending; the allowed ones), with the
// members' frozen weights; members of capacity 0 stay (an Equal member keeps its share).
// The route is the plan's.
[[nodiscard]] PreparedRoute restrict_prepared_route(const Graph& graph, const PreparedRoute& plan, std::span<const NodeId> roots,
                                       std::span<const NodeId> destinations, const Eligibility* eligibility);

// ---------------------------------------------------------------- ShortestPaths creation

// std::lower_bound(a, x) on an ascending `a`, found from position `hint` (the creation walk's
// previous range): a galloping probe toward x (steps of 1, 2, 4, ... between probes), then a
// binary search inside the last step. The answer is std::lower_bound's whatever the hint (the first
// position holding a value >= x, a.size() when none); only the probes differ.
[[nodiscard]] std::size_t lower_bound_from(std::span<const EdgeId> a, std::size_t hint, EdgeId x) noexcept;

// The creation products of a ShortestPaths flow: its record route, its plan and its route cost.
struct PreparedRouteResult {
  std::shared_ptr<const PredecessorDAG> route {};
  PreparedRoute prepared_route {};
  Cost min_route_cost {kUnreachable};
};

// Build a ShortestPaths flow's route, record route DAG, plan and cost in one pass. `support` contains
// all static shortest paths to the sorted `destinations`; start at `source` and apply the route
// choice on `view`. The result matches support_route, the route's construction and prepare_route:
// reuse `support` if it is the whole route, otherwise build the route with
// PredecessorDAGBuilder::of_edges.
// Requires multipath and static membership (all kernels except Weighted with ResidualSnapshot
// weights; a Weighted route edge of weight 0 also returns nullopt).
// While building, check: every route edge can be a member; none leaves a destination; every
// route node reaches a destination; the source is first in topological order; no cycle exists;
// and each tail's members are consecutive. These checks ensure members == route.
// Return nullopt on a failed condition or if `source` reaches no destination; the caller then uses
// the separate routines. No per-call allocation is proportional to the whole graph's size.
[[nodiscard]] std::optional<PreparedRouteResult> prepare_shortest_paths(const Graph& graph,
                                                                 const std::shared_ptr<const PredecessorDAG>& support,
                                                                 NodeId source, std::span<const NodeId> destinations,
                                                                 const RouteChoice& choice, PriorityReader view,
                                                                 RouteAllocation allocation, Admission admission,
                                                                 const WeightSource& weights);

// ---------------------------------------------------------------- allocation

// The cumulative allocation alloc(T) of a Weighted or Equal plan, member-aligned: T
// enters as round robin over the roots (all of it for one root); node by node in topological
// order, a node's inflow (its own injection plus the sum of its in-member allocations) is split
// by round robin or highest averages. With `limit` (member-aligned) it stops at the first member
// above its limit and returns false. A positive share of a root without a member (a share with no
// route) is a logic error: the Lossless searches stay at or below the number of roots before the
// first such root.
bool allocate(const PreparedRoute& plan, Amount total, std::span<Amount> out, std::span<const Amount> limit = {});

// Find the largest T whose allocation fits each entry of member-aligned view_m. The upper
// bound is the sum of the view on the roots' members, and at most the number of roots before the
// first root without a member. Each member's allocation never decreases as
// T grows, so a fitting and a failing T bound the answer. A floating-point estimate chooses
// the first probe only; then double the step size to narrow the bounds and bisect them.
// All fit checks, bounds, the returned value and allocation are exact. `iterations` counts
// at most 62 search probes (enough to bisect the units domain below 2^62); throw std::logic_error
// if that bound is exceeded. Evaluating the upper bound first, or a Lossless step's target,
// is an allocation evaluation and is not counted as a search probe.
struct LosslessBound {
  Amount value {0};
  std::int32_t probes {0};
};
[[nodiscard]] LosslessBound max_lossless_offer(const PreparedRoute& plan, std::span<const Amount> view_m);

// One Lossless step of one flow on an Equal or Weighted plan (the amounts stage): T_new = max(T,
// min(T + request, T_max(V))) on the replacement view V = W + own (member-aligned `view_m`), with
// `own_m` the flow's committed allocation alloc(T) (member-aligned). The step probes T + request
// first and searches T_max (as max_lossless_offer does) only when it does not fit. `request`
// nullopt is +inf. `edge_flows` is the new member-aligned allocation (own_m itself when T is
// unchanged).
struct LosslessStep {
  Amount admitted {0};
  std::vector<Amount> edge_flows {};
  std::int32_t probes {0};  // probes of this step's T_max search (at most 62)
};
[[nodiscard]] LosslessStep lossless_step(const PreparedRoute& plan, Amount total, std::span<const Amount> view_m,
                                         std::span<const Amount> own_m, std::optional<Amount> request);

// The injection of one forwarding: per root (aligned with plan.sources), the amount it injects; or
// (first_hop, the experiment of a Lossy +inf placement) every root's out-members, each injected
// with exactly its own view, while traffic reaching a root from another root is split as usual.
struct Injection {
  std::span<const Amount> per_root {};
  bool first_hop {false};
};

// The per-root part of a forwarding: what each root offered (first hop: injected), lost at the
// source (a root without a member keeps its injection and loses it), delivered, and its (member,
// carried, dropped)
// entries in ascending member order. At a router the members' amounts are attributed to the roots
// by highest averages over the roots' remaining inflows, member by member in member order; on a
// member that drops, its carried amount is attributed by highest averages over the roots'
// arrivals on it, the rest being each root's drop there.
struct Attribution {
  struct Entry {
    std::int32_t edge_position {0};
    Amount carried {0};
    Amount drops {0};
  };
  std::vector<Amount> offered {};
  std::vector<Amount> source_drops_by_source {};
  std::vector<Amount> delivered {};
  std::vector<std::vector<Entry>> entries {};
};

// One forwarding of `in` over the plan on the member-aligned view `w_m` (the working residual;
// empty: no capacity limit, the no-drop forwarding of an admitted total): every node splits its
// inflow (own injection plus carried arrivals) once by the kernel's rule, each member carries
// min(arrivals, view) and drops the rest; with `attr`, the per-root part as well. `offered` is
// the sum of the injections (first_hop: what was injected), `source_drops` the injections of roots
// without a member.
struct Forwarding {
  Amount offered {0};
  Amount network_entered {0};
  Amount source_drops {0};
  Amount delivered {0};
  std::vector<Amount> edge_flows {};  // member-aligned
  std::vector<Amount> edge_drops {};    // member-aligned
};
[[nodiscard]] Forwarding forward(const PreparedRoute& plan, const Injection& in, std::span<const Amount> w_m,
                                 Attribution* attr = nullptr);

} // namespace netgraph::core::detail
