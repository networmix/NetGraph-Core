/* Placement kernels assign one flow's amounts to route edges. Apportionment splits an integer
   total among recipients.

   A plan freezes the route, its split graph and any Weighted weights at creation. The split
   graph keeps edges allowed by the kernel that lie on a source-to-destination walk; it removes
   edges leaving destinations. Its edges are the plan's members. A member-aligned array has one
   entry per member, in the same order as plan.members. Policies combine these allocations; the
   ledger commits them.
*/
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/detail/traversal.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

// ---------------------------------------------------------------- apportionment

// Round robin: give floor(total / k) to each of k recipients, plus one to the first
// total mod k in caller order. Edge shares use ascending EdgeId; shares among a source's
// flows use flow order. Used by EqualBalanced (Lossless and Lossy) at every node and to split a
// source's amount over its flows.
// Visit only the n positive shares (n = min(total, k)); later entries stay 0 in a zeroed output.
// The callback form computes the split once and calls share(i, x) in index order for i in
// [0, n). Stop if it returns false; return the accepted count: n, or that rejected index i.
// The span form writes out[0 .. n), where k = out.size(), and returns n.
template <class Share>
inline std::size_t round_robin(Amount total, std::size_t k, Share&& share) {
  const auto kk = static_cast<Amount>(k);
  const Amount base = total / kk;
  const auto extra = static_cast<std::size_t>(total % kk);
  const std::size_t n = base == 0 ? extra : k;
  for (std::size_t i = 0; i < n; ++i) {
    if (!share(i, base + (i < extra ? 1 : 0))) return i;
  }
  return n;
}

inline std::size_t round_robin(Amount total, std::span<Amount> out) noexcept {
  return round_robin(total, out.size(), [&](std::size_t i, Amount x) noexcept {
    out[i] = x;
    return true;
  });
}

// Highest averages (Jefferson/D'Hondt; Weighted and Lossy Weighted): give each next unit to
// the member with largest w / (a + 1), where w is its weight and a its allocation so far.
// Break ties by lower member index (ascending EdgeId). Start each member at its lower quota
// floor(total * w / W), where W is the weight sum; then assign the fewer than k remaining units
// by those quotients, comparing with exact 128-bit cross products. Weights are positive,
// W < 2^62, total < 2^63. Increasing total never reduces any member's allocation (house
// monotonicity).
void highest_averages(Amount total, std::span<const Amount> weights, std::span<Amount> out);

// ---------------------------------------------------------------- plans

// Throws std::invalid_argument for a value outside its enumeration and for the unsupported
// combinations (Lossy with MaxFlow, Weighted with priority Headroom).
void check_kernel(FlowPlacement placement, Admission admission, EdgePriority priority);
// The combination rule of check_kernel alone, for a caller that has already checked the three
// enumerations (validate_config).
void check_kernel_combination(FlowPlacement placement, Admission admission, EdgePriority priority);

struct KernelPlan {
  FlowPlacement placement {FlowPlacement::MaxFlow};
  Admission admission {Admission::Lossless};
  NodeId src {0};
  std::vector<EdgeId> route {};      // ascending: the flow's route (its record's bundle)
  std::vector<EdgeId> members {};    // ascending: the split graph
  std::vector<Amount> weights {};    // Weighted: per member (positive), frozen; else empty
  // Destinations ascend. Nodes have a topological order: every edge's tail precedes its head,
  // with the source first. Position p's out-members are members[out_begin[p] .. out_end[p]);
  // they are consecutive in ascending EdgeId because EdgeIds sort by source first. Each member
  // also stores its head's position. Empty split graph means empty structure. All views share
  // one block owned by the plan (layout()), so the structure takes one allocation.
  std::span<const NodeId> dsts {};
  std::span<const NodeId> topo {};
  std::span<const std::int32_t> out_begin {};
  std::span<const std::int32_t> out_end {};
  std::span<const std::int32_t> head_pos {};
  std::span<const std::uint8_t> terminal {};
  std::int32_t src_out_begin {0};
  std::int32_t src_out_end {0};

  KernelPlan() = default;
  KernelPlan(const KernelPlan& other);
  KernelPlan(KernelPlan&& other) noexcept;
  KernelPlan& operator=(const KernelPlan& other);
  KernelPlan& operator=(KernelPlan&& other) noexcept;
  ~KernelPlan() = default;

  // The writable views of a new block: the destinations copied in, the structure for `positions`
  // nodes and `num_members` members zeroed. One allocation, made before anything changes (a
  // failure leaves the plan as it was); the previous block is released.
  struct Layout {
    std::span<NodeId> topo;
    std::span<std::int32_t> out_begin;
    std::span<std::int32_t> out_end;
    std::span<std::int32_t> head_pos;
    std::span<std::uint8_t> terminal;
  };
  Layout layout(std::span<const NodeId> destinations, std::size_t positions, std::size_t num_members);

private:
  // Points the views into store_ for these sizes (the layout's order: dsts, topo, out_begin,
  // out_end, head_pos, then terminal as bytes).
  Layout point(std::size_t num_dsts, std::size_t positions, std::size_t num_members) noexcept;
  std::vector<std::int32_t> store_ {};
};

// The edges of `dag` that lie on some walk from `src` to the destination set (edges leaving a
// destination removed: destinations are terminal), ascending.
[[nodiscard]] std::vector<EdgeId> prune(const Graph& g, const PredDAG& dag, NodeId src,
                                        std::span<const NodeId> dsts);

// The minimum walk cost from `src` to the destination set over `route` (a record's bundle:
// its route cost); kUnreachable when the route has no such walk. The split-graph builder's
// cost pass alone: no plan is built.
[[nodiscard]] Cost route_cost(const Graph& g, const PredDAG& route, NodeId src,
                              std::span<const NodeId> dsts);

// The plan of a new flow on `route` (acyclic; the flow's DAG or path, its record's bundle), its
// membership captured on `creation_view` (aligned with the route's edges in ascending EdgeId:
// the view of the search that created the route; read only by a Weighted Residual membership,
// may be empty otherwise): MaxFlow and EqualBalanced every route edge (an edge of capacity 0 is
// a member that carries nothing: EqualBalanced still gives it its round-robin share), Weighted
// the edges of positive weight by `priority` (Residual: creation view >= 1 unit; Capacity: static
// capacity >= 1 unit) with the weights read from the same values; Lossy keeps its kernel's
// membership. The members are pruned to the walks from src to the destination set (none
// leaving a destination) and given their topological structure, in one pass. With `cost`, the
// route's minimum walk cost (route_cost) is computed in the same pass when every route edge is a
// member. `dsts` ascending and disjoint from src.
[[nodiscard]] KernelPlan make_plan(const Graph& g, FlowPlacement placement, Admission admission,
                                   EdgePriority priority, NodeId src, std::span<const NodeId> dsts,
                                   const PredDAG& route, std::span<const Amount> creation_view, Cost* cost = nullptr);

// The usable route of a Pinned bundle as a plan: the bundle's edges that `elig` allows (every
// one when null), whatever their capacity, pruned to the walks from src to the destination set,
// with make_plan's membership: route == members for MaxFlow and EqualBalanced; Weighted with
// Capacity keeps the route edges of positive capacity on the walks that remain (its weights the
// capacities). It is the plan of every kernel whose membership is static; a Weighted Residual
// plan is made from its route at creation (make_plan). With `cost`, the usable route's minimum
// walk cost (kUnreachable when it is empty), computed in the same pass.
[[nodiscard]] KernelPlan make_route_plan(const Graph& g, FlowPlacement placement, Admission admission,
                                         EdgePriority priority, NodeId src, std::span<const NodeId> dsts,
                                         const PredDAG& bundle, const Eligibility* elig, Cost* cost = nullptr);

// ---------------------------------------------------------------- CostOnly creation

// std::lower_bound(a, x) on an ascending `a`, found from position `hint` (the creation walk's
// previous range): a galloping probe toward x (steps of 1, 2, 4, ... between probes), then a
// binary search inside the last step. The answer is std::lower_bound's whatever the hint (the first
// position holding a value >= x, a.size() when none); only the probes differ.
[[nodiscard]] std::size_t lower_bound_from(std::span<const EdgeId> a, std::size_t hint, EdgeId x) noexcept;

// The creation products of a CostOnly flow: its record bundle, its plan and its route cost.
struct CreationProducts {
  std::shared_ptr<const PredDAG> bundle {};
  KernelPlan plan {};
  Cost cost {kUnreachable};
};

// Build a CostOnly flow's route, record bundle, plan and cost in one pass. `support` contains
// all static shortest paths to sorted destinations `dsts`; start at `src` and apply selection
// on `view`. The result matches support_route, bundle construction and make_plan: reuse
// `support` if it is the whole route, otherwise build its bundle with PredDAGBuilder::of_edges.
// Requires multipath and static membership (all kernels except Weighted with Residual).
// While building, check: every route edge can be a member; none leaves a destination; every
// route node reaches a destination; the source is first in topological order; no cycle exists;
// and each tail's members are consecutive. These checks ensure members == route.
// Return nullopt on a failed condition or if src reaches no destination; the caller then uses
// the separate routines. No per-call allocation is proportional to the whole graph's size.
[[nodiscard]] std::optional<CreationProducts> cost_only_creation(const Graph& g,
                                                                 const std::shared_ptr<const PredDAG>& support,
                                                                 NodeId src, std::span<const NodeId> dsts,
                                                                 const PathSelection& selection, PriorityFn view,
                                                                 FlowPlacement placement, Admission admission);

// ---------------------------------------------------------------- allocation

// The cumulative allocation alloc(T) of a Weighted or EqualBalanced plan, member-aligned:
// node by node in topological order, the source's inflow T, every other node's inflow the sum
// of its in-member allocations, split by round robin or highest averages. With `limit`
// (member-aligned) it stops at the first member above its limit and returns false.
bool allocate(const KernelPlan& plan, Amount total, std::span<Amount> out, std::span<const Amount> limit = {});

// Find the largest T whose allocation fits each entry of member-aligned view_m. The upper
// bound is the sum of the view on source members. Each member's allocation never decreases as
// T grows, so a fitting and a failing T bound the answer. A floating-point estimate chooses
// the first probe only; then double the step size to narrow the bounds and bisect them.
// All fit checks, bounds, the returned value and allocation are exact. `iterations` counts
// at most 62 search probes (enough to bisect the units domain below 2^62); throw std::logic_error
// if that bound is exceeded. Evaluating the upper bound first, or a Lossless step's target,
// is an allocation evaluation and is not counted as a search probe.
struct TMax {
  Amount value {0};
  std::int32_t iterations {0};
};
[[nodiscard]] TMax t_max(const KernelPlan& plan, std::span<const Amount> view_m);

// One Lossless step of one flow (the amounts stage): T_new = max(T, min(T + request,
// T_max(V))) on the replacement view V = W + own (member-aligned `view_m`), with `own_m` the
// flow's committed allocation alloc(T) (member-aligned). MaxFlow solves the flow core on the
// split graph with cap T + request (uncapped for +inf) and reports the solve's marginal curve;
// Weighted and EqualBalanced probe T + request first and search T_max (as t_max does) only when
// it does not fit: the result is the same T_new. `request` nullopt is +inf. `alloc` is the new
// member-aligned allocation (own_m itself when T is unchanged on a Weighted or EqualBalanced
// plan).
struct LosslessStep {
  Amount total {0};
  std::vector<Amount> alloc {};
  std::int32_t iterations {0};                                // probes of this step's T_max search (at most 62)
  std::vector<std::pair<Cost, Amount>> cost_distribution {};  // MaxFlow: the solve's marginal curve
};
[[nodiscard]] LosslessStep lossless_step(const Graph& g, const KernelPlan& plan, Amount total,
                                         std::span<const Amount> view_m, std::span<const Amount> own_m,
                                         std::optional<Amount> request);

// One Lossy offer on the working residual (member-aligned `w_m`): split at every node by the
// kernel's rule, each member carrying min(share, residual) and dropping the rest; `offer`
// nullopt is +inf (every source member injected to its residual). `offered` is the offer, or
// the injection for +inf; a source without members injects nothing.
struct LossyStep {
  Amount offered {0};
  Amount injected {0};
  Amount delivered {0};
  std::vector<Amount> carried {};  // member-aligned
  std::vector<Amount> drops {};    // member-aligned
};
[[nodiscard]] LossyStep lossy_step(const KernelPlan& plan, std::optional<Amount> offer, std::span<const Amount> w_m);

} // namespace netgraph::core::detail
