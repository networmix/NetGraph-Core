/* Apportionment and the placement kernels in units (internal; APPROACH 2.1, 7 and 16).

   A kernel computes the pushes of one flow on its route. The plan of a flow is frozen at its
   creation: the route, its split graph (the route's edges that pass the kernel's membership
   rule, pruned to the walks from the source to the destination set, destinations terminal)
   and, for Weighted, the weights. These are the per-flow routines of the placement transaction
   (APPROACH 9.2 stage 3); the policy engine composes them and the ledger commits their output.
*/
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

namespace netgraph::core::detail {

// ---------------------------------------------------------------- apportionment (2.1)

// Round robin (EqualBalanced and Lossy EqualBalanced at every node; the cohort's split over its
// flows and a source's offer over its flows): floor(total / k) to each of the k members, one
// more unit to the first total mod k members (ascending EdgeId). Only the n leading members with
// a positive share (n = total when total < k, else k) are visited: the members from n on carry
// 0, which a caller with a zeroed output need not write or visit. The callback form computes the
// split once and calls share(i, x) for i in [0, n), in index order; share returns false to stop
// early. It returns the number of shares accepted: n, or the index i whose call returned false.
// The span form writes out[0 .. n) of the k = out.size() members and returns n.
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

// Highest averages (Jefferson, D'Hondt; Weighted and Lossy Weighted): the units go one at a
// time to the member maximising w / (a + 1), ties by the lower member index (ascending
// EdgeId). Closed form: every member starts at its lower quota floor(total * w / W), then the
// fewer than k remaining units go by the exact next quotients, compared by 128-bit cross
// products. Weights positive with a sum below 2^62; total below 2^63. House-monotone.
void highest_averages(Amount total, std::span<const Amount> weights, std::span<Amount> out);

// ---------------------------------------------------------------- plans

// Throws std::invalid_argument for a value outside its enumeration and for the combinations the
// contract rejects (Lossy with MaxFlow, Weighted with priority Headroom).
void check_kernel(FlowPlacement placement, Admission admission, EdgePriority priority);

struct KernelPlan {
  FlowPlacement placement {FlowPlacement::MaxFlow};
  Admission admission {Admission::Lossless};
  NodeId src {0};
  std::vector<NodeId> dsts {};       // ascending
  std::vector<EdgeId> route {};      // ascending: the flow's route (its record's bundle)
  std::vector<EdgeId> members {};    // ascending: the split graph
  std::vector<Amount> weights {};    // Weighted: per member (positive), frozen; else empty
  // The split graph's nodes in a topological order (the source first), with per position its
  // out-members members[out_begin[p] .. out_end[p]) (consecutive in ascending EdgeId: B' ids
  // rank edges by source first) and per member the position of its head. Empty when the split
  // graph is empty.
  std::vector<NodeId> topo {};
  std::vector<std::int32_t> out_begin {};
  std::vector<std::int32_t> out_end {};
  std::vector<std::int32_t> head_pos {};
  std::vector<std::uint8_t> terminal {};
  std::int32_t src_out_begin {0};
  std::int32_t src_out_end {0};
};

// The edges of `dag` that lie on some walk from `src` to the destination set (edges leaving a
// destination removed: destinations are terminal), ascending.
[[nodiscard]] std::vector<EdgeId> prune(const Graph& g, const PredDAG& dag, NodeId src,
                                        std::span<const NodeId> dsts);

// The minimum walk cost from `src` to the destination set over `route` (a record's bundle;
// APPROACH 8, rule 12); kUnreachable when the route has no such walk. The split-graph builder's
// cost pass alone: no plan is built.
[[nodiscard]] Cost route_cost(const Graph& g, const PredDAG& route, NodeId src,
                              std::span<const NodeId> dsts);

// The plan of a new flow on `route` (acyclic; the flow's DAG or path, its record's bundle), its
// membership captured on `creation_view` (aligned with the route's edges in ascending EdgeId:
// the view of the search that created the route; read only by a Weighted Residual membership,
// may be empty otherwise): MaxFlow every route edge (its view decides), Weighted by `priority`
// (Residual: creation view >= 1 unit; Capacity: static capacity >= 1 unit) with the weights read
// from the same values, EqualBalanced static capacity >= 1 unit; Lossy keeps its kernel's
// membership. The members are pruned to the walks from src to the destination set (none
// leaving a destination) and given their topological structure, in one pass. With `cost`, the
// route's minimum walk cost (route_cost) is computed in the same pass when every route edge is a
// member. `dsts` ascending and disjoint from src.
[[nodiscard]] KernelPlan make_plan(const Graph& g, FlowPlacement placement, Admission admission,
                                   EdgePriority priority, NodeId src, std::span<const NodeId> dsts,
                                   const PredDAG& route, std::span<const Amount> creation_view, Cost* cost = nullptr);

// The usable route of a Pinned bundle as a plan: of the bundle's edges the ones with static
// capacity >= 1 unit that `elig` allows (every one when null), pruned to the walks from src to
// the destination set; route == members, with the topological structure. It is the plan of
// every kernel whose membership is static (MaxFlow, EqualBalanced, Weighted with Capacity:
// weights the capacities); a Weighted Residual plan is made from its route at creation
// (make_plan). With `cost`, the usable route's minimum walk cost (kUnreachable when it is
// empty), computed in the same pass.
[[nodiscard]] KernelPlan make_route_plan(const Graph& g, FlowPlacement placement, Admission admission,
                                         EdgePriority priority, NodeId src, std::span<const NodeId> dsts,
                                         const PredDAG& bundle, const Eligibility* elig, Cost* cost = nullptr);

// ---------------------------------------------------------------- allocation (16.2, 16.3)

// The cumulative allocation alloc(T) of a Weighted or EqualBalanced plan, member-aligned:
// node by node in topological order, the source's inflow T, every other node's inflow the sum
// of its in-member allocations, split by round robin or highest averages. With `limit`
// (member-aligned) it stops at the first member above its limit and returns false.
bool allocate(const KernelPlan& plan, Amount total, std::span<Amount> out, std::span<const Amount> limit = {});

// The largest T whose allocation fits `view_m` (member-aligned) member by member, by search on
// alloc(T) (exact: alloc is non-decreasing member by member, so fits(T) is monotone); upper
// bound the source members' view. APPROACH 7 (revision 26): approximate arithmetic only
// chooses probes. The first probe is at a floating-point estimate, then a galloping bracket,
// then bisection; every predicate, both bracket ends, the value and the allocation are exact,
// and `iterations` counts the probes, which stay within 62 (the bound of a plain bisection over
// units below 2^62; checked, std::logic_error otherwise). A first evaluation at the upper bound
// (or a Lossless step's target) is another allocation evaluation, not a search probe.
struct TMax {
  Amount value {0};
  std::int32_t iterations {0};
};
[[nodiscard]] TMax t_max(const KernelPlan& plan, std::span<const Amount> view_m);

// One Lossless step of one flow (APPROACH 9.2 stage 3, 16.3): T_new = max(T, min(T + request,
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
