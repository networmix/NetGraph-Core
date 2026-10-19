/* The tunnel models' computations (internal). T1 (distributed_te.cpp): every head sets up its own
   tunnels by CSPF on the reservation view. T2 (centralized_te_tunnels.cpp): a controller lays out the same
   number of tunnels with one common planned bandwidth for every live head. Both programs are built
   by distributed_te_program (program_data.hpp) and forwarded by the forwarding engine like any other
   program. */
#pragma once

#include <optional>
#include <span>
#include <vector>

#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

// T1's setup for the design volume D (< 2^62) over the roster's H heads (ascending)
// with N = ca.tunnel_count.per_head tunnels each: head i's share RR(D, H, i), its tunnel j's request
// RR(share, N, j); the attempts in round order (round j, heads ascending), every one even after
// failures, each on `U` (dense, E entries: the unreserved capacity with this policy's own holds
// released, each entry in [0, capacity]) minus the requests of the tunnels already up in this
// setup. An attempt keeps the allowed edges (`eligibility`) with U >= max(request,
// min_unreserved_bandwidth) (a threshold of 0 keeps every allowed edge), collapses parallel edges
// (ca.parallel_edges = false) and takes the single path in the order of ca.tie_break on static
// costs toward the allowed `destinations` (its scores read U); then the ceiling:
// AbsoluteCost(c) inclusive, CostStretch(f) as cost <= f * c0 exactly (`c0`: the pool baseline, read only
// for CostStretch). Up: the request is subtracted from U on every edge of the path (U is consumed).
// The tunnels come in round order (index = position: round * H + slot).
[[nodiscard]] std::vector<TunnelSetupResult> setup_distributed_te_tunnels(const Graph& graph, const Eligibility* eligibility, std::span<const NodeId> reachable_sources,
                                                std::span<const NodeId> destinations, const DistributedTE& ca, Cost c0, Amount D, std::vector<Amount> U);

// T2's layout on V (dense, E entries: the free capacity with this policy's own
// traffic released) for the roster's H heads (ascending) toward the allowed `destinations`
// under `eligibility`: the live heads (those reaching a destination over links with V >= 1), the count k
// (tunnels.per_head, or the chosen count of "up to"), U (the largest common whole-unit size of k
// tunnels per live head: common), the exact fluid optimum s* = fluid_num / fluid_den (reduced; 0 / 1
// without a live head), the working size u (planned: min(U, max(1, max live ceil(RR(D, H, i) / k)))
// for a finite design D, U for an unlimited one), and at u (U >= 1) the tunnels in round order
// (slot: the roster position, round: the local index, planned u) with their route-cost sum J and
// (with `headroom`) their headroom penalty Phi (t2_objectives). No live head: count N (1 for "up
// to") and nothing else; U = 0 with live heads: no tunnels.
struct T2Layout {
  std::int64_t count {0};
  std::vector<NodeId> live {};
  Amount common {0};
  Amount planned {0};
  Amount fluid_num {0};
  Amount fluid_den {1};
  std::vector<TunnelSetupResult> tunnels {};
  Counter128 route_cost {};
  Counter128 headroom_penalty {};
};
[[nodiscard]] T2Layout t2_layout(const Graph& g, const Eligibility* eligibility, std::span<const NodeId> roster,
                                 std::span<const NodeId> destinations, const TunnelCount& tunnels, bool headroom,
                                 std::optional<Amount> design, std::span<const Amount> V);

// T2's objectives of a count vector x over `edges` (x[i] tunnels on edges[i]) at size u >= 1
// on V (dense, E entries): the route cost J = sum cost x and, with `headroom`, the headroom penalty
// Phi = the sum over links of capacity B > 0 of floor(2^32 (B - V + j u) / B) for j = 1..x (each
// link's sum in closed form). Exact 128-bit: J < K (cost sum) < 2^92; Phi <= 2^32 sum x, where
// sum x counts every link of every tunnel (sum x < 2^62; Phi itself can pass 2^64 for valid
// inputs), so Phi < 2^94. With `headroom`, a link of capacity B > 0 requires V <= B and x u <= V
// (std::logic_error otherwise).
struct T2Objectives {
  Counter128 route_cost {};
  Counter128 headroom_penalty {};
};
[[nodiscard]] T2Objectives t2_objectives(const Graph& g, std::span<const EdgeId> edges, std::span<const Amount> x,
                                         std::span<const Amount> V, Amount u, bool headroom);

// The heads' design shares RR(D, H, i) as (head, share) in roster order (TunnelPlan::design_by_head).
[[nodiscard]] std::vector<std::pair<NodeId, Amount>> design_by_head(std::span<const NodeId> reachable_sources, Amount D);

} // namespace netgraph::core::detail
