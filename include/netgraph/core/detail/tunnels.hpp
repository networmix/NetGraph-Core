/* The tunnel models' computations (internal). T1 (tunnels_t1.cpp): every head sets up its own
   tunnels by CSPF on the reservation view; its program is built by tunnel_program
   (program_data.hpp) and forwarded by the forwarding engine like any other program. */
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

// T1's setup (6.2 to 6.4) for the design volume D (< 2^62) over the roster's H heads (ascending)
// with N = ca.tunnels.count tunnels each: head i's share RR(D, H, i), its tunnel j's request
// RR(share, N, j); the attempts in round order (round j, heads ascending), every one even after
// failures, each on `U` (dense, E entries: the unreserved capacity with this policy's own holds
// released, each entry in [0, capacity]) minus the requests of the tunnels already up in this
// setup. An attempt keeps the allowed edges (`elig`) with U >= max(request,
// min_unreserved_bandwidth) (a threshold of 0 keeps every allowed edge), collapses parallel edges
// (selection.multi_edge = false) and takes the single path in the selection's order on static
// costs toward the allowed destinations `dests` (its priorities read U); then the ceiling:
// Absolute(c) inclusive, Stretch(f) as cost <= f * c0 exactly (`c0`: the pool baseline, read only
// for Stretch). Up: the request is subtracted from U on every edge of the path (U is consumed).
// The tunnels come in round order (index = position: round * H + slot).
[[nodiscard]] std::vector<SetupTunnel> t1_setup(const Graph& g, const Eligibility* elig, std::span<const NodeId> roster,
                                                std::span<const NodeId> dests, const PathSelection& selection,
                                                const CapacityAware& ca, Cost c0, Amount D, std::vector<Amount> U);

// The heads' design shares RR(D, H, i) as (head, share) in roster order (TunnelPlan::head_design).
[[nodiscard]] std::vector<std::pair<NodeId, Amount>> head_design(std::span<const NodeId> roster, Amount D);

} // namespace netgraph::core::detail
