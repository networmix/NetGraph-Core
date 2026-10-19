/*
  T1 setup: every head sets up its own strict tunnels, no controller.

  The design volume D is split RR(D, H, i) over the heads and each head's share RR(d, N, j) over
  its N tunnels, so tunnels of one head differ by at most one unit. The attempts run in round
  order (round j = 0..N-1, the heads ascending within a round), which is also the entry order
  j * H + i of the program; every attempt runs even after failures, on the reservation view U
  minus the requests of the tunnels already up in this setup. An attempt is one CSPF
  (route_search): the allowed edges with U >= max(b, min_unreserved_bandwidth) (no filter at a
  threshold of 0), Collapse of parallel members on the filtered candidates, the canonical single
  path on static costs (cost, then the bottleneck score of a Higher* tie rule, fewer hops, the
  smaller EdgeId sequence), the scores reading U. Then the cost ceiling. An up tunnel holds
  exactly b on every edge of its path; a down one holds nothing. There is no partial reservation,
  retry, growth or count reduction.

  Amounts: U[e] <= capacity < 2^62 and an up tunnel's request is at most U[e] on its edges, so U
  stays in [0, 2^62); the requests of one setup sum to D < 2^62.
*/
#include "netgraph/core/detail/tunnel_setup.hpp"

#include <algorithm>
#include <utility>

#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/route_kernels.hpp"
#include "netgraph/core/detail/search_detail.hpp"

namespace netgraph::core::detail {

std::vector<std::pair<NodeId, Amount>> design_by_head(std::span<const NodeId> reachable_sources, Amount D) {
  std::vector<std::pair<NodeId, Amount>> out;
  out.reserve(reachable_sources.size());
  const auto H = static_cast<std::int64_t>(reachable_sources.size());
  for (std::int64_t i = 0; i < H; ++i) {
    out.emplace_back(reachable_sources[static_cast<std::size_t>(i)], round_robin_share(D, H, i));
  }
  return out;
}

std::vector<TunnelSetupResult> setup_distributed_te_tunnels(const Graph& graph, const Eligibility* eligibility, std::span<const NodeId> reachable_sources,
                                  std::span<const NodeId> destinations, const DistributedTE& ca, Cost c0, Amount D,
                                  std::vector<Amount> U) {
  // A tunnel is one path: the CSPF orders single paths by the model's tie rule.
  const RouteChoice choice {ca.tie_break, false, ca.parallel_edges};
  const auto H = static_cast<std::int64_t>(reachable_sources.size());
  const std::int64_t N = ca.tunnel_count.per_head;
  std::vector<TunnelSetupResult> out;
  out.reserve(static_cast<std::size_t>(H * N));
  for (std::int64_t j = 0; j < N; ++j) {
    for (std::int64_t i = 0; i < H; ++i) {
      TunnelSetupResult t;
      t.roster_index = static_cast<std::int32_t>(i);
      t.local_index = static_cast<std::int32_t>(j);
      t.requested_bandwidth = round_robin_share(round_robin_share(D, H, i), N, j);
      const SearchFilter gate {eligibility, U, std::max(t.requested_bandwidth, ca.min_unreserved_bandwidth)};
      auto found = route_search(graph, gate, choice, reachable_sources[static_cast<std::size_t>(i)], destinations);
      if (!found.has_value()) {
        t.status = TunnelStatus::InsufficientBandwidth;
      } else if (ca.cost_ceiling.has_value() &&
                 (std::holds_alternative<AbsoluteCost>(*ca.cost_ceiling)
                      ? found->cost > std::get<AbsoluteCost>(*ca.cost_ceiling).cost
                      : !stretch_le(found->cost, c0, std::get<CostStretch>(*ca.cost_ceiling).numerator(),
                                  std::get<CostStretch>(*ca.cost_ceiling).denominator()))) {
        t.status = TunnelStatus::CostCeilingExceeded;
      } else {
        t.edges = std::move(found->edges);
        t.path_cost = found->cost;
        for (EdgeId e : t.edges) U[static_cast<std::size_t>(e)] -= t.requested_bandwidth;
      }
      out.push_back(std::move(t));
    }
  }
  return out;
}

} // namespace netgraph::core::detail
