/*
  T1 setup (contract 6.2 to 6.4): every head sets up its own strict tunnels, no controller.

  The design volume D is split RR(D, H, i) over the heads and each head's share RR(d, N, j) over
  its N tunnels, so tunnels of one head differ by at most one unit. The attempts run in round
  order (round j = 0..N-1, the heads ascending within a round), which is also the entry order
  j * H + i of the program; every attempt runs even after failures, on the reservation view U
  minus the requests of the tunnels already up in this setup. An attempt is one CSPF
  (route_search): the allowed edges with U >= max(b, min_unreserved_bandwidth) (no filter at a
  threshold of 0), Collapse of parallel members on the filtered candidates, the canonical single
  path on static costs (cost, then the bottleneck priority under PreferHigherPriority, fewer hops,
  the smaller EdgeId sequence), the priorities reading U. Then the cost ceiling. An up tunnel holds
  exactly b on every edge of its path; a down one holds nothing. There is no partial reservation,
  retry, growth or count reduction.

  Amounts: U[e] <= capacity < 2^62 and an up tunnel's request is at most U[e] on its edges, so U
  stays in [0, 2^62); the requests of one setup sum to D < 2^62.
*/
#include "netgraph/core/detail/tunnels.hpp"

#include <algorithm>
#include <utility>

#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/traversal.hpp"

namespace netgraph::core::detail {

std::vector<std::pair<NodeId, Amount>> head_design(std::span<const NodeId> roster, Amount D) {
  std::vector<std::pair<NodeId, Amount>> out;
  out.reserve(roster.size());
  const auto H = static_cast<std::int64_t>(roster.size());
  for (std::int64_t i = 0; i < H; ++i) {
    out.emplace_back(roster[static_cast<std::size_t>(i)], round_robin_share(D, H, i));
  }
  return out;
}

std::vector<SetupTunnel> t1_setup(const Graph& g, const Eligibility* elig, std::span<const NodeId> roster,
                                  std::span<const NodeId> dests, const PathSelection& selection,
                                  const CapacityAware& ca, Cost c0, Amount D, std::vector<Amount> U) {
  const auto H = static_cast<std::int64_t>(roster.size());
  const std::int64_t N = ca.tunnels.count;
  std::vector<SetupTunnel> out;
  out.reserve(static_cast<std::size_t>(H * N));
  for (std::int64_t j = 0; j < N; ++j) {
    for (std::int64_t i = 0; i < H; ++i) {
      SetupTunnel t;
      t.slot = static_cast<std::int32_t>(i);
      t.round = static_cast<std::int32_t>(j);
      t.requested = round_robin_share(round_robin_share(D, H, i), N, j);
      const Gate gate {elig, U, std::max(t.requested, ca.min_unreserved_bandwidth)};
      auto found = route_search(g, gate, selection, roster[static_cast<std::size_t>(i)], dests);
      if (!found.has_value()) {
        t.status = TunnelStatus::NoBandwidth;
      } else if (ca.max_path_cost.has_value() &&
                 (std::holds_alternative<Absolute>(*ca.max_path_cost)
                      ? found->cost > std::get<Absolute>(*ca.max_path_cost).cost
                      : !stretch_le(found->cost, c0, std::get<Stretch>(*ca.max_path_cost).factor))) {
        t.status = TunnelStatus::CostCeiling;
      } else {
        t.edges = std::move(found->edges);
        t.cost = found->cost;
        for (EdgeId e : t.edges) U[static_cast<std::size_t>(e)] -= t.requested;
      }
      out.push_back(std::move(t));
    }
  }
  return out;
}

} // namespace netgraph::core::detail
