/* Forwarding programs: the routes and integer weights a computation produces and forwarding moves
   an offer over (Tunnel, RouteDefinition, ForwardingProgram). */
#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

namespace detail {
struct ProgramData;
struct ProgramAccess;
struct ProgramBlock;
} // namespace detail

// A tunnel's state: Up; InsufficientBandwidth and CostCeilingExceeded (T1 setup found no path passing the bandwidth
// filter, or the cheapest one exceeds the cost ceiling: the tunnel is down, has no path and weight
// 0, and stays down in every forwarding of its program); and, in one forwarding only (a by_tunnel
// row), Failed (the path touches a node or an edge the forwarding policy's Eligibility excludes:
// the tunnel gets nothing and its head splits over its other entries) and Blackholed (up, but its
// path reaches a black hole of stale forwarding: its traffic is dropped on the first link that is
// one or enters one).
enum class TunnelStatus { Up, InsufficientBandwidth, CostCeilingExceeded, Failed, Blackholed };

// One tunnel (a path entry) of a program. entry_id: the entry identifier (T1 and T2: j * H + i for
// round j and roster position i, T2's dead heads leaving holes; CentralizedTE and ShortestPaths Optimized: the
// path's position in the decomposition order; supplied routes: the route's position in the list),
// never a position in tunnels();
// local_index: its position among its head's entries of both kinds in ascending entry_id;
// destination: where its path ends; edges in path order; path_cost: the sum of their costs (all
// three absent for a tunnel that is down); status: Up, InsufficientBandwidth or CostCeilingExceeded;
// requested_bandwidth and reserved_bandwidth: T1's strict bandwidth (its request, and the request
// held on every link of its path when it is up, else 0; 0 for the other models);
// planned_bandwidth: T2's planned bandwidth u (metadata,
// never held; 0 for the other models); weight: its integer share among the head's entries (T1:
// 1 up, 0 down; T2: 1; a controller path: its solved amount).
struct Tunnel {
  std::int32_t entry_id {0};
  NodeId head {0};
  std::int32_t local_index {0};
  std::optional<NodeId> destination {};
  std::optional<std::vector<EdgeId>> edges {};
  std::optional<Cost> path_cost {};
  TunnelStatus status {TunnelStatus::Up};
  Amount requested_bandwidth {0};
  Amount reserved_bandwidth {0};
  Amount planned_bandwidth {0};
  Amount weight {0};
  friend bool operator==(const Tunnel&, const Tunnel&) = default;
};

// A supplied route of ForwardingProgram::from_routes: a Path with one edge per hop (a tunnel) or
// a PredecessorDAG (a split route, its own split context) starting at `head`, with its integer weight
// among the head's routes (0 allowed: it carries nothing).
struct RouteDefinition {
  NodeId head {0};
  std::variant<std::shared_ptr<const Path>, std::shared_ptr<const PredecessorDAG>> route {};
  Amount weight {1};
};

// The interior split of the PredecessorDAG routes of from_routes and of from_dag's DAG: allocation
// Equal (absent) or Weighted, whose weights come from `weight_basis` (absent: Capacity), `edge_weights` (Explicit: one
// entry per edge) or `residual` (ResidualSnapshot: the free-capacity vector it reads, one entry
// per edge in [0, capacity]; it has no default).
struct RouteOptions {
  std::optional<RouteAllocation> allocation {};
  std::optional<WeightBasis> weight_basis {};
  std::optional<std::vector<Amount>> edge_weights {};
  std::optional<std::vector<Amount>> residual {};
};

// An immutable value: one demand's sources and destinations (normalized sets), its entries and
// its injection rule. Entries are tunnels (listed in tunnels(), ascending index) and split routes
// (a DAG with per-router split tables; not listed). Injection: equal sources (round robin over the
// roster: the sources with an entry, or for T1 and T2 the policy's roster, heads without an up
// tunnel included) or captured source weights (source_weights(): highest
// averages over every roster source's saved supply, zero included); then each head splits its
// share by highest averages over its entries' weights. Identical routes are never merged.
// A program also keeps what a replay under a failure rebuilds its roster from (its basis: the
// topology under the producing policy's exclusions, or each head's original routes). Shared by
// std::shared_ptr<const ForwardingProgram>; every policy that forwards it keeps that reference.
class ForwardingProgram {
public:
  // Validates and freezes supplied routes, std::invalid_argument for the first failure, in this
  // order: the endpoint sets (empty, range, lowest overlap); `routes must not be empty`; when no
  // route is a PredecessorDAG, any option given (`<allocation, weight_basis, edge_weights or residual>
  // applies only to PredecessorDAG routes`); allocation Optimized; the weights' combination;
  // edge_weights' length and entries and residual's (the ledger's baseline messages); then each
  // route in list order: its route's graph, its head a source, a Path's shape (one edge per hop, starting at the head,
  // ending at a destination and meeting none before), a PredecessorDAG's walk from the head to a
  // destination, then its weight in [0, 2^62). Order and duplicates are kept; a route of weight 0 is
  // validated in full. Every split route is its own split context.
  [[nodiscard]] static std::shared_ptr<const ForwardingProgram> from_routes(const Graph& graph,
                                                                           std::span<const NodeId> sources,
                                                                           std::span<const NodeId> destinations,
                                                                           std::span<const RouteDefinition> routes,
                                                                           const RouteOptions& options = {});

  // One shared hop-by-hop program on a supplied DAG of `graph`: its sources allocated equally
  // and pooled at every shared router exactly as ShortestPaths Equal/Weighted (the one pooled
  // builder), the DAG trimmed to the walks from each source to its first destinations.
  // std::invalid_argument for the first failure, in this order: the endpoint sets; `dag belongs to
  // another graph`; options.allocation Optimized (`from_dag allocation must be EQUAL or WEIGHTED`);
  // the weights' combination; edge_weights' and residual's lengths and entries. The roster is the
  // sources that reach a destination in the DAG (possibly none): one split route per roster source,
  // identified by roster position, and no tunnel; a replay rebuilds the roster on the DAG.
  [[nodiscard]] static std::shared_ptr<const ForwardingProgram> from_dag(const Graph& graph, std::span<const NodeId> sources,
                                                                        std::span<const NodeId> destinations,
                                                                        const std::shared_ptr<const PredecessorDAG>& dag,
                                                                        const RouteOptions& options = {});

  ForwardingProgram(const ForwardingProgram&) = delete;
  ForwardingProgram& operator=(const ForwardingProgram&) = delete;
  ~ForwardingProgram();

  [[nodiscard]] std::uint64_t graph_uid() const noexcept;
  [[nodiscard]] std::span<const NodeId> sources() const noexcept;
  [[nodiscard]] std::span<const NodeId> destinations() const noexcept;
  [[nodiscard]] std::span<const Tunnel> tunnels() const noexcept;
  // The captured source weights (roster source, weight) in ascending NodeId, every roster
  // source listed, zeros included: solved supplies for Optimized, declared weights for
  // Weighted CentralizedTE. nullopt for equal sources; an empty captured roster is engaged.
  [[nodiscard]] const std::optional<std::vector<std::pair<NodeId, Amount>>>& source_weights() const noexcept;

private:
  friend struct detail::ProgramAccess;
  friend struct detail::ProgramBlock;
  explicit ForwardingProgram(const detail::ProgramData* data) noexcept;
  // The program's data, in the same allocation as the program (detail::ProgramBlock): the
  // program's owner owns both.
  const detail::ProgramData* data_;
};

} // namespace netgraph::core
