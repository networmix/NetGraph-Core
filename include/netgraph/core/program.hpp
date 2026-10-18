/* Forwarding programs: the routes and integer weights a computation produces and forwarding moves
   an offer over (Tunnel, Route, ForwardingProgram). */
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
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

namespace detail {
struct ProgramData;
struct ProgramAccess;
struct ProgramBlock;
} // namespace detail

// A tunnel's state in one forwarding (a by_tunnel row): Up; Failed (the path touches a node or an
// edge the forwarding policy's Eligibility excludes: the tunnel gets nothing and its head splits
// over its other entries); Blackholed (up, but its path reaches a black hole of stale forwarding:
// its traffic is dropped on the first link that is one or enters one). A program's own tunnels are
// all up.
enum class TunnelStatus { Up, Failed, Blackholed };

// One tunnel (a path entry) of a program. index: the entry identifier (Global and CostOnly
// MaxFlow: the path's position in the decomposition order; supplied routes: the route's position
// in the list), never a position in tunnels(); local_index: its position among its head's entries
// of both kinds in ascending index; destination: where its path ends; edges in path order;
// route_cost: the sum of their costs; weight: its integer share among the head's entries (a
// controller path: its solved amount).
struct Tunnel {
  std::int32_t index {0};
  NodeId head {0};
  std::int32_t local_index {0};
  NodeId destination {0};
  std::vector<EdgeId> edges {};
  Cost route_cost {0};
  Amount weight {0};
  friend bool operator==(const Tunnel&, const Tunnel&) = default;
};

// A supplied route of ForwardingProgram::from_routes: a Path with one edge per hop (a tunnel) or
// a PredDAG (a split route, its own split context) starting at `head`, with its integer weight
// among the head's routes (0 allowed: it carries nothing).
struct Route {
  NodeId head {0};
  std::variant<std::shared_ptr<const Path>, std::shared_ptr<const PredDAG>> bundle {};
  Amount weight {1};
};

// The interior split of the PredDAG routes of from_routes: placement EqualBalanced (absent) or
// Weighted, whose weights come from `weights` (absent: Capacity), `edge_weights` (Explicit: one
// entry per edge) or `residual` (ResidualSnapshot: the free-capacity vector it reads, one entry
// per edge in [0, capacity]; it has no default).
struct RouteOptions {
  std::optional<FlowPlacement> placement {};
  std::optional<WeightBasis> weights {};
  std::optional<std::vector<Amount>> edge_weights {};
  std::optional<std::vector<Amount>> residual {};
};

// An immutable value: one demand's sources and destinations (normalized sets), its entries and
// its injection rule. Entries are tunnels (listed in tunnels(), ascending index) and split routes
// (a DAG with per-router split tables; not listed). Injection: equal sources (round robin over the
// roster, the sources with an entry) or captured source weights (source_weights(): highest
// averages over every roster source's saved supply, zero included); then each head splits its
// share by highest averages over its entries' weights. Identical routes are never merged.
// A program also keeps what a replay under a failure rebuilds its roster from (its basis: the
// topology under the producing policy's exclusions, or each head's original routes). Shared by
// std::shared_ptr<const ForwardingProgram>; every policy that forwards it keeps that reference.
class ForwardingProgram {
public:
  // Validates and freezes supplied routes, std::invalid_argument for the first failure, in this
  // order: the endpoint sets (empty, range, lowest overlap); `routes must not be empty`; when no
  // route is a PredDAG, any option given (`<placement, weights, edge_weights or residual> applies
  // only to PredDAG routes`); placement MaxFlow; the weights' combination; edge_weights' length and
  // entries and residual's (the ledger's baseline messages); then each route in list order: its
  // bundle's graph, its head a source, a Path's shape (one edge per hop, starting at the head,
  // ending at a destination and meeting none before), a PredDAG's walk from the head to a
  // destination, then its weight in [0, 2^62). Order and duplicates are kept; a route of weight 0 is
  // validated in full. Every split route is its own split context.
  [[nodiscard]] static std::shared_ptr<const ForwardingProgram> from_routes(const Graph& g,
                                                                           std::span<const NodeId> sources,
                                                                           std::span<const NodeId> destinations,
                                                                           std::span<const Route> routes,
                                                                           const RouteOptions& options = {});

  ForwardingProgram(const ForwardingProgram&) = delete;
  ForwardingProgram& operator=(const ForwardingProgram&) = delete;
  ~ForwardingProgram();

  [[nodiscard]] std::uint64_t graph_uid() const noexcept;
  [[nodiscard]] std::span<const NodeId> sources() const noexcept;
  [[nodiscard]] std::span<const NodeId> destinations() const noexcept;
  [[nodiscard]] std::span<const Tunnel> tunnels() const noexcept;
  // The captured source weights (roster source, saved supply) in ascending NodeId, every roster
  // source listed; nullopt for equal sources.
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
