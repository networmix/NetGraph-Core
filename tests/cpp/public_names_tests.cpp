// The C++ public names: each renamed declaration is named once here, through the public headers
// alone, so a rename that is undone or a name that is removed fails to compile; the enumerators of
// every public enumeration; and the earlier header names, which no longer exist.

#include <gtest/gtest.h>

#include <initializer_list>
#include <span>
#include <type_traits>
#include <vector>

#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/forwarding_program.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/search.hpp"
#include "netgraph/core/shortest_path_cache.hpp"
#include "netgraph/core/types.hpp"

using namespace netgraph::core;

namespace {

// The earlier headers are gone.
#if __has_include("netgraph/core/dag_cache.hpp") || __has_include("netgraph/core/pred_dag.hpp") || \
    __has_include("netgraph/core/program.hpp")
#error "an earlier public header still exists"
#endif

// Types.
static_assert(std::is_class_v<AbsoluteCost> && std::is_class_v<CostStretch> && std::is_class_v<TunnelCount>);
static_assert(std::is_class_v<ShortestPaths> && std::is_class_v<DistributedTE> && std::is_class_v<CentralizedTE> &&
              std::is_class_v<Pinned>);
static_assert(std::is_class_v<ShortestPathCache> && std::is_base_of_v<std::exception, ShortestPathCacheMiss>);
static_assert(std::is_class_v<PredecessorDAG> && std::is_class_v<EvaluationOptions> && std::is_class_v<RouteDefinition>);
static_assert(std::is_enum_v<RouteTieBreak> && std::is_enum_v<RouteAllocation> && std::is_enum_v<SourceAllocation>);

// Fields.
static_assert(requires(DistributedTE x) { x.cost_ceiling; x.tie_break; x.parallel_edges; });
static_assert(requires(ShortestPaths x) {
  x.allocation;
  x.source_allocation;
  x.weight_basis;
  x.edge_weights;
  x.tie_break;
  x.multipath;
  x.parallel_edges;
});
static_assert(requires(CentralizedTE x) { x.source_allocation; x.support; x.tunnel_count; x.prefer_headroom; });
static_assert(requires(FlowRecord x) { x.route; x.edge_drops; x.edge_flows; x.entry_id; x.policy_id; x.min_route_cost; });
static_assert(requires(KspOptions x) { x.eligibility; x.cost_ceiling; x.tie_break; x.parallel_edges; });
static_assert(requires(SpfOptions x) { x.eligibility; x.tie_break; x.multipath; x.parallel_edges; });
static_assert(requires(OutputOptions x) { x.saturated_boundary_edges; });
static_assert(requires(PlacementResult x) {
  x.saturated_boundary_edges;
  x.call_edge_drops;
  x.call_source_drops;
  x.marginal_cost_curve;
  x.edge_drops;
  x.call_delivered;
  x.call_undelivered;
  x.source_drops;
});
static_assert(requires(EvaluationOptions x) { x.eligibility; x.offer; x.design_demand; });
static_assert(requires(RouteDefinition x) { x.route; x.head; x.weight; });
static_assert(requires(RouteOptions x) { x.allocation; x.weight_basis; x.edge_weights; x.residual; });
static_assert(requires(SensitivityEntry x) { x.delivery_reduction; });
static_assert(requires(SourceFlow x) { x.total_edge_drops; x.source_drops; });
static_assert(requires(Tunnel x) { x.entry_id; x.path_cost; x.planned_bandwidth; });
static_assert(requires(TunnelFlow x) { x.total_edge_drops; x.entry_id; });
static_assert(requires(TunnelCount x) { x.per_head; x.up_to; });

// Functions and methods.
static_assert(requires(const Eligibility& e) { e.allows_edge(EdgeId {0}); e.allows_node(NodeId {0}); });
static_assert(requires(const FlowLedger& f, std::span<const NodeId> s) {
  f.carried();
  f.saturated_boundary_edges(s);
  f.edge_capacity_units();
});
static_assert(requires(const FlowPolicy& p) { p.delivered(); });
static_assert(requires(const Graph& g) {
  g.edge_capacity_units();
  g.edge_costs();
  g.edge_cost_gcd();
  g.edge_destinations();
  g.edge_sources();
  g.incoming_edges();
  g.incoming_offsets();
  g.input_positions();
  g.max_edge_cost();
  g.outgoing_offsets();
});
static_assert(requires(const Path& p) { p.source(); p.destination(); });
static_assert(requires(const Graph& g, std::span<const NodeId> s, const EvaluationOptions& o) {
  evaluate_demand(g, s, s, o);
  evaluate_edge_failures(g, s, s, o, 1);
});
static_assert(requires(const Graph& g, std::span<const NodeId> s, std::shared_ptr<const PredecessorDAG> dag) {
  ForwardingProgram::from_dag(g, s, s, dag);
});

}  // namespace

// Every enumerator of the public enumerations, in declaration order (the order of the stub and
// docs/semantics.md). Only TunnelStatus's and Termination's integer values are promised
// (docs/semantics.md, "Termination" and "Forwarding programs"); for the others only the order is
// checked.
template <typename E>
constexpr bool declared_in_order(std::initializer_list<E> names) {
  const E* prev = nullptr;
  for (const E& e : names) {
    if (prev != nullptr && !(*prev < e)) return false;
    prev = &e;
  }
  return true;
}

TEST(PublicNames, Enumerations) {
  static_assert(declared_in_order({RouteTieBreak::FewerHops, RouteTieBreak::HigherCapacity,
                                   RouteTieBreak::HigherAvailableCapacity, RouteTieBreak::HigherHeadroom}));
  static_assert(declared_in_order({RouteAllocation::Equal, RouteAllocation::Weighted, RouteAllocation::Optimized}));
  static_assert(declared_in_order({SourceAllocation::Equal, SourceAllocation::Optimized}));
  static_assert(declared_in_order({WeightBasis::Capacity, WeightBasis::Explicit, WeightBasis::ResidualSnapshot}));
  static_assert(declared_in_order({Admission::Lossless, Admission::Lossy}));
  static_assert(declared_in_order({Direction::Forward, Direction::Reverse}));
  static_assert(static_cast<int>(TunnelStatus::Up) == 0 && static_cast<int>(TunnelStatus::InsufficientBandwidth) == 1 &&
                static_cast<int>(TunnelStatus::CostCeilingExceeded) == 2 && static_cast<int>(TunnelStatus::Failed) == 3 &&
                static_cast<int>(TunnelStatus::Blackholed) == 4);
  static_assert(static_cast<int>(Termination::Complete) == 0 && static_cast<int>(Termination::NoRoute) == 1 &&
                static_cast<int>(Termination::Exhausted) == 2 && static_cast<int>(Termination::NoActiveWeight) == 3);
  // CentralizedTE names its source allocation: no default construction.
  static_assert(!std::is_default_constructible_v<CentralizedTE>);
  EXPECT_EQ(CentralizedTE(SourceAllocation::Equal).source_allocation, SourceAllocation::Equal);
}
