// Forwarding programs, natively: what every model exports, CentralizedTE with a support, the
// weight bases, ForwardingProgram::from_routes, Pinned over programs, the lossless bound, the loss
// rule inside one program (cycles included), replay under a failure mask, stale forwarding with
// black holes, the operations evaluate, replace and forward with their epochs and flags, the
// per-entry records, the by_tunnel rows and their 128-bit counters, the evaluators' offer, and
// common Pinned setups written as programs. Each test checks hand-derived values on small
// graphs whose EdgeIds equal their input positions (asserted); none is read from a build. Then the storage of the new layouts (a restricted plan, a forwarding view, an epoch)
// against borrowing copies, and a property test of random tunnel programs against the
// independent reference of test_utils.hpp (program_reference).
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "netgraph/core/shortest_path_cache.hpp"
#include "netgraph/core/detail/forwarding.hpp"
#include "netgraph/core/detail/route_kernels.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/forwarding_program.hpp"
#include "state_comparators.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

constexpr RouteAllocation EB = RouteAllocation::Equal;
constexpr RouteAllocation W = RouteAllocation::Weighted;
constexpr RouteAllocation MF = RouteAllocation::Optimized;
constexpr Admission LOSSY = Admission::Lossy;
constexpr Admission LOSSLESS = Admission::Lossless;
constexpr TunnelStatus UP = TunnelStatus::Up;
constexpr TunnelStatus FAILED = TunnelStatus::Failed;
constexpr TunnelStatus BLACKHOLED = TunnelStatus::Blackholed;
const std::optional<Amount> INF;  // +inf

using Spec = std::tuple<NodeId, std::vector<EdgeId>, Amount>;  // (head, path edges, weight)
using Listing = std::vector<Spec>;
using Row = std::tuple<std::int32_t, NodeId, TunnelStatus, Amount, Amount, Amount>;
using Rows = std::vector<Row>;  // (index, head, status, assigned, delivered, dropped)
using Drops = std::vector<std::pair<EdgeId, Amount>>;
using Amounts = std::vector<Amount>;

Counter128 c(Amount x) { return Counter128 {0, static_cast<std::uint64_t>(x)}; }
Amount lo(const Counter128& x) {
  EXPECT_EQ(x.hi, 0U);
  return static_cast<Amount>(x.lo);
}
Amount lo(const std::optional<Counter128>& x) {
  EXPECT_TRUE(x.has_value());
  return x.has_value() ? lo(*x) : -1;
}

// A graph whose rows are already in EdgeId order (source, cost, destination, input).
std::shared_ptr<const Graph> graph(std::int32_t n, const std::vector<InEdge>& rows) {
  auto g = shared(make_graph(n, rows));
  for (std::size_t i = 0; i < rows.size(); ++i) EXPECT_EQ(g->input_positions()[i], static_cast<std::int32_t>(i));
  return g;
}
std::shared_ptr<const Graph> par2(Amount a, Amount b) { return graph(2, {{0, 1, a, 1}, {0, 1, b, 1}}); }

const std::vector<InEdge> kOne = {{0, 1, 10, 1}};
const std::vector<InEdge> kB1 = {{0, 3, 2, 1}, {1, 2, 10, 1}, {2, 3, 10, 2}};
const std::vector<InEdge> kB2 = {{0, 1, 4, 1}, {0, 2, 2, 1}, {1, 2, 10, 1}};
const std::vector<InEdge> kES = {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 2, 1}, {2, 3, 10, 2}};
const std::vector<InEdge> kG1U = {{0, 1, 2, 1}, {0, 2, 10, 1}, {1, 3, 2, 1}, {2, 3, 10, 1}};
const std::vector<InEdge> kPinloss = {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 2, 1}, {2, 3, 10, 1}};
const std::vector<InEdge> kDiam = {{0, 1, 100, 1}, {0, 2, 100, 1}, {1, 3, 100, 1}, {2, 3, 100, 1}};
const std::vector<InEdge> kPool2 = {{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 10, 1}, {2, 3, 10, 1}};
const std::vector<InEdge> kRel = {{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 10, 1}};
const std::vector<InEdge> kTun = {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 10, 1}, {2, 3, 10, 1}};
const std::vector<InEdge> kMB = {{0, 1, 10, 1}, {0, 2, 10, 1}, {0, 3, 10, 1},
                                 {1, 4, 10, 1}, {2, 4, 10, 1}, {3, 4, 10, 1}};
const std::vector<InEdge> kSurv = {{0, 3, 3, 1}, {1, 3, 1, 1}, {2, 3, 5, 1}};
const std::vector<InEdge> kRB = {{0, 2, 5, 1}, {1, 2, 5, 1}, {1, 2, 5, 2}};
const std::vector<InEdge> kInc = {{0, 2, 1, 1}, {1, 2, 3, 1}};
const std::vector<InEdge> kShare = {{0, 1, 4, 1}, {1, 2, 10, 1}, {1, 2, 10, 1}};
const std::vector<InEdge> kHier = {{0, 2, 1000, 1}, {1, 2, 1000, 1}, {2, 3, 1, 1}};

// CYC(x, y): tunnel 0 = [e6, e0, e1, e3, e5] meets X (e0) then Y (e3); tunnel 1 = [e7, e3, e4,
// e0, e2] meets Y then X; 4 -> 5.
std::shared_ptr<const Graph> cyc(Amount x, Amount y) {
  return graph(6, {{0, 1, x, 1},
                   {1, 2, 1000, 1},
                   {1, 5, 1000, 1},
                   {2, 3, y, 1},
                   {3, 0, 1000, 1},
                   {3, 5, 1000, 1},
                   {4, 0, 1000, 1},
                   {4, 2, 1000, 1}});
}
const std::vector<EdgeId> kCycT0 = {6, 0, 1, 3, 5};
const std::vector<EdgeId> kCycT1 = {7, 3, 4, 0, 2};

// CYCS(x, y, z): CYC with Y split in series, e3: 2 -> 6 cap y, e8: 6 -> 3 cap z.
std::shared_ptr<const Graph> cycs(Amount x, Amount y, Amount z) {
  return graph(7, {{0, 1, x, 1},
                   {1, 2, 1000, 1},
                   {1, 5, 1000, 1},
                   {2, 6, y, 1},
                   {3, 0, 1000, 1},
                   {3, 5, 1000, 1},
                   {4, 0, 1000, 1},
                   {4, 2, 1000, 1},
                   {6, 3, z, 1}});
}
const std::vector<EdgeId> kCycsT0 = {6, 0, 1, 3, 8, 5};
const std::vector<EdgeId> kCycsT1 = {7, 3, 8, 4, 0, 2};

std::shared_ptr<const Path> path(const Graph& g, const std::vector<EdgeId>& edges) {
  std::vector<std::vector<EdgeId>> hops;
  for (EdgeId e : edges) hops.push_back({e});
  return std::make_shared<const Path>(g, hops);
}
RouteDefinition tunnel_route(const Graph& g, NodeId head, const std::vector<EdgeId>& edges, Amount weight = 1) {
  return RouteDefinition {head, path(g, edges), weight};
}
RouteDefinition dag_route(const Graph& g, NodeId head, const std::vector<EdgeId>& edges, Amount weight = 1) {
  return RouteDefinition {head, dag_of(g, edges), weight};
}
std::shared_ptr<const ForwardingProgram> program(const Graph& g, const std::vector<NodeId>& s,
                                                 const std::vector<NodeId>& d, const std::vector<RouteDefinition>& routes,
                                                 const RouteOptions& options = {}) {
  return ForwardingProgram::from_routes(g, s, d, routes, options);
}
// A supplied program of tunnels.
std::shared_ptr<const ForwardingProgram> tunnels(const Graph& g, const std::vector<NodeId>& s,
                                                 const std::vector<NodeId>& d, const std::vector<Spec>& specs) {
  std::vector<RouteDefinition> routes;
  for (const auto& [h, es, w] : specs) routes.push_back(tunnel_route(g, h, es, w));
  return program(g, s, d, routes);
}

std::shared_ptr<const Eligibility> no_nodes(const Graph& g, const std::vector<NodeId>& nodes) {
  return std::make_shared<const Eligibility>(g.eligibility(nodes, {}));
}
std::shared_ptr<const Eligibility> no_edges(const Graph& g, const std::vector<EdgeId>& edges) {
  return std::make_shared<const Eligibility>(g.eligibility({}, edges));
}

// A policy with what it points to: the ledger and the Eligibility outlive it (declared first).
struct Pol {
  std::shared_ptr<FlowLedger> ledger;
  std::shared_ptr<const Eligibility> elig;
  std::unique_ptr<FlowPolicy> p;
  FlowPolicy* operator->() const { return p.get(); }
};

Pol make(const std::shared_ptr<const Graph>& g, const std::vector<NodeId>& s, const std::vector<NodeId>& d,
         FlowPolicyConfig cfg, std::shared_ptr<const Eligibility> elig = nullptr,
         std::shared_ptr<FlowLedger> ledger = nullptr) {
  Pol out {ledger ? std::move(ledger) : std::make_shared<FlowLedger>(g), std::move(elig), nullptr};
  out.p = std::make_unique<FlowPolicy>(*out.ledger, s, d, std::move(cfg), out.elig.get());
  return out;
}

FlowPolicyConfig cost_only(RouteAllocation placement = EB, Admission admission = LOSSLESS) {
  return FlowPolicyConfig {.routing = shortest_paths(placement), .admission = admission};
}
FlowPolicyConfig weighted(Admission admission, WeightBasis basis, std::vector<Amount> edge_weights = {}) {
  ShortestPaths sp = shortest_paths(W);
  sp.weight_basis = basis;
  if (basis == WeightBasis::Explicit) sp.edge_weights = std::make_shared<const std::vector<Amount>>(std::move(edge_weights));
  return FlowPolicyConfig {.routing = sp, .admission = admission};
}
FlowPolicyConfig global(std::optional<SourceAllocation> source_allocation = std::nullopt, std::shared_ptr<const PredecessorDAG> support = nullptr) {
  return FlowPolicyConfig {.routing = centralized(source_allocation.value_or(SourceAllocation::Optimized), std::move(support))};
}
FlowPolicyConfig pinned_config(std::shared_ptr<const ForwardingProgram> prog, Admission admission = LOSSLESS,
                               std::shared_ptr<const Eligibility> black_holes = nullptr) {
  return FlowPolicyConfig {.routing = Pinned {std::move(prog), std::move(black_holes)}, .admission = admission};
}
FlowPolicyConfig capacity_aware() { return FlowPolicyConfig {.routing = DistributedTE {TunnelCount {1}}}; }

Pol pinned(const std::shared_ptr<const Graph>& g, const std::shared_ptr<const ForwardingProgram>& prog,
           Admission admission = LOSSLESS, std::shared_ptr<const Eligibility> elig = nullptr,
           std::shared_ptr<const Eligibility> black_holes = nullptr, std::shared_ptr<FlowLedger> ledger = nullptr) {
  return make(g, vec(prog->sources()), vec(prog->destinations()), pinned_config(prog, admission, std::move(black_holes)),
              std::move(elig), std::move(ledger));
}

Rows rows(const PlacementResult& r) {
  Rows out;
  for (const auto& t : r.by_tunnel) out.emplace_back(t.entry_id, t.head, t.status, lo(t.assigned), t.delivered, lo(t.total_edge_drops));
  return out;
}
Listing listing(const ForwardingProgram& p) {
  Listing out;
  for (const auto& t : p.tunnels()) out.emplace_back(t.head, t.edges.value_or(std::vector<EdgeId> {}), t.weight);
  return out;
}
std::vector<std::optional<std::int32_t>> entries(const FlowPolicy& p) {
  std::vector<std::optional<std::int32_t>> out;
  for (const auto& r : p.records()) out.push_back(r.entry_id);
  return out;
}
Amounts delivered(const FlowPolicy& p) {
  Amounts out;
  for (const auto& r : p.records()) out.push_back(r.delivered);
  return out;
}
std::vector<std::int64_t> ids(const FlowPolicy& p) {
  std::vector<std::int64_t> out;
  for (const auto& r : p.records()) out.push_back(r.flow_id.second);
  return out;
}

const OutputOptions kEdges {.edge_flows = true};

EvaluationOptions query(FlowPolicyConfig cfg, std::optional<std::optional<Amount>> offered = std::nullopt) {
  EvaluationOptions q {std::move(cfg)};
  q.offer = offered;
  return q;
}
std::optional<std::optional<Amount>> offer(Amount x) { return std::optional<std::optional<Amount>> {std::in_place, x}; }
std::optional<std::optional<Amount>> offer_inf() { return std::optional<std::optional<Amount>> {std::in_place, std::nullopt}; }

std::vector<std::pair<EdgeId, Amount>> losses(const SensitivityResult& s) {
  std::vector<std::pair<EdgeId, Amount>> out;
  for (const auto& e : s.entries) out.emplace_back(e.edge, e.delivery_reduction);
  return out;
}

const std::string kFinite = "LOSSY forwarding requires a finite offer";
const std::string kNoProgram = "forward requires an installed program";
const std::string kPeriodTwo = "lossy forwarding iteration entered a period-2 cycle";
const std::string kRoundLimit = "lossy forwarding iteration did not reach a fixed point within the round limit";

} // namespace

// ---------------------------------------------------------------- exported programs

// B1, CentralizedTE optimized, place(12).
TEST(ProgramExport, GlobalOptimized) {
  const auto g = graph(4, kB1);
  auto p = make(g, {0, 1}, {3}, global());
  const auto r = p->place(12);
  ASSERT_TRUE(r.program);
  const auto& prog = *r.program;
  EXPECT_EQ(listing(prog), (Listing {{0, {0}, 2}, {1, {1, 2}, 10}}));
  std::vector<std::tuple<std::int32_t, std::int32_t, NodeId, Cost>> fields;
  for (const auto& t : prog.tunnels()) fields.emplace_back(t.entry_id, t.local_index, *t.destination, *t.path_cost);
  EXPECT_EQ(fields, (std::vector<std::tuple<std::int32_t, std::int32_t, NodeId, Cost>> {{0, 0, 3, 1}, {1, 0, 3, 3}}));
  ASSERT_TRUE(prog.source_weights().has_value());
  EXPECT_EQ(*prog.source_weights(), (std::vector<std::pair<NodeId, Amount>> {{0, 2}, {1, 10}}));
  EXPECT_EQ(vec(prog.sources()), v({0, 1}));
  EXPECT_EQ(vec(prog.destinations()), v({3}));
  EXPECT_EQ(prog.graph_uid(), g->uid());
  EXPECT_EQ(entries(*p.p), (std::vector<std::optional<std::int32_t>> {std::nullopt, std::nullopt}));
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 2, 2, 0}, {1, 1, UP, 10, 10, 0}}));
  EXPECT_EQ(p->program(), r.program);
}

// The equal source allocation admits 4 and keeps no source weights.
TEST(ProgramExport, GlobalEqual) {
  const auto g = graph(4, kB1);
  auto p = make(g, {0, 1}, {3}, global(SourceAllocation::Equal));
  const auto r = p->place(12);
  EXPECT_EQ(r.call_delivered, 4);
  EXPECT_EQ(listing(*r.program), (Listing {{0, {0}, 2}, {1, {1, 2}, 2}}));
  EXPECT_FALSE(r.program->source_weights().has_value());
}

// One source is decomposed too.
TEST(ProgramExport, SingleSourceDecomposed) {
  const auto g = graph(4, kES);
  auto p = make(g, {0}, {3}, global());
  const auto r = p->place(12);
  EXPECT_EQ(listing(*r.program), (Listing {{0, {0, 2}, 2}, {0, {1, 3}, 10}}));
  std::vector<std::pair<std::int32_t, Cost>> fields;
  for (const auto& t : r.program->tunnels()) fields.emplace_back(t.local_index, *t.path_cost);
  EXPECT_EQ(fields, (std::vector<std::pair<std::int32_t, Cost>> {{0, 2}, {1, 3}}));
  EXPECT_EQ(r.total_cost, c(34));
  EXPECT_EQ(*r.program->source_weights(), (std::vector<std::pair<NodeId, Amount>> {{0, 12}}));
}

// ShortestPaths Optimized exports its controller's paths.
TEST(ProgramExport, CostOnlyMaxFlow) {
  const auto g = graph(4, kES);
  auto p = make(g, {0}, {3}, cost_only(MF));
  const auto r = p->place(12);
  EXPECT_EQ(r.call_delivered, 2);
  EXPECT_EQ(listing(*r.program), (Listing {{0, {0, 2}, 2}}));
  EXPECT_EQ(*r.program->source_weights(), (std::vector<std::pair<NodeId, Amount>> {{0, 2}}));
}

// A ShortestPaths split program lists no tunnel.
TEST(ProgramExport, CostOnlySplit) {
  const auto g = graph(4, kPool2);
  auto p = make(g, {0, 1}, {3}, cost_only());
  const auto r = p->place(2);
  ASSERT_TRUE(r.program);
  EXPECT_TRUE(r.program->tunnels().empty());
  EXPECT_FALSE(r.program->source_weights().has_value());
  EXPECT_EQ(vec(r.program->sources()), v({0, 1}));
  EXPECT_EQ(vec(r.program->destinations()), v({3}));
  EXPECT_TRUE(r.by_tunnel.empty());
}

// Every roster source keeps its source weight, zero included.
TEST(ProgramExport, ZeroSupplyKept) {
  const auto g = graph(3, {{0, 2, 5, 1}, {1, 2, 5, 1}});
  auto p = make(g, {0, 1}, {2}, global());
  const auto r = p->place(1);
  EXPECT_EQ(listing(*r.program), (Listing {{0, {0}, 1}}));
  EXPECT_EQ(*r.program->source_weights(), (std::vector<std::pair<NodeId, Amount>> {{0, 1}, {1, 0}}));
}

// No program before the first positive call; a zero call reports none.
TEST(ProgramExport, NoneBeforePreparation) {
  const auto g = graph(3, {{0, 2, 5, 1}, {1, 2, 5, 1}});
  for (const auto& cfg : {global(), cost_only()}) {
    auto p = make(g, {0, 1}, {2}, cfg);
    EXPECT_EQ(p->program(), nullptr);
    const auto z = p->place(0);
    EXPECT_EQ(z.program, nullptr);
    EXPECT_TRUE(z.by_tunnel.empty());
    EXPECT_EQ(p->program(), nullptr);
    (void)p->place(3);
    EXPECT_NE(p->program(), nullptr);
  }
}

// B1's optimized program forwarded by Pinned reproduces the solve.
TEST(ProgramExport, ReplayReproducesTheSolve) {
  const auto g = graph(4, kB1);
  auto solver = make(g, {0, 1}, {3}, global());
  const auto prog = solver->place(12).program;
  auto p = pinned(g, prog);
  const auto r = p->place(12, kEdges);
  EXPECT_EQ(r.call_delivered, 12);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_EQ(*r.edge_flows, (Amounts {2, 10, 10}));
  EXPECT_EQ(pinned(g, prog)->place(INF).call_delivered, 12);
}

// The producer's mask is saved: no replay brings source 1 back.
TEST(ProgramExport, MaskSaved) {
  const auto g = graph(4, kB1);
  auto masked = make(g, {0, 1}, {3}, global(), no_nodes(*g, {1}));
  const auto r = masked->place(12);
  EXPECT_EQ(r.call_delivered, 2);
  EXPECT_EQ(listing(*r.program), (Listing {{0, {0}, 2}}));
  const auto lossy = pinned(g, r.program, LOSSY)->replace(12);
  EXPECT_EQ(lossy.call_delivered, 2);
  ASSERT_EQ(lossy.by_source.size(), 2U);
  EXPECT_EQ(lo(lossy.by_source[0].offered), 12);
  EXPECT_EQ(lo(lossy.by_source[1].offered), 0);
}

// An empty roster installs an empty program.
TEST(ProgramExport, EmptyRoster) {
  const auto g = graph(2, kOne);
  auto p = make(g, {0}, {1}, cost_only(), no_nodes(*g, {0}));
  const auto r = p->place(1);
  EXPECT_EQ(r.termination, Termination::NoRoute);
  ASSERT_TRUE(r.program);
  EXPECT_TRUE(r.program->tunnels().empty());
  EXPECT_EQ(vec(r.program->sources()), v({0}));
  EXPECT_FALSE(r.routes_changed);
  EXPECT_FALSE(r.weights_changed);
}

// The policy's program is the configured one.
TEST(PinnedProgram, ProgramIsTheConfiguredOne) {
  const auto g = graph(4, kB1);
  auto solver = make(g, {0, 1}, {3}, global());
  const auto prog = solver->place(12).program;
  auto p = pinned(g, prog);
  const auto r = p->place(1);
  EXPECT_EQ(r.program, prog);
  EXPECT_EQ(p->program(), prog);
  EXPECT_EQ(listing(*r.program), listing(*prog));
}

// ---------------------------------------------------------------- CentralizedTE with a support

// Only the support carries traffic; no cut, no curve.
TEST(GlobalSupport, OnlyTheSupport) {
  const auto g = graph(4, kES);
  auto p = make(g, {0}, {3}, global(std::nullopt, dag_of(*g, {1, 3})));
  const auto r = p->place(12, OutputOptions {.min_cut = true});
  EXPECT_EQ(r.call_delivered, 10);
  EXPECT_EQ(r.total_cost, c(30));
  EXPECT_FALSE(r.min_cut.has_value());
  EXPECT_FALSE(r.marginal_cost_curve.has_value());
  EXPECT_EQ(listing(*r.program), (Listing {{0, {1, 3}, 10}}));
  const auto full = make(g, {0}, {3}, global())->place(12);
  EXPECT_EQ(full.call_delivered, 12);
  EXPECT_EQ(full.total_cost, c(34));
}

// e0 leads nowhere inside the support: trimmed.
TEST(GlobalSupport, Trimmed) {
  const auto g = graph(4, kES);
  const auto r = make(g, {0}, {3}, global(std::nullopt, dag_of(*g, {0, 1, 3})))->place(12, kEdges);
  EXPECT_EQ(r.call_delivered, 10);
  EXPECT_EQ(*r.edge_flows, (Amounts {0, 10, 0, 10}));
}

// The mask removes the support: an empty roster.
TEST(GlobalSupport, MaskRemovesSupport) {
  const auto g = graph(4, kES);
  const auto r = make(g, {0}, {3}, global(std::nullopt, dag_of(*g, {1, 3})), no_edges(*g, {3}))->place(12);
  EXPECT_EQ(r.termination, Termination::NoRoute);
  EXPECT_EQ(r.call_delivered, 0);
  ASSERT_TRUE(r.program);
  EXPECT_TRUE(r.program->tunnels().empty());
}

// The equal source allocation on a support: RR(4, 2) fits, RR(5, 2) does not.
TEST(GlobalSupport, EqualOrigin) {
  const auto g = graph(3, kB2);
  const auto r = make(g, {0, 1}, {2}, global(SourceAllocation::Equal, dag_of(*g, {1, 2})))->place(12);
  EXPECT_EQ(r.call_delivered, 4);
  EXPECT_EQ(r.total_cost, c(4));
}

// Every allowed link of a support is a candidate, no shortcut.
TEST(GlobalSupport, Sensitivity) {
  const auto g = graph(4, kES);
  const auto s = evaluate_edge_failures(*g, v({0}), v({3}), query(global(std::nullopt, dag_of(*g, {1, 3}))), 1);
  EXPECT_EQ(s.baseline.call_delivered, 10);
  EXPECT_EQ(losses(s), (std::vector<std::pair<EdgeId, Amount>> {{0, 0}, {1, 10}, {2, 0}, {3, 10}}));
}

// total_cost 2^70 (high half 64, low half 0).
TEST(GlobalSupport, TotalCostBeyond64Bits) {
  const auto g = graph(2, {{0, 1, Amount {1} << 40, Cost {1} << 30}});
  const auto r = make(g, {0}, {1}, global(std::nullopt, dag_of(*g, {0})))->place(Amount {1} << 40);
  EXPECT_EQ(r.call_delivered, Amount {1} << 40);
  ASSERT_TRUE(r.total_cost.has_value());
  EXPECT_EQ(r.total_cost->hi, 64U);
  EXPECT_EQ(r.total_cost->lo, 0U);
}

// A support of another graph.
TEST(GlobalSupport, Ownership) {
  const auto g = graph(2, {{0, 1, 1, 1}});
  const auto other = graph(2, {{0, 1, 1, 1}});
  FlowLedger fg(g);
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), global(std::nullopt, dag_of(*other, {0}))); }),
            "support belongs to another graph");
}

// ---------------------------------------------------------------- weight bases

// Configured weights; CAPACITY weights give (6, 6); the lossless bound 15.
TEST(WeightBases, Explicit) {
  const auto g = par2(10, 10);
  EXPECT_EQ(*make(g, {0}, {1}, weighted(LOSSY, WeightBasis::Explicit, {1, 2}))->place(12, kEdges).edge_flows,
            (Amounts {4, 8}));
  EXPECT_EQ(*make(g, {0}, {1}, cost_only(W, LOSSY))->place(12, kEdges).edge_flows, (Amounts {6, 6}));
  const auto ll = weighted(LOSSLESS, WeightBasis::Explicit, {1, 2});
  EXPECT_EQ(make(g, {0}, {1}, ll)->place(INF).call_delivered, 15);
  EXPECT_EQ(evaluate_demand(*g, v({0}), v({1}), query(ll)).call_delivered, 15);
}

// Legal weights below 2^62 each may sum to 2^64 or more; the sum is exact. Eight
// equal weights of 2^61 (sum 2^64) and of 2^62 - 1 split 8 units as RR(8, 8) = one each; unequal
// weights (2^62 - 1, 2^61, 2^61 + 1, 1), W = 2^63 + 1: H(9, w) = (5, 2, 2, 0). Through the helper,
// tunnels and explicit WCMP, both admissions.
TEST(WeightBases, WeightSumsBeyond64Bits) {
  for (const Amount w : {Amount {1} << 61, (Amount {1} << 62) - 1}) {
    Amounts out(8, 0);
    detail::highest_averages(8, Amounts(8, w), out);
    EXPECT_EQ(out, Amounts(8, 1)) << w;
    std::vector<InEdge> rows(8, InEdge {0, 1, 100, 1});
    const auto g = graph(2, rows);
    std::vector<std::tuple<NodeId, std::vector<EdgeId>, Amount>> specs;
    for (EdgeId e = 0; e < 8; ++e) specs.emplace_back(0, std::vector<EdgeId> {e}, w);
    for (const auto admission : {LOSSLESS, LOSSY}) {
      auto t = pinned(g, tunnels(*g, {0}, {1}, specs), admission);
      const auto r = t->replace(8, std::nullopt, kEdges);
      EXPECT_EQ(r.call_delivered, 8) << w;
      EXPECT_EQ(*r.edge_flows, Amounts(8, 1)) << w;
      EXPECT_TRUE(r.source_drops.is_zero()) << w;
      const auto x = make(g, {0}, {1}, weighted(admission, WeightBasis::Explicit, Amounts(8, w)))->replace(8, std::nullopt, kEdges);
      EXPECT_EQ(x.call_delivered, 8) << w;
      EXPECT_EQ(*x.edge_flows, Amounts(8, 1)) << w;
    }
  }
  const Amounts unequal = {(Amount {1} << 62) - 1, Amount {1} << 61, (Amount {1} << 61) + 1, 1};
  Amounts out(4, 0);
  detail::highest_averages(9, unequal, out);
  EXPECT_EQ(out, (Amounts {5, 2, 2, 0}));
  // The equal-weight control below 2^64: the same split.
  detail::highest_averages(8, Amounts(8, 3), out = Amounts(8, 0));
  EXPECT_EQ(out, Amounts(8, 1));
}

// Snapshot weights over the baseline (2, 10).
TEST(WeightBases, Snapshot) {
  const auto g = par2(10, 10);
  const Amounts base = {2, 10};
  auto ledger = [&] { return std::make_shared<FlowLedger>(g, base); };
  EXPECT_EQ(make(g, {0}, {1}, cost_only(W), nullptr, ledger())->place(12).call_delivered, 4);
  EXPECT_EQ(make(g, {0}, {1}, weighted(LOSSLESS, WeightBasis::ResidualSnapshot), nullptr, ledger())->place(12).call_delivered,
            12);
}

// The tie rule never chooses weights (CAPACITY is the default).
TEST(WeightBases, TieRuleIsNotAWeight) {
  const auto g = par2(10, 10);
  auto cfg = cost_only(W);
  auto p = make(g, {0}, {1}, cfg, nullptr, std::make_shared<FlowLedger>(g, Amounts {2, 10}));
  EXPECT_EQ(p->place(12).call_delivered, 4);
  EXPECT_EQ(std::get<ShortestPaths>(p->config().routing).weight_basis, WeightBasis::Capacity);
}

// A positive weight on an edge of capacity 0 becomes 0: all 7 on e1, no drop.
TEST(WeightBases, ExplicitOnZeroCapacity) {
  const auto g = par2(0, 10);
  const auto r = make(g, {0}, {1}, weighted(LOSSY, WeightBasis::Explicit, {5, 1}))->place(7);
  EXPECT_EQ(r.call_delivered, 7);
  EXPECT_TRUE(r.call_edge_drops.empty());
}

// A zero weight prunes its edge.
TEST(WeightBases, ExplicitZeroPrunes) {
  const auto g = par2(10, 10);
  EXPECT_EQ(make(g, {0}, {1}, weighted(LOSSLESS, WeightBasis::Explicit, {0, 3}))->place(12).call_delivered, 10);
  const auto r = make(g, {0}, {1}, weighted(LOSSY, WeightBasis::Explicit, {0, 3}))->place(12);
  EXPECT_EQ(r.call_delivered, 10);
  EXPECT_EQ(r.call_edge_drops, (Drops {{1, 2}}));
}

// Every weight 0: NO_ACTIVE_WEIGHT; Lossy loses the offer at the source.
TEST(WeightBases, AllZeroIsNoActiveWeight) {
  const auto g = par2(10, 10);
  for (const auto& [adm, loss] : {std::pair {LOSSY, Amount {5}}, std::pair {LOSSLESS, Amount {0}}}) {
    const auto r = make(g, {0}, {1}, weighted(adm, WeightBasis::Explicit, {0, 0}))->place(5);
    EXPECT_EQ(r.termination, Termination::NoActiveWeight);
    EXPECT_EQ(r.call_delivered, 0);
    EXPECT_EQ(r.source_drops, c(loss));
  }
}

// Weighted with HIGHER_HEADROOM is accepted (ShortestPaths): the single path is e1.
TEST(WeightBases, WeightedWithHeadroom) {
  const auto g = par2(10, 10);
  ShortestPaths headroom = shortest_paths(W);
  headroom.tie_break = RouteTieBreak::HigherHeadroom;
  headroom.multipath = false;
  const FlowPolicyConfig cfg {.routing = headroom};
  auto p = make(g, {0}, {1}, cfg, nullptr, std::make_shared<FlowLedger>(g, Amounts {2, 10}));
  EXPECT_EQ(p->place(12).call_delivered, 10);
}

// The snapshot is frozen for the epoch: call 2 splits 6 as H(6, (2, 10)) = (1, 5), both
// dropped.
TEST(WeightBases, SnapshotFrozenForTheEpoch) {
  const auto g = par2(10, 10);
  auto p = make(g, {0}, {1}, weighted(LOSSY, WeightBasis::ResidualSnapshot), nullptr,
                std::make_shared<FlowLedger>(g, Amounts {2, 10}));
  EXPECT_EQ(*p->place(12, kEdges).edge_flows, (Amounts {2, 10}));
  const auto r = p->place(6);
  EXPECT_EQ(r.call_delivered, 0);
  EXPECT_EQ(r.call_edge_drops, (Drops {{0, 1}, {1, 5}}));
  EXPECT_EQ(r.call_source_drops, 0);
}

// The texts of the weight fields (only ShortestPaths has them).
TEST(WeightBases, MatrixTexts) {
  const auto g = par2(10, 10);
  FlowLedger fg(g);
  auto with = [](FlowPolicyConfig cfg, std::optional<WeightBasis> w, std::optional<std::vector<Amount>> ew) {
    auto& sp = std::get<ShortestPaths>(cfg.routing);
    sp.weight_basis = w;
    if (ew.has_value()) sp.edge_weights = std::make_shared<const std::vector<Amount>>(*ew);
    return cfg;
  };
  const std::vector<std::pair<FlowPolicyConfig, std::string>> cases = {
      {with(cost_only(), WeightBasis::Capacity, std::nullopt), "weight_basis requires allocation=WEIGHTED"},
      {with(cost_only(), std::nullopt, Amounts {1, 1}), "edge_weights require weight_basis=EXPLICIT"},
      {with(cost_only(W), WeightBasis::Explicit, std::nullopt), "weight_basis=EXPLICIT requires edge_weights"},
      {with(cost_only(W), WeightBasis::Capacity, Amounts {1, 1}), "edge_weights require weight_basis=EXPLICIT"},
      {with(cost_only(W), WeightBasis::Explicit, Amounts {1, 2, 3}), "edge_weights: expected 2 entries"},
      {with(cost_only(W), WeightBasis::Explicit, Amounts {1, kL62}), "edge_weights[1] must be an integer in [0, 2^62)"},
      {with(cost_only(W), WeightBasis::Explicit, Amounts {-1, 1}), "edge_weights[0] must be an integer in [0, 2^62)"},
  };
  for (const auto& [cfg, text] : cases) {
    EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), cfg); }), text);
  }
  // The largest entry of the domain is accepted.
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), weighted(LOSSLESS, WeightBasis::Explicit, {1, kL62 - 1})); }),
            "");
}

// One configuration's edge_weights array is shared by every policy built from it.
TEST(WeightBases, SharedArray) {
  const auto g = par2(10, 10);
  const auto cfg = weighted(LOSSLESS, WeightBasis::Explicit, {1, 2});
  FlowLedger fg(g);
  FlowPolicy a(fg, v({0}), v({1}), cfg);
  FlowPolicy b(fg, v({0}), v({1}), cfg);
  const auto weights_of = [](const FlowPolicyConfig& c) { return std::get<ShortestPaths>(c.routing).edge_weights.get(); };
  const auto* buffer = weights_of(cfg);
  EXPECT_EQ(weights_of(detail::PolicyAccess::engine(a).config), buffer);
  EXPECT_EQ(weights_of(detail::PolicyAccess::engine(b).config), buffer);
  EXPECT_EQ(weights_of(a.config()), buffer);
  // Control: an equal array built separately is another buffer.
  FlowPolicy c2(fg, v({0}), v({1}), weighted(LOSSLESS, WeightBasis::Explicit, {1, 2}));
  EXPECT_NE(weights_of(detail::PolicyAccess::engine(c2).config), buffer);
}

// ---------------------------------------------------------------- from_routes

// The validation texts in order (graph PAR2(4, 8), 0 -> 1 unless stated).
TEST(FromRoutes, Texts) {
  const auto g = par2(4, 8);
  const auto p0 = tunnel_route(*g, 0, {0});
  const RouteDefinition d01 = dag_route(*g, 0, {0, 1});
  auto err = [&](const std::vector<NodeId>& s, const std::vector<NodeId>& d, const std::vector<RouteDefinition>& routes,
                 const RouteOptions& o = {}, const Graph* on = nullptr) {
    return error_of([&] { (void)ForwardingProgram::from_routes(on ? *on : *g, s, d, routes, o); });
  };
  EXPECT_EQ(err({}, {1}, {p0}), "sources: the node set is empty");
  EXPECT_EQ(err({0}, {0}, {p0}), "node 0 is in both endpoint sets: sources and destinations must not overlap");
  EXPECT_EQ(err({0}, {1}, {}), "routes must not be empty");
  EXPECT_EQ(err({0}, {1}, {p0}, RouteOptions {.allocation = EB}), "allocation applies only to PredecessorDAG routes");
  EXPECT_EQ(err({0}, {1}, {p0}, RouteOptions {.weight_basis = WeightBasis::Capacity}),
            "weight_basis applies only to PredecessorDAG routes");
  EXPECT_EQ(err({0}, {1}, {p0}, RouteOptions {.edge_weights = Amounts {1, 1}}),
            "edge_weights applies only to PredecessorDAG routes");
  EXPECT_EQ(err({0}, {1}, {p0}, RouteOptions {.residual = Amounts {1, 1}}), "residual applies only to PredecessorDAG routes");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.allocation = MF}),
            "from_routes allocation must be EQUAL or WEIGHTED");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.weight_basis = WeightBasis::Capacity}), "weight_basis requires allocation=WEIGHTED");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.allocation = W, .weight_basis = WeightBasis::Explicit}),
            "weight_basis=EXPLICIT requires edge_weights");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.allocation = W, .edge_weights = Amounts {1, 1}}),
            "edge_weights require weight_basis=EXPLICIT");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.allocation = W, .weight_basis = WeightBasis::ResidualSnapshot}),
            "weight_basis=RESIDUAL_SNAPSHOT requires residual");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.allocation = W, .residual = Amounts {4, 8}}),
            "residual requires weight_basis=RESIDUAL_SNAPSHOT");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.allocation = W, .weight_basis = WeightBasis::Explicit, .edge_weights = Amounts {1}}),
            "edge_weights: expected 2 entries");
  EXPECT_EQ(err({0}, {1}, {d01},
                RouteOptions {.allocation = W, .weight_basis = WeightBasis::Explicit, .edge_weights = Amounts {1, kL62}}),
            "edge_weights[1] must be an integer in [0, 2^62)");
  EXPECT_EQ(err({0}, {1}, {d01},
                RouteOptions {.allocation = W, .weight_basis = WeightBasis::ResidualSnapshot, .residual = Amounts {5, 8}}),
            "residual[0] = 5 is outside [0, edge_capacity_units[0] = 4]");
  const auto other = par2(4, 8);
  EXPECT_EQ(err({0}, {1}, {tunnel_route(*other, 0, {0})}), "routes[0]: route belongs to another graph");
  EXPECT_EQ(err({0}, {1}, {tunnel_route(*g, 1, {0})}), "routes[0]: head 1 is not a source");
  const RouteDefinition hop2 {0, std::make_shared<const Path>(*g, std::vector<std::vector<EdgeId>> {{0, 1}}), 1};
  EXPECT_EQ(err({0}, {1}, {hop2}), "routes[0]: a tunnel has one edge per hop");
  const auto d = graph(4, kDiam);
  EXPECT_EQ(err({0, 1}, {3}, {tunnel_route(*d, 1, {0, 2})}, {}, d.get()),
            "routes[0]: the path starts at 0, not at head 1");
  EXPECT_EQ(err({0}, {3}, {tunnel_route(*d, 0, {0})}, {}, d.get()),
            "routes[0]: the path ends at 1, not at a destination");
  EXPECT_EQ(err({0}, {1, 3}, {tunnel_route(*d, 0, {0, 2})}, {}, d.get()),
            "routes[0]: the path meets destination 1 before its end");
  EXPECT_EQ(err({0}, {3}, {dag_route(*d, 0, {2})}, {}, d.get()), "routes[0]: no walk from head 0 to a destination");
  for (const Amount w : {kL62, Amount {-1}}) {
    EXPECT_EQ(err({0}, {1}, {tunnel_route(*g, 0, {0}, w)}), "routes[0]: weight must be an integer in [0, 2^62)");
  }
  // The last value of the domain is accepted; a later route is named by its position.
  EXPECT_EQ(err({0}, {1}, {tunnel_route(*g, 0, {0}, kL62 - 1)}), "");
  EXPECT_EQ(err({0}, {1}, {p0, tunnel_route(*g, 0, {1}, -1)}), "routes[1]: weight must be an integer in [0, 2^62)");
}

// The routes in list order, each route's checks in order.
TEST(FromRoutes, Order) {
  const auto g = par2(4, 8);
  EXPECT_EQ(error_of([&] {
              (void)program(*g, {0}, {1}, {tunnel_route(*g, 0, {0}, kL62), tunnel_route(*g, 5, {0})});
            }),
            "routes[0]: weight must be an integer in [0, 2^62)");
  // The weight is checked after the shape: route 0's head text comes first.
  EXPECT_EQ(error_of([&] { (void)program(*g, {0}, {1}, {tunnel_route(*g, 1, {0}, kL62)}); }),
            "routes[0]: head 1 is not a source");
}

TEST(FromRoutes, CheckedPlansAreSharedAndOwned) {
  const auto g = par2(4, 8);
  const std::vector<NodeId> sources {0}, destinations {1};
  const std::vector<RouteDefinition> routes {
    tunnel_route(*g, 0, {0}, 0), dag_route(*g, 0, {0, 1}), dag_route(*g, 0, {0, 1}, 0)};
  const std::vector<RouteOptions> options {
    {}, {.allocation = W},
    {.allocation = W, .weight_basis = WeightBasis::Explicit, .edge_weights = Amounts {1, 3}},
    {.allocation = W, .weight_basis = WeightBasis::ResidualSnapshot, .residual = Amounts {0, 8}},
    {.allocation = W, .weight_basis = WeightBasis::Explicit, .edge_weights = Amounts {0, 0}}};
  const std::vector<Amounts> expected_weights {{}, {4, 8}, {1, 3}, {8}, {}};
  const std::vector<Amounts> expected_carried {{2, 2}, {1, 3}, {1, 3}, {0, 4}, {0, 0}};
  for (std::size_t k = 0; k < options.size(); ++k) {
    std::shared_ptr<const ForwardingProgram> prog;
    std::weak_ptr<const detail::PreparedRoute> retained;
    {
      auto checked = detail::check_route_definitions(*g, sources, destinations, routes, options[k]);
      ASSERT_EQ(checked.split_contexts.size(), 2U);
      std::vector<const detail::PreparedRoute*> plans;
      for (const auto& c : checked.split_contexts) plans.push_back(c.prepared_route.get());
      for (const auto& c : checked.split_contexts) {
        ASSERT_TRUE(c.prepared_route);
        EXPECT_EQ(c.prepared_route->support_edges, (std::vector<EdgeId> {0, 1}));
        EXPECT_EQ(c.prepared_route->weights, expected_weights[k]);
        const std::vector<EdgeId> members = k == 4 ? std::vector<EdgeId> {} : k == 3 ? std::vector<EdgeId> {1} : std::vector<EdgeId> {0, 1};
        EXPECT_EQ(c.prepared_route->forwarding_edges, members);
      }
      retained = checked.split_contexts.front().prepared_route;
      prog = detail::ProgramAccess::from_checked_routes(sources, destinations, routes, std::move(checked));
      const auto& data = detail::ProgramAccess::data(*prog);
      ASSERT_EQ(data.split_contexts.size(), 2U);
      EXPECT_EQ(data.entries.size(), 3U);
      EXPECT_EQ(prog->tunnels().size(), 1U);
      for (std::size_t i = 0; i < plans.size(); ++i) {
        EXPECT_EQ(data.split_contexts[i].prepared_route.get(), plans[i]);
        EXPECT_EQ(data.entries[i + 1].min_route_cost, 1);
      }
    }
    ASSERT_FALSE(retained.expired());  // before any read through a borrowed plan
    {
      FlowLedger ledger(g);
      FlowPolicy p(ledger, sources, destinations, {.routing = Pinned {prog, nullptr}});
      const auto result = p.place(4);
      EXPECT_EQ(result.call_delivered, k == 4 ? 0 : 4);
      EXPECT_EQ(ledger.carried(), expected_carried[k]);
    }
    prog.reset();
    EXPECT_TRUE(retained.expired());
  }
}

// A zero-weight route keeps its index and is validated in full.
TEST(FromRoutes, ZeroWeightKept) {
  const auto g = par2(4, 8);
  const auto prog = program(*g, {0}, {1}, {tunnel_route(*g, 0, {0}, 0), tunnel_route(*g, 0, {1})});
  std::vector<std::pair<std::int32_t, Amount>> got;
  for (const auto& t : prog->tunnels()) got.emplace_back(t.entry_id, t.weight);
  EXPECT_EQ(got, (std::vector<std::pair<std::int32_t, Amount>> {{0, 0}, {1, 1}}));
  EXPECT_EQ(error_of([&] { (void)program(*g, {0}, {1}, {tunnel_route(*g, 1, {0}, 0)}); }),
            "routes[0]: head 1 is not a source");
}

// Tunnels listed with indices 0 and 2 (the split route takes 1); local_index 0 and 2.
TEST(FromRoutes, MixedIndices) {
  const auto g = par2(4, 8);
  const auto prog =
      program(*g, {0}, {1}, {tunnel_route(*g, 0, {1}), dag_route(*g, 0, {0, 1}), tunnel_route(*g, 0, {0})});
  std::vector<std::pair<std::int32_t, std::int32_t>> got;
  for (const auto& t : prog->tunnels()) got.emplace_back(t.entry_id, t.local_index);
  EXPECT_EQ(got, (std::vector<std::pair<std::int32_t, std::int32_t>> {{0, 0}, {2, 2}}));
}

// The endpoints are normalized like a policy's (sorted, duplicates removed).
TEST(FromRoutes, Endpoints) {
  const auto d = graph(4, kDiam);
  const auto prog = program(*d, {1, 0, 1}, {3, 3}, {tunnel_route(*d, 0, {0, 2})});
  EXPECT_EQ(vec(prog->sources()), v({0, 1}));
  EXPECT_EQ(vec(prog->destinations()), v({3}));
}

// ---------------------------------------------------------------- Pinned's configuration

// Pinned's row of the matrix: both admissions, no cache (its other options are its program's; the
// model has no field to put one in).
TEST(PinnedProgram, MatrixTexts) {
  const auto g = par2(4, 8);
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}});
  FlowLedger fg(g);
  ShortestPathCache cache(*g);
  for (const auto a : {LOSSLESS, LOSSY}) {
    EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), pinned_config(prog, a)); }), "");
    EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), pinned_config(prog, a), nullptr, &cache); }),
              "Pinned does not use ShortestPathCache");
  }
  fg.reset();
}

// The program's endpoints must be the policy's.
TEST(PinnedProgram, Endpoints) {
  const auto g = graph(3, {{0, 1, 1, 1}, {0, 2, 1, 1}});
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}});
  FlowLedger fg(g);
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1, 2}), pinned_config(prog)); }),
            "program endpoints do not match policy endpoints");
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1, 1}), pinned_config(prog)); }), "");
}

// The ownership texts, before the configuration matrix.
TEST(PinnedProgram, Ownership) {
  const auto g = par2(4, 8);
  const auto other = par2(4, 8);
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}});
  const auto foreign = tunnels(*other, {0}, {1}, {{0, {0}, 1}});
  FlowLedger fg(g);
  ShortestPathCache own_cache(*g);  // a matrix defect after the ownership
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), pinned_config(foreign)); }),
            "program belongs to another graph");
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), pinned_config(foreign), nullptr, &own_cache); }),
            "program belongs to another graph");
  EXPECT_EQ(error_of([&] {
              FlowPolicy p(fg, v({0}), v({1}), pinned_config(prog, LOSSLESS, no_edges(*other, {0})));
            }),
            "black_holes belongs to another graph");
  const auto mask = other->eligibility({}, std::vector<EdgeId> {0});
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), pinned_config(prog), &mask); }),
            "foreign Eligibility: built for another graph");
  ShortestPathCache cache(*other);
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), pinned_config(prog), nullptr, &cache); }),
            "foreign ShortestPathCache: bound to another graph");
  // The mask is checked before the program.
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), pinned_config(foreign), &mask); }),
            "foreign Eligibility: built for another graph");
}

// The endpoint sets come first.
TEST(PinnedProgram, Order) {
  const auto g = par2(4, 8);
  const auto foreign = tunnels(*par2(4, 8), {0}, {1}, {{0, {0}, 1}});
  FlowLedger fg(g);
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, std::vector<NodeId> {}, v({1}), pinned_config(foreign)); }),
            "sources: the node set is empty");
}

// Several sources and destinations.
TEST(PinnedProgram, SeveralEndpoints) {
  const auto g = graph(4, kRel);
  const auto prog = tunnels(*g, {0, 1}, {3}, {{0, {0, 2}, 1}, {1, {1, 2}, 1}});
  EXPECT_EQ(pinned(g, prog)->place(INF).call_delivered, 10);
  const auto two = graph(4, {{0, 2, 3, 1}, {1, 3, 4, 1}});
  const auto many = tunnels(*two, {0, 1}, {2, 3}, {{0, {0}, 1}, {1, {1}, 1}});
  EXPECT_EQ(pinned(two, many)->place(INF).call_delivered, 6);  // RR(6, 2) = (3, 3) fits, RR(7, 2) = (4, 3) does not
}

// ---------------------------------------------------------------- lossless bounds

namespace {
Amount bound(const std::shared_ptr<const Graph>& g, const std::vector<Spec>& specs, const std::vector<NodeId>& s = {0},
             const std::vector<NodeId>& d = {1}) {
  return pinned(g, tunnels(*g, s, d, specs))->place(INF).call_delivered;
}
} // namespace

// Lossless bounds: duplicate routes, deduplicated routes, a static path filled, a zero first or
// last link, a support's Pinned bound.
TEST(PinnedProgram, Duplicates) {
  EXPECT_EQ(bound(par2(4, 8), {{0, {0}, 1}, {0, {1}, 1}, {0, {1}, 1}}), 12);
  EXPECT_EQ(bound(par2(4, 8), {{0, {0}, 1}, {0, {1}, 1}}), 8);
  EXPECT_EQ(bound(par2(6, 4), {{0, {0}, 1}, {0, {1}, 1}}), 9);
  EXPECT_EQ(bound(par2(0, 4), {{0, {0}, 1}, {0, {1}, 1}}), 0);
  EXPECT_EQ(bound(par2(4, 0), {{0, {0}, 1}, {0, {1}, 1}}), 1);
  EXPECT_EQ(bound(graph(4, kPinloss), {{0, {0, 2}, 1}, {0, {1, 3}, 1}}, {0}, {3}), 4);
  // Records per entry: the duplicates' three records, entries 0, 1, 2.
  const auto g = par2(4, 8);
  auto p = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 1}, {0, {1}, 1}, {0, {1}, 1}}));
  (void)p->place(INF);
  EXPECT_EQ(delivered(*p.p), (Amounts {4, 4, 4}));
}

// H(15, (1, 2)) = (5, 10) fits, H(16) = (5, 11) does not.
TEST(PinnedProgram, WeightedTunnels) { EXPECT_EQ(bound(par2(10, 10), {{0, {0}, 1}, {0, {1}, 2}}), 15); }

// A weight-0 tunnel carries nothing; its record and row exist with zeros.
TEST(PinnedProgram, ZeroWeight) {
  const auto g = par2(4, 8);
  auto p = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 0}, {0, {1}, 1}}));
  const auto r = p->place(INF);
  EXPECT_EQ(r.call_delivered, 8);
  EXPECT_EQ(entries(*p.p), (std::vector<std::optional<std::int32_t>> {0, 1}));
  EXPECT_EQ(delivered(*p.p), (Amounts {0, 8}));
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 0, 0, 0}, {1, 0, UP, 8, 8, 0}}));
  // Created at preparation: a zero-delivery record before any positive placement of it.
  auto q = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 0}, {0, {1}, 1}}));
  (void)q->place(1);
  EXPECT_EQ(delivered(*q.p), (Amounts {0, 1}));
}

// Increments, Lossless and Lossy, and a Lossless increment's refused units.
TEST(PinnedProgram, Increments) {
  const auto g = graph(3, kInc);
  const auto prog = tunnels(*g, {0, 1}, {2}, {{0, {0}, 1}, {1, {1}, 1}});
  auto p = pinned(g, prog);
  EXPECT_EQ(p->place(3).call_delivered, 2);
  EXPECT_EQ(p->place(1).call_delivered, 0);
  auto q = pinned(g, prog, LOSSY);
  const auto r1 = q->place(3);
  EXPECT_EQ(r1.call_delivered, 2);
  EXPECT_EQ(r1.call_edge_drops, (Drops {{0, 1}}));
  const auto r2 = q->place(1);  // RR(4) - RR(3) = (0, 1)
  EXPECT_EQ(r2.call_delivered, 1);
  EXPECT_TRUE(r2.call_edge_drops.empty());
  EXPECT_EQ(lo(r2.by_source[0].offered), 2);
  EXPECT_EQ(lo(r2.by_source[1].offered), 2);
  const auto one = graph(2, kOne);
  auto ledger = std::make_shared<FlowLedger>(one);
  auto other = make(one, {0}, {1}, cost_only(), nullptr, ledger);
  (void)other->place(8);
  auto t = pinned(one, tunnels(*one, {0}, {1}, {{0, {0}, 1}}), LOSSLESS, nullptr, nullptr, ledger);
  EXPECT_EQ(t->place(5).call_delivered, 2);
  other->remove();
  (void)other->place(5);
  EXPECT_EQ(t->place(1).call_delivered, 1);
}

// +inf is the bound, EXHAUSTED; Lossy rejects it.
TEST(PinnedProgram, Infinite) {
  const auto g = par2(6, 4);
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}, {0, {1}, 1}});
  const auto r = pinned(g, prog)->place(INF);
  EXPECT_EQ(r.call_delivered, 9);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  EXPECT_FALSE(r.call_undelivered.has_value());
  EXPECT_FALSE(r.requested.has_value());
  auto lossy = pinned(g, prog, LOSSY);
  EXPECT_EQ(error_of([&] { (void)lossy->place(INF); }), kFinite);
}

// A tunnel and a split route side by side: RR(8, 2) = (4, 4).
TEST(PinnedProgram, MixedEntries) {
  const auto g = par2(4, 8);
  auto p = pinned(g, program(*g, {0}, {1}, {tunnel_route(*g, 0, {0}), dag_route(*g, 0, {1})}));
  const auto r = p->place(INF);
  EXPECT_EQ(r.call_delivered, 8);
  EXPECT_EQ(entries(*p.p), (std::vector<std::optional<std::int32_t>> {0, 1}));
  EXPECT_EQ(delivered(*p.p), (Amounts {4, 4}));
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 4, 4, 0}}));
}

// Split contexts: each identical route splits its own unit, also pooled.
TEST(PinnedProgram, SplitContextsAreSeparate) {
  const auto g = par2(1, 1);
  const auto d = dag_of(*g, {0, 1});
  const auto prog = program(*g, {0}, {1}, {RouteDefinition {0, d, 1}, RouteDefinition {0, d, 1}});
  EXPECT_EQ(pinned(g, prog)->place(INF).call_delivered, 1);
  const auto r = pinned(g, prog, LOSSY)->replace(2, std::nullopt, kEdges);
  EXPECT_EQ(r.call_edge_drops, (Drops {{0, 1}}));
  EXPECT_EQ(*r.edge_flows, (Amounts {1, 0}));
  const auto pooled = make(g, {0}, {1}, cost_only())->place(INF, kEdges);
  EXPECT_EQ(pooled.call_delivered, 2);
  EXPECT_EQ(*pooled.edge_flows, (Amounts {1, 1}));
}

// ---------------------------------------------------------------- losses inside one program

namespace {
struct LossyRun {
  PlacementResult r;
  std::vector<std::tuple<Amount, Amount, Amount>> got;  // per tunnel (assigned, delivered, dropped)
};
LossyRun lossy_rows(const std::shared_ptr<const Graph>& g, const std::vector<NodeId>& s, const std::vector<NodeId>& d,
                    const std::vector<Spec>& specs, Amount offer) {
  LossyRun out {pinned(g, tunnels(*g, s, d, specs), LOSSY)->replace(offer), {}};
  for (const auto& t : out.r.by_tunnel) out.got.emplace_back(lo(t.assigned), t.delivered, lo(t.total_edge_drops));
  return out;
}
using Triples = std::vector<std::tuple<Amount, Amount, Amount>>;
Amounts delivered_of(const LossyRun& x) {
  Amounts out;
  for (const auto& [a, d, l] : x.got) out.push_back(d);
  return out;
}
} // namespace

// Losses inside one program: a merge, a tie of the split keys, the rounding, the split hierarchy,
// a link of service 0, subdividing a link before or after, a support's Pinned loss.
TEST(ProgramLossy, Merge) {
  EXPECT_EQ(lossy_rows(graph(2, {{0, 1, 9, 1}}), {0}, {1}, {{0, {0}, 2}, {0, {0}, 1}}, 12).got,
            (Triples {{8, 6, 2}, {4, 3, 1}}));
  EXPECT_EQ(lossy_rows(graph(2, {{0, 1, 1, 1}}), {0}, {1}, {{0, {0}, 1}, {0, {0}, 1}}, 2).got,
            (Triples {{1, 1, 0}, {1, 0, 1}}));
  EXPECT_EQ(lossy_rows(graph(2, {{0, 1, 3, 1}}), {0}, {1}, {{0, {0}, 5}, {0, {0}, 1}, {0, {0}, 1}}, 7).got,
            (Triples {{5, 3, 2}, {1, 0, 1}, {1, 0, 1}}));
  const auto hier = lossy_rows(graph(4, kHier), {0, 1}, {3}, {{0, {0, 2}, 1}, {0, {0, 2}, 1}, {1, {1, 2}, 1}}, 4);
  EXPECT_EQ(delivered_of(hier), (Amounts {1, 0, 0}));
  EXPECT_EQ(hier.got, (Triples {{1, 1, 0}, {1, 0, 1}, {2, 0, 2}}));
  const auto zero = lossy_rows(graph(2, {{0, 1, 0, 1}}), {0}, {1}, {{0, {0}, 1}}, 2);
  EXPECT_EQ(zero.got, (Triples {{2, 0, 2}}));
  EXPECT_EQ(zero.r.call_edge_drops, (Drops {{0, 2}}));
  const auto before = graph(6, {{0, 1, 1, 1}, {1, 2, 1, 1}, {1, 3, 1, 1}, {1, 4, 1, 1},
                                {2, 5, 1, 1}, {3, 5, 1, 1}, {4, 5, 1, 1}});
  EXPECT_EQ(delivered_of(lossy_rows(before, {0}, {5}, {{0, {0, 1, 4}, 1}, {0, {0, 2, 5}, 1}, {0, {0, 3, 6}, 1}}, 3)),
            (Amounts {1, 0, 0}));
  const auto after = graph(7, {{0, 1, 1, 1}, {1, 2, 1, 1}, {1, 3, 1, 1}, {1, 4, 1, 1},
                               {2, 6, 1, 1}, {3, 5, 1, 1}, {4, 5, 1, 1}, {6, 5, 1, 1}});
  EXPECT_EQ(delivered_of(lossy_rows(after, {0}, {5}, {{0, {0, 1, 4, 7}, 1}, {0, {0, 2, 5}, 1}, {0, {0, 3, 6}, 1}}, 3)),
            (Amounts {1, 0, 0}));
  const auto pin = lossy_rows(graph(4, kPinloss), {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}}, 12);
  EXPECT_EQ(pin.got, (Triples {{6, 2, 4}, {6, 6, 0}}));
  EXPECT_EQ(pin.r.call_edge_drops, (Drops {{2, 4}}));
  EXPECT_EQ(pin.r.call_delivered, 8);
}

// No link can drop (U0 = 2 on both barriers' candidates): no rounds.
TEST(ProgramCycles, NoBarrierNoRounds) {
  const auto run = lossy_rows(cyc(2, 2), {4}, {5}, {{4, kCycT0, 1}, {4, kCycT1, 1}}, 2);
  EXPECT_EQ(delivered_of(run), (Amounts {1, 1}));
  EXPECT_TRUE(run.r.call_edge_drops.empty());
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().groups, 0);
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 0);
#endif
}

// X and Y wait on each other; the rounds from zero settle at round 3.
TEST(ProgramCycles, Converges) {
  const auto run = lossy_rows(cyc(1, 1), {4}, {5}, {{4, kCycT0, 1}, {4, kCycT1, 1}}, 2);
  EXPECT_EQ(run.got, (Triples {{1, 1, 0}, {1, 0, 1}}));
  EXPECT_EQ(run.r.call_edge_drops, (Drops {{3, 1}}));
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().groups, 1);
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 3);
#endif
}

// Amounts (1, 2): the rounds alternate between two states; nothing is published.
TEST(ProgramCycles, PeriodTwo) {
  const auto g = cyc(1, 2);
  auto p = pinned(g, tunnels(*g, {4}, {5}, {{4, kCycT0, 1}, {4, kCycT1, 2}}), LOSSY);
  const auto ledger = ledger_state(*p.ledger);
  const auto before = epoch(*p.p);
  EXPECT_EQ(error_of([&] { (void)p->replace(3); }), kPeriodTwo);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  // The error comes at round 3, the first round whose state equals the one two rounds before.
  EXPECT_EQ(detail::thread_lossy_counters().groups, 1);
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 3);
#endif
  expect_same(ledger_state(*p.ledger), ledger, "period two");
  expect_unchanged(epoch(*p.p), before, "period two");
  EXPECT_EQ(p->program(), nullptr);
  EXPECT_EQ(error_of([&] { (void)p->place(3); }), kPeriodTwo);
  EXPECT_EQ(error_of([&] { (void)p->evaluate(3); }), kPeriodTwo);
}

// Three barriers that never settle: the round limit 64 * 3.
TEST(ProgramCycles, RoundLimit) {
  const auto lim = graph(8, {{0, 1, 1, 1}, {1, 4, 1000, 1}, {2, 3, 2, 1}, {3, 0, 1000, 1}, {3, 7, 1000, 1},
                             {4, 5, 1, 1}, {5, 2, 1000, 1}, {5, 7, 1000, 1}, {6, 0, 1000, 1}, {6, 2, 1000, 1}});
  const auto prog = tunnels(*lim, {6}, {7}, {{6, {8, 0, 1, 5, 6, 2, 4}, 2}, {6, {9, 2, 3, 0, 1, 5, 7}, 3}});
  auto p = pinned(lim, prog, LOSSY);
  EXPECT_EQ(error_of([&] { (void)p->replace(3); }), kRoundLimit);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 192);
#endif
}

// Y split in two in series changes nothing (serial normalization).
TEST(ProgramCycles, SubdivisionInvariant) {
  const auto run = lossy_rows(cycs(1, 1, 1), {4}, {5}, {{4, kCycsT0, 1}, {4, kCycsT1, 1}}, 2);
  EXPECT_EQ(run.got, (Triples {{1, 1, 0}, {1, 0, 1}}));
  EXPECT_EQ(run.r.call_edge_drops, (Drops {{3, 1}}));
  const auto g = cycs(1, 2, 2);
  auto p = pinned(g, tunnels(*g, {4}, {5}, {{4, kCycsT0, 1}, {4, kCycsT1, 2}}), LOSSY);
  EXPECT_EQ(error_of([&] { (void)p->replace(3); }), kPeriodTwo);
}

// A scalar call's black holes add no barrier: STALE-black-hole-in-cycle's graph (input rows 0 to
// 18, the example's EdgeIds there) with one head S = 17 feeding its four tunnels (wide rows 19 to
// 22) and wide exits from each tunnel's end to one destination 18 (rows 23 to 25), so every tunnel
// ends at the destination; one Pinned LOSSY program of four tunnels, weights 13, 11, 10, 5, and
// replace(39) (the heads' amounts 13, 11, 10, 5), black holes on rows 2, 6 and 17. Tunnels 0 and 3
// reach a black hole (BLACKHOLED); their arrivals there are dropped without a barrier, so the cyclic
// group keeps its four physical barriers: the round-limit error after 64 * 4 = 256 rounds (the old
// form, a black-holed link at capacity 0 inside the group, made a fifth barrier: 320). Without the
// black holes the same call settles. LOSSLESS keeps its bound: a black-holed link serves 0, so the bound is
// 0 (tunnel 0, on a black hole, takes the first unit).
TEST(ProgramCycles, BlackHolesAddNoBarrier) {
  const std::vector<InEdge> rows = {
      {0, 1, 12, 1},   {1, 6, 1000, 1},  {1, 16, 1000, 1}, {2, 3, 9, 1},     {3, 0, 1000, 1},  {3, 10, 1000, 1},
      {4, 5, 7, 1},    {5, 6, 1000, 1},  {6, 7, 13, 1},    {7, 2, 1000, 1},  {7, 10, 1000, 1}, {8, 9, 1, 1},
      {10, 11, 10, 1}, {11, 2, 1000, 1}, {12, 4, 1000, 1}, {13, 6, 1000, 1}, {14, 0, 1000, 1}, {15, 0, 1000, 1},
      {16, 2, 1000, 1}, {17, 12, 1000, 1}, {17, 13, 1000, 1}, {17, 14, 1000, 1}, {17, 15, 1000, 1},
      {11, 18, 1000, 1}, {1, 18, 1000, 1},  {7, 18, 1000, 1}};
  const auto g = shared(make_graph(19, rows));
  const auto id = id_of_input(*g);
  auto ids = [&](std::vector<std::int32_t> in) {
    std::vector<EdgeId> out;
    for (auto i : in) out.push_back(id[static_cast<std::size_t>(i)]);
    return out;
  };
  const std::vector<Spec> specs = {{17, ids({19, 14, 6, 7, 8, 9, 3, 5, 12, 23}), 13},
                                   {17, ids({20, 15, 8, 10, 12, 13, 3, 4, 0, 24}), 11},
                                   {17, ids({21, 16, 0, 1, 8, 25}), 10},
                                   {17, ids({22, 17, 0, 2, 18, 3, 5, 12, 23}), 5}};
  const auto prog = tunnels(*g, {17}, {18}, specs);
  const auto black = no_edges(*g, ids({2, 6, 17}));
  auto stale = pinned(g, prog, LOSSY, nullptr, black);
  const auto ledger = ledger_state(*stale.ledger);
  EXPECT_EQ(error_of([&] { (void)stale->replace(39); }), kRoundLimit);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().groups, 1);
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 256);
#endif
  expect_same(ledger_state(*stale.ledger), ledger, "the round-limit error");
  EXPECT_EQ(stale->program(), nullptr);
  // The statuses of the same program forwarded where it settles (a smaller offer).
  const auto settled = stale->evaluate(4);
  ASSERT_EQ(settled.by_tunnel.size(), 4U);
  EXPECT_EQ(settled.by_tunnel[0].status, BLACKHOLED);
  EXPECT_EQ(settled.by_tunnel[1].status, UP);
  EXPECT_EQ(settled.by_tunnel[2].status, UP);
  EXPECT_EQ(settled.by_tunnel[3].status, BLACKHOLED);
  // Control: without black holes the call settles.
  EXPECT_NO_THROW((void)pinned(g, prog, LOSSY)->replace(39));
  // LOSSLESS: a black-holed link serves 0.
  EXPECT_EQ(pinned(g, prog, LOSSLESS, nullptr, black)->evaluate(INF).call_delivered, 0);
}

// A tunnel whose share is 0 adds no wait, in a scalar call too: one Pinned LOSSY program, sources
// sA, sB; sA's tunnels A (X, L0, L1; weight 100) and Z (L0, X; weight w), sB's tunnel B (L1, L0).
// replace(4): RR gives each head 2, H(2, [100, 1]) gives Z 0. The period-2 error at round 3 with Z
// at weight 0 (not forwarded) and at weight 1 (forwarded with share 0) alike.
TEST(ProgramCycles, ZeroShareTunnelAddsNoWait) {
  // sA 0, sB 1, n1..n6 2..7, dA 8, dB 9, dZ 10; rows X, L0, L1 first.
  const auto g = shared(make_graph(11, {{2, 3, 1, 1},   {4, 5, 1, 1},   {6, 7, 2, 1},   {0, 2, 100, 1}, {3, 4, 100, 1},
                                        {5, 6, 100, 1}, {7, 8, 100, 1}, {1, 6, 100, 1}, {7, 4, 100, 1}, {5, 9, 100, 1},
                                        {0, 4, 100, 1}, {5, 2, 100, 1}, {3, 10, 100, 1}}));
  const auto id = id_of_input(*g);
  auto ids = [&](std::vector<int> in) {
    std::vector<EdgeId> out;
    for (int i : in) out.push_back(id[static_cast<std::size_t>(i)]);
    return out;
  };
  for (Amount z : {0, 1}) {
    const auto prog = tunnels(*g, {0, 1}, {8, 9, 10},
                              {{0, ids({3, 0, 4, 1, 5, 2, 6}), 100}, {0, ids({10, 1, 11, 0, 12}), z}, {1, ids({7, 2, 8, 1, 9}), 1}});
    auto p = pinned(g, prog, LOSSY);
    EXPECT_EQ(error_of([&] { (void)p->replace(4); }), kPeriodTwo) << "Z weight " << z;
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(detail::thread_lossy_counters().rounds, 3) << "Z weight " << z;
#endif
  }
}

// The LOSSLESS search probes totals up to the first links' capacity: here 3 * 2^60 - 1 over two
// tunnels that share edge 2 (0 -> 1 -> 3 and 0 -> 2 -> 1 -> 3, round robin), so the first probe's
// link sums take the saturating no-drop pass (a total of 2^61 or more on a link shared twice) and
// the later probes the plain one; the bound is the shared edge's 10 either way. LOSSY forwarding of
// a (2^60 and 2^60 - 1 per tunnel): edge 3 (2^40) drops 2^60 - 1 - 2^40 of tunnel 1, edge 2 gives
// its 10 to tunnel 0 (highest averages over 2^60 and 2^40) and drops the rest.
TEST(LosslessBound, WideFirstLinksProbeBothSums) {
  const Amount a = (Amount {1} << 61) - 1;
  const Amount b = Amount {1} << 60;
  const auto g = shared(make_graph(4, {{0, 1, a, 1}, {0, 2, b, 1}, {1, 3, 10, 1}, {2, 1, Amount {1} << 40, 1}}));
  const auto prog = tunnels(*g, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3, 2}, 1}});
  const auto lossless = pinned(g, prog, LOSSLESS)->evaluate(INF);
  EXPECT_EQ(lossless.call_delivered, 10);
  EXPECT_EQ(rows(lossless), (Rows {{0, 0, UP, 5, 5, 0}, {1, 0, UP, 5, 5, 0}}));
  const auto lossy = pinned(g, prog, LOSSY)->replace(a);
  EXPECT_EQ(lossy.call_delivered, 10);
  EXPECT_EQ(lossy.call_edge_drops, (std::vector<std::pair<EdgeId, Amount>> {{2, b - 10 + (Amount {1} << 40)}, {3, b - 1 - (Amount {1} << 40)}}));
  EXPECT_EQ(rows(lossy), (Rows {{0, 0, UP, b, 10, b - 10}, {1, 0, UP, b - 1, 0, b - 1}}));
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
// Through the hook: an iterated result corrupted after its rounds fails the independent check
// with std::runtime_error, and nothing is published.
TEST(ProgramCycles, VerificationCatchesACorruptedResult) {
  const auto g = cyc(1, 1);
  auto p = pinned(g, tunnels(*g, {4}, {5}, {{4, kCycT0, 1}, {4, kCycT1, 1}}), LOSSY);
  const auto ledger = ledger_state(*p.ledger);
  const auto before = epoch(*p.p);
  detail::thread_lossy_counters().corrupt = true;
  std::string what;
  try {
    (void)p->replace(2);
  } catch (const std::runtime_error& e) {
    what = e.what();
  }
  detail::thread_lossy_counters().corrupt = false;
  EXPECT_EQ(what, "lossy forwarding fixed-point verification failed");
  expect_same(ledger_state(*p.ledger), ledger, "corrupted");
  expect_unchanged(epoch(*p.p), before, "corrupted");
  // Control: the same call without the corruption succeeds.
  EXPECT_EQ(p->replace(2).call_delivered, 1);
}
#endif

// ---------------------------------------------------------------- replay under a mask

namespace {
// G1U's optimized program: [(0, [e0, e2], 2), (0, [e1, e3], 10)].
std::shared_ptr<const ForwardingProgram> g1u_program(const std::shared_ptr<const Graph>& g) {
  auto solver = make(g, {0}, {3}, global());
  const auto prog = solver->replace(12).program;
  EXPECT_EQ(listing(*prog), (Listing {{0, {0, 2}, 2}, {0, {1, 3}, 10}}));
  return prog;
}
} // namespace

// The masked-replay call with a CentralizedTE program: tunnel 0 FAILED, the head's 12 goes to tunnel 1.
TEST(Replay, GlobalProgram) {
  const auto g = graph(4, kG1U);
  const auto prog = g1u_program(g);
  const auto r = pinned(g, prog, LOSSY, no_edges(*g, {2}))->replace(12);
  EXPECT_EQ(r.call_delivered, 10);
  EXPECT_EQ(rows(r), (Rows {{0, 0, FAILED, 0, 0, 0}, {1, 0, UP, 12, 10, 2}}));
  EXPECT_EQ(r.call_edge_drops, (Drops {{1, 2}}));
  EXPECT_EQ(pinned(g, prog, LOSSLESS, no_edges(*g, {2}))->place(INF).call_delivered, 10);
}

// The surviving sources' captured weights: H(16, (3, 1)) = (12, 4).
TEST(Replay, SourceSurvivor) {
  const auto g = graph(4, kSurv);
  auto solver = make(g, {0, 1, 2}, {3}, global());
  const auto prog = solver->place(9).program;
  EXPECT_EQ(*prog->source_weights(), (std::vector<std::pair<NodeId, Amount>> {{0, 3}, {1, 1}, {2, 5}}));
  const auto r = pinned(g, prog, LOSSY, no_nodes(*g, {2}))->replace(16);
  EXPECT_EQ(r.call_delivered, 4);
  EXPECT_EQ(lo(r.by_source[0].offered), 12);
  EXPECT_EQ(lo(r.by_source[1].offered), 4);
  EXPECT_EQ(lo(r.by_source[2].offered), 0);
  EXPECT_EQ(r.call_edge_drops, (Drops {{0, 9}, {1, 3}}));
  EXPECT_EQ(std::get<2>(rows(r)[2]), FAILED);
  EXPECT_EQ(pinned(g, prog, LOSSLESS, no_nodes(*g, {2}))->place(INF).call_delivered, 4);
}

// The topology basis keeps head 1 (it still reaches 2 over e2) and it loses its share.
TEST(Replay, TopologyBasisKeepsTheHead) {
  const auto g = graph(3, kRB);
  auto solver = make(g, {0, 1}, {2}, global());
  const auto prog = solver->place(10).program;
  EXPECT_EQ(listing(*prog), (Listing {{0, {0}, 5}, {1, {1}, 5}}));
  const auto r = pinned(g, prog, LOSSY, no_edges(*g, {1}))->replace(10);
  EXPECT_EQ(r.call_delivered, 5);
  EXPECT_EQ(r.source_drops, c(5));
  EXPECT_EQ(pinned(g, prog, LOSSLESS, no_edges(*g, {1}))->place(INF).call_delivered, 1);
}

// The supplied basis drops head 1: head 0 gets all 10.
TEST(Replay, SuppliedBasisDropsTheHead) {
  const auto g = graph(3, kRB);
  const auto prog = tunnels(*g, {0, 1}, {2}, {{0, {0}, 1}, {1, {1}, 1}});
  const auto r = pinned(g, prog, LOSSY, no_edges(*g, {1}))->replace(10);
  EXPECT_EQ(r.call_delivered, 5);
  EXPECT_EQ(r.call_edge_drops, (Drops {{0, 5}}));
  EXPECT_EQ(r.source_drops, c(0));
  EXPECT_EQ(pinned(g, prog, LOSSLESS, no_edges(*g, {1}))->place(INF).call_delivered, 5);
}

// A link of free capacity 0 is not a failure.
TEST(Replay, ZeroCapacityIsNotFailed) {
  const auto g = graph(4, kG1U);
  const auto prog = g1u_program(g);
  const auto r = pinned(g, prog, LOSSY, nullptr, nullptr, std::make_shared<FlowLedger>(g, Amounts {2, 0, 2, 10}))
                     ->replace(12);
  EXPECT_EQ(r.call_delivered, 2);
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 2, 2, 0}, {1, 0, UP, 10, 0, 10}}));
  EXPECT_EQ(r.call_edge_drops, (Drops {{1, 10}}));
}

// The remaining head's captured weight is 0: the whole offer is lost at the sources.
TEST(Replay, AllRemainingWeightsZero) {
  const auto g = graph(3, {{0, 2, 5, 1}, {1, 2, 5, 1}});
  auto solver = make(g, {0, 1}, {2}, global());
  const auto prog = solver->place(1).program;
  const auto r = pinned(g, prog, LOSSY, no_nodes(*g, {0}))->replace(4);
  EXPECT_EQ(r.termination, Termination::NoActiveWeight);
  EXPECT_EQ(r.source_drops, c(4));
  EXPECT_EQ(r.call_delivered, 0);
  EXPECT_EQ(rows(r), (Rows {{0, 0, FAILED, 0, 0, 0}}));
  const auto ll = pinned(g, prog, LOSSLESS, no_nodes(*g, {0}))->replace(4);
  EXPECT_EQ(ll.termination, Termination::NoActiveWeight);
  EXPECT_EQ(ll.call_delivered, 0);
}

// A stale remote tail repaired: the branch 0 -> 2 is pruned.
TEST(Replay, SplitRoutePruned) {
  const auto g = graph(4, kDiam);
  auto producer = make(g, {0}, {3}, cost_only(EB, LOSSY));
  const auto prog = producer->place(1).program;
  const auto r = pinned(g, prog, LOSSY, no_edges(*g, {3}))->replace(100, std::nullopt, kEdges);
  EXPECT_EQ(r.call_delivered, 100);
  EXPECT_EQ(*r.edge_flows, (Amounts {100, 0, 100, 0}));
}

// A tunnel failed by its node.
TEST(Replay, TunnelFailedByNode) {
  const auto g = graph(4, kTun);
  const auto prog = tunnels(*g, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}});
  const auto r = pinned(g, prog, LOSSY, no_nodes(*g, {2}))->replace(10);
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 10, 10, 0}, {1, 0, FAILED, 0, 0, 0}}));
  EXPECT_EQ(r.call_delivered, 10);
  EXPECT_EQ(pinned(g, prog, LOSSLESS, no_nodes(*g, {2}))->place(INF).call_delivered, 10);
}

// Weights frozen at preparation: a member pruned there stays pruned in every replay.
TEST(Replay, FrozenWeights) {
  const auto z = par2(0, 10);
  auto producer = make(z, {0}, {1}, cost_only(W, LOSSY));
  const auto wp = producer->place(1).program;
  const auto r = pinned(z, wp, LOSSY)->replace(5);
  EXPECT_EQ(r.call_delivered, 5);
  EXPECT_TRUE(r.call_edge_drops.empty());
  // A snapshot over the baseline (0, 10) prunes e0; the replay on a fresh ledger at full capacity
  // still sends everything to e1: 10 delivered, 5 dropped there.
  const auto g = par2(10, 10);
  auto snap = make(g, {0}, {1}, weighted(LOSSY, WeightBasis::ResidualSnapshot), nullptr,
                   std::make_shared<FlowLedger>(g, Amounts {0, 10}));
  const auto sp = snap->place(1).program;
  const auto rr = pinned(g, sp, LOSSY)->replace(15, std::nullopt, kEdges);
  EXPECT_EQ(rr.call_delivered, 10);
  EXPECT_EQ(*rr.edge_flows, (Amounts {0, 10}));
  EXPECT_EQ(rr.call_edge_drops, (Drops {{1, 5}}));
}

// ---------------------------------------------------------------- stale forwarding

// Stale forwarding: a remote tail, a WCMP member (also of capacity 0), the relief and its original
// program, a tunnel node, a pooled original and remainder, all weights zero, a black-holed source
// (also of a tunnel), a masked and black-holed element.
TEST(StaleForwarding, RemoteTail) {
  {
    const auto g = graph(4, kDiam);
    auto p = make(g, {0}, {3}, cost_only(EB, LOSSY));
    (void)p->replace(100);
    const auto installed = p->program();
    const auto bh = g->eligibility({}, std::vector<EdgeId> {3});
    const auto r = p->forward(100, &bh);
    EXPECT_EQ(r.call_delivered, 50);
    EXPECT_EQ(r.call_edge_drops, (Drops {{3, 50}}));
    EXPECT_EQ(p->program(), installed);
    EXPECT_FALSE(r.routes_changed || r.weights_changed);
    EXPECT_EQ(p->replace(100).call_delivered, 100);
  }
  {
    const auto w = graph(2, {{0, 1, 10, 1}, {0, 1, 40, 1}});
    auto q = make(w, {0}, {1}, cost_only(W, LOSSY));
    (void)q->replace(40);
    const auto bh = w->eligibility({}, std::vector<EdgeId> {0});
    const auto r = q->forward(40, &bh);
    EXPECT_EQ(r.call_delivered, 32);
    EXPECT_EQ(r.call_edge_drops, (Drops {{0, 8}}));
  }
  {
    const auto wz = graph(2, {{0, 1, 0, 1}, {0, 1, 40, 1}});
    auto q = make(wz, {0}, {1}, cost_only(W, LOSSY));
    (void)q->replace(1);
    EXPECT_EQ(q->forward(40).call_delivered, 40);
  }
  {
    const auto rel = graph(4, kRel);
    const auto prog = tunnels(*rel, {0, 1}, {3}, {{0, {0, 2}, 1}, {1, {1, 2}, 1}});
    EXPECT_EQ(rows(pinned(rel, prog, LOSSY)->replace(20)), (Rows {{0, 0, UP, 10, 5, 5}, {1, 1, UP, 10, 5, 5}}));
    EXPECT_EQ(rows(pinned(rel, prog, LOSSY, nullptr, no_edges(*rel, {0}))->replace(20)),
              (Rows {{0, 0, BLACKHOLED, 10, 0, 10}, {1, 1, UP, 10, 10, 0}}));
  }
  {
    const auto t = graph(4, kTun);
    const auto tp = tunnels(*t, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}});
    const auto r = pinned(t, tp, LOSSY, nullptr, no_nodes(*t, {2}))->replace(10);
    EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 5, 5, 0}, {1, 0, BLACKHOLED, 5, 0, 5}}));
    EXPECT_EQ(r.call_edge_drops, (Drops {{1, 5}}));
    EXPECT_EQ(pinned(t, tp, LOSSLESS, nullptr, no_nodes(*t, {2}))->replace(10).call_delivered, 1);
  }
  {
    const auto pool = graph(4, kPool2);
    auto p = make(pool, {0, 1}, {3}, cost_only(EB, LOSSY));
    (void)p->replace(0);
    const auto r = p->forward(2, nullptr, kEdges);
    EXPECT_EQ(r.call_delivered, 2);
    EXPECT_EQ(*r.edge_flows, (Amounts {1, 1, 1, 1}));
    const auto bh = pool->eligibility({}, std::vector<EdgeId> {0, 3});
    const auto s8 = p->forward(2, &bh);
    EXPECT_EQ(s8.call_delivered, 1);
    EXPECT_EQ(s8.call_edge_drops, (Drops {{0, 1}}));
  }
  {
    const auto zero = graph(2, {{0, 1, 0, 1}});
    auto q = make(zero, {0}, {1}, cost_only(W, LOSSY));
    (void)q->replace(0);
    const auto r = q->forward(2);
    EXPECT_EQ(r.source_drops, c(2));
    EXPECT_EQ(r.termination, Termination::NoActiveWeight);
  }
  {
    const auto one = graph(2, kOne);
    auto q = make(one, {0}, {1}, cost_only(EB, LOSSY));
    (void)q->replace(0);
    const auto bh = one->eligibility(std::vector<NodeId> {0}, {});
    const auto r = q->forward(10, &bh);
    EXPECT_EQ(r.source_drops, c(10));
    EXPECT_TRUE(r.call_edge_drops.empty());
    EXPECT_EQ(r.call_delivered, 0);
    const auto s11 = pinned(one, tunnels(*one, {0}, {1}, {{0, {0}, 1}}), LOSSY, nullptr, no_nodes(*one, {0}))->replace(10);
    EXPECT_EQ(rows(s11), (Rows {{0, 0, BLACKHOLED, 0, 0, 0}}));
    EXPECT_EQ(s11.source_drops, c(10));
  }
  {
    // A masked and black-holed element is masked.
    const auto m = graph(5, kMB);
    const auto mp = tunnels(*m, {0}, {4}, {{0, {0, 3}, 1}, {0, {1, 4}, 1}, {0, {2, 5}, 1}});
    const auto r = pinned(m, mp, LOSSY, no_edges(*m, {3}), no_edges(*m, {3, 4}))->replace(10);
    EXPECT_EQ(rows(r), (Rows {{0, 0, FAILED, 0, 0, 0}, {1, 0, BLACKHOLED, 5, 0, 5}, {2, 0, UP, 5, 5, 0}}));
    EXPECT_EQ(r.call_edge_drops, (Drops {{4, 5}}));
  }
}

// The masked-replay call's CentralizedTE program with e2 black-holed instead of masked.
TEST(StaleForwarding, GlobalProgram) {
  const auto g = graph(4, kG1U);
  const auto prog = g1u_program(g);
  const auto r = pinned(g, prog, LOSSY, nullptr, no_edges(*g, {2}))->replace(12);
  EXPECT_EQ(rows(r), (Rows {{0, 0, BLACKHOLED, 2, 0, 2}, {1, 0, UP, 10, 10, 0}}));
  EXPECT_EQ(r.call_edge_drops, (Drops {{2, 2}}));
  // H(4, (2, 10)) = (0, 4) fits; H(5, (2, 10)) = (1, 4) sends a unit into the black hole.
  EXPECT_EQ(pinned(g, prog, LOSSLESS, nullptr, no_edges(*g, {2}))->replace(12).call_delivered, 4);
}

// An all-allowed Eligibility is no black hole.
TEST(StaleForwarding, AllAllowedIsNone) {
  const auto g = graph(4, kG1U);
  const auto prog = g1u_program(g);
  const auto allowed = std::make_shared<const Eligibility>(g->eligibility());
  const auto with = pinned(g, prog, LOSSY, nullptr, allowed)->replace(12);
  const auto without = pinned(g, prog, LOSSY)->replace(12);
  EXPECT_EQ(rows(with), rows(without));
  EXPECT_EQ(with.call_delivered, 12);
  EXPECT_EQ(with.call_edge_drops, without.call_edge_drops);
  // A forward with an all-allowed value equals a forward without.
  auto p = make(g, {0}, {3}, cost_only(EB, LOSSY));
  (void)p->replace(12);
  const auto all = g->eligibility();
  EXPECT_EQ(p->forward(12, &all).call_delivered, p->forward(12).call_delivered);
}

// Place, evaluate, replace, evaluate_demand and evaluate_edge_failures all see a Pinned policy's black holes
// (TUN, node 2 black-holed).
TEST(StaleForwarding, EveryOperation) {
  const auto t = graph(4, kTun);
  const auto tp = tunnels(*t, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}});
  const auto bh = no_nodes(*t, {2});
  EXPECT_EQ(pinned(t, tp, LOSSY, nullptr, bh)->place(10).call_delivered, 5);
  EXPECT_EQ(pinned(t, tp, LOSSY, nullptr, bh)->evaluate(10).call_delivered, 5);
  EXPECT_EQ(pinned(t, tp, LOSSY, nullptr, bh)->replace(10).call_delivered, 5);
  EXPECT_EQ(evaluate_demand(*t, v({0}), v({3}), query(pinned_config(tp, LOSSY, bh), offer(10))).call_delivered, 5);
  EXPECT_EQ(evaluate_demand(*t, v({0}), v({3}), query(pinned_config(tp, LOSSLESS, bh))).call_delivered, 1);
  EXPECT_EQ(evaluate_edge_failures(*t, v({0}), v({3}), query(pinned_config(tp, LOSSY, bh), offer(10)), 1).baseline.call_delivered, 5);
  // Controls without the black holes.
  EXPECT_EQ(pinned(t, tp, LOSSY)->replace(10).call_delivered, 10);
  EXPECT_EQ(evaluate_demand(*t, v({0}), v({3}), query(pinned_config(tp))).call_delivered, 20);
}

// ---------------------------------------------------------------- operations

// Evaluate publishes nothing: ledger, records and epoch counters unchanged.
TEST(Operations, EvaluatePublishesNothing) {
  const auto g = par2(10, 10);
  auto p = make(g, {0}, {1}, cost_only());
  (void)p->replace(6);
  const auto ledger = ledger_state(*p.ledger);
  const auto before = epoch(*p.p);
  const auto r = p->evaluate(12);
  EXPECT_EQ(r.call_delivered, 12);
  EXPECT_EQ(r.delivery_delta, 6);
  EXPECT_NE(r.program, nullptr);
  expect_same(ledger_state(*p.ledger), ledger, "evaluate");
  expect_unchanged(epoch(*p.p), before, "evaluate");
  EXPECT_EQ(p->place(0).requested, c(6));
}

// Replace retires the old ids and starts a new epoch.
TEST(Operations, ReplaceRetiresIds) {
  const auto g = par2(10, 10);
  auto p = make(g, {0}, {1}, cost_only());
  (void)p->place(4);
  const auto old = p->records().front().flow_id;
  const auto r = p->replace(6);
  EXPECT_EQ(r.requested, c(6));
  EXPECT_EQ(r.admitted, c(6));
  EXPECT_NE(error_of([&] { (void)p.ledger->record(old); }), "");
  EXPECT_NE(p->records().front().flow_id, old);
}

// Forward keeps the frozen weights; replace takes a new snapshot; an identical replace
// sets no flag.
TEST(Operations, ForwardKeepsTheProgram) {
  const auto g = par2(10, 10);
  auto ledger = std::make_shared<FlowLedger>(g, Amounts {2, 10});
  auto a = make(g, {0}, {1}, weighted(LOSSLESS, WeightBasis::ResidualSnapshot), nullptr, ledger);
  EXPECT_EQ(a->replace(6).call_delivered, 6);
  auto b = pinned(g, tunnels(*g, {0}, {1}, {{0, {1}, 1}}), LOSSLESS, nullptr, nullptr, ledger);
  EXPECT_EQ(b->place(5).call_delivered, 5);
  const auto installed = a->program();
  const auto f = a->forward(12);
  EXPECT_EQ(f.call_delivered, 6);
  EXPECT_FALSE(f.routes_changed || f.weights_changed);
  EXPECT_EQ(a->program(), installed);
  const auto r = a->replace(12);
  EXPECT_EQ(r.call_delivered, 7);
  EXPECT_FALSE(r.routes_changed);
  EXPECT_TRUE(r.weights_changed);
}

TEST(Operations, IdenticalReplace) {
  const auto g = par2(10, 10);
  auto ledger = std::make_shared<FlowLedger>(g, Amounts {2, 10});
  auto a = make(g, {0}, {1}, weighted(LOSSLESS, WeightBasis::ResidualSnapshot), nullptr, ledger);
  const auto r1 = a->replace(12);
  const auto ids1 = ids(*a.p);
  const auto r2 = a->replace(12);
  EXPECT_EQ(r2.call_delivered, r1.call_delivered);
  EXPECT_FALSE(r2.routes_changed || r2.weights_changed);
  EXPECT_NE(ids(*a.p), ids1);
  EXPECT_EQ(delivered(*a.p), (Amounts {r1.call_delivered}));
}

// The forward texts.
TEST(Operations, ForwardTexts) {
  const auto g = par2(10, 10);
  auto p = make(g, {0}, {1}, cost_only(EB, LOSSY));
  EXPECT_EQ(error_of([&] { (void)p->forward(1); }), kNoProgram);
  (void)p->place(1);
  for (const auto& call : std::vector<std::function<void()>> {[&] { (void)p->forward(INF); },
                                                              [&] { (void)p->replace(INF); },
                                                              [&] { (void)p->evaluate(INF); }}) {
    EXPECT_EQ(error_of(call), kFinite);
  }
  p->remove();
  EXPECT_EQ(error_of([&] { (void)p->forward(1); }), kNoProgram);
  auto ca = make(g, {0}, {1}, capacity_aware());
  EXPECT_EQ(error_of([&] { (void)ca->forward(1); }), kNoProgram);
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}});
  auto pin = pinned(g, prog, LOSSY);
  EXPECT_EQ(error_of([&] { (void)pin->place(INF); }), kFinite);
  (void)pin->place(1);
  EXPECT_EQ(error_of([&] { (void)pin->forward(1); }), "Pinned forwards with replace");
}

// The order of a call's checks.
TEST(Operations, ForwardOrder) {
  const auto g = par2(10, 10);
  const auto other = par2(1, 1);
  const auto foreign = other->eligibility({}, std::vector<EdgeId> {0});
  auto pin = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 1}}));
  EXPECT_EQ(error_of([&] { (void)pin->forward(INF, &foreign); }), "Pinned forwards with replace");
  auto p = make(g, {0}, {1}, cost_only(EB, LOSSY));
  EXPECT_EQ(error_of([&] { (void)p->forward(-1); }), "offer must be an integer in [0, 2^62) or permitted +inf");
  (void)p->place(1);
  EXPECT_EQ(error_of([&] { (void)p->forward(-1, &foreign); }), "offer must be an integer in [0, 2^62) or permitted +inf");
  EXPECT_EQ(error_of([&] { (void)p->forward(INF, &foreign); }), kFinite);
  EXPECT_EQ(error_of([&] { (void)p->forward(1, &foreign); }), "black_holes belongs to another graph");
  EXPECT_EQ(error_of([&] { (void)p->replace(-1); }), "offer must be an integer in [0, 2^62) or permitted +inf");
  EXPECT_EQ(error_of([&] { (void)p->replace(INF); }), kFinite);
}

// Lossless replace(+inf): the bound, EXHAUSTED.
TEST(Operations, LosslessInfinite) {
  const auto r = make(par2(10, 10), {0}, {1}, cost_only())->replace(INF);
  EXPECT_EQ(r.call_delivered, 20);
  EXPECT_FALSE(r.call_undelivered.has_value());
  EXPECT_EQ(r.termination, Termination::Exhausted);
}

// A zero replace installs a program of two split routes and two records at zero.
TEST(Operations, ZeroReplacePrograms) {
  const auto g = graph(4, kPool2);
  auto p = make(g, {0, 1}, {3}, cost_only());
  const auto r = p->replace(0);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_TRUE(r.routes_changed);
  EXPECT_TRUE(r.weights_changed);
  ASSERT_NE(r.program, nullptr);
  EXPECT_EQ(delivered(*p.p), (Amounts {0, 0}));
}

// A smaller replace: delivery_delta negative.
TEST(Operations, DeliveryDeltaNegative) {
  auto p = make(par2(10, 10), {0}, {1}, cost_only());
  (void)p->replace(10);
  const auto r = p->replace(4);
  EXPECT_EQ(r.call_delivered, 4);
  EXPECT_EQ(r.delivery_delta, -6);
}

// Place continues the replace's epoch (its ids kept).
TEST(Operations, PlaceContinuesTheEpoch) {
  const auto g = par2(10, 10);
  auto p = make(g, {0}, {1}, cost_only());
  (void)p->replace(6);
  const auto before = ids(*p.p);
  const auto r = p->place(3);
  EXPECT_EQ(r.call_delivered, 3);
  EXPECT_EQ(r.requested, c(9));
  EXPECT_EQ(r.admitted, c(9));
  EXPECT_EQ(ids(*p.p), before);
  auto q = make(g, {0}, {1}, cost_only(EB, LOSSY));
  (void)q->replace(6);
  const auto lr = q->place(3);
  EXPECT_EQ(lr.call_delivered, 3);
  EXPECT_EQ(lr.offered, c(9));
  EXPECT_EQ(lr.requested, c(9));
}

// A controller's increment changes the weights.
TEST(Operations, ControllerWeights) {
  const auto g = graph(2, {{0, 1, 20, 1}});
  auto p = make(g, {0}, {1}, global());
  const auto r1 = p->place(5);
  EXPECT_TRUE(r1.routes_changed && r1.weights_changed);
  const auto r2 = p->place(5);
  EXPECT_FALSE(r2.routes_changed);
  EXPECT_TRUE(r2.weights_changed);
  EXPECT_EQ(listing(*r2.program), (Listing {{0, {0}, 10}}));
  const auto r3 = p->replace(10);
  EXPECT_FALSE(r3.routes_changed || r3.weights_changed);
}

// A re-solve that adds a path changes the routes.
TEST(Operations, ControllerRoutes) {
  const auto g = graph(4, kB1);
  auto q = make(g, {0, 1}, {3}, global());
  EXPECT_EQ(listing(*q->place(2).program), (Listing {{0, {0}, 2}}));
  const auto r = q->place(10);
  EXPECT_EQ(listing(*r.program), (Listing {{0, {0}, 2}, {1, {1, 2}, 10}}));
  EXPECT_TRUE(r.routes_changed);
}

// An offer-only change sets no flag.
TEST(Operations, OfferOnly) {
  const auto one = graph(2, {{0, 1, 20, 1}});
  auto pin = pinned(one, tunnels(*one, {0}, {1}, {{0, {0}, 1}}));
  (void)pin->replace(5);
  for (const Amount x : {10, 10}) {
    const auto before = ids(*pin.p);
    const auto r = pin->replace(x);
    EXPECT_FALSE(r.routes_changed || r.weights_changed) << x;
    EXPECT_NE(ids(*pin.p), before);
  }
  auto cost = make(par2(10, 10), {0}, {1}, cost_only(EB, LOSSY));
  (void)cost->place(1);
  for (const Amount x : {5, 7}) {
    const auto r = cost->forward(x);
    EXPECT_FALSE(r.routes_changed || r.weights_changed) << x;
  }
}

// Forward over a controller's program reproduces the solve; above it, EXHAUSTED.
TEST(Operations, ForwardControllerProgram) {
  const auto g = graph(4, kG1U);
  auto p = make(g, {0}, {3}, global());
  (void)p->replace(12);
  EXPECT_EQ(p->forward(12).call_delivered, 12);
  const auto r = p->forward(20);
  EXPECT_EQ(r.call_delivered, 12);
  EXPECT_EQ(r.termination, Termination::Exhausted);
}

// remove() clears the program.
TEST(Operations, RemoveClearsTheProgram) {
  for (const auto& cfg : {global(), cost_only()}) {
    auto p = make(par2(10, 10), {0}, {1}, cfg);
    (void)p->place(3);
    EXPECT_NE(p->program(), nullptr);
    p->remove();
    EXPECT_EQ(p->program(), nullptr);
    const auto z = p->place(0);
    EXPECT_EQ(z.program, nullptr);
    EXPECT_EQ(z.call_delivered, 0);
    EXPECT_EQ(z.delivered, 0);
    EXPECT_TRUE(z.by_tunnel.empty());
  }
}

// ShortestPaths's first positive call sets both flags; the next none.
TEST(Operations, FirstProgramFlags) {
  auto q = make(graph(4, kPool2), {0, 1}, {3}, cost_only());
  const auto r = q->place(2);
  EXPECT_TRUE(r.routes_changed && r.weights_changed);
  const auto r2 = q->place(2);
  EXPECT_FALSE(r2.routes_changed || r2.weights_changed);
}

// evaluate's outputs describe the evaluated state.
TEST(Operations, EvaluateOutputs) {
  auto gl = make(graph(4, kB1), {0, 1}, {3}, global());
  (void)gl->place(2);
  const auto r = gl->evaluate(12, std::nullopt, OutputOptions {.edge_flows = true, .residual = true, .min_cut = true});
  EXPECT_EQ(*r.edge_flows, (Amounts {2, 10, 10}));
  EXPECT_EQ(*r.residual, (Amounts {0, 0, 0}));
  EXPECT_EQ(*r.min_cut, (std::vector<EdgeId> {0, 1}));
  EXPECT_EQ(gl->delivered(), 2);
  EXPECT_EQ(gl.ledger->residual(), (Amounts {0, 10, 10}));
}

// A rejected call of each kind (and an evaluator) changes nothing.
TEST(Operations, RejectedCallsChangeNothing) {
  const auto g = par2(10, 10);
  const auto other = par2(1, 1);
  const auto foreign = other->eligibility({}, std::vector<EdgeId> {0});
  auto p = make(g, {0}, {1}, cost_only(EB, LOSSY));
  (void)p->place(4);
  const auto ledger = ledger_state(*p.ledger);
  const auto before = epoch(*p.p);
  const std::vector<std::pair<std::string, std::function<void()>>> calls = {
      {"place", [&] { (void)p->place(-1); }},
      {"place 2^62", [&] { (void)p->place(kL62); }},
      {"evaluate", [&] { (void)p->evaluate(INF); }},
      {"replace", [&] { (void)p->replace(-1); }},
      {"forward", [&] { (void)p->forward(1, &foreign); }},
      {"capacity", [&] { (void)evaluate_demand(*g, v({0}), v({1}), query(cost_only(EB, LOSSY), offer(-5))); }},
  };
  for (const auto& [name, call] : calls) {
    EXPECT_NE(error_of(call), "") << name;
    expect_same(ledger_state(*p.ledger), ledger, name);
    expect_unchanged(epoch(*p.p), before, name);
  }
  // A rejected Pinned call.
  auto pin = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 1}}), LOSSY);
  (void)pin->place(3);
  const auto pl = ledger_state(*pin.ledger);
  const auto pe = epoch(*pin.p);
  for (const auto& call : std::vector<std::function<void()>> {[&] { (void)pin->place(INF); },
                                                              [&] { (void)pin->forward(1); },
                                                              [&] { (void)pin->replace(-1); }}) {
    EXPECT_NE(error_of(call), "");
    expect_same(ledger_state(*pin.ledger), pl, "pinned");
    expect_unchanged(epoch(*pin.p), pe, "pinned");
  }
}

// ---------------------------------------------------------------- records

// ShortestPaths and CentralizedTE records have no entry; a Pinned record carries its entry identifier.
TEST(ProgramRecords, Fields) {
  const auto g = graph(4, kB1);
  for (const auto& cfg : {global(), cost_only()}) {
    auto p = make(g, {0, 1}, {3}, cfg);
    (void)p->place(4);
    for (const auto& e : entries(*p.p)) EXPECT_FALSE(e.has_value());
  }
  const auto mixed = par2(4, 8);
  auto pin = pinned(mixed, program(*mixed, {0}, {1},
                                   {tunnel_route(*mixed, 0, {1}), dag_route(*mixed, 0, {0, 1}), tunnel_route(*mixed, 0, {0})}));
  (void)pin->place(3);
  EXPECT_EQ(entries(*pin.p), (std::vector<std::optional<std::int32_t>> {0, 1, 2}));
  for (const auto& r : pin->records()) {
    EXPECT_EQ(r.sources, v({0}));
    EXPECT_EQ(r.destinations, v({1}));
  }
}

// One record per entry, duplicates included; a weight-0 entry's record at zero.
TEST(ProgramRecords, PerEntry) {
  const auto g = par2(4, 8);
  auto p = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 1}, {0, {1}, 1}, {0, {1}, 1}}));
  const auto r = p->place(12);
  const auto recs = p->records();
  ASSERT_EQ(recs.size(), 3U);
  EXPECT_EQ(entries(*p.p), (std::vector<std::optional<std::int32_t>> {0, 1, 2}));
  EXPECT_EQ(delivered(*p.p), (Amounts {4, 4, 4}));
  EXPECT_EQ(ascending_edges(*recs[1].route), ascending_edges(*recs[2].route));
  EXPECT_NE(recs[1].flow_id, recs[2].flow_id);
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 4, 4, 0}, {1, 0, UP, 4, 4, 0}, {2, 0, UP, 4, 4, 0}}));
  auto z = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 0}, {0, {1}, 1}}));
  (void)z->place(INF);
  const auto zr = z->records();
  ASSERT_EQ(zr.size(), 2U);
  EXPECT_EQ(zr[0].delivered, 0);
  EXPECT_TRUE(zr[0].edge_flows.empty());
  EXPECT_EQ(ascending_edges(*zr[0].route), (std::vector<EdgeId> {0}));
}

// A FAILED tunnel's record: an empty route, min_route_cost UNREACHABLE.
TEST(ProgramRecords, FailedEntry) {
  const auto g = graph(4, kG1U);
  const auto prog = g1u_program(g);
  auto p = pinned(g, prog, LOSSY, no_edges(*g, {2}));
  (void)p->replace(12);
  const auto recs = p->records();
  ASSERT_EQ(recs.size(), 2U);
  EXPECT_TRUE(ascending_edges(*recs[0].route).empty());
  EXPECT_EQ(recs[0].min_route_cost, kUnreachable);
  EXPECT_EQ(ascending_edges(*recs[1].route), (std::vector<EdgeId> {1, 3}));
  EXPECT_EQ(recs[1].min_route_cost, 2);
}

// A ShortestPaths program in Pinned keeps the pooled attribution (each source's share of a shared
// router's traffic by highest averages).
TEST(ProgramRecords, PooledContext) {
  const auto g = graph(4, kPool2);
  auto producer = make(g, {0, 1}, {3}, cost_only(EB, LOSSY));
  const auto prog = producer->place(2).program;
  auto p = pinned(g, prog, LOSSY);
  (void)p->replace(2);
  std::vector<std::pair<std::optional<std::int32_t>, Drops>> got;
  for (const auto& r : p->records()) got.emplace_back(r.entry_id, r.edge_flows);
  EXPECT_EQ(got, (std::vector<std::pair<std::optional<std::int32_t>, Drops>> {{0, {{0, 1}, {2, 1}}},
                                                                             {1, {{1, 1}, {3, 1}}}}));
}

// ---------------------------------------------------------------- results

// The by_tunnel rows.
TEST(ProgramResults, ByTunnel) {
  const auto pl = graph(4, kPinloss);
  const auto y8 = pinned(pl, tunnels(*pl, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}}), LOSSY)->replace(12);
  EXPECT_EQ(rows(y8), (Rows {{0, 0, UP, 6, 2, 4}, {1, 0, UP, 6, 6, 0}}));
  const auto b1 = graph(4, kB1);
  EXPECT_EQ(rows(make(b1, {0, 1}, {3}, global())->place(12)), (Rows {{0, 0, UP, 2, 2, 0}, {1, 1, UP, 10, 10, 0}}));
  EXPECT_TRUE(make(graph(4, kPool2), {0, 1}, {3}, cost_only())->place(2).by_tunnel.empty());
  // A zero call reports the state's rows.
  auto p = pinned(pl, tunnels(*pl, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}}), LOSSY);
  (void)p->replace(12);
  EXPECT_EQ(rows(p->place(0)), (Rows {{0, 0, UP, 6, 2, 4}, {1, 0, UP, 6, 6, 0}}));
}

// RC02-zero-unprepared extended: no program, no rows.
TEST(ProgramResults, ZeroUnprepared) {
  for (const auto& cfg : {global(), cost_only()}) {
    auto p = make(par2(4, 6), {0}, {1}, cfg);
    const auto z = p->place(0);
    EXPECT_EQ(z.program, nullptr);
    EXPECT_TRUE(z.by_tunnel.empty());
    EXPECT_EQ(z.termination, Termination::Complete);
    (void)p->place(3);
    p->remove();
    const auto after = p->place(0);
    EXPECT_EQ(after.program, nullptr);
    EXPECT_TRUE(after.by_tunnel.empty());
  }
}

// TunnelFlow.assigned and dropped above 2^64 with a carry: 5 (2^62 - 1) =
// 23058430092136939515, high half 1, low half 4611686018427387899.
TEST(WideCounters, TunnelRowsBeyond64Bits) {
  const auto g = graph(2, {{0, 1, 0, 1}});
  auto p = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 1}}), LOSSY);
  const Amount big = kL62 - 1;
  PlacementResult r;
  for (int i = 0; i < 5; ++i) r = p->place(big);
  const Counter128 total {1, 4611686018427387899ULL};
  ASSERT_EQ(r.by_tunnel.size(), 1U);
  EXPECT_EQ(r.by_tunnel[0].assigned, total);
  EXPECT_EQ(r.by_tunnel[0].total_edge_drops, total);
  EXPECT_EQ(r.by_tunnel[0].delivered, 0);
  EXPECT_EQ(r.offered, total);
  EXPECT_EQ(r.requested, total);
  EXPECT_EQ(r.admitted, total);
  EXPECT_EQ(r.edge_drops, (std::vector<std::pair<EdgeId, Counter128>> {{0, total}}));
  EXPECT_TRUE(r.source_drops.is_zero());
  ASSERT_EQ(r.by_source.size(), 1U);
  EXPECT_EQ(r.by_source[0].offered, total);
  EXPECT_EQ(r.by_source[0].admitted, total);
  EXPECT_EQ(r.by_source[0].total_edge_drops, total);
  EXPECT_TRUE(r.by_source[0].source_drops.is_zero());
  const auto rec = p->records().front();
  EXPECT_EQ(rec.offered, total);
  EXPECT_EQ(rec.edge_drops, (std::vector<std::pair<EdgeId, Counter128>> {{0, total}}));
  // The carry: four calls stay below 2^64.
  auto q = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 1}}), LOSSY);
  PlacementResult four;
  for (int i = 0; i < 4; ++i) four = q->place(big);
  EXPECT_EQ(four.by_tunnel[0].assigned, (Counter128 {0, 4 * static_cast<std::uint64_t>(big)}));
}

// ---------------------------------------------------------------- evaluators

// Capacity of a Pinned Lossless program on PAR2(6, 4).
TEST(QueryOptionsEval, PinnedLossless) {
  const auto g = par2(6, 4);
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}, {0, {1}, 1}});
  const auto r = evaluate_demand(*g, v({0}), v({1}), query(pinned_config(prog)));
  EXPECT_EQ(r.call_delivered, 9);
  EXPECT_EQ(r.termination, Termination::Exhausted);
}

// Capacity of PINLOSS's Pinned Lossy with offered 12: 2 + 6.
TEST(QueryOptionsEval, PinnedLossy) {
  const auto pl = graph(4, kPinloss);
  const auto prog = tunnels(*pl, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}});
  EXPECT_EQ(evaluate_demand(*pl, v({0}), v({3}), query(pinned_config(prog, LOSSY), offer(12))).call_delivered, 8);
}

// The offered rules and evaluate_pairs's Pinned text.
TEST(QueryOptionsEval, Texts) {
  const auto g = par2(6, 4);
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}, {0, {1}, 1}});
  auto err = [&](const EvaluationOptions& q) { return error_of([&] { (void)evaluate_demand(*g, v({0}), v({1}), q); }); };
  auto delivered = [&](const EvaluationOptions& q) { return evaluate_demand(*g, v({0}), v({1}), q).call_delivered; };
  const std::string domain = "offer must be an integer in [0, 2^62)";
  EXPECT_EQ(err(query(pinned_config(prog, LOSSY))), "Pinned LOSSY evaluation requires a finite offer");
  EXPECT_EQ(err(query(pinned_config(prog, LOSSY), offer_inf())), domain);
  EXPECT_EQ(err(query(cost_only(), offer_inf())), domain);
  // A finite offer for every model under each admission it supports.
  EXPECT_EQ(delivered(query(cost_only(EB, LOSSY), offer(5))), 5);
  EXPECT_EQ(delivered(query(cost_only(), offer(5))), 5);
  EXPECT_EQ(delivered(query(global(), offer(5))), 5);
  EXPECT_EQ(delivered(query(pinned_config(prog), offer(5))), 5);
  auto designed = query(capacity_aware(), offer(5));
  designed.design_demand = Amount {4};
  EXPECT_EQ(err(designed), "");
  for (const Amount bad : {Amount {-1}, kL62}) {
    EXPECT_EQ(err(query(pinned_config(prog, LOSSY), offer(bad))), domain);
  }
  EXPECT_EQ(err(query(pinned_config(prog, LOSSY), offer(kL62 - 1))), "");
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{{0}, {1}}};
  EXPECT_EQ(error_of([&] { (void)evaluate_pairs(*g, pairs, query(pinned_config(prog)), {}, 1); }),
            "evaluate_pairs does not accept Pinned");
}

// A deletion fails the tunnels through it and the head splits over the others.
TEST(QueryOptionsEval, PinnedSensitivity) {
  const auto pl = graph(4, kPinloss);
  const auto prog = tunnels(*pl, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}});
  const auto s = evaluate_edge_failures(*pl, v({0}), v({3}), query(pinned_config(prog)), 1);
  EXPECT_EQ(s.baseline.call_delivered, 4);
  EXPECT_EQ(losses(s), (std::vector<std::pair<EdgeId, Amount>> {{0, -6}, {1, 2}, {2, -6}, {3, 2}}));
  const auto l = evaluate_edge_failures(*pl, v({0}), v({3}), query(pinned_config(prog, LOSSY), offer(12)), 1);
  EXPECT_EQ(l.baseline.call_delivered, 8);
  EXPECT_EQ(losses(l), (std::vector<std::pair<EdgeId, Amount>> {{0, -2}, {1, 6}, {2, -2}, {3, 6}}));
  const auto w2 = evaluate_edge_failures(*pl, v({0}), v({3}), query(pinned_config(prog, LOSSY), offer(12)), 2);
  EXPECT_EQ(losses(w2), losses(l));
  EXPECT_EQ(w2.baseline.call_delivered, 8);
}

// ---------------------------------------------------------------- common Pinned setups

// Each common Pinned setup and its program form.
TEST(PinnedMapping, SingleBundle) {
  // One EB route DAG (DIAM, the four links).
  const auto d = graph(4, kDiam);
  const auto all = dag_of(*d, {0, 1, 2, 3});
  const auto eb = program(*d, {0}, {3}, {RouteDefinition {0, all, 1}});
  EXPECT_EQ(pinned(d, eb)->place(INF).call_delivered, 200);
  EXPECT_EQ(pinned(d, eb, LOSSY)->place(100).call_delivered, 100);
  // W with capacity weights; with residual weights (the ledger's residual at creation):
  // H over (100, 50) at node 0: H(150) = (100, 50) fits, H(151) = (101, 50) does not.
  const auto wc = program(*d, {0}, {3}, {RouteDefinition {0, all, 1}}, RouteOptions {.allocation = W});
  EXPECT_EQ(pinned(d, wc)->place(INF).call_delivered, 200);
  const auto ws = program(*d, {0}, {3}, {RouteDefinition {0, all, 1}},
                          RouteOptions {.allocation = W, .weight_basis = WeightBasis::ResidualSnapshot,
                                        .residual = Amounts {100, 50, 100, 50}});
  EXPECT_EQ(pinned(d, ws)->place(INF).call_delivered, 150);
  // Two route DAGs on PAR2(6, 4): separate contexts, RR over the entries: 9; W the same 9;
  // an optimized split is CentralizedTE(support) and fills both: 10.
  const auto p = par2(6, 4);
  const auto two = program(*p, {0}, {1}, {dag_route(*p, 0, {0}), dag_route(*p, 0, {1})});
  EXPECT_EQ(pinned(p, two)->place(INF).call_delivered, 9);
  const auto two_w =
      program(*p, {0}, {1}, {dag_route(*p, 0, {0}), dag_route(*p, 0, {1})}, RouteOptions {.allocation = W});
  EXPECT_EQ(pinned(p, two_w)->place(INF).call_delivered, 9);
  EXPECT_EQ(make(p, {0}, {1}, global(std::nullopt, dag_of(*p, {0, 1})))->place(INF).call_delivered, 10);
  // Two route DAGs sharing e0 on SHARE, EB Lossy place(10): both arrive at e0 with 5:
  // H(4, (5, 5)) = (2, 2).
  const auto sh = graph(3, kShare);
  auto shared_e0 = pinned(sh, program(*sh, {0}, {2}, {dag_route(*sh, 0, {0, 1}), dag_route(*sh, 0, {0, 2})}), LOSSY);
  const auto r5 = shared_e0->place(10);
  EXPECT_EQ(r5.call_delivered, 4);
  EXPECT_EQ(delivered(*shared_e0.p), (Amounts {2, 2}));
  // Optimized over one route DAG is CentralizedTE(support): 12, cost 34, one record per source with the
  // used edges, no entry.
  const auto es = graph(4, kES);
  auto m6 = make(es, {0}, {3}, global(std::nullopt, dag_of(*es, {0, 1, 2, 3})));
  const auto r6 = m6->place(INF);
  EXPECT_EQ(r6.call_delivered, 12);
  EXPECT_EQ(r6.total_cost, c(34));
  const auto recs = m6->records();
  ASSERT_EQ(recs.size(), 1U);
  EXPECT_FALSE(recs[0].entry_id.has_value());
  EXPECT_EQ(ascending_edges(*recs[0].route), (std::vector<EdgeId> {0, 1, 2, 3}));
  // Lossy +inf is rejected.
  auto m7 = pinned(d, eb, LOSSY);
  EXPECT_EQ(error_of([&] { (void)m7->place(INF); }), kFinite);
  // Two sources are accepted; duplicates are distinct entries.
  const auto rel = graph(4, kRel);
  EXPECT_EQ(pinned(rel, tunnels(*rel, {0, 1}, {3}, {{0, {0, 2}, 1}, {1, {1, 2}, 1}}))->place(INF).call_delivered, 10);
  const auto par = par2(4, 8);
  EXPECT_EQ(bound(par, {{0, {0}, 1}, {0, {1}, 1}, {0, {1}, 1}}), 12);
  // An empty route list.
  EXPECT_EQ(error_of([&] { (void)program(*par, {0}, {1}, {}); }), "routes must not be empty");
}

// ---------------------------------------------------------------- storage of the new layouts

// A plan's replay under a mask (restrict_prepared_route) is a new plan whose views point into its own block:
// it keeps its values after the source plan is destroyed and plans of the same shape reuse the
// freed memory (a borrowing restriction would read the later plans' roots and destinations).
TEST(ProgramLayouts, RestrictedPlanOwnsItsBlock) {
  const auto g = graph(8, {{0, 1, 5, 1}, {0, 2, 5, 1}, {1, 3, 5, 1}, {2, 3, 5, 1},
                           {4, 5, 5, 1}, {4, 6, 5, 1}, {5, 7, 5, 1}, {6, 7, 5, 1}});
  const auto mask = g->eligibility({}, std::vector<EdgeId> {3, 7});
  auto plan_of = [&](NodeId root, NodeId dst, std::vector<EdgeId> edges) {
    const auto dag = detail::PredecessorDAGBuilder::of_edges(*g, edges);
    return detail::prepare_route(*g, EB, LOSSLESS, detail::WeightSource {}, root, std::vector<NodeId> {dst}, dag);
  };
  std::optional<detail::PreparedRoute> source(plan_of(0, 3, {0, 1, 2, 3}));
  const auto restricted = detail::restrict_prepared_route(*g, *source, v({0}), v({3}), &mask);
  EXPECT_NE(restricted.sources.data(), source->sources.data());
  EXPECT_NE(restricted.destinations.data(), source->destinations.data());
  EXPECT_NE(restricted.topological_order.data(), source->topological_order.data());
  EXPECT_NE(restricted.edge_destination_positions.data(), source->edge_destination_positions.data());
  const auto roots = vec(restricted.sources);
  const auto dsts = vec(restricted.destinations);
  const auto topo = vec(restricted.topological_order);
  const auto members = restricted.forwarding_edges;
  EXPECT_EQ(members, (std::vector<EdgeId> {0, 2}));
  EXPECT_EQ(roots, v({0}));
  EXPECT_EQ(dsts, v({3}));
  source.reset();
  std::vector<detail::PreparedRoute> later;
  for (int i = 0; i < 4; ++i) {
    later.push_back(plan_of(4, 7, {4, 5, 6, 7}));
    later.push_back(detail::restrict_prepared_route(*g, later.back(), v({4}), v({7}), &mask));
  }
  EXPECT_EQ(vec(restricted.sources), roots);
  EXPECT_EQ(vec(restricted.destinations), dsts);
  EXPECT_EQ(vec(restricted.topological_order), topo);
  EXPECT_EQ(restricted.forwarding_edges, members);
  // Its copy owns another block too.
  const detail::PreparedRoute copy(restricted);
  EXPECT_NE(copy.topological_order.data(), restricted.topological_order.data());
  EXPECT_EQ(vec(copy.topological_order), topo);
}

// A program and its data are one allocation (detail::ProgramBlock): the program's data pointer
// reads into the block that every shared reference keeps alive, not into the data it was moved
// from, and stays valid after the producer's locals are gone and programs of the same shape reuse
// the freed memory.
// Every list of a program's data lies in the data's own block, none outside it.
void expect_lists_in_block(const detail::ProgramData& data) {
  const auto block = data.block();
  auto inside = [&](std::span<const std::int32_t> list, const char* name) {
    if (list.empty()) return;
    EXPECT_TRUE(list.data() >= block.data() && list.data() + list.size() <= block.data() + block.size()) << name;
  };
  inside(data.sources, "sources");
  inside(data.destinations, "destinations");
  inside(data.reachable_sources, "roster");
  inside(data.source_entry_offsets, "slot_first");
  inside(data.source_entry_positions, "slot_entries");
  inside(data.tunnel_entry_positions, "tunnel_entry");
  inside(data.excluded_nodes, "excluded_nodes");
  inside(data.excluded_edges, "excluded_edges");
  for (const auto& c : data.split_contexts) inside(c.source_entry_positions, "root_entry");
}

TEST(ProgramLayouts, ProgramOwnsItsDataBlock) {
  const auto g = graph(4, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 10, 1}, {2, 3, 10, 1}});
  std::shared_ptr<const ForwardingProgram> prog;
  {
    auto made = tunnels(*g, {0}, {3}, {{0, {0, 2}, 2}, {0, {1, 3}, 1}});
    prog = made;
  }
  const auto& data = detail::ProgramAccess::data(*prog);
  const auto* block_lo = reinterpret_cast<const char*>(prog.get());
  const auto* data_at = reinterpret_cast<const char*>(&data);
  // The data lies in the same allocation as the program (one block of the two objects).
  EXPECT_LT(std::abs(data_at - block_lo), static_cast<std::ptrdiff_t>(sizeof(detail::ProgramBlock)));
  std::vector<std::shared_ptr<const ForwardingProgram>> later;
  for (int i = 0; i < 8; ++i) later.push_back(tunnels(*g, {0}, {3}, {{0, {1, 3}, 5}}));
  ASSERT_EQ(prog->tunnels().size(), 2U);
  EXPECT_EQ(prog->tunnels()[0].edges, v({0, 2}));
  EXPECT_EQ(prog->tunnels()[1].weight, 1);
  EXPECT_EQ(vec(prog->sources()), v({0}));
  EXPECT_EQ(vec(data.reachable_sources), v({0}));
  EXPECT_EQ(data.entries.size(), 2U);
  EXPECT_EQ(vec(data.entries_for_source(0)), v({0, 1}));
  expect_lists_in_block(data);
  // A copy of the shared reference reads the same block.
  const auto copy = prog;
  prog.reset();
  EXPECT_EQ(&detail::ProgramAccess::data(*copy), &data);
  EXPECT_EQ(copy->tunnels()[1].edges, v({1, 3}));
}

// A controller's program owns its lists: it outlives the policy whose roster and exclusions it
// copied, and reads the same values from its own block.
TEST(ProgramLayouts, ControllerProgramOwnsItsLists) {
  const auto g = graph(5, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 10, 1}, {2, 3, 10, 1}, {4, 3, 10, 1}});
  FlowLedger fg(g);
  const auto el = std::make_shared<const Eligibility>(g->eligibility(std::vector<NodeId> {4}, std::vector<EdgeId> {1}));
  std::shared_ptr<const ForwardingProgram> prog;
  {
    auto p = std::make_unique<FlowPolicy>(fg, v({0, 4}), v({3}), global(std::nullopt), el.get());
    prog = p->place(5).program;
    p.reset();
  }
  // Churn the allocator so a list read from freed policy storage would differ.
  std::vector<std::vector<NodeId>> churn;
  for (int i = 0; i < 64; ++i) churn.emplace_back(8, -7);
  ASSERT_NE(prog, nullptr);
  const auto& data = detail::ProgramAccess::data(*prog);
  EXPECT_EQ(vec(prog->sources()), v({0, 4}));
  EXPECT_EQ(vec(prog->destinations()), v({3}));
  EXPECT_EQ(vec(data.reachable_sources), v({0}));
  EXPECT_EQ(vec(data.excluded_nodes), v({4}));
  EXPECT_EQ(vec(data.excluded_edges), v({1, 4}));
  EXPECT_EQ(vec(data.tunnel_entry_positions), v({0}));
  EXPECT_EQ(vec(data.source_entry_offsets), v({0, 1}));
  EXPECT_EQ(vec(data.source_entry_positions), v({0}));
  expect_lists_in_block(data);
}

// A forwarding view keeps its program alive (its data pointer reads into it), owns its vectors,
// and a view with other black holes shares the program and its structure: both stay valid after
// the source view and every other reference to the program are gone and programs of the same
// shape reuse the memory.
TEST(ProgramLayouts, ViewOwnsItsStorage) {
  const auto g = graph(6, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 2, 10, 1}, {3, 4, 10, 1}, {3, 5, 10, 1}, {4, 5, 10, 1}});
  auto build = [&](NodeId s, std::vector<EdgeId> a, std::vector<EdgeId> b, NodeId dst) {
    return program(*g, {s}, {dst}, {tunnel_route(*g, s, a), dag_route(*g, s, b)});
  };
  const auto dark = g->eligibility({}, std::vector<EdgeId> {2});
  std::optional<detail::ForwardingView> source;
  {
    auto prog = build(0, {1}, {0, 2}, 2);
    source.emplace(detail::make_forwarding_view(*g, prog, nullptr, v({2}), nullptr, true));
  }
  ASSERT_EQ(source->num_sources(), 1U);
  ASSERT_FALSE(source->borrowed);  // a tunnel's record route is made: the view holds its state
  const detail::ForwardingView copy(*source);
  const auto stale = detail::with_black_holes(*source, &dark);
  EXPECT_NE(copy.status_.data(), source->status_.data());
  EXPECT_NE(copy.source_entry_positions_.data(), source->source_entry_positions_.data());
  EXPECT_EQ(copy.program, source->program);
  EXPECT_EQ(stale.program, source->program);
  EXPECT_EQ(stale.black_holes, &dark);
  const auto status = copy.status_;
  const auto entries_of_head = copy.source_entry_positions_;
  const auto plans = copy.prepared_routes_.size();
  source.reset();
  std::vector<std::shared_ptr<const ForwardingProgram>> later;
  for (int i = 0; i < 4; ++i) later.push_back(build(3, {4}, {3, 5}, 5));
  for (const auto* view : {&copy, &stale}) {
    ASSERT_NE(view->program, nullptr);
    EXPECT_EQ(view->data, &detail::ProgramAccess::data(*view->program));
    EXPECT_EQ(view->source(0), 0);
    EXPECT_EQ(vec(view->data->reachable_sources), v({0}));
    EXPECT_EQ(vec(view->data->destinations), v({2}));
    expect_lists_in_block(*view->data);
    EXPECT_EQ(view->status_, status);
    EXPECT_EQ(view->source_entry_positions_, entries_of_head);
    ASSERT_EQ(view->prepared_routes_.size(), plans);
    for (const auto& plan : view->prepared_routes_) EXPECT_EQ(vec(plan->destinations), v({2}));
  }
}

// A borrowed view (no mask, no black holes, no record route to make) holds no state of its own:
// it reads the program, which its shared reference keeps alive after every other reference is
// gone and programs of the same shape reuse the freed memory; a copy reads the same program, and
// the same view with black holes writes the program's state out before changing it.
TEST(ProgramLayouts, BorrowedViewReadsItsProgram) {
  const auto g = graph(4, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 10, 1}, {2, 3, 10, 1}});
  auto build = [&](NodeId s) { return program(*g, {s}, {3}, {dag_route(*g, s, {0, 1, 2, 3})}); };
  const auto dark = g->eligibility({}, std::vector<EdgeId> {2});
  std::optional<detail::ForwardingView> source;
  {
    auto prog = build(0);
    source.emplace(detail::make_forwarding_view(*g, prog, nullptr, v({3}), nullptr, true));
  }
  ASSERT_TRUE(source->borrowed);
  EXPECT_TRUE(source->status_.empty() && source->prepared_routes_.empty() && source->source_entry_positions_.empty());
  const detail::ForwardingView copy(*source);
  const auto stale = detail::with_black_holes(*source, &dark);
  EXPECT_TRUE(copy.borrowed);
  EXPECT_FALSE(stale.borrowed);
  const auto* plan = source->prepared_route(0).get();
  source.reset();
  std::vector<std::shared_ptr<const ForwardingProgram>> later;
  for (int i = 0; i < 4; ++i) later.push_back(build(0));
  for (const auto* view : {&copy, &stale}) {
    ASSERT_NE(view->program, nullptr);
    EXPECT_EQ(view->data, &detail::ProgramAccess::data(*view->program));
    ASSERT_EQ(view->num_sources(), 1U);
    EXPECT_EQ(view->source(0), 0);
    EXPECT_EQ(view->status(0), TunnelStatus::Up);
    EXPECT_EQ(vec(view->entries_for_source(0)), v({0}));
    EXPECT_EQ(view->context_source_index(0), 0);
    EXPECT_EQ(view->prepared_route(0).get(), plan);
    EXPECT_EQ(view->prepared_route(0)->forwarding_edges, v({0, 1, 2, 3}));
    EXPECT_EQ(view->min_route_cost(0), 2);
    expect_lists_in_block(*view->data);
  }
  EXPECT_FALSE(copy.is_blackholed_source(0));
}

// An epoch is a value: its copy keeps the cohorts, the program, the rows and the counters after
// the policy that held the original has moved on and been destroyed.
TEST(ProgramLayouts, EpochCopyOwnsItsStorage) {
  const auto g = graph(4, kPinloss);
  const auto prog = tunnels(*g, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}});
  std::optional<Pol> p(pinned(g, prog, LOSSY));
  (void)(*p)->replace(12);
  const detail::PolicyEpoch copy = detail::PolicyAccess::engine(*p->p).epoch;
  const auto& live = detail::PolicyAccess::engine(*p->p).epoch;
  ASSERT_EQ(copy.record_groups.size(), 2U);
  EXPECT_NE(copy.record_groups.data(), live.record_groups.data());
  EXPECT_NE(copy.by_tunnel.data(), live.by_tunnel.data());
  const auto rows_before = copy.by_tunnel;
  const auto ids_before = std::vector<std::int64_t> {copy.record_groups[0][0].record_counter, copy.record_groups[1][0].record_counter};
  (void)(*p)->replace(4);
  p.reset();
  for (int i = 0; i < 4; ++i) (void)pinned(g, prog, LOSSY)->replace(2);
  EXPECT_EQ(copy.program, prog);
  EXPECT_TRUE(copy.by_tunnel == rows_before);
  EXPECT_EQ(rows_before.size(), 2U);
  EXPECT_EQ(lo(rows_before[0].assigned), 6);
  EXPECT_EQ(rows_before[0].delivered, 2);
  EXPECT_EQ((std::vector<std::int64_t> {copy.record_groups[0][0].record_counter, copy.record_groups[1][0].record_counter}), ids_before);
  EXPECT_EQ(lo(copy.offered), 12);
  EXPECT_EQ(lo(copy.requested), 12);
}

// ---------------------------------------------------------------- property: random tunnel programs

// Random tunnel programs on random DAGs (two heads, one destination, paths drawn from every path
// of each head, weights 0 to 3 with a positive one per head) against program_reference: the
// lossless bound (x fits, x + 1 does not), each tunnel's highest-averages amount and the loads at
// the bound, conservation at every node, a Lossy replace at the bound delivering every share, and
// evaluate above the bound returning the bound.
TEST(ProgramProperty, TunnelProgramsMatchTheReference) {
  namespace ref = program_reference;
  std::mt19937 rng(20261018);
  int checked = 0;
  int two_heads = 0;  // trials whose roster has both heads
  int positive = 0;   // trials with a positive bound
  int shared_links = 0;  // trials where some edge carries two tunnels' amounts
  for (int trial = 0; trial < 60; ++trial) {
    const std::int32_t n = 6;
    std::vector<InEdge> rows;
    for (NodeId u = 0; u < n; ++u) {
      for (NodeId w = u + 1; w < n; ++w) {
        if (u == 0 && w == 1) continue;  // the heads do not feed each other
        const int copies = static_cast<int>(rng() % 3);
        for (int k = 0; k < copies; ++k) rows.push_back({u, w, static_cast<Amount>(rng() % 7), 1});
      }
    }
    const auto g = shared(make_graph(n, rows));
    const NodeId dst = n - 1;
    // Every path from a head to dst (ids of g; destinations are terminal).
    std::vector<std::vector<EdgeId>> found[2];
    for (NodeId h = 0; h < 2; ++h) {
      std::vector<EdgeId> stack;
      std::function<void(NodeId)> walk = [&](NodeId u) {
        if (u == dst) {
          found[h].push_back(stack);
          return;
        }
        for (EdgeId e = 0; e < g->num_edges(); ++e) {
          if (g->edge_sources()[static_cast<std::size_t>(e)] != u) continue;
          stack.push_back(e);
          walk(g->edge_destinations()[static_cast<std::size_t>(e)]);
          stack.pop_back();
        }
      };
      walk(h);
    }
    std::vector<ref::RefTunnel> t;
    for (NodeId h = 0; h < 2; ++h) {
      if (found[h].empty()) continue;
      const int k = 1 + static_cast<int>(rng() % 3);
      for (int i = 0; i < k; ++i) {
        t.push_back({h, found[h][rng() % found[h].size()], static_cast<Amount>(rng() % 4)});
      }
      t.back().weight = std::max<Amount>(t.back().weight, 1);
    }
    if (t.empty()) continue;
    std::vector<Spec> specs;
    for (const auto& x : t) specs.emplace_back(x.head, x.edges, x.weight);
    const auto prog = tunnels(*g, {0, 1}, {dst}, specs);
    const Amount x = ref::lossless_bound(*g, t);
    const std::string where = "trial " + std::to_string(trial) + " bound " + std::to_string(x);
    ASSERT_FALSE(ref::fits(*g, t, x + 1)) << where;
    auto p = pinned(g, prog);
    const auto r = p->place(INF, kEdges);
    ASSERT_EQ(r.call_delivered, x) << where;
    const auto amounts = *ref::tunnel_amounts(x, t);
    two_heads += ref::roster(t).size() == 2 ? 1 : 0;
    positive += x > 0 ? 1 : 0;
    {
      std::vector<int> users(static_cast<std::size_t>(g->num_edges()), 0);
      for (std::size_t i = 0; i < t.size(); ++i) {
        if (amounts[i] > 0) {
          for (EdgeId e : t[i].edges) ++users[static_cast<std::size_t>(e)];
        }
      }
      shared_links += std::any_of(users.begin(), users.end(), [](int u) { return u > 1; }) ? 1 : 0;
    }
    EXPECT_EQ(*r.edge_flows, ref::loads(*g, t, amounts)) << where;
    EXPECT_EQ(delivered(*p.p), amounts) << where;
    ASSERT_EQ(r.by_tunnel.size(), t.size()) << where;
    for (std::size_t i = 0; i < t.size(); ++i) EXPECT_EQ(r.by_tunnel[i].delivered, amounts[i]) << where << " tunnel " << i;
    // Conservation: a transit node passes on what it receives; the heads send their shares.
    std::vector<Amount> in(static_cast<std::size_t>(n), 0), out(static_cast<std::size_t>(n), 0);
    for (EdgeId e = 0; e < g->num_edges(); ++e) {
      out[static_cast<std::size_t>(g->edge_sources()[e])] += (*r.edge_flows)[e];
      in[static_cast<std::size_t>(g->edge_destinations()[e])] += (*r.edge_flows)[e];
    }
    for (NodeId u = 2; u < dst; ++u) EXPECT_EQ(in[u], out[u]) << where << " node " << u;
    EXPECT_EQ(in[dst], x) << where;
    EXPECT_EQ(out[0] + out[1] - in[1], x) << where;
    // Lossy at the bound: every share delivered, nothing dropped.
    const auto l = pinned(g, prog, LOSSY)->replace(x);
    EXPECT_EQ(l.call_delivered, x) << where;
    EXPECT_TRUE(l.call_edge_drops.empty()) << where;
    EXPECT_TRUE(l.source_drops.is_zero()) << where;
    // Lossless above the bound places the bound.
    EXPECT_EQ(pinned(g, prog)->evaluate(x + 1).call_delivered, x) << where;
    ++checked;
  }
  RecordProperty("checked", checked);
  RecordProperty("two_heads", two_heads);
  RecordProperty("positive", positive);
  RecordProperty("shared_links", shared_links);
  EXPECT_GT(checked, 30);  // not vacuous
  EXPECT_GT(two_heads, 10);
  EXPECT_GT(positive, 20);
  EXPECT_GT(shared_links, 5);
}

// The reference's own controls: a known bound and a lost share.
TEST(ProgramProperty, ReferenceControls) {
  namespace ref = program_reference;
  EXPECT_EQ(ref::round_robin(5, 2), (Amounts {3, 2}));
  EXPECT_EQ(ref::highest_averages(15, {1, 2}), (Amounts {5, 10}));
  EXPECT_EQ(ref::highest_averages(2, {1, 1}), (Amounts {1, 1}));
  EXPECT_EQ(ref::highest_averages(4, {0, 3}), (Amounts {0, 4}));
  const auto g = par2(6, 4);
  EXPECT_EQ(ref::lossless_bound(*g, {{0, {0}, 1}, {0, {1}, 1}}), 9);
  EXPECT_EQ(ref::lossless_bound(*g, {{0, {0}, 1}, {0, {1}, 2}}), 6);  // H(6, (1, 2)) = (2, 4); H(7) = (2, 5)
  EXPECT_FALSE(ref::tunnel_amounts(1, {{0, {0}, 0}}).has_value());
}

// ---------------------------------------------------------------- from_dag

// DP: links 0->2 and 1->2 of capacity 10, then two parallel links 2->3 of capacity 1, cost 1;
// sources 0 and 1, destination 3.
const std::vector<InEdge> kDP = {{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 1, 1}, {2, 3, 1, 1}};

std::shared_ptr<const ForwardingProgram> from_dag(const Graph& g, const std::vector<EdgeId>& edges, const RouteOptions& o = {}) {
  return ForwardingProgram::from_dag(g, v({0, 1}), v({3}), dag_of(g, edges), o);
}

// DAG-pooled against DAG-separate: one context splits node 2's two units once; two contexts each
// put their unit on the first link.
TEST(FromDag, Pools) {
  const auto g = graph(4, kDP);
  const auto pooled = from_dag(*g, {0, 1, 2, 3});
  EXPECT_EQ(pinned(g, pooled)->place(INF).call_delivered, 2);
  EXPECT_EQ(pinned(g, pooled, LOSSY)->place(2).call_delivered, 2);
  const auto separate = ForwardingProgram::from_routes(*g, v({0, 1}), v({3}),
                                                       std::vector<RouteDefinition> {dag_route(*g, 0, {0, 1, 2, 3}), dag_route(*g, 1, {0, 1, 2, 3})});
  EXPECT_EQ(pinned(g, separate)->place(INF).call_delivered, 1);
  EXPECT_EQ(pinned(g, separate, LOSSY)->place(2).call_delivered, 1);
  EXPECT_TRUE(pooled->tunnels().empty());
}

// LOSSY 3: node 2 splits RR(3, 2) = (2, 1), e2's 2 units go back by H(2, [2, 1]) = (2, 0).
TEST(FromDag, Attribution) {
  const auto g = graph(4, kDP);
  auto p = pinned(g, from_dag(*g, {0, 1, 2, 3}), LOSSY);
  const auto r = p->place(3);
  EXPECT_EQ(r.call_delivered, 2);
  EXPECT_EQ(r.call_edge_drops, (Drops {{2, 1}}));
  ASSERT_EQ(r.by_source.size(), 2U);
  EXPECT_EQ(std::tuple(lo(r.by_source[0].offered), r.by_source[0].delivered, lo(r.by_source[0].total_edge_drops)),
            std::tuple(Amount {2}, Amount {1}, Amount {1}));
  EXPECT_EQ(std::tuple(lo(r.by_source[1].offered), r.by_source[1].delivered, lo(r.by_source[1].total_edge_drops)),
            std::tuple(Amount {1}, Amount {1}, Amount {0}));
  const auto records = p.ledger->records();
  ASSERT_EQ(records.size(), 2U);
  EXPECT_EQ(records[0].edge_flows, (Drops {{0, 2}, {2, 1}}));
  EXPECT_EQ(records[1].edge_flows, (Drops {{1, 1}, {3, 1}}));
}

// A residual snapshot is read per EdgeId, not per position in the DAG: e0 0->2 is outside the DAG,
// e1 and e2 0->1 (cost 2) are in it; residual {10, 2, 8} weighs them 2 and 8, so an offer of 10
// splits by highest averages as (2, 8) and e0 carries nothing.
TEST(FromDag, SnapshotReadsByEdgeId) {
  const auto g = graph(3, {{0, 2, 10, 1}, {0, 1, 10, 2}, {0, 1, 10, 2}});
  RouteOptions o;
  o.allocation = W;
  o.weight_basis = WeightBasis::ResidualSnapshot;
  o.residual = Amounts {10, 2, 8};
  const auto program = ForwardingProgram::from_dag(*g, v({0}), v({1}), dag_of(*g, {1, 2}), o);
  const auto r = pinned(g, program)->replace(10, std::nullopt, kEdges);
  EXPECT_EQ(r.edge_flows, (Amounts {0, 2, 8}));
}

// ShortestPaths ECMP and Pinned over from_dag of the same DAG: identical delivered, drops and rows
// for Lossy offers 1 to 5 and the Lossless bound; the same under a mask, and WCMP on DP3.
TEST(FromDag, EqualsShortestPaths) {
  const auto check = [](const std::shared_ptr<const Graph>& g, RouteAllocation k, std::shared_ptr<const Eligibility> mask) -> Amount {
    RouteOptions o;
    o.allocation = k;
    const auto program = from_dag(*g, {0, 1, 2, 3}, o);
    for (const Amount offer : {1, 2, 3, 4, 5}) {
      auto a = make(g, v({0, 1}), v({3}), cost_only(k, LOSSY), mask);
      auto b = make(g, v({0, 1}), v({3}), FlowPolicyConfig {.routing = Pinned {program, nullptr}, .admission = LOSSY}, mask);
      const auto x = a->place(offer);
      const auto y = b->place(offer);
      EXPECT_EQ(x.call_delivered, y.call_delivered) << offer;
      EXPECT_EQ(x.call_edge_drops, y.call_edge_drops) << offer;
      EXPECT_EQ(x.by_source.size(), y.by_source.size());
      for (std::size_t i = 0; i < std::min(x.by_source.size(), y.by_source.size()); ++i) {
        EXPECT_EQ(x.by_source[i].delivered, y.by_source[i].delivered) << offer;
        EXPECT_EQ(x.by_source[i].total_edge_drops, y.by_source[i].total_edge_drops) << offer;
      }
    }
    auto a = make(g, v({0, 1}), v({3}), cost_only(k), mask);
    auto b = make(g, v({0, 1}), v({3}), FlowPolicyConfig {.routing = Pinned {program, nullptr}}, mask);
    const Amount bound = a->place(INF).call_delivered;
    EXPECT_EQ(bound, b->place(INF).call_delivered);
    return bound;
  };
  const auto dp = graph(4, kDP);
  EXPECT_EQ(check(dp, EB, nullptr), 2);
  EXPECT_EQ(check(dp, EB, no_edges(*dp, {3})), 1);
  EXPECT_EQ(check(graph(4, {{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 1, 1}, {2, 3, 3, 1}}), W, nullptr), 4);
}

// The validation order of from_dag: the endpoints, the DAG's graph, the allocation, the weights.
TEST(FromDag, Validation) {
  const auto g = graph(4, kDP);
  const auto other = graph(4, kDP);
  const auto whole = dag_of(*g, {0, 1, 2, 3});
  const auto foreign = dag_of(*other, {0, 1, 2, 3});
  RouteOptions optimized;
  optimized.allocation = MF;
  RouteOptions basis;
  basis.weight_basis = WeightBasis::Capacity;
  RouteOptions short_weights;
  short_weights.allocation = W;
  short_weights.weight_basis = WeightBasis::Explicit;
  short_weights.edge_weights = std::vector<Amount> {1, 1, 1};
  const auto err = [&](std::vector<NodeId> s, const std::shared_ptr<const PredecessorDAG>& d, const RouteOptions& o) {
    return error_of([&] { (void)ForwardingProgram::from_dag(*g, s, v({3}), d, o); });
  };
  EXPECT_EQ(err({}, foreign, optimized), "sources: the node set is empty");
  EXPECT_EQ(err({0, 1}, foreign, optimized), "dag belongs to another graph");
  EXPECT_EQ(err({0, 1}, nullptr, {}), "dag belongs to another graph");
  EXPECT_EQ(err({0, 1}, whole, optimized), "from_dag allocation must be EQUAL or WEIGHTED");
  EXPECT_EQ(err({0, 1}, whole, basis), "weight_basis requires allocation=WEIGHTED");
  EXPECT_EQ(err({0, 1}, whole, short_weights), "edge_weights: expected 4 entries");
  // An empty roster is a valid program.
  const auto none = ForwardingProgram::from_dag(*g, v({0, 1}), v({3}), dag_of(*g, {2, 3}));
  EXPECT_TRUE(none->tunnels().empty());
  EXPECT_EQ(pinned(g, none)->place(INF).termination, Termination::NoRoute);
}

// ---------------------------------------------------------------- the forwarding network's size

// Each count of a forwarding network is checked against 2^31 - 1 inclusively, with its exact value
// and its name in the text. 10,000 tunnels of 214,749 hops (2,147,490,000
// occurrences) are refused; exactly 2^31 - 1 passes; 2^31 entries are refused as entries.
TEST(ForwardingSize, EachCountIsCheckedInclusivelyWithItsExactValue) {
  using detail::ForwardingCounter;
  EXPECT_EQ(detail::kForwardingLimit, 2147483647);
  for (const auto c : {ForwardingCounter::Entry, ForwardingCounter::Occurrence, ForwardingCounter::Router,
                       ForwardingCounter::Event, ForwardingCounter::Dependency}) {
    EXPECT_EQ(error_of([&] { detail::check_forwarding_count(c, 2147483647); }), "");
    EXPECT_EQ(error_of([&] { detail::check_forwarding_count(c, 0); }), "");
  }
  EXPECT_EQ(error_of([&] { detail::check_forwarding_count(ForwardingCounter::Occurrence, std::int64_t {10000} * 214749); }),
            "forwarding occurrence limit exceeded: 2147490000 > 2147483647");
  EXPECT_EQ(error_of([&] { detail::check_forwarding_count(ForwardingCounter::Entry, std::int64_t {1} << 31); }),
            "forwarding entry limit exceeded: 2147483648 > 2147483647");
  EXPECT_EQ(error_of([&] { detail::check_forwarding_count(ForwardingCounter::Router, std::int64_t {1} << 31); }),
            "forwarding router limit exceeded: 2147483648 > 2147483647");
  EXPECT_EQ(error_of([&] { detail::check_forwarding_count(ForwardingCounter::Event, std::int64_t {1} << 31); }),
            "forwarding event limit exceeded: 2147483648 > 2147483647");
  EXPECT_EQ(error_of([&] { detail::check_forwarding_count(ForwardingCounter::Dependency, std::int64_t {1} << 31); }),
            "forwarding dependency limit exceeded: 2147483648 > 2147483647");
}

// A call counts its network and checks the entry, occurrence and router counts before the network
// is allocated, the event and dependency counts before the event structures: with one counter's
// limit lowered to its count less one for this thread (a hooks build), the call is refused with
// that count, nothing changed; at the count it runs; with every limit at 1 the entries are
// reported first. PAR2(6, 4) with two one-edge route DAGs (two contexts: the general network),
// LOSSY replace(12): 2 entries, 2 occurrences (one member each), 4 routers (two per route), 6
// units sent on each link; the 4-unit link can drop (one barrier), so 1 + 4 = 5 events and 3
// dependencies (router to router on the 6-unit link; router to event to router on the other).
// LOSSLESS has no events: its entries, occurrences and routers are checked. DIAM's two 2-hop
// tunnels: 2 entries, 4 occurrences.
TEST(ForwardingSize, TheCallChecksEachCountInOrder) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  using detail::ForwardingCounter;
  struct Restore {
    ~Restore() { detail::thread_forwarding_limits() = {}; }
  } restore;
  const auto p2 = par2(6, 4);
  const auto two_dags = program(*p2, {0}, {1}, {dag_route(*p2, 0, {0}), dag_route(*p2, 0, {1})});
  const auto d = graph(4, kDiam);
  const auto two_tunnels = tunnels(*d, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}});
  struct Count {
    ForwardingCounter counter;
    const char* name;
    std::int64_t n;
  };
  struct Case {
    std::shared_ptr<const Graph> g;
    std::shared_ptr<const ForwardingProgram> prog;
    Admission admission;
    Amount delivered;
    std::vector<Count> counts;
  };
  const std::vector<Case> cases {
      {p2, two_dags, LOSSY, 10,
       {{ForwardingCounter::Entry, "entry", 2}, {ForwardingCounter::Occurrence, "occurrence", 2},
        {ForwardingCounter::Router, "router", 4}, {ForwardingCounter::Event, "event", 5},
        {ForwardingCounter::Dependency, "dependency", 3}}},
      {p2, two_dags, LOSSLESS, 9,
       {{ForwardingCounter::Entry, "entry", 2}, {ForwardingCounter::Occurrence, "occurrence", 2},
        {ForwardingCounter::Router, "router", 4}}},
      {d, two_tunnels, LOSSY, 12, {{ForwardingCounter::Entry, "entry", 2}, {ForwardingCounter::Occurrence, "occurrence", 4}}},
      {d, two_tunnels, LOSSLESS, 12, {{ForwardingCounter::Entry, "entry", 2}, {ForwardingCounter::Occurrence, "occurrence", 4}}},
  };
  for (const auto& c : cases) {
    auto pol = pinned(c.g, c.prog, c.admission);
    for (const auto& k : c.counts) {
      const LedgerState before = ledger_state(*pol.ledger);
      auto& limits = detail::thread_forwarding_limits();
      limits = {};
      limits[static_cast<std::size_t>(k.counter)] = k.n - 1;
      EXPECT_EQ(error_of([&] { (void)pol->replace(12); }), std::string("forwarding ") + k.name + " limit exceeded: " +
                                                                std::to_string(k.n) + " > " + std::to_string(k.n - 1));
      expect_same(ledger_state(*pol.ledger), before, k.name);
      EXPECT_TRUE(pol->records().empty());
      limits[static_cast<std::size_t>(k.counter)] = k.n;
      EXPECT_EQ(pol->replace(12).call_delivered, c.delivered) << k.name;
      pol->remove();
    }
    detail::thread_forwarding_limits() = {1, 1, 1, 1, 1};
    EXPECT_EQ(error_of([&] { (void)pol->replace(12); }), "forwarding entry limit exceeded: 2 > 1");
    detail::thread_forwarding_limits() = {};
  }
#else
  GTEST_SKIP() << "a lowered limit needs the private hooks";
#endif
}
