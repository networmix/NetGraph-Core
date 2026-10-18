// Forwarding programs (step 3), natively: what every model exports, Global with a support, the
// weight bases, ForwardingProgram::from_routes, Pinned over programs, the lossless bound, the loss
// rule inside one program (cycles included), replay under a failure mask, stale forwarding with
// black holes, the operations evaluate, replace and forward with their epochs and flags, the
// per-entry records, the by_tunnel rows and their 128-bit counters, the evaluators' offer, and
// today's Pinned mapped onto programs. Each test checks the hand-derived values of the step's
// expected-result rows (LEDGER rows P, G, W, F, C, L, Y, R, S, O, V, D, Z, E, M, named in each
// test) on small graphs whose EdgeIds equal their input positions (asserted); none is read from
// a build. Then the storage of the new layouts (a restricted plan, a forwarding view, an epoch)
// against borrowing copies, and a property test of random tunnel programs against the
// independent reference of test_utils.hpp (program_reference).
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/detail/forwarding.hpp"
#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/program.hpp"
#include "state_comparators.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

constexpr FlowPlacement EB = FlowPlacement::EqualBalanced;
constexpr FlowPlacement W = FlowPlacement::Weighted;
constexpr FlowPlacement MF = FlowPlacement::MaxFlow;
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
  for (std::size_t i = 0; i < rows.size(); ++i) EXPECT_EQ(g->input_index()[i], static_cast<std::int32_t>(i));
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
Route tunnel_route(const Graph& g, NodeId head, const std::vector<EdgeId>& edges, Amount weight = 1) {
  return Route {head, path(g, edges), weight};
}
Route dag_route(const Graph& g, NodeId head, const std::vector<EdgeId>& edges, Amount weight = 1) {
  return Route {head, dag_of(g, edges), weight};
}
std::shared_ptr<const ForwardingProgram> program(const Graph& g, const std::vector<NodeId>& s,
                                                 const std::vector<NodeId>& d, const std::vector<Route>& routes,
                                                 const RouteOptions& options = {}) {
  return ForwardingProgram::from_routes(g, s, d, routes, options);
}
// A supplied program of tunnels.
std::shared_ptr<const ForwardingProgram> tunnels(const Graph& g, const std::vector<NodeId>& s,
                                                 const std::vector<NodeId>& d, const std::vector<Spec>& specs) {
  std::vector<Route> routes;
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

FlowPolicyConfig cost_only(FlowPlacement placement = EB, Admission admission = LOSSLESS) {
  return FlowPolicyConfig {.routing = CostOnly {}, .placement = placement, .admission = admission};
}
FlowPolicyConfig weighted(Admission admission, WeightBasis basis, std::vector<Amount> edge_weights = {}) {
  auto cfg = cost_only(W, admission);
  cfg.weights = basis;
  if (basis == WeightBasis::Explicit) cfg.edge_weights = std::make_shared<const std::vector<Amount>>(std::move(edge_weights));
  return cfg;
}
FlowPolicyConfig global(std::optional<Origin> origin = std::nullopt, std::shared_ptr<const PredDAG> support = nullptr) {
  return FlowPolicyConfig {.routing = Global {std::move(support)}, .origin = origin};
}
FlowPolicyConfig pinned_config(std::shared_ptr<const ForwardingProgram> prog, Admission admission = LOSSLESS,
                               std::shared_ptr<const Eligibility> black_holes = nullptr) {
  return FlowPolicyConfig {.routing = Pinned {std::move(prog), std::move(black_holes)}, .admission = admission};
}
FlowPolicyConfig capacity_aware() {
  return FlowPolicyConfig {.routing = CapacityAware {}, .placement = EB};
}

Pol pinned(const std::shared_ptr<const Graph>& g, const std::shared_ptr<const ForwardingProgram>& prog,
           Admission admission = LOSSLESS, std::shared_ptr<const Eligibility> elig = nullptr,
           std::shared_ptr<const Eligibility> black_holes = nullptr, std::shared_ptr<FlowLedger> ledger = nullptr) {
  return make(g, vec(prog->sources()), vec(prog->destinations()), pinned_config(prog, admission, std::move(black_holes)),
              std::move(elig), std::move(ledger));
}

Rows rows(const PlacementResult& r) {
  Rows out;
  for (const auto& t : r.by_tunnel) out.emplace_back(t.index, t.head, t.status, lo(t.assigned), t.delivered, lo(t.dropped));
  return out;
}
Listing listing(const ForwardingProgram& p) {
  Listing out;
  for (const auto& t : p.tunnels()) out.emplace_back(t.head, t.edges, t.weight);
  return out;
}
std::vector<std::optional<std::int32_t>> entries(const FlowPolicy& p) {
  std::vector<std::optional<std::int32_t>> out;
  for (const auto& r : p.records()) out.push_back(r.entry);
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

QueryOptions query(FlowPolicyConfig cfg, std::optional<std::optional<Amount>> offered = std::nullopt) {
  QueryOptions q {std::move(cfg)};
  q.offered = offered;
  return q;
}
std::optional<std::optional<Amount>> offer(Amount x) { return std::optional<std::optional<Amount>> {std::in_place, x}; }
std::optional<std::optional<Amount>> offer_inf() { return std::optional<std::optional<Amount>> {std::in_place, std::nullopt}; }

std::vector<std::pair<EdgeId, Amount>> losses(const SensitivityResult& s) {
  std::vector<std::pair<EdgeId, Amount>> out;
  for (const auto& e : s.entries) out.emplace_back(e.edge, e.loss);
  return out;
}

const std::string kFinite = "LOSSY forwarding requires a finite offer";
const std::string kNoProgram = "forward requires an installed program";
const std::string kPeriodTwo = "lossy forwarding iteration entered a period-2 cycle";
const std::string kRoundLimit = "lossy forwarding iteration did not reach a fixed point within the round limit";

} // namespace

// ---------------------------------------------------------------- P: exported programs

// P1, Z2: B1, Global optimized, place(12).
TEST(ProgramExport, GlobalOptimized) {
  const auto g = graph(4, kB1);
  auto p = make(g, {0, 1}, {3}, global());
  const auto r = p->place(12);
  ASSERT_TRUE(r.program);
  const auto& prog = *r.program;
  EXPECT_EQ(listing(prog), (Listing {{0, {0}, 2}, {1, {1, 2}, 10}}));
  std::vector<std::tuple<std::int32_t, std::int32_t, NodeId, Cost>> fields;
  for (const auto& t : prog.tunnels()) fields.emplace_back(t.index, t.local_index, t.destination, t.route_cost);
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

// P2: the equal origin admits 4 and keeps no source weights.
TEST(ProgramExport, GlobalEqual) {
  const auto g = graph(4, kB1);
  auto p = make(g, {0, 1}, {3}, global(Origin::Equal));
  const auto r = p->place(12);
  EXPECT_EQ(r.placed, 4);
  EXPECT_EQ(listing(*r.program), (Listing {{0, {0}, 2}, {1, {1, 2}, 2}}));
  EXPECT_FALSE(r.program->source_weights().has_value());
}

// P3: one source is decomposed too (R12).
TEST(ProgramExport, SingleSourceDecomposed) {
  const auto g = graph(4, kES);
  auto p = make(g, {0}, {3}, global());
  const auto r = p->place(12);
  EXPECT_EQ(listing(*r.program), (Listing {{0, {0, 2}, 2}, {0, {1, 3}, 10}}));
  std::vector<std::pair<std::int32_t, Cost>> fields;
  for (const auto& t : r.program->tunnels()) fields.emplace_back(t.local_index, t.route_cost);
  EXPECT_EQ(fields, (std::vector<std::pair<std::int32_t, Cost>> {{0, 2}, {1, 3}}));
  EXPECT_EQ(r.total_cost, c(34));
  EXPECT_EQ(*r.program->source_weights(), (std::vector<std::pair<NodeId, Amount>> {{0, 12}}));
}

// P4: CostOnly MaxFlow exports its controller's paths.
TEST(ProgramExport, CostOnlyMaxFlow) {
  const auto g = graph(4, kES);
  auto p = make(g, {0}, {3}, cost_only(MF));
  const auto r = p->place(12);
  EXPECT_EQ(r.placed, 2);
  EXPECT_EQ(listing(*r.program), (Listing {{0, {0, 2}, 2}}));
  EXPECT_EQ(*r.program->source_weights(), (std::vector<std::pair<NodeId, Amount>> {{0, 2}}));
}

// P5, Z3: a CostOnly split program lists no tunnel.
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

// P7: every roster source keeps its source weight, zero included (R20).
TEST(ProgramExport, ZeroSupplyKept) {
  const auto g = graph(3, {{0, 2, 5, 1}, {1, 2, 5, 1}});
  auto p = make(g, {0, 1}, {2}, global());
  const auto r = p->place(1);
  EXPECT_EQ(listing(*r.program), (Listing {{0, {0}, 1}}));
  EXPECT_EQ(*r.program->source_weights(), (std::vector<std::pair<NodeId, Amount>> {{0, 1}, {1, 0}}));
}

// P8: no program before the first positive call; a zero call reports none.
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

// P9: today's CapacityAware has no program (O1).
TEST(ProgramExport, CapacityAwareHasNone) {
  const auto g = graph(2, kOne);
  auto p = make(g, {0}, {1}, capacity_aware());
  const auto r = p->place(5);
  EXPECT_EQ(r.placed, 5);
  EXPECT_EQ(r.program, nullptr);
  EXPECT_EQ(p->program(), nullptr);
}

// P10: B1's optimized program forwarded by Pinned reproduces the solve.
TEST(ProgramExport, ReplayReproducesTheSolve) {
  const auto g = graph(4, kB1);
  auto solver = make(g, {0, 1}, {3}, global());
  const auto prog = solver->place(12).program;
  auto p = pinned(g, prog);
  const auto r = p->place(12, kEdges);
  EXPECT_EQ(r.placed, 12);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_EQ(*r.edge_flows, (Amounts {2, 10, 10}));
  EXPECT_EQ(pinned(g, prog)->place(INF).placed, 12);
}

// P11: the producer's mask is saved: no replay brings source 1 back (R10).
TEST(ProgramExport, MaskSaved) {
  const auto g = graph(4, kB1);
  auto masked = make(g, {0, 1}, {3}, global(), no_nodes(*g, {1}));
  const auto r = masked->place(12);
  EXPECT_EQ(r.placed, 2);
  EXPECT_EQ(listing(*r.program), (Listing {{0, {0}, 2}}));
  const auto lossy = pinned(g, r.program, LOSSY)->replace(12);
  EXPECT_EQ(lossy.placed, 2);
  ASSERT_EQ(lossy.by_source.size(), 2U);
  EXPECT_EQ(lo(lossy.by_source[0].offered), 12);
  EXPECT_EQ(lo(lossy.by_source[1].offered), 0);
}

// P12: an empty roster installs an empty program (R22).
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

// P6: the policy's program is the configured one (R14).
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

// ---------------------------------------------------------------- G: Global with a support

// GS1: only the support carries traffic; no cut, no curve.
TEST(GlobalSupport, OnlyTheSupport) {
  const auto g = graph(4, kES);
  auto p = make(g, {0}, {3}, global(std::nullopt, dag_of(*g, {1, 3})));
  const auto r = p->place(12, OutputOptions {.min_cut = true});
  EXPECT_EQ(r.placed, 10);
  EXPECT_EQ(r.total_cost, c(30));
  EXPECT_FALSE(r.min_cut.has_value());
  EXPECT_FALSE(r.cost_distribution.has_value());
  EXPECT_EQ(listing(*r.program), (Listing {{0, {1, 3}, 10}}));
  const auto full = make(g, {0}, {3}, global())->place(12);
  EXPECT_EQ(full.placed, 12);
  EXPECT_EQ(full.total_cost, c(34));
}

// GS2: e0 leads nowhere inside the support: trimmed.
TEST(GlobalSupport, Trimmed) {
  const auto g = graph(4, kES);
  const auto r = make(g, {0}, {3}, global(std::nullopt, dag_of(*g, {0, 1, 3})))->place(12, kEdges);
  EXPECT_EQ(r.placed, 10);
  EXPECT_EQ(*r.edge_flows, (Amounts {0, 10, 0, 10}));
}

// GS3: the mask removes the support: an empty roster.
TEST(GlobalSupport, MaskRemovesSupport) {
  const auto g = graph(4, kES);
  const auto r = make(g, {0}, {3}, global(std::nullopt, dag_of(*g, {1, 3})), no_edges(*g, {3}))->place(12);
  EXPECT_EQ(r.termination, Termination::NoRoute);
  EXPECT_EQ(r.placed, 0);
  ASSERT_TRUE(r.program);
  EXPECT_TRUE(r.program->tunnels().empty());
}

// GS4: the equal origin on a support: RR(4, 2) fits, RR(5, 2) does not.
TEST(GlobalSupport, EqualOrigin) {
  const auto g = graph(3, kB2);
  const auto r = make(g, {0, 1}, {2}, global(Origin::Equal, dag_of(*g, {1, 2})))->place(12);
  EXPECT_EQ(r.placed, 4);
  EXPECT_EQ(r.total_cost, c(4));
}

// GS5: every allowed link of a support is a candidate, no shortcut.
TEST(GlobalSupport, Sensitivity) {
  const auto g = graph(4, kES);
  const auto s = sensitivity(*g, v({0}), v({3}), query(global(std::nullopt, dag_of(*g, {1, 3}))), 1);
  EXPECT_EQ(s.baseline.placed, 10);
  EXPECT_EQ(losses(s), (std::vector<std::pair<EdgeId, Amount>> {{0, 0}, {1, 10}, {2, 0}, {3, 10}}));
}

// GS6, Z6: total_cost 2^70 (high half 64, low half 0).
TEST(GlobalSupport, TotalCostBeyond64Bits) {
  const auto g = graph(2, {{0, 1, Amount {1} << 40, Cost {1} << 30}});
  const auto r = make(g, {0}, {1}, global(std::nullopt, dag_of(*g, {0})))->place(Amount {1} << 40);
  EXPECT_EQ(r.placed, Amount {1} << 40);
  ASSERT_TRUE(r.total_cost.has_value());
  EXPECT_EQ(r.total_cost->hi, 64U);
  EXPECT_EQ(r.total_cost->lo, 0U);
}

// GS7: a support of another graph.
TEST(GlobalSupport, Ownership) {
  const auto g = graph(2, {{0, 1, 1, 1}});
  const auto other = graph(2, {{0, 1, 1, 1}});
  FlowLedger fg(g);
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), global(std::nullopt, dag_of(*other, {0}))); }),
            "support belongs to another graph");
}

// ---------------------------------------------------------------- W: weight bases

// W1, W12: E-configured; CAPACITY weights give (6, 6); the lossless bound 15.
TEST(WeightBases, Explicit) {
  const auto g = par2(10, 10);
  EXPECT_EQ(*make(g, {0}, {1}, weighted(LOSSY, WeightBasis::Explicit, {1, 2}))->place(12, kEdges).edge_flows,
            (Amounts {4, 8}));
  EXPECT_EQ(*make(g, {0}, {1}, cost_only(W, LOSSY))->place(12, kEdges).edge_flows, (Amounts {6, 6}));
  const auto ll = weighted(LOSSLESS, WeightBasis::Explicit, {1, 2});
  EXPECT_EQ(make(g, {0}, {1}, ll)->place(INF).placed, 15);
  EXPECT_EQ(capacity(*g, v({0}), v({1}), query(ll)).placed, 15);
}

// S3-01 (review): legal weights below 2^62 each may sum to 2^64 or more; the sum is exact. Eight
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
      const auto r = t->replace(8, kEdges);
      EXPECT_EQ(r.placed, 8) << w;
      EXPECT_EQ(*r.edge_flows, Amounts(8, 1)) << w;
      EXPECT_TRUE(r.source_loss.is_zero()) << w;
      const auto x = make(g, {0}, {1}, weighted(admission, WeightBasis::Explicit, Amounts(8, w)))->replace(8, kEdges);
      EXPECT_EQ(x.placed, 8) << w;
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

// W2: E-snapshot over the baseline (2, 10).
TEST(WeightBases, Snapshot) {
  const auto g = par2(10, 10);
  const Amounts base = {2, 10};
  auto ledger = [&] { return std::make_shared<FlowLedger>(g, base); };
  EXPECT_EQ(make(g, {0}, {1}, cost_only(W), nullptr, ledger())->place(12).placed, 4);
  EXPECT_EQ(make(g, {0}, {1}, weighted(LOSSLESS, WeightBasis::ResidualSnapshot), nullptr, ledger())->place(12).placed,
            12);
}

// W3: the priority no longer chooses weights (CAPACITY is the default).
TEST(WeightBases, PriorityIsNotAWeight) {
  const auto g = par2(10, 10);
  auto cfg = cost_only(W);
  cfg.selection = PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Residual);
  auto p = make(g, {0}, {1}, cfg, nullptr, std::make_shared<FlowLedger>(g, Amounts {2, 10}));
  EXPECT_EQ(p->place(12).placed, 4);
  EXPECT_EQ(p->config().weights, WeightBasis::Capacity);
}

// W4: a positive weight on an edge of capacity 0 becomes 0: all 7 on e1, no drop.
TEST(WeightBases, ExplicitOnZeroCapacity) {
  const auto g = par2(0, 10);
  const auto r = make(g, {0}, {1}, weighted(LOSSY, WeightBasis::Explicit, {5, 1}))->place(7);
  EXPECT_EQ(r.placed, 7);
  EXPECT_TRUE(r.call_drops.empty());
}

// W5: a zero weight prunes its edge.
TEST(WeightBases, ExplicitZeroPrunes) {
  const auto g = par2(10, 10);
  EXPECT_EQ(make(g, {0}, {1}, weighted(LOSSLESS, WeightBasis::Explicit, {0, 3}))->place(12).placed, 10);
  const auto r = make(g, {0}, {1}, weighted(LOSSY, WeightBasis::Explicit, {0, 3}))->place(12);
  EXPECT_EQ(r.placed, 10);
  EXPECT_EQ(r.call_drops, (Drops {{1, 2}}));
}

// W6: every weight 0: NO_ACTIVE_WEIGHT; Lossy loses the offer at the source.
TEST(WeightBases, AllZeroIsNoActiveWeight) {
  const auto g = par2(10, 10);
  for (const auto& [adm, loss] : {std::pair {LOSSY, Amount {5}}, std::pair {LOSSLESS, Amount {0}}}) {
    const auto r = make(g, {0}, {1}, weighted(adm, WeightBasis::Explicit, {0, 0}))->place(5);
    EXPECT_EQ(r.termination, Termination::NoActiveWeight);
    EXPECT_EQ(r.placed, 0);
    EXPECT_EQ(r.source_loss, c(loss));
  }
}

// W7: Weighted with the Headroom priority is accepted (CostOnly): the single path is e1.
TEST(WeightBases, WeightedWithHeadroom) {
  const auto g = par2(10, 10);
  auto cfg = cost_only(W);
  cfg.selection = PathSelection(EdgeTieBreak::PreferHigherPriority, EdgePriority::Headroom, false);
  auto p = make(g, {0}, {1}, cfg, nullptr, std::make_shared<FlowLedger>(g, Amounts {2, 10}));
  EXPECT_EQ(p->place(12).placed, 10);
}

// W8: the snapshot is frozen for the epoch: call 2 splits 6 as H(6, (2, 10)) = (1, 5), both
// dropped.
TEST(WeightBases, SnapshotFrozenForTheEpoch) {
  const auto g = par2(10, 10);
  auto p = make(g, {0}, {1}, weighted(LOSSY, WeightBasis::ResidualSnapshot), nullptr,
                std::make_shared<FlowLedger>(g, Amounts {2, 10}));
  EXPECT_EQ(*p->place(12, kEdges).edge_flows, (Amounts {2, 10}));
  const auto r = p->place(6);
  EXPECT_EQ(r.placed, 0);
  EXPECT_EQ(r.call_drops, (Drops {{0, 1}, {1, 5}}));
  EXPECT_EQ(r.call_source_loss, 0);
}

// W9, W10: the texts of the weight fields (the bool entry is the binding's conversion).
TEST(WeightBases, MatrixTexts) {
  const auto g = par2(10, 10);
  FlowLedger fg(g);
  auto with = [](FlowPolicyConfig cfg, std::optional<WeightBasis> w, std::optional<std::vector<Amount>> ew) {
    cfg.weights = w;
    if (ew.has_value()) cfg.edge_weights = std::make_shared<const std::vector<Amount>>(*ew);
    return cfg;
  };
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}});
  const std::vector<std::pair<FlowPolicyConfig, std::string>> cases = {
      {with(cost_only(), WeightBasis::Capacity, std::nullopt), "weights require CostOnly placement=WEIGHTED"},
      {with(cost_only(), std::nullopt, Amounts {1, 1}), "edge_weights require weights=EXPLICIT"},
      {with(cost_only(W), WeightBasis::Explicit, std::nullopt), "weights=EXPLICIT requires edge_weights"},
      {with(cost_only(W), WeightBasis::Capacity, Amounts {1, 1}), "edge_weights require weights=EXPLICIT"},
      {with(global(), WeightBasis::Capacity, std::nullopt), "Global does not accept weights"},
      {with(pinned_config(prog), WeightBasis::Capacity, std::nullopt), "Pinned does not accept weights"},
      {with(capacity_aware(), WeightBasis::Capacity, std::nullopt), "T1 does not accept weights"},
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

// W11: one configuration's edge_weights array is shared by every policy built from it (R11).
TEST(WeightBases, SharedArray) {
  const auto g = par2(10, 10);
  const auto cfg = weighted(LOSSLESS, WeightBasis::Explicit, {1, 2});
  FlowLedger fg(g);
  FlowPolicy a(fg, v({0}), v({1}), cfg);
  FlowPolicy b(fg, v({0}), v({1}), cfg);
  const auto* buffer = cfg.edge_weights.get();
  EXPECT_EQ(detail::PolicyAccess::engine(a).config.edge_weights.get(), buffer);
  EXPECT_EQ(detail::PolicyAccess::engine(b).config.edge_weights.get(), buffer);
  EXPECT_EQ(a.config().edge_weights.get(), buffer);
  // Control: an equal array built separately is another buffer.
  FlowPolicy c2(fg, v({0}), v({1}), weighted(LOSSLESS, WeightBasis::Explicit, {1, 2}));
  EXPECT_NE(detail::PolicyAccess::engine(c2).config.edge_weights.get(), buffer);
}

// ---------------------------------------------------------------- F: from_routes

// F1 to F22: the validation texts in order (graph PAR2(4, 8), 0 -> 1 unless stated).
TEST(FromRoutes, Texts) {
  const auto g = par2(4, 8);
  const auto p0 = tunnel_route(*g, 0, {0});
  const Route d01 = dag_route(*g, 0, {0, 1});
  auto err = [&](const std::vector<NodeId>& s, const std::vector<NodeId>& d, const std::vector<Route>& routes,
                 const RouteOptions& o = {}, const Graph* on = nullptr) {
    return error_of([&] { (void)ForwardingProgram::from_routes(on ? *on : *g, s, d, routes, o); });
  };
  EXPECT_EQ(err({}, {1}, {p0}), "sources: the node set is empty");
  EXPECT_EQ(err({0}, {0}, {p0}), "node 0 is in both endpoint sets: sources and destinations must not overlap");
  EXPECT_EQ(err({0}, {1}, {}), "routes must not be empty");
  EXPECT_EQ(err({0}, {1}, {p0}, RouteOptions {.placement = EB}), "placement applies only to PredDAG routes");
  EXPECT_EQ(err({0}, {1}, {p0}, RouteOptions {.weights = WeightBasis::Capacity}),
            "weights applies only to PredDAG routes");
  EXPECT_EQ(err({0}, {1}, {p0}, RouteOptions {.edge_weights = Amounts {1, 1}}),
            "edge_weights applies only to PredDAG routes");
  EXPECT_EQ(err({0}, {1}, {p0}, RouteOptions {.residual = Amounts {1, 1}}), "residual applies only to PredDAG routes");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.placement = MF}),
            "from_routes placement must be EQUAL_BALANCED or WEIGHTED");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.weights = WeightBasis::Capacity}), "weights require placement=WEIGHTED");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.placement = W, .weights = WeightBasis::Explicit}),
            "weights=EXPLICIT requires edge_weights");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.placement = W, .edge_weights = Amounts {1, 1}}),
            "edge_weights require weights=EXPLICIT");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.placement = W, .weights = WeightBasis::ResidualSnapshot}),
            "weights=RESIDUAL_SNAPSHOT requires residual");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.placement = W, .residual = Amounts {4, 8}}),
            "residual requires weights=RESIDUAL_SNAPSHOT");
  EXPECT_EQ(err({0}, {1}, {d01}, RouteOptions {.placement = W, .weights = WeightBasis::Explicit, .edge_weights = Amounts {1}}),
            "edge_weights: expected 2 entries");
  EXPECT_EQ(err({0}, {1}, {d01},
                RouteOptions {.placement = W, .weights = WeightBasis::Explicit, .edge_weights = Amounts {1, kL62}}),
            "edge_weights[1] must be an integer in [0, 2^62)");
  EXPECT_EQ(err({0}, {1}, {d01},
                RouteOptions {.placement = W, .weights = WeightBasis::ResidualSnapshot, .residual = Amounts {5, 8}}),
            "residual[0] = 5 is outside [0, capacity_units[0] = 4]");
  const auto other = par2(4, 8);
  EXPECT_EQ(err({0}, {1}, {tunnel_route(*other, 0, {0})}), "routes[0]: bundle belongs to another graph");
  EXPECT_EQ(err({0}, {1}, {tunnel_route(*g, 1, {0})}), "routes[0]: head 1 is not a source");
  const Route hop2 {0, std::make_shared<const Path>(*g, std::vector<std::vector<EdgeId>> {{0, 1}}), 1};
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

// F23: the routes in list order, each route's checks in order.
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

// F24: a zero-weight route keeps its index and is validated in full.
TEST(FromRoutes, ZeroWeightKept) {
  const auto g = par2(4, 8);
  const auto prog = program(*g, {0}, {1}, {tunnel_route(*g, 0, {0}, 0), tunnel_route(*g, 0, {1})});
  std::vector<std::pair<std::int32_t, Amount>> got;
  for (const auto& t : prog->tunnels()) got.emplace_back(t.index, t.weight);
  EXPECT_EQ(got, (std::vector<std::pair<std::int32_t, Amount>> {{0, 0}, {1, 1}}));
  EXPECT_EQ(error_of([&] { (void)program(*g, {0}, {1}, {tunnel_route(*g, 1, {0}, 0)}); }),
            "routes[0]: head 1 is not a source");
}

// F25: tunnels listed with indices 0 and 2 (the split route takes 1); local_index 0 and 2 (R2).
TEST(FromRoutes, MixedIndices) {
  const auto g = par2(4, 8);
  const auto prog =
      program(*g, {0}, {1}, {tunnel_route(*g, 0, {1}), dag_route(*g, 0, {0, 1}), tunnel_route(*g, 0, {0})});
  std::vector<std::pair<std::int32_t, std::int32_t>> got;
  for (const auto& t : prog->tunnels()) got.emplace_back(t.index, t.local_index);
  EXPECT_EQ(got, (std::vector<std::pair<std::int32_t, std::int32_t>> {{0, 0}, {2, 2}}));
}

// F26: the endpoints are normalized like a policy's (sorted, duplicates removed).
TEST(FromRoutes, Endpoints) {
  const auto d = graph(4, kDiam);
  const auto prog = program(*d, {1, 0, 1}, {3, 3}, {tunnel_route(*d, 0, {0, 2})});
  EXPECT_EQ(vec(prog->sources()), v({0, 1}));
  EXPECT_EQ(vec(prog->destinations()), v({3}));
}

// ---------------------------------------------------------------- C: Pinned's configuration

// C1: the fields Pinned does not accept, in column order.
TEST(PinnedProgram, FieldTexts) {
  const auto g = par2(4, 8);
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}});
  FlowLedger fg(g);
  auto base = [&] { return pinned_config(prog); };
  auto placement = base();
  placement.placement = EB;
  auto origin = base();
  origin.origin = Origin::Equal;
  auto selection = base();
  selection.selection = PathSelection();
  auto weights = base();
  weights.weights = WeightBasis::Capacity;
  auto edge_weights = base();
  edge_weights.edge_weights = std::make_shared<const std::vector<Amount>>(Amounts {1, 1});
  auto all = placement;
  all.origin = Origin::Equal;
  all.selection = PathSelection();
  const std::vector<std::pair<FlowPolicyConfig, std::string>> cases = {
      {placement, "Pinned does not accept placement"},   {origin, "Pinned does not accept origin"},
      {selection, "Pinned does not accept selection"},   {weights, "Pinned does not accept weights"},
      {edge_weights, "Pinned does not accept edge_weights"}, {all, "Pinned does not accept placement"},
  };
  for (const auto& [cfg, text] : cases) EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), cfg); }), text);
}

// C2: the program's endpoints must be the policy's.
TEST(PinnedProgram, Endpoints) {
  const auto g = graph(3, {{0, 1, 1, 1}, {0, 2, 1, 1}});
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}});
  FlowLedger fg(g);
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1, 2}), pinned_config(prog)); }),
            "program endpoints do not match policy endpoints");
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1, 1}), pinned_config(prog)); }), "");
}

// C3, C4: the ownership texts, before the configuration matrix (ERR-ownership-texts).
TEST(PinnedProgram, Ownership) {
  const auto g = par2(4, 8);
  const auto other = par2(4, 8);
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}});
  const auto foreign = tunnels(*other, {0}, {1}, {{0, {0}, 1}});
  FlowLedger fg(g);
  auto with_placement = pinned_config(foreign);
  with_placement.placement = EB;
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), pinned_config(foreign)); }),
            "program belongs to another graph");
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), with_placement); }), "program belongs to another graph");
  EXPECT_EQ(error_of([&] {
              FlowPolicy p(fg, v({0}), v({1}), pinned_config(prog, LOSSLESS, no_edges(*other, {0})));
            }),
            "black_holes belongs to another graph");
  const auto mask = other->eligibility({}, std::vector<EdgeId> {0});
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), pinned_config(prog), &mask); }),
            "foreign Eligibility: built for another graph");
  DagCache cache(*other);
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), pinned_config(prog), nullptr, &cache); }),
            "foreign DagCache: bound to another graph");
  // The mask is checked before the program.
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), pinned_config(foreign), &mask); }),
            "foreign Eligibility: built for another graph");
}

// C5: the endpoint sets come first.
TEST(PinnedProgram, Order) {
  const auto g = par2(4, 8);
  const auto foreign = tunnels(*par2(4, 8), {0}, {1}, {{0, {0}, 1}});
  FlowLedger fg(g);
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, std::vector<NodeId> {}, v({1}), pinned_config(foreign)); }),
            "sources: the node set is empty");
}

// C8: several sources and destinations.
TEST(PinnedProgram, SeveralEndpoints) {
  const auto g = graph(4, kRel);
  const auto prog = tunnels(*g, {0, 1}, {3}, {{0, {0, 2}, 1}, {1, {1, 2}, 1}});
  EXPECT_EQ(pinned(g, prog)->place(INF).placed, 10);
  const auto two = graph(4, {{0, 2, 3, 1}, {1, 3, 4, 1}});
  const auto many = tunnels(*two, {0, 1}, {2, 3}, {{0, {0}, 1}, {1, {1}, 1}});
  EXPECT_EQ(pinned(two, many)->place(INF).placed, 6);  // RR(6, 2) = (3, 3) fits, RR(7, 2) = (4, 3) does not
}

// ---------------------------------------------------------------- L: lossless bounds

namespace {
Amount bound(const std::shared_ptr<const Graph>& g, const std::vector<Spec>& specs, const std::vector<NodeId>& s = {0},
             const std::vector<NodeId>& d = {1}) {
  return pinned(g, tunnels(*g, s, d, specs))->place(INF).placed;
}
} // namespace

// L1 to L6: E-duplicates, E-deduplicated, E-static-path-fill, E-zero-first, E-zero-last,
// E-support-pinned-bound.
TEST(PinnedProgram, Duplicates) {
  EXPECT_EQ(bound(par2(4, 8), {{0, {0}, 1}, {0, {1}, 1}, {0, {1}, 1}}), 12);
  EXPECT_EQ(bound(par2(4, 8), {{0, {0}, 1}, {0, {1}, 1}}), 8);
  EXPECT_EQ(bound(par2(6, 4), {{0, {0}, 1}, {0, {1}, 1}}), 9);
  EXPECT_EQ(bound(par2(0, 4), {{0, {0}, 1}, {0, {1}, 1}}), 0);
  EXPECT_EQ(bound(par2(4, 0), {{0, {0}, 1}, {0, {1}, 1}}), 1);
  EXPECT_EQ(bound(graph(4, kPinloss), {{0, {0, 2}, 1}, {0, {1, 3}, 1}}, {0}, {3}), 4);
  // Records per entry: L1's three records, entries 0, 1, 2 (D2 below compares their rows).
  const auto g = par2(4, 8);
  auto p = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 1}, {0, {1}, 1}, {0, {1}, 1}}));
  (void)p->place(INF);
  EXPECT_EQ(delivered(*p.p), (Amounts {4, 4, 4}));
}

// L7: H(15, (1, 2)) = (5, 10) fits, H(16) = (5, 11) does not.
TEST(PinnedProgram, WeightedTunnels) { EXPECT_EQ(bound(par2(10, 10), {{0, {0}, 1}, {0, {1}, 2}}), 15); }

// L8, D3: a weight-0 tunnel carries nothing; its record and row exist with zeros.
TEST(PinnedProgram, ZeroWeight) {
  const auto g = par2(4, 8);
  auto p = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 0}, {0, {1}, 1}}));
  const auto r = p->place(INF);
  EXPECT_EQ(r.placed, 8);
  EXPECT_EQ(entries(*p.p), (std::vector<std::optional<std::int32_t>> {0, 1}));
  EXPECT_EQ(delivered(*p.p), (Amounts {0, 8}));
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 0, 0, 0}, {1, 0, UP, 8, 8, 0}}));
  // Created at preparation: a zero-delivery record before any positive placement of it.
  auto q = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 0}, {0, {1}, 1}}));
  (void)q->place(1);
  EXPECT_EQ(delivered(*q.p), (Amounts {0, 1}));
}

// L9 to L11: INC-lossless, INC-lossy, INC-lossless-refusal.
TEST(PinnedProgram, Increments) {
  const auto g = graph(3, kInc);
  const auto prog = tunnels(*g, {0, 1}, {2}, {{0, {0}, 1}, {1, {1}, 1}});
  auto p = pinned(g, prog);
  EXPECT_EQ(p->place(3).placed, 2);
  EXPECT_EQ(p->place(1).placed, 0);
  auto q = pinned(g, prog, LOSSY);
  const auto r1 = q->place(3);
  EXPECT_EQ(r1.placed, 2);
  EXPECT_EQ(r1.call_drops, (Drops {{0, 1}}));
  const auto r2 = q->place(1);  // RR(4) - RR(3) = (0, 1)
  EXPECT_EQ(r2.placed, 1);
  EXPECT_TRUE(r2.call_drops.empty());
  EXPECT_EQ(lo(r2.by_source[0].offered), 2);
  EXPECT_EQ(lo(r2.by_source[1].offered), 2);
  const auto one = graph(2, kOne);
  auto ledger = std::make_shared<FlowLedger>(one);
  auto other = make(one, {0}, {1}, cost_only(), nullptr, ledger);
  (void)other->place(8);
  auto t = pinned(one, tunnels(*one, {0}, {1}, {{0, {0}, 1}}), LOSSLESS, nullptr, nullptr, ledger);
  EXPECT_EQ(t->place(5).placed, 2);
  other->remove();
  (void)other->place(5);
  EXPECT_EQ(t->place(1).placed, 1);
}

// L12, L13: +inf is the bound, EXHAUSTED; Lossy rejects it.
TEST(PinnedProgram, Infinite) {
  const auto g = par2(6, 4);
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}, {0, {1}, 1}});
  const auto r = pinned(g, prog)->place(INF);
  EXPECT_EQ(r.placed, 9);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  EXPECT_FALSE(r.remaining.has_value());
  EXPECT_FALSE(r.requested.has_value());
  auto lossy = pinned(g, prog, LOSSY);
  EXPECT_EQ(error_of([&] { (void)lossy->place(INF); }), kFinite);
}

// L14: a tunnel and a split route side by side: RR(8, 2) = (4, 4).
TEST(PinnedProgram, MixedEntries) {
  const auto g = par2(4, 8);
  auto p = pinned(g, program(*g, {0}, {1}, {tunnel_route(*g, 0, {0}), dag_route(*g, 0, {1})}));
  const auto r = p->place(INF);
  EXPECT_EQ(r.placed, 8);
  EXPECT_EQ(entries(*p.p), (std::vector<std::optional<std::int32_t>> {0, 1}));
  EXPECT_EQ(delivered(*p.p), (Amounts {4, 4}));
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 4, 4, 0}}));
}

// L15: E-split-contexts (each identical route splits its own unit) and E-split-contexts-pooled.
TEST(PinnedProgram, SplitContextsAreSeparate) {
  const auto g = par2(1, 1);
  const auto d = dag_of(*g, {0, 1});
  const auto prog = program(*g, {0}, {1}, {Route {0, d, 1}, Route {0, d, 1}});
  EXPECT_EQ(pinned(g, prog)->place(INF).placed, 1);
  const auto r = pinned(g, prog, LOSSY)->replace(2, kEdges);
  EXPECT_EQ(r.call_drops, (Drops {{0, 1}}));
  EXPECT_EQ(*r.edge_flows, (Amounts {1, 0}));
  const auto pooled = make(g, {0}, {1}, cost_only())->place(INF, kEdges);
  EXPECT_EQ(pooled.placed, 2);
  EXPECT_EQ(*pooled.edge_flows, (Amounts {1, 1}));
}

// ---------------------------------------------------------------- Y: losses inside one program

namespace {
struct LossyRun {
  PlacementResult r;
  std::vector<std::tuple<Amount, Amount, Amount>> got;  // per tunnel (assigned, delivered, dropped)
};
LossyRun lossy_rows(const std::shared_ptr<const Graph>& g, const std::vector<NodeId>& s, const std::vector<NodeId>& d,
                    const std::vector<Spec>& specs, Amount offer) {
  LossyRun out {pinned(g, tunnels(*g, s, d, specs), LOSSY)->replace(offer), {}};
  for (const auto& t : out.r.by_tunnel) out.got.emplace_back(lo(t.assigned), t.delivered, lo(t.dropped));
  return out;
}
using Triples = std::vector<std::tuple<Amount, Amount, Amount>>;
Amounts delivered_of(const LossyRun& x) {
  Amounts out;
  for (const auto& [a, d, l] : x.got) out.push_back(d);
  return out;
}
} // namespace

// Y1 to Y8: D-merge, D-key-tie, D-rounding, D-hierarchy, D-zero-service, D-subdivision-before
// and -after, E-support-pinned-loss (Z1's rows).
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
  EXPECT_EQ(zero.r.call_drops, (Drops {{0, 2}}));
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
  EXPECT_EQ(pin.r.call_drops, (Drops {{2, 4}}));
  EXPECT_EQ(pin.r.placed, 8);
}

// Y9: no link can drop (U0 = 2 on both barriers' candidates): no rounds.
TEST(ProgramCycles, NoBarrierNoRounds) {
  const auto run = lossy_rows(cyc(2, 2), {4}, {5}, {{4, kCycT0, 1}, {4, kCycT1, 1}}, 2);
  EXPECT_EQ(delivered_of(run), (Amounts {1, 1}));
  EXPECT_TRUE(run.r.call_drops.empty());
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().groups, 0);
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 0);
#endif
}

// Y10: X and Y wait on each other; the rounds from zero settle at round 3.
TEST(ProgramCycles, Converges) {
  const auto run = lossy_rows(cyc(1, 1), {4}, {5}, {{4, kCycT0, 1}, {4, kCycT1, 1}}, 2);
  EXPECT_EQ(run.got, (Triples {{1, 1, 0}, {1, 0, 1}}));
  EXPECT_EQ(run.r.call_drops, (Drops {{3, 1}}));
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().groups, 1);
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 3);
#endif
}

// Y11: amounts (1, 2): the rounds alternate between two states; nothing is published.
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

// Y12: three barriers that never settle: the round limit 64 * 3.
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

// Y13, Y14: Y split in two in series changes nothing (serial normalization).
TEST(ProgramCycles, SubdivisionInvariant) {
  const auto run = lossy_rows(cycs(1, 1, 1), {4}, {5}, {{4, kCycsT0, 1}, {4, kCycsT1, 1}}, 2);
  EXPECT_EQ(run.got, (Triples {{1, 1, 0}, {1, 0, 1}}));
  EXPECT_EQ(run.r.call_drops, (Drops {{3, 1}}));
  const auto g = cycs(1, 2, 2);
  auto p = pinned(g, tunnels(*g, {4}, {5}, {{4, kCycsT0, 1}, {4, kCycsT1, 2}}), LOSSY);
  EXPECT_EQ(error_of([&] { (void)p->replace(3); }), kPeriodTwo);
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
// Y15 through the hook: an iterated result corrupted after its rounds fails the independent check
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
  EXPECT_EQ(p->replace(2).placed, 1);
}
#endif

// ---------------------------------------------------------------- R: replay under a mask

namespace {
// G1U's optimized program: [(0, [e0, e2], 2), (0, [e1, e3], 10)].
std::shared_ptr<const ForwardingProgram> g1u_program(const std::shared_ptr<const Graph>& g) {
  auto solver = make(g, {0}, {3}, global());
  const auto prog = solver->replace(12).program;
  EXPECT_EQ(listing(*prog), (Listing {{0, {0, 2}, 2}, {0, {1, 3}, 10}}));
  return prog;
}
} // namespace

// R1, the 8.5 call with a Global program: tunnel 0 FAILED, the head's 12 goes to tunnel 1.
TEST(Replay, GlobalProgram) {
  const auto g = graph(4, kG1U);
  const auto prog = g1u_program(g);
  const auto r = pinned(g, prog, LOSSY, no_edges(*g, {2}))->replace(12);
  EXPECT_EQ(r.placed, 10);
  EXPECT_EQ(rows(r), (Rows {{0, 0, FAILED, 0, 0, 0}, {1, 0, UP, 12, 10, 2}}));
  EXPECT_EQ(r.call_drops, (Drops {{1, 2}}));
  EXPECT_EQ(pinned(g, prog, LOSSLESS, no_edges(*g, {2}))->place(INF).placed, 10);
}

// R2: E-source-survivor: H(16, (3, 1)) = (12, 4).
TEST(Replay, SourceSurvivor) {
  const auto g = graph(4, kSurv);
  auto solver = make(g, {0, 1, 2}, {3}, global());
  const auto prog = solver->place(9).program;
  EXPECT_EQ(*prog->source_weights(), (std::vector<std::pair<NodeId, Amount>> {{0, 3}, {1, 1}, {2, 5}}));
  const auto r = pinned(g, prog, LOSSY, no_nodes(*g, {2}))->replace(16);
  EXPECT_EQ(r.placed, 4);
  EXPECT_EQ(lo(r.by_source[0].offered), 12);
  EXPECT_EQ(lo(r.by_source[1].offered), 4);
  EXPECT_EQ(lo(r.by_source[2].offered), 0);
  EXPECT_EQ(r.call_drops, (Drops {{0, 9}, {1, 3}}));
  EXPECT_EQ(std::get<2>(rows(r)[2]), FAILED);
  EXPECT_EQ(pinned(g, prog, LOSSLESS, no_nodes(*g, {2}))->place(INF).placed, 4);
}

// R3: the topology basis keeps head 1 (it still reaches 2 over e2) and it loses its share.
TEST(Replay, TopologyBasisKeepsTheHead) {
  const auto g = graph(3, kRB);
  auto solver = make(g, {0, 1}, {2}, global());
  const auto prog = solver->place(10).program;
  EXPECT_EQ(listing(*prog), (Listing {{0, {0}, 5}, {1, {1}, 5}}));
  const auto r = pinned(g, prog, LOSSY, no_edges(*g, {1}))->replace(10);
  EXPECT_EQ(r.placed, 5);
  EXPECT_EQ(r.source_loss, c(5));
  EXPECT_EQ(pinned(g, prog, LOSSLESS, no_edges(*g, {1}))->place(INF).placed, 1);
}

// R4: the supplied basis drops head 1: head 0 gets all 10.
TEST(Replay, SuppliedBasisDropsTheHead) {
  const auto g = graph(3, kRB);
  const auto prog = tunnels(*g, {0, 1}, {2}, {{0, {0}, 1}, {1, {1}, 1}});
  const auto r = pinned(g, prog, LOSSY, no_edges(*g, {1}))->replace(10);
  EXPECT_EQ(r.placed, 5);
  EXPECT_EQ(r.call_drops, (Drops {{0, 5}}));
  EXPECT_EQ(r.source_loss, c(0));
  EXPECT_EQ(pinned(g, prog, LOSSLESS, no_edges(*g, {1}))->place(INF).placed, 5);
}

// R5: a link of free capacity 0 is not a failure.
TEST(Replay, ZeroCapacityIsNotFailed) {
  const auto g = graph(4, kG1U);
  const auto prog = g1u_program(g);
  const auto r = pinned(g, prog, LOSSY, nullptr, nullptr, std::make_shared<FlowLedger>(g, Amounts {2, 0, 2, 10}))
                     ->replace(12);
  EXPECT_EQ(r.placed, 2);
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 2, 2, 0}, {1, 0, UP, 10, 0, 10}}));
  EXPECT_EQ(r.call_drops, (Drops {{1, 10}}));
}

// R6: the remaining head's captured weight is 0: the whole offer is lost at the sources.
TEST(Replay, AllRemainingWeightsZero) {
  const auto g = graph(3, {{0, 2, 5, 1}, {1, 2, 5, 1}});
  auto solver = make(g, {0, 1}, {2}, global());
  const auto prog = solver->place(1).program;
  const auto r = pinned(g, prog, LOSSY, no_nodes(*g, {0}))->replace(4);
  EXPECT_EQ(r.termination, Termination::NoActiveWeight);
  EXPECT_EQ(r.source_loss, c(4));
  EXPECT_EQ(r.placed, 0);
  EXPECT_EQ(rows(r), (Rows {{0, 0, FAILED, 0, 0, 0}}));
  const auto ll = pinned(g, prog, LOSSLESS, no_nodes(*g, {0}))->replace(4);
  EXPECT_EQ(ll.termination, Termination::NoActiveWeight);
  EXPECT_EQ(ll.placed, 0);
}

// R7: STALE-remote-tail-repair: the branch 0 -> 2 is pruned.
TEST(Replay, SplitRoutePruned) {
  const auto g = graph(4, kDiam);
  auto producer = make(g, {0}, {3}, cost_only(EB, LOSSY));
  const auto prog = producer->place(1).program;
  const auto r = pinned(g, prog, LOSSY, no_edges(*g, {3}))->replace(100, kEdges);
  EXPECT_EQ(r.placed, 100);
  EXPECT_EQ(*r.edge_flows, (Amounts {100, 0, 100, 0}));
}

// R8: STALE-tunnel-node-repair.
TEST(Replay, TunnelFailedByNode) {
  const auto g = graph(4, kTun);
  const auto prog = tunnels(*g, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}});
  const auto r = pinned(g, prog, LOSSY, no_nodes(*g, {2}))->replace(10);
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 10, 10, 0}, {1, 0, FAILED, 0, 0, 0}}));
  EXPECT_EQ(r.placed, 10);
  EXPECT_EQ(pinned(g, prog, LOSSLESS, no_nodes(*g, {2}))->place(INF).placed, 10);
}

// R9: weights frozen at preparation: a member pruned there stays pruned in every replay.
TEST(Replay, FrozenWeights) {
  const auto z = par2(0, 10);
  auto producer = make(z, {0}, {1}, cost_only(W, LOSSY));
  const auto wp = producer->place(1).program;
  const auto r = pinned(z, wp, LOSSY)->replace(5);
  EXPECT_EQ(r.placed, 5);
  EXPECT_TRUE(r.call_drops.empty());
  // A snapshot over the baseline (0, 10) prunes e0; the replay on a fresh ledger at full capacity
  // still sends everything to e1: 10 delivered, 5 dropped there.
  const auto g = par2(10, 10);
  auto snap = make(g, {0}, {1}, weighted(LOSSY, WeightBasis::ResidualSnapshot), nullptr,
                   std::make_shared<FlowLedger>(g, Amounts {0, 10}));
  const auto sp = snap->place(1).program;
  const auto rr = pinned(g, sp, LOSSY)->replace(15, kEdges);
  EXPECT_EQ(rr.placed, 10);
  EXPECT_EQ(*rr.edge_flows, (Amounts {0, 10}));
  EXPECT_EQ(rr.call_drops, (Drops {{1, 5}}));
}

// ---------------------------------------------------------------- S: stale forwarding

// S1 to S11: STALE-remote-tail, -wcmp-member, -wcmp-zero-capacity, -relief-original, -relief,
// -tunnel-node, -pooled-original, -pooled-remainder, -all-zero-weight, -black-source,
// -black-source-tunnel, -mask-and-black-hole.
TEST(StaleForwarding, RemoteTail) {
  {
    const auto g = graph(4, kDiam);
    auto p = make(g, {0}, {3}, cost_only(EB, LOSSY));
    (void)p->replace(100);
    const auto installed = p->program();
    const auto bh = g->eligibility({}, std::vector<EdgeId> {3});
    const auto r = p->forward(100, &bh);
    EXPECT_EQ(r.placed, 50);
    EXPECT_EQ(r.call_drops, (Drops {{3, 50}}));
    EXPECT_EQ(p->program(), installed);
    EXPECT_FALSE(r.routes_changed || r.weights_changed);
    EXPECT_EQ(p->replace(100).placed, 100);
  }
  {
    const auto w = graph(2, {{0, 1, 10, 1}, {0, 1, 40, 1}});
    auto q = make(w, {0}, {1}, cost_only(W, LOSSY));
    (void)q->replace(40);
    const auto bh = w->eligibility({}, std::vector<EdgeId> {0});
    const auto r = q->forward(40, &bh);
    EXPECT_EQ(r.placed, 32);
    EXPECT_EQ(r.call_drops, (Drops {{0, 8}}));
  }
  {
    const auto wz = graph(2, {{0, 1, 0, 1}, {0, 1, 40, 1}});
    auto q = make(wz, {0}, {1}, cost_only(W, LOSSY));
    (void)q->replace(1);
    EXPECT_EQ(q->forward(40).placed, 40);
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
    EXPECT_EQ(r.call_drops, (Drops {{1, 5}}));
    EXPECT_EQ(pinned(t, tp, LOSSLESS, nullptr, no_nodes(*t, {2}))->replace(10).placed, 1);
  }
  {
    const auto pool = graph(4, kPool2);
    auto p = make(pool, {0, 1}, {3}, cost_only(EB, LOSSY));
    (void)p->replace(0);
    const auto r = p->forward(2, nullptr, kEdges);
    EXPECT_EQ(r.placed, 2);
    EXPECT_EQ(*r.edge_flows, (Amounts {1, 1, 1, 1}));
    const auto bh = pool->eligibility({}, std::vector<EdgeId> {0, 3});
    const auto s8 = p->forward(2, &bh);
    EXPECT_EQ(s8.placed, 1);
    EXPECT_EQ(s8.call_drops, (Drops {{0, 1}}));
  }
  {
    const auto zero = graph(2, {{0, 1, 0, 1}});
    auto q = make(zero, {0}, {1}, cost_only(W, LOSSY));
    (void)q->replace(0);
    const auto r = q->forward(2);
    EXPECT_EQ(r.source_loss, c(2));
    EXPECT_EQ(r.termination, Termination::NoActiveWeight);
  }
  {
    const auto one = graph(2, kOne);
    auto q = make(one, {0}, {1}, cost_only(EB, LOSSY));
    (void)q->replace(0);
    const auto bh = one->eligibility(std::vector<NodeId> {0}, {});
    const auto r = q->forward(10, &bh);
    EXPECT_EQ(r.source_loss, c(10));
    EXPECT_TRUE(r.call_drops.empty());
    EXPECT_EQ(r.placed, 0);
    const auto s11 = pinned(one, tunnels(*one, {0}, {1}, {{0, {0}, 1}}), LOSSY, nullptr, no_nodes(*one, {0}))->replace(10);
    EXPECT_EQ(rows(s11), (Rows {{0, 0, BLACKHOLED, 0, 0, 0}}));
    EXPECT_EQ(s11.source_loss, c(10));
  }
  {
    // S12: STALE-mask-and-black-hole: a masked and black-holed element is masked.
    const auto m = graph(5, kMB);
    const auto mp = tunnels(*m, {0}, {4}, {{0, {0, 3}, 1}, {0, {1, 4}, 1}, {0, {2, 5}, 1}});
    const auto r = pinned(m, mp, LOSSY, no_edges(*m, {3}), no_edges(*m, {3, 4}))->replace(10);
    EXPECT_EQ(rows(r), (Rows {{0, 0, FAILED, 0, 0, 0}, {1, 0, BLACKHOLED, 5, 0, 5}, {2, 0, UP, 5, 5, 0}}));
    EXPECT_EQ(r.call_drops, (Drops {{4, 5}}));
  }
}

// S13: R1's program with e2 black-holed instead of masked.
TEST(StaleForwarding, GlobalProgram) {
  const auto g = graph(4, kG1U);
  const auto prog = g1u_program(g);
  const auto r = pinned(g, prog, LOSSY, nullptr, no_edges(*g, {2}))->replace(12);
  EXPECT_EQ(rows(r), (Rows {{0, 0, BLACKHOLED, 2, 0, 2}, {1, 0, UP, 10, 10, 0}}));
  EXPECT_EQ(r.call_drops, (Drops {{2, 2}}));
  // H(4, (2, 10)) = (0, 4) fits; H(5, (2, 10)) = (1, 4) sends a unit into the black hole.
  EXPECT_EQ(pinned(g, prog, LOSSLESS, nullptr, no_edges(*g, {2}))->replace(12).placed, 4);
}

// S14: an all-allowed Eligibility is no black hole.
TEST(StaleForwarding, AllAllowedIsNone) {
  const auto g = graph(4, kG1U);
  const auto prog = g1u_program(g);
  const auto allowed = std::make_shared<const Eligibility>(g->eligibility());
  const auto with = pinned(g, prog, LOSSY, nullptr, allowed)->replace(12);
  const auto without = pinned(g, prog, LOSSY)->replace(12);
  EXPECT_EQ(rows(with), rows(without));
  EXPECT_EQ(with.placed, 12);
  EXPECT_EQ(with.call_drops, without.call_drops);
  // A forward with an all-allowed value equals a forward without.
  auto p = make(g, {0}, {3}, cost_only(EB, LOSSY));
  (void)p->replace(12);
  const auto all = g->eligibility();
  EXPECT_EQ(p->forward(12, &all).placed, p->forward(12).placed);
}

// S15: place, evaluate, replace, capacity and sensitivity all see a Pinned policy's black holes
// (S6: TUN, node 2 black-holed).
TEST(StaleForwarding, EveryOperation) {
  const auto t = graph(4, kTun);
  const auto tp = tunnels(*t, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}});
  const auto bh = no_nodes(*t, {2});
  EXPECT_EQ(pinned(t, tp, LOSSY, nullptr, bh)->place(10).placed, 5);
  EXPECT_EQ(pinned(t, tp, LOSSY, nullptr, bh)->evaluate(10).placed, 5);
  EXPECT_EQ(pinned(t, tp, LOSSY, nullptr, bh)->replace(10).placed, 5);
  EXPECT_EQ(capacity(*t, v({0}), v({3}), query(pinned_config(tp, LOSSY, bh), offer(10))).placed, 5);
  EXPECT_EQ(capacity(*t, v({0}), v({3}), query(pinned_config(tp, LOSSLESS, bh))).placed, 1);
  EXPECT_EQ(sensitivity(*t, v({0}), v({3}), query(pinned_config(tp, LOSSY, bh), offer(10)), 1).baseline.placed, 5);
  // Controls without the black holes.
  EXPECT_EQ(pinned(t, tp, LOSSY)->replace(10).placed, 10);
  EXPECT_EQ(capacity(*t, v({0}), v({3}), query(pinned_config(tp))).placed, 20);
}

// ---------------------------------------------------------------- O: operations

// O1: evaluate publishes nothing: ledger, records and epoch counters unchanged.
TEST(Operations, EvaluatePublishesNothing) {
  const auto g = par2(10, 10);
  auto p = make(g, {0}, {1}, cost_only());
  (void)p->replace(6);
  const auto ledger = ledger_state(*p.ledger);
  const auto before = epoch(*p.p);
  const auto r = p->evaluate(12);
  EXPECT_EQ(r.placed, 12);
  EXPECT_EQ(r.delivery_delta, 6);
  EXPECT_NE(r.program, nullptr);
  expect_same(ledger_state(*p.ledger), ledger, "evaluate");
  expect_unchanged(epoch(*p.p), before, "evaluate");
  EXPECT_EQ(p->place(0).requested, c(6));
}

// O2: replace retires the old ids and starts a new epoch.
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

// O3, O10: forward keeps the frozen weights; replace takes a new snapshot; an identical replace
// sets no flag.
TEST(Operations, ForwardKeepsTheProgram) {
  const auto g = par2(10, 10);
  auto ledger = std::make_shared<FlowLedger>(g, Amounts {2, 10});
  auto a = make(g, {0}, {1}, weighted(LOSSLESS, WeightBasis::ResidualSnapshot), nullptr, ledger);
  EXPECT_EQ(a->replace(6).placed, 6);
  auto b = pinned(g, tunnels(*g, {0}, {1}, {{0, {1}, 1}}), LOSSLESS, nullptr, nullptr, ledger);
  EXPECT_EQ(b->place(5).placed, 5);
  const auto installed = a->program();
  const auto f = a->forward(12);
  EXPECT_EQ(f.placed, 6);
  EXPECT_FALSE(f.routes_changed || f.weights_changed);
  EXPECT_EQ(a->program(), installed);
  const auto r = a->replace(12);
  EXPECT_EQ(r.placed, 7);
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
  EXPECT_EQ(r2.placed, r1.placed);
  EXPECT_FALSE(r2.routes_changed || r2.weights_changed);
  EXPECT_NE(ids(*a.p), ids1);
  EXPECT_EQ(delivered(*a.p), (Amounts {r1.placed}));
}

// O4, O5, S16: the forward texts.
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
  (void)ca->place(1);
  EXPECT_EQ(error_of([&] { (void)ca->forward(1); }), kNoProgram);
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}});
  auto pin = pinned(g, prog, LOSSY);
  EXPECT_EQ(error_of([&] { (void)pin->place(INF); }), kFinite);
  (void)pin->place(1);
  EXPECT_EQ(error_of([&] { (void)pin->forward(1); }), "Pinned forwards with replace");
}

// V1 to V4: the order of a call's checks.
TEST(Operations, ForwardOrder) {
  const auto g = par2(10, 10);
  const auto other = par2(1, 1);
  const auto foreign = other->eligibility({}, std::vector<EdgeId> {0});
  auto pin = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 1}}));
  EXPECT_EQ(error_of([&] { (void)pin->forward(INF, &foreign); }), "Pinned forwards with replace");
  auto p = make(g, {0}, {1}, cost_only(EB, LOSSY));
  EXPECT_EQ(error_of([&] { (void)p->forward(-1); }), kNoProgram);
  (void)p->place(1);
  EXPECT_EQ(error_of([&] { (void)p->forward(-1, &foreign); }), "offer must be an integer in [0, 2^62) or permitted +inf");
  EXPECT_EQ(error_of([&] { (void)p->forward(INF, &foreign); }), kFinite);
  EXPECT_EQ(error_of([&] { (void)p->forward(1, &foreign); }), "black_holes belongs to another graph");
  EXPECT_EQ(error_of([&] { (void)p->replace(-1); }), "offer must be an integer in [0, 2^62) or permitted +inf");
  EXPECT_EQ(error_of([&] { (void)p->replace(INF); }), kFinite);
}

// O6: Lossless replace(+inf): the bound, EXHAUSTED.
TEST(Operations, LosslessInfinite) {
  const auto r = make(par2(10, 10), {0}, {1}, cost_only())->replace(INF);
  EXPECT_EQ(r.placed, 20);
  EXPECT_FALSE(r.remaining.has_value());
  EXPECT_EQ(r.termination, Termination::Exhausted);
}

// O7: a zero replace installs a program of two split routes and two records at zero (R21).
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

// O8: a smaller replace: delivery_delta negative.
TEST(Operations, DeliveryDeltaNegative) {
  auto p = make(par2(10, 10), {0}, {1}, cost_only());
  (void)p->replace(10);
  const auto r = p->replace(4);
  EXPECT_EQ(r.placed, 4);
  EXPECT_EQ(r.delivery_delta, -6);
}

// O9: place continues the replace's epoch (its ids kept).
TEST(Operations, PlaceContinuesTheEpoch) {
  const auto g = par2(10, 10);
  auto p = make(g, {0}, {1}, cost_only());
  (void)p->replace(6);
  const auto before = ids(*p.p);
  const auto r = p->place(3);
  EXPECT_EQ(r.placed, 3);
  EXPECT_EQ(r.requested, c(9));
  EXPECT_EQ(r.admitted, c(9));
  EXPECT_EQ(ids(*p.p), before);
  auto q = make(g, {0}, {1}, cost_only(EB, LOSSY));
  (void)q->replace(6);
  const auto lr = q->place(3);
  EXPECT_EQ(lr.placed, 3);
  EXPECT_EQ(lr.offered, c(9));
  EXPECT_EQ(lr.requested, c(9));
}

// O11: FLAGS-controller-increment.
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

// O12: a re-solve that adds a path changes the routes.
TEST(Operations, ControllerRoutes) {
  const auto g = graph(4, kB1);
  auto q = make(g, {0, 1}, {3}, global());
  EXPECT_EQ(listing(*q->place(2).program), (Listing {{0, {0}, 2}}));
  const auto r = q->place(10);
  EXPECT_EQ(listing(*r.program), (Listing {{0, {0}, 2}, {1, {1, 2}, 10}}));
  EXPECT_TRUE(r.routes_changed);
}

// O13: FLAGS-fixed-offer: an offer-only change sets no flag.
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

// O14: forward over a controller's program reproduces the solve; above it, EXHAUSTED (R18).
TEST(Operations, ForwardControllerProgram) {
  const auto g = graph(4, kG1U);
  auto p = make(g, {0}, {3}, global());
  (void)p->replace(12);
  EXPECT_EQ(p->forward(12).placed, 12);
  const auto r = p->forward(20);
  EXPECT_EQ(r.placed, 12);
  EXPECT_EQ(r.termination, Termination::Exhausted);
}

// O15: remove() clears the program.
TEST(Operations, RemoveClearsTheProgram) {
  for (const auto& cfg : {global(), cost_only()}) {
    auto p = make(par2(10, 10), {0}, {1}, cfg);
    (void)p->place(3);
    EXPECT_NE(p->program(), nullptr);
    p->remove();
    EXPECT_EQ(p->program(), nullptr);
    const auto z = p->place(0);
    EXPECT_EQ(z.program, nullptr);
    EXPECT_EQ(z.placed, 0);
    EXPECT_EQ(z.delivered, 0);
    EXPECT_TRUE(z.by_tunnel.empty());
  }
}

// O16: today's CapacityAware: replace is place on a fresh epoch, evaluate publishes nothing.
TEST(Operations, CapacityAwareTransitional) {
  const auto g = graph(2, kOne);
  auto p = make(g, {0}, {1}, capacity_aware());
  (void)p->place(2);
  const auto old = ids(*p.p);
  const auto r = p->replace(5);
  EXPECT_EQ(r.placed, 5);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_NE(ids(*p.p), old);
  const auto state = ledger_state(*p.ledger);
  EXPECT_EQ(p->evaluate(3).placed, 3);
  expect_same(ledger_state(*p.ledger), state, "evaluate");
}

// O17: CostOnly's first positive call sets both flags (R21); the next none.
TEST(Operations, FirstProgramFlags) {
  auto q = make(graph(4, kPool2), {0, 1}, {3}, cost_only());
  const auto r = q->place(2);
  EXPECT_TRUE(r.routes_changed && r.weights_changed);
  const auto r2 = q->place(2);
  EXPECT_FALSE(r2.routes_changed || r2.weights_changed);
}

// O18: evaluate's outputs describe the evaluated state (R8).
TEST(Operations, EvaluateOutputs) {
  auto gl = make(graph(4, kB1), {0, 1}, {3}, global());
  (void)gl->place(2);
  const auto r = gl->evaluate(12, OutputOptions {.edge_flows = true, .residual = true, .min_cut = true});
  EXPECT_EQ(*r.edge_flows, (Amounts {2, 10, 10}));
  EXPECT_EQ(*r.residual, (Amounts {0, 0, 0}));
  EXPECT_EQ(*r.min_cut, (std::vector<EdgeId> {0, 1}));
  EXPECT_EQ(gl->placed_demand(), 2);
  EXPECT_EQ(gl.ledger->residual(), (Amounts {0, 10, 10}));
}

// V5: a rejected call of each kind (and an evaluator) changes nothing.
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
      {"capacity", [&] { (void)capacity(*g, v({0}), v({1}), query(cost_only(EB, LOSSY), offer(5))); }},
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

// ---------------------------------------------------------------- D: records

// D1: CostOnly and Global records have no entry; a Pinned record carries its entry identifier.
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

// D2, D3: one record per entry, duplicates included; a weight-0 entry's record at zero.
TEST(ProgramRecords, PerEntry) {
  const auto g = par2(4, 8);
  auto p = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 1}, {0, {1}, 1}, {0, {1}, 1}}));
  const auto r = p->place(12);
  const auto recs = p->records();
  ASSERT_EQ(recs.size(), 3U);
  EXPECT_EQ(entries(*p.p), (std::vector<std::optional<std::int32_t>> {0, 1, 2}));
  EXPECT_EQ(delivered(*p.p), (Amounts {4, 4, 4}));
  EXPECT_EQ(ascending_edges(*recs[1].bundle), ascending_edges(*recs[2].bundle));
  EXPECT_NE(recs[1].flow_id, recs[2].flow_id);
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 4, 4, 0}, {1, 0, UP, 4, 4, 0}, {2, 0, UP, 4, 4, 0}}));
  auto z = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 0}, {0, {1}, 1}}));
  (void)z->place(INF);
  const auto zr = z->records();
  ASSERT_EQ(zr.size(), 2U);
  EXPECT_EQ(zr[0].delivered, 0);
  EXPECT_TRUE(zr[0].edges.empty());
  EXPECT_EQ(ascending_edges(*zr[0].bundle), (std::vector<EdgeId> {0}));
}

// D4: a FAILED tunnel's record: an empty bundle, route_cost UNREACHABLE.
TEST(ProgramRecords, FailedEntry) {
  const auto g = graph(4, kG1U);
  const auto prog = g1u_program(g);
  auto p = pinned(g, prog, LOSSY, no_edges(*g, {2}));
  (void)p->replace(12);
  const auto recs = p->records();
  ASSERT_EQ(recs.size(), 2U);
  EXPECT_TRUE(ascending_edges(*recs[0].bundle).empty());
  EXPECT_EQ(recs[0].route_cost, kUnreachable);
  EXPECT_EQ(ascending_edges(*recs[1].bundle), (std::vector<EdgeId> {1, 3}));
  EXPECT_EQ(recs[1].route_cost, 2);
}

// D5: a CostOnly program in Pinned keeps the attribution of 5.1.
TEST(ProgramRecords, PooledContext) {
  const auto g = graph(4, kPool2);
  auto producer = make(g, {0, 1}, {3}, cost_only(EB, LOSSY));
  const auto prog = producer->place(2).program;
  auto p = pinned(g, prog, LOSSY);
  (void)p->replace(2);
  std::vector<std::pair<std::optional<std::int32_t>, Drops>> got;
  for (const auto& r : p->records()) got.emplace_back(r.entry, r.edges);
  EXPECT_EQ(got, (std::vector<std::pair<std::optional<std::int32_t>, Drops>> {{0, {{0, 1}, {2, 1}}},
                                                                             {1, {{1, 1}, {3, 1}}}}));
}

// ---------------------------------------------------------------- Z: results

// Z1, Z2, Z3: the by_tunnel rows.
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

// Z4: RC02-zero-unprepared extended: no program, no rows.
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

// Z5: TunnelFlow.assigned and dropped above 2^64 with a carry: 5 (2^62 - 1) =
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
  EXPECT_EQ(r.by_tunnel[0].dropped, total);
  EXPECT_EQ(r.by_tunnel[0].delivered, 0);
  EXPECT_EQ(r.offered, total);
  EXPECT_EQ(r.requested, total);
  EXPECT_EQ(r.admitted, total);
  EXPECT_EQ(r.drops, (std::vector<std::pair<EdgeId, Counter128>> {{0, total}}));
  EXPECT_TRUE(r.source_loss.is_zero());
  ASSERT_EQ(r.by_source.size(), 1U);
  EXPECT_EQ(r.by_source[0].offered, total);
  EXPECT_EQ(r.by_source[0].admitted, total);
  EXPECT_EQ(r.by_source[0].dropped, total);
  EXPECT_TRUE(r.by_source[0].source_loss.is_zero());
  const auto rec = p->records().front();
  EXPECT_EQ(rec.offered, total);
  EXPECT_EQ(rec.drops, (std::vector<std::pair<EdgeId, Counter128>> {{0, total}}));
  // The carry: four calls stay below 2^64.
  auto q = pinned(g, tunnels(*g, {0}, {1}, {{0, {0}, 1}}), LOSSY);
  PlacementResult four;
  for (int i = 0; i < 4; ++i) four = q->place(big);
  EXPECT_EQ(four.by_tunnel[0].assigned, (Counter128 {0, 4 * static_cast<std::uint64_t>(big)}));
}

// ---------------------------------------------------------------- E: evaluators

// E2: capacity of L3's Pinned Lossless.
TEST(QueryOptionsEval, PinnedLossless) {
  const auto g = par2(6, 4);
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}, {0, {1}, 1}});
  const auto r = capacity(*g, v({0}), v({1}), query(pinned_config(prog)));
  EXPECT_EQ(r.placed, 9);
  EXPECT_EQ(r.termination, Termination::Exhausted);
}

// E3: capacity of PINLOSS's Pinned Lossy with offered 12: 2 + 6.
TEST(QueryOptionsEval, PinnedLossy) {
  const auto pl = graph(4, kPinloss);
  const auto prog = tunnels(*pl, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}});
  EXPECT_EQ(capacity(*pl, v({0}), v({3}), query(pinned_config(prog, LOSSY), offer(12))).placed, 8);
}

// E4, E5: the offered rules and batch_capacity's Pinned text.
TEST(QueryOptionsEval, Texts) {
  const auto g = par2(6, 4);
  const auto prog = tunnels(*g, {0}, {1}, {{0, {0}, 1}, {0, {1}, 1}});
  auto err = [&](const QueryOptions& q) { return error_of([&] { (void)capacity(*g, v({0}), v({1}), q); }); };
  EXPECT_EQ(err(query(pinned_config(prog, LOSSY))), "Pinned LOSSY capacity requires finite offered");
  EXPECT_EQ(err(query(pinned_config(prog, LOSSY), offer_inf())), "Pinned LOSSY capacity requires finite offered");
  EXPECT_EQ(err(query(cost_only(EB, LOSSY), offer(5))), "CostOnly LOSSY capacity uses first-hop injection");
  EXPECT_EQ(err(query(cost_only(), offer(5))), "CostOnly LOSSLESS capacity does not accept offered");
  EXPECT_EQ(err(query(global(), offer(5))), "Global LOSSLESS capacity does not accept offered");
  EXPECT_EQ(err(query(pinned_config(prog), offer(5))), "Pinned LOSSLESS capacity does not accept offered");
  EXPECT_EQ(err(query(capacity_aware(), offer(5))), "T1 LOSSLESS capacity does not accept offered");
  EXPECT_EQ(err(query(cost_only(), offer_inf())), "CostOnly LOSSLESS capacity does not accept offered");
  for (const Amount bad : {Amount {-1}, kL62}) {
    EXPECT_EQ(err(query(pinned_config(prog, LOSSY), offer(bad))),
              "offered must be an integer in [0, 2^62) or permitted +inf");
  }
  EXPECT_EQ(err(query(pinned_config(prog, LOSSY), offer(kL62 - 1))), "");
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{{0}, {1}}};
  EXPECT_EQ(error_of([&] { (void)batch_capacity(*g, pairs, query(pinned_config(prog)), {}, 1); }),
            "batch_capacity does not accept Pinned");
}

// E6, E7, E8: a deletion fails the tunnels through it and the head splits over the others.
TEST(QueryOptionsEval, PinnedSensitivity) {
  const auto pl = graph(4, kPinloss);
  const auto prog = tunnels(*pl, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}});
  const auto s = sensitivity(*pl, v({0}), v({3}), query(pinned_config(prog)), 1);
  EXPECT_EQ(s.baseline.placed, 4);
  EXPECT_EQ(losses(s), (std::vector<std::pair<EdgeId, Amount>> {{0, -6}, {1, 2}, {2, -6}, {3, 2}}));
  const auto l = sensitivity(*pl, v({0}), v({3}), query(pinned_config(prog, LOSSY), offer(12)), 1);
  EXPECT_EQ(l.baseline.placed, 8);
  EXPECT_EQ(losses(l), (std::vector<std::pair<EdgeId, Amount>> {{0, -2}, {1, 6}, {2, -2}, {3, 6}}));
  const auto w2 = sensitivity(*pl, v({0}), v({3}), query(pinned_config(prog, LOSSY), offer(12)), 2);
  EXPECT_EQ(losses(w2), losses(l));
  EXPECT_EQ(w2.baseline.placed, 8);
}

// ---------------------------------------------------------------- M: today's Pinned mapped

// M1 to M10: each former Pinned use and its program form.
TEST(PinnedMapping, SingleBundle) {
  // M1: one EB bundle (DIAM, the four links).
  const auto d = graph(4, kDiam);
  const auto all = dag_of(*d, {0, 1, 2, 3});
  const auto eb = program(*d, {0}, {3}, {Route {0, all, 1}});
  EXPECT_EQ(pinned(d, eb)->place(INF).placed, 200);
  EXPECT_EQ(pinned(d, eb, LOSSY)->place(100).placed, 100);
  // M2: W with capacity weights; with residual weights (the ledger's residual at creation):
  // H over (100, 50) at node 0: H(150) = (100, 50) fits, H(151) = (101, 50) does not.
  const auto wc = program(*d, {0}, {3}, {Route {0, all, 1}}, RouteOptions {.placement = W});
  EXPECT_EQ(pinned(d, wc)->place(INF).placed, 200);
  const auto ws = program(*d, {0}, {3}, {Route {0, all, 1}},
                          RouteOptions {.placement = W, .weights = WeightBasis::ResidualSnapshot,
                                        .residual = Amounts {100, 50, 100, 50}});
  EXPECT_EQ(pinned(d, ws)->place(INF).placed, 150);
  // M3, M4: two bundles on PAR2(6, 4): separate contexts, RR over the entries: 9; W the same 9;
  // MaxFlow becomes Global(support) and fills both: 10.
  const auto p = par2(6, 4);
  const auto two = program(*p, {0}, {1}, {dag_route(*p, 0, {0}), dag_route(*p, 0, {1})});
  EXPECT_EQ(pinned(p, two)->place(INF).placed, 9);
  const auto two_w =
      program(*p, {0}, {1}, {dag_route(*p, 0, {0}), dag_route(*p, 0, {1})}, RouteOptions {.placement = W});
  EXPECT_EQ(pinned(p, two_w)->place(INF).placed, 9);
  EXPECT_EQ(make(p, {0}, {1}, global(std::nullopt, dag_of(*p, {0, 1})))->place(INF).placed, 10);
  // M5: two bundles sharing e0 on SHARE, EB Lossy place(10): both arrive at e0 with 5:
  // H(4, (5, 5)) = (2, 2).
  const auto sh = graph(3, kShare);
  auto shared_e0 = pinned(sh, program(*sh, {0}, {2}, {dag_route(*sh, 0, {0, 1}), dag_route(*sh, 0, {0, 2})}), LOSSY);
  const auto r5 = shared_e0->place(10);
  EXPECT_EQ(r5.placed, 4);
  EXPECT_EQ(delivered(*shared_e0.p), (Amounts {2, 2}));
  // M6: MaxFlow over one bundle is Global(support): 12, cost 34, one record per source with the
  // used edges, no entry.
  const auto es = graph(4, kES);
  auto m6 = make(es, {0}, {3}, global(std::nullopt, dag_of(*es, {0, 1, 2, 3})));
  const auto r6 = m6->place(INF);
  EXPECT_EQ(r6.placed, 12);
  EXPECT_EQ(r6.total_cost, c(34));
  const auto recs = m6->records();
  ASSERT_EQ(recs.size(), 1U);
  EXPECT_FALSE(recs[0].entry.has_value());
  EXPECT_EQ(ascending_edges(*recs[0].bundle), (std::vector<EdgeId> {0, 1, 2, 3}));
  // M7: Lossy +inf is rejected.
  auto m7 = pinned(d, eb, LOSSY);
  EXPECT_EQ(error_of([&] { (void)m7->place(INF); }), kFinite);
  // M8: two sources are accepted (C8); M9: duplicates are distinct entries (L1).
  const auto rel = graph(4, kRel);
  EXPECT_EQ(pinned(rel, tunnels(*rel, {0, 1}, {3}, {{0, {0, 2}, 1}, {1, {1, 2}, 1}}))->place(INF).placed, 10);
  const auto par = par2(4, 8);
  EXPECT_EQ(bound(par, {{0, {0}, 1}, {0, {1}, 1}, {0, {1}, 1}}), 12);
  // M10: an empty route list.
  EXPECT_EQ(error_of([&] { (void)program(*par, {0}, {1}, {}); }), "routes must not be empty");
}

// ---------------------------------------------------------------- storage of the new layouts

// A plan's replay under a mask (restrict_plan) is a new plan whose views point into its own block:
// it keeps its values after the source plan is destroyed and plans of the same shape reuse the
// freed memory (a borrowing restriction would read the later plans' roots and destinations).
TEST(ProgramLayouts, RestrictedPlanOwnsItsBlock) {
  const auto g = graph(8, {{0, 1, 5, 1}, {0, 2, 5, 1}, {1, 3, 5, 1}, {2, 3, 5, 1},
                           {4, 5, 5, 1}, {4, 6, 5, 1}, {5, 7, 5, 1}, {6, 7, 5, 1}});
  const auto mask = g->eligibility({}, std::vector<EdgeId> {3, 7});
  auto plan_of = [&](NodeId root, NodeId dst, std::vector<EdgeId> edges) {
    const auto dag = detail::PredDAGBuilder::of_edges(*g, edges);
    return detail::make_plan(*g, EB, LOSSLESS, detail::WeightSource {}, root, std::vector<NodeId> {dst}, dag);
  };
  std::optional<detail::KernelPlan> source(plan_of(0, 3, {0, 1, 2, 3}));
  const auto restricted = detail::restrict_plan(*g, *source, v({0}), v({3}), &mask);
  EXPECT_NE(restricted.roots.data(), source->roots.data());
  EXPECT_NE(restricted.dsts.data(), source->dsts.data());
  EXPECT_NE(restricted.topo.data(), source->topo.data());
  EXPECT_NE(restricted.head_pos.data(), source->head_pos.data());
  const auto roots = vec(restricted.roots);
  const auto dsts = vec(restricted.dsts);
  const auto topo = vec(restricted.topo);
  const auto members = restricted.members;
  EXPECT_EQ(members, (std::vector<EdgeId> {0, 2}));
  EXPECT_EQ(roots, v({0}));
  EXPECT_EQ(dsts, v({3}));
  source.reset();
  std::vector<detail::KernelPlan> later;
  for (int i = 0; i < 4; ++i) {
    later.push_back(plan_of(4, 7, {4, 5, 6, 7}));
    later.push_back(detail::restrict_plan(*g, later.back(), v({4}), v({7}), &mask));
  }
  EXPECT_EQ(vec(restricted.roots), roots);
  EXPECT_EQ(vec(restricted.dsts), dsts);
  EXPECT_EQ(vec(restricted.topo), topo);
  EXPECT_EQ(restricted.members, members);
  // Its copy owns another block too.
  const detail::KernelPlan copy(restricted);
  EXPECT_NE(copy.topo.data(), restricted.topo.data());
  EXPECT_EQ(vec(copy.topo), topo);
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
  inside(data.roster, "roster");
  inside(data.slot_first, "slot_first");
  inside(data.slot_entries, "slot_entries");
  inside(data.tunnel_entry, "tunnel_entry");
  inside(data.excluded_nodes, "excluded_nodes");
  inside(data.excluded_edges, "excluded_edges");
  for (const auto& c : data.contexts) inside(c.root_entry, "root_entry");
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
  EXPECT_EQ(vec(data.roster), v({0}));
  EXPECT_EQ(data.entries.size(), 2U);
  EXPECT_EQ(vec(data.entries_of(0)), v({0, 1}));
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
  EXPECT_EQ(vec(data.roster), v({0}));
  EXPECT_EQ(vec(data.excluded_nodes), v({4}));
  EXPECT_EQ(vec(data.excluded_edges), v({1, 4}));
  EXPECT_EQ(vec(data.tunnel_entry), v({0}));
  EXPECT_EQ(vec(data.slot_first), v({0, 1}));
  EXPECT_EQ(vec(data.slot_entries), v({0}));
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
    source.emplace(detail::make_view(*g, prog, nullptr, v({2}), nullptr, true));
  }
  ASSERT_EQ(source->num_heads(), 1U);
  ASSERT_FALSE(source->borrowed);  // a tunnel's record bundle is made: the view holds its state
  const detail::ForwardingView copy(*source);
  const auto stale = detail::with_black_holes(*source, &dark);
  EXPECT_NE(copy.status_.data(), source->status_.data());
  EXPECT_NE(copy.head_entries_.data(), source->head_entries_.data());
  EXPECT_EQ(copy.program, source->program);
  EXPECT_EQ(stale.program, source->program);
  EXPECT_EQ(stale.dark, &dark);
  const auto status = copy.status_;
  const auto entries_of_head = copy.head_entries_;
  const auto plans = copy.plans_.size();
  source.reset();
  std::vector<std::shared_ptr<const ForwardingProgram>> later;
  for (int i = 0; i < 4; ++i) later.push_back(build(3, {4}, {3, 5}, 5));
  for (const auto* view : {&copy, &stale}) {
    ASSERT_NE(view->program, nullptr);
    EXPECT_EQ(view->data, &detail::ProgramAccess::data(*view->program));
    EXPECT_EQ(view->head(0), 0);
    EXPECT_EQ(vec(view->data->roster), v({0}));
    EXPECT_EQ(vec(view->data->destinations), v({2}));
    expect_lists_in_block(*view->data);
    EXPECT_EQ(view->status_, status);
    EXPECT_EQ(view->head_entries_, entries_of_head);
    ASSERT_EQ(view->plans_.size(), plans);
    for (const auto& plan : view->plans_) EXPECT_EQ(vec(plan->dsts), v({2}));
  }
}

// A borrowed view (no mask, no black holes, no record bundle to make) holds no state of its own:
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
    source.emplace(detail::make_view(*g, prog, nullptr, v({3}), nullptr, true));
  }
  ASSERT_TRUE(source->borrowed);
  EXPECT_TRUE(source->status_.empty() && source->plans_.empty() && source->head_entries_.empty());
  const detail::ForwardingView copy(*source);
  const auto stale = detail::with_black_holes(*source, &dark);
  EXPECT_TRUE(copy.borrowed);
  EXPECT_FALSE(stale.borrowed);
  const auto* plan = source->plan(0).get();
  source.reset();
  std::vector<std::shared_ptr<const ForwardingProgram>> later;
  for (int i = 0; i < 4; ++i) later.push_back(build(0));
  for (const auto* view : {&copy, &stale}) {
    ASSERT_NE(view->program, nullptr);
    EXPECT_EQ(view->data, &detail::ProgramAccess::data(*view->program));
    ASSERT_EQ(view->num_heads(), 1U);
    EXPECT_EQ(view->head(0), 0);
    EXPECT_EQ(view->status(0), TunnelStatus::Up);
    EXPECT_EQ(vec(view->head_entries(0)), v({0}));
    EXPECT_EQ(view->root_of(0), 0);
    EXPECT_EQ(view->plan(0).get(), plan);
    EXPECT_EQ(view->plan(0)->members, v({0, 1, 2, 3}));
    EXPECT_EQ(view->cost(0), 2);
    expect_lists_in_block(*view->data);
  }
  EXPECT_FALSE(copy.dark_head(0));
}

// An epoch is a value: its copy keeps the cohorts, the program, the rows and the counters after
// the policy that held the original has moved on and been destroyed.
TEST(ProgramLayouts, EpochCopyOwnsItsStorage) {
  const auto g = graph(4, kPinloss);
  const auto prog = tunnels(*g, {0}, {3}, {{0, {0, 2}, 1}, {0, {1, 3}, 1}});
  std::optional<Pol> p(pinned(g, prog, LOSSY));
  (void)(*p)->replace(12);
  const detail::Epoch copy = detail::PolicyAccess::engine(*p->p).epoch;
  const auto& live = detail::PolicyAccess::engine(*p->p).epoch;
  ASSERT_EQ(copy.cohorts.size(), 2U);
  EXPECT_NE(copy.cohorts.data(), live.cohorts.data());
  EXPECT_NE(copy.by_tunnel.data(), live.by_tunnel.data());
  const auto rows_before = copy.by_tunnel;
  const auto ids_before = std::vector<std::int64_t> {copy.cohorts[0][0].id, copy.cohorts[1][0].id};
  (void)(*p)->replace(4);
  p.reset();
  for (int i = 0; i < 4; ++i) (void)pinned(g, prog, LOSSY)->replace(2);
  EXPECT_EQ(copy.program, prog);
  EXPECT_TRUE(copy.by_tunnel == rows_before);
  EXPECT_EQ(rows_before.size(), 2U);
  EXPECT_EQ(lo(rows_before[0].assigned), 6);
  EXPECT_EQ(rows_before[0].delivered, 2);
  EXPECT_EQ((std::vector<std::int64_t> {copy.cohorts[0][0].id, copy.cohorts[1][0].id}), ids_before);
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
          if (g->edge_src()[static_cast<std::size_t>(e)] != u) continue;
          stack.push_back(e);
          walk(g->edge_dst()[static_cast<std::size_t>(e)]);
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
    ASSERT_EQ(r.placed, x) << where;
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
      out[static_cast<std::size_t>(g->edge_src()[e])] += (*r.edge_flows)[e];
      in[static_cast<std::size_t>(g->edge_dst()[e])] += (*r.edge_flows)[e];
    }
    for (NodeId u = 2; u < dst; ++u) EXPECT_EQ(in[u], out[u]) << where << " node " << u;
    EXPECT_EQ(in[dst], x) << where;
    EXPECT_EQ(out[0] + out[1] - in[1], x) << where;
    // Lossy at the bound: every share delivered, nothing dropped.
    const auto l = pinned(g, prog, LOSSY)->replace(x);
    EXPECT_EQ(l.placed, x) << where;
    EXPECT_TRUE(l.call_drops.empty()) << where;
    EXPECT_TRUE(l.source_loss.is_zero()) << where;
    // Lossless above the bound places the bound.
    EXPECT_EQ(pinned(g, prog)->evaluate(x + 1).placed, x) << where;
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
