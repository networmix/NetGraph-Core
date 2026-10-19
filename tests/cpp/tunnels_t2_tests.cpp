// T2 (Global(tunnels=Tunnels(N)), Tunnels(L, up_to)): a controller lays out, for every live head,
// the same number of single-path tunnels with one common whole-unit planned bandwidth (contract
// section 7): the contract's T2 examples and calls (7.4 to 7.9, 8.5), sizing, the objectives and
// the canonical tie, count selection, operations, programs, replays and stale forwarding (8.1, 8.4,
// 8.6), records, rows, flags and terminations (10.4, 14.1, 14.2), T2's evaluators (12), the
// numeric ranges of its sums, zero offers and everything masked. Every
// expected value is derived by hand or copied from the contract's examples, never from a build;
// the property test compares random small cases with t2_reference (test_utils.hpp), a brute force
// over path multisets that calls none of the library's flow, search or apportionment code.
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/tunnels.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/program.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

constexpr Admission LOSSLESS = Admission::Lossless;
constexpr Admission LOSSY = Admission::Lossy;
constexpr TunnelStatus UP = TunnelStatus::Up;
constexpr TunnelStatus FAILED = TunnelStatus::Failed;
constexpr TunnelStatus BLACKHOLED = TunnelStatus::Blackholed;
const std::optional<Amount> INF;  // +inf
const std::optional<std::optional<Amount>> UNLIMITED {std::in_place, std::nullopt};  // design +inf

using Amounts = std::vector<Amount>;
using Edges = std::vector<EdgeId>;
using Pairs = std::vector<std::pair<EdgeId, Amount>>;
using SourceRow = std::tuple<NodeId, Amount, Amount, Amount, Amount, Amount>;

std::shared_ptr<const Graph> ordered(std::int32_t n, const std::vector<InEdge>& rows) {
  auto g = shared(make_graph(n, rows));
  for (std::size_t i = 0; i < rows.size(); ++i) EXPECT_EQ(g->input_index()[i], static_cast<std::int32_t>(i));
  return g;
}

FlowPolicyConfig t2(std::int64_t n, Admission a = LOSSLESS, bool headroom = false, bool up_to = false,
                    std::optional<Origin> origin = std::nullopt) {
  return {.routing = Global {.tunnels = Tunnels {n, up_to}, .prefer_headroom = headroom}, .admission = a, .origin = origin};
}

Amount lo(const Counter128& c) {
  EXPECT_EQ(c.hi, 0U);
  return static_cast<Amount>(c.lo);
}
Amount lo(const std::optional<Counter128>& c) { return c.has_value() ? lo(*c) : -1; }
Int128 wide(const Counter128& c) { return int128_of(static_cast<std::int64_t>(c.hi), c.lo); }

std::vector<SourceRow> source_rows(const PlacementResult& r) {
  std::vector<SourceRow> out;
  for (const auto& s : r.by_source) out.emplace_back(s.source, lo(s.offered), lo(s.admitted), s.delivered, lo(s.source_loss), lo(s.dropped));
  return out;
}
std::vector<Edges> routes(const ForwardingProgram& p) {
  std::vector<Edges> out;
  for (const auto& t : p.tunnels()) out.push_back(t.edges.value_or(Edges {}));
  return out;
}
std::vector<std::int32_t> indices(const ForwardingProgram& p) {
  std::vector<std::int32_t> out;
  for (const auto& t : p.tunnels()) out.push_back(t.index);
  return out;
}
Amounts counts(const ForwardingProgram& p, std::size_t m) {
  Amounts x(m, 0);
  for (const auto& t : p.tunnels()) {
    for (EdgeId e : *t.edges) ++x[static_cast<std::size_t>(e)];
  }
  return x;
}
std::vector<TunnelStatus> row_statuses(const PlacementResult& r) {
  std::vector<TunnelStatus> out;
  for (const auto& t : r.by_tunnel) out.push_back(t.status);
  return out;
}

struct Pol {
  std::shared_ptr<FlowLedger> ledger;
  std::unique_ptr<FlowPolicy> p;
  FlowPolicy* operator->() const { return p.get(); }
};
Pol make(const std::shared_ptr<const Graph>& g, std::vector<NodeId> s, std::vector<NodeId> d, const FlowPolicyConfig& c,
         const Eligibility* elig = nullptr, std::shared_ptr<FlowLedger> ledger = nullptr) {
  Pol out;
  out.ledger = ledger ? std::move(ledger) : std::make_shared<FlowLedger>(g);
  out.p = std::make_unique<FlowPolicy>(*out.ledger, s, d, c, elig);
  return out;
}

const std::vector<InEdge> kG1 = {{0, 1, 8192, 1}, {0, 2, 40960, 1}, {1, 3, 8192, 1}, {2, 3, 40960, 1}};
std::shared_ptr<const Graph> g1() { return ordered(4, kG1); }
std::shared_ptr<const Graph> one10() { return ordered(2, {{0, 1, 10, 1}}); }
std::shared_ptr<const Graph> twoh(Amount a, Amount b) { return ordered(3, {{0, 2, a, 1}, {1, 2, b, 1}}); }
std::vector<Edges> g1_paths() {
  std::vector<Edges> out(2, Edges {0, 2});
  out.insert(out.end(), 14, Edges {1, 3});
  return out;
}

// ---------------------------------------------------------------- the contract's examples (X)

struct Example {
  std::string id;
  std::int32_t n;
  std::vector<InEdge> rows;
  std::vector<NodeId> heads;
  NodeId dst;
  std::int64_t N;
  std::optional<Amount> D;  // nullopt: unlimited
  Amounts free;             // empty: the capacities
  Amount U, u;
  std::pair<Amount, Amount> fluid;
  Amount F, total, cost;
  std::optional<Amount> phi;
  std::vector<Edges> routes;
  std::optional<Amount> offer;  // Lossy
  Amount delivery, loss;
};

std::vector<Example> examples() {
  std::vector<InEdge> g4;
  for (std::int32_t i = 1; i <= 15; ++i) g4.push_back({0, i, 4096, 1});
  for (std::int32_t i = 1; i <= 15; ++i) g4.push_back({i, 16, 4096, 1});
  std::vector<Edges> g4_routes(2, Edges {0, 15});
  for (EdgeId i = 1; i <= 14; ++i) g4_routes.push_back({i, 15 + i});
  const std::vector<InEdge> g10 = {
      {0, 1, 40960, 1}, {0, 5, 40960, 1}, {0, 3, 24576, 2}, {1, 0, 40960, 1}, {1, 2, 40960, 1}, {1, 4, 16384, 2},
      {2, 1, 40960, 1}, {2, 3, 40960, 1}, {2, 5, 20480, 2}, {3, 2, 40960, 1}, {3, 4, 40960, 1}, {3, 0, 24576, 2},
      {4, 3, 40960, 1}, {4, 5, 40960, 1}, {4, 1, 16384, 2}, {5, 0, 40960, 1}, {5, 4, 40960, 1}, {5, 2, 20480, 2}};
  auto rep = [](std::vector<std::pair<Edges, int>> parts) {
    std::vector<Edges> out;
    for (const auto& [p, c] : parts) out.insert(out.end(), static_cast<std::size_t>(c), p);
    return out;
  };
  return {
      {"G1 design 49152", 4, kG1, {0}, 3, 16, 49152, {}, 2925, 2925, {20480, 7}, 46812, 46800, 32, 73610035200, g1_paths(), 49152, 47104, 0},
      {"G1 design 8192", 4, kG1, {0}, 3, 16, 8192, {}, 2925, 512, {20480, 7}, 46812, 8192, 32, 12884901876, g1_paths(), {}, 0, 0},
      {"G4 design 61440", 17, g4, {0}, 16, 16, 61440, {}, 2048, 2048, {2048, 1}, 32768, 32768, 32, 73014444032, g4_routes, 61440, 57856, 0},
      {"heads behind 3 and 8", 3, {{0, 2, 3, 1}, {1, 2, 8, 1}}, {0, 1}, 2, 2, 6, {}, 1, 1, {3, 2}, 6, 4, 4, 5905580031,
       {{0}, {1}, {0}, {1}}, {}, 0, 0},
      {"a cheap link of 2", 2, {{0, 1, 2, 1}, {0, 1, 10, 2}}, {0}, 1, 2, 2, {}, 5, 1, {5, 1}, 2, 2, 2, 6442450944, {{0}, {0}}, {}, 0, 0},
      {"G10 head 0", 6, g10, {0}, 3, 16, 49152,
       {40960, 40960, 24576, 40960, 40960, 16384, 40960, 40960, 20480, 40960, 40960, 24576, 40960, 40960, 16384, 40960, 40960, 20480},
       6144, 3072, {6144, 1}, 49160, 49152, 40, 38654705652, rep({{{0, 4, 7}, 4}, {{1, 16, 12}, 4}, {{2}, 8}}), {}, 0, 0},
      {"G10 head 1", 6, g10, {1}, 4, 16, 49152,
       {28672, 28672, 0, 40960, 28672, 16384, 40960, 28672, 20480, 40960, 40960, 24576, 28672, 40960, 16384, 40960, 28672, 20480},
       4096, 3072, {4096, 1}, 52431, 49152, 43, {}, rep({{{3, 1, 16}, 6}, {{4, 7, 10}, 5}, {{5}, 5}}), {}, 0, 0},
      {"G10 head 2", 6, g10, {2}, 5, 16, 49152,
       {28672, 10240, 0, 22528, 13312, 1024, 40960, 13312, 20480, 40960, 25600, 24576, 28672, 40960, 16384, 40960, 10240, 20480},
       2560, 2560, {2560, 1}, 40960, 40960, 40, {}, rep({{{6, 3, 1}, 4}, {{7, 10, 13}, 4}, {{8}, 8}}), {}, 0, 0},
      {"a dead head", 3, {{0, 2, 0, 1}, {1, 2, 10, 1}}, {0, 1}, 2, 2, 10, {}, 5, 3, {5, 1}, 0, 6, 2, 3865470565, {{1}, {1}}, 10, 5, 5},
      {"heads behind 2 and 10", 3, {{0, 2, 2, 1}, {1, 2, 10, 1}}, {0, 1}, 2, 2, {}, {}, 1, 1, {1, 1}, 4, 4, 4, 7730941132,
       {{0}, {1}, {0}, {1}}, {}, 0, 0},
      {"an idle tunnel", 2, {{0, 1, 10, 1}}, {0}, 1, 4, 3, {}, 2, 1, {5, 2}, 10, 4, 4, 4294967294, {{0}, {0}, {0}, {0}}, 3, 3, 0},
      {"a zero design", 2, {{0, 1, 10, 1}}, {0}, 1, 2, 0, {}, 5, 1, {5, 1}, 10, 2, 2, {}, {{0}, {0}}, 4, 4, 0},
      {"the remainder on the first link", 2, {{0, 1, 3, 1}, {0, 1, 2, 1}}, {0}, 1, 2, 5, {}, 2, 2, {2, 1}, 5, 4, 2, {}, {{0}, {1}}, {}, 0, 0},
      {"the remainder on the second link", 2, {{0, 1, 2, 1}, {0, 1, 3, 1}}, {0}, 1, 2, 5, {}, 2, 2, {2, 1}, 4, 4, 2, {}, {{0}, {1}}, {}, 0, 0},
      {"a tiny head", 3, {{0, 2, 1, 1}, {1, 2, 100, 1}}, {0, 1}, 2, 16, {}, {}, 0, 0, {1, 16}, 0, 0, 0, 0, {}, {}, 0, 0},
  };  // clang-format: the table stays one row per example
}

// Every T2 example of the contract (7.3 to 7.6), with prefer_headroom as the contract runs them.
TEST(T2Objectives, ContractExamples) {
  for (const auto& x : examples()) {
    SCOPED_TRACE(x.id);
    const auto g = ordered(x.n, x.rows);
    auto ledger = x.free.empty() ? std::make_shared<FlowLedger>(g) : std::make_shared<FlowLedger>(g, x.free);
    auto p = make(g, x.heads, {x.dst}, t2(x.N, LOSSLESS, true), nullptr, ledger);
    const auto r = p->replace(INF, x.D.has_value() ? std::optional<std::optional<Amount>>(x.D) : UNLIMITED);
    ASSERT_TRUE(r.plan.has_value());
    const auto& pl = *r.plan;
    EXPECT_EQ(pl.common_share, x.U);
    EXPECT_EQ(pl.planned_per_tunnel, x.u);
    EXPECT_EQ(pl.fluid_share, x.fluid);
    EXPECT_EQ(pl.forwarding_bound, x.F);
    EXPECT_EQ(pl.planned_total, x.total);
    EXPECT_EQ(lo(pl.route_cost), x.cost);
    if (x.phi.has_value()) EXPECT_EQ(lo(pl.headroom_penalty), *x.phi);
    EXPECT_EQ(pl.insufficient_resolution, x.U == 0);
    EXPECT_EQ(routes(*r.program), x.routes);
    for (const auto& t : r.program->tunnels()) {
      EXPECT_EQ(t.planned, x.u);
      EXPECT_EQ(t.weight, 1);
      EXPECT_EQ(t.status, UP);
      EXPECT_EQ(t.requested, 0);
      EXPECT_EQ(t.reserved, 0);
    }
    EXPECT_EQ(r.placed, x.F);
    EXPECT_EQ(r.termination, Termination::Exhausted);
    if (!x.offer.has_value()) continue;
    auto lossy_ledger = x.free.empty() ? std::make_shared<FlowLedger>(g) : std::make_shared<FlowLedger>(g, x.free);
    auto q = make(g, x.heads, {x.dst}, t2(x.N, LOSSY, true), nullptr, lossy_ledger);
    const auto s = q->replace(x.offer, x.D.has_value() ? std::optional<std::optional<Amount>>(x.D) : UNLIMITED);
    EXPECT_EQ(s.placed, x.delivery);
    EXPECT_EQ(s.call_source_loss, x.loss);
    EXPECT_EQ(routes(*s.program), x.routes);
  }
}

// The 7.9 call (G1 at q = 2^-12, N = 16, prefer_headroom, Lossy replace(49152)).
TEST(T2Objectives, ContractCall) {
  const auto g = g1();
  auto p = make(g, {0}, {3}, t2(16, LOSSY, true));
  const auto r = p->replace(Amount {49152});
  ASSERT_TRUE(r.plan.has_value());
  const auto& pl = *r.plan;
  EXPECT_EQ(pl.count, 16);
  EXPECT_EQ(pl.design_demand, Amount {49152});
  EXPECT_EQ(pl.head_design, (std::vector<std::pair<NodeId, Amount>> {{0, 49152}}));
  EXPECT_EQ(pl.reserved_total, 0);
  EXPECT_EQ(pl.live_heads, (std::vector<NodeId> {0}));
  EXPECT_EQ(pl.insufficient_resolution, false);
  EXPECT_EQ(pl.fluid_share, (std::pair<Amount, Amount> {20480, 7}));
  EXPECT_EQ(pl.common_share, 2925);
  EXPECT_EQ(pl.planned_per_tunnel, 2925);
  EXPECT_EQ(pl.planned_total, 46800);
  EXPECT_EQ(lo(pl.route_cost), 32);
  EXPECT_EQ(lo(pl.headroom_penalty), 73610035200);
  EXPECT_EQ(pl.forwarding_bound, 46812);
  EXPECT_EQ(r.placed, 47104);
  EXPECT_EQ(routes(*r.program), g1_paths());
  Amounts delivered;
  for (const auto& t : r.by_tunnel) delivered.push_back(t.delivered);
  Amounts want(2, 3072);
  want.insert(want.end(), 10, 2926);
  want.insert(want.end(), 4, 2925);
  EXPECT_EQ(delivered, want);
  EXPECT_EQ(r.call_drops, (Pairs {{1, 2048}}));
  EXPECT_FALSE(r.reservations_changed);
}

// A plan books nothing: no hold, no reserved_edges.
TEST(LedgerProjections, T2BooksNothing) {
  const auto g = g1();
  auto p = make(g, {0}, {3}, t2(16, LOSSY, true));
  (void)p->replace(Amount {49152});
  EXPECT_EQ(p.ledger->reserved(), (Amounts {0, 0, 0, 0}));
  EXPECT_EQ(p.ledger->unreserved(), (Amounts {8192, 40960, 8192, 40960}));
  for (const auto& rec : p->records()) EXPECT_TRUE(rec.reserved_edges.empty());
  EXPECT_EQ(p->records().size(), 16U);
}

// Without prefer_headroom the cost tie goes to the smallest keys.
TEST(T2Objectives, HeadroomOff) {
  const auto g = g1();
  auto p = make(g, {0}, {3}, t2(16, LOSSY));
  const auto r = p->replace(Amount {49152});
  EXPECT_EQ(counts(*r.program, 4), (Amounts {2, 14, 2, 14}));
  EXPECT_EQ(r.plan->forwarding_bound, 46812);
  EXPECT_EQ(r.placed, 47104);
  EXPECT_FALSE(r.plan->headroom_penalty.has_value());
  auto q = make(g, {0}, {3}, t2(16));
  const auto s = q->replace(INF, Amount {8192});
  EXPECT_EQ(s.plan->planned_per_tunnel, 512);
  EXPECT_EQ(counts(*s.program, 4), (Amounts {16, 0, 16, 0}));
  EXPECT_EQ(s.plan->forwarding_bound, 8192);
  const auto par = ordered(2, {{0, 1, 4, 1}, {0, 1, 12, 1}});
  auto h = make(par, {0}, {1}, t2(2, LOSSLESS, true));
  const auto a = h->replace(INF, Amount {2});
  EXPECT_EQ(a.plan->common_share, 6);
  EXPECT_EQ(a.plan->planned_per_tunnel, 1);
  EXPECT_EQ(routes(*a.program), (std::vector<Edges> {{1}, {1}}));
  EXPECT_EQ(a.plan->forwarding_bound, 12);
  EXPECT_EQ(lo(a.plan->headroom_penalty), 1073741823);
  auto c = make(par, {0}, {1}, t2(2));
  const auto b = c->replace(INF, Amount {2});
  EXPECT_EQ(routes(*b.program), (std::vector<Edges> {{0}, {0}}));
  EXPECT_EQ(b.plan->forwarding_bound, 4);
}

// prefer_headroom off and on: two equal-cost parallel links of 10, one head, N = 2, D = 2: U = 10,
// u = 1. Off, the canonical tie puts both tunnels on link 0 (bound 10, Lossy 12 delivers 10); on,
// the second tunnel pays floor(2^32 / 10) on link 1 against floor(2^32 2 / 10) for a second slot on
// link 0 (Phi = 2 floor(2^32 / 10) = 858993458, bound 20, Lossy 12 delivers 12). Both plan 2.
TEST(T2Objectives, PreferHeadroomExample) {
  const auto g = ordered(2, {{0, 1, 10, 1}, {0, 1, 10, 1}});
  for (const bool on : {false, true}) {
    SCOPED_TRACE(on);
    auto p = make(g, {0}, {1}, t2(2, LOSSY, on));
    const auto r = p->replace(Amount {12}, Amount {2});
    EXPECT_EQ(r.plan->common_share, 10);
    EXPECT_EQ(r.plan->planned_per_tunnel, 1);
    EXPECT_EQ(r.plan->planned_total, 2);
    EXPECT_EQ(lo(r.plan->route_cost), 2);
    EXPECT_EQ(counts(*r.program, 2), on ? (Amounts {1, 1}) : (Amounts {2, 0}));
    EXPECT_EQ(r.plan->headroom_penalty, on ? std::optional<Counter128>(Counter128 {0, 858993458}) : std::nullopt);
    EXPECT_EQ(r.plan->forwarding_bound, on ? 20 : 10);
    EXPECT_EQ(r.placed, on ? 12 : 10);
  }
}

// A reverse arc (the second unit cancels the first unit's middle link) and a preload (the slot prices
// read the static capacity B: ONE10 with 6 units free, N = 2, D = 4: U = 3, u = 2, Phi =
// floor(2^32 6 / 10) + floor(2^32 8 / 10)).
TEST(T2Objectives, ReverseArcsAndPreloadedPrices) {
  const auto g = ordered(4, {{0, 1, 1, 1}, {0, 2, 1, 5}, {1, 2, 1, 1}, {1, 3, 1, 5}, {2, 3, 1, 1}});
  auto p = make(g, {0}, {3}, t2(2));
  const auto r = p->replace(INF);
  EXPECT_EQ(r.plan->common_share, 1);
  EXPECT_EQ(routes(*r.program), (std::vector<Edges> {{0, 3}, {1, 4}}));
  EXPECT_EQ(lo(r.plan->route_cost), 12);
  const auto one = one10();
  const Amounts six {6};
  auto q = make(one, {0}, {1}, t2(2, LOSSLESS, true), nullptr, std::make_shared<FlowLedger>(one, six));
  const auto s = q->replace(INF, Amount {4});
  EXPECT_EQ(s.plan->common_share, 3);
  EXPECT_EQ(s.plan->planned_per_tunnel, 2);
  EXPECT_EQ(lo(s.plan->headroom_penalty), 2576980377 + 3435973836);
}

// Heads behind 3 and 8: the tunnels in round order, local indices 0, 0, 1, 1.
TEST(T2Objectives, RoundIndices) {
  const auto g = ordered(3, {{0, 2, 3, 1}, {1, 2, 8, 1}});
  auto p = make(g, {0, 1}, {2}, t2(2));
  const auto r = p->replace(INF, Amount {6});
  std::vector<std::tuple<std::int32_t, NodeId, std::int32_t, Edges>> got;
  for (const auto& t : r.program->tunnels()) got.emplace_back(t.index, t.head, t.local_index, *t.edges);
  EXPECT_EQ(got, (std::vector<std::tuple<std::int32_t, NodeId, std::int32_t, Edges>> {
                     {0, 0, 0, {0}}, {1, 1, 0, {1}}, {2, 0, 1, {0}}, {3, 1, 1, {1}}}));
}

std::vector<std::pair<NodeId, Edges>> head_routes(const ForwardingProgram& p) {
  std::vector<std::pair<NodeId, Edges>> out;
  for (const auto& t : p.tunnels()) out.emplace_back(t.head, t.edges.value_or(Edges {}));
  return out;
}

// The root-arc keys (1, head NodeId, direction): when the heads' first walks tie (price 2, four arcs
// each), the smaller head goes first. Every link capacity 1, cost 1, N = 1, U = u = 1, every slot
// price 2^32, so prefer_headroom changes nothing. (a) e0 0-2, e1 0-4, e2 1-2, e3 1-5, e4 2-3, e5 4-3,
// e6 5-3: head 0 takes e0 e4 (smallest link key at 0); head 1 then takes e3 e6 (four arcs) rather
// than e2, e0 reversed, e1, e5 (six). Heads in descending order would give head 1 e2 e4 and head 0
// e1 e5. (b) e0 0-2, e1 0-5, e2 1-2, e3 1-4, e4 2-3, e5 4-3, e6 5-3: head 0 e0 e4, head 1 e3 e5;
// descending: head 1 e2 e4, head 0 e1 e6.
TEST(T2Objectives, RootKeyTie) {
  const auto a = ordered(6, {{0, 2, 1, 1}, {0, 4, 1, 1}, {1, 2, 1, 1}, {1, 5, 1, 1}, {2, 3, 1, 1}, {4, 3, 1, 1}, {5, 3, 1, 1}});
  const auto b = ordered(6, {{0, 2, 1, 1}, {0, 5, 1, 1}, {1, 2, 1, 1}, {1, 4, 1, 1}, {2, 3, 1, 1}, {4, 3, 1, 1}, {5, 3, 1, 1}});
  for (const bool headroom : {false, true}) {
    auto p = make(a, {0, 1}, {3}, t2(1, LOSSLESS, headroom));
    EXPECT_EQ(head_routes(*p->replace(INF).program), (std::vector<std::pair<NodeId, Edges>> {{0, {0, 4}}, {1, {3, 6}}}));
    auto q = make(b, {0, 1}, {3}, t2(1, LOSSLESS, headroom));
    EXPECT_EQ(head_routes(*q->replace(INF).program), (std::vector<std::pair<NodeId, Edges>> {{0, {0, 4}}, {1, {3, 5}}}));
  }
}

// Fewest arcs before keys (7.3), prefer_headroom off: the direct link ties the longer walk on
// cost and has fewer arcs; the smallest keys alone would take the walk. (a) e0 0-1 (cost 1), e1 0-2
// (cost 2), e2 1-2 (cost 1): e1 (three arcs against four), keys alone e0 e2. (b) e0 0-1, e1 0-3
// (cost 3), e2 1-2, e3 2-3 (cost 1 each): e1, keys alone e0 e2 e3. Capacities 1, N = 1.
TEST(T2Objectives, FewestArcsBeforeKeys) {
  auto p = make(ordered(3, {{0, 1, 1, 1}, {0, 2, 1, 2}, {1, 2, 1, 1}}), {0}, {2}, t2(1));
  EXPECT_EQ(routes(*p->replace(INF).program), (std::vector<Edges> {{1}}));
  auto q = make(ordered(4, {{0, 1, 1, 1}, {0, 3, 1, 3}, {1, 2, 1, 1}, {2, 3, 1, 1}}), {0}, {3}, t2(1));
  EXPECT_EQ(routes(*q->replace(INF).program), (std::vector<Edges> {{1}}));
}

// The cut into unit walks (7.3): live heads ascending, each walk along the smallest positive EdgeId.
// (a) e0 0-1 (capacity 2), e1 1-2 (2), e2 1-2 (3), cost 2 each, heads 0 and 1, N = 2, U = u = 1:
// head 1's two units take e1 (cost 2 against head 0's 4), head 0's e0 then e2; x = (2, 2, 2). Head 0
// cuts first: e0 e1 twice, leaving e2 twice to head 1 (the largest EdgeId would give head 0 e0 e2
// and head 1 e1). (b) e0 0-2, e1 1-2, e2 2-3, e3 2-4, e4 3-5, e5 4-5, capacity 1, cost 1, N = 1:
// x is all ones; head 0 takes e0 e2 e4, head 1 e1 e3 e5 (the largest EdgeId: e0 e3 e5 and e1 e2 e4).
TEST(T2Objectives, CutIntoUnitWalks) {
  auto p = make(ordered(3, {{0, 1, 2, 2}, {1, 2, 2, 2}, {1, 2, 3, 2}}), {0, 1}, {2}, t2(2));
  EXPECT_EQ(head_routes(*p->replace(INF).program),
            (std::vector<std::pair<NodeId, Edges>> {{0, {0, 1}}, {1, {2}}, {0, {0, 1}}, {1, {2}}}));
  auto q = make(ordered(6, {{0, 2, 1, 1}, {1, 2, 1, 1}, {2, 3, 1, 1}, {2, 4, 1, 1}, {3, 5, 1, 1}, {4, 5, 1, 1}}), {0, 1},
                {5}, t2(1));
  EXPECT_EQ(head_routes(*q->replace(INF).program),
            (std::vector<std::pair<NodeId, Edges>> {{0, {0, 2, 4}}, {1, {1, 3, 5}}}));
}

// The route cost above 2^64 (WIDE: one link of 16 at cost 2^61 + 1, N = 16, D = 16).
TEST(T2Objectives, WideRouteCost) {
  const auto g = ordered(2, {{0, 1, 16, (Cost {1} << 61) + 1}});
  auto p = make(g, {0}, {1}, t2(16));
  const auto r = p->replace(INF, Amount {16});
  ASSERT_TRUE(r.plan->route_cost.has_value());
  EXPECT_EQ(wide(*r.plan->route_cost), Int128(16) * Int128((std::int64_t {1} << 61) + 1));
  EXPECT_EQ(r.plan->route_cost->hi, 2U);
}

// A slot price is at most 2^32 (j u <= V): a full link of 10 with 10 tunnels at u = 1
// pays floor(2^32 j / 10) for j = 1..10, the last exactly 2^32: 55 * 2^32 / 10 minus the dropped
// fractions (0.6, 0.2, 0.8, 0.4, 0, 0.6, 0.2, 0.8, 0.4, 0: 4 in all) = 23622320124.
TEST(T2Objectives, PenaltyBound) {
  auto p = make(one10(), {0}, {1}, t2(10, LOSSLESS, true));
  const auto r = p->replace(INF, Amount {10});
  EXPECT_EQ(r.plan->planned_per_tunnel, 1);
  EXPECT_EQ(lo(r.plan->headroom_penalty), 23622320124);
}

// floor(2^32 m / B) for 0 <= m <= B < 2^62, by binary long division in 64 bits (independent of the
// library's 128-bit helpers).
std::uint64_t slot_price_by_division(std::uint64_t m, std::uint64_t B) {
  if (m == B) return std::uint64_t {1} << 32;
  std::uint64_t rem = m;
  std::uint64_t q = 0;
  for (int bit = 0; bit < 32; ++bit) {
    rem <<= 1;  // below 2 B < 2^63
    q <<= 1;
    if (rem >= B) {
      rem -= B;
      q |= 1;
    }
  }
  return q;
}

// The closed-form penalty of t2_objectives equals the slot-by-slot sum of 7.3 on one link (random
// B up to 2^62 - 1, V <= B, u and x with x u <= V, x up to 3000), summed here with its own carry.
TEST(T2Objectives, PenaltySumsEverySlot) {
  std::mt19937_64 rng(4242);
  const std::vector<EdgeId> edges {0};
  for (int trial = 0; trial < 3000; ++trial) {
    const int scale = trial % 4;
    const Amount B = scale == 0   ? Amount(1 + rng() % 20)
                     : scale == 1 ? Amount(1 + rng() % 100000)
                     : scale == 2 ? Amount(1 + rng() % (Amount {1} << 40))
                                  : kL62 - 1 - Amount(rng() % 1000);
    const Amount V = Amount(rng() % static_cast<std::uint64_t>(B + 1));
    if (V == 0) continue;
    const Amount u = trial % 3 == 0 ? 1 : 1 + Amount(rng() % static_cast<std::uint64_t>(V));
    const Amount x = std::min<Amount>(V / u, Amount(rng() % 3001));
    std::uint64_t hi = 0;
    std::uint64_t lo = 0;
    for (Amount j = 1; j <= x; ++j) {
      const std::uint64_t price = slot_price_by_division(static_cast<std::uint64_t>(B - V + j * u), static_cast<std::uint64_t>(B));
      lo += price;
      if (lo < price) ++hi;
    }
    const auto g = ordered(2, {{0, 1, B, 3}});
    const std::vector<Amount> xs {x};
    const std::vector<Amount> free {V};
    const auto obj = detail::t2_objectives(*g, edges, xs, free, u, true);
    ASSERT_EQ(obj.headroom_penalty, (Counter128 {hi, lo})) << "B " << B << " V " << V << " u " << u << " x " << x;
    ASSERT_EQ(obj.route_cost, (Counter128 {0, static_cast<std::uint64_t>(3 * x)}));
    ASSERT_EQ(detail::t2_objectives(*g, edges, xs, free, u, false).headroom_penalty, Counter128 {});
  }
}

// The penalty past 2^64, with a carry into the high word, as a plan publishes it (no layout reaches
// it in a test's time: it takes 2^32 slots). One link B = V = 2^33, u = 1, x = 2^33: slot j pays
// floor(2^32 j / 2^33) = floor(j / 2), and the slots 2i - 1, 2i pay 2i - 1, so the sum is
// (2^32)^2 = 2^64 exactly. Three links B = V = 2^32, u = 1, x = 2^32: slot j pays j, each link
// 2^32 (2^32 + 1) / 2 = 2^63 + 2^31, the three 2^64 + 2^63 + 3 * 2^31; J = 3 * 2^32 (cost 1).
TEST(T2Objectives, PenaltyAboveTwoTo64) {
  const Amount p32 = Amount {1} << 32;
  const auto one = ordered(2, {{0, 1, 2 * p32, 1}});
  const std::vector<EdgeId> e0 {0};
  const auto a = detail::t2_objectives(*one, e0, std::vector<Amount> {2 * p32}, std::vector<Amount> {2 * p32}, 1, true);
  EXPECT_EQ(a.headroom_penalty, (Counter128 {1, 0}));
  EXPECT_EQ(a.route_cost, (Counter128 {0, std::uint64_t {1} << 33}));
  const auto three = ordered(4, {{0, 1, p32, 1}, {1, 2, p32, 1}, {2, 3, p32, 1}});
  const std::vector<EdgeId> e012 {0, 1, 2};
  const std::vector<Amount> x3 {p32, p32, p32};
  const auto b = detail::t2_objectives(*three, e012, x3, x3, 1, true);
  EXPECT_EQ(b.headroom_penalty, (Counter128 {1, (std::uint64_t {1} << 63) + 3 * (std::uint64_t {1} << 31)}));
  EXPECT_EQ(b.route_cost, (Counter128 {0, 3 * static_cast<std::uint64_t>(p32)}));
  // More slots than the free capacity holds is an internal contradiction.
  EXPECT_THROW((void)detail::t2_objectives(*one, e0, std::vector<Amount> {3}, std::vector<Amount> {2}, 1, true),
               std::logic_error);
}

// ---------------------------------------------------------------- up to L (7.7)

// Up to: G1 up to 6, a common size, a dead head, a tiny head: the chosen count,
// the routes, the lossless delivery of D; with every head live, F'(k) is EXACT(k)'s bound.
TEST(T2UpTo, ContractExamples) {
  struct Case {
    std::string id;
    std::int32_t n;
    std::vector<InEdge> rows;
    std::vector<NodeId> heads;
    NodeId dst;
    std::int64_t L;
    std::optional<Amount> D;
    std::int64_t k;
    std::vector<Edges> routes;
    Amount delivered;
    Amounts bounds;
  };
  std::vector<Edges> g1k6 {{0, 2}};
  g1k6.insert(g1k6.end(), 5, Edges {1, 3});
  const std::vector<Case> cases = {
      {"G1 up to 6", 4, kG1, {0}, 3, 6, 49152, 6, g1k6, 49152, {40960, 40960, 40960, 40960, 40960, 49152}},
      {"up to, a common size", 4, {{0, 3, 6, 1}, {1, 3, 3, 1}, {1, 3, 3, 1}}, {0, 1}, 3, 2, 12, 2, {{0}, {1}, {0}, {2}}, 12, {7, 12}},
      {"up to, a dead head", 3, {{0, 2, 0, 1}, {1, 2, 10, 1}}, {0, 1}, 2, 2, 10, 1, {{1}}, 0, {}},
      {"up to, a tiny head", 3, {{0, 2, 1, 1}, {1, 2, 100, 1}}, {0, 1}, 2, 3, {}, 1, {{0}, {1}}, 2, {2, 0, 0}},
  };
  for (const auto& c : cases) {
    SCOPED_TRACE(c.id);
    const auto g = ordered(c.n, c.rows);
    auto p = make(g, c.heads, {c.dst}, t2(c.L, LOSSLESS, true, true));
    const auto r = p->replace(c.D.has_value() ? c.D : INF);
    EXPECT_EQ(r.plan->count, c.k);
    EXPECT_EQ(routes(*r.program), c.routes);
    EXPECT_EQ(r.placed, c.delivered);
    for (std::size_t k = 1; k <= c.bounds.size(); ++k) {
      auto q = make(g, c.heads, {c.dst}, t2(static_cast<std::int64_t>(k), LOSSLESS, true));
      EXPECT_EQ(q->replace(INF).placed, c.bounds[k - 1]) << k;
    }
  }
}

// Three parallel links of 6, 5 and 4, D = 11, up to 3. k = 1: U = 6, F'(1) = 6. k = 2: U = 5 (one
// tunnel on the 6-link, one on the 5-link), the bounds 10 and 11 leave the score open, F'(2) = 11
// (z = 11 gives 6 and 5) = D: chosen. k = 3: U = 4, its lower bound min(11, 3 * 4) already reaches
// D, so a budget F'(2) one unit low would choose 3.
TEST(T2UpTo, ReachesTheDesignExactly) {
  const auto g = ordered(2, {{0, 1, 6, 1}, {0, 1, 5, 1}, {0, 1, 4, 1}});
  auto p = make(g, {0}, {1}, t2(3, LOSSLESS, false, true));
  const auto r = p->replace(Amount {11});
  EXPECT_EQ(r.plan->count, 2);
  EXPECT_EQ(r.plan->common_share, 5);
  EXPECT_EQ(r.plan->planned_per_tunnel, 5);
  EXPECT_EQ(routes(*r.program), (std::vector<Edges> {{0}, {1}}));
  EXPECT_EQ(r.placed, 11);
}

// G1, up to 6, Lossless replace(inf): count 6, placed 49152.
TEST(T2UpTo, Unlimited) {
  auto p = make(g1(), {0}, {3}, t2(6, LOSSLESS, true, true));
  const auto r = p->replace(INF);
  EXPECT_EQ(r.plan->count, 6);
  EXPECT_EQ(r.placed, 49152);
  EXPECT_EQ(r.plan->design_demand, std::nullopt);
  EXPECT_TRUE(r.plan->head_design.empty());
}

// ---------------------------------------------------------------- operations (U)

// An unlimited design (D = offer = +inf): u = U.
TEST(T2Operations, UnlimitedDesign) {
  auto p = make(g1(), {0}, {3}, t2(16, LOSSLESS, true));
  const auto r = p->replace(INF);
  EXPECT_EQ(r.plan->planned_per_tunnel, 2925);
  EXPECT_TRUE(r.plan->head_design.empty());
  EXPECT_EQ(r.plan->design_demand, std::nullopt);
  EXPECT_EQ(r.placed, 46812);
  EXPECT_EQ(r.termination, Termination::Exhausted);
}

// Every head dead: no tunnel, a plan of zeros (count N, 1 for "up to"), the offer lost at
// the heads (Lossy COMPLETE, Lossless EXHAUSTED); "up to" counts 1.
TEST(T2Operations, AllDead) {
  const auto g = twoh(0, 0);
  for (const bool up : {false, true}) {
    for (const auto a : {LOSSY, LOSSLESS}) {
      auto p = make(g, {0, 1}, {2}, t2(up ? 3 : 2, a, false, up));
      const auto r = p->replace(Amount {10}, Amount {10});
      const auto& pl = *r.plan;
      EXPECT_EQ(pl.count, up ? 1 : 2);
      EXPECT_EQ(pl.head_design, (std::vector<std::pair<NodeId, Amount>> {{0, 5}, {1, 5}}));
      EXPECT_EQ(pl.live_heads, std::vector<NodeId> {});
      EXPECT_EQ(pl.insufficient_resolution, false);
      EXPECT_EQ(pl.fluid_share, (std::pair<Amount, Amount> {0, 1}));
      EXPECT_EQ(pl.common_share, 0);
      EXPECT_EQ(pl.planned_per_tunnel, 0);
      EXPECT_EQ(pl.planned_total, 0);
      EXPECT_EQ(lo(pl.route_cost), 0);
      EXPECT_FALSE(pl.headroom_penalty.has_value());
      EXPECT_EQ(pl.forwarding_bound, 0);
      EXPECT_TRUE(r.program->tunnels().empty());
      EXPECT_EQ(r.placed, 0);
      EXPECT_TRUE(p->records().empty());
      if (a == LOSSY) {
        EXPECT_EQ(r.call_source_loss, 10);
        EXPECT_EQ(r.termination, Termination::Complete);
        EXPECT_EQ(source_rows(r), (std::vector<SourceRow> {{0, 5, 5, 0, 5, 0}, {1, 5, 5, 0, 5, 0}}));
      } else {
        EXPECT_EQ(r.termination, Termination::Exhausted);
      }
    }
  }
}

// A tiny head, forwarding without tunnels: U = 0 with live heads: no tunnel,
// insufficient_resolution; an installed program keeping the roster that forward reuses.
TEST(T2Operations, InsufficientResolution) {
  const auto g = twoh(1, 100);
  for (const auto a : {LOSSY, LOSSLESS}) {
    auto p = make(g, {0, 1}, {2}, t2(16, a));
    const auto r = p->replace(Amount {200}, Amount {200});
    EXPECT_EQ(r.plan->insufficient_resolution, true);
    EXPECT_EQ(r.plan->fluid_share, (std::pair<Amount, Amount> {1, 16}));
    EXPECT_EQ(r.plan->common_share, 0);
    EXPECT_EQ(r.plan->live_heads, (std::vector<NodeId> {0, 1}));
    EXPECT_TRUE(r.program->tunnels().empty());
    ASSERT_NE(p->program(), nullptr);
    EXPECT_EQ(r.placed, 0);
    const auto f = p->forward(Amount {10});
    if (a == LOSSY) {
      EXPECT_EQ(source_rows(r), (std::vector<SourceRow> {{0, 100, 100, 0, 100, 0}, {1, 100, 100, 0, 100, 0}}));
      EXPECT_EQ(f.call_source_loss, 10);
      EXPECT_EQ(f.termination, Termination::Complete);
      EXPECT_EQ(source_rows(f), (std::vector<SourceRow> {{0, 5, 5, 0, 5, 0}, {1, 5, 5, 0, 5, 0}}));
    } else {
      EXPECT_EQ(r.termination, Termination::Exhausted);
      EXPECT_EQ(f.placed, 0);
      EXPECT_EQ(f.termination, Termination::Exhausted);
    }
    EXPECT_EQ(f.plan, r.plan);
  }
}

// A dead head: indices 1 and 3; the rows; two records (entries 1, 3).
TEST(T2Operations, DeadHead) {
  auto p = make(twoh(0, 10), {0, 1}, {2}, t2(2, LOSSY));
  const auto r = p->replace(Amount {10}, Amount {10});
  EXPECT_EQ(indices(*r.program), (std::vector<std::int32_t> {1, 3}));
  EXPECT_EQ(source_rows(r), (std::vector<SourceRow> {{0, 5, 5, 0, 5, 0}, {1, 5, 5, 5, 0, 0}}));
  std::vector<std::optional<std::int32_t>> entries;
  for (const auto& rec : p->records()) entries.push_back(rec.entry);
  std::sort(entries.begin(), entries.end());
  EXPECT_EQ(entries, (std::vector<std::optional<std::int32_t>> {1, 3}));
}

// OPTIMIZED with one source is the EQUAL calculation; equal sources.
TEST(T2Operations, OptimizedOneSource) {
  const auto g = g1();
  auto a = make(g, {0}, {3}, t2(16, LOSSY, true));
  auto b = make(g, {0}, {3}, t2(16, LOSSY, true, false, Origin::Optimized));
  const auto ra = a->replace(Amount {49152});
  const auto rb = b->replace(Amount {49152});
  EXPECT_EQ(ra.plan, rb.plan);
  EXPECT_EQ(routes(*ra.program), routes(*rb.program));
  EXPECT_EQ(ra.by_tunnel, rb.by_tunnel);
  EXPECT_EQ(ra.call_drops, rb.call_drops);
  EXPECT_EQ(b->config().origin, Origin::Optimized);
  EXPECT_FALSE(rb.program->source_weights().has_value());
}

// Heads behind 2 and 10 against the global optimum, N = 2: common size 1, lossless bound 4;
// optimized Global delivers 12 on the same network.
TEST(T2Operations, TwoHeadsEqualAgainstGlobal) {
  const auto g = twoh(2, 10);
  auto p = make(g, {0, 1}, {2}, t2(2));
  const auto r = p->replace(INF);
  EXPECT_EQ(r.plan->common_share, 1);
  EXPECT_EQ(r.placed, 4);
  EXPECT_EQ(capacity(*g, v({0, 1}), v({2}), QueryOptions {global_config()}).placed, 12);
}

// An empty roster: NO_ROUTE, an empty program, a plan of zeros.
TEST(T2Operations, EmptyRoster) {
  const auto g = ordered(3, {{1, 2, 5, 1}});
  for (const auto a : {LOSSY, LOSSLESS}) {
    auto p = make(g, {0}, {2}, t2(2, a, true));
    const auto r = p->replace(Amount {5}, Amount {5});
    EXPECT_EQ(r.termination, Termination::NoRoute);
    EXPECT_TRUE(r.program->tunnels().empty());
    EXPECT_EQ(r.plan->count, 2);
    EXPECT_TRUE(r.plan->head_design.empty());
    EXPECT_EQ(r.plan->live_heads, std::vector<NodeId> {});
    EXPECT_EQ(lo(r.plan->headroom_penalty), 0);
    EXPECT_EQ(r.plan->forwarding_bound, 0);
  }
}

// A zero design: D = 0 plans at u = 1.
TEST(T2Operations, ZeroDesign) {
  auto p = make(one10(), {0}, {1}, t2(2));
  const auto r = p->replace(Amount {0});
  EXPECT_EQ(r.plan->planned_per_tunnel, 1);
  EXPECT_EQ(r.program->tunnels().size(), 2U);
  for (const auto& t : r.program->tunnels()) EXPECT_EQ(t.planned, 1);
  EXPECT_EQ(r.placed, 0);
  EXPECT_EQ(r.termination, Termination::Complete);
}

// Forward keeps the program and the plan, sets no flag.
TEST(T2Operations, Forward) {
  const auto g = g1();
  auto p = make(g, {0}, {3}, t2(16, LOSSY, true));
  const auto first = p->replace(Amount {49152});
  const auto f = p->forward(Amount {32768});
  EXPECT_EQ(f.placed, 32768);
  EXPECT_FALSE(f.routes_changed || f.weights_changed || f.reservations_changed);
  EXPECT_EQ(f.program, first.program);
  EXPECT_EQ(f.plan, first.plan);
  auto q = make(g, {0}, {3}, t2(16, LOSSLESS, true));
  (void)q->replace(Amount {49152});
  EXPECT_EQ(q->forward(INF).placed, 46812);
}

// The same routes and weights with another planned size: no flag (a planned size is no
// reservation).
TEST(T2Operations, Flags) {
  auto p = make(g1(), {0}, {3}, t2(16, LOSSY, true));
  const auto a = p->replace(Amount {49152});
  const auto b = p->replace(Amount {49152}, Amount {8192});
  EXPECT_FALSE(b.routes_changed || b.weights_changed || b.reservations_changed);
  EXPECT_EQ(a.plan->planned_per_tunnel, 2925);
  EXPECT_EQ(b.plan->planned_per_tunnel, 512);
}

// ---------------------------------------------------------------- validation (V)

// T2's row of the matrix, in order, at the policy and at an evaluator.
TEST(TunnelValidation, T2Matrix) {
  const auto g = twoh(10, 10);
  FlowLedger fg(g);
  const Global tun {.tunnels = Tunnels {2}};
  auto rejects = [&](FlowPolicyConfig c, const std::string& text, DagCache* cache = nullptr) {
    EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0, 1}), v({2}), c, nullptr, cache); }), text);
    EXPECT_TRUE(detail::LedgerAccess::bound_owners(fg).empty());
    if (cache == nullptr) EXPECT_EQ(error_of([&] { (void)capacity(*g, v({0, 1}), v({2}), QueryOptions {c}); }), text);
  };
  auto with = [&](auto&& set) {
    FlowPolicyConfig c {.routing = tun};
    set(c);
    return c;
  };
  rejects(with([](auto& c) { c.placement = FlowPlacement::EqualBalanced; }), "T2 does not accept placement");
  rejects(with([](auto& c) { c.weights = WeightBasis::Capacity; }), "T2 does not accept weights");
  rejects(with([](auto& c) { c.edge_weights = std::make_shared<const std::vector<Amount>>(Amounts {1, 1}); }),
          "T2 does not accept edge_weights");
  // Every Global rejects any selection, the default and HEADROOM included, before support and
  // tunnels (and so before placement); prefer_headroom without tunnels right after them.
  for (const auto& sel : {PathSelection(), PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Headroom)}) {
    rejects(with([&](auto& c) { c.selection = sel; }), "Global does not accept selection");
    rejects(FlowPolicyConfig {.routing = Global {}, .selection = sel}, "Global does not accept selection");
  }
  const auto dag = dag_of(*g, {0, 1});
  const Global both {.support = dag, .tunnels = Tunnels {2}, .prefer_headroom = true};
  rejects(FlowPolicyConfig {.routing = both, .selection = PathSelection()}, "Global does not accept selection");
  rejects(FlowPolicyConfig {.routing = both}, "Global support and tunnels are mutually exclusive");
  rejects(FlowPolicyConfig {.routing = both, .placement = FlowPlacement::EqualBalanced},
          "Global support and tunnels are mutually exclusive");
  rejects(FlowPolicyConfig {.routing = Global {.prefer_headroom = true}}, "prefer_headroom requires tunnels");
  rejects(FlowPolicyConfig {.routing = Global {.support = dag, .prefer_headroom = true}, .placement = FlowPlacement::EqualBalanced},
          "prefer_headroom requires tunnels");
  rejects(with([](auto& c) {
            c.placement = FlowPlacement::EqualBalanced;
            c.selection = PathSelection();
          }),
          "Global does not accept selection");
  rejects(with([](auto& c) { c.origin = Origin::Optimized; }), "T2 with several sources requires origin=EQUAL");
  DagCache cache(*g);
  rejects(FlowPolicyConfig {.routing = tun}, "T2 does not use DagCache", &cache);
  // The several-sources rule counts the input sources before the mask; one source with OPTIMIZED is accepted.
  const auto mask = g->eligibility(std::vector<NodeId> {1}, {});
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0, 1}), v({2}), with([](auto& c) { c.origin = Origin::Optimized; }), &mask); }),
            "T2 with several sources requires origin=EQUAL");
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0, 0}), v({2}), with([](auto& c) { c.origin = Origin::Optimized; })); }), "");
  // Accepted: both admissions, with and without prefer_headroom; the snapshot resolves the origin
  // and keeps no selection.
  for (const auto a : {LOSSLESS, LOSSY}) {
    for (const bool h : {false, true}) EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0, 1}), v({2}), t2(2, a, h)); }), "");
  }
  FlowPolicy p(fg, v({0, 1}), v({2}), FlowPolicyConfig {.routing = tun});
  EXPECT_EQ(p.config().origin, Origin::Equal);
  EXPECT_EQ(p.config().selection, std::nullopt);
  EXPECT_EQ(std::string(detail::model_name(p.config().routing)), "T2");
}

// The count range for N and L; the tunnel limit with L for "up to".
TEST(TunnelValidation, T2CountAndLimit) {
  const auto g = one10();
  FlowLedger fg(g);
  for (const std::int64_t n : {0, 10001, -1}) {
    for (const bool up : {false, true}) {
      EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), t2(n, LOSSLESS, false, up)); }),
                "tunnel count must be in [1, 10000]");
    }
  }
  const auto fan = [](std::int32_t h) {
    std::vector<InEdge> rows;
    rows.reserve(static_cast<std::size_t>(h));
    for (std::int32_t i = 0; i < h; ++i) rows.push_back({i, h, 1, 1});
    return shared(make_graph(h + 1, rows));
  };
  const auto heads = [](std::int32_t h) {
    std::vector<NodeId> out(static_cast<std::size_t>(h));
    for (std::int32_t i = 0; i < h; ++i) out[static_cast<std::size_t>(i)] = i;
    return out;
  };
  const auto big = fan(107375);
  FlowLedger fb(big);
  for (const bool up : {false, true}) {
    EXPECT_EQ(error_of([&] { FlowPolicy p(fb, heads(107375), v({107375}), t2(10000, LOSSLESS, false, up)); }),
              "tunnel limit exceeded: 107375 heads x 10000 tunnels > 1073741823");
  }
  const auto ok = fan(107374);
  FlowLedger fo(ok);
  EXPECT_EQ(error_of([&] { FlowPolicy p(fo, heads(107374), v({107374}), t2(10000)); }), "");
}

// T2's operations and the design's domain, with their texts.
TEST(TunnelOperations, T2CallValidation) {
  auto p = make(one10(), {0}, {1}, t2(2));
  for (const Amount x : {5, 0}) EXPECT_EQ(error_of([&] { (void)p->place(x); }), "T2 does not support incremental placement");
  EXPECT_EQ(error_of([&] { (void)p->forward(Amount {5}); }), "forward requires an installed program");
  for (const Amount d : {Amount {-1}, kL62}) {
    EXPECT_EQ(error_of([&] { (void)p->replace(Amount {5}, d); }),
              "design_demand must be an integer in [0, 2^62) or permitted +inf");
  }
  EXPECT_EQ(error_of([&] { (void)p->replace(Amount {-1}, Amount {-1}); }),
            "offer must be an integer in [0, 2^62) or permitted +inf");
  auto ly = make(one10(), {0}, {1}, t2(2, LOSSY));
  EXPECT_EQ(error_of([&] { (void)ly->replace(INF); }), "LOSSY forwarding requires a finite offer");
  EXPECT_EQ(error_of([&] { (void)p->replace(INF); }), "");
  EXPECT_EQ(error_of([&] { (void)p->evaluate(INF, UNLIMITED); }), "");
}

// No certificates; the plan only for the tunnel models.
TEST(TunnelOperations, T2NoCertificates) {
  const auto g = g1();
  auto p = make(g, {0}, {3}, t2(16));
  const auto r = p->replace(INF, std::nullopt, OutputOptions {.min_cut = true});
  EXPECT_FALSE(r.total_cost.has_value());
  EXPECT_FALSE(r.min_cut.has_value());
  EXPECT_FALSE(r.cost_distribution.has_value());
  auto pin = make(g, {0}, {3}, FlowPolicyConfig {.routing = Pinned {r.program, nullptr}});
  EXPECT_FALSE(pin->replace(Amount {10}).plan.has_value());
  // An optimized-origin T2 policy whose roster is empty (source 0 reaches nothing) reports no
  // certificate either, never an empty cut or curve: replace(5) is NO_ROUTE, replace(0) COMPLETE.
  auto empty = make(ordered(3, {{1, 2, 10, 1}}), {0}, {2}, t2(2, LOSSLESS, false, false, Origin::Optimized));
  for (const Amount offer : {Amount {5}, Amount {0}}) {
    SCOPED_TRACE(offer);
    const auto e = empty->replace(offer, std::nullopt, OutputOptions {.min_cut = true});
    EXPECT_EQ(e.termination, offer > 0 ? Termination::NoRoute : Termination::Complete);
    EXPECT_EQ(e.placed, 0);
    EXPECT_FALSE(e.total_cost.has_value());
    EXPECT_FALSE(e.min_cut.has_value());
    EXPECT_FALSE(e.cost_distribution.has_value());
    ASSERT_TRUE(e.plan.has_value());
  }
}

// ---------------------------------------------------------------- programs (P)

// The 8.5 call: the 7.9 call's program replayed with e2 masked.
TEST(TunnelPrograms, T2Replay) {
  const auto g = g1();
  auto p = make(g, {0}, {3}, t2(16, LOSSY, true));
  const auto base = p->replace(Amount {49152}).program;
  const auto fail = g->eligibility({}, std::vector<EdgeId> {2});
  auto replay = make(g, {0}, {3}, FlowPolicyConfig {.routing = Pinned {base, nullptr}, .admission = LOSSY}, &fail);
  const auto r = replay->replace(Amount {49152});
  auto st = row_statuses(r);
  st.resize(3);
  EXPECT_EQ(st, (std::vector<TunnelStatus> {FAILED, FAILED, UP}));
  EXPECT_EQ(r.placed, 40960);
  EXPECT_EQ(r.call_drops, (Pairs {{1, 8192}}));
}

// The same program stale: e2 a black hole.
TEST(TunnelPrograms, T2Stale) {
  const auto g = g1();
  auto p = make(g, {0}, {3}, t2(16, LOSSY, true));
  const auto base = p->replace(Amount {49152}).program;
  const auto dark = std::make_shared<const Eligibility>(g->eligibility({}, std::vector<EdgeId> {2}));
  auto stale = make(g, {0}, {3}, FlowPolicyConfig {.routing = Pinned {base, dark}, .admission = LOSSY});
  const auto r = stale->replace(Amount {49152});
  std::vector<TunnelStatus> want(2, BLACKHOLED);
  want.insert(want.end(), 14, UP);
  EXPECT_EQ(row_statuses(r), want);
  EXPECT_EQ(r.call_drops, (Pairs {{1, 2048}, {2, 6144}}));
  EXPECT_EQ(r.placed, 40960);
  for (const auto& t : base->tunnels()) EXPECT_EQ(t.status, UP);
}

// The 7.9 call's tunnels: requested 0, reserved 0, planned 2925, UP, weight 1; equal sources.
TEST(TunnelPrograms, T2Fields) {
  auto p = make(g1(), {0}, {3}, t2(16, LOSSY, true));
  const auto r = p->replace(Amount {49152});
  ASSERT_EQ(r.program->tunnels().size(), 16U);
  for (const auto& t : r.program->tunnels()) {
    EXPECT_EQ(t.requested, 0);
    EXPECT_EQ(t.reserved, 0);
    EXPECT_EQ(t.planned, 2925);
    EXPECT_EQ(t.status, UP);
    EXPECT_EQ(t.weight, 1);
    EXPECT_EQ(t.destination, NodeId {3});
    EXPECT_EQ(t.route_cost, Cost {2});
  }
  EXPECT_FALSE(r.program->source_weights().has_value());
}

// The planned bandwidth is T2's alone: 0 on T1's tunnels, on a controller's paths (Global,
// CostOnly MaxFlow) and on a supplied program's tunnels; Pinned over a T2 program keeps u.
TEST(TunnelPrograms, PlannedOnlyOnT2) {
  const auto g = twoh(2, 10);
  auto planned = [](const PlacementResult& r) {
    Amounts out;
    for (const auto& t : r.program->tunnels()) out.push_back(t.planned);
    return out;
  };
  auto t1 = make(g, {0, 1}, {2}, FlowPolicyConfig {.routing = CapacityAware {.tunnels = Tunnels {2}}, .admission = LOSSY});
  EXPECT_EQ(planned(t1->replace(12, 12)), (Amounts {0, 0, 0, 0}));
  auto global = make(g, {0, 1}, {2}, FlowPolicyConfig {.routing = Global {}});
  EXPECT_EQ(planned(global->place(12)), (Amounts {0, 0}));
  auto max_flow = make(g, {0, 1}, {2}, FlowPolicyConfig {.routing = CostOnly {}, .placement = FlowPlacement::MaxFlow});
  EXPECT_EQ(planned(max_flow->place(12)), (Amounts {0, 0}));
  const std::vector<NodeId> s {0, 1};
  const std::vector<NodeId> d {2};
  const std::vector<Route> routes {{0, std::make_shared<const Path>(*g, std::vector<std::vector<EdgeId>> {{0}}), 1},
                                   {1, std::make_shared<const Path>(*g, std::vector<std::vector<EdgeId>> {{1}}), 1}};
  const auto supplied = ForwardingProgram::from_routes(*g, s, d, routes);
  auto pinned = make(g, {0, 1}, {2}, FlowPolicyConfig {.routing = Pinned {.program = supplied}, .admission = LOSSY});
  EXPECT_EQ(planned(pinned->replace(Amount {4})), (Amounts {0, 0}));
  auto tun = make(g, {0, 1}, {2}, t2(2));
  const auto laid = tun->replace(INF);
  EXPECT_EQ(planned(laid), (Amounts {1, 1, 1, 1}));
  auto replay = make(g, {0, 1}, {2}, FlowPolicyConfig {.routing = Pinned {.program = laid.program}, .admission = LOSSY});
  EXPECT_EQ(planned(replay->replace(Amount {4})), (Amounts {1, 1, 1, 1}));
}

// ---------------------------------------------------------------- evaluators (Q)

// The design rule of 12.1; nothing shared between pairs.
TEST(TunnelEvaluators, T2Capacity) {
  const auto g = g1();
  auto err = [&](const QueryOptions& o) { return error_of([&] { (void)capacity(*g, v({0}), v({3}), o); }); };
  QueryOptions q {t2(16, LOSSLESS, true)};
  EXPECT_EQ(capacity(*g, v({0}), v({3}), q).placed, 46812);
  q.design_demand = Amount {8192};
  EXPECT_EQ(capacity(*g, v({0}), v({3}), q).placed, 46812);
  QueryOptions ly {t2(16, LOSSY, true)};
  ly.design_demand = Amount {49152};
  EXPECT_EQ(capacity(*g, v({0}), v({3}), ly).placed, 47104);
  QueryOptions lo_ {t2(16, LOSSY, true)};
  lo_.offered = Amount {49152};
  EXPECT_EQ(capacity(*g, v({0}), v({3}), lo_).placed, 47104);
  EXPECT_EQ(err(QueryOptions {t2(16, LOSSY, true)}), "T2 LOSSY capacity requires finite offered");
  QueryOptions unlimited {t2(16, LOSSY, true)};
  unlimited.design_demand = std::optional<Amount> {};
  EXPECT_EQ(err(unlimited), "T2 LOSSY capacity requires finite offered");
  QueryOptions offered {t2(16, LOSSLESS, true)};
  offered.offered = Amount {5};
  EXPECT_EQ(err(offered), "T2 LOSSLESS capacity does not accept offered");
  const auto one = one10();
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{v({0}), v({1})}, {v({0}), v({1})}};
  for (const std::int32_t workers : {1, 2}) {
    const auto rs = batch_capacity(*one, pairs, QueryOptions {t2(2)}, {}, workers);
    ASSERT_EQ(rs.size(), 2U);
    EXPECT_EQ(rs[0].placed, 10);
    EXPECT_EQ(rs[1].placed, 10);
  }
}

// T2 LOSSY's offer rule (12.1) in all three evaluators: a finite offer is used as given, an absent
// one is the design if finite, and an explicit +inf is an error whatever the design. One link of 10,
// N = 2: both tunnels share the link, so the Lossy delivery is min(offer, 10).
TEST(TunnelEvaluators, T2LossyOfferRule) {
  const auto g = one10();
  const std::optional<std::optional<Amount>> absent;
  const std::optional<std::optional<Amount>> inf {std::in_place, std::nullopt};
  const auto amount = [](Amount a) { return std::optional<std::optional<Amount>> {std::in_place, a}; };
  const std::string text = "T2 LOSSY capacity requires finite offered";
  struct Row {
    std::optional<std::optional<Amount>> design, offered;
    std::optional<Amount> placed;  // nullopt: the error
  };
  const std::vector<Row> rows = {
      {absent, absent, {}},     {amount(0), absent, 0},   {amount(4), absent, 4},   {inf, absent, {}},
      {absent, inf, {}},        {amount(0), inf, {}},     {amount(4), inf, {}},     {inf, inf, {}},
      {absent, amount(0), 0},   {amount(4), amount(0), 0}, {inf, amount(0), 0},     {absent, amount(7), 7},
      {amount(0), amount(7), 7}, {amount(4), amount(7), 7}, {inf, amount(7), 7},
  };
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{v({0}), v({1})}, {v({0}), v({1})}};
  for (const auto& row : rows) {
    QueryOptions q {t2(2, LOSSY)};
    q.design_demand = row.design;
    q.offered = row.offered;
    const std::string where = "design " + (row.design ? (*row.design ? std::to_string(**row.design) : "inf") : "absent") +
                              ", offered " +
                              (row.offered ? (*row.offered ? std::to_string(**row.offered) : "inf") : "absent");
    SCOPED_TRACE(where);
    if (!row.placed) {
      EXPECT_EQ(error_of([&] { (void)capacity(*g, v({0}), v({1}), q); }), text);
      EXPECT_EQ(error_of([&] { (void)batch_capacity(*g, pairs, q, {}, 2); }), text);
      EXPECT_EQ(error_of([&] { (void)sensitivity(*g, v({0}), v({1}), q, 2); }), text);
      continue;
    }
    EXPECT_EQ(capacity(*g, v({0}), v({1}), q).placed, *row.placed);
    const auto rs = batch_capacity(*g, pairs, q, {}, 2);
    ASSERT_EQ(rs.size(), 2U);
    EXPECT_EQ(rs[0].placed, *row.placed);
    EXPECT_EQ(rs[1].placed, *row.placed);
    EXPECT_EQ(sensitivity(*g, v({0}), v({1}), q, 2).baseline.placed, *row.placed);
  }
}

// Deleting dead head 0's zero link removes it from the roster (loss -10);
// every allowed link is a candidate.
TEST(TunnelEvaluators, T2Sensitivity) {
  QueryOptions q {t2(2)};
  q.design_demand = Amount {10};
  for (const std::int32_t workers : {1, 2}) {
    const auto s = sensitivity(*twoh(0, 10), v({0, 1}), v({2}), q, workers);
    EXPECT_EQ(s.baseline.placed, 0);
    std::vector<std::pair<EdgeId, Amount>> got;
    for (const auto& e : s.entries) got.emplace_back(e.edge, e.loss);
    EXPECT_EQ(got, (std::vector<std::pair<EdgeId, Amount>> {{0, -10}, {1, 0}}));
  }
  const auto g = ordered(3, {{0, 1, 10, 1}, {2, 1, 5, 1}});
  const auto s = sensitivity(*g, v({0}), v({1}), QueryOptions {t2(2)}, 1);
  std::vector<std::pair<EdgeId, Amount>> got;
  for (const auto& e : s.entries) got.emplace_back(e.edge, e.loss);
  EXPECT_EQ(got, (std::vector<std::pair<EdgeId, Amount>> {{0, 10}, {1, 0}}));
}

// ---------------------------------------------------------------- numeric ranges (W)

// ONEMAX, N = 1, prefer_headroom, unlimited: U = u = 2^62 - 1, Phi = 2^32.
TEST(TunnelRanges, T2OneMax) {
  const Amount max = kL62 - 1;
  auto p = make(ordered(2, {{0, 1, max, 1}}), {0}, {1}, t2(1, LOSSLESS, true));
  const auto r = p->replace(INF);
  EXPECT_EQ(r.plan->common_share, max);
  EXPECT_EQ(r.plan->planned_total, max);
  EXPECT_EQ(lo(r.plan->headroom_penalty), Amount {1} << 32);
  EXPECT_EQ(r.plan->forwarding_bound, max);
  EXPECT_EQ(r.placed, max);
  EXPECT_EQ(r.plan->fluid_share, (std::pair<Amount, Amount> {max, 1}));
}

// FAN4MAX, N = 16, unlimited: U = 2^56 - 1, planned total 2^62 - 64, F = 2^62 - 4.
TEST(TunnelRanges, T2ManyTunnels) {
  const Amount c = (Amount {1} << 60) - 1;
  auto p = make(ordered(5, {{0, 4, c, 1}, {1, 4, c, 1}, {2, 4, c, 1}, {3, 4, c, 1}}), {0, 1, 2, 3}, {4}, t2(16));
  const auto r = p->replace(INF);
  EXPECT_EQ(r.plan->common_share, (Amount {1} << 56) - 1);
  EXPECT_EQ(r.plan->planned_total, kL62 - 64);
  EXPECT_EQ(r.plan->forwarding_bound, kL62 - 4);
  EXPECT_EQ(r.plan->fluid_share, (std::pair<Amount, Amount> {c, 16}));
}

// DEAD100: the live budgets F'(1) = F'(2) = 101 2^61 + 100 above 2^63, equal scores
// 2^61: count 1; F = 0; EXACT(1) Lossy replace(101) delivers 1 and loses 100 at the dead heads.
TEST(TunnelRanges, T2DeadBudget) {
  std::vector<InEdge> rows;
  for (std::int32_t h = 0; h < 100; ++h) rows.push_back({h, 101, 0, 1});
  rows.push_back({100, 101, Amount {1} << 61, 1});
  const auto g = ordered(102, rows);
  std::vector<NodeId> heads(101);
  for (std::int32_t h = 0; h < 101; ++h) heads[static_cast<std::size_t>(h)] = h;
  auto p = make(g, heads, {101}, t2(2, LOSSLESS, false, true));
  const auto r = p->replace(INF);
  EXPECT_EQ(r.plan->count, 1);
  EXPECT_EQ(r.plan->forwarding_bound, 0);
  EXPECT_EQ(r.plan->live_heads, (std::vector<NodeId> {100}));
  auto q = make(g, heads, {101}, t2(1, LOSSY));
  const auto s = q->replace(Amount {101});
  EXPECT_EQ(s.placed, 1);
  EXPECT_EQ(s.call_source_loss, 100);
}

// COSTMAX (one link of 10000 at cost 2^62 - 1), N = 10000, prefer_headroom, D = 10000:
// U = u = 1, J = 10000 (2^62 - 1), Phi = sum over j of floor(2^32 j / 10000).
TEST(TunnelRanges, T2CostMax) {
  auto p = make(ordered(2, {{0, 1, 10000, kL62 - 1}}), {0}, {1}, t2(10000, LOSSLESS, true));
  const auto r = p->replace(INF, Amount {10000});
  EXPECT_EQ(r.plan->common_share, 1);
  EXPECT_EQ(r.plan->planned_per_tunnel, 1);
  EXPECT_EQ(wide(*r.plan->route_cost), Int128(10000) * Int128(kL62 - 1));
  EXPECT_EQ(lo(r.plan->headroom_penalty), 21476983958656);
}

// ---------------------------------------------------------------- zero offers and masks (E)

// All dead: zero offers COMPLETE; a positive offer Lossy COMPLETE, Lossless EXHAUSTED.
TEST(TunnelZero, T2AllDead) {
  for (const auto a : {LOSSY, LOSSLESS}) {
    auto p = make(twoh(0, 0), {0, 1}, {2}, t2(2, a));
    EXPECT_EQ(p->replace(Amount {0}, Amount {10}).termination, Termination::Complete);
    EXPECT_EQ(p->evaluate(Amount {0}, Amount {10}).termination, Termination::Complete);
    EXPECT_EQ(p->forward(Amount {0}).termination, Termination::Complete);
    const auto r = p->replace(Amount {10}, Amount {10});
    EXPECT_EQ(r.termination, a == LOSSY ? Termination::Complete : Termination::Exhausted);
    EXPECT_TRUE(r.program->tunnels().empty());
    EXPECT_EQ(vec(r.program->sources()), v({0, 1}));
    EXPECT_TRUE(p->records().empty());
  }
}

// U = 0: replace(0) and forward(0) COMPLETE, forward(10) loses 5 at each head.
TEST(TunnelZero, T2NoTunnels) {
  auto p = make(twoh(1, 100), {0, 1}, {2}, t2(16, LOSSY));
  EXPECT_EQ(p->replace(Amount {0}).termination, Termination::Complete);
  EXPECT_EQ(p->forward(Amount {0}).termination, Termination::Complete);
  const auto f = p->forward(Amount {10});
  EXPECT_EQ(source_rows(f), (std::vector<SourceRow> {{0, 5, 5, 0, 5, 0}, {1, 5, 5, 0, 5, 0}}));
}

// ONE10, N = 2: zero calls COMPLETE; two tunnels planned 1; two zero records; evaluate
// publishes nothing.
TEST(TunnelZero, T2) {
  for (const auto a : {LOSSY, LOSSLESS}) {
    auto p = make(one10(), {0}, {1}, t2(2, a));
    const auto r = p->replace(Amount {0});
    EXPECT_EQ(r.termination, Termination::Complete);
    ASSERT_EQ(r.program->tunnels().size(), 2U);
    std::vector<FlowId> before;
    for (const auto& rec : p->records()) {
      EXPECT_EQ(rec.delivered, 0);
      before.push_back(rec.flow_id);
    }
    ASSERT_EQ(before.size(), 2U);
    EXPECT_EQ(p->evaluate(Amount {0}).termination, Termination::Complete);
    std::vector<FlowId> after;
    for (const auto& rec : p->records()) after.push_back(rec.flow_id);
    EXPECT_EQ(after, before);
    EXPECT_EQ(p->forward(Amount {0}).termination, Termination::Complete);
    EXPECT_EQ(p->records().size(), 2U);
  }
}

// The 7.9 call's program with e0 and e1 masked: no head reaches 3: NO_ROUTE, then COMPLETE;
// sixteen zero records.
TEST(TunnelZero, PinnedT2NoRoster) {
  const auto g = g1();
  auto p = make(g, {0}, {3}, t2(16, LOSSY, true));
  const auto base = p->replace(Amount {49152}).program;
  const auto mask = g->eligibility({}, std::vector<EdgeId> {0, 1});
  auto pin = make(g, {0}, {3}, FlowPolicyConfig {.routing = Pinned {base, nullptr}, .admission = LOSSY}, &mask);
  const auto r = pin->replace(Amount {49152});
  EXPECT_EQ(r.termination, Termination::NoRoute);
  EXPECT_EQ(pin->records().size(), 16U);
  EXPECT_EQ(row_statuses(r), std::vector<TunnelStatus>(16, FAILED));
  EXPECT_EQ(pin->replace(Amount {0}).termination, Termination::Complete);
}

// A zero offered; an Eligibility that excludes every edge.
TEST(TunnelZero, T2Evaluators) {
  const auto g = g1();
  QueryOptions zero {t2(16, LOSSY, true)};
  zero.offered = Amount {0};
  const auto r = capacity(*g, v({0}), v({3}), zero);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_EQ(r.placed, 0);
  for (const auto& e : sensitivity(*g, v({0}), v({3}), zero, 1).entries) EXPECT_EQ(e.loss, 0);
  const auto t = twoh(10, 10);
  const auto none = t->eligibility({}, std::vector<EdgeId> {0, 1});
  const auto c = capacity(*t, v({0, 1}), v({2}), QueryOptions {t2(2), &none});
  EXPECT_EQ(c.termination, Termination::NoRoute);
  EXPECT_EQ(c.placed, 0);
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{v({0, 1}), v({2})}, {v({0}), v({2})}};
  for (const auto& x : batch_capacity(*t, pairs, QueryOptions {t2(2), &none}, {}, 1)) EXPECT_EQ(x.termination, Termination::NoRoute);
}

// ---------------------------------------------------------------- the property test

// Random small graphs (nonuniform costs, masks, zero-capacity links, preloaded residuals, one or two
// heads): U, s* and the minimum (J, Phi) against the brute force; the program's objectives, k
// tunnels per live head in round order with dead holes; F against the independent lossless bound;
// 7.6's remainder bound; "up to" against an unpruned selection over the EXACT(k) programs at U(k).
TEST(T2Property, AgainstTheReference) {
  std::mt19937 rng(20261019);
  int checked = 0;
  int optimal = 0;
  int selections = 0;
  for (int trial = 0; trial < 250; ++trial) {
    SCOPED_TRACE(trial);
    const int n = 3 + static_cast<int>(rng() % 3);
    std::vector<InEdge> rows;
    const Amount capacities[] = {0, 1, 2, 3, 4, 6, 8};
    for (int i = 0; i + 1 < n; ++i) rows.push_back({i, i + 1, capacities[rng() % 7], 1 + static_cast<Cost>(rng() % 4)});
    for (int x = 0, extra = 1 + static_cast<int>(rng() % 4); x < extra; ++x) {
      const int a = static_cast<int>(rng() % static_cast<unsigned>(n - 1));
      const int b = 1 + static_cast<int>(rng() % static_cast<unsigned>(n - 1));
      if (a != b) rows.push_back({a, b, capacities[rng() % 7], 1 + static_cast<Cost>(rng() % 4)});
    }
    const auto g = shared(make_graph(n, rows));
    std::vector<NodeId> heads {0};
    if (n > 3 && rng() % 2 == 0) heads.push_back(1);
    const std::vector<NodeId> dests {n - 1};
    std::optional<Eligibility> mask;
    if (rng() % 4 == 0) mask = g->eligibility({}, std::vector<EdgeId> {static_cast<EdgeId>(rng() % rows.size())});
    const Eligibility* elig = mask.has_value() ? &*mask : nullptr;
    Amounts V(static_cast<std::size_t>(g->num_edges()));
    for (std::size_t e = 0; e < V.size(); ++e) {
      const Amount c = g->capacity_units()[e];
      V[e] = rng() % 3 == 0 ? static_cast<Amount>(rng() % static_cast<unsigned>(c + 1)) : c;
    }
    const bool hr = rng() % 2 == 0;
    const std::int64_t N = 1 + static_cast<std::int64_t>(rng() % 3);
    const std::optional<Amount> D = rng() % 3 == 0 ? std::nullopt : std::optional<Amount>(static_cast<Amount>(rng() % 15));
    const auto design = D.has_value() ? std::optional<std::optional<Amount>>(D) : UNLIMITED;
    auto p = make(g, heads, dests, t2(N, LOSSLESS, hr), elig, std::make_shared<FlowLedger>(g, V));
    const auto r = p->replace(INF, design);
    const auto& pl = *r.plan;
    const auto& roster = detail::PolicyAccess::engine(*p.p).roster;
    const auto H = static_cast<std::int64_t>(roster.size());
    if (H == 0) continue;
    const auto u_of = [&](Amount U) {
      if (!D.has_value()) return U;
      Amount need = 1;
      for (NodeId h : pl.live_heads.value_or(std::vector<NodeId> {})) {
        const auto i = std::find(roster.begin(), roster.end(), h) - roster.begin();
        const Amount share = t1_reference::rr(*D, H, i);
        need = std::max(need, (share + N - 1) / N);
      }
      return std::min(U, need);
    };
    const auto want = t2_reference::solve(*g, elig, roster, dests, N, V, u_of, hr);
    EXPECT_EQ(pl.live_heads, want.live);
    EXPECT_EQ(pl.common_share, want.U);
    if (!want.live.empty()) EXPECT_EQ(pl.fluid_share, want.fluid);
    ++checked;
    if (want.U == 0) {
      EXPECT_TRUE(r.program->tunnels().empty());
      continue;
    }
    ASSERT_TRUE(want.best.has_value());
    // The objectives of the program's own counts are the plan's and the brute-force minimum.
    Int128 J = 0;
    Int128 phi = 0;
    Amounts x(V.size(), 0);
    for (const auto& t : r.program->tunnels()) {
      for (EdgeId e : *t.edges) ++x[static_cast<std::size_t>(e)];
    }
    for (std::size_t e = 0; e < x.size(); ++e) {
      J = J + Int128(g->cost()[e]) * x[e];
      for (Amount j = 1; hr && g->capacity_units()[e] > 0 && j <= x[e]; ++j) {
        phi = phi + t2_reference::slot_price(g->capacity_units()[e], V[e], j, *pl.planned_per_tunnel);
      }
      EXPECT_LE(x[e] * *pl.planned_per_tunnel, V[e]);
    }
    EXPECT_EQ(J, wide(*pl.route_cost));
    EXPECT_EQ((std::pair {J, phi}), *want.best);
    if (hr) EXPECT_EQ(phi, wide(*pl.headroom_penalty));
    ++optimal;
    // k tunnels per live head, round order j * H + i with holes.
    std::vector<std::int32_t> want_index;
    for (std::int64_t j = 0; j < N; ++j) {
      for (std::int64_t i = 0; i < H; ++i) {
        if (std::count(want.live.begin(), want.live.end(), roster[static_cast<std::size_t>(i)]) != 0) {
          want_index.push_back(static_cast<std::int32_t>(j * H + i));
        }
      }
    }
    EXPECT_EQ(indices(*r.program), want_index);
    // F: the independent lossless bound of the program on V.
    std::vector<t1_reference::RefTunnel> ref;
    for (const auto& t : r.program->tunnels()) {
      t1_reference::RefTunnel rt;
      rt.head = t.head;
      rt.slot = static_cast<std::size_t>(std::find(roster.begin(), roster.end(), t.head) - roster.begin());
      rt.path = *t.edges;
      ref.push_back(rt);
    }
    EXPECT_EQ(pl.forwarding_bound, t1_reference::lossless_bound(roster.size(), ref, V));
    const auto K = static_cast<Amount>(want.live.size()) * N;
    if (pl.planned_per_tunnel == pl.common_share && want.live.size() == roster.size()) {
      EXPECT_GE(*pl.forwarding_bound - K * want.U, 0);
      EXPECT_LT(*pl.forwarding_bound - K * want.U, K);
    }
    // Up to N: the unpruned selection over the EXACT(k) programs at U(k) (unlimited design).
    std::vector<Int128> scores;
    for (std::int64_t k = 1; k <= N; ++k) {
      auto e = make(g, heads, dests, t2(k, LOSSLESS, hr), elig, std::make_shared<FlowLedger>(g, V));
      const auto rk = e->replace(INF);
      if (*rk.plan->common_share == 0) {
        scores.push_back(0);
        continue;
      }
      const std::int64_t M = H * k;
      auto fits = [&](Int128 z) {
        Amounts load(V.size(), 0);
        for (const auto& t : rk.program->tunnels()) {
          const Int128 amount = z / M + (Int128(t.index) < z % M ? 1 : 0);
          for (EdgeId x2 : *t.edges) load[static_cast<std::size_t>(x2)] += static_cast<Amount>(amount);
        }
        for (std::size_t x2 = 0; x2 < V.size(); ++x2) {
          if (load[x2] > V[x2]) return false;
        }
        return true;
      };
      Int128 z = 0;
      while (fits(z + 1)) z = z + 1;
      const Int128 capped = D.has_value() && Int128(*D) < z ? Int128(*D) : z;
      Int128 g_part = 0;
      for (std::int64_t i = 0; i < H; ++i) {
        if (std::count(want.live.begin(), want.live.end(), roster[static_cast<std::size_t>(i)]) == 0) continue;
        g_part = g_part + capped / H + (Int128(i) < capped % H ? 1 : 0);
      }
      scores.push_back(g_part);
    }
    Int128 target = -1;
    if (D.has_value()) {
      target = 0;
      for (std::int64_t i = 0; i < H; ++i) {
        if (std::count(want.live.begin(), want.live.end(), roster[static_cast<std::size_t>(i)]) != 0) {
          target = target + t1_reference::rr(*D, H, i);
        }
      }
    }
    std::int64_t chosen = 0;
    for (std::int64_t k = 1; k <= N && chosen == 0; ++k) {
      if (scores[static_cast<std::size_t>(k - 1)] == target) chosen = k;
    }
    if (chosen == 0) {
      chosen = 1;
      for (std::int64_t k = 2; k <= N; ++k) {
        if (scores[static_cast<std::size_t>(k - 1)] > scores[static_cast<std::size_t>(chosen - 1)]) chosen = k;
      }
    }
    auto up = make(g, heads, dests, t2(N, LOSSLESS, hr, true), elig, std::make_shared<FlowLedger>(g, V));
    EXPECT_EQ(up->replace(INF, design).plan->count, chosen);
    ++selections;
  }
  EXPECT_GT(checked, 200);
  EXPECT_GT(optimal, 100);
  EXPECT_GT(selections, 100);
}

} // namespace
