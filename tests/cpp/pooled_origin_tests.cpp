// One demand, several sources (docs/semantics.md, section 7): the configuration table and its
// error order, the pooled single-path rule, pooled CostOnly splitting (one split per shared
// router, attribution by highest averages, the equal origin's entry that cannot pass anything,
// Lossy forwarding through a source, weights frozen per epoch, the first-hop experiment),
// increments and zero calls, the controllers (Global and CostOnly MaxFlow) under both origins,
// min_cut as its own output, the result's accounting and per-source rows, the terminations and
// the evaluators. Each test checks the hand-derived values of the step's expected-result rows
// (LEDGER rows C, Q, P, I, M, G, S, T, V) on small graphs whose EdgeIds equal their input
// positions (asserted). The property tests check the controllers' certificates (conservation,
// maximality under the origin, no negative residual cycle) with an independent augmenting-path
// maximum flow and Bellman-Ford, and pooled records against their own routes; the pooled
// forwarding itself is compared with an independent reference in tests/py/test_pooled.py.
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

using Drops = std::vector<std::pair<EdgeId, Amount>>;
using Counters = std::vector<std::pair<EdgeId, Counter128>>;
using Row = std::array<Amount, 6>;  // source, offered (-1: none), admitted, delivered, source loss, dropped

constexpr Amount kNone = -1;

Counter128 c(Amount x) { return Counter128 {0, static_cast<std::uint64_t>(x)}; }
Amount lo(const Counter128& x) {
  EXPECT_EQ(x.hi, 0U);
  return static_cast<Amount>(x.lo);
}
Amount lo(const std::optional<Counter128>& x) { return x.has_value() ? lo(*x) : kNone; }

std::vector<Row> rows(const PlacementResult& r) {
  std::vector<Row> out;
  for (const auto& s : r.by_source) {
    out.push_back({s.source, lo(s.offered), lo(s.admitted), s.delivered, lo(s.source_loss), lo(s.dropped)});
  }
  return out;
}

// A graph whose rows are already in EdgeId order (source, cost, destination, input).
std::shared_ptr<const Graph> ordered(std::int32_t n, const std::vector<InEdge>& rows) {
  auto g = shared(make_graph(n, rows));
  for (std::size_t i = 0; i < rows.size(); ++i) EXPECT_EQ(g->input_index()[i], static_cast<std::int32_t>(i));
  return g;
}

// The ledger's graphs: (u, v, capacity, cost) in EdgeId order.
const std::vector<InEdge> kPool = {{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 1, 1}, {2, 3, 1, 1}};  // {0, 1} -> 3
const std::vector<InEdge> kB1 = {{0, 3, 2, 1}, {1, 2, 10, 1}, {2, 3, 10, 2}};                // {0, 1} -> 3
const std::vector<InEdge> kB2 = {{0, 1, 4, 1}, {0, 2, 2, 1}, {1, 2, 10, 1}};                 // {0, 1} -> 2
const std::vector<InEdge> kB2S = {{0, 2, 2, 1}, {1, 2, 10, 1}};                              // {0, 1} -> 2
const std::vector<InEdge> kB3Z = {{0, 2, 0, 1}, {1, 2, 10, 1}};                              // {0, 1} -> 2
const std::vector<InEdge> kB3P = {{0, 2, 10, 1}, {1, 2, 0, 1}};                              // {0, 1} -> 2
const std::vector<InEdge> kES = {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 2, 1}, {2, 3, 10, 2}};  // 0 -> 3
const std::vector<InEdge> kR4 = {{0, 2, 2, 9}, {1, 2, 10, 1}};                               // {0, 1} -> 2
const std::vector<InEdge> kInc = {{0, 2, 1, 1}, {1, 2, 3, 1}};                               // {0, 1} -> 2
const std::vector<InEdge> kOne = {{0, 1, 10, 1}};
const std::vector<InEdge> kRc = {{0, 1, 1, 1}, {0, 1, 1, 1}};
const std::vector<InEdge> kPar = {{0, 1, 5, 1}, {0, 1, 8, 1}};
const std::vector<InEdge> kChain = {{0, 2, 3, 1}, {1, 0, 10, 1}};  // {0, 1} -> 2, source 1 passes source 0
const std::vector<InEdge> kSnap = {{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 10, 1}, {2, 3, 10, 1}};  // {0, 1} -> 3
const std::vector<InEdge> kZero = {{0, 1, 0, 1}};
const std::vector<InEdge> kSpf = {{0, 1, 0, 1}, {0, 2, 5, 1}, {2, 1, 5, 1}};                  // 0 -> 1
const std::vector<InEdge> kRows = {{0, 3, 5, 1}, {1, 2, 5, 1}, {4, 3, 5, 1}};                 // {0, 1, 4} -> 3
const std::vector<InEdge> kPriv = {{0, 2, 4, 1}, {1, 2, 6, 1}};                               // {0, 1} -> 2
const std::vector<InEdge> kWit = {{0, 1, 1, 1}, {0, 2, 1, 1}, {1, 0, 1, 1}};                  // 0 -> 2

std::shared_ptr<const Graph> graph(const std::vector<InEdge>& edges) {
  NodeId n = 0;
  for (const auto& e : edges) n = std::max({n, e.u + 1, e.v + 1});
  return ordered(n, edges);
}

PathSelection sel(EdgeTieBreak t = EdgeTieBreak::Deterministic, EdgePriority p = EdgePriority::Capacity,
                  bool multipath = true) {
  return PathSelection(t, p, multipath, true);
}

FlowPolicyConfig cost_only(FlowPlacement k = FlowPlacement::EqualBalanced, Admission a = Admission::Lossless,
                           std::optional<PathSelection> s = std::nullopt, std::optional<Origin> o = std::nullopt) {
  return FlowPolicyConfig {CostOnly {}, k, a, s, o};
}
FlowPolicyConfig glob(std::optional<Origin> o = std::nullopt) {
  return FlowPolicyConfig {.routing = Global {}, .origin = o};
}

const PathSelection kPreferPath = sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Capacity, false);
const std::string kQ12 = "CostOnly with multiple input sources rejects multipath=false with PREFER_HIGHER_PRIORITY";
const OutputOptions kEdges {true, false, false, false};
const OutputOptions kCut {false, false, false, true};
const OutputOptions kAll {true, true, true, true};

// A policy on a fresh ledger (an optional baseline residual) of its own graph.
struct Case {
  std::shared_ptr<const Graph> g;
  FlowLedger ledger;
  std::optional<Eligibility> mask;
  std::unique_ptr<FlowPolicy> p;
  Case(const std::vector<InEdge>& edges, std::vector<NodeId> s, std::vector<NodeId> t, const FlowPolicyConfig& cfg,
       std::vector<NodeId> xn = {}, std::vector<Amount> residual = {})
      : g(graph(edges)), ledger(residual.empty() ? FlowLedger(g) : FlowLedger(g, residual)) {
    if (!xn.empty()) mask = g->eligibility(xn, {});
    p = std::make_unique<FlowPolicy>(ledger, s, t, cfg, mask ? &*mask : nullptr);
  }
  PlacementResult place(std::optional<Amount> v, const OutputOptions& o = {}) { return p->place(v, o); }
  std::vector<std::pair<std::vector<std::pair<EdgeId, Amount>>, Amount>> records() const {
    std::vector<std::pair<std::vector<std::pair<EdgeId, Amount>>, Amount>> out;
    for (const auto& r : p->records()) out.emplace_back(r.edges, r.delivered);
    return out;
  }
};

std::string construct_error(const std::vector<InEdge>& edges, std::vector<NodeId> s, std::vector<NodeId> t,
                            const FlowPolicyConfig& cfg, std::vector<NodeId> xn = {}) {
  auto g = graph(edges);
  FlowLedger fg(g);
  const auto mask = g->eligibility(xn, {});
  return error_of([&] { FlowPolicy p(fg, s, t, cfg, xn.empty() ? nullptr : &mask); });
}

}  // namespace

// ---------------------------------------------------------------- configuration (C rows)

// C4 to C12: each defect's text, at a policy and at capacity().
TEST(ConfigOrigin, MatrixTexts) {
  const std::vector<std::pair<FlowPolicyConfig, std::string>> rows = {
      {FlowPolicyConfig {.routing = CostOnly {}}, "CostOnly requires placement"},
      {FlowPolicyConfig {.routing = Global {}, .placement = FlowPlacement::MaxFlow}, "Global does not accept placement"},
      {FlowPolicyConfig {.routing = Global {}, .selection = PathSelection()}, "Global does not accept selection"},
      {cost_only(FlowPlacement::EqualBalanced, Admission::Lossless, std::nullopt, Origin::Optimized),
       "CostOnly with EQUAL_BALANCED or WEIGHTED requires origin=EQUAL"},
      {cost_only(FlowPlacement::Weighted, Admission::Lossless, std::nullopt, Origin::Optimized),
       "CostOnly with EQUAL_BALANCED or WEIGHTED requires origin=EQUAL"},
      {cost_only(FlowPlacement::MaxFlow, Admission::Lossy), "MAX_FLOW requires admission=LOSSLESS"},
      {FlowPolicyConfig {.routing = Global {}, .admission = Admission::Lossy}, "MAX_FLOW requires admission=LOSSLESS"},
      {FlowPolicyConfig {.routing = CapacityAware {}, .placement = FlowPlacement::MaxFlow, .origin = Origin::Optimized},
       "T1 requires origin=EQUAL"},
      {FlowPolicyConfig {.routing = CapacityAware {}}, "CapacityAware requires placement"},
      {FlowPolicyConfig {.routing = CapacityAware {}, .placement = FlowPlacement::EqualBalanced,
                         .admission = Admission::Lossy},
       "CapacityAware requires admission=LOSSLESS"},
  };
  const auto g = graph(kB1);
  for (const auto& [cfg, text] : rows) {
    EXPECT_EQ(construct_error(kB1, {0, 1}, {3}, cfg), text);
    EXPECT_EQ(error_of([&] { (void)capacity(*g, v({0, 1}), v({3}), QueryOptions {cfg}); }), text);
  }
  // C10, C12 and C8 for Pinned.
  Pinned pin;
  pin.bundles = {dag_of(*g, {0})};
  EXPECT_EQ(construct_error(kB1, {0}, {3}, FlowPolicyConfig {pin, FlowPlacement::MaxFlow, Admission::Lossless,
                                                             std::nullopt, Origin::Equal}),
            "Pinned does not accept origin");
  EXPECT_EQ(construct_error(kB1, {0}, {3}, FlowPolicyConfig {.routing = pin}), "Pinned requires placement");
  EXPECT_EQ(construct_error(kB1, {0}, {3}, FlowPolicyConfig {pin, FlowPlacement::MaxFlow, Admission::Lossy}),
            "MAX_FLOW requires admission=LOSSLESS");
  // C11: only CostOnly reads a cache (policy and evaluators).
  DagCache cache(*g);
  FlowLedger fg(g);
  const std::vector<std::tuple<FlowPolicyConfig, std::vector<NodeId>, std::string>> cached = {
      {glob(), {0, 1}, "Global does not use DagCache"},
      {FlowPolicyConfig {.routing = CapacityAware {}, .placement = FlowPlacement::MaxFlow}, {0, 1},
       "T1 does not use DagCache"},
      {FlowPolicyConfig {.routing = pin, .placement = FlowPlacement::MaxFlow}, {0}, "Pinned does not use DagCache"},
  };
  for (const auto& [cfg, s, text] : cached) {
    EXPECT_EQ(error_of([&] { FlowPolicy p(fg, s, v({3}), cfg, nullptr, &cache); }), text);
    EXPECT_EQ(error_of([&] { (void)capacity(*g, s, v({3}), QueryOptions {cfg, nullptr, &cache}); }), text);
    EXPECT_EQ(error_of([&] { (void)sensitivity(*g, s, v({3}), QueryOptions {cfg, nullptr, &cache}, 1); }), text);
  }
  EXPECT_EQ(cache.stats().misses + cache.stats().hits, 0);
}

// C13, C14, C16, Q6: the row, then the columns placement, origin, admission, selection; the table
// before the cache check.
TEST(ConfigOrigin, MatrixOrder) {
  EXPECT_EQ(construct_error(kB1, {0, 1}, {3},
                            FlowPolicyConfig {Global {}, FlowPlacement::MaxFlow, Admission::Lossy, PathSelection()}),
            "Global does not accept placement");
  EXPECT_EQ(construct_error(kB1, {0, 1}, {3},
                            FlowPolicyConfig {Global {}, std::nullopt, Admission::Lossy, PathSelection(), Origin::Equal}),
            "MAX_FLOW requires admission=LOSSLESS");
  EXPECT_EQ(construct_error(kPool, {0, 1}, {3},
                            cost_only(FlowPlacement::EqualBalanced, Admission::Lossless, kPreferPath, Origin::Optimized)),
            "CostOnly with EQUAL_BALANCED or WEIGHTED requires origin=EQUAL");
  const auto g = graph(kB1);
  DagCache cache(*g);
  FlowLedger fg(g);
  EXPECT_EQ(error_of([&] {
              FlowPolicy p(fg, v({0, 1}), v({3}), FlowPolicyConfig {.routing = Global {}, .placement = FlowPlacement::MaxFlow},
                           nullptr, &cache);
            }),
            "Global does not accept placement");
}

// C15: the endpoints come before the ownership of the mask.
TEST(ConfigOrigin, EndpointsBeforeOwnership) {
  const auto g = graph(kB1);
  const auto other = graph(kB1);
  const auto foreign = other->eligibility(v({0}), {});
  FlowLedger fg(g);
  const auto text = error_of([&] { FlowPolicy p(fg, std::vector<NodeId> {}, v({3}), glob(), &foreign); });
  EXPECT_FALSE(text.empty());
  EXPECT_FALSE(contains(text, "foreign")) << text;
}

// C19: finite after infinite, place(0) included, until remove().
TEST(ConfigOrigin, FiniteAfterInfinite) {
  for (const auto& cfg : {cost_only(), glob(), cost_only(FlowPlacement::EqualBalanced, Admission::Lossy)}) {
    Case k(kOne, {0}, {1}, cfg);
    (void)k.place(std::nullopt);
    for (Amount x : {Amount {1}, Amount {0}}) {
      EXPECT_EQ(error_of([&] { (void)k.place(x); }), "finite placement after infinite placement requires a new epoch");
    }
    k.p->remove();
    EXPECT_EQ(k.place(1).placed, 1);
  }
}

// ---------------------------------------------------------------- the pooled single-path rule (Q rows)

// Q1 to Q4: two normalized input sources (also when a mask leaves one) reject PreferHigherPriority
// single paths for every kernel; one source and multipath are accepted.
TEST(PooledSinglePath, Rejected) {
  for (auto k : {FlowPlacement::EqualBalanced, FlowPlacement::Weighted, FlowPlacement::MaxFlow}) {
    EXPECT_EQ(construct_error(kPool, {0, 1}, {3}, cost_only(k, Admission::Lossless, kPreferPath)), kQ12);
    EXPECT_EQ(construct_error(kPool, {0, 1}, {3}, cost_only(k, Admission::Lossless, kPreferPath), {1}), kQ12);
    EXPECT_EQ(construct_error(kPool, {0, 0}, {3}, cost_only(k, Admission::Lossless, kPreferPath)), "");
    EXPECT_EQ(construct_error(kPool, {0, 1}, {3},
                              cost_only(k, Admission::Lossless, sel(EdgeTieBreak::PreferHigherPriority))),
              "");
  }
}

// Q5: DETERMINISTIC single paths pool: 0: [e0, e2], 1: [e1, e2] (equal cost and hops: e2).
TEST(PooledSinglePath, DeterministicPaths) {
  const auto det = sel(EdgeTieBreak::Deterministic, EdgePriority::Capacity, false);
  Case k(kPool, {0, 1}, {3}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossless, det));
  const auto r = k.place(std::nullopt);
  EXPECT_EQ(r.placed, 1);
  EXPECT_EQ(rows(r), (std::vector<Row> {{0, kNone, 1, 1, 0, 0}, {1, kNone, 0, 0, 0, 0}}));
  Case lossy(kPool, {0, 1}, {3}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy, det));
  const auto l = lossy.place(2);
  EXPECT_EQ(l.placed, 1);
  EXPECT_EQ(l.call_drops, (Drops {{2, 1}}));
  EXPECT_EQ(rows(l), (std::vector<Row> {{0, 1, 1, 1, 0, 0}, {1, 1, 1, 0, 0, 1}}));
  const auto recs = lossy.p->records();
  ASSERT_EQ(recs.size(), 2U);
  EXPECT_EQ(ascending_edges(*recs[0].bundle), (std::vector<EdgeId> {0, 2}));
  EXPECT_EQ(ascending_edges(*recs[1].bundle), (std::vector<EdgeId> {1, 2}));
}

// ---------------------------------------------------------------- pooled CostOnly (P rows)

// P1 to P3 (B-pooled): one split at the shared router, attributed member by member.
TEST(Pooled, SharedRouterSplitsOnce) {
  Case k(kPool, {0, 1}, {3}, cost_only());
  const auto r = k.place(2, kEdges);
  EXPECT_EQ(r.placed, 2);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_EQ(*r.edge_flows, (std::vector<Amount> {1, 1, 1, 1}));
  EXPECT_EQ(rows(r), (std::vector<Row> {{0, kNone, 1, 1, 0, 0}, {1, kNone, 1, 1, 0, 0}}));
  EXPECT_EQ(k.records(), (decltype(k.records()) {{{{0, 1}, {2, 1}}, 1}, {{{1, 1}, {3, 1}}, 1}}));
  Case inf(kPool, {0, 1}, {3}, cost_only());
  const auto i = inf.place(std::nullopt);
  EXPECT_EQ(i.placed, 2);
  EXPECT_EQ(i.termination, Termination::Exhausted);
  EXPECT_FALSE(i.requested.has_value());
  EXPECT_FALSE(i.unadmitted.has_value());
  EXPECT_FALSE(i.remaining.has_value());
  Case w(kPool, {0, 1}, {3}, cost_only(FlowPlacement::Weighted));
  EXPECT_EQ(w.place(2).placed, 2);
}

// P4: each member's arrivals divided by what each source brought.
TEST(Pooled, LossyAttributionByArrivals) {
  Case k(kPool, {0, 1}, {3}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy));
  const auto r = k.place(4);
  EXPECT_EQ(r.placed, 2);
  EXPECT_EQ(r.call_drops, (Drops {{2, 1}, {3, 1}}));
  EXPECT_EQ(rows(r), (std::vector<Row> {{0, 2, 2, 2, 0, 0}, {1, 2, 2, 0, 0, 2}}));
  for (const auto& rec : k.p->records()) expect_conserved(*k.g, rec, "P4");
}

// P5 to P7, P9 (B3-zero-first, B3-positive-first): Lossless admits only totals that give the
// source behind the zero edge nothing; Lossy drops its share on that edge.
TEST(Pooled, EqualEntryThatCannotPass) {
  for (const auto& [edges, want] : {std::pair {kB3Z, Amount {0}}, std::pair {kB3P, Amount {1}}}) {
    Case eb(edges, {0, 1}, {2}, cost_only());
    const auto r = eb.place(10);
    EXPECT_EQ(r.placed, want);
    EXPECT_EQ(r.termination, Termination::Exhausted);
    EXPECT_EQ(eb.p->records().size(), 2U);
    Case w(edges, {0, 1}, {2}, cost_only(FlowPlacement::Weighted));
    EXPECT_EQ(w.place(10).placed, want);
  }
  Case lossy(kB3Z, {0, 1}, {2}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy));
  const auto l = lossy.place(10);
  EXPECT_EQ(l.placed, 5);
  EXPECT_EQ(l.call_drops, (Drops {{0, 5}}));
  EXPECT_EQ(l.call_source_loss, 0);
  EXPECT_EQ(l.termination, Termination::Complete);
}

// P8: Weighted gives the zero edge weight 0; source 0's share has no member and is source loss.
TEST(Pooled, WeightedZeroRouteIsSourceLoss) {
  Case k(kB3Z, {0, 1}, {2}, cost_only(FlowPlacement::Weighted, Admission::Lossy));
  const auto r = k.place(10);
  EXPECT_EQ(r.placed, 5);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_TRUE(r.call_drops.empty());
  EXPECT_EQ(r.call_source_loss, 5);
  EXPECT_EQ(lo(r.source_loss), 5);
  EXPECT_EQ(rows(r), (std::vector<Row> {{0, 5, 5, 0, 5, 0}, {1, 5, 5, 5, 0, 0}}));
}

// P10, P11 on B1.
TEST(Pooled, EqualOriginLimitsTheTotal) {
  Case k(kB1, {0, 1}, {3}, cost_only());
  const auto r = k.place(12);
  EXPECT_EQ(r.placed, 4);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  EXPECT_EQ(rows(r), (std::vector<Row> {{0, kNone, 2, 2, 0, 0}, {1, kNone, 2, 2, 0, 0}}));
  Case lossy(kB1, {0, 1}, {3}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy));
  const auto l = lossy.place(12);
  EXPECT_EQ(l.placed, 8);
  EXPECT_EQ(l.call_drops, (Drops {{0, 4}}));
}

// P12, P13, P17 on CHAIN: source 1's traffic passes source 0's router.
TEST(Pooled, TransitThroughASource) {
  Case k(kChain, {0, 1}, {2}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy));
  const auto r = k.place(10);
  EXPECT_EQ(r.placed, 3);
  EXPECT_EQ(r.call_drops, (Drops {{0, 7}}));
  EXPECT_EQ(rows(r), (std::vector<Row> {{0, 5, 5, 2, 0, 3}, {1, 5, 5, 1, 0, 4}}));
  EXPECT_EQ(k.records(), (decltype(k.records()) {{{{0, 2}}, 2}, {{{0, 1}, {1, 5}}, 1}}));
  for (const auto& rec : k.p->records()) expect_conserved(*k.g, rec, "P12");
  Case lossless(kChain, {0, 1}, {2}, cost_only());
  const auto l = lossless.place(std::nullopt);
  EXPECT_EQ(l.placed, 3);
  EXPECT_EQ(rows(l), (std::vector<Row> {{0, kNone, 2, 2, 0, 0}, {1, kNone, 1, 1, 0, 0}}));
  Case hop(kChain, {0, 1}, {2}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy));
  const auto h = hop.place(std::nullopt);
  EXPECT_EQ(h.placed, 3);
  EXPECT_EQ(h.termination, Termination::Exhausted);
  EXPECT_EQ(lo(h.offered), 13);
  EXPECT_EQ(h.drops, (Counters {{0, c(10)}}));
  EXPECT_EQ(rows(h), (std::vector<Row> {{0, 3, 3, 0, 0, 3}, {1, 10, 10, 3, 0, 7}}));
}

// P14: Residual weights read once, on the residual at the first positive call's entry.
TEST(Pooled, ResidualWeightsOnTheCommonView) {
  Case k(kSnap, {0, 1}, {3}, cost_only(FlowPlacement::Weighted, Admission::Lossless, sel(EdgeTieBreak::Deterministic, EdgePriority::Residual)),
         {}, {10, 10, 2, 10});
  const auto r = k.place(std::nullopt);
  EXPECT_EQ(r.placed, 12);
  EXPECT_EQ(k.records(), (decltype(k.records()) {{{{0, 6}, {2, 1}, {3, 5}}, 6}, {{{1, 6}, {2, 1}, {3, 5}}, 6}}));
}

// P15: the weights stay frozen for the epoch after the other owner leaves.
TEST(Pooled, WeightsFrozenForTheEpoch) {
  const auto g = graph(kSnap);
  FlowLedger fg(g);
  Pinned pin;
  pin.bundles = {dag_of(*g, {2})};
  auto other = std::make_unique<FlowPolicy>(fg, v({2}), v({3}), FlowPolicyConfig {.routing = pin, .placement = FlowPlacement::EqualBalanced});
  EXPECT_EQ(other->place(8).placed, 8);
  FlowPolicy p(fg, v({0, 1}), v({3}),
               cost_only(FlowPlacement::Weighted, Admission::Lossless, sel(EdgeTieBreak::Deterministic, EdgePriority::Residual)));
  EXPECT_EQ(p.place(2).placed, 2);
  other->remove();
  const auto r = p.place(10, kEdges);
  EXPECT_EQ(r.placed, 10);
  EXPECT_EQ(*r.edge_flows, (std::vector<Amount> {6, 6, 2, 10}));
}

// P16: the first-hop experiment.
TEST(Pooled, FirstHopExperiment) {
  Case k(kPriv, {0, 1}, {2}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy));
  const auto r = k.place(std::nullopt);
  EXPECT_EQ(r.placed, 10);
  EXPECT_EQ(lo(r.offered), 10);
  EXPECT_FALSE(r.requested.has_value());
  EXPECT_EQ(r.termination, Termination::Exhausted);
}

// P18, P19 (properties): on random graphs with several sources, masks, zero capacities and dispersed
// costs, every pooled record (EB and W, both admissions, multipath or deterministic single paths)
// carries and drops only on its own route, conserves its traffic, and the rows sum to the totals
// with the accounting identities; the ledger stays within the capacities.
TEST(Pooled, RecordsStayOnTheirRoutesAndAccountExactly) {
  std::mt19937 rng(20261018);
  auto pick = [&](int n) { return static_cast<int>(rng() % static_cast<unsigned>(n)); };
  int multi = 0, carrying = 0;
  for (int trial = 0; trial < 300; ++trial) {
    const int n = 4 + pick(5);
    std::vector<InEdge> es;
    for (int k = 0, m = n + pick(2 * n); k < m; ++k) {
      const NodeId u = pick(n);
      const NodeId w = pick(3) == 0 ? n - 1 : pick(n);
      if (u != w) es.push_back({u, w, pick(4) == 0 ? 0 : 1 + pick(6), static_cast<Cost>(1 + pick(3))});
    }
    if (es.empty()) es.push_back({0, n - 1, 1, 1});
    auto g = shared(make_graph(n, es));
    std::vector<NodeId> s = {0, 1};
    if (pick(2) == 0) s.push_back(2);
    const std::vector<NodeId> t = {static_cast<NodeId>(n - 1)};
    const auto mask = g->eligibility(pick(3) == 0 ? v({static_cast<NodeId>(3)}) : std::vector<NodeId> {},
                                     pick(3) == 0 ? std::vector<EdgeId> {static_cast<EdgeId>(pick(g->num_edges()))}
                                                  : std::vector<EdgeId> {});
    const auto k = pick(2) == 0 ? FlowPlacement::EqualBalanced : FlowPlacement::Weighted;
    const auto a = pick(2) == 0 ? Admission::Lossless : Admission::Lossy;
    const auto prio = pick(2) == 0 ? EdgePriority::Capacity : EdgePriority::Residual;
    const auto selection = sel(EdgeTieBreak::Deterministic, prio, pick(3) != 0);
    FlowLedger fg(g);
    FlowPolicy p(fg, s, t, cost_only(k, a, selection), &mask);
    const std::string where = "trial " + std::to_string(trial);
    for (int call = 0, calls = 1 + pick(3); call < calls; ++call) {
      const std::optional<Amount> vol = pick(5) == 0 && call == calls - 1 ? std::nullopt : std::optional<Amount> {pick(9)};
      const auto r = p.place(vol);
      Amount admitted = 0, delivered = 0, loss = 0, dropped = 0;
      for (const auto& row : r.by_source) {
        admitted += lo(row.admitted);
        delivered += row.delivered;
        loss += lo(row.source_loss);
        dropped += lo(row.dropped);
      }
      Amount drops = 0;
      for (const auto& [e, x] : r.drops) drops += lo(x);
      const bool empty = fg.records().empty();
      if (!empty || a == Admission::Lossless) {
        EXPECT_EQ(admitted, lo(r.admitted)) << where;
        EXPECT_EQ(loss, lo(r.source_loss)) << where;
      }
      EXPECT_EQ(delivered, r.delivered) << where;
      EXPECT_EQ(dropped, drops) << where;
      if (a == Admission::Lossy) {
        EXPECT_EQ(lo(r.offered), lo(r.source_loss) + drops + r.delivered) << where;
        EXPECT_EQ(lo(r.unadmitted), 0) << where;
      } else if (r.requested.has_value()) {
        EXPECT_EQ(lo(r.requested), lo(r.unadmitted) + r.delivered) << where;
      }
      EXPECT_EQ(r.delivery_delta, r.placed) << where;
      for (const auto& rec : p.records()) {
        const auto route = ascending_edges(*rec.bundle);
        for (const auto& [e, x] : rec.edges) {
          EXPECT_TRUE(std::binary_search(route.begin(), route.end(), e)) << where << " edge " << e;
        }
        for (const auto& [e, x] : rec.drops) {
          EXPECT_TRUE(std::binary_search(route.begin(), route.end(), e)) << where << " drop " << e;
        }
        if (expect_conserved(*g, rec, where)) ++carrying;
      }
      const auto alloc = fg.allocated();
      for (EdgeId e = 0; e < g->num_edges(); ++e) {
        EXPECT_LE(alloc[static_cast<std::size_t>(e)], g->capacity_units()[static_cast<std::size_t>(e)]) << where;
      }
      if (!vol.has_value()) break;
    }
    if (p.records().size() > 1) ++multi;
  }
  EXPECT_GT(multi, 100);
  EXPECT_GT(carrying, 200);
}

// ---------------------------------------------------------------- increments and zero calls (I rows)

// I1 to I4.
TEST(Increments, LosslessRefusedUnitsStayRefused) {
  Case k(kInc, {0, 1}, {2}, cost_only());
  const auto a = k.place(3);
  const auto b = k.place(1);
  EXPECT_EQ(a.placed, 2);
  EXPECT_EQ(b.placed, 0);
  EXPECT_EQ(b.termination, Termination::Exhausted);
  EXPECT_EQ(lo(b.requested), 4);
  EXPECT_EQ(lo(b.admitted), 2);
  EXPECT_EQ(lo(b.unadmitted), 2);
  EXPECT_FALSE(b.offered.has_value());
  Case lossy(kInc, {0, 1}, {2}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy));
  const auto la = lossy.place(3);
  const auto lb = lossy.place(1);
  EXPECT_EQ(la.placed, 2);
  EXPECT_EQ(la.call_drops, (Drops {{0, 1}}));
  EXPECT_EQ(lb.placed, 1);
  EXPECT_TRUE(lb.call_drops.empty());
  EXPECT_EQ(lb.drops, (Counters {{0, c(1)}}));
  EXPECT_EQ(lo(lb.offered), 4);
  EXPECT_EQ(lb.delivered, 3);
  EXPECT_EQ(rows(lb), (std::vector<Row> {{0, 2, 2, 1, 0, 1}, {1, 2, 2, 2, 0, 0}}));
  // I3: another owner holds 8, leaves, places 5.
  const auto g = graph(kOne);
  FlowLedger fg(g);
  FlowPolicy other(fg, v({0}), v({1}), cost_only());
  (void)other.place(8);
  FlowPolicy p(fg, v({0}), v({1}), cost_only());
  EXPECT_EQ(p.place(5).placed, 2);
  other.remove();
  (void)other.place(5);
  const auto r = p.place(1);
  EXPECT_EQ(r.placed, 1);
  EXPECT_EQ(lo(r.requested), 6);
  EXPECT_EQ(lo(r.unadmitted), 3);
  // I4.
  Case rc(kRc, {0}, {1}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy));
  EXPECT_EQ(rc.place(1).placed, 1);
  EXPECT_EQ(rc.place(1).placed, 0);
  Case rc2(kRc, {0}, {1}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy));
  EXPECT_EQ(rc2.place(2).placed, 2);
  Case rc3(kRc, {0}, {1}, cost_only());
  EXPECT_EQ(rc3.place(1).placed, 1);
  EXPECT_EQ(*rc3.place(1, kEdges).edge_flows, (std::vector<Amount> {1, 1}));
}

// I5, I6 (RC02-zero-unprepared): a zero call before preparation, after a positive call, after
// remove().
TEST(Epochs, ZeroCallBeforePreparation) {
  for (const bool global : {false, true}) {
    Case k(kPar, {0}, {1}, global ? glob() : cost_only(), {}, {4, 6});
    auto unprepared = [&](const PlacementResult& r) {
      EXPECT_EQ(r.placed, 0);
      EXPECT_EQ(r.remaining, 0);
      EXPECT_EQ(r.delivery_delta, 0);
      EXPECT_TRUE(r.call_drops.empty());
      EXPECT_EQ(r.call_source_loss, 0);
      EXPECT_EQ(r.termination, Termination::Complete);
      EXPECT_FALSE(r.routes_changed || r.weights_changed);
      EXPECT_EQ(*r.edge_flows, (std::vector<Amount> {0, 0}));
      EXPECT_EQ(*r.residual, (std::vector<Amount> {4, 6}));
      EXPECT_FALSE(r.min_cut.has_value());
      EXPECT_FALSE(r.cost_distribution.has_value());
      EXPECT_EQ(lo(r.requested), 0);
      EXPECT_EQ(rows(r), (std::vector<Row> {{0, kNone, 0, 0, 0, 0}}));
      EXPECT_EQ(r.total_cost.has_value(), global);
      if (global) EXPECT_EQ(*r.total_cost, c(0));
      EXPECT_TRUE(k.p->records().empty());
    };
    unprepared(k.place(0, kAll));
    (void)k.place(3);
    const auto z = k.place(0, kAll);
    EXPECT_EQ(z.delivered, 3);
    EXPECT_FALSE(z.routes_changed);
    if (global) EXPECT_EQ(*z.cost_distribution, (std::vector<std::pair<Cost, Amount>> {{1, 3}}));
    k.p->remove();
    unprepared(k.place(0, kAll));
  }
}

// ---------------------------------------------------------------- the controllers (M, G rows)

// M1 to M7: one controller on the union of the sources' shortest routes.
TEST(CostOnlyMaxFlow, OptimizedPoolsTheRoster) {
  struct Want {
    std::vector<InEdge> edges;
    std::vector<NodeId> s;
    NodeId t;
    std::optional<Origin> origin;
    Amount request, placed, cost;
    std::vector<std::vector<std::pair<EdgeId, Amount>>> records;
  };
  const std::vector<Want> wants = {
      {kB1, {0, 1}, 3, std::nullopt, 12, 12, 32, {{{0, 2}}, {{1, 10}, {2, 10}}}},
      {kB1, {0, 1}, 3, Origin::Equal, 12, 4, 8, {{{0, 2}}, {{1, 2}, {2, 2}}}},
      {kB2, {0, 1}, 2, Origin::Equal, 12, 4, 4, {{{1, 2}}, {{2, 2}}}},
      {kB2, {0, 1}, 2, Origin::Optimized, 12, 12, 12, {{{1, 2}}, {{2, 10}}}},
      {kES, {0}, 3, std::nullopt, 12, 2, 4, {{{0, 2}, {2, 2}}}},
      {kSpf, {0}, 1, std::nullopt, 5, 0, 0, {{}}},
  };
  for (const auto& w : wants) {
    Case k(w.edges, w.s, {w.t}, cost_only(FlowPlacement::MaxFlow, Admission::Lossless, std::nullopt, w.origin));
    const auto r = k.place(w.request, kCut);
    EXPECT_EQ(r.placed, w.placed);
    EXPECT_EQ(r.termination, w.placed == w.request ? Termination::Complete : Termination::Exhausted);
    EXPECT_EQ(*r.total_cost, c(w.cost));
    EXPECT_FALSE(r.cost_distribution.has_value());
    EXPECT_FALSE(r.min_cut.has_value());
    std::vector<std::vector<std::pair<EdgeId, Amount>>> got;
    for (const auto& rec : k.p->records()) got.push_back(rec.edges);
    EXPECT_EQ(got, w.records);
  }
  // M6.
  Case r4(kR4, {0, 1}, {2}, cost_only(FlowPlacement::MaxFlow, Admission::Lossless, std::nullopt, Origin::Equal));
  EXPECT_EQ(r4.place(3).placed, 3);
  EXPECT_EQ(r4.place(3).placed, 1);
  EXPECT_EQ(*r4.place(0).total_cost, c(20));
}

// G1 to G7, G11.
TEST(GlobalEqual, Examples) {
  struct Want {
    std::vector<InEdge> edges;
    std::vector<NodeId> s;
    NodeId t;
    Origin origin;
    Amount request, placed, cost;
    std::vector<Amount> flows;
  };
  const std::vector<Want> wants = {
      {kB1, {0, 1}, 3, Origin::Equal, 12, 4, 8, {2, 2, 2}},
      {kB1, {0, 1}, 3, Origin::Optimized, 12, 12, 32, {2, 10, 10}},
      {kB2, {0, 1}, 2, Origin::Equal, 12, 12, 16, {4, 2, 10}},
      {kB2, {0, 1}, 2, Origin::Optimized, 12, 12, 12, {0, 2, 10}},
      {kB2S, {0, 1}, 2, Origin::Equal, 12, 4, 4, {2, 2}},
      {kB2S, {0, 1}, 2, Origin::Optimized, 12, 12, 12, {2, 10}},
      {kB3Z, {0, 1}, 2, Origin::Equal, 10, 0, 0, {0, 0}},
      {kB3P, {0, 1}, 2, Origin::Equal, 10, 1, 1, {1, 0}},
      {kB3Z, {0, 1}, 2, Origin::Optimized, 10, 10, 10, {0, 10}},
      {kES, {0}, 3, Origin::Optimized, 12, 12, 34, {2, 10, 2, 10}},
      {kWit, {0}, 2, Origin::Optimized, 1, 1, 1, {0, 1, 0}},
  };
  for (const auto& w : wants) {
    Case k(w.edges, w.s, {w.t}, glob(w.origin));
    const auto r = k.place(w.request, kEdges);
    EXPECT_EQ(r.placed, w.placed);
    EXPECT_EQ(r.remaining, w.request - w.placed);
    EXPECT_EQ(r.termination, w.placed == w.request ? Termination::Complete : Termination::Exhausted);
    EXPECT_EQ(*r.total_cost, c(w.cost));
    EXPECT_EQ(*r.edge_flows, w.flows);
    EXPECT_EQ(r.cost_distribution.has_value(), w.origin == Origin::Optimized);
  }
  // G3: source 0's 6 is 2 direct and 4 through node 1.
  Case b2(kB2, {0, 1}, {2}, glob(Origin::Equal));
  (void)b2.place(12);
  EXPECT_EQ(b2.records(), (decltype(b2.records()) {{{{0, 4}, {1, 2}, {2, 4}}, 6}, {{{2, 6}}, 6}}));
  // G6 B3-delete-zero.
  Case masked(kB3Z, {0, 1}, {2}, glob(Origin::Equal), {0});
  EXPECT_EQ(masked.place(10).placed, 10);
}

// G8 to G10.
TEST(GlobalEqual, Increments) {
  Case k(kR4, {0, 1}, {2}, glob(Origin::Equal));
  const auto a = k.place(3);
  const auto b = k.place(3, kEdges);
  EXPECT_EQ(a.placed, 3);
  EXPECT_EQ(a.termination, Termination::Complete);
  EXPECT_EQ(b.placed, 1);
  EXPECT_EQ(b.termination, Termination::Exhausted);
  EXPECT_EQ(*b.edge_flows, (std::vector<Amount> {2, 2}));
  EXPECT_EQ(*b.total_cost, c(20));
  EXPECT_EQ(lo(b.requested), 6);
  EXPECT_EQ(lo(b.unadmitted), 2);
  const auto g = graph(kOne);
  FlowLedger fg(g);
  FlowPolicy other(fg, v({0}), v({1}), cost_only());
  (void)other.place(8);
  FlowPolicy p(fg, v({0}), v({1}), glob());
  EXPECT_EQ(p.place(5).placed, 2);
  other.remove();
  (void)other.place(5);
  EXPECT_EQ(p.place(1).placed, 1);
  Case one(kB2, {0, 1}, {2}, glob(Origin::Equal));
  const auto whole = one.place(12, kEdges);
  Case two(kB2, {0, 1}, {2}, glob(Origin::Equal));
  (void)two.place(6);
  const auto split = two.place(6, kEdges);
  EXPECT_EQ(*split.edge_flows, *whole.edge_flows);
  EXPECT_EQ(*split.total_cost, *whole.total_cost);
  EXPECT_EQ(two.records(), one.records());
}

namespace {

// Independent maximum flow (augmenting paths found by BFS) from a virtual root with an arc of
// capacity caps[i] (nullopt: unbounded) to roots[i], to a virtual sink fed by every destination,
// over the edges with `allowed`.
Amount max_flow(const Graph& g, const std::vector<Amount>& view, const std::vector<char>& allowed,
                const std::vector<NodeId>& roots, const std::vector<std::optional<Amount>>& caps,
                const std::vector<NodeId>& dsts) {
  struct Arc {
    std::size_t u, v;
    Amount cap, flow;
  };
  const auto n = static_cast<std::size_t>(g.num_nodes());
  const std::size_t root = n, sink = n + 1;
  Amount big = 1;
  for (Amount x : view) big += x;
  std::vector<Arc> arcs;
  for (std::size_t e = 0; e < view.size(); ++e) {
    if (allowed[e]) arcs.push_back({static_cast<std::size_t>(g.edge_src()[e]), static_cast<std::size_t>(g.edge_dst()[e]), view[e], 0});
  }
  for (std::size_t i = 0; i < roots.size(); ++i) arcs.push_back({root, static_cast<std::size_t>(roots[i]), caps[i].value_or(big), 0});
  for (NodeId d : dsts) arcs.push_back({static_cast<std::size_t>(d), sink, big, 0});
  Amount total = 0;
  for (;;) {
    std::vector<std::pair<std::size_t, int>> parent(n + 2, {SIZE_MAX, 0});
    std::vector<char> seen(n + 2, 0);
    std::vector<std::size_t> queue = {root};
    seen[root] = 1;
    for (std::size_t qi = 0; qi < queue.size(); ++qi) {
      const std::size_t u = queue[qi];
      for (std::size_t i = 0; i < arcs.size(); ++i) {
        if (arcs[i].u == u && arcs[i].flow < arcs[i].cap && !seen[arcs[i].v]) {
          seen[arcs[i].v] = 1;
          parent[arcs[i].v] = {i, 1};
          queue.push_back(arcs[i].v);
        }
        if (arcs[i].v == u && arcs[i].flow > 0 && !seen[arcs[i].u]) {
          seen[arcs[i].u] = 1;
          parent[arcs[i].u] = {i, -1};
          queue.push_back(arcs[i].u);
        }
      }
    }
    if (!seen[sink]) return total;
    Amount amount = big;
    for (std::size_t x = sink; x != root;) {
      const auto [i, d] = parent[x];
      amount = std::min(amount, d == 1 ? arcs[i].cap - arcs[i].flow : arcs[i].flow);
      x = d == 1 ? arcs[i].u : arcs[i].v;
    }
    for (std::size_t x = sink; x != root;) {
      const auto [i, d] = parent[x];
      arcs[i].flow += d * amount;
      x = d == 1 ? arcs[i].u : arcs[i].v;
    }
    total += amount;
  }
}

// Whether the residual graph of `flows` on the allowed edges, with a virtual sink fed by every
// destination (and arcs back to a destination that receives) and, when `free_supply`, a virtual
// root with arcs to every root (and back from a root that sends), has a negative-cost cycle
// (Bellman-Ford from a virtual start joined to every node).
bool negative_cycle(const Graph& g, const std::vector<Amount>& view, const std::vector<char>& allowed,
                    const std::vector<Amount>& flows, const std::vector<NodeId>& roots,
                    const std::vector<Amount>& supply, const std::vector<NodeId>& dsts,
                    const std::vector<Amount>& inflow, bool free_supply) {
  struct Arc {
    std::size_t u, v;
    Int128 w;
  };
  const auto n = static_cast<std::size_t>(g.num_nodes());
  const std::size_t root = n, sink = n + 1;
  std::vector<Arc> arcs;
  for (std::size_t e = 0; e < view.size(); ++e) {
    if (!allowed[e]) continue;
    const auto u = static_cast<std::size_t>(g.edge_src()[e]);
    const auto w = static_cast<std::size_t>(g.edge_dst()[e]);
    if (flows[e] < view[e]) arcs.push_back({u, w, g.cost()[e]});
    if (flows[e] > 0) arcs.push_back({w, u, -static_cast<Int128>(g.cost()[e])});
  }
  for (std::size_t i = 0; i < dsts.size(); ++i) {
    arcs.push_back({static_cast<std::size_t>(dsts[i]), sink, 0});
    if (inflow[i] > 0) arcs.push_back({sink, static_cast<std::size_t>(dsts[i]), 0});
  }
  if (free_supply) {
    for (std::size_t i = 0; i < roots.size(); ++i) {
      arcs.push_back({root, static_cast<std::size_t>(roots[i]), 0});
      if (supply[i] > 0) arcs.push_back({static_cast<std::size_t>(roots[i]), root, 0});
    }
  }
  std::vector<Int128> dist(n + 2, 0);
  for (std::size_t round = 0; round <= n + 2; ++round) {
    bool changed = false;
    for (const auto& a : arcs) {
      if (dist[a.u] + a.w < dist[a.v]) {
        dist[a.v] = dist[a.u] + a.w;
        changed = true;
      }
    }
    if (!changed) return false;
  }
  return true;
}

// The roster (sources that reach a destination over allowed nodes and edges, capacity-blind) and
// the static distance to the destinations, by Bellman-Ford.
std::vector<std::optional<Int128>> distances(const Graph& g, const std::vector<NodeId>& dsts, const Eligibility* elig) {
  const auto n = static_cast<std::size_t>(g.num_nodes());
  std::vector<std::optional<Int128>> d(n);
  auto node_ok = [&](NodeId x) { return elig == nullptr || elig->node_ok(x); };
  const std::set<NodeId> D(dsts.begin(), dsts.end());
  for (NodeId t : dsts) {
    if (node_ok(t)) d[static_cast<std::size_t>(t)] = 0;
  }
  for (std::size_t round = 0; round < n; ++round) {
    for (EdgeId e = 0; e < g.num_edges(); ++e) {
      const NodeId u = g.edge_src()[static_cast<std::size_t>(e)];
      const NodeId w = g.edge_dst()[static_cast<std::size_t>(e)];
      if ((elig != nullptr && !elig->edge_ok(e)) || !node_ok(u) || !node_ok(w) || D.count(u)) continue;
      const auto& dw = d[static_cast<std::size_t>(w)];
      auto& du = d[static_cast<std::size_t>(u)];
      if (dw && (!du || *dw + g.cost()[static_cast<std::size_t>(e)] < *du)) du = *dw + g.cost()[static_cast<std::size_t>(e)];
    }
  }
  return d;
}

// The controllers' certificates on random graphs with masks, zero capacities and dispersed costs:
// Global (every allowed edge) or CostOnly MaxFlow (the tight edges the roster reaches) under one
// origin. Flow only on the support, conservation and bounds; supplies RR(T) under the equal
// origin; maximality (equal: RR(T + 1) does not fit; optimized: the maximum flow); minimum cost (no
// negative residual cycle, with free supplies only under the optimized origin).
void certify_controllers(bool global) {
  std::mt19937 rng(global ? 20261019 : 20261020);
  auto pick = [&](int n) { return static_cast<int>(rng() % static_cast<unsigned>(n)); };
  int positive = 0, differ = 0;
  for (int trial = 0; trial < 200; ++trial) {
    const int n = 4 + pick(5);
    std::vector<InEdge> es;
    for (int k = 0, m = n + pick(2 * n); k < m; ++k) {
      const NodeId u = pick(n);
      const NodeId w = pick(3) == 0 ? n - 1 : pick(n);
      if (u != w) es.push_back({u, w, pick(5), static_cast<Cost>(1 + pick(4))});
    }
    if (es.empty()) es.push_back({0, n - 1, 1, 1});
    auto g = shared(make_graph(n, es));
    std::vector<NodeId> s = {0, 1};
    if (pick(2) == 0) s.push_back(2);
    std::vector<NodeId> t = {static_cast<NodeId>(n - 1)};
    if (pick(3) == 0 && n - 2 > 2) t.push_back(static_cast<NodeId>(n - 2));
    const auto mask = g->eligibility(pick(3) == 0 && n > 5 ? v({static_cast<NodeId>(3)}) : std::vector<NodeId> {},
                                     pick(3) == 0 ? std::vector<EdgeId> {static_cast<EdgeId>(pick(g->num_edges()))}
                                                  : std::vector<EdgeId> {});
    const auto m = static_cast<std::size_t>(g->num_edges());
    const auto dist = distances(*g, t, &mask);
    std::vector<NodeId> roster;
    for (NodeId x : s) {
      if (mask.node_ok(x) && dist[static_cast<std::size_t>(x)]) roster.push_back(x);
    }
    std::vector<char> allowed(m, 0);
    for (std::size_t e = 0; e < m; ++e) {
      const NodeId u = g->edge_src()[e], w = g->edge_dst()[e];
      allowed[e] = mask.edge_ok(static_cast<EdgeId>(e)) && mask.node_ok(u) && mask.node_ok(w);
    }
    if (!global) {
      // The tight edges the roster reaches, not out of a destination.
      std::vector<char> reach(static_cast<std::size_t>(n), 0), support(m, 0);
      std::vector<NodeId> stack = roster;
      for (NodeId x : roster) reach[static_cast<std::size_t>(x)] = 1;
      while (!stack.empty()) {
        const NodeId u = stack.back();
        stack.pop_back();
        if (std::find(t.begin(), t.end(), u) != t.end()) continue;
        for (std::size_t e = 0; e < m; ++e) {
          const NodeId w = g->edge_src()[e] == u ? g->edge_dst()[e] : -1;
          if (w < 0 || !allowed[e]) continue;
          const auto& du = dist[static_cast<std::size_t>(u)];
          const auto& dw = dist[static_cast<std::size_t>(w)];
          if (du && dw && *du == *dw + g->cost()[e]) {
            support[e] = 1;
            if (!reach[static_cast<std::size_t>(w)]) {
              reach[static_cast<std::size_t>(w)] = 1;
              stack.push_back(w);
            }
          }
        }
      }
      allowed = support;
    }
    std::vector<Amount> view(g->capacity_units().begin(), g->capacity_units().end());
    std::vector<Amount> placed_by_origin;
    for (const auto origin : {Origin::Equal, Origin::Optimized}) {
      const auto cfg = global ? glob(origin) : cost_only(FlowPlacement::MaxFlow, Admission::Lossless, std::nullopt, origin);
      const std::string where = "trial " + std::to_string(trial) + (origin == Origin::Equal ? " equal" : " optimized");
      const auto r = capacity(*g, s, t, QueryOptions {cfg, &mask, nullptr, kEdges});
      const auto& flows = *r.edge_flows;
      placed_by_origin.push_back(r.placed);
      std::vector<std::string> fails;
      const FlowBalance bal = flow_clauses(*g, view, flows, roster, t, &mask, r.placed, fails);
      EXPECT_TRUE(fails.empty()) << where << ": " << (fails.empty() ? "" : fails.front());
      for (std::size_t e = 0; e < m; ++e) {
        if (!allowed[e]) EXPECT_EQ(flows[e], 0) << where << " edge " << e;
      }
      std::vector<Amount> supply, inflow;
      for (NodeId x : roster) supply.push_back(static_cast<Amount>(bal.supply[static_cast<std::size_t>(x)]));
      for (NodeId x : t) inflow.push_back(static_cast<Amount>(bal.inflow[static_cast<std::size_t>(x)]));
      Int128 cost = 0;
      for (std::size_t e = 0; e < m; ++e) cost += static_cast<Int128>(g->cost()[e]) * flows[e];
      EXPECT_EQ(*r.total_cost, (Counter128 {0, static_cast<std::uint64_t>(cost)})) << where;
      const auto k = roster.size();
      if (origin == Origin::Equal) {
        if (k > 0) {
          // The round-robin quotas written out (not the library's round_robin): T / k each, one more
          // for the first T % k roster sources.
          auto quotas = [&](Amount total) {
            std::vector<Amount> q(k);
            for (std::size_t i = 0; i < k; ++i) {
              q[i] = total / static_cast<Amount>(k) + (static_cast<Amount>(i) < total % static_cast<Amount>(k) ? 1 : 0);
            }
            return q;
          };
          EXPECT_EQ(supply, quotas(r.placed)) << where;
          const auto next = quotas(r.placed + 1);
          std::vector<std::optional<Amount>> caps(next.begin(), next.end());
          EXPECT_LT(max_flow(*g, view, allowed, roster, caps, t), r.placed + 1) << where;
        }
        EXPECT_FALSE(r.cost_distribution.has_value()) << where;
      } else {
        const std::vector<std::optional<Amount>> free(k);
        EXPECT_EQ(max_flow(*g, view, allowed, roster, free, t), r.placed) << where;
      }
      EXPECT_FALSE(negative_cycle(*g, view, allowed, flows, roster, supply, t, inflow, origin == Origin::Optimized)) << where;
      // A row for every input source; the supplies are the rows' admitted amounts.
      ASSERT_EQ(r.by_source.size(), std::set<NodeId>(s.begin(), s.end()).size()) << where;
      for (const auto& row : r.by_source) {
        const auto it = std::find(roster.begin(), roster.end(), row.source);
        EXPECT_EQ(lo(row.admitted), it == roster.end() ? 0 : supply[static_cast<std::size_t>(it - roster.begin())]) << where;
      }
      if (r.placed > 0) ++positive;
    }
    if (placed_by_origin[0] != placed_by_origin[1]) ++differ;
  }
  EXPECT_GT(positive, 150);
  EXPECT_GT(differ, 15);
}

}  // namespace

// G16 (property).
TEST(GlobalEqual, Certificates) { certify_controllers(true); }

// M8 (property).
TEST(CostOnlyMaxFlow, Certificates) { certify_controllers(false); }

// total_cost is exact beyond 64 bits, for both controllers under both origins (the values as
// hexadecimal halves {hi, lo}). One link of capacity = cost = 2^40: 2^80. Two sources, 0 -> 2
// (capacity 2^40, cost 2^40) and 1 -> 2 (capacity 2^39, cost 2^41): optimized carries both links
// full, 2^80 + 2^80 = 2^81; equal stops at T = 2^40 + 1 (shares 2^39 + 1 and 2^39; one more unit
// would give source 1 more than 2^39): (2^39 + 1) * 2^40 + 2^39 * 2^41 = 3 * 2^79 + 2^40.
// A chain 0 -> 1 -> 2 of two links of capacity 2^32 - 1 and cost 2^32 + 1: each product is
// 2^64 - 1, so the second one carries into the high half: 2^65 - 2.
TEST(ResultShape, TotalCostBeyond64Bits) {
  const Amount p40 = Amount {1} << 40, p39 = Amount {1} << 39, p41 = Amount {1} << 41;
  const Amount lo32 = (Amount {1} << 32) - 1, hi32 = (Amount {1} << 32) + 1;
  struct Want {
    std::vector<InEdge> edges;
    std::vector<NodeId> s;
    NodeId t;
    Counter128 optimized, equal;
  };
  const std::vector<Want> wants = {
      {{{0, 1, p40, p40}}, {0}, 1, {0x10000, 0}, {0x10000, 0}},
      {{{0, 2, p40, p40}, {1, 2, p39, p41}}, {0, 1}, 2, {0x20000, 0}, {0x18000, Amount {1} << 40}},
      {{{0, 1, lo32, hi32}, {1, 2, lo32, hi32}}, {0}, 2, {1, 0xfffffffffffffffeULL}, {1, 0xfffffffffffffffeULL}},
  };
  for (const auto& w : wants) {
    const auto g = graph(w.edges);
    for (const auto origin : {Origin::Optimized, Origin::Equal}) {
      const Counter128 want = origin == Origin::Optimized ? w.optimized : w.equal;
      for (const auto& cfg :
           {glob(origin), cost_only(FlowPlacement::MaxFlow, Admission::Lossless, std::nullopt, origin)}) {
        Case k(w.edges, w.s, {w.t}, cfg);
        EXPECT_EQ(*k.place(std::nullopt).total_cost, want);
        EXPECT_EQ(*capacity(*g, w.s, v({w.t}), QueryOptions {cfg}).total_cost, want);
        const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{w.s, {w.t}}, {w.s, {w.t}}};
        for (const auto& r : batch_capacity(*g, pairs, QueryOptions {cfg}, {}, 2)) EXPECT_EQ(*r.total_cost, want);
      }
    }
  }
}

// G12 to G15: min_cut is its own output, for the optimized Global only.
TEST(GlobalMinCut, Flag) {
  for (const auto& [edges, t] : {std::pair {kB1, NodeId {3}}, std::pair {kB3Z, NodeId {2}}}) {
    Case k(edges, {0, 1}, {t}, glob());
    EXPECT_EQ(*k.place(std::nullopt, kCut).min_cut, (std::vector<EdgeId> {0, 1}));
  }
  Case zero(kZero, {0}, {1}, glob());
  const auto z = zero.place(5, kCut);
  EXPECT_EQ(z.termination, Termination::Exhausted);
  EXPECT_EQ(*z.min_cut, (std::vector<EdgeId> {0}));
  Case bottleneck(kB1, {0, 1}, {3}, glob());
  EXPECT_FALSE(bottleneck.place(std::nullopt, OutputOptions {false, false, true, false}).min_cut.has_value());
  Case equal(kB1, {0, 1}, {3}, glob(Origin::Equal));
  EXPECT_FALSE(equal.place(std::nullopt, kCut).min_cut.has_value());
  Case below(kB1, {0, 1}, {3}, glob());
  EXPECT_FALSE(below.place(5, kCut).min_cut.has_value());
  EXPECT_FALSE(below.place(0, kCut).min_cut.has_value());
  Case empty({{1, 0, 5, 1}}, {0}, {1}, glob());
  const auto e = empty.place(5, kCut);
  EXPECT_EQ(e.termination, Termination::NoRoute);
  EXPECT_EQ(*e.min_cut, (std::vector<EdgeId> {}));
  EXPECT_EQ(*e.cost_distribution, (std::vector<std::pair<Cost, Amount>> {}));
  EXPECT_EQ(*e.total_cost, c(0));
}

// ---------------------------------------------------------------- terminations and the result (T, S rows)

// T1 to T5.
TEST(Terminations, Order) {
  for (const auto& cfg : {glob(), cost_only(), cost_only(FlowPlacement::EqualBalanced, Admission::Lossy)}) {
    Case k({{1, 0, 5, 1}}, {0}, {1}, cfg);
    EXPECT_EQ(k.place(0).termination, Termination::Complete);
    const auto r = k.place(5);
    EXPECT_EQ(r.termination, Termination::NoRoute);
    if (cfg.admission == Admission::Lossy) {
      // O2: nothing enters: the offer is all source loss, every row zero.
      EXPECT_EQ(lo(r.offered), 5);
      EXPECT_EQ(lo(r.source_loss), 5);
      EXPECT_EQ(r.call_source_loss, 5);
      EXPECT_EQ(rows(r), (std::vector<Row> {{0, 0, 0, 0, 0, 0}}));
    } else {
      EXPECT_EQ(lo(r.unadmitted), 5);
    }
  }
  Case w(kZero, {0}, {1}, cost_only(FlowPlacement::Weighted, Admission::Lossy));
  const auto wr = w.place(2);
  EXPECT_EQ(wr.termination, Termination::NoActiveWeight);
  EXPECT_EQ(wr.call_source_loss, 2);
  EXPECT_EQ(lo(wr.offered), 2);
  Case wl(kZero, {0}, {1}, cost_only(FlowPlacement::Weighted));
  const auto wlr = wl.place(2);
  EXPECT_EQ(wlr.termination, Termination::NoActiveWeight);
  EXPECT_EQ(lo(wlr.unadmitted), 2);
  Case lossy(kB1, {0, 1}, {3}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy));
  EXPECT_EQ(lossy.place(12).termination, Termination::Complete);
  Case unlimited(kPool, {0, 1}, {3}, cost_only());
  EXPECT_EQ(unlimited.place(std::nullopt).termination, Termination::Exhausted);
}

// T6: a controller whose maximum is 0 over existing routes is exhausted, not without a route.
TEST(Terminations, ControllerAtZero) {
  Case g(kZero, {0}, {1}, glob());
  EXPECT_EQ(g.place(5).termination, Termination::Exhausted);
  Case m(kSpf, {0}, {1}, cost_only(FlowPlacement::MaxFlow));
  EXPECT_EQ(m.place(5).termination, Termination::Exhausted);
}

// S2: Lossless accounting (G1).
TEST(ResultShape, LosslessAccounting) {
  Case k(kB1, {0, 1}, {3}, glob(Origin::Equal));
  const auto r = k.place(12);
  EXPECT_EQ(lo(r.requested), 12);
  EXPECT_FALSE(r.offered.has_value());
  EXPECT_EQ(lo(r.admitted), 4);
  EXPECT_EQ(r.delivered, 4);
  EXPECT_EQ(lo(r.unadmitted), 8);
  EXPECT_EQ(lo(r.source_loss), 0);
  EXPECT_TRUE(r.drops.empty() && r.call_drops.empty());
  EXPECT_EQ(r.call_source_loss, 0);
  EXPECT_EQ(r.delivery_delta, 4);
  EXPECT_EQ(rows(r), (std::vector<Row> {{0, kNone, 2, 2, 0, 0}, {1, kNone, 2, 2, 0, 0}}));
}

// S3, S4: Lossy accounting (P12, I2): state versus call.
TEST(ResultShape, LossyAccounting) {
  Case k(kChain, {0, 1}, {2}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy));
  const auto r = k.place(10);
  EXPECT_EQ(lo(r.requested), 10);
  EXPECT_EQ(lo(r.offered), 10);
  EXPECT_EQ(lo(r.admitted), 10);
  EXPECT_EQ(r.delivered, 3);
  EXPECT_EQ(lo(r.unadmitted), 0);
  EXPECT_EQ(lo(r.source_loss), 0);
  EXPECT_EQ(r.drops, (Counters {{0, c(7)}}));
  Case inc(kInc, {0, 1}, {2}, cost_only(FlowPlacement::EqualBalanced, Admission::Lossy));
  (void)inc.place(3);
  const auto b = inc.place(1);
  EXPECT_TRUE(b.call_drops.empty());
  EXPECT_EQ(b.drops, (Counters {{0, c(1)}}));
}

// S5, S6: a row for every input source; total_cost only for the controllers.
TEST(ResultShape, RowsForEveryInputSource) {
  Case k(kRows, {0, 1, 4}, {3}, cost_only(), {4});
  EXPECT_EQ(rows(k.place(3)),
            (std::vector<Row> {{0, kNone, 3, 3, 0, 0}, {1, kNone, 0, 0, 0, 0}, {4, kNone, 0, 0, 0, 0}}));
  const auto g = graph(kB1);
  Pinned pin;
  pin.bundles = {dag_of(*g, {0})};
  for (const auto& cfg : {cost_only(), cost_only(FlowPlacement::Weighted),
                          FlowPolicyConfig {.routing = CapacityAware {}, .placement = FlowPlacement::EqualBalanced}}) {
    FlowLedger fg(g);
    EXPECT_FALSE(FlowPolicy(fg, v({0, 1}), v({3}), cfg).place(1).total_cost.has_value());
  }
  FlowLedger fg(g);
  EXPECT_FALSE(FlowPolicy(fg, v({0}), v({3}), FlowPolicyConfig {.routing = pin, .placement = FlowPlacement::MaxFlow})
                   .place(1)
                   .total_cost.has_value());
  Case mf(kB1, {0, 1}, {3}, cost_only(FlowPlacement::MaxFlow));
  EXPECT_EQ(*mf.place(12).total_cost, c(32));
}

// S7: CapacityAware keeps its per-source computation (offers 6 and 6).
TEST(ResultShape, CapacityAwareKeepsItsComputation) {
  Case k(kB1, {0, 1}, {3}, FlowPolicyConfig {.routing = CapacityAware {}, .placement = FlowPlacement::EqualBalanced});
  const auto r = k.place(12);
  EXPECT_EQ(r.placed, 8);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  EXPECT_EQ(lo(r.requested), 12);
  EXPECT_EQ(lo(r.unadmitted), 4);
  EXPECT_EQ(rows(r), (std::vector<Row> {{0, kNone, 2, 2, 0, 0}, {1, kNone, 6, 6, 0, 0}}));
}

// S8, S9 (D-zero-service): Pinned over a zero edge in the new shape.
TEST(ResultShape, PinnedInTheNewShape) {
  const auto g = graph(kZero);
  Pinned pin;
  pin.bundles = {dag_of(*g, {0})};
  FlowLedger fg(g);
  FlowPolicy w(fg, v({0}), v({1}), FlowPolicyConfig {pin, FlowPlacement::Weighted, Admission::Lossy});
  const auto r = w.place(2);
  EXPECT_EQ(r.placed, 0);
  EXPECT_EQ(r.termination, Termination::NoActiveWeight);
  EXPECT_EQ(r.call_source_loss, 2);
  EXPECT_EQ(lo(r.source_loss), 2);
  EXPECT_EQ(lo(r.offered), 2);
  EXPECT_EQ(rows(r), (std::vector<Row> {{0, 2, 2, 0, 2, 0}}));
  FlowLedger fg2(g);
  FlowPolicy eb(fg2, v({0}), v({1}), FlowPolicyConfig {pin, FlowPlacement::EqualBalanced, Admission::Lossy});
  const auto e = eb.place(2);
  EXPECT_EQ(e.drops, (Counters {{0, c(2)}}));
  EXPECT_EQ(e.termination, Termination::Complete);
}

// ---------------------------------------------------------------- evaluators (V rows, Q7)

// V1 to V4.
TEST(Evaluators, CapacityUnderBothOrigins) {
  const auto g = graph(kB1);
  const std::vector<std::pair<FlowPolicyConfig, Amount>> wants = {
      {glob(Origin::Equal), 4},
      {glob(), 12},
      {cost_only(FlowPlacement::MaxFlow), 12},
      {cost_only(FlowPlacement::MaxFlow, Admission::Lossless, std::nullopt, Origin::Equal), 4},
      {cost_only(), 4},
  };
  for (const auto& [cfg, placed] : wants) EXPECT_EQ(capacity(*g, v({0, 1}), v({3}), QueryOptions {cfg}).placed, placed);
  EXPECT_EQ(capacity(*graph(kPool), v({0, 1}), v({3}), QueryOptions {cost_only()}).placed, 2);
  const auto r = capacity(*graph(kPriv), v({0, 1}), v({2}),
                          QueryOptions {cost_only(FlowPlacement::EqualBalanced, Admission::Lossy)});
  EXPECT_EQ(r.delivered, 10);
  EXPECT_EQ(lo(r.offered), 10);
  EXPECT_FALSE(r.requested.has_value());
  EXPECT_EQ(r.termination, Termination::Exhausted);
}

// V5.
TEST(Evaluators, BatchEqual) {
  const auto g = graph(kB1);
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{{0, 1}, {3}}, {{0}, {3}}, {{1}, {3}}};
  std::vector<Amount> got;
  for (const auto& r : batch_capacity(*g, pairs, QueryOptions {glob(Origin::Equal)}, {}, 2)) got.push_back(r.placed);
  EXPECT_EQ(got, (std::vector<Amount> {4, 2, 10}));
}

// Q7: the pooled single-path rule per pair, the lowest failing pair.
TEST(Evaluators, BatchQ12PerPair) {
  const auto g = graph(kPool);
  const QueryOptions opts {cost_only(FlowPlacement::EqualBalanced, Admission::Lossless, kPreferPath)};
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> bad = {{{0}, {3}}, {{0, 1}, {3}}};
  EXPECT_EQ(error_of([&] { (void)batch_capacity(*g, bad, opts, {}, 1); }), kQ12);
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> good = {{{0}, {3}}, {{1}, {3}}};
  EXPECT_EQ(batch_capacity(*g, good, opts, {}, 1).size(), 2U);
}

// V6 (B3-delete-zero): every allowed edge is a candidate under the equal origin.
TEST(Evaluators, SensitivityEqualGlobal) {
  for (int workers : {1, 2}) {
    const auto r = sensitivity(*graph(kB3Z), v({0, 1}), v({2}), QueryOptions {glob(Origin::Equal)}, workers);
    EXPECT_EQ(r.baseline.placed, 0);
    std::vector<std::pair<EdgeId, Amount>> got;
    for (const auto& e : r.entries) got.emplace_back(e.edge, e.loss);
    EXPECT_EQ(got, (std::vector<std::pair<EdgeId, Amount>> {{0, -10}, {1, 0}}));
  }
}

// V7: pooled CostOnly: deleting e0 takes source 0 out of the roster.
TEST(Evaluators, SensitivityPooled) {
  const auto r = sensitivity(*graph(kB1), v({0, 1}), v({3}), QueryOptions {cost_only()}, 1);
  EXPECT_EQ(r.baseline.placed, 4);
  std::vector<std::pair<EdgeId, Amount>> got;
  for (const auto& e : r.entries) got.emplace_back(e.edge, e.loss);
  EXPECT_EQ(got, (std::vector<std::pair<EdgeId, Amount>> {{0, -6}, {1, 2}, {2, 2}}));
}
