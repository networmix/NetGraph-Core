// Edges of capacity 0 stay in routing (docs/semantics.md): searches, cached supports, rosters and
// Stretch baselines read no capacity (only an explicit residual filters), an EqualBalanced split
// still gives a zero-capacity member its round-robin share (dropped under Lossy, blocking under
// Lossless), a Weighted split gives it weight 0 and never routes around it, sensitivity lists it,
// and Global's minimum cut lists it when it leaves the sources' residual closure. Each test checks
// hand-derived values on small graphs whose EdgeIds equal their input positions (asserted). The
// graphs: S (the contract's search example, whose only shortest route crosses the zero edge): 0:
// 0 -> 1 cap 0, 1: 0 -> 2 cap 5, 2: 2 -> 1 cap 5, costs 1;
// B (two sources): 0: 0 -> 2 cap 0, 1: 1 -> 2 cap 10.
#include <gtest/gtest.h>

#include <algorithm>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/search.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

// A graph whose rows are already in EdgeId order (source, cost, destination, input).
std::shared_ptr<const Graph> ordered(std::int32_t n, const std::vector<InEdge>& rows) {
  auto g = shared(make_graph(n, rows));
  for (std::size_t i = 0; i < rows.size(); ++i) EXPECT_EQ(g->input_index()[i], static_cast<std::int32_t>(i));
  return g;
}

std::shared_ptr<const Graph> graph_s() { return ordered(3, {{0, 1, 0, 1}, {0, 2, 5, 1}, {2, 1, 5, 1}}); }
std::shared_ptr<const Graph> graph_b(Amount c0, Amount c1) { return ordered(3, {{0, 2, c0, 1}, {1, 2, c1, 1}}); }

PathSelection sel(EdgeTieBreak t = EdgeTieBreak::Deterministic, EdgePriority p = EdgePriority::Capacity,
                  bool multipath = true, bool multi_edge = true) {
  return PathSelection(t, p, multipath, multi_edge);
}

// A configuration; Global takes neither a placement nor a selection, so `k` and `s` are left out.
FlowPolicyConfig cfg(Routing r, FlowPlacement k, Admission a = Admission::Lossless, PathSelection s = sel()) {
  if (std::holds_alternative<Global>(r)) return FlowPolicyConfig {.routing = std::move(r), .admission = a};
  return FlowPolicyConfig {std::move(r), k, a, s};
}

SpfResult search(const Graph& g, std::vector<NodeId> roots, Direction d, std::optional<std::vector<NodeId>> targets,
                 const PathSelection& s = sel(), const std::vector<Amount>* residual = nullptr,
                 const Eligibility* elig = nullptr) {
  SpfOptions o(d, s);
  o.targets = std::move(targets);
  if (residual != nullptr) o.residual = std::span<const Amount>(*residual);
  o.elig = elig;
  return spf(g, roots, o);
}

std::vector<EdgeId> dag_edges(const PredDAG& dag) { return ascending_edges(dag); }

std::vector<std::vector<EdgeId>> walks(const std::vector<std::pair<Cost, Path>>& r) {
  std::vector<std::vector<EdgeId>> out;
  for (const auto& [c, p] : r) {
    std::vector<EdgeId> w;
    for (const auto& hop : p.hops()) w.push_back(hop.front());
    out.push_back(w);
  }
  return out;
}

using Drops = std::vector<std::pair<EdgeId, Amount>>;

// One placement on a fresh ledger of g.
struct Placed {
  PlacementResult r;
  std::vector<FlowRecord> records;
};
Placed place_once(const std::shared_ptr<const Graph>& g, std::vector<NodeId> s, std::vector<NodeId> d,
                  const FlowPolicyConfig& c, std::optional<Amount> v, const Eligibility* elig = nullptr,
                  OutputOptions out = {}) {
  FlowLedger fg(g);
  FlowPolicy p(fg, s, d, c, elig);
  Placed x;
  x.r = p.place(v, out);
  x.records = p.records();
  return x;
}

// T1 with N tunnels per head: one replace(offer, design) on a fresh ledger of g.
FlowPolicyConfig t1(std::int64_t n, Admission a = Admission::Lossless, std::optional<CostCeiling> ceiling = std::nullopt) {
  return FlowPolicyConfig {.routing = CapacityAware {Tunnels {n}, ceiling}, .admission = a};
}
Placed replace_once(const std::shared_ptr<const Graph>& g, std::vector<NodeId> s, std::vector<NodeId> d,
                    const FlowPolicyConfig& c, Amount offer, Amount design) {
  FlowLedger fg(g);
  FlowPolicy p(fg, s, d, c);
  Placed x;
  x.r = p.replace(offer, design);
  x.records = p.records();
  return x;
}

std::vector<std::pair<EdgeId, Amount>> entries(const SensitivityResult& s) {
  std::vector<std::pair<EdgeId, Amount>> out;
  for (const auto& e : s.entries) out.emplace_back(e.edge, e.loss);
  return out;
}

} // namespace

// SPF over the zero edge and its variants on S.
TEST(ZeroCapacity, SpfSearchesTheZeroEdgeAndOnlyAResidualFiltersIt) {
  const auto g = graph_s();
  const auto t = search(*g, {0}, Direction::Forward, std::vector<NodeId> {1});
  EXPECT_EQ(t.target_cost, 1);
  EXPECT_EQ(dag_edges(t.dag), (std::vector<EdgeId>{0}));
  const std::vector<Amount> res = {0, 5, 5};  // a residual that filters the zero edge
  const auto gated = search(*g, {0}, Direction::Forward, std::vector<NodeId> {1}, sel(), &res);
  EXPECT_EQ(gated.target_cost, 2);
  EXPECT_EQ(dag_edges(gated.dag), (std::vector<EdgeId>{1, 2}));
  const auto full = search(*g, {0}, Direction::Forward, std::nullopt);
  EXPECT_EQ(*full.costs, (std::vector<Cost>{0, 1, 1}));
  EXPECT_EQ(dag_edges(full.dag), (std::vector<EdgeId>{0, 1}));
  const auto rev = search(*g, {1}, Direction::Reverse, std::nullopt);
  EXPECT_EQ(*rev.costs, (std::vector<Cost>{1, 0, 1}));
  EXPECT_EQ(dag_edges(rev.dag), (std::vector<EdgeId>{0, 2}));
  const auto single = search(*g, {0}, Direction::Forward, std::vector<NodeId> {1}, sel(EdgeTieBreak::Deterministic,
                                                                                     EdgePriority::Capacity, false));
  EXPECT_EQ(single.target_cost, 1);
  EXPECT_EQ(dag_edges(single.dag), (std::vector<EdgeId>{0}));
  const auto mask = g->eligibility({}, std::vector<EdgeId> {0});
  const auto masked = search(*g, {0}, Direction::Forward, std::vector<NodeId> {1}, sel(), nullptr, &mask);
  EXPECT_EQ(masked.target_cost, 2);
  EXPECT_EQ(dag_edges(masked.dag), (std::vector<EdgeId>{1, 2}));
}

// KSP over the zero edge, with and without a residual (the Stretch case is in ksp_tests.cpp).
TEST(ZeroCapacity, KspListsTheZeroEdgeWalk) {
  const auto g = graph_s();
  auto r = ksp(*g, 0, 1, KspOptions(2, sel()));
  EXPECT_EQ(walks(r), (std::vector<std::vector<EdgeId>>{{0}, {1, 2}}));
  EXPECT_EQ(r[0].first, 1);
  EXPECT_EQ(r[1].first, 2);
  const std::vector<Amount> res = {0, 5, 5};
  KspOptions o(2, sel());
  o.residual = std::span<const Amount>(res);
  EXPECT_EQ(walks(ksp(*g, 0, 1, o)), (std::vector<std::vector<EdgeId>>{{1, 2}}));
}

// The parallel-member rule and the priorities without a residual.
TEST(ZeroCapacity, CollapseAndPrioritiesWithoutAResidual) {
  const auto pair = ordered(2, {{0, 1, 0, 1}, {0, 1, 3, 1}});
  auto collapse = [&](EdgeTieBreak t, EdgePriority p) {
    return dag_edges(search(*pair, {0}, Direction::Forward, std::vector<NodeId> {1}, sel(t, p, true, false)).dag);
  };
  EXPECT_EQ(collapse(EdgeTieBreak::Deterministic, EdgePriority::Capacity), (std::vector<EdgeId>{0}));
  EXPECT_EQ(collapse(EdgeTieBreak::PreferHigherPriority, EdgePriority::Capacity), (std::vector<EdgeId>{1}));
  EXPECT_EQ(collapse(EdgeTieBreak::PreferHigherPriority, EdgePriority::Headroom), (std::vector<EdgeId>{1}));
  EXPECT_EQ(collapse(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual), (std::vector<EdgeId>{1}));
  // On positive capacities 1 and 5, Headroom reads 1 everywhere (EdgeId order
  // decides); Capacity and Residual read the capacities.
  const auto d = ordered(4, {{0, 1, 1, 1}, {0, 2, 5, 1}, {1, 3, 1, 1}, {2, 3, 5, 1}});
  auto path = [&](const std::shared_ptr<const Graph>& h, EdgeTieBreak t, EdgePriority p) {
    return dag_edges(search(*h, {0}, Direction::Forward, std::vector<NodeId> {3}, sel(t, p, false)).dag);
  };
  EXPECT_EQ(path(d, EdgeTieBreak::PreferHigherPriority, EdgePriority::Headroom), (std::vector<EdgeId>{0, 2}));
  EXPECT_EQ(path(d, EdgeTieBreak::PreferHigherPriority, EdgePriority::Capacity), (std::vector<EdgeId>{1, 3}));
  EXPECT_EQ(path(d, EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual), (std::vector<EdgeId>{1, 3}));
  // Headroom 0 on the zero-capacity edge, 1 elsewhere.
  const auto z = ordered(4, {{0, 1, 0, 1}, {0, 2, 5, 1}, {1, 3, 5, 1}, {2, 3, 5, 1}});
  EXPECT_EQ(path(z, EdgeTieBreak::PreferHigherPriority, EdgePriority::Headroom), (std::vector<EdgeId>{1, 3}));
  EXPECT_EQ(path(z, EdgeTieBreak::Deterministic, EdgePriority::Headroom), (std::vector<EdgeId>{0, 2}));
  KspOptions o(2, sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Headroom));
  EXPECT_EQ(walks(ksp(*z, 0, 3, o)), (std::vector<std::vector<EdgeId>>{{1, 3}, {0, 2}}));
}

// The cached support is capacity-blind, and a policy reading it places 0 on S.
TEST(ZeroCapacity, DagCacheValuesAndTheirConsumers) {
  const auto g = graph_s();
  DagCache cache(*g);
  const auto v1 = cache.lookup(*g, v({1}), Direction::Reverse, nullptr);
  EXPECT_EQ(dag_edges(*v1), (std::vector<EdgeId>{0, 2}));
  EXPECT_EQ(cache.stats().misses, 1);
  (void)cache.lookup(*g, v({1}), Direction::Reverse, nullptr);
  EXPECT_EQ(cache.stats().hits, 1);
  FlowLedger fg(g);
  FlowPolicy p(fg, v({0}), v({1}), cost_only_config(), nullptr, &cache);
  EXPECT_EQ(p.place(5).placed, 0);
  cache.freeze();
  EXPECT_EQ(capacity(*g, v({0}), v({1}), QueryOptions {cost_only_config(), nullptr, &cache}).placed, 0);
}

// S under every model: its only shortest route crosses the zero edge.
TEST(ZeroCapacity, G8ForEveryModel) {
  const auto g = graph_s();
  const auto eb = place_once(g, {0}, {1}, cost_only_config(), 5);
  EXPECT_EQ(eb.r.placed, 0);
  EXPECT_EQ(eb.r.remaining, 5);
  EXPECT_EQ(eb.r.termination, Termination::Exhausted);
  EXPECT_TRUE(eb.r.routes_changed);
  ASSERT_EQ(eb.records.size(), 1U);
  EXPECT_EQ(eb.records[0].delivered, 0);
  EXPECT_EQ(dag_edges(*eb.records[0].bundle), (std::vector<EdgeId>{0}));
  EXPECT_EQ(eb.records[0].route_cost, 1);
  EXPECT_TRUE(eb.records[0].edges.empty());
  const auto lossy = place_once(g, {0}, {1}, cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossy), 5);
  EXPECT_EQ(lossy.r.placed, 0);
  EXPECT_EQ(lossy.r.call_drops, (Drops {{0, 5}}));
  EXPECT_EQ(lossy.r.termination, Termination::Complete);
  ASSERT_EQ(lossy.records.size(), 1U);
  EXPECT_EQ(lossy.records[0].offered, (Counter128 {0, 5}));
  EXPECT_EQ(lossy.records[0].delivered, 0);
  EXPECT_EQ(lossy.records[0].drops, (std::vector<std::pair<EdgeId, Counter128>> {{0, Counter128 {0, 5}}}));
  for (const auto a : {Admission::Lossless, Admission::Lossy}) {
    const auto w = place_once(g, {0}, {1}, cfg(CostOnly {}, FlowPlacement::Weighted, a), 5);
    EXPECT_EQ(w.r.placed, 0);
    EXPECT_TRUE(w.r.call_drops.empty());  // no member: the offer is lost at the source, never bypassed
    EXPECT_EQ(w.r.call_source_loss, a == Admission::Lossy ? 5 : 0);
    EXPECT_EQ(w.r.termination, Termination::NoActiveWeight);
    ASSERT_EQ(w.records.size(), 1U);
    EXPECT_EQ(dag_edges(*w.records[0].bundle), (std::vector<EdgeId>{0}));
    EXPECT_EQ(w.records[0].route_cost, 1);
    EXPECT_TRUE(w.records[0].edges.empty());
    if (a == Admission::Lossy) EXPECT_EQ(w.records[0].offered, (Counter128 {0, 5}));
  }
  const auto mf = place_once(g, {0}, {1}, cfg(CostOnly {}, FlowPlacement::MaxFlow), 5);
  EXPECT_EQ(mf.r.placed, 0);
  EXPECT_EQ(mf.r.termination, Termination::Exhausted);
  const auto single =
      place_once(g, {0}, {1}, cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless,
                                  sel(EdgeTieBreak::Deterministic, EdgePriority::Capacity, false)), 5);
  EXPECT_EQ(single.r.placed, 0);
  EXPECT_EQ(place_once(g, {0}, {1}, global_config(), 5).r.placed, 5);
  // T1: the request of 5 filters the zero link out (6.4): the tunnel comes up over 0 -> 2 -> 1.
  const auto tunnel = replace_once(g, {0}, {1}, t1(1), 5, 5);
  EXPECT_EQ(tunnel.r.placed, 5);
  ASSERT_EQ(tunnel.r.program->tunnels().size(), 1U);
  EXPECT_EQ(tunnel.r.program->tunnels()[0].edges, (std::optional<std::vector<EdgeId>> {{1, 2}}));
}

// A zero-capacity member first, last and deeper, and the parallel-member choice.
TEST(ZeroCapacity, EqualSplitsKeepTheZeroMemberWeightedSplitsPruneIt) {
  struct Want {
    Amount eb_lossless, eb_lossy;
    Drops drops;
    Amount weighted, maxflow;
  };
  auto check = [](const std::shared_ptr<const Graph>& g, NodeId dst, const Want& w, const std::string& where) {
    EXPECT_EQ(place_once(g, {0}, {dst}, cost_only_config(), 10).r.placed, w.eb_lossless) << where;
    const auto lossy = place_once(g, {0}, {dst}, cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossy), 10);
    EXPECT_EQ(lossy.r.placed, w.eb_lossy) << where;
    EXPECT_EQ(lossy.r.call_drops, w.drops) << where;
    for (const auto a : {Admission::Lossless, Admission::Lossy}) {
      const auto wt = place_once(g, {0}, {dst}, cfg(CostOnly {}, FlowPlacement::Weighted, a), 10);
      EXPECT_EQ(wt.r.placed, w.weighted) << where;
      EXPECT_TRUE(wt.r.call_drops.empty()) << where;
    }
    EXPECT_EQ(place_once(g, {0}, {dst}, cfg(CostOnly {}, FlowPlacement::MaxFlow), 10).r.placed, w.maxflow) << where;
  };
  const auto first = ordered(2, {{0, 1, 0, 1}, {0, 1, 10, 1}});
  check(first, 1, {0, 5, {{0, 5}}, 10, 10}, "zero first");
  check(ordered(2, {{0, 1, 10, 1}, {0, 1, 0, 1}}), 1, {1, 5, {{1, 5}}, 10, 10}, "zero last");
  const auto deep = ordered(4, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 10, 1}, {2, 3, 0, 1}});
  check(deep, 3, {1, 5, {{3, 5}}, 10, 10}, "zero deeper");
  const auto w = place_once(deep, {0}, {3}, cfg(CostOnly {}, FlowPlacement::Weighted), 10);
  ASSERT_EQ(w.records.size(), 1U);
  EXPECT_EQ(w.records[0].edges, (std::vector<std::pair<EdgeId, Amount>> {{0, 10}, {2, 10}}));
  // One member per pair; Headroom reads the working residual, 0/0 on the zero link being 0
  // (were it compared equal to every value, the lower EdgeId 0 would stay).
  auto collapsed = [&](EdgeTieBreak t, EdgePriority p) {
    return place_once(first, {0}, {1},
                      cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless, sel(t, p, true, false)), 10)
        .r.placed;
  };
  EXPECT_EQ(collapsed(EdgeTieBreak::Deterministic, EdgePriority::Capacity), 0);
  EXPECT_EQ(collapsed(EdgeTieBreak::PreferHigherPriority, EdgePriority::Capacity), 10);
  EXPECT_EQ(collapsed(EdgeTieBreak::PreferHigherPriority, EdgePriority::Headroom), 10);
  EXPECT_EQ(collapsed(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual), 10);
  const auto head = place_once(first, {0}, {1},
                               cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless,
                                   sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Headroom, false)), 10);
  EXPECT_EQ(head.r.placed, 10);  // the canonical single path by bottleneck: Headroom 1 over 0
}

// The first-hop experiment, a stale WCMP member of capacity 0, all weights zero and a link of
// service 0.
TEST(ZeroCapacity, LossyInjectionAndAllZeroWeights) {
  const auto lossy_eb = cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossy);
  const auto first = ordered(2, {{0, 1, 0, 1}, {0, 1, 10, 1}});
  EXPECT_EQ(capacity(*first, v({0}), v({1}), QueryOptions {lossy_eb}).placed, 10);
  const auto s = capacity(*graph_s(), v({0}), v({1}), QueryOptions {lossy_eb});
  EXPECT_EQ(s.placed, 0);
  EXPECT_EQ(s.termination, Termination::Exhausted);
  const auto wcmp = ordered(2, {{0, 1, 0, 1}, {0, 1, 40, 1}});
  for (const auto a : {Admission::Lossless, Admission::Lossy}) {
    EXPECT_EQ(place_once(wcmp, {0}, {1}, cfg(CostOnly {}, FlowPlacement::Weighted, a), 40).r.placed, 40);
  }
  const auto one = ordered(2, {{0, 1, 0, 1}});
  const auto w = place_once(one, {0}, {1}, cfg(CostOnly {}, FlowPlacement::Weighted, Admission::Lossy), 2);
  EXPECT_EQ(w.r.placed, 0);
  EXPECT_EQ(w.r.remaining, 2);
  EXPECT_TRUE(w.r.call_drops.empty());
  EXPECT_EQ(w.r.call_source_loss, 2);
  EXPECT_EQ(w.r.termination, Termination::NoActiveWeight);
  ASSERT_EQ(w.records.size(), 1U);
  EXPECT_EQ(w.records[0].offered, (Counter128 {0, 2}));  // 2 offered, 0 delivered, no drop: lost at the source
  EXPECT_EQ(w.records[0].delivered, 0);
  const auto d = place_once(one, {0}, {1}, lossy_eb, 2);
  EXPECT_EQ(d.r.placed, 0);
  EXPECT_EQ(d.r.call_drops, (Drops {{0, 2}}));
  EXPECT_EQ(d.r.termination, Termination::Complete);
}

// A zero-capacity first or second source, and its deletion, under the equal origin: a
// source behind a zero-capacity link stays in the roster and receives its share, which Lossless
// admits only when it is 0 (0 when it is the first source, 1 when the second) and Lossy loses.
TEST(ZeroCapacity, ZeroCapacitySourcesStayInTheRoster) {
  const auto b = graph_b(0, 10);
  const auto ll = place_once(b, {0, 1}, {2}, cost_only_config(), 10);
  EXPECT_EQ(ll.r.placed, 0);
  EXPECT_EQ(ll.r.termination, Termination::Exhausted);
  ASSERT_EQ(ll.records.size(), 2U);
  EXPECT_EQ(ll.records[0].sources, v({0}));
  EXPECT_EQ(ll.records[0].delivered, 0);
  EXPECT_EQ(dag_edges(*ll.records[0].bundle), (std::vector<EdgeId>{0}));
  EXPECT_EQ(ll.records[1].delivered, 0);
  const auto lossy_eb = cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossy);
  const auto ly = place_once(b, {0, 1}, {2}, lossy_eb, 10);
  EXPECT_EQ(ly.r.placed, 5);
  EXPECT_EQ(ly.r.call_drops, (Drops {{0, 5}}));
  EXPECT_EQ(ly.r.termination, Termination::Complete);
  const auto bp = graph_b(10, 0);
  const auto pl = place_once(bp, {0, 1}, {2}, cost_only_config(), 10);
  EXPECT_EQ(pl.r.placed, 1);
  const auto py = place_once(bp, {0, 1}, {2}, lossy_eb, 10);
  EXPECT_EQ(py.r.placed, 5);
  EXPECT_EQ(py.r.call_drops, (Drops {{1, 5}}));
  // The source or its link excluded: the roster is [1].
  const auto no_node = b->eligibility(std::vector<NodeId> {0}, {});
  const auto no_link = b->eligibility({}, std::vector<EdgeId> {0});
  for (const Eligibility* el : {&no_node, &no_link}) {
    EXPECT_EQ(place_once(b, {0, 1}, {2}, cost_only_config(), 10, el).r.placed, 10);
    EXPECT_EQ(place_once(b, {0, 1}, {2}, global_config(), 10, el).r.placed, 10);
  }
  // Global keeps a zero-delivery record of source 0 with an empty bundle.
  const auto gl = place_once(b, {0, 1}, {2}, global_config(), 10);
  EXPECT_EQ(gl.r.placed, 10);
  EXPECT_EQ(gl.r.termination, Termination::Complete);
  ASSERT_EQ(gl.records.size(), 2U);
  EXPECT_EQ(gl.records[0].delivered, 0);
  EXPECT_TRUE(dag_edges(*gl.records[0].bundle).empty());
  EXPECT_EQ(gl.records[0].route_cost, kUnreachable);
  // T1 (one tunnel per head, design 10: requests 5 and 5): head 0's link has no unreserved
  // capacity, so its tunnel is down (NoBandwidth) and the head stays in the roster with its share
  // (3.4): Lossless admits only totals that give it nothing, 0 as it is the first head; Lossy
  // delivers head 1's 5 and loses head 0's 5 at the head. One record per entry, the down one too.
  const auto ll_t1 = replace_once(b, {0, 1}, {2}, t1(1), 10, 10);
  EXPECT_EQ(ll_t1.r.placed, 0);
  EXPECT_EQ(ll_t1.r.termination, Termination::Exhausted);
  ASSERT_EQ(ll_t1.records.size(), 2U);
  EXPECT_EQ(ll_t1.records[0].entry, std::optional<std::int32_t> {0});
  EXPECT_EQ(ll_t1.r.by_tunnel[0].status, TunnelStatus::NoBandwidth);
  const auto ly_t1 = replace_once(b, {0, 1}, {2}, t1(1, Admission::Lossy), 10, 10);
  EXPECT_EQ(ly_t1.r.placed, 5);
  EXPECT_EQ(ly_t1.r.call_source_loss, 5);
  EXPECT_EQ(ly_t1.r.termination, Termination::Complete);
}

// T1's setup on X8 and the pool baseline of Stretch(1.2) (6.4).
TEST(ZeroCapacity, T1RosterAndStretchBaseline) {
  // X8: D = 4, N = 1: head 0's link has capacity 0 (NoBandwidth), head 1's is up; Lossy delivers
  // head 1's 2 and loses head 0's 2.
  const auto x8 = ordered(3, {{0, 2, 0, 1}, {1, 2, 8, 3}});
  const auto r = replace_once(x8, {0, 1}, {2}, t1(1, Admission::Lossy), 4, 4);
  EXPECT_EQ(r.r.placed, 2);
  EXPECT_EQ(r.r.call_source_loss, 2);
  EXPECT_EQ(r.r.by_tunnel[0].status, TunnelStatus::NoBandwidth);
  EXPECT_EQ(r.r.by_tunnel[1].status, TunnelStatus::Up);
  // Stretch(1.2): c0 = 5 (the zero-capacity link counts); the stored 1.2 times 5 is just
  // below 6, so the only feasible path (cost 6) exceeds the ceiling: CostCeiling, nothing
  // delivered; a factor just above 1.2 admits it.
  const auto h = ordered(2, {{0, 1, 0, 5}, {0, 1, 1, 6}});
  const auto stretched = [&](double f) { return replace_once(h, {0}, {1}, t1(1, Admission::Lossy, Stretch {f}), 1, 1); };
  const auto exact = stretched(1.2);
  EXPECT_EQ(exact.r.placed, 0);
  EXPECT_EQ(exact.r.by_tunnel[0].status, TunnelStatus::CostCeiling);
  EXPECT_EQ(exact.r.termination, Termination::Complete);  // every tunnel down: not rule (3)
  EXPECT_EQ(exact.records.size(), 1U);
  EXPECT_EQ(stretched(1.2 + 1e-9).r.placed, 1);
}

// Pinned split routes over a zero-capacity link (first or last): each
// route is its own split context, the head's share split over them by round robin. EqualBalanced
// keeps the zero link as a member, so the first unit (route 0) blocks a Lossless total: 0 when the
// zero route comes first, 1 when it comes last; Lossy offer 4 sends 2 down each route and loses
// the zero route's 2 on its link. A Weighted route over the zero link has no member (weight 0), so
// a share it receives has no route: Lossless 0 again. Global with the two links as its support (a
// controller) carries 4.
TEST(ZeroCapacity, PinnedRoutesKeepTheirZeroLinks) {
  const auto p0 = ordered(2, {{0, 1, 0, 1}, {0, 1, 4, 1}});
  const auto p1 = ordered(2, {{0, 1, 4, 1}, {0, 1, 0, 1}});
  const auto two = [&](const std::shared_ptr<const Graph>& g, std::optional<FlowPlacement> k,
                       Admission a = Admission::Lossless) {
    return FlowPolicyConfig {.routing = pinned_dags(*g, 0, 1, {dag_of(*g, {0}), dag_of(*g, {1})}, k), .admission = a};
  };
  const auto first = place_once(p0, {0}, {1}, two(p0, std::nullopt), std::nullopt);
  EXPECT_EQ(first.r.placed, 0);
  EXPECT_EQ(first.records.size(), 2U);
  EXPECT_EQ(place_once(p1, {0}, {1}, two(p1, std::nullopt), std::nullopt).r.placed, 1);
  const auto lossy = place_once(p0, {0}, {1}, two(p0, std::nullopt, Admission::Lossy), 4);
  EXPECT_EQ(lossy.r.placed, 2);
  EXPECT_EQ(lossy.r.call_drops, (Drops {{0, 2}}));
  EXPECT_EQ(place_once(p0, {0}, {1}, two(p0, FlowPlacement::Weighted), std::nullopt).r.placed, 0);
  const auto weighted_lossy = place_once(p0, {0}, {1}, two(p0, FlowPlacement::Weighted, Admission::Lossy), 4);
  EXPECT_EQ(weighted_lossy.r.placed, 2);
  EXPECT_EQ(weighted_lossy.r.call_source_loss, 2);
  EXPECT_EQ(place_once(p0, {0}, {1}, FlowPolicyConfig {.routing = Global {dag_of(*p0, {0, 1})}}, std::nullopt).r.placed,
            4);
  const auto alone =
      place_once(p0, {0}, {1}, FlowPolicyConfig {.routing = pinned_dags(*p0, 0, 1, {dag_of(*p0, {0})})}, 3);
  EXPECT_EQ(alone.r.placed, 0);
  EXPECT_EQ(alone.r.termination, Termination::Exhausted);
  ASSERT_EQ(alone.records.size(), 1U);
  EXPECT_EQ(alone.records[0].delivered, 0);
}

// The plan of a Weighted split route with capacity weights. A route without a
// zero-capacity edge is its own member set; a route with one keeps only the positive edges on the
// walks that remain, its route and route cost unchanged.
TEST(ZeroCapacity, WeightedRoutePlansKeepOnlyPositiveEdges) {
  // 0: 0 -> 1 cap 0, 1: 0 -> 1 cap 4, 2: 0 -> 2 cap 3, 3: 1 -> 3 cap 5, 4: 2 -> 3 cap 0 (costs 1).
  const auto g = ordered(4, {{0, 1, 0, 1}, {0, 1, 4, 1}, {0, 2, 3, 1}, {1, 3, 5, 1}, {2, 3, 0, 1}});
  const NodeId dst[1] = {3};
  auto plan = [&](std::vector<EdgeId> edges, Cost& cost) {
    return detail::make_plan(*g, FlowPlacement::Weighted, Admission::Lossless, detail::WeightSource {}, 0, dst,
                             *dag_of(*g, std::move(edges)), &cost);
  };
  Cost cost = -1;
  const auto positive = plan({1, 3}, cost);  // no zero-capacity edge: route == members
  EXPECT_EQ(positive.route, (std::vector<EdgeId>{1, 3}));
  EXPECT_EQ(positive.members, (std::vector<EdgeId>{1, 3}));
  EXPECT_EQ(positive.weights, (std::vector<Amount>{4, 5}));
  EXPECT_EQ(cost, 2);
  const auto mixed = plan({0, 1, 2, 3, 4}, cost);  // the zero links 0 and 4 drop; so does 2 (dead end)
  EXPECT_EQ(mixed.route, (std::vector<EdgeId>{0, 1, 2, 3, 4}));
  EXPECT_EQ(mixed.members, (std::vector<EdgeId>{1, 3}));
  EXPECT_EQ(mixed.weights, (std::vector<Amount>{4, 5}));
  EXPECT_EQ(cost, 2);
  const auto zero = plan({0, 3}, cost);  // every walk crosses a zero link: no member
  EXPECT_EQ(zero.route, (std::vector<EdgeId>{0, 3}));
  EXPECT_TRUE(zero.members.empty());
  EXPECT_TRUE(zero.weights.empty());
  EXPECT_EQ(cost, 2);
}

// The evaluators.
TEST(ZeroCapacity, EvaluatorsCountZeroCapacityRoutesAndCandidates) {
  const auto s = graph_s();
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{{0}, {1}}, {{2}, {1}}};
  const auto batch = batch_capacity(*s, pairs, QueryOptions {cost_only_config()}, {}, 1);
  ASSERT_EQ(batch.size(), 2U);
  EXPECT_EQ(batch[0].placed, 0);
  EXPECT_EQ(batch[1].placed, 5);
  // CostOnly's negative loss is in flow_policy_tests.cpp; Global keeps its
  // positive-flow candidates.
  EXPECT_EQ(entries(sensitivity(*s, v({0}), v({1}), QueryOptions {global_config()}, 1)),
            (std::vector<std::pair<EdgeId, Amount>> {{1, 5}, {2, 5}}));
  const auto b = graph_b(0, 10);
  // CostOnly EB under the equal origin (the zero-capacity source first): the baseline is 0; deleting the zero link
  // takes source 0 out of the roster and raises the total to 10 (loss -10); deleting link 1 leaves
  // source 0 alone (0).
  const auto eb = sensitivity(*b, v({0, 1}), v({2}), QueryOptions {cost_only_config()}, 1);
  EXPECT_EQ(eb.baseline.placed, 0);
  EXPECT_EQ(entries(eb), (std::vector<std::pair<EdgeId, Amount>> {{0, -10}, {1, 0}}));
  // T1 with the design 10 (X8's shape): head 0 is down and first, so the
  // baseline admits 0; deleting the zero link takes head 0 out of the roster and head 1's one
  // tunnel carries 10 (loss -10); deleting link 1 leaves head 0 alone, down (0).
  QueryOptions tq {t1(1)};
  tq.design_demand = Amount {10};
  const auto ca = sensitivity(*b, v({0, 1}), v({2}), tq, 2);
  EXPECT_EQ(ca.baseline.placed, 0);
  EXPECT_EQ(entries(ca), (std::vector<std::pair<EdgeId, Amount>> {{0, -10}, {1, 0}}));
  const auto p0 = ordered(2, {{0, 1, 0, 1}, {0, 1, 4, 1}});
  const Pinned pin = pinned_dags(*p0, 0, 1, {dag_of(*p0, {0}), dag_of(*p0, {1})});
  const auto ps = sensitivity(*p0, v({0}), v({1}), QueryOptions {{.routing = pin}}, 1);
  EXPECT_EQ(ps.baseline.placed, 0);
  EXPECT_EQ(entries(ps), (std::vector<std::pair<EdgeId, Amount>> {{0, -4}, {1, 0}}));
}

// Global's minimum cut lists the zero-capacity link leaving the closure (value unchanged);
// an excluded link is never listed.
TEST(ZeroCapacity, GlobalMinCutListsZeroCapacityLinksLeavingTheClosure) {
  const OutputOptions want {false, false, true, true};
  const auto b = graph_b(0, 10);
  const auto rb = capacity(*b, v({0, 1}), v({2}), QueryOptions {global_config(), nullptr, nullptr, want});
  EXPECT_EQ(rb.placed, 10);
  EXPECT_EQ(rb.min_cut, (std::optional<std::vector<EdgeId>> {{0, 1}}));
  const auto s = graph_s();
  const auto rs = capacity(*s, v({0}), v({1}), QueryOptions {global_config(), nullptr, nullptr, want});
  EXPECT_EQ(rs.placed, 5);
  EXPECT_EQ(rs.min_cut, (std::optional<std::vector<EdgeId>> {{0, 1}}));
  const auto no_zero = s->eligibility({}, std::vector<EdgeId> {0});
  const auto rm = capacity(*s, v({0}), v({1}), QueryOptions {global_config(), &no_zero, nullptr, want});
  EXPECT_EQ(rm.min_cut, (std::optional<std::vector<EdgeId>> {{1}}));
}
