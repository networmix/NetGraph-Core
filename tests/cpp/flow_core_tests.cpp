// The integer minimum-cost flow core through place_global_units: a first solve on a fresh graph,
// the witness graphs (hand-built examples; expected values derived by hand and from the reference flow core,
// tests/reference/policy/p1_amount/flowcore.py), the certificate in its two branches, input
// validation, determinism, and randomized graphs certified against an exhaustive optimum (tier a)
// and an independent solver (tier b), which is itself checked against the exhaustive optimum.
#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "netgraph/core/detail/flow_core.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

constexpr Amount kL62m1 = (Amount{1} << 62) - 1;

// The counters fixed by the per-tier order and the augmenting-path sequence, compared
// with the reference flow core (tests/reference/policy/p1_amount/flowcore.py): dijkstra_runs,
// dijkstra_pops, settled, dinic_phases, augmentations. arc_scans
// is the implementation's own diagnostic (the level search here first marks the vertices that
// reach t* and stops at t*'s level, so it examines a different number of arcs than the
// reference's full BFS); it is checked for determinism and presence only.
using Counters = std::vector<std::int64_t>;
Counters counters_of(const FlowCoreResult& r) {
  const auto& c = r.counters;
  EXPECT_GE(c.arc_scans, c.dijkstra_runs > 0 ? 1 : 0);
  return {c.dijkstra_runs, c.dijkstra_pops, c.settled, c.dinic_phases, c.augmentations};
}

std::vector<InEdge> chain_plus_bypasses(Cost factor = 1) {
  // Chain 0-1-2-3 of cost 1 with bypasses 0->2 and 1->3 of cost 3, unit capacities.
  return {{0, 1, 1, factor}, {1, 2, 1, factor}, {2, 3, 1, factor}, {0, 2, 1, 3 * factor}, {1, 3, 1, 3 * factor}};
}

void expect_certified(const Graph& g, const std::vector<NodeId>& s, const std::vector<NodeId>& d,
                      std::optional<Amount> cap, const FlowCoreResult& r, const std::vector<Amount>* view = nullptr,
                      const Eligibility* elig = nullptr, bool tier_a = true) {
  std::vector<Amount> full(g.capacity_units().begin(), g.capacity_units().end());
  const auto& v = view ? *view : full;
  const auto fails = certify(g, v, s, d, cap, elig, r, tier_a);
  EXPECT_TRUE(fails.empty()) << (fails.empty() ? "" : fails.front());
}

} // namespace

TEST(FlowCore, ConstructAndSolve) {
  // 0 -> 1 (cost 10), 1 -> 2 (cost 20), 0 -> 2 (cost 100); capacities 10 units each.
  const auto g = make_graph(3, {{0, 1, 10, 10}, {1, 2, 10, 20}, {0, 2, 10, 100}});
  EXPECT_EQ(g.num_nodes(), 3);
  EXPECT_EQ(g.num_edges(), 3);
  const auto r = solve(g, {0}, {2}, std::nullopt);
  EXPECT_EQ(r.placed, 20);
  ASSERT_EQ(r.cost_distribution.size(), 2u);
  EXPECT_EQ(r.cost_distribution[0], (std::pair<Cost, Amount> {30, 10}));   // 0 -> 1 -> 2 first
  EXPECT_EQ(r.cost_distribution[1], (std::pair<Cost, Amount> {100, 10}));  // then the direct edge
  expect_certified(g, {0}, {2}, std::nullopt, r);
}

TEST(FlowCore, QuantisedCapacitiesFamily) {
  // Physical 2.5q is 3 units and 0.5q is 1 unit (the quantised network).
  const double q = std::ldexp(1.0, -12);
  const Amount a = to_units(2.5 * q, q);
  const Amount b = to_units(0.5 * q, q);
  ASSERT_EQ(a, 3);
  ASSERT_EQ(b, 1);
  auto g = make_graph(3, {{0, 1, a, 1}, {1, 2, a, 1}, {0, 2, b, 3}}, q);
  auto r = solve(g, {0}, {2}, std::nullopt);
  EXPECT_EQ(r.placed, 4);
  EXPECT_EQ(r.cost_distribution, curve({{2, 3}, {3, 1}}));
  ASSERT_TRUE(r.min_cut.has_value());
  EXPECT_EQ(ids_to_input(g, *r.min_cut), (std::vector<std::int32_t>{0, 2}));
  EXPECT_EQ(r.potentials, (std::vector<Cost>{0, 2, 3}));
  EXPECT_EQ(counters_of(r), (Counters{3, 13, 12, 4, 2}));
  expect_certified(g, {0}, {2}, std::nullopt, r);
}

TEST(FlowCore, ChainPlusBypassesCancellation) {
  // A chain with bypasses: value 2 needs the middle edge cancelled.
  auto g = make_graph(4, chain_plus_bypasses());
  auto one = solve(g, {0}, {3}, 1);
  EXPECT_EQ(one.placed, 1);
  EXPECT_EQ(one.cost_distribution, curve({{3, 1}}));
  EXPECT_EQ(to_input(g, one.edge_flows), (std::vector<Amount>{1, 1, 1, 0, 0}));
  EXPECT_FALSE(one.min_cut.has_value());  // value 1 of a maximum 2: no cut
  EXPECT_EQ(counters_of(one), (Counters{1, 8, 6, 1, 1}));
  expect_certified(g, {0}, {3}, 1, one);
  auto two = solve(g, {0}, {3}, 2);
  EXPECT_EQ(two.placed, 2);
  EXPECT_EQ(two.cost_distribution, curve({{3, 1}, {5, 1}}));  // C(1) = 3, C(2) = 8
  EXPECT_EQ(to_input(g, two.edge_flows), (std::vector<Amount>{1, 0, 1, 1, 1}));
  ASSERT_TRUE(two.min_cut.has_value());
  EXPECT_EQ(ids_to_input(g, *two.min_cut), (std::vector<std::int32_t>{0, 3}));
  EXPECT_EQ(two.potentials, (std::vector<Cost>{0, 2, 3, 5}));
  EXPECT_EQ(counters_of(two), (Counters{2, 14, 12, 3, 2}));
  expect_certified(g, {0}, {3}, 2, two);
  auto inf = solve(g, {0}, {3}, std::nullopt);
  EXPECT_EQ(inf.placed, 2);
  EXPECT_EQ(inf.cost_distribution, two.cost_distribution);
  EXPECT_EQ(inf.edge_flows, two.edge_flows);
  EXPECT_EQ(inf.min_cut, two.min_cut);
  EXPECT_EQ(counters_of(inf), (Counters{3, 16, 14, 4, 2}));  // the final failed Dijkstra
  expect_certified(g, {0}, {3}, std::nullopt, inf);
}

TEST(FlowCore, ChainPlusBypassesMasksAndCostFactors) {
  // The middle edge disabled by a zero view or excluded by Eligibility: curve [(4, 2)].
  auto g = make_graph(4, chain_plus_bypasses());
  const auto inv = id_of_input(g);
  std::vector<Amount> view(g.capacity_units().begin(), g.capacity_units().end());
  view[static_cast<std::size_t>(inv[1])] = 0;
  auto zero = solve(g, {0}, {3}, std::nullopt, &view);
  EXPECT_EQ(zero.cost_distribution, curve({{4, 2}}));
  EXPECT_EQ(to_input(g, zero.edge_flows), (std::vector<Amount>{1, 0, 1, 1, 1}));
  EXPECT_EQ(zero.potentials, (std::vector<Cost>{0, 1, 3, 4}));
  expect_certified(g, {0}, {3}, std::nullopt, zero, &view);
  const std::vector<EdgeId> excluded = {inv[1]};
  auto el = g.eligibility({}, excluded);
  auto masked = solve(g, {0}, {3}, std::nullopt, nullptr, &el);
  EXPECT_EQ(masked.cost_distribution, curve({{4, 2}}));
  EXPECT_EQ(masked.edge_flows, zero.edge_flows);
  EXPECT_EQ(masked.min_cut, zero.min_cut);
  expect_certified(g, {0}, {3}, std::nullopt, masked, nullptr, &el);
  // Excluding node 1 leaves 0 -> 2 -> 3 only.
  const std::vector<NodeId> xn = {1};
  auto el_node = g.eligibility(xn, {});
  auto node_masked = solve(g, {0}, {3}, std::nullopt, nullptr, &el_node);
  EXPECT_EQ(node_masked.cost_distribution, curve({{4, 1}}));
  EXPECT_EQ(node_masked.potentials[1], kNeverReached);
  expect_certified(g, {0}, {3}, std::nullopt, node_masked, nullptr, &el_node);
  // Scaling every cost by a factor scales the curve and the potentials; the flows, the cut and
  // the operation counters are unchanged.
  for (Cost factor : {2, 7, 1000003}) {
    auto gf = make_graph(4, chain_plus_bypasses(factor));
    auto base = solve(g, {0}, {3}, std::nullopt);
    auto scaled = solve(gf, {0}, {3}, std::nullopt);
    ASSERT_EQ(scaled.cost_distribution.size(), base.cost_distribution.size());
    for (std::size_t i = 0; i < base.cost_distribution.size(); ++i) {
      EXPECT_EQ(scaled.cost_distribution[i].first, factor * base.cost_distribution[i].first);
      EXPECT_EQ(scaled.cost_distribution[i].second, base.cost_distribution[i].second);
    }
    for (std::size_t v = 0; v < base.potentials.size(); ++v) {
      EXPECT_EQ(scaled.potentials[v], factor * base.potentials[v]);
    }
    EXPECT_EQ(scaled.edge_flows, base.edge_flows);
    EXPECT_EQ(scaled.min_cut, base.min_cut);
    EXPECT_EQ(counters_of(scaled), counters_of(base));
    expect_certified(gf, {0}, {3}, std::nullopt, scaled);
  }
}

TEST(FlowCore, EqualCostTierNeedsTwoDinicPhases) {
  // Equal-cost walks of several hop counts in one tier:
  // C = {0: 0, 1: 3, 2: 6} in one tier of cost 3 drained by two successful level builds (plus the
  // tier's last unsuccessful one).
  auto g = make_graph(4, {{0, 3, 1, 3}, {0, 1, 1, 1}, {1, 2, 1, 1}, {2, 3, 1, 1}});
  auto r = solve(g, {0}, {3}, std::nullopt);
  EXPECT_EQ(r.placed, 2);
  EXPECT_EQ(r.cost_distribution, curve({{3, 2}}));
  EXPECT_EQ(r.counters.dijkstra_runs, 2);
  EXPECT_EQ(r.counters.dinic_phases, 3);
  EXPECT_EQ(counters_of(r), (Counters{2, 8, 8, 3, 2}));
  expect_certified(g, {0}, {3}, std::nullopt, r);
}

TEST(FlowCore, CyclicAdmissibleNeedsUpdateBeforeDrain) {
  // Two parallel cost-1 links. Potentials are updated by
  // the first Dijkstra before the drain (on the zero potentials no arc has reduced cost 0), and
  // the zero-reduced-cost residual cycle 0 -> 1 -> 0 left after the first push is drained in
  // the same tier.
  auto g = make_graph(2, {{0, 1, 1, 1}, {0, 1, 1, 1}});
  auto r = solve(g, {0}, {1}, 2);
  EXPECT_EQ(r.placed, 2);
  EXPECT_EQ(r.cost_distribution, curve({{1, 2}}));
  EXPECT_EQ(counters_of(r), (Counters{1, 4, 4, 1, 2}));
  expect_certified(g, {0}, {1}, 2, r);
  // A partial drain: cap 1 stops inside the tier; the destinations stay connected.
  auto p = solve(g, {0}, {1}, 1);
  EXPECT_EQ(p.placed, 1);
  EXPECT_EQ(p.cost_distribution, curve({{1, 1}}));
  EXPECT_FALSE(p.min_cut.has_value());
  EXPECT_EQ(counters_of(p), (Counters{1, 4, 4, 1, 1}));
  expect_certified(g, {0}, {1}, 1, p);
}

TEST(FlowCore, FailedFinalRunRaisesPotentials) {
  // A failed final run: the side branch 3 is reached at 8 by the first run and at 10 by
  // the failed final run; dual feasibility is claimed inside the closure {0, 3} only.
  auto g = make_graph(4, {{0, 1, 1, 1}, {1, 2, 1, 1}, {1, 3, 1, 7}, {0, 3, 1, 10}});
  auto r = solve(g, {0}, {2}, std::nullopt);
  EXPECT_EQ(r.placed, 1);
  EXPECT_EQ(r.cost_distribution, curve({{2, 1}}));
  EXPECT_EQ(r.potentials, (std::vector<Cost>{0, 1, 2, 10}));
  ASSERT_TRUE(r.min_cut.has_value());
  EXPECT_EQ(ids_to_input(g, *r.min_cut), (std::vector<std::int32_t>{0}));
  EXPECT_EQ(counters_of(r), (Counters{2, 10, 9, 2, 1}));
  expect_certified(g, {0}, {2}, std::nullopt, r);
}

TEST(FlowCore, FreshSolveUsesReleasedCheapEdge) {
  // A solve on the view [1, 1, 1] (as a re-placement sees it: the residual plus the owner's own
  // allocation) at value 1 takes the cheap link (path cost 2), not the cost-3 link (path cost 4).
  auto g = make_graph(3, {{0, 1, 1, 1}, {0, 1, 2, 3}, {1, 2, 2, 1}});
  const std::vector<Amount> view = to_ids(g, std::vector<Amount>{1, 1, 1});
  auto r = solve(g, {0}, {2}, 1, &view);
  EXPECT_EQ(to_input(g, r.edge_flows), (std::vector<Amount>{1, 0, 1}));
  EXPECT_EQ(r.cost_distribution, curve({{2, 1}}));
  std::vector<Amount> residual(view.size());
  for (std::size_t e = 0; e < view.size(); ++e) residual[e] = view[e] - r.edge_flows[e];
  EXPECT_EQ(to_input(g, residual), (std::vector<Amount>{0, 1, 0}));
  expect_certified(g, {0}, {2}, 1, r, &view);
  EXPECT_TRUE(same_payload(r, solve(g, {0}, {2}, 1, &view)));
}

TEST(FlowCore, AttainedValueUnderCapAndInfinity) {
  // The attained value under a cap, without one and at the internal cap, on one edge of capacity 3 and cost 7.
  auto g = make_graph(2, {{0, 1, 3, 7}});
  auto c2 = solve(g, {0}, {1}, 2);
  EXPECT_EQ(c2.placed, 2);
  EXPECT_FALSE(c2.min_cut.has_value());
  EXPECT_EQ(counters_of(c2), (Counters{1, 4, 4, 1, 1}));
  expect_certified(g, {0}, {1}, 2, c2);
  auto c3 = solve(g, {0}, {1}, 3);
  EXPECT_EQ(c3.placed, 3);
  ASSERT_TRUE(c3.min_cut.has_value());
  EXPECT_EQ(*c3.min_cut, (std::vector<EdgeId>{0}));
  EXPECT_EQ(counters_of(c3), (Counters{1, 4, 4, 1, 1}));  // the cap stops before a failed run
  expect_certified(g, {0}, {1}, 3, c3);
  auto inf = solve(g, {0}, {1}, std::nullopt);
  EXPECT_EQ(inf.placed, 3);
  EXPECT_EQ(inf.cost_distribution, curve({{7, 3}}));
  EXPECT_EQ(counters_of(inf), (Counters{2, 6, 6, 2, 1}));
  expect_certified(g, {0}, {1}, std::nullopt, inf);
  // The internal cap 2^63 - 2: the solve runs to disconnection (no clamp).
  auto big = solve(g, {0}, {1}, kMaxFlowCap);
  EXPECT_TRUE(same_payload(big, inf));
  expect_certified(g, {0}, {1}, kMaxFlowCap, big);
}

TEST(FlowCore, ZeroCapIsThePrescribedZeroPayload) {
  // A zero cap: no flow, no cut, every potential the never-reached sentinel, no work counted.
  auto g = make_graph(4, {{0, 1, 2, 1}, {1, 2, 2, 2}, {3, 2, 1, 1}});
  auto r = solve(g, {0}, {2}, 0);
  EXPECT_EQ(r.placed, 0);
  EXPECT_TRUE(r.cost_distribution.empty());
  EXPECT_EQ(r.edge_flows, (std::vector<Amount>{0, 0, 0}));
  EXPECT_FALSE(r.min_cut.has_value());
  EXPECT_EQ(r.potentials, (std::vector<Cost>(4, kNeverReached)));
  EXPECT_EQ(counters_of(r), (Counters{0, 0, 0, 0, 0}));
  expect_certified(g, {0}, {2}, 0, r);
}

TEST(FlowCore, PositiveCapOnZeroView) {
  // A positive cap runs one failed Dijkstra even when
  // the view gives the sources no out-capacity; roots at 0, min_cut [].
  auto g = make_graph(3, {{0, 1, 5, 1}, {1, 2, 5, 1}});
  const std::vector<Amount> view = to_ids(g, std::vector<Amount>{0, 5});
  auto r = solve(g, {0}, {2}, 4, &view);
  EXPECT_EQ(r.placed, 0);
  EXPECT_TRUE(r.cost_distribution.empty());
  ASSERT_TRUE(r.min_cut.has_value());
  EXPECT_TRUE(r.min_cut->empty());
  EXPECT_EQ(r.potentials, (std::vector<Cost>{0, kNeverReached, kNeverReached}));
  EXPECT_EQ(counters_of(r), (Counters{1, 2, 2, 0, 0}));
  expect_certified(g, {0}, {2}, 4, r, &view);
  // An all-zero view behaves the same way.
  const std::vector<Amount> zeros(2, 0);
  auto z = solve(g, {0}, {2}, 1, &zeros);
  EXPECT_EQ(z.placed, 0);
  ASSERT_TRUE(z.min_cut.has_value());
  EXPECT_TRUE(z.min_cut->empty());
  EXPECT_EQ(z.potentials[0], 0);
  expect_certified(g, {0}, {2}, 1, z, &zeros);
}

TEST(FlowCore, LegalMaximumUnits) {
  // Terminal saturation: one edge at 2^62 - 1 units; the closure is seeded at the real
  // sources, so the cut holds the real edge (virtual arcs never saturate).
  auto g = make_graph(2, {{0, 1, kL62m1, 1}});
  auto r = solve(g, {0}, {1}, std::nullopt);
  EXPECT_EQ(r.placed, kL62m1);
  EXPECT_EQ(r.cost_distribution, curve({{1, kL62m1}}));
  ASSERT_TRUE(r.min_cut.has_value());
  EXPECT_EQ(*r.min_cut, (std::vector<EdgeId>{0}));
  EXPECT_EQ(r.potentials, (std::vector<Cost>{0, 1}));
  expect_certified(g, {0}, {1}, std::nullopt, r, nullptr, nullptr, /*tier_a=*/false);
  // A sum of capacities at 2^62 - 1 split over parallel paths with large costs.
  auto h = make_graph(3, {{0, 1, kL62m1 - 10, (Amount{1} << 60)}, {0, 2, 5, 3}, {2, 1, 5, (Amount{1} << 60)}});
  auto rh = solve(h, {0}, {1}, std::nullopt);
  EXPECT_EQ(rh.placed, kL62m1 - 5);
  EXPECT_EQ(rh.cost_distribution, curve({{(Amount{1} << 60), kL62m1 - 10}, {(Amount{1} << 60) + 3, 5}}));
  expect_certified(h, {0}, {1}, std::nullopt, rh, nullptr, nullptr, /*tier_a=*/false);
  auto capped = solve(h, {0}, {1}, kL62m1 - 6);
  EXPECT_EQ(capped.placed, kL62m1 - 6);
  EXPECT_FALSE(capped.min_cut.has_value());
  expect_certified(h, {0}, {1}, kL62m1 - 6, capped, nullptr, nullptr, /*tier_a=*/false);
}

TEST(FlowCore, EndpointSetsAndVirtualTerminals) {
  // Several destinations through the virtual terminal: both destinations have
  // inflow, so both carry pi(t*) = 2.
  auto g = make_graph(3, {{0, 1, 1, 1}, {0, 2, 1, 2}});
  auto r = solve(g, {0}, {1, 2}, std::nullopt);
  EXPECT_EQ(r.placed, 2);
  EXPECT_EQ(r.cost_distribution, curve({{1, 1}, {2, 1}}));
  EXPECT_EQ(r.potentials, (std::vector<Cost>{0, 2, 2}));
  expect_certified(g, {0}, {1, 2}, std::nullopt, r);
  // Source transit is dominated: source 0 supplies nothing.
  auto t = make_graph(3, {{0, 1, 1, 1}, {1, 2, 2, 1}});
  auto rt = solve(t, {0, 1}, {2}, std::nullopt);
  EXPECT_EQ(rt.placed, 2);
  EXPECT_EQ(to_input(t, rt.edge_flows), (std::vector<Amount>{0, 2}));
  EXPECT_EQ(rt.potentials, (std::vector<Cost>{0, 0, 1}));
  expect_certified(t, {0, 1}, {2}, std::nullopt, rt);
  // Duplicates are removed.
  EXPECT_TRUE(same_payload(solve(t, {1, 0, 1}, {2, 2}, std::nullopt), rt));
}

TEST(FlowCore, ClosureAndSentinels) {
  // The cut: the closure follows reverse arcs, so the forward-only boundary is wrong.
  auto g = make_graph(5, {{0, 1, 1, 1}, {0, 2, 1, 1}, {1, 3, 1, 1}, {2, 3, 1, 1}, {3, 4, 1, 1}});
  auto r = solve(g, {0}, {4}, std::nullopt);
  EXPECT_EQ(r.placed, 1);
  ASSERT_TRUE(r.min_cut.has_value());
  EXPECT_EQ(ids_to_input(g, *r.min_cut), (std::vector<std::int32_t>{4}));
  expect_certified(g, {0}, {4}, std::nullopt, r);
  // Node 3 is never settled when the cap stops the solve with the destination connected.
  auto h = make_graph(4, {{0, 1, 2, 1}, {1, 2, 2, 2}, {3, 2, 1, 1}});
  auto rh = solve(h, {0}, {2}, 1);
  EXPECT_FALSE(rh.min_cut.has_value());
  EXPECT_EQ(rh.potentials, (std::vector<Cost>{0, 1, 3, kNeverReached}));
  expect_certified(h, {0}, {2}, 1, rh);
}

TEST(FlowCore, ExcludedEndpointsAreSkipped) {
  auto g = make_graph(3, {{0, 1, 2, 1}, {1, 2, 2, 1}});
  const std::vector<NodeId> x0 = {0};
  auto no_source = g.eligibility(x0, {});
  auto r = solve(g, {0}, {2}, std::nullopt, nullptr, &no_source);
  EXPECT_EQ(r.placed, 0);
  ASSERT_TRUE(r.min_cut.has_value());
  EXPECT_TRUE(r.min_cut->empty());
  EXPECT_EQ(r.potentials, (std::vector<Cost>(3, kNeverReached)));
  expect_certified(g, {0}, {2}, std::nullopt, r, nullptr, &no_source);
  const std::vector<NodeId> x2 = {2};
  auto no_dest = g.eligibility(x2, {});
  auto rd = solve(g, {0}, {2}, std::nullopt, nullptr, &no_dest);
  EXPECT_EQ(rd.placed, 0);
  EXPECT_EQ(rd.potentials, (std::vector<Cost>{0, 1, kNeverReached}));
  expect_certified(g, {0}, {2}, std::nullopt, rd, nullptr, &no_dest);
  // One allowed source of two.
  auto rs = solve(g, {0, 1}, {2}, std::nullopt, nullptr, &no_source);
  EXPECT_EQ(rs.placed, 2);
  expect_certified(g, {0, 1}, {2}, std::nullopt, rs, nullptr, &no_source);
}

TEST(FlowCore, MaskedSolvesReachTheValueOfTheirCuts) {
  // Sources {0, 1}, destinations {3, 4} under three Eligibilities.
  auto g = make_graph(5, {{0, 2, 2, 2}, {1, 2, 1, 1}, {2, 3, 2, 3}, {2, 4, 1, 5}});
  const auto inv = id_of_input(g);
  for (auto [mask, expect] : std::vector<std::pair<std::vector<std::int32_t>, Amount>>{{{}, 3}, {{1}, 2}, {{2}, 1}}) {
    std::vector<EdgeId> excluded;
    for (auto i : mask) excluded.push_back(inv[static_cast<std::size_t>(i)]);
    auto el = g.eligibility({}, excluded);
    auto r = solve(g, {0, 1}, {3, 4}, std::nullopt, nullptr, &el);
    EXPECT_EQ(r.placed, expect);
    expect_certified(g, {0, 1}, {3, 4}, std::nullopt, r, nullptr, &el);
  }
}

TEST(FlowCore, InputValidationBeforeAnyWork) {
  auto g = make_graph(3, {{0, 1, 2, 1}, {1, 2, 2, 1}});
  auto h = make_graph(3, {{0, 1, 2, 1}, {1, 2, 2, 1}});
  const std::vector<Amount> view = {2, 2};
  const std::vector<NodeId> s = {0}, d = {2};
  auto call = [&](std::vector<Amount> v, std::vector<NodeId> src, std::vector<NodeId> dst, std::optional<Amount> cap,
                  const Eligibility* el = nullptr) {
    return error_of([&] { (void)place_global_units(g, v, src, dst, el, cap); });
  };
  EXPECT_EQ(call({2}, s, d, 1), "view has length 1, expected 2 (one entry per edge)");
  EXPECT_EQ(call({2, 3}, s, d, 1), "view[1] = 3 is outside [0, capacity_units[1] = 2]");
  EXPECT_EQ(call({-1, 2}, s, d, 1), "view[0] = -1 is outside [0, capacity_units[0] = 2]");
  EXPECT_EQ(call(view, {}, d, 1), "sources: the node set is empty");
  EXPECT_EQ(call(view, s, {}, 1), "destinations: the node set is empty");
  EXPECT_EQ(call(view, {3}, d, 1), "sources: node 3 out of node range [0, 3)");
  EXPECT_EQ(call(view, s, {-1}, 1), "destinations: node -1 out of node range [0, 3)");
  EXPECT_EQ(call(view, {0, 2}, {2}, 1), "node 2 is in both endpoint sets: sources and destinations must not overlap");
  EXPECT_EQ(call(view, s, d, -1), "cap -1 is outside [0, 2^63 - 2]");
  EXPECT_EQ(call(view, s, d, kMaxFlowCap + 1), "cap 9223372036854775807 is outside [0, 2^63 - 2]");
  auto foreign = h.eligibility();
  EXPECT_EQ(call(view, s, d, 1, &foreign), "foreign Eligibility: built for another graph");
  // The overlap is checked on the full sets, also when an Eligibility excludes the vertex.
  const std::vector<NodeId> x0 = {0};
  auto el = g.eligibility(x0, {});
  EXPECT_EQ(call(view, {0}, {0}, 1, &el), "node 0 is in both endpoint sets: sources and destinations must not overlap");
  // Order: the view before the endpoints, the sources before the destinations, the destinations
  // before the overlap, the overlap before the cap, the cap before the Eligibility.
  EXPECT_EQ(call({2}, {}, d, -1), "view has length 1, expected 2 (one entry per edge)");
  EXPECT_EQ(call(view, {}, d, -1), "sources: the node set is empty");
  EXPECT_EQ(call(view, {}, {}, -1, &foreign), "sources: the node set is empty");
  EXPECT_EQ(call(view, {0, 2}, {2, 9}, -1, &foreign), "destinations: node 9 out of node range [0, 3)");
  EXPECT_EQ(call(view, {0, 2}, {2}, -1, &foreign),
            "node 2 is in both endpoint sets: sources and destinations must not overlap");
  EXPECT_EQ(call(view, s, d, -1, &foreign), "cap -1 is outside [0, 2^63 - 2]");
}

// A view with several defects at once (one error, from the shared edge-amount check, labelled "view"):
// its length before any entry, then its first entry outside [0, capacity] in EdgeId order, before
// the endpoints, the cap and the Eligibility.
TEST(FlowCore, ViewDefectsReportTheFirstInOrder) {
  const auto g = edge_amount_graph();
  const auto other = make_graph(4, {{0, 1, 2, 1}});
  const auto foreign = other.eligibility();
  for (const auto& c : edge_amount_defects()) {
    EXPECT_EQ(error_of([&] { (void)place_global_units(g, c.amounts, {}, {}, &foreign, -1); }), "view" + c.error);
  }
  const auto bounds = edge_amounts_at_the_bounds();
  EXPECT_EQ(place_global_units(g, bounds, v({0}), v({3}), nullptr, std::nullopt).placed, 0);
}

// A solve numbers its network in int32: the N real vertices plus s* and t*, and two arc slots per
// edge, per source and per destination. Every solver takes its vertex count from
// network_vertices(), called here with the counts alone, so no network of the boundary sizes is
// allocated.
TEST(FlowCore, NetworkSizeLimits) {
  constexpr std::int32_t kMax = std::numeric_limits<std::int32_t>::max();
  constexpr std::size_t kPairs = (std::size_t {1} << 30) - 1;  // twice it, 2^31 - 2, fits int32; 2^31 does not
  constexpr std::size_t kHuge = std::numeric_limits<std::size_t>::max();
  const std::string nodes =
      "flow core: the number of nodes exceeds INT32_MAX - 2 (a solve adds a virtual source and sink)";
  const std::string slots =
      "flow core: the edges, sources and destinations of a solve number more than 2^30 - 1 together (each takes "
      "two int32 arc slots)";
  auto error = [](std::int32_t n, std::size_t e, std::size_t s, std::size_t d) {
    return error_of([&] { (void)detail::network_vertices(n, e, s, d); });
  };
  // The vertices: N + 2, up to N = INT32_MAX - 2.
  EXPECT_EQ(detail::network_vertices(0, 0, 0, 0), 2);
  EXPECT_EQ(detail::network_vertices(kMax - 2, 0, 1, 1), kMax);
  EXPECT_EQ(error(kMax - 1, 0, 1, 1), nodes);
  EXPECT_EQ(error(kMax, 0, 1, 1), nodes);
  // The slots: edges + sources + destinations up to 2^30 - 1, every term counted.
  EXPECT_EQ(detail::network_vertices(2, kPairs - 2, 1, 1), 4);
  EXPECT_EQ(error(2, kPairs - 1, 1, 1), slots);
  EXPECT_EQ(error(2, kPairs - 1, 2, 0), slots);
  EXPECT_EQ(error(2, kPairs - 1, 0, 2), slots);
  EXPECT_EQ(detail::network_vertices(2, kPairs, 0, 0), 4);
  EXPECT_EQ(detail::network_vertices(2, 0, kPairs, 0), 4);
  EXPECT_EQ(detail::network_vertices(2, 0, 0, kPairs), 4);
  EXPECT_EQ(error(2, kPairs + 1, 0, 0), slots);
  EXPECT_EQ(error(2, 0, kPairs + 1, 0), slots);
  EXPECT_EQ(error(2, 0, 0, kPairs + 1), slots);
  // Counts whose plain sum would wrap around are refused, not wrapped into range.
  EXPECT_EQ(error(2, kHuge, 2, 0), slots);
  EXPECT_EQ(error(2, 0, kHuge, 2), slots);
  EXPECT_EQ(error(2, 2, 0, kHuge), slots);
  // The node limit is reported first.
  EXPECT_EQ(error(kMax, kPairs + 1, 0, 0), nodes);
}

TEST(FlowCore, IdenticalCallsIdenticalPayloads) {
  std::mt19937 rng(7);
  int positive = 0;
  for (int trial = 0; trial < 20; ++trial) {
    std::vector<InEdge> edges;
    for (int e = 0; e < 40; ++e) {
      edges.push_back({static_cast<NodeId>(rng() % 12), static_cast<NodeId>(rng() % 12),
                       static_cast<Amount>(rng() % 9), static_cast<Cost>(1 + rng() % 6)});
    }
    auto g = make_graph(12, edges);
    auto a = solve(g, {0, 1}, {10, 11}, std::nullopt);
    auto b = solve(g, {0, 1}, {10, 11}, std::nullopt);
    EXPECT_TRUE(same_payload(a, b));
    positive += a.placed > 0;
  }
  EXPECT_GT(positive, 10);  // two empty payloads would be identical too
}

namespace {

struct RandomCase {
  Graph g;
  std::vector<NodeId> s, d;
  std::vector<Amount> view;
  std::optional<Amount> cap;
  std::vector<NodeId> xn;
  std::vector<EdgeId> xe;
};

RandomCase random_case(std::mt19937& rng, int n, int m, Amount max_cap, Cost max_cost, bool masks) {
  std::vector<InEdge> edges;
  for (int e = 0; e < m; ++e) {
    edges.push_back({static_cast<NodeId>(rng() % n), static_cast<NodeId>(rng() % n),
                     static_cast<Amount>(rng() % (max_cap + 1)), static_cast<Cost>(1 + rng() % max_cost)});
  }
  auto g = make_graph(n, edges);
  std::vector<NodeId> nodes(static_cast<std::size_t>(n));
  std::iota(nodes.begin(), nodes.end(), 0);
  std::shuffle(nodes.begin(), nodes.end(), rng);
  const int ns = 1 + static_cast<int>(rng() % static_cast<unsigned>(std::min(2, n - 1)));
  const int nd = 1 + static_cast<int>(rng() % static_cast<unsigned>(std::min(2, n - ns)));
  std::vector<NodeId> s(nodes.begin(), nodes.begin() + ns);
  std::vector<NodeId> d(nodes.begin() + ns, nodes.begin() + std::min(n, ns + nd));
  std::vector<Amount> view;
  for (Amount c : g.capacity_units()) view.push_back(rng() % 4 == 0 ? static_cast<Amount>(rng() % (c + 1)) : c);
  std::optional<Amount> cap;
  switch (rng() % 5) {
    case 0: cap = std::nullopt; break;
    case 1: cap = 0; break;
    case 2: cap = 1; break;
    case 3: cap = static_cast<Amount>(rng() % (2 * max_cap + 2)); break;
    default: cap = kMaxFlowCap; break;
  }
  std::vector<NodeId> xn;
  std::vector<EdgeId> xe;
  if (masks && rng() % 2 == 0) {
    if (rng() % 3 == 0) xn.push_back(static_cast<NodeId>(rng() % n));
    for (int k = 0; k < 2; ++k) if (m > 0 && rng() % 2 == 0) xe.push_back(static_cast<EdgeId>(rng() % m));
  }
  return {std::move(g), s, d, view, cap, xn, xe};
}

// The kinds of random trials that reach a certificate branch only by chance, counted so a
// population that stopped drawing them fails: masked (an Eligibility), stopped by a positive cap
// below the maximum (no min_cut: the destinations stay connected) and capped at zero.
struct Population {
  int masked = 0;
  int capped = 0;
  int zero_cap = 0;
  void count(const RandomCase& c, const Eligibility* elig, const FlowCoreResult& r) {
    masked += elig != nullptr ? 1 : 0;
    capped += c.cap.has_value() && *c.cap > 0 && !r.min_cut.has_value() ? 1 : 0;
    zero_cap += c.cap == Amount {0} ? 1 : 0;
  }
  void expect_at_least(int m, int c, int z) const {
    EXPECT_GE(masked, m);
    EXPECT_GE(capped, c);
    EXPECT_GE(zero_cap, z);
  }
};

} // namespace

TEST(FlowCore, RandomGraphsTierAExhaustive) {
  // Nonuniform costs, parallel edges, self-loops, partial views, endpoint sets, masks and caps;
  // every payload certified, the curve against the minimum cost of every integer value.
  std::mt19937 rng(20261004);
  int positive = 0;
  Population kinds;
  for (int trial = 0; trial < 2000; ++trial) {
    const int n = 2 + static_cast<int>(rng() % 4);
    const int m = 1 + static_cast<int>(rng() % 6);
    auto c = random_case(rng, n, m, 2, 5, true);
    std::optional<Eligibility> el;
    if (!c.xn.empty() || !c.xe.empty()) el = c.g.eligibility(c.xn, c.xe);
    const Eligibility* elp = el ? &*el : nullptr;
    auto r = place_global_units(c.g, c.view, c.s, c.d, elp, c.cap);
    const auto fails = certify(c.g, c.view, c.s, c.d, c.cap, elp, r, true);
    ASSERT_TRUE(fails.empty()) << "trial " << trial << ": " << fails.front();
    positive += r.placed > 0;
    kinds.count(c, elp, r);
  }
  EXPECT_GT(positive, 250);  // the population is not vacuous
  // 859, 68 and 443 at this seed with libc++: random_case's std::shuffle differs between standard
  // libraries, so the floors stay well below the counts.
  kinds.expect_at_least(400, 30, 200);
}

TEST(FlowCore, RandomGraphsTierBIndependent) {
  std::mt19937 rng(1816);
  int multi_tier = 0;
  Population kinds;
  for (int trial = 0; trial < 150; ++trial) {
    const int n = 6 + static_cast<int>(rng() % 20);
    const int m = n + static_cast<int>(rng() % (4 * n));
    auto c = random_case(rng, n, m, 1000, 40, true);
    std::optional<Eligibility> el;
    if (!c.xn.empty() || !c.xe.empty()) el = c.g.eligibility(c.xn, c.xe);
    const Eligibility* elp = el ? &*el : nullptr;
    auto r = place_global_units(c.g, c.view, c.s, c.d, elp, c.cap);
    const auto fails = certify(c.g, c.view, c.s, c.d, c.cap, elp, r, false);
    ASSERT_TRUE(fails.empty()) << "trial " << trial << ": " << fails.front();
    multi_tier += r.cost_distribution.size() > 1;
    kinds.count(c, elp, r);
  }
  EXPECT_GT(multi_tier, 30);
  kinds.expect_at_least(30, 15, 10);  // 64, 35 and 23 at this seed with libc++ (as above)
}

// The independent solver is tier b's only optimum (and the deletion and Global tests'), so it is
// itself checked against the exhaustive optimum: on small random graphs with partial views,
// endpoint sets and masks, its curve has the same maximum and equals the minimum cost of every
// integer value, many of them after two augmentations or more.
TEST(FlowCore, IndependentSolverEqualsTheExhaustiveOptimum) {
  // The chain with bypasses needs a cancellation: the second unit reroutes the first.
  const auto chain = make_graph(4, chain_plus_bypasses());
  const std::vector<Amount> unit(5, 1);
  EXPECT_EQ(independent_curve(chain, unit, {0}, {3}, nullptr),
            (std::vector<std::pair<Amount, Int128>> {{0, 0}, {1, 3}, {2, 8}}));
  EXPECT_EQ(exhaustive_costs(chain, unit, {0}, {3}, nullptr), (std::map<Amount, Int128> {{0, 0}, {1, 3}, {2, 8}}));
  std::mt19937 rng(2028);
  int several = 0;
  int masked = 0;
  for (int trial = 0; trial < 1000; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 3);
    const int m = 4 + static_cast<int>(rng() % 5);
    const auto c = random_case(rng, n, m, 2, 5, true);
    std::optional<Eligibility> el;
    if (!c.xn.empty() || !c.xe.empty()) el = c.g.eligibility(c.xn, c.xe);
    const Eligibility* elp = el ? &*el : nullptr;
    const auto best = exhaustive_costs(c.g, c.view, c.s, c.d, elp);
    const auto ref = independent_curve(c.g, c.view, c.s, c.d, elp);
    ASSERT_EQ(ref.back().first, best.rbegin()->first) << "trial " << trial;
    for (const auto& [x, cost] : best) ASSERT_EQ(curve_at(ref, x), cost) << "trial " << trial << " value " << x;
    several += ref.size() > 2;
    masked += elp != nullptr && ref.size() > 1;
  }
  // Not vacuous: curves of two augmentations or more, and masked cases with a positive flow (91
  // and 116 at this seed with libc++, whose std::shuffle in random_case another standard library
  // does not reproduce).
  EXPECT_GT(several, 60);
  EXPECT_GT(masked, 90);
}

TEST(FlowCore, RandomLargeUnitsTierB) {
  // Capacities near the domain bound (sum below 2^62) and costs near 2^40.
  std::mt19937_64 rng(62);
  int positive = 0;
  for (int trial = 0; trial < 40; ++trial) {
    const int n = 4 + static_cast<int>(rng() % 6);
    const int m = 3 + static_cast<int>(rng() % 12);
    std::vector<InEdge> edges;
    const Amount per = ((Amount{1} << 62) - 1) / m;
    for (int e = 0; e < m; ++e) {
      edges.push_back({static_cast<NodeId>(rng() % n), static_cast<NodeId>(rng() % n),
                       static_cast<Amount>(rng() % static_cast<std::uint64_t>(per)),
                       static_cast<Cost>(1 + rng() % (std::uint64_t{1} << 40))});
    }
    auto g = make_graph(n, edges);
    std::vector<Amount> view(g.capacity_units().begin(), g.capacity_units().end());
    auto r = place_global_units(g, view, std::vector<NodeId>{0}, std::vector<NodeId>{static_cast<NodeId>(n - 1)},
                                nullptr, std::nullopt);
    const auto fails = certify(g, view, {0}, {n - 1}, std::nullopt, nullptr, r, false);
    ASSERT_TRUE(fails.empty()) << "trial " << trial << ": " << fails.front();
    positive += r.placed > 0;
  }
  // Only a positive flow reaches the arithmetic near the bound, and many of these sparse random
  // graphs place nothing, so the trials that place something are counted.
  EXPECT_GT(positive, 10);
}

TEST(FlowCore, CertificateRejectsMutatedPayloads) {
  // The test oracle itself is not vacuous: each negative control below is rejected by the clause
  // it breaks. Every clean payload certifies first, so the rejection comes from the mutation.
  auto one = make_graph(2, {{0, 1, 3, 7}});
  std::vector<Amount> v1 = {3};
  auto under = solve(one, {0}, {1}, 1);
  ASSERT_TRUE(certify(one, v1, {0}, {1}, 1, nullptr, under, true).empty());
  // Placed 1 under cap 2 with the edge still open.
  EXPECT_TRUE(has_failure(certify(one, v1, {0}, {1}, 2, nullptr, under, true), "attained value: early stop"));
  // The same payload without a cap, in both branches: an uncapped solve must be maximal, so the
  // two clauses that prove it fail and nothing else does (the curve up to the placed unit is
  // optimal): the destination is still reachable, and 1 is not F_max = 3.
  for (const bool tier_a : {true, false}) {
    EXPECT_EQ(certify(one, v1, {0}, {1}, std::nullopt, nullptr, under, tier_a),
              (std::vector<std::string> {"uncapped solve left the destinations connected",
                                         "placed differs from min(cap, F_max)"}))
        << tier_a;
  }
  auto slope = solve(one, {0}, {1}, 2);
  ASSERT_TRUE(certify(one, v1, {0}, {1}, 2, nullptr, slope, true).empty());
  slope.cost_distribution = curve({{6, 2}});  // moment 12, the flow costs 14
  EXPECT_TRUE(has_failure(certify(one, v1, {0}, {1}, 2, nullptr, slope, true), "curve moment differs"));
  // The curve-shape clauses, each the only failure in both branches: an entry of volume 0 before
  // the real one (an empty curve segment, which the curve comparison must skip, not divide by its
  // zero length), and the one tier split in two at the same cost. Total, moment and pi(t*) hold.
  auto full = solve(one, {0}, {1}, std::nullopt);
  ASSERT_EQ(full.cost_distribution, curve({{7, 3}}));
  for (const bool tier_a : {true, false}) {
    ASSERT_TRUE(certify(one, v1, {0}, {1}, std::nullopt, nullptr, full, tier_a).empty()) << tier_a;
    auto empty_entry = full;
    empty_entry.cost_distribution = curve({{6, 0}, {7, 3}});
    EXPECT_EQ(certify(one, v1, {0}, {1}, std::nullopt, nullptr, empty_entry, tier_a),
              std::vector<std::string> {"curve volume not positive"})
        << tier_a;
    auto split = full;
    split.cost_distribution = curve({{7, 1}, {7, 2}});
    EXPECT_EQ(certify(one, v1, {0}, {1}, std::nullopt, nullptr, split, tier_a),
              std::vector<std::string> {"curve costs not strictly increasing"})
        << tier_a;
  }
  auto cut = solve(one, {0}, {1}, std::nullopt);
  ASSERT_TRUE(certify(one, v1, {0}, {1}, std::nullopt, nullptr, cut, true).empty());
  cut.min_cut = std::vector<EdgeId>{};
  EXPECT_TRUE(has_failure(certify(one, v1, {0}, {1}, std::nullopt, nullptr, cut, true), "min_cut differs"));
  auto par = make_graph(2, {{0, 1, 1, 1}, {0, 1, 1, 3}});
  std::vector<Amount> vp = {1, 1};
  auto omitted = solve(par, {0}, {1}, std::nullopt);
  ASSERT_TRUE(certify(par, vp, {0}, {1}, std::nullopt, nullptr, omitted, true).empty());
  omitted.cost_distribution = curve({{2, 2}});  // same moment, omitted breakpoint at 1
  EXPECT_TRUE(has_failure(certify(par, vp, {0}, {1}, std::nullopt, nullptr, omitted, true),
                          "curve differs from the exhaustive optimum at 1"));
  EXPECT_TRUE(has_failure(certify(par, vp, {0}, {1}, std::nullopt, nullptr, omitted, false),
                          "curve differs from the independent solver at 1"));
  auto any = make_graph(3, {{0, 1, 1, 1}, {0, 2, 1, 2}});
  std::vector<Amount> va = {1, 1};
  auto pot = solve(any, {0}, {1, 2}, std::nullopt);
  ASSERT_TRUE(certify(any, va, {0}, {1, 2}, std::nullopt, nullptr, pot, true).empty());
  pot.potentials[1] = 1;  // an inflow destination below pi(t*)
  EXPECT_TRUE(has_failure(certify(any, va, {0}, {1, 2}, std::nullopt, nullptr, pot, true),
                          "inflow destination potential differs"));
  auto sentinel = solve(one, {0}, {1}, 1);
  ASSERT_TRUE(certify(one, v1, {0}, {1}, 1, nullptr, sentinel, true).empty());
  sentinel.potentials = {kNeverReached, kNeverReached};
  EXPECT_TRUE(has_failure(certify(one, v1, {0}, {1}, 1, nullptr, sentinel, true), "closure vertex holds NEVER_REACHED"));
  auto zero = solve(one, {0}, {1}, 0);
  ASSERT_TRUE(certify(one, v1, {0}, {1}, 0, nullptr, zero, true).empty());
  zero.counters.dijkstra_runs = 1;
  EXPECT_TRUE(has_failure(certify(one, v1, {0}, {1}, 0, nullptr, zero, true), "zero cap: the prescribed zero payload"));

  // The flow clauses. Two parallel edges of capacity 1 and cost 1 carry one unit each; both
  // units on EdgeId 0 keep the value, conservation and the cost but exceed its view of 1.
  auto twin = make_graph(2, {{0, 1, 1, 1}, {0, 1, 1, 1}});
  std::vector<Amount> vt = {1, 1};
  auto over = solve(twin, {0}, {1}, std::nullopt);
  ASSERT_TRUE(certify(twin, vt, {0}, {1}, std::nullopt, nullptr, over, true).empty());
  ASSERT_EQ(over.edge_flows, (std::vector<Amount>{1, 1}));
  over.edge_flows = {2, 0};
  EXPECT_TRUE(has_failure(certify(twin, vt, {0}, {1}, std::nullopt, nullptr, over, true), "feasibility at edge 0"));
  // EdgeId 1 masked: the solve uses EdgeId 0; moving its unit to EdgeId 1 puts flow on the
  // ineligible edge at the same cost.
  const std::vector<EdgeId> second = {1};
  const auto mask = twin.eligibility({}, second);
  auto masked = solve(twin, {0}, {1}, std::nullopt, nullptr, &mask);
  ASSERT_TRUE(certify(twin, vt, {0}, {1}, std::nullopt, &mask, masked, true).empty());
  ASSERT_EQ(masked.edge_flows, (std::vector<Amount>{1, 0}));
  masked.edge_flows = {0, 1};
  EXPECT_TRUE(has_failure(certify(twin, vt, {0}, {1}, std::nullopt, &mask, masked, true),
                          "flow on an ineligible edge 1"));
  // The chain 0 -> 1 -> 2 carries 2 units; one unit less out of node 1 leaves it unbalanced.
  auto chain = make_graph(3, {{0, 1, 2, 1}, {1, 2, 2, 1}});
  std::vector<Amount> vc = {2, 2};
  auto leak = solve(chain, {0}, {2}, std::nullopt);
  ASSERT_TRUE(certify(chain, vc, {0}, {2}, std::nullopt, nullptr, leak, true).empty());
  const auto out_of_1 = static_cast<std::size_t>(chain.edge_src()[0] == 1 ? 0 : 1);
  leak.edge_flows[out_of_1] = 1;
  EXPECT_TRUE(has_failure(certify(chain, vc, {0}, {2}, std::nullopt, nullptr, leak, true), "conservation at node 1"));
  // The same chain with capacities 5, capped at 1: every vertex is in the closure with
  // potentials (0, 1, 2). Raising pi(1) to 3 makes the reduced costs of 0 -> 1 (1 + 0 - 3) and of
  // the backward arc 2 -> 1 (-1 + 2 - 3) negative; pi(t*) = 2 is still the last curve cost, so
  // these two are the only failures.
  auto wide_chain = make_graph(3, {{0, 1, 5, 1}, {1, 2, 5, 1}});
  std::vector<Amount> vw = {5, 5};
  auto reduced = solve(wide_chain, {0}, {2}, 1);
  ASSERT_TRUE(certify(wide_chain, vw, {0}, {2}, 1, nullptr, reduced, true).empty());
  ASSERT_EQ(reduced.potentials, (std::vector<Cost>{0, 1, 2}));
  reduced.potentials[1] = 3;
  EXPECT_EQ(certify(wide_chain, vw, {0}, {2}, 1, nullptr, reduced, true),
            (std::vector<std::string>{"negative reduced cost inside the closure",
                                      "negative reduced cost inside the closure"}));
}

TEST(FlowCore, CertificateRejectsEveryOtherClause) {
  // The clauses (and each condition of a clause with several) that the test above does not reach,
  // one control each. `mutated` solves on the full capacities, certifies the clean payload, then
  // returns the failed clauses of the mutated copy. A control expects the exact list where the
  // mutation can break one clause alone, else the clause it targets among those it also breaks.
  using Mutation = std::function<void(FlowCoreResult&)>;
  auto mutated = [](const Graph& g, const std::vector<NodeId>& s, const std::vector<NodeId>& d,
                    std::optional<Amount> cap, const Eligibility* el, const Mutation& mutate) {
    std::vector<Amount> view(g.capacity_units().begin(), g.capacity_units().end());
    auto r = place_global_units(g, view, s, d, el, cap);
    EXPECT_TRUE(certify(g, view, s, d, cap, el, r, true).empty());
    mutate(r);
    return certify(g, view, s, d, cap, el, r, true);
  };
  using Fails = std::vector<std::string>;
  // One edge 0 -> 1 of capacity 3 and cost 7. Uncapped, the edge is saturated, so the closure is
  // {0} and the destination's potential (7) is pinned only by the last curve cost; under cap 1 the
  // destination stays in the closure.
  const auto one = make_graph(2, {{0, 1, 3, 7}});
  const auto s = v({0});
  const auto d = v({1});
  // The shape: either array one entry too long.
  EXPECT_EQ(mutated(one, s, d, std::nullopt, nullptr, [](auto& r) { r.edge_flows.push_back(0); }), Fails{"shape"});
  EXPECT_EQ(mutated(one, s, d, std::nullopt, nullptr, [](auto& r) { r.potentials.push_back(0); }), Fails{"shape"});
  // The zero-cap payload, each of its fields alone (the counters are the test above's control).
  const std::vector<Mutation> nonzero = {
      [](auto& r) { r.placed = 1; },
      [](auto& r) { r.cost_distribution = curve({{7, 1}}); },
      [](auto& r) { r.min_cut = std::vector<EdgeId>{}; },
      [](auto& r) { r.edge_flows[0] = 1; },
      [](auto& r) { r.potentials[0] = 0; },
  };
  for (std::size_t k = 0; k < nonzero.size(); ++k) {
    EXPECT_EQ(mutated(one, s, d, 0, nullptr, nonzero[k]), Fails{"zero cap: the prescribed zero payload"}) << k;
  }
  // A negative edge flow (the test above sends one above the view).
  EXPECT_TRUE(has_failure(mutated(one, s, d, std::nullopt, nullptr, [](auto& r) { r.edge_flows[0] = -1; }),
                          "feasibility at edge 0"));
  // Sources {0, 1}, destination 2, edges 0 -> 1 and 1 -> 2: source 1 sends its unit directly.
  // Moving it to 0 -> 1 makes source 1 absorb a unit.
  const auto two_sources = make_graph(3, {{0, 1, 1, 1}, {1, 2, 1, 1}});
  const auto ts = id_of_input(two_sources);
  EXPECT_TRUE(has_failure(mutated(two_sources, v({0, 1}), v({2}), std::nullopt, nullptr,
                                  [&](auto& r) {
                                    r.edge_flows[static_cast<std::size_t>(ts[0])] = 1;
                                    r.edge_flows[static_cast<std::size_t>(ts[1])] = 0;
                                  }),
                          "source absorbs flow"));
  // The same edges with destinations {1, 2} and source 0: the unit stops at destination 1;
  // moving it to 1 -> 2 makes destination 1 emit a unit.
  EXPECT_TRUE(has_failure(mutated(two_sources, v({0}), v({1, 2}), std::nullopt, nullptr,
                                  [&](auto& r) {
                                    r.edge_flows[static_cast<std::size_t>(ts[0])] = 0;
                                    r.edge_flows[static_cast<std::size_t>(ts[1])] = 1;
                                  }),
                          "destination emits flow"));
  // Excluded endpoints: 0 -> 2 and 1 -> 2 with sources {0, 1} and source 0 excluded (and with it
  // the edge 0 -> 2); a unit on 0 -> 2 is supplied by the excluded source. Mirrored, 0 -> 1 and
  // 0 -> 2 with destinations {1, 2} and destination 2 excluded.
  const auto into_2 = make_graph(3, {{0, 2, 1, 1}, {1, 2, 1, 1}});
  const auto x0 = into_2.eligibility(v({0}), {});
  EXPECT_TRUE(has_failure(
      mutated(into_2, v({0, 1}), v({2}), std::nullopt, &x0,
              [&](auto& r) { r.edge_flows[static_cast<std::size_t>(id_of_input(into_2)[0])] = 1; }),
      "excluded source supplies"));
  const auto out_of_0 = make_graph(3, {{0, 1, 1, 1}, {0, 2, 1, 1}});
  const auto x2 = out_of_0.eligibility(v({2}), {});
  EXPECT_TRUE(has_failure(
      mutated(out_of_0, v({0}), v({1, 2}), std::nullopt, &x2,
              [&](auto& r) { r.edge_flows[static_cast<std::size_t>(id_of_input(out_of_0)[1])] = 1; }),
      "excluded destination receives"));
  // The value against the supply alone and the inflow alone: on the chain 0 -> 1 -> 2 carrying
  // 2, one unit less into node 1 lowers the supply only, one unit less out of it the inflow only.
  const auto chain = make_graph(3, {{0, 1, 2, 1}, {1, 2, 2, 1}});
  for (const int edge : {0, 1}) {
    EXPECT_TRUE(has_failure(
        mutated(chain, s, v({2}), std::nullopt, nullptr,
                [&](auto& r) { r.edge_flows[static_cast<std::size_t>(id_of_input(chain)[edge])] = 1; }),
        "value differs from supply or inflow"))
        << edge;
  }
  // placed outside [0, cap]: above the cap, and negative without one.
  EXPECT_TRUE(has_failure(mutated(one, s, d, 1, nullptr, [](auto& r) { r.placed = 2; }), "placed outside [0, cap]"));
  EXPECT_TRUE(
      has_failure(mutated(one, s, d, std::nullopt, nullptr, [](auto& r) { r.placed = -1; }), "placed outside [0, cap]"));
  // The cut: present while cap 1 leaves the destination reachable, missing on the maximal solve,
  // and its value (3) against a lowered placed (which also breaks the value and the curve total).
  EXPECT_EQ(mutated(one, s, d, 1, nullptr, [](auto& r) { r.min_cut = std::vector<EdgeId>{}; }),
            Fails{"min_cut present on a connected solve"});
  EXPECT_EQ(mutated(one, s, d, std::nullopt, nullptr, [](auto& r) { r.min_cut.reset(); }),
            Fails{"min_cut missing on a maximal solve"});
  EXPECT_TRUE(
      has_failure(mutated(one, s, d, std::nullopt, nullptr, [](auto& r) { r.placed = 2; }), "cut value differs from placed"));
  // The potentials: node 1's (outside the closure) below the sentinel (also below pi(t*)), the
  // root's (the closure), pi(t*) against the last curve cost, and no destination potential at all
  // although 3 units arrive (node 1 is reached at zero flow, so the sentinel also breaks that
  // clause).
  EXPECT_TRUE(has_failure(mutated(one, s, d, std::nullopt, nullptr, [](auto& r) { r.potentials[1] = kNeverReached - 1; }),
                          "potential below the sentinel"));
  EXPECT_EQ(mutated(one, s, d, std::nullopt, nullptr, [](auto& r) { r.potentials[0] = 1; }),
            Fails{"allowed root with potential other than 0"});
  EXPECT_EQ(mutated(one, s, d, std::nullopt, nullptr, [](auto& r) { r.potentials[1] = 8; }),
            Fails{"pi(t*) differs from the last curve cost"});
  EXPECT_EQ(mutated(one, s, d, std::nullopt, nullptr, [](auto& r) { r.potentials[1] = kNeverReached; }),
            (Fails{"NEVER_REACHED differs from reachability at zero flow at node 1",
                   "flow delivered but no destination reached"}));
  // The sentinel against reachability at zero flow, each direction alone. On 0 -> 1 and 2 -> 1
  // nothing reaches node 2, which must hold NEVER_REACHED, not 0 or another value. On 0 -> 1
  // (capacity 3, cost 7) and 1 -> 2, uncapped, the first search settles node 2 through node 1 and
  // the saturated edge then leaves both outside the closure: node 2 must not hold NEVER_REACHED.
  const auto unreached = make_graph(3, {{0, 1, 1, 1}, {2, 1, 1, 1}});
  for (const Cost p : {Cost {0}, Cost {123456}}) {
    EXPECT_EQ(mutated(unreached, s, d, std::nullopt, nullptr, [&](auto& r) { r.potentials[2] = p; }),
              Fails{"NEVER_REACHED differs from reachability at zero flow at node 2"})
        << p;
  }
  const auto beyond = make_graph(3, {{0, 1, 3, 7}, {1, 2, 1, 1}});
  EXPECT_EQ(mutated(beyond, s, d, std::nullopt, nullptr, [](auto& r) { r.potentials[2] = kNeverReached; }),
            Fails{"NEVER_REACHED differs from reachability at zero flow at node 2"});
  // The curve total: a free tier of 1 before the real one keeps the moment (21) and pi(t*), and
  // the curve then also differs from the optimum at 1.
  EXPECT_EQ(mutated(one, s, d, std::nullopt, nullptr, [](auto& r) { r.cost_distribution = curve({{0, 1}, {7, 3}}); }),
            (Fails{"curve total differs from placed", "curve differs from the exhaustive optimum at 1"}));
}
