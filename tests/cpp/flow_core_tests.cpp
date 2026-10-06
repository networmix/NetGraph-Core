// The flow core of APPROACH 2.2 through place_global_units: a first solve on a fresh graph, the
// review witnesses (expected
// values derived by hand and from prototype 1, `prototypes/p1_amount/flowcore.py`), the
// certificate in its two branches, input validation, determinism, and randomized graphs
// certified against an exhaustive optimum (tier a) and an independent solver (tier b).
#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "netgraph/core/detail/flow_core.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

constexpr Amount kL62m1 = (Amount{1} << 62) - 1;

// The counters fixed by the normative per-tier order and the augmenting-path sequence, compared
// with prototype 1: dijkstra_runs, dijkstra_pops, settled, dinic_phases, augmentations. arc_scans
// is the implementation's own diagnostic (the level search here first marks the vertices that
// reach t* and stops at t*'s level, so it examines a different number of arcs than the
// prototype's full BFS); it is checked for determinism and presence only.
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

std::string error_of(const std::function<void()>& fn) {
  try {
    fn();
  } catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
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
}

TEST(FlowCore, QuantisedCapacitiesFamily) {
  // Physical 2.5q is 3 units and 0.5q is 1 unit (review-2 R1 at 1a: the quantised network).
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
  // review-14 flow.chain, review-1 C7: value 2 needs the middle edge cancelled.
  auto g = make_graph(4, chain_plus_bypasses());
  auto one = solve(g, {0}, {3}, 1);
  EXPECT_EQ(one.placed, 1);
  EXPECT_EQ(one.cost_distribution, curve({{3, 1}}));
  EXPECT_EQ(to_input(g, one.edge_flows), (std::vector<Amount>{1, 1, 1, 0, 0}));
  EXPECT_FALSE(one.min_cut.has_value());  // review-14 flow.finite_cap: value 1 of a maximum 2
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
  // review-14 flow.equal_cost_multiple_hop_tiers, review-16 solver.multiple_dinic_phases:
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
  // review-16 solver.cyclic_admissible: two parallel cost-1 links. Potentials are updated by
  // the first Dijkstra before the drain (on the zero potentials no arc has reduced cost 0), and
  // the zero-reduced-cost residual cycle 0 -> 1 -> 0 left after the first push is drained in
  // the same tier.
  auto g = make_graph(2, {{0, 1, 1, 1}, {0, 1, 1, 1}});
  auto r = solve(g, {0}, {1}, 2);
  EXPECT_EQ(r.placed, 2);
  EXPECT_EQ(r.cost_distribution, curve({{1, 2}}));
  EXPECT_EQ(counters_of(r), (Counters{1, 4, 4, 1, 2}));
  expect_certified(g, {0}, {1}, 2, r);
  // review-17 partial_drain: cap 1 stops inside the tier; the destinations stay connected.
  auto p = solve(g, {0}, {1}, 1);
  EXPECT_EQ(p.placed, 1);
  EXPECT_EQ(p.cost_distribution, curve({{1, 1}}));
  EXPECT_FALSE(p.min_cut.has_value());
  EXPECT_EQ(counters_of(p), (Counters{1, 4, 4, 1, 1}));
  expect_certified(g, {0}, {1}, 1, p);
}

TEST(FlowCore, FailedFinalRunRaisesPotentials) {
  // review-17 failed_final_run: the side branch 3 is reached at 8 by the first run and at 10 by
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
  // review-14 flow.release_no_growth, review-15 fresh.owner_release_atomic: on the view W + own
  // = [1, 1, 1] the solve at value 1 takes the cheap link (cost 2), not the old cost-4 one.
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
  // review-15 C1 and review-16 certificate.internal_cap on one edge of capacity 3 and cost 7.
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
  auto inf = solve(g, {0}, {1}, std::nullopt);
  EXPECT_EQ(inf.placed, 3);
  EXPECT_EQ(inf.cost_distribution, curve({{7, 3}}));
  EXPECT_EQ(counters_of(inf), (Counters{2, 6, 6, 2, 1}));
  expect_certified(g, {0}, {1}, std::nullopt, inf);
  // The internal cap 2^63 - 2: the solve runs to disconnection (no clamp, APPROACH 2.2 r25).
  auto big = solve(g, {0}, {1}, kMaxFlowCap);
  EXPECT_TRUE(same_payload(big, inf));
  expect_certified(g, {0}, {1}, kMaxFlowCap, big);
}

TEST(FlowCore, ZeroCapIsThePrescribedZeroPayload) {
  // review-16 certificate.zero_cap, review-17 C2.zero_cap_branch.
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
  // review-17 C2.positive_cap_no_out_capacity: a positive cap runs one failed Dijkstra even when
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
}

TEST(FlowCore, LegalMaximumUnits) {
  // review-17 terminal_saturation: one edge at 2^62 - 1 units; the closure is seeded at the real
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
  // review-17 C2.reverse_virtual and review-16 global.all_destinations: both destinations have
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
  // review-14 flow.cut: the closure follows reverse arcs, so the forward-only boundary is wrong.
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

TEST(FlowCore, Review13MasksEqualTheirCuts) {
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
  EXPECT_NE(call({2}, s, d, 1).find("view has length 1, expected 2"), std::string::npos);
  EXPECT_NE(call({2, 3}, s, d, 1).find("capacity_units"), std::string::npos);
  EXPECT_NE(call({-1, 2}, s, d, 1).find("outside"), std::string::npos);
  EXPECT_NE(call(view, {}, d, 1).find("sources: the node set is empty"), std::string::npos);
  EXPECT_NE(call(view, s, {}, 1).find("destinations: the node set is empty"), std::string::npos);
  EXPECT_NE(call(view, {3}, d, 1).find("out of node range"), std::string::npos);
  EXPECT_NE(call(view, s, {-1}, 1).find("out of node range"), std::string::npos);
  EXPECT_NE(call(view, {0, 2}, {2}, 1).find("both endpoint sets"), std::string::npos);
  EXPECT_NE(call(view, s, d, -1).find("cap"), std::string::npos);
  EXPECT_NE(call(view, s, d, kMaxFlowCap + 1).find("cap"), std::string::npos);
  auto foreign = h.eligibility();
  EXPECT_NE(call(view, s, d, 1, &foreign).find("foreign Eligibility"), std::string::npos);
  // The overlap is checked on the full sets, also when an Eligibility excludes the vertex.
  const std::vector<NodeId> x0 = {0};
  auto el = g.eligibility(x0, {});
  EXPECT_NE(call(view, {0}, {0}, 1, &el).find("both endpoint sets"), std::string::npos);
  // Order: the view before the endpoints, the endpoints before the cap.
  EXPECT_NE(call({2}, {}, d, -1).find("view"), std::string::npos);
  EXPECT_NE(call(view, {}, d, -1).find("sources"), std::string::npos);
}

TEST(FlowCore, IdenticalCallsIdenticalPayloads) {
  std::mt19937 rng(7);
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
  }
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

} // namespace

TEST(FlowCore, RandomGraphsTierAExhaustive) {
  // Nonuniform costs, parallel edges, self-loops, partial views, endpoint sets, masks and caps;
  // every payload certified, the curve against the minimum cost of every integer value.
  std::mt19937 rng(20261004);
  int positive = 0;
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
  }
  EXPECT_GT(positive, 250);  // the population is not vacuous
}

TEST(FlowCore, RandomGraphsTierBIndependent) {
  std::mt19937 rng(1816);
  int multi_tier = 0;
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
  }
  EXPECT_GT(multi_tier, 30);
}

TEST(FlowCore, RandomLargeUnitsTierB) {
  // Capacities near the domain bound (sum below 2^62) and costs near 2^40.
  std::mt19937_64 rng(62);
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
  }
}

TEST(FlowCore, CertificateRejectsMutatedPayloads) {
  // The test oracle itself is not vacuous: the review-17 negative controls are rejected.
  auto one = make_graph(2, {{0, 1, 3, 7}});
  std::vector<Amount> v1 = {3};
  auto under = solve(one, {0}, {1}, 1);
  EXPECT_FALSE(certify(one, v1, {0}, {1}, 2, nullptr, under, true).empty());  // underfill under cap 2
  auto slope = solve(one, {0}, {1}, 2);
  slope.cost_distribution = curve({{6, 2}});
  EXPECT_FALSE(certify(one, v1, {0}, {1}, 2, nullptr, slope, true).empty());
  auto cut = solve(one, {0}, {1}, std::nullopt);
  cut.min_cut = std::vector<EdgeId>{};
  EXPECT_FALSE(certify(one, v1, {0}, {1}, std::nullopt, nullptr, cut, true).empty());
  auto par = make_graph(2, {{0, 1, 1, 1}, {0, 1, 1, 3}});
  std::vector<Amount> vp = {1, 1};
  auto omitted = solve(par, {0}, {1}, std::nullopt);
  omitted.cost_distribution = curve({{2, 2}});  // same moment, omitted breakpoint at 1
  EXPECT_FALSE(certify(par, vp, {0}, {1}, std::nullopt, nullptr, omitted, true).empty());
  EXPECT_FALSE(certify(par, vp, {0}, {1}, std::nullopt, nullptr, omitted, false).empty());
  auto any = make_graph(3, {{0, 1, 1, 1}, {0, 2, 1, 2}});
  std::vector<Amount> va = {1, 1};
  auto pot = solve(any, {0}, {1, 2}, std::nullopt);
  pot.potentials[1] = 1;  // an inflow destination below pi(t*)
  EXPECT_FALSE(certify(any, va, {0}, {1, 2}, std::nullopt, nullptr, pot, true).empty());
  auto sentinel = solve(one, {0}, {1}, 1);
  sentinel.potentials = {kNeverReached, kNeverReached};
  EXPECT_FALSE(certify(one, v1, {0}, {1}, 1, nullptr, sentinel, true).empty());
  auto zero = solve(one, {0}, {1}, 0);
  zero.counters.dijkstra_runs = 1;
  EXPECT_FALSE(certify(one, v1, {0}, {1}, 0, nullptr, zero, true).empty());
}
