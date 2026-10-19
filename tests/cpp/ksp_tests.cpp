// KSP and resolve_to_paths (docs/semantics.md#k-shortest-paths): the witness graphs (hand-built
// examples: the diamond and the three-hop graph, whose first path is the canonical single path
// of the forward and the reverse search; a shared bottleneck that leaves the fewest hops first; the prefix-capped spur; a
// candidate's bottleneck including its prefix; parallel members under Expand and Collapse;
// accepted paths that never return after a shared prefix; zero-capacity edges (searched like any
// other; a residual filters them; Headroom 0 on them); the ceilings on
// chain-plus-bypasses and the exact CostStretch ceiling at large costs and factors), the defects in the documented order, and agreement with exhaustive
// enumeration (search_reference.hpp) on random graphs; the resolver's order, splitting and
// truncation against an enumeration of every walk of a DAG, from either search direction and on
// long chains.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <functional>
#include <map>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "netgraph/core/detail/route_builders.hpp"
#include "netgraph/core/search.hpp"
#include "search_reference.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

// The route choice of ksp: parallel_edges and tie_break (always single paths).
struct KspChoice {
  bool parallel_edges {true};
  RouteTieBreak tie_break {RouteTieBreak::FewerHops};
};
KspChoice sel(bool parallel_edges = true, RouteTieBreak tie = RouteTieBreak::FewerHops) { return {parallel_edges, tie}; }

// A single-path spf with targets under `tie`.
SpfOptions single_path(Direction dir, RouteTieBreak tie, NodeId target) {
  SpfOptions o(dir);
  o.tie_break = tie;
  o.multipath = false;
  o.targets = std::vector<NodeId> {target};
  return o;
}

// A graph from rows (u, v, capacity, cost), stable-sorted by (source, cost, destination) first: that is the
// EdgeId order, so each EdgeId is the row's position in the sorted list.
Graph bgraph(std::int32_t n, std::vector<InEdge> rows) {
  std::stable_sort(rows.begin(), rows.end(), [](const InEdge& a, const InEdge& b) {
    return std::tie(a.u, a.cost, a.v) < std::tie(b.u, b.cost, b.v);
  });
  return make_graph(n, rows);
}

std::vector<std::pair<Cost, Path>> run(const Graph& g, NodeId s, NodeId t, std::int64_t k, const KspChoice& s1,
                                       std::optional<CostCeiling> ceiling = std::nullopt,
                                       const std::vector<Amount>* residual = nullptr, const Eligibility* elig = nullptr) {
  KspOptions o(k);
  o.tie_break = s1.tie_break;
  o.parallel_edges = s1.parallel_edges;
  o.cost_ceiling = ceiling;
  if (residual != nullptr) o.residual = std::span<const Amount>(*residual);
  o.eligibility = elig;
  return ksp(g, s, t, o);
}

std::vector<EdgeId> walk_of(const Path& p) {
  std::vector<EdgeId> out;
  for (const auto& hop : p.hops()) {
    EXPECT_EQ(hop.size(), 1U);
    out.push_back(hop.front());
  }
  return out;
}

std::vector<NodeId> nodes_of(const Graph& g, const Path& p) {
  std::vector<NodeId> out = {p.source()};
  for (const auto& hop : p.hops()) out.push_back(g.edge_destinations()[static_cast<std::size_t>(hop.front())]);
  return out;
}

std::vector<Cost> costs_of(const std::vector<std::pair<Cost, Path>>& r) {
  std::vector<Cost> out;
  for (const auto& [c, p] : r) {
    EXPECT_EQ(c, p.cost());
    out.push_back(c);
  }
  return out;
}

} // namespace

TEST(Ksp, FirstPathIsTheCanonicalSinglePathOfBothDirections) {
  auto diamond = bgraph(4, {{0, 1, 1, 9}, {0, 2, 1, 1}, {1, 3, 1, 1}, {2, 3, 1, 9}});
  auto three = bgraph(6, {{0, 2, 1, 1}, {2, 3, 1, 1}, {3, 5, 1, 1}, {0, 1, 1, 1}, {1, 4, 1, 1}, {4, 5, 1, 1}});
  struct Case {
    const Graph* g;
    NodeId t;
    std::vector<NodeId> want;
  };
  for (const auto& c : {Case {&diamond, 3, {0, 2, 3}}, Case {&three, 5, {0, 1, 4, 5}}}) {
    for (auto tie : {RouteTieBreak::FewerHops, RouteTieBreak::HigherAvailableCapacity}) {
      auto got = run(*c.g, 0, c.t, 2, sel(true, tie));
      ASSERT_EQ(got.size(), 2U);
      EXPECT_EQ(nodes_of(*c.g, got[0].second), c.want);
      // The first path is the canonical single path of the forward and the reverse search.
      const SpfOptions f = single_path(Direction::Forward, tie, c.t);
      const SpfOptions r = single_path(Direction::Reverse, tie, 0);
      auto fw = spf(*c.g, std::vector<NodeId> {0}, f);
      auto rv = spf(*c.g, std::vector<NodeId> {c.t}, r);
      auto w = walk_of(got[0].second);
      std::sort(w.begin(), w.end());
      auto fe = std::vector<EdgeId>(fw.dag.via_edges().begin(), fw.dag.via_edges().end());
      auto re = std::vector<EdgeId>(rv.dag.via_edges().begin(), rv.dag.via_edges().end());
      std::sort(fe.begin(), fe.end());
      std::sort(re.begin(), re.end());
      EXPECT_EQ(fe, w);
      EXPECT_EQ(re, w);
    }
  }
}

TEST(Ksp, WithASharedBottleneckTheFewestHopsComeFirst) {
  auto g = bgraph(5, {{0, 1, 1, 1}, {1, 4, 1, 3}, {1, 2, 2, 1}, {2, 3, 2, 1}, {3, 4, 2, 1}});
  auto got = run(g, 0, 4, 1, sel(true, RouteTieBreak::HigherCapacity));
  ASSERT_EQ(got.size(), 1U);
  EXPECT_EQ(nodes_of(g, got[0].second), (std::vector<NodeId>{0, 1, 4}));
}

TEST(Ksp, PrefixCappedSpurOrder) {
  // After 0-1-4 the spur at 1 takes 1-5-4 (bottleneck 1, the prefix's) before 1-2-3-4: with the
  // spur's B* capped by the prefix, the 3-hop path comes second (uncapped it would come third).
  auto g = bgraph(6, {{0, 1, 1, 1}, {1, 4, 1, 3}, {1, 5, 1, 1}, {5, 4, 1, 2}, {1, 2, 2, 1}, {2, 3, 2, 1}, {3, 4, 2, 1}});
  auto got = run(g, 0, 4, 3, sel(true, RouteTieBreak::HigherCapacity));
  ASSERT_EQ(got.size(), 3U);
  EXPECT_EQ(nodes_of(g, got[0].second), (std::vector<NodeId>{0, 1, 4}));
  EXPECT_EQ(nodes_of(g, got[1].second), (std::vector<NodeId>{0, 1, 5, 4}));
  EXPECT_EQ(nodes_of(g, got[2].second), (std::vector<NodeId>{0, 1, 2, 3, 4}));
}

TEST(Ksp, CandidateBottleneckIncludesThePrefix) {
  // The first path 0-1-3 (cost 2) has bottleneck 1 (edge 0-1). Its deviations cost 3: at 0,
  // 0-2-3 (bottleneck 5); at 1, the prefix 0-1 then the second 1-3 edge (the spur alone has
  // bottleneck 10, the candidate min(1, 10) = 1). HIGHER_CAPACITY orders 0-2-3 first; a
  // candidate keyed by its spur's bottleneck alone would put the prefix's deviation first.
  // FEWER_HOPS orders the two by EdgeIds: 0-1-3 (ids 0, 3) before 0-2-3 (ids 1, 4).
  auto g = bgraph(4, {{0, 1, 1, 1}, {1, 3, 10, 1}, {0, 2, 5, 1}, {2, 3, 5, 2}, {1, 3, 10, 2}});
  auto prefer = run(g, 0, 3, 8, sel(true, RouteTieBreak::HigherCapacity));
  ASSERT_EQ(costs_of(prefer), (std::vector<Cost> {2, 3, 3}));
  EXPECT_EQ(nodes_of(g, prefer[0].second), (std::vector<NodeId> {0, 1, 3}));
  EXPECT_EQ(nodes_of(g, prefer[1].second), (std::vector<NodeId> {0, 2, 3}));
  EXPECT_EQ(nodes_of(g, prefer[2].second), (std::vector<NodeId> {0, 1, 3}));
  EXPECT_EQ(walk_of(prefer[2].second), (std::vector<EdgeId> {0, 3}));
  auto det = run(g, 0, 3, 8, sel(true, RouteTieBreak::FewerHops));
  ASSERT_EQ(costs_of(det), (std::vector<Cost> {2, 3, 3}));
  EXPECT_EQ(walk_of(det[1].second), (std::vector<EdgeId> {0, 3}));
  EXPECT_EQ(walk_of(det[2].second), (std::vector<EdgeId> {1, 4}));
}

TEST(Ksp, ParallelMembersUnderExpandAndCollapse) {
  // Two hops of two parallel members: four equal paths under Expand, one under Collapse.
  auto pr = bgraph(3, {{0, 1, 1, 1}, {0, 1, 1, 1}, {1, 2, 1, 1}, {1, 2, 1, 1}});
  EXPECT_EQ(run(pr, 0, 2, 8, sel()).size(), 4U);
  EXPECT_EQ(run(pr, 0, 2, 8, sel(false)).size(), 1U);
  auto k3 = run(pr, 0, 2, 3, sel());
  ASSERT_EQ(k3.size(), 3U);
  EXPECT_EQ(walk_of(k3[0].second), (std::vector<EdgeId>{0, 2}));
  EXPECT_EQ(walk_of(k3[1].second), (std::vector<EdgeId>{0, 3}));
  EXPECT_EQ(walk_of(k3[2].second), (std::vector<EdgeId>{1, 2}));
  // Costs 1/3 then 2/3: all four simple paths under Expand in cost order; Collapse keeps one.
  auto fe = bgraph(3, {{0, 1, 1, 1}, {0, 1, 1, 3}, {1, 2, 1, 2}, {1, 2, 1, 3}});
  EXPECT_EQ(costs_of(run(fe, 0, 2, 8, sel())), (std::vector<Cost>{3, 4, 5, 6}));
  EXPECT_EQ(costs_of(run(fe, 0, 2, 8, sel(false))), (std::vector<Cost>{3}));
}

TEST(Ksp, AcceptedPathsNeverReturnAfterASharedPrefix) {
  // A candidate goes into the pending set without a search of the accepted
  // paths, because each accepted path with the spur's prefix has its next edge (Expand) or its
  // next transition (Collapse) banned at the spur. Spurs here share their prefix with one or two
  // accepted paths at once, over parallel members; every path comes once, in the derived order.
  // EdgeIds: 0, 1: 0->1 (1); 2: 0->2 (3); 3: 1->2 (1); 4, 5: 1->3 (1); 6: 2->3 (1).
  auto g = bgraph(4, {{0, 1, 1, 1}, {0, 1, 1, 1}, {0, 2, 1, 3}, {1, 2, 1, 1}, {1, 3, 1, 1}, {1, 3, 1, 1}, {2, 3, 1, 1}});
  auto walks = [](const std::vector<std::pair<Cost, Path>>& r) {
    std::vector<std::vector<EdgeId>> out;
    for (const auto& [c, p] : r) out.push_back(walk_of(p));
    return out;
  };
  // Expand, derived by hand from the deviations: [0, 4] first; its spurs give [1, 4] and [0, 5].
  // After [0, 5] the spur at node 1 bans 4 and 5 (both accepted paths share the prefix [0]) and
  // finds [0, 3, 6]; the spur at node 0 bans 0 (twice) and finds [1, 4] again (pending: kept
  // once). [1, 4] adds [2, 6] and [1, 5]; [1, 5] adds [1, 3, 6]; the rest add nothing new.
  const auto expand = run(g, 0, 3, 16, sel());
  EXPECT_EQ(walks(expand), (std::vector<std::vector<EdgeId>> {{0, 4}, {0, 5}, {1, 4}, {1, 5}, {0, 3, 6}, {1, 3, 6}, {2, 6}}));
  EXPECT_EQ(costs_of(expand), (std::vector<Cost> {2, 2, 2, 2, 3, 3, 4}));
  // Collapse: [0, 4] first; at node 1 the transition to 3 is banned, so 5 never follows [0]; at
  // node 0 the transition to 1 is banned. Then [0, 3, 6] (node 1 bans the heads 3 and 2) and
  // [2, 6].
  const auto collapse = run(g, 0, 3, 16, sel(false));
  EXPECT_EQ(walks(collapse), (std::vector<std::vector<EdgeId>> {{0, 4}, {0, 3, 6}, {2, 6}}));
  EXPECT_EQ(costs_of(collapse), (std::vector<Cost> {2, 3, 4}));
  // Every k gives the prefix of the full answer.
  const auto all_e = walks(expand);
  const auto all_c = walks(collapse);
  for (std::int64_t k = 1; k <= 8; ++k) {
    EXPECT_EQ(walks(run(g, 0, 3, k, sel())), decltype(all_e)(all_e.begin(), all_e.begin() + std::min<std::int64_t>(k, 7)));
    EXPECT_EQ(walks(run(g, 0, 3, k, sel(false))), decltype(all_c)(all_c.begin(), all_c.begin() + std::min<std::int64_t>(k, 3)));
  }
}

TEST(Ksp, EndpointDefectsReportTheFirstInOrder) {
  // The scalar endpoints are checked as one-element sets (range, then the
  // overlap), with every other argument also wrong where the order allows it; the first defect
  // in the validation order is reported, for ksp and resolve_to_paths.
  auto g = bgraph(3, {{0, 1, 1, 1}, {1, 2, 1, 1}});
  auto other = bgraph(3, {{0, 1, 1, 1}});
  auto foreign = other.eligibility({}, std::vector<EdgeId> {0});
  const std::vector<Amount> bad = {1, 2};
  const auto all_wrong = [&](NodeId s, NodeId t) {
    return error_of([&] { (void)run(g, s, t, 0, sel(), AbsoluteCost {-1}, &bad, &foreign); });
  };
  EXPECT_EQ(all_wrong(5, -1), "source: node 5 out of node range [0, 3)");
  EXPECT_EQ(all_wrong(-1, -1), "source: node -1 out of node range [0, 3)");
  EXPECT_EQ(all_wrong(3, 3), "source: node 3 out of node range [0, 3)");  // range before the overlap
  EXPECT_EQ(all_wrong(0, 3), "destination: node 3 out of node range [0, 3)");
  EXPECT_EQ(all_wrong(2, -4), "destination: node -4 out of node range [0, 3)");
  EXPECT_EQ(all_wrong(1, 1), "node 1 is in both endpoint sets: source and destination must not overlap");
  EXPECT_EQ(all_wrong(0, 2), "k must be >= 1 (got 0)");
  KspOptions bad_tie(0);
  bad_tie.tie_break = static_cast<RouteTieBreak>(9);
  EXPECT_EQ(error_of([&] { (void)ksp(g, 7, 7, bad_tie); }), "tie_break: 9 is not a value of RouteTieBreak");
  const auto dag = detail::PredecessorDAGBuilder::of_edges(g, std::vector<EdgeId> {0, 1});
  const auto dag_other = detail::PredecessorDAGBuilder::of_edges(other, std::vector<EdgeId> {0});
  EXPECT_EQ(error_of([&] { (void)resolve_to_paths(g, dag_other, 9, 9, false, 0); }), "foreign PredecessorDAG: built for another graph");
  EXPECT_EQ(error_of([&] { (void)resolve_to_paths(g, dag, 9, -9, false, 0); }), "source: node 9 out of node range [0, 3)");
  EXPECT_EQ(error_of([&] { (void)resolve_to_paths(g, dag, 3, 3, false, 0); }), "source: node 3 out of node range [0, 3)");
  EXPECT_EQ(error_of([&] { (void)resolve_to_paths(g, dag, 0, 3, false, 0); }), "destination: node 3 out of node range [0, 3)");
  EXPECT_EQ(error_of([&] { (void)resolve_to_paths(g, dag, 2, 2, false, 0); }),
            "node 2 is in both endpoint sets: source and destination must not overlap");
  EXPECT_EQ(error_of([&] { (void)resolve_to_paths(g, dag, 0, 2, false, 0); }), "max_paths must be >= 1 (got 0)");
  // Valid endpoints at both ends of the range.
  EXPECT_EQ(run(g, 0, 2, 1, sel()).size(), 1U);
  EXPECT_EQ(resolve_to_paths(g, dag, 0, 2).size(), 1U);
  EXPECT_TRUE(run(g, 2, 0, 1, sel()).empty());
}

TEST(Ksp, ZeroCapacityEdgesAreSearchedAndOnlyAResidualFiltersThem) {
  // Ids: 0: 0->1, 1: 0->2 (cost 1), 2: 0->3 cap 0, 3: 0->3 cap 1 (cost 2), 4: 1->3, 5: 2->3 (cost 1).
  auto g = bgraph(4, {{0, 1, 1, 1}, {1, 3, 1, 1}, {0, 2, 1, 1}, {2, 3, 1, 1}, {0, 3, 0, 2}, {0, 3, 1, 2}});
  ASSERT_EQ(g.edge_capacity_units()[2], 0);
  auto walks = [](const std::vector<std::pair<Cost, Path>>& r) {
    std::vector<std::vector<EdgeId>> out;
    for (const auto& [c, p] : r) out.push_back(walk_of(p));
    return out;
  };
  // Four walks of cost 2, the zero-capacity edge searched like any other: the one-hop walks
  // first ([2] before [3] by EdgeId), then the two-hop walks.
  auto all = run(g, 0, 3, 8, sel());
  EXPECT_EQ(costs_of(all), (std::vector<Cost>{2, 2, 2, 2}));
  EXPECT_EQ(walks(all), (std::vector<std::vector<EdgeId>>{{2}, {3}, {0, 4}, {1, 5}}));
  // A residual equal to the capacities filters the zero-capacity edge (it has no unit free), and
  // so does an Eligibility excluding it.
  const std::vector<EdgeId> zero = {2};
  auto el = g.eligibility({}, zero);
  const std::vector<std::vector<EdgeId>> without = {{3}, {0, 4}, {1, 5}};
  EXPECT_EQ(walks(run(g, 0, 3, 8, sel(), std::nullopt, nullptr, &el)), without);
  const std::vector<Amount> full(g.edge_capacity_units().begin(), g.edge_capacity_units().end());
  EXPECT_EQ(walks(run(g, 0, 3, 8, sel(), std::nullopt, &full)), without);
}

TEST(Ksp, ZeroEdgeExamplesOfTheContract) {
  // KSP over the zero edge, with and without a residual (ids 0: 0->1 cap 0, 1: 0->2 cap 5, 2: 2->1 cap 5,
  // costs 1): k = 2 lists [0] (cost 1) and [1, 2] (cost 2); the residual (0, 5, 5) leaves [1, 2].
  auto g = bgraph(3, {{0, 1, 0, 1}, {0, 2, 5, 1}, {2, 1, 5, 1}});
  auto r = run(g, 0, 1, 2, sel());
  ASSERT_EQ(r.size(), 2U);
  EXPECT_EQ(costs_of(r), (std::vector<Cost>{1, 2}));
  EXPECT_EQ(walk_of(r[0].second), (std::vector<EdgeId>{0}));
  EXPECT_EQ(walk_of(r[1].second), (std::vector<EdgeId>{1, 2}));
  const std::vector<Amount> res = {0, 5, 5};
  auto gated = run(g, 0, 1, 2, sel(), std::nullopt, &res);
  ASSERT_EQ(gated.size(), 1U);
  EXPECT_EQ(costs_of(gated), (std::vector<Cost>{2}));
  EXPECT_EQ(walk_of(gated[0].second), (std::vector<EdgeId>{1, 2}));
  // The ceiling CostStretch(6/5): ids 0: 0->1 cap 0 cost 5, 1: 0->1 cap 1 cost 6. The static
  // baseline counts the zero-capacity edge: 5; 6 * 5 <= 6 * 5 exactly, so both walks fit (the
  // binary value of the float 1.2 would keep only the first); at 11/10 the cost-6 walk does not.
  auto h = bgraph(2, {{0, 1, 0, 5}, {0, 1, 1, 6}});
  EXPECT_EQ(costs_of(run(h, 0, 1, 2, sel(), CostStretch(6, 5))), (std::vector<Cost>{5, 6}));
  auto ceiled = run(h, 0, 1, 2, sel(), CostStretch(11, 10));
  ASSERT_EQ(ceiled.size(), 1U);
  EXPECT_EQ(costs_of(ceiled), (std::vector<Cost>{5}));
  EXPECT_EQ(walk_of(ceiled[0].second), (std::vector<EdgeId>{0}));
}

TEST(Ksp, HeadroomWithoutAResidualReadsOneOrZero) {
  // Two equal-cost walks: [0, 2] over 0->1 (cap 0) and 1->3, [1, 3] over 0->2 and 2->3 (cap 5).
  // HEADROOM without a residual reads capacity / capacity: 0 on the zero-capacity edge, else 1.
  // HIGHER_HEADROOM ranks the bottleneck-1 walk first; FEWER_HOPS the EdgeId order.
  auto g = bgraph(4, {{0, 1, 0, 1}, {0, 2, 5, 1}, {1, 3, 5, 1}, {2, 3, 5, 1}});
  const auto head = sel(true, RouteTieBreak::HigherHeadroom);
  auto r = run(g, 0, 3, 2, head);
  ASSERT_EQ(r.size(), 2U);
  EXPECT_EQ(walk_of(r[0].second), (std::vector<EdgeId>{1, 3}));
  EXPECT_EQ(walk_of(r[1].second), (std::vector<EdgeId>{0, 2}));
  auto d = run(g, 0, 3, 2, sel());
  ASSERT_EQ(d.size(), 2U);
  EXPECT_EQ(walk_of(d[0].second), (std::vector<EdgeId>{0, 2}));
  // On positive capacities 1 and 5 every edge reads 1: the bottlenecks tie and the EdgeId order
  // decides (CAPACITY would rank the capacity-5 walk first).
  auto q = bgraph(4, {{0, 1, 1, 1}, {0, 2, 5, 1}, {1, 3, 1, 1}, {2, 3, 5, 1}});
  EXPECT_EQ(walk_of(run(q, 0, 3, 1, head)[0].second), (std::vector<EdgeId>{0, 2}));
  const auto cap = sel(true, RouteTieBreak::HigherCapacity);
  EXPECT_EQ(walk_of(run(q, 0, 3, 1, cap)[0].second), (std::vector<EdgeId>{1, 3}));
}

TEST(Ksp, CeilingsOnChainPlusBypasses) {
  // Chain 0-1-2-3 (cost 3) with bypasses 0->2 and 1->3 (cost 3 each): static baseline 3.
  auto g = bgraph(4, {{0, 1, 1, 1}, {1, 2, 1, 1}, {2, 3, 1, 1}, {0, 2, 1, 3}, {1, 3, 1, 3}});
  EXPECT_EQ(costs_of(run(g, 0, 3, 8, sel())), (std::vector<Cost>{3, 4, 4}));
  EXPECT_EQ(costs_of(run(g, 0, 3, 8, sel(), CostStretch(1))), (std::vector<Cost>{3}));
  EXPECT_EQ(costs_of(run(g, 0, 3, 8, sel(), CostStretch(3, 2))), (std::vector<Cost>{3, 4, 4}));
  EXPECT_EQ(costs_of(run(g, 0, 3, 8, sel(), CostStretch(5, 4))), (std::vector<Cost>{3}));  // 4 * 4 > 5 * 3
  EXPECT_EQ(costs_of(run(g, 0, 3, 8, sel(), CostStretch(4, 3))), (std::vector<Cost>{3, 4, 4}));  // 4 * 3 = 4 * 3
  EXPECT_EQ(costs_of(run(g, 0, 3, 8, sel(), AbsoluteCost {3})), (std::vector<Cost>{3}));
  EXPECT_EQ(costs_of(run(g, 0, 3, 8, sel(), AbsoluteCost {4})), (std::vector<Cost>{3, 4, 4}));
  EXPECT_TRUE(run(g, 0, 3, 8, sel(), AbsoluteCost {0}).empty());  // a valid empty budget
  // The CostStretch baseline is static: a residual gating the middle edge removes the cost-3 path but
  // not the baseline, so CostStretch(1) keeps nothing (the residual paths cost 4).
  std::vector<Amount> res(5, 1);
  const auto inv = id_of_input(g);
  res[static_cast<std::size_t>(inv[2])] = 0;  // the middle edge 1 -> 2
  EXPECT_TRUE(run(g, 0, 3, 8, sel(), CostStretch(1), &res).empty());
  EXPECT_EQ(costs_of(run(g, 0, 3, 8, sel(), std::nullopt, &res)), (std::vector<Cost>{4, 4}));
  // No static route under a CostStretch ceiling: [].
  const std::vector<NodeId> xn = {1, 2};
  auto el = g.eligibility(xn, {});
  EXPECT_TRUE(run(g, 0, 3, 8, sel(), CostStretch(2), nullptr, &el).empty());
}

// The CostStretch ceiling is exact where a 64-bit product overflows: cost * d <= n * base in 128
// bits. Paths 0->1 of cost b and 0->2->1 of cost b + 1 (b + 2 on the second graph). b = 2^60: b + 1
// rounds to b in a double, so CostStretch(1) admitting both paths would be a double comparison;
// (2^60 + 1) / 2^60 admits b + 1 at equality and (2^61 + 1) / 2^61 does not ((2^60 + 1) * 2^61 =
// 2^121 + 2^61 > 2^121 + 2^60); the largest numerator admits both. b = 2^52: (2^52 + 1) / 2^52
// admits b + 1, not b + 2.
TEST(Ksp, StretchCeilingIsExactForLargeCostsAndFactors) {
  const Cost big = Cost {1} << 60;
  auto g = bgraph(3, {{0, 1, 1, big}, {0, 2, 1, 1}, {2, 1, 1, big}});
  EXPECT_EQ(costs_of(run(g, 0, 1, 8, sel(), CostStretch(1))), (std::vector<Cost> {big}));
  EXPECT_EQ(costs_of(run(g, 0, 1, 8, sel(), CostStretch(big + 1, big))), (std::vector<Cost> {big, big + 1}));
  EXPECT_EQ(costs_of(run(g, 0, 1, 8, sel(), CostStretch(2 * big + 1, 2 * big))), (std::vector<Cost> {big}));
  EXPECT_EQ(costs_of(run(g, 0, 1, 8, sel(), CostStretch(std::numeric_limits<std::int64_t>::max()))),
            (std::vector<Cost> {big, big + 1}));
  const Cost b52 = Cost {1} << 52;
  auto h = bgraph(4, {{0, 1, 1, b52}, {0, 2, 1, 1}, {2, 1, 1, b52}, {0, 3, 1, 2}, {3, 1, 1, b52}});
  EXPECT_EQ(costs_of(run(h, 0, 1, 8, sel(), CostStretch(1))), (std::vector<Cost> {b52}));
  EXPECT_EQ(costs_of(run(h, 0, 1, 8, sel(), CostStretch(b52 + 1, b52))), (std::vector<Cost> {b52, b52 + 1}));
  EXPECT_EQ(costs_of(run(h, 0, 1, 8, sel(), CostStretch(2))), (std::vector<Cost> {b52, b52 + 1, b52 + 2}));
  // Components at 2^63 - 1 on a cost of 2^62 - 1 (c = c0): the products stay below 2^125.
  const Cost top = (Cost {1} << 62) - 1;
  constexpr std::int64_t kMax = std::numeric_limits<std::int64_t>::max();
  EXPECT_EQ(costs_of(run(bgraph(2, {{0, 1, 1, top}}), 0, 1, 1, sel(), CostStretch(kMax, kMax - 1))),
            (std::vector<Cost> {top}));
}

TEST(Ksp, ValidationInTheContractOrder) {
  auto g = bgraph(3, {{0, 1, 1, 1}, {1, 2, 1, 1}});
  auto other = bgraph(3, {{0, 1, 1, 1}});
  auto foreign = other.eligibility({}, std::vector<EdgeId> {0});
  const std::vector<Amount> bad = {1, 2};
  // The endpoints and k, alone and with every later argument wrong: the test above.
  EXPECT_EQ(error_of([&] { (void)run(g, 0, 2, 1, sel(), AbsoluteCost {-1}); }), "cost_ceiling AbsoluteCost cost must be >= 0 (got -1)");
  EXPECT_EQ(error_of([&] { (void)run(g, 0, 2, 1, sel(), CostStretch(1, 2)); }),
            "cost_ceiling CostStretch factor must be >= 1 (got 1/2)");
  EXPECT_EQ(error_of([&] { (void)run(g, 0, 2, 1, sel(), CostStretch(std::numeric_limits<std::int64_t>::min())); }),
            "cost_ceiling CostStretch factor must be >= 1 (got -9223372036854775808/1)");
  EXPECT_EQ(error_of([&] { (void)run(g, 0, 2, 1, sel(), CostStretch(1, 2), &bad); }),
            "cost_ceiling CostStretch factor must be >= 1 (got 1/2)");  // before the residual
  EXPECT_EQ(error_of([&] { (void)run(g, 0, 2, 1, sel(), std::nullopt, &bad); }), "residual[1] = 2 is outside [0, edge_capacity_units[1] = 1]");
  EXPECT_EQ(error_of([&] { (void)run(g, 0, 2, 1, sel(), std::nullopt, nullptr, &foreign); }),
            "foreign Eligibility: built for another graph");
  // Order: the ceiling before the residual (the endpoints before k before the ceiling: the test
  // above).
  EXPECT_EQ(error_of([&] { (void)run(g, 0, 2, 1, sel(), AbsoluteCost {-1}, &bad); }), "cost_ceiling AbsoluteCost cost must be >= 0 (got -1)");
  KspOptions o(1);
  o.tie_break = static_cast<RouteTieBreak>(4);
  EXPECT_EQ(error_of([&] { (void)ksp(g, 0, 2, o); }), "tie_break: 4 is not a value of RouteTieBreak");
  // An excluded endpoint, the source or the destination, gives [], never an error.
  for (const NodeId x : {0, 2}) {
    const std::vector<NodeId> xn = {x};
    auto el = g.eligibility(xn, {});
    EXPECT_TRUE(run(g, 0, 2, 3, sel(), std::nullopt, nullptr, &el).empty()) << x;
    EXPECT_TRUE(run(g, 0, 2, 3, sel(), CostStretch(2), nullptr, &el).empty()) << x;
  }
}

// A residual array with several defects at once (one error, from the shared edge-amount check): its
// length before any entry, then its first entry outside [0, capacity] in EdgeId order, before
// the Eligibility.
TEST(Ksp, ResidualDefectsReportTheFirstInOrder) {
  const auto g = edge_amount_graph();
  const auto other = make_graph(4, {{0, 1, 2, 1}});
  const auto foreign = other.eligibility({}, std::vector<EdgeId> {0});
  for (const auto& c : edge_amount_defects()) {
    EXPECT_EQ(error_of([&] { (void)run(g, 0, 3, 2, sel(), std::nullopt, &c.amounts, &foreign); }), "residual" + c.error);
  }
  const auto bounds = edge_amounts_at_the_bounds();
  EXPECT_EQ(error_of([&] { (void)run(g, 0, 3, 2, sel(), std::nullopt, &bounds); }), "");
}

namespace {

// The reference KSP list: every simple source -> destination path of the member graph sorted by the key,
// truncated at the inclusive ceiling and at k.
std::vector<std::vector<EdgeId>> reference_ksp(const Graph& g, const ref::Query& q, NodeId s, NodeId t, std::int64_t k,
                                               std::optional<CostCeiling> ceiling) {
  const auto mem = ref::members(q);
  std::vector<char> is_end(static_cast<std::size_t>(g.num_nodes()), 0);
  is_end[static_cast<std::size_t>(t)] = 1;
  std::vector<std::vector<EdgeId>> all;
  std::vector<EdgeId> walk;
  std::vector<char> on(static_cast<std::size_t>(g.num_nodes()), 0);
  on[static_cast<std::size_t>(s)] = 1;
  ref::walks(g, mem, s, is_end, walk, on, all);
  std::sort(all.begin(), all.end(), [&](const auto& a, const auto& b) { return ref::key_less(q, a, b); });
  // The static baseline: the member graph without the residual, Expand.
  ref::Query st = q;
  st.residual = nullptr;
  st.parallel_edges = true;
  const auto d = ref::distances(g, ref::members(st), {s}, false);
  std::vector<std::vector<EdgeId>> out;
  for (const auto& p : all) {
    if (static_cast<std::int64_t>(out.size()) == k) break;
    Cost c = 0;
    for (EdgeId e : p) c += g.edge_costs()[static_cast<std::size_t>(e)];
    if (ceiling.has_value()) {
      if (const auto* a = std::get_if<AbsoluteCost>(&*ceiling)) {
        if (c > a->cost) break;
      } else {
        const Cost base = d[static_cast<std::size_t>(t)];
        if (base == ref::kNone) return {};
        // c <= (n / d) * base exactly: c * d <= n * base.
        const auto& f = std::get<CostStretch>(*ceiling);
        if (static_cast<Int128>(c) * f.denominator() > static_cast<Int128>(f.numerator()) * base) break;
      }
    }
    out.push_back(p);
  }
  return out;
}

} // namespace

TEST(Ksp, RandomGraphsAgreeWithExhaustiveEnumeration) {
  // Random graphs (seed 6006) checked against exhaustive enumeration: parallel members, residuals,
  // masks, every tie_break, Expand and Collapse, k in {1, 3, 8}, AbsoluteCost and CostStretch
  // ceilings.
  std::mt19937 rng(6006);
  int nonempty = 0;
  int deep = 0;
  for (int gi = 0; gi < 300; ++gi) {
    const int n = 3 + static_cast<int>(rng() % 5);
    const int m = n + static_cast<int>(rng() % static_cast<unsigned>(2 * n + 2));
    std::vector<InEdge> rows;
    for (int i = 0; i < m; ++i) {
      rows.push_back({static_cast<NodeId>(rng() % static_cast<unsigned>(n)), static_cast<NodeId>(rng() % static_cast<unsigned>(n)),
                      static_cast<Amount>(rng() % 5), static_cast<Cost>(1 + rng() % 5)});
    }
    if (rng() % 2 == 0) {
      const auto u = static_cast<NodeId>(rng() % static_cast<unsigned>(n));
      const auto v = static_cast<NodeId>(rng() % static_cast<unsigned>(n));
      for (int i = 0; i < 2; ++i) rows.push_back({u, v, static_cast<Amount>(1 + rng() % 4), static_cast<Cost>(1 + rng() % 3)});
    }
    rows.push_back({0, static_cast<NodeId>(n - 1), static_cast<Amount>(1 + rng() % 3), static_cast<Cost>(3 + rng() % 6)});
    auto g = make_graph(n, rows);
    std::vector<Amount> residual(static_cast<std::size_t>(g.num_edges()));
    for (std::size_t e = 0; e < residual.size(); ++e) {
      residual[e] = static_cast<Amount>(rng() % static_cast<unsigned>(g.edge_capacity_units()[e] + 1));
    }
    std::vector<NodeId> xn;
    std::vector<EdgeId> xe;
    for (NodeId v = 1; v + 1 < n; ++v) {
      if (rng() % 10 == 0) xn.push_back(v);
    }
    for (EdgeId e = 0; e < g.num_edges(); ++e) {
      if (rng() % 10 == 0) xe.push_back(e);
    }
    const auto el = g.eligibility(xn, xe);
    for (int c = 0; c < 24; ++c) {
      const bool use_res = rng() % 2 == 0;
      const bool masks = rng() % 10 < 3;
      const RouteTieBreak higher[3] = {RouteTieBreak::HigherAvailableCapacity, RouteTieBreak::HigherCapacity,
                                       RouteTieBreak::HigherHeadroom};
      const bool parallel = rng() % 2 == 0;
      const bool prefer = rng() % 2 != 0;
      const RouteTieBreak tie = higher[rng() % 3];
      const auto s1 = sel(parallel, prefer ? tie : RouteTieBreak::FewerHops);
      const std::int64_t ks[3] = {1, 3, 8};
      const std::int64_t k = ks[rng() % 3];
      std::optional<CostCeiling> ceiling;
      const auto pick = rng() % 4;
      if (pick == 2) ceiling = AbsoluteCost {static_cast<Cost>(3 + rng() % 9)};
      if (pick == 3) {
        const std::pair<std::int64_t, std::int64_t> fs[4] = {{1, 1}, {5, 4}, {3, 2}, {2, 1}};
        const auto [n, d] = fs[rng() % 4];
        ceiling = CostStretch(n, d);
      }
      ref::Query q;
      q.g = &g;
      q.elig = masks ? &el : nullptr;
      q.residual = use_res ? &residual : nullptr;
      q.parallel_edges = s1.parallel_edges;
      q.tie_break = s1.tie_break;
      const auto want = reference_ksp(g, q, 0, static_cast<NodeId>(n - 1), k, ceiling);
      auto got = run(g, 0, static_cast<NodeId>(n - 1), k, s1, ceiling, use_res ? &residual : nullptr, masks ? &el : nullptr);
      std::vector<std::vector<EdgeId>> walks;
      for (const auto& [cost, p] : got) {
        walks.push_back(walk_of(p));
        EXPECT_EQ(cost, p.cost());
      }
      ASSERT_EQ(walks, want) << "graph " << gi << " case " << c;
      if (!want.empty()) ++nonempty;
      if (want.size() >= 3) ++deep;
    }
  }
  EXPECT_GT(nonempty, 3000);
  EXPECT_GT(deep, 300);
}

namespace {

// Every source -> destination walk of a DAG given by its edges, as hops (grouped or split), sorted by the
// resolver's key (cost, hops, edge-id sequence with a grouped hop keyed by its lowest member).
std::vector<std::vector<std::vector<EdgeId>>> reference_walks(const Graph& g, const std::vector<EdgeId>& dag, NodeId s,
                                                              NodeId t, bool split) {
  std::map<std::pair<NodeId, NodeId>, std::vector<EdgeId>> groups;
  for (EdgeId e : dag) groups[{g.edge_sources()[static_cast<std::size_t>(e)], g.edge_destinations()[static_cast<std::size_t>(e)]}].push_back(e);
  for (auto& [k, v] : groups) std::sort(v.begin(), v.end());
  std::vector<std::vector<std::vector<EdgeId>>> out;
  std::function<void(NodeId, std::vector<std::vector<EdgeId>>&)> walk = [&](NodeId u, std::vector<std::vector<EdgeId>>& hops) {
    if (u == t) {
      if (!split) {
        out.push_back(hops);
        return;
      }
      std::vector<std::vector<std::vector<EdgeId>>> combos(1);
      for (const auto& h : hops) {
        std::vector<std::vector<std::vector<EdgeId>>> next;
        for (const auto& c : combos) {
          for (EdgeId e : h) {
            auto cc = c;
            cc.push_back({e});
            next.push_back(cc);
          }
        }
        combos = next;
      }
      out.insert(out.end(), combos.begin(), combos.end());
      return;
    }
    for (const auto& [k, v] : groups) {
      if (k.first != u) continue;
      hops.push_back(v);
      walk(k.second, hops);
      hops.pop_back();
    }
  };
  std::vector<std::vector<EdgeId>> hops;
  walk(s, hops);
  auto key = [&](const std::vector<std::vector<EdgeId>>& p) {
    Cost c = 0;
    std::vector<EdgeId> seq;
    for (const auto& h : p) {
      Cost best = g.edge_costs()[static_cast<std::size_t>(h.front())];
      for (EdgeId e : h) best = std::min(best, g.edge_costs()[static_cast<std::size_t>(e)]);
      c += best;
      seq.push_back(h.front());
    }
    return std::make_tuple(c, p.size(), seq);
  };
  std::stable_sort(out.begin(), out.end(), [&](const auto& a, const auto& b) { return key(a) < key(b); });
  return out;
}

} // namespace

TEST(ResolveToPaths, OrderSplitAndTruncation) {
  // A DAG with a parallel pair, two routes of equal cost and a dearer one:
  // EdgeIds: 0: 0->1 (1), 1: 0->1 (1), 2: 0->2 (1), 3: 0->3 (5), 4: 1->3 (1), 5: 2->3 (1).
  auto g = bgraph(4, {{0, 1, 1, 1}, {0, 1, 1, 1}, {0, 2, 1, 1}, {0, 3, 1, 5}, {1, 3, 1, 1}, {2, 3, 1, 1}});
  const auto dag = detail::PredecessorDAGBuilder::of_edges(g, std::vector<EdgeId> {0, 1, 2, 3, 4, 5});
  auto hops_of = [](const std::vector<Path>& ps) {
    std::vector<std::vector<std::vector<EdgeId>>> out;
    for (const auto& p : ps) out.push_back(p.hops());
    return out;
  };
  // Grouped: [[0, 1], [4]] (cost 2, key 0), [[2], [5]] (cost 2, key 2), [[3]] (cost 5).
  auto grouped = resolve_to_paths(g, dag, 0, 3);
  EXPECT_EQ(hops_of(grouped), (std::vector<std::vector<std::vector<EdgeId>>> {{{0, 1}, {4}}, {{2}, {5}}, {{3}}}));
  EXPECT_EQ(grouped[0].cost(), 2);
  EXPECT_EQ(grouped[2].cost(), 5);
  // Split: one path per member combination, sorted by its own (cost, hops, edge-id sequence).
  auto split = resolve_to_paths(g, dag, 0, 3, true);
  EXPECT_EQ(hops_of(split), (std::vector<std::vector<std::vector<EdgeId>>> {{{0}, {4}}, {{1}, {4}}, {{2}, {5}}, {{3}}}));
  EXPECT_EQ(resolve_to_paths(g, dag, 0, 3, true, 2).size(), 2U);
  EXPECT_EQ(hops_of(resolve_to_paths(g, dag, 0, 3, false, 1)), (std::vector<std::vector<std::vector<EdgeId>>> {{{0, 1}, {4}}}));
  // An endpoint the DAG does not join: [].
  EXPECT_TRUE(resolve_to_paths(g, dag, 3, 0).empty());
  EXPECT_TRUE(resolve_to_paths(g, detail::PredecessorDAGBuilder::of_edges(g, std::vector<EdgeId> {4}), 0, 3).empty());
  // Unequal parallel costs: expanding the first grouped walk in place would put cost 4 before
  // the direct cost-3 path, and truncating there would discard the wrong path.
  const auto mixed = bgraph(3, {{0, 1, 1, 1}, {0, 1, 1, 3}, {1, 2, 1, 1}, {0, 2, 1, 3}});
  const auto all = detail::PredecessorDAGBuilder::of_edges(mixed, all_ids(mixed));
  const auto grouped_mixed = resolve_to_paths(mixed, all, 0, 2);
  EXPECT_EQ(hops_of(grouped_mixed), (std::vector<std::vector<std::vector<EdgeId>>> {{{0, 1}, {3}}, {{2}}}));
  const auto split_mixed = resolve_to_paths(mixed, all, 0, 2, true);
  EXPECT_EQ(hops_of(split_mixed), (std::vector<std::vector<std::vector<EdgeId>>> {{{0}, {3}}, {{2}}, {{1}, {3}}}));
  ASSERT_EQ(split_mixed.size(), 3U);
  EXPECT_EQ((std::vector<Cost> {split_mixed[0].cost(), split_mixed[1].cost(), split_mixed[2].cost()}),
            (std::vector<Cost> {2, 3, 4}));
  EXPECT_EQ(hops_of(resolve_to_paths(mixed, all, 0, 2, true, 2)),
            (std::vector<std::vector<std::vector<EdgeId>>> {{{0}, {3}}, {{2}}}));
  // Validation: Ksp.EndpointDefectsReportTheFirstInOrder.
}

TEST(ResolveToPaths, RandomDagsAgreeWithEnumeration) {
  // Random DAGs (edges u -> v with u < v, parallel members, mixed costs, not necessarily tight):
  // the order, both forms and the truncation against the enumeration of every walk.
  std::mt19937 rng(55);
  int nonempty = 0;   // checks with at least one walk
  int truncated = 0;  // checks whose max_paths cut the walks short
  for (int trial = 0; trial < 600; ++trial) {
    const int n = 2 + static_cast<int>(rng() % 7);
    std::vector<InEdge> rows;
    const int m = static_cast<int>(rng() % 18);
    for (int i = 0; i < m; ++i) {
      const auto u = static_cast<NodeId>(rng() % static_cast<unsigned>(n - 1));
      const auto v = static_cast<NodeId>(u + 1 + rng() % static_cast<unsigned>(n - 1 - u));
      rows.push_back({u, v, 1, static_cast<Cost>(1 + rng() % 4)});
    }
    auto g = make_graph(n, rows);
    std::vector<EdgeId> edges;
    for (EdgeId e = 0; e < g.num_edges(); ++e) {
      if (rng() % 5 != 0) edges.push_back(e);
    }
    const auto dag = detail::PredecessorDAGBuilder::of_edges(g, edges);
    const NodeId s = 0;
    const auto t = static_cast<NodeId>(n - 1);
    for (bool split : {false, true}) {
      const auto want = reference_walks(g, edges, s, t, split);
      for (std::optional<std::int64_t> mp : {std::optional<std::int64_t> {}, std::optional<std::int64_t> {1},
                                            std::optional<std::int64_t> {3}}) {
        const auto got = resolve_to_paths(g, dag, s, t, split, mp);
        const std::size_t expect_n = mp.has_value() ? std::min<std::size_t>(want.size(), static_cast<std::size_t>(*mp)) : want.size();
        ASSERT_EQ(got.size(), expect_n) << trial;
        for (std::size_t i = 0; i < got.size(); ++i) EXPECT_EQ(got[i].hops(), want[i]) << trial << " path " << i;
        nonempty += want.empty() ? 0 : 1;
        truncated += expect_n < want.size() ? 1 : 0;
      }
    }
  }
  RecordProperty("nonempty", nonempty);
  RecordProperty("truncated", truncated);
  EXPECT_GT(nonempty, 1800);  // not vacuous: most DAGs hold a walk
  EXPECT_GT(truncated, 400);
}

TEST(ResolveToPaths, ResolvesASearchDagFromEitherDirection) {
  // The DAG of a reverse search is forward-oriented: resolving it from a source gives the same
  // walks as resolving the forward search's DAG.
  auto g = bgraph(5, {{0, 1, 1, 1}, {0, 2, 1, 1}, {1, 3, 1, 1}, {2, 3, 1, 1}, {3, 4, 1, 1}, {0, 4, 1, 3}});
  SpfOptions f(Direction::Forward);
  SpfOptions r(Direction::Reverse);
  auto fw = spf(g, std::vector<NodeId> {0}, f);
  auto rv = spf(g, std::vector<NodeId> {4}, r);
  auto a = resolve_to_paths(g, fw.dag, 0, 4);
  auto b = resolve_to_paths(g, rv.dag, 0, 4);
  ASSERT_EQ(a.size(), 3U);
  ASSERT_EQ(b.size(), 3U);
  for (std::size_t i = 0; i < a.size(); ++i) EXPECT_TRUE(a[i] == b[i]);
  EXPECT_EQ(a[0].cost(), 3);
}

TEST(ResolveToPaths, LongChainsNeedNoRecursion) {
  // A 200,000-node chain with one parallel pair: the enumeration keeps an explicit stack.
  const int n = 200000;
  std::vector<InEdge> rows;
  for (int v = 0; v + 1 < n; ++v) rows.push_back({v, v + 1, 1, 1});
  rows.push_back({n / 2, n / 2 + 1, 1, 1});
  auto g = make_graph(n, rows);
  const auto dag = detail::PredecessorDAGBuilder::of_edges(g, all_ids(g));
  EXPECT_EQ(resolve_to_paths(g, dag, 0, n - 1, true).size(), 2U);
  EXPECT_EQ(resolve_to_paths(g, dag, 0, n - 1).size(), 1U);
}
