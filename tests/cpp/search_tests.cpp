// The search through the public spf and the traversal core's internal entries: the
// member graph and replaces() (the parallel-member choice), the unique DAG in both directions and the values that reach
// it only through the member graph, anycast pruning, the staged
// canonical single path on its witness graphs (hand-built examples), the canonical predecessor order, both queues and
// their dispatch, and agreement with an independent reference (search_reference.hpp) on random
// graphs, one population of them built so that every rule of the single path's key decides
// some trials.
#include <gtest/gtest.h>

#include <algorithm>
#include <functional>
#include <numeric>
#include <optional>
#include <queue>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "netgraph/core/shortest_path_cache.hpp"
#include "netgraph/core/detail/route_builders.hpp"
#include "netgraph/core/detail/search_detail.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/search.hpp"
#include "search_reference.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

// A route choice of spf: multipath, parallel_edges, tie_break.
detail::RouteChoice sel(bool multipath = true, bool parallel_edges = true, RouteTieBreak tie = RouteTieBreak::FewerHops) {
  return detail::RouteChoice {tie, multipath, parallel_edges};
}

constexpr RouteTieBreak kTieBreaks[] = {RouteTieBreak::FewerHops, RouteTieBreak::HigherCapacity,
                                        RouteTieBreak::HigherAvailableCapacity, RouteTieBreak::HigherHeadroom};

SpfResult run(const Graph& g, std::vector<NodeId> roots, Direction dir, std::optional<std::vector<NodeId>> targets,
              const detail::RouteChoice& s, const std::vector<Amount>* residual = nullptr, Amount min_residual = 1,
              const Eligibility* elig = nullptr, SpfQueue queue = SpfQueue::Auto) {
  SpfOptions o(dir);
  o.tie_break = s.tie_break;
  o.multipath = s.multipath;
  o.parallel_edges = s.parallel_edges;
  o.targets = std::move(targets);
  if (residual != nullptr) o.residual = std::span<const Amount>(*residual);
  o.min_residual = min_residual;
  o.eligibility = elig;
  return spf(g, roots, o, queue);
}

std::vector<EdgeId> edges_of(const PredecessorDAG& dag) {
  std::vector<EdgeId> v(dag.via_edges().begin(), dag.via_edges().end());
  std::sort(v.begin(), v.end());
  return v;
}

std::vector<EdgeId> sorted(std::vector<EdgeId> v) {
  std::sort(v.begin(), v.end());
  return v;
}

// The DAG's canonical form: entries grouped by child, strictly ascending via within a child,
// every via joining (parent, child); the internal ascending list equals the sorted via.
void expect_canonical(const Graph& g, const PredecessorDAG& dag) {
  const auto off = dag.parent_offsets();
  ASSERT_EQ(off.size(), static_cast<std::size_t>(g.num_nodes()) + 1);
  for (std::size_t v = 0; v + 1 < off.size(); ++v) {
    for (auto k = off[v]; k < off[v + 1]; ++k) {
      const auto e = static_cast<std::size_t>(dag.via_edges()[static_cast<std::size_t>(k)]);
      EXPECT_EQ(g.edge_destinations()[e], static_cast<NodeId>(v));
      EXPECT_EQ(g.edge_sources()[e], dag.parents()[static_cast<std::size_t>(k)]);
      if (k > off[v]) EXPECT_LT(dag.via_edges()[static_cast<std::size_t>(k) - 1], dag.via_edges()[static_cast<std::size_t>(k)]);
    }
  }
  const auto asc = detail::PredecessorDAGBuilder::ascending(dag);
  EXPECT_EQ(std::vector<EdgeId>(asc.begin(), asc.end()), edges_of(dag));
}

} // namespace

TEST(Search, MemberGraphIsEligibilityAndResidualNeverCapacity) {
  // ids: 0: 0->1 cap 0, 1: 0->1 cap 3, 2: 0->1 cap 5 (cost 2), 3: 1->2 cap 4. A full search keeps
  // the tight members, the zero-capacity id 0 included (capacity is never a condition): the
  // cost-2 member is never tight while a cost-1 member exists.
  auto g = make_graph(3, {{0, 1, 0, 1}, {0, 1, 3, 1}, {0, 1, 5, 2}, {1, 2, 4, 1}});
  EXPECT_EQ(edges_of(run(g, {0}, Direction::Forward, std::nullopt, sel()).dag), (std::vector<EdgeId>{0, 1, 3}));
  // A residual gates at min_residual: residual [0, 1, 5, 4] at 2 units removes id 1, so the
  // cost-2 member carries node 1 (distance 2).
  const std::vector<Amount> res = {0, 1, 5, 4};
  auto gated = run(g, {0}, Direction::Forward, std::nullopt, sel(), &res, 2);
  EXPECT_EQ(*gated.costs, (std::vector<Cost>{0, 2, 3}));
  EXPECT_EQ(edges_of(gated.dag), (std::vector<EdgeId>{2, 3}));
  // min_residual without a residual gates nothing.
  auto static_only = run(g, {0}, Direction::Forward, std::nullopt, sel(), nullptr, 4);
  EXPECT_EQ(*static_only.costs, (std::vector<Cost>{0, 1, 2}));
  const std::vector<EdgeId> xe = {3};
  auto el = g.eligibility({}, xe);
  auto masked = run(g, {0}, Direction::Forward, std::nullopt, sel(), nullptr, 1, &el);
  EXPECT_EQ(*masked.costs, (std::vector<Cost>{0, 1, kUnreachable}));
}

TEST(Search, CollapseBySelectAdjacency) {
  // Parallel 0->1: id 0 cost 1 cap 2, id 1 cost 1 cap 6, id 2 cost 1 cap 4, id 3 cost 2 cap 9.
  auto g = make_graph(2, {{0, 1, 2, 1}, {0, 1, 6, 1}, {0, 1, 4, 1}, {0, 1, 9, 2}});
  auto member = [&](const detail::RouteChoice& s, const std::vector<Amount>* res = nullptr) {
    return edges_of(run(g, {0}, Direction::Forward, std::nullopt, s, res).dag);
  };
  // Expand keeps every tight member; FEWER_HOPS Collapse the lowest EdgeId among the cheapest.
  EXPECT_EQ(member(sel()), (std::vector<EdgeId>{0, 1, 2}));
  // multipath and parallel_edges default to true: a full search keeping every tight member.
  EXPECT_EQ(edges_of(spf(g, std::vector<NodeId> {0}, SpfOptions(Direction::Forward)).dag), (std::vector<EdgeId>{0, 1, 2}));
  EXPECT_EQ(member(sel(true, false)), (std::vector<EdgeId>{0}));
  // HIGHER_CAPACITY: the highest capacity among the cheapest.
  EXPECT_EQ(member(sel(true, false, RouteTieBreak::HigherCapacity)), (std::vector<EdgeId>{1}));
  // HIGHER_AVAILABLE_CAPACITY reads the residual; an equal residual keeps the lowest id.
  const auto resid = sel(true, false, RouteTieBreak::HigherAvailableCapacity);
  const std::vector<Amount> r1 = {2, 1, 2, 9};
  const std::vector<Amount> r2 = {1, 1, 3, 9};
  EXPECT_EQ(member(resid, &r1), (std::vector<EdgeId>{0}));
  EXPECT_EQ(member(resid, &r2), (std::vector<EdgeId>{2}));
  // HIGHER_HEADROOM compares residual / capacity exactly: 1/2 < 4/6 = 2/3 > 2/4.
  const auto head = sel(true, false, RouteTieBreak::HigherHeadroom);
  const std::vector<Amount> h1 = {1, 4, 2, 9};
  const std::vector<Amount> h2 = {1, 1, 2, 9};
  EXPECT_EQ(member(head, &h1), (std::vector<EdgeId>{1}));
  // Equal rational headroom (2/4 == 1/2): the lowest EdgeId.
  EXPECT_EQ(member(head, &h2), (std::vector<EdgeId>{0}));
  // Without a residual HIGHER_AVAILABLE_CAPACITY reads the static capacities.
  EXPECT_EQ(member(resid), (std::vector<EdgeId>{1}));
}

TEST(Search, HeadroomRationalsBeyondDoubles) {
  // Headroom ratios a double cannot separate: (a - 1) / a < a / (a + 1) for a = 2^54, while both round to the
  // same double. Parallel 0->1 members: id 0 residual a - 1 of capacity a, id 1 residual a of
  // capacity a + 1; Headroom Collapse must keep id 1.
  const Amount a = Amount {1} << 54;
  auto g = make_graph(2, {{0, 1, a, 1}, {0, 1, a + 1, 1}});
  const std::vector<Amount> res = {a - 1, a};
  auto r = run(g, {0}, Direction::Forward, std::nullopt,
               sel(true, false, RouteTieBreak::HigherHeadroom), &res);
  EXPECT_EQ(edges_of(r.dag), (std::vector<EdgeId>{1}));
}

TEST(Search, FullDagForwardAndReverse) {
  // Diamond with a parallel member and a longer branch: 0->1, 0->1, 0->2 (cost 2), 1->3, 2->3.
  auto g = make_graph(4, {{0, 1, 1, 1}, {0, 1, 1, 1}, {0, 2, 1, 2}, {1, 3, 1, 1}, {2, 3, 1, 1}});
  auto f = run(g, {0}, Direction::Forward, std::nullopt, sel());
  ASSERT_TRUE(f.costs.has_value());
  EXPECT_FALSE(f.target_cost.has_value());
  EXPECT_EQ(*f.costs, (std::vector<Cost>{0, 1, 2, 2}));
  // Tight edges: both parallel members, 0->2 (d(2) = 2) and 1->3; 2->3 is not tight (2 + 1 > 2).
  EXPECT_EQ(edges_of(f.dag), (std::vector<EdgeId>{0, 1, 2, 3}));
  expect_canonical(g, f.dag);
  auto r = run(g, {3}, Direction::Reverse, std::nullopt, sel());
  EXPECT_EQ(*r.costs, (std::vector<Cost>{2, 1, 1, 0}));
  // The reverse DAG toward 3 is the forward-oriented set of tight edges toward 3. It differs from
  // the forward DAG from 0 although 3 is a node that search reaches: each search keeps the edges
  // tight in its own direction (0->2 forward only, 2->3 reverse only).
  EXPECT_EQ(edges_of(r.dag), (std::vector<EdgeId>{0, 1, 3, 4}));
  expect_canonical(g, r.dag);
  // On the same root the directions differ too: no edge enters 0, so the reverse DAG is empty.
  EXPECT_TRUE(run(g, {0}, Direction::Reverse, std::nullopt, sel()).dag.via_edges().empty());
  // Unreached nodes keep kUnreachable.
  auto u = run(g, {2}, Direction::Forward, std::nullopt, sel());
  EXPECT_EQ(*u.costs, (std::vector<Cost>{kUnreachable, kUnreachable, 0, 1}));
}

TEST(Search, PredecessorOrderIsAscendingViaIndependentlyOfDiscovery) {
  // A four-node witness: ids 0: 0->1 (1), 1: 0->2 (2), 2: 1->3 (2), 3: 2->3 (1).
  // Node 3's tight parents are discovered via 2 (from node 1, settled first) and then via 3; the
  // list holds them newest first, [3, 2]; the DAG must list [2, 3] whatever the queue.
  auto g = make_graph(4, {{0, 1, 1, 1}, {0, 2, 1, 2}, {1, 3, 1, 2}, {2, 3, 1, 1}});
  for (auto q : {SpfQueue::Auto, SpfQueue::Heap}) {
    for (bool targets : {false, true}) {
      auto r = run(g, {0}, Direction::Forward, targets ? std::optional<std::vector<NodeId>>(std::vector<NodeId> {3})
                                                       : std::nullopt,
                   sel(), nullptr, 1, nullptr, q);
      const auto off = r.dag.parent_offsets();
      ASSERT_EQ(off[4] - off[3], 2);
      EXPECT_EQ(r.dag.via_edges()[static_cast<std::size_t>(off[3])], 2);
      EXPECT_EQ(r.dag.via_edges()[static_cast<std::size_t>(off[3]) + 1], 3);
      expect_canonical(g, r.dag);
    }
  }
  // The same through the reverse search, whose lists hang on the tails.
  auto rev = run(g, {3}, Direction::Reverse, std::nullopt, sel());
  expect_canonical(g, rev.dag);
  EXPECT_EQ(edges_of(rev.dag), (std::vector<EdgeId>{0, 1, 2, 3}));
}

// The canonical-form check reports each defect of a hand-built DAG: a descending via within a
// child, a wrong parent, an entry under a child its edge does not enter, an ascending list that
// is not the sorted via; the correct DAG passes.
TEST(Search, CanonicalFormCheckReportsEachDefect) {
  // ids: 0: 0->1, 1: 0->1, 2: 0->2 (cost 2), 3: 1->2.
  auto g = make_graph(3, {{0, 1, 1, 1}, {0, 1, 1, 1}, {0, 2, 1, 2}, {1, 2, 1, 1}});
  const auto failures = [&](std::vector<NodeId> nodes, std::vector<std::int32_t> offsets, std::vector<NodeId> parents,
                            std::vector<EdgeId> via, std::vector<EdgeId> ascending) {
    const auto dag = detail::PredecessorDAGBuilder::of_parts(g, std::move(nodes), std::move(offsets), std::move(parents),
                                                      std::move(via), std::move(ascending));
    return failures_of([&] { expect_canonical(g, dag); });
  };
  EXPECT_EQ(failures({1, 2}, {0, 2, 3}, {0, 0, 1}, {0, 1, 3}, {0, 1, 3}), 0);
  EXPECT_GT(failures({1, 2}, {0, 2, 3}, {0, 0, 1}, {1, 0, 3}, {0, 1, 3}), 0) << "descending via";
  EXPECT_GT(failures({1, 2}, {0, 2, 3}, {0, 0, 0}, {0, 1, 3}, {0, 1, 3}), 0) << "wrong parent";
  EXPECT_GT(failures({1}, {0, 3}, {0, 0, 1}, {0, 1, 3}, {0, 1, 3}), 0) << "wrong child";
  EXPECT_GT(failures({1, 2}, {0, 2, 3}, {0, 0, 1}, {0, 1, 3}, {0, 1, 2}), 0) << "ascending list";
}

TEST(Search, AnycastPruningAndDirectionIndependence) {
  // 0->1 (1), 1->2 (1), 0->3 (2), 3->2 (5), 1->4 (5): targets {2, 4}: nearest is 2 at cost 2.
  auto g = make_graph(5, {{0, 1, 1, 1}, {1, 2, 1, 1}, {0, 3, 1, 2}, {3, 2, 1, 5}, {1, 4, 1, 5}});
  const auto inv = id_of_input(g);
  auto f = run(g, {0}, Direction::Forward, std::vector<NodeId> {2, 4}, sel());
  ASSERT_TRUE(f.target_cost.has_value());
  EXPECT_EQ(*f.target_cost, 2);
  EXPECT_FALSE(f.costs.has_value());
  EXPECT_EQ(edges_of(f.dag), sorted({inv[0], inv[1]}));
  // A reverse search rooted at the targets with the source as its target: the same edges.
  auto r = run(g, {2, 4}, Direction::Reverse, std::vector<NodeId> {0}, sel());
  EXPECT_EQ(edges_of(r.dag), edges_of(f.dag));
  EXPECT_EQ(r.target_cost, f.target_cost);
  // No reachable target: kUnreachable, no edges.
  auto none = run(g, {2}, Direction::Forward, std::vector<NodeId> {0}, sel());
  EXPECT_EQ(none.target_cost, kUnreachable);
  EXPECT_TRUE(none.dag.via_edges().empty());
}

TEST(Search, AnycastKeepsEveryEquallyNearTarget) {
  // Two targets at the same distance through different branches, a third farther: the search
  // stops at the first target popped but keeps the walks to both nearest ones.
  // ids by input: 0: 0->1, 1: 0->2, 2: 1->3, 3: 2->4, 4: 0->5 (cost 5).
  auto g = make_graph(6, {{0, 1, 1, 1}, {0, 2, 1, 1}, {1, 3, 1, 1}, {2, 4, 1, 1}, {0, 5, 1, 5}});
  const auto inv = id_of_input(g);
  for (auto q : {SpfQueue::Auto, SpfQueue::Heap}) {
    auto r = run(g, {0}, Direction::Forward, std::vector<NodeId> {3, 4, 5}, sel(), nullptr, 1, nullptr, q);
    EXPECT_EQ(r.target_cost, 2);
    EXPECT_EQ(edges_of(r.dag), sorted({inv[0], inv[1], inv[2], inv[3]}));
  }
}

TEST(Search, ExcludedRootsAndTargetsAreSkipped) {
  auto g = make_graph(3, {{0, 1, 1, 1}, {1, 2, 1, 1}});
  const std::vector<NodeId> xn = {0};
  auto el = g.eligibility(xn, {});
  auto r = run(g, {0}, Direction::Forward, std::vector<NodeId> {2}, sel(), nullptr, 1, &el);
  EXPECT_EQ(r.target_cost, kUnreachable);
  EXPECT_TRUE(r.dag.via_edges().empty());
  auto full = run(g, {0}, Direction::Forward, std::nullopt, sel(), nullptr, 1, &el);
  EXPECT_EQ(*full.costs, (std::vector<Cost>{kUnreachable, kUnreachable, kUnreachable}));
  // An excluded root in a set is skipped; the others search.
  auto two = run(g, {0, 1}, Direction::Forward, std::nullopt, sel(), nullptr, 1, &el);
  EXPECT_EQ(*two.costs, (std::vector<Cost>{kUnreachable, 0, 1}));
  // An excluded target is unreachable; another target still counts.
  const std::vector<NodeId> x2 = {2};
  auto el2 = g.eligibility(x2, {});
  auto t = run(g, {0}, Direction::Forward, std::vector<NodeId> {1, 2}, sel(), nullptr, 1, &el2);
  EXPECT_EQ(t.target_cost, 1);
}

TEST(Search, ValidationBeforeAnyWorkInTheContractOrder) {
  auto g = make_graph(3, {{0, 1, 2, 1}, {1, 2, 2, 1}});
  auto other = make_graph(3, {{0, 1, 2, 1}});
  auto foreign = other.eligibility({}, std::vector<EdgeId> {0});
  const std::vector<Amount> short_res = {1};
  const std::vector<Amount> big_res = {1, 3};
  const std::vector<Amount> ok_res = {1, 1};
  EXPECT_EQ(error_of([&] { (void)run(g, {}, Direction::Forward, std::nullopt, sel()); }), "roots: the node set is empty");
  EXPECT_EQ(error_of([&] { (void)run(g, {3}, Direction::Forward, std::nullopt, sel()); }),
            "roots: node 3 out of node range [0, 3)");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::vector<NodeId> {}, sel()); }),
            "targets: the node set is empty");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::vector<NodeId> {-1}, sel()); }),
            "targets: node -1 out of node range [0, 3)");
  EXPECT_EQ(error_of([&] { (void)run(g, {0, 1}, Direction::Forward, std::vector<NodeId> {1}, sel()); }),
            "node 1 is in both endpoint sets: roots and targets must not overlap");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(false)); }),
            "multipath = false requires targets (the single path joins a root to a target)");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(), &short_res); }),
            "residual has length 1, expected 2 (one entry per edge)");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(), &big_res); }),
            "residual[1] = 3 is outside [0, edge_capacity_units[1] = 2]");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(), nullptr, 0); }),
            "min_residual must be at least 1 unit (got 0)");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(), nullptr, 1, &foreign); }),
            "foreign Eligibility: built for another graph");
  // The order when several defects apply: roots before targets before the overlap before the
  // single path before the ineffective tie before the residual before min_residual before the
  // Eligibility.
  EXPECT_EQ(error_of([&] { (void)run(g, {}, Direction::Forward, std::vector<NodeId> {}, sel()); }),
            "roots: the node set is empty");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::vector<NodeId> {0, 7}, sel()); }),
            "targets: node 7 out of node range [0, 3)");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::vector<NodeId> {0}, sel(), &short_res); }),
            "node 0 is in both endpoint sets: roots and targets must not overlap");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(false), &short_res, 0); }),
            "multipath = false requires targets (the single path joins a root to a target)");
  EXPECT_EQ(error_of([&] {
              (void)run(g, {0}, Direction::Forward, std::nullopt, sel(true, true, RouteTieBreak::HigherCapacity),
                        &short_res, 0);
            }),
            "spf tie_break requires multipath=False or parallel_edges=False");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(), &big_res, 0, &foreign); }),
            "residual[1] = 3 is outside [0, edge_capacity_units[1] = 2]");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(), &ok_res, 0, &foreign); }),
            "min_residual must be at least 1 unit (got 0)");
  // Every root excluded with valid sets is an empty result, never an error; an excluded root
  // does not excuse a defect.
  const std::vector<NodeId> all = {0, 1, 2};
  auto none = g.eligibility(all, {});
  EXPECT_EQ(run(g, {0}, Direction::Forward, std::vector<NodeId> {2}, sel(), nullptr, 1, &none).target_cost, kUnreachable);
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::vector<NodeId> {0}, sel(), nullptr, 1, &none); }),
            "node 0 is in both endpoint sets: roots and targets must not overlap");
  // Enumeration values outside their enumeration, at construction and at the entry.
  EXPECT_EQ(error_of([&] { SpfOptions o(static_cast<Direction>(9)); }), "direction: 9 is not a value of Direction");
  SpfOptions bad(Direction::Forward);
  bad.direction = static_cast<Direction>(7);
  EXPECT_EQ(error_of([&] { (void)spf(g, std::vector<NodeId> {0}, bad); }), "direction: 7 is not a value of Direction");
  SpfOptions bad_tie(Direction::Forward);
  bad_tie.tie_break = static_cast<RouteTieBreak>(5);
  EXPECT_EQ(error_of([&] { (void)spf(g, std::vector<NodeId> {0}, bad_tie); }),
            "tie_break: 5 is not a value of RouteTieBreak");
  SpfOptions good(Direction::Forward);
  EXPECT_EQ(error_of([&] { (void)spf(g, std::vector<NodeId> {0}, good, static_cast<SpfQueue>(3)); }),
            "queue: 3 is not a value of SpfQueue");
  // The enumerations come before every other check, in the order direction, tie_break, the queue.
  // Each call has every later defect too (empty roots, multipath = false without targets, a short
  // residual, min_residual 0, a foreign Eligibility); repairing the reported defect moves the
  // error to the next one.
  SpfOptions every(Direction::Forward);
  every.multipath = false;
  every.residual = std::span<const Amount>(short_res);
  every.min_residual = 0;
  every.eligibility = &foreign;
  every.direction = static_cast<Direction>(7);
  every.tie_break = static_cast<RouteTieBreak>(6);
  const auto bad_queue = static_cast<SpfQueue>(3);
  const auto staircase_error = [&] { return error_of([&] { (void)spf(g, std::vector<NodeId> {}, every, bad_queue); }); };
  EXPECT_EQ(staircase_error(), "direction: 7 is not a value of Direction");
  every.direction = Direction::Forward;
  EXPECT_EQ(staircase_error(), "tie_break: 6 is not a value of RouteTieBreak");
  every.tie_break = RouteTieBreak::FewerHops;
  EXPECT_EQ(staircase_error(), "queue: 3 is not a value of SpfQueue");
  EXPECT_EQ(error_of([&] { (void)spf(g, std::vector<NodeId> {}, every, SpfQueue::Auto); }), "roots: the node set is empty");
}

// A residual array with several defects at once (one error, from the shared edge-amount check): its
// length before any entry, then its first entry outside [0, capacity] in EdgeId order, before
// min_residual and the Eligibility; in both directions, with and without targets.
TEST(Search, ResidualDefectsReportTheFirstInOrder) {
  const auto g = edge_amount_graph();
  const auto other = make_graph(4, {{0, 1, 2, 1}});
  const auto foreign = other.eligibility({}, std::vector<EdgeId> {0});
  for (const auto& c : edge_amount_defects()) {
    for (auto dir : {Direction::Forward, Direction::Reverse}) {
      EXPECT_EQ(error_of([&] { (void)run(g, {0}, dir, std::nullopt, sel(), &c.amounts, 0, &foreign); }),
                "residual" + c.error);
      EXPECT_EQ(error_of([&] { (void)run(g, {0}, dir, std::vector<NodeId> {3}, sel(false), &c.amounts, 0, &foreign); }),
                "residual" + c.error);
    }
  }
  const auto bounds = edge_amounts_at_the_bounds();
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(), &bounds); }), "");
}

TEST(Search, StagedPathWithASharedBottleneckTakesTheFewestHops) {
  // s=0 -> v=1 cost 1 prio 1; v -> t=4 cost 3 prio 1; v -> a=2 -> b=3 -> t costs 1+1+1 prio 2.
  // Capacities carry the score. Under HIGHER_CAPACITY the bottleneck of every walk is 1 (the shared
  // s->v), so B* = 1 and the canonical walk is the two-hop s, v, t.
  auto g = make_graph(5, {{0, 1, 1, 1}, {1, 4, 1, 3}, {1, 2, 2, 1}, {2, 3, 2, 1}, {3, 4, 2, 1}});
  const auto inv = id_of_input(g);
  for (auto tie : {RouteTieBreak::FewerHops, RouteTieBreak::HigherCapacity}) {
    auto r = run(g, {0}, Direction::Forward, std::vector<NodeId> {4}, sel(false, true, tie));
    EXPECT_EQ(edges_of(r.dag), sorted({inv[0], inv[1]}));
    EXPECT_EQ(r.target_cost, 4);
    // The reverse call with the roles swapped returns the same walk.
    auto rev = run(g, {4}, Direction::Reverse, std::vector<NodeId> {0}, sel(false, true, tie));
    EXPECT_EQ(edges_of(rev.dag), edges_of(r.dag));
  }
}

TEST(Search, StagedPathBottleneckThenHopsThenIds) {
  // Two cost-2 walks 0->1->3 (caps 5, 1) and 0->2->3 (caps 3, 3), and a cost-2 direct edge
  // 0->3 of capacity 2. FEWER_HOPS: fewest hops, the direct edge. HIGHER_CAPACITY: B* = 3
  // (0->2->3), so the two-hop walk.
  auto g = make_graph(4, {{0, 1, 5, 1}, {1, 3, 1, 1}, {0, 2, 3, 1}, {2, 3, 3, 1}, {0, 3, 2, 2}});
  const auto inv = id_of_input(g);
  auto det = run(g, {0}, Direction::Forward, std::vector<NodeId> {3}, sel(false));
  EXPECT_EQ(edges_of(det.dag), (std::vector<EdgeId>{inv[4]}));
  auto pri = run(g, {0}, Direction::Forward, std::vector<NodeId> {3},
                 sel(false, true, RouteTieBreak::HigherCapacity));
  EXPECT_EQ(edges_of(pri.dag), sorted({inv[2], inv[3]}));
  // HIGHER_AVAILABLE_CAPACITY reads the residual: on capacities 5 everywhere and residuals by input
  // [5, 4, 1, 3, 2], 0->1->3 has bottleneck 4, 0->2->3 has 1, the direct edge 2: B* = 4.
  auto g5 = make_graph(4, {{0, 1, 5, 1}, {1, 3, 5, 1}, {0, 2, 5, 1}, {2, 3, 5, 1}, {0, 3, 5, 2}});
  const auto inv5 = id_of_input(g5);
  const auto residual = to_ids(g5, std::vector<Amount> {5, 4, 1, 3, 2});
  auto res = run(g5, {0}, Direction::Forward, std::vector<NodeId> {3},
                 sel(false, true, RouteTieBreak::HigherAvailableCapacity), &residual);
  EXPECT_EQ(edges_of(res.dag), sorted({inv5[0], inv5[1]}));
  // HIGHER_CAPACITY on the same graph: every walk has bottleneck 5, so hops decide.
  auto cap5 = run(g5, {0}, Direction::Forward, std::vector<NodeId> {3},
                  sel(false, true, RouteTieBreak::HigherCapacity), &residual);
  EXPECT_EQ(edges_of(cap5.dag), (std::vector<EdgeId>{inv5[4]}));
  // FEWER_HOPS ignores capacities and residuals entirely.
  auto det_res = run(g5, {0}, Direction::Forward, std::vector<NodeId> {3}, sel(false), &residual);
  EXPECT_EQ(edges_of(det_res.dag), (std::vector<EdgeId>{inv5[4]}));
}

TEST(Search, StagedPathPrefixTrapNeedsTheFullKey) {
  // A cost-2 direct prefix 0->2 with bottleneck 5 against a two-edge
  // cost-2 prefix 0->1->2 with bottleneck 10, then a cost-1 edge of priority 1 to the target 3.
  // Both complete walks cost 3 with bottleneck 1, so hops decide: [0->2, 2->3]. A search that kept
  // only the wider prefix would return the three-hop walk.
  auto g = make_graph(4, {{0, 2, 5, 2}, {0, 1, 10, 1}, {1, 2, 10, 1}, {2, 3, 1, 1}});
  const auto inv = id_of_input(g);
  for (auto dir : {Direction::Forward, Direction::Reverse}) {
    auto r = dir == Direction::Forward
                 ? run(g, {0}, dir, std::vector<NodeId> {3}, sel(false, true, RouteTieBreak::HigherCapacity))
                 : run(g, {3}, dir, std::vector<NodeId> {0}, sel(false, true, RouteTieBreak::HigherCapacity));
    EXPECT_EQ(edges_of(r.dag), sorted({inv[0], inv[3]}));
  }
}

TEST(Search, MultiRootSinglePathStartsAtTheNearestRootByHops) {
  // Roots {0, 1}, both at cost 0; 0 -> 3 in one hop of cost 2 and 1 -> 2 -> 3 in two hops of
  // cost 1 each: the walk of fewest hops (from root 0).
  auto g = make_graph(4, {{0, 3, 1, 2}, {1, 2, 1, 1}, {2, 3, 1, 1}});
  const auto inv = id_of_input(g);
  auto r = run(g, {0, 1}, Direction::Forward, std::vector<NodeId> {3}, sel(false));
  EXPECT_EQ(edges_of(r.dag), (std::vector<EdgeId>{inv[0]}));
}

TEST(Search, StagedPathTiesGoToTheSmallerEdgeIdSequence) {
  // The diamond (0-1-3 and 0-2-3 both cost 10) takes 0-2-3,
  // whose edge-id sequence is smaller; the equal three-hop graph takes 0-1-4-5. Every tie_break
  // (all capacities 1: score ties), both directions.
  auto diamond = make_graph(4, {{0, 1, 1, 9}, {0, 2, 1, 1}, {1, 3, 1, 1}, {2, 3, 1, 9}});
  auto three = make_graph(6, {{0, 2, 1, 1}, {2, 3, 1, 1}, {3, 5, 1, 1}, {0, 1, 1, 1}, {1, 4, 1, 1}, {4, 5, 1, 1}});
  auto nodes_of = [](const Graph& g, const PredecessorDAG& dag) {
    std::set<NodeId> out;
    for (EdgeId e : edges_of(dag)) {
      out.insert(g.edge_sources()[static_cast<std::size_t>(e)]);
      out.insert(g.edge_destinations()[static_cast<std::size_t>(e)]);
    }
    return out;
  };
  for (auto tie : kTieBreaks) {
    {
      const auto s = sel(false, true, tie);
      auto a = run(diamond, {0}, Direction::Forward, std::vector<NodeId> {3}, s);
      EXPECT_EQ(nodes_of(diamond, a.dag), (std::set<NodeId>{0, 2, 3}));
      auto ar = run(diamond, {3}, Direction::Reverse, std::vector<NodeId> {0}, s);
      EXPECT_EQ(edges_of(ar.dag), edges_of(a.dag));
      auto b = run(three, {0}, Direction::Forward, std::vector<NodeId> {5}, s);
      EXPECT_EQ(nodes_of(three, b.dag), (std::set<NodeId>{0, 1, 4, 5}));
      auto br = run(three, {5}, Direction::Reverse, std::vector<NodeId> {0}, s);
      EXPECT_EQ(edges_of(br.dag), edges_of(b.dag));
    }
  }
}

namespace {

// One random query against the reference (search_reference.hpp), on both queues.
struct RandomCase {
  Graph g;
  std::vector<Amount> residual;
  bool use_res;
  Amount min_residual;
  std::optional<Eligibility> elig;
  detail::RouteChoice choice;
  Direction direction;
  std::vector<NodeId> roots;
  std::optional<std::vector<NodeId>> targets;
};

// The expected result of a query by the reference: costs, the DAG's edges, target_cost.
struct Expected {
  std::optional<std::vector<Cost>> costs;
  std::vector<EdgeId> edges;
  std::optional<Cost> target_cost;
};

ref::Query query_of(const RandomCase& c) {
  ref::Query q;
  q.g = &c.g;
  q.elig = c.elig ? &*c.elig : nullptr;
  q.residual = c.use_res ? &c.residual : nullptr;
  q.min_residual = c.min_residual;
  q.parallel_edges = c.choice.parallel_edges;
  q.tie_break = c.choice.tie_break;
  return q;
}

Expected expect_of(const RandomCase& c) {
  const Graph& g = c.g;
  const ref::Query q = query_of(c);
  const auto mem = ref::members(q);
  const auto roots = allowed_nodes(c.roots, q.elig);
  const bool reverse = c.direction == Direction::Reverse;
  Expected out;
  if (!c.targets.has_value()) {
    auto d = ref::distances(g, mem, roots, reverse);
    const auto dag = ref::full_dag(g, mem, d, reverse);
    out.costs = d;
    out.edges.assign(dag.begin(), dag.end());
    return out;
  }
  const auto targets = allowed_nodes(*c.targets, q.elig);
  out.target_cost = kUnreachable;
  if (roots.empty() || targets.empty()) return out;
  const auto& starts = reverse ? targets : roots;
  const auto& ends = reverse ? roots : targets;
  auto [dstar, edges] = ref::pruned(g, mem, starts, ends);
  if (dstar == ref::kNone) return out;
  out.target_cost = dstar;
  if (!c.choice.single_path()) {
    out.edges.assign(edges.begin(), edges.end());
  } else {
    auto walk = ref::canonical(q, mem, starts, ends);
    out.edges = sorted(*walk);
  }
  return out;
}

} // namespace

// A HIGHER_* tie_break with multipath and parallel_edges both true is rejected (nothing to choose);
// such a draw checks the rejection, then runs with FEWER_HOPS.
TEST(Search, RandomQueriesAgreeWithTheReferenceOnBothQueues) {
  std::mt19937 rng(20261006);
  int full = 0;         // full searches that keep an edge
  int multipath = 0;    // multipath target searches that reach a target
  int single = 0;       // single-path target searches that reach a target
  int ineffective = 0;  // draws of the rejected combination
  for (int trial = 0; trial < 3000; ++trial) {
    const int n = 2 + static_cast<int>(rng() % 7);
    const int m = static_cast<int>(rng() % 16);
    const Cost max_edge_cost = (trial % 5 == 0) ? 40 : 4;
    const Cost scale = (trial % 7 == 0) ? 3 : 1;  // a common gcd above 1
    RandomCase c {random_graph(rng, n, m, max_edge_cost, scale), {}, (rng() % 2) == 0, 1 + static_cast<Amount>(rng() % 3),
                  std::nullopt, sel(),
                  (rng() % 2) == 0 ? Direction::Forward : Direction::Reverse, {}, std::nullopt};
    c.residual.resize(static_cast<std::size_t>(c.g.num_edges()));
    for (std::size_t e = 0; e < c.residual.size(); ++e) {
      const Amount cap = c.g.edge_capacity_units()[e];
      c.residual[e] = cap == 0 ? 0 : static_cast<Amount>(rng() % static_cast<std::uint32_t>(cap + 1));
    }
    if (rng() % 3 == 0) {
      std::vector<NodeId> xn;
      std::vector<EdgeId> xe;
      if (rng() % 2 == 0) xn.push_back(static_cast<NodeId>(rng() % static_cast<unsigned>(n)));
      if (c.g.num_edges() > 0 && rng() % 2 == 0) xe.push_back(static_cast<EdgeId>(rng() % static_cast<unsigned>(c.g.num_edges())));
      c.elig = c.g.eligibility(xn, xe);
    }
    const RouteTieBreak higher[3] = {RouteTieBreak::HigherAvailableCapacity, RouteTieBreak::HigherCapacity,
                                     RouteTieBreak::HigherHeadroom};
    const bool with_targets = (rng() % 3) != 0;
    const bool prefer = (rng() % 2) != 0;
    const RouteTieBreak tie = prefer ? higher[rng() % 3] : RouteTieBreak::FewerHops;
    const bool mp = !with_targets || (rng() % 2) == 0;
    c.choice = sel(mp, (rng() % 2) == 0, tie);
    // Disjoint random root and target sets (duplicates allowed in the input).
    std::vector<NodeId> perm(static_cast<std::size_t>(n));
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    const std::size_t nr = 1 + rng() % static_cast<unsigned>(std::max(1, n / 2));
    c.roots.assign(perm.begin(), perm.begin() + static_cast<std::ptrdiff_t>(nr));
    if (rng() % 4 == 0) c.roots.push_back(c.roots.front());
    if (with_targets) {  // never empty: nr <= n / 2 < n
      const std::size_t nt = 1 + rng() % static_cast<unsigned>(std::max<std::size_t>(1, perm.size() - nr));
      c.targets = std::vector<NodeId>(perm.begin() + static_cast<std::ptrdiff_t>(nr),
                                      perm.begin() + static_cast<std::ptrdiff_t>(std::min(perm.size(), nr + nt)));
    }
    if (c.choice.prefers_score() && c.choice.multipath && c.choice.parallel_edges) {
      ASSERT_EQ(error_of([&] {
                  (void)run(c.g, c.roots, c.direction, c.targets, c.choice, c.use_res ? &c.residual : nullptr,
                            c.min_residual, c.elig ? &*c.elig : nullptr);
                }),
                "spf tie_break requires multipath=False or parallel_edges=False")
          << "trial " << trial;
      ++ineffective;
      c.choice.tie_break = RouteTieBreak::FewerHops;
    }
    const auto want = expect_of(c);
    for (auto queue : {SpfQueue::Auto, SpfQueue::Heap}) {
      auto got = run(c.g, c.roots, c.direction, c.targets, c.choice, c.use_res ? &c.residual : nullptr, c.min_residual,
                     c.elig ? &*c.elig : nullptr, queue);
      ASSERT_EQ(got.costs, want.costs) << "trial " << trial;
      ASSERT_EQ(edges_of(got.dag), want.edges) << "trial " << trial << " queue " << static_cast<int>(queue);
      ASSERT_EQ(got.target_cost, want.target_cost) << "trial " << trial;
      expect_canonical(c.g, got.dag);
    }
    full += !c.targets.has_value() && !want.edges.empty() ? 1 : 0;
    const bool reached = c.targets.has_value() && want.target_cost != kUnreachable;
    multipath += reached && !c.choice.single_path() ? 1 : 0;
    single += reached && c.choice.single_path() ? 1 : 0;
  }
  EXPECT_GT(full, 300);
  EXPECT_GT(multipath, 300);
  EXPECT_GT(single, 300);
  EXPECT_GT(ineffective, 300);
}

namespace {

// Whether the hop count decides a single-path query's canonical walk: another walk ties it on
// cost (and bottleneck under a HIGHER_* tie_break) with more hops but a smaller edge-id sequence.
bool hops_decide(const RandomCase& c) {
  const ref::Query q = query_of(c);
  const auto mem = ref::members(q);
  const bool reverse = c.direction == Direction::Reverse;
  const auto roots = allowed_nodes(c.roots, q.elig);
  const auto targets = allowed_nodes(*c.targets, q.elig);
  const auto& starts = reverse ? targets : roots;
  const auto& ends = reverse ? roots : targets;
  const auto canonical = ref::canonical(q, mem, starts, ends);
  if (!canonical) return false;
  const std::vector<EdgeId>& walk = *canonical;
  for (const auto& w : ref::all_walks(c.g, mem, starts, ends)) {
    if (ref::walk_cost(c.g, w) == ref::walk_cost(c.g, walk) &&
        (!q.prefer() || ref::cmp(ref::bottleneck(q, w), ref::bottleneck(q, walk)) == 0) && w.size() > walk.size() &&
        w < walk) {
      return true;
    }
  }
  return false;
}

} // namespace

// The single-path rules that the random comparison above reaches only now and then, on a
// population built for them: a few nodes joined by many parallel edges of cost 1 or 2, a residual
// on every query, one root and one target. Each rule is guarded by the trials in which it decides
// the reference's answer, against a twin query without it: the HIGHER_HEADROOM ratio (the twin
// reads the residual alone), HIGHER_CAPACITY (the twin reads the residual), the bottleneck key of a
// HIGHER_* tie_break (the twin takes FEWER_HOPS, on Expand member graphs, where the tie_break
// decides nothing else), and the hop count (hops_decide). Both queues agree with the reference on
// every trial.
TEST(Search, RandomSinglePathsExerciseEveryKeyRule) {
  std::mt19937 rng(20261010);
  int headroom = 0, capacity = 0, bottleneck = 0, hops = 0;
  for (int trial = 0; trial < 8000; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 4);
    const int m = 6 + static_cast<int>(rng() % 14);
    const RouteTieBreak higher[3] = {RouteTieBreak::HigherAvailableCapacity, RouteTieBreak::HigherCapacity,
                                     RouteTieBreak::HigherHeadroom};
    const RouteTieBreak prio = higher[rng() % 3];
    const bool prefer = rng() % 4 != 0;
    const bool expand = rng() % 2 == 0;
    RandomCase c {random_graph(rng, n, m, 2), {}, true, 1, std::nullopt,
                  sel(false, expand, prefer ? prio : RouteTieBreak::FewerHops),
                  (rng() % 2) == 0 ? Direction::Forward : Direction::Reverse, {}, std::nullopt};
    c.residual.resize(static_cast<std::size_t>(c.g.num_edges()));
    for (std::size_t e = 0; e < c.residual.size(); ++e) {
      const Amount cap = c.g.edge_capacity_units()[e];
      c.residual[e] = cap == 0 ? 0 : static_cast<Amount>(rng() % static_cast<std::uint32_t>(cap + 1));
    }
    const auto root = static_cast<NodeId>(rng() % static_cast<unsigned>(n));
    c.roots = {root};
    c.targets = std::vector<NodeId> {static_cast<NodeId>((root + 1 + rng() % static_cast<unsigned>(n - 1)) % n)};
    const auto want = expect_of(c);
    for (auto queue : {SpfQueue::Auto, SpfQueue::Heap}) {
      auto got = run(c.g, c.roots, c.direction, c.targets, c.choice, &c.residual, c.min_residual, nullptr, queue);
      ASSERT_EQ(edges_of(got.dag), want.edges) << "trial " << trial << " queue " << static_cast<int>(queue);
      ASSERT_EQ(got.target_cost, want.target_cost) << "trial " << trial;
    }
    if (want.target_cost == kUnreachable) continue;
    auto differs = [&](const std::function<void(RandomCase&)>& change) {
      RandomCase twin = c;
      change(twin);
      return expect_of(twin).edges != want.edges ? 1 : 0;
    };
    const auto residual_only = [](RandomCase& t) { t.choice.tie_break = RouteTieBreak::HigherAvailableCapacity; };
    if (prefer && prio == RouteTieBreak::HigherHeadroom) headroom += differs(residual_only);
    if (prefer && prio == RouteTieBreak::HigherCapacity) capacity += differs(residual_only);
    if (prefer && expand) {
      bottleneck += differs([](RandomCase& t) { t.choice.tie_break = RouteTieBreak::FewerHops; });
    }
    hops += hops_decide(c) ? 1 : 0;
  }
  // Each rule decides 42 to 66 of the 8,000 trials with this seed: the floors leave room for a
  // generator change while staying far above the one or two trials of the comparison above.
  EXPECT_GT(headroom, 20);
  EXPECT_GT(capacity, 20);
  EXPECT_GT(bottleneck, 20);
  EXPECT_GT(hops, 20);
}

TEST(Search, DirectionsWithExchangedEndpointsDescribeTheSameWalks) {
  // A forward search from S with targets T and a reverse search rooted at T with targets S give
  // the same pruned DAG, target cost and single path, under FEWER_HOPS and HIGHER_CAPACITY with and
  // without multipath and parallel edges (HIGHER_CAPACITY not with both, where it is rejected; the
  // other tie rules are compared with the reference in RandomQueriesAgreeWithTheReferenceOnBothQueues).
  std::mt19937 rng(7);
  int reached = 0;  // comparisons where a target is reached (a non-empty DAG)
  for (int trial = 0; trial < 500; ++trial) {
    const int n = 4 + static_cast<int>(rng() % 8);
    auto g = random_graph(rng, n, 4 + static_cast<int>(rng() % 20), 3);
    std::vector<NodeId> perm(static_cast<std::size_t>(n));
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    std::vector<NodeId> s(perm.begin(), perm.begin() + 1 + static_cast<std::ptrdiff_t>(rng() % 2));
    std::vector<NodeId> t(perm.end() - 1 - static_cast<std::ptrdiff_t>(rng() % 2), perm.end());
    for (auto tie : {RouteTieBreak::FewerHops, RouteTieBreak::HigherCapacity}) {
      for (bool mp : {true, false}) {
        for (bool me : {true, false}) {
          if (tie != RouteTieBreak::FewerHops && mp && me) continue;
          const auto s1 = sel(mp, me, tie);
          auto f = run(g, s, Direction::Forward, t, s1);
          auto r = run(g, t, Direction::Reverse, s, s1);
          ASSERT_EQ(f.target_cost, r.target_cost) << trial;
          ASSERT_EQ(edges_of(f.dag), edges_of(r.dag)) << trial;
          reached += f.target_cost != kUnreachable ? 1 : 0;
        }
      }
    }
  }
  RecordProperty("reached", reached);
  EXPECT_GT(reached, 1500);  // not vacuous: of 3500 comparisons
}

TEST(Search, TieBreakWithEveryRouteAndEveryParallelEdgeIsRejected) {
  // With multipath and parallel_edges both true the DAG keeps every tight member and every
  // equal-cost walk, so a tie rule has nothing to choose: every HIGHER_* tie_break is rejected,
  // before the residual, and FEWER_HOPS (never a hop-count prune) gives that DAG. Half of the
  // graphs, with and without targets, have a random residual (0 to the capacity).
  std::mt19937 rng(11);
  int nonempty = 0;
  int with_residual = 0;  // nonempty DAGs under a residual
  for (int trial = 0; trial < 300; ++trial) {
    auto g = random_graph(rng, 6, 14, 2);
    std::optional<std::vector<NodeId>> t;
    if (trial % 2 == 0) t = std::vector<NodeId> {5};
    std::vector<Amount> res;
    if (trial % 4 < 2) {
      res.resize(static_cast<std::size_t>(g.num_edges()));
      for (std::size_t e = 0; e < res.size(); ++e) {
        const Amount cap = g.edge_capacity_units()[e];
        res[e] = cap == 0 ? 0 : static_cast<Amount>(rng() % static_cast<std::uint32_t>(cap + 1));
      }
    }
    const std::vector<Amount>* r = res.empty() ? nullptr : &res;
    const auto base = edges_of(run(g, {0}, Direction::Forward, t, sel(), r).dag);
    nonempty += base.empty() ? 0 : 1;
    with_residual += !base.empty() && r != nullptr ? 1 : 0;
    const RandomCase c {g, res, r != nullptr, 1, std::nullopt, sel(), Direction::Forward, {0}, t};
    EXPECT_EQ(base, expect_of(c).edges) << trial;
    for (auto tie : kTieBreaks) {
      if (tie == RouteTieBreak::FewerHops) continue;
      EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, t, sel(true, true, tie), r); }),
                "spf tie_break requires multipath=False or parallel_edges=False")
          << trial;
    }
  }
  RecordProperty("nonempty", nonempty);
  RecordProperty("with_residual", with_residual);
  EXPECT_GT(nonempty, 150);  // not vacuous: of 300 graphs
  EXPECT_GT(with_residual, 50);  // of the 150 graphs with a residual (this seed: 76)
}

TEST(Search, CapacityAndResidualValuesActOnlyThroughTheMemberGraph) {
  // Two graphs with the same edges, in EdgeId order so the ids agree, where the second raises
  // every positive capacity by 1 to 5 units (a zero capacity stays zero). On half of the trials
  // each gets a residual that gates the same edges at min_residual (1 or 2 units) with other
  // values. The member graph is then the same, so the DAG must be: in both directions, with and
  // without targets, under the FEWER_HOPS rules, which read no capacity or residual value: every
  // member and every walk, one member per (u, v) pair, and with targets the single path with and
  // without parallel edges.
  std::mt19937 rng(23);
  int nonempty = 0;
  int gated = 0;  // nonempty DAGs where the residual excluded an edge of at least 1 unit
  for (int trial = 0; trial < 300; ++trial) {
    auto g = random_graph(rng, 6, 14, 2);
    const auto m = static_cast<std::size_t>(g.num_edges());
    std::vector<InEdge> rows;
    for (std::size_t e = 0; e < m; ++e) {
      const Amount cap = g.edge_capacity_units()[e];
      rows.push_back({g.edge_sources()[e], g.edge_destinations()[e], cap == 0 ? 0 : cap + 1 + static_cast<Amount>(rng() % 5),
                      g.edge_costs()[e]});
    }
    auto h = make_graph(6, rows);
    for (std::size_t e = 0; e < m; ++e) {
      ASSERT_EQ(h.edge_sources()[e], g.edge_sources()[e]);
      ASSERT_EQ(h.edge_destinations()[e], g.edge_destinations()[e]);
      ASSERT_EQ(h.edge_costs()[e], g.edge_costs()[e]);
      const Amount cg = g.edge_capacity_units()[e];
      ASSERT_TRUE(cg == 0 ? h.edge_capacity_units()[e] == 0 : h.edge_capacity_units()[e] > cg);
    }
    std::vector<Amount> rg;
    std::vector<Amount> rh;
    const Amount minr = 1 + static_cast<Amount>(rng() % 2);
    bool excluded = false;
    if (trial % 2 == 0) {
      for (std::size_t e = 0; e < m; ++e) {
        const Amount cg = g.edge_capacity_units()[e];
        const Amount ch = h.edge_capacity_units()[e];
        rg.push_back(cg == 0 ? 0 : static_cast<Amount>(rng() % static_cast<std::uint32_t>(cg + 1)));
        if (rg.back() >= minr) {
          rh.push_back(minr + static_cast<Amount>(rng() % static_cast<std::uint32_t>(ch - minr + 1)));
        } else {
          rh.push_back(static_cast<Amount>(rng() % static_cast<std::uint32_t>(std::min(minr - 1, ch) + 1)));
          excluded = excluded || cg >= 1;
        }
      }
    }
    const std::vector<Amount>* res_g = rg.empty() ? nullptr : &rg;
    const std::vector<Amount>* res_h = rh.empty() ? nullptr : &rh;
    std::optional<std::vector<NodeId>> t;
    if (trial % 4 < 2) t = std::vector<NodeId> {5};
    std::vector<detail::RouteChoice> sels {sel(true, true), sel(true, false)};
    if (t.has_value()) {
      sels.push_back(sel(false, true));
      sels.push_back(sel(false, false));
    }
    for (const Direction dir : {Direction::Forward, Direction::Reverse}) {
      for (const auto& s : sels) {
        const auto want = edges_of(run(g, {0}, dir, t, s, res_g, minr).dag);
        EXPECT_EQ(edges_of(run(h, {0}, dir, t, s, res_h, minr).dag), want) << trial;
        nonempty += want.empty() ? 0 : 1;
        gated += !want.empty() && excluded ? 1 : 0;
      }
    }
  }
  RecordProperty("nonempty", nonempty);
  RecordProperty("gated", gated);
  EXPECT_GT(nonempty, 700);  // not vacuous: of 150 x 2 x 4 + 150 x 2 x 2 = 1800 comparisons (this seed: 1018)
  EXPECT_GT(gated, 200);     // of the 900 comparisons under a residual (this seed: 314)
}

TEST(Search, QueueDispatchAtTheWidthLimit) {
  // W = max_edge_cost / gcd + 1: 65,535 / 1 + 1 = 65,536 buckets (the limit) and 65,536 + 1 above it;
  // a large common gcd makes the width 2 big / big + 1 = 3; a graph without edges has width 1.
  // The widths and the dispatch rule (Auto takes the bucket queue up to the limit, Heap never)
  // hold in every build; the queue each search takes is counted by a private hook.
  auto at = make_graph(3, {{0, 1, 1, 65535}, {1, 0, 1, 1}});
  auto above = make_graph(3, {{0, 1, 1, 65536}, {1, 0, 1, 1}});
  const Cost big = Cost {1} << 40;
  auto gcd = make_graph(3, {{0, 1, 1, big}, {1, 0, 1, 2 * big}});
  auto bare = make_graph(2, {});
  EXPECT_EQ(detail::bucket_width(at), 65536);
  EXPECT_EQ(detail::bucket_width(above), 65537);
  EXPECT_EQ(detail::bucket_width(gcd), 3);
  EXPECT_EQ(detail::bucket_width(bare), 1);
  EXPECT_TRUE(detail::uses_buckets(at, SpfQueue::Auto));
  EXPECT_FALSE(detail::uses_buckets(above, SpfQueue::Auto));
  EXPECT_TRUE(detail::uses_buckets(gcd, SpfQueue::Auto));
  EXPECT_TRUE(detail::uses_buckets(bare, SpfQueue::Auto));
  EXPECT_FALSE(detail::uses_buckets(at, SpfQueue::Heap));
  EXPECT_FALSE(detail::uses_buckets(gcd, SpfQueue::Heap));
#if NETGRAPH_CORE_PRIVATE_HOOKS
  auto& counts = detail::thread_queue_counts();
  auto queue_of = [&](const Graph& g, SpfQueue q) {
    const auto before = counts;
    auto r = run(g, {0}, Direction::Forward, std::nullopt, sel(), nullptr, 1, nullptr, q);
    EXPECT_EQ(*r.costs, (std::vector<Cost>{0, g.edge_costs()[0], kUnreachable}));
    return counts.buckets > before.buckets ? std::string("buckets") : std::string("heap");
  };
  EXPECT_EQ(queue_of(at, SpfQueue::Auto), "buckets");
  EXPECT_EQ(queue_of(above, SpfQueue::Auto), "heap");
  EXPECT_EQ(queue_of(gcd, SpfQueue::Auto), "buckets");
  EXPECT_EQ(queue_of(at, SpfQueue::Heap), "heap");
  EXPECT_EQ(queue_of(gcd, SpfQueue::Heap), "heap");
#else
  GTEST_SKIP() << "the queue counters are a private hook";
#endif
}

namespace {

// A textbook Dijkstra over every edge, whatever its capacity (the test's own, for large graphs
// where Bellman-Ford is too slow): the costs and the tight-edge set of a full forward search.
std::pair<std::vector<Cost>, std::vector<EdgeId>> dijkstra(const Graph& g, const std::vector<NodeId>& roots) {
  std::vector<std::vector<std::size_t>> out(static_cast<std::size_t>(g.num_nodes()));
  for (std::size_t e = 0; e < static_cast<std::size_t>(g.num_edges()); ++e) {
    out[static_cast<std::size_t>(g.edge_sources()[e])].push_back(e);
  }
  std::vector<Cost> d(static_cast<std::size_t>(g.num_nodes()), kUnreachable);
  std::priority_queue<std::pair<Cost, NodeId>, std::vector<std::pair<Cost, NodeId>>, std::greater<>> pq;
  for (NodeId r : roots) {
    d[static_cast<std::size_t>(r)] = 0;
    pq.emplace(0, r);
  }
  while (!pq.empty()) {
    auto [du, u] = pq.top();
    pq.pop();
    if (du != d[static_cast<std::size_t>(u)]) continue;
    for (std::size_t e : out[static_cast<std::size_t>(u)]) {
      const auto v = static_cast<std::size_t>(g.edge_destinations()[e]);
      if (du + g.edge_costs()[e] < d[v]) {
        d[v] = du + g.edge_costs()[e];
        pq.emplace(d[v], static_cast<NodeId>(v));
      }
    }
  }
  std::vector<EdgeId> tight;
  for (std::size_t e = 0; e < static_cast<std::size_t>(g.num_edges()); ++e) {
    const Cost du = d[static_cast<std::size_t>(g.edge_sources()[e])];
    const Cost dv = d[static_cast<std::size_t>(g.edge_destinations()[e])];
    if (du != kUnreachable && dv != kUnreachable && du + g.edge_costs()[e] == dv) {
      tight.push_back(static_cast<EdgeId>(e));
    }
  }
  return {d, tight};
}

} // namespace

TEST(Search, BucketWraparoundOnRings) {
  // Small widths with distances many times the width: the cyclic buckets wrap around repeatedly.
  // A plain ring, searched from node 0, runs like a chain with a frontier of one item (the inline
  // item and its promotion); chords, shortcuts and parallel edges on the other three rings
  // promote it. Every result equals the heap's and the textbook Dijkstra's.
  std::mt19937 rng(5);
  for (int shape = 0; shape < 4; ++shape) {
    const int n = 400 + 37 * shape;
    std::vector<InEdge> edges;
    for (int v = 0; v < n; ++v) {
      edges.push_back({v, (v + 1) % n, 1, 1 + static_cast<Cost>(rng() % 3)});  // the ring
      if (shape >= 1 && rng() % 5 == 0) edges.push_back({v, (v + 7) % n, 1, 3});  // chords
      if (shape >= 2 && rng() % 3 == 0) edges.push_back({v, (v + 2) % n, 1, 2});
      if (shape == 3 && rng() % 4 == 0) edges.push_back({v, (v + 1) % n, 1, 1});  // parallel
    }
    auto g = make_graph(n, edges);
    EXPECT_LE(detail::bucket_width(g), 4);
    const auto [d, tight] = dijkstra(g, {0});
    for (auto q : {SpfQueue::Auto, SpfQueue::Heap}) {
      auto r = run(g, {0}, Direction::Forward, std::nullopt, sel(), nullptr, 1, nullptr, q);
      EXPECT_EQ(*r.costs, d) << "shape " << shape;
      EXPECT_EQ(edges_of(r.dag), tight) << "shape " << shape;
    }
  }
}

TEST(Search, EarlyStopsLeaveNothingForTheNextQuery) {
  // A target search stops with pending entries in the queue and touched tentative nodes; the next
  // searches on this thread (other graphs, sizes, directions and queues) must see a clean
  // workspace. Alternating graph sizes also exercise the arrays' growth and shrink.
  std::mt19937 rng(3);
  auto big = random_graph(rng, 3000, 24000, 30);
  auto small = random_graph(rng, 9, 30, 3);
  auto huge = random_graph(rng, 20000, 40000, 5);
  const auto [d_big, t_big] = dijkstra(big, {0});
  const auto [d_small, t_small] = dijkstra(small, {1});
  // The huge reverse search, first on a fresh thread (an empty workspace): each repetition below
  // must give the same.
  std::optional<SpfResult> h_want;
  std::thread([&] { h_want.emplace(run(huge, {5}, Direction::Reverse, std::vector<NodeId> {7}, sel(false))); }).join();
  ASSERT_NE(h_want->target_cost, kUnreachable);
  int stopped = 0;  // target searches that reached their target (and stopped there)
  for (int round = 0; round < 6; ++round) {
    for (auto q : {SpfQueue::Auto, SpfQueue::Heap}) {
      auto stop = run(big, {0}, Direction::Forward, std::vector<NodeId> {static_cast<NodeId>(1 + round)}, sel(), nullptr, 1,
                      nullptr, q);
      EXPECT_EQ(stop.target_cost, d_big[static_cast<std::size_t>(1 + round)]);
      stopped += stop.target_cost != kUnreachable ? 1 : 0;
      auto s = run(small, {1}, Direction::Forward, std::nullopt, sel(), nullptr, 1, nullptr, q);
      EXPECT_EQ(*s.costs, d_small);
      EXPECT_EQ(edges_of(s.dag), t_small);
      if (round % 2 == 0) {
        const auto h = run(huge, {5}, Direction::Reverse, std::vector<NodeId> {7}, sel(false), nullptr, 1, nullptr, q);
        EXPECT_EQ(h.target_cost, h_want->target_cost);
        EXPECT_EQ(edges_of(h.dag), edges_of(h_want->dag));
      }
      auto b = run(big, {0}, Direction::Forward, std::nullopt, sel(), nullptr, 1, nullptr, q);
      EXPECT_EQ(*b.costs, d_big);
      EXPECT_EQ(edges_of(b.dag), t_big);
    }
  }
  EXPECT_EQ(stopped, 12);
}

TEST(Search, LargeGraphsBothQueuesIdentical) {
  // Every SpfResult field identical on both queues, on graphs where Auto takes the buckets: in
  // both directions, the full multipath search and the target searches (multipath, and single
  // path under FEWER_HOPS and HIGHER_CAPACITY; a single-path search needs targets), with and
  // without parallel edges across the shapes.
  std::mt19937 rng(17);
  int reached = 0;
  for (int shape = 0; shape < 3; ++shape) {
    auto g = random_graph(rng, 2000, 16000, shape == 2 ? 1 : 31);
    ASSERT_TRUE(detail::uses_buckets(g, SpfQueue::Auto));
    int k = 0;
    for (const Direction dir : {Direction::Forward, Direction::Reverse}) {
      for (int kind = 0; kind < 4; ++kind, ++k) {
        const auto tie = kind == 3 ? RouteTieBreak::HigherCapacity : RouteTieBreak::FewerHops;
        const auto s = sel(kind < 2, (kind + shape) % 2 == 0, tie);
        std::optional<std::vector<NodeId>> t;
        if (kind > 0) t = std::vector<NodeId> {static_cast<NodeId>(100 + k), static_cast<NodeId>(900 + k)};
        auto a = run(g, {static_cast<NodeId>(k)}, dir, t, s, nullptr, 1, nullptr, SpfQueue::Auto);
        auto h = run(g, {static_cast<NodeId>(k)}, dir, t, s, nullptr, 1, nullptr, SpfQueue::Heap);
        EXPECT_EQ(a.costs, h.costs);
        EXPECT_EQ(a.target_cost, h.target_cost);
        EXPECT_TRUE(a.dag == h.dag);
        reached += t.has_value() && a.target_cost.value_or(kUnreachable) != kUnreachable;
      }
    }
  }
  EXPECT_GT(reached, 12);  // of the 18 target searches: not comparisons of unreached targets
}

TEST(Search, InternalEntriesMatchThePublicSearch) {
  // route_search (T1's CSPF: the single path, in walk order), support_route (ShortestPaths on the cached
  // support) and
  // static_costs (rosters) against spf on the same questions; support_route under HIGHER_* scores
  // read on a residual view against the independent reference (search_reference.hpp). A HIGHER_*
  // tie_break never meets multipath and parallel_edges both true (the configuration rejects it).
  std::mt19937 rng(23);
  int routed = 0;        // gated route searches that reach a destination
  int projected = 0;     // static projections that reach a destination
  int residual_views = 0;
  int view_changes = 0;  // the residual view selected another route than the static one
  for (int trial = 0; trial < 400; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 8);
    auto g = random_graph(rng, n, 3 + static_cast<int>(rng() % 22), 3);
    const NodeId src = 0;
    std::vector<NodeId> dests = {static_cast<NodeId>(n - 1)};
    if (n > 4 && rng() % 2 == 0) dests.insert(dests.begin(), static_cast<NodeId>(n - 2));
    std::vector<Amount> residual(static_cast<std::size_t>(g.num_edges()));
    for (std::size_t e = 0; e < residual.size(); ++e) {
      const Amount cap = g.edge_capacity_units()[e];
      residual[e] = cap == 0 ? 0 : static_cast<Amount>(rng() % static_cast<std::uint32_t>(cap + 1));
    }
    for (bool mp : {true, false}) {
      for (bool me : {true, false}) {
        for (bool prefer : {false, true}) {
          if (prefer && mp && me) continue;
          const auto s = sel(mp, me, prefer ? RouteTieBreak::HigherAvailableCapacity : RouteTieBreak::FewerHops);
          // T1's CSPF: gated at 2 units on the residual, the single path (it reads no multipath)
          // in walk order from src.
          if (!mp) {
            auto want = run(g, {src}, Direction::Forward, dests, s, &residual, 2);
            auto got = detail::route_search(g, detail::SearchFilter {nullptr, residual, 2}, s, src, dests);
            if (want.target_cost == kUnreachable) {
              EXPECT_FALSE(got.has_value());
            } else {
              ASSERT_TRUE(got.has_value());
              EXPECT_EQ(got->cost, *want.target_cost);
              NodeId at = src;
              for (EdgeId e : got->edges) {
                EXPECT_EQ(g.edge_sources()[static_cast<std::size_t>(e)], at);
                at = g.edge_destinations()[static_cast<std::size_t>(e)];
              }
              auto sorted = got->edges;
              std::sort(sorted.begin(), sorted.end());
              EXPECT_EQ(sorted, edges_of(want.dag));
              ++routed;
            }
          }
          // ShortestPaths: the static support toward dests projected with static scores equals the
          // static search from src (no residual: HIGHER_AVAILABLE_CAPACITY reads the capacities).
          const auto support = detail::static_support(g, dests, Direction::Reverse, nullptr);
          const detail::PriorityValues values {&g, s.tie_break, {}};
          auto proj = detail::support_route(g, detail::PredecessorDAGBuilder::ascending(support), src, dests, s,
                                            detail::PriorityReader::of(values));
          auto stat = run(g, {src}, Direction::Forward, dests, s);
          if (stat.target_cost == kUnreachable) {
            EXPECT_TRUE(proj.edges.empty());
          } else {
            EXPECT_EQ(proj.edges, edges_of(stat.dag)) << trial;
            EXPECT_EQ(proj.cost, *stat.target_cost);
            ++projected;
          }
          // ShortestPaths creation reads its scores on a residual view without gating on it: the
          // projection of the static support under each tie_break read on `residual` equals the
          // independent reference's walks (or canonical path) over the static member graph with
          // that tie_break, gated at no residual (min_residual 0).
          const std::vector<RouteTieBreak> ties =
              prefer ? std::vector<RouteTieBreak> {RouteTieBreak::HigherAvailableCapacity, RouteTieBreak::HigherHeadroom,
                                                   RouteTieBreak::HigherCapacity}
                     : std::vector<RouteTieBreak> {RouteTieBreak::FewerHops};
          for (auto prio : ties) {
            const auto sp = sel(mp, me, prio);
            const detail::PriorityValues on_residual {&g, prio, residual};
            auto rproj = detail::support_route(g, detail::PredecessorDAGBuilder::ascending(support), src, dests, sp,
                                               detail::PriorityReader::of(on_residual));
            ref::Query q;
            q.g = &g;
            q.residual = &residual;
            q.min_residual = 0;
            q.parallel_edges = me;
            q.tie_break = prio;
            const auto mem = ref::members(q);
            const auto [dstar, walks] = ref::pruned(g, mem, {src}, dests);
            if (dstar == ref::kNone) {
              EXPECT_TRUE(rproj.edges.empty()) << trial;
              continue;
            }
            std::vector<EdgeId> want(walks.begin(), walks.end());
            if (!mp) want = sorted(*ref::canonical(q, mem, {src}, dests));
            EXPECT_EQ(rproj.edges, want) << "trial " << trial << " tie_break " << static_cast<int>(prio);
            EXPECT_EQ(rproj.cost, dstar) << trial;
            ++residual_views;
            if (rproj.edges != proj.edges) ++view_changes;
          }
        }
      }
    }
    const auto costs = detail::static_costs(g, nullptr, dests, std::vector<NodeId> {0, 1});
    auto rev = run(g, dests, Direction::Reverse, std::nullopt, sel());
    EXPECT_EQ(costs, (std::vector<Cost>{(*rev.costs)[0], (*rev.costs)[1]}));
  }
  EXPECT_GT(routed, 300);           // not vacuous: many queries reach a destination
  EXPECT_GT(projected, 1000);
  EXPECT_GT(residual_views, 3000);
  EXPECT_GT(view_changes, 30);      // and the residual view decides some of them
}

TEST(Search, ReverseSupportsOrderEveryTailsEdgesAtAnyDegree) {
  // A reverse full search emits each tail's tight out-edges in ascending EdgeId, by one
  // bit per out-edge id up to 64 out-edges and by sorting above: hubs with 1, 63, 64, 65 and 129
  // out-edges, tight (cost 1, parallel copies, ties; a fifth of capacity 0, tight like any other)
  // and not (cost 2), and hubs of
  // 64, 65 and 66 out-edges that are all tight, against the reference's full DAG, with and
  // without an Eligibility.
  // Ids rank a node's costlier edges last, so the all-tight variant makes the last ids tight too.
  std::mt19937 rng(64);
  for (const auto& [degree, all_tight] : std::vector<std::pair<int, bool>> {
           {1, false}, {63, false}, {64, false}, {65, false}, {129, false}, {64, true}, {65, true}, {66, true}}) {
    const NodeId root = 0;
    const NodeId hub = 1;
    const int spokes = 40;  // nodes 2 .. 41, each one hop from the root
    std::vector<InEdge> edges;
    for (int k = 0; k < degree; ++k) {
      const auto to = static_cast<NodeId>(2 + rng() % spokes);
      const Cost cost = !all_tight && rng() % 4 == 0 ? 2 : 1;
      edges.push_back({hub, to, static_cast<Amount>(!all_tight && rng() % 5 == 0 ? 0 : 3), cost});
    }
    for (NodeId s = 2; s < 2 + spokes; ++s) edges.push_back({s, root, 4, 1});
    const Graph g = make_graph(2 + spokes, edges);
    ASSERT_EQ(g.outgoing_offsets()[2] - g.outgoing_offsets()[1], degree);
    for (bool with_elig : {false, true}) {
      std::optional<Eligibility> elig;
      if (with_elig) elig = g.eligibility(std::vector<NodeId> {5}, std::vector<EdgeId> {g.outgoing_offsets()[1]});
      ref::Query q;
      q.g = &g;
      q.elig = elig ? &*elig : nullptr;
      const auto mem = ref::members(q);
      const auto d = ref::distances(g, mem, {root}, true);
      const auto want = ref::full_dag(g, mem, d, true);
      const auto got = detail::static_support(g, std::vector<NodeId> {root}, Direction::Reverse, q.elig);
      const auto asc = detail::PredecessorDAGBuilder::ascending(got);
      EXPECT_EQ(std::vector<EdgeId>(asc.begin(), asc.end()), std::vector<EdgeId>(want.begin(), want.end()))
          << "degree " << degree << " all tight " << all_tight << " elig " << with_elig;
      if (all_tight && !with_elig) EXPECT_EQ(asc.size(), static_cast<std::size_t>(degree + spokes));
      EXPECT_TRUE(std::is_sorted(asc.begin(), asc.end()));
      // The hub's own entries, whose order is the subject here: at least half its edges (none when
      // the Eligibility masks the degree-1 hub's only edge).
      const auto hub_entries = std::count_if(asc.begin(), asc.end(), [&](EdgeId e) { return g.edge_sources()[e] == hub; });
      EXPECT_GE(hub_entries, with_elig && degree == 1 ? 0 : (degree + 1) / 2)
          << "degree " << degree << " all tight " << all_tight << " elig " << with_elig;
    }
  }
}

TEST(Search, StaticSupportTakesNormalisedRootsAndAllowsNone) {
  // static_support and static_costs take root sets their callers normalised
  // (ascending, no duplicates, in range) and filter the Eligibility themselves. On random graphs
  // with masks they equal the public search on the same roots: the full support in both
  // directions, the reverse costs at every node. With no allowed root the support is the empty
  // DAG and every cost kUnreachable, as the public search gives (and static_costs of no
  // destination at all). The public entries keep normalising: spf and a cache lookup on the
  // roots unsorted and repeated give the same support.
  std::mt19937 rng(813);
  int masked_roots = 0;  // supports with an edge although the Eligibility removes a root
  for (int trial = 0; trial < 300; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 9);
    const auto g = random_graph(rng, n, 2 + static_cast<int>(rng() % 30), 4);
    std::vector<NodeId> roots;
    for (NodeId v = 0; v < n; ++v) {
      if (rng() % 3 == 0) roots.push_back(v);
    }
    if (roots.empty()) roots.push_back(static_cast<NodeId>(rng() % static_cast<unsigned>(n)));
    std::vector<NodeId> xn;
    for (NodeId v = 0; v < n; ++v) {
      if (rng() % 4 == 0) xn.push_back(v);
    }
    if (trial % 10 == 0) xn = roots;  // no root allowed
    std::vector<EdgeId> xe;
    for (EdgeId e = 0; e < g.num_edges(); ++e) {
      if (rng() % 6 == 0) xe.push_back(e);
    }
    const auto el = g.eligibility(xn, xe);
    const Eligibility* elig = trial % 2 == 0 ? &el : nullptr;
    const bool any_allowed =
        std::any_of(roots.begin(), roots.end(), [&](NodeId r) { return elig == nullptr || elig->allows_node(r); });
    const bool masked =
        elig != nullptr && std::any_of(roots.begin(), roots.end(), [&](NodeId r) { return !elig->allows_node(r); });
    std::vector<NodeId> all(static_cast<std::size_t>(n));
    std::iota(all.begin(), all.end(), 0);
    for (Direction dir : {Direction::Forward, Direction::Reverse}) {
      const auto got = detail::static_support(g, roots, dir, elig);
      const auto [want_costs, want_dag, want_cost] = run(g, roots, dir, std::nullopt, sel(), nullptr, 1, elig);
      EXPECT_TRUE(got == want_dag) << "trial " << trial << " direction " << static_cast<int>(dir);
      if (!any_allowed) EXPECT_TRUE(got == detail::PredecessorDAGBuilder::empty(g)) << trial;
      if (dir == Direction::Forward && masked && !detail::PredecessorDAGBuilder::ascending(got).empty()) ++masked_roots;
    }
    const auto costs = detail::static_costs(g, elig, roots, all);
    const auto [rev_costs, rev_dag, rev_cost] = run(g, roots, Direction::Reverse, std::nullopt, sel(), nullptr, 1, elig);
    EXPECT_EQ(costs, *rev_costs) << trial;
    if (!any_allowed) EXPECT_EQ(costs, std::vector<Cost>(all.size(), kUnreachable)) << trial;  // every tenth trial, and an even one whose exclusions remove every root
    EXPECT_EQ(detail::static_costs(g, elig, std::vector<NodeId> {}, all), std::vector<Cost>(all.size(), kUnreachable));
    // The public entries on the same roots reversed and repeated.
    std::vector<NodeId> messy(roots.rbegin(), roots.rend());
    messy.push_back(roots.front());
    ShortestPathCache cache(g);
    const auto cached = cache.lookup(g, messy, Direction::Reverse, elig);
    EXPECT_TRUE(*cached == detail::static_support(g, roots, Direction::Reverse, elig)) << trial;
    EXPECT_TRUE(cache.lookup(g, roots, Direction::Reverse, elig) == cached) << trial;  // the same entry
    EXPECT_TRUE(rev_dag == *cached) << trial;
    const auto [messy_costs, messy_dag, messy_cost] = run(g, messy, Direction::Reverse, std::nullopt, sel(), nullptr, 1, elig);
    EXPECT_TRUE(messy_dag == *cached) << trial;
  }
  EXPECT_GT(masked_roots, 10);  // not vacuous: the Eligibility removes some roots of real supports
}

// The route choice. TRI: three parallel cost-1 links 0 -> 1 of capacities 4, 10,
// 100 and a cost-2 link of 50, availability 4, 8, 5, 50. The single path takes e0 under
// FEWER_HOPS, e2 under HIGHER_CAPACITY, e1 under HIGHER_AVAILABLE_CAPACITY and e0 under
// HIGHER_HEADROOM (1 against 4/5 and 1/20); e3 never.
TEST(SpfTie, Parallel) {
  const auto g = make_graph(2, {{0, 1, 4, 1}, {0, 1, 10, 1}, {0, 1, 100, 1}, {0, 1, 50, 2}});
  const std::vector<Amount> available = {4, 8, 5, 50};
  for (const auto& [tie, want] : {std::pair {RouteTieBreak::FewerHops, EdgeId {0}},
                                  std::pair {RouteTieBreak::HigherCapacity, EdgeId {2}},
                                  std::pair {RouteTieBreak::HigherAvailableCapacity, EdgeId {1}},
                                  std::pair {RouteTieBreak::HigherHeadroom, EdgeId {0}}}) {
    const auto single = run(g, {0}, Direction::Forward, std::vector<NodeId> {1}, sel(false, true, tie), &available);
    EXPECT_EQ(edges_of(single.dag), (std::vector<EdgeId> {want})) << static_cast<int>(tie);
    // The same rule keeps one member of the pair in a multipath search without parallel edges.
    const auto collapsed = run(g, {0}, Direction::Forward, std::nullopt, sel(true, false, tie), &available);
    EXPECT_EQ(edges_of(collapsed.dag), (std::vector<EdgeId> {want})) << static_cast<int>(tie);
  }
}

// TIE-routes: the direct link e1 (capacity 3, cost 2) against e0, e2 (capacity 10, cost 1 each).
// Without a residual HIGHER_HEADROOM reads 1 on every positive link, so hops decide. TIE-
// bottleneck: e2 at capacity 1 makes the two-hop bottleneck 1, and every rule takes e1 although
// the two-hop route's first link is the larger one (the paths are compared whole).
TEST(SpfTie, Routes) {
  for (const Amount e2 : {Amount {10}, Amount {1}}) {
    const auto g = make_graph(3, {{0, 1, 10, 1}, {0, 2, 3, 2}, {1, 2, e2, 1}});
    const std::vector<EdgeId> two_hops = e2 == 10 ? std::vector<EdgeId> {0, 2} : std::vector<EdgeId> {1};
    for (const auto& [tie, want] : {std::pair {RouteTieBreak::FewerHops, std::vector<EdgeId> {1}},
                                    std::pair {RouteTieBreak::HigherCapacity, two_hops},
                                    std::pair {RouteTieBreak::HigherAvailableCapacity, two_hops},
                                    std::pair {RouteTieBreak::HigherHeadroom, std::vector<EdgeId> {1}}}) {
      const auto r = run(g, {0}, Direction::Forward, std::vector<NodeId> {2}, sel(false, true, tie));
      EXPECT_EQ(edges_of(r.dag), want) << e2 << " " << static_cast<int>(tie);
    }
  }
}

// spf's entry: the ineffective tie rule after the single-path rule, before the residual.
TEST(SpfValidation, TieRule) {
  const auto g = make_graph(3, {{0, 1, 2, 1}, {1, 2, 2, 1}});
  const std::string ineffective = "spf tie_break requires multipath=False or parallel_edges=False";
  const std::vector<Amount> short_residual = {1};
  const auto hc = RouteTieBreak::HigherCapacity;
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(true, true, hc)); }), ineffective);
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(true, true, hc), &short_residual); }),
            ineffective);
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(false, true, hc)); }),
            "multipath = false requires targets (the single path joins a root to a target)");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(true, false, hc)); }), "");
  EXPECT_EQ(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel()); }), "");
}
