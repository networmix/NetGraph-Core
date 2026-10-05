// The internal search of APPROACH 6 (sub-PR 1b; public spf in 1d): the member graph and
// select_adjacency, the unique DAG in both directions, anycast pruning, the staged canonical
// single path on the review witnesses and against exhaustive enumeration, validation.
#include <gtest/gtest.h>

#include <algorithm>
#include <functional>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "netgraph/core/detail/search.hpp"
#include "netgraph/core/flow.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using netgraph::core::detail::Direction;

namespace {

PathSelection sel(bool multipath = true, bool multi_edge = true,
                  EdgeTieBreak tie = EdgeTieBreak::Deterministic, EdgePriority prio = EdgePriority::Residual) {
  return PathSelection(tie, prio, multipath, multi_edge);
}

std::vector<EdgeId> members_of(const std::vector<std::uint8_t>& mask) {
  std::vector<EdgeId> out;
  for (std::size_t e = 0; e < mask.size(); ++e) {
    if (mask[e] != 0) out.push_back(static_cast<EdgeId>(e));
  }
  return out;
}

std::string error_of(const std::function<void()>& fn) {
  try {
    fn();
  } catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
}

detail::SpfOut run(const StrictMultiDiGraph& g, std::vector<NodeId> roots, Direction dir,
                   std::optional<std::vector<NodeId>> targets, const PathSelection& s,
                   std::optional<std::vector<Amount>> residual = std::nullopt, Amount min_residual = 1,
                   const Eligibility* elig = nullptr) {
  std::optional<std::span<const NodeId>> t;
  if (targets.has_value()) t = std::span<const NodeId>(*targets);
  std::optional<std::span<const Amount>> r;
  if (residual.has_value()) r = std::span<const Amount>(*residual);
  return detail::spf(g, roots, dir, t, s, r, min_residual, elig);
}

} // namespace

TEST(Search, MemberGraphStaticEligibilityCapacityAndResidual) {
  // ids: 0: 0->1 cap 0, 1: 0->1 cap 3, 2: 0->1 cap 5 (cost 2), 3: 1->2 cap 4
  auto g = make_graph(3, {{0, 1, 0, 1}, {0, 1, 3, 1}, {0, 1, 5, 2}, {1, 2, 4, 1}});
  auto m = detail::member_mask(g, nullptr, sel(), {}, 1);
  EXPECT_EQ(members_of(m), (std::vector<EdgeId>{1, 2, 3}));  // capacity < 1 unit never routable
  const std::vector<Amount> res = {0, 1, 5, 4};
  EXPECT_EQ(members_of(detail::member_mask(g, nullptr, sel(), res, 1)), (std::vector<EdgeId>{1, 2, 3}));
  EXPECT_EQ(members_of(detail::member_mask(g, nullptr, sel(), res, 2)), (std::vector<EdgeId>{2, 3}));
  const std::vector<EdgeId> xe = {3};
  auto el = g.eligibility({}, xe);
  EXPECT_EQ(members_of(detail::member_mask(g, &el, sel(), {}, 1)), (std::vector<EdgeId>{1, 2}));
}

TEST(Search, CollapseBySelectAdjacency) {
  // Parallel 0->1: id 0 cost 1 cap 2, id 1 cost 1 cap 6, id 2 cost 1 cap 4, id 3 cost 2 cap 9.
  auto g = make_graph(2, {{0, 1, 2, 1}, {0, 1, 6, 1}, {0, 1, 4, 1}, {0, 1, 9, 2}});
  // Deterministic: minimum cost, then the lowest EdgeId.
  EXPECT_EQ(members_of(detail::member_mask(g, nullptr, sel(true, false), {}, 1)), (std::vector<EdgeId>{0}));
  // PreferHigherPriority on Capacity: the highest capacity among the cheapest.
  auto cap = sel(true, false, EdgeTieBreak::PreferHigherPriority, EdgePriority::Capacity);
  EXPECT_EQ(members_of(detail::member_mask(g, nullptr, cap, {}, 1)), (std::vector<EdgeId>{1}));
  // Residual priority reads the residual; equal residual keeps the lowest id.
  auto resid = sel(true, false, EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual);
  EXPECT_EQ(members_of(detail::member_mask(g, nullptr, resid, std::vector<Amount> {2, 1, 2, 9}, 1)),
            (std::vector<EdgeId>{0}));
  EXPECT_EQ(members_of(detail::member_mask(g, nullptr, resid, std::vector<Amount> {1, 1, 3, 9}, 1)),
            (std::vector<EdgeId>{2}));
  // Headroom compares residual / capacity exactly: 1/2 < 4/6 = 2/3 > 2/4.
  auto head = sel(true, false, EdgeTieBreak::PreferHigherPriority, EdgePriority::Headroom);
  EXPECT_EQ(members_of(detail::member_mask(g, nullptr, head, std::vector<Amount> {1, 4, 2, 9}, 1)),
            (std::vector<EdgeId>{1}));
  // Equal rational headroom (2/4 == 1/2): the lowest EdgeId.
  EXPECT_EQ(members_of(detail::member_mask(g, nullptr, head, std::vector<Amount> {1, 1, 2, 9}, 1)),
            (std::vector<EdgeId>{0}));
  // Without a residual every priority reads static capacities.
  EXPECT_EQ(members_of(detail::member_mask(g, nullptr, resid, {}, 1)), (std::vector<EdgeId>{1}));
}

TEST(Search, FullDagForwardAndReverse) {
  // Diamond with a parallel member and a longer branch: 0->1, 0->1, 0->2 (cost 2), 1->3, 2->3.
  auto g = make_graph(4, {{0, 1, 1, 1}, {0, 1, 1, 1}, {0, 2, 1, 2}, {1, 3, 1, 1}, {2, 3, 1, 1}});
  auto f = run(g, {0}, Direction::Forward, std::nullopt, sel());
  ASSERT_TRUE(f.costs.has_value());
  EXPECT_EQ(*f.costs, (std::vector<Cost>{0, 1, 2, 2}));
  // Tight edges: both parallel members, 0->2 (d(2) = 2) and 1->3; 2->3 is not tight (2 + 1 > 2).
  EXPECT_EQ(f.edges, (std::vector<EdgeId>{0, 1, 2, 3}));
  auto r = run(g, {3}, Direction::Reverse, std::nullopt, sel());
  EXPECT_EQ(*r.costs, (std::vector<Cost>{2, 1, 1, 0}));
  // The reverse DAG toward 3 is the forward-oriented set of tight edges toward 3.
  EXPECT_EQ(r.edges, (std::vector<EdgeId>{0, 1, 3, 4}));
  // Unreached nodes keep kUnreachable.
  auto u = run(g, {2}, Direction::Forward, std::nullopt, sel());
  EXPECT_EQ(*u.costs, (std::vector<Cost>{kUnreachable, kUnreachable, 0, 1}));
}

TEST(Search, AnycastPruningAndDirectionIndependence) {
  // 0->1 (1), 1->2 (1), 0->3 (2), 3->2 (5), 1->4 (5): targets {2, 4}: nearest is 2 at cost 2.
  auto g = make_graph(5, {{0, 1, 1, 1}, {1, 2, 1, 1}, {0, 3, 1, 2}, {3, 2, 1, 5}, {1, 4, 1, 5}});
  const auto inv = id_of_input(g);
  auto f = run(g, {0}, Direction::Forward, std::vector<NodeId> {2, 4}, sel());
  ASSERT_TRUE(f.target_cost.has_value());
  EXPECT_EQ(*f.target_cost, 2);
  EXPECT_FALSE(f.costs.has_value());
  std::vector<EdgeId> want = {inv[0], inv[1]};
  std::sort(want.begin(), want.end());
  EXPECT_EQ(f.edges, want);
  // A reverse search rooted at the targets with the source as its target: the same edges.
  auto r = run(g, {2, 4}, Direction::Reverse, std::vector<NodeId> {0}, sel());
  EXPECT_EQ(r.edges, f.edges);
  EXPECT_EQ(r.target_cost, f.target_cost);
  // No reachable target: nullopt, no edges.
  auto none = run(g, {2}, Direction::Forward, std::vector<NodeId> {0}, sel());
  EXPECT_FALSE(none.target_cost.has_value());
  EXPECT_TRUE(none.edges.empty());
}

TEST(Search, ExcludedRootsAndTargetsAreSkipped) {
  auto g = make_graph(3, {{0, 1, 1, 1}, {1, 2, 1, 1}});
  const std::vector<NodeId> xn = {0};
  auto el = g.eligibility(xn, {});
  auto r = run(g, {0}, Direction::Forward, std::vector<NodeId> {2}, sel(), std::nullopt, 1, &el);
  EXPECT_FALSE(r.target_cost.has_value());
  EXPECT_TRUE(r.edges.empty());
  auto full = run(g, {0}, Direction::Forward, std::nullopt, sel(), std::nullopt, 1, &el);
  EXPECT_EQ(*full.costs, (std::vector<Cost>{kUnreachable, kUnreachable, kUnreachable}));
}

TEST(Search, ValidationBeforeAnyWork) {
  auto g = make_graph(3, {{0, 1, 2, 1}, {1, 2, 2, 1}});
  auto other = make_graph(3, {{0, 1, 2, 1}});
  auto foreign = other.eligibility({}, std::vector<EdgeId> {0});
  EXPECT_NE(error_of([&] { (void)run(g, {}, Direction::Forward, std::nullopt, sel()); }).find("roots"), std::string::npos);
  EXPECT_NE(error_of([&] { (void)run(g, {3}, Direction::Forward, std::nullopt, sel()); }).find("out of node range"),
            std::string::npos);
  EXPECT_NE(error_of([&] { (void)run(g, {0}, Direction::Forward, std::vector<NodeId> {}, sel()); }).find("targets"),
            std::string::npos);
  EXPECT_NE(error_of([&] { (void)run(g, {0, 1}, Direction::Forward, std::vector<NodeId> {1}, sel()); })
                .find("both endpoint sets"),
            std::string::npos);
  EXPECT_NE(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(false)); }).find("requires targets"),
            std::string::npos);
  EXPECT_NE(error_of([&] {
              (void)run(g, {0}, Direction::Forward, std::nullopt, sel(), std::vector<Amount> {1});
            }).find("residual has length"),
            std::string::npos);
  EXPECT_NE(error_of([&] {
              (void)run(g, {0}, Direction::Forward, std::nullopt, sel(), std::vector<Amount> {1, 3});
            }).find("outside [0, capacity_units"),
            std::string::npos);
  EXPECT_NE(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(), std::nullopt, 0); })
                .find("min_residual"),
            std::string::npos);
  EXPECT_NE(error_of([&] { (void)run(g, {0}, Direction::Forward, std::nullopt, sel(), std::nullopt, 1, &foreign); })
                .find("foreign Eligibility"),
            std::string::npos);
}

TEST(Search, StagedPathReview13D3Witness) {
  // s=0 -> v=1 cost 1 prio 1; v -> t=4 cost 3 prio 1; v -> a=2 -> b=3 -> t costs 1+1+1 prio 2.
  // Capacities carry the priority. Under PreferHigherPriority the bottleneck of every walk is 1
  // (the shared s->v), so B* = 1 and the canonical walk is the two-hop s, v, t.
  auto g = make_graph(5, {{0, 1, 1, 1}, {1, 4, 1, 3}, {1, 2, 2, 1}, {2, 3, 2, 1}, {3, 4, 2, 1}});
  const auto inv = id_of_input(g);
  for (auto tie : {EdgeTieBreak::Deterministic, EdgeTieBreak::PreferHigherPriority}) {
    auto r = run(g, {0}, Direction::Forward, std::vector<NodeId> {4},
                 sel(false, true, tie, EdgePriority::Capacity));
    ASSERT_TRUE(r.path.has_value());
    EXPECT_EQ(r.path->walk, (std::vector<EdgeId>{inv[0], inv[1]}));
    EXPECT_EQ(r.path->cost, 4);
    // The reverse call with the roles swapped returns the same walk.
    auto rev = run(g, {4}, Direction::Reverse, std::vector<NodeId> {0}, sel(false, true, tie, EdgePriority::Capacity));
    ASSERT_TRUE(rev.path.has_value());
    EXPECT_EQ(rev.path->walk, r.path->walk);
  }
}

TEST(Search, StagedPathBottleneckThenHopsThenIds) {
  // Two cost-2 walks 0->1->3 (caps 5, 1) and 0->2->3 (caps 3, 3), and a cost-2 direct edge
  // 0->3 of capacity 2. Deterministic: fewest hops, the direct edge. PreferHigherPriority on
  // Capacity: B* = 3 (0->2->3), so the two-hop walk.
  auto g = make_graph(4, {{0, 1, 5, 1}, {1, 3, 1, 1}, {0, 2, 3, 1}, {2, 3, 3, 1}, {0, 3, 2, 2}});
  const auto inv = id_of_input(g);
  auto det = run(g, {0}, Direction::Forward, std::vector<NodeId> {3}, sel(false));
  EXPECT_EQ(det.path->walk, (std::vector<EdgeId>{inv[4]}));
  auto pri = run(g, {0}, Direction::Forward, std::vector<NodeId> {3},
                 sel(false, true, EdgeTieBreak::PreferHigherPriority, EdgePriority::Capacity));
  EXPECT_EQ(pri.path->walk, (std::vector<EdgeId>{inv[2], inv[3]}));
  ASSERT_TRUE(pri.path->b_star.has_value());
  EXPECT_EQ(pri.path->b_star->num, 3);
  // Residual priority reads the residual: on capacities 5 everywhere and residuals by input
  // [5, 4, 1, 3, 2], 0->1->3 has bottleneck 4, 0->2->3 has 1, the direct edge 2: B* = 4.
  auto g5 = make_graph(4, {{0, 1, 5, 1}, {1, 3, 5, 1}, {0, 2, 5, 1}, {2, 3, 5, 1}, {0, 3, 5, 2}});
  const auto inv5 = id_of_input(g5);
  const auto residual = to_ids(g5, std::vector<Amount> {5, 4, 1, 3, 2});
  auto res = run(g5, {0}, Direction::Forward, std::vector<NodeId> {3},
                 sel(false, true, EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual), residual);
  EXPECT_EQ(res.path->walk, (std::vector<EdgeId>{inv5[0], inv5[1]}));
  // Capacity priority on the same graph: every walk has bottleneck 5, so hops decide.
  auto cap5 = run(g5, {0}, Direction::Forward, std::vector<NodeId> {3},
                  sel(false, true, EdgeTieBreak::PreferHigherPriority, EdgePriority::Capacity), residual);
  EXPECT_EQ(cap5.path->walk, (std::vector<EdgeId>{inv5[4]}));
  // Deterministic ignores capacities and residuals entirely.
  auto det_res = run(g5, {0}, Direction::Forward, std::vector<NodeId> {3}, sel(false), residual);
  EXPECT_EQ(det_res.path->walk, (std::vector<EdgeId>{inv5[4]}));
}

namespace {

// Every simple path from a source to a destination over the members (tiny graphs).
void enumerate(const StrictMultiDiGraph& g, const std::vector<std::uint8_t>& member, NodeId u,
               const std::vector<NodeId>& dsts, std::vector<EdgeId>& walk, std::vector<std::uint8_t>& on,
               std::vector<std::vector<EdgeId>>& out) {
  if (std::find(dsts.begin(), dsts.end(), u) != dsts.end()) {
    out.push_back(walk);
    return;
  }
  const auto row = g.row_offsets();
  for (auto e = row[static_cast<std::size_t>(u)]; e < row[static_cast<std::size_t>(u) + 1]; ++e) {
    if (member[static_cast<std::size_t>(e)] == 0) continue;
    const NodeId v = g.edge_dst()[static_cast<std::size_t>(e)];
    if (on[static_cast<std::size_t>(v)] != 0) continue;
    on[static_cast<std::size_t>(v)] = 1;
    walk.push_back(e);
    enumerate(g, member, v, dsts, walk, on, out);
    walk.pop_back();
    on[static_cast<std::size_t>(v)] = 0;
  }
}

} // namespace

TEST(Search, StagedPathAgainstEnumeration) {
  // The canonical key (cost, -bottleneck under PreferHigherPriority, hops, edge-id sequence)
  // over every simple walk of the member graph, against the staged construction.
  std::mt19937 rng(20261004);
  int checked = 0;
  for (int trial = 0; trial < 2000; ++trial) {
    const int n = 2 + static_cast<int>(rng() % 6);
    const int m = 1 + static_cast<int>(rng() % 14);
    std::vector<InEdge> edges;
    for (int i = 0; i < m; ++i) {
      edges.push_back({static_cast<NodeId>(rng() % n), static_cast<NodeId>(rng() % n),
                       static_cast<Amount>(rng() % 5), static_cast<Cost>(1 + rng() % 3)});
    }
    auto g = make_graph(n, edges);
    std::vector<Amount> residual(static_cast<std::size_t>(g.num_edges()));
    for (std::size_t e = 0; e < residual.size(); ++e) {
      const Amount c = g.capacity_units()[e];
      residual[e] = c == 0 ? 0 : static_cast<Amount>(rng() % static_cast<std::uint32_t>(c + 1));
    }
    const bool use_res = (rng() % 2) == 0;
    const auto tie = (rng() % 2) == 0 ? EdgeTieBreak::Deterministic : EdgeTieBreak::PreferHigherPriority;
    const EdgePriority prios[3] = {EdgePriority::Residual, EdgePriority::Capacity, EdgePriority::Headroom};
    const auto prio = prios[rng() % 3];
    const bool multi_edge = (rng() % 2) == 0;
    const NodeId s = static_cast<NodeId>(rng() % n);
    NodeId t = static_cast<NodeId>(rng() % n);
    if (t == s) t = static_cast<NodeId>((s + 1) % n);
    const auto selection = sel(false, multi_edge, tie, prio);
    std::optional<std::vector<Amount>> res;
    if (use_res) res = residual;
    auto got = run(g, {s}, Direction::Forward, std::vector<NodeId> {t}, selection, res);
    const auto member = detail::member_mask(g, nullptr, selection, use_res ? std::span<const Amount>(residual)
                                                                           : std::span<const Amount> {}, 1);
    std::vector<std::vector<EdgeId>> walks;
    std::vector<EdgeId> walk;
    std::vector<std::uint8_t> on(static_cast<std::size_t>(n), 0);
    on[static_cast<std::size_t>(s)] = 1;
    enumerate(g, member, s, {t}, walk, on, walks);
    if (walks.empty()) {
      EXPECT_FALSE(got.path.has_value());
      continue;
    }
    const detail::PriorityValues pv {&g, prio, use_res ? std::span<const Amount>(residual) : std::span<const Amount> {}};
    auto key_less = [&](const std::vector<EdgeId>& a, const std::vector<EdgeId>& b) {
      Cost ca = 0, cb = 0;
      for (EdgeId e : a) ca += g.cost()[static_cast<std::size_t>(e)];
      for (EdgeId e : b) cb += g.cost()[static_cast<std::size_t>(e)];
      if (ca != cb) return ca < cb;
      if (tie == EdgeTieBreak::PreferHigherPriority) {
        detail::Priority ba {0, 1, true}, bb {0, 1, true};
        for (EdgeId e : a) if (detail::compare(pv.of(e), ba) < 0) ba = pv.of(e);
        for (EdgeId e : b) if (detail::compare(pv.of(e), bb) < 0) bb = pv.of(e);
        const int c = detail::compare(ba, bb);
        if (c != 0) return c > 0;
      }
      if (a.size() != b.size()) return a.size() < b.size();
      return a < b;
    };
    const auto best = *std::min_element(walks.begin(), walks.end(), key_less);
    ASSERT_TRUE(got.path.has_value());
    EXPECT_EQ(got.path->walk, best) << "trial " << trial;
    ++checked;
  }
  EXPECT_GT(checked, 600);
}

TEST(Search, MultiRootSinglePathStartsAtTheNearestRootByHops) {
  // Roots {0, 1}, both at cost 0; 0 -> 3 in one hop of cost 2 and 1 -> 2 -> 3 in two hops of
  // cost 1 each: the walk of fewest hops (from root 0).
  auto g = make_graph(4, {{0, 3, 1, 2}, {1, 2, 1, 1}, {2, 3, 1, 1}});
  const auto inv = id_of_input(g);
  auto r = run(g, {0, 1}, Direction::Forward, std::vector<NodeId> {3}, sel(false));
  ASSERT_TRUE(r.path.has_value());
  EXPECT_EQ(r.path->walk, (std::vector<EdgeId>{inv[0]}));
}

TEST(Search, StagedPathReview13CanonicalWitnesses) {
  // review-13 contract_probes.py: the C2 diamond (0-1-3 and 0-2-3 both cost 10) takes 0-2-3,
  // whose edge-id sequence is smaller; the equal three-hop graph takes 0-1-4-5. Both tie modes
  // (all capacities 1: priority ties), both directions.
  auto c2 = make_graph(4, {{0, 1, 1, 9}, {0, 2, 1, 1}, {1, 3, 1, 1}, {2, 3, 1, 9}});
  auto three = make_graph(6, {{0, 2, 1, 1}, {2, 3, 1, 1}, {3, 5, 1, 1}, {0, 1, 1, 1}, {1, 4, 1, 1}, {4, 5, 1, 1}});
  auto nodes_of = [](const StrictMultiDiGraph& g, const std::vector<EdgeId>& walk) {
    std::vector<NodeId> out = {g.edge_src()[static_cast<std::size_t>(walk.front())]};
    for (EdgeId e : walk) out.push_back(g.edge_dst()[static_cast<std::size_t>(e)]);
    return out;
  };
  for (auto tie : {EdgeTieBreak::Deterministic, EdgeTieBreak::PreferHigherPriority}) {
    for (auto prio : {EdgePriority::Residual, EdgePriority::Capacity, EdgePriority::Headroom}) {
      const auto s = sel(false, true, tie, prio);
      auto a = run(c2, {0}, Direction::Forward, std::vector<NodeId> {3}, s);
      EXPECT_EQ(nodes_of(c2, a.path->walk), (std::vector<NodeId>{0, 2, 3}));
      auto ar = run(c2, {3}, Direction::Reverse, std::vector<NodeId> {0}, s);
      EXPECT_EQ(ar.path->walk, a.path->walk);
      auto b = run(three, {0}, Direction::Forward, std::vector<NodeId> {5}, s);
      EXPECT_EQ(nodes_of(three, b.path->walk), (std::vector<NodeId>{0, 1, 4, 5}));
      auto br = run(three, {5}, Direction::Reverse, std::vector<NodeId> {0}, s);
      EXPECT_EQ(br.path->walk, b.path->walk);
    }
  }
}
