// Eligibility as the graph builds it: dense node and edge flags from excluded ids,
// graph-bound, uid 0 for the all-allowed value.
#include <gtest/gtest.h>

#include <string>
#include <type_traits>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

std::vector<int> flags(std::span<const std::uint8_t> s) { return std::vector<int>(s.begin(), s.end()); }

} // namespace

TEST(Eligibility, AllAllowedHasUidZero) {
  // Built only by Graph::eligibility and intersect: no default construction.
  static_assert(!std::is_default_constructible_v<Eligibility>);
  auto g = make_graph(3, {{0, 1, 1, 1}, {1, 2, 1, 1}});
  auto el = g.eligibility();
  EXPECT_EQ(el.uid(), 0u);
  EXPECT_EQ(el.graph_uid(), g.uid());
  EXPECT_EQ(flags(el.node_allowed()), (std::vector<int>{1, 1, 1}));
  EXPECT_EQ(flags(el.edge_allowed()), (std::vector<int>{1, 1}));
}

TEST(Eligibility, ExcludedNodeExcludesIncidentEdges) {
  // In-edges, out-edges and a self-loop of the excluded node.
  auto g = make_graph(4, {{0, 1, 1, 1}, {1, 2, 1, 1}, {2, 3, 1, 1}, {1, 1, 1, 2}, {3, 0, 1, 1}});
  const auto inv = id_of_input(g);
  const std::vector<NodeId> xn = {1};
  auto el = g.eligibility(xn, {});
  EXPECT_NE(el.uid(), 0u);
  EXPECT_EQ(flags(el.node_allowed()), (std::vector<int>{1, 0, 1, 1}));
  auto allowed_in = to_input(g, std::vector<std::uint8_t>(el.edge_allowed().begin(), el.edge_allowed().end()));
  EXPECT_EQ(std::vector<int>(allowed_in.begin(), allowed_in.end()), (std::vector<int>{0, 0, 1, 0, 1}));
  EXPECT_FALSE(el.allows_edge(inv[0]));
  EXPECT_TRUE(el.allows_edge(inv[2]));
}

TEST(Eligibility, ExcludedEdgesDuplicatesAndOrder) {
  // EdgeIds 0: 0 -> 1, 1: 0 -> 2, 2: 1 -> 2, 3: 3 -> 0 (the one edge left allowed).
  auto g = make_graph(4, {{0, 1, 1, 1}, {1, 2, 1, 1}, {0, 2, 1, 5}, {3, 0, 1, 1}});
  const std::vector<EdgeId> xe = {2, 0, 2};
  const std::vector<NodeId> xn = {2, 2};
  auto el = g.eligibility(xn, xe);
  EXPECT_EQ(flags(el.node_allowed()), (std::vector<int>{1, 1, 0, 1}));
  int allowed = 0;
  for (EdgeId e = 0; e < 4; ++e) {
    const bool excluded = e == 0 || e == 2 || g.edge_destinations()[static_cast<std::size_t>(e)] == 2;
    EXPECT_EQ(el.allows_edge(e), !excluded) << e;
    allowed += excluded ? 0 : 1;
  }
  EXPECT_EQ(allowed, 1);  // both answers are expected
}

TEST(Eligibility, DistinctUidsForEqualLists) {
  auto g = make_graph(2, {{0, 1, 1, 1}});
  const std::vector<EdgeId> xe = {0};
  auto a = g.eligibility({}, xe);
  auto b = g.eligibility({}, xe);
  EXPECT_NE(a.uid(), b.uid());
  EXPECT_NE(a.uid(), 0u);
  EXPECT_EQ(flags(a.edge_allowed()), flags(b.edge_allowed()));
}

TEST(Eligibility, OutOfRangeIdsRejected) {
  auto g = make_graph(2, {{0, 1, 1, 1}});
  const std::vector<NodeId> bad_node = {2};
  const std::vector<EdgeId> bad_edge = {1};
  const std::vector<NodeId> neg_node = {-1};
  EXPECT_EQ(error_of([&] { (void)g.eligibility(bad_node, {}); }), "excluded node 2 out of node range [0, 2)");
  EXPECT_EQ(error_of([&] { (void)g.eligibility({}, bad_edge); }), "excluded edge 1 out of edge range [0, 1)");
  EXPECT_EQ(error_of([&] { (void)g.eligibility(neg_node, {}); }), "excluded node -1 out of node range [0, 2)");
}

// Eligibility.intersect: an operand when the result equals it (its uid kept, a
// first), otherwise a new value with a fresh uid; values of different graphs rejected.
TEST(Eligibility, IntersectReturnsAnOperandWhenEqual) {
  auto g = make_graph(3, {{0, 1, 1, 1}, {1, 2, 1, 1}, {0, 2, 1, 5}});
  const auto all = g.eligibility();
  const std::vector<EdgeId> e0 = {0};
  const std::vector<EdgeId> e1 = {1};
  const std::vector<EdgeId> e01 = {0, 1};
  const auto a = g.eligibility({}, e0);
  const auto b = g.eligibility({}, e1);
  // all-allowed with a: a (its uid); a with all-allowed: a.
  EXPECT_EQ(Eligibility::intersect(all, a).uid(), a.uid());
  EXPECT_EQ(Eligibility::intersect(a, all).uid(), a.uid());
  EXPECT_EQ(Eligibility::intersect(all, all).uid(), 0u);
  EXPECT_EQ(Eligibility::intersect(a, a).uid(), a.uid());
  {
    // Two equal values built separately (distinct uids): the first operand wins the tie.
    const auto x = g.eligibility({}, e0);
    const auto y = g.eligibility({}, e0);
    ASSERT_NE(x.uid(), y.uid());
    EXPECT_EQ(Eligibility::intersect(x, y).uid(), x.uid());
    EXPECT_EQ(Eligibility::intersect(y, x).uid(), y.uid());
  }
  // a subset relation: the smaller value.
  const auto ab_list = g.eligibility({}, e01);
  EXPECT_EQ(Eligibility::intersect(a, ab_list).uid(), ab_list.uid());
  // A new value: a fresh uid, flags are the conjunction; equal results built twice differ.
  const auto x = Eligibility::intersect(a, b);
  const auto y = Eligibility::intersect(a, b);
  EXPECT_NE(x.uid(), a.uid());
  EXPECT_NE(x.uid(), b.uid());
  EXPECT_NE(x.uid(), 0u);
  EXPECT_NE(x.uid(), y.uid());
  EXPECT_EQ(flags(x.edge_allowed()), flags(ab_list.edge_allowed()));
  EXPECT_EQ(x.graph_uid(), g.uid());
  // Node exclusions intersect too (and keep excluding incident edges).
  const std::vector<NodeId> n2 = {2};
  const auto c = g.eligibility(n2, {});
  const auto z = Eligibility::intersect(a, c);
  EXPECT_EQ(flags(z.node_allowed()), (std::vector<int>{1, 1, 0}));
  EXPECT_EQ(flags(z.edge_allowed()), (std::vector<int>{0, 0, 0}));
}

TEST(Eligibility, IntersectRejectsValuesOfDifferentGraphs) {
  auto g = make_graph(2, {{0, 1, 1, 1}});
  auto h = make_graph(2, {{0, 1, 1, 1}});
  EXPECT_EQ(error_of([&] { (void)Eligibility::intersect(g.eligibility(), h.eligibility()); }),
            "Eligibility.intersect: the two values belong to different graphs");
}
