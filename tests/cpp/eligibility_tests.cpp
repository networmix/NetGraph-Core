// Eligibility (APPROACH 4) as the graph builds it: dense node and edge flags from excluded ids,
// graph-bound, uid 0 for the all-allowed value.
#include <gtest/gtest.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

std::vector<int> flags(std::span<const std::uint8_t> s) { return std::vector<int>(s.begin(), s.end()); }

} // namespace

TEST(Eligibility, AllAllowedHasUidZero) {
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
  EXPECT_FALSE(el.edge_ok(inv[0]));
  EXPECT_TRUE(el.edge_ok(inv[2]));
}

TEST(Eligibility, ExcludedEdgesDuplicatesAndOrder) {
  auto g = make_graph(3, {{0, 1, 1, 1}, {1, 2, 1, 1}, {0, 2, 1, 5}});
  const std::vector<EdgeId> xe = {2, 0, 2};
  const std::vector<NodeId> xn = {2, 2};
  auto el = g.eligibility(xn, xe);
  EXPECT_EQ(flags(el.node_allowed()), (std::vector<int>{1, 1, 0}));
  for (EdgeId e = 0; e < 3; ++e) {
    const bool excluded = e == 0 || e == 2 || g.edge_dst()[static_cast<std::size_t>(e)] == 2;
    EXPECT_EQ(el.edge_ok(e), !excluded) << e;
  }
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
  EXPECT_THROW((void)g.eligibility(bad_node, {}), std::invalid_argument);
  EXPECT_THROW((void)g.eligibility({}, bad_edge), std::invalid_argument);
  EXPECT_THROW((void)g.eligibility(neg_node, {}), std::invalid_argument);
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
  EXPECT_THROW((void)Eligibility::intersect(g.eligibility(), h.eligibility()), std::invalid_argument);
}
