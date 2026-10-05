#include <gtest/gtest.h>

#include <optional>
#include <span>
#include <vector>

#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

using namespace netgraph::core;

TEST(GraphSmoke, ConstructAndSolve) {
  // 0 -> 1 (cost 10), 1 -> 2 (cost 20), 0 -> 2 (cost 100); capacities 10 units each.
  const std::vector<NodeId> src = {0, 1, 0};
  const std::vector<NodeId> dst = {1, 2, 2};
  const std::vector<Amount> cap = {10, 10, 10};
  const std::vector<Cost> cost = {10, 20, 100};
  StrictMultiDiGraph g(3, src, dst, cap, cost, 1.0);
  EXPECT_EQ(g.num_nodes(), 3);
  EXPECT_EQ(g.num_edges(), 3);

  const std::vector<NodeId> s = {0};
  const std::vector<NodeId> t = {2};
  auto r = place_global_units(g, g.capacity_units(), s, t, nullptr, std::nullopt);
  EXPECT_EQ(r.placed, 20);
  ASSERT_EQ(r.cost_distribution.size(), 2u);
  EXPECT_EQ(r.cost_distribution[0], (std::pair<Cost, Amount>{30, 10}));   // 0 -> 1 -> 2 first
  EXPECT_EQ(r.cost_distribution[1], (std::pair<Cost, Amount>{100, 10}));  // then the direct edge
}
