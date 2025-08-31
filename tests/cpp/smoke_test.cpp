#include <gtest/gtest.h>
#include <span>
#include "netgraph/core/strict_multidigraph.hpp"

using namespace netgraph::core;

TEST(GraphSmoke, ConstructFromArrays) {
  const std::int32_t N = 3;
  std::int32_t src[2] = {0, 1};
  std::int32_t dst[2] = {1, 2};
  double cap[2] = {1.0, 2.0};
  double cost[2] = {0.5, 1.5};
  auto g = StrictMultiDiGraph::from_arrays(N,
      std::span<const std::int32_t>(src, 2),
      std::span<const std::int32_t>(dst, 2),
      std::span<const double>(cap, 2),
      std::span<const double>(cost, 2),
      std::span<const std::int64_t>(), false);
  EXPECT_EQ(g.num_nodes(), N);
  EXPECT_EQ(g.num_edges(), 2);
}
