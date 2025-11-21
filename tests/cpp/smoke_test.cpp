#include <gtest/gtest.h>
#include <span>
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/types.hpp"

using namespace netgraph::core;

TEST(GraphSmoke, ConstructAndRunSPF) {
  const std::int32_t N = 3;
  // 0 -> 1 (cost 10)
  // 1 -> 2 (cost 20)
  // 0 -> 2 (cost 100)
  std::int32_t src[3] = {0, 1, 0};
  std::int32_t dst[3] = {1, 2, 2};
  double cap[3] = {10.0, 10.0, 10.0};
  std::int64_t cost[3] = {10, 20, 100};

  auto g = StrictMultiDiGraph::from_arrays(N,
      std::span<const std::int32_t>(src, 3),
      std::span<const std::int32_t>(dst, 3),
      std::span<const double>(cap, 3),
      std::span<const std::int64_t>(cost, 3));

  EXPECT_EQ(g.num_nodes(), N);
  EXPECT_EQ(g.num_edges(), 3);

  // Run SPF from 0
  EdgeSelection sel;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel);

  // Check distances
  EXPECT_EQ(dist[0], 0);
  EXPECT_EQ(dist[1], 10);
  EXPECT_EQ(dist[2], 30); // Path 0->1->2 is cost 10+20=30, better than 0->2 (100)
}
