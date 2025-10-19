#include <gtest/gtest.h>
#include "netgraph/core/k_shortest_paths.hpp"
#include "netgraph/core/backend.hpp"
#include "netgraph/core/algorithms.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

TEST(KShortestPaths, FindsKPaths) {
  auto g = make_square_graph(1);  // Has at least 2 paths

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  KspOptions opts;
  opts.k = 2;
  opts.unique = true;
  opts.max_cost_factor = std::nullopt;

  auto paths = algs.ksp(gh, 0, 2, opts);

  // Should find at least 1 path, up to 2
  EXPECT_GT(paths.size(), 0);
  EXPECT_LE(paths.size(), 2);
}

TEST(KShortestPaths, KEqualsZeroReturnsEmpty) {
  auto g = make_square_graph(1);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);
  KspOptions opts; opts.k = 0; opts.unique = true;
  auto paths = algs.ksp(gh, 0, 2, opts);
  EXPECT_TRUE(paths.empty());
}

TEST(KShortestPaths, NodeMaskBlocksPaths) {
  // Verify that node masks prevent path discovery in KSP
  auto g = make_line_graph(3);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);
  KspOptions opts; opts.k = 5; opts.unique = true;
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[1] = false;  // Block middle node
  opts.node_mask = std::span<const bool>(node_mask.get(), static_cast<std::size_t>(g.num_nodes()));
  auto paths = algs.ksp(gh, 0, 2, opts);
  EXPECT_TRUE(paths.empty());
}

TEST(KShortestPaths, UniqueModeFiltering) {
  auto g = make_square_graph(1);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  KspOptions opts1;
  opts1.k = 10;
  opts1.unique = true;
  auto paths_unique = algs.ksp(gh, 0, 2, opts1);

  KspOptions opts2;
  opts2.k = 10;
  opts2.unique = false;
  auto paths_non_unique = algs.ksp(gh, 0, 2, opts2);

  // Non-unique mode should find at least as many paths as unique mode
  EXPECT_GE(paths_non_unique.size(), paths_unique.size());
}

TEST(KShortestPaths, MaxCostFactorLimit) {
  auto g = make_square_graph(1);  // Has paths of cost 2 and 4

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  KspOptions opts1;
  opts1.k = 10;
  opts1.max_cost_factor = 1.5;  // Allow up to 1.5x shortest
  auto paths_limited = algs.ksp(gh, 0, 2, opts1);

  KspOptions opts2;
  opts2.k = 10;
  opts2.max_cost_factor = std::nullopt;  // No limit
  auto paths_unlimited = algs.ksp(gh, 0, 2, opts2);

  // Unlimited should find at least as many paths
  EXPECT_GE(paths_unlimited.size(), paths_limited.size());

  // All limited paths should have reasonable cost (check via distance to destination)
  if (!paths_limited.empty()) {
    Cost min_cost = paths_limited[0].first[2];  // Cost to destination node 2

    for (const auto& [dist, dag] : paths_limited) {
      Cost path_cost = dist[2];  // Cost to destination
      EXPECT_LE(path_cost, min_cost * 1.5 + 1e-9);
    }
  }
}

TEST(KShortestPaths, PathsSortedByCost) {
  auto g = make_square_graph(1);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  KspOptions opts;
  opts.k = 5;
  opts.unique = true;

  auto paths = algs.ksp(gh, 0, 2, opts);

  if (paths.size() > 1) {
    // Extract costs to destination (node 2)
    std::vector<Cost> costs;
    for (const auto& [dist, dag] : paths) {
      costs.push_back(dist[2]);  // Cost to destination node 2
    }

    // Check sorted
    for (std::size_t i = 0; i < costs.size() - 1; ++i) {
      EXPECT_LE(costs[i], costs[i + 1] + 1e-9) << "Paths not sorted by cost";
    }
  }
}

TEST(KShortestPaths, DisconnectedReturnsEmpty) {
  // Create disconnected graph: 0-1 and 2-3
  std::int32_t src[2] = {0, 2};
  std::int32_t dst[2] = {1, 3};
  double cap[2] = {1.0, 1.0};
  std::int64_t cost[2] = {1, 1};

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src, 2), std::span(dst, 2),
    std::span(cap, 2), std::span(cost, 2));

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  KspOptions opts;
  opts.k = 5;

  // Try to find path from 0 to 3 (disconnected)
  auto paths = algs.ksp(gh, 0, 3, opts);

  EXPECT_EQ(paths.size(), 0);
}

TEST(KShortestPaths, LooplessPaths) {
  // Create a graph with cycles: 0->1->2->3 (forward path) and 1->0, 2->0 (back edges)
  std::int32_t src[6] = {0, 1, 2, 1, 0, 2};
  std::int32_t dst[6] = {1, 2, 3, 0, 2, 0};  // Has cycles
  double cap[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  std::int64_t cost[6] = {1, 1, 1, 1, 2, 2};

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src, 6), std::span(dst, 6),
    std::span(cap, 6), std::span(cost, 6));

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  KspOptions opts;
  opts.k = 10;
  opts.unique = true;

  auto paths = algs.ksp(gh, 0, 3, opts);

  // KSP should return loop-free paths (Yen's algorithm property)
  // Verify by checking that destination is reachable and has finite cost
  for (const auto& [dist, dag] : paths) {
    // Destination (node 3) should be reachable with finite cost
    EXPECT_LT(dist[3], std::numeric_limits<Cost>::max())
      << "Path to destination should have finite cost (no infinite loops)";

    // Source (node 0) should have distance 0
    EXPECT_DOUBLE_EQ(dist[0], 0.0) << "Source distance should be 0";

    // Path cost should be positive (since all edge costs are positive)
    EXPECT_GT(dist[3], 0.0) << "Path cost should be positive";
  }

  // Should find at least one path (0->1->2->3 with cost 3)
  EXPECT_GT(paths.size(), 0) << "Should find at least one loop-free path";

  // First path should be shortest (cost 3)
  if (!paths.empty()) {
    EXPECT_DOUBLE_EQ(paths[0].first[3], 3.0) << "Shortest path should have cost 3";
  }
}
