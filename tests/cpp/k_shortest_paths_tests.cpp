#include <gtest/gtest.h>
#include <limits>
#include <set>
#include "netgraph/core/k_shortest_paths.hpp"
#include "netgraph/core/shortest_paths.hpp"
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

TEST(KShortestPaths, OutOfOrderNodePredDAGSemantics) {
  // Graph topology that forces paths visiting nodes out of numerical order.
  // Path 1 (shortest): 0 -> 1 -> 2 (cost 2, in-order)
  // Path 2 (longer):   0 -> 3 -> 2 (cost 4, visits node 3 before node 2)
  // This was the exact bug scenario: PredDAG fill for path [0,3,2] swapped
  // parent entries between nodes 2 and 3 when using a linear index.
  auto g = make_square_graph(1);  // 0->1->2 (cost 2) vs 0->3->2 (cost 4)

  auto items = k_shortest_paths(g, 0, 2, 5, std::nullopt, true);

  ASSERT_GE(items.size(), 2) << "Should find at least 2 paths";

  // Verify all PredDAGs have correct semantic structure
  for (std::size_t p = 0; p < items.size(); ++p) {
    const auto& [dist, dag] = items[p];
    SCOPED_TRACE("Path " + std::to_string(p));

    expect_pred_dag_valid(dag, g.num_nodes());
    expect_pred_dag_semantically_valid(g, dag, dist);

    // Source must have distance 0 and destination must be reachable
    EXPECT_EQ(dist[0], 0);
    EXPECT_LT(dist[2], std::numeric_limits<Cost>::max());

    // Verify path can be reconstructed via resolve_to_paths
    auto concrete = resolve_to_paths(dag, 0, 2, /*split_parallel_edges=*/true);
    EXPECT_GE(concrete.size(), 1) << "resolve_to_paths should reconstruct at least one path";

    // Verify reconstructed path is valid: consecutive nodes connected by claimed edges
    for (const auto& path : concrete) {
      ASSERT_GE(path.size(), 2);  // at least src and dst
      EXPECT_EQ(path.front().first, 0) << "Path should start at source";
      EXPECT_EQ(path.back().first, 2) << "Path should end at destination";
      // Check edge connectivity
      for (std::size_t i = 0; i + 1 < path.size(); ++i) {
        auto from_node = path[i].first;
        auto to_node = path[i + 1].first;
        const auto& edges = path[i].second;
        for (auto eid : edges) {
          EXPECT_EQ(g.edge_src_view()[static_cast<std::size_t>(eid)], from_node);
          EXPECT_EQ(g.edge_dst_view()[static_cast<std::size_t>(eid)], to_node);
        }
      }
    }
  }
}

TEST(KShortestPaths, LargerOutOfOrderTopology) {
  // Larger graph where multiple k-shortest paths visit nodes out of numerical order.
  // Topology (6 nodes):
  //   0 -> 5 -> 3 -> 1 (cost 3, severely out-of-order: 0,5,3,1)
  //   0 -> 2 -> 4 -> 1 (cost 6, also out-of-order: 0,2,4,1)
  //   0 -> 1           (cost 10, direct)
  std::int32_t src_arr[7] = {0, 5, 3, 0, 2, 4, 0};
  std::int32_t dst_arr[7] = {5, 3, 1, 2, 4, 1, 1};
  double cap_arr[7] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  std::int64_t cost_arr[7] = {1, 1, 1, 2, 2, 2, 10};

  auto g = StrictMultiDiGraph::from_arrays(6,
    std::span(src_arr, 7), std::span(dst_arr, 7),
    std::span(cap_arr, 7), std::span(cost_arr, 7));

  auto items = k_shortest_paths(g, 0, 1, 5, std::nullopt, true);

  ASSERT_GE(items.size(), 2) << "Should find at least 2 paths";

  // Verify costs are non-decreasing
  for (std::size_t i = 0; i + 1 < items.size(); ++i) {
    EXPECT_LE(items[i].first[1], items[i + 1].first[1])
        << "Paths should be sorted by cost to destination";
  }

  // The first path should be 0->5->3->1 (cost 3)
  EXPECT_EQ(items[0].first[1], 3) << "Shortest path cost should be 3";

  // Verify all PredDAGs are semantically correct and can be resolved
  for (std::size_t p = 0; p < items.size(); ++p) {
    const auto& [dist, dag] = items[p];
    SCOPED_TRACE("Path " + std::to_string(p));

    expect_pred_dag_valid(dag, g.num_nodes());
    expect_pred_dag_semantically_valid(g, dag, dist);

    // Source must have distance 0 and destination must be reachable
    EXPECT_EQ(dist[0], 0);
    EXPECT_LT(dist[1], std::numeric_limits<Cost>::max());

    auto concrete = resolve_to_paths(dag, 0, 1, true);
    EXPECT_GE(concrete.size(), 1) << "Should reconstruct at least one concrete path";

    // Verify node connectivity in reconstructed paths
    for (const auto& path : concrete) {
      ASSERT_GE(path.size(), 2);
      EXPECT_EQ(path.front().first, 0);
      EXPECT_EQ(path.back().first, 1);
      // Verify no repeated nodes (simple path)
      std::set<NodeId> visited;
      for (const auto& [node, edges] : path) {
        EXPECT_TRUE(visited.insert(node).second)
            << "Node " << node << " appears twice in path (cycle)";
      }
    }
  }
}
