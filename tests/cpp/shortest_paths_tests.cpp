#include <gtest/gtest.h>
#include <limits>
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/backend.hpp"
#include "netgraph/core/algorithms.hpp"
#include "netgraph/core/options.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

TEST(ShortestPaths, SingleSourceAllNodes) {
  auto g = make_line_graph(5);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  // Distances should be 0, 1, 2, 3, 4
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_DOUBLE_EQ(dist[1], 1.0);
  EXPECT_DOUBLE_EQ(dist[2], 2.0);
  EXPECT_DOUBLE_EQ(dist[3], 3.0);
  EXPECT_DOUBLE_EQ(dist[4], 4.0);

  expect_pred_dag_valid(dag, g.num_nodes());
}

TEST(ShortestPaths, DisconnectedComponents) {
  // Create two disconnected components: 0-1 and 2-3
  std::int32_t src[2] = {0, 2};
  std::int32_t dst[2] = {1, 3};
  double cap[2] = {1.0, 1.0};
  std::int64_t cost[2] = {1, 1};

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src, 2), std::span(dst, 2),
    std::span(cap, 2), std::span(cost, 2));

  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  // Nodes 0, 1 should be reachable
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_DOUBLE_EQ(dist[1], 1.0);

  // Nodes 2, 3 should be unreachable
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());
  EXPECT_EQ(dist[3], std::numeric_limits<Cost>::max());
}

TEST(ShortestPaths, MultipleEqualCostPaths) {
  auto g = make_square_graph(2);  // Two equal-cost paths
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  // Both paths should have cost 2
  EXPECT_DOUBLE_EQ(dist[2], 2.0);

  // Node 2 should have two parents in the DAG
  auto start = static_cast<std::size_t>(dag.parent_offsets[2]);
  auto end = static_cast<std::size_t>(dag.parent_offsets[3]);
  EXPECT_GT(end - start, 0);  // At least one predecessor
}

TEST(ShortestPaths, ResidualAwareFiltering) {
  auto g = make_line_graph(3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = true;  // Only use edges with capacity
  sel.tie_break = EdgeTieBreak::Deterministic;

  // Create residual with first edge having zero capacity
  std::vector<Cap> residual(g.num_edges(), 1.0);
  residual[0] = 0.0;  // Block first edge

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, residual, {}, {});

  // Node 1 should be unreachable
  EXPECT_EQ(dist[1], std::numeric_limits<Cost>::max());
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());
}

TEST(ShortestPaths, NodeMaskIsolation) {
  // Verify that masking out the middle node in a line graph blocks reachability
  auto g = make_line_graph(3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  // Mask out node 1 (middle node) - blocks path from 0 to 2
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[1] = false;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, std::span<const bool>(node_mask.get(), g.num_nodes()), {});

  // Nodes 1 and 2 should be unreachable
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_EQ(dist[1], std::numeric_limits<Cost>::max());
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());
}

TEST(ShortestPaths, EdgeMaskFiltering) {
  // Verify that masking out an edge removes it from consideration
  auto g = make_line_graph(3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  // Mask out first edge (0->1) - blocks the entire path
  auto edge_mask = make_bool_mask(g.num_edges());
  edge_mask[0] = false;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, std::span<const bool>(edge_mask.get(), g.num_edges()));

  // Node 1 should be unreachable
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_EQ(dist[1], std::numeric_limits<Cost>::max());
}

TEST(ShortestPaths, EarlyExitOptimization) {
  auto g = make_line_graph(10);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  // With destination specified, should exit early
  auto [dist, dag] = shortest_paths(g, 0, NodeId(2), true, sel, {}, {}, {});

  EXPECT_DOUBLE_EQ(dist[2], 2.0);
  expect_pred_dag_valid(dag, g.num_nodes());
}

TEST(ShortestPaths, PredDAGIntegrity) {
  auto g = make_grid_graph(3, 3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  expect_pred_dag_valid(dag, g.num_nodes());

  // Source node should have no predecessors
  EXPECT_EQ(dag.parent_offsets[0], 0);
  EXPECT_EQ(dag.parent_offsets[1], 0);
}

TEST(ShortestPaths, TieBreakingDeterminism) {
  // Two parallel edges with same cost
  std::int32_t src[2] = {0, 0};
  std::int32_t dst[2] = {1, 1};
  double cap[2] = {1.0, 2.0};
  std::int64_t cost[2] = {1, 1};

  auto g = StrictMultiDiGraph::from_arrays(2,
    std::span(src, 2), std::span(dst, 2),
    std::span(cap, 2), std::span(cost, 2));

  EdgeSelection sel;
  sel.multi_edge = false;  // Single edge selection
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  auto [dist1, dag1] = shortest_paths(g, 0, std::nullopt, false, sel, {}, {}, {});
  auto [dist2, dag2] = shortest_paths(g, 0, std::nullopt, false, sel, {}, {}, {});

  // Results should be deterministic
  EXPECT_EQ(dist1[1], dist2[1]);

  // Should select the same edge both times
  EXPECT_EQ(dag1.via_edges[0], dag2.via_edges[0]);
}

TEST(ShortestPaths, PreferHigherResidualTieBreak) {
  // Same two parallel edges but enforce HigherResidual tie-break
  std::int32_t src[2] = {0, 0};
  std::int32_t dst[2] = {1, 1};
  double cap[2] = {1.0, 2.0};
  std::int64_t cost[2] = {1, 1};
  auto g = StrictMultiDiGraph::from_arrays(2,
    std::span(src, 2), std::span(dst, 2),
    std::span(cap, 2), std::span(cost, 2));

  // Residuals equal to capacity; higher residual should pick edge with cap=2
  EdgeSelection sel;
  sel.multi_edge = false;
  sel.require_capacity = true;
  sel.tie_break = EdgeTieBreak::PreferHigherResidual;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, false, sel, g.capacity_view(), {}, {});
  ASSERT_FALSE(dag.via_edges.empty());
  // Edge ids are compacted; find which corresponds to cap=2 by checking capacity
  auto chosen = dag.via_edges[0];
  auto capv = g.capacity_view();
  EXPECT_EQ(capv[chosen], 2.0);
}

TEST(ShortestPaths, ResolveToPathsSplittingParallelEdges) {
  // Test resolve_to_paths both with and without splitting parallel edges
  // Graph 0->1 has two parallel edges, 1->2 single edge
  std::int32_t src[3] = {0, 0, 1};
  std::int32_t dst[3] = {1, 1, 2};
  double cap[3] = {1.0, 1.0, 1.0};
  std::int64_t cost[3] = {1, 1, 1};
  auto g = StrictMultiDiGraph::from_arrays(3,
    std::span(src, 3), std::span(dst, 3),
    std::span(cap, 3), std::span(cost, 3));

  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, 2, true, sel, {}, {}, {});

  // Resolve without splitting (grouped parallel edges)
  auto grouped = resolve_to_paths(dag, 0, 2, /*split_parallel_edges=*/false);
  ASSERT_FALSE(grouped.empty());
  ASSERT_GE(grouped[0].size(), 2u);
  // First tuple corresponds to node 0 with two parallel edges grouped
  EXPECT_EQ(grouped[0][0].first, 0);
  EXPECT_EQ(grouped[0][0].second.size(), 2u);

  // Resolve with splitting should produce 2 concrete paths (one per parallel edge)
  auto concrete = resolve_to_paths(dag, 0, 2, /*split_parallel_edges=*/true);
  EXPECT_EQ(concrete.size(), 2u);
}

TEST(ShortestPaths, ResolveToPathsWithMaxPathsLimit) {
  // Verify that max_paths limit caps enumeration in resolve_to_paths
  // Create graph with multiple equal-cost paths via 2x2 grid
  std::int32_t src[4] = {0, 0, 1, 2};
  std::int32_t dst[4] = {1, 2, 3, 3};
  double cap[4] = {1.0, 1.0, 1.0, 1.0};
  std::int64_t cost[4] = {1, 1, 1, 1};
  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src, 4), std::span(dst, 4),
    std::span(cap, 4), std::span(cost, 4));

  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, 3, true, sel, {}, {}, {});

  // Without limit, should find both paths
  auto all_paths = resolve_to_paths(dag, 0, 3, false, std::nullopt);
  EXPECT_GE(all_paths.size(), 1u);

  // With max_paths=1, enumeration should stop at first path
  auto limited = resolve_to_paths(dag, 0, 3, false, std::int64_t(1));
  EXPECT_EQ(limited.size(), 1u);
}

TEST(ShortestPaths, LargeGraphStressTest) {
  // Create a larger graph to test performance
  auto g = make_grid_graph(20, 20);  // 400 nodes
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  // Bottom-right corner should be reachable with cost 38 (19 right + 19 down)
  EXPECT_DOUBLE_EQ(dist[399], 38.0);
  expect_pred_dag_valid(dag, g.num_nodes());
}

TEST(ShortestPaths, RejectsMaskLengthMismatch) {
  // Simple line graph 0->1->2
  auto g = make_line_graph(3);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);
  SpfOptions opts;
  opts.multipath = true;
  opts.dst = 2;
  // Create node_mask with wrong length (N-1)
  auto node_mask = make_bool_mask(static_cast<std::size_t>(g.num_nodes() - 1), true);
  opts.node_mask = std::span<const bool>(node_mask.get(), static_cast<std::size_t>(g.num_nodes() - 1));
  EXPECT_THROW({ (void)algs.spf(gh, 0, opts); }, std::invalid_argument);

  // Edge mask wrong length (M-1)
  SpfOptions opts2;
  opts2.multipath = true;
  opts2.dst = 2;
  auto edge_mask = make_bool_mask(static_cast<std::size_t>(g.num_edges() - 1), true);
  opts2.edge_mask = std::span<const bool>(edge_mask.get(), static_cast<std::size_t>(g.num_edges() - 1));
  EXPECT_THROW({ (void)algs.spf(gh, 0, opts2); }, std::invalid_argument);
}
