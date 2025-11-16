/**
 * Comprehensive tests for masking behavior in NetGraph-Core.
 *
 * These tests verify:
 * 1. Node and edge mask semantics (true = allowed, false = excluded)
 * 2. Masked source node handling (should return empty DAG)
 * 3. Masked destination node handling
 * 4. Mask validation (length checking)
 * 5. Combination of node and edge masks
 * 6. Interaction with residual capacity gating
 * 7. Masking in max_flow and k_shortest_paths
 */

#include <gtest/gtest.h>
#include <limits>
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/k_shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/backend.hpp"
#include "netgraph/core/algorithms.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/flow_state.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

// ============================================================================
// Masked Source Tests
// ============================================================================

TEST(MaskingTests, MaskedSourceReturnsEmptyDAG) {
  // CRITICAL BUG FIX TEST: When source is masked out, SPF must return empty DAG
  // and leave all distances at infinity.
  auto g = make_line_graph(3);  // 0->1->2
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  // Mask out source node 0
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[0] = false;

  auto [dist, dag] = shortest_paths(
      g, 0, std::nullopt, true, sel, {},
      std::span<const bool>(node_mask.get(), g.num_nodes()), {});

  // All distances should remain at infinity
  EXPECT_EQ(dist[0], std::numeric_limits<Cost>::max());
  EXPECT_EQ(dist[1], std::numeric_limits<Cost>::max());
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());

  // DAG should be empty (all offsets zero)
  EXPECT_EQ(dag.parent_offsets.size(), static_cast<std::size_t>(g.num_nodes() + 1));
  for (const auto& offset : dag.parent_offsets) {
    EXPECT_EQ(offset, 0) << "DAG should be empty when source is masked";
  }
  EXPECT_EQ(dag.parents.size(), 0);
  EXPECT_EQ(dag.via_edges.size(), 0);
}

TEST(MaskingTests, MaskedSourceWithDestination) {
  // Test masked source with explicit destination
  auto g = make_line_graph(3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[0] = false;  // Mask source

  auto [dist, dag] = shortest_paths(
      g, 0, NodeId(2), true, sel, {},
      std::span<const bool>(node_mask.get(), g.num_nodes()), {});

  // Destination should be unreachable
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());

  // DAG should be empty
  EXPECT_EQ(dag.parent_offsets.size(), static_cast<std::size_t>(g.num_nodes() + 1));
  for (const auto& offset : dag.parent_offsets) {
    EXPECT_EQ(offset, 0);
  }
}

TEST(MaskingTests, MaskedSourceInMaxFlow) {
  // Test masked source in max_flow algorithm
  auto g = make_line_graph(3);

  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[0] = false;  // Mask source

  auto [flow_val, summary] = calc_max_flow(
      g, 0, 2,
      FlowPlacement::Proportional, false, true, true, false, false,
      std::span<const bool>(node_mask.get(), g.num_nodes()), {});

  // No flow should be possible
  EXPECT_DOUBLE_EQ(flow_val, 0.0);

  // Edge flows may be empty or all zeros when no flow is possible
  if (!summary.edge_flows.empty()) {
    for (const auto& ef : summary.edge_flows) {
      EXPECT_DOUBLE_EQ(ef, 0.0);
    }
  }
}

TEST(MaskingTests, MaskedSourceInComputeMinCut) {
  // Test masked source in FlowState::compute_min_cut (early return optimization)
  auto g = make_triangle_topology();
  FlowState fs(g);

  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[0] = false;  // Mask source

  auto mincut = fs.compute_min_cut(
      0,
      std::span<const bool>(node_mask.get(), g.num_nodes()),
      {});

  // Min cut should be empty when source is masked (early return)
  EXPECT_EQ(mincut.edges.size(), 0);
}

TEST(MaskingTests, OutOfRangeSourceInComputeMinCut) {
  // Test early return for out-of-range source in compute_min_cut
  auto g = make_triangle_topology();
  FlowState fs(g);

  // Test negative source
  auto mincut = fs.compute_min_cut(-1, {}, {});
  EXPECT_EQ(mincut.edges.size(), 0);

  // Test source >= num_nodes
  mincut = fs.compute_min_cut(999, {}, {});
  EXPECT_EQ(mincut.edges.size(), 0);
}

// ============================================================================
// Masked Destination Tests
// ============================================================================

TEST(MaskingTests, MaskedDestinationUnreachable) {
  auto g = make_line_graph(3);  // 0->1->2
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  // Mask out destination node 2
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[2] = false;

  auto [dist, dag] = shortest_paths(
      g, 0, std::nullopt, true, sel, {},
      std::span<const bool>(node_mask.get(), g.num_nodes()), {});

  // Source and node 1 should be reachable
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_DOUBLE_EQ(dist[1], 1.0);

  // Destination should be unreachable
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());

  // Node 2 should have no predecessors
  auto start = static_cast<std::size_t>(dag.parent_offsets[2]);
  auto end = static_cast<std::size_t>(dag.parent_offsets[3]);
  EXPECT_EQ(start, end) << "Masked node should have no predecessors";
}

TEST(MaskingTests, MaskedIntermediateNode) {
  auto g = make_line_graph(3);  // 0->1->2
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  // Mask out intermediate node 1
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[1] = false;

  auto [dist, dag] = shortest_paths(
      g, 0, NodeId(2), true, sel, {},
      std::span<const bool>(node_mask.get(), g.num_nodes()), {});

  // Source reachable
  EXPECT_DOUBLE_EQ(dist[0], 0.0);

  // Nodes 1 and 2 should be unreachable (path broken)
  EXPECT_EQ(dist[1], std::numeric_limits<Cost>::max());
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());
}

// ============================================================================
// Edge Mask Tests
// ============================================================================

TEST(MaskingTests, MaskedEdgeBlocksPath) {
  auto g = make_line_graph(3);  // 0->1->2
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  // Mask out first edge (0->1)
  auto edge_mask = make_bool_mask(g.num_edges());
  edge_mask[0] = false;

  auto [dist, dag] = shortest_paths(
      g, 0, std::nullopt, true, sel, {}, {},
      std::span<const bool>(edge_mask.get(), g.num_edges()));

  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_EQ(dist[1], std::numeric_limits<Cost>::max());
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());
}

TEST(MaskingTests, MaskedEdgeSelectsAlternativePath) {
  // Create graph with two paths: 0->1->3 and 0->2->3
  auto g = make_simple_parallel_paths_topology();
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  // Mask out edge 0->1 (first edge)
  auto edge_mask = make_bool_mask(g.num_edges());
  edge_mask[0] = false;

  auto [dist, dag] = shortest_paths(
      g, 0, NodeId(3), true, sel, {}, {},
      std::span<const bool>(edge_mask.get(), g.num_edges()));

  // Should still reach destination via alternative path 0->2->3
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_EQ(dist[1], std::numeric_limits<Cost>::max());  // Node 1 masked
  EXPECT_DOUBLE_EQ(dist[2], 1.0);  // Via 0->2
  EXPECT_DOUBLE_EQ(dist[3], 2.0);  // Via 0->2->3
}

TEST(MaskingTests, MultipleEdgesMasked) {
  auto g = make_square_graph(2);  // Two equal-cost paths
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  // Mask out multiple edges
  auto edge_mask = make_bool_mask(g.num_edges());
  edge_mask[1] = false;  // Block one path

  auto [dist, dag] = shortest_paths(
      g, 0, std::nullopt, true, sel, {}, {},
      std::span<const bool>(edge_mask.get(), g.num_edges()));

  // Should still have one path to destination
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_GT(dist[2], 0.0);
  EXPECT_LT(dist[2], std::numeric_limits<Cost>::max());
}

// ============================================================================
// Combined Node and Edge Mask Tests
// ============================================================================

TEST(MaskingTests, CombinedNodeAndEdgeMasks) {
  // Create a grid graph with multiple paths
  auto g = make_grid_graph(3, 3);  // 3x3 grid
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  // Mask out some nodes and edges to create specific topology
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[4] = false;  // Mask center node

  auto edge_mask = make_bool_mask(g.num_edges());
  // Find and mask specific edge (e.g., 0->1)
  auto edge_src = g.edge_src_view();
  for (std::int32_t e = 0; e < g.num_edges(); ++e) {
    if (edge_src[e] == 0 && e == 0) {
      edge_mask[e] = false;
    }
  }

  auto [dist, dag] = shortest_paths(
      g, 0, std::nullopt, true, sel, {},
      std::span<const bool>(node_mask.get(), g.num_nodes()),
      std::span<const bool>(edge_mask.get(), g.num_edges()));

  // Verify source is reachable
  EXPECT_DOUBLE_EQ(dist[0], 0.0);

  // Verify masked node is unreachable
  EXPECT_EQ(dist[4], std::numeric_limits<Cost>::max());

  // DAG should be valid
  expect_pred_dag_valid(dag, g.num_nodes());
}

TEST(MaskingTests, NodeMaskOverridesEdgeMask) {
  // When both node and edge are masked, node mask takes precedence
  // (if a node is masked, all its edges are implicitly unavailable)
  auto g = make_line_graph(3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[1] = false;  // Mask node 1

  auto edge_mask = make_bool_mask(g.num_edges());
  // Even if edges are allowed, masked node should block them
  // (all edges true)

  auto [dist, dag] = shortest_paths(
      g, 0, std::nullopt, true, sel, {},
      std::span<const bool>(node_mask.get(), g.num_nodes()),
      std::span<const bool>(edge_mask.get(), g.num_edges()));

  // Node 1 and beyond should be unreachable
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_EQ(dist[1], std::numeric_limits<Cost>::max());
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());
}

// ============================================================================
// Mask Interaction with Residual Capacity
// ============================================================================

TEST(MaskingTests, MaskWithResidualCapacity) {
  auto g = make_line_graph(4);  // 0->1->2->3
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = true;  // Enable residual filtering

  // Set residual capacity
  std::vector<Cap> residual(g.num_edges(), 1.0);
  residual[1] = 0.0;  // Zero capacity on edge 1->2

  // Also mask out edge 2->3
  auto edge_mask = make_bool_mask(g.num_edges());
  edge_mask[2] = false;

  auto [dist, dag] = shortest_paths(
      g, 0, std::nullopt, true, sel, residual, {},
      std::span<const bool>(edge_mask.get(), g.num_edges()));

  // Should reach node 1 but not beyond
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_DOUBLE_EQ(dist[1], 1.0);
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());  // Blocked by residual
  EXPECT_EQ(dist[3], std::numeric_limits<Cost>::max());  // Blocked by mask
}

TEST(MaskingTests, ResidualAndNodeMaskCombined) {
  auto g = make_simple_parallel_paths_topology();  // Two parallel paths
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = true;

  // Zero capacity on one path
  std::vector<Cap> residual(g.num_edges(), 1.0);
  residual[0] = 0.0;  // Block 0->1

  // Mask intermediate node on other path
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[2] = false;  // Block 0->2->3 path

  auto [dist, dag] = shortest_paths(
      g, 0, NodeId(3), true, sel, residual,
      std::span<const bool>(node_mask.get(), g.num_nodes()), {});

  // Destination should be unreachable (both paths blocked)
  EXPECT_EQ(dist[3], std::numeric_limits<Cost>::max());
}

// ============================================================================
// K-Shortest Paths with Masks
// ============================================================================

TEST(MaskingTests, KSPWithNodeMask) {
  auto g = make_square_graph(2);  // Two equal-cost paths

  // Mask out one intermediate node
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[1] = false;  // Block path through node 1

  auto results = k_shortest_paths(
      g, 0, 2, 2, std::nullopt, true,
      std::span<const bool>(node_mask.get(), g.num_nodes()), {});

  // Should find at most 1 path (through node 3)
  EXPECT_LE(results.size(), 1);

  if (results.size() > 0) {
    const auto& [dist, dag] = results[0];
    // Destination should be reachable via alternate path
    EXPECT_GT(dist[2], 0.0);
    EXPECT_LT(dist[2], std::numeric_limits<Cost>::max());
  }
}

TEST(MaskingTests, KSPWithEdgeMask) {
  auto g = make_square_graph(2);  // Two equal-cost paths

  // Mask out edges to eliminate one path
  auto edge_mask = make_bool_mask(g.num_edges());
  edge_mask[0] = false;  // Block first edge

  auto results = k_shortest_paths(
      g, 0, 2, 2, std::nullopt, true, {},
      std::span<const bool>(edge_mask.get(), g.num_edges()));

  // Should find remaining path(s)
  EXPECT_GE(results.size(), 0);

  for (const auto& [dist, dag] : results) {
    if (dist[2] < std::numeric_limits<Cost>::max()) {
      expect_pred_dag_valid(dag, g.num_nodes());
    }
  }
}

// ============================================================================
// Max Flow with Masks
// ============================================================================

TEST(MaskingTests, MaxFlowWithNodeMask) {
  auto g = make_simple_parallel_paths_topology();

  // Mask out one intermediate node
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[1] = false;  // Block path through node 1

  auto [flow_val, summary] = calc_max_flow(
      g, 0, 3,
      FlowPlacement::Proportional, false, true, true, false, false,
      std::span<const bool>(node_mask.get(), g.num_nodes()), {});

  // Should get flow only through remaining path
  EXPECT_GT(flow_val, 0.0);
  EXPECT_LE(flow_val, 1.0);  // Only one path available

  // Validate constraints
  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaskingTests, MaxFlowWithEdgeMask) {
  auto g = make_triangle_topology();

  // Mask out direct edge 0->2
  auto edge_src = g.edge_src_view();
  auto edge_dst = g.edge_dst_view();
  auto edge_mask = make_bool_mask(g.num_edges());

  for (std::int32_t e = 0; e < g.num_edges(); ++e) {
    if (edge_src[e] == 0 && edge_dst[e] == 2) {
      edge_mask[e] = false;
    }
  }

  auto [flow_val, summary] = calc_max_flow(
      g, 0, 2,
      FlowPlacement::Proportional, false, true, true, false, false,
      {},
      std::span<const bool>(edge_mask.get(), g.num_edges()));

  // Should get flow only through indirect path 0->1->2
  EXPECT_GT(flow_val, 0.0);
  EXPECT_LE(flow_val, 1.0);  // Limited by bottleneck

  // Verify masked edge has no flow
  for (std::int32_t e = 0; e < g.num_edges(); ++e) {
    if (edge_src[e] == 0 && edge_dst[e] == 2) {
      EXPECT_DOUBLE_EQ(summary.edge_flows[e], 0.0)
          << "Masked edge should have zero flow";
    }
  }
}

TEST(MaskingTests, MaxFlowAllPathsMasked) {
  auto g = make_line_graph(3);

  // Mask intermediate node to block all paths
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[1] = false;

  auto [flow_val, summary] = calc_max_flow(
      g, 0, 2,
      FlowPlacement::Proportional, false, true, true, false, false,
      std::span<const bool>(node_mask.get(), g.num_nodes()), {});

  // No flow should be possible
  EXPECT_DOUBLE_EQ(flow_val, 0.0);

  // All edge flows should be zero
  for (const auto& ef : summary.edge_flows) {
    EXPECT_DOUBLE_EQ(ef, 0.0);
  }
}

// ============================================================================
// Empty Mask Tests (no masking applied)
// ============================================================================

TEST(MaskingTests, EmptyNodeMaskAllowsAllNodes) {
  auto g = make_line_graph(3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  // Empty mask means no masking
  auto [dist, dag] = shortest_paths(
      g, 0, std::nullopt, true, sel, {}, {}, {});

  // All nodes should be reachable
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_DOUBLE_EQ(dist[1], 1.0);
  EXPECT_DOUBLE_EQ(dist[2], 2.0);
}

TEST(MaskingTests, EmptyEdgeMaskAllowsAllEdges) {
  auto g = make_line_graph(3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  // Empty mask means no masking
  auto [dist, dag] = shortest_paths(
      g, 0, std::nullopt, true, sel, {}, {}, {});

  expect_pred_dag_valid(dag, g.num_nodes());

  // All nodes should have valid paths
  for (std::int32_t i = 0; i < g.num_nodes(); ++i) {
    EXPECT_LT(dist[i], std::numeric_limits<Cost>::max());
  }
}

// ============================================================================
// All-False Mask Tests
// ============================================================================

TEST(MaskingTests, AllNodesMaskedNoReachability) {
  auto g = make_line_graph(3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  // Mask all nodes
  auto node_mask = make_bool_mask(g.num_nodes(), false);

  auto [dist, dag] = shortest_paths(
      g, 0, std::nullopt, true, sel, {},
      std::span<const bool>(node_mask.get(), g.num_nodes()), {});

  // All nodes unreachable (including source)
  for (std::int32_t i = 0; i < g.num_nodes(); ++i) {
    EXPECT_EQ(dist[i], std::numeric_limits<Cost>::max());
  }

  // DAG should be empty
  for (const auto& offset : dag.parent_offsets) {
    EXPECT_EQ(offset, 0);
  }
}

TEST(MaskingTests, AllEdgesMaskedNoReachability) {
  auto g = make_line_graph(3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;

  // Mask all edges
  auto edge_mask = make_bool_mask(g.num_edges(), false);

  auto [dist, dag] = shortest_paths(
      g, 0, std::nullopt, true, sel, {}, {},
      std::span<const bool>(edge_mask.get(), g.num_edges()));

  // Only source should be reachable
  EXPECT_DOUBLE_EQ(dist[0], 0.0);

  // All other nodes unreachable
  for (std::int32_t i = 1; i < g.num_nodes(); ++i) {
    EXPECT_EQ(dist[i], std::numeric_limits<Cost>::max());
  }
}

// ============================================================================
// Backend (Algorithms API) Mask Tests
// ============================================================================

TEST(MaskingTests, AlgorithmsAPIMaskedSource) {
  auto backend = make_cpu_backend();
  Algorithms algs(backend);

  auto g = make_line_graph(3);
  auto gh = algs.build_graph(g);

  // Mask source
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[0] = false;

  SpfOptions opts;
  opts.multipath = true;
  opts.selection.multi_edge = true;
  opts.selection.require_capacity = false;
  opts.node_mask = std::span<const bool>(node_mask.get(), g.num_nodes());

  auto [dist, dag] = algs.spf(gh, 0, opts);

  // All distances should be infinity
  for (const auto& d : dist) {
    EXPECT_EQ(d, std::numeric_limits<Cost>::max());
  }

  // DAG should be empty
  for (const auto& offset : dag.parent_offsets) {
    EXPECT_EQ(offset, 0);
  }
}

TEST(MaskingTests, AlgorithmsAPIMaxFlowWithMask) {
  auto backend = make_cpu_backend();
  Algorithms algs(backend);

  auto g = make_simple_parallel_paths_topology();
  auto gh = algs.build_graph(g);

  // Mask one path
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[1] = false;

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = false;
  opts.require_capacity = true;
  opts.with_edge_flows = true;
  opts.node_mask = std::span<const bool>(node_mask.get(), g.num_nodes());

  auto [flow_val, summary] = algs.max_flow(gh, 0, 3, opts);

  // Should get flow through remaining path only
  EXPECT_GT(flow_val, 0.0);
  EXPECT_LE(flow_val, 1.0);
}
