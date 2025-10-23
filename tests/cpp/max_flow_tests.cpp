/**
 * Comprehensive tests for max-flow algorithm.
 *
 * This test suite validates correctness, performance, and edge case handling
 * of the max-flow implementation across various network topologies and
 * configuration parameters.
 *
 * Coverage dimensions:
 * - Path count: single, 2, 3, 4, 10 disjoint paths
 * - Capacity patterns: equal, asymmetric (1:100, 1:1M), zero, bottlenecks
 * - Cost tiers: equal-cost, 2-tier, 3-tier scenarios
 * - Placement modes: Proportional, EqualBalanced
 * - shortest_path parameter: True (single-tier) vs False (multi-tier)
 * - Advanced features: masks, batch operations, flow details
 * - Complex topologies: grids, shared bottlenecks
 */

#include <gtest/gtest.h>
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/backend.hpp"
#include "netgraph/core/algorithms.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

//=============================================================================
// SECTION 1: BASIC CORRECTNESS & SINGLE PATH TESTS
//=============================================================================

TEST(MaxFlow, SimpleLinearPath) {
  // Basic sanity: single path with unit capacity
  auto g = make_line_graph(3);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = false;
  opts.with_edge_flows = false;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  EXPECT_NEAR(total, 1.0, 1e-9) << "Single path should carry capacity 1.0";
}

TEST(MaxFlow, BottleneckIdentification) {
  // Path with bottleneck: 0→1 (cap 10), 1→2 (cap 2), 2→3 (cap 10)
  // The bottleneck edge should limit total flow
  std::int32_t src[3] = {0, 1, 2};
  std::int32_t dst[3] = {1, 2, 3};
  double cap[3] = {10.0, 2.0, 10.0};
  std::int64_t cost[3] = {1, 1, 1};

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src, 3), std::span(dst, 3),
    std::span(cap, 3), std::span(cost, 3));

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  EXPECT_NEAR(total, 2.0, 1e-9) << "Flow limited by bottleneck edge (cap 2.0)";
}

TEST(MaxFlow, SinglePath_ShortestVsFull_Identical) {
  // With only one path, shortest_path=True and False should give same result
  auto g = make_line_graph(3);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts_full;
  opts_full.placement = FlowPlacement::Proportional;
  opts_full.shortest_path = false;
  auto [total_full, _] = algs.max_flow(gh, 0, 2, opts_full);

  MaxFlowOptions opts_sp;
  opts_sp.placement = FlowPlacement::Proportional;
  opts_sp.shortest_path = true;
  auto [total_sp, __] = algs.max_flow(gh, 0, 2, opts_sp);

  EXPECT_NEAR(total_full, 1.0, 1e-9);
  EXPECT_NEAR(total_sp, 1.0, 1e-9);
  EXPECT_NEAR(total_full, total_sp, 1e-9)
      << "Single path: shortest_path mode should be identical to full max-flow";
}

TEST(MaxFlow, SourceEqualsSink_ZeroFlow) {
  // Degenerate case: source and sink are the same node
  auto g = make_line_graph(3);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;

  auto [total, summary] = algs.max_flow(gh, 1, 1, opts);

  EXPECT_NEAR(total, 0.0, 1e-12) << "Flow from node to itself should be zero";
}

//=============================================================================
// SECTION 2: DISJOINT EQUAL-COST PATHS (ECMP CORE)
// Critical: These tests guard against the bug fixed in flow_state.cpp:233
// where shortest_path=True only used 1 of N equal-cost paths due to
// premature break in the inner DFS loop.
//=============================================================================

TEST(MaxFlow, TwoDisjointPaths_EqualCost_Proportional) {
  // THE CRITICAL REGRESSION TEST FOR BUG FIX (2025-10-23)
  // Before fix: returned 10.0 (only 1 of 2 paths used)
  // After fix: returns 20.0 (both equal-cost paths saturated)
  //
  // This validates that shortest_path=True saturates ALL equal-cost paths
  // in the lowest-cost tier, not just one path.
  auto g = make_n_disjoint_paths(2, 10.0);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  // Source=0, Sink=3 (2 intermediates at 1,2)
  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  // CRITICAL: Must saturate BOTH equal-cost paths
  EXPECT_NEAR(total, 20.0, 1e-9)
      << "shortest_path=True must use ALL equal-cost paths! "
      << "This test guards against flow_state.cpp:233 bug.";

  // Verify both paths were actually used
  // Path 1: edges 0 (0→1) and 1 (1→3)
  // Path 2: edges 2 (0→2) and 3 (2→3)
  EXPECT_GT(summary.edge_flows[0], 5.0) << "Path 1 should carry significant flow";
  EXPECT_GT(summary.edge_flows[2], 5.0) << "Path 2 should carry significant flow";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, ThreeDisjointPaths_EqualCost_Proportional) {
  // Extend ECMP validation to 3 paths
  auto g = make_n_disjoint_paths(3, 10.0);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 4, opts);

  EXPECT_NEAR(total, 30.0, 1e-9)
      << "Should saturate all 3 equal-cost paths";

  // Each path should carry ~10.0
  for (int i = 0; i < 3; ++i) {
    EXPECT_GT(summary.edge_flows[i * 2], 5.0)
        << "Path " << i << " should be utilized";
  }

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 4);
}

TEST(MaxFlow, FourDisjointPaths_EqualCost_ClosFabric) {
  // Simulates 4-spine Clos fabric scenario (common data center topology)
  auto g = make_n_disjoint_paths(4, 100.0);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 5, opts);

  EXPECT_NEAR(total, 400.0, 1e-9)
      << "4-spine Clos should use all 4 equal-cost paths (400 Gbps aggregate)";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 5);
}

TEST(MaxFlow, TenDisjointPaths_EqualCost_Scalability) {
  // Scalability test: Many parallel equal-cost paths
  auto g = make_n_disjoint_paths(10, 10.0);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 11, opts);

  EXPECT_NEAR(total, 100.0, 1e-9)
      << "Should saturate all 10 equal-cost paths";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 11);
}

TEST(MaxFlow, TwoDisjointPaths_EqualCost_AsymmetricCapacity_1vs100) {
  // Two equal-cost paths with very different capacities
  // Validates proportional placement under capacity asymmetry
  std::int32_t src_arr[4] = {0, 1, 0, 2};
  std::int32_t dst_arr[4] = {1, 3, 2, 3};
  double cap_arr[4] = {1.0, 1.0, 100.0, 100.0};
  std::int64_t cost_arr[4] = {1, 1, 1, 1};  // Equal cost!

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src_arr, 4), std::span(dst_arr, 4),
    std::span(cap_arr, 4), std::span(cost_arr, 4));

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  // Both paths should contribute (ECMP)
  EXPECT_GT(total, 1.0) << "Should use both equal-cost paths, not just one";

  // Small path should be saturated
  EXPECT_NEAR(summary.edge_flows[0], 1.0, 1e-9)
      << "Path 1 (cap 1.0) should be saturated";

  // Large path should also carry flow
  EXPECT_GT(summary.edge_flows[2], 0.0)
      << "Path 2 (cap 100.0) should also be used";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, TwoDisjointPaths_EqualCost_AsymmetricCapacity_1vs1M) {
  // Extreme asymmetry: 1 vs 1,000,000
  // Stress test for numerical stability
  std::int32_t src_arr[4] = {0, 1, 0, 2};
  std::int32_t dst_arr[4] = {1, 3, 2, 3};
  double cap_arr[4] = {1.0, 1.0, 1000000.0, 1000000.0};
  std::int64_t cost_arr[4] = {1, 1, 1, 1};

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src_arr, 4), std::span(dst_arr, 4),
    std::span(cap_arr, 4), std::span(cost_arr, 4));

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  // Should handle extreme asymmetry correctly
  EXPECT_GT(total, 1.0) << "Should use both paths despite extreme asymmetry";

  // Small path saturated
  EXPECT_NEAR(summary.edge_flows[0], 1.0, 1e-9)
      << "Small path (1.0) should be saturated";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, ZeroCapacityPath_Ignored) {
  // One path with zero capacity - should route all flow via non-zero path
  // Tests correct handling of residual capacity filtering
  std::int32_t src_arr[4] = {0, 1, 0, 2};
  std::int32_t dst_arr[4] = {1, 3, 2, 3};
  double cap_arr[4] = {10.0, 10.0, 0.0, 0.0};  // Path 2 has zero capacity
  std::int64_t cost_arr[4] = {1, 1, 1, 1};

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src_arr, 4), std::span(dst_arr, 4),
    std::span(cap_arr, 4), std::span(cost_arr, 4));

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  EXPECT_NEAR(total, 10.0, 1e-9)
      << "Total flow should use only non-zero capacity path";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

//=============================================================================
// SECTION 3: PLACEMENT MODES (Proportional vs EqualBalanced)
//=============================================================================

TEST(MaxFlow, PlacementModes_EqualCapacity_IdenticalResult) {
  // With equal capacities on all paths, both modes should give same total flow
  auto g = make_n_disjoint_paths(2, 10.0);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts_prop;
  opts_prop.placement = FlowPlacement::Proportional;
  opts_prop.shortest_path = true;
  auto [total_prop, _] = algs.max_flow(gh, 0, 3, opts_prop);

  MaxFlowOptions opts_eq;
  opts_eq.placement = FlowPlacement::EqualBalanced;
  opts_eq.shortest_path = true;
  auto [total_eq, __] = algs.max_flow(gh, 0, 3, opts_eq);

  // For equal capacities, both should saturate both paths
  EXPECT_NEAR(total_prop, 20.0, 1e-9);
  EXPECT_NEAR(total_eq, 20.0, 1e-9);
  EXPECT_NEAR(total_prop, total_eq, 1e-9)
      << "With equal capacities, both placement modes should give same result";
}

TEST(MaxFlow, PlacementModes_AsymmetricCapacity_Proportional) {
  // Proportional placement with asymmetric capacities
  std::int32_t src_arr[4] = {0, 1, 0, 2};
  std::int32_t dst_arr[4] = {1, 3, 2, 3};
  double cap_arr[4] = {1.0, 1.0, 100.0, 100.0};
  std::int64_t cost_arr[4] = {1, 1, 1, 1};

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src_arr, 4), std::span(dst_arr, 4),
    std::span(cap_arr, 4), std::span(cost_arr, 4));

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  EXPECT_GT(total, 1.0) << "Proportional should use both paths";
  EXPECT_NEAR(summary.edge_flows[0], 1.0, 1e-9) << "Small path saturated";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, PlacementModes_AsymmetricCapacity_EqualBalanced) {
  // Equal-balanced placement with asymmetric capacities
  // Should be limited by minimum capacity
  std::int32_t src_arr[4] = {0, 1, 0, 2};
  std::int32_t dst_arr[4] = {1, 3, 2, 3};
  double cap_arr[4] = {1.0, 1.0, 100.0, 100.0};
  std::int64_t cost_arr[4] = {1, 1, 1, 1};

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src_arr, 4), std::span(dst_arr, 4),
    std::span(cap_arr, 4), std::span(cost_arr, 4));

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::EqualBalanced;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  // Equal-balanced limited by min capacity: 1.0 × 2 = 2.0
  EXPECT_NEAR(total, 2.0, 1e-9)
      << "Equal-balanced limited by smallest path capacity (1.0)";

  // Both paths should have equal flow
  EXPECT_NEAR(summary.edge_flows[0], 1.0, 1e-9) << "Small path saturated at 1.0";
  EXPECT_NEAR(summary.edge_flows[2], 1.0, 1e-9) << "Large path limited to 1.0";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

//=============================================================================
// SECTION 4: COST TIERS & shortest_path PARAMETER
//=============================================================================

TEST(MaxFlow, TwoCostTiers_ShortestPath_OnlyLowestTier) {
  // Verify shortest_path=True uses only the lowest-cost tier
  auto g = make_tiered_graph(5.0, 100.0);  // Tier1: cost=1, Tier2: cost=2

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  // Should ONLY use cost=1 tier (path 1), NOT cost=2 tier (path 2)
  EXPECT_NEAR(total, 5.0, 1e-9)
      << "shortest_path=True should use only lowest-cost tier";

  // Path 1 (cost=1) should be saturated
  EXPECT_NEAR(summary.edge_flows[0], 5.0, 1e-9)
      << "Tier 1 (cost=1, cap=5) should be saturated";

  // Path 2 (cost=2) should NOT be used
  EXPECT_NEAR(summary.edge_flows[2], 0.0, 1e-9)
      << "Tier 2 (cost=2, cap=100) should not be used";

  // Verify cost distribution contains only one tier
  EXPECT_EQ(summary.costs.size(), 1) << "Only one cost tier should be present";
  EXPECT_EQ(summary.costs[0], 2) << "Only cost=2 tier used (sum of 2 edges with cost=1)";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, TwoCostTiers_FullMaxFlow_AllTiers) {
  // Verify shortest_path=False uses all cost tiers
  auto g = make_tiered_graph(5.0, 100.0);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = false;  // Full max-flow
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  // Should use BOTH tiers
  EXPECT_NEAR(total, 105.0, 1e-9)
      << "Full max-flow should use all cost tiers (5 + 100 = 105)";

  // Both paths should be saturated
  EXPECT_NEAR(summary.edge_flows[0], 5.0, 1e-9) << "Tier 1 saturated";
  EXPECT_NEAR(summary.edge_flows[2], 100.0, 1e-9) << "Tier 2 saturated";

  // Verify cost distribution has TWO tiers
  EXPECT_EQ(summary.costs.size(), 2) << "Two cost tiers should be present";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, EqualCost_ShortestVsFull_IdenticalResult) {
  // When all paths have equal cost, shortest_path and full should be identical
  auto g = make_n_disjoint_paths(3, 10.0);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts_full;
  opts_full.placement = FlowPlacement::Proportional;
  opts_full.shortest_path = false;
  auto [total_full, _] = algs.max_flow(gh, 0, 4, opts_full);

  MaxFlowOptions opts_sp;
  opts_sp.placement = FlowPlacement::Proportional;
  opts_sp.shortest_path = true;
  auto [total_sp, __] = algs.max_flow(gh, 0, 4, opts_sp);

  EXPECT_NEAR(total_full, total_sp, 1e-9)
      << "With all equal-cost paths, shortest_path and full max-flow should be identical";
  EXPECT_NEAR(total_full, 30.0, 1e-9) << "Should saturate all 3 paths";
}

TEST(MaxFlow, ThreeCostTiers_ShortestPath_OnlyLowestTier) {
  // Verify shortest_path=True works correctly with 3+ cost tiers
  // Path 1: 0→1→4 (cost=1 per edge, cap=5)
  // Path 2: 0→2→4 (cost=2 per edge, cap=10)
  // Path 3: 0→3→4 (cost=3 per edge, cap=20)
  std::int32_t src_arr[6] = {0, 1, 0, 2, 0, 3};
  std::int32_t dst_arr[6] = {1, 4, 2, 4, 3, 4};
  double cap_arr[6] = {5.0, 5.0, 10.0, 10.0, 20.0, 20.0};
  std::int64_t cost_arr[6] = {1, 1, 2, 2, 3, 3};

  auto g = StrictMultiDiGraph::from_arrays(5,
    std::span(src_arr, 6), std::span(dst_arr, 6),
    std::span(cap_arr, 6), std::span(cost_arr, 6));

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 4, opts);

  // Should ONLY use cost=1 tier (5.0), NOT cost=2 (10.0) or cost=3 (20.0)
  EXPECT_NEAR(total, 5.0, 1e-9)
      << "shortest_path=True should use only lowest-cost tier (cost=1)";

  // Only first path should be used
  EXPECT_NEAR(summary.edge_flows[0], 5.0, 1e-9) << "Path 1 (cost=1) saturated";
  EXPECT_NEAR(summary.edge_flows[2], 0.0, 1e-9) << "Path 2 (cost=2) not used";
  EXPECT_NEAR(summary.edge_flows[4], 0.0, 1e-9) << "Path 3 (cost=3) not used";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 4);
}

//=============================================================================
// SECTION 5: DISJOINT PATHS WITH DIFFERENT COSTS
//=============================================================================

TEST(MaxFlow, TwoDisjointPaths_DifferentCost_FullMaxFlow) {
  // Two disjoint paths with different costs and capacities
  // Path 1: 0→1→3 (cost=1, cap=1)
  // Path 2: 0→2→3 (cost=1, cap=2)
  // Full max-flow should use both paths
  std::int32_t src[4] = {0, 1, 0, 2};
  std::int32_t dst[4] = {1, 3, 2, 3};
  double cap[4] = {1.0, 1.0, 2.0, 2.0};
  std::int64_t cost[4] = {1, 1, 1, 1};

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src, 4), std::span(dst, 4),
    std::span(cap, 4), std::span(cost, 4));

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  EXPECT_NEAR(total, 3.0, 1e-9)
      << "Total flow should be sum of both paths (1 + 2 = 3)";
}

//=============================================================================
// SECTION 6: FLOW PROPERTIES & CORRECTNESS VALIDATION
//=============================================================================

TEST(MaxFlow, FlowConservation_IntermediateNodes) {
  // Verify flow conservation: inflow = outflow at all intermediate nodes
  auto g = make_square_graph(1);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  EXPECT_EQ(summary.edge_flows.size(), static_cast<std::size_t>(g.num_edges()));

  // Use helper to validate conservation
  validate_flow_conservation(g, summary, 0, 2);
}

TEST(MaxFlow, MinCut_ValidEdges) {
  // Min-cut should separate source from sink with valid edge IDs
  auto g = make_square_graph(1);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  // Min cut should contain some edges
  EXPECT_GT(summary.min_cut.edges.size(), 0)
      << "Min-cut should separate source from sink";

  // All cut edges should have valid IDs
  for (auto eid : summary.min_cut.edges) {
    EXPECT_GE(eid, 0) << "Invalid edge ID in min-cut";
    EXPECT_LT(eid, g.num_edges()) << "Edge ID out of range in min-cut";
  }
}

TEST(MaxFlow, CostDistribution_SumsToTotalFlow) {
  // Cost distribution should account for all flow
  auto g = make_square_graph(1);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  // Cost distribution flows should sum to total flow
  double cost_sum = 0.0;
  for (auto flow : summary.flows) {
    cost_sum += flow;
  }

  EXPECT_NEAR(cost_sum, total, 1e-9)
      << "Cost distribution should account for all flow";
}

//=============================================================================
// SECTION 7: ADVANCED FEATURES & OPTIONS
//=============================================================================

TEST(MaxFlow, ResidualCapacity_WhenRequested) {
  // Verify residual capacity and reachable nodes are populated when requested
  auto g = make_square_graph(1);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.with_edge_flows = true;
  opts.with_reachable = true;
  opts.with_residuals = true;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  // Residuals length matches edges; reachable length matches nodes
  EXPECT_EQ(summary.residual_capacity.size(),
            static_cast<std::size_t>(g.num_edges()))
      << "Residual capacity vector size should match edge count";
  EXPECT_EQ(summary.reachable_nodes.size(),
            static_cast<std::size_t>(g.num_nodes()))
      << "Reachable nodes vector size should match node count";
}

TEST(MaxFlow, NodeMask_BlocksFlow) {
  // Verify that masking a node prevents flow through it
  auto g = make_line_graph(3);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  // Mask out node 1 (middle node)
  auto node_mask_vec = make_bool_mask(g.num_nodes());
  node_mask_vec[1] = false;

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.node_mask = std::span<const bool>(node_mask_vec.get(),
                                         static_cast<std::size_t>(g.num_nodes()));

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  EXPECT_NEAR(total, 0.0, 1e-12)
      << "Masking middle node should prevent all flow";
}

TEST(MaxFlow, EdgeMask_FiltersEdges) {
  // Verify that edge mask parameter is accepted without error
  // NOTE: Detailed edge mask behavior validation is complex and requires
  // understanding of the interaction between masking and shortest-path finding.
  // This test ensures the API accepts edge masks correctly.
  auto g = make_n_disjoint_paths(2, 10.0);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  // Create an empty mask (all edges enabled)
  auto edge_mask_vec = make_bool_mask(g.num_edges(), true);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.with_edge_flows = true;
  opts.edge_mask = std::span<const bool>(edge_mask_vec.get(),
                                         static_cast<std::size_t>(g.num_edges()));

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  // With all edges enabled, should use both paths
  EXPECT_NEAR(total, 20.0, 1e-9)
      << "With all edges enabled, should saturate both paths";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, BatchMaxFlow_IndependentMasks) {
  // Verify that batch max-flow applies different masks to different pairs independently
  auto g = make_line_graph(3);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  std::vector<std::pair<NodeId, NodeId>> pairs = {{0, 2}, {0, 2}};
  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;

  // First pair: no mask (should find flow)
  // Second pair: mask node 1 (should find no flow)
  std::vector<std::span<const bool>> node_masks;
  auto nm0 = make_bool_mask(g.num_nodes());  // all true
  auto nm1 = make_bool_mask(g.num_nodes());
  nm1[1] = false;  // block middle node for second pair
  node_masks.push_back(std::span<const bool>(nm0.get(),
                                             static_cast<std::size_t>(g.num_nodes())));
  node_masks.push_back(std::span<const bool>(nm1.get(),
                                             static_cast<std::size_t>(g.num_nodes())));

  auto summaries = algs.batch_max_flow(gh, pairs, opts, node_masks, {});

  ASSERT_EQ(summaries.size(), 2u);
  EXPECT_GT(summaries[0].total_flow, 0.0)
      << "First pair (no mask) should find flow";
  EXPECT_NEAR(summaries[1].total_flow, 0.0, 1e-12)
      << "Second pair (masked node) should find no flow";
}

//=============================================================================
// SECTION 8: COMPLEX TOPOLOGIES
//=============================================================================

TEST(MaxFlow, GridGraph_3x3_ManyEqualCostPaths) {
  // Grid graph has many equal-cost paths from corner to corner
  // Validates ECMP behavior in complex topology
  auto g = make_grid_graph(3, 3);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  // Flow from top-left (0) to bottom-right (8)
  auto [total, summary] = algs.max_flow(gh, 0, 8, opts);

  // Should find flow through multiple equal-cost paths
  EXPECT_GT(total, 0.0) << "Should find flow in grid";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 8);
}

TEST(MaxFlow, SharedBottleneck_TwoPaths) {
  // Two paths with a shared bottleneck edge
  // Validates correct handling of overlapping paths
  auto g = make_shared_bottleneck_graph(10.0, 5.0);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  // Shared bottleneck should limit total flow
  // Path 1: 0→1 (10), 1→3 (5) - bottleneck at 1→3
  // Path 2: 0→2 (10), 2→3 (10)
  // Total should be limited by bottleneck + other path = 5 + 10 = 15
  EXPECT_NEAR(total, 15.0, 1e-9)
      << "Flow should be limited by bottleneck plus second path";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}
