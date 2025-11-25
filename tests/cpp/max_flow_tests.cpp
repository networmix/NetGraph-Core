/**
 * Comprehensive Max-Flow Algorithm Test Suite
 *
 * This suite systematically validates max-flow correctness across all parameter
 * combinations and network topologies. The test matrix ensures complete coverage
 * of algorithm behavior under various configurations.
 *
 * TEST MATRIX COVERAGE:
 *
 * Dimension 1: Path Structure
 *   - Single path
 *   - Multiple disjoint equal-cost paths (2, 3, 4, 10)
 *   - Multiple cost tiers (2-tier, 3-tier)
 *   - Complex topologies (grids, shared bottlenecks)
 *
 * Dimension 2: Capacity Patterns
 *   - Equal capacities across paths
 *   - Asymmetric capacities (1:100, 1:1M ratios)
 *   - Zero capacity edges
 *   - Bottleneck scenarios
 *
 * Dimension 3: Flow Placement
 *   - PROPORTIONAL: Capacity-weighted distribution
 *   - EQUAL_BALANCED: Uniform distribution across paths
 *
 * Dimension 4: Path Selection
 *   - shortest_path=False: Multi-tier max-flow (use all cost tiers)
 *   - shortest_path=True: Single-tier flow (lowest cost tier only)
 *
 * Dimension 5: Optional Features
 *   - Node/edge masking
 *   - Batch operations
 *   - Flow details (edge flows, residuals, min-cut)
 *
 * INVARIANT PROPERTIES VALIDATED:
 *   - Flow conservation at intermediate nodes
 *   - Capacity constraints on all edges
 *   - Cost distribution consistency
 *   - Min-cut correctness
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
// SECTION 1: BASIC CORRECTNESS - SINGLE PATH SCENARIOS
//=============================================================================

TEST(MaxFlow, SinglePath_BasicLinear) {
  // Validates basic max-flow computation on simplest topology
  auto g = make_line_graph(3);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = false;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  EXPECT_NEAR(total, 1.0, 1e-9) << "Single path with unit capacity";
}

TEST(MaxFlow, SinglePath_BottleneckConstraint) {
  // Validates that bottleneck edge correctly limits total flow
  // Topology: 0→1 (cap 10), 1→2 (cap 2), 2→3 (cap 10)
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

  EXPECT_NEAR(total, 2.0, 1e-9) << "Flow limited by bottleneck edge capacity";
}

TEST(MaxFlow, SinglePath_ShortestPathEquivalence) {
  // For single-path networks, shortest_path=True/False should be identical
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
      << "shortest_path flag should not affect single-path networks";
}

TEST(MaxFlow, DegenerateCase_SourceEqualsSink) {
  // Validates handling of degenerate case where source = sink
  auto g = make_line_graph(3);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;

  auto [total, summary] = algs.max_flow(gh, 1, 1, opts);

  EXPECT_NEAR(total, 0.0, 1e-12) << "Flow from node to itself must be zero";
}

//=============================================================================
// SECTION 2: EQUAL-COST MULTIPATH (ECMP) - CORE VALIDATION
//
// This section systematically validates that when multiple paths have equal
// cost, the algorithm correctly utilizes ALL paths to maximize throughput.
// This is critical for data center fabrics (Clos, spine-leaf) where ECMP
// load balancing across equal-cost paths is the primary traffic distribution
// mechanism.
//=============================================================================

TEST(MaxFlow, ECMP_TwoPaths_Proportional_ShortestPath) {
  // Validates ECMP with 2 parallel equal-cost paths
  // Topology: Two disjoint paths from 0→3, each with capacity 10.0
  // shortest_path=True should saturate both paths (not just one)
  auto g = make_n_disjoint_paths(2, 10.0);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  EXPECT_NEAR(total, 20.0, 1e-9)
      << "shortest_path=True must saturate all equal-cost paths";

  // Verify both paths carry significant flow
  EXPECT_GT(summary.edge_flows[0], 5.0) << "Path 1 should be utilized";
  EXPECT_GT(summary.edge_flows[2], 5.0) << "Path 2 should be utilized";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, ECMP_TwoPaths_Proportional_FullMaxFlow) {
  // Same topology, but shortest_path=False for comparison
  auto g = make_n_disjoint_paths(2, 10.0);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = false;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  EXPECT_NEAR(total, 20.0, 1e-9) << "Full max-flow should use both paths";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, ECMP_TwoPaths_EqualBalanced_ShortestPath) {
  // Same topology with EQUAL_BALANCED placement
  // Should distribute flow uniformly across both paths
  auto g = make_n_disjoint_paths(2, 10.0);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::EqualBalanced;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  EXPECT_NEAR(total, 20.0, 1e-9) << "Equal-balanced should use both paths";

  // With equal capacities, both paths should have equal flow
  EXPECT_NEAR(summary.edge_flows[0], 10.0, 1e-9) << "Path 1 should be saturated";
  EXPECT_NEAR(summary.edge_flows[2], 10.0, 1e-9) << "Path 2 should be saturated";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, ECMP_ThreePaths_Proportional_ShortestPath) {
  // Extends ECMP validation to 3 equal-cost paths
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

  // Each path should carry significant flow (~10.0)
  for (int i = 0; i < 3; ++i) {
    EXPECT_GT(summary.edge_flows[i * 2], 5.0)
        << "Path " << i << " should be utilized";
  }

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 4);
}

TEST(MaxFlow, ECMP_FourPaths_ClosFabric_Simulation) {
  // Simulates 4-spine Clos fabric (common data center topology)
  // Each spine provides an equal-cost path with 100 Gbps capacity
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
      << "4-spine Clos: aggregate bandwidth should be 400 Gbps";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 5);
}

TEST(MaxFlow, ECMP_TenPaths_ScalabilityValidation) {
  // Validates algorithm scalability with many parallel equal-cost paths
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

//=============================================================================
// SECTION 3: ASYMMETRIC CAPACITY PATTERNS
//
// Real networks often have paths with very different capacities due to
// heterogeneous equipment, upgrades, or failures. These tests validate
// correct behavior under capacity asymmetry.
//=============================================================================

TEST(MaxFlow, AsymmetricCapacity_1to100_Proportional) {
  // Two equal-cost paths with 100x capacity difference
  // Path 1: capacity 1.0, Path 2: capacity 100.0
  std::int32_t src_arr[4] = {0, 1, 0, 2};
  std::int32_t dst_arr[4] = {1, 3, 2, 3};
  double cap_arr[4] = {1.0, 1.0, 100.0, 100.0};
  std::int64_t cost_arr[4] = {1, 1, 1, 1};  // Equal cost

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

  // Should use both paths despite capacity asymmetry
  EXPECT_NEAR(total, 101.0, 1e-9) << "Total flow should be 1.0 + 100.0";

  // Validate per-edge flows by (u,v) lookup instead of relying on edge indices
  auto row = g.row_offsets_view();
  auto col = g.col_indices_view();
  auto aei = g.adj_edge_index_view();
  auto eid_uv = [&](int u, int v) -> std::size_t {
    auto s = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto e = static_cast<std::size_t>(row[static_cast<std::size_t>(u) + 1]);
    for (std::size_t j = s; j < e; ++j) {
      if (static_cast<int>(col[j]) == v) return static_cast<std::size_t>(aei[j]);
    }
    ADD_FAILURE() << "Edge (" << u << "->" << v << ") not found";
    return static_cast<std::size_t>(-1);
  };
  auto e01 = eid_uv(0, 1);
  auto e13 = eid_uv(1, 3);
  auto e02 = eid_uv(0, 2);
  auto e23 = eid_uv(2, 3);
  ASSERT_NE(e01, static_cast<std::size_t>(-1));
  ASSERT_NE(e13, static_cast<std::size_t>(-1));
  ASSERT_NE(e02, static_cast<std::size_t>(-1));
  ASSERT_NE(e23, static_cast<std::size_t>(-1));
  // Small path saturated
  EXPECT_NEAR(summary.edge_flows[e01], 1.0, 1e-9)
      << "Low-capacity first hop should be fully utilized";
  EXPECT_NEAR(summary.edge_flows[e13], 1.0, 1e-9)
      << "Low-capacity second hop should be fully utilized";
  // Large path saturated
  EXPECT_NEAR(summary.edge_flows[e02], 100.0, 1e-9)
      << "High-capacity first hop should be fully utilized";
  EXPECT_NEAR(summary.edge_flows[e23], 100.0, 1e-9)
      << "High-capacity second hop should be fully utilized";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, AsymmetricCapacity_1to1M_NumericalStability) {
  // Extreme capacity asymmetry: 1M:1 ratio
  // Tests numerical stability of proportional placement
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

  // Should handle extreme asymmetry without numerical issues
  EXPECT_NEAR(total, 1000001.0, 1e-6) << "Should use both paths";

  // Small path saturated
  EXPECT_NEAR(summary.edge_flows[0], 1.0, 1e-9)
      << "Small path should be saturated";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, AsymmetricCapacity_1to100_EqualBalanced) {
  // Equal-balanced placement with asymmetric capacities
  // Should be limited by minimum capacity across paths
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

  // Equal-balanced: both paths get same flow, limited by min capacity
  EXPECT_NEAR(total, 2.0, 1e-9)
      << "Equal-balanced: 1.0 per path × 2 paths = 2.0";

  // Both paths should have equal flow
  EXPECT_NEAR(summary.edge_flows[0], 1.0, 1e-9) << "Path 1 saturated at min";
  EXPECT_NEAR(summary.edge_flows[2], 1.0, 1e-9) << "Path 2 limited to min";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, EqualBalanced_NodeMask_ReducesECMPWidth) {
  // Two equal-cost paths 0->1->3 and 0->2->3 with cap=5 each.
  // Mask out node 2; EB should place only on path via 1, limited by 5.
  std::int32_t src_arr[4] = {0, 1, 0, 2};
  std::int32_t dst_arr[4] = {1, 3, 2, 3};
  double cap_arr[4] = {5.0, 5.0, 5.0, 5.0};
  std::int64_t cost_arr[4] = {1, 1, 1, 1};
  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src_arr, 4), std::span(dst_arr, 4),
    std::span(cap_arr, 4), std::span(cost_arr, 4));
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);
  MaxFlowOptions opts;
  opts.placement = FlowPlacement::EqualBalanced;
  opts.shortest_path = false;
  opts.with_edge_flows = true;
  // Node mask: disable node 2
  auto node_mask_vec = make_bool_mask(g.num_nodes());
  node_mask_vec[2] = false;
  opts.node_mask = std::span<const bool>(node_mask_vec.get(), static_cast<std::size_t>(g.num_nodes()));
  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);
  EXPECT_NEAR(total, 5.0, 1e-9) << "Only unmasked path should carry flow";
  // Validate that only path via 1 carries 5
  auto row = g.row_offsets_view();
  auto col = g.col_indices_view();
  auto aei = g.adj_edge_index_view();
  auto eid_uv = [&](int u, int v) -> std::size_t {
    auto s = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto e = static_cast<std::size_t>(row[static_cast<std::size_t>(u) + 1]);
    for (std::size_t j = s; j < e; ++j) {
      if (static_cast<int>(col[j]) == v) return static_cast<std::size_t>(aei[j]);
    }
    return static_cast<std::size_t>(-1);
  };
  auto e01 = eid_uv(0, 1);
  auto e13 = eid_uv(1, 3);
  auto e02 = eid_uv(0, 2);
  auto e23 = eid_uv(2, 3);
  ASSERT_NE(e01, static_cast<std::size_t>(-1));
  ASSERT_NE(e13, static_cast<std::size_t>(-1));
  ASSERT_NE(e02, static_cast<std::size_t>(-1));
  ASSERT_NE(e23, static_cast<std::size_t>(-1));
  EXPECT_NEAR(summary.edge_flows[e01], 5.0, 1e-9);
  EXPECT_NEAR(summary.edge_flows[e13], 5.0, 1e-9);
  EXPECT_NEAR(summary.edge_flows[e02], 0.0, 1e-9);
  EXPECT_NEAR(summary.edge_flows[e23], 0.0, 1e-9);
  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}
TEST(MaxFlow, ZeroCapacityPath_CorrectFiltering) {
  // One path with zero capacity should be ignored
  // Tests correct residual capacity filtering
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
      << "Should route all flow via non-zero capacity path";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

//=============================================================================
// SECTION 4: COST TIER VALIDATION
//
// Networks often have paths with different costs (distance, latency, monetary).
// The shortest_path parameter controls whether to use only lowest-cost paths
// or to utilize higher-cost paths for additional capacity.
//=============================================================================

TEST(MaxFlow, TwoCostTiers_ShortestPath_SingleTierOnly) {
  // Validates that shortest_path=True uses only lowest-cost tier
  // Tier 1: cost=1, capacity=5
  // Tier 2: cost=2, capacity=100
  auto g = make_tiered_graph(5.0, 100.0);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  // Should use ONLY cost=1 tier, NOT cost=2 tier
  EXPECT_NEAR(total, 5.0, 1e-9)
      << "shortest_path=True: use only lowest-cost tier";

  // Tier 1 (cost=1) should be saturated
  EXPECT_NEAR(summary.edge_flows[0], 5.0, 1e-9)
      << "Lowest-cost path should be saturated";

  // Tier 2 (cost=2) should NOT be used
  EXPECT_NEAR(summary.edge_flows[2], 0.0, 1e-9)
      << "Higher-cost tier should not be used";

  // Verify cost distribution contains only one tier
  EXPECT_EQ(summary.costs.size(), 1) << "Only one cost tier used";
  EXPECT_EQ(summary.costs[0], 2) << "Total cost = 2 (sum of 2 edges @ cost 1)";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, TwoCostTiers_FullMaxFlow_AllTiers) {
  // Validates that shortest_path=False uses all available tiers
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
      << "Full max-flow should use all tiers (5 + 100 = 105)";

  // Both tiers should be saturated
  EXPECT_NEAR(summary.edge_flows[0], 5.0, 1e-9) << "Tier 1 saturated";
  EXPECT_NEAR(summary.edge_flows[2], 100.0, 1e-9) << "Tier 2 saturated";

  // Verify cost distribution has TWO tiers
  EXPECT_EQ(summary.costs.size(), 2) << "Two cost tiers used";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, EqualCost_ShortestVsFull_Equivalent) {
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
      << "With all equal-cost paths, both modes should be equivalent";
  EXPECT_NEAR(total_full, 30.0, 1e-9) << "Should saturate all 3 paths";
}

TEST(MaxFlow, ThreeCostTiers_ShortestPath_LowestTierOnly) {
  // Validates shortest_path with 3+ cost tiers
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

  // Should use ONLY cost=1 tier
  EXPECT_NEAR(total, 5.0, 1e-9)
      << "Should use only lowest-cost tier (cost=1)";

  // Only first path should have flow
  EXPECT_NEAR(summary.edge_flows[0], 5.0, 1e-9) << "Path 1 (cost=1) saturated";
  EXPECT_NEAR(summary.edge_flows[2], 0.0, 1e-9) << "Path 2 (cost=2) unused";
  EXPECT_NEAR(summary.edge_flows[4], 0.0, 1e-9) << "Path 3 (cost=3) unused";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 4);
}

//=============================================================================
// SECTION 5: PLACEMENT MODE COMPARISON
//
// Direct comparison of PROPORTIONAL vs EQUAL_BALANCED behavior under
// various conditions.
//=============================================================================

TEST(MaxFlow, PlacementModes_EqualCapacity_IdenticalResult) {
  // With equal capacities, both modes should give same total flow
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

  // Both should saturate both paths
  EXPECT_NEAR(total_prop, 20.0, 1e-9);
  EXPECT_NEAR(total_eq, 20.0, 1e-9);
  EXPECT_NEAR(total_prop, total_eq, 1e-9)
      << "With equal capacities, both modes give same total";
}

//=============================================================================
// SECTION 6: FLOW PROPERTIES & INVARIANTS
//
// Validates fundamental flow properties: conservation, capacity constraints,
// min-cut correctness, cost distribution consistency.
//=============================================================================

TEST(MaxFlow, FlowConservation_IntermediateNodes) {
  // Validates flow conservation: inflow = outflow at intermediate nodes
  auto g = make_square_graph(1);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  EXPECT_EQ(summary.edge_flows.size(), static_cast<std::size_t>(g.num_edges()));

  // Use helper to validate conservation at all intermediate nodes
  validate_flow_conservation(g, summary, 0, 2);
}

TEST(MaxFlow, MinCut_ValidEdgeIdentification) {
  // Min-cut should separate source from sink with valid edge IDs
  auto g = make_square_graph(1);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  // Min cut should contain edges
  EXPECT_GT(summary.min_cut.edges.size(), 0)
      << "Min-cut must separate source from sink";

  // All cut edges should have valid IDs
  for (auto eid : summary.min_cut.edges) {
    EXPECT_GE(eid, 0) << "Edge ID must be non-negative";
    EXPECT_LT(eid, g.num_edges()) << "Edge ID must be in valid range";
  }
}

TEST(MaxFlow, CostDistribution_AccountsForAllFlow) {
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
      << "Cost distribution must account for all flow";
}

//=============================================================================
// SECTION 7: ADVANCED FEATURES
//=============================================================================

TEST(MaxFlow, OptionalOutputs_ResidualAndReachable) {
  // Validates optional residual capacity and reachable nodes outputs
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

  // Verify output dimensions
  EXPECT_EQ(summary.residual_capacity.size(),
            static_cast<std::size_t>(g.num_edges()))
      << "Residual capacity vector size must match edge count";
  EXPECT_EQ(summary.reachable_nodes.size(),
            static_cast<std::size_t>(g.num_nodes()))
      << "Reachable nodes vector size must match node count";
}

TEST(MaxFlow, NodeMask_BlocksFlowThroughNode) {
  // Validates that node masking prevents flow through masked nodes
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

TEST(MaxFlow, EdgeMask_SelectiveEdgeFiltering) {
  // Validates edge mask functionality
  auto g = make_n_disjoint_paths(2, 10.0);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  // Create mask with all edges enabled
  auto edge_mask_vec = make_bool_mask(g.num_edges(), true);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.with_edge_flows = true;
  opts.edge_mask = std::span<const bool>(edge_mask_vec.get(),
                                        static_cast<std::size_t>(g.num_edges()));

  auto [total, summary] = algs.max_flow(gh, 0, 3, opts);

  // With all edges enabled, should use both paths
  EXPECT_NEAR(total, 20.0, 1e-9)
      << "With all edges enabled, both paths should be used";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

TEST(MaxFlow, BatchMaxFlow_IndependentMaskApplication) {
  // Validates batch max-flow with different masks per pair
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
      << "First pair (unmasked) should find flow";
  EXPECT_NEAR(summaries[1].total_flow, 0.0, 1e-12)
      << "Second pair (masked) should find no flow";
}

//=============================================================================
// SECTION 8: COMPLEX TOPOLOGIES
//
// Validates algorithm correctness on realistic network topologies beyond
// simple test graphs.
//=============================================================================

TEST(MaxFlow, GridGraph_3x3_MultipleEqualCostPaths) {
  // Grid graphs have many equal-cost paths between corners
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
  EXPECT_GT(total, 0.0) << "Should find flow in grid topology";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 8);
}

TEST(MaxFlow, SharedBottleneck_MultiplePathsMerge) {
  // Two paths that merge at an intermediate node before reaching sink
  // Validates correct handling of shared edges
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
  // Total limited by bottleneck + other path = 5 + 10 = 15
  EXPECT_NEAR(total, 15.0, 1e-9)
      << "Flow limited by shared bottleneck plus alternate path";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 3);
}

//=============================================================================
// SECTION 10: REAL-WORLD TOPOLOGY VALIDATION (NetGraph Integration Tests)
//=============================================================================
//
// These tests validate that NetGraph-Core produces the same maxflow results
// as the original NetGraph repository for important real-world topologies.
// Reference: NetGraph tests/integration/test_scenario_3.py and test_maxflow_api.py

TEST(MaxFlow, Clos_3Tier_Topology_Proportional_Shortest) {
  /**
   * Validates 3-tier Clos fabric maxflow matches NetGraph scenario_3 expectations.
   *
   * Topology: 2 interconnected 3-tier Clos fabrics
   * Source: Clos1 T1 nodes (0-3, 8-11) = 8 nodes
   * Sink: Clos2 T1 nodes (32-35, 40-43) = 8 nodes
   *
   * Expected: 3200.0 Gbps (from NetGraph integration test expectations.py)
   * Bottleneck: 8 active inter-clos spine links × 400 Gbps = 3200 Gbps
   */
  auto g = make_3tier_clos_graph();
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  // Test with PROPORTIONAL placement (matches NetGraph scenario_3 "capacity_analysis_forward")
  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;  // shortest_path: true in scenario_3.yaml
  opts.with_edge_flows = true;

  // Compute maxflow from all Clos1 T1 nodes to all Clos2 T1 nodes
  // We'll compute from representative nodes and verify total capacity
  // For simplicity, test one representative path: node 0 (clos1/b1/t1-1) to node 32 (clos2/b1/t1-1)
  auto [total, summary] = algs.max_flow(gh, 0, 32, opts);

  // Expected: Each T1 node can push 400 Gbps through the network
  // (4 T2 paths × 100 Gbps = 400 Gbps from T1 to T2 layer)
  // But bottleneck at inter-clos spines limits per-T1 flow
  EXPECT_NEAR(total, 400.0, 1.0)
      << "Single T1-to-T1 maxflow should be ~400 Gbps";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 32);
}

TEST(MaxFlow, Clos_3Tier_Topology_EqualBalanced_Shortest) {
  /**
   * Validates 3-tier Clos with EQUAL_BALANCED placement.
   * Reference: NetGraph scenario_3 "capacity_analysis_forward_balanced"
   */
  auto g = make_3tier_clos_graph();
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::EqualBalanced;
  opts.shortest_path = true;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 32, opts);

  // With EQUAL_BALANCED, flow is split equally across paths
  // Expected similar capacity as PROPORTIONAL for this balanced topology
  EXPECT_GT(total, 0.0) << "Should find positive flow";
  EXPECT_LE(total, 400.0 + 1.0) << "Should not exceed topology capacity";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 32);
}

TEST(MaxFlow, Triangle_Topology_Combine_Mode) {
  /**
   * Triangle topology from NetGraph solver tests.
   * Reference: test_maxflow_api.py::test_max_flow_combine_basic()
   *
   * Expected maxflow A->C: 2.0 (1 direct + 1 via B)
   */
  auto g = make_triangle_topology();
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = false;  // Use all paths
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  EXPECT_NEAR(total, 2.0, 1e-9)
      << "Triangle A->C maxflow should be 2.0 (1 direct + 1 via B)";

  validate_capacity_constraints(g, summary);
  validate_flow_conservation(g, summary, 0, 2);
}

TEST(MaxFlow, SimpleParallelPaths_Topology) {
  /**
   * Two disjoint parallel paths topology.
   * Reference: test_maxflow_api.py::test_max_flow_combine_basic() using _simple_network()
   *
   * Expected maxflow S->T: 2.0 (both paths saturated)
   */
  auto g = make_simple_parallel_paths_topology();
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  // Test with shortest_path=True (should use both equal-cost paths)
  MaxFlowOptions opts_sp;
  opts_sp.placement = FlowPlacement::Proportional;
  opts_sp.shortest_path = true;
  opts_sp.with_edge_flows = true;

  auto [total_sp, summary_sp] = algs.max_flow(gh, 0, 3, opts_sp);

  EXPECT_NEAR(total_sp, 2.0, 1e-9)
      << "Should saturate both equal-cost parallel paths with shortest_path=True";

  validate_capacity_constraints(g, summary_sp);
  validate_flow_conservation(g, summary_sp, 0, 3);

  // Test with shortest_path=False (should also be 2.0 since all paths are equal cost)
  MaxFlowOptions opts_full;
  opts_full.placement = FlowPlacement::Proportional;
  opts_full.shortest_path = false;
  opts_full.with_edge_flows = true;

  auto [total_full, summary_full] = algs.max_flow(gh, 0, 3, opts_full);

  EXPECT_NEAR(total_full, 2.0, 1e-9)
      << "Should saturate both paths with shortest_path=False";

  validate_capacity_constraints(g, summary_full);
  validate_flow_conservation(g, summary_full, 0, 3);
}

TEST(MaxFlow, ShortestPath_MustSaturateAllEqualCostPaths) {
  /**
   * Validates that shortest_path=True correctly saturates ALL equal-cost paths,
   * not just a single path.
   * Reference: test_maxflow_api.py::test_shortest_path_vs_full_max_flow()
   *
   * When multiple parallel paths have equal cost, the algorithm must utilize
   * all of them to maximize throughput, even with shortest_path=True.
   */
  auto g = make_simple_parallel_paths_topology();
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = true;

  auto [total, _] = algs.max_flow(gh, 0, 3, opts);

  // Two equal-cost paths with capacity 1.0 each should yield total flow of 2.0
  EXPECT_NEAR(total, 2.0, 1e-9)
      << "shortest_path=True must saturate all equal-cost paths, not just one";
}

//=============================================================================
// SECTION 11: SENSITIVITY ANALYSIS
//=============================================================================

TEST(MaxFlow, Sensitivity_SinglePath) {
  auto g = make_line_graph(3);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;

  auto results = algs.sensitivity_analysis(gh, 0, 2, opts);

  // Path: 0->1->2. Both edges saturated and critical.
  // Capacity is 1.0. Flow is 1.0.
  // Removing any edge drops flow to 0. Delta = 1.0.
  EXPECT_EQ(results.size(), 2u);
  for (const auto& p : results) {
    EXPECT_NEAR(p.second, 1.0, 1e-9);
  }
}

TEST(MaxFlow, Sensitivity_TwoParallelPaths) {
  // Two paths, capacity 10 each. Total 20.
  auto g = make_n_disjoint_paths(2, 10.0);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;

  auto results = algs.sensitivity_analysis(gh, 0, 3, opts);

  // All 4 edges are saturated.
  // Removing any edge kills one path (flow 20 -> 10). Delta = 10.
  EXPECT_EQ(results.size(), 4u);
  for (const auto& p : results) {
    EXPECT_NEAR(p.second, 10.0, 1e-9);
  }
}

TEST(MaxFlow, Sensitivity_PartialSaturation) {
  // Topology: S->A (cap 10), S->B (cap 5)
  //           A->T (cap 5),  B->T (cap 10)
  // Flow: S->A->T limited by A->T (5). S->B->T limited by S->B (5). Total 10.
  // Saturated: A->T (edge 2), S->B (edge 1).
  // Non-Saturated: S->A (edge 0), B->T (edge 3).

  std::int32_t src[4] = {0, 0, 1, 2};
  std::int32_t dst[4] = {1, 2, 3, 3};
  double cap[4] = {10.0, 5.0, 5.0, 10.0};
  std::int64_t cost[4] = {1, 1, 1, 1};

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src, 4), std::span(dst, 4),
    std::span(cap, 4), std::span(cost, 4));

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = false;

  auto results = algs.sensitivity_analysis(gh, 0, 3, opts);

  // Should only report saturated edges: S->B and A->T.
  // (We assume edge indices follow insertion order: 0, 1, 2, 3)
  // S->B is index 1. A->T is index 2.

  // Check results
  int count = 0;
  for (const auto& p : results) {
    if (p.first == 1 || p.first == 2) {
      EXPECT_NEAR(p.second, 5.0, 1e-9);
      count++;
    } else {
      // If non-saturated edges are reported (which they shouldn't be based on logic), fail.
      FAIL() << "Reported sensitivity for non-saturated edge " << p.first;
    }
  }
  EXPECT_EQ(count, 2);
}

TEST(MaxFlow, Sensitivity_ShortestPathVsMaxFlow) {
  // Tests that shortest_path mode produces different sensitivity results
  // than full max-flow mode when there are paths of different costs.
  //
  // Topology: S(0) -> A(1) -> T(2) [cost 1+1=2, cap 10 each]
  //           S(0) -> B(3) -> T(2) [cost 2+2=4, cap 5 each]
  //
  // With shortest_path=false (full max-flow):
  //   - Uses both paths: S->A->T (10) + S->B->T (5) = 15 total
  //   - All 4 edges are saturated and critical
  //
  // With shortest_path=true (single-pass, IP/IGP mode):
  //   - Only uses cheapest path: S->A->T (10)
  //   - Edges 2,3 (S->B->T path) are NOT used, so NOT critical

  std::int32_t src[4] = {0, 1, 0, 3};
  std::int32_t dst[4] = {1, 2, 3, 2};
  double cap[4] = {10.0, 10.0, 5.0, 5.0};
  std::int64_t cost[4] = {1, 1, 2, 2};  // S->A->T costs 2, S->B->T costs 4

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src, 4), std::span(dst, 4),
    std::span(cap, 4), std::span(cost, 4));

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  // Step 1: Verify baseline flow values with max_flow
  // This validates the routing semantics before testing sensitivity
  {
    MaxFlowOptions opts;
    opts.placement = FlowPlacement::Proportional;
    opts.shortest_path = false;
    auto [flow_full, _] = algs.max_flow(gh, 0, 2, opts);
    EXPECT_NEAR(flow_full, 15.0, 1e-9) << "Full max-flow should be 15";
  }
  {
    MaxFlowOptions opts;
    opts.placement = FlowPlacement::Proportional;
    opts.shortest_path = true;
    auto [flow_sp, _] = algs.max_flow(gh, 0, 2, opts);
    EXPECT_NEAR(flow_sp, 10.0, 1e-9) << "Shortest-path flow should be 10";
  }

  // Step 2: Test sensitivity with shortest_path=false (full max-flow)
  {
    MaxFlowOptions opts;
    opts.placement = FlowPlacement::Proportional;
    opts.shortest_path = false;

    auto results = algs.sensitivity_analysis(gh, 0, 2, opts);

    // All 4 edges should be saturated and critical
    EXPECT_EQ(results.size(), 4u) << "Full max-flow should report all 4 edges as critical";

    // Edges 0,1 (S->A->T path) have delta 10 when removed
    // Edges 2,3 (S->B->T path) have delta 5 when removed
    for (const auto& p : results) {
      if (p.first == 0 || p.first == 1) {
        EXPECT_NEAR(p.second, 10.0, 1e-9);
      } else {
        EXPECT_NEAR(p.second, 5.0, 1e-9);
      }
    }
  }

  // Step 3: Test sensitivity with shortest_path=true (IP/IGP mode)
  {
    MaxFlowOptions opts;
    opts.placement = FlowPlacement::Proportional;
    opts.shortest_path = true;

    auto results = algs.sensitivity_analysis(gh, 0, 2, opts);

    // Only edges 0,1 (S->A->T path) should be critical
    // Edges 2,3 (S->B->T path) are not used under shortest-path routing
    EXPECT_EQ(results.size(), 2u) << "Shortest-path mode should only report 2 edges as critical";

    for (const auto& p : results) {
      EXPECT_TRUE(p.first == 0 || p.first == 1)
          << "Edge " << p.first << " should not be critical in shortest-path mode";
      // Baseline=10, removing either S->A or A->T forces traffic to S->B->T (cap 5)
      // Delta = 10 - 5 = 5
      EXPECT_NEAR(p.second, 5.0, 1e-9);
    }
  }
}
