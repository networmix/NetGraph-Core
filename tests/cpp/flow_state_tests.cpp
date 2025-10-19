#include <gtest/gtest.h>
#include "netgraph/core/flow_state.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/constants.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

TEST(FlowState, InitialResidualEqualsCapacity) {
  auto g = make_line_graph(3);
  FlowState fs(g);

  auto residual = fs.residual_view();
  auto capacity = g.capacity_view();

  for (std::size_t i = 0; i < residual.size(); ++i) {
    EXPECT_DOUBLE_EQ(residual[i], capacity[i]);
  }
}

TEST(FlowState, PlacementReducesResidual) {
  auto g = make_line_graph(3);
  FlowState fs(g);

  // Get shortest path DAG
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  auto initial_residual = fs.residual_view()[0];

  // Place 0.5 units of flow
  Flow placed = fs.place_on_dag(0, 2, dag, 0.5, FlowPlacement::Proportional, false);

  EXPECT_GT(placed, 0.0);

  // Residual should be reduced
  auto new_residual = fs.residual_view()[0];
  EXPECT_LT(new_residual, initial_residual);
}

TEST(FlowState, ProportionalDistribution) {
  // Create graph with two parallel edges of different capacities
  std::int32_t src[3] = {0, 1, 1};
  std::int32_t dst[3] = {1, 2, 2};
  double cap[3] = {10.0, 3.0, 7.0};  // Total 10 capacity B->C
  std::int64_t cost[3] = {1, 1, 1};

  auto g = StrictMultiDiGraph::from_arrays(3,
    std::span(src, 3), std::span(dst, 3),
    std::span(cap, 3), std::span(cost, 3));

  FlowState fs(g);

  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  // Place 10 units (should use both B->C edges proportionally)
  Flow placed = fs.place_on_dag(0, 2, dag, 10.0, FlowPlacement::Proportional, false);

  EXPECT_NEAR(placed, 10.0, 1e-9);

  auto flows = fs.edge_flow_view();
  // First edge (A->B) should have 10
  EXPECT_NEAR(flows[0], 10.0, 1e-9);

  // Parallel edges should have flow proportional to capacity: 3 and 7
  EXPECT_NEAR(flows[1], 3.0, 1e-9);
  EXPECT_NEAR(flows[2], 7.0, 1e-9);
}

TEST(FlowState, EqualBalancedDistribution) {
  // Same graph as proportional test
  std::int32_t src[3] = {0, 1, 1};
  std::int32_t dst[3] = {1, 2, 2};
  double cap[3] = {10.0, 3.0, 7.0};
  std::int64_t cost[3] = {1, 1, 1};

  auto g = StrictMultiDiGraph::from_arrays(3,
    std::span(src, 3), std::span(dst, 3),
    std::span(cap, 3), std::span(cost, 3));

  FlowState fs(g);

  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  // Place 6 units with equal-balanced (limited by min capacity * 2 = 6)
  Flow placed = fs.place_on_dag(0, 2, dag, 10.0, FlowPlacement::EqualBalanced, false);

  EXPECT_NEAR(placed, 6.0, 1e-9);

  auto flows = fs.edge_flow_view();
  // Parallel edges should each get 3 (equal split)
  EXPECT_NEAR(flows[1], 3.0, 1e-9);
  EXPECT_NEAR(flows[2], 3.0, 1e-9);
}

TEST(FlowState, ZeroCapacityHandling) {
  std::int32_t src[1] = {0};
  std::int32_t dst[1] = {1};
  double cap[1] = {0.0};
  std::int64_t cost[1] = {1};

  auto g = StrictMultiDiGraph::from_arrays(2,
    std::span(src, 1), std::span(dst, 1),
    std::span(cap, 1), std::span(cost, 1));

  FlowState fs(g);

  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = true;  // Should filter zero capacity
  sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, fs.residual_view(), {}, {});

  // No path should be found due to zero capacity
  EXPECT_EQ(dag.parent_offsets[1], dag.parent_offsets[2]);
}

TEST(FlowState, ApplyDeltasAdditive) {
  auto g = make_line_graph(3);
  FlowState fs(g);

  std::vector<std::pair<EdgeId, Flow>> deltas = {{0, 0.5}, {1, 0.3}};

  auto initial_flow_0 = fs.edge_flow_view()[0];
  auto initial_residual_0 = fs.residual_view()[0];

  fs.apply_deltas(deltas, true);  // Add

  EXPECT_NEAR(fs.edge_flow_view()[0], initial_flow_0 + 0.5, 1e-9);
  EXPECT_NEAR(fs.residual_view()[0], initial_residual_0 - 0.5, 1e-9);
}

TEST(FlowState, ApplyDeltasSubtractive) {
  auto g = make_line_graph(3);
  FlowState fs(g);

  // First add some flow
  std::vector<std::pair<EdgeId, Flow>> deltas = {{0, 0.5}};
  fs.apply_deltas(deltas, true);

  auto flow_after_add = fs.edge_flow_view()[0];
  auto residual_after_add = fs.residual_view()[0];

  // Now subtract
  fs.apply_deltas(deltas, false);

  EXPECT_NEAR(fs.edge_flow_view()[0], flow_after_add - 0.5, 1e-9);
  EXPECT_NEAR(fs.residual_view()[0], residual_after_add + 0.5, 1e-9);
}

TEST(FlowState, MinCutComputation) {
  // Test min-cut computation after saturating max-flow on square graph.
  // Square has two paths: 0->1->2 (cap 1, cost 2) and 0->3->2 (cap 2, cost 4).
  // Max flow is 3.0, and min-cut should separate source from sink with saturated edges.
  auto g = make_square_graph(1);
  FlowState fs(g);

  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = true;
  sel.tie_break = EdgeTieBreak::Deterministic;

  // Saturate max-flow: iteratively place flow until exhausted
  Flow total_placed = 0.0;
  for (int iter = 0; iter < 10; ++iter) {
    auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, fs.residual_view(), {}, {});
    // Check if destination is reachable
    if (dag.parent_offsets[2] == dag.parent_offsets[3]) break; // No path to node 2
    Flow placed = fs.place_on_dag(0, 2, dag, 10.0, FlowPlacement::Proportional, false);
    if (placed < kMinFlow) break;
    total_placed += placed;
  }

  // Verify flow was actually placed
  ASSERT_GT(total_placed, 0.0) << "No flow placed, test setup is broken";

  // Compute min-cut from source
  auto min_cut = fs.compute_min_cut(0, {}, {});

  // Min-cut must have at least one edge (separating source from sink)
  EXPECT_GT(min_cut.edges.size(), 0) << "Min-cut should contain edges after flow saturation";

  // All cut edges should be saturated (residual <= kMinCap)
  auto residual = fs.residual_view();
  for (auto eid : min_cut.edges) {
    EXPECT_LE(residual[eid], kMinCap)
      << "Cut edge " << eid << " not saturated (residual=" << residual[eid] << ")";
  }

  // Verify min-cut capacity equals max-flow (by max-flow min-cut theorem)
  double cut_capacity = 0.0;
  auto capacity = g.capacity_view();
  for (auto eid : min_cut.edges) {
    cut_capacity += capacity[eid];
  }
  EXPECT_NEAR(cut_capacity, total_placed, 1e-9)
    << "Min-cut capacity should equal max-flow by max-flow min-cut theorem";
}

TEST(FlowState, ResetRestoresCapacity) {
  auto g = make_line_graph(3);
  FlowState fs(g);

  // Place some flow
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});
  Flow placed = fs.place_on_dag(0, 2, dag, 0.5, FlowPlacement::Proportional, false);
  EXPECT_GT(placed, 0.0);

  EXPECT_GT(fs.edge_flow_view()[0], 0.0);

  // Reset
  fs.reset();

  // Flow should be zero
  EXPECT_DOUBLE_EQ(fs.edge_flow_view()[0], 0.0);

  // Residual should equal capacity
  auto residual = fs.residual_view();
  auto capacity = g.capacity_view();
  for (std::size_t i = 0; i < residual.size(); ++i) {
    EXPECT_DOUBLE_EQ(residual[i], capacity[i]);
  }
}

TEST(FlowState, CustomInitialResidual) {
  auto g = make_line_graph(3);

  // Create custom residual (half capacity)
  std::vector<Cap> custom_residual(g.num_edges());
  auto capacity = g.capacity_view();
  for (std::size_t i = 0; i < custom_residual.size(); ++i) {
    custom_residual[i] = capacity[i] * 0.5;
  }

  FlowState fs(g, custom_residual);

  auto residual = fs.residual_view();
  for (std::size_t i = 0; i < residual.size(); ++i) {
    EXPECT_DOUBLE_EQ(residual[i], custom_residual[i]);
  }
}

TEST(FlowState, InvalidResidualInitThrows) {
  auto g = make_line_graph(3);
  // residual_init wrong length
  std::vector<Cap> bad(static_cast<std::size_t>(g.num_edges() + 1), 1.0);
  EXPECT_THROW({ FlowState fs_bad(g, bad); (void)fs_bad; }, std::invalid_argument);
}

TEST(FlowState, ResetWithInvalidResidualThrows) {
  auto g = make_line_graph(3);
  FlowState fs(g);
  std::vector<Cap> bad(static_cast<std::size_t>(g.num_edges() + 2), 1.0);
  EXPECT_THROW({ fs.reset(bad); }, std::invalid_argument);
}
