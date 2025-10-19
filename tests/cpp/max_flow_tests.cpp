#include <gtest/gtest.h>
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/backend.hpp"
#include "netgraph/core/algorithms.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

TEST(MaxFlow, SimpleLinearPath) {
  auto g = make_line_graph(3);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.shortest_path = false;
  opts.with_edge_flows = false;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  // Single path with capacity 1.0
  EXPECT_NEAR(total, 1.0, 1e-9);
}

TEST(MaxFlow, MultipleDisjointPaths) {
  // Create two disjoint paths: 0->1->3 and 0->2->3
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

  // Total flow is sum of both paths: 1 + 2 = 3
  EXPECT_NEAR(total, 3.0, 1e-9);
}

TEST(MaxFlow, BottleneckIdentification) {
  // Path with bottleneck: 0->1 (cap 10), 1->2 (cap 2), 2->3 (cap 10)
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

  // Bottleneck is edge 1->2 with capacity 2
  EXPECT_NEAR(total, 2.0, 1e-9);
}

TEST(MaxFlow, MinCutCorrectness) {
  auto g = make_square_graph(1);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  // Min cut should separate source from sink
  EXPECT_GT(summary.min_cut.edges.size(), 0);

  // All cut edges should be valid
  for (auto eid : summary.min_cut.edges) {
    EXPECT_GE(eid, 0);
    EXPECT_LT(eid, g.num_edges());
  }
}

TEST(MaxFlow, EdgeFlowConservation) {
  auto g = make_square_graph(1);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  EXPECT_EQ(summary.edge_flows.size(), static_cast<std::size_t>(g.num_edges()));

  // Flow conservation: inflow = outflow for intermediate nodes
  auto row = g.row_offsets_view();
  auto col = g.col_indices_view();
  auto aei = g.adj_edge_index_view();
  auto in_row = g.in_row_offsets_view();
  auto in_col = g.in_col_indices_view();
  auto in_aei = g.in_adj_edge_index_view();

  for (std::int32_t u = 0; u < g.num_nodes(); ++u) {
    if (u == 0 || u == 2) continue;  // Skip source and sink

    double outflow = 0.0;
    for (std::size_t j = row[u]; j < row[u + 1]; ++j) {
      outflow += summary.edge_flows[aei[j]];
    }

    double inflow = 0.0;
    for (std::size_t j = in_row[u]; j < in_row[u + 1]; ++j) {
      inflow += summary.edge_flows[in_aei[j]];
    }

    EXPECT_NEAR(inflow, outflow, 1e-9) << "Flow not conserved at node " << u;
  }
}

TEST(MaxFlow, ProportionalVsEqualBalanced) {
  auto g = make_square_graph(1);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts1;
  opts1.placement = FlowPlacement::Proportional;
  auto [total1, summary1] = algs.max_flow(gh, 0, 2, opts1);

  MaxFlowOptions opts2;
  opts2.placement = FlowPlacement::EqualBalanced;
  auto [total2, summary2] = algs.max_flow(gh, 0, 2, opts2);

  // Both should find the same max flow value
  EXPECT_NEAR(total1, total2, 1e-9);
}

TEST(MaxFlow, ShortestPathOnlyMode) {
  auto g = make_square_graph(1);

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

  // Shortest path only should be less than or equal to full max flow
  EXPECT_LE(total_sp, total_full + 1e-9);
}

TEST(MaxFlow, CostDistributionAccuracy) {
  auto g = make_square_graph(1);

  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.with_edge_flows = true;

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);

  // Cost distribution should sum to total flow
  double cost_sum = 0.0;
  for (auto flow : summary.flows) {
    cost_sum += flow;
  }

  EXPECT_NEAR(cost_sum, total, 1e-9);
}

TEST(MaxFlow, ReachableAndResidualsFlags) {
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
  EXPECT_EQ(summary.residual_capacity.size(), static_cast<std::size_t>(g.num_edges()));
  EXPECT_EQ(summary.reachable_nodes.size(), static_cast<std::size_t>(g.num_nodes()));
}

TEST(MaxFlow, NodeMaskBlocksFlow) {
  // Verify that masking the middle node in a line graph prevents all flow
  auto g = make_line_graph(3);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  // Mask out node 1 (middle node)
  auto node_mask_vec = make_bool_mask(g.num_nodes());
  node_mask_vec[1] = false;

  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;
  opts.node_mask = std::span<const bool>(node_mask_vec.get(), static_cast<std::size_t>(g.num_nodes()));

  auto [total, summary] = algs.max_flow(gh, 0, 2, opts);
  EXPECT_NEAR(total, 0.0, 1e-12);
}

TEST(MaxFlow, BatchMaxFlowRespectsMasks) {
  // Verify that batch max-flow applies different masks to different pairs independently
  auto g = make_line_graph(3);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);

  std::vector<std::pair<NodeId,NodeId>> pairs = {{0,2}, {0,2}};
  MaxFlowOptions opts;
  opts.placement = FlowPlacement::Proportional;

  // First pair: no mask (should find flow)
  // Second pair: mask node 1 (should find no flow)
  std::vector<std::span<const bool>> node_masks;
  auto nm0 = make_bool_mask(g.num_nodes());  // all true
  auto nm1 = make_bool_mask(g.num_nodes());
  nm1[1] = false; // block middle node for second pair
  node_masks.push_back(std::span<const bool>(nm0.get(), static_cast<std::size_t>(g.num_nodes())));
  node_masks.push_back(std::span<const bool>(nm1.get(), static_cast<std::size_t>(g.num_nodes())));

  auto summaries = algs.batch_max_flow(gh, pairs, opts, node_masks, {});
  ASSERT_EQ(summaries.size(), 2u);
  EXPECT_GT(summaries[0].total_flow, 0.0);
  EXPECT_NEAR(summaries[1].total_flow, 0.0, 1e-12);
}
