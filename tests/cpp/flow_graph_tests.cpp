#include <gtest/gtest.h>
#include "netgraph/core/flow_graph.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

TEST(FlowGraph, EmptyLedgerInitially) {
  auto g = make_line_graph(3);
  FlowGraph fg(g);

  FlowIndex idx{0, 2, 0, 0};
  auto edges = fg.get_flow_edges(idx);

  EXPECT_TRUE(edges.empty());
}

TEST(FlowGraph, PlacementCreatesEntry) {
  auto g = make_line_graph(3);
  FlowGraph fg(g);

  // Get shortest path
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  FlowIndex idx{0, 2, 0, 0};
  Flow placed = fg.place(idx, 0, 2, dag, 0.5, FlowPlacement::Proportional);

  EXPECT_GT(placed, 0.0);

  auto edges = fg.get_flow_edges(idx);
  EXPECT_FALSE(edges.empty());
}

TEST(FlowGraph, RemovalDeletesEntry) {
  auto g = make_line_graph(3);
  FlowGraph fg(g);

  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  FlowIndex idx{0, 2, 0, 0};
  Flow placed = fg.place(idx, 0, 2, dag, 0.5, FlowPlacement::Proportional);
  EXPECT_GT(placed, 0.0);

  EXPECT_FALSE(fg.get_flow_edges(idx).empty());

  fg.remove(idx);

  EXPECT_TRUE(fg.get_flow_edges(idx).empty());
}

TEST(FlowGraph, RemoveByClassFiltersCorrectly) {
  auto g = make_line_graph(3);
  FlowGraph fg(g);

  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  // Create flows with different classes
  FlowIndex idx1{0, 2, 10, 0};
  FlowIndex idx2{0, 2, 20, 1};
  FlowIndex idx3{0, 2, 10, 2};

  (void)fg.place(idx1, 0, 2, dag, 0.5, FlowPlacement::Proportional);
  (void)fg.place(idx2, 0, 2, dag, 0.5, FlowPlacement::Proportional);
  (void)fg.place(idx3, 0, 2, dag, 0.5, FlowPlacement::Proportional);

  // Remove flows of class 10
  fg.remove_by_class(10);

  // idx1 and idx3 should be removed (class 10)
  EXPECT_TRUE(fg.get_flow_edges(idx1).empty());
  EXPECT_TRUE(fg.get_flow_edges(idx3).empty());

  // idx2 should remain (class 20)
  EXPECT_FALSE(fg.get_flow_edges(idx2).empty());
}

TEST(FlowGraph, GetFlowEdgesReturnsDeltas) {
  auto g = make_line_graph(3);
  FlowGraph fg(g);

  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  FlowIndex idx{0, 2, 0, 0};
  Flow placed = fg.place(idx, 0, 2, dag, 0.5, FlowPlacement::Proportional);

  auto edges = fg.get_flow_edges(idx);

  // Should have entries for edges used in the path
  EXPECT_GT(edges.size(), 0);

  // Total flow across all edges should match placed amount times path length
  double total = 0.0;
  for (const auto& [eid, flow] : edges) {
    EXPECT_GE(flow, 0.0);
    total += flow;
  }
  EXPECT_NEAR(total, placed * 2, 1e-9);  // Two edges in path
}

TEST(FlowGraph, GetFlowPathForLinearFlow) {
  auto g = make_line_graph(3);
  FlowGraph fg(g);

  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  FlowIndex idx{0, 2, 0, 0};
  Flow placed = fg.place(idx, 0, 2, dag, 0.5, FlowPlacement::Proportional);
  EXPECT_GT(placed, 0.0);

  auto path = fg.get_flow_path(idx);

  // Should find a simple path with 2 edges
  EXPECT_EQ(path.size(), 2);
}

TEST(FlowGraph, GetFlowPathForDAGReturnsEmpty) {
  // Create a graph with two disjoint equal-cost paths: 0->1->2 and 0->3->2
  // When flow splits proportionally across both paths, get_flow_path should
  // return empty (since there's no single simple path)
  std::int32_t src[4] = {0, 1, 0, 3};
  std::int32_t dst[4] = {1, 2, 3, 2};
  double cap[4] = {1.0, 1.0, 1.0, 1.0};
  std::int64_t cost[4] = {1, 1, 1, 1};

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src, 4), std::span(dst, 4),
    std::span(cap, 4), std::span(cost, 4));

  FlowGraph fg(g);

  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  FlowIndex idx{0, 2, 0, 0};
  // Place 2.0 units - should split across both paths (1.0 each)
  Flow placed = fg.place(idx, 0, 2, dag, 2.0, FlowPlacement::Proportional);
  EXPECT_GT(placed, 0.0);

  auto path = fg.get_flow_path(idx);

  // get_flow_path should return empty when flow is split across multiple paths
  // (it only returns a simple path if one exists)
  EXPECT_TRUE(path.empty()) << "Expected empty path for split flow, got " << path.size() << " edges";
}

TEST(FlowGraph, CoalescesDuplicateEdgeIds) {
  // Create graph with parallel edges
  std::int32_t src[3] = {0, 1, 1};
  std::int32_t dst[3] = {1, 2, 2};
  double cap[3] = {10.0, 1.0, 1.0};
  std::int64_t cost[3] = {1, 1, 1};

  auto g = StrictMultiDiGraph::from_arrays(3,
    std::span(src, 3), std::span(dst, 3),
    std::span(cap, 3), std::span(cost, 3));

  FlowGraph fg(g);

  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  FlowIndex idx{0, 2, 0, 0};
  (void)fg.place(idx, 0, 2, dag, 2.0, FlowPlacement::Proportional);

  auto edges = fg.get_flow_edges(idx);

  // Check that no EdgeId appears more than once
  std::set<EdgeId> seen;
  for (const auto& [eid, flow] : edges) {
    EXPECT_EQ(seen.count(eid), 0) << "Duplicate EdgeId " << eid;
    seen.insert(eid);
  }
}

TEST(FlowGraph, LedgerHistory_RemovalRestoresResidual) {
  // edges: 0->1, 1->2 (both cap=1.0)
  auto g = make_line_graph(3);
  FlowGraph fg(g);

  // Build DAG once
  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
  auto sp = shortest_paths(g, 0, 2, /*multipath=*/true, sel, {}, {}, {});
  const PredDAG& dag = sp.second;

  FlowIndex idx{0, 2, 42, 7};

  // Two placements on same flow
  Flow p1 = fg.place(idx, 0, 2, dag, 0.4, FlowPlacement::Proportional);
  Flow p2 = fg.place(idx, 0, 2, dag, 0.5, FlowPlacement::Proportional);
  ASSERT_GT(p1, 0.0);
  ASSERT_GT(p2, 0.0);

  // Remove the flow (should restore both p1 and p2 if ledger retained history)
  fg.remove(idx);

  auto cap = g.capacity_view();
  auto res = fg.residual_view();
  ASSERT_EQ(res.size(), static_cast<std::size_t>(g.num_edges()));

  // Expect residual to be fully restored to capacity if removal undid entire history
  for (std::size_t i = 0; i < res.size(); ++i) {
    EXPECT_NEAR(res[i], cap[i], 1e-9) << "Residual not fully restored at edge " << i;
  }
}

TEST(FlowGraph, LedgerMicroFlows_RemovalRestoresResidual) {
  // Split 1.0 unit equally over many parallel edges so per-edge share is < kMinFlow.
  const int parallel = 5000; // 1.0 / 5000 < kMinFlow (1/4096)
  std::vector<std::int32_t> src(parallel), dst(parallel);
  std::vector<double> cap(parallel, 1.0);
  std::vector<std::int64_t> cost(parallel, 1);
  for (int i = 0; i < parallel; ++i) { src[i] = 0; dst[i] = 1; }
  auto g = StrictMultiDiGraph::from_arrays(2, src, dst, cap, cost);

  FlowGraph fg(g);
  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
  auto sp = shortest_paths(g, 0, 1, /*multipath=*/true, sel, {}, {}, {});
  const PredDAG& dag = sp.second;

  FlowIndex idx{0, 1, 1, 1};

  // Place 1.0 evenly across many parallels; per-edge contribution < kMinFlow
  Flow placed = fg.place(idx, 0, 1, dag, 1.0, FlowPlacement::EqualBalanced);
  ASSERT_NEAR(placed, 1.0, 1e-9);

  // Ledger must record all deltas, including micro-flows smaller than kMinFlow
  auto edges = fg.get_flow_edges(idx);
  EXPECT_FALSE(edges.empty()) << "Ledger should record micro-deltas even when per-edge < kMinFlow";

  // Removing the flow should fully restore residual capacity with no drift
  auto capv = g.capacity_view();
  fg.remove(idx);
  auto res_after = fg.residual_view();
  ASSERT_EQ(res_after.size(), static_cast<std::size_t>(g.num_edges()));
  for (std::size_t i = 0; i < res_after.size(); ++i) {
    EXPECT_NEAR(res_after[i], capv[i], 1e-9) << "Residual not restored at edge " << i;
  }
}
