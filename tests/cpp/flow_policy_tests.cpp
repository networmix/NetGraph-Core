#include <gtest/gtest.h>

#include "netgraph/core/flow_graph.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/backend.hpp"
#include "netgraph/core/algorithms.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

using namespace netgraph::core;

namespace {

StrictMultiDiGraph make_square1() {
  // Nodes: 0:A, 1:B, 2:C, 3:D
  // Edges:
  // 0: A->B cap=1 cost=1
  // 1: B->C cap=1 cost=1
  // 2: A->D cap=2 cost=2
  // 3: D->C cap=2 cost=2
  std::int32_t num_nodes = 4;
  std::int32_t src_arr[4] = {0, 1, 0, 3};
  std::int32_t dst_arr[4] = {1, 2, 3, 2};
  double cap_arr[4] = {1.0, 1.0, 2.0, 2.0};
  std::int64_t cost_arr[4] = {1, 1, 2, 2};
  std::int64_t ext_arr[4] = {0, 1, 2, 3};
  std::span<const std::int32_t> src(src_arr, 4);
  std::span<const std::int32_t> dst(dst_arr, 4);
  std::span<const double> cap(cap_arr, 4);
  std::span<const std::int64_t> cost(cost_arr, 4);
  std::span<const std::int64_t> ext(ext_arr, 4);
  return StrictMultiDiGraph::from_arrays(num_nodes, src, dst, cap, cost, ext);
}

StrictMultiDiGraph make_line1() {
  // Nodes: 0:A, 1:B, 2:C
  // Edges forward:
  // A->B cap=5 cost=1
  // B->C cap=1 cost=1
  // B->C cap=3 cost=1
  // B->C cap=7 cost=2
  std::int32_t num_nodes = 3;
  std::int32_t src_arr[4] = {0, 1, 1, 1};
  std::int32_t dst_arr[4] = {1, 2, 2, 2};
  double cap_arr[4] = {5.0, 1.0, 3.0, 7.0};
  std::int64_t cost_arr[4] = {1, 1, 1, 2};
  std::int64_t ext_arr[4] = {0, 1, 2, 3};
  std::span<const std::int32_t> src(src_arr, 4);
  std::span<const std::int32_t> dst(dst_arr, 4);
  std::span<const double> cap(cap_arr, 4);
  std::span<const std::int64_t> cost(cost_arr, 4);
  std::span<const std::int64_t> ext(ext_arr, 4);
  return StrictMultiDiGraph::from_arrays(num_nodes, src, dst, cap, cost, ext);
}

StrictMultiDiGraph make_square3() {
  // Nodes: 0:A, 1:B, 2:C, 3:D
  // A->B 100@1, B->C 125@1, A->D 75@1, D->C 50@1
  std::int32_t num_nodes = 4;
  std::int32_t src_arr[4] = {0, 1, 0, 3};
  std::int32_t dst_arr[4] = {1, 2, 3, 2};
  double cap_arr[4] = {100.0, 125.0, 75.0, 50.0};
  std::int64_t cost_arr[4] = {1, 1, 1, 1};
  std::int64_t ext_arr[4] = {0, 1, 2, 3};
  std::span<const std::int32_t> src(src_arr, 4);
  std::span<const std::int32_t> dst(dst_arr, 4);
  std::span<const double> cap(cap_arr, 4);
  std::span<const std::int64_t> cost(cost_arr, 4);
  std::span<const std::int64_t> ext(ext_arr, 4);
  return StrictMultiDiGraph::from_arrays(num_nodes, src, dst, cap, cost, ext);
}

void expect_edge_flows_by_uv(const FlowGraph& fg, std::initializer_list<std::tuple<int,int,double>> exp) {
  const auto& g = fg.graph();
  auto row = g.row_offsets_view();
  auto col = g.col_indices_view();
  auto aei = g.adj_edge_index_view();
  // Build map EdgeId->flow
  auto ef = fg.edge_flow_view();
  // Check each expected (u,v,flow)
  // When a (u,v) pair appears multiple times, advance through successive edges
  // in the adjacency row to validate each parallel edge in order.
  std::map<std::pair<int,int>, std::size_t> cursor;
  for (auto [u, v, expected] : exp) {
    bool found = false;
    // Scan adjacency row for (u->v)
    auto s = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto e = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
    std::size_t start = s;
    auto key = std::make_pair(u, v);
    if (auto it = cursor.find(key); it != cursor.end()) start = it->second;
    for (std::size_t j = start; j < e; ++j) {
      if (static_cast<int>(col[j]) == v) {
        auto eid = static_cast<std::size_t>(aei[j]);
        EXPECT_NEAR(static_cast<double>(ef[eid]), expected, 1e-9) << "edge (" << u << "," << v << ")";
        cursor[key] = j + 1; // subsequent matches continue after this position
        found = true;
        break;
      }
    }
    ASSERT_TRUE(found) << "edge (" << u << "," << v << ") not found";
  }
}

}

TEST(FlowPolicyCore, Proportional_SingleDemand_UsesShortestPath) {
  auto g = make_square1();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  auto be = make_cpu_backend(); auto algs = std::make_shared<Algorithms>(be); auto gh = algs->build_graph(g);
  ExecutionContext ctx(algs, gh);
  FlowPolicyConfig cfg;
  cfg.flow_placement = FlowPlacement::Proportional;
  cfg.selection = sel;
  cfg.require_capacity = true;
  FlowPolicy policy(ctx, cfg);
  auto res = policy.place_demand(fg, /*src=*/0, /*dst=*/2, /*flowClass=*/0, /*volume=*/1.0);
  EXPECT_NEAR(res.first, 1.0, 1e-9);
  EXPECT_NEAR(res.second, 0.0, 1e-9);
  expect_edge_flows_by_uv(fg, {{0,1,1.0}, {1,2,1.0}, {0,3,0.0}, {3,2,0.0}});
}

TEST(FlowPolicyCore, Proportional_SingleDemand_SplitsAcrossTwoEqualCostPaths) {
  auto g = make_square1();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  auto be = make_cpu_backend(); auto algs = std::make_shared<Algorithms>(be); auto gh = algs->build_graph(g);
  ExecutionContext ctx(algs, gh);
  FlowPolicyConfig cfg;
  cfg.flow_placement = FlowPlacement::Proportional;
  cfg.selection = sel;
  cfg.require_capacity = true;
  FlowPolicy policy(ctx, cfg);
  auto res = policy.place_demand(fg, 0, 2, 0, 2.0);
  EXPECT_NEAR(res.first, 2.0, 1e-9);
  EXPECT_NEAR(res.second, 0.0, 1e-9);
  expect_edge_flows_by_uv(fg, {{0,1,1.0}, {1,2,1.0}, {0,3,1.0}, {3,2,1.0}});
}

TEST(FlowPolicyCore, Proportional_SingleDemand_MaxFlowCount1_UsesSinglePath) {
  auto g = make_square1();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  auto be = make_cpu_backend(); auto algs = std::make_shared<Algorithms>(be); auto gh = algs->build_graph(g);
  ExecutionContext ctx(algs, gh);
  FlowPolicyConfig cfg;
  cfg.flow_placement = FlowPlacement::Proportional;
  cfg.selection = sel;
  cfg.require_capacity = true;
  cfg.multipath = false;
  cfg.min_flow_count = 1;
  cfg.max_flow_count = 1;
  FlowPolicy policy(ctx, cfg);
  auto res = policy.place_demand(fg, 0, 2, 0, 2.0);
  EXPECT_NEAR(res.first, 2.0, 1e-9);
  EXPECT_NEAR(res.second, 0.0, 1e-9);
  expect_edge_flows_by_uv(fg, {{0,1,0.0}, {1,2,0.0}, {0,3,2.0}, {3,2,2.0}});
}

TEST(FlowPolicyCore, Proportional_SingleDemand_PartialPlacement_InsufficientCapacity) {
  auto g = make_square1();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  auto be = make_cpu_backend(); auto algs = std::make_shared<Algorithms>(be); auto gh = algs->build_graph(g);
  ExecutionContext ctx(algs, gh);
  FlowPolicyConfig cfg;
  cfg.flow_placement = FlowPlacement::Proportional;
  cfg.selection = sel;
  cfg.require_capacity = true;
  FlowPolicy policy(ctx, cfg);
  auto res = policy.place_demand(fg, 0, 2, 0, 5.0);
  EXPECT_NEAR(res.first, 3.0, 1e-9);
  EXPECT_NEAR(res.second, 2.0, 1e-9);
  expect_edge_flows_by_uv(fg, {{0,1,1.0}, {1,2,1.0}, {0,3,2.0}, {3,2,2.0}});
}

TEST(FlowPolicyCore, EqualBalanced_SingleDemand_FlowCount2_BalancedDistribution) {
  auto g = make_line1();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multi_edge = false; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  auto be = make_cpu_backend(); auto algs = std::make_shared<Algorithms>(be); auto gh = algs->build_graph(g);
  ExecutionContext ctx(algs, gh);
  FlowPolicyConfig cfg;
  cfg.flow_placement = FlowPlacement::EqualBalanced;
  cfg.selection = sel;
  cfg.require_capacity = true;
  cfg.min_flow_count = 2;
  cfg.max_flow_count = 2;
  FlowPolicy policy(ctx, cfg);
  auto res = policy.place_demand(fg, 0, 2, 0, 7.0);
  EXPECT_NEAR(res.first, 4.0, 1e-9);
  EXPECT_NEAR(res.second, 3.0, 1e-9);
  // A->B: 4, B->C: cap1 gets 1.0, cap3 gets 3.0, cap7(cost2) unused
  expect_edge_flows_by_uv(fg, {{0,1,4.0}, {1,2,1.0}, {1,2,3.0}, {1,2,0.0}});
}

TEST(FlowPolicyCore, EqualBalanced_SingleDemand_FlowCount3_SaturatesBothPaths) {
  // Graph: A->B (cap 100), B->C (cap 125), A->D (cap 75), D->C (cap 50)
  // Two equal-cost paths exist with capacities 100 (A->B->C) and 50 (A->D->C).
  // EqualBalanced with 3 flows should saturate both paths for a total of 150,
  // leaving 50 unplaced out of the 200 requested.
  auto g = make_square3();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multi_edge = false; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  auto be = make_cpu_backend(); auto algs = std::make_shared<Algorithms>(be); auto gh = algs->build_graph(g);
  ExecutionContext ctx(algs, gh);
  FlowPolicyConfig cfg;
  cfg.flow_placement = FlowPlacement::EqualBalanced;
  cfg.selection = sel;
  cfg.require_capacity = true;
  cfg.min_flow_count = 3;
  cfg.max_flow_count = 3;
  FlowPolicy policy(ctx, cfg);
  auto res = policy.place_demand(fg, 0, 2, 0, 200.0);
  EXPECT_NEAR(res.first, 150.0, 1e-9);
  EXPECT_NEAR(res.second, 50.0, 1e-9);
  // Path A->B->C carries 100, A->D->C carries 50
  expect_edge_flows_by_uv(fg, {{0,1,100.0}, {1,2,100.0}, {0,3,50.0}, {3,2,50.0}});
}

TEST(FlowPolicyCore, DiminishingReturnsCutoff) {
  auto g = make_line1();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
  auto be = make_cpu_backend(); auto algs = std::make_shared<Algorithms>(be); auto gh = algs->build_graph(g);
  ExecutionContext ctx(algs, gh);
  FlowPolicyConfig cfg;
  cfg.flow_placement = FlowPlacement::EqualBalanced;
  cfg.selection = sel;
  cfg.require_capacity = true;
  cfg.multipath = true;
  cfg.min_flow_count = 1;
  cfg.max_flow_count = 1000000;
  cfg.shortest_path = false;
  cfg.reoptimize_flows_on_each_placement = false;
  cfg.max_no_progress_iterations = 100;
  cfg.max_total_iterations = 10000;
  cfg.diminishing_returns_enabled = true;
  cfg.diminishing_returns_window = 8;
  cfg.diminishing_returns_epsilon_frac = 1e-3;
  FlowPolicy policy(ctx, cfg);
  auto res = policy.place_demand(fg, 0, 2, 0, 7.0);
  EXPECT_GE(res.first, 0.0);
  EXPECT_GE(res.second, 0.0);
  EXPECT_GT(res.second, 0.0);
}

TEST(FlowPolicyCore, EqualBalanced_ShortestPath_ECMP3_SaturatesAllEqualCostPaths) {
  // Graph: 0->{1,2,3}->4, all caps=5, all costs=1 (ECMP width 3)
  // ECMP-like config: EqualBalanced, multi_edge=true, require_capacity=false, shortest_path=true, single flow
  std::int32_t num_nodes = 5;
  std::int32_t src_arr[6] = {0, 0, 0, 1, 2, 3};
  std::int32_t dst_arr[6] = {1, 2, 3, 4, 4, 4};
  double cap_arr[6] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
  std::int64_t cost_arr[6] = {1, 1, 1, 1, 1, 1};
  std::int64_t ext_arr[6] = {0, 1, 2, 3, 4, 5};
  auto g = StrictMultiDiGraph::from_arrays(num_nodes,
                                           std::span<const std::int32_t>(src_arr, 6),
                                           std::span<const std::int32_t>(dst_arr, 6),
                                           std::span<const double>(cap_arr, 6),
                                           std::span<const std::int64_t>(cost_arr, 6),
                                           std::span<const std::int64_t>(ext_arr, 6));
  FlowGraph fg(g);
  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
  auto be = make_cpu_backend(); auto algs = std::make_shared<Algorithms>(be); auto gh = algs->build_graph(g);
  ExecutionContext ctx(algs, gh);
  FlowPolicyConfig cfg;
  cfg.flow_placement = FlowPlacement::EqualBalanced;
  cfg.selection = sel;
  cfg.require_capacity = false;
  cfg.min_flow_count = 1;
  cfg.max_flow_count = 1;
  cfg.shortest_path = true;
  FlowPolicy policy(ctx, cfg);
  auto res = policy.place_demand(fg, 0, 4, 0, 100.0);
  // Should place sum of shortest-tier capacities = 15.0, and stop
  EXPECT_NEAR(res.first, 15.0, 1e-9);
  EXPECT_NEAR(res.second, 85.0, 1e-9);
  // Equal-balanced split across the three equal-cost paths
  expect_edge_flows_by_uv(fg, {{0,1,5.0}, {1,4,5.0}, {0,2,5.0}, {2,4,5.0}, {0,3,5.0}, {3,4,5.0}});
}

TEST(FlowPolicyCore, EqualBalanced_ShortestPath_DownstreamSplit_BalancesBottleneck) {
  // Graph: 0->1 (cap 10) -> {2,3} (cap 5 each) -> 4, all costs=1
  // Expect placement limited to 10, split 5/5 on downstream equal-cost arcs.
  std::int32_t num_nodes = 5;
  std::int32_t src_arr[5] = {0, 1, 1, 2, 3};
  std::int32_t dst_arr[5] = {1, 2, 3, 4, 4};
  double cap_arr[5] = {10.0, 5.0, 5.0, 10.0, 10.0};
  std::int64_t cost_arr[5] = {1, 1, 1, 1, 1};
  std::int64_t ext_arr[5] = {0, 1, 2, 3, 4};
  auto g = StrictMultiDiGraph::from_arrays(num_nodes,
                                           std::span<const std::int32_t>(src_arr, 5),
                                           std::span<const std::int32_t>(dst_arr, 5),
                                           std::span<const double>(cap_arr, 5),
                                           std::span<const std::int64_t>(cost_arr, 5),
                                           std::span<const std::int64_t>(ext_arr, 5));
  FlowGraph fg(g);
  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
  auto be = make_cpu_backend(); auto algs = std::make_shared<Algorithms>(be); auto gh = algs->build_graph(g);
  ExecutionContext ctx(algs, gh);
  FlowPolicyConfig cfg;
  cfg.flow_placement = FlowPlacement::EqualBalanced;
  cfg.selection = sel;
  cfg.require_capacity = false;
  cfg.min_flow_count = 1;
  cfg.max_flow_count = 1;
  cfg.shortest_path = true;
  FlowPolicy policy(ctx, cfg);
  auto res = policy.place_demand(fg, 0, 4, 0, 100.0);
  EXPECT_NEAR(res.first, 10.0, 1e-9);
  EXPECT_NEAR(res.second, 90.0, 1e-9);
  expect_edge_flows_by_uv(fg, {{0,1,10.0}, {1,2,5.0}, {1,3,5.0}, {2,4,5.0}, {3,4,5.0}});
}

TEST(FlowPolicyCore, EqualBalanced_ShortestPath_IgnoresHigherCostTier) {
  // Graph:
  //   Shortest tier: 0->1->4, caps 10, costs=1
  //   Higher tier:   0->2->4, caps 50, costs=2
  // With shortest_path=True, placement must saturate shortest tier only (10) and stop.
  std::int32_t num_nodes = 5;
  std::int32_t src_arr[4] = {0, 1, 0, 2};
  std::int32_t dst_arr[4] = {1, 4, 2, 4};
  double cap_arr[4] = {10.0, 10.0, 50.0, 50.0};
  std::int64_t cost_arr[4] = {1, 1, 2, 2};
  std::int64_t ext_arr[4] = {0, 1, 2, 3};
  auto g = StrictMultiDiGraph::from_arrays(num_nodes,
                                           std::span<const std::int32_t>(src_arr, 4),
                                           std::span<const std::int32_t>(dst_arr, 4),
                                           std::span<const double>(cap_arr, 4),
                                           std::span<const std::int64_t>(cost_arr, 4),
                                           std::span<const std::int64_t>(ext_arr, 4));
  FlowGraph fg(g);
  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
  auto be = make_cpu_backend(); auto algs = std::make_shared<Algorithms>(be); auto gh = algs->build_graph(g);
  ExecutionContext ctx(algs, gh);
  FlowPolicyConfig cfg;
  cfg.flow_placement = FlowPlacement::EqualBalanced;
  cfg.selection = sel;
  cfg.require_capacity = false;
  cfg.min_flow_count = 1;
  cfg.max_flow_count = 1;
  cfg.shortest_path = true;
  FlowPolicy policy(ctx, cfg);
  auto res = policy.place_demand(fg, 0, 4, 0, 100.0);
  EXPECT_NEAR(res.first, 10.0, 1e-9);
  EXPECT_NEAR(res.second, 90.0, 1e-9);
  // Only shortest tier edges should carry flow
  expect_edge_flows_by_uv(fg, {{0,1,10.0}, {1,4,10.0}, {0,2,0.0}, {2,4,0.0}});
}
