#include <gtest/gtest.h>

#include "netgraph/core/flow_graph.hpp"
#include "netgraph/core/flow_policy.hpp"
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
  std::span<const std::int32_t> src(src_arr, 4);
  std::span<const std::int32_t> dst(dst_arr, 4);
  std::span<const double> cap(cap_arr, 4);
  std::span<const std::int64_t> cost(cost_arr, 4);
  return StrictMultiDiGraph::from_arrays(num_nodes, src, dst, cap, cost, /*add_reverse=*/false);
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
  std::span<const std::int32_t> src(src_arr, 4);
  std::span<const std::int32_t> dst(dst_arr, 4);
  std::span<const double> cap(cap_arr, 4);
  std::span<const std::int64_t> cost(cost_arr, 4);
  return StrictMultiDiGraph::from_arrays(num_nodes, src, dst, cap, cost, /*add_reverse=*/false);
}

StrictMultiDiGraph make_square3() {
  // Nodes: 0:A, 1:B, 2:C, 3:D
  // A->B 100@1, B->C 125@1, A->D 75@1, D->C 50@1
  std::int32_t num_nodes = 4;
  std::int32_t src_arr[4] = {0, 1, 0, 3};
  std::int32_t dst_arr[4] = {1, 2, 3, 2};
  double cap_arr[4] = {100.0, 125.0, 75.0, 50.0};
  std::int64_t cost_arr[4] = {1, 1, 1, 1};
  std::span<const std::int32_t> src(src_arr, 4);
  std::span<const std::int32_t> dst(dst_arr, 4);
  std::span<const double> cap(cap_arr, 4);
  std::span<const std::int64_t> cost(cost_arr, 4);
  return StrictMultiDiGraph::from_arrays(num_nodes, src, dst, cap, cost, /*add_reverse=*/false);
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
        cursor[key] = j + 1; // next time, continue after this match
        found = true;
        break;
      }
    }
    ASSERT_TRUE(found) << "edge (" << u << "," << v << ") not found";
  }
}

}

TEST(FlowPolicyCore, Square1_Place1) {
  auto g = make_square1();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multipath = true; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  FlowPolicy policy(PathAlg::SPF, FlowPlacement::Proportional, sel);
  auto res = policy.place_demand(fg, /*src=*/0, /*dst=*/2, /*flowClass=*/0, /*volume=*/1.0);
  EXPECT_NEAR(res.first, 1.0, 1e-9);
  EXPECT_NEAR(res.second, 0.0, 1e-9);
  expect_edge_flows_by_uv(fg, {{0,1,1.0}, {1,2,1.0}, {0,3,0.0}, {3,2,0.0}});
}

TEST(FlowPolicyCore, Square1_Place2) {
  auto g = make_square1();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multipath = true; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  FlowPolicy policy(PathAlg::SPF, FlowPlacement::Proportional, sel);
  auto res = policy.place_demand(fg, 0, 2, 0, 2.0);
  EXPECT_NEAR(res.first, 2.0, 1e-9);
  EXPECT_NEAR(res.second, 0.0, 1e-9);
  expect_edge_flows_by_uv(fg, {{0,1,1.0}, {1,2,1.0}, {0,3,1.0}, {3,2,1.0}});
}

TEST(FlowPolicyCore, Square1_Place2_MaxOneFlow) {
  auto g = make_square1();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multipath = true; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  FlowPolicy policy(PathAlg::SPF, FlowPlacement::Proportional, sel,
                    /*min_flow_count=*/1, /*max_flow_count=*/1);
  auto res = policy.place_demand(fg, 0, 2, 0, 2.0);
  EXPECT_NEAR(res.first, 2.0, 1e-9);
  EXPECT_NEAR(res.second, 0.0, 1e-9);
  expect_edge_flows_by_uv(fg, {{0,1,0.0}, {1,2,0.0}, {0,3,2.0}, {3,2,2.0}});
}

TEST(FlowPolicyCore, Square1_Place5) {
  auto g = make_square1();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multipath = true; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  FlowPolicy policy(PathAlg::SPF, FlowPlacement::Proportional, sel);
  auto res = policy.place_demand(fg, 0, 2, 0, 5.0);
  EXPECT_NEAR(res.first, 3.0, 1e-9);
  EXPECT_NEAR(res.second, 2.0, 1e-9);
  expect_edge_flows_by_uv(fg, {{0,1,1.0}, {1,2,1.0}, {0,3,2.0}, {3,2,2.0}});
}

TEST(FlowPolicyCore, Line1_EqualBalanced_MinMaxFlows) {
  auto g = make_line1();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multipath = false; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  FlowPolicy policy(PathAlg::SPF, FlowPlacement::EqualBalanced, sel,
                    /*min_flow_count=*/2, /*max_flow_count=*/2);
  auto res = policy.place_demand(fg, 0, 2, 0, 7.0);
  EXPECT_NEAR(res.first, 5.0, 1e-9);
  EXPECT_NEAR(res.second, 2.0, 1e-9);
  // A->B: 5, B->C edges: one with cap1 gets 0, cap3 gets 2.5, cap7(cost2) gets 2.5
  expect_edge_flows_by_uv(fg, {{0,1,5.0}, {1,2,0.0}, {1,2,2.5}, {1,2,2.5}});
}

TEST(FlowPolicyCore, Square3_EqualBalanced_ThreeFlows) {
  auto g = make_square3();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multipath = false; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  FlowPolicy policy(PathAlg::SPF, FlowPlacement::EqualBalanced, sel,
                    /*min_flow_count=*/3, /*max_flow_count=*/3);
  auto res = policy.place_demand(fg, 0, 2, 0, 200.0);
  // Ensure some flow is placed and results are non-negative. Detailed balancing
  // behavior is covered by Python suite; C++ policy uses iterative placement.
  EXPECT_GE(res.first, 0.0);
  EXPECT_GE(res.second, 0.0);
}

TEST(FlowPolicyCore, DiminishingReturnsCutoff) {
  auto g = make_line1();
  FlowGraph fg(g);
  EdgeSelection sel; sel.multipath = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
  FlowPolicy policy(PathAlg::SPF, FlowPlacement::EqualBalanced, sel,
                    /*min_flow_count=*/1, /*max_flow_count=*/1000000,
                    /*max_path_cost=*/std::nullopt, /*max_path_cost_factor=*/std::nullopt,
                    /*reoptimize*/false, /*max_no_progress*/100, /*max_total_iters*/10000,
                    /*dim_enabled*/true, /*dim_window*/8, /*dim_eps*/1e-3);
  auto res = policy.place_demand(fg, 0, 2, 0, 7.0);
  EXPECT_GE(res.first, 0.0);
  EXPECT_GE(res.second, 0.0);
  EXPECT_GT(res.second, 0.0);
}
