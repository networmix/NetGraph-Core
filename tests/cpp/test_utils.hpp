#pragma once

#include <gtest/gtest.h>
#include <span>
#include <vector>
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/flow_graph.hpp"
#include "netgraph/core/max_flow.hpp"

namespace netgraph::core::test {

// Helper to create bool arrays for masks (avoids std::vector<bool> issues with .data())
inline std::unique_ptr<bool[]> make_bool_mask(std::size_t n, bool default_val = true) {
  auto mask = std::make_unique<bool[]>(n);
  std::fill_n(mask.get(), n, default_val);
  return mask;
}

// Graph builders matching Python fixtures
inline StrictMultiDiGraph make_line_graph(int n) {
  // Create a simple line graph 0->1->2->...->n-1
  if (n <= 1) {
    return StrictMultiDiGraph::from_arrays(n, {}, {}, {}, {});
  }
  std::vector<std::int32_t> src, dst;
  std::vector<double> cap;
  std::vector<std::int64_t> cost;
  for (int i = 0; i < n - 1; ++i) {
    src.push_back(i);
    dst.push_back(i + 1);
    cap.push_back(1.0);
    cost.push_back(1);
  }
  return StrictMultiDiGraph::from_arrays(n, src, dst, cap, cost);
}

inline StrictMultiDiGraph make_square_graph(int type = 1) {
  // Type 1: one shortest route and one longer alternative
  // 0->1->2 (cost 2) vs 0->3->2 (cost 4)
  if (type == 1) {
    std::int32_t src_arr[4] = {0, 1, 0, 3};
    std::int32_t dst_arr[4] = {1, 2, 3, 2};
    double cap_arr[4] = {1.0, 1.0, 2.0, 2.0};
    std::int64_t cost_arr[4] = {1, 1, 2, 2};
    return StrictMultiDiGraph::from_arrays(4,
      std::span(src_arr, 4), std::span(dst_arr, 4),
      std::span(cap_arr, 4), std::span(cost_arr, 4));
  }
  // Type 2: two equal-cost shortest routes
  std::int32_t src_arr[4] = {0, 1, 0, 3};
  std::int32_t dst_arr[4] = {1, 2, 3, 2};
  double cap_arr[4] = {1.0, 1.0, 2.0, 2.0};
  std::int64_t cost_arr[4] = {1, 1, 1, 1};
  return StrictMultiDiGraph::from_arrays(4,
    std::span(src_arr, 4), std::span(dst_arr, 4),
    std::span(cap_arr, 4), std::span(cost_arr, 4));
}

inline StrictMultiDiGraph make_grid_graph(int rows, int cols) {
  // Create a grid graph with edges going right and down
  int n = rows * cols;
  std::vector<std::int32_t> src, dst;
  std::vector<double> cap;
  std::vector<std::int64_t> cost;

  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      int node = r * cols + c;
      // Right edge
      if (c < cols - 1) {
        src.push_back(node);
        dst.push_back(node + 1);
        cap.push_back(1.0);
        cost.push_back(1);
      }
      // Down edge
      if (r < rows - 1) {
        src.push_back(node);
        dst.push_back(node + cols);
        cap.push_back(1.0);
        cost.push_back(1);
      }
    }
  }
  return StrictMultiDiGraph::from_arrays(n, src, dst, cap, cost);
}

inline StrictMultiDiGraph make_n_disjoint_paths(int n, double capacity, std::int64_t cost = 1) {
  // Create N disjoint equal-cost paths from source to sink
  // Topology: 0 → [1..n] → (n+1)
  // Each path has 2 hops with specified capacity and cost
  int num_nodes = n + 2;  // source, n intermediates, sink
  std::vector<std::int32_t> src, dst;
  std::vector<double> cap;
  std::vector<std::int64_t> cost_vec;

  for (int i = 0; i < n; ++i) {
    // First hop: 0 → (i+1)
    src.push_back(0);
    dst.push_back(i + 1);
    cap.push_back(capacity);
    cost_vec.push_back(cost);

    // Second hop: (i+1) → (n+1)
    src.push_back(i + 1);
    dst.push_back(n + 1);
    cap.push_back(capacity);
    cost_vec.push_back(cost);
  }

  return StrictMultiDiGraph::from_arrays(num_nodes, src, dst, cap, cost_vec);
}

inline StrictMultiDiGraph make_tiered_graph(double cap_tier1, double cap_tier2,
                                             std::int64_t cost_tier1 = 1,
                                             std::int64_t cost_tier2 = 2) {
  // Create topology with two different-cost paths
  // Path 1: 0→1→3 (tier1 cost and capacity)
  // Path 2: 0→2→3 (tier2 cost and capacity)
  std::int32_t src_arr[4] = {0, 1, 0, 2};
  std::int32_t dst_arr[4] = {1, 3, 2, 3};
  double cap_arr[4] = {cap_tier1, cap_tier1, cap_tier2, cap_tier2};
  std::int64_t cost_arr[4] = {cost_tier1, cost_tier1, cost_tier2, cost_tier2};

  return StrictMultiDiGraph::from_arrays(4,
    std::span(src_arr, 4), std::span(dst_arr, 4),
    std::span(cap_arr, 4), std::span(cost_arr, 4));
}

inline StrictMultiDiGraph make_shared_bottleneck_graph(double path_cap, double bottleneck_cap) {
  // Create topology with shared bottleneck
  // Two paths: 0→1→3 and 0→2→3, where edge 1→3 is shared
  // Path 1: 0→1 (path_cap), 1→3 (bottleneck_cap)
  // Path 2: 0→2 (path_cap), 2→3 (path_cap)
  // Bottleneck: 1→3 limits total flow
  std::int32_t src_arr[4] = {0, 1, 0, 2};
  std::int32_t dst_arr[4] = {1, 3, 2, 3};
  double cap_arr[4] = {path_cap, bottleneck_cap, path_cap, path_cap};
  std::int64_t cost_arr[4] = {1, 1, 1, 1};

  return StrictMultiDiGraph::from_arrays(4,
    std::span(src_arr, 4), std::span(dst_arr, 4),
    std::span(cap_arr, 4), std::span(cost_arr, 4));
}

inline StrictMultiDiGraph make_3tier_clos_graph() {
  /**
   * Build a 3-tier Clos topology matching NetGraph scenario_3.yaml structure.
   *
   * Structure (replicates scenario_3 from NetGraph integration tests):
   * - 2 Clos instances (clos1, clos2)
   * - Each Clos has:
   *   - 2 bricks (b1, b2)
   *   - Each brick has:
   *     - 4 T1 nodes (tier1/leaf layer)
   *     - 4 T2 nodes (tier2/aggregation layer)
   *     - Mesh: all T1->T2 edges within brick (100 Gbps, cost 1)
   *   - 16 spine nodes (T3)
   *   - Connections: Each brick's 4 T2 nodes connect one-to-one to 4 spines (400 Gbps, cost 1)
   * - Inter-Clos: 16 spine-to-spine links (400 Gbps, cost 1)
   *
   * Node numbering:
   * - Clos1:
   *   - b1/t1: nodes 0-3
   *   - b1/t2: nodes 4-7
   *   - b2/t1: nodes 8-11
   *   - b2/t2: nodes 12-15
   *   - spine: nodes 16-31
   * - Clos2:
   *   - b1/t1: nodes 32-35
   *   - b1/t2: nodes 36-39
   *   - b2/t1: nodes 40-43
   *   - b2/t2: nodes 44-47
   *   - spine: nodes 48-63
   *
   * Total nodes: 64
   *
   * Expected max flow from clos1 T1 nodes to clos2 T1 nodes: 3200.0 Gbps
   * (8 inter-clos spine links × 400 Gbps = 3200 Gbps bottleneck)
   */

  constexpr int num_nodes = 64;
  std::vector<std::int32_t> src, dst;
  std::vector<double> cap;
  std::vector<std::int64_t> cost;

  // Helper to add mesh connections between two groups
  auto add_mesh = [&](int group1_start, int group1_count, int group2_start, int group2_count, double capacity) {
    for (int i = 0; i < group1_count; ++i) {
      for (int j = 0; j < group2_count; ++j) {
        src.push_back(group1_start + i);
        dst.push_back(group2_start + j);
        cap.push_back(capacity);
        cost.push_back(1);
      }
    }
  };

  // Helper to add one-to-one connections
  auto add_one_to_one = [&](int group1_start, int group1_count, int group2_start, double capacity) {
    for (int i = 0; i < group1_count; ++i) {
      src.push_back(group1_start + i);
      dst.push_back(group2_start + i);
      cap.push_back(capacity);
      cost.push_back(1);
    }
  };

  // Clos1 brick1: T1 (0-3) mesh to T2 (4-7)
  add_mesh(0, 4, 4, 4, 100.0);

  // Clos1 brick2: T1 (8-11) mesh to T2 (12-15)
  add_mesh(8, 4, 12, 4, 100.0);

  // Clos1 b1/T2 (4-7) one-to-one to spine (16-19)
  add_one_to_one(4, 4, 16, 400.0);

  // Clos1 b2/T2 (12-15) one-to-one to spine (20-23)
  add_one_to_one(12, 4, 20, 400.0);

  // Inter-Clos: Clos1 spine (16-31) one-to-one to Clos2 spine (48-63)
  add_one_to_one(16, 16, 48, 400.0);

  // Clos2 spine (48-63) to b1/T2 (36-39) - only first 4 spines connect
  add_one_to_one(48, 4, 36, 400.0);

  // Clos2 spine (52-55) to b2/T2 (44-47)
  add_one_to_one(52, 4, 44, 400.0);

  // Clos2 brick1: T2 (36-39) mesh to T1 (32-35)
  add_mesh(36, 4, 32, 4, 100.0);

  // Clos2 brick2: T2 (44-47) mesh to T1 (40-43)
  add_mesh(44, 4, 40, 4, 100.0);

  return StrictMultiDiGraph::from_arrays(num_nodes, src, dst, cap, cost);
}

inline StrictMultiDiGraph make_triangle_topology() {
  /**
   * Triangle topology from NetGraph solver tests.
   * Used in test_maxflow_api.py::_triangle_network()
   *
   * Topology:
   *   A(0) -> B(1) (cap 2)
   *   B(1) -> C(2) (cap 1)
   *   A(0) -> C(2) (cap 1)
   *
   * Expected max flow A->C: 2.0 (1 direct + 1 via B)
   */
  std::int32_t src_arr[3] = {0, 1, 0};
  std::int32_t dst_arr[3] = {1, 2, 2};
  double cap_arr[3] = {2.0, 1.0, 1.0};
  std::int64_t cost_arr[3] = {1, 1, 1};

  return StrictMultiDiGraph::from_arrays(3,
    std::span(src_arr, 3), std::span(dst_arr, 3),
    std::span(cap_arr, 3), std::span(cost_arr, 3));
}

inline StrictMultiDiGraph make_simple_parallel_paths_topology() {
  /**
   * Simple network with two disjoint parallel paths.
   * From NetGraph solver tests: test_maxflow_api.py::_simple_network()
   *
   * Topology:
   *   S(0) -> A(1) (cap 1) -> T(3) (cap 1)
   *   S(0) -> B(2) (cap 1) -> T(3) (cap 1)
   *
   * Expected max flow S->T: 2.0
   */
  std::int32_t src_arr[4] = {0, 1, 0, 2};
  std::int32_t dst_arr[4] = {1, 3, 2, 3};
  double cap_arr[4] = {1.0, 1.0, 1.0, 1.0};
  std::int64_t cost_arr[4] = {1, 1, 1, 1};

  return StrictMultiDiGraph::from_arrays(4,
    std::span(src_arr, 4), std::span(dst_arr, 4),
    std::span(cap_arr, 4), std::span(cost_arr, 4));
}

// Assertion helpers
inline void expect_csr_valid(const StrictMultiDiGraph& g) {
  auto row = g.row_offsets_view();
  auto col = g.col_indices_view();
  auto aei = g.adj_edge_index_view();

  // Offsets should be monotonic
  EXPECT_EQ(row.size(), static_cast<std::size_t>(g.num_nodes() + 1));
  for (std::size_t i = 0; i < row.size() - 1; ++i) {
    EXPECT_LE(row[i], row[i + 1]) << "Row offsets not monotonic at " << i;
  }

  // Total edges match
  EXPECT_EQ(row[row.size() - 1], g.num_edges());
  EXPECT_EQ(col.size(), static_cast<std::size_t>(g.num_edges()));
  EXPECT_EQ(aei.size(), static_cast<std::size_t>(g.num_edges()));

  // All column indices and edge indices in range
  for (std::size_t i = 0; i < col.size(); ++i) {
    EXPECT_GE(col[i], 0) << "Invalid column index at " << i;
    EXPECT_LT(col[i], g.num_nodes()) << "Column index out of range at " << i;
    EXPECT_GE(aei[i], 0) << "Invalid edge index at " << i;
    EXPECT_LT(aei[i], g.num_edges()) << "Edge index out of range at " << i;
  }
}

inline void expect_pred_dag_valid(const PredDAG& dag, int num_nodes) {
  EXPECT_EQ(dag.parent_offsets.size(), static_cast<std::size_t>(num_nodes + 1));

  // Offsets monotonic
  for (std::size_t i = 0; i < dag.parent_offsets.size() - 1; ++i) {
    EXPECT_LE(dag.parent_offsets[i], dag.parent_offsets[i + 1]);
  }

  auto total = static_cast<std::size_t>(dag.parent_offsets.back());
  EXPECT_EQ(dag.parents.size(), total);
  EXPECT_EQ(dag.via_edges.size(), total);

  // Node IDs in range
  for (auto p : dag.parents) {
    EXPECT_GE(p, 0);
    EXPECT_LT(p, num_nodes);
  }
}

inline void expect_flow_conservation(const FlowGraph& fg, NodeId src, NodeId dst) {
  const auto& g = fg.graph();
  auto row = g.row_offsets_view();
  auto col = g.col_indices_view();
  auto aei = g.adj_edge_index_view();
  auto in_row = g.in_row_offsets_view();
  auto in_col = g.in_col_indices_view();
  auto in_aei = g.in_adj_edge_index_view();
  auto flows = fg.edge_flow_view();

  // For each intermediate node, inflow should equal outflow
  for (std::int32_t u = 0; u < g.num_nodes(); ++u) {
    if (u == src || u == dst) continue;

    double outflow = 0.0;
    auto s = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto e = static_cast<std::size_t>(row[static_cast<std::size_t>(u) + 1]);
    for (std::size_t j = s; j < e; ++j) {
      outflow += flows[static_cast<std::size_t>(aei[j])];
    }

    double inflow = 0.0;
    auto is = static_cast<std::size_t>(in_row[static_cast<std::size_t>(u)]);
    auto ie = static_cast<std::size_t>(in_row[static_cast<std::size_t>(u) + 1]);
    for (std::size_t j = is; j < ie; ++j) {
      inflow += flows[static_cast<std::size_t>(in_aei[j])];
    }

    EXPECT_NEAR(inflow, outflow, 1e-9) << "Flow not conserved at node " << u;
  }
}

// FlowSummary validation helpers
inline void validate_capacity_constraints(const StrictMultiDiGraph& g, const FlowSummary& summary) {
  ASSERT_EQ(summary.edge_flows.size(), static_cast<std::size_t>(g.num_edges()))
      << "Edge flow vector size mismatch";

  auto caps = g.capacity_view();
  for (std::size_t i = 0; i < summary.edge_flows.size(); ++i) {
    EXPECT_GE(summary.edge_flows[i], 0.0)
        << "Edge " << i << " has negative flow: " << summary.edge_flows[i];
    EXPECT_LE(summary.edge_flows[i], caps[i] + 1e-9)
        << "Edge " << i << " exceeds capacity: flow=" << summary.edge_flows[i]
        << " cap=" << caps[i];
  }
}

inline void validate_flow_conservation(const StrictMultiDiGraph& g, const FlowSummary& summary,
                                       NodeId src, NodeId dst) {
  auto row = g.row_offsets_view();
  auto aei = g.adj_edge_index_view();
  auto in_row = g.in_row_offsets_view();
  auto in_aei = g.in_adj_edge_index_view();

  for (std::int32_t u = 0; u < g.num_nodes(); ++u) {
    if (u == src || u == dst) continue;  // Skip source and sink

    double outflow = 0.0;
    for (std::size_t j = row[u]; j < row[u + 1]; ++j) {
      outflow += summary.edge_flows[aei[j]];
    }

    double inflow = 0.0;
    for (std::size_t j = in_row[u]; j < in_row[u + 1]; ++j) {
      inflow += summary.edge_flows[in_aei[j]];
    }

    EXPECT_NEAR(inflow, outflow, 1e-9)
        << "Flow not conserved at node " << u
        << ": inflow=" << inflow << " outflow=" << outflow;
  }
}

} // namespace netgraph::core::test
