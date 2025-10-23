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
