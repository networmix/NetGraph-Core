#include <gtest/gtest.h>
#include <limits>
#include "netgraph/core/strict_multidigraph.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

TEST(StrictMultiDiGraph, EmptyGraph) {
  auto g = StrictMultiDiGraph::from_arrays(10, {}, {}, {}, {});
  EXPECT_EQ(g.num_nodes(), 10);
  EXPECT_EQ(g.num_edges(), 0);
  expect_csr_valid(g);

  // All row offsets should be 0
  auto row = g.row_offsets_view();
  for (auto offset : row) {
    EXPECT_EQ(offset, 0);
  }
}

TEST(StrictMultiDiGraph, SingleNode) {
  auto g = StrictMultiDiGraph::from_arrays(1, {}, {}, {}, {});
  EXPECT_EQ(g.num_nodes(), 1);
  EXPECT_EQ(g.num_edges(), 0);
  expect_csr_valid(g);
}

TEST(StrictMultiDiGraph, SingleEdge) {
  std::int32_t src[1] = {0};
  std::int32_t dst[1] = {1};
  double cap[1] = {5.0};
  std::int64_t cost[1] = {10};

  auto g = StrictMultiDiGraph::from_arrays(2,
    std::span(src, 1), std::span(dst, 1),
    std::span(cap, 1), std::span(cost, 1));

  EXPECT_EQ(g.num_nodes(), 2);
  EXPECT_EQ(g.num_edges(), 1);
  expect_csr_valid(g);

  // Check edge properties
  auto cap_view = g.capacity_view();
  auto cost_view = g.cost_view();
  auto src_view = g.edge_src_view();
  auto dst_view = g.edge_dst_view();

  EXPECT_EQ(cap_view[0], 5.0);
  EXPECT_EQ(cost_view[0], 10);
  EXPECT_EQ(src_view[0], 0);
  EXPECT_EQ(dst_view[0], 1);
}

TEST(StrictMultiDiGraph, ParallelEdges) {
  // Three parallel edges from 0 to 1
  std::int32_t src[3] = {0, 0, 0};
  std::int32_t dst[3] = {1, 1, 1};
  double cap[3] = {1.0, 2.0, 3.0};
  std::int64_t cost[3] = {1, 1, 2};

  auto g = StrictMultiDiGraph::from_arrays(2,
    std::span(src, 3), std::span(dst, 3),
    std::span(cap, 3), std::span(cost, 3));

  EXPECT_EQ(g.num_nodes(), 2);
  EXPECT_EQ(g.num_edges(), 3);
  expect_csr_valid(g);

  // All three edges should be in node 0's adjacency list
  auto row = g.row_offsets_view();
  EXPECT_EQ(row[1] - row[0], 3);
}

TEST(StrictMultiDiGraph, EdgeSortingByCost) {
  // Edges with different costs should be sorted
  std::int32_t src[3] = {0, 0, 0};
  std::int32_t dst[3] = {1, 1, 1};
  double cap[3] = {1.0, 2.0, 3.0};
  std::int64_t cost[3] = {3, 1, 2};  // Unsorted

  auto g = StrictMultiDiGraph::from_arrays(2,
    std::span(src, 3), std::span(dst, 3),
    std::span(cap, 3), std::span(cost, 3));

  // Edges should be reordered by cost
  auto cost_view = g.cost_view();
  EXPECT_EQ(cost_view[0], 1);  // Lowest cost first
  EXPECT_EQ(cost_view[1], 2);
  EXPECT_EQ(cost_view[2], 3);  // Highest cost last
}

TEST(StrictMultiDiGraph, ExternalEdgeIdsPreserved) {
  std::int32_t src[3] = {0, 0, 0};
  std::int32_t dst[3] = {1, 1, 1};
  double cap[3] = {1.0, 2.0, 3.0};
  std::int64_t cost[3] = {3, 1, 2};
  std::int64_t ext[3] = {100, 200, 300};

  auto g = StrictMultiDiGraph::from_arrays(2,
    std::span(src, 3), std::span(dst, 3),
    std::span(cap, 3), std::span(cost, 3),
    std::span(ext, 3));

  // External IDs should be reordered alongside edges
  auto ext_view = g.ext_edge_ids_view();
  EXPECT_EQ(ext_view.size(), 3);

  // Edge with cost 1 (input index 1) should have ext_id 200
  auto cost_view = g.cost_view();
  for (std::size_t i = 0; i < 3; ++i) {
    if (cost_view[i] == 1) {
      EXPECT_EQ(ext_view[i], 200);
    } else if (cost_view[i] == 2) {
      EXPECT_EQ(ext_view[i], 300);
    } else if (cost_view[i] == 3) {
      EXPECT_EQ(ext_view[i], 100);
    }
  }
}

TEST(StrictMultiDiGraph, CSROffsetsMonotonicity) {
  auto g = make_grid_graph(3, 3);
  expect_csr_valid(g);

  auto row = g.row_offsets_view();
  for (std::size_t i = 0; i < row.size() - 1; ++i) {
    EXPECT_LE(row[i], row[i + 1]) << "CSR offsets not monotonic";
  }
}

TEST(StrictMultiDiGraph, CSRReverseAdjacency) {
  // Create a simple graph: 0->1, 1->2, 0->2
  std::int32_t src[3] = {0, 1, 0};
  std::int32_t dst[3] = {1, 2, 2};
  double cap[3] = {1.0, 1.0, 1.0};
  std::int64_t cost[3] = {1, 1, 1};

  auto g = StrictMultiDiGraph::from_arrays(3,
    std::span(src, 3), std::span(dst, 3),
    std::span(cap, 3), std::span(cost, 3));

  // Check reverse adjacency (incoming edges)
  auto in_row = g.in_row_offsets_view();
  auto in_col = g.in_col_indices_view();

  // Node 0 has no incoming edges
  EXPECT_EQ(in_row[1] - in_row[0], 0);

  // Node 1 has one incoming edge (from 0)
  EXPECT_EQ(in_row[2] - in_row[1], 1);
  EXPECT_EQ(in_col[in_row[1]], 0);

  // Node 2 has two incoming edges (from 0 and 1)
  EXPECT_EQ(in_row[3] - in_row[2], 2);
}

TEST(StrictMultiDiGraph, MaxCapacityEdges) {
  std::int32_t src[1] = {0};
  std::int32_t dst[1] = {1};
  double cap[1] = {std::numeric_limits<double>::max()};
  std::int64_t cost[1] = {1};

  auto g = StrictMultiDiGraph::from_arrays(2,
    std::span(src, 1), std::span(dst, 1),
    std::span(cap, 1), std::span(cost, 1));

  auto cap_view = g.capacity_view();
  EXPECT_EQ(cap_view[0], std::numeric_limits<double>::max());
}

TEST(StrictMultiDiGraph, ZeroCostEdges) {
  std::int32_t src[2] = {0, 0};
  std::int32_t dst[2] = {1, 1};
  double cap[2] = {1.0, 2.0};
  std::int64_t cost[2] = {0, 1};

  auto g = StrictMultiDiGraph::from_arrays(2,
    std::span(src, 2), std::span(dst, 2),
    std::span(cap, 2), std::span(cost, 2));

  // Zero cost edge should be first in sorted order
  auto cost_view = g.cost_view();
  EXPECT_EQ(cost_view[0], 0);
  EXPECT_EQ(cost_view[1], 1);
}

TEST(StrictMultiDiGraph, InvalidNumNodesThrows) {
  // Negative node count
  EXPECT_THROW({
    auto g = StrictMultiDiGraph::from_arrays(-1, {}, {}, {}, {});
    (void)g;
  }, std::invalid_argument);
}

TEST(StrictMultiDiGraph, MismatchedArrayLengthsThrow) {
  std::int32_t src[1] = {0};
  std::int32_t dst[1] = {0};
  double cap[1] = {1.0};
  std::int64_t cost[2] = {1, 2};
  EXPECT_THROW({
    auto g = StrictMultiDiGraph::from_arrays(1,
      std::span(src, 1), std::span(dst, 1),
      std::span(cap, 1), std::span(cost, 2));
    (void)g;
  }, std::invalid_argument);
}

TEST(StrictMultiDiGraph, ExtEdgeIdsWrongLengthThrows) {
  std::int32_t src[1] = {0};
  std::int32_t dst[1] = {0};
  double cap[1] = {1.0};
  std::int64_t cost[1] = {1};
  std::int64_t ext[2] = {10, 20};
  EXPECT_THROW({
    auto g = StrictMultiDiGraph::from_arrays(1,
      std::span(src, 1), std::span(dst, 1),
      std::span(cap, 1), std::span(cost, 1),
      std::span(ext, 2));
    (void)g;
  }, std::invalid_argument);
}

TEST(StrictMultiDiGraph, OutOfRangeNodeThrows) {
  std::int32_t src[1] = {0};
  std::int32_t dst[1] = {2}; // out of range for num_nodes=2 (0..1)
  double cap[1] = {1.0};
  std::int64_t cost[1] = {1};
  EXPECT_THROW({
    auto g = StrictMultiDiGraph::from_arrays(2,
      std::span(src, 1), std::span(dst, 1),
      std::span(cap, 1), std::span(cost, 1));
    (void)g;
  }, std::out_of_range);
}

TEST(StrictMultiDiGraph, NegativeCapacityThrows) {
  std::int32_t src[1] = {0};
  std::int32_t dst[1] = {1};
  double cap[1] = {-1.0};
  std::int64_t cost[1] = {1};
  EXPECT_THROW({
    auto g = StrictMultiDiGraph::from_arrays(2,
      std::span(src, 1), std::span(dst, 1),
      std::span(cap, 1), std::span(cost, 1));
    (void)g;
  }, std::invalid_argument);
}

TEST(StrictMultiDiGraph, NegativeCostThrows) {
  std::int32_t src[1] = {0};
  std::int32_t dst[1] = {1};
  double cap[1] = {1.0};
  std::int64_t cost[1] = {-1};
  EXPECT_THROW({
    auto g = StrictMultiDiGraph::from_arrays(2,
      std::span(src, 1), std::span(dst, 1),
      std::span(cap, 1), std::span(cost, 1));
    (void)g;
  }, std::invalid_argument);
}

TEST(StrictMultiDiGraph, SelfLoopAllowed) {
  // Self-loops (edges from node to itself) should be valid
  std::int32_t src[1] = {0};
  std::int32_t dst[1] = {0};  // Self-loop
  double cap[1] = {5.0};
  std::int64_t cost[1] = {1};

  auto g = StrictMultiDiGraph::from_arrays(1,
    std::span(src, 1), std::span(dst, 1),
    std::span(cap, 1), std::span(cost, 1));

  EXPECT_EQ(g.num_nodes(), 1);
  EXPECT_EQ(g.num_edges(), 1);

  // Verify self-loop appears in both forward and reverse adjacency
  auto row = g.row_offsets_view();
  EXPECT_EQ(row[1] - row[0], 1);  // One outgoing edge

  auto in_row = g.in_row_offsets_view();
  EXPECT_EQ(in_row[1] - in_row[0], 1);  // One incoming edge
}

TEST(StrictMultiDiGraph, NegativeNodeIdThrows) {
  // Negative node IDs should be rejected
  std::int32_t src[1] = {-1};  // Invalid negative ID
  std::int32_t dst[1] = {0};
  double cap[1] = {1.0};
  std::int64_t cost[1] = {1};
  EXPECT_THROW({
    auto g = StrictMultiDiGraph::from_arrays(1,
      std::span(src, 1), std::span(dst, 1),
      std::span(cap, 1), std::span(cost, 1));
    (void)g;
  }, std::out_of_range);
}

TEST(StrictMultiDiGraph, MultipleCostTiers) {
  // Verify edges are sorted by cost into distinct tiers
  std::int32_t src[6] = {0, 0, 0, 0, 0, 0};
  std::int32_t dst[6] = {1, 1, 1, 1, 1, 1};
  double cap[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  std::int64_t cost[6] = {10, 5, 10, 1, 5, 1};  // Three tiers: 1, 5, 10

  auto g = StrictMultiDiGraph::from_arrays(2,
    std::span(src, 6), std::span(dst, 6),
    std::span(cap, 6), std::span(cost, 6));

  auto cost_view = g.cost_view();
  // Sorted costs should be [1, 1, 5, 5, 10, 10]
  EXPECT_EQ(cost_view[0], 1);
  EXPECT_EQ(cost_view[1], 1);
  EXPECT_EQ(cost_view[2], 5);
  EXPECT_EQ(cost_view[3], 5);
  EXPECT_EQ(cost_view[4], 10);
  EXPECT_EQ(cost_view[5], 10);
}
