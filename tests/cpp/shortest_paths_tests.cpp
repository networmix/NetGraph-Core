#include <gtest/gtest.h>
#include <limits>
#include <memory>
#include <random>
#include "netgraph/core/constants.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/backend.hpp"
#include "netgraph/core/algorithms.hpp"
#include "netgraph/core/options.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

TEST(ShortestPaths, SingleSourceAllNodes) {
  auto g = make_line_graph(5);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  // Distances should be 0, 1, 2, 3, 4
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_DOUBLE_EQ(dist[1], 1.0);
  EXPECT_DOUBLE_EQ(dist[2], 2.0);
  EXPECT_DOUBLE_EQ(dist[3], 3.0);
  EXPECT_DOUBLE_EQ(dist[4], 4.0);

  expect_pred_dag_valid(dag, g.num_nodes());
}

TEST(ShortestPaths, DisconnectedComponents) {
  // Create two disconnected components: 0-1 and 2-3
  std::int32_t src[2] = {0, 2};
  std::int32_t dst[2] = {1, 3};
  double cap[2] = {1.0, 1.0};
  std::int64_t cost[2] = {1, 1};

  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src, 2), std::span(dst, 2),
    std::span(cap, 2), std::span(cost, 2));

  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  // Nodes 0, 1 should be reachable
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_DOUBLE_EQ(dist[1], 1.0);

  // Nodes 2, 3 should be unreachable
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());
  EXPECT_EQ(dist[3], std::numeric_limits<Cost>::max());
}

TEST(ShortestPaths, MultipleEqualCostPaths) {
  auto g = make_square_graph(2);  // Two equal-cost paths
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  // Both paths should have cost 2
  EXPECT_DOUBLE_EQ(dist[2], 2.0);

  // Node 2 should have two parents in the DAG
  auto start = static_cast<std::size_t>(dag.parent_offsets[2]);
  auto end = static_cast<std::size_t>(dag.parent_offsets[3]);
  EXPECT_GT(end - start, 0);  // At least one predecessor
}

TEST(ShortestPaths, ResidualAwareFiltering) {
  auto g = make_line_graph(3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = true;  // Only use edges with capacity
  sel.tie_break = EdgeTieBreak::Deterministic;

  // Create residual with first edge having zero capacity
  std::vector<Cap> residual(g.num_edges(), 1.0);
  residual[0] = 0.0;  // Block first edge

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, residual, {}, {});

  // Node 1 should be unreachable
  EXPECT_EQ(dist[1], std::numeric_limits<Cost>::max());
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());
}

TEST(ShortestPaths, NodeMaskIsolation) {
  // Verify that masking out the middle node in a line graph blocks reachability
  auto g = make_line_graph(3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  // Mask out node 1 (middle node) - blocks path from 0 to 2
  auto node_mask = make_bool_mask(g.num_nodes());
  node_mask[1] = false;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, std::span<const bool>(node_mask.get(), g.num_nodes()), {});

  // Nodes 1 and 2 should be unreachable
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_EQ(dist[1], std::numeric_limits<Cost>::max());
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());
}

TEST(ShortestPaths, EdgeMaskFiltering) {
  // Verify that masking out an edge removes it from consideration
  auto g = make_line_graph(3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  // Mask out first edge (0->1) - blocks the entire path
  auto edge_mask = make_bool_mask(g.num_edges());
  edge_mask[0] = false;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, std::span<const bool>(edge_mask.get(), g.num_edges()));

  // Node 1 should be unreachable
  EXPECT_DOUBLE_EQ(dist[0], 0.0);
  EXPECT_EQ(dist[1], std::numeric_limits<Cost>::max());
}

TEST(ShortestPaths, EarlyExitOptimization) {
  auto g = make_line_graph(10);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  // With destination specified, should exit early
  auto [dist, dag] = shortest_paths(g, 0, NodeId(2), true, sel, {}, {}, {});

  EXPECT_DOUBLE_EQ(dist[2], 2.0);
  expect_pred_dag_valid(dag, g.num_nodes());
}

TEST(ShortestPaths, PredDAGIntegrity) {
  auto g = make_grid_graph(3, 3);
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  expect_pred_dag_valid(dag, g.num_nodes());

  // Source node should have no predecessors
  EXPECT_EQ(dag.parent_offsets[0], 0);
  EXPECT_EQ(dag.parent_offsets[1], 0);
}

TEST(ShortestPaths, TieBreakingDeterminism) {
  // Two parallel edges with same cost
  std::int32_t src[2] = {0, 0};
  std::int32_t dst[2] = {1, 1};
  double cap[2] = {1.0, 2.0};
  std::int64_t cost[2] = {1, 1};

  auto g = StrictMultiDiGraph::from_arrays(2,
    std::span(src, 2), std::span(dst, 2),
    std::span(cap, 2), std::span(cost, 2));

  EdgeSelection sel;
  sel.multi_edge = false;  // Single edge selection
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  auto [dist1, dag1] = shortest_paths(g, 0, std::nullopt, false, sel, {}, {}, {});
  auto [dist2, dag2] = shortest_paths(g, 0, std::nullopt, false, sel, {}, {}, {});

  // Results should be deterministic
  EXPECT_EQ(dist1[1], dist2[1]);

  // Should select the same edge both times
  EXPECT_EQ(dag1.via_edges[0], dag2.via_edges[0]);
}

TEST(ShortestPaths, PreferHigherResidualTieBreak) {
  // Same two parallel edges but enforce HigherResidual tie-break
  std::int32_t src[2] = {0, 0};
  std::int32_t dst[2] = {1, 1};
  double cap[2] = {1.0, 2.0};
  std::int64_t cost[2] = {1, 1};
  auto g = StrictMultiDiGraph::from_arrays(2,
    std::span(src, 2), std::span(dst, 2),
    std::span(cap, 2), std::span(cost, 2));

  // Residuals equal to capacity; higher residual should pick edge with cap=2
  EdgeSelection sel;
  sel.multi_edge = false;
  sel.require_capacity = true;
  sel.tie_break = EdgeTieBreak::PreferHigherResidual;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, false, sel, g.capacity_view(), {}, {});
  ASSERT_FALSE(dag.via_edges.empty());
  // Edge ids are compacted; find which corresponds to cap=2 by checking capacity
  auto chosen = dag.via_edges[0];
  auto capv = g.capacity_view();
  EXPECT_EQ(capv[chosen], 2.0);
}

TEST(ShortestPaths, ResolveToPathsSplittingParallelEdges) {
  // Test resolve_to_paths both with and without splitting parallel edges
  // Graph 0->1 has two parallel edges, 1->2 single edge
  std::int32_t src[3] = {0, 0, 1};
  std::int32_t dst[3] = {1, 1, 2};
  double cap[3] = {1.0, 1.0, 1.0};
  std::int64_t cost[3] = {1, 1, 1};
  auto g = StrictMultiDiGraph::from_arrays(3,
    std::span(src, 3), std::span(dst, 3),
    std::span(cap, 3), std::span(cost, 3));

  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, 2, true, sel, {}, {}, {});

  // Resolve without splitting (grouped parallel edges)
  auto grouped = resolve_to_paths(dag, 0, 2, /*split_parallel_edges=*/false);
  ASSERT_FALSE(grouped.empty());
  ASSERT_GE(grouped[0].size(), 2u);
  // First tuple corresponds to node 0 with two parallel edges grouped
  EXPECT_EQ(grouped[0][0].first, 0);
  EXPECT_EQ(grouped[0][0].second.size(), 2u);

  // Resolve with splitting should produce 2 concrete paths (one per parallel edge)
  auto concrete = resolve_to_paths(dag, 0, 2, /*split_parallel_edges=*/true);
  EXPECT_EQ(concrete.size(), 2u);
}

TEST(ShortestPaths, ResolveToPathsWithMaxPathsLimit) {
  // Verify that max_paths limit caps enumeration in resolve_to_paths
  // Create graph with multiple equal-cost paths via 2x2 grid
  std::int32_t src[4] = {0, 0, 1, 2};
  std::int32_t dst[4] = {1, 2, 3, 3};
  double cap[4] = {1.0, 1.0, 1.0, 1.0};
  std::int64_t cost[4] = {1, 1, 1, 1};
  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src, 4), std::span(dst, 4),
    std::span(cap, 4), std::span(cost, 4));

  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g, 0, 3, true, sel, {}, {}, {});

  // Without limit, should find both paths
  auto all_paths = resolve_to_paths(dag, 0, 3, false, std::nullopt);
  EXPECT_GE(all_paths.size(), 1u);

  // With max_paths=1, enumeration should stop at first path
  auto limited = resolve_to_paths(dag, 0, 3, false, std::int64_t(1));
  EXPECT_EQ(limited.size(), 1u);
}

TEST(ShortestPaths, LargeGraphStressTest) {
  // Create a larger graph to test performance
  auto g = make_grid_graph(20, 20);  // 400 nodes
  EdgeSelection sel;
  sel.multi_edge = true;
  sel.require_capacity = false;
  sel.tie_break = EdgeTieBreak::Deterministic;

  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {});

  // Bottom-right corner should be reachable with cost 38 (19 right + 19 down)
  EXPECT_DOUBLE_EQ(dist[399], 38.0);
  expect_pred_dag_valid(dag, g.num_nodes());
}

TEST(ShortestPaths, RejectsMaskLengthMismatch) {
  // Simple line graph 0->1->2
  auto g = make_line_graph(3);
  auto be = make_cpu_backend();
  Algorithms algs(be);
  auto gh = algs.build_graph(g);
  SpfOptions opts;
  opts.multipath = true;
  opts.dst = 2;
  // Create node_mask with wrong length (N-1)
  auto node_mask = make_bool_mask(static_cast<std::size_t>(g.num_nodes() - 1), true);
  opts.node_mask = std::span<const bool>(node_mask.get(), static_cast<std::size_t>(g.num_nodes() - 1));
  EXPECT_THROW({ (void)algs.spf(gh, 0, opts); }, std::invalid_argument);

  // Edge mask wrong length (M-1)
  SpfOptions opts2;
  opts2.multipath = true;
  opts2.dst = 2;
  auto edge_mask = make_bool_mask(static_cast<std::size_t>(g.num_edges() - 1), true);
  opts2.edge_mask = std::span<const bool>(edge_mask.get(), static_cast<std::size_t>(g.num_edges() - 1));
  EXPECT_THROW({ (void)algs.spf(gh, 0, opts2); }, std::invalid_argument);
}

// ============================================================================
// Regression: zero-cost edges must not create cycles in the PredDAG
// ============================================================================

// A zero-cost pair 1<->2 previously recorded each node as the other's parent
// (both relaxations see an "equal cost" path), producing a cyclic PredDAG.
// EqualBalanced placement then stalled (returned 0 flow) and path enumeration
// walked the cycle forever. Equal-cost predecessors are now only accepted while
// the child is unsettled, which keeps the DAG acyclic by settle order.
TEST(ShortestPaths, ZeroCostEdges_PredDAGIsAcyclic) {
  std::int32_t src_arr[4] = {0, 1, 2, 2};
  std::int32_t dst_arr[4] = {1, 2, 1, 3};
  double cap_arr[4] = {10.0, 10.0, 10.0, 10.0};
  std::int64_t cost_arr[4] = {1, 0, 0, 1};
  auto g = StrictMultiDiGraph::from_arrays(4,
    std::span(src_arr, 4), std::span(dst_arr, 4),
    std::span(cap_arr, 4), std::span(cost_arr, 4));

  EdgeSelection sel;
  sel.multi_edge = true;
  auto [dist, dag] = shortest_paths(g, 0, std::nullopt, /*multipath=*/true, sel);

  EXPECT_EQ(dist[3], 2);

  // No 2-cycles: v's parent u must not also have v as a parent.
  for (std::int32_t v = 0; v < 4; ++v) {
    for (auto i = dag.parent_offsets[static_cast<std::size_t>(v)];
         i < dag.parent_offsets[static_cast<std::size_t>(v) + 1]; ++i) {
      auto u = dag.parents[static_cast<std::size_t>(i)];
      for (auto j = dag.parent_offsets[static_cast<std::size_t>(u)];
           j < dag.parent_offsets[static_cast<std::size_t>(u) + 1]; ++j) {
        EXPECT_NE(dag.parents[static_cast<std::size_t>(j)], v)
            << "cycle " << u << "<->" << v << " in PredDAG";
      }
    }
  }

  // Path enumeration must terminate and yield the single simple path.
  auto paths = resolve_to_paths(dag, 0, 3);
  ASSERT_EQ(paths.size(), 1u);
  EXPECT_EQ(paths[0].size(), 4u);  // 0 -> 1 -> 2 -> 3
}

// ---------------------------------------------------------------------------
// shortest_paths_to: reverse SPF toward one destination, with forced fan-out.
// ---------------------------------------------------------------------------
#include "netgraph/core/flow_state.hpp"

namespace {
// Nodes: 0=P (pseudo source), 1=S1, 2=S2, 3=M, 4=T.
// Edges: 0: P->S1 cost 0 cap 1e6; 1: P->S2 cost 0 cap 1e6;
//        2: S1->T cost 1 cap 100;  3: S2->M cost 1 cap 20; 4: M->T cost 1 cap 60.
// S1 is one hop from T, S2 two hops, so an SPF from P keeps only S1.
StrictMultiDiGraph make_fanout_graph() {
  std::int32_t src[5]  = {0, 0, 1, 2, 3};
  std::int32_t dst[5]  = {1, 2, 4, 3, 4};
  double       cap[5]  = {1e6, 1e6, 100.0, 20.0, 60.0};
  std::int64_t cost[5] = {0, 0, 1, 1, 1};
  return StrictMultiDiGraph::from_arrays(5,
    std::span(src, 5), std::span(dst, 5), std::span(cap, 5), std::span(cost, 5));
}
EdgeSelection cost_only_sel() {
  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = false; sel.tie_break = EdgeTieBreak::Deterministic;
  return sel;
}
} // namespace

TEST(ShortestPathsTo, DistancesMatchForwardDistancesToDst) {
  auto g = make_grid_graph(3, 4);
  const NodeId t = g.num_nodes() - 1;
  auto [dist_to, dag] = shortest_paths_to(g, t, /*multipath=*/true, cost_only_sel());
  expect_pred_dag_valid(dag, g.num_nodes());
  for (NodeId s = 0; s < g.num_nodes(); ++s) {
    auto [dist_from_s, fwd] = shortest_paths(g, s, t, /*multipath=*/true, cost_only_sel());
    EXPECT_EQ(dist_to[static_cast<std::size_t>(s)], dist_from_s[static_cast<std::size_t>(t)]) << "node " << s;
  }
  // Every DAG entry u -> v via e lies on a shortest u -> t walk.
  const auto esrc = g.edge_src_view(); const auto edst = g.edge_dst_view(); const auto cost = g.cost_view();
  for (NodeId v = 0; v < g.num_nodes(); ++v) {
    for (auto i = dag.parent_offsets[static_cast<std::size_t>(v)]; i < dag.parent_offsets[static_cast<std::size_t>(v)+1]; ++i) {
      const auto e = static_cast<std::size_t>(dag.via_edges[static_cast<std::size_t>(i)]);
      const NodeId u = dag.parents[static_cast<std::size_t>(i)];
      EXPECT_EQ(esrc[e], u); EXPECT_EQ(edst[e], v);
      EXPECT_EQ(dist_to[static_cast<std::size_t>(u)], cost[e] + dist_to[static_cast<std::size_t>(v)]);
    }
  }
}

TEST(ShortestPathsTo, ForwardSpfFromPseudoSourceKeepsOnlyNearestSource) {
  auto g = make_fanout_graph();
  auto [dist, dag] = shortest_paths(g, 0, 4, /*multipath=*/true, cost_only_sel());
  // T (node 4) is reached only via S1 -> T (edge 2): S2's branch costs 2 and is
  // not on a shortest P -> T path, so a placement from P never uses S2.
  ASSERT_EQ(dag.parent_offsets[5] - dag.parent_offsets[4], 1);
  EXPECT_EQ(dag.via_edges[static_cast<std::size_t>(dag.parent_offsets[4])], 2);
}

TEST(ShortestPathsTo, FanoutEdgesForceEverySourceIntoTheDag) {
  auto g = make_fanout_graph();
  EdgeId fan[2] = {0, 1};
  auto [dist, dag] = shortest_paths_to(g, 4, /*multipath=*/true, cost_only_sel(), {}, {}, {}, std::span<const EdgeId>(fan, 2));
  expect_pred_dag_valid(dag, g.num_nodes());
  EXPECT_EQ(dist[1], 1); EXPECT_EQ(dist[2], 2); EXPECT_EQ(dist[0], 1) << "pseudo source keeps the SPF distance via S1";
  // Both S1 and S2 now have P as parent via their attachment edge.
  ASSERT_EQ(dag.parent_offsets[2] - dag.parent_offsets[1], 1); EXPECT_EQ(dag.via_edges[static_cast<std::size_t>(dag.parent_offsets[1])], 0);
  ASSERT_EQ(dag.parent_offsets[3] - dag.parent_offsets[2], 1); EXPECT_EQ(dag.via_edges[static_cast<std::size_t>(dag.parent_offsets[2])], 1);
  // Entry the SPF already recorded (P->S1) is not duplicated.
  int p_entries = 0;
  for (auto v : dag.parents) if (v == 0) ++p_entries;
  EXPECT_EQ(p_entries, 2);

  // Lossless equal-balanced admission over the fan-out: shares 50/50 of 100;
  // S2's branch admits 20 of 50 (link 3) so the whole demand scales to 0.4.
  FlowState fs(g);
  EXPECT_NEAR(fs.place_on_dag(0, 4, dag, 100.0, FlowPlacement::EqualBalanced), 40.0, 1e-9);
  EXPECT_NEAR(fs.edge_flow_view()[2], 20.0, 1e-9);
  EXPECT_NEAR(fs.edge_flow_view()[3], 20.0, 1e-9);
}

TEST(ShortestPathsTo, FanoutSkipsUnreachableAndMaskedHeads) {
  auto g = make_fanout_graph();
  EdgeId fan[2] = {0, 1};
  auto edge_mask = make_bool_mask(static_cast<std::size_t>(g.num_edges()), true);
  edge_mask[3] = false;  // S2 -> M down: S2 cannot reach T
  auto [dist, dag] = shortest_paths_to(g, 4, true, cost_only_sel(), {}, {},
                                       std::span<const bool>(edge_mask.get(), static_cast<std::size_t>(g.num_edges())),
                                       std::span<const EdgeId>(fan, 2));
  EXPECT_EQ(dist[2], std::numeric_limits<Cost>::max());
  EXPECT_EQ(dag.parent_offsets[3] - dag.parent_offsets[2], 0) << "no fan-out entry to an unreachable source";
  EXPECT_EQ(dag.parent_offsets[2] - dag.parent_offsets[1], 1);
}

TEST(ShortestPathsTo, FanoutFromInteriorNodeIsRejected) {
  auto g = make_fanout_graph();
  EdgeId fan[1] = {4};  // M -> T, but M already has the incoming entry S2 -> M
  EXPECT_THROW((void)shortest_paths_to(g, 4, true, cost_only_sel(), {}, {}, {}, std::span<const EdgeId>(fan, 1)),
               std::invalid_argument);
  EdgeId bad[1] = {99};
  EXPECT_THROW((void)shortest_paths_to(g, 4, true, cost_only_sel(), {}, {}, {}, std::span<const EdgeId>(bad, 1)),
               std::invalid_argument);
}

TEST(ShortestPathsTo, SinglePathModeKeepsOneSuccessorPerNode) {
  auto g = make_n_disjoint_paths(3, 10.0);
  const NodeId t = g.num_nodes() - 1;
  auto [dist, dag] = shortest_paths_to(g, t, /*multipath=*/false, cost_only_sel());
  expect_pred_dag_valid(dag, g.num_nodes());
  // Count DAG entries leaving node 0: exactly one successor.
  int leaving_src = 0;
  for (auto p : dag.parents) if (p == 0) ++leaving_src;
  EXPECT_EQ(leaving_src, 1);
}

// ---------------------------------------------------------------------------
// Queue differential: the bucket queue must reproduce the reference heap's
// results bit for bit (distances, parent_offsets, parents, via_edges) across
// the whole option matrix, including zero-cost edges, masks, residuals,
// single-path capacity ties, destination early exit and workspace reuse.
// ---------------------------------------------------------------------------
namespace {

using SpfResult = std::pair<std::vector<Cost>, PredDAG>;

struct DiffGraph {
  StrictMultiDiGraph g;
  std::unique_ptr<bool[]> node_mask, edge_mask;
  std::vector<Cap> residual;
  std::span<const bool> nm() const { return {node_mask.get(), static_cast<std::size_t>(g.num_nodes())}; }
  std::span<const bool> em() const { return {edge_mask.get(), static_cast<std::size_t>(g.num_edges())}; }
};

struct DiffSpec {
  unsigned seed; int n; int degree; bool masked; Cost max_cost; int zero_pct; bool int64_offset;
};

DiffGraph make_diff_graph(const DiffSpec& s) {
  std::mt19937 rng(s.seed);
  std::vector<std::int32_t> src, dst;
  std::vector<Cap> cap;
  std::vector<Cost> cost;
  auto draw_cost = [&]() -> Cost {
    Cost c = (static_cast<int>(rng() % 100) < s.zero_pct) ? 0 : 1 + static_cast<Cost>(rng() % static_cast<unsigned>(s.max_cost));
    return c + (s.int64_offset ? (Cost(1) << 34) : 0);
  };
  for (int u = 0; u < s.n; ++u) {
    src.push_back(u); dst.push_back((u + 1) % s.n); cost.push_back(draw_cost()); cap.push_back(1.0 + double(rng() % 100) / 7.0);
    for (int k = 1; k < s.degree; ++k) {
      src.push_back(u); dst.push_back(static_cast<int>(rng() % s.n));   // self-loops and parallel edges allowed
      cost.push_back(draw_cost()); cap.push_back(1.0 + double(rng() % 100) / 7.0);
    }
  }
  // A root node with out-edges only (never has an incoming DAG entry): its
  // out-edges serve as legal fanout edges for the reverse search.
  const int root = s.n;
  for (int k = 0; k < 3; ++k) { src.push_back(root); dst.push_back(static_cast<int>(rng() % s.n)); cost.push_back(draw_cost()); cap.push_back(2.0); }
  DiffGraph d{StrictMultiDiGraph::from_arrays(s.n + 1, src, dst, cap, cost), nullptr, nullptr, {}};
  const int N = d.g.num_nodes(); const int E = d.g.num_edges();
  d.node_mask = std::make_unique<bool[]>(static_cast<std::size_t>(N));
  d.edge_mask = std::make_unique<bool[]>(static_cast<std::size_t>(E));
  d.residual.assign(d.g.capacity_view().begin(), d.g.capacity_view().end());
  for (int i = 0; i < N; ++i) d.node_mask[static_cast<std::size_t>(i)] = !s.masked || rng() % 23 != 0;
  for (int i = 0; i < E; ++i) {
    d.edge_mask[static_cast<std::size_t>(i)] = !s.masked || rng() % 11 != 0;
    if (s.masked && rng() % 13 == 0) d.residual[static_cast<std::size_t>(i)] = 0;
    if (s.masked && rng() % 17 == 0) d.residual[static_cast<std::size_t>(i)] = kMinCap * static_cast<Cap>(rng() % 3);
  }
  return d;
}

bool same_result(const SpfResult& a, const SpfResult& b) {
  return a.first == b.first && a.second.parent_offsets == b.second.parent_offsets &&
         a.second.parents == b.second.parents && a.second.via_edges == b.second.via_edges;
}

const std::vector<DiffSpec> kDiffSpecs = {
  {1, 40, 5, true, 31, 0, false},   {2, 43, 6, true, 31, 0, false},   {3, 46, 6, true, 31, 25, false},
  {4, 60, 4, false, 1, 0, false},   {5, 64, 5, true, 1, 30, false},   {6, 80, 6, true, 1000, 10, false},
  {7, 50, 6, true, 31, 0, true},    {8, 120, 8, false, 31, 0, false}, {9, 90, 3, true, 5, 40, false},
  {10, 30, 12, true, 3, 20, false},
  // Ring graphs (degree 1): a frontier of one node exercises the bucket
  // queue's inline item, and zero-cost edges its promotion into the buckets.
  {11, 300, 1, false, 31, 0, false}, {12, 300, 1, true, 5, 30, false},
};

template <class Fn>
void for_each_option(const DiffGraph& d, Fn&& fn) {
  for (int mp = 0; mp < 2; ++mp) for (int me = 0; me < 2; ++me) for (int tb = 0; tb < 2; ++tb)
  for (int rc = 0; rc < 2; ++rc) for (int ur = 0; ur < 2; ++ur) for (int um = 0; um < 2; ++um) {
    EdgeSelection sel;
    sel.multi_edge = me != 0;
    sel.require_capacity = rc != 0;
    sel.tie_break = tb ? EdgeTieBreak::PreferHigherResidual : EdgeTieBreak::Deterministic;
    auto res = ur ? std::span<const Cap>(d.residual) : std::span<const Cap>{};
    auto nm = um ? d.nm() : std::span<const bool>{};
    auto em = um ? d.em() : std::span<const bool>{};
    fn(mp != 0, sel, res, nm, em);
  }
}

} // namespace

TEST(ShortestPaths, QueueDifferential_Forward) {
  std::size_t checked = 0;
  for (const auto& spec : kDiffSpecs) {
    auto d = make_diff_graph(spec);
    const int n = d.g.num_nodes();
    for (NodeId s : {NodeId(0), NodeId(n / 3), NodeId(n - 2), NodeId(n - 1)}) {
      for_each_option(d, [&](bool mp, const EdgeSelection& sel, auto res, auto nm, auto em) {
        auto full = shortest_paths(d.g, s, std::nullopt, mp, sel, res, nm, em, SpfQueue::Heap);
        NodeId far = -1; Cost fd = -1;
        for (int v = 0; v < n; ++v) if (full.first[static_cast<std::size_t>(v)] != std::numeric_limits<Cost>::max() && full.first[static_cast<std::size_t>(v)] > fd) { fd = full.first[static_cast<std::size_t>(v)]; far = v; }
        const auto row = d.g.row_offsets_view(); const auto col = d.g.col_indices_view();
        std::vector<std::optional<NodeId>> dsts = {std::nullopt, s, NodeId((s + n / 2) % n), far};
        if (row[static_cast<std::size_t>(s) + 1] > row[static_cast<std::size_t>(s)]) dsts.push_back(col[static_cast<std::size_t>(row[static_cast<std::size_t>(s)])]);
        for (auto dst : dsts) {
          auto expected = shortest_paths(d.g, s, dst, mp, sel, res, nm, em, SpfQueue::Heap);
          for (int rep = 0; rep < 2; ++rep) {
            auto actual = shortest_paths(d.g, s, dst, mp, sel, res, nm, em, SpfQueue::Bucket);
            ++checked;
            ASSERT_TRUE(same_result(expected, actual))
                << "seed=" << spec.seed << " src=" << s << " dst=" << (dst ? *dst : -1)
                << " mp=" << mp << " me=" << sel.multi_edge << " rc=" << sel.require_capacity;
          }
          // Auto must equal one of the two explicit routes, and both equal each other.
          auto automatic = shortest_paths(d.g, s, dst, mp, sel, res, nm, em);
          ASSERT_TRUE(same_result(expected, automatic));
        }
      });
    }
  }
  EXPECT_GT(checked, 20000u);
}

TEST(ShortestPaths, QueueDifferential_Reverse) {
  std::size_t checked = 0;
  for (const auto& spec : kDiffSpecs) {
    auto d = make_diff_graph(spec);
    const int n = d.g.num_nodes();
    const NodeId root = n - 1;
    const auto row = d.g.row_offsets_view(); const auto aei = d.g.adj_edge_index_view();
    std::vector<EdgeId> fanout(aei.begin() + row[static_cast<std::size_t>(root)], aei.begin() + row[static_cast<std::size_t>(root) + 1]);
    for (NodeId t : {NodeId(0), NodeId(n / 3), NodeId(n - 2)}) {
      for_each_option(d, [&](bool mp, const EdgeSelection& sel, auto res, auto nm, auto em) {
        for (int with_fanout = 0; with_fanout < 2; ++with_fanout) {
          auto fo = with_fanout ? std::span<const EdgeId>(fanout) : std::span<const EdgeId>{};
          auto expected = shortest_paths_to(d.g, t, mp, sel, res, nm, em, fo, SpfQueue::Heap);
          for (int rep = 0; rep < 2; ++rep) {
            auto actual = shortest_paths_to(d.g, t, mp, sel, res, nm, em, fo, SpfQueue::Bucket);
            ++checked;
            ASSERT_TRUE(same_result(expected, actual))
                << "seed=" << spec.seed << " dst=" << t << " mp=" << mp << " fanout=" << with_fanout;
          }
        }
      });
    }
  }
  EXPECT_GT(checked, 5000u);
}

TEST(ShortestPaths, QueueDifferential_ComparisonHasTeeth) {
  // Guard: the comparison must reject a result whose predecessor order differs.
  auto d = make_diff_graph(kDiffSpecs[7]);
  EdgeSelection sel;
  auto a = shortest_paths(d.g, 0, std::nullopt, true, sel, {}, {}, {}, SpfQueue::Heap);
  auto b = a;
  bool permuted = false;
  for (std::size_t v = 0; v + 1 < b.second.parent_offsets.size() && !permuted; ++v) {
    auto lo = static_cast<std::size_t>(b.second.parent_offsets[v]);
    auto hi = static_cast<std::size_t>(b.second.parent_offsets[v + 1]);
    if (hi - lo >= 2 && b.second.parents[lo] != b.second.parents[hi - 1]) {
      std::swap(b.second.parents[lo], b.second.parents[hi - 1]);
      std::swap(b.second.via_edges[lo], b.second.via_edges[hi - 1]);
      permuted = true;
    }
  }
  ASSERT_TRUE(permuted) << "fixture has no multi-parent node";
  EXPECT_FALSE(same_result(a, b));
  EXPECT_TRUE(same_result(a, a));
}

TEST(ShortestPaths, BucketQueueKeepsCapacityOrderedPredecessors) {
  // 0 -> {1,2,3} at cost 1 with bottlenecks 1, 3, 2; each -> 4 at cost 1.
  // Equal-cost nodes settle in descending bottleneck order, so 4's parents are
  // recorded as 2, 3, 1 (not in node-id order). Both queues must agree.
  std::vector<std::int32_t> src = {0, 0, 0, 1, 2, 3};
  std::vector<std::int32_t> dst = {1, 2, 3, 4, 4, 4};
  std::vector<Cap> cap = {1.0, 3.0, 2.0, 5.0, 5.0, 5.0};
  std::vector<Cost> cost = {1, 1, 1, 1, 1, 1};
  auto g = StrictMultiDiGraph::from_arrays(5, src, dst, cap, cost);
  EdgeSelection sel;
  for (SpfQueue q : {SpfQueue::Heap, SpfQueue::Bucket, SpfQueue::Auto}) {
    auto [dist, dag] = shortest_paths(g, 0, std::nullopt, true, sel, {}, {}, {}, q);
    std::vector<NodeId> parents(dag.parents.begin() + dag.parent_offsets[4], dag.parents.begin() + dag.parent_offsets[5]);
    EXPECT_EQ(parents, (std::vector<NodeId>{2, 3, 1}));
    EXPECT_EQ(dist[4], 2);
  }
}

TEST(ShortestPaths, WorkspaceIsResetAfterFanoutThrow) {
  auto d = make_diff_graph(kDiffSpecs[1]);
  const int n = d.g.num_nodes();
  EdgeSelection sel;
  auto before = shortest_paths_to(d.g, 0, true, sel, d.residual, d.nm(), d.em(), {}, SpfQueue::Bucket);
  auto before_fwd = shortest_paths(d.g, 3, std::nullopt, true, sel, d.residual, d.nm(), d.em(), SpfQueue::Bucket);
  // A fanout edge leaving a node that already has an incoming DAG entry throws
  // after the search has filled the workspace.
  EdgeId offending = -1;
  const auto esrc = d.g.edge_src_view();
  for (std::size_t v = 0; v + 1 < before.second.parent_offsets.size() && offending < 0; ++v) {
    if (before.second.parent_offsets[v + 1] > before.second.parent_offsets[v]) {
      const auto row = d.g.row_offsets_view(); const auto aei = d.g.adj_edge_index_view();
      if (row[v + 1] > row[v]) offending = aei[static_cast<std::size_t>(row[v])];
    }
  }
  ASSERT_GE(offending, 0);
  ASSERT_TRUE(before.second.parent_offsets[static_cast<std::size_t>(esrc[static_cast<std::size_t>(offending)]) + 1] >
              before.second.parent_offsets[static_cast<std::size_t>(esrc[static_cast<std::size_t>(offending)])]);
  std::vector<EdgeId> bad = {offending};
  for (SpfQueue q : {SpfQueue::Bucket, SpfQueue::Heap}) {
    EXPECT_THROW((void)shortest_paths_to(d.g, 0, true, sel, d.residual, d.nm(), d.em(), bad, q), std::invalid_argument);
    EXPECT_THROW((void)shortest_paths_to(d.g, 0, true, sel, d.residual, d.nm(), d.em(), std::vector<EdgeId>{n * 10}, q), std::invalid_argument);
    auto after = shortest_paths_to(d.g, 0, true, sel, d.residual, d.nm(), d.em(), {}, q);
    EXPECT_TRUE(same_result(before, after));
    auto after_fwd = shortest_paths(d.g, 3, std::nullopt, true, sel, d.residual, d.nm(), d.em(), q);
    EXPECT_TRUE(same_result(before_fwd, after_fwd));
  }
}

TEST(ShortestPaths, WorkspaceReuseAcrossGraphSizesAndEmptyGraphs) {
  auto small = make_diff_graph(kDiffSpecs[0]);
  auto large = make_diff_graph(kDiffSpecs[7]);
  auto empty = StrictMultiDiGraph::from_arrays(0, {}, {}, {}, {});
  EdgeSelection sel;
  for (SpfQueue q : {SpfQueue::Bucket, SpfQueue::Heap}) {
    auto s1 = shortest_paths(small.g, 1, std::nullopt, true, sel, small.residual, small.nm(), small.em(), q);
    auto l1 = shortest_paths(large.g, 5, NodeId(7), false, sel, {}, {}, {}, q);
    auto r1 = shortest_paths_to(large.g, 9, true, sel, large.residual, {}, {}, {}, q);
    auto e1 = shortest_paths(empty, 0, std::nullopt, true, sel, {}, {}, {}, q);
    EXPECT_TRUE(e1.first.empty());
    EXPECT_EQ(e1.second.parent_offsets.size(), 1u);
    auto e2 = shortest_paths_to(empty, 0, true, sel, {}, {}, {}, {}, q);
    EXPECT_TRUE(e2.first.empty());
    auto s2 = shortest_paths(small.g, 1, std::nullopt, true, sel, small.residual, small.nm(), small.em(), q);
    auto l2 = shortest_paths(large.g, 5, NodeId(7), false, sel, {}, {}, {}, q);
    auto r2 = shortest_paths_to(large.g, 9, true, sel, large.residual, {}, {}, {}, q);
    EXPECT_TRUE(same_result(s1, s2));
    EXPECT_TRUE(same_result(l1, l2));
    EXPECT_TRUE(same_result(r1, r2));
    // A masked-out source still returns an all-INF result without touching the workspace.
    auto masked = make_bool_mask(static_cast<std::size_t>(small.g.num_nodes()), true);
    masked[1] = false;
    auto m = shortest_paths(small.g, 1, std::nullopt, true, sel, {}, std::span<const bool>(masked.get(), static_cast<std::size_t>(small.g.num_nodes())), {}, q);
    for (auto dv : m.first) EXPECT_EQ(dv, std::numeric_limits<Cost>::max());
    auto s3 = shortest_paths(small.g, 1, std::nullopt, true, sel, small.residual, small.nm(), small.em(), q);
    EXPECT_TRUE(same_result(s1, s3));
  }
}
