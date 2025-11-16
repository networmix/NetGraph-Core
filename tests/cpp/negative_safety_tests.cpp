#include <gtest/gtest.h>
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/flow_graph.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

/**
 * Documents unsafe behaviors that lack runtime guards.
 *
 * These tests are SKIPPED by default and only document behaviors that SHOULD crash
 * when runtime validation is missing. They are kept separate from normative
 * correctness tests to avoid encoding unsafe behaviors as acceptable.
 *
 * Note: Corresponding Python-level lifetime tests are in test_lifetime_safety.py.
 * These C++ tests specifically document raw pointer dangers that manifest at C++ level.
 */

TEST(UnsafeBehaviorsDeathTest, GraphMismatch_DAGFromDifferentGraphCausesCrash) {
  GTEST_SKIP() << "Unsafe behavior doc test; enable under sanitizers or UB checks only.";
#if GTEST_HAS_DEATH_TEST
  // Policy graph with two parallel edges (EdgeIds 0 and 1)
  std::int32_t src1[2] = {0, 0};
  std::int32_t dst1[2] = {1, 1};
  double cap1[2] = {1.0, 1.0};
  std::int64_t cost1[2] = {1, 1};
  auto g_policy = StrictMultiDiGraph::from_arrays(2,
    std::span(src1, 2), std::span(dst1, 2),
    std::span(cap1, 2), std::span(cost1, 2));

  // FlowGraph over a different graph with only one edge (EdgeId 0 only)
  std::int32_t src2[1] = {0};
  std::int32_t dst2[1] = {1};
  double cap2[1] = {1.0};
  std::int64_t cost2[1] = {1};
  auto g_flow = StrictMultiDiGraph::from_arrays(2,
    std::span(src2, 1), std::span(dst2, 1),
    std::span(cap2, 1), std::span(cost2, 1));

  // Build DAG from policy graph (contains via edge 1 as valid index there)
  EdgeSelection sel; sel.multi_edge = true; sel.require_capacity = true; sel.tie_break = EdgeTieBreak::Deterministic;
  auto [dist, dag] = shortest_paths(g_policy, 0, 1, /*multipath=*/true, sel, {}, {}, {});

  FlowGraph fg(g_flow);
  FlowIndex idx{0, 1, 0, 0};

  // Expect a crash due to OOB when interpreting dag.via_edges against fg's graph
  EXPECT_DEATH({
    (void)fg.place(idx, 0, 1, dag, 1.0, FlowPlacement::Proportional);
  }, "");
#endif
}

TEST(UnsafeBehaviorsDeathTest, DanglingGraphPointer_CapacityViewAfterGraphDestructionCrashes) {
  GTEST_SKIP() << "Unsafe behavior doc test; enable under sanitizers or UB checks only.";
#if GTEST_HAS_DEATH_TEST
  // Allocate graph on heap to control lifetime explicitly
  auto* gp = new StrictMultiDiGraph(make_line_graph(3));
  {
    FlowGraph fg(*gp);
    // Destroy the underlying graph; fg holds a dangling pointer
    delete gp; gp = nullptr;
    // Any access into fg that consults the graph pointer is undefined
    EXPECT_DEATH({
      auto s = fg.capacity_view();
      (void)s;
    }, "");
  }
#endif
}
