// Allocation-failure atomicity of the graph's construction and copy, and the allocation-free
// move (APPROACH 13; task 1e). Runs in the fault-injection executable (fault_injection.hpp), whose
// replacement operator new fails the chosen allocation and counts the live allocations
// (g_live, in flow_ledger_fault_tests.cpp's replacement functions).
//
// Every allocation of every selected construction path fails once: the constructor raises
// std::bad_alloc, the caller's input arrays are unchanged, no allocation stays live (leaks are
// counted, not left to ASan, whose leak detection is off on macOS), and a clean retry gives the
// reference arrays and passes the complete layout checker. The paths: an empty graph, isolated
// nodes, the hash path from creation through growth and collisions, the K-limit transition into
// the fallback, the fallback from the start, the probe-budget fallbacks during collection and
// during a table growth, and the copy constructor.
#include <gtest/gtest.h>

#include <atomic>
#include <cstdlib>
#include <string>
#include <vector>

#include "fault_injection.hpp"
#include "graph_reference_1d.hpp"
#include "layout_checker.hpp"
#include "netgraph/core/detail/graph_access.hpp"
#include "netgraph/core/graph.hpp"



using namespace netgraph::core;
using namespace netgraph::core::test;
using detail::GraphAccess;

namespace {

struct Input {
  std::string name;
  std::int32_t n {0};
  std::vector<NodeId> src {}, dst {};
  std::vector<Amount> cap {};
  std::vector<Cost> cost {};
};

Input exact_k(std::int64_t k, std::size_t extra, std::int32_t n) {
  Input in {"K=" + std::to_string(k), n};
  for (std::int64_t c = 1; c <= k; ++c) in.cost.push_back(c * 5 + 2);
  for (std::size_t j = 0; j < extra; ++j) in.cost.push_back(in.cost[(j * 7919) % in.cost.size()]);
  for (std::size_t i = 0; i < in.cost.size(); ++i) {
    in.src.push_back(static_cast<NodeId>((i * 31) % static_cast<std::size_t>(n)));
    in.dst.push_back(static_cast<NodeId>((i * 17 + 3) % static_cast<std::size_t>(n)));
    in.cap.push_back(static_cast<Amount>(i % 9));
  }
  return in;
}

Input colliding(std::size_t count) {
  Input in {"collisions", 8};
  for (Cost c = 1; in.cost.size() < count; ++c) {
    if ((detail::cost_hash(c) & 63) == 9) in.cost.push_back(c);
  }
  for (std::size_t i = 0; i < in.cost.size(); ++i) {
    in.src.push_back(static_cast<NodeId>(i % 8));
    in.dst.push_back(static_cast<NodeId>((i + 3) % 8));
    in.cap.push_back(1);
  }
  return in;
}

// The probe-budget fallbacks (closed-form shapes of GraphLayout.ProbeBudgetFallbackInCollectionAndGrowth):
// 64 slot-zero costs plus D duplicates alternating the 63rd and 64th, optionally a 65th one.
Input budget_shape(const std::string& name, std::size_t dups, bool tail) {
  std::vector<Cost> z;
  for (Cost c = 1; z.size() < 65; ++c) {
    if ((detail::cost_hash(c) & 2047) == 0) z.push_back(c);
  }
  Input in {name, 7};
  std::vector<Cost> c(z.begin(), z.begin() + 64);
  for (std::size_t i = 0; i < dups; ++i) c.push_back(z[i % 2 == 0 ? 62 : 63]);
  if (tail) c.push_back(z[64]);
  for (std::size_t i = 0; i < c.size(); ++i) {
    in.src.push_back(static_cast<NodeId>(i % 7));
    in.dst.push_back(static_cast<NodeId>((i * 3) % 7));
    in.cap.push_back(static_cast<Amount>(i % 11));
    in.cost.push_back(c[i]);
  }
  return in;
}

std::vector<Input> paths() {
  return {Input {"empty", 0},
          Input {"isolated", 64},
          exact_k(1, 10, 4),
          exact_k(31, 200, 16),
          exact_k(500, 100, 32),  // table growths
          colliding(40),
          exact_k(65537, 50, 4096),  // the K-limit transition after a full table
          exact_k(70000, 0, 4096),
          budget_shape("budget exhausted in collection", 10000, false),
          budget_shape("budget exhausted in growth", 9380, true)};
}

void expect_built(const Graph& g, const Input& in, const std::string& where) {
  const auto r = reference_1d(in.n, in.src, in.dst, in.cap, in.cost, 1.0);
  EXPECT_EQ(std::vector<std::int32_t>(g.input_index().begin(), g.input_index().end()), r.input_index) << where;
  EXPECT_EQ(std::vector<EdgeId>(g.in_adj_edge_index().begin(), g.in_adj_edge_index().end()), r.in_adj_edge_index)
      << where;
  const auto bad = check_layout({g.num_nodes(), g.num_edges(), g.capacity_units(), g.cost(), g.edge_src(), g.edge_dst(),
                                 g.input_index(), g.row_offsets(), g.in_row_offsets(), g.in_adj_edge_index(),
                                 GraphAccess::in_col_indices(g), g.max_cost(), g.cost_gcd()},
                                {in.n, in.src, in.dst, in.cap, in.cost});
  EXPECT_TRUE(bad.empty()) << where << ": " << (bad.empty() ? "" : bad.front());
}

} // namespace

TEST(GraphFaults, EveryConstructionAllocationFailsCleanly) {
  for (const auto& in : paths()) {
    const Input saved = in;
    long k = 0;
    for (;; ++k) {
      const long live = g_live.load();
      const bool failed = fails_at(k, [&] { const Graph g(in.n, in.src, in.dst, in.cap, in.cost, 1.0); });
      EXPECT_EQ(g_live.load(), live) << in.name << " allocation " << k << ": a leak";
      EXPECT_TRUE(in.src == saved.src && in.dst == saved.dst && in.cap == saved.cap && in.cost == saved.cost)
          << in.name << " allocation " << k << ": the input changed";
      if (!failed) break;
      const Graph retry(in.n, in.src, in.dst, in.cap, in.cost, 1.0);
      expect_built(retry, in, in.name + " retry after allocation " + std::to_string(k));
    }
    // The allocation count of each path is recorded for NOTES.md (the observed list).
    RecordProperty(in.name + " allocations", static_cast<int>(k));
    // An edgeless graph allocates its two offset arrays only; with edges, the nine retained
    // arrays at least.
    if (in.src.empty()) {
      EXPECT_EQ(k, 2) << in.name;
    } else {
      EXPECT_GE(k, 9) << in.name;
    }
  }
}

TEST(GraphFaults, EveryCopyAllocationFailsCleanly) {
  for (const auto& in : paths()) {
    const Graph g(in.n, in.src, in.dst, in.cap, in.cost, 1.0);
    for (long k = 0;; ++k) {
      const long live = g_live.load();
      const bool failed = fails_at(k, [&] { const Graph c(g); });
      EXPECT_EQ(g_live.load(), live) << in.name << " allocation " << k << ": a leak";
      expect_built(g, in, in.name + " copied-from after allocation " + std::to_string(k));
      if (!failed) break;
      const Graph retry(g);
      expect_built(retry, in, in.name + " copy retry");
      EXPECT_NE(retry.uid(), g.uid());
    }
  }
}

// The declared noexcept move allocates nothing. The armed injection fails the first allocation
// anywhere in the child: an allocating move would throw inside noexcept and terminate the child,
// which the parent reports as a failure (death-test subprocess, so the suite survives).
TEST(GraphFaults, MoveIsAllocationFreeUnderAnArmedFailure) {
  GTEST_FLAG_SET(death_test_style, "threadsafe");
  const auto in = exact_k(31, 200, 16);
  EXPECT_EXIT(
      {
        Graph g(in.n, in.src, in.dst, in.cap, in.cost, 1.0);
        const auto uid = g.uid();
        const auto edges = g.num_edges();
        const long live = g_live.load();
        g_countdown = 0;  // the next allocation, in any thread, throws
        Graph m(std::move(g));
        const bool untouched = g_countdown.load() == 0;  // no allocation was attempted
        g_countdown = -1;
        const bool ok = untouched && g_live.load() == live && m.uid() == uid && m.num_edges() == edges &&
                        g.num_edges() == 0 && g.row_offsets().size() == 1 &&  // NOLINT(bugprone-use-after-move)
                        g.resolution() == m.resolution();
        std::_Exit(ok ? 0 : 3);
      },
      ::testing::ExitedWithCode(0), "");
}
