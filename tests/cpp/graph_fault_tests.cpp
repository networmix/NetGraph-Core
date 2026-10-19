// Allocation-failure atomicity of the graph's construction and copy, the allocation-free move,
// and the calls that build a new Eligibility, PredecessorDAG or Path from const inputs. Runs in the
// fault-injection executable (fault_injection.hpp), whose
// replacement operator new fails the chosen allocation and counts the live allocations
// (g_live, in fault_injection.cpp's replacement functions).
//
// Every allocation of every selected construction path fails once: the constructor raises
// std::bad_alloc and constructs nothing, no allocation stays live (leaks are counted, not left to
// ASan, whose leak detection is off on macOS), and a clean retry, like the construction without a
// failure, gives the reference arrays and passes the complete layout checker. The inputs are built for these paths
// (graph_inputs.hpp): an empty graph, isolated nodes, the hash path from creation through growth
// and collisions, the K-limit transition into the fallback, the probe-budget fallbacks during
// collection and during a table growth; and the copy constructor. Each factory call (the
// Eligibility of a graph and of an intersection, PredecessorDAG::from_csr and parent_offsets, a Path and
// its PredecessorDAG) fails at every allocation without building anything or leaking, and the completing
// call equals the reference.
// No branch counter is read here: graph_layout_tests.cpp observes the probe-budget paths of the
// same shapes.
#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <optional>
#include <string>
#include <vector>

#include "fault_injection.hpp"
#include "graph_inputs.hpp"
#include "graph_reference.hpp"
#include "layout_checker.hpp"
#include "netgraph/core/detail/graph_access.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/route_structures.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

std::vector<Input> paths() {
  return {Input {"empty", 0},
          Input {"isolated", 64},
          exact_k_strided(1, 10, 4),
          exact_k_strided(31, 200, 16),
          exact_k_strided(500, 100, 32),  // table growths
          colliding(40),
          exact_k_strided(65537, 50, 4096),  // the K-limit transition, then repeats
          budget_shape("budget exhausted in collection", 10000, false),
          budget_shape("budget exhausted in growth", 9380, true)};
}

void expect_built(const Graph& g, const Input& in, const std::string& where) {
  const auto r = reference_graph(in.n, in.src, in.dst, in.cap, in.cost, 1.0);
  EXPECT_EQ(std::vector<std::int32_t>(g.input_positions().begin(), g.input_positions().end()), r.input_positions) << where;
  EXPECT_EQ(std::vector<EdgeId>(g.incoming_edges().begin(), g.incoming_edges().end()), r.incoming_edges)
      << where;
  const auto bad = check_layout(spans_of(g), view(in));
  EXPECT_TRUE(bad.empty()) << where << ": " << (bad.empty() ? "" : bad.front());
}

} // namespace

TEST(GraphFaults, EveryConstructionAllocationFailsCleanly) {
  for (const auto& in : paths()) {
    long k = 0;
    for (;; ++k) {
      std::optional<Graph> g;
      const long live = g_live.load();
      const bool failed = fails_at(k, [&] { g.emplace(in.n, in.src, in.dst, in.cap, in.cost, 1.0); });
      if (!failed) {
        expect_built(*g, in, in.name + " without a failure");
        // The graph holds its retained arrays and nothing else: nine with edges, the two offset
        // arrays without (counted by the allocator, not by the graph's own list of its arrays).
        const long held = g_live.load() - live;
        EXPECT_EQ(held, in.src.empty() ? 2 : 9) << in.name << ": allocations the graph holds";
        // With edges every retained array is allocated, so the graph's own list of its arrays
        // (detail::GraphAccess::retained, which the payload checks of graph_layout_tests.cpp
        // read) must name every allocation it holds.
        if (!in.src.empty()) {
          EXPECT_EQ(static_cast<long>(detail::GraphAccess::retained(*g).arrays), held)
              << in.name << ": an array left off the list";
        }
        break;
      }
      EXPECT_FALSE(g.has_value()) << in.name << " allocation " << k;
      EXPECT_EQ(g_live.load(), live) << in.name << " allocation " << k << ": a leak";
      const Graph retry = build(in);
      expect_built(retry, in, in.name + " retry after allocation " + std::to_string(k));
    }
    // The allocation count of each path is recorded as a test property.
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
    const Graph g = build(in);
    long k = 0;
    for (;; ++k) {
      std::optional<Graph> c;
      const long live = g_live.load();
      const bool failed = fails_at(k, [&] { c.emplace(g); });
      if (!failed) {
        expect_built(*c, in, in.name + " copy without a failure");
        EXPECT_EQ(g_live.load() - live, in.src.empty() ? 2 : 9) << in.name << ": allocations the copy holds";
        EXPECT_NE(c->uid(), g.uid());
        break;
      }
      EXPECT_FALSE(c.has_value()) << in.name << " allocation " << k;
      EXPECT_EQ(g_live.load(), live) << in.name << " allocation " << k << ": a leak";
      const Graph retry(g);
      expect_built(retry, in, in.name + " copy retry");
      EXPECT_NE(retry.uid(), g.uid());
    }
    // k allocations failed before the copy succeeded: as in the construction sweep, the two offset
    // arrays of an edgeless graph, at least the nine retained arrays with edges.
    if (in.src.empty()) {
      EXPECT_EQ(k, 2) << in.name;
    } else {
      EXPECT_GE(k, 9) << in.name;
    }
  }
}

// The declared noexcept move allocates nothing. The armed injection fails the first allocation
// anywhere in the child: an allocating move would throw inside noexcept and terminate the child,
// which the parent reports as a failure (death-test subprocess, so the suite survives).
TEST(GraphFaults, MoveIsAllocationFreeUnderAnArmedFailure) {
  GTEST_FLAG_SET(death_test_style, "threadsafe");
  const auto in = exact_k_strided(31, 200, 16);
  EXPECT_EXIT(
      {
        Graph g = build(in);
        const auto uid = g.uid();
        const auto edges = g.num_edges();
        const long live = g_live.load();
        g_countdown = 0;  // the next allocation, in any thread, throws
        Graph m(std::move(g));
        const bool untouched = g_countdown.load() == 0;  // no allocation was attempted
        g_countdown = -1;
        const bool ok = untouched && g_live.load() == live && m.uid() == uid && m.num_edges() == edges &&
                        g.num_edges() == 0 && g.outgoing_offsets().size() == 1 &&  // NOLINT(bugprone-use-after-move)
                        g.resolution() == m.resolution();
        std::_Exit(ok ? 0 : 3);
      },
      ::testing::ExitedWithCode(0), "");
}

namespace {

// Every allocation of make() fails once: nothing is built and nothing stays allocated; the
// completing call equals the reference (the first call, which also sizes any per-thread scratch).
template <class T, class Make, class Same>
void sweep_factory(const std::string& name, Make make, Same same) {
  const T want = make();
  long k = 0;
  for (;; ++k) {
    std::optional<T> got;
    const long live = g_live.load();
    if (!fails_at(k, [&] { got.emplace(make()); })) {
      ASSERT_TRUE(got.has_value()) << name;
      EXPECT_TRUE(same(*got, want)) << name;
      break;
    }
    EXPECT_FALSE(got.has_value()) << name << " allocation " << k;
    EXPECT_EQ(g_live.load(), live) << name << " allocation " << k << ": a leak";
  }
  EXPECT_GT(k, 0) << name;  // not vacuous: every factory allocates
}

} // namespace

TEST(GraphFaults, EveryFactoryAllocationFailsWithoutALeak) {
  // 0 -> 1 -> 2 -> 3 with a parallel pair 1 -> 2 and a bypass 0 -> 3.
  const std::vector<NodeId> src = {0, 1, 1, 2, 0};
  const std::vector<NodeId> dst = {1, 2, 2, 3, 3};
  const std::vector<Amount> cap = {5, 5, 3, 5, 1};
  const std::vector<Cost> cost = {1, 1, 2, 1, 9};
  const Graph g(4, src, dst, cap, cost, 1.0);
  const auto edge = [&](NodeId u, NodeId v, Cost c) {
    for (EdgeId e = 0; e < g.num_edges(); ++e) {
      const auto i = static_cast<std::size_t>(e);
      if (g.edge_sources()[i] == u && g.edge_destinations()[i] == v && g.edge_costs()[i] == c) return e;
    }
    ADD_FAILURE() << "no edge " << u << " -> " << v;
    return EdgeId {0};
  };
  const EdgeId bypass[] = {edge(0, 3, 9)};
  const EdgeId costly[] = {edge(1, 2, 2)};
  const NodeId none[] = {0};
  const auto same_elig = [](const Eligibility& a, const Eligibility& b) {
    return std::ranges::equal(a.node_allowed(), b.node_allowed()) && std::ranges::equal(a.edge_allowed(), b.edge_allowed());
  };
  // Controls: an isolated node changes only the node flags; an edge changes only edge flags.
  const auto isolated = make_graph(3, {{0, 1, 1, 1}});
  const auto all = isolated.eligibility();
  EXPECT_TRUE(same_elig(all, all));
  EXPECT_FALSE(same_elig(all, isolated.eligibility(std::vector<NodeId> {2}, {})));
  EXPECT_FALSE(same_elig(all, isolated.eligibility({}, std::vector<EdgeId> {0})));
  sweep_factory<Eligibility>("Graph::eligibility", [&] { return g.eligibility(none, bypass); }, same_elig);
  const Eligibility no_bypass = g.eligibility({}, bypass);
  const Eligibility no_costly = g.eligibility({}, costly);
  sweep_factory<Eligibility>("Eligibility::intersect", [&] { return Eligibility::intersect(no_bypass, no_costly); },
                             same_elig);
  const std::vector<std::vector<EdgeId>> hops = {{edge(0, 1, 1)}, {edge(1, 2, 1), edge(1, 2, 2)}, {edge(2, 3, 1)}};
  const auto same = [](const auto& a, const auto& b) { return a == b; };
  sweep_factory<Path>("Path", [&] { return Path(g, hops); }, same);
  const Path path(g, hops);
  sweep_factory<PredecessorDAG>("Path::to_dag", [&] { return path.to_dag(); }, same);
  const PredecessorDAG dag = path.to_dag();
  sweep_factory<std::vector<std::int32_t>>("PredecessorDAG::parent_offsets", [&] { return dag.parent_offsets(); }, same);
  const auto offsets = dag.parent_offsets();
  sweep_factory<PredecessorDAG>("PredecessorDAG::from_csr", [&] { return PredecessorDAG::from_csr(g, offsets, dag.parents(), dag.via_edges()); },
                         same);
}

TEST(GraphFaults, ConstructionComparisonRejectsDifferentInputs) {
  const Input in {"comparison", 4, {0, 1, 0, 2}, {1, 3, 2, 3}, {2, 3, 4, 5}, {2, 4, 6, 8}};
  const auto g = build(in);
  EXPECT_EQ(failures_of([&] { expect_built(g, in, "same"); }), 0);
  for (int kind = 0; kind < 3; ++kind) {
    auto bad = in;
    if (kind == 0) std::swap(bad.cost[0], bad.cost[2]);
    if (kind == 1) bad.dst[3] = 0;
    if (kind == 2) ++bad.cap.back();
    EXPECT_GT(failures_of([&] { expect_built(g, bad, "different"); }), 0) << kind;
  }
}
