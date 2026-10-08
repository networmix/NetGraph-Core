// Allocation-failure atomicity of the search (APPROACH 13): every allocation an spf call makes
// (the workspace's acquisition and growth, the tight-edge pool, both queues, the pruning and the
// extraction scratch, the emitted DAG and the result) fails once; the call raises std::bad_alloc,
// and a clean retry and a different query afterwards give their references, so the thread's
// workspace is restored on every exit (no touched node, mark, list or queue entry survives).
// Runs in the fault-injection executable (fault_injection.hpp).
#include <gtest/gtest.h>

#include <exception>
#include <functional>
#include <new>
#include <optional>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "fault_injection.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/traversal.hpp"
#include "netgraph/core/search.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

struct Query {
  std::string name;
  std::vector<NodeId> roots;
  Direction direction;
  std::optional<std::vector<NodeId>> targets;
  PathSelection selection;
  bool residual;
  SpfQueue queue;
};

Graph random_graph(int n, int m, Cost max_cost, unsigned seed) {
  std::mt19937 rng(seed);
  std::vector<InEdge> edges;
  for (int i = 0; i < m; ++i) {
    edges.push_back({static_cast<NodeId>(rng() % static_cast<unsigned>(n)), static_cast<NodeId>(rng() % static_cast<unsigned>(n)),
                     static_cast<Amount>(1 + rng() % 4), static_cast<Cost>(1 + rng() % static_cast<unsigned>(max_cost))});
  }
  return make_graph(n, edges);
}

SpfResult run(const Graph& g, const Query& q, const std::vector<Amount>& residual) {
  SpfOptions o(q.direction, q.selection);
  o.targets = q.targets;
  if (q.residual) o.residual = std::span<const Amount>(residual);
  return spf(g, q.roots, o, q.queue);
}

void expect_same(const SpfResult& a, const SpfResult& b, const std::string& where) {
  EXPECT_EQ(a.costs, b.costs) << where;
  EXPECT_EQ(a.target_cost, b.target_cost) << where;
  EXPECT_TRUE(a.dag == b.dag) << where;
}

// Runs `attempt` in a new thread, whose workspace starts empty (every growth path of the first
// search is then reached), and on a std::bad_alloc runs `after` in the same thread (the injection
// is spent by then) before the exception reaches the caller's thread.
void in_fresh_thread(const std::function<void()>& attempt, const std::function<void()>& after) {
  std::exception_ptr error;
  std::thread worker([&] {
    try {
      attempt();
    } catch (const std::bad_alloc&) {
      error = std::current_exception();
      after();
    } catch (...) {
      error = std::current_exception();
    }
  });
  worker.join();
  if (error) std::rethrow_exception(error);
}

} // namespace

TEST(SearchFaults, EveryAllocationOfSpfFailsCleanly) {
  // A graph large enough that the pool and the queues grow past their first capacities, with
  // lists of several entries (unit costs on a dense part). The tiny search before every attempt
  // runs on the same workspace and reuses it, sized for the large graph (Lease reallocates only
  // below n or above 4n + 4096 nodes, so it never shrinks here): each attempt starts from a
  // workspace another graph used last. The acquisition of the node arrays is swept below, from
  // fresh threads.
  const auto big = random_graph(400, 2400, 3, 99);
  const auto tiny = make_graph(3, {{0, 1, 1, 1}, {1, 2, 1, 1}});
  std::vector<Amount> residual(static_cast<std::size_t>(big.num_edges()));
  for (std::size_t e = 0; e < residual.size(); ++e) residual[e] = big.capacity_units()[e] - static_cast<Amount>(e % 2);
  const PathSelection det(EdgeTieBreak::Deterministic, EdgePriority::Capacity);
  const PathSelection collapse(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, true, false);
  const PathSelection single(EdgeTieBreak::PreferHigherPriority, EdgePriority::Headroom, false, true);
  const std::vector<Query> queries = {
      {"full forward buckets", {0}, Direction::Forward, std::nullopt, det, false, SpfQueue::Auto},
      {"full forward heap", {0, 7}, Direction::Forward, std::nullopt, det, false, SpfQueue::Heap},
      {"full reverse collapse", {3}, Direction::Reverse, std::nullopt, collapse, true, SpfQueue::Auto},
      {"anycast forward", {1}, Direction::Forward, std::vector<NodeId> {250, 399, 17}, det, true, SpfQueue::Auto},
      {"anycast reverse heap", {250, 399}, Direction::Reverse, std::vector<NodeId> {1, 2}, collapse, false, SpfQueue::Heap},
      {"single path forward", {5}, Direction::Forward, std::vector<NodeId> {300}, single, true, SpfQueue::Auto},
      {"single path reverse", {300}, Direction::Reverse, std::vector<NodeId> {5, 6}, single, true, SpfQueue::Heap},
  };
  const Query other {"other", {2}, Direction::Forward, std::vector<NodeId> {123}, single, true, SpfQueue::Auto};
  std::vector<SpfResult> want;
  for (const auto& q : queries) want.push_back(run(big, q, residual));
  const auto want_other = run(big, other, residual);
  const auto want_tiny = spf(tiny, std::vector<NodeId> {0}, SpfOptions(Direction::Forward, det));
  for (std::size_t i = 0; i < queries.size(); ++i) {
    const auto& q = queries[i];
    long failures = 0;
    for (long k = 0;; ++k) {
      (void)spf(tiny, std::vector<NodeId> {0}, SpfOptions(Direction::Forward, det));  // reuses the workspace
      std::optional<SpfResult> got;
      const bool hit = fails_at(k, [&] { got.emplace(run(big, q, residual)); });
      if (!hit) {
        ASSERT_TRUE(got.has_value()) << q.name;
        expect_same(*got, want[i], q.name + " without a failure");
        break;
      }
      ++failures;
      EXPECT_FALSE(got.has_value()) << q.name << " at " << k;
      // The retry, another query and a tiny search on this thread see a clean workspace.
      expect_same(run(big, q, residual), want[i], q.name + " retry after " + std::to_string(k));
      expect_same(run(big, other, residual), want_other, q.name + " other after " + std::to_string(k));
      expect_same(spf(tiny, std::vector<NodeId> {0}, SpfOptions(Direction::Forward, det)), want_tiny,
                  q.name + " tiny after " + std::to_string(k));
    }
    EXPECT_GE(failures, 8) << q.name;  // the sweep reached the search's allocations, not only the result
    // The same sweep from an empty workspace (a new thread per attempt): the first acquisition of
    // the node arrays, the bucket arrays and the queues' storage fail too.
    long fresh = 0;
    for (long k = 0;; ++k) {
      std::optional<SpfResult> got;
      const bool hit = fails_at(k, [&] {
        in_fresh_thread([&] { got.emplace(run(big, q, residual)); },
                        [&] {
                          expect_same(run(big, q, residual), want[i], q.name + " fresh retry after " + std::to_string(k));
                          expect_same(run(big, other, residual), want_other, q.name + " fresh other after " + std::to_string(k));
                        });
      });
      if (!hit) {
        ASSERT_TRUE(got.has_value()) << q.name;
        expect_same(*got, want[i], q.name + " fresh without a failure");
        break;
      }
      ++fresh;
    }
    EXPECT_GT(fresh, failures) << q.name;
  }
}

TEST(SearchFaults, InternalRouteSearchesFailCleanly) {
  // The policies' entries (route_search, support_route, static_support, static_costs) under the
  // same sweep: a failure leaves nothing behind for the next search.
  const auto g = random_graph(200, 1200, 2, 7);
  std::vector<Amount> view(static_cast<std::size_t>(g.num_edges()));
  for (std::size_t e = 0; e < view.size(); ++e) view[e] = g.capacity_units()[e];
  const PathSelection single(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, false, false);
  const std::vector<NodeId> dests = {150, 199};
  const auto want_route = detail::route_search(g, detail::Gate {nullptr, view, 1}, single, 3, dests);
  const auto support = detail::static_support(g, dests, Direction::Reverse, nullptr);
  const detail::PriorityValues values {&g, EdgePriority::Residual, view};
  const auto prio = detail::PriorityFn::of(values);
  const auto want_support = detail::support_route(g, detail::PredDAGBuilder::ascending(support), 3, dests, single, prio);
  const auto want_costs = detail::static_costs(g, nullptr, dests, std::vector<NodeId> {0, 1, 2, 3});
  std::vector<std::function<void()>> calls = {
      [&] {
        auto r = detail::route_search(g, detail::Gate {nullptr, view, 1}, single, 3, dests);
        EXPECT_EQ(r->edges, want_route->edges);
      },
      [&] {
        auto r = detail::support_route(g, detail::PredDAGBuilder::ascending(support), 3, dests, single, prio);
        EXPECT_EQ(r.edges, want_support.edges);
      },
      [&] { EXPECT_TRUE(detail::static_support(g, dests, Direction::Reverse, nullptr) == support); },
      [&] { EXPECT_EQ(detail::static_costs(g, nullptr, dests, std::vector<NodeId> {0, 1, 2, 3}), want_costs); },
  };
  for (std::size_t c = 0; c < calls.size(); ++c) {
    long failures = 0;
    for (long k = 0;; ++k) {
      const bool hit = fails_at(k, [&] {
        in_fresh_thread(calls[c], [&] {
          for (auto& call : calls) call();  // every entry clean afterwards, on the same thread
        });
      });
      if (!hit) break;
      ++failures;
    }
    EXPECT_GE(failures, 6) << "call " << c;
  }
}

TEST(SearchFaults, KspAndTheResolverFailCleanly) {
  // ksp (the first path, every spur search, the candidate set, the Paths) and resolve_to_paths
  // (the walk back, the hops, the per-node lists and heaps, the frame stack, the Paths): every
  // allocation fails once, in a new thread each time (an empty workspace), and the thread's next
  // calls give their references.
  const auto g = random_graph(60, 300, 3, 21);
  const PathSelection expand(EdgeTieBreak::PreferHigherPriority, EdgePriority::Capacity, true, true);
  const PathSelection collapse(EdgeTieBreak::Deterministic, EdgePriority::Capacity, true, false);
  auto ksp_of = [&](const PathSelection& s) {
    KspOptions o(6, s);
    o.max_path_cost = Stretch {3.0};
    return ksp(g, 0, 59, o);
  };
  const auto want_expand = ksp_of(expand);
  const auto want_collapse = ksp_of(collapse);
  ASSERT_GE(want_expand.size(), 3U);
  // Layers 0 -> {1..5} -> {6..10} -> 59 with a parallel pair: 26 equal shortest walks.
  std::vector<InEdge> layered;
  for (NodeId a = 1; a <= 5; ++a) layered.push_back({0, a, 1, 1});
  for (NodeId a = 1; a <= 5; ++a) {
    for (NodeId b = 6; b <= 10; ++b) layered.push_back({a, b, 1, 1});
  }
  for (NodeId b = 6; b <= 10; ++b) layered.push_back({b, 59, 1, 1});
  layered.push_back({10, 59, 1, 1});
  const auto unit = make_graph(60, layered);
  SpfOptions full(Direction::Forward, PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity));
  const auto dag = spf(unit, std::vector<NodeId> {0}, full).dag;
  const auto want_paths = resolve_to_paths(unit, dag, 0, 59, true, 7);
  ASSERT_GE(want_paths.size(), 2U);
  auto same = [](const std::vector<std::pair<Cost, Path>>& a, const std::vector<std::pair<Cost, Path>>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (a[i].first != b[i].first || !(a[i].second == b[i].second)) return false;
    }
    return true;
  };
  std::vector<std::function<void()>> calls = {
      [&] { EXPECT_TRUE(same(ksp_of(expand), want_expand)); },
      [&] { EXPECT_TRUE(same(ksp_of(collapse), want_collapse)); },
      [&] {
        const auto got = resolve_to_paths(unit, dag, 0, 59, true, 7);
        ASSERT_EQ(got.size(), want_paths.size());
        for (std::size_t i = 0; i < got.size(); ++i) EXPECT_TRUE(got[i] == want_paths[i]);
      },
  };
  for (std::size_t c = 0; c < calls.size(); ++c) {
    long failures = 0;
    for (long k = 0;; ++k) {
      const bool hit = fails_at(k, [&] {
        in_fresh_thread(calls[c], [&] {
          for (auto& call : calls) call();
        });
      });
      if (!hit) break;
      ++failures;
    }
    EXPECT_GE(failures, 10) << "call " << c;
  }
}
