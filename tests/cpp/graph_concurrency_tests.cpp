// The graph under concurrency (APPROACH 13 threading; task 1e), run in the ThreadSanitizer suite:
// independent constructions at the same time (no shared hash or counting scratch), and full
// forward and reverse searches through both queues on one published immutable graph and on its
// copy, each thread with its own options, compared with a serial reference; every array of the
// graph is unchanged after. One batch path builds its own reverse support in its workers (no
// prepopulated frozen cache), so the reverse-source copy is read by worker threads.
#include <gtest/gtest.h>

#include <cstdint>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "graph_reference_1d.hpp"
#include "netgraph/core/detail/graph_access.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/search.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

struct Input {
  std::int32_t n {0};
  std::vector<NodeId> src {}, dst {};
  std::vector<Amount> cap {};
  std::vector<Cost> cost {};
};

Input random_input(std::int32_t n, std::size_t m, Cost max_cost, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  Input in {n};
  for (std::size_t i = 0; i < m; ++i) {
    in.src.push_back(static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n)));
    in.dst.push_back(static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n)));
    in.cap.push_back(static_cast<Amount>(1 + rng() % 4));
    in.cost.push_back(static_cast<Cost>(1 + rng() % static_cast<std::uint64_t>(max_cost)));
  }
  return in;
}

// Every array of a graph, for the before/after comparison.
struct Snapshot {
  std::vector<Amount> capacity;
  std::vector<Cost> cost;
  std::vector<NodeId> src, dst, in_col;
  std::vector<std::int32_t> input_index, row_offsets, in_row_offsets;
  std::vector<EdgeId> in_adj;
  bool operator==(const Snapshot&) const = default;
};

Snapshot snapshot(const Graph& g) {
  const auto col = detail::GraphAccess::in_col_indices(g);
  return {{g.capacity_units().begin(), g.capacity_units().end()},
          {g.cost().begin(), g.cost().end()},
          {g.edge_src().begin(), g.edge_src().end()},
          {g.edge_dst().begin(), g.edge_dst().end()},
          {col.begin(), col.end()},
          {g.input_index().begin(), g.input_index().end()},
          {g.row_offsets().begin(), g.row_offsets().end()},
          {g.in_row_offsets().begin(), g.in_row_offsets().end()},
          {g.in_adj_edge_index().begin(), g.in_adj_edge_index().end()}};
}

} // namespace

TEST(GraphConcurrency, IndependentConstructionsAtTheSameTime) {
  constexpr int kThreads = 8;
  std::vector<Input> inputs;
  for (int t = 0; t < kThreads; ++t) {
    // Alternate the hash path and the fallback (K > 65,536) across threads.
    inputs.push_back(t % 2 == 0 ? random_input(2000, 20000, 31, static_cast<std::uint64_t>(t))
                                : random_input(4096, 70000, Cost{1} << 40, static_cast<std::uint64_t>(t)));
  }
  std::vector<std::vector<std::int32_t>> got(kThreads), got_rev(kThreads);
  std::vector<std::thread> threads;
  for (int t = 0; t < kThreads; ++t) {
    threads.emplace_back([&, t] {
      const auto& in = inputs[static_cast<std::size_t>(t)];
      const Graph g(in.n, in.src, in.dst, in.cap, in.cost, 1.0);
      got[static_cast<std::size_t>(t)].assign(g.input_index().begin(), g.input_index().end());
      got_rev[static_cast<std::size_t>(t)].assign(g.in_adj_edge_index().begin(), g.in_adj_edge_index().end());
    });
  }
  for (auto& th : threads) th.join();
  for (int t = 0; t < kThreads; ++t) {
    const auto& in = inputs[static_cast<std::size_t>(t)];
    const auto r = reference_1d(in.n, in.src, in.dst, in.cap, in.cost, 1.0);
    EXPECT_EQ(got[static_cast<std::size_t>(t)], r.input_index) << t;
    EXPECT_EQ(got_rev[static_cast<std::size_t>(t)], r.in_adj_edge_index) << t;
  }
}

TEST(GraphConcurrency, SharedGraphSearchesBothDirectionsBothQueues) {
  const auto in = random_input(3000, 30000, 20, 99);
  const Graph g(in.n, in.src, in.dst, in.cap, in.cost, 1.0);
  const Graph copy(g);
  const auto before = snapshot(g);
  struct Job {
    const Graph* graph;
    Direction direction;
    SpfQueue queue;
    NodeId root;
  };
  std::vector<Job> jobs;
  for (const Graph* x : {&g, &copy}) {
    for (Direction d : {Direction::Forward, Direction::Reverse}) {
      for (SpfQueue q : {SpfQueue::Auto, SpfQueue::Heap}) {
        for (NodeId r : {0, 17, 2999}) jobs.push_back({x, d, q, r});
      }
    }
  }
  auto run = [](const Job& j) {
    const std::vector<NodeId> roots {j.root};
    // A full search keeping every tight edge (multipath, every parallel member).
    const SpfOptions o(j.direction, PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity));
    return spf(*j.graph, roots, o, j.queue);
  };
  std::vector<SpfResult> serial;
  for (const auto& j : jobs) serial.push_back(run(j));
  std::vector<std::optional<SpfResult>> parallel(jobs.size());
  std::vector<std::thread> threads;
  for (std::size_t i = 0; i < jobs.size(); ++i) threads.emplace_back([&, i] { parallel[i] = run(jobs[i]); });
  for (auto& th : threads) th.join();
  for (std::size_t i = 0; i < jobs.size(); ++i) {
    ASSERT_TRUE(parallel[i].has_value());
    EXPECT_EQ(parallel[i]->costs, serial[i].costs) << i;
    EXPECT_TRUE(parallel[i]->dag == serial[i].dag) << i;
  }
  EXPECT_TRUE(snapshot(g) == before);
  EXPECT_TRUE(snapshot(copy) == before);
}

// Workers of batch_capacity without a supplied cache each build their own reverse static support
// (the CostOnly policy's DAG toward the destinations) on the shared graph.
TEST(GraphConcurrency, BatchWorkersBuildTheirOwnReverseSupport) {
  const auto in = random_input(400, 4000, 9, 5);
  const Graph g(in.n, in.src, in.dst, in.cap, in.cost, 1.0);
  const auto before = snapshot(g);
  std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs;
  for (NodeId s = 0; s < 32; ++s) pairs.push_back({{s}, {static_cast<NodeId>(399 - s)}});
  // CostOnly routing: each policy's DAG is the reverse static support toward its destination.
  const CapacityOptions opts {FlowPolicyConfig(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless,
                                               PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity))};
  const auto serial = batch_capacity(g, pairs, opts, {}, 1);
  const auto parallel = batch_capacity(g, pairs, opts, {}, 4);
  ASSERT_EQ(serial.size(), parallel.size());
  for (std::size_t i = 0; i < serial.size(); ++i) {
    EXPECT_EQ(serial[i].placed, parallel[i].placed) << i;
    EXPECT_EQ(serial[i].edge_flows, parallel[i].edge_flows) << i;
  }
  EXPECT_TRUE(snapshot(g) == before);
}
