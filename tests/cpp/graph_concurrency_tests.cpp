// The graph under concurrency, run in the ThreadSanitizer suite:
// independent constructions at the same time (no shared hash or counting scratch; distinct uids
// from the one counter), and full forward and reverse searches through both queues on one
// published immutable graph and on its copy, each thread with its own options, compared with a
// serial reference; every array of the graph is unchanged after. One batch path builds its own
// reverse support in its workers (no prepopulated frozen cache), so the reverse-source copy is read
// by worker threads. One PredecessorDAG and one Path are read by every thread. The tests' own threads
// start their work together once all of them run (a start latch), so the calls overlap as far as
// the scheduler allows; ThreadSanitizer reports a race between unordered accesses whether or not
// they happened to overlap in time.
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <latch>
#include <optional>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "graph_inputs.hpp"
#include "graph_reference.hpp"
#include "netgraph/core/detail/route_builders.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/search.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

// m random edges of 1 to 4 units (no zero capacity).
Input random_input(std::int32_t n, std::size_t m, Cost max_edge_cost, std::uint64_t seed) {
  return test::random_input("", n, m, max_edge_cost, seed, 1, 4);
}

} // namespace

TEST(GraphConcurrency, IndependentConstructionsAtTheSameTime) {
  constexpr int kThreads = 8;
  std::vector<Input> inputs;
  for (int t = 0; t < kThreads; ++t) {
    // Alternate the hash path and the fallback (K > 65,536) across threads.
    inputs.push_back(t % 2 == 0 ? random_input(2000, 20000, 31, static_cast<std::uint64_t>(t))
                                : random_input(4096, 70000, Cost{1} << 40, static_cast<std::uint64_t>(t)));
    if (t % 2 == 0) {
      ASSERT_LE(distinct(inputs.back().cost), 65536) << t;
    } else {
      ASSERT_GT(distinct(inputs.back().cost), 65536) << t;
    }
  }
  std::vector<std::vector<std::int32_t>> got(kThreads), got_rev(kThreads);
  std::vector<std::uint64_t> uids(kThreads);
  std::vector<std::thread> threads;
  std::latch start(kThreads);  // the threads start together, once all of them are running
  for (int t = 0; t < kThreads; ++t) {
    threads.emplace_back([&, t] {
      const auto& in = inputs[static_cast<std::size_t>(t)];
      start.arrive_and_wait();
      const Graph g = build(in);
      got[static_cast<std::size_t>(t)].assign(g.input_positions().begin(), g.input_positions().end());
      got_rev[static_cast<std::size_t>(t)].assign(g.incoming_edges().begin(), g.incoming_edges().end());
      uids[static_cast<std::size_t>(t)] = g.uid();
    });
  }
  for (auto& th : threads) th.join();
  // One process-wide uid counter: the graphs of different threads never share a uid, even after
  // their threads have destroyed them.
  std::sort(uids.begin(), uids.end());
  EXPECT_GT(uids.front(), 0U);
  EXPECT_EQ(std::adjacent_find(uids.begin(), uids.end()), uids.end());
  for (int t = 0; t < kThreads; ++t) {
    const auto& in = inputs[static_cast<std::size_t>(t)];
    const auto r = reference_graph(in.n, in.src, in.dst, in.cap, in.cost, 1.0);
    EXPECT_EQ(got[static_cast<std::size_t>(t)], r.input_positions) << t;
    EXPECT_EQ(got_rev[static_cast<std::size_t>(t)], r.incoming_edges) << t;
  }
}

TEST(GraphConcurrency, SharedGraphSearchesBothDirectionsBothQueues) {
  const auto in = random_input(3000, 30000, 20, 99);
  const Graph g = build(in);
  const Graph copy(g);
  const auto before = arrays_of(g);
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
    const SpfOptions o(j.direction);
    return spf(*j.graph, roots, o, j.queue);
  };
  std::vector<SpfResult> serial;
  for (const auto& j : jobs) serial.push_back(run(j));
  // Not vacuous: every search reaches beyond its root.
  for (std::size_t i = 0; i < jobs.size(); ++i) EXPECT_FALSE(serial[i].dag.via_edges().empty()) << i;
  std::vector<std::optional<SpfResult>> parallel(jobs.size());
  std::vector<std::thread> threads;
  // The searches start together, once all threads run.
  std::latch start(static_cast<std::ptrdiff_t>(jobs.size()));
  for (std::size_t i = 0; i < jobs.size(); ++i) {
    threads.emplace_back([&, i] {
      start.arrive_and_wait();
      parallel[i] = run(jobs[i]);
    });
  }
  for (auto& th : threads) th.join();
  for (std::size_t i = 0; i < jobs.size(); ++i) {
    ASSERT_TRUE(parallel[i].has_value());
    EXPECT_EQ(parallel[i]->costs, serial[i].costs) << i;
    EXPECT_TRUE(parallel[i]->dag == serial[i].dag) << i;
  }
  EXPECT_TRUE(arrays_of(g) == before);
  EXPECT_TRUE(arrays_of(copy) == before);
}

// The immutable route values shared by threads: every thread reads every member of one PredecessorDAG
// (its sparse storage and its ascending edge list too, through the internal read access) and of
// one Path (its DAG included) 20 times, and compares each read with the values read before the
// threads started (the DAG also with a copy taken then). ThreadSanitizer reports any write racing
// with these reads (a lazily cached member, for example).
TEST(GraphConcurrency, SharedRouteValuesAreReadByEveryThread) {
  const auto in = random_input(200, 2000, 9, 7);
  const Graph g = build(in);
  const std::vector<NodeId> roots {0};
  const auto search = spf(g, roots, SpfOptions(Direction::Forward));
  const PredecessorDAG& dag = search.dag;
  const PredecessorDAG dag_before = dag;
  ASSERT_FALSE(dag.via_edges().empty());
  // A walk of single-member hops from the root: the first entry of each node, followed back.
  NodeId v = static_cast<NodeId>(g.edge_destinations()[static_cast<std::size_t>(dag.via_edges().back())]);
  std::vector<std::vector<EdgeId>> hops;
  const auto offsets = dag.parent_offsets();
  while (v != 0) {
    const auto k = static_cast<std::size_t>(offsets[static_cast<std::size_t>(v)]);
    hops.insert(hops.begin(), std::vector<EdgeId> {dag.via_edges()[k]});
    v = dag.parents()[k];
  }
  ASSERT_FALSE(hops.empty());
  const Path path(g, hops);
  using Builder = detail::PredecessorDAGBuilder;
  struct Read {
    std::uint64_t graph_uid;
    std::vector<std::int32_t> offsets;
    std::vector<NodeId> parents;
    std::vector<EdgeId> via;
    std::vector<NodeId> member_nodes;
    std::vector<std::int32_t> member_offsets;
    std::vector<EdgeId> ascending;
    std::uint64_t path_graph_uid;
    std::vector<std::vector<EdgeId>> hops;
    NodeId src;
    NodeId dst;
    Cost cost;
    bool dag_equal;
    bool path_dag_equal;
    bool operator==(const Read&) const = default;
  };
  const PredecessorDAG path_dag = path.to_dag();
  auto read = [&] {
    const auto nodes = Builder::member_nodes(dag);
    const auto node_offsets = Builder::member_offsets(dag);
    const auto ascending = Builder::ascending(dag);
    return Read {dag.graph_uid(),
                 dag.parent_offsets(),
                 {dag.parents().begin(), dag.parents().end()},
                 {dag.via_edges().begin(), dag.via_edges().end()},
                 {nodes.begin(), nodes.end()},
                 {node_offsets.begin(), node_offsets.end()},
                 {ascending.begin(), ascending.end()},
                 path.graph_uid(),
                 path.hops(),
                 path.source(),
                 path.destination(),
                 path.cost(),
                 dag == dag_before,
                 path.to_dag() == path_dag};
  };
  const Read serial = read();
  EXPECT_TRUE(serial.dag_equal && serial.path_dag_equal);
  constexpr int kThreads = 8;
  constexpr int kRounds = 20;
  std::vector<int> wrong(kThreads, 0);  // per thread, the reads that differ from the serial one
  std::latch start(kThreads);
  std::vector<std::thread> threads;
  for (int t = 0; t < kThreads; ++t) {
    threads.emplace_back([&, t] {
      start.arrive_and_wait();
      for (int i = 0; i < kRounds; ++i) {
        if (!(read() == serial)) ++wrong[static_cast<std::size_t>(t)];
      }
    });
  }
  for (auto& th : threads) th.join();
  for (int t = 0; t < kThreads; ++t) EXPECT_EQ(wrong[static_cast<std::size_t>(t)], 0) << t;
}

// Workers of evaluate_pairs without a supplied cache each build their own reverse static support
// (the ShortestPaths policy's DAG toward the destinations) on the shared graph.
TEST(GraphConcurrency, BatchWorkersBuildTheirOwnReverseSupport) {
  const auto in = random_input(400, 4000, 9, 5);
  const Graph g = build(in);
  const auto before = arrays_of(g);
  std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs;
  for (NodeId s = 0; s < 32; ++s) pairs.push_back({{s}, {static_cast<NodeId>(399 - s)}});
  // ShortestPaths routing: each policy's DAG is the reverse static support toward its destination.
  // Every output is requested, so the per-edge flows and the residuals are compared too.
  const EvaluationOptions opts {shortest_paths_config(), nullptr, nullptr, OutputOptions {true, true, true}};
  const auto serial = evaluate_pairs(g, pairs, opts, {}, 1);
  const auto parallel = evaluate_pairs(g, pairs, opts, {}, 4);
  ASSERT_EQ(serial.size(), pairs.size());
  ASSERT_EQ(parallel.size(), pairs.size());
  int placing = 0;
  for (std::size_t i = 0; i < serial.size(); ++i) {
    ASSERT_TRUE(serial[i].edge_flows.has_value() && parallel[i].edge_flows.has_value()) << i;
    ASSERT_TRUE(serial[i].residual.has_value() && parallel[i].residual.has_value()) << i;
    EXPECT_EQ(serial[i].call_delivered, parallel[i].call_delivered) << i;
    EXPECT_EQ(serial[i].edge_flows, parallel[i].edge_flows) << i;
    EXPECT_EQ(serial[i].residual, parallel[i].residual) << i;
    EXPECT_EQ(serial[i].saturated_boundary_edges, parallel[i].saturated_boundary_edges) << i;
    placing += serial[i].call_delivered > 0 ? 1 : 0;
  }
  EXPECT_GT(placing, 16);  // not vacuous: most pairs place something
  EXPECT_TRUE(arrays_of(g) == before);
}
