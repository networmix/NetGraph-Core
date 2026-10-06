// Allocation-failure atomicity of the deletion re-solve's paths (APPROACH 2.2, 13): the capture
// of a solve's final state (on cold and growing threads), the copy and the copy assignment of a
// state (between states of different sizes, on 1, 2 and 4 threads), a resolver's construction
// and every re-solve are swept over every allocation they make, and Global sensitivity over every
// allocation on 1, 2 and 4 workers. After each injected std::bad_alloc nothing is captured or
// constructed, an assigned state keeps its old value, the saved state is unchanged, and a clean
// retry (the same resolver included) gives the reference result. Runs in the fault-injection
// executable (fault_injection.hpp).
#include <gtest/gtest.h>

#include <cstdint>
#include <exception>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "fault_injection.hpp"
#include "netgraph/core/detail/flow_core.hpp"
#include "netgraph/core/evaluators.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using detail::DeletionCapture;
using detail::DeletionResolver;
using detail::DeletionResult;
using detail::DeletionState;

namespace {

// Two sources and two destinations; deletions with full, partial and no reroute.
std::shared_ptr<const Graph> graph() {
  return std::make_shared<const Graph>(
      make_graph(8, {{0, 2, 2, 1}, {1, 2, 1, 2}, {2, 4, 2, 1}, {2, 5, 2, 3}, {1, 5, 1, 7}, {0, 4, 1, 9},
                     {4, 5, 1, 1}, {3, 6, 4, 1}, {6, 7, 4, 1}, {5, 2, 1, 1}, {2, 3, 1, 2}, {6, 4, 2, 1}}));
}

const std::vector<NodeId> kSources {0, 1};
const std::vector<NodeId> kDestinations {4, 5, 7};

FlowCoreResult capture(const Graph& g, std::shared_ptr<const Graph> owner,
                       DeletionCapture& req) {
  req.graph = std::move(owner);
  const std::vector<Amount> view(g.capacity_units().begin(), g.capacity_units().end());
  return place_global_units(g, view, kSources, kDestinations, nullptr, std::nullopt, &req);
}

bool same_result(const DeletionResult& a, const DeletionResult& b) {
  const auto& x = a.counters;
  const auto& y = b.counters;
  return a.value == b.value && a.potentials == b.potentials && a.cut == b.cut && a.edge_flows == b.edge_flows &&
         x.dijkstra_runs == y.dijkstra_runs && x.dijkstra_pops == y.dijkstra_pops && x.settled == y.settled &&
         x.arc_scans == y.arc_scans && x.dinic_phases == y.dinic_phases && x.augmentations == y.augmentations;
}

bool same_state(const DeletionState& a, const DeletionState& b) {
  return std::vector<Amount>(a.residual().begin(), a.residual().end()) ==
             std::vector<Amount>(b.residual().begin(), b.residual().end()) &&
         std::vector<detail::I128>(a.potentials().begin(), a.potentials().end()) ==
             std::vector<detail::I128>(b.potentials().begin(), b.potentials().end()) &&
         a.value() == b.value() && a.maximum() == b.maximum();
}

// Every field of two populated states: the graph and the network by identity (copies share
// them), the rest by value.
bool identical(const DeletionState& a, const DeletionState& b) {
  const auto& x = a.counters();
  const auto& y = b.counters();
  return &a.graph() == &b.graph() && &a.network() == &b.network() && same_state(a, b) &&
         x.dijkstra_runs == y.dijkstra_runs && x.dijkstra_pops == y.dijkstra_pops && x.settled == y.settled &&
         x.arc_scans == y.arc_scans && x.dinic_phases == y.dinic_phases && x.augmentations == y.augmentations;
}

// A state that never held a solve: no fields, and no resolver accepts it.
bool empty_state(const DeletionState& s) {
  if (!s.residual().empty() || !s.potentials().empty() || s.value() != 0 || s.maximum()) return false;
  try {
    DeletionResolver rs(s);
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

// A captured state with its re-solve of every positive-flow edge.
struct Saved {
  DeletionState state;
  std::vector<std::pair<EdgeId, DeletionResult>> results;
};

Saved saved(std::shared_ptr<const Graph> g, const std::vector<NodeId>& sources,
            const std::vector<NodeId>& destinations) {
  DeletionCapture req;
  req.graph = g;
  (void)place_global_units(*g, g->capacity_units(), sources, destinations, nullptr, std::nullopt, &req);
  Saved out {*req.state, {}};
  const auto flows = out.state.edge_flows();
  for (EdgeId e = 0; e < g->num_edges(); ++e) {
    if (flows[static_cast<std::size_t>(e)] > 0) out.results.emplace_back(e, DeletionResolver(out.state).resolve(e));
  }
  return out;
}

// The path 0 -> 1 -> ... -> n - 1 (cost 1 per edge) at capacity c.
Saved chain(std::int32_t n, Amount c) {
  std::vector<InEdge> edges;
  for (NodeId v = 0; v + 1 < n; ++v) edges.push_back({v, v + 1, c, 1});
  return saved(std::make_shared<const Graph>(make_graph(n, edges)), {0}, {n - 1});
}

// The state is identical to `want`'s and serves its re-solves: a resolver on it gives every
// reference result (never run on a state that is not identical: a hybrid could crash it).
bool holds(const DeletionState& s, const Saved& want) {
  if (!identical(s, want.state)) return false;
  DeletionResolver rs(s);
  for (const auto& [e, r] : want.results) {
    if (!same_result(rs.resolve(e), r)) return false;
  }
  return !want.results.empty();
}

// Whether a copy of `from` needs more storage than `to` holds: then every assignment allocates.
bool grows(const DeletionState& to, const DeletionState& from) {
  return from.residual().size() > to.residual().size() || from.potentials().size() > to.potentials().size();
}

// States of different coordinate sizes (N + 2 vertices and slots: 4 and 6 for the 2-node path,
// 8 and 14 for the 6-node paths, 10 and 34 for the 8-node graph, 13 and 24 for the 11-node path),
// and two of equal sizes with different residuals and values.
std::vector<Saved> assignment_states() {
  std::vector<Saved> out;
  out.push_back(chain(2, 2));
  out.push_back(chain(6, 2));
  out.push_back(chain(6, 3));
  const auto g = graph();
  out.push_back(saved(g, kSources, kDestinations));
  out.push_back(chain(11, 1));
  return out;
}

} // namespace

TEST(DeletionFaults, CaptureFailsWithoutAStateOrAChangedPayload) {
  const auto g = graph();
  DeletionCapture want;
  const auto payload = capture(*g, g, want);
  ASSERT_TRUE(want.state.has_value());
  int failures = 0;
  for (long at = 0;; ++at) {
    DeletionCapture req;
    FlowCoreResult got;
    if (!fails_at(at, [&] { got = capture(*g, g, req); })) {
      EXPECT_TRUE(same_payload(got, payload));
      ASSERT_TRUE(req.state.has_value());
      EXPECT_TRUE(same_state(*req.state, *want.state));
      break;
    }
    ++failures;
    EXPECT_FALSE(req.state.has_value()) << "allocation " << at;
    DeletionCapture retry;
    EXPECT_TRUE(same_payload(capture(*g, g, retry), payload)) << "retry after allocation " << at;
  }
  EXPECT_GT(failures, 10);
}

TEST(DeletionFaults, CopyAndResolverConstructionFailAtomically) {
  const auto g = graph();
  DeletionCapture req;
  (void)capture(*g, g, req);
  const DeletionState saved = *req.state;
  int failures = 0;
  for (long at = 0;; ++at) {
    std::optional<DeletionState> copy;
    if (!fails_at(at, [&] { copy.emplace(*req.state); })) break;
    ++failures;
    EXPECT_FALSE(copy.has_value()) << "copy, allocation " << at;
    EXPECT_TRUE(same_state(*req.state, saved)) << "copy, allocation " << at;
  }
  for (long at = 0;; ++at) {
    std::unique_ptr<DeletionResolver> rs;
    if (!fails_at(at, [&] { rs = std::make_unique<DeletionResolver>(*req.state); })) break;
    ++failures;
    EXPECT_EQ(rs, nullptr) << "resolver, allocation " << at;
    EXPECT_TRUE(same_state(*req.state, saved)) << "resolver, allocation " << at;
  }
  EXPECT_GT(failures, 10);
}

// Every positive-flow edge, every allocation of a re-solve on a fresh resolver (its working
// storage grows) and on a warm one; the failed resolver serves the next call correctly.
TEST(DeletionFaults, EveryReSolveFailsWithoutPoisoningItsResolver) {
  const auto g = graph();
  DeletionCapture req;
  (void)capture(*g, g, req);
  const DeletionState saved = *req.state;
  const auto flows = saved.edge_flows();
  int failures = 0;
  int candidates = 0;
  for (EdgeId e = 0; e < g->num_edges(); ++e) {
    if (flows[static_cast<std::size_t>(e)] <= 0) continue;
    ++candidates;
    const auto want = DeletionResolver(saved).resolve(e);
    for (bool warm : {false, true}) {
      DeletionResolver kept(saved);
      if (warm) (void)kept.resolve(e);
      for (long at = 0;; ++at) {
        DeletionResolver fresh(saved);
        DeletionResolver& rs = warm ? kept : fresh;
        DeletionResult got;
        if (!fails_at(at, [&] { got = rs.resolve(e); })) {
          EXPECT_TRUE(same_result(got, want)) << "edge " << e;
          break;
        }
        ++failures;
        EXPECT_TRUE(same_state(*req.state, saved)) << "edge " << e << " allocation " << at;
        EXPECT_TRUE(same_result(rs.resolve(e), want)) << "edge " << e << " retry after allocation " << at;
      }
    }
  }
  EXPECT_GE(candidates, 5);
  EXPECT_GT(failures, 50);
}

// Review-coherence RC1: a copy assignment between populated states of different coordinate sizes
// (and onto itself) fails at every allocation without changing its target; the reviewer's
// reproducer (the 6-node path assigned to the 2-node path) is one of the pairs. A failed target is
// identical to its old value and serves its re-solves; a finished one is identical to the source.
TEST(DeletionFaults, CopyAssignmentFailsWithoutChangingItsTarget) {
  const auto states = assignment_states();
  int sweeps = 0;
  int failures = 0;
  for (std::size_t i = 0; i < states.size(); ++i) {
    for (std::size_t j = 0; j < states.size(); ++j) {
      int failed_here = 0;
      for (long at = 0;; ++at) {
        DeletionState target = states[i].state;
        const std::string where = "target " + std::to_string(i) + " source " + std::to_string(j) + " allocation " +
                                  std::to_string(at);
        const bool failed = fails_at(at, [&] {
          if (i == j) {
            const DeletionState& self = target;
            target = self;
          } else {
            target = states[j].state;
          }
        });
        if (!failed) {
          EXPECT_TRUE(holds(target, states[j])) << where;
          break;
        }
        ++failed_here;
        EXPECT_TRUE(holds(target, states[i])) << where;
        ASSERT_FALSE(::testing::Test::HasFailure()) << where;
      }
      ++sweeps;
      failures += failed_here;
      if (grows(states[i].state, states[j].state)) EXPECT_GE(failed_here, 1) << "target " << i << " source " << j;
    }
  }
  EXPECT_EQ(sweeps, 25);
  EXPECT_GE(failures, 2 * sweeps);  // the copy built aside allocates the residual and the potentials
}

// An empty state (no solve) as the target keeps nothing after a failure and accepts no resolver;
// as the source it empties the target without allocating.
TEST(DeletionFaults, CopyAssignmentWithAnEmptyStateIsAllOrNothing) {
  const auto states = assignment_states();
  for (std::size_t j = 0; j < states.size(); ++j) {
    for (long at = 0;; ++at) {
      DeletionState target;
      if (!fails_at(at, [&] { target = states[j].state; })) {
        EXPECT_TRUE(holds(target, states[j])) << "source " << j;
        break;
      }
      EXPECT_TRUE(empty_state(target)) << "source " << j << " allocation " << at;
    }
    DeletionState target = states[j].state;
    const DeletionState none;
    EXPECT_FALSE(fails_at(0, [&] { target = none; })) << "target " << j;
    EXPECT_TRUE(empty_state(target)) << "target " << j;
  }
}

// The moves of a state, of a resolver and of a capture request neither allocate nor throw (the
// exchange of the copy assignment); a moved resolver serves its state.
TEST(DeletionFaults, MovesNeverAllocate) {
  const auto states = assignment_states();
  for (std::size_t i = 0; i + 1 < states.size(); ++i) {
    DeletionState target = states[i].state;
    DeletionState source = states[i + 1].state;
    EXPECT_FALSE(fails_at(0, [&] { target = std::move(source); })) << i;
    EXPECT_TRUE(identical(target, states[i + 1].state)) << i;
    EXPECT_FALSE(fails_at(0, [&] { DeletionState moved(std::move(target)); })) << i;
    DeletionResolver first(states[i].state);
    DeletionResolver second(states[i + 1].state);
    EXPECT_FALSE(fails_at(0, [&] { first = std::move(second); })) << i;
    const auto& [e, want] = states[i + 1].results.front();
    EXPECT_TRUE(same_result(first.resolve(e), want)) << i;
    DeletionCapture req;
    req.state = states[i].state;
    DeletionCapture other;
    EXPECT_FALSE(fails_at(0, [&] { other = std::move(req); })) << i;
    ASSERT_TRUE(other.state.has_value());
    EXPECT_TRUE(identical(*other.state, states[i].state)) << i;
  }
}

// Copy assignments on 1, 2 and 4 threads at once, each thread into its own target from a source
// whose graph and network the other threads' sources share; every allocation of the run fails
// once (any thread's, or the workers' only). A target whose assignment failed or never started
// keeps its old value, every other one holds its source.
TEST(DeletionFaults, CopyAssignmentsOnWorkerThreadsAreAllOrNothing) {
  const auto states = assignment_states();
  const auto k = states.size();
  for (std::size_t workers : {1, 2, 4}) {
    for (bool workers_only : {false, true}) {
      if (workers == 1 && workers_only) continue;
      int failures = 0;
      for (long at = 0;; ++at) {
        std::vector<DeletionState> targets;
        for (std::size_t t = 0; t < workers; ++t) targets.push_back(states[t % k].state);
        std::vector<std::uint8_t> done(workers, 0);
        std::vector<std::exception_ptr> errors(workers);
        auto body = [&](std::size_t t) noexcept {
          try {
            targets[t] = states[(t + 1) % k].state;
            done[t] = 1;
          } catch (...) {
            errors[t] = std::current_exception();
          }
        };
        auto run = [&] {
          if (workers == 1) {
            body(0);
          } else {
            std::vector<std::thread> pool;
            pool.reserve(workers);
            try {
              for (std::size_t t = 0; t < workers; ++t) pool.emplace_back(body, t);
            } catch (...) {
              for (auto& th : pool) th.join();
              throw;
            }
            for (auto& th : pool) th.join();
          }
          for (const auto& e : errors) {
            if (e) std::rethrow_exception(e);
          }
        };
        const bool failed = workers_only ? fails_in_workers_at(at, run) : fails_at(at, run);
        for (std::size_t t = 0; t < workers; ++t) {
          const auto& want = done[t] != 0 ? states[(t + 1) % k] : states[t % k];
          EXPECT_TRUE(holds(targets[t], want)) << "workers " << workers << " thread " << t << " allocation " << at;
        }
        if (!failed) {
          for (std::size_t t = 0; t < workers; ++t) EXPECT_EQ(done[t], 1) << "workers " << workers << " thread " << t;
          break;
        }
        ++failures;
        ASSERT_FALSE(::testing::Test::HasFailure());
        ASSERT_LT(at, 1000);
      }
      // Thread 0 assigns the 6-node path to the 2-node path: it allocates in every run.
      EXPECT_GE(failures, 1) << "workers " << workers << " workers_only " << workers_only;
    }
  }
}

// Review-coherence's strengthened sweep, kept as regression evidence: every allocation of the
// capture on a thread whose solver is cold, and on one whose solver first served a smaller graph
// (its storage grows), leaves no state and an unchanged payload; the retry on the same thread
// gives the reference payload and state.
TEST(DeletionFaults, CaptureFailsAtEveryAllocationOnAColdOrGrowingThread) {
  const auto g = graph();
  DeletionCapture ref;
  const auto payload = capture(*g, g, ref);
  for (bool grow : {false, true}) {
    int failures = 0;
    for (long at = 0;; ++at) {
      bool failed = false;
      std::thread worker([&] {
        if (grow) (void)chain(2, 2);
        DeletionCapture req;
        FlowCoreResult got;
        failed = fails_at(at, [&] { got = capture(*g, g, req); });
        if (failed) {
          EXPECT_FALSE(req.state.has_value()) << "allocation " << at;
          got = capture(*g, g, req);
        }
        EXPECT_TRUE(same_payload(got, payload)) << "allocation " << at;
        ASSERT_TRUE(req.state.has_value());
        EXPECT_TRUE(same_state(*req.state, *ref.state)) << "allocation " << at;
      });
      worker.join();
      if (!failed) break;
      ++failures;
      ASSERT_FALSE(::testing::Test::HasFailure());
      ASSERT_LT(at, 1000);
    }
    EXPECT_GT(failures, 30) << "grow " << grow;
  }
}

// Review-coherence's strengthened sweep, kept as regression evidence: Global sensitivity at every
// allocation (any thread's, and the workers' only) on 1, 2 and 4 workers; after each failure a
// retry gives the reference baseline and every entry.
TEST(DeletionFaults, GlobalSensitivityFailsAtEveryAllocation) {
  const auto g = std::make_shared<const Graph>(
      make_graph(5, {{0, 1, 2, 1}, {1, 2, 2, 1}, {1, 3, 1, 2}, {3, 2, 1, 2}, {2, 4, 2, 1}, {0, 3, 1, 5}}));
  const FlowPolicyConfig cfg {Global {}, FlowPlacement::MaxFlow, Admission::Lossless,
                              PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity)};
  const SensitivityOptions opts {cfg, nullptr, nullptr, OutputOptions {true, true, true}};
  const std::vector<NodeId> sources {0};
  const std::vector<NodeId> destinations {4};
  const auto want = sensitivity(*g, sources, destinations, opts, 1);
  ASSERT_EQ(want.entries.size(), 3U);  // the cheapest path 0 -> 1 -> 2 -> 4 carries the maximum 2
  const auto same = [&](const SensitivityResult& got, const std::string& where) {
    expect_same(got.baseline, want.baseline, where);
    ASSERT_EQ(got.entries.size(), want.entries.size()) << where;
    for (std::size_t i = 0; i < want.entries.size(); ++i) {
      EXPECT_EQ(got.entries[i].edge, want.entries[i].edge) << where;
      EXPECT_EQ(got.entries[i].loss, want.entries[i].loss) << where;
      EXPECT_EQ(got.entries[i].status, want.entries[i].status) << where;
      EXPECT_EQ(got.entries[i].reason, want.entries[i].reason) << where;
    }
  };
  for (std::int32_t workers : {1, 2, 4}) {
    for (bool workers_only : {false, true}) {
      if (workers == 1 && workers_only) continue;
      int failures = 0;
      for (long at = 0;; ++at) {
        SensitivityResult got;
        const auto call = [&] { got = sensitivity(*g, sources, destinations, opts, workers); };
        const std::string where = "workers " + std::to_string(workers) + (workers_only ? " (workers only)" : "") +
                                  " allocation " + std::to_string(at);
        if (!(workers_only ? fails_in_workers_at(at, call) : fails_at(at, call))) {
          same(got, where);
          break;
        }
        ++failures;
        same(sensitivity(*g, sources, destinations, opts, workers), "retry after " + where);
        ASSERT_FALSE(::testing::Test::HasFailure());
        ASSERT_LT(at, 10000);
      }
      EXPECT_GT(failures, workers_only ? 10 : 100) << "workers " << workers << " workers_only " << workers_only;
    }
  }
}
