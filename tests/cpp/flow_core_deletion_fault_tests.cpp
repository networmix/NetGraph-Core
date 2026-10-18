// Allocation-failure atomicity of the deletion re-solve's paths: the capture
// of a solve's final state (on cold and growing threads), the copy and the copy assignment of a
// state (between states of different sizes, on 1, 2 and 4 threads), a resolver's construction
// and every re-solve are swept over every allocation they make, and Global sensitivity over every
// allocation on 1, 2 and 4 workers. After each injected std::bad_alloc nothing is captured or
// constructed, no allocation is left live once the attempt's objects are gone, an
// assigned state keeps its old value, the network every state and resolver shares is unchanged,
// and a clean retry (the same resolver included) gives the reference result. A state's and a
// resolver's moves allocate nothing. The state comparators of these checks have negative
// controls field by field (the last test). Runs in the fault-injection executable
// (fault_injection.hpp).
#include <gtest/gtest.h>

#include <cstdint>
#include <exception>
#include <functional>
#include <memory>
#include <optional>
#include <span>
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

// Two sources and three destinations; deletions with full and with no reroute.
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

// The comparators below read a state through its accessors only, so they apply alike to the
// StateFields stand-in of their negative controls (the last test), whose fields can be changed
// one at a time.

// The values of two states: residual, potentials, value and the maximum flag.
template <typename State>
bool same_state(const State& a, const State& b) {
  return std::vector<Amount>(a.residual().begin(), a.residual().end()) ==
             std::vector<Amount>(b.residual().begin(), b.residual().end()) &&
         std::vector<detail::I128>(a.potentials().begin(), a.potentials().end()) ==
             std::vector<detail::I128>(b.potentials().begin(), b.potentials().end()) &&
         a.value() == b.value() && a.maximum() == b.maximum();
}

// Every field of two populated states: the graph and the network by identity (copies share
// them), the rest by value.
template <typename State>
bool identical(const State& a, const State& b) {
  return &a.graph() == &b.graph() && &a.network() == &b.network() && same_state(a, b) &&
         same_counters(a.counters(), b.counters());
}

// A state without fields: no residual, no potentials, value 0, not maximum.
template <typename State>
bool no_fields(const State& s) {
  return s.residual().empty() && s.potentials().empty() && s.value() == 0 && !s.maximum();
}

// A state that never held a solve: no fields, and no resolver accepts it.
bool empty_state(const DeletionState& s) {
  if (!no_fields(s)) return false;
  try {
    DeletionResolver rs(s);
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

// A state's fields, read through DeletionState's accessors and writable one at a time.
struct StateFields {
  const Graph* graph_ {nullptr};
  const detail::Network* network_ {nullptr};
  std::vector<Amount> residual_ {};
  std::vector<detail::I128> potentials_ {};
  Amount value_ {0};
  bool maximum_ {false};
  detail::FlowCoreCounters counters_ {};
  [[nodiscard]] const Graph& graph() const { return *graph_; }
  [[nodiscard]] const detail::Network& network() const { return *network_; }
  [[nodiscard]] std::span<const Amount> residual() const { return residual_; }
  [[nodiscard]] std::span<const detail::I128> potentials() const { return potentials_; }
  [[nodiscard]] Amount value() const { return value_; }
  [[nodiscard]] bool maximum() const { return maximum_; }
  [[nodiscard]] const detail::FlowCoreCounters& counters() const { return counters_; }
};

StateFields fields_of(const DeletionState& s) {
  return {&s.graph(),
          &s.network(),
          {s.residual().begin(), s.residual().end()},
          {s.potentials().begin(), s.potentials().end()},
          s.value(),
          s.maximum(),
          s.counters()};
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
Saved chain(std::int32_t n, Amount c) { return saved(chain_graph(n, c), {0}, {n - 1}); }

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
  // The same solve without a capture, its view allocated alike: the sweep below must also fail
  // the allocations of the captured state.
  long plain = 0;
  while (fails_at(plain, [&] {
    const std::vector<Amount> view(g->capacity_units().begin(), g->capacity_units().end());
    (void)place_global_units(*g, view, kSources, kDestinations, nullptr, std::nullopt);
  })) {
    ++plain;
  }
  int failures = 0;
  for (long at = 0;; ++at) {
    // Everything an attempt allocates is freed when it ends (the reference capture above grew
    // this thread's solver), so an allocation left live is a leak of the failed capture.
    const long live = g_live.load();
    bool failed = false;
    {
      DeletionCapture req;
      FlowCoreResult got;
      failed = fails_at(at, [&] { got = capture(*g, g, req); });
      if (!failed) {
        EXPECT_TRUE(same_payload(got, payload));
        ASSERT_TRUE(req.state.has_value());
        EXPECT_TRUE(same_state(*req.state, *want.state));
      } else {
        ++failures;
        EXPECT_FALSE(req.state.has_value()) << "allocation " << at;
        DeletionCapture retry;
        EXPECT_TRUE(same_payload(capture(*g, g, retry), payload)) << "retry after allocation " << at;
      }
    }
    EXPECT_EQ(g_live.load(), live) << "allocation " << at << ": a leak";
    if (!failed) break;
  }
  EXPECT_GT(failures, plain);
}

// A failed copy or resolver construction leaves nothing constructed and nothing allocated live;
// each sweep fails every allocation of its call (counted in the test below).
TEST(DeletionFaults, CopyAndResolverConstructionFailAtomically) {
  const auto g = graph();
  DeletionCapture req;
  (void)capture(*g, g, req);
  int copy_failures = 0;
  for (long at = 0;; ++at) {
    std::optional<DeletionState> copy;
    const long live = g_live.load();
    if (!fails_at(at, [&] { copy.emplace(*req.state); })) break;
    ++copy_failures;
    EXPECT_FALSE(copy.has_value()) << "copy, allocation " << at;
    EXPECT_EQ(g_live.load(), live) << "copy, allocation " << at << ": a leak";
  }
  int resolver_failures = 0;
  for (long at = 0;; ++at) {
    std::unique_ptr<DeletionResolver> rs;
    const long live = g_live.load();
    if (!fails_at(at, [&] { rs = std::make_unique<DeletionResolver>(*req.state); })) break;
    ++resolver_failures;
    EXPECT_EQ(rs, nullptr) << "resolver, allocation " << at;
    EXPECT_EQ(g_live.load(), live) << "resolver, allocation " << at << ": a leak";
  }
  EXPECT_EQ(copy_failures, 2);           // the residual and the potentials
  EXPECT_EQ(resolver_failures, 1 + 10);  // the resolver itself, then its construction (below)
}

// A resolver reads its state's network in place and copies none of its 14
// arrays. Its construction makes exactly the state copy's allocations (the residual and the
// potentials) plus 8 of its own: the solver, and its working arrays (levels, marks, current
// arcs, residual, wide potentials, wide distances, stamps). Counted on states of 4 to 13
// vertices; a failure at each of them constructs nothing (the sweep above).
TEST(DeletionFaults, AResolverAllocatesNoNetwork) {
  auto allocations = [](const std::function<void()>& fn) {
    long at = 0;
    while (fails_at(at, fn)) ++at;
    return at;
  };
  for (const auto& s : assignment_states()) {
    const DeletionState& st = s.state;
    const long copy = allocations([&] { DeletionState c(st); });
    const long resolver = allocations([&] { DeletionResolver rs(st); });
    EXPECT_EQ(copy, 2) << st.network().num_nodes;
    EXPECT_EQ(resolver, copy + 8) << st.network().num_nodes;
  }
}

// Every positive-flow edge, every allocation of a re-solve on a fresh resolver (its working
// storage grows) and on a warm one; the failed resolver serves the next call correctly.
TEST(DeletionFaults, EveryReSolveFailsWithoutPoisoningItsResolver) {
  const auto g = graph();
  DeletionCapture req;
  (void)capture(*g, g, req);
  const DeletionState saved = *req.state;
  const detail::Network network = saved.network();  // the shared network, as captured
  const auto flows = saved.edge_flows();
  int candidates = 0;
  for (EdgeId e = 0; e < g->num_edges(); ++e) {
    if (flows[static_cast<std::size_t>(e)] <= 0) continue;
    ++candidates;
    const auto want = DeletionResolver(saved).resolve(e);
    for (bool warm : {false, true}) {
      DeletionResolver kept(saved);
      if (warm) (void)kept.resolve(e);
      int failures = 0;
      for (long at = 0;; ++at) {
        // A fresh resolver ends with the attempt; the warm one already holds its grown storage.
        // So an allocation left live after the attempt is a leak of the failed re-solve.
        const long live = g_live.load();
        bool failed = false;
        {
          DeletionResolver fresh(saved);
          DeletionResolver& rs = warm ? kept : fresh;
          DeletionResult got;
          failed = fails_at(at, [&] { got = rs.resolve(e); });
          if (!failed) {
            EXPECT_TRUE(same_result(got, want)) << "edge " << e;
          } else {
            ++failures;
            EXPECT_TRUE(same_network(saved.network(), network)) << "edge " << e << " allocation " << at;
            EXPECT_TRUE(same_result(rs.resolve(e), want)) << "edge " << e << " retry after allocation " << at;
          }
        }
        EXPECT_EQ(g_live.load(), live) << "edge " << e << " warm " << warm << " allocation " << at << ": a leak";
        if (!failed) break;
      }
      // Every re-solve allocates (its result's arrays at least), warm or fresh.
      EXPECT_GE(failures, 3) << "edge " << e << " warm " << warm;
    }
  }
  EXPECT_GE(candidates, 5);
}

// A copy assignment between populated states of different coordinate sizes
// (and onto itself) fails at every allocation without changing its target; the pairs include a
// larger state assigned to a smaller one (the 6-node path to the 2-node path). A failed target is
// identical to its old value and serves its re-solves; a finished one is identical to the source.
TEST(DeletionFaults, CopyAssignmentFailsWithoutChangingItsTarget) {
  const auto states = assignment_states();
  for (std::size_t i = 0; i < states.size(); ++i) {
    for (std::size_t j = 0; j < states.size(); ++j) {
      int failed_here = 0;
      for (long at = 0;; ++at) {
        const std::string where = "target " + std::to_string(i) + " source " + std::to_string(j) + " allocation " +
                                  std::to_string(at);
        const long live = g_live.load();
        bool failed = false;
        {
          DeletionState target = states[i].state;
          failed = fails_at(at, [&] {
            if (i == j) {
              const DeletionState& self = target;
              target = self;
            } else {
              target = states[j].state;
            }
          });
          if (!failed) {
            EXPECT_TRUE(holds(target, states[j])) << where;
          } else {
            ++failed_here;
            EXPECT_TRUE(holds(target, states[i])) << where;
          }
        }
        EXPECT_EQ(g_live.load(), live) << where << ": a leak";  // the copy built aside is gone
        if (!failed) break;
        ASSERT_FALSE(::testing::Test::HasFailure()) << where;
      }
      // The copy built aside allocates the residual and the potentials, whatever the sizes.
      EXPECT_EQ(failed_here, 2) << "target " << i << " source " << j;
    }
  }
}

// An empty state (no solve) as the target keeps nothing after a failure and accepts no resolver;
// as the source it empties the target without allocating.
TEST(DeletionFaults, CopyAssignmentWithAnEmptyStateIsAllOrNothing) {
  const auto states = assignment_states();
  for (std::size_t j = 0; j < states.size(); ++j) {
    long failed = 0;
    for (long at = 0;; ++at) {
      const long live = g_live.load();
      bool done = false;
      {
        DeletionState target;
        done = !fails_at(at, [&] { target = states[j].state; });
        if (done) {
          EXPECT_TRUE(holds(target, states[j])) << "source " << j;
        } else {
          ++failed;
          EXPECT_TRUE(empty_state(target)) << "source " << j << " allocation " << at;
        }
      }
      EXPECT_EQ(g_live.load(), live) << "source " << j << " allocation " << at << ": a leak";
      if (done) break;
    }
    // The copy built aside allocates the residual and the potentials (as a copy construction).
    EXPECT_EQ(failed, 2) << "source " << j;
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
// whose graph and network the next thread's target, when there is one, shares (a copy of the same
// state); every allocation of the run fails once (any thread's, or the workers' only). A target
// whose assignment failed or never started keeps its old value, every other one holds its source;
// the runs that failed are those of the threads' own allocations plus two per assignment.
TEST(DeletionFaults, CopyAssignmentsOnWorkerThreadsAreAllOrNothing) {
  const auto states = assignment_states();
  const auto k = states.size();
  // Runs body(t) on `workers` threads (on the caller for one), then raises the first error.
  const auto on_threads = [](std::size_t workers, const std::function<void(std::size_t)>& body,
                             const std::vector<std::exception_ptr>& errors) {
    if (workers == 1) {
      body(0);
    } else {
      std::vector<std::thread> pool;
      pool.reserve(workers);
      try {
        for (std::size_t t = 0; t < workers; ++t) pool.emplace_back([&body, t] { body(t); });  // not a copy of body
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
  for (std::size_t workers : {1, 2, 4}) {
    for (bool workers_only : {false, true}) {
      if (workers == 1 && workers_only) continue;
      // The allocations of the threads alone (their bodies empty).
      const std::vector<std::exception_ptr> none(workers);
      const std::function<void(std::size_t)> empty = [](std::size_t) {};
      const auto idle = [&] { on_threads(workers, empty, none); };
      long threads_only = 0;
      while (workers_only ? fails_in_workers_at(threads_only, idle) : fails_at(threads_only, idle)) ++threads_only;
      long failures = 0;
      for (long at = 0;; ++at) {
        // The targets, the threads and the copies built aside end with the attempt: an
        // allocation left live is a leak of a failed assignment.
        const long live = g_live.load();
        bool failed = false;
        {
          std::vector<DeletionState> targets;
          for (std::size_t t = 0; t < workers; ++t) targets.push_back(states[t % k].state);
          std::vector<std::uint8_t> done(workers, 0);
          std::vector<std::exception_ptr> errors(workers);
          // Wrapped before the sweep's window (the wrapper may allocate).
          const std::function<void(std::size_t)> body = [&](std::size_t t) noexcept {
            try {
              targets[t] = states[(t + 1) % k].state;
              done[t] = 1;
            } catch (...) {
              errors[t] = std::current_exception();
            }
          };
          auto run = [&] { on_threads(workers, body, errors); };
          failed = workers_only ? fails_in_workers_at(at, run) : fails_at(at, run);
          for (std::size_t t = 0; t < workers; ++t) {
            const auto& want = done[t] != 0 ? states[(t + 1) % k] : states[t % k];
            EXPECT_TRUE(holds(targets[t], want)) << "workers " << workers << " thread " << t << " allocation " << at;
          }
          if (!failed) {
            for (std::size_t t = 0; t < workers; ++t) EXPECT_EQ(done[t], 1) << "workers " << workers << " thread " << t;
          }
        }
        EXPECT_EQ(g_live.load(), live) << "workers " << workers << " allocation " << at << ": a leak";
        if (!failed) break;
        ++failures;
        ASSERT_FALSE(::testing::Test::HasFailure());
        ASSERT_LT(at, 1000);
      }
      // Each assignment allocates twice (the copy built aside), and every allocation failed once.
      EXPECT_EQ(failures, threads_only + 2 * static_cast<long>(workers))
          << "workers " << workers << " workers_only " << workers_only;
    }
  }
}

// Every allocation of the capture on a thread whose solver is cold, and on one whose solver first
// served a smaller graph (its storage grows), leaves no state and an unchanged payload; the retry
// on the same thread gives the reference payload and state. The sweep's failures outnumber the
// allocations of a capture on a thread whose solver already served this graph: the sweep reached
// the solver's growth (sweep_scratch_growth).
TEST(DeletionFaults, CaptureFailsAtEveryAllocationOnAColdOrGrowingThread) {
  const auto g = graph();
  DeletionCapture ref;
  const auto payload = capture(*g, g, ref);
  for (bool grow : {false, true}) {
    sweep_scratch_growth(
        grow ? "growing thread" : "cold thread",
        [&] {
          if (grow) (void)chain(2, 2);
        },
        [&](long at, const std::string& w) {
          DeletionCapture req;
          FlowCoreResult got;
          const bool failed = fails_at(at, [&] { got = capture(*g, g, req); });
          if (failed) {
            EXPECT_FALSE(req.state.has_value()) << w;
            got = capture(*g, g, req);
          }
          EXPECT_TRUE(same_payload(got, payload)) << w;
          EXPECT_TRUE(req.state.has_value() && same_state(*req.state, *ref.state)) << w;
          return failed;
        },
        [&] {
          DeletionCapture first;
          (void)capture(*g, g, first);  // the solver's storage grows to this graph
          long warm = 0;
          while (fails_at(warm, [&] {
            DeletionCapture again;
            (void)capture(*g, g, again);
          })) {
            ++warm;
          }
          return warm;
        });
  }
}

// Global sensitivity at every allocation (any thread's, and the workers' only) on 1, 2 and 4
// workers; after each failure a retry gives the reference baseline and every entry.
TEST(DeletionFaults, GlobalSensitivityFailsAtEveryAllocation) {
  const auto g = std::make_shared<const Graph>(
      make_graph(5, {{0, 1, 2, 1}, {1, 2, 2, 1}, {1, 3, 1, 2}, {3, 2, 1, 2}, {2, 4, 2, 1}, {0, 3, 1, 5}}));
  const QueryOptions opts {global_config(), nullptr, nullptr, OutputOptions {true, true, true}};
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
    }
  };
  for (std::int32_t workers : {1, 2, 4}) {
    for (bool workers_only : {false, true}) {
      if (workers == 1 && workers_only) continue;
      int failures = 0;
      for (long at = 0;; ++at) {
        const std::string where = "workers " + std::to_string(workers) + (workers_only ? " (workers only)" : "") +
                                  " allocation " + std::to_string(at);
        // The result and the retry end with the attempt, the workers' scratch with their threads,
        // and the reference above grew this thread's: an allocation left live is a leak.
        const long live = g_live.load();
        bool failed = false;
        {
          SensitivityResult got;
          const auto call = [&] { got = sensitivity(*g, sources, destinations, opts, workers); };
          failed = workers_only ? fails_in_workers_at(at, call) : fails_at(at, call);
          if (!failed) {
            same(got, where);
          } else {
            same(sensitivity(*g, sources, destinations, opts, workers), "retry after " + where);
          }
        }
        EXPECT_EQ(g_live.load(), live) << where << ": a leak";
        if (!failed) break;
        ++failures;
        ASSERT_FALSE(::testing::Test::HasFailure());
        ASSERT_LT(at, 10000);
      }
      EXPECT_GT(failures, workers_only ? 10 : 100) << "workers " << workers << " workers_only " << workers_only;
    }
  }
}

// The negative controls of the comparators above: a copy of a captured state's fields with any
// one field changed is not identical, and with any value field changed not the same state
// either; a vector field is changed both by one more entry and by its last entry's value (so a
// comparator of sizes only, or of a prefix only, fails). A state without fields has none, and
// each field set makes it one with a field. A state serves its own re-solves; another state
// (also one of equal values on another graph), a changed reference result or an empty list of
// results does not.
TEST(DeletionFaults, TheStateComparatorsCompareEveryField) {
  const Saved six = chain(6, 2);
  const Saved other = chain(6, 3);
  const StateFields st = fields_of(six.state);
  const StateFields elsewhere = fields_of(other.state);
  ASSERT_TRUE(identical(st, st));
  ASSERT_TRUE(st.maximum_ && st.value_ > 0 && !st.residual_.empty() && !st.potentials_.empty());
  ASSERT_NE(st.graph_, elsewhere.graph_);
  ASSERT_NE(st.network_, elsewhere.network_);
  const std::vector<std::function<void(StateFields&)>> values = {
      [](auto& x) { x.residual_.push_back(0); },
      [](auto& x) { ++x.residual_.back(); },
      [](auto& x) { x.potentials_.push_back(detail::I128 {0}); },
      [](auto& x) { x.potentials_.back() += detail::I128 {1}; },
      [](auto& x) { ++x.value_; },
      [](auto& x) { x.maximum_ = !x.maximum_; },
  };
  for (std::size_t k = 0; k < values.size(); ++k) {
    StateFields x = st;
    values[k](x);
    EXPECT_FALSE(same_state(x, st)) << "value field " << k;
    EXPECT_FALSE(identical(x, st)) << "value field " << k;
  }
  const std::vector<std::function<void(StateFields&)>> others = {
      [&](auto& x) { x.graph_ = elsewhere.graph_; },
      [&](auto& x) { x.network_ = elsewhere.network_; },
      [](auto& x) { ++x.counters_.dijkstra_runs; },
      [](auto& x) { ++x.counters_.augmentations; },
  };
  for (std::size_t k = 0; k < others.size(); ++k) {
    StateFields x = st;
    others[k](x);
    EXPECT_TRUE(same_state(x, st)) << "identity or counter field " << k;
    EXPECT_FALSE(identical(x, st)) << "identity or counter field " << k;
  }
  const StateFields none {};
  EXPECT_TRUE(no_fields(none));
  for (std::size_t k = 0; k < values.size(); ++k) {
    if (k == 1 || k == 3) continue;  // the last entry's value: none has no entry
    StateFields x = none;
    values[k](x);
    EXPECT_FALSE(no_fields(x)) << "value field " << k;
  }
  EXPECT_TRUE(empty_state(DeletionState {}));
  EXPECT_FALSE(empty_state(six.state));
  EXPECT_TRUE(holds(six.state, six));
  EXPECT_FALSE(holds(other.state, six));
  const Saved twin = chain(6, 2);  // equal values on another graph and network
  ASSERT_TRUE(same_state(twin.state, six.state));
  EXPECT_FALSE(holds(twin.state, six));
  Saved changed = six;
  ASSERT_FALSE(changed.results.empty());
  ++changed.results.back().second.value;
  EXPECT_FALSE(holds(six.state, changed));
  Saved no_results = six;
  no_results.results.clear();
  EXPECT_FALSE(holds(six.state, no_results));
}
