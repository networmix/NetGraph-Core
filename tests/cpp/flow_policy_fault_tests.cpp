// Allocation-failure atomicity of the policy engine, the DagCache and the evaluators (APPROACH 8,
// 9.2, 10, 11, 13). Every path that grows a buffer or publishes state is swept over every
// allocation it makes: after each injected std::bad_alloc the ledger (aggregate, residual, every
// record field, owners, next FlowId) and the policy's epoch state (requested, +inf, cohorts with
// their ids and plans) equal their state before the call; a clean retry gives the reference
// result (no poisoned per-thread scratch); and remove() restores the ledger. A supplied building
// cache is a memo (APPROACH 13, limit 2): it may keep the complete entries of the lookups that
// succeeded before the failure, never a partial one. Runs in the fault-injection executable
// (fault_injection.hpp).
#include <gtest/gtest.h>

#include <algorithm>
#include <functional>
#include <memory>
#include <new>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "fault_injection.hpp"
#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/traversal.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using detail::LedgerAccess;

namespace {

// The shared comparators of fault_injection.hpp, beside this file's own overloads.
using test::expect_same;

// The policy's epoch state: the call state and per roster source its flows (ids, routes, members).
struct Epoch {
  Counter128 requested;
  bool infinite {false};
  std::vector<std::vector<std::pair<std::int64_t, std::vector<EdgeId>>>> cohorts;
};

Epoch epoch(const FlowPolicy& p) {
  const auto& pe = detail::PolicyAccess::engine(p);
  Epoch e {pe.requested, pe.infinite, {}};
  for (const auto& c : pe.cohorts) {
    auto& out = e.cohorts.emplace_back();
    for (const auto& f : c) {
      if (f.plan) {
        out.emplace_back(f.id, f.plan->members);
      } else {
        const auto edges = detail::PredDAGBuilder::ascending(*f.bundle);
        out.emplace_back(f.id, std::vector<EdgeId>(edges.begin(), edges.end()));
      }
    }
  }
  return e;
}

void expect_same(const Epoch& a, const Epoch& b, const std::string& where) {
  EXPECT_EQ(a.requested, b.requested) << where;
  EXPECT_EQ(a.infinite, b.infinite) << where;
  EXPECT_EQ(a.cohorts, b.cohorts) << where;
}

PathSelection sel(EdgeTieBreak t, EdgePriority p, bool multipath = true, bool multi_edge = true) {
  return PathSelection(t, p, multipath, multi_edge);
}

// Sources 0 and 1, destinations 5 and 6: three equal-cost routes with a parallel pair and a
// shared middle, a costlier bypass, and a second destination.
std::shared_ptr<const Graph> graph() {
  return std::make_shared<const Graph>(
      make_graph(7, {{0, 2, 6, 1}, {0, 3, 4, 1}, {1, 3, 5, 1}, {1, 4, 3, 1}, {2, 5, 3, 1}, {2, 5, 2, 1},
                     {3, 5, 4, 1}, {4, 6, 9, 1}, {3, 6, 2, 2}, {0, 6, 3, 4}, {2, 4, 2, 1}}));
}

struct Mode {
  std::string name;
  FlowPolicyConfig cfg;
  std::vector<NodeId> sources;
  std::vector<NodeId> destinations;
};

std::vector<Mode> modes(const Graph& g) {
  std::vector<Mode> out;
  const auto det = sel(EdgeTieBreak::Deterministic, EdgePriority::Capacity);
  const auto res = sel(EdgeTieBreak::Deterministic, EdgePriority::Residual);
  const auto tunnel = sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, false, false);
  const std::vector<NodeId> S = {0, 1};
  const std::vector<NodeId> D = {5, 6};
  out.push_back({"costonly maxflow", {CostOnly {}, FlowPlacement::MaxFlow, Admission::Lossless, det}, S, D});
  out.push_back({"costonly eb", {CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless, det}, S, D});
  out.push_back({"costonly weighted residual", {CostOnly {}, FlowPlacement::Weighted, Admission::Lossless, res}, S, D});
  out.push_back({"costonly lossy", {CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossy, det}, S, D});
  out.push_back({"costonly single path collapse",
                 {CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless, tunnel}, S, D});
  CapacityAware grow;
  grow.flow_count = 1;
  grow.max_flow_count = 4;
  out.push_back({"capacityaware eb growth", {grow, FlowPlacement::EqualBalanced, Admission::Lossless, tunnel}, S, D});
  CapacityAware reopt;
  reopt.flow_count = 2;
  reopt.reoptimize = true;
  reopt.max_path_cost = Stretch {2.0};
  out.push_back({"capacityaware maxflow reopt", {reopt, FlowPlacement::MaxFlow, Admission::Lossless, det}, S, D});
  CapacityAware dim;
  dim.flow_count = 1;
  dim.diminishing = DiminishingReturns {1, 0.5};
  dim.min_flow_bandwidth = 1;
  out.push_back({"capacityaware weighted diminishing", {dim, FlowPlacement::Weighted, Admission::Lossless, res}, S, D});
  out.push_back({"global", {Global {}, FlowPlacement::MaxFlow, Admission::Lossless, det}, S, D});
  Pinned two;
  const auto inv = id_of_input(g);
  two.bundles = {std::make_shared<const PredDAG>(detail::PredDAGBuilder::of_edges(g, std::vector<EdgeId> {inv[0], inv[4], inv[5]})),
                 std::make_shared<const PredDAG>(detail::PredDAGBuilder::of_edges(g, std::vector<EdgeId> {inv[1], inv[6]}))};
  out.push_back({"pinned eb cohort", {two, FlowPlacement::EqualBalanced, Admission::Lossless, det}, {0}, {5}});
  out.push_back({"pinned lossy", {two, FlowPlacement::Weighted, Admission::Lossy, det}, {0}, {5}});
  // CostOnly creation in one pass (task 1d-msd): its Collapse branch (multipath, multi_edge false,
  // on the creation view) and its Weighted weights (Capacity priority).
  out.push_back({"costonly eb multipath collapse",
                 {CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless,
                  sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, true, false)},
                 S, D});
  out.push_back({"costonly weighted capacity", {CostOnly {}, FlowPlacement::Weighted, Admission::Lossless, det}, S, D});
  return out;
}

struct Scenario {
  std::vector<std::optional<Amount>> before;
  std::optional<Amount> call;
  OutputOptions outputs;
  std::string name;
};

std::vector<Scenario> scenarios() {
  return {
      {{}, 5, {}, "create"},
      {{3}, 4, {}, "update"},
      {{2}, std::nullopt, {}, "update to +inf"},
      {{}, 7, {true, true, true}, "create with outputs"},
      {{4}, 0, {true, true, true}, "zero with outputs"},
  };
}

// A third party holding one unit on the first route's edges.
std::unique_ptr<FlowPolicy> third_party(FlowLedger& fg, const Graph& g) {
  const auto inv = id_of_input(g);
  Pinned p;
  p.bundles = {std::make_shared<const PredDAG>(detail::PredDAGBuilder::of_edges(g, std::vector<EdgeId> {inv[0], inv[4]}))};
  auto t = std::make_unique<FlowPolicy>(fg, std::vector<NodeId> {0}, std::vector<NodeId> {5},
                                        FlowPolicyConfig {p, FlowPlacement::MaxFlow, Admission::Lossless,
                                                          PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity)});
  (void)t->place(1);
  return t;
}

} // namespace

TEST(PolicyFaults, EveryAllocationOfEveryRoutingModeFailsAtomically) {
  const auto g = graph();
  int sweeps = 0;
  long injected = 0;
  for (const auto& m : modes(*g)) {
    for (const auto& sc : scenarios()) {
      const std::string where = m.name + " / " + sc.name;
      FlowLedger ref(g);
      auto ref_third = third_party(ref, *g);
      FlowPolicy pref(ref, m.sources, m.destinations, m.cfg);
      for (auto v : sc.before) (void)pref.place(v);
      const LedgerState initial = ledger_state(ref);
      const Epoch initial_epoch = epoch(pref);
      const PlacementResult want = pref.place(sc.call, sc.outputs);
      const LedgerState after = ledger_state(ref);
      const Epoch want_epoch = epoch(pref);
      pref.remove();
      const LedgerState removed = ledger_state(ref);
      ++sweeps;
      for (long at = 0;; ++at) {
        const std::string w = where + " allocation " + std::to_string(at);
        FlowLedger fg(g);
        auto third = third_party(fg, *g);
        FlowPolicy p(fg, m.sources, m.destinations, m.cfg);
        for (auto v : sc.before) (void)p.place(v);
        expect_same(ledger_state(fg), initial, w + " (setup)");
        PlacementResult got;
        if (!fails_at(at, [&] { got = p.place(sc.call, sc.outputs); })) {
          expect_same(got, want, w + " (no failure)");
          expect_same(ledger_state(fg), after, w + " (no failure)");
          break;
        }
        ++injected;
        expect_same(ledger_state(fg), initial, w + " (failed call)");
        expect_same(epoch(p), initial_epoch, w + " (failed call)");
        expect_balanced(fg, w + " (failed call)");
        expect_same(p.place(sc.call, sc.outputs), want, w + " (retry)");
        expect_same(ledger_state(fg), after, w + " (retry)");
        expect_same(epoch(p), want_epoch, w + " (retry)");
        p.remove();
        expect_same(ledger_state(fg), removed, w + " (remove)");
        if (::testing::Test::HasFailure()) return;
      }
    }
  }
  EXPECT_EQ(sweeps, 13 * 5);
  EXPECT_GT(injected, sweeps * 3);  // not vacuous
}

// The calls that place nothing publish only their outputs and, for a positive request, the
// epoch's request: a zero request, and a positive or +inf request when no source reaches a
// destination (an empty roster; the endpoints swapped: 5 and 6 have no out-edges). A failure
// changes neither the ledger nor the epoch state.
TEST(PolicyFaults, ZeroAndEmptyRosterCallsFailAtomically) {
  const auto g = graph();
  const OutputOptions all {true, true, true};
  long injected = 0;
  struct Case {
    bool swap;
    std::optional<Amount> volume;
    std::string name;
  };
  const std::vector<Case> cases = {{false, 0, "zero"}, {true, 4, "empty roster"}, {true, std::nullopt, "empty roster +inf"}};
  for (const auto& m : modes(*g)) {
    if (std::holds_alternative<Pinned>(m.cfg.routing)) continue;  // bundles bind their endpoints
    for (const auto& c : cases) {
      const auto& s = c.swap ? m.destinations : m.sources;
      const auto& d = c.swap ? m.sources : m.destinations;
      FlowLedger ref(g);
      auto ref_third = third_party(ref, *g);
      FlowPolicy pref(ref, s, d, m.cfg);
      const PlacementResult want = pref.place(c.volume, all);
      const Epoch want_epoch = epoch(pref);
      if (c.swap) EXPECT_TRUE(detail::PolicyAccess::engine(pref).roster.empty()) << m.name;
      for (long at = 0;; ++at) {
        const std::string w = m.name + " / " + c.name + " allocation " + std::to_string(at);
        FlowLedger fg(g);
        auto third = third_party(fg, *g);
        FlowPolicy p(fg, s, d, m.cfg);
        const LedgerState before = ledger_state(fg);
        const Epoch before_epoch = epoch(p);
        PlacementResult got;
        if (!fails_at(at, [&] { got = p.place(c.volume, all); })) {
          expect_same(got, want, w + " (no failure)");
          break;
        }
        ++injected;
        expect_same(ledger_state(fg), before, w + " (failed call)");
        expect_same(epoch(p), before_epoch, w + " (failed call)");
        expect_same(p.place(c.volume, all), want, w + " (retry)");
        expect_same(epoch(p), want_epoch, w + " (retry)");
        if (::testing::Test::HasFailure()) return;
      }
    }
  }
  EXPECT_GT(injected, 27);  // every mode and case fails at least once
}

TEST(PolicyFaults, ConstructionFailsWithoutBindingAndTheCacheStaysConsistent) {
  const auto g = graph();
  for (const auto& m : modes(*g)) {
    for (const bool with_cache : {false, true}) {
      for (long at = 0;; ++at) {
        FlowLedger fg(g);
        DagCache cache(*g);
        const LedgerState before = ledger_state(fg);
        std::optional<FlowPolicy> p;
        if (!fails_at(at, [&] { p.emplace(fg, m.sources, m.destinations, m.cfg, nullptr, with_cache ? &cache : nullptr); })) {
          EXPECT_EQ(LedgerAccess::bound_owners(fg).size(), 1U);
          break;
        }
        const std::string w = m.name + (with_cache ? " cached" : "") + " construction allocation " + std::to_string(at);
        EXPECT_FALSE(p.has_value()) << w;
        expect_same(ledger_state(fg), before, w);
        // The cache holds whole entries only: its statistics agree with a fresh lookup.
        const auto s = cache.stats();
        EXPECT_EQ(s.live_values, s.entries) << w;
        EXPECT_LE(s.entries, 1) << w;
        // A clean construction works and places like the reference.
        p.emplace(fg, m.sources, m.destinations, m.cfg, nullptr, with_cache ? &cache : nullptr);
        EXPECT_GE(p->place(2).placed, 0) << w;
        expect_balanced(fg, w);
      }
    }
  }
}

TEST(PolicyFaults, DagCacheInsertionIsAllOrNothing) {
  const auto g = graph();
  const std::vector<NodeId> roots = {5, 6};
  for (long at = 0;; ++at) {
    DagCache cache(*g);
    (void)cache.lookup(*g, std::vector<NodeId> {5}, Direction::Reverse, nullptr);  // one entry before
    const auto before = cache.stats();
    std::shared_ptr<const PredDAG> got;
    if (!fails_at(at, [&] { got = cache.lookup(*g, roots, Direction::Reverse, nullptr); })) {
      EXPECT_EQ(cache.stats().entries, 2);
      break;
    }
    const auto s = cache.stats();
    const std::string w = "lookup allocation " + std::to_string(at);
    EXPECT_EQ(s.entries, before.entries) << w;
    EXPECT_EQ(s.bytes, before.bytes) << w;
    EXPECT_EQ(s.live_values, before.live_values) << w;
    EXPECT_EQ(s.hits, before.hits) << w;
    // The key is the lookup's first allocation: its failure counts nothing; a failed search (or
    // insertion) after it counts its miss.
    EXPECT_EQ(s.misses, before.misses + (at == 0 ? 0 : 1)) << w;
    const auto v = cache.lookup(*g, roots, Direction::Reverse, nullptr);  // a clean retry inserts
    EXPECT_EQ(cache.stats().entries, 2) << w;
    EXPECT_TRUE(*v == detail::static_support(*g, roots, Direction::Reverse, nullptr)) << w;
  }
}

namespace {

// The evaluators under a failure on `workers` threads: anywhere in the process (every `step`-th
// allocation: the evaluators allocate a lot) or, with several workers, in the worker threads only
// (their setup and their tasks; every third worker allocation, from the first). The failure
// propagates as std::bad_alloc after every worker has joined, never terminating the process
// (review-1c F1: a worker's cache allocation escaped its thread), and a clean retry gives the
// reference. Returns the number of injected failures.
long sweep_evaluator(const std::string& name, std::int32_t workers, bool in_workers, long step,
                     const std::function<void(std::int32_t)>& call, const std::function<void(std::int32_t)>& retry) {
  long failures = 0;
  for (long at = 0;; at += in_workers ? 3 : step) {
    const auto fn = [&] { call(workers); };
    if (!(in_workers ? fails_in_workers_at(at, fn) : fails_at(at, fn))) break;
    ++failures;
    retry(workers);
    if (::testing::Test::HasFailure()) {
      ADD_FAILURE() << name << " workers " << workers << (in_workers ? " (worker threads)" : "") << " allocation " << at;
      break;
    }
  }
  return failures;
}

void sweep_batch(const Mode& m, const Graph& g, std::span<const std::int32_t> workers_set) {
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {m.sources, m.destinations}, {{1}, {6}}, {{0}, {5}}};
  const auto want = batch_capacity(g, pairs, CapacityOptions {m.cfg}, {}, 1);
  for (const std::int32_t workers : workers_set) {
    for (const bool in_workers : {false, true}) {
      if (in_workers && workers == 1) continue;  // one worker runs on the calling thread
      const long failures = sweep_evaluator(
          m.name + " batch", workers, in_workers, 7,
          [&](std::int32_t w) { (void)batch_capacity(g, pairs, CapacityOptions {m.cfg}, {}, w); },
          [&](std::int32_t w) {
            const auto b = batch_capacity(g, pairs, CapacityOptions {m.cfg}, {}, w);
            ASSERT_EQ(b.size(), want.size());
            for (std::size_t i = 0; i < b.size(); ++i) expect_same(b[i], want[i], m.name + " batch retry");
          });
      EXPECT_GT(failures, 2) << m.name << " batch workers " << workers;  // not vacuous
    }
  }
}

void sweep_sensitivity(const Mode& m, const Graph& g, std::span<const std::int32_t> workers_set) {
  const SensitivityOptions opts {m.cfg};
  const auto want = sensitivity(g, m.sources, m.destinations, opts, 1);
  for (const std::int32_t workers : workers_set) {
    for (const bool in_workers : {false, true}) {
      if (in_workers && workers == 1) continue;
      const long failures = sweep_evaluator(
          m.name + " sensitivity", workers, in_workers, 11,
          [&](std::int32_t w) { (void)sensitivity(g, m.sources, m.destinations, opts, w); },
          [&](std::int32_t w) {
            const auto s = sensitivity(g, m.sources, m.destinations, opts, w);
            ASSERT_EQ(s.entries.size(), want.entries.size());
            for (std::size_t i = 0; i < s.entries.size(); ++i) {
              EXPECT_EQ(s.entries[i].edge, want.entries[i].edge) << m.name;
              EXPECT_EQ(s.entries[i].loss, want.entries[i].loss) << m.name;
            }
          });
      EXPECT_GT(failures, 2) << m.name << " sensitivity workers " << workers;
    }
  }
}

} // namespace

TEST(PolicyFaults, EvaluatorsLeaveNoPoisonedScratch) {
  const auto g = graph();
  const std::int32_t one[] = {1};
  for (const auto& m : modes(*g)) {
    const CapacityOptions copts {m.cfg, nullptr, nullptr, OutputOptions {true, true, true}};
    const auto want = capacity(*g, m.sources, m.destinations, copts);
    for (long at = 0;; ++at) {
      PlacementResult got;
      if (!fails_at(at, [&] { got = capacity(*g, m.sources, m.destinations, copts); })) {
        expect_same(got, want, m.name + " capacity (no failure)");
        break;
      }
      expect_same(capacity(*g, m.sources, m.destinations, copts), want, m.name + " capacity retry " + std::to_string(at));
    }
    if (std::holds_alternative<Pinned>(m.cfg.routing)) continue;
    sweep_batch(m, *g, one);
    sweep_sensitivity(m, *g, one);
    if (::testing::Test::HasFailure()) return;
  }
}

TEST(PolicyFaults, BatchWorkersFailWithoutTerminating) {
  const auto g = graph();
  const std::int32_t several[] = {2, 4};
  for (const auto& m : modes(*g)) {
    if (std::holds_alternative<Pinned>(m.cfg.routing)) continue;
    sweep_batch(m, *g, several);
    if (::testing::Test::HasFailure()) return;
  }
}

TEST(PolicyFaults, SensitivityWorkersFailWithoutTerminating) {
  const auto g = graph();
  const std::int32_t several[] = {2, 4};
  for (const auto& m : modes(*g)) {
    if (std::holds_alternative<Pinned>(m.cfg.routing)) continue;
    sweep_sensitivity(m, *g, several);
    if (::testing::Test::HasFailure()) return;
  }
}

TEST(PolicyFaults, RemovalAndReleaseNeverAllocate) {
  const auto g = graph();
  for (const auto& m : modes(*g)) {
    FlowLedger fg(g);
    auto p = std::make_unique<FlowPolicy>(fg, m.sources, m.destinations, m.cfg);
    (void)p->place(5);
    EXPECT_FALSE(fails_at(0, [&] { p->remove(); })) << m.name;
    EXPECT_TRUE(fg.records().empty()) << m.name;
    (void)p->place(3);
    const auto held = fg.allocated();
    EXPECT_FALSE(fails_at(0, [&] { p.reset(); })) << m.name;
    EXPECT_EQ(fg.allocated(), held) << m.name;
    EXPECT_TRUE(LedgerAccess::bound_owners(fg).empty()) << m.name;
  }
}

// The working-residual overlay and the per-source scratch grow with the graph: every position
// fails on a graph larger than any before, so the failure lands while they grow; the failed call
// changes nothing and a clean retry on the same graph works (ASan reports any access outside a
// partly grown array).
TEST(PolicyFaults, ScratchGrowthFailureLeavesNothingHalfGrown) {
  int failures = 0;
  std::int32_t n = 400;
  for (long at = 0;; ++at, n += 300) {
    std::vector<InEdge> chain;
    for (NodeId u = 0; u + 1 < n; ++u) chain.push_back({u, u + 1, 5, 1});
    auto g = std::make_shared<const Graph>(make_graph(n, chain));
    FlowLedger fg(g);
    FlowPolicy p(fg, std::vector<NodeId> {0}, std::vector<NodeId> {n - 1},
                 FlowPolicyConfig {CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless,
                                   PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity)});
    const LedgerState before = ledger_state(fg);
    PlacementResult got;
    const bool failed = fails_at(at, [&] { got = p.place(2, OutputOptions {true, true, true}); });
    if (failed) {
      ++failures;
      expect_same(ledger_state(fg), before, "growth allocation " + std::to_string(at));
      got = p.place(2, OutputOptions {true, true, true});
    }
    EXPECT_EQ(got.placed, 2) << "growth allocation " << at;
    expect_balanced(fg, "growth allocation " + std::to_string(at));
    if (!failed) break;
  }
  EXPECT_GT(failures, 4);
}

TEST(PolicyFaults, FrozenLookupsFailWithoutCounting) {
  // A frozen hit allocates only the key of roots that are not already a sorted set (roots that are
  // one are borrowed: nothing; task 1d-msd), a frozen miss its key and the DagCacheMiss message,
  // all before counting: a failure there counts nothing and changes no entry; a clean retry counts.
  const auto g = graph();
  DagCache cache(*g);
  const std::vector<NodeId> unsorted = {6, 5};
  const std::vector<NodeId> sorted = {5, 6};
  const std::vector<NodeId> absent = {5};
  const auto value = cache.lookup(*g, unsorted, Direction::Reverse, nullptr);
  cache.freeze();
  for (const int kind : {0, 1, 2}) {  // a hit on unsorted roots, a hit on a sorted set, a miss
    const bool hit = kind != 2;
    const auto& present = kind == 0 ? unsorted : sorted;
    long failures = 0;
    for (long at = 0;; ++at) {
      const auto before = cache.stats();
      bool missed = false;
      const auto call = [&] {
        try {
          if (cache.lookup(*g, hit ? present : absent, Direction::Reverse, nullptr).get() != value.get()) ADD_FAILURE();
        } catch (const DagCacheMiss&) {
          missed = true;
        }
      };
      const bool failed = fails_at(at, call);
      const auto s = cache.stats();
      const std::string w = std::string(hit ? "hit" : "miss") + " allocation " + std::to_string(at);
      EXPECT_EQ(s.entries, 1) << w;
      EXPECT_EQ(s.bytes, before.bytes) << w;
      EXPECT_EQ(s.live_values, before.live_values) << w;
      if (failed) {
        ++failures;
        EXPECT_EQ(s.hits, before.hits) << w;
        EXPECT_EQ(s.misses, before.misses) << w;
        continue;
      }
      EXPECT_EQ(missed, !hit) << w;
      EXPECT_EQ(s.hits - before.hits, hit ? 1 : 0) << w;
      EXPECT_EQ(s.misses - before.misses, hit ? 0 : 1) << w;
      break;
    }
    if (kind == 1) {
      EXPECT_EQ(failures, 0) << "hit on a sorted set";  // borrowed: no allocation at all
    } else {
      EXPECT_GE(failures, hit ? 1 : 2) << (hit ? "hit" : "miss");  // the key; the message (and a key)
    }
  }
}

namespace {

void expect_stats(const DagCache::Stats& got, const DagCache::Stats& want, const std::string& where) {
  EXPECT_EQ(got.hits, want.hits) << where;
  EXPECT_EQ(got.misses, want.misses) << where;
  EXPECT_EQ(got.entries, want.entries) << where;
  EXPECT_EQ(got.bytes, want.bytes) << where;
  EXPECT_EQ(got.live_values, want.live_values) << where;
}

// One public call reading a supplied frozen cache, under every injected failure (every `step`-th
// allocation in any thread, or every third in the worker threads only): a call failing with
// std::bad_alloc leaves all five statistics as they were; the call that completes adds exactly
// `hits` hits, or, when its outcome is the DagCacheMiss (`miss`), exactly one miss and no hit
// (NOTES R27). Returns the number of injected failures.
long sweep_frozen_stats(const std::string& name, const DagCache& cache, bool in_workers, long step, std::int64_t hits,
                        bool miss, const std::function<void()>& call) {
  long failures = 0;
  for (long at = 0;; at += in_workers ? 3 : step) {
    const auto before = cache.stats();
    bool missed = false;
    const auto fn = [&] {
      try {
        call();
      } catch (const DagCacheMiss&) {
        missed = true;
      }
    };
    const bool failed = in_workers ? fails_in_workers_at(at, fn) : fails_at(at, fn);
    const std::string w = name + " allocation " + std::to_string(at) + (in_workers ? " (worker threads)" : "");
    auto want = before;
    if (failed) {
      ++failures;
      EXPECT_FALSE(missed) << w;
    } else {
      EXPECT_EQ(missed, miss) << w;
      if (miss) {
        want.misses += 1;
      } else {
        want.hits += hits;
      }
    }
    expect_stats(cache.stats(), want, w);
    if (!failed || ::testing::Test::HasFailure()) break;
  }
  return failures;
}

} // namespace

TEST(PolicyFaults, FrozenStatisticsFollowThePublicCall) {
  // review-1d F1: a frozen cache's statistics are observable state (APPROACH 13; the building
  // cache's limit 2 does not cover them). Every public call that can read a supplied frozen cache
  // (the policy constructor in every mode, capacity in every mode, batch_capacity and sensitivity
  // at 1, 2 and 4 workers, in any thread and in the worker threads only) is swept: hits, misses,
  // entries, bytes and live values after a failed call equal those before it.
  const auto g = graph();
  const auto all = modes(*g);
  DagCache cache(*g);
  for (const auto& t : {std::vector<NodeId> {5, 6}, std::vector<NodeId> {6}, std::vector<NodeId> {5}}) {
    (void)cache.lookup(*g, t, Direction::Reverse, nullptr);
  }
  cache.freeze();
  const auto reads = [](const Mode& m) -> std::int64_t { return std::holds_alternative<CostOnly>(m.cfg.routing) ? 1 : 0; };
  for (const auto& m : all) {
    const long built = sweep_frozen_stats(m.name + " frozen construction", cache, false, 1, reads(m), false, [&] {
      FlowLedger fg(g);
      FlowPolicy p(fg, m.sources, m.destinations, m.cfg, nullptr, &cache);
    });
    EXPECT_GT(built, 2) << m.name;
    const CapacityOptions opts {m.cfg, nullptr, &cache, OutputOptions {true, true, true}};
    const long placed = sweep_frozen_stats(m.name + " frozen capacity", cache, false, 3, reads(m), false,
                                           [&] { (void)capacity(*g, m.sources, m.destinations, opts); });
    EXPECT_GT(placed, 2) << m.name;
    if (::testing::Test::HasFailure()) return;
  }
  // A CostOnly policy whose key is absent: the constructor's outcome is the miss.
  const auto& eb = all[1];  // costonly eb
  ASSERT_TRUE(std::holds_alternative<CostOnly>(eb.cfg.routing));
  const std::vector<NodeId> absent = {4};
  EXPECT_GT(sweep_frozen_stats("frozen construction miss", cache, false, 1, 0, true,
                               [&] {
                                 FlowLedger fg(g);
                                 FlowPolicy p(fg, eb.sources, absent, eb.cfg, nullptr, &cache);
                               }),
            1);
  const CapacityOptions shared {eb.cfg, nullptr, &cache};
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {eb.sources, eb.destinations}, {{1}, {6}}, {{0}, {5}}};
  // Pair 0 hits, pair 1 misses: a lower pair's hit is not kept by the failed batch.
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> with_miss = {{eb.sources, eb.destinations},
                                                                                     {{0}, absent}};
  for (const std::int32_t workers : {1, 2, 4}) {
    for (const bool in_workers : {false, true}) {
      if (in_workers && workers == 1) continue;  // one worker runs on the calling thread
      const std::string w = " workers " + std::to_string(workers);
      EXPECT_GT(sweep_frozen_stats("frozen batch" + w, cache, in_workers, 5, 3, false,
                                   [&] { (void)batch_capacity(*g, pairs, shared, {}, workers); }),
                2)
          << w;
      EXPECT_GT(sweep_frozen_stats("frozen batch with a miss" + w, cache, in_workers, 3, 0, true,
                                   [&] { (void)batch_capacity(*g, with_miss, shared, {}, workers); }),
                1)
          << w;
      if (::testing::Test::HasFailure()) return;
    }
    // CostOnly sensitivity from a source that reaches no destination: one lookup (the baseline's
    // hit), no placed deletion, and the call succeeds; its failures happen after the hit.
    EXPECT_GT(sweep_frozen_stats("frozen costonly sensitivity" + std::string(" workers ") + std::to_string(workers),
                                 cache, false, 3, 1, false,
                                 [&] {
                                   (void)sensitivity(*g, std::vector<NodeId> {5}, std::vector<NodeId> {6},
                                                     SensitivityOptions {eb.cfg, nullptr, &cache}, workers);
                                 }),
              2);
  }
  // CostOnly sensitivity whose deletions miss (a fresh uid per deletion): the baseline's hit is
  // dropped with the failed call, and the miss counts once.
  EXPECT_GT(sweep_frozen_stats("frozen costonly sensitivity miss", cache, false, 1, 0, true,
                               [&] {
                                 (void)sensitivity(*g, eb.sources, eb.destinations,
                                                   SensitivityOptions {eb.cfg, nullptr, &cache}, 1);
                               }),
            2);
  EXPECT_EQ(cache.stats().entries, 3);
}

TEST(PolicyFaults, WorkersSharingAFrozenCacheFailCleanly) {
  // batch_capacity and CostOnly sensitivity with a frozen cache shared by two and four workers:
  // a failure anywhere, or in the worker threads only, propagates as std::bad_alloc after the
  // join, the cache keeps its entries, and a clean retry gives the reference.
  const auto g = graph();
  const auto all = modes(*g);
  const auto& m = all[1];  // costonly eb: reads the cache
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {m.sources, m.destinations}, {{1}, {6}}, {{0}, {5}}};
  DagCache cache(*g);
  for (const auto& [s, t] : pairs) (void)cache.lookup(*g, t, Direction::Reverse, nullptr);
  const auto entries = cache.stats().entries;
  cache.freeze();
  const auto want = batch_capacity(*g, pairs, CapacityOptions {m.cfg}, {}, 1);
  const CapacityOptions shared {m.cfg, nullptr, &cache};
  for (const std::int32_t workers : {2, 4}) {
    for (const bool in_workers : {false, true}) {
      const long failures = sweep_evaluator(
          "frozen batch", workers, in_workers, 5,
          [&](std::int32_t w) { (void)batch_capacity(*g, pairs, shared, {}, w); },
          [&](std::int32_t w) {
            EXPECT_EQ(cache.stats().entries, entries);
            const auto b = batch_capacity(*g, pairs, shared, {}, w);
            ASSERT_EQ(b.size(), want.size());
            for (std::size_t i = 0; i < b.size(); ++i) expect_same(b[i], want[i], "frozen batch retry");
          });
      EXPECT_GT(failures, 2) << "frozen batch workers " << workers;
    }
  }
  // CostOnly sensitivity through a frozen cache: every re-evaluated deletion misses (a fresh
  // uid), so the reference outcome is the DagCacheMiss itself (its key; the uid is fresh per
  // call). With one worker an injected failure before that miss gives std::bad_alloc, and the
  // retry raises the same miss. (With several workers a lower candidate's miss may win over a
  // worker's std::bad_alloc: the lowest index's exception is raised.)
  const auto key_of = [](const DagCacheMiss& e) {
    const std::string w = e.what();
    return w.substr(0, w.find(", Eligibility uid "));
  };
  const auto costonly_sensitivity = [&](std::int32_t w) {
    (void)sensitivity(*g, m.sources, m.destinations, SensitivityOptions {m.cfg, nullptr, &cache}, w);
  };
  std::string want_miss;
  try {
    costonly_sensitivity(1);
  } catch (const DagCacheMiss& e) {
    want_miss = key_of(e);
  }
  EXPECT_EQ(want_miss, "DagCache: the frozen cache has no entry for roots {5, 6}, direction Reverse");
  const long failures = sweep_evaluator(
      "frozen sensitivity", 1, false, 1,
      [&](std::int32_t w) {
        try {
          costonly_sensitivity(w);
        } catch (const DagCacheMiss&) {
        }
      },
      [&](std::int32_t w) {
        EXPECT_EQ(cache.stats().entries, entries);
        try {
          costonly_sensitivity(w);
          ADD_FAILURE() << "no DagCacheMiss";
        } catch (const DagCacheMiss& e) {
          EXPECT_EQ(key_of(e), want_miss);
        }
      });
  EXPECT_GT(failures, 2) << "frozen sensitivity";
  // Global (deletion-mode) sensitivity with the frozen cache shared by every worker: it reads no
  // cache, and fails and retries as without one.
  const auto git = std::find_if(all.begin(), all.end(), [](const Mode& x) { return x.name == "global"; });
  ASSERT_NE(git, all.end());
  const auto& gm = *git;
  const auto want_global = sensitivity(*g, gm.sources, gm.destinations, SensitivityOptions {gm.cfg}, 1);
  const SensitivityOptions gopts {gm.cfg, nullptr, &cache};
  for (const std::int32_t workers : {2, 4}) {
    for (const bool in_workers : {false, true}) {
      const long global_failures = sweep_evaluator(
          "frozen global sensitivity", workers, in_workers, 11,
          [&](std::int32_t w) { (void)sensitivity(*g, gm.sources, gm.destinations, gopts, w); },
          [&](std::int32_t w) {
            const auto r = sensitivity(*g, gm.sources, gm.destinations, gopts, w);
            ASSERT_EQ(r.entries.size(), want_global.entries.size());
            for (std::size_t i = 0; i < r.entries.size(); ++i) {
              EXPECT_EQ(r.entries[i].edge, want_global.entries[i].edge);
              EXPECT_EQ(r.entries[i].loss, want_global.entries[i].loss);
            }
          });
      EXPECT_GT(global_failures, 2) << "frozen global sensitivity workers " << workers;
    }
  }
  EXPECT_EQ(cache.stats().entries, entries);
}
