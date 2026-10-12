// Allocation-failure atomicity of the policy engine, the DagCache and the evaluators.
// Exhaustive sweeps (each allocation fails once) cover the placement calls of the 13 routing and
// kernel modes below, with and without a static Eligibility; the zero and empty-roster calls; the
// policy construction, with and without a supplied building cache and an Eligibility; the
// DagCache construction and insertion; capacity() and batch_capacity() on one worker in every
// mode (batch_capacity: every mode but Pinned), without options, with an Eligibility and, in the
// CostOnly modes, with a supplied building cache; sensitivity() on one worker in one mode per
// routing model with the same options; and the frozen lookups. Sampled sweeps cover the other
// modes' sensitivity on one worker (every 11th allocation), the evaluators with 2 and 4 workers
// and the frozen-cache statistics: every 7th (batch_capacity; every 5th on a frozen cache), 11th
// (sensitivity) or step-th allocation in the process, and every third allocation of the worker
// threads only. The scratch-growth test sweeps every allocation of a placement on a 3,000-node
// graph after warming up on a 3-node graph, on a fresh thread for each attempt. After each injected
// std::bad_alloc the ledger (aggregate, residual, every record field, owners, next FlowId) and
// the policy's epoch state (requested, +inf, Global's cumulative curve, and per flow its id, route
// cost and plan: route, members, weights) equal their state before the call; a clean retry gives
// the reference result (no poisoned per-thread scratch); and remove() leaves the ledger as the
// reference's own remove() does. No failed call leaves an allocation live once its attempt's
// objects are gone (a supplied building cache with its kept entries included). remove(), a
// policy's destruction and DagCache::clear() allocate nothing. A supplied building cache is a
// memo: it keeps the complete entries of the lookups that succeeded before the failure (with
// their misses), never a partial one. The hook's worker-only mode, which the sweeps with several
// workers use, is checked against the caller's own allocations. Runs in the fault-injection
// executable (fault_injection.hpp).
#include <gtest/gtest.h>

#include <algorithm>
#include <exception>
#include <functional>
#include <memory>
#include <new>
#include <optional>
#include <span>
#include <string>
#include <thread>
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
  const auto det = PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity);
  const auto res = PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Residual);
  const auto tunnel = PathSelection(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, false, false);
  const std::vector<NodeId> S = {0, 1};
  const std::vector<NodeId> D = {5, 6};
  out.push_back({"costonly maxflow", {CostOnly {}, FlowPlacement::MaxFlow, Admission::Lossless, det}, S, D});
  out.push_back({"costonly eb", cost_only_config(), S, D});
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
  out.push_back({"global", global_config(), S, D});
  Pinned two;
  const auto inv = id_of_input(g);
  two.bundles = {dag_of(g, {inv[0], inv[4], inv[5]}), dag_of(g, {inv[1], inv[6]})};
  out.push_back({"pinned eb cohort", {two, FlowPlacement::EqualBalanced, Admission::Lossless, det}, {0}, {5}});
  out.push_back({"pinned lossy", {two, FlowPlacement::Weighted, Admission::Lossy, det}, {0}, {5}});
  // CostOnly creation in one pass: its Collapse branch (multipath, multi_edge false,
  // on the creation view) and its Weighted weights (Capacity priority).
  out.push_back({"costonly eb multipath collapse",
                 {CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless,
                  PathSelection(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, true, false)},
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

// The static Eligibility of g without the edge at input position `position`.
Eligibility without_input_edge(const Graph& g, std::size_t position) {
  const EdgeId excluded[] = {id_of_input(g)[position]};
  return g.eligibility({}, excluded);
}

// graph() without the costly bypass 0 -> 6 (input 9), and without the cross edge 2 -> 4 (input
// 10: batch_capacity's per-pair value, whose intersection with the first is a new value).
Eligibility without_bypass(const Graph& g) { return without_input_edge(g, 9); }
Eligibility without_cross(const Graph& g) { return without_input_edge(g, 10); }

// A third party holding one unit on the first route's edges.
std::unique_ptr<FlowPolicy> third_party(FlowLedger& fg, const Graph& g) {
  const auto inv = id_of_input(g);
  Pinned p;
  p.bundles = {dag_of(g, {inv[0], inv[4]})};
  auto t = std::make_unique<FlowPolicy>(fg, std::vector<NodeId> {0}, std::vector<NodeId> {5},
                                        FlowPolicyConfig {p, FlowPlacement::MaxFlow, Admission::Lossless,
                                                          PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity)});
  (void)t->place(1);
  return t;
}

} // namespace

TEST(PolicyFaults, EveryAllocationOfEveryRoutingModeFailsAtomically) {
  const auto g = graph();
  const Eligibility no_bypass = without_bypass(*g);
  for (const auto& m : modes(*g)) {
    for (const auto& sc : scenarios()) {
      for (const Eligibility* elig : {static_cast<const Eligibility*>(nullptr), &no_bypass}) {
        const std::string where = m.name + " / " + sc.name + (elig != nullptr ? " / eligibility" : "");
        // The third party, then the policy's earlier calls.
        const auto setup = [&](FlowLedger& fg) {
          auto third = std::shared_ptr<FlowPolicy>(third_party(fg, *g));
          auto p = std::make_shared<FlowPolicy>(fg, m.sources, m.destinations, m.cfg, elig);
          for (auto v : sc.before) (void)p->place(v);
          FlowPolicy* policy = p.get();
          return Staged {std::make_shared<std::pair<std::shared_ptr<FlowPolicy>, std::shared_ptr<FlowPolicy>>>(
                             std::move(third), std::move(p)),
                         policy};
        };
        const long injected =
            sweep_placement(g, where, setup, [&](FlowPolicy& p) { return p.place(sc.call, sc.outputs); });
        if (::testing::Test::HasFailure()) return;
        EXPECT_GT(injected, 2) << where;  // not vacuous: every mode and scenario allocates
      }
    }
  }
}

// The calls that place nothing publish only their outputs and, for a positive request, the
// epoch's request: a zero request, and a positive or +inf request when no source reaches a
// destination (an empty roster; the endpoints swapped: 5 and 6 have no out-edges). Each runs
// through the shared placement sweep: a failure changes neither the ledger nor the epoch state and
// leaves no allocation live; the completing call, the retry and remove() leave the ledger as the
// reference's do.
TEST(PolicyFaults, ZeroAndEmptyRosterCallsFailAtomically) {
  const auto g = graph();
  const OutputOptions all {true, true, true};
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
      const std::string where = m.name + " / " + c.name;
      const auto setup = [&](FlowLedger& fg) {
        auto third = std::shared_ptr<FlowPolicy>(third_party(fg, *g));
        auto p = std::make_shared<FlowPolicy>(fg, s, d, m.cfg);
        if (c.swap) EXPECT_TRUE(detail::PolicyAccess::engine(*p).roster.empty()) << where;
        FlowPolicy* policy = p.get();
        return Staged {std::make_shared<std::pair<std::shared_ptr<FlowPolicy>, std::shared_ptr<FlowPolicy>>>(
                           std::move(third), std::move(p)),
                       policy};
      };
      const long injected = sweep_placement(g, where, setup, [&](FlowPolicy& p) { return p.place(c.volume, all); });
      if (::testing::Test::HasFailure()) return;
      EXPECT_GT(injected, 0) << where;  // not vacuous: every mode and case fails
    }
  }
}

TEST(PolicyFaults, ConstructionFailsWithoutBindingAndTheCacheStaysConsistent) {
  const auto g = graph();
  const Eligibility no_bypass = without_bypass(*g);
  for (const auto& m : modes(*g)) {
    for (const int v : {0, 1, 2, 3}) {
      const bool with_cache = (v & 1) != 0;
      const Eligibility* elig = (v & 2) != 0 ? &no_bypass : nullptr;
      const std::string variant = std::string(with_cache ? " cached" : "") + (elig != nullptr ? " eligibility" : "");
      // The reference: a clean construction on a fresh ledger and cache, then place(2).
      FlowLedger rf(g);
      DagCache rc(*g);
      FlowPolicy pr(rf, m.sources, m.destinations, m.cfg, elig, with_cache ? &rc : nullptr);
      const DagCache::Stats fresh = rc.stats();
      const PlacementResult want = pr.place(2);
      const LedgerState placed = ledger_state(rf);
      long injected = 0;
      int empty = 0;  // failures that left the supplied cache empty, and that left its whole entry
      int whole = 0;
      for (long at = 0;; ++at) {
        const std::string w = m.name + variant + " construction allocation " + std::to_string(at);
        // Everything an attempt allocates (its ledger, cache, failed construction, then the clean
        // construction and its placement) is freed when the attempt ends: the thread's scratch is
        // already grown by the reference, so an allocation left live is a leak of the failed
        // construction (a cache's kept entry is freed with the cache).
        const long live = g_live.load();
        {
          FlowLedger fg(g);
          DagCache cache(*g);
          const LedgerState before = ledger_state(fg);
          std::optional<FlowPolicy> p;
          if (!fails_at(at, [&] { p.emplace(fg, m.sources, m.destinations, m.cfg, elig, with_cache ? &cache : nullptr); })) {
            EXPECT_EQ(LedgerAccess::bound_owners(fg).size(), 1U);
            break;
          }
          ++injected;
          EXPECT_FALSE(p.has_value()) << w;
          expect_same(ledger_state(fg), before, w);
          // A supplied cache holds whole entries only: none, or exactly the entries and bytes a
          // clean construction leaves in a fresh cache; no policy holds a value any more.
          const auto s = cache.stats();
          EXPECT_LE(fresh.entries, 1) << w;
          if (s.entries == 0) {
            EXPECT_EQ(s.bytes, 0) << w;
            ++empty;
          } else {
            EXPECT_EQ(s.entries, fresh.entries) << w;
            EXPECT_EQ(s.bytes, fresh.bytes) << w;
            ++whole;
          }
          EXPECT_EQ(s.live_values, s.entries) << w;
          // A clean construction works, leaves the fresh cache's entries and places like the reference.
          p.emplace(fg, m.sources, m.destinations, m.cfg, elig, with_cache ? &cache : nullptr);
          EXPECT_EQ(cache.stats().entries, fresh.entries) << w;
          EXPECT_EQ(cache.stats().bytes, fresh.bytes) << w;
          expect_same(p->place(2), want, w);
          expect_same(ledger_state(fg), placed, w);
          expect_balanced(fg, w);
        }
        EXPECT_EQ(g_live.load(), live) << w << ": a leak";
        if (::testing::Test::HasFailure()) return;
      }
      EXPECT_GT(injected, 0) << m.name << variant;  // not vacuous
      // A mode that reads the supplied cache fails both before its entry exists and after.
      if (fresh.entries == 1) {
        EXPECT_GT(empty, 0) << m.name;
        EXPECT_GT(whole, 0) << m.name;
      }
    }
  }
  // Every CostOnly mode reads the supplied cache (an entry after a clean construction), the
  // other modes none: the memo outcomes above are checked for every CostOnly mode.
  for (const auto& m : modes(*g)) {
    FlowLedger fg(g);
    DagCache cache(*g);
    const FlowPolicy p(fg, m.sources, m.destinations, m.cfg, nullptr, &cache);
    EXPECT_EQ(cache.stats().entries, std::holds_alternative<CostOnly>(m.cfg.routing) ? 1 : 0) << m.name;
  }
}

TEST(PolicyFaults, DagCacheInsertionIsAllOrNothing) {
  const auto g = graph();
  const std::vector<NodeId> roots = {5, 6};
  {  // clean searches first grow this thread's search scratch, which outlives the attempts
    DagCache warm(*g);
    (void)warm.lookup(*g, std::vector<NodeId> {5}, Direction::Reverse, nullptr);
    (void)warm.lookup(*g, roots, Direction::Reverse, nullptr);
  }
  long injected = 0;
  for (long at = 0;; ++at) {
    const std::string w = "lookup allocation " + std::to_string(at);
    // Everything an attempt allocates (its cache and values, the retry) is freed when it ends, so
    // an allocation left live is a leak of the failed lookup.
    const long live = g_live.load();
    bool failed = false;
    {
      DagCache cache(*g);
      (void)cache.lookup(*g, std::vector<NodeId> {5}, Direction::Reverse, nullptr);  // one entry before
      const auto before = cache.stats();
      std::shared_ptr<const PredDAG> got;
      failed = fails_at(at, [&] { got = cache.lookup(*g, roots, Direction::Reverse, nullptr); });
      if (!failed) {
        EXPECT_EQ(cache.stats().entries, 2);
      } else {
        ++injected;
        const auto s = cache.stats();
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
    EXPECT_EQ(g_live.load(), live) << w << ": a leak";
    if (!failed) break;
  }
  EXPECT_GE(injected, 2);  // the key (allocation 0) and at least one allocation after it
}

// A failed DagCache construction (its allocations: the counter of live values that every value's
// deleter shares, and on MSVC the map's own) builds nothing and leaks nothing; clear() and
// freeze()/thaw() make no allocation.
TEST(PolicyFaults, DagCacheConstructionFailsWithoutALeakAndClearNeverAllocates) {
  const auto g = graph();
  long injected = 0;
  for (long at = 0;; ++at) {
    std::optional<DagCache> cache;
    const long live = g_live.load();
    if (!fails_at(at, [&] { cache.emplace(*g); })) {
      ASSERT_TRUE(cache.has_value());
      EXPECT_EQ(cache->stats().entries, 0);
      break;
    }
    ++injected;
    EXPECT_FALSE(cache.has_value()) << "allocation " << at;
    EXPECT_EQ(g_live.load(), live) << "allocation " << at << ": a leak";
  }
  EXPECT_GE(injected, 1);  // the counter (MSVC's map adds its own)
  DagCache cache(*g);
  const auto held = cache.lookup(*g, std::vector<NodeId> {5}, Direction::Reverse, nullptr);
  (void)cache.lookup(*g, std::vector<NodeId> {5, 6}, Direction::Reverse, nullptr);
  ASSERT_EQ(cache.stats().entries, 2);
  EXPECT_FALSE(fails_at(0, [&] {
    cache.freeze();
    cache.thaw();
    cache.clear();
  }));
  const auto s = cache.stats();
  EXPECT_EQ(s.entries, 0);
  EXPECT_EQ(s.bytes, 0);
  EXPECT_EQ(s.live_values, 1);  // the value held here outlives its entry
}

namespace {

// The evaluators under a failure on `workers` threads: anywhere in the process (every `step`-th
// allocation: 1 sweeps them all) or, with several workers, in the worker threads only (their
// setup and their tasks; every third worker allocation, from the first). `setup` runs before each
// attempt and `teardown` after it, both outside the failure; teardown frees what setup made (a
// supplied building cache, whose kept entries are no leak). The failure propagates as
// std::bad_alloc after every worker has joined, never terminating the process, and a clean retry
// gives the reference. Each attempt leaves no allocation live: the workers' scratch ends with
// their threads, and the caller computed the reference on its own thread first, which grew that
// thread's scratch. Returns the number of injected failures.
long sweep_evaluator(const std::string& name, std::int32_t workers, bool in_workers, long step,
                     const std::function<void(std::int32_t)>& call, const std::function<void(std::int32_t)>& retry,
                     const std::function<void()>& setup = [] {}, const std::function<void()>& teardown = [] {}) {
  long failures = 0;
  for (long at = 0;; at += in_workers ? 3 : step) {
    const std::string w = name + " workers " + std::to_string(workers) + (in_workers ? " (worker threads)" : "") +
                          " allocation " + std::to_string(at);
    const long live = g_live.load();
    setup();
    const auto fn = [&] { call(workers); };
    const bool failed = in_workers ? fails_in_workers_at(at, fn) : fails_at(at, fn);
    if (failed) {
      ++failures;
      retry(workers);
    }
    teardown();
    EXPECT_EQ(g_live.load(), live) << w << ": a leak";
    if (!failed) break;
    if (::testing::Test::HasFailure()) {
      ADD_FAILURE() << w;
      break;
    }
  }
  return failures;
}

// The options of an evaluator sweep: none, a static Eligibility (for batch_capacity also a
// per-pair value, so each pair's intersection is a new value), or a supplied building cache (one
// worker; a fresh cache before each attempt). A building cache is a memo: after a failed call it
// holds whole entries only (no value alive outside it), some failed calls keep the entries of the
// lookups that completed before the failure, and the retry leaves the entries of a clean call.
struct Options {
  std::string name;
  const Eligibility* elig {nullptr};
  const Eligibility* per_pair {nullptr};
  bool cache {false};
};

std::vector<Options> option_sets(const Eligibility& no_bypass, const Eligibility& no_cross) {
  return {{"", nullptr, nullptr, false}, {" eligibility", &no_bypass, &no_cross, false}, {" building cache", nullptr, nullptr, true}};
}

// A supplied building cache after a failed call: whole entries only, no more than a clean call
// leaves.
void expect_whole_entries(const DagCache::Stats& s, std::int64_t clean, const std::string& where) {
  EXPECT_LE(s.entries, clean) << where;
  EXPECT_EQ(s.live_values, s.entries) << where;
  if (s.entries == 0) EXPECT_EQ(s.bytes, 0) << where;
}

void expect_whole_entries(const std::optional<DagCache>& cache, std::int64_t clean, const std::string& where) {
  if (cache) expect_whole_entries(cache->stats(), clean, where);
}

void sweep_batch(const Mode& m, const Graph& g, std::span<const std::int32_t> workers_set, long step,
                 const Options& o = {}) {
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {m.sources, m.destinations}, {{1}, {6}}, {{0}, {5}}};
  const std::vector<const Eligibility*> per_pair(o.per_pair != nullptr ? pairs.size() : 0, o.per_pair);
  std::optional<DagCache> cache;
  const auto opts = [&] { return CapacityOptions {m.cfg, o.elig, cache ? &*cache : nullptr}; };
  if (o.cache) cache.emplace(g);
  const auto want = batch_capacity(g, pairs, opts(), per_pair, 1);
  const std::int64_t clean = cache ? cache->stats().entries : 0;
  cache.reset();  // each attempt gets a fresh one (setup)
  const std::string name = m.name + " batch" + o.name;
  long kept = 0;  // failed calls whose cache kept the entries of the lookups before the failure
  for (const std::int32_t workers : workers_set) {
    for (const bool in_workers : {false, true}) {
      if (in_workers && workers == 1) continue;  // one worker runs on the calling thread
      const long failures = sweep_evaluator(
          name, workers, in_workers, step, [&](std::int32_t w) { (void)batch_capacity(g, pairs, opts(), per_pair, w); },
          [&](std::int32_t w) {
            expect_whole_entries(cache, clean, name + " failed call");
            if (cache && cache->stats().entries > 0) ++kept;
            const auto b = batch_capacity(g, pairs, opts(), per_pair, w);
            ASSERT_EQ(b.size(), want.size());
            for (std::size_t i = 0; i < b.size(); ++i) expect_same(b[i], want[i], name + " retry");
            if (cache) EXPECT_EQ(cache->stats().entries, clean) << name << " retry";
          },
          [&] {
            if (o.cache) cache.emplace(g);
          },
          [&] { cache.reset(); });
      EXPECT_GT(failures, 2) << name << " workers " << workers;  // not vacuous
    }
  }
  if (o.cache) EXPECT_GT(kept, 0) << name;
}

void sweep_sensitivity(const Mode& m, const Graph& g, std::span<const std::int32_t> workers_set, long step,
                       const Options& o = {}) {
  std::optional<DagCache> cache;
  const auto opts = [&] { return SensitivityOptions {m.cfg, o.elig, cache ? &*cache : nullptr}; };
  if (o.cache) cache.emplace(g);
  const auto want = sensitivity(g, m.sources, m.destinations, opts(), 1);
  const std::int64_t clean = cache ? cache->stats().entries : 0;
  cache.reset();  // each attempt gets a fresh one (setup)
  const std::string name = m.name + " sensitivity" + o.name;
  long kept = 0;  // failed calls whose cache kept the entries of the lookups before the failure
  for (const std::int32_t workers : workers_set) {
    for (const bool in_workers : {false, true}) {
      if (in_workers && workers == 1) continue;
      const long failures = sweep_evaluator(
          name, workers, in_workers, step,
          [&](std::int32_t w) { (void)sensitivity(g, m.sources, m.destinations, opts(), w); },
          [&](std::int32_t w) {
            expect_whole_entries(cache, clean, name + " failed call");
            if (cache && cache->stats().entries > 0) ++kept;
            const auto s = sensitivity(g, m.sources, m.destinations, opts(), w);
            expect_same(s.baseline, want.baseline, name + " retry");
            ASSERT_EQ(s.entries.size(), want.entries.size());
            for (std::size_t i = 0; i < s.entries.size(); ++i) {
              EXPECT_EQ(s.entries[i].edge, want.entries[i].edge) << name;
              EXPECT_EQ(s.entries[i].loss, want.entries[i].loss) << name;
            }
          },
          [&] {
            if (o.cache) cache.emplace(g);
          },
          [&] { cache.reset(); });
      EXPECT_GT(failures, 2) << name << " workers " << workers;
    }
  }
  if (o.cache) EXPECT_GT(kept, 0) << name;
}

} // namespace

// The hook's worker-only mode, which every sweep of the evaluators with several workers relies on:
// it counts and fails only the allocations of threads other than the caller's, so an allocation
// of the calling thread neither fails nor advances the countdown (were the thread filter lost,
// those sweeps would quietly fail any thread's allocations). reaches_in_workers_at accepts an
// injected failure that reaches the caller or loses to another exception that fn handles, and
// reports a std::bad_alloc that it did not inject.
TEST(FaultHook, WorkerOnlyModeFailsOnlyWorkerAllocations) {
  const auto caller_allocates = [] { ::operator delete(::operator new(8)); };
  // Three allocations of the calling thread, then one of a worker, whose exception the caller
  // rethrows after the join (as the evaluators do).
  const auto worker_allocates = [](bool handle) {
    for (int i = 0; i < 3; ++i) ::operator delete(::operator new(8));
    std::exception_ptr error;
    std::thread([&] {
      try {
        ::operator delete(::operator new(8));
      } catch (...) {
        error = std::current_exception();
      }
    }).join();
    if (error && !handle) std::rethrow_exception(error);
  };
  EXPECT_EQ(failures_of([&] { EXPECT_FALSE(fails_in_workers_at(0, caller_allocates)); }), 0);
  EXPECT_EQ(failures_of([&] { EXPECT_TRUE(fails_in_workers_at(0, [&] { worker_allocates(false); })); }), 0);
  EXPECT_EQ(failures_of([&] { EXPECT_FALSE(fails_in_workers_at(1, [&] { worker_allocates(false); })); }), 0);
  bool bad_alloc = false;
  EXPECT_EQ(failures_of([&] { EXPECT_TRUE(reaches_in_workers_at(0, [&] { worker_allocates(false); }, bad_alloc)); }), 0);
  EXPECT_TRUE(bad_alloc);
  EXPECT_EQ(failures_of([&] { EXPECT_TRUE(reaches_in_workers_at(0, [&] { worker_allocates(true); }, bad_alloc)); }), 0);
  EXPECT_FALSE(bad_alloc);  // lost to an outcome fn handled
  EXPECT_EQ(failures_of([&] { EXPECT_FALSE(reaches_in_workers_at(0, caller_allocates, bad_alloc)); }), 0);
  // A std::bad_alloc of fn's own while the armed worker failure is never reached.
  EXPECT_EQ(failures_of([&] { (void)reaches_in_workers_at(1000, [] { throw std::bad_alloc(); }, bad_alloc); }), 1);
}

// Every allocation of capacity and batch_capacity on one worker fails once, for every mode, without
// options, with a static Eligibility and (for CostOnly, the one routing model that reads a cache)
// with a supplied building cache. sensitivity repeats capacity's work once per candidate edge, so
// its exhaustive sweep runs on one mode per routing model with the same options (every allocation
// site of its own code: the deletion values, the intersections, the CostOnly skip, the per-
// candidate cache clear); the other modes fail at every 11th allocation.
TEST(PolicyFaults, EvaluatorsLeaveNoPoisonedScratch) {
  const auto g = graph();
  const Eligibility no_bypass = without_bypass(*g);
  const Eligibility no_cross = without_cross(*g);
  const std::int32_t one[] = {1};
  const std::vector<std::string> exhaustive_sensitivity = {"costonly eb", "capacityaware weighted diminishing", "global",
                                                           "pinned eb cohort"};
  int swept_sensitivity = 0;
  for (const auto& m : modes(*g)) {
    const bool reads_cache = std::holds_alternative<CostOnly>(m.cfg.routing);
    const bool every_sensitivity_allocation =
        std::find(exhaustive_sensitivity.begin(), exhaustive_sensitivity.end(), m.name) != exhaustive_sensitivity.end();
    for (const auto& o : option_sets(no_bypass, no_cross)) {
      if (o.cache && !reads_cache) continue;
      const std::string name = m.name + " capacity" + o.name;
      std::optional<DagCache> cache;
      const auto copts = [&] {
        return CapacityOptions {m.cfg, o.elig, cache ? &*cache : nullptr, OutputOptions {true, true, true}};
      };
      if (o.cache) cache.emplace(*g);
      const auto want = capacity(*g, m.sources, m.destinations, copts());
      const std::int64_t clean = cache ? cache->stats().entries : 0;
      if (o.cache) EXPECT_GT(clean, 0) << name;  // the cache is read
      cache.reset();  // each attempt supplies a fresh one
      long injected = 0;
      long kept = 0;  // failed calls whose cache kept the entry of the lookup before the failure
      long empty = 0;
      for (long at = 0;; ++at) {
        const std::string w = name + " allocation " + std::to_string(at);
        // Everything an attempt allocates is freed when it ends, the supplied cache with it (its
        // kept entries are no leak), so an allocation left live is a leak of the failed call.
        const long live = g_live.load();
        if (o.cache) cache.emplace(*g);
        bool failed = false;
        {
          PlacementResult got;
          failed = fails_at(at, [&] { got = capacity(*g, m.sources, m.destinations, copts()); });
          if (!failed) {
            expect_same(got, want, name + " (no failure)");
          } else {
            ++injected;
            expect_whole_entries(cache, clean, w);
            if (cache) {
              // The memo keeps a completed lookup's entry and its statistics: one miss per entry.
              const auto s = cache->stats();
              ++(s.entries > 0 ? kept : empty);
              if (s.entries > 0) {
                EXPECT_EQ(s.misses, s.entries) << w;
                EXPECT_EQ(s.hits, 0) << w;
              }
            }
            expect_same(capacity(*g, m.sources, m.destinations, copts()), want, w + " (retry)");
            if (cache) EXPECT_EQ(cache->stats().entries, clean) << w;
          }
        }
        cache.reset();
        EXPECT_EQ(g_live.load(), live) << w << ": a leak";
        if (!failed) break;
      }
      EXPECT_GT(injected, 0) << name;  // not vacuous
      if (o.cache) {
        EXPECT_GT(kept, 0) << name;  // failures after the lookup, and before it
        EXPECT_GT(empty, 0) << name;
      }
      if (!std::holds_alternative<Pinned>(m.cfg.routing)) sweep_batch(m, *g, one, 1, o);  // batch rejects Pinned
      if (every_sensitivity_allocation) {
        sweep_sensitivity(m, *g, one, 1, o);
        ++swept_sensitivity;
      } else if (o.elig == nullptr && !o.cache) {
        sweep_sensitivity(m, *g, one, 11, o);
      }
      if (::testing::Test::HasFailure()) return;
    }
  }
  EXPECT_EQ(swept_sensitivity, 4 * 2 + 1);  // four modes with two option sets, CostOnly's with three
}

TEST(PolicyFaults, BatchWorkersFailWithoutTerminating) {
  const auto g = graph();
  const std::int32_t several[] = {2, 4};
  for (const auto& m : modes(*g)) {
    if (std::holds_alternative<Pinned>(m.cfg.routing)) continue;
    sweep_batch(m, *g, several, 7);
    if (::testing::Test::HasFailure()) return;
  }
}

TEST(PolicyFaults, SensitivityWorkersFailWithoutTerminating) {
  const auto g = graph();
  const std::int32_t several[] = {2, 4};
  for (const auto& m : modes(*g)) {
    sweep_sensitivity(m, *g, several, 11);
    if (::testing::Test::HasFailure()) return;
  }
}

TEST(PolicyFaults, RemovalAndReleaseNeverAllocate) {
  const auto g = graph();
  for (const auto& m : modes(*g)) {
    FlowLedger fg(g);
    auto p = std::make_unique<FlowPolicy>(fg, m.sources, m.destinations, m.cfg);
    expect_removal_and_release_allocate_nothing(fg, p, m.name);
  }
}

// The working-residual overlay and the per-source scratch grow with the graph. Each attempt runs
// on a new thread (its scratch empty) that first places on a 3-node graph, then the swept call on
// a larger one, so every allocation of that call fails once while its scratch grows from the
// small size; the failed call changes nothing, and a clean retry on the same graph (or the call
// itself, once no failure is reached) gives the reference call's result, every output included,
// ledger and epoch (ASan reports any access outside a partly grown array). The same call on a new
// thread whose scratch has first grown to the large graph allocates less: the sweep reached the
// growth.
TEST(PolicyFaults, ScratchGrowthFailureLeavesNothingHalfGrown) {
  const FlowPolicyConfig cfg = cost_only_config();
  const OutputOptions all {true, true, true};
  const auto small = chain_graph(3, 5);
  const auto large = chain_graph(3000, 5);
  const std::int32_t n = large->num_nodes();
  // The reference: the call on a fresh ledger without a failure.
  PlacementResult want;
  LedgerState want_state;
  Epoch want_epoch;
  {
    FlowLedger fg(large);
    FlowPolicy p(fg, std::vector<NodeId> {0}, std::vector<NodeId> {n - 1}, cfg);
    want = p.place(2, all);
    want_state = ledger_state(fg);
    want_epoch = epoch(p);
  }
  ASSERT_EQ(want.placed, 2);  // the whole chain carries the request: the comparison is not of zeros
  sweep_scratch_growth(
      "placement",
      [&] {
        FlowLedger fg(small);
        FlowPolicy p(fg, std::vector<NodeId> {0}, std::vector<NodeId> {2}, cfg);
        EXPECT_EQ(p.place(2, all).placed, 2);
      },
      [&](long at, const std::string& w) {
        FlowLedger fg(large);
        FlowPolicy p(fg, std::vector<NodeId> {0}, std::vector<NodeId> {n - 1}, cfg);
        const LedgerState before = ledger_state(fg);
        const Epoch before_epoch = epoch(p);
        PlacementResult got;
        const bool failed = fails_at(at, [&] { got = p.place(2, all); });
        if (failed) {
          expect_same(ledger_state(fg), before, w);
          expect_same(epoch(p), before_epoch, w);
          got = p.place(2, all);
        }
        expect_same(got, want, w);
        expect_same(ledger_state(fg), want_state, w);
        expect_same(epoch(p), want_epoch, w);
        expect_balanced(fg, w);
        return failed;
      },
      [&] {
        FlowLedger fg(large);
        FlowPolicy p(fg, std::vector<NodeId> {0}, std::vector<NodeId> {n - 1}, cfg);
        (void)p.place(2, all);  // the scratch grows to the large graph
        p.remove();
        long warm = 0;
        while (fails_at(warm, [&] {
          (void)p.place(2, all);
          p.remove();
        })) {
          ++warm;
        }
        return warm;
      });
}

TEST(PolicyFaults, FrozenLookupsFailWithoutCounting) {
  // A frozen lookup allocates a key only for roots that are not already a sorted set (a sorted set
  // is borrowed: nothing), and a frozen miss builds the DagCacheMiss message, all before counting:
  // a failure there counts nothing and changes no entry; a clean retry counts.
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
      const std::string w = std::string(hit ? "hit" : "miss") + " allocation " + std::to_string(at);
      const long live = g_live.load();
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
      EXPECT_EQ(g_live.load(), live) << w << ": a leak";  // the key and the miss's message are gone
      const auto s = cache.stats();
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
      EXPECT_GE(failures, hit ? 1 : 2) << (hit ? "hit" : "miss");  // the key; the message (sorted roots)
    }
  }
}

namespace {

void expect_stats(const DagCache::Stats& got, const DagCache::Stats& want, const std::string& where) {
  EXPECT_EQ(stats_of(got), stats_of(want)) << where << " (hits, misses, entries, bytes, live values)";
}

// One public call reading a supplied frozen cache, under every injected failure (every `step`-th
// allocation in any thread, or every third in the worker threads only): a call failing with
// std::bad_alloc leaves all five statistics as they were and no allocation live; the call that
// completes adds exactly `hits` hits, or, when its outcome is the DagCacheMiss (`miss`), exactly
// one miss and no hit. A clean call before the sweep grows the calling thread's scratch, which
// outlives the attempts (the workers' scratch ends with their threads). Returns the number of
// injected failures.
long sweep_frozen_stats(const std::string& name, const DagCache& cache, bool in_workers, long step, std::int64_t hits,
                        bool miss, const std::function<void()>& call) {
  try {
    call();
  } catch (const DagCacheMiss&) {
  }
  long failures = 0;
  for (long at = 0;; at += in_workers ? 3 : step) {
    const std::string w = name + " allocation " + std::to_string(at) + (in_workers ? " (worker threads)" : "");
    const long live = g_live.load();
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
    EXPECT_EQ(g_live.load(), live) << w << ": a leak";
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

// The cache comparisons of the sweeps report each defect they exist for: expect_stats any one of
// the five statistics changed; expect_whole_entries more entries than a clean call leaves, a
// value alive outside the cache, and payload bytes without an entry.
TEST(PolicyFaults, TheCacheComparatorsReportEveryDefect) {
  const DagCache::Stats s {3, 2, 1, 40, 1};
  EXPECT_EQ(failures_of([&] { expect_stats(s, s, "same"); }), 0);
  const std::vector<std::function<void(DagCache::Stats&)>> fields = {
      [](auto& x) { ++x.hits; },    [](auto& x) { ++x.misses; },      [](auto& x) { ++x.entries; },
      [](auto& x) { ++x.bytes; },   [](auto& x) { ++x.live_values; },
  };
  for (std::size_t k = 0; k < fields.size(); ++k) {
    DagCache::Stats x = s;
    fields[k](x);
    EXPECT_GT(failures_of([&] { expect_stats(x, s, "stats"); }), 0) << "statistic " << k;
  }
  EXPECT_EQ(failures_of([&] { expect_whole_entries(s, 1, "whole"); }), 0);
  EXPECT_EQ(failures_of([&] { expect_whole_entries(DagCache::Stats {0, 0, 0, 0, 0}, 1, "empty"); }), 0);
  EXPECT_GT(failures_of([&] { expect_whole_entries(s, 0, "more entries than a clean call"); }), 0);
  EXPECT_GT(failures_of([&] { expect_whole_entries(DagCache::Stats {3, 2, 1, 40, 2}, 1, "a value outside"); }), 0);
  EXPECT_GT(failures_of([&] { expect_whole_entries(DagCache::Stats {3, 2, 0, 40, 0}, 1, "bytes, no entry"); }), 0);
}

TEST(PolicyFaults, FrozenStatisticsFollowThePublicCall) {
  // A frozen cache's statistics are observable state (the memo allowance of a building cache
  // does not cover them). Every public call that can read a supplied frozen cache is swept: the
  // policy constructor and capacity in every mode; batch_capacity at 1, 2 and 4 workers, in any
  // thread and in the worker threads only, with and without a miss; CostOnly sensitivity without
  // a re-evaluated deletion at 1, 2 and 4 workers in any thread; and CostOnly sensitivity whose
  // deletions miss at 1 worker in any thread and at 2 and 4 workers in the worker threads only.
  // Hits, misses, entries, bytes and live values after a failed call equal those before it.
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
  const auto costonly_miss = [&](std::int32_t workers) {
    (void)sensitivity(*g, eb.sources, eb.destinations, SensitivityOptions {eb.cfg, nullptr, &cache}, workers);
  };
  EXPECT_GT(sweep_frozen_stats("frozen costonly sensitivity miss", cache, false, 1, 0, true, [&] { costonly_miss(1); }),
            2);
  // The same with the deletions evaluated by 2 and 4 workers, every third allocation of the
  // worker threads failing. The tasks start in increasing index and the lowest failing index's
  // exception is raised, so an injected failure either is raised (nothing counted) or loses to a
  // lower deletion's miss (the miss counted once, no hit).
  for (const std::int32_t workers : {2, 4}) {
    const std::string w = "frozen costonly sensitivity miss, worker threads, workers " + std::to_string(workers);
    long reached = 0, surfaced = 0;
    for (long at = 0;; at += 3) {
      const std::string where = w + " allocation " + std::to_string(at);
      const long live = g_live.load();
      const auto before = cache.stats();
      bool missed = false;
      bool bad_alloc = false;
      const bool hit = reaches_in_workers_at(
          at,
          [&] {
            try {
              costonly_miss(workers);
            } catch (const DagCacheMiss&) {
              missed = true;
            }
          },
          bad_alloc);
      EXPECT_EQ(g_live.load(), live) << where << ": a leak";
      EXPECT_NE(missed, bad_alloc) << where;  // one outcome or the other
      surfaced += hit && bad_alloc ? 1 : 0;
      auto want = before;
      if (missed) want.misses += 1;
      expect_stats(cache.stats(), want, where);
      if (!hit || ::testing::Test::HasFailure()) break;
      ++reached;
    }
    EXPECT_GT(reached, 2) << w;
    EXPECT_GT(surfaced, 0) << w;
  }
  EXPECT_EQ(cache.stats().entries, 3);
}

TEST(PolicyFaults, WorkersSharingAFrozenCacheFailCleanly) {
  // batch_capacity with a frozen cache shared by two and four workers, a failure anywhere or in
  // the worker threads only, then CostOnly sensitivity with that cache on one worker (its
  // statistics with two and four workers: the test above), and Global sensitivity, which reads no
  // cache, on two and four workers: the failure propagates as std::bad_alloc after the join, the
  // cache keeps its entries, and a clean retry gives the reference.
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
  // Global sensitivity (each deletion repaired from the baseline's state) with the frozen cache
  // shared by every worker: it reads no cache, and fails and retries as without one.
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
