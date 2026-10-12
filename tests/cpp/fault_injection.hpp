// The allocation-failure hook of the fault-injection executable (netgraph_core_fault_tests), with
// the state its sweeps compare (state_comparators.hpp): fault_injection.cpp replaces the global
// operator new so that the chosen allocation throws std::bad_alloc (the nothrow form never fails
// and does not advance the countdown; its allocations count as live like the others); every test
// file of the executable sweeps its operations through fails_at(). The counter is shared by every thread (the
// evaluators' workers allocate too): the k-th allocation in the process fails, whichever thread
// makes it.
#pragma once

#include <gtest/gtest.h>

#include <atomic>
#include <functional>
#include <memory>
#include <string>
#include <thread>

#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "state_comparators.hpp"
#include "test_utils.hpp"

// Allocations left before the injected failure; negative: no injection.
extern std::atomic<long> g_countdown;
// Live allocations of the replacement operator new (leak accounting).
extern std::atomic<long> g_live;

// Runs fn with the k-th allocation (from 0) failing; returns whether that failure happened
// (and expects that a reached failure propagated as std::bad_alloc).
bool fails_at(long k, const std::function<void()>& fn);

// As fails_at, counting only the allocations of threads other than the caller's (the
// evaluators' worker threads): the k-th worker allocation fails.
bool fails_in_workers_at(long k, const std::function<void()>& fn);

// As fails_in_workers_at, for a call whose lowest-index task exception is the one raised (the
// evaluators): the injected std::bad_alloc may lose to another task's exception, which fn
// handles. Returns whether the k-th worker allocation was reached, and sets `bad_alloc` when fn
// raised std::bad_alloc (which must then be the injected one).
bool reaches_in_workers_at(long k, const std::function<void()>& fn, bool& bad_alloc);

namespace netgraph::core::test {

// remove() and a policy's destruction allocate nothing: after place(3), remove() takes the
// policy's records away; after place(2), the destruction leaves them as orphans, the aggregate as
// it was and the policy no longer bound.
template <class Policy>
void expect_removal_and_release_allocate_nothing(const FlowLedger& fg, std::unique_ptr<Policy>& p,
                                                 const std::string& where) {
  (void)p->place(3);
  ASSERT_FALSE(fg.records().empty()) << where;  // something to remove
  EXPECT_FALSE(fails_at(0, [&] { p->remove(); })) << where;
  EXPECT_TRUE(fg.records().empty()) << where;
  (void)p->place(2);
  const LedgerState held = ledger_state(fg);
  ASSERT_FALSE(held.records.empty()) << where;  // something to orphan
  ASSERT_EQ(held.bound.size(), 1U) << where;
  EXPECT_FALSE(fails_at(0, [&] { p.reset(); })) << where;
  const LedgerState orphaned = ledger_state(fg);
  EXPECT_EQ(orphaned.allocated, held.allocated) << where;
  EXPECT_EQ(orphaned.records.size(), held.records.size()) << where;
  EXPECT_TRUE(orphaned.bound.empty()) << where;
  EXPECT_EQ(orphaned.orphans, held.bound) << where;
}

// What a placement sweep sets up on a fresh ledger before the call under test: the policy it
// calls, and what must stay alive with it (third parties, a wrapper of the policy).
struct Staged {
  std::shared_ptr<void> keep;
  FlowPolicy* policy {nullptr};
};

// The placement-atomicity sweep of one call. setup(ledger) stages everything before the call on a
// fresh ledger of g (third parties, the policy and its earlier calls); call(policy) is the call
// under test. The reference stages and calls once without a failure, then removes the policy.
// Attempt `at` (from 0) stages afresh (the same ledger state and epoch as the reference's) and
// makes the at-th allocation of the call fail; after the failure the ledger and the policy's epoch
// are as before the call and the aggregate is balanced, a clean retry gives the reference's result,
// ledger and epoch, and remove() leaves the ledger the reference's remove() left. Everything an
// attempt allocates is freed when it ends (the reference has already grown this thread's scratch),
// so a live allocation left after an attempt is a leak of the failed call. The sweep ends at the
// first attempt that completes, whose result, ledger and epoch must be the reference's, or at the first
// reported difference. Returns the failures injected.
inline long sweep_placement(const std::shared_ptr<const Graph>& g, const std::string& where,
                            const std::function<Staged(FlowLedger&)>& setup,
                            const std::function<PlacementResult(FlowPolicy&)>& call) {
  FlowLedger ref(g);
  const Staged staged_ref = setup(ref);
  const LedgerState initial = ledger_state(ref);
  const Epoch initial_epoch = epoch(*staged_ref.policy);
  const PlacementResult want = call(*staged_ref.policy);
  const LedgerState after = ledger_state(ref);
  const Epoch want_epoch = epoch(*staged_ref.policy);
  staged_ref.policy->remove();
  const LedgerState removed = ledger_state(ref);
  long injected = 0;
  for (long at = 0;; ++at) {
    const std::string w = where + " allocation " + std::to_string(at);
    const long live = g_live.load();
    bool completed = false;
    {
      FlowLedger fg(g);
      const Staged staged = setup(fg);
      FlowPolicy& p = *staged.policy;
      expect_same(ledger_state(fg), initial, w + " (setup)");
      expect_same(epoch(p), initial_epoch, w + " (setup)");
      PlacementResult got;
      if (!fails_at(at, [&] { got = call(p); })) {
        expect_same(got, want, w + " (no failure)");
        expect_same(ledger_state(fg), after, w + " (no failure)");
        expect_same(epoch(p), want_epoch, w + " (no failure)");
        completed = true;
      } else {
        ++injected;
        // Failed: nothing changed, ledger and policy.
        expect_same(ledger_state(fg), initial, w + " (failed call)");
        expect_same(epoch(p), initial_epoch, w + " (failed call)");
        expect_balanced(fg, w + " (failed call)");
        // A clean retry is the reference call (no poisoned per-thread scratch).
        expect_same(call(p), want, w + " (retry)");
        expect_same(ledger_state(fg), after, w + " (retry)");
        expect_same(epoch(p), want_epoch, w + " (retry)");
        expect_balanced(fg, w + " (retry)");
        p.remove();
        expect_same(ledger_state(fg), removed, w + " (remove)");
      }
    }
    EXPECT_EQ(g_live.load(), live) << w << ": a leak";
    if (completed || ::testing::Test::HasFailure()) break;
  }
  return injected;
}

// The scratch-growth sweep. Attempt `at` (from 0) runs on a new thread, whose per-thread scratch
// starts empty: warm_up() (a call on a small graph), then attempt(at, where), which makes the
// at-th allocation of the swept call on a larger graph fail, checks the failed call, compares its
// retry (or the call itself, when the failure is not reached) with a reference call made
// beforehand without a failure, and returns whether the failure was reached; the sweep ends at the
// first attempt that completes. Each attempt's thread frees its scratch when it ends, so a live
// allocation left after the join is a leak of the failed call. Then, on one more new thread,
// warm_up() and warm(), which grows the scratch to the larger graph and returns how many
// allocations the swept call still makes: fewer than the sweep's failures, so the sweep reached
// the growth. Returns the sweep's failures.
inline long sweep_scratch_growth(const std::string& what, const std::function<void()>& warm_up,
                                 const std::function<bool(long, const std::string&)>& attempt,
                                 const std::function<long()>& warm) {
  long failures = 0;
  for (long at = 0;; ++at) {
    bool failed = false;
    const std::string w = what + " allocation " + std::to_string(at);
    const long live = g_live.load();
    std::thread([&] {
      warm_up();
      failed = attempt(at, w);
    }).join();
    EXPECT_EQ(g_live.load(), live) << w << ": a leak";  // the thread's scratch is gone with it
    if (!failed) break;
    ++failures;
    if (::testing::Test::HasFailure() || at >= 1000) {
      ADD_FAILURE() << what << ": stopped at allocation " << at;
      return failures;
    }
  }
  long warm_count = 0;
  std::thread([&] {
    warm_up();
    warm_count = warm();
  }).join();
  EXPECT_GT(failures, warm_count) << what;
  return failures;
}

} // namespace netgraph::core::test
