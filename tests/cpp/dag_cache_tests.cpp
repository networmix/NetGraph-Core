// The DagCache's modes (APPROACH 10): the transition table (building and frozen lookups, freeze,
// clear, thaw, destruction), the statistics (a frozen miss counts one miss), DagCacheMiss as a
// dedicated exception, frozen hits from several threads, and the evaluators' sharing rules
// (section 11): a frozen cache shared by every worker, a building one rejected with two or more,
// DagCacheMiss propagated unchanged through batch_capacity and sensitivity, the reuse of a
// prepopulated Eligibility through an unchanged intersection, and the deletion-mode sensitivity
// with a supplied frozen cache. The thread cases are named DagCacheFrozen* (the TSan filter).
#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <vector>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

static_assert(!std::is_base_of_v<std::runtime_error, DagCacheMiss>);
static_assert(!std::is_base_of_v<std::logic_error, DagCacheMiss>);
static_assert(std::is_base_of_v<std::exception, DagCacheMiss>);
static_assert(std::is_nothrow_copy_constructible_v<DagCacheMiss>);

namespace {

// A small two-tier fabric: leaves 0..5, spines 6..8, every leaf joined to every spine both ways.
std::shared_ptr<const Graph> fabric() {
  std::vector<InEdge> e;
  for (NodeId l = 0; l < 6; ++l) {
    for (NodeId s = 6; s < 9; ++s) {
      e.push_back({l, s, 4, 1});
      e.push_back({s, l, 4, 1});
    }
  }
  return std::make_shared<const Graph>(make_graph(9, e));
}

std::vector<NodeId> v(std::initializer_list<NodeId> x) { return x; }

FlowPolicyConfig costonly() {
  return FlowPolicyConfig {CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless,
                           PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity)};
}

FlowPolicyConfig global() {
  return FlowPolicyConfig {Global {}, FlowPlacement::MaxFlow, Admission::Lossless,
                           PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity)};
}

std::string miss_of(const std::function<void()>& fn) {
  try {
    fn();
  } catch (const DagCacheMiss& e) {
    return e.what();
  }
  return "";
}

} // namespace

TEST(DagCacheModes, TransitionTableAndStatistics) {
  const auto g = fabric();
  DagCache cache(*g);
  // Building: a miss searches and inserts, a hit returns the same value.
  const auto a = cache.lookup(*g, v({1, 0}), Direction::Reverse, nullptr);
  const auto a2 = cache.lookup(*g, v({0, 1}), Direction::Reverse, nullptr);
  EXPECT_EQ(a.get(), a2.get());
  auto st = cache.stats();
  EXPECT_EQ(st.misses, 1);
  EXPECT_EQ(st.hits, 1);
  EXPECT_EQ(st.entries, 1);
  // Frozen: hits return the value; a miss counts one miss, inserts nothing, throws DagCacheMiss.
  cache.freeze();
  cache.freeze();  // idempotent
  EXPECT_EQ(cache.lookup(*g, v({0, 1}), Direction::Reverse, nullptr).get(), a.get());
  const auto msg = miss_of([&] { (void)cache.lookup(*g, v({2}), Direction::Forward, nullptr); });
  EXPECT_EQ(msg, "DagCache: the frozen cache has no entry for roots {2}, direction Forward, Eligibility uid 0 (a frozen "
                 "cache never searches)");
  st = cache.stats();
  EXPECT_EQ(st.misses, 2);
  EXPECT_EQ(st.hits, 2);
  EXPECT_EQ(st.entries, 1);
  // Invalid queries throw std::invalid_argument and count nothing, in frozen mode too.
  const auto other = fabric();
  EXPECT_THROW((void)cache.lookup(*other, v({0}), Direction::Reverse, nullptr), std::invalid_argument);
  EXPECT_THROW((void)cache.lookup(*g, v({}), Direction::Reverse, nullptr), std::invalid_argument);
  EXPECT_THROW((void)cache.lookup(*g, v({9}), Direction::Reverse, nullptr), std::invalid_argument);
  EXPECT_EQ(cache.stats().misses, 2);
  EXPECT_EQ(cache.stats().hits, 2);
  // clear() throws when frozen and keeps the entries.
  try {
    cache.clear();
    FAIL() << "clear() on a frozen cache";
  } catch (const std::invalid_argument& e) {
    EXPECT_EQ(std::string(e.what()),
              "DagCache.clear is rejected while the cache is frozen (join every reader, then thaw() and clear())");
  }
  EXPECT_EQ(cache.stats().entries, 1);
  // thaw(): building again; a miss searches and inserts, clear() empties (counters persist).
  cache.thaw();
  cache.thaw();  // no effect when building
  const auto b = cache.lookup(*g, v({2}), Direction::Forward, nullptr);
  EXPECT_EQ(cache.stats().entries, 2);
  cache.clear();
  st = cache.stats();
  EXPECT_EQ(st.entries, 0);
  EXPECT_EQ(st.bytes, 0);
  EXPECT_EQ(st.misses, 3);
  EXPECT_EQ(st.hits, 2);
  EXPECT_EQ(st.live_values, 2);  // a and b are held by the test
  // Held values stay valid after clear() and after a frozen cache is destroyed.
  EXPECT_EQ(b->via_edges().size(), 18U);  // roots {2} forward: three 2 -> spine edges, 15 spine -> leaf
  EXPECT_EQ(a2->via_edges().size(), 18U);
  std::shared_ptr<const PredDAG> c;
  {
    DagCache frozen(*g);
    c = frozen.lookup(*g, v({3}), Direction::Reverse, nullptr);
    frozen.freeze();
  }
  // Roots {3}: the three spine -> 3 edges and the five other leaves' 15 leaf -> spine edges;
  // roots {0, 1}: six spine -> root edges and the four other leaves' 12.
  EXPECT_EQ(c->via_edges().size(), 18U);
  EXPECT_EQ(a->via_edges().size(), 18U);
}

TEST(DagCacheModes, RootsInAnyOrderShareOneEntryAndOneValidation) {
  // Task 1d-msd: a lookup borrows roots that are already a sorted set (no copy on a hit) and
  // normalises any other; either way the key is the sorted set, ordered by roots, then direction,
  // then Eligibility uid, and the checks and their messages are node_set's.
  const auto g = fabric();
  DagCache cache(*g);
  const auto elig = g->eligibility(v({}), std::vector<EdgeId> {0});
  const auto a = cache.lookup(*g, v({0, 2}), Direction::Reverse, nullptr);
  for (const auto& roots : {v({2, 0}), v({0, 0, 2}), v({2, 2, 0, 2}), v({0, 2})}) {
    EXPECT_EQ(cache.lookup(*g, roots, Direction::Reverse, nullptr).get(), a.get());
  }
  EXPECT_EQ(cache.stats().misses, 1);
  EXPECT_EQ(cache.stats().hits, 4);
  // Keys that differ in a prefix, the length, the direction or the Eligibility uid are distinct
  // entries, each found again by its own key in any order.
  const std::vector<std::tuple<std::vector<NodeId>, Direction, const Eligibility*>> keys = {
      {v({0}), Direction::Reverse, nullptr},        {v({0, 1}), Direction::Reverse, nullptr},
      {v({0, 2, 3}), Direction::Reverse, nullptr},  {v({1}), Direction::Reverse, nullptr},
      {v({0, 2}), Direction::Forward, nullptr},     {v({0, 2}), Direction::Reverse, &elig},
      {v({8}), Direction::Forward, &elig},          {v({0, 1, 2, 3, 4, 5}), Direction::Reverse, nullptr}};
  std::vector<std::shared_ptr<const PredDAG>> values;
  for (const auto& [roots, dir, el] : keys) values.push_back(cache.lookup(*g, roots, dir, el));
  EXPECT_EQ(cache.stats().entries, 9);
  EXPECT_EQ(cache.stats().misses, 9);
  for (std::size_t i = 0; i < keys.size(); ++i) {
    EXPECT_NE(values[i].get(), a.get());
    for (std::size_t j = 0; j < i; ++j) EXPECT_NE(values[i].get(), values[j].get());
  }
  for (std::size_t i = keys.size(); i-- > 0;) {
    auto roots = std::get<0>(keys[i]);
    std::reverse(roots.begin(), roots.end());
    EXPECT_EQ(cache.lookup(*g, roots, std::get<1>(keys[i]), std::get<2>(keys[i])).get(), values[i].get());
  }
  EXPECT_EQ(cache.stats().misses, 9);
  EXPECT_EQ(cache.stats().hits, 4 + static_cast<std::int64_t>(keys.size()));
  // The checks in node_set's order, sorted or not: empty, then the first id out of range in input
  // order; nothing counted.
  auto error_of = [&](std::vector<NodeId> roots) {
    try {
      (void)cache.lookup(*g, roots, Direction::Reverse, nullptr);
    } catch (const std::invalid_argument& e) {
      return std::string(e.what());
    }
    return std::string();
  };
  const auto before = cache.stats();
  EXPECT_EQ(error_of(v({9})), "roots: node 9 out of node range [0, 9)");
  EXPECT_EQ(error_of(v({3, 9})), "roots: node 9 out of node range [0, 9)");
  EXPECT_EQ(error_of(v({9, 3})), "roots: node 9 out of node range [0, 9)");
  EXPECT_EQ(error_of(v({-1, 9})), "roots: node -1 out of node range [0, 9)");
  EXPECT_EQ(error_of(v({9, -1})), "roots: node 9 out of node range [0, 9)");
  EXPECT_EQ(error_of(v({3, 3, 12})), "roots: node 12 out of node range [0, 9)");
  EXPECT_FALSE(error_of(v({})).empty());
  const auto after = cache.stats();
  EXPECT_EQ(after.misses, before.misses);
  EXPECT_EQ(after.hits, before.hits);
  EXPECT_EQ(after.entries, before.entries);
  // A frozen miss names the sorted set, borrowed or normalised.
  cache.freeze();
  const std::string want = "DagCache: the frozen cache has no entry for roots {1, 4}, direction Reverse, Eligibility uid "
                           "0 (a frozen cache never searches)";
  EXPECT_EQ(miss_of([&] { (void)cache.lookup(*g, v({1, 4}), Direction::Reverse, nullptr); }), want);
  EXPECT_EQ(miss_of([&] { (void)cache.lookup(*g, v({4, 1, 4}), Direction::Reverse, nullptr); }), want);
  EXPECT_EQ(cache.stats().misses, 11);
}

TEST(DagCacheModes, MissIsADedicatedException) {
  const auto g = fabric();
  DagCache cache(*g);
  cache.freeze();
  bool as_miss = false;
  try {
    (void)cache.lookup(*g, v({0}), Direction::Reverse, nullptr);
  } catch (const std::runtime_error&) {
    FAIL() << "DagCacheMiss caught as std::runtime_error";
  } catch (const std::logic_error&) {
    FAIL() << "DagCacheMiss caught as std::logic_error";
  } catch (const DagCacheMiss&) {
    as_miss = true;
  }
  EXPECT_TRUE(as_miss);
}

namespace {

// Frozen hits from `threads` threads at once over prepopulated keys; every value is the
// prepopulated one and every hit is counted.
void frozen_readers(int threads) {
  const auto g = fabric();
  DagCache cache(*g);
  std::vector<std::shared_ptr<const PredDAG>> want;
  for (NodeId r = 0; r < 9; ++r) want.push_back(cache.lookup(*g, v({r}), Direction::Reverse, nullptr));
  cache.freeze();
  const auto before = cache.stats();
  std::atomic<int> wrong {0};
  std::vector<std::thread> pool;
  constexpr int kRounds = 200;
  for (int t = 0; t < threads; ++t) {
    pool.emplace_back([&, t] {
      for (int round = 0; round < kRounds; ++round) {
        const NodeId r = static_cast<NodeId>((round + t) % 9);
        const NodeId roots[1] = {r};
        if (cache.lookup(*g, roots, Direction::Reverse, nullptr).get() != want[static_cast<std::size_t>(r)].get()) ++wrong;
      }
    });
  }
  for (auto& th : pool) th.join();
  EXPECT_EQ(wrong.load(), 0);
  const auto after = cache.stats();
  EXPECT_EQ(after.hits - before.hits, static_cast<std::int64_t>(threads) * kRounds);
  EXPECT_EQ(after.misses, before.misses);
  EXPECT_EQ(after.entries, 9);
  // Joined: thaw and build again.
  cache.thaw();
  (void)cache.lookup(*g, v({0, 1}), Direction::Reverse, nullptr);
  EXPECT_EQ(cache.stats().entries, 10);
}

} // namespace

TEST(DagCacheFrozen, HitsFromOneThread) { frozen_readers(1); }
TEST(DagCacheFrozen, HitsFromTwoThreads) { frozen_readers(2); }
TEST(DagCacheFrozen, HitsFromFourThreads) { frozen_readers(4); }

TEST(DagCacheFrozen, SharedByEveryWorkerOfABatch) {
  const auto g = fabric();
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {v({0}), v({3})}, {v({1}), v({4})}, {v({2}), v({5})}, {v({0, 1}), v({3, 4})}};
  const auto want = batch_capacity(*g, pairs, CapacityOptions {costonly()}, {}, 1);
  for (int workers : {1, 2, 4}) {
    DagCache cache(*g);
    for (const auto& [s, t] : pairs) (void)cache.lookup(*g, t, Direction::Reverse, nullptr);
    // A building cache with two or more workers is rejected; frozen, it is shared.
    if (workers > 1) {
      EXPECT_THROW((void)batch_capacity(*g, pairs, CapacityOptions {costonly(), nullptr, &cache}, {}, workers),
                   std::invalid_argument);
    }
    cache.freeze();
    const auto before = cache.stats();
    const auto got = batch_capacity(*g, pairs, CapacityOptions {costonly(), nullptr, &cache}, {}, workers);
    ASSERT_EQ(got.size(), want.size());
    for (std::size_t i = 0; i < got.size(); ++i) EXPECT_EQ(got[i].placed, want[i].placed) << workers;
    EXPECT_EQ(cache.stats().hits - before.hits, static_cast<std::int64_t>(pairs.size()));
    EXPECT_EQ(cache.stats().misses, before.misses);
  }
}

TEST(DagCacheFrozen, MissPropagatesThroughBatchCapacity) {
  // An empty frozen cache: every pair misses; the lowest index's DagCacheMiss reaches the caller
  // unchanged, with one worker and with several.
  const auto g = fabric();
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{v({0}), v({3})}, {v({1}), v({4})},
                                                                                   {v({2}), v({5})}};
  for (int workers : {1, 2, 4}) {
    DagCache cache(*g);
    (void)cache.lookup(*g, v({4}), Direction::Reverse, nullptr);  // only pair 1's key
    cache.freeze();
    const auto msg = miss_of([&] { (void)batch_capacity(*g, pairs, CapacityOptions {costonly(), nullptr, &cache}, {}, workers); });
    EXPECT_NE(msg.find("roots {3}"), std::string::npos) << workers << ": " << msg;  // pair 0, the lowest index
  }
}

TEST(DagCacheFrozen, PerPairIntersectionsAndPrepopulatedUids) {
  // APPROACH 10: an intersection equal to an operand returns that operand (its uid), so the
  // prepopulated key hits; one that differs from both operands gets a fresh uid and misses.
  const auto g = fabric();
  const EdgeId e0[1] = {0};
  const EdgeId e01[2] = {0, 1};
  const EdgeId e1[1] = {1};
  const auto common = g->eligibility({}, e0);    // excludes edge 0
  const auto narrow = g->eligibility({}, e01);   // excludes 0 and 1: contained in common
  const auto other = g->eligibility({}, e1);     // excludes 1: neither contains the other
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{v({1}), v({4})}, {v({2}), v({5})}};
  DagCache cache(*g);
  (void)cache.lookup(*g, v({4}), Direction::Reverse, &narrow);
  (void)cache.lookup(*g, v({5}), Direction::Reverse, &common);
  (void)cache.lookup(*g, v({4}), Direction::Reverse, &common);
  cache.freeze();
  const auto hits = [&] { return cache.stats().hits; };
  // A common value containing the per-pair one: intersect returns the per-pair operand.
  auto before = hits();
  const std::vector<const Eligibility*> contained = {&narrow, nullptr};
  EXPECT_EQ(batch_capacity(*g, pairs, CapacityOptions {costonly(), &common, &cache}, contained, 2).size(), 2U);
  EXPECT_EQ(hits() - before, 2);
  // An unchanged intersection (the operand itself) and already-effective values with a null
  // common value.
  before = hits();
  const std::vector<const Eligibility*> same = {&common, &common};
  EXPECT_EQ(batch_capacity(*g, pairs, CapacityOptions {costonly(), &common, &cache}, same, 2).size(), 2U);
  const std::vector<const Eligibility*> own = {&narrow, &common};
  EXPECT_EQ(batch_capacity(*g, pairs, CapacityOptions {costonly(), nullptr, &cache}, own, 2).size(), 2U);
  EXPECT_EQ(hits() - before, 4);
  // An intersection that differs from both operands: a fresh uid, absent from the cache.
  const auto misses = cache.stats().misses;
  const std::vector<const Eligibility*> changed = {&other, nullptr};
  const auto msg = miss_of([&] { (void)batch_capacity(*g, pairs, CapacityOptions {costonly(), &common, &cache}, changed, 2); });
  EXPECT_NE(msg.find("roots {4}"), std::string::npos) << msg;
  EXPECT_EQ(msg.find("uid " + std::to_string(common.uid()) + " "), std::string::npos) << msg;
  EXPECT_EQ(msg.find("uid " + std::to_string(other.uid()) + " "), std::string::npos) << msg;
  EXPECT_EQ(cache.stats().misses - misses, 1);
  EXPECT_EQ(cache.stats().entries, 3);
}

TEST(DagCacheFrozen, MissPropagatesThroughSensitivityNeverAsUnevaluable) {
  // CostOnly sensitivity evaluates each candidate under a new Eligibility (a new uid): with a
  // frozen cache holding only the baseline's key, the first re-evaluated deletion misses and the
  // miss reaches the caller (it is not a std::runtime_error, so never `unevaluable`).
  const auto g = fabric();
  for (int workers : {1, 2, 4}) {
    DagCache cache(*g);
    (void)cache.lookup(*g, v({4}), Direction::Reverse, nullptr);
    cache.freeze();
    const auto msg = miss_of([&] { (void)sensitivity(*g, v({1}), v({4}), SensitivityOptions {costonly(), nullptr, &cache}, workers); });
    EXPECT_NE(msg.find("no entry for roots {4}"), std::string::npos) << workers << ": " << msg;
    EXPECT_EQ(cache.stats().entries, 1);  // nothing inserted
  }
}

TEST(DagCacheFrozen, DeletionModeSensitivityWithASuppliedFrozenCache) {
  // Global sensitivity repairs the baseline's state per deletion and reads no cache: a supplied
  // frozen cache, shared by two or four workers, gives the result without a cache.
  const auto g = fabric();
  const auto want = sensitivity(*g, v({0, 1}), v({4, 5}), SensitivityOptions {global()}, 1);
  for (int workers : {1, 2, 4}) {
    DagCache cache(*g);
    cache.freeze();
    const auto got = sensitivity(*g, v({0, 1}), v({4, 5}), SensitivityOptions {global(), nullptr, &cache}, workers);
    EXPECT_EQ(got.baseline.placed, want.baseline.placed);
    ASSERT_EQ(got.entries.size(), want.entries.size());
    for (std::size_t i = 0; i < got.entries.size(); ++i) {
      EXPECT_EQ(got.entries[i].edge, want.entries[i].edge);
      EXPECT_EQ(got.entries[i].loss, want.entries[i].loss);
    }
    EXPECT_EQ(cache.stats().hits + cache.stats().misses, 0);
  }
}

namespace {

std::tuple<std::int64_t, std::int64_t, std::int64_t, std::int64_t, std::int64_t> all_of(const DagCache::Stats& s) {
  return {s.hits, s.misses, s.entries, s.bytes, s.live_values};
}

} // namespace

TEST(DagCacheFrozen, AFailedCallChangesNoStatisticButItsMiss) {
  // review-1d F1, NOTES R27: a public call reading a frozen cache publishes its hits only when it
  // succeeds. A batch whose pair 0 hits and pair 1 misses raises the miss and keeps no hit; the
  // miss counts once, with one worker and with several (the boundary counts it, not the workers);
  // the policy constructor and CostOnly sensitivity do the same.
  const auto g = fabric();
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {v({0}), v({3})}, {v({1}), v({4})}, {v({2}), v({5})}, {v({0}), v({3})}};
  for (int workers : {1, 2, 4}) {
    DagCache cache(*g);
    (void)cache.lookup(*g, v({3}), Direction::Reverse, nullptr);
    cache.freeze();
    auto before = all_of(cache.stats());
    // Pairs 1 and 2 miss; pairs 0 and 3 hit.
    EXPECT_NE(miss_of([&] { (void)batch_capacity(*g, pairs, CapacityOptions {costonly(), nullptr, &cache}, {}, workers); })
                  .find("roots {4}"),
              std::string::npos);
    auto want = before;
    std::get<1>(want) += 1;
    EXPECT_EQ(all_of(cache.stats()), want) << workers;
    before = all_of(cache.stats());
    EXPECT_NE(miss_of([&] {
                (void)sensitivity(*g, v({0}), v({3}), SensitivityOptions {costonly(), nullptr, &cache}, workers);
              }).find("roots {3}"),
              std::string::npos);
    want = before;
    std::get<1>(want) += 1;
    EXPECT_EQ(all_of(cache.stats()), want) << workers;  // the baseline's hit is dropped
  }
  DagCache cache(*g);
  cache.freeze();
  FlowLedger fg(g);
  EXPECT_FALSE(miss_of([&] { FlowPolicy p(fg, v({0}), v({3}), costonly(), nullptr, &cache); }).empty());
  EXPECT_EQ(all_of(cache.stats()), std::make_tuple(0, 1, 0, 0, 0));
}

TEST(DagCacheFrozen, ConcurrentCallsKeepEverySuccessfulHit) {
  // Calls on one frozen cache from four threads at once, each a successful batch (two hits), a
  // failing batch (its hit dropped, one miss) and a successful capacity (one hit): every
  // successful call's hits are kept, no failed call's, and every failed call's miss is counted.
  const auto g = fabric();
  DagCache cache(*g);
  for (NodeId d : {3, 4}) (void)cache.lookup(*g, v({d}), Direction::Reverse, nullptr);
  cache.freeze();
  const auto before = cache.stats();
  const CapacityOptions opts {costonly(), nullptr, &cache};
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> good = {{v({0}), v({3})}, {v({1}), v({4})}};
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> bad = {{v({0}), v({3})}, {v({1}), v({5})}};
  constexpr int kThreads = 4;
  constexpr int kRounds = 25;
  std::atomic<int> wrong {0};
  std::vector<std::thread> pool;
  for (int t = 0; t < kThreads; ++t) {
    pool.emplace_back([&] {
      for (int round = 0; round < kRounds; ++round) {
        if (batch_capacity(*g, good, opts, {}, 2).size() != 2) ++wrong;
        try {
          (void)batch_capacity(*g, bad, opts, {}, 2);
          ++wrong;
        } catch (const DagCacheMiss&) {
        }
        if (capacity(*g, v({2}), v({3}), opts).placed <= 0) ++wrong;
      }
    });
  }
  for (auto& th : pool) th.join();
  EXPECT_EQ(wrong.load(), 0);
  const auto after = cache.stats();
  EXPECT_EQ(after.hits - before.hits, static_cast<std::int64_t>(kThreads) * kRounds * 3);
  EXPECT_EQ(after.misses - before.misses, static_cast<std::int64_t>(kThreads) * kRounds);
  EXPECT_EQ(after.entries, 2);
  EXPECT_EQ(after.bytes, before.bytes);
  EXPECT_EQ(after.live_values, before.live_values);
}

TEST(DagCacheFrozen, ConcurrentPublicationsLoseNoIncrement) {
  // Many short public calls at once (a CostOnly policy constructor: one frozen hit, published as
  // its last step), so that publications from different threads overlap often: the hits are
  // added atomically, none is lost.
  const auto g = fabric();
  DagCache cache(*g);
  (void)cache.lookup(*g, v({3}), Direction::Reverse, nullptr);
  cache.freeze();
  const auto before = cache.stats().hits;
  constexpr int kThreads = 8;
  constexpr int kRounds = 1500;
  std::vector<std::thread> pool;
  for (int t = 0; t < kThreads; ++t) {
    pool.emplace_back([&] {
      FlowLedger fg(g);
      for (int round = 0; round < kRounds; ++round) FlowPolicy p(fg, v({0}), v({3}), costonly(), nullptr, &cache);
    });
  }
  for (auto& th : pool) th.join();
  EXPECT_EQ(cache.stats().hits - before, static_cast<std::int64_t>(kThreads) * kRounds);
}
