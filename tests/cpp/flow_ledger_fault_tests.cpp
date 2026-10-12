// Allocation-failure atomicity of the ledger and of a one-bundle Pinned policy's transaction
// (APPROACH 8, 9.2; review-1b D1). This executable replaces the global operator new so that one chosen
// allocation throws std::bad_alloc. Every operation is swept over every allocation it makes:
// after each injected failure the ledger (aggregate, residual, every record field, bound and
// orphaned owners, next FlowId) and the policy (epoch state, flow id, plan) must equal their
// state before the call, the per-thread scratch must not be poisoned (a clean retry gives the
// reference result), and remove() must restore the state before the policy placed anything.
#include <gtest/gtest.h>

#include <cstdlib>
#include <functional>
#include <memory>
#include <new>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/node_map.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "fault_injection.hpp"
#include "pinned_one.hpp"
#include "test_utils.hpp"

// Allocations left before the injected failure; negative: no injection (fault_injection.hpp).
// Shared by every thread; with g_workers_only the calling test thread's allocations do not count.
std::atomic<long> g_countdown {-1};
// Live allocations of the replacement operator new (a successful allocation counts up, the
// delete of a non-null pointer counts down): the sweeps compare it before and after an operation
// to count leaks, which ASan's leak detection (off on macOS) would not.
std::atomic<long> g_live {0};
static std::atomic<bool> g_workers_only {false};
static thread_local bool t_test_thread = false;

void* operator new(std::size_t n) {
  if (g_countdown.load(std::memory_order_relaxed) >= 0 &&
      !(t_test_thread && g_workers_only.load(std::memory_order_relaxed)) &&
      g_countdown.fetch_sub(1, std::memory_order_relaxed) == 0) {
    throw std::bad_alloc();
  }
  if (void* p = std::malloc(n == 0 ? 1 : n)) {
    g_live.fetch_add(1, std::memory_order_relaxed);
    return p;
  }
  throw std::bad_alloc();
}
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept {
  if (p != nullptr) g_live.fetch_sub(1, std::memory_order_relaxed);
  std::free(p);
}
void operator delete[](void* p) noexcept { ::operator delete(p); }
void operator delete(void* p, std::size_t) noexcept { ::operator delete(p); }
void operator delete[](void* p, std::size_t) noexcept { ::operator delete(p); }

using namespace netgraph::core;
using namespace netgraph::core::test;
using detail::LedgerAccess;

namespace {

// The shared comparators of fault_injection.hpp, beside this file's own overloads.
using test::expect_same;

bool inject(long k, bool workers_only, const std::function<void()>& fn) {
  t_test_thread = true;
  g_workers_only = workers_only;
  g_countdown = k;
  bool failed = false;
  try {
    fn();
  } catch (const std::bad_alloc&) {
    failed = true;
  }
  const bool reached = g_countdown.load() < 0;
  g_countdown = -1;
  g_workers_only = false;
  EXPECT_EQ(failed, reached) << "a std::bad_alloc was swallowed or replaced at allocation " << k
                             << (workers_only ? " (worker threads)" : "");
  return failed;
}

} // namespace

// Runs fn with the k-th allocation (from 0) failing; returns whether that failure happened.
bool fails_at(long k, const std::function<void()>& fn) { return inject(k, false, fn); }

bool fails_in_workers_at(long k, const std::function<void()>& fn) { return inject(k, true, fn); }

namespace {

struct Epoch {
  Counter128 requested;
  bool infinite {false};
  std::optional<FlowId> flow_id;
  bool planned {false};
  Amount placed {0};
};

Epoch epoch(const PinnedOne& p) {
  auto id = p.flow_id();
  if (id.has_value()) id->first = 0;
  return {p.requested(), p.infinite(), id, p.plan() != nullptr, p.placed_demand()};
}

void expect_same(const Epoch& a, const Epoch& b, const std::string& where) {
  EXPECT_EQ(a.requested, b.requested) << where;
  EXPECT_EQ(a.infinite, b.infinite) << where;
  EXPECT_EQ(a.flow_id, b.flow_id) << where;
  EXPECT_EQ(a.planned, b.planned) << where;
  EXPECT_EQ(a.placed, b.placed) << where;
}

struct Kernel {
  FlowPlacement placement;
  Admission admission;
  EdgePriority priority;
};

const std::vector<Kernel>& kernels() {
  static const std::vector<Kernel> all = {
      {FlowPlacement::MaxFlow, Admission::Lossless, EdgePriority::Residual},
      {FlowPlacement::Weighted, Admission::Lossless, EdgePriority::Residual},
      {FlowPlacement::Weighted, Admission::Lossless, EdgePriority::Capacity},
      {FlowPlacement::EqualBalanced, Admission::Lossless, EdgePriority::Residual},
      {FlowPlacement::EqualBalanced, Admission::Lossy, EdgePriority::Residual},
      {FlowPlacement::Weighted, Admission::Lossy, EdgePriority::Capacity},
  };
  return all;
}

struct Case {
  std::string name;
  std::shared_ptr<const Graph> g;
  NodeId src;
  NodeId dst;
};

std::vector<Case> cases() {
  auto shared = [](Graph g) { return std::make_shared<const Graph>(std::move(g)); };
  return {
      {"edge", shared(make_graph(2, {{0, 1, 9, 1}})), 0, 1},
      {"chain", shared(make_graph(3, {{0, 1, 9, 1}, {1, 2, 9, 1}})), 0, 2},
      // Split, reconvergence, a parallel pair and a costlier branch.
      {"diamond", shared(make_graph(5, {{0, 1, 7, 1}, {0, 2, 5, 1}, {1, 3, 4, 1}, {1, 3, 3, 1}, {2, 3, 9, 1},
                                        {3, 4, 20, 1}, {0, 4, 2, 5}})),
       0, 4},
  };
}

std::shared_ptr<const PredDAG> all_edges(const Graph& g) {
  std::vector<EdgeId> e(static_cast<std::size_t>(g.num_edges()));
  for (std::size_t i = 0; i < e.size(); ++i) e[i] = static_cast<EdgeId>(i);
  return std::make_shared<const PredDAG>(detail::PredDAGBuilder::of_edges(g, e));
}

// One scenario: a third party holding 1 unit on every edge, then the policy's earlier calls,
// then the call under test.
struct Scenario {
  std::vector<std::optional<Amount>> before;
  std::optional<Amount> call;
  std::string name;
};

const std::vector<Scenario>& scenarios() {
  static const std::vector<Scenario> all = {
      {{}, 2, "create"},
      {{1}, 2, "update"},
      {{1}, std::nullopt, "update to +inf"},
      {{}, std::nullopt, "create at +inf"},
  };
  return all;
}

} // namespace

TEST(LedgerFaults, EveryAllocationOfAPlacementFailsAtomically) {
  int sweeps = 0;
  long injected = 0;
  for (const auto& c : cases()) {
    const auto dag = all_edges(*c.g);
    for (const auto& k : kernels()) {
      for (const auto& sc : scenarios()) {
        const std::string where = c.name + " kernel " + std::to_string(static_cast<int>(k.placement)) + "/" +
                                  std::to_string(static_cast<int>(k.admission)) + "/" +
                                  std::to_string(static_cast<int>(k.priority)) + " " + sc.name;
        // The reference: the same sequence without a failure.
        FlowLedger ref(c.g);
        PinnedOne third_ref(ref, c.src, c.dst, dag, FlowPlacement::EqualBalanced, Admission::Lossless,
                             EdgePriority::Residual);
        (void)third_ref.place(1);
        PinnedOne pref(ref, c.src, c.dst, dag, k.placement, k.admission, k.priority);
        for (auto v : sc.before) (void)pref.place(v);
        const LedgerState initial = ledger_state(ref);
        const PlacementResult want = pref.place(sc.call);
        const LedgerState after = ledger_state(ref);
        const Epoch want_epoch = epoch(pref);
        pref.remove();
        const LedgerState removed = ledger_state(ref);
        ++sweeps;
        for (long at = 0;; ++at) {
          const std::string w = where + " allocation " + std::to_string(at);
          FlowLedger fg(c.g);
          PinnedOne third(fg, c.src, c.dst, dag, FlowPlacement::EqualBalanced, Admission::Lossless,
                           EdgePriority::Residual);
          (void)third.place(1);
          PinnedOne p(fg, c.src, c.dst, dag, k.placement, k.admission, k.priority);
          for (auto v : sc.before) (void)p.place(v);
          const LedgerState before = ledger_state(fg);
          const Epoch before_epoch = epoch(p);
          expect_same(before, initial, w + " (setup)");
          PlacementResult got;
          if (!fails_at(at, [&] { got = p.place(sc.call); })) {
            expect_same(got, want, w + " (no failure)");
            expect_same(ledger_state(fg), after, w + " (no failure)");
            break;
          }
          ++injected;
          // Failed: nothing changed, ledger and policy.
          expect_same(ledger_state(fg), before, w + " (failed call)");
          expect_same(epoch(p), before_epoch, w + " (failed call)");
          expect_balanced(fg, w + " (failed call)");
          // A clean retry is the reference call (no poisoned per-thread scratch).
          expect_same(p.place(sc.call), want, w + " (retry)");
          expect_same(ledger_state(fg), after, w + " (retry)");
          expect_same(epoch(p), want_epoch, w + " (retry)");
          expect_balanced(fg, w + " (retry)");
          p.remove();
          expect_same(ledger_state(fg), removed, w + " (remove)");
          if (::testing::Test::HasFailure()) return;
        }
      }
    }
  }
  EXPECT_EQ(sweeps, 3 * 6 * 4);
  EXPECT_GT(injected, sweeps * 5);  // not vacuous: every placement allocates
}

TEST(LedgerFaults, PolicyConstructionFailsWithoutBinding) {
  for (const auto& c : cases()) {
    const auto dag = all_edges(*c.g);
    for (const auto& k : kernels()) {
      for (long at = 0;; ++at) {
        FlowLedger fg(c.g);
        {
          PinnedOne other(fg, c.src, c.dst, dag, FlowPlacement::EqualBalanced, Admission::Lossless,
                           EdgePriority::Residual);
          (void)other.place(1);
        }  // an orphan
        const LedgerState before = ledger_state(fg);
        std::optional<PinnedOne> p;
        if (!fails_at(at, [&] { p.emplace(fg, c.src, c.dst, dag, k.placement, k.admission, k.priority); })) {
          EXPECT_EQ(LedgerAccess::bound_owners(fg).size(), 1U);
          break;
        }
        expect_same(ledger_state(fg), before, c.name + " construction allocation " + std::to_string(at));
        EXPECT_FALSE(p.has_value());
      }
    }
  }
}

// review-1b-fixes N1: the per-thread scratch of the plan builders and of the commit grows with
// the graph. Every position fails on a graph larger than any before it in this process, so the
// failure lands while the scratch grows; the failed call changes nothing and a clean retry on the
// same graph works (ASan/UBSan builds report any access outside a partly grown array).
TEST(LedgerFaults, ScratchGrowthFailureLeavesNothingHalfGrown) {
  int failures = 0;
  std::int32_t n = 200;
  // Construction (the route plan's per-node scratch): one edge between the last two nodes.
  for (long at = 0;; ++at, n += 150) {
    auto g = std::make_shared<const Graph>(make_graph(n, {{n - 2, n - 1, 3, 1}}));
    const auto dag = all_edges(*g);
    FlowLedger fg(g);
    std::optional<PinnedOne> p;
    const bool failed = fails_at(at, [&] {
      p.emplace(fg, n - 2, n - 1, dag, FlowPlacement::EqualBalanced, Admission::Lossless, EdgePriority::Residual);
    });
    if (failed) {
      ++failures;
      EXPECT_FALSE(p.has_value());
      EXPECT_TRUE(LedgerAccess::bound_owners(fg).empty());
      p.emplace(fg, n - 2, n - 1, dag, FlowPlacement::EqualBalanced, Admission::Lossless, EdgePriority::Residual);
    }
    EXPECT_EQ(p->place(2).placed, 2) << "construction allocation " << at;
    expect_balanced(fg, "construction allocation " + std::to_string(at));
    p->remove();
    EXPECT_EQ(fg.allocated()[0], 0);
    if (!failed) break;
  }
  // Placement (the commit's per-edge scratch and the kernels' buffers): a chain of n - 1 edges.
  for (long at = 0;; ++at, n += 150) {
    std::vector<InEdge> chain;
    for (NodeId v = 0; v + 1 < n; ++v) chain.push_back({v, v + 1, 5, 1});
    auto g = std::make_shared<const Graph>(make_graph(n, chain));
    const auto dag = all_edges(*g);
    FlowLedger fg(g);
    PinnedOne p(fg, 0, n - 1, dag, FlowPlacement::EqualBalanced, Admission::Lossless, EdgePriority::Residual);
    const LedgerState before = ledger_state(fg);
    PlacementResult got;
    const bool failed = fails_at(at, [&] { got = p.place(2); });
    if (failed) {
      ++failures;
      expect_same(ledger_state(fg), before, "placement allocation " + std::to_string(at));
      got = p.place(2);
    }
    EXPECT_EQ(got.placed, 2) << "placement allocation " << at;
    expect_balanced(fg, "placement allocation " + std::to_string(at));
    p.remove();  // ids are not reused: only the aggregate and the records return
    EXPECT_EQ(fg.allocated(), before.allocated) << "placement allocation " << at;
    EXPECT_TRUE(fg.records().empty());
    if (!failed) break;
  }
  EXPECT_GT(failures, 4);  // not vacuous
}

// The bimodality guard (NOTES "Bimodality"): NodeMap::add runs once per edge endpoint in the
// flow core's subset solve. It must stay allocation-free (reset() sizes both arrays) so that it
// stays a few inlined instructions: the push_back version was an out-of-line call whose register
// save/restore stalled in some process layouts.
static_assert(noexcept(std::declval<detail::NodeMap&>().add(NodeId {0})));

TEST(LedgerFaults, NodeMapAddNeverAllocates) {
  auto& node_map = detail::NodeMap::local();
  for (std::size_t n : {std::size_t {10}, std::size_t {5000}, std::size_t {300}}) {
    node_map.reset(n);
    EXPECT_FALSE(fails_at(0, [&] {
      for (std::size_t v = n; v-- > 0;) (void)node_map.add(static_cast<NodeId>(v));
      for (std::size_t v = 0; v < n; ++v) (void)node_map.add(static_cast<NodeId>(v));
    }));
    ASSERT_EQ(node_map.size(), n);
    EXPECT_EQ(node_map.of(static_cast<NodeId>(n - 1)), 0);  // first-touch order
    EXPECT_EQ(node_map.of(0), static_cast<std::int32_t>(n - 1));
  }
  node_map.reset(0);
}

TEST(LedgerFaults, RemovalAndReleaseNeverAllocate) {
  for (const auto& c : cases()) {
    const auto dag = all_edges(*c.g);
    for (const auto& k : kernels()) {
      FlowLedger fg(c.g);
      auto p = std::make_unique<PinnedOne>(fg, c.src, c.dst, dag, k.placement, k.admission, k.priority);
      (void)p->place(3);
      // remove() and the destructor (the records become orphans) make no allocation.
      EXPECT_FALSE(fails_at(0, [&] { p->remove(); }));
      (void)p->place(2);
      const LedgerState held = ledger_state(fg);
      EXPECT_FALSE(fails_at(0, [&] { p.reset(); }));
      const LedgerState orphaned = ledger_state(fg);
      EXPECT_EQ(orphaned.allocated, held.allocated);
      EXPECT_TRUE(orphaned.bound.empty());
      EXPECT_EQ(orphaned.orphans.size(), 1U);
    }
  }
}
