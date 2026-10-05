// Allocation-failure atomicity of the ledger and the single-bundle transaction (APPROACH 8,
// 9.2; review-1b D1). This executable replaces the global operator new so that one chosen
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
#include "netgraph/core/detail/pinned_flow.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/flow_graph.hpp"
#include "test_utils.hpp"

namespace {
// Allocations left before the injected failure; negative: no injection.
long g_countdown = -1;
} // namespace

void* operator new(std::size_t n) {
  if (g_countdown >= 0 && g_countdown-- == 0) throw std::bad_alloc();
  if (void* p = std::malloc(n == 0 ? 1 : n)) return p;
  throw std::bad_alloc();
}
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

using namespace netgraph::core;
using namespace netgraph::core::test;
using detail::LedgerAccess;
using detail::PinnedFlow;
using detail::PinnedResult;

namespace {

// Runs fn with the k-th allocation (from 0) failing; returns whether that failure happened.
bool fails_at(long k, const std::function<void()>& fn) {
  g_countdown = k;
  bool failed = false;
  try {
    fn();
  } catch (const std::bad_alloc&) {
    failed = true;
  }
  const bool reached = g_countdown < 0;
  g_countdown = -1;
  EXPECT_EQ(failed, reached) << "a std::bad_alloc was swallowed or replaced at allocation " << k;
  return failed;
}

struct Snapshot {
  std::vector<Amount> allocated;
  std::vector<Amount> residual;
  std::vector<FlowRecord> records;
  std::int64_t next_id {0};
  std::vector<std::int64_t> bound;
  std::vector<std::int64_t> orphans;
};

// Records keep their local ids only: the reference sequence runs on another ledger (uid).
Snapshot snapshot(const FlowGraph& fg) {
  Snapshot s {fg.allocated(), fg.residual(), fg.records(), LedgerAccess::next_flow_id(fg),
              LedgerAccess::bound_owners(fg), LedgerAccess::orphan_owners(fg)};
  for (auto& r : s.records) {
    EXPECT_EQ(r.flow_id.first, fg.uid());
    r.flow_id.first = 0;
  }
  return s;
}

bool same(const FlowRecord& a, const FlowRecord& b) {
  return a.flow_id == b.flow_id && a.owner == b.owner && a.sources == b.sources && a.destinations == b.destinations &&
         a.placement == b.placement && a.admission == b.admission && a.bundle && b.bundle && *a.bundle == *b.bundle &&
         a.route_cost == b.route_cost && a.offered == b.offered && a.delivered == b.delivered && a.edges == b.edges &&
         a.drops == b.drops;
}

void expect_same(const Snapshot& a, const Snapshot& b, const std::string& where) {
  EXPECT_EQ(a.allocated, b.allocated) << where;
  EXPECT_EQ(a.residual, b.residual) << where;
  EXPECT_EQ(a.next_id, b.next_id) << where;
  EXPECT_EQ(a.bound, b.bound) << where;
  EXPECT_EQ(a.orphans, b.orphans) << where;
  ASSERT_EQ(a.records.size(), b.records.size()) << where;
  for (std::size_t i = 0; i < a.records.size(); ++i) EXPECT_TRUE(same(a.records[i], b.records[i])) << where << " record " << i;
}

// The ledger invariant: the aggregate is the sum of the records, within [0, baseline].
void expect_balanced(const FlowGraph& fg, const std::string& where) {
  std::vector<Amount> sum(fg.allocated().size(), 0);
  for (const auto& r : fg.records()) {
    for (const auto& [e, x] : r.edges) sum[static_cast<std::size_t>(e)] += x;
  }
  EXPECT_EQ(sum, fg.allocated()) << where;
  for (Amount x : fg.residual()) EXPECT_GE(x, 0) << where;
}

struct Epoch {
  Counter128 requested;
  bool infinite {false};
  std::optional<FlowId> flow_id;
  bool planned {false};
  Amount placed {0};
};

Epoch epoch(const PinnedFlow& p) {
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

void expect_same(const PinnedResult& a, const PinnedResult& b, const std::string& where) {
  EXPECT_EQ(a.placed, b.placed) << where;
  EXPECT_EQ(a.remaining, b.remaining) << where;
  EXPECT_EQ(a.termination, b.termination) << where;
  EXPECT_EQ(a.routes_changed, b.routes_changed) << where;
  EXPECT_EQ(a.drops, b.drops) << where;
  EXPECT_EQ(a.iterations, b.iterations) << where;
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
  std::shared_ptr<const StrictMultiDiGraph> g;
  NodeId src;
  NodeId dst;
};

std::vector<Case> cases() {
  auto shared = [](StrictMultiDiGraph g) { return std::make_shared<const StrictMultiDiGraph>(std::move(g)); };
  return {
      {"edge", shared(make_graph(2, {{0, 1, 9, 1}})), 0, 1},
      {"chain", shared(make_graph(3, {{0, 1, 9, 1}, {1, 2, 9, 1}})), 0, 2},
      // Split, reconvergence, a parallel pair and a costlier branch.
      {"diamond", shared(make_graph(5, {{0, 1, 7, 1}, {0, 2, 5, 1}, {1, 3, 4, 1}, {1, 3, 3, 1}, {2, 3, 9, 1},
                                        {3, 4, 20, 1}, {0, 4, 2, 5}})),
       0, 4},
  };
}

std::shared_ptr<const PredDAG> all_edges(const StrictMultiDiGraph& g) {
  std::vector<EdgeId> e(static_cast<std::size_t>(g.num_edges()));
  for (std::size_t i = 0; i < e.size(); ++i) e[i] = static_cast<EdgeId>(i);
  return std::make_shared<const PredDAG>(detail::PredDAGBuilder::from_edges(g, e));
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
        FlowGraph ref(c.g);
        PinnedFlow third_ref(ref, c.src, c.dst, dag, FlowPlacement::EqualBalanced, Admission::Lossless,
                             EdgePriority::Residual);
        (void)third_ref.place(1);
        PinnedFlow pref(ref, c.src, c.dst, dag, k.placement, k.admission, k.priority);
        for (auto v : sc.before) (void)pref.place(v);
        const Snapshot initial = snapshot(ref);
        const PinnedResult want = pref.place(sc.call);
        const Snapshot after = snapshot(ref);
        const Epoch want_epoch = epoch(pref);
        pref.remove();
        const Snapshot removed = snapshot(ref);
        ++sweeps;
        for (long at = 0;; ++at) {
          const std::string w = where + " allocation " + std::to_string(at);
          FlowGraph fg(c.g);
          PinnedFlow third(fg, c.src, c.dst, dag, FlowPlacement::EqualBalanced, Admission::Lossless,
                           EdgePriority::Residual);
          (void)third.place(1);
          PinnedFlow p(fg, c.src, c.dst, dag, k.placement, k.admission, k.priority);
          for (auto v : sc.before) (void)p.place(v);
          const Snapshot before = snapshot(fg);
          const Epoch before_epoch = epoch(p);
          expect_same(before, initial, w + " (setup)");
          PinnedResult got;
          if (!fails_at(at, [&] { got = p.place(sc.call); })) {
            expect_same(got, want, w + " (no failure)");
            expect_same(snapshot(fg), after, w + " (no failure)");
            break;
          }
          ++injected;
          // Failed: nothing changed, ledger and policy.
          expect_same(snapshot(fg), before, w + " (failed call)");
          expect_same(epoch(p), before_epoch, w + " (failed call)");
          expect_balanced(fg, w + " (failed call)");
          // A clean retry is the reference call (no poisoned per-thread scratch).
          expect_same(p.place(sc.call), want, w + " (retry)");
          expect_same(snapshot(fg), after, w + " (retry)");
          expect_same(epoch(p), want_epoch, w + " (retry)");
          expect_balanced(fg, w + " (retry)");
          p.remove();
          expect_same(snapshot(fg), removed, w + " (remove)");
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
        FlowGraph fg(c.g);
        {
          PinnedFlow other(fg, c.src, c.dst, dag, FlowPlacement::EqualBalanced, Admission::Lossless,
                           EdgePriority::Residual);
          (void)other.place(1);
        }  // an orphan
        const Snapshot before = snapshot(fg);
        std::optional<PinnedFlow> p;
        if (!fails_at(at, [&] { p.emplace(fg, c.src, c.dst, dag, k.placement, k.admission, k.priority); })) {
          EXPECT_EQ(LedgerAccess::bound_owners(fg).size(), 1U);
          break;
        }
        expect_same(snapshot(fg), before, c.name + " construction allocation " + std::to_string(at));
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
    auto g = std::make_shared<const StrictMultiDiGraph>(make_graph(n, {{n - 2, n - 1, 3, 1}}));
    const auto dag = all_edges(*g);
    FlowGraph fg(g);
    std::optional<PinnedFlow> p;
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
    auto g = std::make_shared<const StrictMultiDiGraph>(make_graph(n, chain));
    const auto dag = all_edges(*g);
    FlowGraph fg(g);
    PinnedFlow p(fg, 0, n - 1, dag, FlowPlacement::EqualBalanced, Admission::Lossless, EdgePriority::Residual);
    const Snapshot before = snapshot(fg);
    PinnedResult got;
    const bool failed = fails_at(at, [&] { got = p.place(2); });
    if (failed) {
      ++failures;
      expect_same(snapshot(fg), before, "placement allocation " + std::to_string(at));
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
// flow core's subset solve and in the plan builder. It must stay allocation-free (reset() sizes
// both arrays) so that it stays a few inlined instructions: the push_back version was an
// out-of-line call whose register save/restore stalled in some process layouts.
static_assert(noexcept(std::declval<detail::NodeMap&>().add(NodeId {0})));

TEST(LedgerFaults, NodeMapAddNeverAllocates) {
  auto& node_map = detail::NodeMap::local();
  // Called through a member pointer: the inventory scanner attributes this member call on a
  // receiver it does not list to a deleted 1a owner of the same member name.
  constexpr auto reset = &detail::NodeMap::reset;
  for (std::size_t n : {std::size_t {10}, std::size_t {5000}, std::size_t {300}}) {
    std::invoke(reset, node_map, n);
    EXPECT_FALSE(fails_at(0, [&] {
      for (std::size_t v = n; v-- > 0;) (void)node_map.add(static_cast<NodeId>(v));
      for (std::size_t v = 0; v < n; ++v) (void)node_map.add(static_cast<NodeId>(v));
    }));
    ASSERT_EQ(node_map.size(), n);
    EXPECT_EQ(node_map.of(static_cast<NodeId>(n - 1)), 0);  // first-touch order
    EXPECT_EQ(node_map.nodes()[0], static_cast<NodeId>(n - 1));
  }
  std::invoke(reset, node_map, std::size_t {0});
}

TEST(LedgerFaults, RemovalAndReleaseNeverAllocate) {
  for (const auto& c : cases()) {
    const auto dag = all_edges(*c.g);
    for (const auto& k : kernels()) {
      FlowGraph fg(c.g);
      auto p = std::make_unique<PinnedFlow>(fg, c.src, c.dst, dag, k.placement, k.admission, k.priority);
      (void)p->place(3);
      // remove() and the destructor's release() make no allocation (release() is noexcept).
      EXPECT_FALSE(fails_at(0, [&] { p->remove(); }));
      (void)p->place(2);
      const Snapshot held = snapshot(fg);
      EXPECT_FALSE(fails_at(0, [&] { p.reset(); }));
      const Snapshot orphaned = snapshot(fg);
      EXPECT_EQ(orphaned.allocated, held.allocated);
      EXPECT_TRUE(orphaned.bound.empty());
      EXPECT_EQ(orphaned.orphans.size(), 1U);
    }
  }
}
