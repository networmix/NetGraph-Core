// The integer ledger FlowLedger (APPROACH 8): the baseline, snapshots, FlowIds, bottleneck, and
// every row and transition of the lifetime table, driven through a FlowPolicy with Pinned
// routing and one bundle (test::PinnedOne, pinned_one.hpp).
//
// | State             | place                         | remove()           | destroyed  | reset    | copy     |
// | bound, no flows   | place(0) no-op; positive      | resets the epoch   | none       | rejected | rejected |
// |                   | creates (no route: nothing)   |                    |            |          |          |
// | bound, flows      | adds or replaces              | reverts, retires   | orphaned   | rejected | rejected |
// | orphaned records  | n/a                           | n/a                | n/a        | clears   | copied   |
#include <gtest/gtest.h>

#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "pinned_one.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using detail::LedgerAccess;
using netgraph::core::Termination;

namespace {

std::shared_ptr<const Graph> shared(Graph g) {
  return std::make_shared<const Graph>(std::move(g));
}

std::shared_ptr<const PredDAG> bundle(const Graph& g, std::vector<EdgeId> edges) {
  return std::make_shared<const PredDAG>(detail::PredDAGBuilder::of_edges(g, edges));
}

std::shared_ptr<const PredDAG> all_edges(const Graph& g) {
  std::vector<EdgeId> e(static_cast<std::size_t>(g.num_edges()));
  for (std::size_t i = 0; i < e.size(); ++i) e[i] = static_cast<EdgeId>(i);
  return bundle(g, e);
}

std::string error_of(const std::function<void()>& fn) {
  try {
    fn();
  } catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
}

std::vector<Amount> alloc_in(const FlowLedger& fg) { return to_input(*fg.graph(), fg.allocated()); }

} // namespace

TEST(FlowLedger, BaselineValidationAndSnapshots) {
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  EXPECT_THROW(FlowLedger(g, std::vector<Amount> {6}), std::invalid_argument);
  EXPECT_THROW(FlowLedger(g, std::vector<Amount> {-1}), std::invalid_argument);
  EXPECT_THROW(FlowLedger(g, std::vector<Amount> {1, 1}), std::invalid_argument);
  EXPECT_THROW(FlowLedger(nullptr), std::invalid_argument);
  FlowLedger fg(g, std::vector<Amount> {3});
  EXPECT_EQ(fg.residual(), (std::vector<Amount>{3}));
  EXPECT_EQ(fg.allocated(), (std::vector<Amount>{0}));
  EXPECT_EQ(std::vector<Amount>(fg.capacity_units().begin(), fg.capacity_units().end()), (std::vector<Amount>{5}));
  EXPECT_TRUE(fg.records().empty());
  EXPECT_EQ(fg.graph().get(), g.get());
  FlowLedger other(g);
  EXPECT_NE(fg.uid(), other.uid());
}

TEST(FlowLedger, BoundNoFlowsRow) {
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  PinnedOne p(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless, EdgePriority::Residual);
  // place(0): a no-op that creates nothing.
  auto r = p.place(0);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_FALSE(r.routes_changed);
  EXPECT_TRUE(fg.records().empty());
  EXPECT_EQ(LedgerAccess::next_flow_id(fg), 0);
  // reset and copy rejected while a policy is bound; move rejected too.
  EXPECT_NE(error_of([&] { fg.reset(); }).find("bound"), std::string::npos);
  EXPECT_THROW(FlowLedger copy(fg), std::invalid_argument);
  EXPECT_THROW(FlowLedger moved(std::move(fg)), std::invalid_argument);
  // remove() with no flows: no allocation to revert; the epoch state resets.
  p.remove();
  EXPECT_TRUE(fg.records().empty());
  EXPECT_FALSE(p.infinite());
  EXPECT_TRUE(p.requested().is_zero());
  // A positive request creates the cohort (new id), routes_changed true.
  r = p.place(2);
  EXPECT_EQ(r.placed, 2);
  EXPECT_TRUE(r.routes_changed);
  ASSERT_TRUE(p.flow_id().has_value());
  EXPECT_EQ(p.flow_id()->second, 0);
}

TEST(FlowLedger, NoRouteCreatesNothing) {
  // A bundle whose only walk is excluded by the Eligibility: usable_bundles() == 0.
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  auto el = std::make_shared<const Eligibility>(g->eligibility({}, std::vector<EdgeId> {0}));
  PinnedOne p(fg, 0, 1, all_edges(*g), FlowPlacement::MaxFlow, Admission::Lossless, EdgePriority::Residual, el);
  EXPECT_EQ(p.usable_bundles(), 0);
  auto r = p.place(3);
  EXPECT_EQ(r.termination, Termination::NoRoute);
  EXPECT_EQ(r.placed, 0);
  ASSERT_TRUE(r.remaining.has_value());
  EXPECT_EQ(*r.remaining, 3);
  EXPECT_FALSE(r.routes_changed);
  EXPECT_TRUE(fg.records().empty());
  EXPECT_EQ(LedgerAccess::next_flow_id(fg), 0);
  // A bundle without any walk before pruning is malformed.
  auto g2 = shared(make_graph(3, {{0, 1, 5, 1}, {2, 1, 5, 1}}));
  FlowLedger fg2(g2);
  EXPECT_NE(error_of([&] {
              PinnedOne bad(fg2, 0, 2, all_edges(*g2), FlowPlacement::MaxFlow, Admission::Lossless,
                             EdgePriority::Residual);
            }).find("malformed bundle"),
            std::string::npos);
}

TEST(FlowLedger, BoundWithFlowsRowAndEpochs) {
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  PinnedOne p(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless, EdgePriority::Residual);
  EXPECT_EQ(p.place(2).placed, 2);
  // Adds: T grows within the epoch.
  auto r = p.place(std::nullopt);
  EXPECT_EQ(r.placed, 3);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  EXPECT_FALSE(r.remaining.has_value());
  EXPECT_TRUE(p.infinite());
  // A finite request after +inf is rejected within the epoch, place(0) included (validation
  // precedes the zero shortcut); the state is unchanged.
  EXPECT_NE(error_of([&] { (void)p.place(1); }).find("after +inf"), std::string::npos);
  EXPECT_NE(error_of([&] { (void)p.place(0); }).find("after +inf"), std::string::npos);
  EXPECT_EQ(alloc_in(fg), (std::vector<Amount>{5}));
  EXPECT_NE(error_of([&] { fg.reset(); }).find("bound"), std::string::npos);
  EXPECT_THROW(FlowLedger copy(fg), std::invalid_argument);
  // remove(): allocations reverted, record gone, id retired.
  const auto old = *p.flow_id();
  p.remove();
  EXPECT_EQ(alloc_in(fg), (std::vector<Amount>{0}));
  EXPECT_TRUE(fg.records().empty());
  EXPECT_NE(error_of([&] { (void)fg.record(old); }).find("retired or unknown"), std::string::npos);
  // place(+inf); remove(); place(1) is accepted: a new epoch with a new id.
  r = p.place(1);
  EXPECT_EQ(r.placed, 1);
  EXPECT_TRUE(r.routes_changed);
  EXPECT_EQ(p.flow_id()->second, 1);
  p.remove();
  p.remove();  // on an empty cohort: resets the epoch
  EXPECT_TRUE(p.requested().is_zero());
  EXPECT_FALSE(p.infinite());
  EXPECT_EQ(LedgerAccess::next_flow_id(fg), 2);
}

TEST(FlowLedger, OrphanedRecordsResetAndCopy) {
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  auto p = std::make_unique<PinnedOne>(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless,
                                        EdgePriority::Residual);
  {
    PinnedOne q(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless, EdgePriority::Residual);
    EXPECT_EQ(q.place(2).placed, 2);
  }  // destroyed without remove(): the record stays as an orphan
  ASSERT_EQ(fg.records().size(), 1u);
  EXPECT_EQ(fg.records()[0].delivered, 2);
  EXPECT_EQ(LedgerAccess::orphan_owners(fg), (std::vector<std::int64_t>{1}));
  EXPECT_EQ(LedgerAccess::bound_owners(fg), (std::vector<std::int64_t>{0}));
  // Another policy sees the orphan's allocation (it is never changed again).
  EXPECT_EQ(p->place(5).placed, 3);
  EXPECT_NE(error_of([&] { fg.reset(); }).find("bound"), std::string::npos);
  p->remove();
  p.reset();  // destroyed: nothing bound
  // Copy of an unbound ledger: new uid, same ids under the new uid.
  FlowLedger c(fg);
  EXPECT_NE(c.uid(), fg.uid());
  ASSERT_EQ(c.records().size(), 1u);
  EXPECT_EQ(c.records()[0].flow_id, (FlowId {c.uid(), 0}));
  EXPECT_EQ(c.records()[0].delivered, 2);
  EXPECT_EQ(c.allocated(), fg.allocated());
  // reset: orphans and their allocations cleared, ids retired (never reused), new baseline.
  fg.reset(std::vector<Amount> {4});
  EXPECT_TRUE(fg.records().empty());
  EXPECT_EQ(fg.residual(), (std::vector<Amount>{4}));
  EXPECT_TRUE(LedgerAccess::orphan_owners(fg).empty());
  EXPECT_EQ(LedgerAccess::next_flow_id(fg), 2);
  PinnedOne r(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless, EdgePriority::Residual);
  EXPECT_EQ(r.place(9).placed, 4);
  EXPECT_EQ(r.flow_id()->second, 2);
  // The copy is independent of the original.
  EXPECT_EQ(c.records()[0].delivered, 2);
}

TEST(FlowLedger, MovedFromLedgerIsUnusable) {
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  FlowLedger moved(std::move(fg));
  EXPECT_THROW((void)fg.uid(), std::logic_error);
  EXPECT_THROW((void)fg.records(), std::logic_error);
  EXPECT_THROW((void)fg.residual(), std::logic_error);
  EXPECT_THROW(fg.reset(), std::logic_error);
  EXPECT_THROW((void)LedgerAccess::bind(fg), std::logic_error);
  EXPECT_EQ(moved.residual(), (std::vector<Amount>{5}));
}

TEST(FlowLedger, FlowIdsOfOtherLedgersAndUnknownIds) {
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger a(g), b(g);
  PinnedOne p(a, 0, 1, all_edges(*g), FlowPlacement::MaxFlow, Admission::Lossless, EdgePriority::Residual);
  (void)p.place(1);
  const auto id = *p.flow_id();
  EXPECT_EQ(a.record(id).delivered, 1);
  EXPECT_NE(error_of([&] { (void)b.record(id); }).find("another ledger"), std::string::npos);
  EXPECT_NE(error_of([&] { (void)a.record(FlowId {a.uid(), 7}); }).find("retired or unknown"), std::string::npos);
}

TEST(FlowLedger, RecordsCarryTheContractFields) {
  // Lossy EqualBalanced on the diamond 10/10/1/10 (prototype 1 lossy_diamond): 10 then 1.
  auto g = shared(make_graph(4, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 1, 1}, {2, 3, 10, 1}}));
  const auto inv = id_of_input(*g);
  FlowLedger fg(g);
  PinnedOne p(fg, 0, 3, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossy, EdgePriority::Residual);
  auto r = p.place(10);
  EXPECT_EQ(r.placed, 6);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_EQ(r.drops, (std::vector<std::pair<EdgeId, Amount>>{{inv[2], 4}}));
  r = p.place(1);
  EXPECT_EQ(r.placed, 0);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_EQ(*r.remaining, 1);
  EXPECT_EQ(alloc_in(fg), (std::vector<Amount>{6, 5, 1, 5}));
  const auto rec = fg.records().at(0);
  EXPECT_EQ(rec.flow_id, (FlowId {fg.uid(), 0}));
  EXPECT_EQ(rec.owner, 0);
  EXPECT_EQ(rec.sources, (std::vector<NodeId>{0}));
  EXPECT_EQ(rec.destinations, (std::vector<NodeId>{3}));
  EXPECT_EQ(rec.placement, FlowPlacement::EqualBalanced);
  EXPECT_EQ(rec.admission, Admission::Lossy);
  ASSERT_TRUE(rec.bundle);
  EXPECT_EQ(rec.bundle->via_edges().size(), 4u);
  EXPECT_TRUE(*rec.bundle == *all_edges(*g));
  EXPECT_EQ(rec.route_cost, 2);
  ASSERT_TRUE(rec.offered.has_value());
  EXPECT_EQ(rec.offered->lo, 11u);
  EXPECT_EQ(rec.delivered, 6);
  ASSERT_EQ(rec.drops.size(), 1u);
  EXPECT_EQ(rec.drops[0].first, inv[2]);
  EXPECT_EQ(rec.drops[0].second.lo, 5u);
  // +inf on a fresh ledger: offered is the injection (20), delivered 11, drops 9.
  FlowLedger f2(g);
  PinnedOne q(f2, 0, 3, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossy, EdgePriority::Residual);
  r = q.place(std::nullopt);
  EXPECT_EQ(r.placed, 11);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  EXPECT_EQ(f2.records()[0].offered->lo, 20u);
  // Lossless records carry no offered.
  FlowLedger f3(g);
  PinnedOne m(f3, 0, 3, all_edges(*g), FlowPlacement::MaxFlow, Admission::Lossless, EdgePriority::Residual);
  (void)m.place(3);
  EXPECT_FALSE(f3.records()[0].offered.has_value());
  EXPECT_TRUE(f3.records()[0].drops.empty());
}

TEST(FlowLedger, CountersBeyondSixtyFourBits) {
  // Prototype 1 lossy_oversized_and_saturated_offers: every unit of a 2^62 - 1 offer dropped on
  // a saturated link, many times: offered and drops grow beyond 2^64.
  auto g = shared(make_graph(2, {{0, 1, 1, 1}}));
  FlowLedger fg(g);
  PinnedOne p(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossy, EdgePriority::Residual);
  const Amount big = kDomainLimit - 1;
  for (int i = 0; i < 9; ++i) (void)p.place(big);
  const auto rec = fg.records().at(0);
  // offered = 9 (2^62 - 1) = 2^65 + 2^62 - 9
  EXPECT_EQ(rec.offered->hi, 2u);
  EXPECT_EQ(rec.offered->lo, (std::uint64_t {1} << 62) - 9);
  EXPECT_EQ(rec.delivered, 1);
  EXPECT_EQ(rec.drops[0].second.hi, 2u);
  EXPECT_EQ(rec.drops[0].second.lo, (std::uint64_t {1} << 62) - 10);
  EXPECT_EQ(p.requested().hi, 2u);
}

// The public counter is self-contained (flow_record.hpp): its own hi and lo, no base class from
// the private arithmetic, the same carry and the same failure at 2^127, which changes nothing.
static_assert(std::is_aggregate_v<Counter128> && std::is_trivially_copyable_v<Counter128>);
static_assert(std::is_standard_layout_v<Counter128> && sizeof(Counter128) == 2 * sizeof(std::uint64_t));
static_assert(!std::is_base_of_v<netgraph::core::detail::U128, Counter128> &&
              !std::is_convertible_v<Counter128, netgraph::core::detail::U128>);

TEST(FlowLedger, Counter128CarriesAndFailsAt2To127) {
  const std::uint64_t max64 = std::numeric_limits<std::uint64_t>::max();
  Counter128 c;
  EXPECT_TRUE(c.is_zero());
  EXPECT_EQ(c, Counter128 {});
  // A carry from lo into hi.
  c.lo = max64 - 2;
  EXPECT_TRUE(c.add(5));
  EXPECT_EQ(c.hi, 1u);
  EXPECT_EQ(c.lo, 2u);
  EXPECT_FALSE(c.is_zero());
  EXPECT_EQ(c, (Counter128 {1, 2}));
  EXPECT_NE(c, (Counter128 {0, 2}));
  EXPECT_NE(c, (Counter128 {1, 3}));
  // 2^127 - 1 is the largest value; adding 0 keeps it, any positive amount fails unchanged,
  // also through the carry.
  Counter128 top {(std::uint64_t {1} << 63) - 1, max64};
  EXPECT_TRUE(top.add(0));
  EXPECT_FALSE(top.add(1));
  EXPECT_EQ(top, (Counter128 {(std::uint64_t {1} << 63) - 1, max64}));
  Counter128 near {(std::uint64_t {1} << 63) - 1, max64 - 4};
  EXPECT_TRUE(near.add(4));
  EXPECT_EQ(near, top);
  Counter128 below {(std::uint64_t {1} << 63) - 1, max64 - 4};
  EXPECT_FALSE(below.add(kDomainLimit - 1));
  EXPECT_EQ(below, (Counter128 {(std::uint64_t {1} << 63) - 1, max64 - 4}));
}

TEST(FlowLedger, RequestDomainAndFailedCallsLeaveState) {
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  PinnedOne p(fg, 0, 1, all_edges(*g), FlowPlacement::Weighted, Admission::Lossless, EdgePriority::Capacity);
  (void)p.place(2);
  const auto before = fg.records();
  EXPECT_THROW((void)p.place(kDomainLimit), std::invalid_argument);
  EXPECT_THROW((void)p.place(-1), std::invalid_argument);
  EXPECT_EQ(fg.records().size(), before.size());
  EXPECT_EQ(fg.allocated(), (std::vector<Amount>{2}));
  // Constructor validation.
  auto other = make_graph(2, {{0, 1, 5, 1}});
  EXPECT_NE(error_of([&] {
              PinnedOne x(fg, 0, 1, bundle(other, {0}), FlowPlacement::MaxFlow, Admission::Lossless,
                           EdgePriority::Residual);
            }).find("foreign PredDAG"),
            std::string::npos);
  EXPECT_NE(error_of([&] {
              PinnedOne x(fg, 0, 0, all_edges(*g), FlowPlacement::MaxFlow, Admission::Lossless, EdgePriority::Residual);
            }).find("both endpoint sets"),
            std::string::npos);
  EXPECT_NE(error_of([&] {
              PinnedOne x(fg, 0, 1, all_edges(*g), FlowPlacement::MaxFlow, Admission::Lossy, EdgePriority::Residual);
            }).find("Lossy"),
            std::string::npos);
  EXPECT_NE(error_of([&] {
              PinnedOne x(fg, 0, 1, all_edges(*g), FlowPlacement::Weighted, Admission::Lossless, EdgePriority::Headroom);
            }).find("Headroom"),
            std::string::npos);
  // A rejected construction binds nothing.
  EXPECT_EQ(LedgerAccess::bound_owners(fg), (std::vector<std::int64_t>{0}));
}

TEST(FlowLedger, BottleneckVersusMinCut) {
  // Prototype 1 ledger_bottleneck_vs_min_cut: (0->1 cap 1 cost 1), (0->2 cap 2 cost 3),
  // (1->2 cap 1), (2->3 cap 2). After the MaxFlow +inf [1, 1, 1, 2] the bottleneck (forward
  // residual) is {0->1, 2->3}, not the cut {2->3}.
  auto g = shared(make_graph(4, {{0, 1, 1, 1}, {0, 2, 2, 3}, {1, 2, 1, 1}, {2, 3, 2, 1}}));
  const auto inv = id_of_input(*g);
  FlowLedger fg(g);
  PinnedOne p(fg, 0, 3, all_edges(*g), FlowPlacement::MaxFlow, Admission::Lossless, EdgePriority::Residual);
  EXPECT_EQ(p.place(std::nullopt).placed, 2);
  EXPECT_EQ(alloc_in(fg), (std::vector<Amount>{1, 1, 1, 2}));
  std::vector<EdgeId> want = {inv[0], inv[3]};
  std::sort(want.begin(), want.end());
  EXPECT_EQ(fg.bottleneck(std::vector<NodeId> {0}), want);
  // A zero-capacity edge at the frontier is not listed; an excluded edge is not crossed.
  auto g2 = shared(make_graph(3, {{0, 1, 0, 1}, {0, 2, 1, 1}, {2, 1, 1, 1}}));
  FlowLedger f2(g2);
  EXPECT_TRUE(f2.bottleneck(std::vector<NodeId> {0}).empty());
  const auto inv2 = id_of_input(*g2);
  auto el = g2->eligibility({}, std::vector<EdgeId> {inv2[2]});
  EXPECT_TRUE(f2.bottleneck(std::vector<NodeId> {0}, &el).empty());  // 2->1 excluded: not listed
  // An excluded source contributes nothing; a foreign Eligibility throws.
  const std::vector<NodeId> xn = {0};
  auto ex = g2->eligibility(xn, {});
  EXPECT_TRUE(f2.bottleneck(std::vector<NodeId> {0}, &ex).empty());
  auto other = make_graph(3, {{0, 1, 1, 1}});
  auto foreign = other.eligibility({}, std::vector<EdgeId> {0});
  EXPECT_THROW((void)f2.bottleneck(std::vector<NodeId> {0}, &foreign), std::invalid_argument);
  EXPECT_THROW((void)f2.bottleneck(std::vector<NodeId> {}), std::invalid_argument);
}

TEST(FlowLedger, TwoOwnersNeverCancelEachOther) {
  // Prototype 1 maxflow_two_owner_no_cross_cancellation: A holds the chain 0-1-2-3; B's bundle
  // (chain plus bypasses) could reach 2 only by cancelling A's middle edge: B places 0 (record
  // at zero delivery, exhausted) and A is untouched; after A's remove() B places 2.
  auto g = shared(make_graph(4, {{0, 1, 1, 1}, {1, 2, 1, 1}, {2, 3, 1, 1}, {0, 2, 1, 3}, {1, 3, 1, 3}}));
  const auto inv = id_of_input(*g);
  FlowLedger fg(g);
  PinnedOne a(fg, 0, 3, bundle(*g, {inv[0], inv[1], inv[2]}), FlowPlacement::MaxFlow, Admission::Lossless,
               EdgePriority::Residual);
  EXPECT_EQ(a.place(1).placed, 1);
  PinnedOne b(fg, 0, 3, all_edges(*g), FlowPlacement::MaxFlow, Admission::Lossless, EdgePriority::Residual);
  auto r = b.place(2);
  EXPECT_EQ(r.placed, 0);
  EXPECT_TRUE(r.routes_changed);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  EXPECT_EQ(fg.record(*b.flow_id()).delivered, 0);
  EXPECT_EQ(alloc_in(fg), (std::vector<Amount>{1, 1, 1, 0, 0}));
  a.remove();
  r = b.place(2);
  EXPECT_EQ(r.placed, 2);
  EXPECT_EQ(alloc_in(fg), (std::vector<Amount>{1, 0, 1, 1, 1}));
}

TEST(FlowLedger, AggregatesEqualTheSumOfRecords) {
  // Random sequences of placements, removals and destructions on one ledger: allocated() is the
  // exact sum of the records' edges and residual() = baseline - allocated() within [0, baseline].
  std::mt19937 rng(5);
  for (int trial = 0; trial < 60; ++trial) {
    const int n = 4 + static_cast<int>(rng() % 3);
    std::vector<InEdge> edges;
    for (int i = 0; i < 3 * n; ++i) {
      const NodeId u = static_cast<NodeId>(rng() % (n - 1));
      const NodeId v = static_cast<NodeId>(u + 1 + rng() % static_cast<std::uint32_t>(n - 1 - u));
      edges.push_back({u, v, static_cast<Amount>(rng() % 7), 1 + static_cast<Cost>(rng() % 3)});
    }
    auto g = shared(make_graph(n, edges));
    FlowLedger fg(g);
    std::vector<std::unique_ptr<PinnedOne>> pols;
    const FlowPlacement kinds[3] = {FlowPlacement::MaxFlow, FlowPlacement::Weighted, FlowPlacement::EqualBalanced};
    for (int k = 0; k < 4; ++k) {
      const auto kind = kinds[rng() % 3];
      const auto adm = kind != FlowPlacement::MaxFlow && rng() % 3 == 0 ? Admission::Lossy : Admission::Lossless;
      try {
        pols.push_back(std::make_unique<PinnedOne>(fg, 0, n - 1, all_edges(*g), kind, adm, EdgePriority::Capacity));
      } catch (const std::invalid_argument&) {
        break;  // no walk at all
      }
    }
    for (int step = 0; step < 12 && !pols.empty(); ++step) {
      auto& p = pols[rng() % pols.size()];
      if (!p) continue;
      const int op = static_cast<int>(rng() % 6);
      if (op == 0) {
        p->remove();
      } else if (op == 1) {
        p.reset();
      } else if (!p->infinite()) {
        (void)p->place(op == 2 ? std::optional<Amount>() : std::optional<Amount>(rng() % 9));
      }
      std::vector<Amount> sum(static_cast<std::size_t>(g->num_edges()), 0);
      for (const auto& rec : fg.records()) {
        for (const auto& [e, x] : rec.edges) sum[static_cast<std::size_t>(e)] += x;
      }
      ASSERT_EQ(fg.allocated(), sum);
      const auto res = fg.residual();
      for (std::size_t e = 0; e < res.size(); ++e) {
        ASSERT_GE(res[e], 0);
        ASSERT_EQ(res[e] + sum[e], g->capacity_units()[e]);
      }
    }
  }
}

TEST(FlowLedger, NoRouteCallsRecordTheEpochState) {
  // A positive request with no usable bundle creates nothing, but the epoch records it: after
  // +inf a finite request (place(0) included) is rejected until remove() (APPROACH 8, 9.2;
  // prototype 2's reading; prototype 1 kept no state on such a call).
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  auto el = std::make_shared<const Eligibility>(g->eligibility({}, std::vector<EdgeId> {0}));
  PinnedOne p(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless, EdgePriority::Residual,
               el);
  EXPECT_EQ(p.place(3).termination, Termination::NoRoute);
  EXPECT_EQ(p.requested().lo, 3u);
  auto r = p.place(std::nullopt);
  EXPECT_EQ(r.termination, Termination::NoRoute);
  EXPECT_FALSE(r.remaining.has_value());
  EXPECT_TRUE(p.infinite());
  EXPECT_NE(error_of([&] { (void)p.place(1); }).find("after +inf"), std::string::npos);
  EXPECT_NE(error_of([&] { (void)p.place(0); }).find("after +inf"), std::string::npos);
  p.remove();
  EXPECT_FALSE(p.infinite());
  EXPECT_EQ(p.place(1).termination, Termination::NoRoute);
  EXPECT_TRUE(fg.records().empty());
  EXPECT_EQ(LedgerAccess::next_flow_id(fg), 0);
}

TEST(FlowLedger, RouteCostIsTheBundlesMinimumWalk) {
  // Weighted Residual on a route whose only edge another owner fills: the split graph is
  // empty (record at zero delivery), route_cost is still the bundle's walk cost 3.
  auto g = shared(make_graph(2, {{0, 1, 2, 3}}));
  FlowLedger fg(g);
  PinnedOne holder(fg, 0, 1, all_edges(*g), FlowPlacement::MaxFlow, Admission::Lossless, EdgePriority::Residual);
  EXPECT_EQ(holder.place(2).placed, 2);
  PinnedOne w(fg, 0, 1, all_edges(*g), FlowPlacement::Weighted, Admission::Lossless, EdgePriority::Residual);
  auto r = w.place(1);
  EXPECT_EQ(r.placed, 0);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  EXPECT_TRUE(r.routes_changed);
  ASSERT_NE(w.plan(), nullptr);
  EXPECT_TRUE(w.plan()->members.empty());
  const auto rec = fg.record(*w.flow_id());
  EXPECT_EQ(rec.route_cost, 3);
  EXPECT_EQ(rec.delivered, 0);
  // A mixed-cost bundle: the cheapest walk (1 + 1 + 1) is the route cost, not the cost-4 edge.
  auto g2 = shared(make_graph(4, {{0, 1, 1, 1}, {1, 2, 1, 1}, {2, 3, 1, 1}, {0, 3, 1, 4}}));
  FlowLedger f2(g2);
  PinnedOne m(f2, 0, 3, all_edges(*g2), FlowPlacement::MaxFlow, Admission::Lossless, EdgePriority::Residual);
  (void)m.place(1);
  EXPECT_EQ(f2.records()[0].route_cost, 3);
}

TEST(FlowLedger, RemovalRestoresTheBaselineAndLeavesOtherOwnersIntact) {
  // review-2 R7: removing a flow restores the baseline residual (5 on a capacity-10 link), not
  // the capacity; review-check check_review6: another owner's allocation stays intact.
  auto g = shared(make_graph(2, {{0, 1, 10, 1}}));
  FlowLedger fg(g, std::vector<Amount> {5});
  PinnedOne a(fg, 0, 1, all_edges(*g), FlowPlacement::MaxFlow, Admission::Lossless, EdgePriority::Residual);
  PinnedOne b(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless, EdgePriority::Residual);
  EXPECT_EQ(a.place(3).placed, 3);
  EXPECT_EQ(b.place(1).placed, 1);
  EXPECT_EQ(fg.residual(), (std::vector<Amount>{1}));
  a.remove();
  EXPECT_EQ(fg.residual(), (std::vector<Amount>{4}));
  EXPECT_EQ(fg.record(*b.flow_id()).edges, (std::vector<std::pair<EdgeId, Amount>>{{0, 1}}));
  b.remove();
  EXPECT_EQ(fg.residual(), (std::vector<Amount>{5}));
  // A baseline one unit above the quantised capacity is rejected (no clamp at this boundary).
  EXPECT_THROW(FlowLedger(g, std::vector<Amount> {11}), std::invalid_argument);
}

TEST(FlowLedger, LossyOffersArePerCall) {
  // review-13 lossy / rule 19: on the diamond, two calls of 1 offer each unit to the lowest
  // EdgeId first (both on the first branch) while one call of 2 splits [1, 1]: not additive.
  auto g = shared(make_graph(4, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 1, 1}, {2, 3, 10, 1}}));
  FlowLedger f1(g), f2(g);
  PinnedOne p1(f1, 0, 3, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossy, EdgePriority::Residual);
  PinnedOne p2(f2, 0, 3, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossy, EdgePriority::Residual);
  auto r = p1.place(1);
  EXPECT_EQ(r.placed, 1);
  r = p1.place(1);
  EXPECT_EQ(r.placed, 0);  // the tail 1->3 is full: the unit is dropped there
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_EQ(alloc_in(f1), (std::vector<Amount>{2, 0, 1, 0}));
  r = p2.place(2);
  EXPECT_EQ(r.placed, 2);
  EXPECT_EQ(alloc_in(f2), (std::vector<Amount>{1, 1, 1, 1}));
}

// review-1b-fixes N2: a commit carrying several records of one owner touches an edge from more
// than one record; its delta may pass through zero while it accumulates. Each edge's delta is
// published once.
namespace {

detail::Draft draft_on(const std::shared_ptr<const PredDAG>& dag, std::vector<std::pair<EdgeId, Amount>> edges,
                       std::optional<std::int64_t> id = std::nullopt) {
  detail::Draft d;
  d.id = id;
  FlowRecord& r = d.record;
  r.sources = {0};
  r.destinations = {1};
  r.bundle = dag;
  r.route_cost = 1;
  for (const auto& [e, x] : edges) r.delivered += x;
  r.edges = std::move(edges);
  return d;
}

Amount record_sum(const FlowLedger& fg, EdgeId e) {
  Amount sum = 0;
  for (const auto& r : fg.records()) {
    for (const auto& [f, x] : r.edges) {
      if (f == e) sum += x;
    }
  }
  return sum;
}

} // namespace

TEST(FlowLedger, MultiRecordCommitsPublishEachEdgeOnce) {
  auto g = shared(make_graph(2, {{0, 1, 3, 1}}));
  const auto dag = all_edges(*g);
  FlowLedger fg(g);
  const auto owner = LedgerAccess::bind(fg);
  std::vector<detail::Draft> first {draft_on(dag, {{0, 1}})};
  const auto ids = LedgerAccess::commit(fg, owner, first);
  // The record kept unchanged (-1, then +1: the delta passes through 0) plus a new 2-unit record.
  std::vector<detail::Draft> next {draft_on(dag, {{0, 1}}, ids[0].second), draft_on(dag, {{0, 2}})};
  EXPECT_EQ(LedgerAccess::commit(fg, owner, next).size(), 1U);
  EXPECT_EQ(fg.allocated(), (std::vector<Amount>{3}));
  EXPECT_EQ(fg.residual(), (std::vector<Amount>{0}));
  EXPECT_EQ(record_sum(fg, 0), 3);
  // Over capacity through the same pattern: rejected, nothing changed.
  const auto held = fg.records();
  const auto live = LedgerAccess::bound_owners(fg);
  std::vector<detail::Draft> over {draft_on(dag, {{0, 1}}, ids[0].second), draft_on(dag, {{0, 2}}, ids[0].second + 1),
                                   draft_on(dag, {{0, 1}})};
  EXPECT_THROW((void)LedgerAccess::commit(fg, owner, over), std::logic_error);
  EXPECT_EQ(fg.allocated(), (std::vector<Amount>{3}));
  EXPECT_EQ(fg.records().size(), held.size());
  EXPECT_EQ(LedgerAccess::next_flow_id(fg), 2);
  // A kept record gives its unit up while a new record takes the edge: published once.
  std::vector<detail::Draft> moved {draft_on(dag, {}, ids[0].second), draft_on(dag, {{0, 2}}, ids[0].second + 1),
                                    draft_on(dag, {{0, 1}})};
  EXPECT_EQ(LedgerAccess::commit(fg, owner, moved).size(), 1U);
  EXPECT_EQ(fg.allocated(), (std::vector<Amount>{3}));
  EXPECT_EQ(record_sum(fg, 0), 3);
  EXPECT_EQ(LedgerAccess::detach(fg, owner).size(), 3U);
  EXPECT_EQ(fg.allocated(), (std::vector<Amount>{0}));
  LedgerAccess::release(fg, owner);
}

TEST(FlowLedger, RandomMultiOwnerCommitsMatchAModel) {
  // Three owners on four edges; every commit carries all of its owner's records (each kept,
  // changed, emptied or moved) plus up to two new ones, so deltas cancel and reappear on shared
  // edges; detaches between commits. The ledger must equal an independent model after every call,
  // and a rejected commit must change nothing.
  std::mt19937 rng(11);
  auto g = shared(make_graph(4, {{0, 1, 6, 1}, {1, 2, 4, 1}, {2, 3, 5, 1}, {0, 3, 3, 1}}));
  const auto dag = all_edges(*g);
  const auto m = static_cast<std::size_t>(g->num_edges());
  FlowLedger fg(g);
  std::vector<std::int64_t> owners;
  for (int i = 0; i < 3; ++i) owners.push_back(LedgerAccess::bind(fg));
  // The model: per owner, its live records (id -> per-edge amounts).
  std::vector<std::map<std::int64_t, std::vector<Amount>>> model(owners.size());
  std::int64_t next_id = 0;
  int accepted = 0;
  int rejected = 0;
  auto random_edges = [&] {
    std::vector<Amount> a(m, 0);
    for (std::size_t e = 0; e < m; ++e) a[e] = (rng() % 3) == 0 ? static_cast<Amount>(rng() % 3) : 0;
    return a;
  };
  auto as_pairs = [](const std::vector<Amount>& a) {
    std::vector<std::pair<EdgeId, Amount>> out;
    for (std::size_t e = 0; e < a.size(); ++e) {
      if (a[e] != 0) out.emplace_back(static_cast<EdgeId>(e), a[e]);
    }
    return out;
  };
  for (int step = 0; step < 3000; ++step) {
    const std::size_t o = rng() % owners.size();
    if ((rng() % 10) == 0) {
      (void)LedgerAccess::detach(fg, owners[o]);
      model[o].clear();
    } else {
      std::vector<detail::Draft> drafts;
      std::map<std::int64_t, std::vector<Amount>> after;
      for (const auto& [id, a] : model[o]) {
        const auto kind = rng() % 4;  // kept, changed, emptied, moved to random edges
        std::vector<Amount> b = a;
        if (kind == 1) {
          b[rng() % m] += 1;
        } else if (kind == 2) {
          b.assign(m, 0);
        } else if (kind == 3) {
          b = random_edges();
        }
        drafts.push_back(draft_on(dag, as_pairs(b), id));
        after[id] = b;
      }
      const int fresh = static_cast<int>(rng() % 3);
      for (int k = 0; k < fresh; ++k) {
        const auto b = random_edges();
        drafts.push_back(draft_on(dag, as_pairs(b)));
        after[next_id + k] = b;
      }
      // Feasible exactly when every edge's new total stays within its capacity.
      std::vector<Amount> total(m, 0);
      for (std::size_t q = 0; q < owners.size(); ++q) {
        for (const auto& [id, a] : q == o ? after : model[q]) {
          for (std::size_t e = 0; e < m; ++e) total[e] += a[e];
        }
      }
      bool feasible = true;
      for (std::size_t e = 0; e < m; ++e) feasible = feasible && total[e] <= g->capacity_units()[e];
      const auto before = fg.allocated();
      const auto before_records = fg.records().size();
      if (feasible) {
        EXPECT_EQ(LedgerAccess::commit(fg, owners[o], drafts).size(), static_cast<std::size_t>(fresh));
        model[o] = after;
        next_id += fresh;
        ++accepted;
      } else {
        EXPECT_THROW((void)LedgerAccess::commit(fg, owners[o], drafts), std::logic_error);
        EXPECT_EQ(fg.allocated(), before);
        EXPECT_EQ(fg.records().size(), before_records);
        ++rejected;
      }
    }
    std::vector<Amount> want(m, 0);
    for (const auto& recs : model) {
      for (const auto& [id, a] : recs) {
        for (std::size_t e = 0; e < m; ++e) want[e] += a[e];
      }
    }
    ASSERT_EQ(fg.allocated(), want) << "step " << step;
    for (std::size_t e = 0; e < m; ++e) ASSERT_EQ(record_sum(fg, static_cast<EdgeId>(e)), want[e]);
    ASSERT_EQ(LedgerAccess::next_flow_id(fg), next_id);
  }
  EXPECT_GT(accepted, 1000);
  EXPECT_GT(rejected, 100);
  for (auto owner : owners) LedgerAccess::release(fg, owner);
}
