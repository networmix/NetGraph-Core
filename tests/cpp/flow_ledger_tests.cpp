// The integer ledger FlowLedger: the baseline, snapshots, FlowIds, bottleneck, and
// every row and transition of the lifetime table, driven through a policy on one route (Global
// with the route as its support, or Pinned over it as a supplied split route: test::RoutePolicy,
// route_policy.hpp); the 128-bit counters of the records;
// and the commit rule itself (each edge published once, the multi-owner model, the fresh-record
// path), driven through detail::LedgerAccess directly.
//
// | State             | place                         | remove()           | destroyed  | reset    | copy     |
// | bound, no flows   | place(0) no-op; positive      | resets the epoch   | none       | rejected | rejected |
// |                   | creates (no route: nothing)   |                    |            |          |          |
// | bound, flows      | adds or replaces              | reverts, retires   | orphaned   | rejected | rejected |
// | orphaned records  | n/a                           | n/a                | n/a        | clears   | copied   |
#include <gtest/gtest.h>

#include <algorithm>
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
#include <typeinfo>
#include <vector>

#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "route_policy.hpp"
#include "state_comparators.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using detail::LedgerAccess;
using netgraph::core::Termination;

namespace {

std::vector<Amount> alloc_in(const FlowLedger& fg) { return to_input(*fg.graph(), fg.allocated()); }

// What fn raises: the message when the exception's dynamic type is exactly std::logic_error (the
// ledger's refusals of a moved-from ledger and of a malformed commit), or a description of anything
// else (std::invalid_argument and the other standard errors derived from std::logic_error are other
// checks' refusals and must not pass for these).
std::string logic_error_of(const std::function<void()>& fn) {
  try {
    fn();
  } catch (const std::exception& e) {
    if (typeid(e) == typeid(std::logic_error)) return e.what();
    return std::string("another exception type (") + typeid(e).name() + "): " + e.what();
  } catch (...) {
    return "an exception not derived from std::exception";
  }
  return "no exception";
}

} // namespace

TEST(FlowLedger, BaselineValidationAndSnapshots) {
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  EXPECT_THROW(FlowLedger(g, std::vector<Amount> {6}), std::invalid_argument);
  EXPECT_THROW(FlowLedger(g, std::vector<Amount> {-1}), std::invalid_argument);
  EXPECT_THROW(FlowLedger(g, std::vector<Amount> {1, 1}), std::invalid_argument);
  EXPECT_EQ(error_of([] { FlowLedger fg(nullptr); }), "FlowLedger: graph is null");
  EXPECT_EQ(error_of([] { FlowLedger fg(nullptr, std::vector<Amount> {}); }), "FlowLedger: graph is null");
  FlowLedger fg(g, std::vector<Amount> {3});
  EXPECT_EQ(fg.residual(), (std::vector<Amount>{3}));
  EXPECT_EQ(fg.allocated(), (std::vector<Amount>{0}));
  EXPECT_EQ(std::vector<Amount>(fg.capacity_units().begin(), fg.capacity_units().end()), (std::vector<Amount>{5}));
  EXPECT_TRUE(fg.records().empty());
  EXPECT_EQ(fg.graph().get(), g.get());
  FlowLedger other(g);
  EXPECT_NE(fg.uid(), other.uid());
  EXPECT_EQ(other.residual(), (std::vector<Amount>{5}));  // no baseline: the capacities
  // Either constructor shares the ownership of the graph: it lives while the ledger does.
  for (const bool with_baseline : {false, true}) {
    auto owned = shared(make_graph(2, {{0, 1, 5, 1}}));
    const std::weak_ptr<const Graph> watch = owned;
    const FlowLedger held = with_baseline ? FlowLedger(owned, std::vector<Amount> {3}) : FlowLedger(owned);
    owned.reset();
    EXPECT_FALSE(watch.expired()) << with_baseline;
    EXPECT_EQ(held.graph()->num_edges(), 1) << with_baseline;
  }
}

// A baseline with one or several defects (one error, from the shared edge-amount check): its length
// before any entry, then its first entry outside [0, capacity] in EdgeId order, at construction
// and at reset; a bound policy is reported before the baseline, and a rejected reset changes
// nothing (records, allocation, baseline, the next FlowId).
TEST(FlowLedger, BaselineDefectsReportTheFirstInOrder) {
  const auto g = shared(edge_amount_graph());
  for (const auto& c : edge_amount_defects()) {
    EXPECT_EQ(error_of([&] { FlowLedger fg(g, c.amounts); }), "residual" + c.error);
  }
  FlowLedger fg(g);
  {
    RoutePolicy q(fg, 0, 3, split_route(fg, 0, 3, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless));
    ASSERT_EQ(q.place(1).placed, 1);
    for (const auto& c : edge_amount_defects()) {
      EXPECT_EQ(error_of([&] { fg.reset(c.amounts); }), "FlowLedger.reset is rejected while a policy is bound (1 bound)");
    }
  }  // destroyed without remove(): the record stays as an orphan
  const auto records = fg.records();
  const auto allocated = fg.allocated();
  const auto residual = fg.residual();
  ASSERT_EQ(records.size(), 1U);
  for (const auto& c : edge_amount_defects()) {
    EXPECT_EQ(error_of([&] { fg.reset(c.amounts); }), "residual" + c.error);
    ASSERT_EQ(fg.records().size(), 1U);
    EXPECT_EQ(fg.records()[0].flow_id, records[0].flow_id);
    EXPECT_EQ(fg.allocated(), allocated);
    EXPECT_EQ(fg.residual(), residual);
    EXPECT_EQ(LedgerAccess::next_flow_id(fg), 1);
  }
  fg.reset(edge_amounts_at_the_bounds());
  EXPECT_TRUE(fg.records().empty());
  EXPECT_EQ(fg.residual(), edge_amounts_at_the_bounds());
  fg.reset();  // no baseline: back to the capacities
  EXPECT_EQ(fg.residual(), vec(g->capacity_units()));
}

TEST(FlowLedger, BoundNoFlowsRow) {
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  RoutePolicy p(fg, 0, 1, split_route(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless));
  // place(0): a no-op that creates nothing.
  auto r = p.place(0);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_FALSE(r.routes_changed);
  EXPECT_TRUE(fg.records().empty());
  EXPECT_EQ(LedgerAccess::next_flow_id(fg), 0);
  // reset and copy rejected while a policy is bound; move rejected too.
  EXPECT_EQ(error_of([&] { fg.reset(); }), "FlowLedger.reset is rejected while a policy is bound (1 bound)");
  EXPECT_EQ(error_of([&] { FlowLedger copy(fg); }), "FlowLedger: a ledger with bound policies cannot be copied");
  EXPECT_EQ(error_of([&] { FlowLedger moved(std::move(fg)); }),
            "FlowLedger: a ledger with bound policies cannot be moved");
  // remove() with no flows and no request: nothing to revert, the epoch stays empty (the reset of a
  // recorded request without flows: NoRouteCallsRecordTheEpochState).
  p.remove();
  EXPECT_TRUE(fg.records().empty());
  EXPECT_FALSE(p.infinite());
  EXPECT_TRUE(p.requested().is_zero());
  // A positive request creates the flow (new id), routes_changed true.
  r = p.place(2);
  EXPECT_EQ(r.placed, 2);
  EXPECT_TRUE(r.routes_changed);
  ASSERT_TRUE(p.flow_id().has_value());
  EXPECT_EQ(p.flow_id()->second, 0);
}

TEST(FlowLedger, NoRouteCreatesNothing) {
  // A support whose only walk is excluded by the Eligibility: the roster is empty.
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  auto el = std::make_shared<const Eligibility>(g->eligibility({}, std::vector<EdgeId> {0}));
  RoutePolicy p(fg, 0, 1, on_support(all_edges(*g)), el);
  EXPECT_TRUE(detail::PolicyAccess::engine(p.policy()).roster.empty());
  auto r = p.place(3);
  EXPECT_EQ(r.termination, Termination::NoRoute);
  EXPECT_EQ(r.placed, 0);
  ASSERT_TRUE(r.remaining.has_value());
  EXPECT_EQ(*r.remaining, 3);
  EXPECT_FALSE(r.routes_changed);
  EXPECT_TRUE(fg.records().empty());
  EXPECT_EQ(LedgerAccess::next_flow_id(fg), 0);
  // A route without any walk: a supplied split route is rejected; as a support it leaves the
  // roster empty.
  auto g2 = shared(make_graph(3, {{0, 1, 5, 1}, {2, 1, 5, 1}}));
  FlowLedger fg2(g2);
  EXPECT_EQ(error_of([&] {
              RoutePolicy bad(fg2, 0, 2, split_route(fg2, 0, 2, all_edges(*g2), FlowPlacement::EqualBalanced, Admission::Lossless));
            }),
            "routes[0]: no walk from head 0 to a destination");
  RoutePolicy none(fg2, 0, 2, on_support(all_edges(*g2)));
  EXPECT_TRUE(detail::PolicyAccess::engine(none.policy()).roster.empty());
  EXPECT_EQ(none.place(1).termination, Termination::NoRoute);
}

TEST(FlowLedger, BoundWithFlowsRowAndEpochs) {
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  RoutePolicy p(fg, 0, 1, split_route(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless));
  EXPECT_EQ(p.place(2).placed, 2);
  // Adds: the requested total grows within the epoch.
  auto r = p.place(std::nullopt);
  EXPECT_EQ(r.placed, 3);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  EXPECT_FALSE(r.remaining.has_value());
  EXPECT_TRUE(p.infinite());
  // A finite request after +inf is rejected within the epoch, place(0) included (validation
  // precedes the zero shortcut); the requested total, the +inf flag, the flow and the
  // allocations are unchanged.
  const Counter128 requested = p.requested();
  const auto id = p.flow_id();
  EXPECT_NE(error_of([&] { (void)p.place(1); }).find("after infinite placement"), std::string::npos);
  EXPECT_NE(error_of([&] { (void)p.place(0); }).find("after infinite placement"), std::string::npos);
  EXPECT_EQ(p.requested(), requested);
  EXPECT_TRUE(p.infinite());
  EXPECT_EQ(p.flow_id(), id);
  EXPECT_EQ(alloc_in(fg), (std::vector<Amount>{5}));
  // A rejected reset, copy or move of a ledger holding a record leaves the ledger as it was (the
  // record, the allocations and the binding): remove() below still reverts the flow.
  const LedgerState held = ledger_state(fg);
  EXPECT_EQ(error_of([&] { fg.reset(); }), "FlowLedger.reset is rejected while a policy is bound (1 bound)");
  EXPECT_EQ(error_of([&] { FlowLedger copy(fg); }), "FlowLedger: a ledger with bound policies cannot be copied");
  EXPECT_EQ(error_of([&] { FlowLedger moved(std::move(fg)); }),
            "FlowLedger: a ledger with bound policies cannot be moved");
  expect_same(ledger_state(fg), held, "the rejected reset, copy and move");
  EXPECT_EQ(fg.records().at(0).flow_id, *p.flow_id());
  EXPECT_EQ(alloc_in(fg), (std::vector<Amount>{5}));
  // remove(): allocations reverted, record gone, id retired.
  const auto old = *p.flow_id();
  p.remove();
  EXPECT_EQ(alloc_in(fg), (std::vector<Amount>{0}));
  EXPECT_TRUE(fg.records().empty());
  EXPECT_EQ(error_of([&] { (void)fg.record(old); }),
            "FlowId (" + std::to_string(old.first) + ", " + std::to_string(old.second) + ") is retired or unknown");
  // place(+inf); remove(); place(1) is accepted: a new epoch with a new id.
  r = p.place(1);
  EXPECT_EQ(r.placed, 1);
  EXPECT_TRUE(r.routes_changed);
  EXPECT_EQ(p.flow_id()->second, 1);
  p.remove();  // the epoch resets: the request of 1 is gone
  EXPECT_TRUE(p.requested().is_zero());
  EXPECT_FALSE(p.infinite());
  p.remove();  // again, with no flow: nothing changes
  EXPECT_TRUE(fg.records().empty());
  EXPECT_TRUE(p.requested().is_zero());
  EXPECT_FALSE(p.infinite());
  EXPECT_EQ(LedgerAccess::next_flow_id(fg), 2);
}

TEST(FlowLedger, OrphanedRecordsResetAndCopy) {
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  auto p = std::make_unique<RoutePolicy>(fg, 0, 1, split_route(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless));
  {
    RoutePolicy q(fg, 0, 1, split_route(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless));
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
  RoutePolicy r(fg, 0, 1, split_route(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless));
  EXPECT_EQ(r.place(9).placed, 4);
  EXPECT_EQ(r.flow_id()->second, 2);
  // The copy is independent of the original.
  EXPECT_EQ(c.records()[0].delivered, 2);
}

// The copy and the move constructors carry everything the source holds: the baseline, the
// aggregate, every record field, the next FlowId and the owners (the orphaned ones and the id the
// next policy gets). The source holds values a new or reset ledger would not have, so a member
// left behind shows: a baseline below the capacities, a retired FlowId (the next FlowId is past
// every live one), a Lossy record with an offer and drops beside a Lossless one, and an orphaned
// owner without records. The copy has a new uid and leaves its source as it was; the move keeps
// the uid and leaves its source unusable.
TEST(FlowLedger, CopyAndMoveCarryTheWholeLedger) {
  auto g = shared(make_graph(4, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 1, 1}, {2, 3, 10, 1}}));
  const auto inv = id_of_input(*g);
  std::vector<Amount> baseline(4);
  const Amount by_input[4] = {9, 8, 1, 9};
  for (std::size_t i = 0; i < 4; ++i) baseline[static_cast<std::size_t>(inv[i])] = by_input[i];
  // A reservation baseline below the capacities: 7, 7, 1, 7 by input.
  std::vector<Amount> reservation(4);
  const Amount reserve_by_input[4] = {7, 7, 1, 7};
  for (std::size_t i = 0; i < 4; ++i) reservation[static_cast<std::size_t>(inv[i])] = reserve_by_input[i];
  FlowLedger fg(g, baseline, reservation);
  {
    RoutePolicy lossy(fg, 0, 3, split_route(fg, 0, 3, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossy));
    (void)lossy.place(10);  // owner 0, FlowId 0: delivers 6, drops 4 on the 1 -> 3 edge
    RoutePolicy removed(fg, 0, 3, on_support(all_edges(*g)));
    (void)removed.place(1);  // owner 1, FlowId 1, retired by remove()
    removed.remove();
    RoutePolicy lossless(fg, 0, 3, on_support(all_edges(*g)));
    (void)lossless.place(2);  // owner 2, FlowId 2
    // owner 3, FlowId 3: one T1 tunnel over 0 -> 2 -> 3 holding 5 (the request filters out the
    // 1-unit link 1 -> 3), carrying nothing.
    FlowPolicy holds(fg, v({0}), v({3}), FlowPolicyConfig {.routing = CapacityAware {Tunnels {1}}});
    (void)holds.replace(0, 5);
  }  // four orphaned owners, three records: a copy or a move needs a ledger without bound policies
  const LedgerState source = ledger_state(fg);
  const std::vector<Amount> source_baseline = vec(LedgerAccess::baseline(fg));
  const std::vector<Amount> source_reservation = vec(LedgerAccess::reservation_baseline(fg));
  // The source is what the comment above says, so that each part below can differ.
  ASSERT_EQ(source.records.size(), 3u);
  EXPECT_EQ(source.records[2].reserved_edges, (std::vector<std::pair<EdgeId, Amount>> {{inv[1], 5}, {inv[3], 5}}));
  EXPECT_EQ(source.reserved[static_cast<std::size_t>(inv[1])], 5);
  EXPECT_EQ(source_reservation, reservation);
  EXPECT_TRUE(source.records[0].offered.has_value());
  EXPECT_EQ(source.records[0].delivered, 6);
  EXPECT_EQ(source.records[0].drops, (std::vector<std::pair<EdgeId, Counter128>>{{inv[2], Counter128 {0, 4}}}));
  EXPECT_FALSE(source.records[1].offered.has_value());
  EXPECT_EQ(source.records[1].flow_id.second, 2);
  EXPECT_EQ(source.next_id, 4);
  EXPECT_TRUE(source.bound.empty());
  EXPECT_EQ(source.orphans, (std::vector<std::int64_t>{0, 1, 2, 3}));
  EXPECT_NE(source_baseline, vec(g->capacity_units()));
  EXPECT_NE(source.allocated, std::vector<Amount>(4, 0));

  // Everything the source held (both baselines and both projections), then the next policy: it
  // gets the next owner id (4) and the next FlowId.
  auto expect_whole = [&](FlowLedger& x, const std::string& what) {
    expect_same(ledger_state(x), source, what);
    EXPECT_EQ(vec(LedgerAccess::baseline(x)), source_baseline) << what;
    EXPECT_EQ(vec(LedgerAccess::reservation_baseline(x)), source_reservation) << what;
    RoutePolicy next(x, 0, 3, on_support(all_edges(*g)));
    EXPECT_EQ(next.place(1).placed, 1) << what;
    ASSERT_TRUE(next.flow_id().has_value()) << what;
    EXPECT_EQ(next.flow_id()->second, source.next_id) << what;
    EXPECT_EQ(x.record(*next.flow_id()).owner, 4) << what;
    EXPECT_EQ(LedgerAccess::bound_owners(x), (std::vector<std::int64_t>{4})) << what;
    next.remove();
  };

  FlowLedger copy(fg);
  EXPECT_NE(copy.uid(), fg.uid());
  expect_same(ledger_state(fg), source, "the source after the copy");
  expect_whole(copy, "the copy");
  expect_same(ledger_state(fg), source, "the source after a placement on the copy");

  const auto uid = fg.uid();
  FlowLedger moved(std::move(fg));
  EXPECT_EQ(moved.uid(), uid);
  expect_whole(moved, "the move");
  // Every call on the moved-from ledger raises the moved-from refusal (the copy and the move of
  // it included), never a result or another check's error.
  // NOLINTBEGIN(bugprone-use-after-move): the documented moved-from state.
  const FlowId some {uid, 0};
  const std::vector<NodeId> sources {0};
  const std::vector<Amount> baseline_again(4, 1);
  const std::vector<std::pair<std::string, std::function<void()>>> calls = {
      {"uid", [&] { (void)fg.uid(); }},
      {"graph", [&] { (void)fg.graph(); }},
      {"capacity_units", [&] { (void)fg.capacity_units(); }},
      {"residual", [&] { (void)fg.residual(); }},
      {"allocated", [&] { (void)fg.allocated(); }},
      {"reserved", [&] { (void)fg.reserved(); }},
      {"unreserved", [&] { (void)fg.unreserved(); }},
      {"records", [&] { (void)fg.records(); }},
      {"record", [&] { (void)fg.record(some); }},
      {"bottleneck", [&] { (void)fg.bottleneck(sources); }},
      {"reset", [&] { fg.reset(); }},
      {"reset(residual)", [&] { fg.reset(baseline_again); }},
      {"reset(residual, unreserved)", [&] { fg.reset(baseline_again, baseline_again); }},
      {"copy", [&] { FlowLedger again(fg); }},
      {"move", [&] { FlowLedger again(std::move(fg)); }},
      {"bind", [&] { (void)LedgerAccess::bind(fg); }},
  };
  for (const auto& [name, call] : calls) {
    EXPECT_EQ(logic_error_of(call), "FlowLedger: a moved-from ledger is empty and unusable") << name;
  }
  // NOLINTEND(bugprone-use-after-move)
}

// The commit checks holds against the reservation baseline exactly as it checks carried load
// against the service baseline, in its own pass after the carried one, and publishes both or
// neither: a hold above R0, a negative hold, and a draft whose load is infeasible (its holds
// fine) each throw with nothing changed; the exchange releases the owner's old holds before the
// check; detach releases both projections.
TEST(FlowLedger, HoldsAreCheckedAgainstTheReservationBaseline) {
  auto g = shared(make_graph(2, {{0, 1, 10, 1}}));
  FlowLedger fg(g, std::vector<Amount> {10}, std::vector<Amount> {6});
  const auto owner = LedgerAccess::bind(fg);
  auto draft = [&](Amount load, Amount hold) {
    detail::Draft d;
    d.record.sources = {0};
    d.record.destinations = {1};
    d.record.bundle = all_edges(*g);
    if (load != 0) d.record.edges = {{0, load}};
    if (hold != 0) d.record.reserved_edges = {{0, hold}};
    return d;
  };
  const LedgerState empty = ledger_state(fg);
  for (const auto& [load, hold, text] :
       {std::tuple {Amount {2}, Amount {7}, std::string("ledger: infeasible reservation on edge 0 (7 of 6)")},
        std::tuple {Amount {2}, Amount {-1}, std::string("ledger: a negative hold in a commit")},
        std::tuple {Amount {11}, Amount {6}, std::string("ledger: infeasible commit on edge 0 (11 of 10)")}}) {
    std::vector<detail::Draft> ds {draft(load, hold)};
    EXPECT_EQ(logic_error_of([&] { (void)LedgerAccess::commit(fg, owner, ds); }), text);
    expect_same(ledger_state(fg), empty, text);
  }
  std::vector<detail::Draft> ds {draft(2, 6)};
  const auto ids = LedgerAccess::commit(fg, owner, ds);
  EXPECT_EQ(fg.reserved(), (std::vector<Amount> {6}));
  EXPECT_EQ(fg.unreserved(), (std::vector<Amount> {0}));
  EXPECT_EQ(fg.allocated(), (std::vector<Amount> {2}));
  // An exchange retires the old record and its hold before checking the new one's: 6 again fits.
  std::vector<detail::Draft> again {draft(3, 6)};
  (void)LedgerAccess::exchange(fg, owner, again);
  EXPECT_EQ(fg.reserved(), (std::vector<Amount> {6}));
  EXPECT_EQ(fg.allocated(), (std::vector<Amount> {3}));
  EXPECT_EQ(LedgerAccess::live(fg, ids[0].second + 1).reserved_edges, (std::vector<std::pair<EdgeId, Amount>> {{0, 6}}));
  (void)LedgerAccess::detach(fg, owner);
  EXPECT_EQ(fg.reserved(), (std::vector<Amount> {0}));
  EXPECT_EQ(fg.allocated(), (std::vector<Amount> {0}));
  expect_balanced(fg, "after detach");
  LedgerAccess::release(fg, owner);
}

TEST(FlowLedger, FlowIdsOfOtherLedgersAndUnknownIds) {
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger a(g), b(g);
  RoutePolicy p(a, 0, 1, on_support(all_edges(*g)));
  (void)p.place(1);
  const auto id = *p.flow_id();
  EXPECT_EQ(a.record(id).delivered, 1);
  EXPECT_EQ(error_of([&] { (void)b.record(id); }), "FlowId (" + std::to_string(a.uid()) + ", " + std::to_string(id.second) +
                                                     ") belongs to another ledger (uid " + std::to_string(b.uid()) + ")");
  EXPECT_EQ(error_of([&] { (void)a.record(FlowId {a.uid(), 7}); }),
            "FlowId (" + std::to_string(a.uid()) + ", 7) is retired or unknown");
  // A policy's record() reads its own records only: another policy's record of the same ledger is
  // refused, and the ledger's own checks come first.
  RoutePolicy q(a, 0, 1, on_support(all_edges(*g)));
  EXPECT_EQ(p.policy().record(id).delivered, 1);
  EXPECT_EQ(error_of([&] { (void)q.policy().record(id); }),
            "FlowId (" + std::to_string(a.uid()) + ", " + std::to_string(id.second) + ") is not a record of this policy");
  EXPECT_EQ(error_of([&] { (void)q.policy().record(FlowId {b.uid(), id.second}); }),
            "FlowId (" + std::to_string(b.uid()) + ", " + std::to_string(id.second) + ") belongs to another ledger (uid " +
                std::to_string(a.uid()) + ")");
}

TEST(FlowLedger, RecordsCarryTheContractFields) {
  // Lossy EqualBalanced on the diamond 10/10/1/10: 10 then 1 (the reference fixture
  // lossy_diamond_10_then_1, tests/reference/policy/p2_transaction/fixtures.py).
  auto g = shared(make_graph(4, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 1, 1}, {2, 3, 10, 1}}));
  const auto inv = id_of_input(*g);
  FlowLedger fg(g);
  RoutePolicy p(fg, 0, 3, split_route(fg, 0, 3, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossy));
  auto r = p.place(10);
  EXPECT_EQ(r.placed, 6);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_EQ(r.call_drops, (std::vector<std::pair<EdgeId, Amount>>{{inv[2], 4}}));
  r = p.place(1);
  EXPECT_EQ(r.placed, 0);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_EQ(*r.remaining, 1);
  // The call's own drops; the record keeps the cumulative 4 + 1 below.
  EXPECT_EQ(r.call_drops, (std::vector<std::pair<EdgeId, Amount>>{{inv[2], 1}}));
  EXPECT_EQ(alloc_in(fg), (std::vector<Amount>{6, 5, 1, 5}));
  const auto rec = fg.records().at(0);
  EXPECT_EQ(rec.flow_id, (FlowId {fg.uid(), 0}));
  EXPECT_EQ(rec.owner, 0);
  EXPECT_EQ(rec.sources, (std::vector<NodeId>{0}));
  EXPECT_EQ(rec.destinations, (std::vector<NodeId>{3}));
  EXPECT_EQ(rec.admission, Admission::Lossy);
  ASSERT_TRUE(rec.bundle);
  EXPECT_EQ(rec.bundle->via_edges().size(), 4u);
  EXPECT_TRUE(*rec.bundle == *all_edges(*g));
  EXPECT_EQ(rec.route_cost, 2);
  ASSERT_TRUE(rec.offered.has_value());
  EXPECT_EQ(*rec.offered, (Counter128 {0, 11}));
  EXPECT_EQ(rec.delivered, 6);
  ASSERT_EQ(rec.drops.size(), 1u);
  EXPECT_EQ(rec.drops[0].first, inv[2]);
  EXPECT_EQ(rec.drops[0].second, (Counter128 {0, 5}));
  // +inf on a fresh ledger (CostOnly's first-hop experiment over the same DAG): offered is the
  // injection (20), delivered 11, drops 9.
  FlowLedger f2(g);
  FlowPolicy q(f2, std::vector<NodeId> {0}, std::vector<NodeId> {3},
               FlowPolicyConfig {CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossy});
  r = q.place(std::nullopt);
  EXPECT_EQ(r.placed, 11);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  EXPECT_EQ(*f2.records()[0].offered, (Counter128 {0, 20}));
  // Lossless records carry no offered.
  FlowLedger f3(g);
  RoutePolicy m(f3, 0, 3, on_support(all_edges(*g)));
  (void)m.place(3);
  EXPECT_FALSE(f3.records()[0].offered.has_value());
  EXPECT_TRUE(f3.records()[0].drops.empty());
}

TEST(FlowLedger, CountersBeyondSixtyFourBits) {
  // Oversized offers on a link of capacity 1: of nine 2^62 - 1 offers one unit is delivered and
  // every other unit dropped: offered and drops grow beyond 2^64.
  auto g = shared(make_graph(2, {{0, 1, 1, 1}}));
  FlowLedger fg(g);
  RoutePolicy p(fg, 0, 1, split_route(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossy));
  const Amount big = kDomainLimit - 1;
  for (int i = 0; i < 9; ++i) (void)p.place(big);
  const auto rec = fg.records().at(0);
  // offered = 9 (2^62 - 1) = 2^65 + 2^62 - 9
  EXPECT_EQ(*rec.offered, (Counter128 {2, (std::uint64_t {1} << 62) - 9}));
  EXPECT_EQ(rec.delivered, 1);
  EXPECT_EQ(rec.drops[0].second, (Counter128 {2, (std::uint64_t {1} << 62) - 10}));
  EXPECT_EQ(p.requested(), (Counter128 {2, (std::uint64_t {1} << 62) - 9}));  // the offers' total
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
  RoutePolicy p(fg, 0, 1, split_route(fg, 0, 1, all_edges(*g), FlowPlacement::Weighted, Admission::Lossless));
  (void)p.place(2);
  const LedgerState before = ledger_state(fg);
  const auto id = p.flow_id();
  EXPECT_EQ(error_of([&] { (void)p.place(kDomainLimit); }),
            "offer must be an integer in [0, 2^62) or permitted +inf");
  EXPECT_EQ(error_of([&] { (void)p.place(-1); }),
            "offer must be an integer in [0, 2^62) or permitted +inf");
  // The failed calls changed nothing: the requested total, the +inf flag, the flow and the
  // ledger (every record field, the aggregate, the residual, the next FlowId and the owners).
  EXPECT_EQ(p.requested(), (Counter128 {0, 2}));
  EXPECT_FALSE(p.infinite());
  EXPECT_EQ(p.flow_id(), id);
  expect_same(ledger_state(fg), before, "after the rejected requests");
  EXPECT_EQ(fg.allocated(), (std::vector<Amount>{2}));
  EXPECT_EQ(fg.residual(), (std::vector<Amount>{3}));
  // Constructor validation.
  auto other = make_graph(2, {{0, 1, 5, 1}});
  EXPECT_EQ(error_of([&] {
              RoutePolicy x(fg, 0, 1, on_support(dag_of(other, {0})));
            }),
            "support belongs to another graph");
  EXPECT_EQ(error_of([&] {
              RoutePolicy x(fg, 0, 1, split_route(fg, 0, 1, dag_of(other, {0}), FlowPlacement::EqualBalanced, Admission::Lossless));
            }),
            "routes[0]: bundle belongs to another graph");
  EXPECT_EQ(error_of([&] {
              RoutePolicy x(fg, 0, 0, on_support(all_edges(*g)));
            }),
            "node 0 is in both endpoint sets: sources and destinations must not overlap");
  EXPECT_EQ(error_of([&] {
              RoutePolicy x(fg, 0, 1, on_support(all_edges(*g), Admission::Lossy));
            }),
            kLossyKernelError);
  // A rejected construction binds nothing (the accepted one above was destroyed: an orphan owner).
  EXPECT_EQ(LedgerAccess::bound_owners(fg), (std::vector<std::int64_t>{0}));
}

TEST(FlowLedger, BottleneckVersusMinCut) {
  // The ledger's bottleneck against the minimum cut: (0->1 cap 1 cost 1), (0->2 cap 2 cost 3),
  // (1->2 cap 1), (2->3 cap 2). After the MaxFlow +inf [1, 1, 1, 2] the bottleneck (forward
  // residual) is {0->1, 2->3}, not the cut {2->3}.
  auto g = shared(make_graph(4, {{0, 1, 1, 1}, {0, 2, 2, 3}, {1, 2, 1, 1}, {2, 3, 2, 1}}));
  const auto inv = id_of_input(*g);
  FlowLedger fg(g);
  RoutePolicy p(fg, 0, 3, on_support(all_edges(*g)));
  EXPECT_EQ(p.place(std::nullopt).placed, 2);
  EXPECT_EQ(alloc_in(fg), (std::vector<Amount>{1, 1, 1, 2}));
  std::vector<EdgeId> want = {inv[0], inv[3]};
  std::sort(want.begin(), want.end());
  EXPECT_EQ(fg.bottleneck(std::vector<NodeId> {0}), want);
  // The masks, on a graph whose unmasked bottleneck is not empty: (0->1 cap 1), (1->2 cap 1),
  // (0->3 cap 1), (0->4 cap 0). Owner a fills 1->2 and owner b 0->3, so only 0->1 keeps a
  // residual. Unmasked, the sources reach {0, 1}; the edges leaving that set are 1->2 and 0->3
  // (saturated, listed) and 0->4 (zero capacity, not listed).
  auto g2 = shared(make_graph(5, {{0, 1, 1, 1}, {1, 2, 1, 1}, {0, 3, 1, 1}, {0, 4, 0, 1}}));
  const auto inv2 = id_of_input(*g2);
  FlowLedger f2(g2);
  RoutePolicy a(f2, 1, 2, on_support(dag_of(*g2, {inv2[1]})));
  RoutePolicy b(f2, 0, 3, on_support(dag_of(*g2, {inv2[2]})));
  ASSERT_EQ(a.place(1).placed, 1);
  ASSERT_EQ(b.place(1).placed, 1);
  std::vector<EdgeId> want2 = {inv2[1], inv2[2]};
  std::sort(want2.begin(), want2.end());
  EXPECT_EQ(f2.bottleneck(std::vector<NodeId> {0}), want2);
  // 0->1 excluded: neither crossed (1->2 would be listed) nor listed itself; the sources reach
  // only {0}, whose leaving edges are 0->1 (excluded), 0->3 (listed) and 0->4 (zero capacity).
  auto el = g2->eligibility({}, std::vector<EdgeId> {inv2[0]});
  EXPECT_EQ(f2.bottleneck(std::vector<NodeId> {0}, &el), (std::vector<EdgeId> {inv2[2]}));
  // An excluded source contributes nothing; a foreign Eligibility throws.
  const std::vector<NodeId> xn = {0};
  auto ex = g2->eligibility(xn, {});
  EXPECT_TRUE(f2.bottleneck(std::vector<NodeId> {0}, &ex).empty());
  auto other = make_graph(3, {{0, 1, 1, 1}});
  auto foreign = other.eligibility({}, std::vector<EdgeId> {0});
  EXPECT_EQ(error_of([&] { (void)f2.bottleneck(std::vector<NodeId> {0}, &foreign); }),
            "foreign Eligibility: built for another graph");
  EXPECT_EQ(error_of([&] { (void)f2.bottleneck(std::vector<NodeId> {}); }), "sources: the node set is empty");
}

TEST(FlowLedger, TwoOwnersNeverCancelEachOther) {
  // Two MaxFlow owners never cancel each other's flow: A holds the chain 0-1-2-3; B's bundle
  // (chain plus bypasses) could reach 2 only by cancelling A's middle edge: B places 0 (record
  // at zero delivery, exhausted) and A is untouched; after A's remove() B places 2.
  auto g = shared(make_graph(4, {{0, 1, 1, 1}, {1, 2, 1, 1}, {2, 3, 1, 1}, {0, 2, 1, 3}, {1, 3, 1, 3}}));
  const auto inv = id_of_input(*g);
  FlowLedger fg(g);
  RoutePolicy a(fg, 0, 3, on_support(dag_of(*g, {inv[0], inv[1], inv[2]})));
  EXPECT_EQ(a.place(1).placed, 1);
  RoutePolicy b(fg, 0, 3, on_support(all_edges(*g)));
  auto r = b.place(2);
  EXPECT_EQ(r.placed, 0);
  EXPECT_FALSE(r.routes_changed);  // a controller program without a path: no program change
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
  // exact sum of the records' edges and residual() = baseline - allocated() within [0, baseline],
  // and every record is a flow of its own (expect_conserved, test_utils.hpp).
  std::mt19937 rng(5);
  int delivering = 0;   // placements that delivered units
  int with_records = 0;  // trials whose ledger held a record after some step
  int removed = 0;       // removals that dropped records
  int orphaned = 0;      // destructions that left records as orphans
  int lossless_carrying = 0;  // record checks of a Lossless record that carries something
  int lossy_dropping = 0;     // record checks of a Lossy record with drops
  for (int trial = 0; trial < 60; ++trial) {
    bool recorded = false;
    const int n = 4 + static_cast<int>(rng() % 3);
    std::vector<InEdge> edges;
    for (int i = 0; i < 3 * n; ++i) {
      const NodeId u = static_cast<NodeId>(rng() % (n - 1));
      const NodeId v = static_cast<NodeId>(u + 1 + rng() % static_cast<std::uint32_t>(n - 1 - u));
      edges.push_back({u, v, static_cast<Amount>(rng() % 7), 1 + static_cast<Cost>(rng() % 3)});
    }
    auto g = shared(make_graph(n, edges));
    FlowLedger fg(g);
    std::vector<std::unique_ptr<RoutePolicy>> pols;
    for (int k = 0; k < 4; ++k) {
      // Global on the route, or Pinned over it as a split route (Weighted or EqualBalanced, either
      // admission).
      const int kind = static_cast<int>(rng() % 3);
      const auto adm = kind != 0 && rng() % 3 == 0 ? Admission::Lossy : Admission::Lossless;
      try {
        auto cfg = kind == 0 ? on_support(all_edges(*g))
                             : split_route(fg, 0, n - 1, all_edges(*g),
                                           kind == 1 ? FlowPlacement::Weighted : FlowPlacement::EqualBalanced, adm);
        pols.push_back(std::make_unique<RoutePolicy>(fg, 0, n - 1, std::move(cfg)));
      } catch (const std::invalid_argument&) {
        break;  // no walk at all
      }
    }
    for (int step = 0; step < 12 && !pols.empty(); ++step) {
      auto& p = pols[rng() % pols.size()];
      if (!p) continue;
      const int op = static_cast<int>(rng() % 6);
      const auto records = fg.records().size();
      const auto orphans = LedgerAccess::orphan_owners(fg).size();
      if (op == 0) {
        p->remove();
        removed += fg.records().size() < records ? 1 : 0;
      } else if (op == 1) {
        p.reset();
        orphaned += LedgerAccess::orphan_owners(fg).size() > orphans ? 1 : 0;
      } else if (!p->infinite()) {
        // +inf only Lossless (a Lossy forwarding of a program needs a finite offer).
        const bool inf = op == 2 && !p->lossy();
        delivering += p->place(inf ? std::optional<Amount>() : std::optional<Amount>(rng() % 9)).placed > 0;
      }
      recorded = recorded || !fg.records().empty();
      const std::string where = "trial " + std::to_string(trial) + " step " + std::to_string(step);
      expect_balanced(fg, where);
      for (const auto& r : fg.records()) {
        const bool carries = expect_conserved(*g, r, where + " record " + std::to_string(r.flow_id.second));
        if (r.admission == Admission::Lossless) lossless_carrying += carries ? 1 : 0;
        if (r.admission == Admission::Lossy) lossy_dropping += r.drops.empty() ? 0 : 1;
      }
      const auto res = fg.residual();
      const auto allocated = fg.allocated();
      for (std::size_t e = 0; e < res.size(); ++e) {
        EXPECT_EQ(res[e] + allocated[e], g->capacity_units()[e]) << where << " edge " << e;
      }
      ASSERT_FALSE(::testing::Test::HasFailure()) << where;
    }
    with_records += recorded;
  }
  // The trials did real work: a generator change that left every bundle without a walk, or every
  // placement empty, would otherwise pass with nothing compared.
  EXPECT_GT(delivering, 100);
  EXPECT_GT(with_records, 50);
  EXPECT_GT(removed, 25);
  EXPECT_GT(orphaned, 40);
  EXPECT_GT(lossless_carrying, 100);
  EXPECT_GT(lossy_dropping, 20);
}

// expect_conserved's clauses, each broken in a copy of a real record by a change that breaks no
// other clause: the check reports exactly that clause's violations, so a check that lost the
// clause would report fewer. The Lossless record is the max flow of the diamond 10/10/1/10 (1 unit
// over 0-1-3, 10 over 0-2-3); the Lossy one is RecordsCarryTheContractFields' record (offered 11,
// 6 delivered, 5 dropped on 1 -> 3, so offered covers what was delivered and dropped exactly). A
// drop at the source is outside the transit balance (the source's own drops are what offered
// leaves over).
TEST(FlowLedger, TheRecordCheckReportsEachBrokenClause) {
  auto g = shared(make_graph(4, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 1, 1}, {2, 3, 10, 1}}));
  const auto inv = id_of_input(*g);
  FlowLedger f1(g);
  FlowLedger f2(g);
  RoutePolicy max_flow(f1, 0, 3, on_support(all_edges(*g)));
  ASSERT_EQ(max_flow.place(std::nullopt).placed, 11);
  RoutePolicy lossy(f2, 0, 3, split_route(f2, 0, 3, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossy));
  ASSERT_EQ(lossy.place(10).placed, 6);
  ASSERT_EQ(lossy.place(1).placed, 0);
  const FlowRecord lossless_record = f1.records().at(0);
  const FlowRecord lossy_record = f2.records().at(0);
  ASSERT_EQ(lossless_record.edges,
            (std::vector<std::pair<EdgeId, Amount>> {{inv[0], 1}, {inv[1], 10}, {inv[2], 1}, {inv[3], 10}}));
  ASSERT_EQ(lossless_record.delivered, 11);
  ASSERT_EQ(lossy_record.edges,
            (std::vector<std::pair<EdgeId, Amount>> {{inv[0], 6}, {inv[1], 5}, {inv[2], 1}, {inv[3], 5}}));
  ASSERT_EQ(lossy_record.drops, (std::vector<std::pair<EdgeId, Counter128>> {{inv[2], Counter128 {0, 5}}}));
  ASSERT_EQ(lossy_record.offered, (std::optional<Counter128> {Counter128 {0, 11}}));
  ASSERT_EQ(lossy_record.delivered, 6);
  // The amount on the edge of input position i.
  auto on = [&](FlowRecord& r, std::size_t i) -> Amount& {
    const auto it = std::find_if(r.edges.begin(), r.edges.end(), [&](const auto& x) { return x.first == inv[i]; });
    if (it == r.edges.end()) throw std::logic_error("the record has no amount on that edge");
    return it->second;
  };
  bool carries = false;
  EXPECT_EQ(failures_of([&] { carries = expect_conserved(*g, lossless_record, "lossless"); }), 0);
  EXPECT_TRUE(carries);
  EXPECT_EQ(failures_of([&] { carries = expect_conserved(*g, lossy_record, "lossy"); }), 0);
  EXPECT_TRUE(carries);
  FlowRecord empty = lossless_record;
  empty.edges.clear();
  empty.delivered = 0;
  EXPECT_EQ(failures_of([&] { carries = expect_conserved(*g, empty, "empty"); }), 0);
  EXPECT_FALSE(carries);

  struct Broken {
    std::string clause;
    std::function<void(FlowRecord&)> change;
    int failures;
  };
  const std::vector<Broken> lossless_cases = {
      {"positive edge amounts (the path 0-1-3 emptied, delivered with it)",
       [&](FlowRecord& r) {
         on(r, 0) = 0;
         on(r, 2) = 0;
         r.delivered -= 1;
       },
       2},
      {"edges in ascending EdgeId", [](FlowRecord& r) { std::swap(r.edges[0], r.edges[1]); }, 1},
      {"offered only on a Lossy record", [](FlowRecord& r) { r.offered = Counter128 {0, 11}; }, 1},
      {"drops only on a Lossy record", [&](FlowRecord& r) { r.drops.emplace_back(inv[0], Counter128 {0, 1}); }, 1},
      {"transit conservation (one unit moved from 0 -> 2 to 0 -> 1)",
       [&](FlowRecord& r) {
         ++on(r, 0);
         --on(r, 1);
       },
       2},
      // Lossless: delivered is both the destinations' net inflow and the sources' net outflow.
      {"delivered equal to the destinations' inflow and the sources' outflow", [](FlowRecord& r) { ++r.delivered; },
       2},
  };
  const std::vector<Broken> lossy_cases = {
      {"drops below 2^64", [](FlowRecord& r) { r.drops[0].second.hi = 1; }, 1},
      {"positive drops (the transit node's inflow lowered with it)",
       [&](FlowRecord& r) {
         r.drops[0].second.lo = 0;
         on(r, 0) = 1;
       },
       1},
      {"drops in ascending EdgeId (a second drop at the source, offered raised with it)",
       [&](FlowRecord& r) {
         r.drops.emplace_back(inv[0], Counter128 {0, 1});
         r.offered->lo += 1;
       },
       1},
      {"transit conservation net of drops (offered raised with the drop)",
       [](FlowRecord& r) {
         r.drops[0].second.lo += 1;
         r.offered->lo += 1;
       },
       1},
      {"delivered equal to the destinations' inflow (offered raised with it)",
       [](FlowRecord& r) {
         ++r.delivered;
         r.offered->lo += 1;
       },
       1},
      {"offered covers what was delivered and dropped", [](FlowRecord& r) { r.offered->lo -= 1; }, 1},
      {"offered below 2^64", [](FlowRecord& r) { r.offered->hi = 1; }, 1},
      // Without offered both the presence clause and the comparison with offered report.
      {"offered on a Lossy record", [](FlowRecord& r) { r.offered.reset(); }, 2},
  };
  for (const auto& [base, cases] : {std::pair {&lossless_record, &lossless_cases}, std::pair {&lossy_record, &lossy_cases}}) {
    for (const auto& c : *cases) {
      FlowRecord x = *base;
      c.change(x);
      EXPECT_EQ(failures_of([&] { (void)expect_conserved(*g, x, c.clause); }), c.failures) << c.clause;
    }
  }
}

TEST(FlowLedger, NoRouteCallsRecordTheEpochState) {
  // A positive request with an empty roster creates nothing, but the epoch records it: after
  // +inf a finite request (place(0) included) is rejected until remove() (the epoch state:
  // the reading of the reference in tests/reference/policy/p2_transaction/). Global on a
  // support: Pinned keeps one record per program entry even when every entry is masked
  // (below, and test_programs.py::test_all_failed_entries_still_have_records).
  auto g = shared(make_graph(2, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  auto el = std::make_shared<const Eligibility>(g->eligibility({}, std::vector<EdgeId> {0}));
  RoutePolicy p(fg, 0, 1, on_support(all_edges(*g)), el);
  EXPECT_EQ(p.place(3).termination, Termination::NoRoute);
  EXPECT_EQ(p.requested(), (Counter128 {0, 3}));
  auto r = p.place(std::nullopt);
  EXPECT_EQ(r.termination, Termination::NoRoute);
  EXPECT_FALSE(r.remaining.has_value());
  EXPECT_TRUE(p.infinite());
  EXPECT_NE(error_of([&] { (void)p.place(1); }).find("after infinite placement"), std::string::npos);
  EXPECT_NE(error_of([&] { (void)p.place(0); }).find("after infinite placement"), std::string::npos);
  EXPECT_EQ(p.requested(), (Counter128 {0, 3}));  // the rejected calls added nothing
  EXPECT_TRUE(p.infinite());
  p.remove();
  EXPECT_FALSE(p.infinite());
  EXPECT_TRUE(p.requested().is_zero());
  EXPECT_EQ(p.place(1).termination, Termination::NoRoute);
  EXPECT_TRUE(fg.records().empty());
  EXPECT_EQ(LedgerAccess::next_flow_id(fg), 0);
  // Pinned with its one entry masked: NoRoute, and the entry's record at zero.
  FlowLedger fg2(g);
  RoutePolicy q(fg2, 0, 1, split_route(fg2, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless),
                el);
  auto rq = q.place(3);
  EXPECT_EQ(rq.termination, Termination::NoRoute);
  EXPECT_EQ(rq.placed, 0);
  const auto kept = fg2.records();
  ASSERT_EQ(kept.size(), 1U);
  EXPECT_EQ(kept[0].entry, std::optional<std::int32_t> {0});
  EXPECT_EQ(kept[0].delivered, 0);
  EXPECT_TRUE(kept[0].edges.empty());
  EXPECT_EQ(kept[0].route_cost, kUnreachable);
}

TEST(FlowLedger, RouteCostIsTheBundlesMinimumWalk) {
  // A split route weighted by its residual snapshot, made while another owner fills its only
  // edge: no positive weight (record at zero delivery; NoActiveWeight), route_cost is still the
  // route's walk cost 3.
  auto g = shared(make_graph(2, {{0, 1, 2, 3}}));
  FlowLedger fg(g);
  RoutePolicy holder(fg, 0, 1, on_support(all_edges(*g)));
  EXPECT_EQ(holder.place(2).placed, 2);
  RoutePolicy w(fg, 0, 1, split_route(fg, 0, 1, all_edges(*g), FlowPlacement::Weighted, Admission::Lossless, WeightBasis::ResidualSnapshot));
  auto r = w.place(1);
  EXPECT_EQ(r.placed, 0);
  EXPECT_EQ(r.termination, Termination::NoActiveWeight);
  EXPECT_TRUE(r.routes_changed);
  ASSERT_NE(w.plan(), nullptr);
  EXPECT_TRUE(w.plan()->members.empty());
  const auto rec = fg.record(*w.flow_id());
  EXPECT_EQ(rec.route_cost, 3);
  EXPECT_EQ(rec.delivered, 0);
  // A mixed-cost bundle: the cheapest walk (1 + 1 + 1) is the route cost, not the cost-4 edge.
  auto g2 = shared(make_graph(4, {{0, 1, 1, 1}, {1, 2, 1, 1}, {2, 3, 1, 1}, {0, 3, 1, 4}}));
  FlowLedger f2(g2);
  RoutePolicy m(f2, 0, 3, on_support(all_edges(*g2)));
  (void)m.place(1);
  EXPECT_EQ(f2.records()[0].route_cost, 3);
}

TEST(FlowLedger, RemovalRestoresTheBaselineAndLeavesOtherOwnersIntact) {
  // Removing a flow restores the baseline residual (5 on a capacity-10 link), not
  // the capacity; another owner's allocation stays intact.
  auto g = shared(make_graph(2, {{0, 1, 10, 1}}));
  FlowLedger fg(g, std::vector<Amount> {5});
  RoutePolicy a(fg, 0, 1, on_support(all_edges(*g)));
  RoutePolicy b(fg, 0, 1, split_route(fg, 0, 1, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossless));
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
  // On the diamond, two calls of 1 offer each unit to the lowest
  // EdgeId first (both on the first branch) while one call of 2 splits [1, 1]: not additive.
  auto g = shared(make_graph(4, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 1, 1}, {2, 3, 10, 1}}));
  FlowLedger f1(g), f2(g);
  RoutePolicy p1(f1, 0, 3, split_route(f1, 0, 3, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossy));
  RoutePolicy p2(f2, 0, 3, split_route(f2, 0, 3, all_edges(*g), FlowPlacement::EqualBalanced, Admission::Lossy));
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

// A commit carrying several records of one owner touches an edge from more
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
  const LedgerState held = ledger_state(fg);
  std::vector<detail::Draft> over {draft_on(dag, {{0, 1}}, ids[0].second), draft_on(dag, {{0, 2}}, ids[0].second + 1),
                                   draft_on(dag, {{0, 1}})};
  EXPECT_EQ(logic_error_of([&] { (void)LedgerAccess::commit(fg, owner, over); }),
            "ledger: infeasible commit on edge 0 (4 of 3)");
  expect_same(ledger_state(fg), held, "the refused commit");
  EXPECT_EQ(fg.allocated(), (std::vector<Amount>{3}));
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
  // edges; detaches between commits. After every call the ledger must equal an independent model
  // (the aggregate, every live record's owner and edge amounts, the next id), and a rejected
  // commit must name an edge over its capacity and change nothing.
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
  int moved = 0;     // accepted commits that changed the aggregate
  int rejected = 0;  // commits the ledger rejected
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
      const LedgerState before_state = ledger_state(fg);
      if (feasible) {
        // The new records take the next ids in draft order.
        std::vector<FlowId> want_ids;
        for (int k = 0; k < fresh; ++k) want_ids.emplace_back(fg.uid(), next_id + k);
        EXPECT_EQ(LedgerAccess::commit(fg, owners[o], drafts), want_ids) << "step " << step;
        model[o] = after;
        next_id += fresh;
        moved += fg.allocated() != before ? 1 : 0;
      } else {
        // The refusal names one edge over its capacity, with its new total.
        std::vector<std::string> refusals;
        for (std::size_t e = 0; e < m; ++e) {
          if (total[e] <= g->capacity_units()[e]) continue;
          refusals.push_back("ledger: infeasible commit on edge " + std::to_string(e) + " (" + std::to_string(total[e]) +
                             " of " + std::to_string(g->capacity_units()[e]) + ")");
        }
        const std::string got = logic_error_of([&] { (void)LedgerAccess::commit(fg, owners[o], drafts); });
        EXPECT_TRUE(std::find(refusals.begin(), refusals.end(), got) != refusals.end()) << "step " << step << ": " << got;
        ++rejected;
        expect_same(ledger_state(fg), before_state, "step " + std::to_string(step));
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
    // Every live record by id: its owner and its edge amounts are the model's.
    std::map<std::int64_t, std::pair<std::int64_t, std::vector<std::pair<EdgeId, Amount>>>> want_records;
    for (std::size_t q = 0; q < owners.size(); ++q) {
      for (const auto& [id, a] : model[q]) want_records[id] = {owners[q], as_pairs(a)};
    }
    const auto records = fg.records();
    ASSERT_EQ(records.size(), want_records.size()) << "step " << step;
    for (const auto& r : records) {
      const auto it = want_records.find(r.flow_id.second);
      ASSERT_NE(it, want_records.end()) << "step " << step << ": record " << r.flow_id.second;
      EXPECT_EQ(r.owner, it->second.first) << "step " << step << ": record " << r.flow_id.second;
      EXPECT_EQ(r.edges, it->second.second) << "step " << step << ": record " << r.flow_id.second;
    }
  }
  EXPECT_GT(moved, 500);
  EXPECT_GT(rejected, 500);
  for (auto owner : owners) LedgerAccess::release(fg, owner);
}

TEST(FlowLedger, FreshSingleRecordCommitsMatchTheGeneralRule) {
  // A commit of one new record by an owner without live records takes a shortcut
  // when its edges strictly ascend (no aggregation). Against an independent model of the general
  // rule, for drafts that ascend, repeat an edge, run backwards, hold zeros or negatives, and
  // overflow at any position (another owner holding part of the capacity): the same outcome, the
  // same message (negative amounts first, then the first infeasible edge in first-touch order with
  // its aggregated total), and nothing changed by a rejected commit.
  std::mt19937 rng(1007);
  auto g = shared(make_graph(6, {{0, 1, 4, 1}, {1, 2, 3, 1}, {2, 3, 5, 1}, {0, 3, 2, 1}, {3, 4, 4, 1}, {4, 5, 3, 1}}));
  const auto dag = all_edges(*g);
  const auto m = static_cast<std::size_t>(g->num_edges());
  FlowLedger fg(g);
  const auto other = LedgerAccess::bind(fg);
  std::vector<detail::Draft> held {draft_on(dag, {{0, 1}, {2, 2}, {5, 1}})};
  (void)LedgerAccess::commit(fg, other, held);
  int shortcut = 0;    // accepted non-empty drafts whose edges strictly ascend
  int general = 0;     // accepted drafts that repeat or reorder an edge
  int negative = 0;    // rejected for a negative amount
  int infeasible = 0;  // rejected for an infeasible edge
  for (int step = 0; step < 4000; ++step) {
    const auto owner = LedgerAccess::bind(fg);
    std::vector<std::pair<EdgeId, Amount>> edges;
    const auto k = rng() % 5;
    for (std::size_t i = 0; i < k; ++i) {
      edges.emplace_back(static_cast<EdgeId>(rng() % m), rng() % 9 == 0 ? Amount {-1} : static_cast<Amount>(rng() % 4));
    }
    const int shape = static_cast<int>(rng() % 3);  // ascending, sorted with repeats kept, as drawn
    if (shape < 2) std::sort(edges.begin(), edges.end());
    if (shape == 0) {
      edges.erase(std::unique(edges.begin(), edges.end(), [](const auto& a, const auto& b) { return a.first == b.first; }),
                  edges.end());
    }
    const bool strictly = std::adjacent_find(edges.begin(), edges.end(), [](const auto& a, const auto& b) {
                            return a.first >= b.first;
                          }) == edges.end();
    // The general rule: negatives first; then each touched edge in first-touch order with its sum.
    std::string want;
    for (const auto& [e, x] : edges) {
      if (x < 0) want = "ledger: a negative allocation in a commit";
    }
    std::vector<Amount> delta(m, 0);
    std::vector<EdgeId> touched;
    for (const auto& [e, x] : edges) {
      if (std::find(touched.begin(), touched.end(), e) == touched.end()) touched.push_back(e);
      delta[static_cast<std::size_t>(e)] += x;
    }
    const auto before = fg.allocated();
    const LedgerState before_state = ledger_state(fg);
    const auto next_id = LedgerAccess::next_flow_id(fg);
    if (want.empty()) {
      for (EdgeId e : touched) {
        const auto ei = static_cast<std::size_t>(e);
        const Amount next = before[ei] + delta[ei];
        if (next < 0 || next > g->capacity_units()[ei]) {
          want = "ledger: infeasible commit on edge " + std::to_string(e) + " (" + std::to_string(next) + " of " +
                 std::to_string(g->capacity_units()[ei]) + ")";
          break;
        }
      }
    }
    std::vector<detail::Draft> drafts {draft_on(dag, edges)};
    if (!want.empty()) {
      try {
        (void)LedgerAccess::commit(fg, owner, drafts);
        ADD_FAILURE() << "step " << step << ": accepted, want " << want;
      } catch (const std::logic_error& e) {
        EXPECT_EQ(std::string(e.what()), want) << "step " << step;
        ++(std::string(e.what()).find("negative") != std::string::npos ? negative : infeasible);
      }
      expect_same(ledger_state(fg), before_state, "step " + std::to_string(step));
    } else {
      const auto ids = LedgerAccess::commit(fg, owner, drafts);
      ASSERT_EQ(ids.size(), 1U);
      EXPECT_EQ(ids[0], (FlowId {fg.uid(), next_id}));
      const auto rec = fg.record(ids[0]);
      std::vector<std::pair<EdgeId, Amount>> kept;
      for (const auto& p : edges) {
        if (p.second != 0) kept.push_back(p);
      }
      EXPECT_EQ(rec.edges, kept) << "step " << step;
      EXPECT_EQ(rec.owner, owner);
      for (std::size_t e = 0; e < m; ++e) EXPECT_EQ(fg.allocated()[e], before[e] + delta[e]) << "step " << step;
      if (strictly && !edges.empty()) ++shortcut;
      if (!strictly) ++general;
      if (rng() % 2 == 0) (void)LedgerAccess::detach(fg, owner);
    }
    LedgerAccess::release(fg, owner);
    if (::testing::Test::HasFailure()) return;
  }
  EXPECT_GT(shortcut, 140);
  EXPECT_GT(general, 10);
  EXPECT_GT(negative, 300);
  EXPECT_GT(infeasible, 300);
  LedgerAccess::release(fg, other);
  // An owner with a live record, and a draft naming an id the owner does not hold, keep the
  // general path and its rule (a fresh ledger).
  FlowLedger fg2(g);
  const auto owner = LedgerAccess::bind(fg2);
  std::vector<detail::Draft> one {draft_on(dag, {{3, 1}})};
  const auto ids = LedgerAccess::commit(fg2, owner, one);
  std::vector<detail::Draft> without_it {draft_on(dag, {{4, 1}})};
  const std::string must_carry = "ledger: a commit must carry every live record of its policy";
  EXPECT_EQ(logic_error_of([&] { (void)LedgerAccess::commit(fg2, owner, without_it); }), must_carry);
  const auto fresh = LedgerAccess::bind(fg2);
  std::vector<detail::Draft> claims {draft_on(dag, {{4, 1}}, ids[0].second)};
  EXPECT_EQ(logic_error_of([&] { (void)LedgerAccess::commit(fg2, fresh, claims); }), must_carry);
  EXPECT_EQ(fg2.records().size(), 1U);
  LedgerAccess::release(fg2, owner);
  LedgerAccess::release(fg2, fresh);
}
