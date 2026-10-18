// Allocation-failure atomicity of the ledger and of a one-bundle Pinned policy's transaction.
// The fault-injection executable replaces the global operator new so that one chosen allocation
// throws std::bad_alloc (fault_injection.cpp). Every operation is swept over every allocation it
// makes: after each injected failure the ledger (aggregate, residual, every record field, bound
// and orphaned owners, next FlowId) and the policy's epoch state (requested, +inf, placed demand,
// and its flow's id, route cost and plan: route, members, weights) must equal their state before
// the call, the per-thread scratch must not be poisoned (a clean retry gives the reference
// result), and remove() must leave the state the reference's own remove() leaves (the policy's
// records gone, their FlowIds retired). The ledger's own calls are swept the same way: reset() and
// reset(residual) change nothing and leak nothing when they fail; a failed construction or copy
// of a ledger, and a failed construction of a policy, build nothing and leak nothing; a move
// changes its source only once it can no longer fail; the completing copy and move hold
// everything the source held; remove() and a policy's destruction allocate nothing. One check
// belongs to the flow core: its node map never allocates when it adds a node (the last test but
// one). The first test checks the comparators of state_comparators.hpp themselves: each must
// report a copy with any one field changed; the next two check the hook of fault_injection.cpp: it
// reports a swallowed or a foreign std::bad_alloc, its nothrow form never fails, and any other
// exception of a swept call leaves it disarmed.
#include <gtest/gtest.h>

#include <algorithm>
#include <functional>
#include <memory>
#include <new>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/node_map.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "fault_injection.hpp"
#include "route_policy.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using detail::LedgerAccess;

namespace {

// The policies on one route the sweeps run: Global with the route as its support, and Pinned over
// the route as a supplied split route under each kernel, weight basis and admission.
struct Kernel {
  std::string name;
  Admission admission;
  std::function<FlowPolicyConfig(const FlowLedger&, NodeId, NodeId, std::shared_ptr<const PredDAG>)> config;
};

const std::vector<Kernel>& kernels() {
  using D = std::shared_ptr<const PredDAG>;
  auto split = [](FlowPlacement placement, Admission admission, WeightBasis weights) {
    return [=](const FlowLedger& fg, NodeId s, NodeId d, D dag) {
      return split_route(fg, s, d, std::move(dag), placement, admission, weights);
    };
  };
  static const std::vector<Kernel> all = {
      {"global support", Admission::Lossless, [](const FlowLedger&, NodeId, NodeId, D dag) { return on_support(std::move(dag)); }},
      {"weighted snapshot", Admission::Lossless, split(FlowPlacement::Weighted, Admission::Lossless, WeightBasis::ResidualSnapshot)},
      {"weighted capacity", Admission::Lossless, split(FlowPlacement::Weighted, Admission::Lossless, WeightBasis::Capacity)},
      {"equal balanced", Admission::Lossless, split(FlowPlacement::EqualBalanced, Admission::Lossless, WeightBasis::Capacity)},
      {"equal balanced lossy", Admission::Lossy, split(FlowPlacement::EqualBalanced, Admission::Lossy, WeightBasis::Capacity)},
      {"weighted capacity lossy", Admission::Lossy, split(FlowPlacement::Weighted, Admission::Lossy, WeightBasis::Capacity)},
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
  return {
      {"edge", shared(make_graph(2, {{0, 1, 9, 1}})), 0, 1},
      {"chain", shared(make_graph(3, {{0, 1, 9, 1}, {1, 2, 9, 1}})), 0, 2},
      // Split, reconvergence, a parallel pair and a costlier branch.
      {"diamond", shared(make_graph(5, {{0, 1, 7, 1}, {0, 2, 5, 1}, {1, 3, 4, 1}, {1, 3, 3, 1}, {2, 3, 9, 1},
                                        {3, 4, 20, 1}, {0, 4, 2, 5}})),
       0, 4},
  };
}

// One scenario: a third party's 1-unit EqualBalanced placement on the same bundle, then the
// policy's earlier calls, then the call under test.
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

// The sweeps' only "nothing changed" comparisons compare every field: same_record and the
// expect_same overloads of the ledger, a result and a policy's epoch report a copy with any one
// field changed, and expect_balanced reports an aggregate that leaves the sum of the records or
// [0, baseline] on one edge.
TEST(LedgerFaults, TheSweepComparatorsCompareEveryField) {
  const auto g = shared(make_graph(3, {{0, 1, 9, 1}, {1, 2, 9, 1}, {0, 2, 4, 3}}));
  FlowLedger fg(g);
  auto p = std::make_unique<FlowPolicy>(fg, v({0}), v({2}), cost_only_config());
  const PlacementResult r = p->place(5, OutputOptions {true, true, true});
  const LedgerState st = ledger_state(fg);
  ASSERT_EQ(st.records.size(), 1U);
  const FlowRecord& rec = st.records.front();
  ASSERT_TRUE(rec.bundle != nullptr);
  ASSERT_TRUE(same_record(rec, rec));
  const auto other_bundle = all_edges(*g);
  ASSERT_FALSE(*other_bundle == *rec.bundle);
  const std::vector<std::function<void(FlowRecord&)>> records = {
      [](auto& x) { ++x.flow_id.second; },
      [](auto& x) { ++x.owner; },
      [](auto& x) { x.sources.push_back(1); },
      [](auto& x) { x.destinations.push_back(1); },
      [](auto& x) { x.admission = Admission::Lossy; },
      [&](auto& x) { x.bundle = other_bundle; },
      [](auto& x) { ++x.route_cost; },
      [](auto& x) { x.offered = Counter128 {}; },
      [](auto& x) { ++x.delivered; },
      [](auto& x) { x.edges.emplace_back(2, 1); },
      [](auto& x) { x.drops.emplace_back(2, Counter128 {}); },
  };
  ASSERT_EQ(rec.admission, Admission::Lossless);
  ASSERT_FALSE(rec.offered.has_value());
  for (std::size_t k = 0; k < records.size(); ++k) {
    FlowRecord x = rec;
    records[k](x);
    EXPECT_FALSE(same_record(x, rec)) << "record field " << k;
  }
  // The uid clause of ledger_state: the records' own uid passes, another fails.
  {
    auto own = fg.records();
    EXPECT_EQ(failures_of([&] { expect_uid_then_clear(own, fg.uid()); }), 0);
    EXPECT_EQ(own.front().flow_id.first, 0U);
    auto foreign = fg.records();
    EXPECT_EQ(failures_of([&] { expect_uid_then_clear(foreign, fg.uid() + 1); }), 1);
  }
  EXPECT_EQ(failures_of([&] { expect_same(st, st, "same"); }), 0);
  const std::vector<std::function<void(LedgerState&)>> states = {
      [](auto& x) { ++x.allocated[0]; },
      [](auto& x) { ++x.residual[0]; },
      [](auto& x) { x.records.emplace_back(); },
      [](auto& x) { ++x.records[0].delivered; },
      [](auto& x) { ++x.next_id; },
      [](auto& x) { x.bound.push_back(7); },
      [](auto& x) { x.orphans.push_back(7); },
  };
  for (std::size_t k = 0; k < states.size(); ++k) {
    LedgerState x = st;
    states[k](x);
    EXPECT_GT(failures_of([&] { expect_same(x, st, "state"); }), 0) << "ledger field " << k;
  }
  ASSERT_TRUE(r.edge_flows.has_value() && r.residual.has_value() && r.bottleneck.has_value());
  EXPECT_EQ(failures_of([&] { expect_same(r, r, "same"); }), 0);
  const std::vector<std::function<void(PlacementResult&)>> results = {
      [](auto& x) { ++x.placed; },
      [](auto& x) { x.remaining = Amount {1}; },
      [](auto& x) { ++x.delivery_delta; },
      [](auto& x) { x.call_drops.emplace_back(2, 1); },
      [](auto& x) { ++x.call_source_loss; },
      [](auto& x) { x.termination = Termination::NoRoute; },
      [](auto& x) { x.routes_changed = !x.routes_changed; },
      [](auto& x) { x.weights_changed = !x.weights_changed; },
      [](auto& x) { x.requested.reset(); },
      [](auto& x) { x.offered = Counter128 {0, 1}; },
      [](auto& x) { ++x.admitted.lo; },
      [](auto& x) { ++x.delivered; },
      [](auto& x) { x.unadmitted.reset(); },
      [](auto& x) { ++x.source_loss.lo; },
      [](auto& x) { x.drops.emplace_back(2, Counter128 {0, 1}); },
      [](auto& x) { x.by_source.emplace_back(); },
      [](auto& x) { ++x.by_source[0].delivered; },
      [](auto& x) { x.total_cost = Counter128 {0, 1}; },
      [](auto& x) { x.cost_distribution = std::vector<std::pair<Cost, Amount>> {{9, 1}}; },
      [](auto& x) { x.edge_flows->push_back(0); },
      [](auto& x) { x.residual->push_back(0); },
      [](auto& x) { x.bottleneck->push_back(0); },
      [](auto& x) { x.min_cut = std::vector<EdgeId> {0}; },
  };
  ASSERT_NE(r.remaining, std::optional<Amount> {1});
  ASSERT_NE(r.termination, Termination::NoRoute);
  ASSERT_FALSE(r.cost_distribution.has_value());  // CostOnly: Global's curve only
  ASSERT_FALSE(r.min_cut.has_value());
  ASSERT_FALSE(r.offered.has_value() || r.total_cost.has_value());  // Lossless CostOnly EqualBalanced
  ASSERT_TRUE(r.requested.has_value() && r.unadmitted.has_value() && !r.by_source.empty());
  for (std::size_t k = 0; k < results.size(); ++k) {
    PlacementResult x = r;
    results[k](x);
    EXPECT_GT(failures_of([&] { expect_same(x, r, "result"); }), 0) << "result field " << k;
  }
  // The epoch: every field of the call state and of a flow.
  const Epoch ep = epoch(*p);
  ASSERT_FALSE(ep.cohorts.empty() || ep.cohorts.front().empty());
  EXPECT_EQ(failures_of([&] { expect_same(ep, ep, "same"); }), 0);
  const std::vector<std::function<void(Epoch&)>> epochs = {
      [](auto& x) { ++x.requested.hi; },
      [](auto& x) { ++x.requested.lo; },
      [](auto& x) { ++x.offered.lo; },
      [](auto& x) { x.infinite = !x.infinite; },
      [](auto& x) { x.prepared = !x.prepared; },
      [](auto& x) { ++x.placed; },
      [](auto& x) {
        if (!x.curve) x.curve.emplace();
        x.curve->emplace_back(9, 1);
      },
      [](auto& x) { x.program.present = !x.program.present; },
      [](auto& x) { x.program.tunnels.emplace_back(); },
      [](auto& x) { x.program.members.emplace_back(); },
      [](auto& x) { x.support.push_back(0); },
      [](auto& x) { x.head_lost.emplace_back(); },
      [](auto& x) { x.by_tunnel.emplace_back(); },
      [](auto& x) { x.cohorts.emplace_back(); },
      [](auto& x) { ++x.cohorts[0][0].id; },
      [](auto& x) { ++x.cohorts[0][0].cost; },
      [](auto& x) { x.cohorts[0][0].route.push_back(0); },
      [](auto& x) { x.cohorts[0][0].members.push_back(0); },
      [](auto& x) { x.cohorts[0][0].weights.push_back(1); },
  };
  for (std::size_t k = 0; k < epochs.size(); ++k) {
    Epoch x = ep;
    epochs[k](x);
    EXPECT_GT(failures_of([&] { expect_same(x, ep, "epoch"); }), 0) << "epoch field " << k;
  }
  // expect_balanced: one clause broken at a time on a copy of the ledger (its record an orphan:
  // a ledger with a bound policy cannot be copied), through the copy's stored aggregate and
  // record (the copy is not const, so writing through the borrowed views is defined).
  p.reset();
  EXPECT_EQ(failures_of([&] { expect_balanced(fg, "balanced"); }), 0);
  const EdgeId e = rec.edges.front().first;
  const auto ei = static_cast<std::size_t>(e);
  const std::vector<std::pair<std::string, std::function<void(FlowLedger&)>>> ledgers = {
      {"aggregate above the sum of the records",
       [&](FlowLedger& x) { ++const_cast<Amount&>(detail::LedgerAccess::allocated(x)[ei]); }},
      {"record and aggregate above the baseline",
       [&](FlowLedger& x) {
         const Amount over = detail::LedgerAccess::baseline(x)[ei] + 1 - detail::LedgerAccess::allocated(x)[ei];
         const_cast<Amount&>(detail::LedgerAccess::allocated(x)[ei]) += over;
         auto& live = const_cast<FlowRecord&>(detail::LedgerAccess::live(x, rec.flow_id.second));
         for (auto& [edge, amount] : live.edges) amount += edge == e ? over : 0;
       }},
      {"record and aggregate below zero",
       [&](FlowLedger& x) {
         const Amount under = detail::LedgerAccess::allocated(x)[ei] + 1;
         const_cast<Amount&>(detail::LedgerAccess::allocated(x)[ei]) -= under;
         auto& live = const_cast<FlowRecord&>(detail::LedgerAccess::live(x, rec.flow_id.second));
         for (auto& [edge, amount] : live.edges) amount -= edge == e ? under : 0;
       }},
  };
  for (const auto& [name, breaks] : ledgers) {
    FlowLedger x(fg);
    breaks(x);
    EXPECT_GT(failures_of([&] { expect_balanced(x, name); }), 0) << name;
  }
}

// The hook itself. fails_at's check that an injected std::bad_alloc propagated is every sweep's
// only check of it, so it must report a failure that fn swallows and a std::bad_alloc it did not
// inject; a call that allocates nothing reaches no failure. The nothrow operator new never fails
// and counts as live like the throwing one, so a call that asks for memory it can do without (the
// temporary buffer of std::stable_partition) is neither a reported failure nor a leak.
TEST(FaultHook, ReportsSwallowedAndForeignFailuresAndNeverFailsNothrow) {
  EXPECT_EQ(failures_of([] { EXPECT_FALSE(fails_at(0, [] {})); }), 0);
  // Reached and swallowed: the injected failure does not reach fails_at.
  EXPECT_EQ(failures_of([] {
              (void)fails_at(0, [] {
                try {
                  ::operator delete(::operator new(8));
                } catch (const std::bad_alloc&) {
                }
              });
            }),
            1);
  // A std::bad_alloc of fn's own while the armed failure is never reached.
  EXPECT_EQ(failures_of([] { (void)fails_at(1000, [] { throw std::bad_alloc(); }); }), 1);
  // A standard algorithm's nothrow temporary buffer: whether the library takes it through the
  // nothrow form (not a failure point) or the throwing one (a failure that propagates), the hook
  // reports nothing, and the completing call partitions.
  std::vector<std::vector<int>> v(6);
  for (int i = 0; i < 6; ++i) v[static_cast<std::size_t>(i)].push_back(i);
  const auto even = [](const std::vector<int>& x) { return x[0] % 2 == 0; };
  EXPECT_EQ(failures_of([&] { (void)fails_at(0, [&] { (void)std::stable_partition(v.begin(), v.end(), even); }); }),
            0);
  (void)std::stable_partition(v.begin(), v.end(), even);
  EXPECT_EQ(v, (std::vector<std::vector<int>> {{0}, {2}, {4}, {1}, {3}, {5}}));
  // The nothrow form under an armed failure: it allocates, counts and leaves the countdown alone.
  const long live = g_live.load();
  void* q = nullptr;
  EXPECT_EQ(failures_of([&] { EXPECT_FALSE(fails_at(0, [&] { q = ::operator new(16, std::nothrow); })); }), 0);
  ASSERT_NE(q, nullptr);
  EXPECT_EQ(g_live.load(), live + 1);
  ::operator delete(q);
  EXPECT_EQ(g_live.load(), live);
}

// A swept call that throws something other than std::bad_alloc (a validation error, a failed
// check of the library) leaves the hook with that exception, unchanged, and disarmed: an armed
// countdown would fail an allocation of whatever runs next (the test framework's report of the
// exception, a later test), so one real failure would show up as unrelated ones. Each check reads
// the countdown and disarms it in one step, so a failure here does not spread either.
TEST(FaultHook, AnyOtherExceptionLeavesTheHookDisarmed) {
  EXPECT_THROW((void)fails_at(1000, [] { throw std::invalid_argument("a validation error"); }), std::invalid_argument);
  EXPECT_LT(g_countdown.exchange(-1), 0);
  EXPECT_THROW((void)fails_in_workers_at(1000, [] { throw std::logic_error("a failed check"); }), std::logic_error);
  EXPECT_LT(g_countdown.exchange(-1), 0);
  bool bad_alloc = true;
  EXPECT_THROW((void)reaches_in_workers_at(1000, [] { throw std::invalid_argument("x"); }, bad_alloc),
               std::invalid_argument);
  EXPECT_LT(g_countdown.exchange(-1), 0);
  EXPECT_FALSE(bad_alloc);
}

TEST(LedgerFaults, EveryAllocationOfAPlacementFailsAtomically) {
  for (const auto& c : cases()) {
    const auto dag = all_edges(*c.g);
    for (const auto& k : kernels()) {
      for (const auto& sc : scenarios()) {
        // A Lossy program forwarding needs a finite offer (no first-hop experiment).
        if (k.admission == Admission::Lossy && !sc.call.has_value()) continue;
        const std::string where = c.name + " " + k.name + " " + sc.name;
        // A third party's 1-unit EqualBalanced placement on the same bundle, then the policy's
        // earlier calls.
        const auto setup = [&](FlowLedger& fg) {
          auto third = std::make_shared<RoutePolicy>(fg, c.src, c.dst, split_route(fg, c.src, c.dst, dag, FlowPlacement::EqualBalanced, Admission::Lossless));
          (void)third->place(1);
          auto p = std::make_shared<RoutePolicy>(fg, c.src, c.dst, k.config(fg, c.src, c.dst, dag));
          for (auto v : sc.before) (void)p->place(v);
          FlowPolicy* policy = &p->policy();
          return Staged {std::make_shared<std::pair<std::shared_ptr<RoutePolicy>, std::shared_ptr<RoutePolicy>>>(
                             std::move(third), std::move(p)),
                         policy};
        };
        const long injected =
            sweep_placement(c.g, where, setup, [&](FlowPolicy& p) { return p.place(sc.call); });
        if (::testing::Test::HasFailure()) return;
        EXPECT_GT(injected, 2) << where;  // not vacuous: every placement allocates
      }
    }
  }
}

// A failed policy construction binds nothing, changes nothing and leaks nothing (the thread's
// scratch is already grown by a clean construction of the same policy, so every allocation the
// failed one leaves live is a leak).
TEST(LedgerFaults, PolicyConstructionFailsWithoutBinding) {
  for (const auto& c : cases()) {
    const auto dag = all_edges(*c.g);
    for (const auto& k : kernels()) {
      {
        FlowLedger warm(c.g);
        (void)RoutePolicy(warm, c.src, c.dst, k.config(warm, c.src, c.dst, dag)).place(1);
      }
      long injected = 0;
      for (long at = 0;; ++at) {
        const std::string w = c.name + " construction allocation " + std::to_string(at);
        FlowLedger fg(c.g);
        {
          RoutePolicy other(fg, c.src, c.dst, split_route(fg, c.src, c.dst, dag, FlowPlacement::EqualBalanced, Admission::Lossless));
          (void)other.place(1);
        }  // an orphan
        const LedgerState before = ledger_state(fg);
        std::optional<RoutePolicy> p;
        const long live = g_live.load();
        if (!fails_at(at, [&] { p.emplace(fg, c.src, c.dst, k.config(fg, c.src, c.dst, dag)); })) {
          EXPECT_EQ(LedgerAccess::bound_owners(fg).size(), 1U);
          break;
        }
        ++injected;
        expect_same(ledger_state(fg), before, w);
        EXPECT_FALSE(p.has_value()) << w;
        EXPECT_EQ(g_live.load(), live) << w << ": a leak";
        if (::testing::Test::HasFailure()) return;
      }
      EXPECT_GT(injected, 0) << c.name;  // not vacuous
    }
  }
}

// reset() and reset(residual) build the new baseline before they change anything: a failure at
// any of their allocations leaves the orphaned records, the aggregate, the baseline and the next
// FlowId as they were, and the completing call gives the reference reset state.
TEST(LedgerFaults, EveryAllocationOfAResetFailsAtomically) {
  for (const auto& c : cases()) {
    const auto dag = all_edges(*c.g);
    const auto with_orphan = [&] {
      auto fg = std::make_unique<FlowLedger>(c.g);
      // The policy is destroyed without remove(): its record is an orphan.
      (void)RoutePolicy(*fg, c.src, c.dst, split_route(*fg, c.src, c.dst, dag, FlowPlacement::EqualBalanced, Admission::Lossless))
          .place(2);
      return fg;
    };
    const std::vector<Amount> half = [&] {
      std::vector<Amount> r(c.g->capacity_units().begin(), c.g->capacity_units().end());
      for (auto& x : r) x /= 2;
      return r;
    }();
    for (const bool to_capacities : {true, false}) {
      const std::string where = c.name + (to_capacities ? " reset()" : " reset(residual)");
      const auto reset = [&](FlowLedger& fg) {
        if (to_capacities) {
          fg.reset();
        } else {
          fg.reset(half);
        }
      };
      const auto ref = with_orphan();
      const LedgerState initial = ledger_state(*ref);
      ASSERT_FALSE(initial.records.empty()) << where;  // something to clear
      reset(*ref);
      const LedgerState want = ledger_state(*ref);
      EXPECT_EQ(want.residual, to_capacities ? std::vector<Amount>(c.g->capacity_units().begin(),
                                                                   c.g->capacity_units().end())
                                             : half)
          << where;
      long injected = 0;
      for (long at = 0;; ++at) {
        const std::string w = where + " allocation " + std::to_string(at);
        // The attempt's ledger ends with it (the reference above grew this thread's scratch), so
        // an allocation left live is a leak of the failed reset.
        const long live = g_live.load();
        bool failed = false;
        {
          const auto fg = with_orphan();
          expect_same(ledger_state(*fg), initial, w + " (setup)");
          failed = fails_at(at, [&] { reset(*fg); });
          if (!failed) {
            expect_same(ledger_state(*fg), want, w + " (no failure)");
          } else {
            ++injected;
            expect_same(ledger_state(*fg), initial, w + " (failed call)");
            // A clean retry is the reference reset.
            reset(*fg);
            expect_same(ledger_state(*fg), want, w + " (retry)");
          }
        }
        EXPECT_EQ(g_live.load(), live) << w << ": a leak";
        if (!failed) break;
        if (::testing::Test::HasFailure()) return;
      }
      EXPECT_GT(injected, 0) << where;  // not vacuous: the new baseline is allocated
    }
  }
}

// A failed construction of a ledger, from the capacities or from a residual, builds nothing and
// leaks nothing; the completing one holds the expected baseline and nothing else.
TEST(LedgerFaults, EveryAllocationOfALedgerConstructionFailsWithoutALeak) {
  for (const auto& c : cases()) {
    const std::vector<Amount> capacities(c.g->capacity_units().begin(), c.g->capacity_units().end());
    std::vector<Amount> half = capacities;
    for (auto& x : half) x /= 2;
    for (const bool from_residual : {false, true}) {
      const std::string where = c.name + (from_residual ? " FlowLedger(graph, residual)" : " FlowLedger(graph)");
      long injected = 0;
      for (long at = 0;; ++at) {
        const std::string w = where + " allocation " + std::to_string(at);
        std::optional<FlowLedger> fg;
        const long live = g_live.load();
        if (!fails_at(at, [&] {
              if (from_residual) {
                fg.emplace(c.g, half);
              } else {
                fg.emplace(c.g);
              }
            })) {
          ASSERT_TRUE(fg.has_value()) << w;
          EXPECT_EQ(fg->residual(), from_residual ? half : capacities) << w;
          EXPECT_EQ(fg->allocated(), std::vector<Amount>(capacities.size(), 0)) << w;
          EXPECT_TRUE(fg->records().empty()) << w;
          break;
        }
        ++injected;
        EXPECT_FALSE(fg.has_value()) << w;
        EXPECT_EQ(g_live.load(), live) << w << ": a leak";
        if (::testing::Test::HasFailure()) return;
      }
      EXPECT_GE(injected, 2) << where;  // the baseline and the aggregate (MSVC's maps add their own)
    }
  }
}

// A failed copy of a ledger constructs nothing, leaks nothing and leaves the source as it was;
// the completing copy holds the source's records, aggregate, baseline, next FlowId and owners,
// and so does the moved-to ledger. A move moves every member, so it allocates nothing with libc++
// and libstdc++ (MSVC's maps allocate a new empty node when constructed): any allocation it makes
// comes before it changes the source.
TEST(LedgerFaults, EveryAllocationOfALedgerCopyOrMoveFailsWithoutALeak) {
  for (const auto& c : cases()) {
    const auto dag = all_edges(*c.g);
    // A baseline below the capacities: a ledger that rebuilt its baseline from them would differ.
    std::vector<Amount> half(c.g->capacity_units().begin(), c.g->capacity_units().end());
    for (auto& x : half) x /= 2;
    FlowLedger fg(c.g, half);
    {
      RoutePolicy first(fg, c.src, c.dst, split_route(fg, c.src, c.dst, dag, FlowPlacement::EqualBalanced, Admission::Lossless));
      (void)first.place(1);
      RoutePolicy second(fg, c.src, c.dst, split_route(fg, c.src, c.dst, dag, FlowPlacement::Weighted, Admission::Lossy));
      (void)second.place(2);
    }  // two orphans: a copy may only be taken of a ledger without bound policies
    const LedgerState source = ledger_state(fg);
    const auto baseline_of = [](const FlowLedger& x) {
      const auto b = LedgerAccess::baseline(x);
      return std::vector<Amount>(b.begin(), b.end());
    };
    const std::vector<Amount> source_baseline = baseline_of(fg);
    const auto fg_uid = fg.uid();
    ASSERT_EQ(source.records.size(), 2U) << c.name;
    long injected = 0;
    for (long at = 0;; ++at) {
      const std::string w = c.name + " allocation " + std::to_string(at);
      std::optional<FlowLedger> copy;
      const long live = g_live.load();
      if (!fails_at(at, [&] { copy.emplace(fg); })) {
        ASSERT_TRUE(copy.has_value()) << w;
        EXPECT_NE(copy->uid(), fg.uid()) << w;
        expect_same(ledger_state(*copy), source, w + " (copy)");
        EXPECT_EQ(baseline_of(*copy), source_baseline) << w << " (copy)";
        break;
      }
      ++injected;
      EXPECT_FALSE(copy.has_value()) << w;
      EXPECT_EQ(g_live.load(), live) << w << ": a leak";
      expect_same(ledger_state(fg), source, w + " (source)");
      if (::testing::Test::HasFailure()) return;
    }
    EXPECT_GT(injected, 2) << c.name;  // not vacuous: the baseline, the aggregate and the records
    for (long at = 0;; ++at) {
      const std::string w = c.name + " move, allocation " + std::to_string(at);
      std::optional<FlowLedger> moved;
      const long live = g_live.load();
      if (!fails_at(at, [&] { moved.emplace(std::move(fg)); })) {
        ASSERT_TRUE(moved.has_value()) << w;
        EXPECT_EQ(moved->uid(), fg_uid) << w;
        EXPECT_THROW((void)fg.records(), std::logic_error) << w;  // a moved-from ledger is unusable
        // Everything the source held: the records, the aggregate, the baseline, the next FlowId
        // and the owners (the uid is kept, so the records' FlowIds compare as they are).
        expect_same(ledger_state(*moved), source, w + " (moved)");
        EXPECT_EQ(baseline_of(*moved), source_baseline) << w << " (moved)";
        moved->reset();  // the moved ledger works
        break;
      }
      EXPECT_FALSE(moved.has_value()) << w;
      EXPECT_EQ(g_live.load(), live) << w << ": a leak";
      expect_same(ledger_state(fg), source, w + " (source)");
      if (::testing::Test::HasFailure()) return;
    }
  }
}

// The per-thread scratch of the plan builders and of the commit grows with the graph. Each
// attempt runs on a new thread (its scratch empty) that first serves a 3-node graph, then the
// swept call on a larger one, so every allocation of that call fails once while its scratch grows
// from the small size; the failed call changes nothing, and a clean retry on the same graph (or
// the call itself, once no failure is reached) and then a placement give the reference's ledger,
// epoch and result (ASan/UBSan builds report any access outside a partly grown array). The same
// call on a new thread whose scratch has first grown to the large graph allocates less: the sweep
// reached the growth.
TEST(LedgerFaults, ScratchGrowthFailureLeavesNothingHalfGrown) {
  const auto small = chain_graph(3, 5);
  const auto small_dag = all_edges(*small);
  const auto large = chain_graph(3000, 5);
  const auto dag = all_edges(*large);
  const std::int32_t n = large->num_nodes();
  const auto policy = [&](FlowLedger& fg) {
    return RoutePolicy(fg, 0, n - 1, split_route(fg, 0, n - 1, dag, FlowPlacement::EqualBalanced, Admission::Lossless));
  };
  // The reference, without a failure: the empty ledger, the ledger and epoch right after the
  // construction, and the result, ledger and epoch of place(2).
  LedgerState fresh;
  LedgerState built;
  Epoch built_epoch;
  PlacementResult want;
  LedgerState want_state;
  Epoch want_epoch;
  {
    FlowLedger fg(large);
    fresh = ledger_state(fg);
    RoutePolicy p = policy(fg);
    built = ledger_state(fg);
    built_epoch = epoch(p.policy());
    want = p.place(2);
    want_state = ledger_state(fg);
    want_epoch = epoch(p.policy());
  }
  ASSERT_EQ(want.placed, 2);  // the whole chain carries the request: the comparison is not of zeros
  ASSERT_EQ(want_state.records.size(), 1U);
  // Construction (the route plan's per-node scratch) and a placement (the commit's per-edge
  // scratch and the kernels' buffers) on the small graph.
  const auto warm_up = [&] {
    FlowLedger fg(small);
    RoutePolicy p(fg, 0, 2, split_route(fg, 0, 2, small_dag, FlowPlacement::EqualBalanced, Admission::Lossless));
    EXPECT_EQ(p.place(2).placed, 2);
  };
  sweep_scratch_growth(
      "construction", warm_up,
      [&](long at, const std::string& w) {
        FlowLedger fg(large);
        std::optional<RoutePolicy> p;
        const bool failed = fails_at(at, [&] { p.emplace(policy(fg)); });
        if (failed) {
          EXPECT_FALSE(p.has_value()) << w;
          expect_same(ledger_state(fg), fresh, w);  // nothing bound
          p.emplace(policy(fg));
        }
        expect_same(ledger_state(fg), built, w);
        expect_same(epoch(p->policy()), built_epoch, w);
        expect_same(p->place(2), want, w);
        expect_same(ledger_state(fg), want_state, w);
        expect_same(epoch(p->policy()), want_epoch, w);
        expect_balanced(fg, w);
        return failed;
      },
      [&] {
        FlowLedger fg(large);
        std::optional<RoutePolicy> p(policy(fg));
        (void)p->place(2);  // the scratch grows to the large graph
        p->remove();
        p.reset();
        long warm = 0;
        while (fails_at(warm, [&] { p.emplace(policy(fg)); })) ++warm;
        return warm;
      });
  sweep_scratch_growth(
      "placement", warm_up,
      [&](long at, const std::string& w) {
        FlowLedger fg(large);
        RoutePolicy p = policy(fg);
        const LedgerState before = ledger_state(fg);
        const Epoch before_epoch = epoch(p.policy());
        PlacementResult got;
        const bool failed = fails_at(at, [&] { got = p.place(2); });
        if (failed) {
          expect_same(ledger_state(fg), before, w);
          expect_same(epoch(p.policy()), before_epoch, w);
          got = p.place(2);
        }
        expect_same(got, want, w);
        expect_same(ledger_state(fg), want_state, w);
        expect_same(epoch(p.policy()), want_epoch, w);
        expect_balanced(fg, w);
        return failed;
      },
      [&] {
        FlowLedger fg(large);
        RoutePolicy p = policy(fg);
        (void)p.place(2);  // the scratch grows to the large graph
        p.remove();
        long warm = 0;
        while (fails_at(warm, [&] { (void)p.place(2); })) ++warm;
        return warm;
      });
}

// The flow core's NodeMap::add runs once per edge endpoint in its subset solve. reset() sizes both
// of its arrays, so add() never allocates and is noexcept: it stays a few inlined instructions (an
// out-of-line call there made the solve's timing switch between two levels from one build to the
// next).
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
      auto p = std::make_unique<RoutePolicy>(fg, c.src, c.dst, k.config(fg, c.src, c.dst, dag));
      expect_removal_and_release_allocate_nothing(fg, p, c.name + " " + k.name);
    }
  }
}

TEST(LedgerFaults, PlacementSweepRejectsAPartialCommit) {
  const auto g = shared(make_graph(3, {{0, 1, 10, 1}, {1, 2, 10, 1}}));
  const auto setup = [](FlowLedger& fg) {
    auto p = std::make_shared<FlowPolicy>(fg, std::vector<NodeId> {0}, std::vector<NodeId> {2}, cost_only_config());
    return Staged {p, p.get()};
  };
  long reached = 0;
  EXPECT_EQ(failures_of([&] {
    reached = sweep_placement(g, "one call", setup, [](FlowPolicy& p) { return p.place(2); });
  }), 0);
  EXPECT_GT(reached, 0);
  EXPECT_GT(failures_of([&] {
    (void)sweep_placement(g, "partial commit", setup, [](FlowPolicy& p) {
      (void)p.place(1);
      return p.place(1);
    });
  }), 0);
}
