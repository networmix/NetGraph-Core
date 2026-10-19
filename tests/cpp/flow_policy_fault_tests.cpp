// Allocation-failure atomicity of the policy engine, the DagCache and the evaluators.
// Exhaustive sweeps (each allocation fails once) cover the placement calls of the routing and
// kernel modes below (supplied programs of split routes and of crossing tunnels, Global with a
// support among them), with and without a static Eligibility; the full-offer calls evaluate,
// replace and forward (with black holes) of every mode; the zero and empty-roster calls; the
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
// the policy's epoch state (requested, +inf, Global's cumulative curve, the installed program, the
// by_tunnel rows, and per flow its id, route cost and plan: route, members, weights) equal their
// state before the call; a clean retry gives
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
#include <numeric>
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
// shared middle, a costlier bypass, and a second destination. Two links of capacity 0 (inputs 11
// and 12): 1 -> 6 ties source 1's other routes, so its shortest support (and every split there)
// holds a zero-capacity member; 0 -> 5 is costlier than 0's routes and serves as an all-zero
// Pinned bundle (a Weighted plan without a member).
std::shared_ptr<const Graph> graph() {
  return std::make_shared<const Graph>(
      make_graph(7, {{0, 2, 6, 1}, {0, 3, 4, 1}, {1, 3, 5, 1}, {1, 4, 3, 1}, {2, 5, 3, 1}, {2, 5, 2, 1},
                     {3, 5, 4, 1}, {4, 6, 9, 1}, {3, 6, 2, 2}, {0, 6, 3, 4}, {2, 4, 2, 1}, {1, 6, 0, 2},
                     {0, 5, 0, 3}}));
}

struct Mode {
  std::string name;
  FlowPolicyConfig cfg;
  std::vector<NodeId> sources;
  std::vector<NodeId> destinations;
  // The evaluators' offer (a T1 or Pinned Lossy evaluation requires a finite one).
  std::optional<std::optional<Amount>> offered {};
  // T1's design: the evaluators' design_demand and the design of its calls.
  std::optional<std::optional<Amount>> design {};
};

// The evaluators' options of a mode.
QueryOptions query_of(const Mode& m, const Eligibility* elig = nullptr, DagCache* cache = nullptr,
                      OutputOptions outputs = {}) {
  return QueryOptions {m.cfg, elig, cache, outputs, m.design, m.offered};
}

// A Pinned Lossy forwarding needs a finite offer (no first-hop experiment).
bool finite_only(const Mode& m) {
  return std::holds_alternative<Pinned>(m.cfg.routing) && m.cfg.admission == Admission::Lossy;
}

std::vector<Mode> modes(const Graph& g) {
  std::vector<Mode> out;
  const auto det = PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity);
  const auto res = PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Residual);
  const std::vector<NodeId> S = {0, 1};
  const std::vector<NodeId> D = {5, 6};
  out.push_back({"costonly maxflow", {CostOnly {}, FlowPlacement::MaxFlow, Admission::Lossless, det}, S, D});
  out.push_back({"costonly eb", cost_only_config(), S, D});
  out.push_back({"costonly weighted residual",
                 {CostOnly {}, FlowPlacement::Weighted, Admission::Lossless, res, std::nullopt, WeightBasis::ResidualSnapshot},
                 S, D});
  out.push_back({"costonly lossy", {CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossy, det}, S, D});
  // Several sources: the deterministic single path (the pooled single-path rule rejects
  // PreferHigherPriority with multipath false).
  out.push_back({"costonly single path collapse",
                 {CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless,
                  PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Residual, false, false)},
                 S, D});
  out.push_back({"global", global_config(), S, D});
  // The equal origin of the controllers: the bisection's probe solves and the final solve.
  out.push_back({"global equal", {.routing = Global {}, .origin = Origin::Equal}, S, D});
  out.push_back({"costonly maxflow equal",
                 {CostOnly {}, FlowPlacement::MaxFlow, Admission::Lossless, det, Origin::Equal}, S, D});
  // Pooled Weighted Lossy: source 1's zero-capacity member has weight 0 (pruned), the shares
  // attributed across the shared routers.
  out.push_back({"costonly weighted lossy", {CostOnly {}, FlowPlacement::Weighted, Admission::Lossy, det}, S, D});
  // Pinned over supplied programs: three split routes of source 0 (one all-zero: a Weighted
  // route without a member), Lossless EqualBalanced and Lossy Weighted; and two sources' crossing
  // tunnels, Lossy, which meet on shared links (the general loss rule).
  const auto inv = id_of_input(g);
  const NodeId s0[1] = {0};
  const NodeId d5[1] = {5};
  const std::vector<Route> split = {Route {0, dag_of(g, {inv[0], inv[4], inv[5]}), 1},
                                    Route {0, dag_of(g, {inv[1], inv[6]}), 2}, Route {0, dag_of(g, {inv[12]}), 1}};
  RouteOptions weighted;
  weighted.placement = FlowPlacement::Weighted;
  out.push_back({"pinned eb split routes",
                 {.routing = Pinned {ForwardingProgram::from_routes(g, s0, d5, split), nullptr}}, {0}, {5}});
  out.push_back({"pinned lossy",
                 {.routing = Pinned {ForwardingProgram::from_routes(g, s0, d5, split, weighted), nullptr},
                  .admission = Admission::Lossy},
                 {0}, {5}, std::optional<Amount> {6}});
  auto walk = [&](std::vector<std::size_t> inputs) {
    std::vector<std::vector<EdgeId>> hops;
    for (auto i : inputs) hops.push_back({inv[i]});
    return std::make_shared<const Path>(g, hops);
  };
  const std::vector<Route> crossing = {Route {0, walk({0, 10, 7}), 2}, Route {0, walk({1, 8}), 1},
                                       Route {1, walk({3, 7}), 1}, Route {1, walk({2, 8}), 1}};
  out.push_back({"pinned tunnels lossy",
                 {.routing = Pinned {ForwardingProgram::from_routes(g, S, std::vector<NodeId> {6}, crossing), nullptr},
                  .admission = Admission::Lossy},
                 S, {6}, std::optional<Amount> {6}});
  out.push_back({"global support", {.routing = Global {dag_of(g, {inv[0], inv[1], inv[4], inv[6], inv[2]})}}, S, D});
  // CostOnly creation in one pass: its Collapse branch (multipath, multi_edge false,
  // on the creation view) and its Weighted weights (Capacity priority).
  out.push_back({"costonly eb multipath collapse",
                 {CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless,
                  PathSelection(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, true, false)},
                 S, D});
  out.push_back({"costonly weighted capacity", {CostOnly {}, FlowPlacement::Weighted, Admission::Lossless, det}, S, D});
  out.push_back({"costonly weighted explicit",
                 {CostOnly {}, FlowPlacement::Weighted, Admission::Lossy, det, std::nullopt, WeightBasis::Explicit,
                  std::make_shared<const std::vector<Amount>>(std::vector<Amount> {3, 1, 2, 2, 1, 4, 2, 1, 1, 5, 2, 7, 3})},
                 S, D});
  return out;
}

// T1's modes: they take no place, so their sweeps are their own (T1 below). Lossless with two
// tunnels per head; Lossy with three, PreferHigherPriority on the unreserved view without parallel
// links, a Stretch ceiling and a search minimum.
std::vector<Mode> t1_modes() {
  const std::vector<NodeId> S = {0, 1};
  const std::vector<NodeId> D = {5, 6};
  std::vector<Mode> out;
  out.push_back({"t1 lossless", {.routing = CapacityAware {Tunnels {2}}}, S, D, std::nullopt, Amount {8}});
  out.push_back({"t1 lossy",
                 {.routing = CapacityAware {Tunnels {3}, Stretch {2.0}, 1},
                  .admission = Admission::Lossy,
                  .selection = PathSelection(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, false, false)},
                 S, D, Amount {6}, Amount {9}});
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

// A third party holding one unit on the first route's edges (Global on that route as its support).
std::unique_ptr<FlowPolicy> third_party(FlowLedger& fg, const Graph& g) {
  const auto inv = id_of_input(g);
  auto t = std::make_unique<FlowPolicy>(fg, std::vector<NodeId> {0}, std::vector<NodeId> {5},
                                        FlowPolicyConfig {.routing = Global {dag_of(g, {inv[0], inv[4]})}});
  (void)t->place(1);
  return t;
}

} // namespace

// The sweeps below reach the zero-capacity paths (not vacuous): source 1's CostOnly EqualBalanced
// route holds the zero link 1 -> 6 as a member, and its sensitivity lists it; the all-zero Pinned
// bundle is usable, with no member under the Weighted kernel.
TEST(PolicyFaults, ZeroCapacityLinksTakePartInTheSweeps) {
  const auto g = graph();
  const auto inv = id_of_input(*g);
  for (const auto& m : modes(*g)) {
    FlowLedger fg(g);
    FlowPolicy p(fg, m.sources, m.destinations, m.cfg);
    (void)p.place(5);
    const auto& pe = detail::PolicyAccess::engine(p);
    if (m.name == "costonly eb") {
      // The pooled program of both sources holds the zero link as a member.
      ASSERT_TRUE(pe.epoch.program != nullptr);
      const auto& plan = *detail::ProgramAccess::data(*pe.epoch.program).contexts.front().plan;
      EXPECT_EQ(vec(plan.roots), (std::vector<NodeId> {0, 1}));
      EXPECT_TRUE(std::find(plan.members.begin(), plan.members.end(), inv[11]) != plan.members.end());
      const auto s = sensitivity(*g, m.sources, m.destinations, query_of(m), 1);
      EXPECT_TRUE(std::any_of(s.entries.begin(), s.entries.end(), [&](const auto& e) { return e.edge == inv[11]; }));
    }
    if (m.name == "pinned eb split routes" || m.name == "pinned lossy") {
      // Three split routes, three records; the all-zero route keeps its EqualBalanced member and
      // has none under the Weighted kernel.
      ASSERT_EQ(pe.view.data->contexts.size(), 3U) << m.name;
      ASSERT_EQ(pe.epoch.cohorts.size(), 3U) << m.name;
      EXPECT_EQ(pe.view.plan(2)->members.empty(), m.name == "pinned lossy") << m.name;
    }
  }
}

TEST(PolicyFaults, EveryAllocationOfEveryRoutingModeFailsAtomically) {
  const auto g = graph();
  const Eligibility no_bypass = without_bypass(*g);
  for (const auto& m : modes(*g)) {
    for (const auto& sc : scenarios()) {
      if (finite_only(m) && !sc.call.has_value()) continue;
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

// The full-offer calls on a policy that has placed (its program installed): evaluate (nothing
// published, the epoch untouched), replace (the old records retired and new ones created in one
// exchange of the ledger; to zero, and to +inf when Lossless) and, for the computing models,
// forward over the installed program, with and without black holes. Every allocation of the call
// fails once, through the shared placement sweep: a failure changes neither the ledger nor the
// epoch state and leaves no allocation live; the completing call, the retry and remove() leave the
// ledger as the reference's do.
TEST(PolicyFaults, EveryFullOfferCallFailsAtomically) {
  const auto g = graph();
  const auto inv = id_of_input(*g);
  const EdgeId dark_edges[] = {inv[4], inv[6]};
  const Eligibility dark = g->eligibility({}, dark_edges);
  struct Op {
    std::string name;
    std::function<PlacementResult(FlowPolicy&)> call;
  };
  for (const auto& m : modes(*g)) {
    const bool lossy = m.cfg.admission == Admission::Lossy;
    const bool pinned = std::holds_alternative<Pinned>(m.cfg.routing);
    std::vector<Op> ops = {{"evaluate", [](FlowPolicy& p) { return p.evaluate(6, std::nullopt, {true, true, true}); }},
                           {"replace", [](FlowPolicy& p) { return p.replace(6, std::nullopt, {true, true, true}); }},
                           {"replace to zero", [](FlowPolicy& p) { return p.replace(0); }}};
    if (!lossy) ops.push_back({"replace +inf", [](FlowPolicy& p) { return p.replace(std::nullopt); }});
    if (!pinned) {
      ops.push_back({"forward", [](FlowPolicy& p) { return p.forward(6, nullptr, {true, true, true}); }});
      ops.push_back({"forward black holes", [&](FlowPolicy& p) { return p.forward(6, &dark); }});
    }
    for (const auto& op : ops) {
      const std::string where = m.name + " / " + op.name;
      const auto setup = [&](FlowLedger& fg) {
        auto third = std::shared_ptr<FlowPolicy>(third_party(fg, *g));
        auto p = std::make_shared<FlowPolicy>(fg, m.sources, m.destinations, m.cfg);
        (void)p->place(3);
        FlowPolicy* policy = p.get();
        return Staged {std::make_shared<std::pair<std::shared_ptr<FlowPolicy>, std::shared_ptr<FlowPolicy>>>(
                           std::move(third), std::move(p)),
                       policy};
      };
      const long injected = sweep_placement(g, where, setup, op.call);
      if (::testing::Test::HasFailure()) return;
      EXPECT_GT(injected, 0) << where;  // not vacuous
    }
  }
}

// Lossy forwarding with links that wait on each other in a cycle (contract 9.4): the rounds'
// scratch, the independent check and the error path. CYC(x, y): tunnel 0 meets X (input edge 0)
// before Y (input edge 3), tunnel 1 Y before X. With x = y = 1 the rounds settle: every
// allocation of a replace(2) fails atomically. With weights 1, 2 and y = 2 the rounds alternate:
// the period-two ValueError and every allocation failure before it publish nothing and leave no
// allocation live.
TEST(PolicyFaults, CyclicLossyForwardingFailsAtomically) {
  auto cyc = [](Amount x, Amount y) {
    return shared(make_graph(6, {{0, 1, x, 1}, {1, 2, 1000, 1}, {1, 5, 1000, 1}, {2, 3, y, 1}, {3, 0, 1000, 1},
                                 {3, 5, 1000, 1}, {4, 0, 1000, 1}, {4, 2, 1000, 1}}));
  };
  auto program = [](const Graph& g, Amount w1) {
    const auto inv = id_of_input(g);
    auto walk = [&](std::vector<std::size_t> inputs) {
      std::vector<std::vector<EdgeId>> hops;
      for (auto i : inputs) hops.push_back({inv[i]});
      return std::make_shared<const Path>(g, hops);
    };
    const NodeId s[1] = {4};
    const NodeId d[1] = {5};
    const std::vector<Route> routes = {Route {4, walk({6, 0, 1, 3, 5}), 1}, Route {4, walk({7, 3, 4, 0, 2}), w1}};
    return ForwardingProgram::from_routes(g, s, d, routes);
  };
  auto config = [](std::shared_ptr<const ForwardingProgram> prog) {
    return FlowPolicyConfig {.routing = Pinned {std::move(prog), nullptr}, .admission = Admission::Lossy};
  };
  const std::vector<NodeId> S = {4};
  const std::vector<NodeId> D = {5};
  {
    const auto g = cyc(1, 1);
    const auto prog = program(*g, 1);
    const auto setup = [&](FlowLedger& fg) {
      auto p = std::make_shared<FlowPolicy>(fg, S, D, config(prog));
      FlowPolicy* policy = p.get();
      return Staged {std::move(p), policy};
    };
    const long injected = sweep_placement(g, "converging cycle", setup, [](FlowPolicy& p) { return p.replace(2); });
    EXPECT_GT(injected, 0);  // not vacuous
  }
  const auto g = cyc(1, 2);
  const auto prog = program(*g, 2);
  long injected = 0;
  for (long at = 0;; ++at) {
    const std::string w = "period two allocation " + std::to_string(at);
    const long live = g_live.load();
    bool completed = false;
    {
      FlowLedger fg(g);
      FlowPolicy p(fg, S, D, config(prog));
      const LedgerState initial = ledger_state(fg);
      const Epoch initial_epoch = epoch(p);
      bool raised = false;
      if (fails_at(at, [&] {
            try {
              (void)p.replace(3);
            } catch (const std::invalid_argument& e) {
              raised = std::string(e.what()) == "lossy forwarding iteration entered a period-2 cycle";
            }
          })) {
        ++injected;
      } else {
        EXPECT_TRUE(raised) << w;
        completed = true;
      }
      expect_same(ledger_state(fg), initial, w);
      expect_same(epoch(p), initial_epoch, w);
    }
    EXPECT_EQ(g_live.load(), live) << w << ": a leak";
    if (completed || ::testing::Test::HasFailure()) break;
  }
  EXPECT_GT(injected, 0);  // not vacuous
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
    if (std::holds_alternative<Pinned>(m.cfg.routing)) continue;  // a program binds its endpoints
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
      if (with_cache && !std::holds_alternative<CostOnly>(m.cfg.routing)) continue;  // rejected (below)
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
  // Every CostOnly mode reads the supplied cache (an entry after a clean construction), and the
  // memo outcomes above are checked for every CostOnly mode; the other modes reject a cache before
  // binding, without touching it.
  for (const auto& m : modes(*g)) {
    FlowLedger fg(g);
    DagCache cache(*g);
    if (std::holds_alternative<CostOnly>(m.cfg.routing)) {
      const FlowPolicy p(fg, m.sources, m.destinations, m.cfg, nullptr, &cache);
      EXPECT_EQ(cache.stats().entries, 1) << m.name;
    } else {
      EXPECT_THROW(FlowPolicy(fg, m.sources, m.destinations, m.cfg, nullptr, &cache), std::invalid_argument) << m.name;
      EXPECT_EQ(cache.stats().entries, 0) << m.name;
      EXPECT_TRUE(LedgerAccess::bound_owners(fg).empty()) << m.name;
    }
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

long sweep_batch(const Mode& m, const Graph& g, std::span<const std::int32_t> workers_set, long step,
                 const Options& o = {}) {
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {m.sources, m.destinations}, {{1}, {6}}, {{0}, {5}}};
  const std::vector<const Eligibility*> per_pair(o.per_pair != nullptr ? pairs.size() : 0, o.per_pair);
  std::optional<DagCache> cache;
  const auto opts = [&] { return query_of(m, o.elig, cache ? &*cache : nullptr); };
  if (o.cache) cache.emplace(g);
  const auto want = batch_capacity(g, pairs, opts(), per_pair, 1);
  const std::int64_t clean = cache ? cache->stats().entries : 0;
  cache.reset();  // each attempt gets a fresh one (setup)
  const std::string name = m.name + " batch" + o.name;
  long kept = 0;  // failed calls whose cache kept the entries of the lookups before the failure
  long attempts = 0;
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
      attempts += failures + 1;
      EXPECT_GT(failures, 2) << name << " workers " << workers;  // not vacuous
    }
  }
  if (o.cache) EXPECT_GT(kept, 0) << name;
  return attempts;
}

long sweep_sensitivity(const Mode& m, const Graph& g, std::span<const std::int32_t> workers_set, long step,
                       const Options& o = {}) {
  std::optional<DagCache> cache;
  const auto opts = [&] { return query_of(m, o.elig, cache ? &*cache : nullptr); };
  if (o.cache) cache.emplace(g);
  const auto want = sensitivity(g, m.sources, m.destinations, opts(), 1);
  const std::int64_t clean = cache ? cache->stats().entries : 0;
  cache.reset();  // each attempt gets a fresh one (setup)
  const std::string name = m.name + " sensitivity" + o.name;
  long kept = 0;  // failed calls whose cache kept the entries of the lookups before the failure
  long attempts = 0;
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
      attempts += failures + 1;
      EXPECT_GT(failures, 2) << name << " workers " << workers;
    }
  }
  if (o.cache) EXPECT_GT(kept, 0) << name;
  return attempts;
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
  const std::vector<std::string> exhaustive_sensitivity = {"costonly eb", "global", "pinned eb split routes",
                                                           "pinned tunnels lossy"};
  int swept_sensitivity = 0;
  long attempts = 0;
  for (const auto& m : modes(*g)) {
    const bool reads_cache = std::holds_alternative<CostOnly>(m.cfg.routing);
    const bool every_sensitivity_allocation =
        std::find(exhaustive_sensitivity.begin(), exhaustive_sensitivity.end(), m.name) != exhaustive_sensitivity.end();
    for (const auto& o : option_sets(no_bypass, no_cross)) {
      if (o.cache && !reads_cache) continue;
      const std::string name = m.name + " capacity" + o.name;
      std::optional<DagCache> cache;
      const auto copts = [&] {
        return query_of(m, o.elig, cache ? &*cache : nullptr, OutputOptions {true, true, true});
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
        ++attempts;
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
      if (!std::holds_alternative<Pinned>(m.cfg.routing)) attempts += sweep_batch(m, *g, one, 1, o);  // batch rejects Pinned
      if (every_sensitivity_allocation) {
        attempts += sweep_sensitivity(m, *g, one, 1, o);
        ++swept_sensitivity;
      } else if (o.elig == nullptr && !o.cache) {
        attempts += sweep_sensitivity(m, *g, one, 11, o);
      }
      if (::testing::Test::HasFailure()) return;
    }
  }
  RecordProperty("allocation_attempts", std::to_string(attempts));
  EXPECT_EQ(swept_sensitivity, 3 * 2 + 3);  // three modes with two option sets, CostOnly's with three
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
  long attempts = 0;
  for (const auto& m : modes(*g)) {
    attempts += sweep_sensitivity(m, *g, several, 11);
    if (::testing::Test::HasFailure()) return;
  }
  RecordProperty("allocation_attempts", std::to_string(attempts));
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
// small size; the failed call changes nothing, and a clean retry on the same thread and graph (or
// the call itself, once no failure is reached) gives the reference call's result, every output
// included, ledger and epoch (ASan reports any access outside a partly grown array). The same call
// on a new thread whose scratch has first grown to the large graph allocates less: the sweep
// reached the growth. `config` gives the policy's configuration on a graph.
void expect_scratch_growth_failure_leaves_nothing_half_grown(
    const std::string& name, const std::function<FlowPolicyConfig(const Graph&)>& config) {
  SCOPED_TRACE(name);
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
    FlowPolicy p(fg, std::vector<NodeId> {0}, std::vector<NodeId> {n - 1}, config(*large));
    want = p.place(2, all);
    want_state = ledger_state(fg);
    want_epoch = epoch(p);
  }
  ASSERT_EQ(want.placed, 2);  // the whole chain carries the request: the comparison is not of zeros
  sweep_scratch_growth(
      name + " placement",
      [&] {
        FlowLedger fg(small);
        FlowPolicy p(fg, std::vector<NodeId> {0}, std::vector<NodeId> {2}, config(*small));
        EXPECT_EQ(p.place(2, all).placed, 2);
      },
      [&](long at, const std::string& w) {
        FlowLedger fg(large);
        FlowPolicy p(fg, std::vector<NodeId> {0}, std::vector<NodeId> {n - 1}, config(*large));
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
        FlowPolicy p(fg, std::vector<NodeId> {0}, std::vector<NodeId> {n - 1}, config(*large));
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

// A split kernel's placement (CostOnly EqualBalanced): the working residual and the kernels'
// scratch.
TEST(PolicyFaults, ScratchGrowthFailureLeavesNothingHalfGrown) {
  expect_scratch_growth_failure_leaves_nothing_half_grown("costonly eb",
                                                          [](const Graph&) { return cost_only_config(); });
}

// A controller's placement also grows the decomposition's per-thread scratch (cursors, range ends,
// destination and slot marks, each grown on its own): Global over every edge (positions are
// EdgeIds) and Global on a support (positions index the support).
TEST(PolicyFaults, ControllerScratchGrowthFailureLeavesNothingHalfGrown) {
  expect_scratch_growth_failure_leaves_nothing_half_grown("global", [](const Graph&) { return global_config(); });
  if (::testing::Test::HasFailure()) return;
  expect_scratch_growth_failure_leaves_nothing_half_grown("global support", [](const Graph& g) {
    return FlowPolicyConfig {.routing = Global {all_edges(g)}};
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
  // policy constructor and capacity in every CostOnly mode (the other models reject a cache before
  // reading it, statistics unchanged); batch_capacity at 1, 2 and 4 workers, in any
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
    if (!std::holds_alternative<CostOnly>(m.cfg.routing)) {
      const auto before = cache.stats();
      FlowLedger fg(g);
      EXPECT_THROW(FlowPolicy(fg, m.sources, m.destinations, m.cfg, nullptr, &cache), std::invalid_argument) << m.name;
      EXPECT_THROW((void)capacity(*g, m.sources, m.destinations, query_of(m, nullptr, &cache)),
                   std::invalid_argument)
          << m.name;
      expect_stats(cache.stats(), before, m.name + " rejected");
      continue;
    }
    const long built = sweep_frozen_stats(m.name + " frozen construction", cache, false, 1, reads(m), false, [&] {
      FlowLedger fg(g);
      FlowPolicy p(fg, m.sources, m.destinations, m.cfg, nullptr, &cache);
    });
    EXPECT_GT(built, 2) << m.name;
    const QueryOptions opts {m.cfg, nullptr, &cache, OutputOptions {true, true, true}};
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
  const QueryOptions shared {eb.cfg, nullptr, &cache};
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
                                                     QueryOptions {eb.cfg, nullptr, &cache}, workers);
                                 }),
              2);
  }
  // CostOnly sensitivity whose deletions miss (a fresh uid per deletion): the baseline's hit is
  // dropped with the failed call, and the miss counts once.
  const auto costonly_miss = [&](std::int32_t workers) {
    (void)sensitivity(*g, eb.sources, eb.destinations, QueryOptions {eb.cfg, nullptr, &cache}, workers);
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
  // statistics with two and four workers: the test above): the failure propagates as
  // std::bad_alloc after the join, the cache keeps its entries, and a clean retry gives the
  // reference. Global sensitivity, which reads no cache, rejects it before any worker starts.
  const auto g = graph();
  const auto all = modes(*g);
  const auto& m = all[1];  // costonly eb: reads the cache
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {m.sources, m.destinations}, {{1}, {6}}, {{0}, {5}}};
  DagCache cache(*g);
  for (const auto& [s, t] : pairs) (void)cache.lookup(*g, t, Direction::Reverse, nullptr);
  const auto entries = cache.stats().entries;
  cache.freeze();
  const auto want = batch_capacity(*g, pairs, query_of(m), {}, 1);
  const QueryOptions shared {m.cfg, nullptr, &cache};
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
    (void)sensitivity(*g, m.sources, m.destinations, query_of(m, nullptr, &cache), w);
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
  // Global sensitivity reads no cache: with the frozen cache it is rejected before any worker
  // starts, the cache unchanged.
  const auto git = std::find_if(all.begin(), all.end(), [](const Mode& x) { return x.name == "global"; });
  ASSERT_NE(git, all.end());
  EXPECT_THROW((void)sensitivity(*g, git->sources, git->destinations, QueryOptions {git->cfg, nullptr, &cache}, 2),
               std::invalid_argument);
  EXPECT_EQ(cache.stats().entries, entries);
}

// ---------------------------------------------------------------- T1

// T1's calls (replace, evaluate, replace(0), replace(+inf), forward with and without black holes,
// and two whose holds change: a replace with a new design and the first replace of a policy)
// with a third party's load and another T1 owner's holds on the ledger: every allocation of the
// call fails once, through the shared placement sweep, which compares both projections and every
// record field (reserved_edges included): a failure changes neither the ledger nor the epoch and
// leaves no allocation live; the completing call, the retry and remove() leave the ledger as the
// reference's do.
TEST(PolicyFaults, T1CallsFailAtomically) {
  const auto g = graph();
  const auto inv = id_of_input(*g);
  const EdgeId dark_edges[] = {inv[4], inv[6]};
  const Eligibility dark = g->eligibility({}, dark_edges);
  const OutputOptions all {.edge_flows = true, .residual = true, .bottleneck = true, .reservation_edge_flows = true};
  struct Op {
    std::string name;
    std::function<PlacementResult(FlowPolicy&)> call;
    bool fresh {false};  // on a policy without a call: every hold is new
  };
  for (const auto& m : t1_modes()) {
    const bool lossy = m.cfg.admission == Admission::Lossy;
    const auto design = *m.design;
    std::vector<Op> ops = {{"evaluate", [&](FlowPolicy& p) { return p.evaluate(6, design, all); }},
                           {"replace", [&](FlowPolicy& p) { return p.replace(6, design, all); }},
                           {"replace to zero", [&](FlowPolicy& p) { return p.replace(0, design); }},
                           {"forward", [](FlowPolicy& p) { return p.forward(6, nullptr, {true, true, true}); }},
                           {"forward black holes", [&](FlowPolicy& p) { return p.forward(6, &dark); }},
                           // The holds change: a failed commit must leave the old ones.
                           {"replace with a new design",
                            [&](FlowPolicy& p) { return p.replace(6, std::optional<Amount> {*design + 2}, all); }},
                           {"first replace", [&](FlowPolicy& p) { return p.replace(6, design, all); }, true}};
    if (!lossy) ops.push_back({"replace +inf", [&](FlowPolicy& p) { return p.replace(std::nullopt, design); }});
    for (const auto& op : ops) {
      const std::string where = m.name + " / " + op.name;
      const auto setup = [&](FlowLedger& fg) {
        auto third = std::shared_ptr<FlowPolicy>(third_party(fg, *g));
        auto holder = std::make_shared<FlowPolicy>(fg, std::vector<NodeId> {1}, std::vector<NodeId> {6},
                                                   FlowPolicyConfig {.routing = CapacityAware {Tunnels {1}}});
        (void)holder->replace(0, 2);
        auto p = std::make_shared<FlowPolicy>(fg, m.sources, m.destinations, m.cfg);
        if (!op.fresh) (void)p->replace(3, design);
        FlowPolicy* policy = p.get();
        auto keep = std::make_shared<std::vector<std::shared_ptr<FlowPolicy>>>(
            std::vector<std::shared_ptr<FlowPolicy>> {std::move(third), std::move(holder), std::move(p)});
        return Staged {keep, policy};
      };
      const long injected = sweep_placement(g, where, setup, op.call);
      if (::testing::Test::HasFailure()) return;
      EXPECT_GT(injected, 0) << where;  // not vacuous
    }
  }
}

// A T1 policy's construction (the roster, the pool baseline, the replay basis of its exclusions)
// fails at every allocation without binding; a clean construction then sets up like the
// reference.
TEST(PolicyFaults, T1ConstructionFailsWithoutBinding) {
  const auto g = graph();
  const Eligibility no_bypass = without_bypass(*g);
  for (const auto& m : t1_modes()) {
    for (const Eligibility* elig : {static_cast<const Eligibility*>(nullptr), &no_bypass}) {
      const std::string variant = m.name + (elig != nullptr ? " eligibility" : "");
      FlowLedger rf(g);
      FlowPolicy pr(rf, m.sources, m.destinations, m.cfg, elig);
      const PlacementResult want = pr.replace(2, *m.design);
      const LedgerState placed = ledger_state(rf);
      long injected = 0;
      for (long at = 0;; ++at) {
        const std::string w = variant + " construction allocation " + std::to_string(at);
        const long live = g_live.load();
        {
          FlowLedger fg(g);
          const LedgerState before = ledger_state(fg);
          std::optional<FlowPolicy> p;
          if (!fails_at(at, [&] { p.emplace(fg, m.sources, m.destinations, m.cfg, elig); })) {
            EXPECT_EQ(LedgerAccess::bound_owners(fg).size(), 1U);
            break;
          }
          ++injected;
          EXPECT_FALSE(p.has_value()) << w;
          expect_same(ledger_state(fg), before, w);
          p.emplace(fg, m.sources, m.destinations, m.cfg, elig);
          expect_same(p->replace(2, *m.design), want, w);
          expect_same(ledger_state(fg), placed, w);
          expect_balanced(fg, w);
        }
        EXPECT_EQ(g_live.load(), live) << w << ": a leak";
        if (::testing::Test::HasFailure()) return;
      }
      EXPECT_GT(injected, 0) << variant;
    }
  }
}

// T1's evaluators (capacity: setup and forwarding on a fresh ledger; batch_capacity; sensitivity:
// setup redone under every deletion) under the shared evaluator sweeps, one worker exhaustively
// and several workers in their threads.
TEST(PolicyFaults, T1EvaluatorsLeaveNoPoisonedScratch) {
  const auto g = graph();
  const Eligibility no_bypass = without_bypass(*g);
  const Eligibility no_cross = without_cross(*g);
  const std::int32_t one[] = {1};
  const std::int32_t several[] = {2, 3};
  long attempts = 0;
  for (const auto& m : t1_modes()) {
    for (const auto& o : option_sets(no_bypass, no_cross)) {
      if (o.cache) continue;  // T1 rejects a cache
      const std::string name = m.name + " capacity" + o.name;
      const auto want = capacity(*g, m.sources, m.destinations, query_of(m, o.elig, nullptr, {true, true, true}));
      long injected = 0;
      for (long at = 0;; ++at) {
        const long live = g_live.load();
        bool failed = false;
        {
          PlacementResult got;
          failed = fails_at(at, [&] { got = capacity(*g, m.sources, m.destinations, query_of(m, o.elig, nullptr, {true, true, true})); });
          if (failed) {
            ++injected;
            expect_same(capacity(*g, m.sources, m.destinations, query_of(m, o.elig, nullptr, {true, true, true})), want,
                        name + " retry");
          } else {
            expect_same(got, want, name);
          }
        }
        EXPECT_EQ(g_live.load(), live) << name << " allocation " << at << ": a leak";
        if (!failed || ::testing::Test::HasFailure()) break;
      }
      EXPECT_GT(injected, 0) << name;
      attempts += injected + sweep_batch(m, *g, one, 1, o) + sweep_sensitivity(m, *g, one, 1, o);
      attempts += sweep_sensitivity(m, *g, several, 7, o);
      if (::testing::Test::HasFailure()) return;
    }
  }
  RecordProperty("allocation_attempts", std::to_string(attempts));
}

// remove() and the destruction of a T1 policy allocate nothing: remove() takes the records and the
// holds away; the destruction leaves both projections as orphans and the policy no longer bound.
TEST(PolicyFaults, T1RemovalAndReleaseNeverAllocate) {
  const auto g = graph();
  for (const auto& m : t1_modes()) {
    FlowLedger fg(g);
    auto p = std::make_unique<FlowPolicy>(fg, m.sources, m.destinations, m.cfg);
    (void)p->replace(4, *m.design);
    ASSERT_FALSE(fg.records().empty()) << m.name;
    EXPECT_FALSE(fails_at(0, [&] { p->remove(); })) << m.name;
    EXPECT_TRUE(fg.records().empty()) << m.name;
    EXPECT_EQ(fg.reserved(), std::vector<Amount>(fg.reserved().size(), 0)) << m.name;
    (void)p->replace(2, *m.design);
    const LedgerState held = ledger_state(fg);
    ASSERT_GT(std::accumulate(held.reserved.begin(), held.reserved.end(), Amount {0}), 0) << m.name;
    EXPECT_FALSE(fails_at(0, [&] { p.reset(); })) << m.name;
    const LedgerState orphaned = ledger_state(fg);
    EXPECT_EQ(orphaned.allocated, held.allocated) << m.name;
    EXPECT_EQ(orphaned.reserved, held.reserved) << m.name;
    EXPECT_EQ(orphaned.records.size(), held.records.size()) << m.name;
    EXPECT_TRUE(orphaned.bound.empty()) << m.name;
  }
}

// T1's commit books the holds through the reservation projection's own per-thread scratch, so the
// shared placement sweep (whose reference call has already grown this thread's scratch) does not
// reach its growth. Here each attempt runs on a new thread whose scratch is cold, or has first
// served a T1 replace on a 3-node chain (small to large); the swept replace on a 3000-node chain
// fails at its at-th allocation and changes neither projection, no record and not the epoch, and a
// retry on the same thread gives the reference's result, ledger and epoch. With `prior`, the policy
// already holds an earlier design (set up on another thread, so the attempt's scratch stays as it
// was) and the swept call replaces those holds. The same call on a thread whose scratch has grown to
// the large chain allocates less: the sweep reached the growth.
void expect_t1_hold_scratch_growth(const std::string& name, bool from_small, bool prior) {
  SCOPED_TRACE(name);
  const auto small = chain_graph(3, 5);
  const auto large = chain_graph(3000, 5);
  const NodeId last = large->num_nodes() - 1;
  const FlowPolicyConfig cfg {.routing = CapacityAware {Tunnels {2}}};
  const auto make = [&](FlowLedger& fg) {
    auto p = std::make_unique<FlowPolicy>(fg, std::vector<NodeId> {0}, std::vector<NodeId> {last}, cfg);
    if (prior) std::thread([&] { (void)p->replace(2, 2); }).join();
    return p;
  };
  const OutputOptions holds {.reservation_edge_flows = true};
  const auto call = [&](FlowPolicy& p) { return p.replace(4, 4, holds); };
  PlacementResult want;
  LedgerState want_state;
  Epoch want_epoch;
  {
    FlowLedger fg(large);
    auto p = make(fg);
    want = call(*p);
    want_state = ledger_state(fg);
    want_epoch = epoch(*p);
  }
  // Both tunnels up, each holding 2 on every link: the comparison is not of zeros.
  ASSERT_EQ(want_state.reserved, std::vector<Amount>(want_state.reserved.size(), 4));
  ASSERT_EQ(want.placed, 4);
  sweep_scratch_growth(
      name,
      [&] {
        if (!from_small) return;
        FlowLedger fg(small);
        FlowPolicy p(fg, std::vector<NodeId> {0}, std::vector<NodeId> {2}, cfg);
        EXPECT_EQ(p.replace(4, 4).placed, 4);
      },
      [&](long at, const std::string& w) {
        FlowLedger fg(large);
        auto p = make(fg);
        const LedgerState before = ledger_state(fg);
        const Epoch before_epoch = epoch(*p);
        PlacementResult got;
        const bool failed = fails_at(at, [&] { got = call(*p); });
        if (failed) {
          expect_same(ledger_state(fg), before, w + " (failed call)");
          expect_same(epoch(*p), before_epoch, w + " (failed call)");
          expect_balanced(fg, w + " (failed call)");
          got = call(*p);
        }
        expect_same(got, want, w);
        expect_same(ledger_state(fg), want_state, w);
        expect_same(epoch(*p), want_epoch, w);
        expect_balanced(fg, w);
        return failed;
      },
      [&] {
        {
          FlowLedger fg(large);
          auto p = make(fg);
          (void)call(*p);  // the scratch grows to the large chain
        }
        for (long warm = 0;; ++warm) {
          FlowLedger fg(large);
          auto p = make(fg);
          if (!fails_at(warm, [&] { (void)call(*p); })) return warm;
        }
      });
}

TEST(PolicyFaults, T1HoldScratchGrowthFailureLeavesNothingHalfGrown) {
  for (const bool prior : {false, true}) {
    for (const bool from_small : {false, true}) {
      expect_t1_hold_scratch_growth(std::string(prior ? "new design" : "first replace") +
                                        (from_small ? " from small" : " from cold"),
                                    from_small, prior);
      if (::testing::Test::HasFailure()) return;
    }
  }
}
