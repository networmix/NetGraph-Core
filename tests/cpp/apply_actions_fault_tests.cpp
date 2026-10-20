// apply_actions under allocation failures (the fault-injection executable, fault_injection.cpp):
// every allocation of an ORDERED list mixing every model and of a SIMULTANEOUS list (a pooled
// ShortestPaths replace, a Pinned tunnel program, a T1 forward with holds and black holes, a cyclic
// group across demands) fails once in turn; after each failure the ledger and every listed policy's
// epoch are as before the call, the aggregate is balanced, nothing is left allocated, and a clean
// retry gives the reference's results, ledger and epochs. Because the sweep reaches every allocation
// of the call, a publication that allocated after the ledger's exchange would show as a changed
// state here: the publication allocates nothing. The batch exchange is swept on its own (owners that
// retire and keep records, holds), and both modes are swept on a fresh thread whose scratch grows
// from a 3-node to a 3,000-node graph, the retry on that thread equal to the reference.
#include <gtest/gtest.h>

#include <algorithm>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/forwarding_program.hpp"
#include "fault_injection.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using detail::LedgerAccess;

namespace {

using Actions = std::vector<Action>;

const OutputOptions kAll {.edge_flows = true, .residual = true, .saturated_boundary_edges = true, .edge_reservations = true};

// What a list sweep stages on a fresh ledger: the policies (kept alive), their earlier calls, the
// list and its mode.
struct ListScenario {
  std::vector<std::shared_ptr<const Eligibility>> masks;
  std::vector<std::unique_ptr<FlowPolicy>> policies;
  Actions actions;
  BatchMode mode {BatchMode::Ordered};
};

EdgeId edge(const Graph& g, NodeId u, NodeId v) {
  for (EdgeId e = 0; e < g.num_edges(); ++e) {
    if (g.edge_sources()[static_cast<std::size_t>(e)] == u && g.edge_destinations()[static_cast<std::size_t>(e)] == v) return e;
  }
  ADD_FAILURE() << "no edge " << u << " -> " << v;
  return 0;
}

std::shared_ptr<const ForwardingProgram> tunnels(const Graph& g, const std::vector<NodeId>& s, const std::vector<NodeId>& d,
                                                 const std::vector<std::pair<NodeId, std::vector<NodeId>>>& walks) {
  std::vector<RouteDefinition> routes;
  for (const auto& [h, nodes] : walks) {
    std::vector<std::vector<EdgeId>> hops;
    for (std::size_t i = 0; i + 1 < nodes.size(); ++i) hops.push_back({edge(g, nodes[i], nodes[i + 1])});
    routes.push_back(RouteDefinition {h, std::make_shared<const Path>(g, hops), 1});
  }
  return ForwardingProgram::from_routes(g, s, d, routes);
}

std::vector<Epoch> epochs(const ListScenario& s) {
  std::vector<Epoch> out;
  for (const auto& p : s.policies) out.push_back(epoch(*p));
  return out;
}

void expect_same_epochs(const std::vector<Epoch>& a, const std::vector<Epoch>& b, const std::string& where) {
  ASSERT_EQ(a.size(), b.size()) << where;
  for (std::size_t i = 0; i < a.size(); ++i) expect_same(a[i], b[i], where + ", policy " + std::to_string(i));
}

// The list sweep (as sweep_placement for one call): returns the failures injected.
long sweep_list(const std::shared_ptr<const Graph>& g, const std::string& where, const std::function<ListScenario(FlowLedger&)>& setup) {
  FlowLedger ref(g);
  ListScenario sref = setup(ref);
  const LedgerState initial = ledger_state(ref);
  const auto initial_epochs = epochs(sref);
  const auto want = apply_actions(ref, sref.actions, sref.mode, kAll);
  const LedgerState after = ledger_state(ref);
  const auto want_epochs = epochs(sref);
  long injected = 0;
  for (long at = 0;; ++at) {
    const std::string w = where + " allocation " + std::to_string(at);
    const long live = g_live.load();
    bool completed = false;
    {
      FlowLedger fg(g);
      ListScenario s = setup(fg);
      expect_same(ledger_state(fg), initial, w + " (setup)");
      expect_same_epochs(epochs(s), initial_epochs, w + " (setup)");
      std::vector<PlacementResult> got;
      if (!fails_at(at, [&] { got = apply_actions(fg, s.actions, s.mode, kAll); })) {
        EXPECT_EQ(got.size(), want.size()) << w;
        for (std::size_t k = 0; k < std::min(got.size(), want.size()); ++k) expect_same(got[k], want[k], w + " entry " + std::to_string(k));
        expect_same(ledger_state(fg), after, w + " (no failure)");
        expect_same_epochs(epochs(s), want_epochs, w + " (no failure)");
        completed = true;
      } else {
        ++injected;
        expect_same(ledger_state(fg), initial, w + " (failed list)");
        expect_same_epochs(epochs(s), initial_epochs, w + " (failed list)");
        expect_balanced(fg, w + " (failed list)");
        const auto retry = apply_actions(fg, s.actions, s.mode, kAll);
        for (std::size_t k = 0; k < std::min(retry.size(), want.size()); ++k) expect_same(retry[k], want[k], w + " retry " + std::to_string(k));
        expect_same(ledger_state(fg), after, w + " (retry)");
        expect_same_epochs(epochs(s), want_epochs, w + " (retry)");
      }
    }
    EXPECT_EQ(g_live.load(), live) << w << ": a leak";
    if (completed || ::testing::Test::HasFailure()) break;
  }
  return injected;
}

// Nodes 0..5 (sources 0 and 1, destination 5) and a separate link 6 -> 7.
std::shared_ptr<const Graph> mixed_graph() {
  return shared(make_graph(8, {{0, 2, 6, 1}, {0, 3, 4, 2}, {1, 2, 5, 1}, {1, 3, 3, 1}, {2, 4, 5, 1}, {3, 4, 6, 2}, {2, 5, 4, 3},
                               {4, 5, 8, 1}, {3, 5, 2, 1}, {6, 7, 10, 1}}));
}

const std::vector<NodeId> kS {0, 1};
const std::vector<NodeId> kD {5};

FlowPolicyConfig ecmp(Admission a = Admission::Lossy) { return FlowPolicyConfig {.routing = shortest_paths(), .admission = a}; }
FlowPolicyConfig pinned_config(std::shared_ptr<const ForwardingProgram> prog, std::shared_ptr<const Eligibility> black = nullptr) {
  return FlowPolicyConfig {.routing = Pinned {std::move(prog), std::move(black)}, .admission = Admission::Lossy};
}

// An ORDERED list over every model: a ShortestPaths place, a LOSSLESS replace(+inf), a
// controller's place, T1 and T2 replaces, a Pinned replace, a forward with black holes of an
// installed program, and a zero place.
TEST(PolicyFaults, EveryAllocationOfAnOrderedListFailsAtomically) {
  const auto g = mixed_graph();
  const auto prog = tunnels(*g, kS, kD, {{0, {0, 2, 4, 5}}, {1, {1, 3, 5}}});
  const auto injected = sweep_list(g, "ordered", [&](FlowLedger& fg) {
    ListScenario s;
    auto add = [&](FlowPolicyConfig cfg, std::vector<NodeId> src = kS) -> FlowPolicy& {
      s.policies.push_back(std::make_unique<FlowPolicy>(fg, src, kD, std::move(cfg)));
      return *s.policies.back();
    };
    auto& sp = add(ecmp());
    auto& ll = add(ecmp(Admission::Lossless), {0});
    auto& ctl = add(FlowPolicyConfig {.routing = centralized()});
    auto& t1 = add(FlowPolicyConfig {.routing = DistributedTE {TunnelCount {1}}, .admission = Admission::Lossy});
    auto& t2 = add(FlowPolicyConfig {.routing = centralized(SourceAllocation::Equal, nullptr, TunnelCount {1}), .admission = Admission::Lossy});
    auto& pin = add(pinned_config(prog));
    auto& fw = add(ecmp());
    auto& zero = add(ecmp());
    (void)sp.place(2);
    (void)fw.replace(2);
    (void)t1.replace(2, std::optional<Amount> {2});
    const std::vector<EdgeId> black {edge(*g, 2, 4)};
    s.masks.push_back(std::make_shared<const Eligibility>(g->eligibility({}, black)));
    s.actions = {Place {sp, 4},
                 Replace {ll, std::nullopt},
                 Place {ctl, 3},
                 Replace {t1, 3, std::optional<Amount> {3}},
                 Replace {t2, 2},
                 Replace {pin, 4},
                 Forward {fw, 3, s.masks.back().get()},
                 Place {zero, 0}};
    s.mode = BatchMode::Ordered;
    return s;
  });
  EXPECT_GT(injected, 100);
}

// A SIMULTANEOUS list: a pooled ShortestPaths replace, a Pinned tunnel program, a T1 forward with
// holds and black holes, and (on its own graph) a cyclic group across two demands.
TEST(PolicyFaults, EveryAllocationOfASimultaneousListFailsAtomically) {
  const auto g = mixed_graph();
  const auto prog = tunnels(*g, kS, kD, {{0, {0, 2, 5}}, {1, {1, 2, 4, 5}}});
  const auto injected = sweep_list(g, "simultaneous", [&](FlowLedger& fg) {
    ListScenario s;
    auto add = [&](FlowPolicyConfig cfg) -> FlowPolicy& {
      s.policies.push_back(std::make_unique<FlowPolicy>(fg, kS, kD, std::move(cfg)));
      return *s.policies.back();
    };
    auto& sp = add(ecmp());
    auto& pin = add(pinned_config(prog));
    auto& t1 = add(FlowPolicyConfig {.routing = DistributedTE {TunnelCount {1}}, .admission = Admission::Lossy});
    (void)sp.place(3);
    (void)t1.replace(3, std::optional<Amount> {3});
    const std::vector<EdgeId> black {edge(*g, 3, 5)};
    s.masks.push_back(std::make_shared<const Eligibility>(g->eligibility({}, black)));
    s.actions = {Replace {sp, 5}, Replace {pin, 6}, Forward {t1, 4, s.masks.back().get()}};
    s.mode = BatchMode::Simultaneous;
    return s;
  });
  EXPECT_GT(injected, 100);
  // The cyclic group: two demands crossing X (e0) and Y (e3) in opposite orders, settling at round 3.
  const auto cyc = shared(make_graph(6, {{0, 1, 1, 1}, {1, 2, 1000, 1}, {1, 5, 1000, 1}, {2, 3, 1, 1}, {3, 0, 1000, 1},
                                         {3, 5, 1000, 1}, {4, 0, 1000, 1}, {4, 2, 1000, 1}}));
  const auto a = tunnels(*cyc, {4}, {5}, {{4, {4, 0, 1, 2, 3, 5}}});
  const auto b = tunnels(*cyc, {4}, {5}, {{4, {4, 2, 3, 0, 1, 5}}});
  const auto cyclic = sweep_list(cyc, "cyclic", [&](FlowLedger& fg) {
    ListScenario s;
    s.policies.push_back(std::make_unique<FlowPolicy>(fg, std::vector<NodeId> {4}, std::vector<NodeId> {5}, pinned_config(a)));
    s.policies.push_back(std::make_unique<FlowPolicy>(fg, std::vector<NodeId> {4}, std::vector<NodeId> {5}, pinned_config(b)));
    s.actions = {Replace {*s.policies[0], 1}, Replace {*s.policies[1], 1}};
    s.mode = BatchMode::Simultaneous;
    return s;
  });
  EXPECT_GT(cyclic, 50);
}

// Both modes on a fresh thread, the scratch grown from a 3-node to a 3,000-node chain, each
// allocation failing once; the retry on that thread equals the reference.
TEST(PolicyFaults, ListScratchGrowthFailureLeavesNothingHalfGrown) {
  for (BatchMode mode : {BatchMode::Ordered, BatchMode::Simultaneous}) {
    const std::string name = mode == BatchMode::Ordered ? "ordered" : "simultaneous";
    SCOPED_TRACE(name);
    const auto small = chain_graph(3, 5);
    const auto large = chain_graph(3000, 5);
    auto chain = [](const Graph& g) {
      std::vector<std::vector<EdgeId>> hops;
      for (EdgeId e = 0; e < g.num_edges(); ++e) hops.push_back({e});
      const std::vector<RouteDefinition> routes {RouteDefinition {0, std::make_shared<const Path>(g, hops), 1}};
      return ForwardingProgram::from_routes(g, std::vector<NodeId> {0}, std::vector<NodeId> {g.num_nodes() - 1}, routes);
    };
    struct Run {
      std::vector<PlacementResult> results;
      LedgerState state;
    };
    auto list_on = [&](const std::shared_ptr<const Graph>& g, FlowLedger& fg, std::unique_ptr<FlowPolicy>& a,
                       std::unique_ptr<FlowPolicy>& b) {
      const std::vector<NodeId> s {0};
      const std::vector<NodeId> d {g->num_nodes() - 1};
      a = std::make_unique<FlowPolicy>(fg, s, d, ecmp());
      b = std::make_unique<FlowPolicy>(fg, s, d, pinned_config(chain(*g)));
      return Actions {Replace {*a, 2}, Replace {*b, 1}};
    };
    Run want;
    {
      FlowLedger fg(large);
      std::unique_ptr<FlowPolicy> a;
      std::unique_ptr<FlowPolicy> b;
      const auto list = list_on(large, fg, a, b);
      want.results = apply_actions(fg, list, mode, kAll);
      want.state = ledger_state(fg);
    }
    ASSERT_EQ(want.results[0].call_delivered, 2);
    sweep_scratch_growth(
        name,
        [&] {
          FlowLedger fg(small);
          std::unique_ptr<FlowPolicy> a;
          std::unique_ptr<FlowPolicy> b;
          const auto list = list_on(small, fg, a, b);
          EXPECT_EQ(apply_actions(fg, list, mode, kAll).size(), 2U);
        },
        [&](long at, const std::string& w) {
          FlowLedger fg(large);
          std::unique_ptr<FlowPolicy> a;
          std::unique_ptr<FlowPolicy> b;
          const auto list = list_on(large, fg, a, b);
          const LedgerState before = ledger_state(fg);
          std::vector<PlacementResult> got;
          const bool failed = fails_at(at, [&] { got = apply_actions(fg, list, mode, kAll); });
          if (failed) {
            expect_same(ledger_state(fg), before, w);
            got = apply_actions(fg, list, mode, kAll);
          }
          for (std::size_t k = 0; k < got.size(); ++k) expect_same(got[k], want.results[k], w);
          expect_same(ledger_state(fg), want.state, w);
          return failed;
        },
        [&] {
          FlowLedger fg(large);
          std::unique_ptr<FlowPolicy> a;
          std::unique_ptr<FlowPolicy> b;
          const auto list = list_on(large, fg, a, b);
          (void)apply_actions(fg, list, mode, kAll);  // the scratch grows to the large graph
          long warm = 0;
          while (fails_at(warm, [&] { (void)apply_actions(fg, list, mode, kAll); })) ++warm;
          return warm;
        });
    if (::testing::Test::HasFailure()) return;
  }
}

// The batch exchange: every allocation of a two-owner and a three-owner commit_batch (an owner
// that retires its records, one that keeps them, holds) fails once; the ledger is unchanged after
// each failure, and the completing exchange equals the reference.
TEST(LedgerFaults, EveryAllocationOfABatchExchangeFailsAtomically) {
  const auto g = shared(make_graph(3, {{0, 1, 10, 1}, {1, 2, 10, 1}}));
  const auto dag = all_edges(*g);
  auto draft = [&](std::vector<std::pair<EdgeId, Amount>> flows, std::optional<std::int64_t> id = std::nullopt, Amount hold = 0) {
    detail::RecordDraft d;
    d.record_counter = id;
    d.record.sources = {0};
    d.record.destinations = {2};
    d.record.route = dag;
    d.record.min_route_cost = 2;
    for (const auto& [e, x] : flows) d.record.delivered = std::max(d.record.delivered, x);
    d.record.edge_flows = std::move(flows);
    if (hold != 0) d.record.edge_reservations = {{0, hold}};
    return d;
  };
  for (int owners : {2, 3}) {
    const std::string where = std::to_string(owners) + " owners";
    struct Staged {
      std::vector<std::int64_t> ids;
      std::vector<std::vector<detail::RecordDraft>> drafts;
      std::vector<bool> retire;
    };
    auto stage = [&](FlowLedger& fg) {
      Staged s;
      for (int i = 0; i < owners; ++i) s.ids.push_back(LedgerAccess::bind(fg));
      std::vector<detail::RecordDraft> first {draft({{0, 2}, {1, 2}}, std::nullopt, 3)};
      const auto kept = LedgerAccess::commit(fg, s.ids[0], first);
      std::vector<detail::RecordDraft> second {draft({{0, 1}})};
      (void)LedgerAccess::commit(fg, s.ids[1], second);
      s.drafts.push_back({draft({{0, 3}, {1, 3}}, kept[0].second, 2), draft({{1, 1}})});  // keeps and adds
      s.drafts.push_back({draft({{0, 2}}, std::nullopt, 4)});                                // retires and adds a hold
      s.retire = {false, true};
      if (owners == 3) {
        s.drafts.push_back({draft({{0, 1}, {1, 1}}), draft({{1, 2}})});
        s.retire.push_back(true);
      }
      return s;
    };
    auto commit = [&](FlowLedger& fg, Staged& s) {
      std::vector<detail::OwnerDrafts> parts;
      for (std::size_t i = 0; i < s.ids.size(); ++i) parts.push_back({s.ids[i], s.drafts[i], s.retire[i]});
      std::vector<std::vector<FlowId>> created(s.ids.size());
      LedgerAccess::commit_batch(fg, parts, created);
      return created;
    };
    FlowLedger ref(g);
    Staged sref = stage(ref);
    const LedgerState initial = ledger_state(ref);
    const auto want_ids = commit(ref, sref);
    const LedgerState after = ledger_state(ref);
    long injected = 0;
    for (long at = 0;; ++at) {
      const std::string w = where + " allocation " + std::to_string(at);
      const long live = g_live.load();
      bool completed = false;
      {
        FlowLedger fg(g);
        Staged s = stage(fg);
        expect_same(ledger_state(fg), initial, w + " (setup)");
        std::vector<std::vector<FlowId>> got;
        if (!fails_at(at, [&] { got = commit(fg, s); })) {
          EXPECT_EQ(got.size(), want_ids.size()) << w;
          expect_same(ledger_state(fg), after, w + " (no failure)");
          completed = true;
        } else {
          ++injected;
          expect_same(ledger_state(fg), initial, w + " (failed)");
          expect_balanced(fg, w + " (failed)");
        }
      }
      EXPECT_EQ(g_live.load(), live) << w << ": a leak";
      if (completed || ::testing::Test::HasFailure()) break;
    }
    EXPECT_GT(injected, 5) << where;
  }
}

} // namespace
