#include <gtest/gtest.h>

#include "fault_injection.hpp"
#include "netgraph/core/shortest_paths_reuse.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

TEST(ShortestPathsReuseFault, CaptureDeriveAndConstructionRetryWithoutPublication) {
  const auto g = shared(make_graph(4, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 10, 1}, {2, 3, 10, 1}}));
  FlowLedger baseline(g), ledger(g);
  FlowPolicy base(baseline, v({0}), v({3}), {ShortestPaths {}});
  (void)base.replace(100);
  FlowPolicy* owners[] = {&base};
  const ShortestPathsReuse reuse(owners);
  const auto mask = g->eligibility({}, v({2}));
  const auto hints = reuse.prepare(&mask);
  const auto original = ledger_state(baseline);
  for (int operation = 0; operation < 3; ++operation) {
    for (long at = 0; at < 2048; ++at) {
      const auto before = ledger_state(ledger);
      const long live = g_live.load();
      const auto call = [&] {
        if (operation == 0) { const ShortestPathsReuse captured(owners); }
        else if (operation == 1) { const auto prepared = reuse.prepare(&mask); }
        else { FlowPolicy p(ledger, v({0}), v({3}), {ShortestPaths {}}, &mask, nullptr, hints.front().get()); }
      };
      const bool failed = fails_at(at, call);
      if (failed) {
        expect_same(ledger_state(baseline), original, "failed immutable preparation");
        expect_same(ledger_state(ledger), before, "failed construction");
        call();
      }
      EXPECT_EQ(g_live.load(), live) << operation << " allocation " << at;
      if (!failed) break;
      ASSERT_LT(at, 2047);
    }
  }
}

TEST(ShortestPathsReuseFault, PlacementFailuresDoNotRetainPrecommitProgram) {
  const auto g = shared(make_graph(2, {{0, 1, 10, 1}, {0, 1, 10, 1}}));
  FlowLedger baseline(g);
  FlowPolicy base(baseline, v({0}), v({1}), {ShortestPaths {}});
  (void)base.replace(100);
  FlowPolicy* owners[] = {&base};
  const ShortestPathsReuse reuse(owners);
  const auto mask = g->eligibility({}, v({0}));
  const auto hint = reuse.prepare(&mask).front();
  for (bool evaluate : {false, true}) {
    const auto injected = sweep_placement(g, "prepared placement", [&](FlowLedger& ledger) {
      auto p = std::make_shared<FlowPolicy>(ledger, v({0}), v({1}), FlowPolicyConfig {ShortestPaths {}}, &mask, nullptr, hint.get());
      return Staged {p, p.get()};
    }, [&](FlowPolicy& p) {
      const OutputOptions out {.edge_flows = true, .residual = true, .program = true, .by_source = true};
      return evaluate ? p.evaluate(7, {}, out) : p.replace(7, {}, out);
    });
    EXPECT_GT(injected, 0);
  }
}
