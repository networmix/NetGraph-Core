#include <gtest/gtest.h>

#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/shortest_paths_reuse.hpp"
#include "state_comparators.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

TEST(ShortestPathsReuse, PrecommitIdentityAndFreshCaseValues) {
  const auto g = shared(make_graph(2, {{0, 1, 10, 1}, {0, 1, 10, 1}}));
  FlowLedger baseline(g), ledger(g), reference(g);
  FlowPolicy base(baseline, v({0}), v({1}), {ShortestPaths {}});
  (void)base.replace(100);
  FlowPolicy* owners[] = {&base};
  const ShortestPathsReuse reuse(owners);
  const auto mask = g->eligibility({}, v({0}));
  const auto hints = reuse.prepare(&mask);
  ASSERT_TRUE(hints.front());
  FlowPolicy p(ledger, v({0}), v({1}), {ShortestPaths {}}, &mask, nullptr, hints.front().get());
  FlowPolicy ref(reference, v({0}), v({1}), {ShortestPaths {}}, &mask);
  EXPECT_TRUE(detail::PolicyAccess::engine(p).preparation_consumed);
  EXPECT_FALSE(detail::PolicyAccess::engine(p).kept_program);
  const OutputOptions out {.edge_flows = true, .residual = true, .saturated_boundary_edges = true,
      .min_cut = true, .edge_reservations = true, .program = true, .by_source = true};
  const auto a = p.evaluate(7, {}, out), b = p.evaluate(7, {}, out);
  ASSERT_TRUE(a.program);
  EXPECT_NE(a.program.get(), b.program.get());
  EXPECT_FALSE(p.program());
  EXPECT_FALSE(detail::PolicyAccess::engine(p).kept_program);
  const auto got = p.replace(7, {}, out), want = ref.replace(7, {}, out);
  EXPECT_NE(a.program.get(), got.program.get());
  EXPECT_EQ(got.program.get(), p.program().get());
  expect_same(got, want, "restricted case");
  expect_same(ledger_state(ledger), ledger_state(reference), "restricted case state");
  p.remove();
  EXPECT_EQ(p.replace(8, {}, out).program.get(), got.program.get());
}

TEST(ShortestPathsReuse, EffectiveContainmentAndFailedSources) {
  const auto g = shared(make_graph(3, {{0, 1, 2, 1}, {0, 2, 2, 1}}));
  const auto base_mask = g->eligibility(v({1}), v({0}));
  FlowLedger baseline(g);
  FlowPolicy base(baseline, v({0}), v({1, 2}), {ShortestPaths {}}, &base_mask);
  (void)base.replace(100);
  FlowPolicy* owners[] = {&base};
  const ShortestPathsReuse reuse(owners);
  const auto alias = g->eligibility(v({1}), {});
  EXPECT_TRUE(reuse.prepare(&alias).front());
  EXPECT_FALSE(reuse.prepare().front());
  const auto failed_source = g->eligibility(v({0, 1}), {});
  const auto hints = reuse.prepare(&failed_source);
  ASSERT_TRUE(hints.front());
  FlowLedger ledger(g);
  FlowPolicy p(ledger, v({0}), v({1, 2}), {ShortestPaths {}}, &failed_source, nullptr, hints.front().get());
  const auto r = p.replace(100);
  EXPECT_EQ(r.termination, Termination::NoRoute);
  EXPECT_FALSE(r.routes_changed);
  EXPECT_FALSE(r.weights_changed);
}

TEST(ShortestPathsReuse, CaptureOwnsItsProgramAfterOriginalPoliciesDie) {
  const auto g = shared(make_graph(2, {{0, 1, 10, 1}, {0, 1, 10, 1}}));
  const auto reuse = [&] {
    FlowLedger baseline(g);
    FlowPolicy base(baseline, v({0}), v({1}), {ShortestPaths {}});
    (void)base.replace(100);
    FlowPolicy* owners[] = {&base};
    return ShortestPathsReuse(owners);
  }();
  const auto mask = g->eligibility({}, v({0}));
  const auto hint = reuse.prepare(&mask).front();
  ASSERT_TRUE(hint);
  FlowLedger ledger(g);
  FlowPolicy p(ledger, v({0}), v({1}), {ShortestPaths {}}, &mask, nullptr, hint.get());
  const auto result = p.replace(100, {}, OutputOptions {.program = true, .by_source = true});
  EXPECT_EQ(result.delivered, 10);
  EXPECT_EQ(p.records().front().min_route_cost, 1);
}

TEST(ShortestPathsReuse, SurvivingSourcesGetFreshEntryIdentifiers) {
  const auto g = shared(make_graph(3, {{0, 2, 3, 2}, {1, 2, 2, 1}}));
  FlowLedger baseline(g);
  FlowPolicy base(baseline, v({0, 1}), v({2}), {ShortestPaths {}});
  (void)base.replace(100);
  FlowPolicy* owners[] = {&base};
  const ShortestPathsReuse reuse(owners);
  const auto mask = g->eligibility(v({0}), {});
  const auto hint = reuse.prepare(&mask).front();
  ASSERT_TRUE(hint);
  FlowLedger ledger(g);
  FlowPolicy p(ledger, v({0, 1}), v({2}), {ShortestPaths {}}, &mask, nullptr, hint.get());
  const auto r = p.replace(100, {}, OutputOptions {.program = true, .by_source = true});
  ASSERT_EQ(r.delivered, 2);
  const auto& data = detail::ProgramAccess::data(*r.program);
  ASSERT_EQ(data.entries.size(), 1U);
  EXPECT_EQ(data.entries.front().entry_id, 0);
  EXPECT_EQ(data.entries.front().head, 1);
  EXPECT_EQ(data.entries.front().roster_index, 0);
  EXPECT_EQ(data.entries.front().min_route_cost, 1);
}
