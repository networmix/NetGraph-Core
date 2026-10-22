// Static preparation ownership, value equivalence, every epoch entry point and
// dynamic exclusions. The fresh control discards only its private retained owner.
#include <gtest/gtest.h>

#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/detail/trace.hpp"
#include "state_comparators.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {
const OutputOptions all {true, true, true, true, true};
std::shared_ptr<const Graph> epoch_graph() {
  return std::make_shared<const Graph>(make_graph(6, {{0, 2, 20, 1}, {1, 2, 20, 1}, {2, 3, 4, 1},
      {2, 3, 10, 1}, {2, 4, 20, 1}, {3, 5, 20, 1}, {4, 5, 20, 1}}));
}
std::vector<FlowPolicyConfig> static_configs(const Graph& g, bool pooled) {
  std::vector<FlowPolicyConfig> out;
  for (auto allocation : {RouteAllocation::Equal, RouteAllocation::Weighted})
    for (auto tie : {RouteTieBreak::FewerHops, RouteTieBreak::HigherCapacity})
      for (bool multipath : {false, true}) for (bool parallel : {false, true}) {
        if (tie == RouteTieBreak::HigherCapacity && ((multipath && parallel) || (pooled && !multipath))) continue;
        for (bool explicit_weights : {false, true}) {
          if (allocation == RouteAllocation::Equal && explicit_weights) continue;
          ShortestPaths sp;
          sp.allocation = allocation; sp.tie_break = tie; sp.multipath = multipath; sp.parallel_edges = parallel;
          if (explicit_weights) {
            sp.weight_basis = WeightBasis::Explicit;
            auto weights = std::make_shared<std::vector<Amount>>(static_cast<std::size_t>(g.num_edges()));
            for (std::size_t i = 0; i < weights->size(); ++i) (*weights)[i] = static_cast<Amount>(i % 3);
            sp.edge_weights = weights;
          }
          for (auto admission : {Admission::Lossless, Admission::Lossy}) out.push_back({sp, admission});
        }
      }
  return out;
}
void discard_kept(FlowPolicy& p) { detail::PolicyAccess::engine(p).kept_program.reset(); }
void compare(const FlowPolicy& p, const FlowPolicy& ref, const FlowLedger& ledger, const FlowLedger& fresh,
             const PlacementResult& got, const PlacementResult& want) {
  expect_same(got, want, "fresh preparation result");
  expect_same(ledger_state(ledger), ledger_state(fresh), "fresh preparation ledger");
  expect_same(epoch(p), epoch(ref), "fresh preparation epoch");
  expect_balanced(ledger, "retained preparation");
}
}  // namespace

TEST(EpochStatic, EveryStaticCombinationReusesAcrossScalarAndListEpochsWithFreshValues) {
  const auto g = epoch_graph();
  for (bool pooled : {false, true}) for (const auto& cfg : static_configs(*g, pooled)) {
    const std::vector<NodeId> sources = pooled ? std::vector<NodeId>{0, 1} : std::vector<NodeId>{0};
    const std::vector<NodeId> targets {5};
    FlowLedger ledger(g), fresh(g);
    FlowPolicy p(ledger, sources, targets, cfg), ref(fresh, sources, targets, cfg);
    auto& engine = detail::PolicyAccess::engine(p);
    ASSERT_TRUE(engine.static_preparation);
    EXPECT_FALSE(engine.kept_program);
    (void)p.evaluate(0, {}, all);
    EXPECT_FALSE(engine.kept_program);
    EXPECT_FALSE(p.program());
    (void)p.place(0);
    EXPECT_FALSE(engine.kept_program);
    auto got = p.place(6, all), want = ref.place(6, all);
    compare(p, ref, ledger, fresh, got, want);
    const auto kept = got.program;
    ASSERT_TRUE(kept);
    EXPECT_EQ(engine.kept_program.get(), kept.get());
    const auto before = epoch(p);
    const auto before_ledger = ledger_state(ledger);
    discard_kept(ref);
    expect_same(p.evaluate(8, {}, all), ref.evaluate(8, {}, all), "warm evaluate");
    EXPECT_EQ(p.evaluate(8).program.get(), kept.get());
    expect_same(epoch(p), before, "evaluate read only");
    expect_same(ledger_state(ledger), before_ledger, "evaluate read only");
    p.remove(); ref.remove();
    EXPECT_FALSE(p.program());
    EXPECT_EQ(engine.kept_program.get(), kept.get());
    EXPECT_FALSE(p.place(0).program);
    EXPECT_THROW((void)p.forward(1), std::invalid_argument);
    discard_kept(ref);
    got = p.place(7, all); want = ref.place(7, all);
    EXPECT_EQ(got.program.get(), kept.get());
    compare(p, ref, ledger, fresh, got, want);
    const auto old_id = p.records().front().flow_id.second;
    discard_kept(ref);
    got = p.replace(3, {}, all); want = ref.replace(3, {}, all);
    EXPECT_EQ(got.program.get(), kept.get());
    EXPECT_GT(p.records().front().flow_id.second, old_id);
    compare(p, ref, ledger, fresh, got, want);
    for (bool replace : {false, true}) {
      if (!replace) { p.remove(); ref.remove(); }
      discard_kept(ref);
      const Action a = replace ? Action(Replace {p, 9}) : Action(Place {p, 9});
      const Action b = replace ? Action(Replace {ref, 9}) : Action(Place {ref, 9});
      got = apply_actions(ledger, std::span(&a, 1), BatchMode::Ordered, all).front();
      want = apply_actions(fresh, std::span(&b, 1), BatchMode::Ordered, all).front();
      EXPECT_EQ(got.program.get(), kept.get());
      compare(p, ref, ledger, fresh, got, want);
    }
    if (cfg.admission == Admission::Lossy) {
      discard_kept(ref);
      const Action a = Replace {p, 4}, b = Replace {ref, 4};
      got = apply_actions(ledger, std::span(&a, 1), BatchMode::Simultaneous, all).front();
      want = apply_actions(fresh, std::span(&b, 1), BatchMode::Simultaneous, all).front();
      EXPECT_EQ(got.program.get(), kept.get());
      compare(p, ref, ledger, fresh, got, want);
    }
    const auto holes = g->eligibility({}, std::vector<EdgeId>{2});
    got = p.forward(8, &holes, all); want = ref.forward(8, &holes, all);
    compare(p, ref, ledger, fresh, got, want);
    discard_kept(ref);
    got = p.replace(5, {}, all); want = ref.replace(5, {}, all);
    EXPECT_EQ(got.program.get(), kept.get());
    compare(p, ref, ledger, fresh, got, want);
  }
}

TEST(EpochStatic, DynamicChoicesNeverRetainAndChangeOnTheNewResidual) {
  for (int mode = 0; mode < 4; ++mode) {
    const auto g = std::make_shared<const Graph>(make_graph(2, {{0, 1, 10, 1}, {0, 1, mode == 1 ? 20 : 10, 1}}));
    const std::vector<NodeId> src {0}, dst {1};
    ShortestPaths sp;
    if (mode < 2) {
      sp.tie_break = mode == 0 ? RouteTieBreak::HigherAvailableCapacity : RouteTieBreak::HigherHeadroom;
      sp.parallel_edges = false;
    } else if (mode == 2) {
      sp.allocation = RouteAllocation::Weighted;
      sp.weight_basis = WeightBasis::ResidualSnapshot;
    } else sp.allocation = RouteAllocation::Optimized;
    FlowLedger ledger(g);
    FlowPolicy p(ledger, src, dst, {sp, mode < 3 ? Admission::Lossy : Admission::Lossless});
    ASSERT_FALSE(detail::PolicyAccess::engine(p).static_preparation);
    const auto first = p.place(2);
    const auto original = program_snapshot(first.program.get());
    // Outside traffic is pinned to e0 so its own route choice cannot hide the change.
    const auto route = std::make_shared<const PredecessorDAG>(*dag_of(*g, {0}));
    RouteDefinition def {0, route};
    const auto outside_program = ForwardingProgram::from_routes(*g, src, dst, std::span(&def, 1));
    FlowPolicy outside(ledger, src, dst, {Pinned {outside_program}, Admission::Lossless});
    ASSERT_EQ(outside.place(mode == 1 ? 1 : 8).delivered, mode == 1 ? 1 : 8);
    const auto next = p.replace(mode < 2 ? 1 : 12);
    EXPECT_NE(next.program.get(), first.program.get());
    EXPECT_FALSE(detail::PolicyAccess::engine(p).kept_program);
    EXPECT_FALSE(program_snapshot(next.program.get()) == original);
    if (mode < 2) EXPECT_EQ(ascending_edges(*p.records().front().route), std::vector<EdgeId> {1});
    else EXPECT_EQ(next.delivered, 12);
    p.remove();
    EXPECT_NE(p.place(1).program.get(), first.program.get());
    for (bool replace : {false, true}) {
      if (!replace) p.remove();
      const auto before_program = p.program();
      const Action action = replace ? Action(Replace {p, 1}) : Action(Place {p, 1});
      const auto result = apply_actions(ledger, std::span(&action, 1)).front();
      EXPECT_NE(result.program.get(), first.program.get());
      if (before_program) EXPECT_NE(result.program.get(), before_program.get());
      EXPECT_FALSE(detail::PolicyAccess::engine(p).kept_program);
    }
    if (mode < 3) {
      const auto previous = p.program();
      const Action action = Replace {p, 1};
      const auto result = apply_actions(ledger, std::span(&action, 1), BatchMode::Simultaneous).front();
      EXPECT_NE(result.program.get(), previous.get());
      EXPECT_FALSE(detail::PolicyAccess::engine(p).kept_program);
    }

  }
}

TEST(EpochStatic, KeptProgramOwnsUnprunedRoutesAndLifetimeWithoutOldRecordIds) {
  const auto g = std::make_shared<const Graph>(make_graph(2, {{0, 1, 0, 1}, {0, 1, 2, 1}}));
  const std::vector<NodeId> s {0}, t {1};
  FlowLedger ledger(g);
  ShortestPaths sp; sp.allocation = RouteAllocation::Weighted;
  auto p = std::make_unique<FlowPolicy>(ledger, s, t, FlowPolicyConfig {sp, Admission::Lossy});
  std::weak_ptr<const ForwardingProgram> weak;
  {
    const auto result = p->place(2);
    ASSERT_EQ(result.delivered, 2);
    weak = result.program;
    const auto& owner = detail::PolicyAccess::engine(*p).kept_program;
    ASSERT_FALSE(owner.owner_before(result.program));
    ASSERT_FALSE(result.program.owner_before(owner));

    EXPECT_EQ(ascending_edges(*p->records().front().route), (std::vector<EdgeId>{0, 1}));
    EXPECT_EQ(detail::ProgramAccess::data(*result.program).split_contexts[0].prepared_route->forwarding_edges,
              (std::vector<EdgeId>{1}));
  }
  const auto id = p->records().front().flow_id;
  p->remove();
  EXPECT_FALSE(p->program()); EXPECT_FALSE(weak.expired());
  auto next = p->place(1);
  EXPECT_EQ(next.program.get(), weak.lock().get());
  EXPECT_GT(p->records().front().flow_id.second, id.second);
  EXPECT_THROW((void)ledger.record(id), std::invalid_argument);
  p->remove(); p.reset();
  EXPECT_FALSE(weak.expired());  // independently returned immutable value owns it
  next.program.reset();
  EXPECT_TRUE(weak.expired());
}

TEST(EpochStatic, ProgramComparatorDetectsWrongRetainedRouteCostWeightAndStructure) {
  const auto g = epoch_graph();
  FlowLedger ledger(g);
  const std::vector<NodeId> sources {0, 1}, dst {5};
  ShortestPaths sp; sp.allocation = RouteAllocation::Weighted;
  FlowPolicy p(ledger, sources, dst, {sp, Admission::Lossy});
  const auto r = p.place(3);
  const auto good = program_snapshot(r.program.get());
  auto wrong = good;
  ASSERT_FALSE(wrong.weights.empty()); ASSERT_FALSE(wrong.weights.front().empty());
  ++wrong.weights.front().front();
  EXPECT_FALSE(good == wrong);
  wrong = good; ASSERT_FALSE(wrong.structure.empty()); ++wrong.structure.back();
  EXPECT_FALSE(good == wrong);
  wrong = good; ++wrong.members.front().front();
  EXPECT_FALSE(good == wrong);
  // A writable, independently owned test value: never write through the real
  // program's const ownership. initialize_storage rebuilds its own block/views.
  const auto& original = detail::ProgramAccess::data(*r.program);
  detail::ProgramData writable;
  writable.graph_uid = original.graph_uid;
  writable.tunnels = original.tunnels;
  writable.source_weights = original.source_weights;
  writable.source_weight_values = original.source_weight_values;
  writable.entries = original.entries;
  writable.split_contexts = original.split_contexts;
  writable.uses_topology_basis = original.uses_topology_basis;
  writable.common_replay_support = original.common_replay_support;
  writable.source_supports = original.source_supports;
  writable.pools_source_inflows = original.pools_source_inflows;
  writable.initialize_storage(original.sources, original.destinations, original.reachable_sources,
                              original.excluded_nodes, original.excluded_edges);
  auto block = std::make_shared<detail::ProgramBlock>(std::move(writable));
  const std::shared_ptr<const ForwardingProgram> poisoned(block, &block->program);
  ASSERT_TRUE(good == program_snapshot(poisoned.get()));
  ++block->data.entries.front().min_route_cost;
  EXPECT_FALSE(good == program_snapshot(poisoned.get()));
  --block->data.entries.front().min_route_cost;
  block->data.entries.front().route = dag_of(*g, {0});
  EXPECT_FALSE(good == program_snapshot(poisoned.get()));
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
namespace {
bool equal_trace(const detail::Trace& a, const detail::Trace& b) {
  if (a.kind != b.kind || a.i != b.i || a.s != b.s || a.items.size() != b.items.size() || a.fields.size() != b.fields.size()) return false;
  for (std::size_t i = 0; i < a.items.size(); ++i) if (!equal_trace(a.items[i], b.items[i])) return false;
  for (std::size_t i = 0; i < a.fields.size(); ++i)
    if (a.fields[i].first != b.fields[i].first || !equal_trace(a.fields[i].second, b.fields[i].second)) return false;
  return true;
}
}
TEST(EpochStatic, TraceHasTheFreshPrepareMapOnReuseSingleAndPooled) {
  const auto g = epoch_graph();
  for (bool pooled : {false, true}) {
    FlowLedger ledger(g);
    const std::vector<NodeId> src = pooled ? std::vector<NodeId>{0, 1} : std::vector<NodeId>{0};
    const std::vector<NodeId> dst {5};
    FlowPolicy p(ledger, src, dst, {});
    detail::Trace first, repeated, increment;
    (void)detail::place_traced(p, 3, all, first);
    p.remove();
    (void)detail::place_traced(p, 3, all, repeated);
    ASSERT_TRUE(equal_trace(first, repeated));
    bool found = false;
    for (const auto& field : repeated.fields) if (field.first == "prepare") {
      found = true; EXPECT_EQ(field.second.kind, detail::Trace::Kind::Map);
    }
    EXPECT_TRUE(found);
    (void)detail::place_traced(p, 1, all, increment);
    for (const auto& field : increment.fields) if (field.first == "prepare") EXPECT_EQ(field.second.kind, detail::Trace::Kind::Null);
  }
}
#endif

TEST(EpochStatic, MaskedBranchJoinAndParallelPopulationsEqualFreshOnChangedResiduals) {
  std::mt19937 random(8371);
  for (int trial = 0; trial < 160; ++trial) {
    std::vector<InEdge> edges;
    for (NodeId u = 0; u < 6; ++u) {
      edges.push_back({u, u + 1, 1 + Amount(random() % 8), 2});
      for (NodeId v = u + 1; v < 8; ++v) if (random() % 3 == 0)
        for (int copy = 0; copy < 2; ++copy) edges.push_back({u, v, Amount(random() % 8), 1 + Cost(random() % 7)});
    }
    edges.push_back({6, 7, 3, 2});
    const auto g = std::make_shared<const Graph>(make_graph(8, edges));
    std::vector<NodeId> excluded;
    if (trial % 5 == 0) excluded = {0, 1};
    else if (trial % 7 == 0) excluded = {6, 7};
    else if (trial % 3 == 0) excluded = {3};
    const auto mask = g->eligibility(excluded, std::vector<EdgeId>{EdgeId(random() % g->num_edges())});
    ShortestPaths sp;
    sp.allocation = trial % 3 == 0 ? RouteAllocation::Equal : RouteAllocation::Weighted;
    if (trial % 2 == 0) { sp.tie_break = RouteTieBreak::HigherCapacity; sp.parallel_edges = false; }
    if (trial % 3 == 2) {
      sp.weight_basis = WeightBasis::Explicit;
      auto weights = std::make_shared<std::vector<Amount>>(static_cast<std::size_t>(g->num_edges()));
      for (auto& w : *weights) w = Amount(random() % 4);
      sp.edge_weights = weights;
    }
    const FlowPolicyConfig cfg {sp, trial % 2 == 0 ? Admission::Lossless : Admission::Lossy};
    FlowLedger ledger(g), fresh(g);
    const std::vector<NodeId> sources {0, 1}, destinations {6, 7};
    FlowPolicy p(ledger, sources, destinations, cfg, &mask), ref(fresh, sources, destinations, cfg, &mask);
    auto got = p.place(3, all), want = ref.place(3, all);
    compare(p, ref, ledger, fresh, got, want);
    p.remove(); ref.remove();
    const auto chosen = EdgeId(random() % g->num_edges());
    const std::vector<NodeId> from {g->edge_sources()[static_cast<std::size_t>(chosen)]};
    const std::vector<NodeId> to {g->edge_destinations()[static_cast<std::size_t>(chosen)]};
    const RouteDefinition route {from[0], dag_of(*g, {chosen})};
    const auto outside = ForwardingProgram::from_routes(*g, from, to, std::span(&route, 1));
    const FlowPolicyConfig pinned {Pinned {outside}, Admission::Lossless};
    FlowPolicy q(ledger, from, to, pinned), q_ref(fresh, from, to, pinned);
    (void)q.place(3); (void)q_ref.place(3);
    discard_kept(ref);
    got = p.place(9, all); want = ref.place(9, all);
    compare(p, ref, ledger, fresh, got, want);
  }
}

TEST(EpochStatic, WideCountersAndLargestLegalRouteCostAreFreshAfterRemove) {
  constexpr Amount m = kL62 - 1;
  const auto g = std::make_shared<const Graph>(make_graph(3, {{0, 1, 0, kL62 / 2 - 1}, {1, 2, 0, kL62 / 2}}));
  FlowLedger ledger(g);
  const std::vector<NodeId> source {0}, destination {2};
  FlowPolicy p(ledger, source, destination, {ShortestPaths {}, Admission::Lossy});
  PlacementResult r;
  for (int call = 0; call < 5; ++call) r = p.place(m);
  ASSERT_TRUE(r.offered); ASSERT_TRUE(r.requested);
  EXPECT_EQ(r.offered->hi, 1U); EXPECT_EQ(r.offered->lo, static_cast<std::uint64_t>(kL62 - 5));
  EXPECT_EQ(*r.requested, *r.offered);
  EXPECT_EQ(p.records().front().min_route_cost, m);
  const auto kept = r.program;
  p.remove();
  r = p.place(1);
  EXPECT_EQ(r.program.get(), kept.get());
  EXPECT_EQ(r.offered->hi, 0U); EXPECT_EQ(r.offered->lo, 1U);
  EXPECT_EQ(p.records().front().min_route_cost, m);
}

TEST(EpochStatic, ColdListPreparationSeedsTheOwnerBeforeLaterEpochs) {
  const auto g = epoch_graph();
  for (bool pooled : {false, true}) for (const auto& cfg : static_configs(*g, pooled))
    for (auto mode : {BatchMode::Ordered, BatchMode::Simultaneous}) {
      if (mode == BatchMode::Simultaneous && cfg.admission != Admission::Lossy) continue;
      FlowLedger ledger(g);
      const std::vector<NodeId> sources = pooled ? std::vector<NodeId>{0, 1} : std::vector<NodeId>{0};
      FlowPolicy p(ledger, sources, std::vector<NodeId>{5}, cfg);
      const Action action = Replace {p, 3};
      const auto first = apply_actions(ledger, std::span(&action, 1), mode, all).front();
      EXPECT_EQ(detail::PolicyAccess::engine(p).kept_program.get(), first.program.get());
      p.remove();
      const auto repeated = apply_actions(ledger, std::span(&action, 1), mode, all).front();
      EXPECT_EQ(repeated.program.get(), first.program.get());
      expect_same(repeated, first, "cold-list then retained list");
      p.remove();
      EXPECT_EQ(p.place(3).program.get(), first.program.get());
    }
}

TEST(EpochStatic, SuppliedCacheStatisticsFreezeAndClearKeepTheirOriginalBehavior) {
  const auto g = std::make_shared<const Graph>(make_graph(2, {{0, 1, 4, 1}, {0, 1, 10, 1}}));
  ShortestPathCache cache(*g);
  {
    FlowLedger ledger(g);
    const std::vector<NodeId> source {0}, destination {1};
    FlowPolicy p(ledger, source, destination, {}, nullptr, &cache);
    const auto before = cache.stats();
    ASSERT_EQ(before.entries, 1); ASSERT_EQ(before.misses, 1); ASSERT_EQ(before.live_values, 1);
    cache.freeze();
    const auto original = p.place(3).program;
    for (int epoch = 0; epoch < 3; ++epoch) {
      p.remove();
      EXPECT_EQ(p.place(3).program.get(), original.get());
      EXPECT_EQ(p.evaluate(3).program.get(), original.get());
      const auto now = cache.stats();
      EXPECT_EQ(now.hits, before.hits); EXPECT_EQ(now.misses, before.misses);
      EXPECT_EQ(now.entries, before.entries); EXPECT_EQ(now.bytes, before.bytes); EXPECT_EQ(now.live_values, before.live_values);
    }
    EXPECT_THROW(cache.clear(), std::invalid_argument);
    cache.thaw(); cache.clear();
    EXPECT_EQ(cache.stats().entries, 0); EXPECT_EQ(cache.stats().bytes, 0); EXPECT_EQ(cache.stats().live_values, 1);
    p.remove();
    EXPECT_EQ(p.place(1).program.get(), original.get());
    EXPECT_EQ(cache.stats().misses, before.misses);
  }
  EXPECT_EQ(cache.stats().live_values, 0);
}

namespace {
void check_separate_programs(bool explicit_weights) {
  const auto g = std::make_shared<const Graph>(make_graph(2, {{0, 1, 12, 1}, {0, 1, 12, 1}}));
  const auto mask = g->eligibility({}, std::vector<EdgeId>{0});
  const std::vector<NodeId> source {0}, destination {1};
  ShortestPaths first_choice, second_choice;
  if (explicit_weights) {
    first_choice.allocation = second_choice.allocation = RouteAllocation::Weighted;
    first_choice.weight_basis = second_choice.weight_basis = WeightBasis::Explicit;
    first_choice.edge_weights = std::make_shared<const std::vector<Amount>>(std::vector<Amount>{1, 2});
    second_choice.edge_weights = std::make_shared<const std::vector<Amount>>(std::vector<Amount>{2, 1});
  }
  const FlowPolicyConfig first_config {first_choice, Admission::Lossless};
  const FlowPolicyConfig second_config {second_choice, Admission::Lossless};
  const Eligibility* second_mask = explicit_weights ? nullptr : &mask;
  FlowLedger ledger(g), fresh(g);
  FlowPolicy first(ledger, source, destination, first_config);
  FlowPolicy second(ledger, source, destination, second_config, second_mask);
  FlowPolicy first_ref(fresh, source, destination, first_config);
  FlowPolicy second_ref(fresh, source, destination, second_config, second_mask);
  std::shared_ptr<const ForwardingProgram> kept_first, kept_second;
  auto check = [&](FlowPolicy& policy, FlowPolicy& ref, std::shared_ptr<const ForwardingProgram>& kept,
                   const std::shared_ptr<const ForwardingProgram>& other, Amount offer,
                   const std::vector<Amount>& carried, const std::vector<EdgeId>& edges) {
    discard_kept(ref);
    const auto got = policy.place(offer, all), want = ref.place(offer, all);
    compare(policy, ref, ledger, fresh, got, want);
    ASSERT_TRUE(got.program);
    if (!kept) kept = got.program;
    EXPECT_EQ(got.program.get(), kept.get());
    if (other) EXPECT_NE(got.program.get(), other.get());
    EXPECT_EQ(detail::PolicyAccess::engine(policy).kept_program.get(), kept.get());
    EXPECT_EQ(got.delivered, offer);
    EXPECT_EQ(ledger.carried(), carried);
    const auto records = policy.records();
    ASSERT_EQ(records.size(), 1U);
    EXPECT_EQ(records.front().delivered, offer);
    EXPECT_EQ(records.front().min_route_cost, 1);
    EXPECT_EQ(ascending_edges(*records.front().route), edges);
    std::vector<std::pair<EdgeId, Amount>> edge_flows;
    for (EdgeId e = 0; e < 2; ++e)
      if (carried[static_cast<std::size_t>(e)]) edge_flows.emplace_back(e, carried[static_cast<std::size_t>(e)]);
    EXPECT_EQ(records.front().edge_flows, edge_flows);
    policy.remove(); ref.remove();
    EXPECT_FALSE(policy.program());
    EXPECT_EQ(detail::PolicyAccess::engine(policy).kept_program.get(), kept.get());
    expect_same(ledger_state(ledger), ledger_state(fresh), "separate-policy removal");
    EXPECT_TRUE(ledger.records().empty());
  };
  for (Amount offer : {3, 6, 9, 6}) {
    SCOPED_TRACE(offer);
    // Equal RR gives [(offer+1)/2, offer/2]; masking e0 gives [0, offer].
    // For offers divisible by 3, weights [1,2] and [2,1] give [n,2n] and [2n,n].
    check(first, first_ref, kept_first, kept_second, offer,
          explicit_weights ? std::vector<Amount>{offer / 3, 2 * (offer / 3)}
                           : std::vector<Amount>{(offer + 1) / 2, offer / 2}, {0, 1});
    check(second, second_ref, kept_second, kept_first, offer,
          explicit_weights ? std::vector<Amount>{2 * (offer / 3), offer / 3}
                           : std::vector<Amount>{0, offer},
          explicit_weights ? std::vector<EdgeId>{0, 1} : std::vector<EdgeId>{1});
    ASSERT_TRUE(kept_first); ASSERT_TRUE(kept_second);
    EXPECT_NE(kept_first.get(), kept_second.get());
    EXPECT_FALSE(program_snapshot(kept_first.get()) == program_snapshot(kept_second.get()));
  }
}
}  // namespace

TEST(EpochStatic, SeparateMaskedPoliciesKeepOwnProgramsAcrossEpochs) {
  check_separate_programs(false);
}

TEST(EpochStatic, SeparateExplicitWeightPoliciesKeepOwnProgramsAcrossEpochs) {
  check_separate_programs(true);
}
