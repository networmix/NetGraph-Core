// Prescribed source quotas and replay. Expected examples and small exhaustive
// flows are independent of production apportionment and flow-core helpers.
#include <gtest/gtest.h>

#include <functional>
#include <future>
#include <numeric>

#include "netgraph/core/detail/forwarding.hpp"
#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/evaluators.hpp"
#include "state_comparators.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {
const OutputOptions all {true, true, true, true, true};
FlowPolicyConfig weighted(std::vector<Amount> w, std::shared_ptr<const PredecessorDAG> support = {}) {
  CentralizedTE model(SourceAllocation::Weighted);
  model.source_weights = std::make_shared<const std::vector<Amount>>(std::move(w));
  model.support = std::move(support);
  return {.routing = std::move(model)};
}
std::shared_ptr<const Graph> witness(bool coarse = false) {
  const Amount q = coarse ? 1 : 8;
  return std::make_shared<const Graph>(make_graph(6, {{0,2,2*q,1}, {1,2,q,1}, {1,3,q,2},
      {2,3,q,1}, {2,4,q,1}, {3,5,q,1}, {4,5,q,1}}));
}
std::vector<Amount> sources(const PlacementResult& r) {
  std::vector<Amount> out;
  for (const auto& row : r.by_source) out.push_back(row.delivered);
  return out;
}
std::string refused(const std::function<void()>& call) {
  try { call(); } catch (const std::invalid_argument& e) { return e.what(); }
  return "no error";
}
// One-unit-at-a-time priority oracle for small totals; test integer arithmetic
// is independent of the library's portable multiplication/division.
std::vector<Amount> sequential(Amount t, const std::vector<Amount>& w) {
  std::vector<Amount> a(w.size());
  if (std::none_of(w.begin(), w.end(), [](Amount x) { return x != 0; })) return a;
  for (Amount unit = 0; unit < t; ++unit) {
    std::size_t best = 0;
    for (std::size_t j = 1; j < w.size(); ++j) {
      if (UInt128(w[j]) * UInt128(a[best] + 1) > UInt128(w[best]) * UInt128(a[j] + 1)) best = j;
    }
    ++a[best];
  }
  return a;
}
} // namespace

TEST(WeightedSources, WitnessAndDeclaredReplayOnDenseAndSupportPaths) {
  for (bool support : {false, true}) {
    const auto g = witness();
    const auto cfg = weighted({2,1,0,0,0,0}, support ? dag_of(*g, {0,1,2,3,4,5,6}) : nullptr);
    FlowLedger ledger(g);
    FlowPolicy p(ledger, v({0,1}), v({5}), cfg);
    const auto before = ledger_state(ledger);
    const auto r = p.evaluate(std::nullopt, {}, all);
    expect_same(ledger_state(ledger), before, "evaluate is read only");
    EXPECT_FALSE(p.program());
    expect_same(p.replace(std::nullopt, {}, all), r, "evaluate/replace");
    EXPECT_EQ(r.delivered, 16);
    EXPECT_FALSE(r.call_undelivered);
    EXPECT_EQ(r.by_tunnel, (std::vector<TunnelFlow>{{0,0,TunnelStatus::Up,{0,3},3,{}},
        {1,0,TunnelStatus::Up,{0,8},8,{}},{2,1,TunnelStatus::Up,{0,5},5,{}}}));
    EXPECT_TRUE(r.routes_changed); EXPECT_TRUE(r.weights_changed); EXPECT_FALSE(r.reservations_changed);
    EXPECT_EQ(sources(r), (std::vector<Amount>{11,5}));
    EXPECT_EQ(r.total_cost, (Counter128{0,48}));
    EXPECT_EQ(r.edge_flows, (std::optional<std::vector<Amount>>{{11,0,5,3,8,8,8}}));
    ASSERT_TRUE(r.program->source_weights());
    EXPECT_EQ(*r.program->source_weights(), (std::vector<std::pair<NodeId,Amount>>{{0,2},{1,1}}));
    ASSERT_EQ(r.program->tunnels().size(), 3U);
    const std::vector<std::vector<EdgeId>> paths {{0,3,5},{0,4,6},{2,5}};
    const Amount amounts[] = {3,8,5};
    for (std::size_t i = 0; i < paths.size(); ++i) {
      const auto& t = r.program->tunnels()[i];
      EXPECT_EQ(t.entry_id, static_cast<std::int32_t>(i));
      EXPECT_EQ(t.head, i == 2 ? 1 : 0);
      EXPECT_EQ(t.edges, paths[i]);
      EXPECT_EQ(t.weight, amounts[i]);
      EXPECT_EQ(t.path_cost, 3);
    }
    EXPECT_FALSE(r.plan); EXPECT_FALSE(r.min_cut); EXPECT_FALSE(r.marginal_cost_curve);
    EXPECT_FALSE(r.offered); EXPECT_FALSE(r.unadmitted);
    EXPECT_EQ(r.source_drops, (Counter128{})); EXPECT_TRUE(r.edge_drops.empty());
    EXPECT_EQ(r.edge_reservations, (std::optional<std::vector<Amount>>(std::vector<Amount>(7,0))));
    EXPECT_FALSE(detail::PolicyAccess::engine(p).static_preparation);
    EXPECT_FALSE(detail::PolicyAccess::engine(p).kept_program);
    for (Admission admission : {Admission::Lossless, Admission::Lossy}) {
      FlowLedger replay(g);
      FlowPolicy q(replay, v({0,1}), v({5}), {.routing=Pinned{r.program}, .admission=admission});
      const auto at15 = q.replace(15, {}, all);
      EXPECT_EQ(sources(at15), (std::vector<Amount>{10,5}));
      EXPECT_EQ(at15.edge_flows, (std::optional<std::vector<Amount>>{{10,0,5,3,7,8,7}}));
      EXPECT_EQ(at15.termination, Termination::Complete);
      EXPECT_EQ(q.replace(16).delivered, 16);
      EXPECT_EQ(refused([&] { (void)q.forward(15); }), "Pinned forwards with replace");
    }
    const auto forwarded=p.forward(15);
    EXPECT_EQ(sources(forwarded), (std::vector<Amount>{10,5}));
    EXPECT_FALSE(forwarded.routes_changed);EXPECT_FALSE(forwarded.weights_changed);EXPECT_FALSE(forwarded.reservations_changed);
    expect_balanced(ledger, "controller replay");
  }
}

TEST(WeightedSources, CoarsePositiveWeightWithoutAPathKeepsItsShare) {
  const auto g = witness(true);
  FlowLedger ledger(g);
  FlowPolicy p(ledger, v({0,1}), v({5}), weighted({2,1,0,0,0,0}));
  const auto r = p.replace(std::nullopt, {}, all);
  EXPECT_EQ(r.delivered, 2); EXPECT_EQ(r.total_cost, (Counter128{0,6}));
  EXPECT_EQ(sources(r), (std::vector<Amount>{2,0}));
  EXPECT_EQ(r.edge_flows, (std::optional<std::vector<Amount>>{{2,0,0,1,1,1,1}}));
  ASSERT_TRUE(r.program->source_weights());
  EXPECT_EQ(*r.program->source_weights(), (std::vector<std::pair<NodeId,Amount>>{{0,2},{1,1}}));
  ASSERT_EQ(r.program->tunnels().size(), 2U);
  for (const auto& t : r.program->tunnels()) { EXPECT_EQ(t.head,0); EXPECT_EQ(t.weight,1); }
  for (Admission admission : {Admission::Lossless, Admission::Lossy}) {
    FlowLedger replay(g);
    FlowPolicy q(replay, v({0,1}), v({5}), {.routing=Pinned{r.program}, .admission=admission});
    const auto at3 = q.replace(3);
    EXPECT_EQ(at3.delivered, 2);
    EXPECT_EQ(at3.source_drops, (Counter128{0,admission==Admission::Lossy ? 1U : 0U}));
    EXPECT_EQ(at3.termination, admission==Admission::Lossy ? Termination::Complete : Termination::Exhausted);
  }
}

TEST(WeightedSources, RosterZeroWeightsAndEmptyReplayKeepDifferentOutcomes) {
  const auto g = std::make_shared<const Graph>(make_graph(4, {{0,3,4,1},{1,3,4,1},{2,3,4,1}}));
  for (bool masked : {false,true}) {
    const auto mask = g->eligibility(masked ? v({2}) : v({}), {});
    FlowLedger ledger(g);
    FlowPolicy p(ledger, v({2,0,1,2}), v({3}), weighted({1,1,2,0}), &mask);
    const auto r = p.replace(std::nullopt, {}, all);
    EXPECT_EQ(r.delivered,8); EXPECT_EQ(r.total_cost,(Counter128{0,8}));
    EXPECT_EQ(sources(r), (masked ? std::vector<Amount>{4,4,0} : std::vector<Amount>{2,2,4}));
  }
  for (bool support : {false,true}) for (bool empty : {false,true}) {
    const auto mask = g->eligibility(empty ? v({0,1,2}) : v({2}), {});
    FlowLedger ledger(g);
    FlowPolicy p(ledger, v({0,1,2}), v({3}), weighted({0,0,5,0}, support ? dag_of(*g,{0,1,2}) : nullptr), &mask);
    const auto r = p.replace(6, {}, all);
    EXPECT_EQ(r.termination,empty ? Termination::NoRoute : Termination::NoActiveWeight);
    EXPECT_EQ(r.delivered,0); EXPECT_EQ(r.unadmitted,(Counter128{0,6}));
    EXPECT_EQ(sources(r),(std::vector<Amount>{0,0,0}));
    ASSERT_TRUE(r.program->source_weights());
    EXPECT_EQ(r.program->source_weights()->size(),empty ? 0U : 2U);
    EXPECT_TRUE(r.program->tunnels().empty());
    const auto unlimited=p.evaluate(std::nullopt);
    EXPECT_EQ(unlimited.termination,empty ? Termination::NoRoute : Termination::NoActiveWeight);
    EXPECT_FALSE(unlimited.requested);EXPECT_FALSE(unlimited.unadmitted);EXPECT_FALSE(unlimited.call_undelivered);
    EXPECT_EQ(p.forward(6).termination,empty ? Termination::NoRoute : Termination::Exhausted);
    for (Admission admission : {Admission::Lossless,Admission::Lossy}) {
      FlowLedger replay(g);
      FlowPolicy q(replay,v({0,1,2}),v({3}),{.routing=Pinned{r.program},.admission=admission},&mask);
      EXPECT_EQ(q.replace(0).termination,Termination::Complete);
      EXPECT_EQ(q.evaluate(6).termination,empty ? Termination::NoRoute : admission==Admission::Lossy ? Termination::Complete : Termination::Exhausted);
      const auto at6=q.place(6);
      EXPECT_EQ(at6.source_drops,(Counter128{0,admission==Admission::Lossy ? 6U : 0U}));
      EXPECT_EQ(at6.termination,empty ? Termination::NoRoute : admission==Admission::Lossy ? Termination::Complete : Termination::Exhausted);
      EXPECT_EQ(refused([&]{(void)q.forward(6);}),"Pinned forwards with replace");
    }
  }
}

TEST(WeightedSources, EqualPositiveWeightsAndZeroWeightHead) {
  const auto g=std::make_shared<const Graph>(make_graph(4,{{0,3,2,1},{1,2,10,1},{2,3,10,2}}));
  auto r=evaluate_demand(*g,v({0,1}),v({3}),EvaluationOptions{weighted({7,7,0,0}),nullptr,nullptr,all,{},Amount{12}});
  EXPECT_EQ(r.delivered,4); EXPECT_EQ(r.total_cost,(Counter128{0,8}));
  EXPECT_EQ(r.edge_flows,(std::optional<std::vector<Amount>>{{2,2,2}}));
  EXPECT_EQ(sources(r),(std::vector<Amount>{2,2}));
  const auto h=std::make_shared<const Graph>(make_graph(3,{{0,2,5,1},{1,2,5,1}}));
  r=evaluate_demand(*h,v({0,1}),v({2}),EvaluationOptions{weighted({0,3,0}),nullptr,nullptr,all});
  EXPECT_EQ(r.delivered,5); EXPECT_EQ(r.total_cost,(Counter128{0,5}));
  EXPECT_EQ(sources(r),(std::vector<Amount>{0,5}));
  EXPECT_EQ(r.edge_flows,(std::optional<std::vector<Amount>>{{0,5}}));
  ASSERT_TRUE(r.program->source_weights());
  EXPECT_EQ(*r.program->source_weights(),(std::vector<std::pair<NodeId,Amount>>{{0,0},{1,3}}));
}

TEST(WeightedSources, RefusalReliefUsesDeliveredRatherThanRequestedInScalarAndOrdered) {
  const auto g=std::make_shared<const Graph>(make_graph(4,{{0,2,10,1},{1,2,10,1},{2,3,10,1}}));
  for(bool ordered:{false,true}) {
    FlowLedger ledger(g);
    FlowPolicy other(ledger,v({2}),v({3}),{}),p(ledger,v({0,1}),v({3}),weighted({2,1,0,0}));
    (void)other.replace(8);
    const auto call=[&](Amount offer) {
      if(!ordered)return p.place(offer,all);
      const std::vector<Action> actions{Place{p,offer}};
      return apply_actions(ledger,actions,BatchMode::Ordered,all).front();
    };
    const auto initial=call(5);
    EXPECT_EQ(initial.delivered,2); EXPECT_EQ(initial.call_undelivered,3);
    EXPECT_EQ(initial.by_tunnel,(std::vector<TunnelFlow>{{0,0,TunnelStatus::Up,{0,2},2,{}}}));
    (void)other.replace(5);
    const auto r=call(1);
    EXPECT_EQ(r.delivered,3); EXPECT_EQ(r.call_delivered,1); EXPECT_EQ(r.total_cost,(Counter128{0,6}));
    EXPECT_EQ(r.requested,(Counter128{0,6})); EXPECT_EQ(r.unadmitted,(Counter128{0,3}));
    EXPECT_EQ(sources(r),(std::vector<Amount>{2,1}));
    EXPECT_EQ(r.call_undelivered,0);
    EXPECT_EQ(r.by_tunnel,(std::vector<TunnelFlow>{{0,0,TunnelStatus::Up,{0,2},2,{}},{1,1,TunnelStatus::Up,{0,1},1,{}}}));
    EXPECT_EQ(r.edge_flows,(std::optional<std::vector<Amount>>{{2,1,3}}));
    const auto state=epoch(p); (void)p.evaluate(4); expect_same(epoch(p),state,"evaluate after refusal");
    p.remove(); EXPECT_FALSE(p.program());
  }
}

TEST(WeightedSources, SimultaneousAwardsFollowCreationOrderInBothActionOrders) {
  const auto g=std::make_shared<const Graph>(make_graph(2,{{0,1,3,1}}));
  const auto program=evaluate_demand(*g,v({0}),v({1}),EvaluationOptions{weighted({1,0})}).program;
  for(bool reverse:{false,true}) for(bool reverse_creation:{false,true}) {
    FlowLedger ledger(g);
    FlowPolicy first(ledger,v({0}),v({1}),{.routing=Pinned{program},.admission=Admission::Lossy});
    FlowPolicy second(ledger,v({0}),v({1}),{.routing=Pinned{program},.admission=Admission::Lossy});
    FlowPolicy* a=reverse_creation ? &second : &first;
    FlowPolicy* b=reverse_creation ? &first : &second;
    if(reverse)std::swap(a,b);
    const std::vector<Action> actions{Replace{*a,3},Replace{*b,3}};
    const auto r=apply_actions(ledger,actions,BatchMode::Simultaneous,all);
    EXPECT_EQ(r[0].delivered,a==&first ? 2 : 1); EXPECT_EQ(r[1].delivered,b==&first ? 2 : 1);
    EXPECT_EQ(first.delivered(),2); EXPECT_EQ(second.delivered(),1);
    first.remove(); second.remove();
    const auto ordered=apply_actions(ledger,actions,BatchMode::Ordered,all);
    EXPECT_EQ(ordered[0].delivered,3); EXPECT_EQ(ordered[1].delivered,0);
  }
}

TEST(WeightedSources, VectorMatrixAndValidationOrderHaveSuccessfulNeighbours) {
  const auto g=std::make_shared<const Graph>(make_graph(3,{{0,2,2,1},{1,2,2,1}}));
  const auto check=[&](FlowPolicyConfig cfg) {FlowLedger ledger(g);return refused([&]{FlowPolicy p(ledger,v({0}),v({2}),cfg);});};
  auto cfg=weighted({1,0,0});
  EXPECT_EQ(check(cfg),"no error");
  std::get<CentralizedTE>(cfg.routing).tunnel_count=TunnelCount{1};
  EXPECT_EQ(check(cfg),"source_allocation=WEIGHTED requires CentralizedTE without tunnel_count");
  cfg.routing=CentralizedTE{SourceAllocation::Weighted};
  EXPECT_EQ(check(cfg),"source_allocation=WEIGHTED requires source_weights");
  cfg.admission=Admission::Lossy;
  EXPECT_EQ(check(cfg),"CentralizedTE without tunnel_count requires admission=LOSSLESS");
  for(auto allocation:{RouteAllocation::Equal,RouteAllocation::Weighted,RouteAllocation::Optimized}) {
    ShortestPaths sp;sp.allocation=allocation;sp.source_allocation=SourceAllocation::Weighted;
    EXPECT_EQ(check({sp}),"source_allocation=WEIGHTED requires CentralizedTE without tunnel_count");
  }
  for(auto allocation:{SourceAllocation::Equal,SourceAllocation::Optimized}) {
    cfg=weighted({1,0,0}); std::get<CentralizedTE>(cfg.routing).source_allocation=allocation;
    EXPECT_EQ(check(cfg),"source_weights require source_allocation=WEIGHTED");
  }
  for(const auto& w:std::vector<std::vector<Amount>>{{1,0},{-1,0},{0,kDomainLimit},{1,0,0,kDomainLimit}}) {
    EXPECT_EQ(check(weighted(w)),"source_weights: expected 3 entries");
    EXPECT_EQ(refused([&]{(void)evaluate_demand(*g,v({0}),v({2}),EvaluationOptions{weighted(w),nullptr,nullptr,all,{},Amount{0}});}),
        "source_weights: expected 3 entries");
  }
  for(Amount bad:{Amount{-1},kDomainLimit}) EXPECT_EQ(check(weighted({1,0,bad})),"source_weights[2] must be an integer in [0, 2^62)");
  EXPECT_EQ(check(weighted({kDomainLimit-1,0,0})),"no error");
  cfg=weighted({1,0});
  ShortestPathCache cache(*g);FlowLedger ledger(g);
  EXPECT_EQ(refused([&]{FlowPolicy p(ledger,v({0}),v({2}),cfg,nullptr,&cache);}),"CentralizedTE does not use ShortestPathCache");
}

TEST(WeightedSources, FailureStudyRecomputesUnusedEdgesAndEveryWorkerMatches) {
  const auto g=std::make_shared<const Graph>(make_graph(3,{{0,2,0,1},{1,2,10,1}}));
  for(int workers:{1,2,4}) {
    const EvaluationOptions options{weighted({2,1,0}),nullptr,nullptr,all};
    const auto result=evaluate_edge_failures(*g,v({0,1}),v({2}),options,workers);
    EXPECT_EQ(result.baseline.delivered,0);ASSERT_EQ(result.entries.size(),2U);
    EXPECT_EQ(result.entries[0].edge,0);EXPECT_EQ(result.entries[0].delivery_reduction,-10);
    EXPECT_EQ(result.entries[1].edge,1);EXPECT_EQ(result.entries[1].delivery_reduction,0);
    const std::vector<std::pair<std::vector<NodeId>,std::vector<NodeId>>> pairs(5,{{1,0,1},{2}});
    const auto batch=evaluate_pairs(*g,pairs,options,{},workers);
    for(const auto& r:batch)expect_same(r,result.baseline,"pair-local source weights");
  }
}

TEST(WeightedSources, SmallFlowEnumerationConfirmsThePrescribedMaximumAndMinimumCost) {
  EXPECT_EQ(sequential(2,{2,1}),(std::vector<Amount>{2,0}));
  EXPECT_EQ(sequential(50,{100,1,1,1}),(std::vector<Amount>{50,0,0,0}));
  {
    const auto g=make_graph(3,{{0,2,2,9},{1,2,10,1}});
    const auto r=evaluate_demand(g,v({0,1}),v({2}),EvaluationOptions{weighted({2,1,0})});
    EXPECT_EQ(r.delivered,3);EXPECT_EQ(r.total_cost,(Counter128{0,19}));EXPECT_EQ(sources(r),(std::vector<Amount>{2,1}));
    const auto h=make_graph(3,{{0,1,2,1},{0,2,2,5},{1,2,2,1}});
    const auto cheap=evaluate_demand(h,v({0}),v({2}),EvaluationOptions{weighted({1,0,0}),nullptr,nullptr,all,{},Amount{2}});
    EXPECT_EQ(cheap.delivered,2);EXPECT_EQ(cheap.total_cost,(Counter128{0,4}));
    EXPECT_EQ(cheap.edge_flows,(std::optional<std::vector<Amount>>{{2,0,2}}));
  }
  std::mt19937 rng(20261022);
  for(int trial=0;trial<80;++trial) {
    const std::vector<Amount> w{static_cast<Amount>(rng()%4),static_cast<Amount>(rng()%4)};
    const auto g=std::make_shared<const Graph>(make_graph(4,{{0,2,static_cast<Amount>(rng()%3),1},
      {0,3,static_cast<Amount>(rng()%3),7},{1,2,static_cast<Amount>(rng()%3),3},
      {1,3,static_cast<Amount>(rng()%3),1},{2,3,static_cast<Amount>(rng()%3),1}}));
    const auto u=g->edge_sources(),d=g->edge_destinations();const auto cap=g->edge_capacity_units(),cost=g->edge_costs();
    std::vector<Amount> flow(5);std::map<Amount,Amount> best{{0,0}};
    std::function<void(std::size_t)> enumerate=[&](std::size_t e) {
      if(e<flow.size()){for(Amount x=0;x<=cap[e];++x){flow[e]=x;enumerate(e+1);}return;}
      Amount bal[4]={0,0,0,0};Amount c=0;
      for(std::size_t j=0;j<flow.size();++j){bal[u[j]]+=flow[j];bal[d[j]]-=flow[j];c+=flow[j]*cost[j];}
      if(bal[2]!=0)return;
      const Amount t=bal[0]+bal[1];const auto share=sequential(t,w);
      if(bal[0]!=share[0]||bal[1]!=share[1])return;
      auto [it,inserted]=best.emplace(t,c);if(!inserted)it->second=std::min(it->second,c);
    };
    enumerate(0);
    const Amount maximum=best.rbegin()->first;
    for(Amount offer:{Amount{0},Amount{1},Amount{9}}) {
      const auto r=evaluate_demand(*g,v({0,1}),v({3}),EvaluationOptions{weighted({w[0],w[1],0,0}),nullptr,nullptr,all,{},offer});
      const Amount t=std::min(offer,maximum);
      EXPECT_EQ(r.delivered,t);EXPECT_EQ(r.total_cost,(Counter128{0,static_cast<std::uint64_t>(best.at(t))}));
      EXPECT_EQ(sources(r),sequential(t,w));
      Amount bal[4]={0,0,0,0};
      for(std::size_t e=0;e<flow.size();++e){const Amount x=(*r.edge_flows)[e];EXPECT_GE(x,0);EXPECT_LE(x,cap[e]);bal[u[e]]+=x;bal[d[e]]-=x;}
      EXPECT_EQ(bal[2],0);EXPECT_EQ(bal[0],sources(r)[0]);EXPECT_EQ(bal[1],sources(r)[1]);EXPECT_EQ(bal[3],-t);
    }
  }
}

namespace {
// Compare two nonnegative rational numbers by Euclidean division, without any
// cross-product or production arithmetic; denominators can use 127 bits.
int rational_order(UInt128 a,UInt128 b,UInt128 c,UInt128 d) {
  bool reverse=false;
  for(;;) {
    const auto x=a/b,y=c/d;
    if(x!=y)return (x<y ? -1 : 1)*(reverse ? -1 : 1);
    a=a%b;c=c%d;
    if(a==0||c==0)return (a==c ? 0 : a==0 ? -1 : 1)*(reverse ? -1 : 1);
    std::swap(a,b);std::swap(c,d);reverse=!reverse;
  }
}
UInt128 test_value(detail::U128 a) {return (UInt128(a.hi)<<64)|UInt128(a.lo);}
}

TEST(WeightedSources, WideApportionmentHasExactQuotientOrderAndIncrementCarries) {
  for(int count:{1,2,3,4,5,16,17,40})for(bool unequal:{false,true}) {
    std::vector<Amount> weights(static_cast<std::size_t>(count));
    for(int i=0;i<count;++i)weights[static_cast<std::size_t>(i)]=kDomainLimit-1-(unequal ? i : 0);
    for(detail::U128 x:{detail::U128{0,16},detail::U128{1,1},detail::U128{(std::uint64_t{1}<<62)+1,3}}) {
      std::vector<detail::U128> a(weights.size());detail::highest_averages(x,weights,a);
      UInt128 sum=0;for(auto z:a)sum+=test_value(z);EXPECT_EQ(sum,test_value(x));
      for(std::size_t i=0;i<a.size();++i)for(std::size_t j=0;j<a.size();++j)if(test_value(a[i])!=0) {
        const int cmp=rational_order(UInt128(weights[i]),test_value(a[i]),UInt128(weights[j]),test_value(a[j])+1);
        EXPECT_TRUE(cmp>0||(cmp==0&&i<=j));
      }
    }
    std::vector<InEdge> edges;
    std::vector<NodeId> heads;
    for(int i=0;i<count;++i){edges.push_back({i,count,5,1});heads.push_back(i);}
    const auto g=std::make_shared<const Graph>(make_graph(count+1,edges));
    auto declared=weights;declared.push_back(0);
    const auto program=evaluate_demand(*g,heads,v({count}),EvaluationOptions{weighted(declared)}).program;
    const auto view=detail::make_forwarding_view(*g,program,nullptr,v({count}),nullptr);
    for(Counter128 x:{Counter128{1,1},Counter128{(std::uint64_t{1}<<62)+1,3}}) {
      std::vector<Amount> delta(weights.size());detail::source_offer_increments(view,x,17,delta);
      EXPECT_EQ(std::accumulate(delta.begin(),delta.end(),Amount{0}),17);
      if(count==5&&!unequal&&x.hi==1)EXPECT_EQ(delta,(std::vector<Amount>{3,3,4,4,3}));
      if(count==5&&unequal&&x.hi>1)EXPECT_EQ(delta,(std::vector<Amount>{4,3,3,4,3}));
    }
    FlowLedger ledger(g);
    FlowPolicy p(ledger,heads,v({count}),{.routing=Pinned{program},.admission=Admission::Lossy});
    UInt128 total=0;
    for(Amount offer:{kDomainLimit-1,kDomainLimit-1,kDomainLimit-1,Amount{7},kDomainLimit-1,Amount{5},Amount{17}}) {
      const auto r=p.place(offer);
      total+=UInt128(offer);
      ASSERT_TRUE(r.requested);
      EXPECT_EQ((UInt128(r.requested->hi)<<64)|UInt128(r.requested->lo),total);
      UInt128 offered=0;for(const auto& row:r.by_source){ASSERT_TRUE(row.offered);offered+=(UInt128(row.offered->hi)<<64)|UInt128(row.offered->lo);}
      EXPECT_EQ(offered,total);
    }
  }
}

TEST(WeightedSources, WideCostKeepsBothWordsAndSharedInputOwnership) {
  const Amount m=kDomainLimit-1;
  const auto g=std::make_shared<const Graph>(make_graph(2,{{0,1,m,m}}));
  auto owner=std::make_shared<const std::vector<Amount>>(std::vector<Amount>{m,0});
  std::weak_ptr<const std::vector<Amount>> weak=owner;
  CentralizedTE cfg(SourceAllocation::Weighted);cfg.source_weights=owner;
  FlowLedger ledger(g);std::optional<FlowPolicy> p;
  p.emplace(ledger,v({0}),v({1}),FlowPolicyConfig{cfg});
  cfg.source_weights.reset();owner.reset();ASSERT_FALSE(weak.expired());
  const auto r=p->replace(m,{},all);
  const UInt128 square=UInt128(m)*UInt128(m);
  ASSERT_TRUE(r.total_cost);
  EXPECT_EQ((UInt128(r.total_cost->hi)<<64)|UInt128(r.total_cost->lo),square);
  EXPECT_EQ(r.delivered,m);EXPECT_EQ(sources(r),(std::vector<Amount>{m}));
  auto& requested=detail::PolicyAccess::engine(*p).epoch.requested;
  requested=Counter128{(std::uint64_t{1}<<63)-1,std::numeric_limits<std::uint64_t>::max()};
  const auto before=ledger_state(ledger);const auto old=epoch(*p);
  EXPECT_EQ(refused([&]{(void)p->place(1);}),"requested limit exceeded");
  expect_same(ledger_state(ledger),before,"counter limit");expect_same(epoch(*p),old,"counter limit");
  p.reset();EXPECT_TRUE(weak.expired());
}

TEST(WeightedSources, MasksChangeTheRosterButBlackHolesKeepDeclaredShares) {
  const auto g=std::make_shared<const Graph>(make_graph(3,{{0,2,10,1},{1,2,10,1}}));
  const auto program=evaluate_demand(*g,v({0,1}),v({2}),EvaluationOptions{weighted({2,1,0}),nullptr,nullptr,{}, {},Amount{3}}).program;
  const auto mask=g->eligibility(v({0}),{});
  FlowLedger a(g),b(g);
  FlowPolicy repaired(a,v({0,1}),v({2}),{.routing=Pinned{program},.admission=Admission::Lossy},&mask);
  FlowPolicy stale(b,v({0,1}),v({2}),{.routing=Pinned{program,std::make_shared<const Eligibility>(mask)},.admission=Admission::Lossy});
  EXPECT_EQ(sources(repaired.replace(3)),(std::vector<Amount>{0,3}));
  const auto r=stale.replace(3);
  EXPECT_EQ(sources(r),(std::vector<Amount>{0,1}));
  EXPECT_EQ(r.source_drops,(Counter128{0,2}));
  const auto supported=evaluate_demand(*g,v({0,1}),v({2}),EvaluationOptions{weighted({2,1,0},dag_of(*g,{1})),nullptr,nullptr,all});
  EXPECT_EQ(sources(supported),(std::vector<Amount>{0,10}));
  EXPECT_EQ(supported.edge_flows,(std::optional<std::vector<Amount>>{{0,10}}));
}

TEST(WeightedSources, ZeroInjectionChangesNoOtherLossOrForwardingStructure) {
  const auto g=witness();
  const auto empty=evaluate_demand(*g,v({0,1}),v({5}),EvaluationOptions{weighted({0,0,0,0,0,0})}).program;
  const auto positive=evaluate_demand(*g,v({0,1}),v({5}),EvaluationOptions{weighted({2,1,0,0,0,0})}).program;
  const auto dag=dag_of(*g,{0,1,2,3,4,5,6});
  const auto pooled=ForwardingProgram::from_dag(*g,v({0,1}),v({5}),dag);
  const std::vector<RouteDefinition> definitions{{0,dag,1},{1,dag,1}};
  const auto split=ForwardingProgram::from_routes(*g,v({0,1}),v({5}),definitions);
  const std::vector<Routing> models{Pinned{positive},ShortestPaths{},Pinned{pooled},Pinned{split}};
  for(const auto& model:models)for(Amount zero_offer:{Amount{0},Amount{7}}) {
    FlowLedger ledger(g);
    FlowPolicy p(ledger,v({0,1}),v({5}),FlowPolicyConfig{model,Admission::Lossy});
    FlowPolicy zero(ledger,v({0,1}),v({5}),{.routing=Pinned{empty},.admission=Admission::Lossy});
    const std::vector<Action> baseline{Replace{p,19}};
    const auto want=apply_actions(ledger,baseline,BatchMode::Simultaneous,all).front();
#if NETGRAPH_CORE_PRIVATE_HOOKS
    const auto old=detail::thread_lossy_counters();
#endif
    p.remove();
    const std::vector<Action> with_zero{Replace{zero,zero_offer},Replace{p,19}};
    const auto actual=apply_actions(ledger,with_zero,BatchMode::Simultaneous,all);
    expect_same(actual[1],want,"zero source injection");
#if NETGRAPH_CORE_PRIVATE_HOOKS
    const auto& now=detail::thread_lossy_counters();
    EXPECT_EQ(now.groups,old.groups);EXPECT_EQ(now.rounds,old.rounds);EXPECT_EQ(now.barriers,old.barriers);
    EXPECT_EQ(now.dependencies,old.dependencies);EXPECT_EQ(now.event_structures,old.event_structures);
#endif
  }
}

TEST(WeightedSources, TransitThroughAnotherSourceIsNotItsInjection) {
  const auto g=std::make_shared<const Graph>(make_graph(3,{{0,1,5,1},{1,2,10,1}}));
  FlowLedger ledger(g);FlowPolicy p(ledger,v({0,1}),v({2}),weighted({1,1,0}));
  const auto r=p.replace(std::nullopt,{},all);
  EXPECT_EQ(r.delivered,10);EXPECT_EQ(r.total_cost,(Counter128{0,15}));
  EXPECT_EQ(sources(r),(std::vector<Amount>{5,5}));
  EXPECT_EQ(r.edge_flows,(std::optional<std::vector<Amount>>{{5,10}}));
  const auto records=p.records();ASSERT_EQ(records.size(),2U);
  EXPECT_EQ(records[0].delivered,5);EXPECT_EQ(records[1].delivered,5);
}

TEST(WeightedSources, LiveTunnelRowsAndUndeliveredAreExactAcrossReaders) {
  const auto g=std::make_shared<const Graph>(make_graph(2,{{0,1,5,1}}));
  const auto cfg=weighted({1,0});
  const auto check=[](const PlacementResult& r,Amount offer) {
    const Amount carried=std::min(offer,Amount{5});
    const Counter128 amount{0,static_cast<std::uint64_t>(carried)};
    EXPECT_EQ(r.call_delivered,carried);EXPECT_EQ(r.call_undelivered,offer-carried);
    EXPECT_EQ(r.delivery_delta,carried);EXPECT_TRUE(r.call_edge_drops.empty());EXPECT_EQ(r.call_source_drops,0);
    EXPECT_EQ(r.termination,offer<=5 ? Termination::Complete : Termination::Exhausted);
    EXPECT_TRUE(r.routes_changed);EXPECT_TRUE(r.weights_changed);EXPECT_FALSE(r.reservations_changed);
    ASSERT_TRUE(r.program);EXPECT_FALSE(r.plan);EXPECT_FALSE(r.offered);
    EXPECT_EQ(r.requested,(Counter128{0,static_cast<std::uint64_t>(offer)}));
    EXPECT_EQ(r.admitted,amount);EXPECT_EQ(r.delivered,carried);
    EXPECT_EQ(r.unadmitted,(Counter128{0,static_cast<std::uint64_t>(offer-carried)}));
    EXPECT_EQ(r.source_drops,(Counter128{}));EXPECT_TRUE(r.edge_drops.empty());EXPECT_EQ(r.total_cost,amount);
    EXPECT_EQ(r.by_source,(std::vector<SourceFlow>{{0,{},amount,carried,{},{}}}));
    EXPECT_FALSE(r.min_cut);EXPECT_FALSE(r.marginal_cost_curve);
    EXPECT_EQ(r.edge_flows,(std::optional<std::vector<Amount>>{{carried}}));
    EXPECT_EQ(r.residual,(std::optional<std::vector<Amount>>{{5-carried}}));
    EXPECT_EQ(r.edge_reservations,(std::optional<std::vector<Amount>>{{0}}));
    EXPECT_EQ(r.saturated_boundary_edges,(std::optional<std::vector<EdgeId>>{carried==5 ? v({0}) : v({})}));
    EXPECT_EQ(r.by_tunnel,(std::vector<TunnelFlow>{{0,0,TunnelStatus::Up,{0,static_cast<std::uint64_t>(carried)},carried,{}}}));
  };
  for(Amount offer:{Amount{3},Amount{7}}) {
    for(int operation=0;operation<3;++operation)for(bool ordered:{false,true}) {
      FlowLedger ledger(g);FlowPolicy p(ledger,v({0}),v({1}),cfg);
      if(operation==2)check(p.evaluate(offer,{},all),offer);
      else if(!ordered)check(operation==0 ? p.place(offer,all) : p.replace(offer,{},all),offer);
      else {
        const std::vector<Action> actions=operation==0 ? std::vector<Action>{Place{p,offer}} : std::vector<Action>{Replace{p,offer}};
        check(apply_actions(ledger,actions,BatchMode::Ordered,all).front(),offer);
      }
    }
    const EvaluationOptions options{cfg,nullptr,nullptr,all,{},offer};
    check(evaluate_demand(*g,v({0}),v({1}),options),offer);
    for(int workers:{1,2,4}) {
      const std::vector<std::pair<std::vector<NodeId>,std::vector<NodeId>>> pairs(5,{{0},{1}});
      for(const auto& r:evaluate_pairs(*g,pairs,options,{},workers))check(r,offer);
      check(evaluate_edge_failures(*g,v({0}),v({1}),options,workers).baseline,offer);
    }
  }
}

TEST(WeightedSources, RemainingZeroCapturedWeightsWithAProgramPathHaveNoActiveWeight) {
  const auto g=std::make_shared<const Graph>(make_graph(3,{{0,2,10,1},{1,2,10,1}}));
  const auto program=evaluate_demand(*g,v({0,1}),v({2}),EvaluationOptions{weighted({2,0,0}),nullptr,nullptr,all,{},Amount{5}}).program;
  ASSERT_EQ(program->tunnels().size(),1U);
  const auto mask=g->eligibility(v({0}),{});
  for(Admission admission:{Admission::Lossy,Admission::Lossless}) {
    FlowLedger ledger(g);FlowPolicy p(ledger,v({0,1}),v({2}),{.routing=Pinned{program},.admission=admission},&mask);
    for(int operation=0;operation<3;++operation) {
      p.remove();
      const auto r=operation==0 ? p.replace(4,{},all) : operation==1 ? p.evaluate(4,{},all) : p.place(4,all);
      EXPECT_EQ(r.termination,Termination::NoActiveWeight);EXPECT_EQ(r.delivered,0);
      EXPECT_EQ(r.call_undelivered,4);EXPECT_EQ(r.source_drops,(Counter128{0,admission==Admission::Lossy ? 4U : 0U}));
      EXPECT_EQ(r.unadmitted,(std::optional<Counter128>{{0,admission==Admission::Lossy ? 0U : 4U}}));
      EXPECT_EQ(sources(r),(std::vector<Amount>{0,0}));
      for(const auto& row:r.by_source) EXPECT_EQ(row.source_drops,(Counter128{}));
    }
  }
}

TEST(WeightedSources, SmallWeightSumsUseWideCumulativeArithmetic) {
  const Amount m=kDomainLimit-1;
  for(int count:{2,3,4})for(bool unequal:{false,true}) {
    std::vector<InEdge> edges;std::vector<NodeId> heads;std::vector<Amount> weights;
    for(int i=0;i<count;++i){edges.push_back({i,count,5,1});heads.push_back(i);weights.push_back(m-(unequal ? i : 0));}
    weights.push_back(0);
    const auto g=std::make_shared<const Graph>(make_graph(count+1,edges));
    const auto program=evaluate_demand(*g,heads,v({count}),EvaluationOptions{weighted(weights)}).program;
    const auto view=detail::make_forwarding_view(*g,program,nullptr,v({count}),nullptr);
    // The sum of four M weights is 2^64-4; X=3M+7 reaches the old product overflow.
    std::vector<Amount> delta(static_cast<std::size_t>(count));
    detail::source_offer_increments(view,{0,3*static_cast<std::uint64_t>(m)+7},m,delta);
    EXPECT_TRUE(std::all_of(delta.begin(),delta.end(),[](Amount a){return a>=0;}));
    UInt128 sum=0;for(auto a:delta)sum+=UInt128(a);EXPECT_EQ(sum,UInt128(m));
    if(count==4&&!unequal) EXPECT_EQ(delta,(std::vector<Amount>{m/4+1,m/4+1,m/4+1,m/4}));
  }
}
