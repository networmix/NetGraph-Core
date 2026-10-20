// apply_actions, natively: the list's validation order (every refusal before anything is computed,
// nothing changed), ORDERED as the scalar calls on a private working state (earlier actions first, FlowIds in
// list order, outputs and holds after each action, normal outcomes, zero places, an abort
// mid-list), SIMULTANEOUS on one view with the loss rule across demands (the view, preparation on
// it, permutation, creation-order ties, pooled sources, the demand level, each participant's own
// black holes, T1 holds, empty rosters, zero offers, repeats, orphans, iteration errors), link sums
// beyond 2^63, counters beyond 2^64, per-entry outputs and every field of an entry. Hand-derived
// values on small graphs (EdgeIds equal input positions unless mapped); the generators compare an
// ORDERED list with the same scalar calls and a SIMULTANEOUS participant with its scalar call.
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "netgraph/core/detail/forwarding.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/forwarding_program.hpp"
#include "state_comparators.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using detail::LedgerAccess;

namespace {

using Amounts = std::vector<Amount>;
using Pairs = std::vector<std::pair<EdgeId, Amount>>;
using Actions = std::vector<Action>;

constexpr Admission LOSSY = Admission::Lossy;
constexpr Admission LOSSLESS = Admission::Lossless;
constexpr BatchMode ORDERED = BatchMode::Ordered;
constexpr BatchMode SIMULTANEOUS = BatchMode::Simultaneous;
const std::optional<Amount> INF;

const std::string kPeriodTwo = "lossy forwarding iteration entered a period-2 cycle";
const std::string kRoundLimit = "lossy forwarding iteration did not reach a fixed point within the round limit";
const std::string kOffer = "offer must be an integer in [0, 2^62) or permitted +inf";

Amount lo(const Counter128& x) {
  EXPECT_EQ(x.hi, 0U);
  return static_cast<Amount>(x.lo);
}

std::shared_ptr<const Graph> graph(std::int32_t n, const std::vector<InEdge>& rows) {
  auto g = shared(make_graph(n, rows));
  for (std::size_t i = 0; i < rows.size(); ++i) EXPECT_EQ(g->input_positions()[i], static_cast<std::int32_t>(i));
  return g;
}

FlowPolicyConfig ecmp(Admission admission = LOSSY) { return FlowPolicyConfig {.routing = shortest_paths(), .admission = admission}; }
FlowPolicyConfig t1(Admission admission = LOSSY) { return FlowPolicyConfig {.routing = DistributedTE {TunnelCount {1}}, .admission = admission}; }
FlowPolicyConfig t2(Admission admission = LOSSY) {
  return FlowPolicyConfig {.routing = centralized(SourceAllocation::Equal, nullptr, TunnelCount {1}), .admission = admission};
}
FlowPolicyConfig global(Admission admission = LOSSLESS) { return FlowPolicyConfig {.routing = centralized(), .admission = admission}; }
FlowPolicyConfig pinned(std::shared_ptr<const ForwardingProgram> prog, Admission admission = LOSSY,
                        std::shared_ptr<const Eligibility> black_holes = nullptr) {
  return FlowPolicyConfig {.routing = Pinned {std::move(prog), std::move(black_holes)}, .admission = admission};
}

// A supplied program of tunnels: (head, EdgeIds in path order, weight).
std::shared_ptr<const ForwardingProgram> tunnels(const Graph& g, const std::vector<NodeId>& s, const std::vector<NodeId>& d,
                                                 const std::vector<std::tuple<NodeId, std::vector<EdgeId>, Amount>>& specs) {
  std::vector<RouteDefinition> routes;
  for (const auto& [h, es, w] : specs) {
    std::vector<std::vector<EdgeId>> hops;
    for (EdgeId e : es) hops.push_back({e});
    routes.push_back(RouteDefinition {h, std::make_shared<const Path>(g, hops), w});
  }
  return ForwardingProgram::from_routes(g, s, d, routes);
}
std::shared_ptr<const ForwardingProgram> tunnel(const Graph& g, NodeId s, NodeId d, const std::vector<EdgeId>& path, Amount w = 1) {
  return tunnels(g, {s}, {d}, {{s, path, w}});
}

std::shared_ptr<const Eligibility> no_edges(const Graph& g, const std::vector<EdgeId>& edges) {
  return std::make_shared<const Eligibility>(g.eligibility({}, edges));
}
std::shared_ptr<const Eligibility> no_nodes(const Graph& g, const std::vector<NodeId>& nodes) {
  return std::make_shared<const Eligibility>(g.eligibility(nodes, {}));
}

// A ledger with its policies (created in order) and what they point to (declared first).
struct World {
  std::shared_ptr<const Graph> g;
  std::unique_ptr<FlowLedger> ledger;
  std::vector<std::shared_ptr<const Eligibility>> keep;
  std::vector<std::unique_ptr<FlowPolicy>> p;
  explicit World(std::shared_ptr<const Graph> graph) : g(std::move(graph)), ledger(std::make_unique<FlowLedger>(g)) {}
  World(std::shared_ptr<const Graph> graph, std::vector<Amount> service, std::vector<Amount> reservation)
      : g(std::move(graph)), ledger(std::make_unique<FlowLedger>(g, service, reservation)) {}
  FlowPolicy& add(const std::vector<NodeId>& s, const std::vector<NodeId>& d, FlowPolicyConfig cfg,
                  std::shared_ptr<const Eligibility> elig = nullptr) {
    keep.push_back(elig);
    p.push_back(std::make_unique<FlowPolicy>(*ledger, s, d, std::move(cfg), elig.get()));
    return *p.back();
  }
  FlowPolicy& pin(const std::shared_ptr<const ForwardingProgram>& prog, Admission admission = LOSSY,
                  std::shared_ptr<const Eligibility> black_holes = nullptr, std::shared_ptr<const Eligibility> elig = nullptr) {
    keep.push_back(black_holes);
    return add(vec(prog->sources()), vec(prog->destinations()), pinned(prog, admission, black_holes), std::move(elig));
  }
  // The ledger and every policy's epoch.
  struct State {
    LedgerState ledger;
    std::vector<Epoch> epochs;
  };
  [[nodiscard]] State state() const {
    State s {ledger_state(*ledger), {}};
    for (const auto& x : p) s.epochs.push_back(epoch(*x));
    return s;
  }
};

void expect_unchanged(const World& w, const World::State& before, const std::string& where) {
  expect_same(ledger_state(*w.ledger), before.ledger, where);
  ASSERT_EQ(w.p.size(), before.epochs.size());
  for (std::size_t i = 0; i < w.p.size(); ++i) expect_unchanged(epoch(*w.p[i]), before.epochs[i], where + ", policy " + std::to_string(i));
}

// The message of the std::invalid_argument the list raises ("" for none), and that nothing changed.
std::string refused(World& w, const Actions& actions, BatchMode mode, const OutputOptions& outputs = {}) {
  const auto before = w.state();
  const std::string what = error_of([&] { (void)apply_actions(*w.ledger, actions, mode, outputs); });
  expect_unchanged(w, before, "refused: " + what);
  return what;
}

std::vector<PlacementResult> run(World& w, const Actions& actions, BatchMode mode, const OutputOptions& outputs = {}) {
  auto out = apply_actions(*w.ledger, actions, mode, outputs);
  expect_balanced(*w.ledger, "after the list");
  return out;
}

Amounts delivered(const std::vector<PlacementResult>& rs) {
  Amounts out;
  for (const auto& r : rs) out.push_back(r.call_delivered);
  return out;
}

const OutputOptions kAll {.edge_flows = true, .residual = true, .saturated_boundary_edges = true, .edge_reservations = true};

// The step-3 cycle Pc: CYC(1, 2), two tunnels crossing X (e0) and Y (e3) in opposite orders,
// weights 1 and 2; its replace(3) enters a period-2 cycle.
std::shared_ptr<const Graph> cyc12() {
  return graph(6, {{0, 1, 1, 1}, {1, 2, 1000, 1}, {1, 5, 1000, 1}, {2, 3, 2, 1}, {3, 0, 1000, 1}, {3, 5, 1000, 1}, {4, 0, 1000, 1},
                   {4, 2, 1000, 1}});
}
std::shared_ptr<const ForwardingProgram> pc_program(const Graph& g) {
  return tunnels(g, {4}, {5}, {{4, {6, 0, 1, 3, 5}, 1}, {4, {7, 3, 4, 0, 2}, 2}});
}

// ---------------------------------------------------------------- validation

// An empty list returns nothing in both modes; a mode outside BatchMode is refused even then.
TEST(ApplyActions, EmptyList) {
  World w(graph(2, {{0, 1, 10, 1}}));
  auto& e = w.add({0}, {1}, ecmp());
  (void)e.place(3);
  const auto before = w.state();
  EXPECT_TRUE(run(w, {}, ORDERED).empty());
  EXPECT_TRUE(run(w, {}, SIMULTANEOUS).empty());
  expect_unchanged(w, before, "empty");
  EXPECT_EQ(refused(w, {}, static_cast<BatchMode>(7)), "mode: 7 is not a value of BatchMode");
  EXPECT_EQ(refused(w, {Place {e, 1}}, static_cast<BatchMode>(7)), "mode: 7 is not a value of BatchMode");
}

// Ownership, then the duplicate rule, per action in list order.
TEST(ApplyActions, Duplicate) {
  World w(graph(2, {{0, 1, 10, 1}}));
  auto& e = w.add({0}, {1}, ecmp());
  auto& e2 = w.add({0}, {1}, ecmp());
  World other(graph(2, {{0, 1, 10, 1}}));
  auto& x = other.add({0}, {1}, ecmp());
  (void)e.place(2);
  for (BatchMode mode : {ORDERED}) {
    EXPECT_EQ(refused(w, {Place {e, 1}, Place {e, 2}}, mode), "duplicate policy in list");
    EXPECT_EQ(refused(w, {Place {x, 1}}, mode), "list policy belongs to another ledger");
    EXPECT_EQ(refused(w, {Place {x, 1}, Place {x, 1}}, mode), "list policy belongs to another ledger");
    EXPECT_EQ(refused(w, {Place {e, 1}, Place {e, -1}}, mode), "duplicate policy in list");
    EXPECT_EQ(refused(w, {Place {e, -1}, Place {e2, 1}}, mode), kOffer);
    EXPECT_EQ(refused(w, {Place {e, -1}, Place {x, 1}}, mode), kOffer);
  }
  EXPECT_EQ(refused(w, {Replace {e, 1}, Forward {e, 1}}, SIMULTANEOUS), "duplicate policy in list");
  EXPECT_EQ(refused(w, {Replace {e, 1}, Forward {e, 1}}, ORDERED), "duplicate policy in list");
}

// The SIMULTANEOUS rules in their order (kind, the Replace model, admission,
// finiteness), then the scalar call's checks.
TEST(ApplyActions, SimultaneousRefusals) {
  World w(graph(2, {{0, 1, 10, 1}}));
  auto& e = w.add({0}, {1}, ecmp());
  auto& e2 = w.add({0}, {1}, ecmp());
  auto& f = w.add({0}, {1}, ecmp(LOSSLESS));
  auto& t = w.add({0}, {1}, t1());
  auto& k = w.add({0}, {1}, global());
  auto& opt = w.add({0}, {1}, FlowPolicyConfig {.routing = shortest_paths(RouteAllocation::Optimized), .admission = LOSSLESS});
  auto& tt = w.add({0}, {1}, t2());
  auto& p = w.pin(tunnel(*w.g, 0, 1, {0}));
  auto& pll = w.pin(tunnel(*w.g, 0, 1, {0}), LOSSLESS);
  (void)k.place(1);
  EXPECT_EQ(refused(w, {Place {e, 1}}, SIMULTANEOUS), "SIMULTANEOUS does not accept Place");
  const std::string model = "SIMULTANEOUS Replace requires ShortestPaths EQUAL or WEIGHTED, or Pinned";
  EXPECT_EQ(refused(w, {Replace {t, 1}}, SIMULTANEOUS), model);
  EXPECT_EQ(refused(w, {Replace {k, 1}}, SIMULTANEOUS), model);
  EXPECT_EQ(refused(w, {Replace {opt, 1}}, SIMULTANEOUS), model);
  EXPECT_EQ(refused(w, {Replace {tt, 1}}, SIMULTANEOUS), model);
  const std::string lossy = "SIMULTANEOUS requires admission=LOSSY";
  EXPECT_EQ(refused(w, {Replace {f, 1}}, SIMULTANEOUS), lossy);
  EXPECT_EQ(refused(w, {Replace {pll, 1}}, SIMULTANEOUS), lossy);
  EXPECT_EQ(refused(w, {Forward {k, 1}}, SIMULTANEOUS), lossy);
  EXPECT_EQ(refused(w, {Replace {e, INF}}, SIMULTANEOUS), "SIMULTANEOUS requires finite offers");
  EXPECT_EQ(refused(w, {Forward {p, 1}}, SIMULTANEOUS), "Pinned forwards with replace");
  EXPECT_EQ(refused(w, {Forward {e2, 1}}, SIMULTANEOUS), "forward requires an installed program");
  EXPECT_EQ(refused(w, {Replace {e, -1}}, SIMULTANEOUS), kOffer);
  EXPECT_EQ(refused(w, {Replace {e, 1, std::optional<Amount> {1}}}, SIMULTANEOUS), "ShortestPaths does not accept design_demand");
  // Several invalid parts in one action.
  EXPECT_EQ(refused(w, {Forward {e2, INF}}, SIMULTANEOUS), "SIMULTANEOUS requires finite offers");
  EXPECT_EQ(refused(w, {Forward {k, INF}}, SIMULTANEOUS), lossy);
  EXPECT_EQ(refused(w, {Replace {t, -1}}, SIMULTANEOUS), model);
}

// ORDERED runs the scalar checks.
TEST(ApplyActions, OrderedScalarChecks) {
  World w(graph(2, {{0, 1, 10, 1}}));
  World other(graph(2, {{0, 1, 10, 1}}));
  auto& e = w.add({0}, {1}, ecmp());
  auto& t = w.add({0}, {1}, t1());
  auto& tt = w.add({0}, {1}, t2());
  auto& inf = w.add({0}, {1}, ecmp());
  (void)inf.place(INF);
  (void)e.place(1);
  EXPECT_EQ(refused(w, {Place {t, 1}}, ORDERED), "DistributedTE does not support incremental placement");
  EXPECT_EQ(refused(w, {Place {tt, 0}}, ORDERED), "CentralizedTE with tunnel_count does not support incremental placement");
  const auto foreign = no_edges(*other.g, {0});
  EXPECT_EQ(refused(w, {Forward {e, 1, foreign.get()}}, ORDERED), "black_holes belongs to another graph");
  EXPECT_EQ(refused(w, {Place {inf, 1}}, ORDERED), "finite placement after infinite placement requires a new epoch");
  EXPECT_EQ(refused(w, {Place {t, -1}}, ORDERED), "DistributedTE does not support incremental placement");
}

// An invalid last action of a long list: refused before anything is computed (not the
// period-2 error the first action's computation raises).
TEST(ApplyActions, InvalidLastActionChangesNothing) {
  World w(cyc12());
  auto& pc = w.pin(pc_program(*w.g));
  auto& t = w.add({4}, {5}, t1());
  std::vector<FlowPolicy*> fresh;
  for (int i = 0; i < 48; ++i) fresh.push_back(&w.add({4}, {5}, ecmp()));
  auto& e = w.add({4}, {5}, ecmp());
  Actions ordered {Replace {pc, 3}};
  for (auto* f : fresh) ordered.push_back(Place {*f, 1});
  ordered.push_back(Place {t, 1});
  EXPECT_EQ(refused(w, ordered, ORDERED), "DistributedTE does not support incremental placement");
  Actions simultaneous {Replace {pc, 3}};
  for (auto* f : fresh) simultaneous.push_back(Replace {*f, 1});
  simultaneous.push_back(Place {e, 1});
  EXPECT_EQ(refused(w, simultaneous, SIMULTANEOUS), "SIMULTANEOUS does not accept Place");
  // Controls: the same lists without their last action reach the computation's error.
  ordered.pop_back();
  simultaneous.pop_back();
  EXPECT_EQ(refused(w, ordered, ORDERED), kPeriodTwo);
  EXPECT_EQ(refused(w, simultaneous, SIMULTANEOUS), kPeriodTwo);
}

// ---------------------------------------------------------------- ORDERED

// Earlier actions are served first; a reservation protects nothing.
TEST(ApplyActions, OrderIsPriority) {
  for (bool tunnel_first : {false, true}) {
    World w(graph(2, {{0, 1, 10, 1}}));
    auto& e = w.add({0}, {1}, ecmp());
    auto& t = w.add({0}, {1}, t1());
    const Actions list = tunnel_first ? Actions {Replace {t, 10, std::optional<Amount> {10}}, Replace {e, 10}}
                                      : Actions {Replace {e, 10}, Replace {t, 10, std::optional<Amount> {10}}};
    const auto r = run(w, list, ORDERED);
    const auto& er = tunnel_first ? r[1] : r[0];
    const auto& tr = tunnel_first ? r[0] : r[1];
    if (!tunnel_first) {
      EXPECT_EQ(er.call_delivered, 10);
      EXPECT_EQ(tr.call_delivered, 0);
      EXPECT_EQ(tr.call_edge_drops, (Pairs {{0, 10}}));
      ASSERT_EQ(tr.by_tunnel.size(), 1U);
      EXPECT_EQ(tr.by_tunnel[0].status, TunnelStatus::Up);
    } else {
      EXPECT_EQ(tr.call_delivered, 10);
      EXPECT_EQ(er.call_delivered, 0);
      EXPECT_EQ(er.call_edge_drops, (Pairs {{0, 10}}));
    }
    EXPECT_EQ(w.ledger->carried(), (Amounts {10}));
    EXPECT_EQ(w.ledger->reserved(), (Amounts {10}));
    EXPECT_EQ(w.ledger->residual(), (Amounts {0}));
    EXPECT_EQ(w.ledger->unreserved(), (Amounts {0}));
  }
}

// New FlowIds follow list order (as the scalar sequence B, A), not creation order.
TEST(ApplyActions, FlowIdsInListOrder) {
  World w(graph(2, {{0, 1, 10, 1}}));
  auto& a = w.add({0}, {1}, ecmp());
  auto& b = w.add({0}, {1}, ecmp());
  (void)run(w, {Replace {b, 1}, Replace {a, 1}}, ORDERED);
  ASSERT_EQ(b.records().size(), 1U);
  ASSERT_EQ(a.records().size(), 1U);
  EXPECT_EQ(b.records()[0].flow_id, (FlowId {w.ledger->uid(), 0}));
  EXPECT_EQ(a.records()[0].flow_id, (FlowId {w.ledger->uid(), 1}));
}

// LIST-outputs-ordered: each action's outputs read the state after it.
TEST(ApplyActions, OrderedOutputsAfterEachAction) {
  World w(graph(2, {{0, 1, 10, 1}}));
  auto& a = w.add({0}, {1}, ecmp());
  auto& b = w.add({0}, {1}, ecmp());
  const OutputOptions o {.edge_flows = true, .residual = true, .saturated_boundary_edges = true};
  const auto r = run(w, {Replace {a, 7}, Replace {b, 7}}, ORDERED, o);
  EXPECT_EQ(r[0].call_delivered, 7);
  EXPECT_EQ(*r[0].residual, (Amounts {3}));
  EXPECT_EQ(*r[0].saturated_boundary_edges, (std::vector<EdgeId> {}));
  EXPECT_EQ(*r[0].edge_flows, (Amounts {7}));
  EXPECT_EQ(r[1].call_delivered, 3);
  EXPECT_EQ(r[1].call_edge_drops, (Pairs {{0, 4}}));
  EXPECT_EQ(*r[1].residual, (Amounts {0}));
  EXPECT_EQ(*r[1].saturated_boundary_edges, (std::vector<EdgeId> {0}));
  EXPECT_EQ(*r[1].edge_flows, (Amounts {3}));
}

// Holds after each action: the second T1 sets up on the first's hold.
TEST(ApplyActions, OrderedHoldsAfterEachAction) {
  for (bool reversed : {false, true}) {
    World w(graph(2, {{0, 1, 10, 1}}));
    auto& t = w.add({0}, {1}, t1());
    auto& u = w.add({0}, {1}, t1());
    const OutputOptions o {.edge_reservations = true};
    const Replace rt {t, 10, std::optional<Amount> {10}};
    const Replace ru {u, 10, std::optional<Amount> {10}};
    const auto r = run(w, reversed ? Actions {ru, rt} : Actions {rt, ru}, ORDERED, o);
    const auto& first = r[0];
    const auto& second = r[1];
    ASSERT_EQ(first.by_tunnel.size(), 1U);
    EXPECT_EQ(first.by_tunnel[0].status, TunnelStatus::Up);
    EXPECT_EQ(first.call_delivered, 10);
    EXPECT_EQ(*first.edge_reservations, (Amounts {10}));
    ASSERT_EQ(second.by_tunnel.size(), 1U);
    EXPECT_EQ(second.by_tunnel[0].status, TunnelStatus::InsufficientBandwidth);
    EXPECT_EQ(second.call_source_drops, 10);
    EXPECT_EQ(second.call_delivered, 0);
    EXPECT_EQ(*second.edge_reservations, (Amounts {0}));
    EXPECT_EQ(w.ledger->reserved(), (Amounts {10}));
  }
}

// Normal outcomes commit; a zero place publishes nothing.
TEST(ApplyActions, NormalOutcomesCommit) {
  World w(graph(2, {{0, 1, 10, 1}}));
  auto& e = w.add({0}, {1}, ecmp());
  auto& f = w.add({0}, {1}, ecmp(LOSSLESS));
  auto& e2 = w.add({0}, {1}, ecmp());
  const auto r = run(w, {Replace {e, 10}, Replace {f, 5}}, ORDERED);
  EXPECT_EQ(r[0].call_delivered, 10);
  EXPECT_EQ(r[1].call_delivered, 0);
  EXPECT_EQ(r[1].termination, Termination::Exhausted);
  EXPECT_EQ(r[1].call_undelivered, 5);
  ASSERT_EQ(f.records().size(), 1U);
  EXPECT_EQ(f.records()[0].delivered, 0);
  const auto before = w.state();
  const auto z = run(w, {Place {e2, 0}}, ORDERED, OutputOptions {.residual = true});
  EXPECT_EQ(z[0].call_delivered, 0);
  EXPECT_EQ(z[0].program, nullptr);
  EXPECT_FALSE(z[0].routes_changed);
  EXPECT_EQ(*z[0].residual, (Amounts {0}));
  expect_unchanged(w, before, "a zero place");
}

// An error mid-list publishes nothing (A, on its own link 6 -> 7, keeps its old record and id).
TEST(ApplyActions, ErrorLeavesEverythingUnchanged) {
  World w(graph(8, {{0, 1, 1, 1}, {1, 2, 1000, 1}, {1, 5, 1000, 1}, {2, 3, 2, 1}, {3, 0, 1000, 1}, {3, 5, 1000, 1}, {4, 0, 1000, 1},
                    {4, 2, 1000, 1}, {6, 7, 10, 1}}));
  auto& a = w.add({6}, {7}, ecmp());
  auto& pc = w.pin(pc_program(*w.g));
  (void)a.place(2);
  const auto id = a.records()[0].flow_id;
  EXPECT_EQ(refused(w, {Replace {a, 5}, Replace {pc, 3}}, ORDERED), kPeriodTwo);
  EXPECT_EQ(a.delivered(), 2);
  EXPECT_EQ(a.records()[0].flow_id, id);
}

// Successful infinities in ORDERED lists and the first-hop experiment.
TEST(ApplyActions, InfiniteOrderedActions) {
  {
    World w(graph(2, {{0, 1, 10, 1}}));
    auto& f = w.add({0}, {1}, ecmp(LOSSLESS));
    auto& e = w.add({0}, {1}, ecmp());
    const auto r = run(w, {Replace {f, INF}, Replace {e, 1}}, ORDERED);
    EXPECT_EQ(r[0].call_delivered, 10);
    EXPECT_FALSE(r[0].requested.has_value());
    EXPECT_FALSE(r[0].offered.has_value());
    EXPECT_FALSE(r[0].unadmitted.has_value());
    EXPECT_FALSE(r[0].call_undelivered.has_value());
    EXPECT_EQ(r[0].termination, Termination::Exhausted);
    EXPECT_EQ(r[1].call_delivered, 0);
    EXPECT_EQ(r[1].call_edge_drops, (Pairs {{0, 1}}));
    EXPECT_EQ(r[1].termination, Termination::Complete);
  }
  const auto par = [] { return graph(2, {{0, 1, 4, 1}, {0, 1, 10, 1}}); };
  {
    World w(par());
    auto& e = w.add({0}, {1}, ecmp());
    const auto r = run(w, {Place {e, INF}}, ORDERED);
    EXPECT_FALSE(r[0].requested.has_value());
    EXPECT_EQ(lo(*r[0].offered), 14);
    EXPECT_EQ(r[0].call_delivered, 14);
    EXPECT_TRUE(r[0].call_edge_drops.empty());
    EXPECT_EQ(r[0].termination, Termination::Exhausted);
    EXPECT_EQ(refused(w, {Place {e, 0}}, ORDERED), "finite placement after infinite placement requires a new epoch");
  }
  {
    World w(par());
    auto& e = w.add({0}, {1}, ecmp());
    const auto r = run(w, {Replace {e, 12}}, ORDERED);
    EXPECT_EQ(r[0].call_delivered, 10);
    EXPECT_EQ(r[0].call_edge_drops, (Pairs {{0, 2}}));
    const auto before = w.state();
    (void)run(w, {Place {e, 0}}, ORDERED);
    expect_unchanged(w, before, "a zero place after a finite replace");
  }
  {
    World w(par());
    auto& e = w.add({0}, {1}, ecmp());
    (void)run(w, {Place {e, INF}}, ORDERED);
    const auto r = run(w, {Replace {e, 12}}, ORDERED);
    EXPECT_EQ(r[0].call_delivered, 10);
    EXPECT_EQ(r[0].call_edge_drops, (Pairs {{0, 2}}));
    EXPECT_NO_THROW((void)run(w, {Place {e, 0}}, ORDERED));
  }
}

// A returned counter crosses 2^64 (a place keeps its record's id).
TEST(ApplyActions, CountersCrossTwoToThe64) {
  const Amount w62 = kL62 - 1;
  World w(graph(2, {{0, 1, 10, 1}}));
  auto& e = w.add({0}, {1}, ecmp());
  for (int i = 0; i < 4; ++i) (void)e.place(w62);
  const auto id = e.records()[0].flow_id;
  const auto r = run(w, {Place {e, w62}}, ORDERED);
  // 5w = 5 * (2^62 - 1) = 2^64 + 2^62 - 5.
  const Counter128 want {1, (std::uint64_t {1} << 62) - 5};
  EXPECT_EQ(r[0].requested, want);
  EXPECT_EQ(r[0].offered, want);
  EXPECT_EQ(r[0].admitted, want);
  EXPECT_EQ(r[0].delivered, 10);
  ASSERT_EQ(r[0].edge_drops.size(), 1U);
  EXPECT_EQ(r[0].edge_drops[0].second, (Counter128 {1, (std::uint64_t {1} << 62) - 15}));
  EXPECT_EQ(r[0].call_delivered, 0);
  EXPECT_EQ(r[0].call_undelivered, w62);
  EXPECT_EQ(r[0].call_edge_drops, (Pairs {{0, w62}}));
  EXPECT_EQ(e.records()[0].flow_id, id);
}

// ---------------------------------------------------------------- SIMULTANEOUS

// D-simultaneous with its ORDERED contrast.
TEST(ApplyActions, OrderedAgainstSimultaneous) {
  for (BatchMode mode : {SIMULTANEOUS, ORDERED}) {
    World w(graph(2, {{0, 1, 10, 1}}));
    auto& a = w.add({0}, {1}, ecmp());
    auto& b = w.add({0}, {1}, ecmp());
    const auto r = run(w, {Replace {a, 10}, Replace {b, 10}}, mode);
    if (mode == SIMULTANEOUS) {
      EXPECT_EQ(delivered(r), (Amounts {5, 5}));
      EXPECT_EQ(r[0].call_edge_drops, (Pairs {{0, 5}}));
      EXPECT_EQ(r[1].call_edge_drops, (Pairs {{0, 5}}));
    } else {
      EXPECT_EQ(delivered(r), (Amounts {10, 0}));
      EXPECT_TRUE(r[0].call_edge_drops.empty());
      EXPECT_EQ(r[1].call_edge_drops, (Pairs {{0, 10}}));
    }
  }
}

// The participants' old load is released together; permutation, repeat, removal.
TEST(ApplyActions, SimultaneousViewReleasesTogether) {
  auto setup = [](World& w) {
    auto& o = w.add({0}, {1}, ecmp());
    auto& a = w.add({0}, {1}, ecmp());
    auto& b = w.add({0}, {1}, ecmp());
    (void)o.place(4);
    (void)a.place(1);
    (void)b.place(5);
    return std::tuple<FlowPolicy*, FlowPolicy*, FlowPolicy*> {&o, &a, &b};
  };
  World w(graph(2, {{0, 1, 10, 1}}));
  auto [o, a, b] = setup(w);
  const auto o_before = epoch(*o);
  const auto r = run(w, {Replace {*a, 10}, Replace {*b, 10}}, SIMULTANEOUS);
  EXPECT_EQ(delivered(r), (Amounts {3, 3}));
  EXPECT_EQ(r[0].call_edge_drops, (Pairs {{0, 7}}));
  EXPECT_EQ(r[0].delivery_delta, 2);
  EXPECT_EQ(r[1].delivery_delta, -2);
  EXPECT_EQ(w.ledger->carried(), (Amounts {10}));
  expect_unchanged(epoch(*o), o_before, "the outside owner");
  // Permuted, the same per-policy results in list order.
  World w2(graph(2, {{0, 1, 10, 1}}));
  auto [o2, a2, b2] = setup(w2);
  const auto r2 = run(w2, {Replace {*b2, 10}, Replace {*a2, 10}}, SIMULTANEOUS);
  expect_same(r2[0], r[1], "permuted B");
  expect_same(r2[1], r[0], "permuted A");
  // The same loads; the new records' FlowIds follow list order.
  EXPECT_EQ(w2.ledger->carried(), w.ledger->carried());
  EXPECT_EQ(a2->delivered(), a->delivered());
  EXPECT_EQ(b2->delivered(), b->delivered());
  EXPECT_EQ(b2->records()[0].flow_id.second + 1, a2->records()[0].flow_id.second);
  // Repeated: the same amounts, delivery_delta 0, new FlowIds.
  const auto ids = a->records()[0].flow_id;
  const auto again = run(w, {Replace {*a, 10}, Replace {*b, 10}}, SIMULTANEOUS);
  EXPECT_EQ(delivered(again), (Amounts {3, 3}));
  EXPECT_EQ(again[0].delivery_delta, 0);
  EXPECT_NE(a->records()[0].flow_id, ids);
  // Removing A later frees only its traffic.
  a->remove();
  EXPECT_EQ(w.ledger->carried(), (Amounts {7}));
  EXPECT_EQ(b->delivered(), 3);
  // ORDERED the same list (the contrast): A on 0 + 1, B on 0 + 5.
  World w3(graph(2, {{0, 1, 10, 1}}));
  auto [o3, a3, b3] = setup(w3);
  EXPECT_EQ(delivered(run(w3, {Replace {*a3, 10}, Replace {*b3, 10}}, ORDERED)), (Amounts {1, 5}));
}

// A ShortestPaths replace is prepared on the list's view (a residual snapshot reads it).
TEST(ApplyActions, SimultaneousPreparationReadsTheView) {
  for (BatchMode mode : {SIMULTANEOUS, ORDERED}) {
    World w(graph(2, {{0, 1, 10, 1}, {0, 1, 10, 1}}));
    ShortestPaths snap = shortest_paths(RouteAllocation::Weighted);
    snap.weight_basis = WeightBasis::ResidualSnapshot;
    auto& a = w.add({0}, {1}, FlowPolicyConfig {.routing = snap, .admission = LOSSY});
    (void)a.replace(2);
    auto& o = w.pin(tunnel(*w.g, 0, 1, {0}));
    (void)o.place(8);
    auto& b = w.add({0}, {1}, ecmp());
    const auto r = run(w, {Replace {a, 12}, Replace {b, 4}}, mode);
    if (mode == SIMULTANEOUS) {
      EXPECT_TRUE(r[0].weights_changed);
      EXPECT_FALSE(r[0].routes_changed);
      EXPECT_EQ(r[0].call_delivered, 10);
      EXPECT_EQ(r[0].call_edge_drops, (Pairs {{0, 1}, {1, 1}}));
      EXPECT_EQ(r[1].call_delivered, 2);
      EXPECT_EQ(r[1].call_edge_drops, (Pairs {{0, 1}, {1, 1}}));
      EXPECT_EQ(w.ledger->carried(), (Amounts {10, 10}));
    } else {
      EXPECT_EQ(r[0].call_delivered, 12);
      EXPECT_EQ(r[1].call_delivered, 0);
      EXPECT_EQ(r[1].call_edge_drops, (Pairs {{0, 2}, {1, 2}}));
    }
  }
  // The other participant's old load is released before the preparation: B (a tunnel on edge 0)
  // carries 8, the list replaces A (12) and B (8). On the list's view A's snapshot is (10, 10):
  // A splits (6, 6), edge 0 shares 10 over A 6 and B 8 by highest averages, (4, 6) (A created
  // first). On A's own view (B's 8 counted) it would be (2, 10), and nothing would drop. A single
  // path by HIGHER_AVAILABLE_CAPACITY ties on the list's view and takes edge 0 (on A's own view,
  // edge 1): edge 0 shares 10 over A 12 and B 8, (6, 4).
  ShortestPaths snap = shortest_paths(RouteAllocation::Weighted);
  snap.weight_basis = WeightBasis::ResidualSnapshot;
  ShortestPaths single = shortest_paths();
  single.tie_break = RouteTieBreak::HigherAvailableCapacity;
  single.multipath = false;
  single.parallel_edges = false;
  for (const auto& [model, a_flows, want_delivered, want_carried] :
       std::vector<std::tuple<ShortestPaths, Pairs, Amounts, Amounts>> {{snap, {{0, 4}, {1, 6}}, {10, 6}, {10, 6}},
                                                                        {single, {{0, 6}}, {6, 4}, {10, 0}}}) {
    World w(graph(2, {{0, 1, 10, 1}, {0, 1, 10, 1}}));
    auto& a = w.add({0}, {1}, FlowPolicyConfig {.routing = model, .admission = LOSSY});
    auto& b = w.pin(tunnel(*w.g, 0, 1, {0}));
    (void)b.replace(8);
    const auto r = run(w, {Replace {a, 12}, Replace {b, 8}}, SIMULTANEOUS);
    ASSERT_EQ(a.records().size(), 1U);
    EXPECT_EQ(a.records()[0].edge_flows, a_flows);
    EXPECT_EQ(delivered(r), want_delivered);
    EXPECT_EQ(w.ledger->carried(), want_carried);
  }
}

// Ties follow creation order, never list position.
TEST(ApplyActions, TiesFollowCreationOrder) {
  World w(graph(2, {{0, 1, 1, 1}}));
  const auto prog = tunnel(*w.g, 0, 1, {0});
  auto& p0 = w.pin(prog);
  auto& p1 = w.pin(prog);
  const auto r = run(w, {Replace {p1, 1}, Replace {p0, 1}}, SIMULTANEOUS);
  EXPECT_EQ(delivered(r), (Amounts {0, 1}));
  EXPECT_EQ(r[0].call_edge_drops, (Pairs {{0, 1}}));
  World w2(graph(2, {{0, 1, 1, 1}}));
  const auto prog2 = tunnel(*w2.g, 0, 1, {0});
  auto& q1 = w2.pin(prog2);
  auto& q0 = w2.pin(prog2);
  EXPECT_EQ(delivered(run(w2, {Replace {q1, 1}, Replace {q0, 1}}, SIMULTANEOUS)), (Amounts {1, 0}));
}

// A pooled demand shares a link; every field of both entries.
TEST(ApplyActions, PooledDemandSharesALink) {
  World w(graph(5, {{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 1, 1}, {2, 3, 1, 1}, {4, 2, 10, 1}}));
  auto& a = w.add({0, 1}, {3}, ecmp());
  auto& b = w.pin(tunnel(*w.g, 4, 3, {4, 2}));
  const auto r = run(w, {Replace {a, 2}, Replace {b, 2}}, SIMULTANEOUS, kAll);
  const auto& ra = r[0];
  EXPECT_EQ(ra.call_delivered, 1);
  EXPECT_EQ(ra.call_undelivered, 1);
  EXPECT_EQ(ra.delivery_delta, 1);
  EXPECT_EQ(ra.call_edge_drops, (Pairs {{2, 1}}));
  EXPECT_EQ(ra.call_source_drops, 0);
  EXPECT_EQ(ra.termination, Termination::Complete);
  EXPECT_TRUE(ra.routes_changed);
  EXPECT_TRUE(ra.weights_changed);
  EXPECT_FALSE(ra.reservations_changed);
  ASSERT_NE(ra.program, nullptr);
  EXPECT_TRUE(ra.program->tunnels().empty());
  EXPECT_FALSE(ra.plan.has_value());
  EXPECT_EQ(lo(*ra.requested), 2);
  EXPECT_EQ(lo(*ra.offered), 2);
  EXPECT_EQ(lo(ra.admitted), 2);
  EXPECT_EQ(ra.delivered, 1);
  EXPECT_EQ(lo(*ra.unadmitted), 0);
  EXPECT_EQ(lo(ra.source_drops), 0);
  ASSERT_EQ(ra.edge_drops.size(), 1U);
  EXPECT_EQ(ra.edge_drops[0].first, 2);
  EXPECT_EQ(lo(ra.edge_drops[0].second), 1);
  ASSERT_EQ(ra.by_source.size(), 2U);
  EXPECT_EQ(ra.by_source[0].source, 0);
  EXPECT_EQ(lo(*ra.by_source[0].offered), 1);
  EXPECT_EQ(ra.by_source[0].delivered, 0);
  EXPECT_EQ(lo(ra.by_source[0].total_edge_drops), 1);
  EXPECT_EQ(ra.by_source[1].source, 1);
  EXPECT_EQ(lo(*ra.by_source[1].offered), 1);
  EXPECT_EQ(ra.by_source[1].delivered, 1);
  EXPECT_TRUE(ra.by_tunnel.empty());
  EXPECT_FALSE(ra.total_cost.has_value());
  EXPECT_FALSE(ra.marginal_cost_curve.has_value());
  EXPECT_EQ(*ra.edge_flows, (Amounts {1, 1, 0, 1, 0}));
  EXPECT_EQ(*ra.residual, (Amounts {9, 9, 0, 0, 8}));
  EXPECT_EQ(*ra.saturated_boundary_edges, (std::vector<EdgeId> {2, 3}));
  EXPECT_EQ(*ra.edge_reservations, (Amounts {0, 0, 0, 0, 0}));
  EXPECT_FALSE(ra.min_cut.has_value());
  const auto recs = a.records();
  ASSERT_EQ(recs.size(), 2U);
  EXPECT_EQ(recs[0].edge_flows, (Pairs {{0, 1}}));
  ASSERT_EQ(recs[0].edge_drops.size(), 1U);
  EXPECT_EQ(recs[1].edge_flows, (Pairs {{1, 1}, {3, 1}}));
  EXPECT_TRUE(recs[1].edge_drops.empty());
  const auto& rb = r[1];
  EXPECT_EQ(rb.call_delivered, 1);
  EXPECT_EQ(rb.call_undelivered, 1);
  EXPECT_EQ(rb.delivery_delta, 1);
  EXPECT_EQ(rb.call_edge_drops, (Pairs {{2, 1}}));
  EXPECT_EQ(rb.termination, Termination::Complete);
  EXPECT_TRUE(rb.routes_changed);
  EXPECT_TRUE(rb.weights_changed);
  ASSERT_EQ(rb.by_tunnel.size(), 1U);
  EXPECT_EQ(rb.by_tunnel[0].head, 4);
  EXPECT_EQ(rb.by_tunnel[0].status, TunnelStatus::Up);
  EXPECT_EQ(lo(rb.by_tunnel[0].assigned), 2);
  EXPECT_EQ(rb.by_tunnel[0].delivered, 1);
  EXPECT_EQ(lo(rb.by_tunnel[0].total_edge_drops), 1);
  EXPECT_EQ(*rb.edge_flows, (Amounts {0, 0, 1, 0, 2}));
  EXPECT_EQ(*rb.residual, (Amounts {9, 9, 0, 0, 8}));
  EXPECT_EQ(w.ledger->carried(), (Amounts {1, 1, 1, 1, 2}));
}

// The list itself: n entries in input order.
TEST(ApplyActions, ResultsInInputOrder) {
  World w(graph(4, {{0, 1, 10, 1}, {2, 3, 10, 1}}));
  auto& a = w.pin(tunnel(*w.g, 0, 1, {0}));
  auto& b = w.pin(tunnel(*w.g, 2, 3, {1}));
  for (BatchMode mode : {ORDERED, SIMULTANEOUS}) {
    const auto r = run(w, {Replace {b, 7}, Replace {a, 3}}, mode);
    ASSERT_EQ(r.size(), 2U);
    EXPECT_EQ(delivered(r), (Amounts {7, 3}));
    EXPECT_EQ(r[0].program->tunnels()[0].head, 2);
  }
}

// The demand level above tunnels and sources.
TEST(ApplyActions, DemandLevelAboveTunnels) {
  World w(graph(2, {{0, 1, 1, 1}}));
  auto& a = w.pin(tunnels(*w.g, {0}, {1}, {{0, {0}, 1}, {0, {0}, 1}}));
  auto& b = w.pin(tunnel(*w.g, 0, 1, {0}));
  const auto r = run(w, {Replace {a, 2}, Replace {b, 2}}, SIMULTANEOUS);
  EXPECT_EQ(r[0].call_delivered, 1);
  EXPECT_EQ(r[0].by_tunnel[0].delivered, 1);
  EXPECT_EQ(r[0].by_tunnel[1].delivered, 0);
  EXPECT_EQ(r[1].call_delivered, 0);
  World w2(graph(2, {{0, 1, 1, 1}}));
  auto& b2 = w2.pin(tunnel(*w2.g, 0, 1, {0}));
  auto& a2 = w2.pin(tunnels(*w2.g, {0}, {1}, {{0, {0}, 1}, {0, {0}, 1}}));
  EXPECT_EQ(delivered(run(w2, {Replace {a2, 2}, Replace {b2, 2}}, SIMULTANEOUS)), (Amounts {0, 1}));
}

TEST(ApplyActions, DemandLevelAboveSources) {
  World w(graph(6, {{1, 3, 100, 1}, {2, 3, 100, 1}, {3, 4, 1, 1}, {5, 3, 100, 1}}));
  auto& a = w.pin(tunnels(*w.g, {1, 2}, {4}, {{1, {0, 2}, 1}, {2, {1, 2}, 1}}));
  auto& b = w.pin(tunnel(*w.g, 5, 4, {3, 2}));
  const auto r = run(w, {Replace {a, 2}, Replace {b, 2}}, SIMULTANEOUS);
  EXPECT_EQ(r[0].call_delivered, 1);
  EXPECT_EQ(r[0].by_tunnel[0].delivered, 1);
  EXPECT_EQ(r[1].call_delivered, 0);
}

// LIST-black-hole-own and a Forward with black holes: each participant's black holes are its own.
TEST(ApplyActions, BlackHolesOfOneParticipant) {
  {
    World w(graph(2, {{0, 1, 10, 1}}));
    auto& a = w.pin(tunnel(*w.g, 0, 1, {0}), LOSSY, no_edges(*w.g, {0}));
    auto& b = w.add({0}, {1}, ecmp());
    const auto r = run(w, {Replace {a, 10}, Replace {b, 10}}, SIMULTANEOUS, OutputOptions {.residual = true, .saturated_boundary_edges = true});
    EXPECT_EQ(r[0].call_delivered, 0);
    EXPECT_EQ(r[0].call_edge_drops, (Pairs {{0, 10}}));
    EXPECT_EQ(r[0].by_tunnel[0].status, TunnelStatus::Blackholed);
    EXPECT_EQ(r[1].call_delivered, 10);
    EXPECT_TRUE(r[1].call_edge_drops.empty());
    EXPECT_EQ(*r[0].residual, (Amounts {0}));
    EXPECT_EQ(*r[1].residual, (Amounts {0}));
    EXPECT_EQ(*r[0].saturated_boundary_edges, (std::vector<EdgeId> {0}));
  }
  {
    World w(graph(2, {{0, 1, 10, 1}, {0, 1, 10, 1}}));
    auto& a = w.add({0}, {1}, ecmp());
    auto& b = w.add({0}, {1}, ecmp());
    (void)a.replace(2);
    const auto bh = no_edges(*w.g, {1});
    const auto r = run(w, {Forward {a, 10, bh.get()}, Replace {b, 10}}, SIMULTANEOUS);
    EXPECT_EQ(r[0].call_delivered, 5);
    EXPECT_EQ(r[0].call_edge_drops, (Pairs {{1, 5}}));
    EXPECT_EQ(r[1].call_delivered, 10);
    EXPECT_EQ(w.ledger->carried(), (Amounts {10, 5}));
  }
}

// LIST-black-hole-upstream: a black hole downstream of a shared link; the two links are no
// transparent pair.
TEST(ApplyActions, BlackHoleDownstreamOfASharedLink) {
  World w(graph(3, {{0, 1, 10, 1}, {1, 2, 10, 1}}));
  auto& a = w.pin(tunnel(*w.g, 0, 2, {0, 1}), LOSSY, no_edges(*w.g, {1}));
  auto& b = w.pin(tunnel(*w.g, 0, 2, {0, 1}));
  const auto r = run(w, {Replace {a, 10}, Replace {b, 10}}, SIMULTANEOUS, OutputOptions {.residual = true, .saturated_boundary_edges = true});
  EXPECT_EQ(r[0].call_edge_drops, (Pairs {{0, 5}, {1, 5}}));
  EXPECT_EQ(r[0].by_tunnel[0].status, TunnelStatus::Blackholed);
  EXPECT_EQ(r[1].call_delivered, 5);
  EXPECT_EQ(r[1].call_edge_drops, (Pairs {{0, 5}}));
  EXPECT_EQ(*r[0].residual, (Amounts {0, 5}));
  EXPECT_EQ(*r[0].saturated_boundary_edges, (std::vector<EdgeId> {0}));
  EXPECT_EQ(*r[1].saturated_boundary_edges, (std::vector<EdgeId> {0}));
}

// A T1 Forward keeps its holds and plan.
TEST(ApplyActions, ForwardKeepsHolds) {
  World w(graph(2, {{0, 1, 10, 1}}));
  auto& t = w.add({0}, {1}, t1());
  auto& e = w.add({0}, {1}, ecmp());
  (void)t.replace(5, std::optional<Amount> {5});
  const auto r = run(w, {Forward {t, 10}, Replace {e, 10}}, SIMULTANEOUS);
  EXPECT_EQ(r[0].call_delivered, 5);
  EXPECT_EQ(r[0].call_edge_drops, (Pairs {{0, 5}}));
  EXPECT_FALSE(r[0].routes_changed);
  EXPECT_FALSE(r[0].reservations_changed);
  ASSERT_TRUE(r[0].plan.has_value());
  EXPECT_EQ(r[0].plan->total_reserved_bandwidth, 5);
  ASSERT_EQ(t.records().size(), 1U);
  EXPECT_EQ(t.records()[0].edge_reservations, (Pairs {{0, 5}}));
  EXPECT_EQ(r[1].call_delivered, 5);
  EXPECT_EQ(w.ledger->carried(), (Amounts {10}));
  EXPECT_EQ(w.ledger->reserved(), (Amounts {5}));
}

// Empty rosters and masked participants.
TEST(ApplyActions, EmptyRosterParticipants) {
  World w(graph(2, {{0, 1, 10, 1}}));
  const auto mask = no_edges(*w.g, {0});
  auto& m = w.add({0}, {1}, ecmp(), mask);
  auto& e = w.add({0}, {1}, ecmp());
  auto& p = w.pin(tunnel(*w.g, 0, 1, {0}), LOSSY, nullptr, mask);
  const auto r = run(w, {Replace {m, 5}, Replace {e, 10}, Replace {p, 5}}, SIMULTANEOUS);
  EXPECT_EQ(r[0].termination, Termination::NoRoute);
  EXPECT_EQ(r[0].call_source_drops, 5);
  EXPECT_EQ(lo(*r[0].offered), 5);
  EXPECT_EQ(lo(r[0].source_drops), 5);
  EXPECT_TRUE(m.records().empty());
  EXPECT_EQ(r[1].call_delivered, 10);
  EXPECT_EQ(r[2].termination, Termination::NoRoute);
  EXPECT_EQ(r[2].call_source_drops, 5);
  ASSERT_EQ(p.records().size(), 1U);
  EXPECT_EQ(p.records()[0].delivered, 0);
}

TEST(ApplyActions, EverythingMasked) {
  for (BatchMode mode : {ORDERED, SIMULTANEOUS}) {
    World w(graph(2, {{0, 1, 10, 1}}));
    const auto mask = no_edges(*w.g, {0});
    auto& m1 = w.add({0}, {1}, ecmp(), mask);
    auto& m2 = w.add({0}, {1}, ecmp(), mask);
    auto& p = w.pin(tunnel(*w.g, 0, 1, {0}), LOSSY, nullptr, mask);
    const auto r = run(w, {Replace {m1, 5}, Replace {m2, 5}, Replace {p, 5}}, mode);
    for (const auto& x : r) {
      EXPECT_EQ(x.termination, Termination::NoRoute);
      EXPECT_EQ(x.call_source_drops, 5);
    }
    EXPECT_TRUE(m1.records().empty());
    EXPECT_TRUE(m2.records().empty());
    ASSERT_EQ(p.records().size(), 1U);
    EXPECT_EQ(w.ledger->carried(), (Amounts {0}));
  }
}

// Zero offers prepare and record; an ORDERED zero place does not.
TEST(ApplyActions, ZeroOffers) {
  World w(graph(2, {{0, 1, 10, 1}}));
  auto& a = w.add({0}, {1}, ecmp());
  auto& b = w.add({0}, {1}, ecmp());
  const auto r = run(w, {Replace {a, 0}, Replace {b, 0}}, SIMULTANEOUS);
  for (const auto& x : r) {
    EXPECT_EQ(x.termination, Termination::Complete);
    EXPECT_EQ(x.call_delivered, 0);
    EXPECT_NE(x.program, nullptr);
    EXPECT_TRUE(x.routes_changed);
    EXPECT_TRUE(x.weights_changed);
  }
  EXPECT_EQ(a.records().size(), 1U);
  EXPECT_EQ(b.records().size(), 1U);
  EXPECT_EQ(w.ledger->carried(), (Amounts {0}));
}

// Holds never enter the view; orphans are fixed load.
TEST(ApplyActions, HoldsAndOrphansStayOutside) {
  {
    World w(graph(2, {{0, 1, 10, 1}}));
    auto& t = w.add({0}, {1}, t1());
    auto& e = w.add({0}, {1}, ecmp());
    (void)t.replace(0, std::optional<Amount> {10});
    EXPECT_EQ(run(w, {Replace {e, 10}}, SIMULTANEOUS)[0].call_delivered, 10);
  }
  {
    World w(graph(2, {{0, 1, 10, 1}}));
    w.add({0}, {1}, ecmp()).place(4);
    w.p.pop_back();  // orphaned
    auto& a = w.add({0}, {1}, ecmp());
    auto& b = w.add({0}, {1}, ecmp());
    EXPECT_EQ(delivered(run(w, {Replace {a, 10}, Replace {b, 10}}, SIMULTANEOUS)), (Amounts {3, 3}));
    EXPECT_EQ(w.ledger->carried(), (Amounts {10}));
  }
}

// An iteration error, and a failed check, publish nothing.
TEST(ApplyActions, IterationErrorsChangeNothing) {
  // D-C3 as three Pinned tunnels: A 0, 1 (1); B 1, 2 (1); C 2, 0 (2); V 1, 1, 2 (CYC-like layout:
  // each link its own nodes, wide connectors).
  World w(shared(make_graph(8, {{0, 1, 1, 1}, {2, 3, 1, 1}, {4, 5, 2, 1}, {1, 2, 1000, 1}, {3, 4, 1000, 1}, {5, 0, 1000, 1},
                                {6, 0, 1000, 1}, {6, 2, 1000, 1}, {6, 4, 1000, 1}, {3, 7, 1000, 1}, {5, 7, 1000, 1}, {1, 7, 1000, 1}})));
  // EdgeIds sort by tail: map the input rows.
  const auto id = id_of_input(*w.g);
  auto ids = [&](std::vector<int> in) {
    std::vector<EdgeId> out;
    for (int i : in) out.push_back(id[static_cast<std::size_t>(i)]);
    return out;
  };
  auto& a = w.pin(tunnel(*w.g, 6, 7, ids({6, 0, 3, 1, 9})));
  auto& b = w.pin(tunnel(*w.g, 6, 7, ids({7, 1, 4, 2, 10})));
  auto& c = w.pin(tunnel(*w.g, 6, 7, ids({8, 2, 5, 0, 11})));
  (void)a.replace(1);
  (void)b.replace(1);
  (void)c.replace(1);
  EXPECT_EQ(refused(w, {Replace {a, 1}, Replace {b, 1}, Replace {c, 2}}, SIMULTANEOUS), kPeriodTwo);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 3);
  // The independent check across demands, through the hook (a settling cyclic list).
  World s(graph(6, {{0, 1, 1, 1}, {1, 2, 1000, 1}, {1, 5, 1000, 1}, {2, 3, 1, 1}, {3, 0, 1000, 1}, {3, 5, 1000, 1},
                    {4, 0, 1000, 1}, {4, 2, 1000, 1}}));
  auto& x = s.pin(tunnel(*s.g, 4, 5, {6, 0, 1, 3, 5}));
  auto& y = s.pin(tunnel(*s.g, 4, 5, {7, 3, 4, 0, 2}));
  const auto before = s.state();
  detail::thread_lossy_counters().corrupt = true;
  std::string what;
  try {
    (void)apply_actions(*s.ledger, Actions {Replace {x, 1}, Replace {y, 1}}, SIMULTANEOUS);
  } catch (const std::runtime_error& e) {
    what = e.what();
  }
  detail::thread_lossy_counters().corrupt = false;
  EXPECT_EQ(what, "lossy forwarding fixed-point verification failed");
  expect_unchanged(s, before, "the failed check");
  // One unit moved between the two occurrences of one barrier before the final forwarding: every
  // continuation stays consistent, only that link's loss rule is wrong.
  const auto before_shift = s.state();
  detail::thread_lossy_counters().shift = true;
  what.clear();
  try {
    (void)apply_actions(*s.ledger, Actions {Replace {x, 1}, Replace {y, 1}}, SIMULTANEOUS);
  } catch (const std::runtime_error& e) {
    what = e.what();
  }
  detail::thread_lossy_counters().shift = false;
  EXPECT_EQ(what, "lossy forwarding fixed-point verification failed");
  expect_unchanged(s, before_shift, "the failed check of a moved unit");
  EXPECT_EQ(delivered(run(s, {Replace {x, 1}, Replace {y, 1}}, SIMULTANEOUS)), (Amounts {1, 0}));
#endif
}

// Link sums beyond 2^63 and 2^64; disjoint maximal participants.
TEST(ApplyActions, MaximalOffersShareALink) {
  const Amount w62 = kL62 - 1;
  for (const auto& [n, cap, want] : std::vector<std::tuple<int, Amount, Amounts>> {
           {2, 10, {5, 5}}, {3, 10, {4, 3, 3}}, {5, 10, {2, 2, 2, 2, 2}}, {5, 11, {3, 2, 2, 2, 2}}}) {
    World w(graph(2, {{0, 1, cap, 1}}));
    const auto prog = tunnel(*w.g, 0, 1, {0});
    Actions list;
    for (int i = 0; i < n; ++i) list.push_back(Replace {w.pin(prog), w62});
    const auto r = run(w, list, SIMULTANEOUS);
    EXPECT_EQ(delivered(r), want) << n << " demands";
    for (int i = 0; i < n; ++i) EXPECT_EQ(r[static_cast<std::size_t>(i)].call_edge_drops, (Pairs {{0, w62 - want[static_cast<std::size_t>(i)]}}));
  }
}

TEST(ApplyActions, DisjointMaximalParticipants) {
  const Amount w62 = (kL62 - 1) / 5;  // the graph's capacity sum stays below 2^62
  std::vector<InEdge> rows;
  for (int i = 0; i < 5; ++i) rows.push_back({0, 1, w62, 1});
  World w(graph(2, rows));
  Actions list;
  for (EdgeId e = 0; e < 5; ++e) list.push_back(Replace {w.pin(tunnel(*w.g, 0, 1, {e})), w62});
  const auto r = run(w, list, SIMULTANEOUS);
  EXPECT_EQ(delivered(r), Amounts(5, w62));
  EXPECT_EQ(w.ledger->carried(), Amounts(5, w62));
}

// An imported T1 program beside its live owner, unequal baselines; the orphan variant.
TEST(ApplyActions, ImportedTunnelsBookNoHolds) {
  {
    World w(graph(2, {{0, 1, 10, 1}}), {8}, {6});
    auto& t = w.add({0}, {1}, t1());
    (void)t.replace(0, std::optional<Amount> {6});
    auto& p = w.pin(t.program());
    const auto r = run(w, {Forward {t, 8}, Replace {p, 8}}, SIMULTANEOUS, OutputOptions {.edge_reservations = true});
    EXPECT_EQ(delivered(r), (Amounts {4, 4}));
    EXPECT_EQ(w.ledger->carried(), (Amounts {8}));
    EXPECT_EQ(w.ledger->reserved(), (Amounts {6}));
    EXPECT_EQ(*r[1].edge_reservations, (Amounts {0}));
    EXPECT_FALSE(r[1].plan.has_value());
    ASSERT_TRUE(r[0].plan.has_value());
  }
  {
    World w(graph(2, {{0, 1, 10, 1}}), {8}, {6});
    w.add({0}, {1}, t1()).replace(0, std::optional<Amount> {6});
    const auto prog = w.p.back()->program();
    w.p.pop_back();
    auto& p = w.pin(prog);
    EXPECT_EQ(run(w, {Replace {p, 8}}, SIMULTANEOUS)[0].call_delivered, 8);
    EXPECT_EQ(w.ledger->reserved(), (Amounts {6}));
  }
}

// A zero Forward keeps holds, plan and program, retires its old id.
TEST(ApplyActions, ZeroForward) {
  World w(graph(2, {{0, 1, 10, 1}}));
  auto& t = w.add({0}, {1}, t1());
  auto& e = w.add({0}, {1}, ecmp());
  (void)t.replace(6, std::optional<Amount> {6});
  const auto old = t.records()[0].flow_id;
  const auto prog = t.program();
  const auto r = run(w, {Forward {t, 0}, Replace {e, 10}}, SIMULTANEOUS);
  EXPECT_EQ(r[0].termination, Termination::Complete);
  EXPECT_EQ(r[0].call_delivered, 0);
  EXPECT_EQ(r[0].delivery_delta, -6);
  EXPECT_EQ(lo(r[0].by_tunnel[0].assigned), 0);
  EXPECT_EQ(t.program(), prog);
  ASSERT_EQ(t.records().size(), 1U);
  EXPECT_NE(t.records()[0].flow_id, old);
  EXPECT_EQ(t.records()[0].edge_reservations, (Pairs {{0, 6}}));
  EXPECT_EQ(r[1].call_delivered, 10);
  const auto f = t.forward(8);
  EXPECT_EQ(f.call_delivered, 0);
  EXPECT_EQ(f.call_edge_drops, (Pairs {{0, 8}}));
  EXPECT_EQ(w.ledger->reserved(), (Amounts {6}));
}

// Black-holed sources and nodes; a repair mask first.
TEST(ApplyActions, BlackHoledSourcesAndNodes) {
  World w(graph(5, {{0, 2, 10, 1}, {1, 2, 10, 1}, {3, 4, 10, 1}}));
  auto& a = w.add({0, 1}, {2}, ecmp());
  (void)a.replace(2);
  auto& z = w.pin(tunnel(*w.g, 3, 4, {2}));
  const auto bh = no_nodes(*w.g, {0});
  const auto r = run(w, {Forward {a, 4, bh.get()}, Replace {z, 3}}, SIMULTANEOUS);
  EXPECT_EQ(r[0].call_source_drops, 2);
  EXPECT_TRUE(r[0].call_edge_drops.empty());
  EXPECT_EQ(r[0].call_delivered, 2);
  EXPECT_EQ(r[1].call_delivered, 3);
  // A Pinned over A's program under a repair mask excluding node 0: source 1 takes the whole offer.
  auto& q = w.pin(a.program(), LOSSY, nullptr, no_nodes(*w.g, {0}));
  const auto rq = run(w, {Replace {q, 4}, Replace {z, 3}}, SIMULTANEOUS);
  EXPECT_EQ(rq[0].call_delivered, 4);
  EXPECT_EQ(rq[0].call_source_drops, 0);
  // An internal black-holed node.
  World w2(graph(3, {{0, 1, 10, 1}, {1, 2, 10, 1}}));
  auto& n = w2.pin(tunnel(*w2.g, 0, 2, {0, 1}), LOSSY, no_nodes(*w2.g, {1}));
  const auto rn = run(w2, {Replace {n, 6}}, SIMULTANEOUS);
  EXPECT_EQ(rn[0].by_tunnel[0].status, TunnelStatus::Blackholed);
  EXPECT_EQ(rn[0].call_edge_drops, (Pairs {{0, 6}}));
}

// Zero service against zero weight, beside a positive participant.
TEST(ApplyActions, ZeroServiceAndZeroWeight) {
  World w(graph(4, {{0, 1, 0, 1}, {2, 3, 10, 1}}));
  auto& z = w.add({0}, {1}, ecmp());
  auto& pos = w.pin(tunnel(*w.g, 2, 3, {1}));
  const auto r = run(w, {Replace {z, 2}, Replace {pos, 4}}, SIMULTANEOUS);
  EXPECT_EQ(r[0].termination, Termination::Complete);
  EXPECT_EQ(r[0].call_edge_drops, (Pairs {{0, 2}}));
  EXPECT_EQ(r[1].call_delivered, 4);
  World w2(graph(4, {{0, 1, 10, 1}, {2, 3, 10, 1}}));
  auto& zw = w2.pin(tunnel(*w2.g, 0, 1, {0}, 0));
  auto& pos2 = w2.pin(tunnel(*w2.g, 2, 3, {1}));
  const auto r2 = run(w2, {Replace {zw, 2}, Replace {pos2, 4}}, SIMULTANEOUS);
  EXPECT_EQ(r2[0].termination, Termination::NoActiveWeight);
  EXPECT_EQ(r2[0].call_source_drops, 2);
  EXPECT_EQ(zw.records().size(), 1U);
  EXPECT_EQ(r2[1].call_delivered, 4);
}

// Every root of a pooled demand is an injection (alone, and beside a disjoint participant),
// equal to the scalar replace.
TEST(ApplyActions, PooledRootsBlockSerialPairs) {
  auto make = [] { return graph(6, {{0, 1, 2, 1}, {1, 2, 2, 1}, {2, 3, 2, 1}, {4, 5, 10, 1}}); };
  World scalar(make());
  const auto want = scalar.add({0, 1, 2}, {3}, ecmp()).replace(6);
  for (bool beside : {false, true}) {
    World w(make());
    auto& a = w.add({0, 1, 2}, {3}, ecmp());
    auto& z = w.pin(tunnel(*w.g, 4, 5, {3}));
    const auto r = run(w, beside ? Actions {Replace {a, 6}, Replace {z, 1}} : Actions {Replace {a, 6}}, SIMULTANEOUS);
    expect_same(r[0], want, beside ? "beside" : "alone");
    EXPECT_EQ(r[0].call_delivered, 2);
    EXPECT_EQ(r[0].by_source[0].delivered, 1);
    EXPECT_EQ(r[0].by_source[1].delivered, 0);
    EXPECT_EQ(r[0].by_source[2].delivered, 1);
  }
}

// A pooled cycle keeps per-source state (solved at round 3).
TEST(ApplyActions, PooledCycleKeepsSourceState) {
  World w(shared(make_graph(9, {{0, 1, 1, 1}, {2, 3, 1, 1}, {4, 0, 100, 1}, {1, 2, 100, 1}, {5, 2, 100, 1}, {3, 6, 100, 1},
                                {7, 2, 100, 1}, {3, 0, 100, 1}, {1, 8, 100, 1}})));
  const auto id = id_of_input(*w.g);
  auto& a = w.add({4, 5}, {6}, ecmp());
  auto& b = w.pin(tunnel(*w.g, 7, 8, {id[6], id[1], id[7], id[0], id[8]}));
  const auto r = run(w, {Replace {a, 2}, Replace {b, 1}}, SIMULTANEOUS);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().groups, 1);
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 3);
#endif
  EXPECT_EQ(r[0].call_delivered, 1);
  EXPECT_EQ(r[0].by_source[0].delivered, 1);
  EXPECT_EQ(r[0].by_source[1].delivered, 0);
  EXPECT_EQ(r[1].call_delivered, 0);
  EXPECT_EQ(r[1].call_edge_drops, (Pairs {{id[1], 1}}));
}

// Dark arrivals add no barrier: the round-limit example's graph plus 16 -> 6, its three demands
// (heads 12, 13, 14, offers 13, 11, 10; black holes {5}, used by head 12 only) and a fourth demand D
// from 16 to 1 (offer 1) whose black holes {16, 11} drop its offer on its first link. The cyclic group
// keeps its four physical links (0, 2, 7, 11): the round-limit error after 256 rounds, nothing
// published (a dark copy of D's link 11 as a fifth barrier would allow 320).
TEST(ApplyActions, DarkArrivalsAddNoBarrier) {
  World w(shared(make_graph(17, {{0, 1, 12, 1},    {1, 6, 1000, 1},  {2, 3, 9, 1},     {3, 0, 1000, 1},  {3, 10, 1000, 1},
                                 {4, 5, 7, 1},     {5, 6, 1000, 1},  {6, 7, 13, 1},    {7, 2, 1000, 1},  {7, 10, 1000, 1},
                                 {8, 9, 1, 1},     {10, 11, 10, 1},  {11, 2, 1000, 1}, {12, 4, 1000, 1}, {13, 6, 1000, 1},
                                 {14, 0, 1000, 1}, {16, 6, 1000, 1}})));
  const auto id = id_of_input(*w.g);
  auto ids = [&](std::vector<int> in) {
    std::vector<EdgeId> out;
    for (int i : in) out.push_back(id[static_cast<std::size_t>(i)]);
    return out;
  };
  auto& a = w.pin(tunnel(*w.g, 12, 11, ids({13, 5, 6, 7, 8, 2, 4, 11})), LOSSY, no_edges(*w.g, ids({5})));
  auto& b = w.pin(tunnel(*w.g, 13, 1, ids({14, 7, 9, 11, 12, 2, 3, 0})));
  auto& c = w.pin(tunnel(*w.g, 14, 7, ids({15, 0, 1, 7})));
  const Actions three {Replace {a, 13}, Replace {b, 11}, Replace {c, 10}};
  EXPECT_EQ(refused(w, three, SIMULTANEOUS), kRoundLimit);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 256);
#endif
  auto& d = w.pin(tunnel(*w.g, 16, 1, ids({16, 7, 9, 11, 12, 2, 3, 0})), LOSSY, no_edges(*w.g, ids({16, 11})));
  Actions four = three;
  four.push_back(Replace {d, 1});
  EXPECT_EQ(refused(w, four, SIMULTANEOUS), kRoundLimit);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().groups, 1);
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 256);
#endif
  // D alone: its offer dies on its first link.
  const auto r = run(w, {Replace {d, 1}}, SIMULTANEOUS);
  EXPECT_EQ(r[0].call_edge_drops, (Pairs {{id[16], 1}}));
  EXPECT_EQ(r[0].by_tunnel[0].status, TunnelStatus::Blackholed);
  // A dark arrival on a link that otherwise fits: E from 16 over a wide feeder to node 7, black-holed
  // on the cycle's connector 7 -> 2 (1000 units, row 8), offers 1000. The connector carries only head
  // 12's arrivals (at most 13): it stays no barrier, and the group fails after 256 rounds (counting
  // E's 1000 in its no-drop load would make it a fifth barrier: 320).
  World w2(shared(make_graph(17, {{0, 1, 12, 1},    {1, 6, 1000, 1},  {2, 3, 9, 1},     {3, 0, 1000, 1},  {3, 10, 1000, 1},
                                  {4, 5, 7, 1},     {5, 6, 1000, 1},  {6, 7, 13, 1},    {7, 2, 1000, 1},  {7, 10, 1000, 1},
                                  {8, 9, 1, 1},     {10, 11, 10, 1},  {11, 2, 1000, 1}, {12, 4, 1000, 1}, {13, 6, 1000, 1},
                                  {14, 0, 1000, 1}, {16, 7, 2000, 1}})));
  const auto id2 = id_of_input(*w2.g);
  auto ids2 = [&](std::vector<int> in) {
    std::vector<EdgeId> out;
    for (int i : in) out.push_back(id2[static_cast<std::size_t>(i)]);
    return out;
  };
  auto& a2 = w2.pin(tunnel(*w2.g, 12, 11, ids2({13, 5, 6, 7, 8, 2, 4, 11})), LOSSY, no_edges(*w2.g, ids2({5})));
  auto& b2 = w2.pin(tunnel(*w2.g, 13, 1, ids2({14, 7, 9, 11, 12, 2, 3, 0})));
  auto& c2 = w2.pin(tunnel(*w2.g, 14, 7, ids2({15, 0, 1, 7})));
  auto& e2 = w2.pin(tunnel(*w2.g, 16, 2, ids2({16, 8})), LOSSY, no_edges(*w2.g, ids2({8})));
  EXPECT_EQ(refused(w2, {Replace {a2, 13}, Replace {b2, 11}, Replace {c2, 10}, Replace {e2, 1000}}, SIMULTANEOUS), kRoundLimit);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 256);
#endif
}

// An entry whose share is 0 can carry nothing, so it makes no link wait and joins no group. Two
// step-3 cycles (CYC(1, 2) twice, 4 -> 5 and 16 -> 17, offers 1 and 2 each) bridged by two Pinned
// tunnels of offer 0 (4 -> 17 over the first cycle's X and the second's links, 16 -> 5 the other
// way): the period-2 error at round 3 with or without the bridges (with arcs through them the two
// groups merged and failed at the round limit). And a link X before a period-2 pair: A over X, L0,
// L1, B over L1, L0, offers 2; Z over L0 then X offers 0, as a Pinned tunnel or a ShortestPaths
// demand: the period-2 error at round 3 with or without Z.
TEST(ApplyActions, ZeroSharesAddNoArc) {
  auto ids_of = [](const std::vector<EdgeId>& id, std::vector<int> in) {
    std::vector<EdgeId> out;
    for (int i : in) out.push_back(id[static_cast<std::size_t>(i)]);
    return out;
  };
  const std::vector<InEdge> bridge_rows {{0, 1, 1, 1},       {1, 2, 1000, 1},   {1, 5, 1000, 1},   {2, 3, 2, 1},       {3, 0, 1000, 1},
                                         {3, 5, 1000, 1},    {4, 0, 1000, 1},   {4, 2, 1000, 1},   {10, 11, 1, 1},     {11, 14, 1000, 1},
                                         {12, 13, 2, 1},     {13, 10, 1000, 1}, {13, 17, 1000, 1}, {14, 15, 1, 1},     {15, 12, 1000, 1},
                                         {15, 17, 1000, 1},  {16, 10, 1000, 1}, {16, 12, 1000, 1}, {1, 10, 1000, 1},   {11, 0, 1000, 1}};
  const std::vector<std::tuple<NodeId, NodeId, std::vector<int>, Amount>> bridged {
      {4, 5, {6, 0, 1, 3, 5}, 1},           {4, 5, {7, 3, 4, 0, 2}, 2},       {16, 17, {16, 8, 9, 13, 14, 10, 12}, 1},
      {16, 17, {17, 10, 11, 8, 9, 13, 15}, 2}, {4, 17, {6, 0, 18, 8, 9, 13, 15}, 0}, {16, 5, {16, 8, 19, 0, 2}, 0}};
  for (std::size_t count : {std::size_t {4}, std::size_t {6}}) {
    World w(shared(make_graph(18, bridge_rows)));
    const auto id = id_of_input(*w.g);
    Actions list;
    for (std::size_t i = 0; i < count; ++i) {
      const auto& [s, d, path, offer] = bridged[i];
      auto& p = w.pin(tunnel(*w.g, s, d, ids_of(id, path)));
      (void)p.replace(1);
      list.push_back(Replace {p, offer});
    }
    EXPECT_EQ(refused(w, list, SIMULTANEOUS), kPeriodTwo) << count << " demands";
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(detail::thread_lossy_counters().rounds, 3) << count << " demands";
#endif
  }
  // sA 0, sB 1, sZ 2, n1..n6 3..8, dA 9, dB 10, dZ 11; rows X, L0, L1 first.
  const std::vector<InEdge> rows {{3, 4, 1, 1},   {5, 6, 1, 1},   {7, 8, 2, 1},    {0, 3, 100, 1}, {4, 5, 100, 1},
                                  {6, 7, 100, 1}, {8, 9, 100, 1}, {1, 7, 100, 1},  {8, 5, 100, 1}, {6, 10, 100, 1},
                                  {2, 5, 100, 1}, {6, 3, 100, 1}, {4, 11, 100, 1}};
  for (int z : {0, 1, 2}) {  // without Z, Z as a tunnel, Z as a ShortestPaths demand
    World w(shared(make_graph(12, rows)));
    const auto id = id_of_input(*w.g);
    auto& a = w.pin(tunnel(*w.g, 0, 9, ids_of(id, {3, 0, 4, 1, 5, 2, 6})));
    auto& b = w.pin(tunnel(*w.g, 1, 10, ids_of(id, {7, 2, 8, 1, 9})));
    Actions list {Replace {a, 2}, Replace {b, 2}};
    if (z == 1) list.push_back(Replace {w.pin(tunnel(*w.g, 2, 11, ids_of(id, {10, 1, 11, 0, 12}))), 0});
    if (z == 2) list.push_back(Replace {w.add({2}, {11}, ecmp()), 0});
    EXPECT_EQ(refused(w, list, SIMULTANEOUS), kPeriodTwo) << "Z form " << z;
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(detail::thread_lossy_counters().rounds, 3) << "Z form " << z;
#endif
  }
  // A member that can carry nothing in a program that is not pooled: Z's program from sZ has a
  // supplied DAG route (sZ's exit to dZ0 below its L0 branch, L0 then X) and a tunnel over that
  // exit; offer 2 gives each entry 1 and the DAG route's unit takes the exit (round robin at sZ,
  // the lower EdgeId first), so the L0 branch never carries anything: the period-2 error at round 3,
  // as without Z. sA 0, sB 1, sZ 2, dZ0 3, n1..n6 4..9, dA 10, dB 11, dZ 12; rows X, L0, L1 first.
  {
    World w(shared(make_graph(13, {{4, 5, 1, 1},   {6, 7, 1, 1},   {8, 9, 2, 1},   {0, 4, 100, 1}, {5, 6, 100, 1},
                                   {7, 8, 100, 1}, {9, 10, 100, 1}, {1, 8, 100, 1}, {9, 6, 100, 1}, {7, 11, 100, 1},
                                   {2, 6, 100, 1}, {7, 4, 100, 1}, {5, 12, 100, 1}, {2, 3, 100, 1}})));
    const auto id = id_of_input(*w.g);
    ASSERT_LT(id[13], id[10]);
    auto& a = w.pin(tunnel(*w.g, 0, 10, ids_of(id, {3, 0, 4, 1, 5, 2, 6})));
    auto& b = w.pin(tunnel(*w.g, 1, 11, ids_of(id, {7, 2, 8, 1, 9})));
    const std::vector<RouteDefinition> routes {
        RouteDefinition {2, dag_of(*w.g, ids_of(id, {10, 1, 11, 0, 12, 13})), 1},
        RouteDefinition {2, std::make_shared<const Path>(*w.g, std::vector<std::vector<EdgeId>> {{id[13]}}), 1}};
    auto& z = w.pin(ForwardingProgram::from_routes(*w.g, std::vector<NodeId> {2}, std::vector<NodeId> {3, 12}, routes));
    EXPECT_EQ(refused(w, {Replace {a, 2}, Replace {b, 2}, Replace {z, 2}}, SIMULTANEOUS), kPeriodTwo);
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(detail::thread_lossy_counters().rounds, 3);
#endif
  }
}

// What can carry nothing changes nothing, a router's dead merges and branches included.
// D-cycle-subdivision-after (A over links 1, 0, 2, 3 with offer 2; B over 2, 3, 0 with offer 1; C over
// 2, 3, 0, 1 with offer 4; links of 2, 4, 2, 2; links 2 and 3 joined directly at node 10) delivers 0,
// 0, 1 at round 4. B becomes a DAG program with a member that carries nothing at node 10:
// - a second pooled source, node 12, whose own link (capacity 9) merges into node 10: with offer 1
//   its share is 0; with offer 2 each source gets 1 and the new link is B's own black hole (it keeps
//   its own router arc: one more arc);
// - a branch from node 10 to a destination of its own, node 13: round robin sends B's unit to link 3
//   (the lower EdgeId), so the branch carries nothing;
// - a merge into node 10 from B's source over node 12: round robin sends B's unit to link 2's side
//   (the lower EdgeId), so that path carries nothing;
// the last two pooled (from_dag) and not (the DAG route beside a tunnel of weight 0 over B's links).
// Links 2 and 3 stay a pair, and the outcome and the structure are those of B without the member.
// sA 0, sB 1, sC 2, d 3; link e from 4 + 2e (link 2 into node 10).
TEST(ApplyActions, MembersThatCarryNothingChangeNothing) {
  const std::vector<InEdge> rows {{4, 5, 2, 1},    {6, 7, 4, 1},    {8, 10, 2, 1},   {10, 11, 2, 1},  {0, 6, 100, 1},
                                  {7, 4, 100, 1},  {5, 8, 100, 1},  {11, 3, 100, 1}, {1, 8, 100, 1},  {11, 4, 100, 1},
                                  {5, 3, 100, 1},  {2, 8, 100, 1},  {5, 6, 100, 1},  {7, 3, 100, 1},  {12, 10, 9, 1},
                                  {10, 13, 100, 1}, {1, 12, 100, 1}};
  struct Out {
    std::string error;
    Amounts delivered;
    std::vector<Amount> carried;
    std::int64_t rounds {-1};
    std::int64_t groups {-1};
    std::vector<EdgeId> barriers;
    std::int64_t dependencies {-1};
  };
  enum Case { kWithout, kZeroShare, kBlackHoled, kDeadBranch, kDeadMerge };
  auto call = [&](Case c, bool pooled) {
    World w(shared(make_graph(14, rows)));
    const auto id = id_of_input(*w.g);
    auto ids = [&](std::vector<int> in) {
      std::vector<EdgeId> out;
      for (int i : in) out.push_back(id[static_cast<std::size_t>(i)]);
      return out;
    };
    EXPECT_LT(id[3], id[15]);
    EXPECT_LT(id[8], id[16]);
    auto& a = w.pin(tunnel(*w.g, 0, 3, ids({4, 1, 5, 0, 6, 2, 3, 7})));
    const auto path = ids({8, 2, 3, 9, 0, 10});
    auto members = path;
    std::vector<NodeId> sources {1};
    std::vector<NodeId> dests {3};
    if (c == kZeroShare || c == kBlackHoled) {
      members.push_back(id[14]);
      sources.push_back(12);
    }
    if (c == kDeadBranch) {
      members.push_back(id[15]);
      dests.push_back(13);
    }
    if (c == kDeadMerge) {
      members.push_back(id[16]);
      members.push_back(id[14]);
    }
    std::shared_ptr<const ForwardingProgram> prog;
    if (pooled) {
      prog = ForwardingProgram::from_dag(*w.g, sources, dests, dag_of(*w.g, members), RouteOptions {});
    } else {
      std::vector<std::vector<EdgeId>> hops;
      for (EdgeId e : path) hops.push_back({e});
      prog = ForwardingProgram::from_routes(*w.g, sources, dests,
                                            std::vector<RouteDefinition> {RouteDefinition {1, dag_of(*w.g, members), 1},
                                                                          RouteDefinition {1, std::make_shared<const Path>(*w.g, hops), 0}});
    }
    auto& b = w.pin(prog, LOSSY, c == kBlackHoled ? no_edges(*w.g, {id[14]}) : nullptr);
    auto& cc = w.pin(tunnel(*w.g, 2, 3, ids({11, 2, 3, 9, 0, 12, 1, 13})));
    Out out;
    std::vector<PlacementResult> res;
    out.error = error_of([&] { res = run(w, Actions {Replace {a, 2}, Replace {b, c == kBlackHoled ? 2 : 1}, Replace {cc, 4}}, SIMULTANEOUS); });
#if NETGRAPH_CORE_PRIVATE_HOOKS
    const auto& k = detail::thread_lossy_counters();
    out.rounds = k.rounds;
    out.groups = k.groups;
    out.barriers = k.barriers;
    out.dependencies = k.dependencies;
#endif
    if (out.error.empty()) {
      out.delivered = delivered(res);
      for (int e : {0, 1, 2, 3}) out.carried.push_back(w.ledger->carried()[static_cast<std::size_t>(id[static_cast<std::size_t>(e)])]);
    }
    return out;
  };
  struct Row {
    Case c;
    bool pooled;
    const char* where;
  };
  for (const auto& row : {Row {kZeroShare, true, "zero-share source"}, Row {kBlackHoled, true, "black-holed source"},
                          Row {kDeadBranch, true, "dead branch, pooled"}, Row {kDeadBranch, false, "dead branch, not pooled"},
                          Row {kDeadMerge, true, "dead merge, pooled"}, Row {kDeadMerge, false, "dead merge, not pooled"}}) {
    const auto without = call(kWithout, row.pooled);
    EXPECT_EQ(without.error, "") << row.where;
    EXPECT_EQ(without.delivered, (Amounts {0, 0, 1})) << row.where;
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(without.rounds, 4) << row.where;
#endif
    const auto with = call(row.c, row.pooled);
    EXPECT_EQ(with.error, without.error) << row.where;
    EXPECT_EQ(with.delivered, without.delivered) << row.where;
    EXPECT_EQ(with.carried, without.carried) << row.where;
    EXPECT_EQ(with.rounds, without.rounds) << row.where;
    EXPECT_EQ(with.groups, without.groups) << row.where;
    EXPECT_EQ(with.barriers, without.barriers) << row.where;
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(with.dependencies, without.dependencies + (row.c == kBlackHoled ? 1 : 0)) << row.where;
#endif
  }
}

// Traffic behind a demand's own black hole can carry nothing: D-dead-tail (A over links 2, 0 with
// offer 1; B over 1, 0, 2 with offer 2; links of 2, 2, 1, 9) gives the period-2 error at round 3
// without Z, and with Z (its own black hole on link 3) over links 3, 0, 1 with offer 1 as a tunnel
// or as a supplied DAG route (a pooled context), or with offer 2 as that DAG route beside a tunnel
// on its own exit (a split route of a program that is not pooled; each entry gets 1): what reaches
// the black hole goes no further down the plan. sA 0, sB 1, sZ 2, d 3; link e from 4 + 2e to
// 5 + 2e (rows 0 to 3), joined by links of 100; Z's exit is row 15.
TEST(ApplyActions, DeadTailsTakeNoPart) {
  for (int z : {0, 1, 2, 3}) {  // without Z, a tunnel, a DAG route, a DAG route beside a tunnel
    World w(shared(make_graph(12, {{4, 5, 2, 1},
                                   {6, 7, 2, 1},
                                   {8, 9, 1, 1},
                                   {10, 11, 9, 1},
                                   {0, 8, 100, 1},
                                   {9, 4, 100, 1},
                                   {5, 3, 100, 1},
                                   {1, 6, 100, 1},
                                   {7, 4, 100, 1},
                                   {5, 8, 100, 1},
                                   {9, 3, 100, 1},
                                   {2, 10, 100, 1},
                                   {11, 4, 100, 1},
                                   {5, 6, 100, 1},
                                   {7, 3, 100, 1},
                                   {2, 3, 100, 1}})));
    const auto id = id_of_input(*w.g);
    auto ids = [&](std::vector<int> in) {
      std::vector<EdgeId> out;
      for (int i : in) out.push_back(id[static_cast<std::size_t>(i)]);
      return out;
    };
    auto& a = w.pin(tunnel(*w.g, 0, 3, ids({4, 2, 5, 0, 6})));
    auto& b = w.pin(tunnel(*w.g, 1, 3, ids({7, 1, 8, 0, 9, 2, 10})));
    Actions list {Replace {a, 1}, Replace {b, 2}};
    const auto zr = ids({11, 3, 12, 0, 13, 1, 14});
    if (z == 1) list.push_back(Replace {w.pin(tunnel(*w.g, 2, 3, zr), LOSSY, no_edges(*w.g, ids({3}))), 1});
    if (z >= 2) {
      std::vector<RouteDefinition> routes {RouteDefinition {2, dag_of(*w.g, zr), 1}};
      if (z == 3) routes.push_back(RouteDefinition {2, std::make_shared<const Path>(*w.g, std::vector<std::vector<EdgeId>> {{id[15]}}), 1});
      const auto prog = ForwardingProgram::from_routes(*w.g, std::vector<NodeId> {2}, std::vector<NodeId> {3}, routes);
      ASSERT_EQ(detail::ProgramAccess::data(*prog).pools_source_inflows, z == 2);
      list.push_back(Replace {w.pin(prog, LOSSY, no_edges(*w.g, ids({3}))), z == 3 ? 2 : 1});
    }
    EXPECT_EQ(refused(w, list, SIMULTANEOUS), kPeriodTwo) << "Z form " << z;
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(detail::thread_lossy_counters().rounds, 3) << "Z form " << z;
#endif
  }
}

// Per-entry boundaries on the one residual.
TEST(ApplyActions, OutputsPerEntry) {
  World w(graph(4, {{0, 2, 10, 1}, {1, 2, 1, 1}, {2, 3, 10, 1}}));
  auto& a = w.pin(tunnel(*w.g, 0, 3, {0, 2}));
  auto& b = w.pin(tunnel(*w.g, 1, 3, {1, 2}));
  const auto r = run(w, {Replace {a, 5}, Replace {b, 1}}, SIMULTANEOUS, OutputOptions {.residual = true, .saturated_boundary_edges = true});
  EXPECT_EQ(*r[0].residual, (Amounts {5, 0, 4}));
  EXPECT_EQ(*r[1].residual, (Amounts {5, 0, 4}));
  EXPECT_EQ(*r[0].saturated_boundary_edges, (std::vector<EdgeId> {}));
  EXPECT_EQ(*r[1].saturated_boundary_edges, (std::vector<EdgeId> {1}));
}

// A controller program's captured weights inject in a list.
TEST(ApplyActions, CapturedWeightsInAList) {
  auto g = graph(3, {{0, 2, 3, 1}, {1, 2, 1, 1}});
  World solver(g);
  const auto prog = solver.add({0, 1}, {2}, global()).replace(4).program;
  World w(g);
  auto& a = w.pin(prog);
  auto& b = w.pin(tunnel(*w.g, 0, 2, {0}));
  const auto r = run(w, {Replace {a, 4}, Replace {b, 1}}, SIMULTANEOUS);
  EXPECT_EQ(delivered(r), (Amounts {4, 0}));
  EXPECT_EQ(r[1].call_edge_drops, (Pairs {{0, 1}}));
}

// The record identities of the whole list are counted first.
TEST(ApplyActions, BatchRecordIdentities) {
  World w(graph(2, {{0, 1, 10, 1}}));
  auto& a = w.add({0}, {1}, ecmp());
  auto& b = w.add({0}, {1}, ecmp());
  auto& c = w.add({0}, {1}, ecmp());
  LedgerAccess::set_next_ids(*w.ledger, std::numeric_limits<std::int64_t>::max() - 2, 3);
  const auto before = w.state();
  std::string what;
  try {
    (void)apply_actions(*w.ledger, Actions {Replace {a, 1}, Replace {b, 1}, Replace {c, 1}}, ORDERED);
  } catch (const std::runtime_error& e) {
    what = e.what();
  }
  EXPECT_EQ(what, "record identity limit exceeded");
  expect_unchanged(w, before, "three new records, two ids left");
  EXPECT_EQ(run(w, {Replace {a, 1}, Replace {b, 1}}, SIMULTANEOUS).size(), 2U);
  EXPECT_EQ(LedgerAccess::next_record_counter(*w.ledger), std::numeric_limits<std::int64_t>::max());
}


// A list's joint network counts every participant: a Pinned program of two one-edge routes on
// PAR(6, 4) (2 entries, 2 occurrences, 4 routers; the 4-unit link drops: 1 barrier, 3 dependencies),
// a Pinned program of two 2-hop tunnels on a diamond (2 entries, 4 occurrences) and a pooled
// ShortestPaths demand of sources 6 and 7 over 6 -> 8, 7 -> 8, 8 -> 9 (2 entries; each source reaches
// two members: 4 occurrences; 4 routers, counted once; 3 dependencies): 6 entries, 10 occurrences,
// 8 routers, 9 events, 6 dependencies. With one limit lowered to its count less one the list is
// refused with the joint count, nothing changed (each participant alone stays within); at the count
// it runs; with every limit at 1 the entries come first. An event or dependency refusal comes before
// the event structures are sized (the hooks' counter reads 0).
TEST(ApplyActions, ForwardingSizesCountEveryParticipant) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  using detail::ForwardingCounter;
  struct Restore {
    ~Restore() { detail::thread_forwarding_limits() = {}; }
  } restore;
  World w(graph(10, {{0, 1, 6, 1}, {0, 1, 4, 1}, {2, 3, 100, 1}, {2, 4, 100, 1}, {3, 5, 100, 1}, {4, 5, 100, 1}, {6, 8, 100, 1},
                     {7, 8, 100, 1}, {8, 9, 100, 1}}));
  const std::vector<RouteDefinition> dags {RouteDefinition {0, dag_of(*w.g, {0}), 1}, RouteDefinition {0, dag_of(*w.g, {1}), 1}};
  auto& a = w.pin(ForwardingProgram::from_routes(*w.g, std::vector<NodeId> {0}, std::vector<NodeId> {1}, dags));
  auto& b = w.pin(tunnels(*w.g, {2}, {5}, {{2, {2, 4}, 1}, {2, {3, 5}, 1}}));
  auto& c = w.add({6, 7}, {9}, ecmp());
  const Actions list {Replace {a, 12}, Replace {b, 2}, Replace {c, 2}};
  struct Count {
    ForwardingCounter counter;
    const char* name;
    std::int64_t n;
  };
  for (const auto& k : std::vector<Count> {{ForwardingCounter::Entry, "entry", 6}, {ForwardingCounter::Occurrence, "occurrence", 10},
                                           {ForwardingCounter::Router, "router", 8}, {ForwardingCounter::Event, "event", 9},
                                           {ForwardingCounter::Dependency, "dependency", 6}}) {
    auto& limits = detail::thread_forwarding_limits();
    limits = {};
    limits[static_cast<std::size_t>(k.counter)] = k.n - 1;
    EXPECT_EQ(refused(w, list, SIMULTANEOUS), std::string("forwarding ") + k.name + " limit exceeded: " + std::to_string(k.n) + " > " +
                                                 std::to_string(k.n - 1));
    if (k.counter == ForwardingCounter::Event || k.counter == ForwardingCounter::Dependency) {
      EXPECT_EQ(detail::thread_lossy_counters().event_structures, 0) << k.name;
    }
    // Each participant alone stays within.
    for (const auto& one : list) EXPECT_NO_THROW((void)apply_actions(*w.ledger, Actions {one}, SIMULTANEOUS)) << k.name;
    limits[static_cast<std::size_t>(k.counter)] = k.n;
    EXPECT_NO_THROW((void)run(w, list, SIMULTANEOUS)) << k.name;
    EXPECT_EQ(detail::thread_lossy_counters().event_structures, 1) << k.name;
  }
  detail::thread_forwarding_limits() = {1, 1, 1, 1, 1};
  EXPECT_EQ(refused(w, list, SIMULTANEOUS), "forwarding entry limit exceeded: 6 > 1");
#else
  GTEST_SKIP() << "a lowered limit needs the private hooks";
#endif
}

// The timing rows' ring: routers 0..63, ring links i -> i + 1 of 3 units, heads 64..127 with
// wide feeders, wide exits to 128; tunnel i over ring links i and i + 1. Every link carries two
// tunnels of 3: 64 barriers in one cyclic group. As one 64-head program (replace(192): 3 per head)
// and as 64 one-tunnel policies in one SIMULTANEOUS list: both settle at round 6, every head
// delivering 1, every link carrying 3.
TEST(ApplyActions, RingSettlesAtRoundSix) {
  constexpr int n = 64;
  const Amount wide = Amount {1} << 40;
  std::vector<InEdge> rows;
  for (int i = 0; i < n; ++i) rows.push_back({i, (i + 1) % n, 3, 1});
  for (int i = 0; i < n; ++i) rows.push_back({i, 2 * n, wide, 1});
  for (int i = 0; i < n; ++i) rows.push_back({n + i, i, wide, 1});
  const auto g = shared(make_graph(2 * n + 1, rows));
  const auto id = id_of_input(*g);
  auto tunnel_of = [&](int i) {
    return std::tuple<NodeId, std::vector<EdgeId>, Amount> {
        n + i, std::vector<EdgeId> {id[static_cast<std::size_t>(2 * n + i)], id[static_cast<std::size_t>(i)], id[static_cast<std::size_t>((i + 1) % n)],
                                    id[static_cast<std::size_t>(n + (i + 2) % n)]},
        1};
  };
  std::vector<NodeId> heads;
  std::vector<std::tuple<NodeId, std::vector<EdgeId>, Amount>> all;
  for (int i = 0; i < n; ++i) heads.push_back(n + i), all.push_back(tunnel_of(i));
  World one(g);
  auto& program = one.pin(tunnels(*g, heads, {2 * n}, all));
  const auto r = program.replace(192);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().groups, 1);
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 6);
#endif
  Amounts per_head;
  for (const auto& t : r.by_tunnel) per_head.push_back(t.delivered);
  EXPECT_EQ(per_head, Amounts(n, 1));
  World many(g);
  Actions list;
  for (int i = 0; i < n; ++i) list.push_back(Replace {many.pin(tunnels(*g, {n + i}, {2 * n}, {tunnel_of(i)})), 3});
  const auto res = run(many, list, SIMULTANEOUS);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().groups, 1);
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 6);
#endif
  EXPECT_EQ(delivered(res), Amounts(n, 1));
  EXPECT_EQ(many.ledger->carried(), one.ledger->carried());
  for (int i = 0; i < n; ++i) EXPECT_EQ(one.ledger->carried()[static_cast<std::size_t>(id[static_cast<std::size_t>(i)])], 3);
}

// ---------------------------------------------------------------- generators

// A world built from a seed (built twice, identically): a random DAG on nodes 0..5 (sources 0 and 1,
// destination 5; nonuniform costs, capacities 0 to 8) plus a separate link 6 -> 7, and one policy of
// each kind: ShortestPaths EQUAL LY (two sources), EQUAL LL, WEIGHTED LY with capacity, explicit and
// residual-snapshot weights, EQUAL LY under a mask; OPTIMIZED LL; CentralizedTE LL; T1 LY (two tunnels
// per head); T2 LY; Pinned LY over tunnels, over split routes, over from_dag under a mask, over a
// controller program with black holes; Pinned LY on the separate link. Then scalar calls, an orphan
// and an outside owner. Empty when a construction fails (both copies fail alike).
struct Generated {
  std::unique_ptr<World> w;
  std::vector<std::size_t> lossy_replace;  // policies a SIMULTANEOUS Replace accepts
  std::vector<std::size_t> forwardable;    // LOSSY computing policies (Forward after an install)
  std::size_t separate {0};                // the Pinned policy on 6 -> 7
};

// `graph`: the first copy's graph, shared by the second (records compare their routes by graph).
Generated generate(std::uint32_t seed, std::shared_ptr<const Graph> graph = nullptr) {
  std::mt19937 rng(seed);
  std::vector<InEdge> rows;
  for (NodeId u = 0; u < 6; ++u) {
    for (NodeId v = u + 1; v < 6; ++v) {
      if (rng() % 3 != 0) rows.push_back({u, v, static_cast<Amount>(rng() % 9), 1 + static_cast<Cost>(rng() % 3)});
    }
  }
  rows.push_back({6, 7, 10, 1});
  Generated out;
  try {
    out.w = std::make_unique<World>(graph ? std::move(graph) : shared(make_graph(8, rows)));
    World& w = *out.w;
    const Graph& g = *w.g;
    const auto E = static_cast<std::size_t>(g.num_edges());
    const EdgeId sep = g.num_edges() - 1;  // 6 -> 7: the largest tail
    auto random_edge = [&] { return static_cast<EdgeId>(rng() % (E - 1)); };
    // A random path from h to 5 (empty when none).
    auto path_from = [&](NodeId h) {
      std::vector<EdgeId> path;
      NodeId u = h;
      while (u != 5) {
        std::vector<EdgeId> outs;
        for (EdgeId e = 0; e < sep; ++e) {
          if (g.edge_sources()[static_cast<std::size_t>(e)] == u) outs.push_back(e);
        }
        if (outs.empty()) return std::vector<EdgeId> {};
        const EdgeId e = outs[rng() % outs.size()];
        path.push_back(e);
        u = g.edge_destinations()[static_cast<std::size_t>(e)];
      }
      return path;
    };
    const std::vector<NodeId> S {0, 1};
    const std::vector<NodeId> D {5};
    w.add(S, D, ecmp());
    out.lossy_replace.push_back(0);
    out.forwardable.push_back(0);
    w.add({0}, D, ecmp(LOSSLESS));
    {
      ShortestPaths cap = shortest_paths(RouteAllocation::Weighted);
      w.add(S, D, FlowPolicyConfig {.routing = cap, .admission = LOSSY});
      out.lossy_replace.push_back(2);
      out.forwardable.push_back(2);
      ShortestPaths ex = shortest_paths(RouteAllocation::Weighted);
      ex.weight_basis = WeightBasis::Explicit;
      std::vector<Amount> weights(E);
      for (auto& x : weights) x = static_cast<Amount>(rng() % 4);
      ex.edge_weights = std::make_shared<const std::vector<Amount>>(weights);
      w.add(S, D, FlowPolicyConfig {.routing = ex, .admission = LOSSY});
      out.lossy_replace.push_back(3);
      out.forwardable.push_back(3);
      ShortestPaths snap = shortest_paths(RouteAllocation::Weighted);
      snap.weight_basis = WeightBasis::ResidualSnapshot;
      w.add({0}, D, FlowPolicyConfig {.routing = snap, .admission = LOSSY});
      out.lossy_replace.push_back(4);
      out.forwardable.push_back(4);
    }
    w.add(S, D, ecmp(), no_edges(g, {random_edge()}));
    out.lossy_replace.push_back(5);
    out.forwardable.push_back(5);
    w.add(S, D, FlowPolicyConfig {.routing = shortest_paths(RouteAllocation::Optimized), .admission = LOSSLESS});
    w.add(S, D, global());
    w.add(S, D, FlowPolicyConfig {.routing = DistributedTE {TunnelCount {2}}, .admission = LOSSY});
    out.forwardable.push_back(8);
    w.add(S, D, t2());
    out.forwardable.push_back(9);
    // Pinned programs.
    std::vector<std::tuple<NodeId, std::vector<EdgeId>, Amount>> specs;
    for (NodeId h : S) {
      for (int k = 0; k < 2; ++k) {
        auto path = path_from(h);
        if (!path.empty()) specs.emplace_back(h, path, 1 + static_cast<Amount>(rng() % 3));
      }
    }
    if (specs.empty()) return {};
    std::vector<NodeId> heads;
    for (const auto& x : specs) heads.push_back(std::get<0>(x));
    std::sort(heads.begin(), heads.end());
    heads.erase(std::unique(heads.begin(), heads.end()), heads.end());
    w.pin(tunnels(g, heads, D, specs));
    out.lossy_replace.push_back(10);
    std::vector<RouteDefinition> split;
    for (NodeId h : heads) split.push_back(RouteDefinition {h, all_edges(g), 1 + static_cast<Amount>(rng() % 2)});
    w.pin(ForwardingProgram::from_routes(g, heads, D, split));
    out.lossy_replace.push_back(11);
    w.pin(ForwardingProgram::from_dag(g, S, D, all_edges(g), RouteOptions {}), LOSSY, nullptr, no_edges(g, {random_edge()}));
    out.lossy_replace.push_back(12);
    {
      World solver(w.g);
      const auto prog = solver.add(S, D, global()).replace(6).program;
      w.pin(prog, LOSSY, no_edges(g, {random_edge()}));
      out.lossy_replace.push_back(13);
    }
    w.pin(tunnel(g, 6, 7, {sep}));
    out.separate = 14;
    // State before the list: scalar calls, an orphan, an outside owner.
    (void)w.p[0]->place(static_cast<Amount>(rng() % 6));
    (void)w.p[2]->replace(static_cast<Amount>(rng() % 6));
    (void)w.p[6]->place(static_cast<Amount>(rng() % 6));
    (void)w.p[8]->replace(2, std::optional<Amount> {2});
    (void)w.p[9]->replace(2);
    (void)w.p[10]->replace(static_cast<Amount>(rng() % 6));
    w.add(S, D, ecmp()).place(1 + static_cast<Amount>(rng() % 3));  // the outside owner (index 15)
    auto& orphan = w.add({0}, D, ecmp());
    (void)orphan.place(1 + static_cast<Amount>(rng() % 3));
    w.p.pop_back();
  } catch (const std::invalid_argument&) {
    return {};
  }
  return out;
}

// One random valid ORDERED action on policy i of a generated world (a fresh draw of the generator).
Action ordered_action(std::mt19937& rng, World& w, std::size_t i) {
  FlowPolicy& p = *w.p[i];
  const auto& cfg = p.config();
  const bool tunnel_model = std::holds_alternative<DistributedTE>(cfg.routing) ||
                            (std::holds_alternative<CentralizedTE>(cfg.routing) && std::get<CentralizedTE>(cfg.routing).tunnel_count);
  const bool pinned_model = std::holds_alternative<Pinned>(cfg.routing);
  const bool lossy = cfg.admission == LOSSY;
  const Amount v = static_cast<Amount>(rng() % 13);
  const auto kind = rng() % 3;
  if (kind == 2 && !pinned_model && p.program() != nullptr) {
    const Eligibility* bh = nullptr;
    if (rng() % 3 == 0) {
      w.keep.push_back(no_edges(*w.g, {static_cast<EdgeId>(rng() % (w.g->num_edges() - 1))}));
      bh = w.keep.back().get();
    }
    return Forward {p, (!lossy && rng() % 4 == 0) ? INF : std::optional<Amount> {v}, bh};
  }
  if (kind == 0 && !tunnel_model && !detail::PolicyAccess::engine(p).epoch.infinite) return Place {p, v};
  if (std::holds_alternative<DistributedTE>(cfg.routing)) return Replace {p, v, std::optional<Amount> {v + static_cast<Amount>(rng() % 3)}};
  return Replace {p, (!lossy && rng() % 4 == 0) ? INF : std::optional<Amount> {v}};
}

// The same action as a scalar call.
PlacementResult scalar(const Action& a, const OutputOptions& o) {
  return std::visit(
      [&](const auto& x) -> PlacementResult {
        using T = std::decay_t<decltype(x)>;
        if constexpr (std::is_same_v<T, Place>) return x.policy.get().place(x.offer, o);
        if constexpr (std::is_same_v<T, Replace>) return x.policy.get().replace(x.offer, x.design_demand, o);
        if constexpr (std::is_same_v<T, Forward>) return x.policy.get().forward(x.offer, x.black_holes, o);
      },
      a);
}

// The action `a` of world `from` aimed at the same policy (and black holes) of world `to`.
Action retarget(const Action& a, World& from, World& to) {
  auto index = [&](const FlowPolicy& p) {
    for (std::size_t i = 0; i < from.p.size(); ++i) {
      if (from.p[i].get() == &p) return i;
    }
    return from.p.size();
  };
  return std::visit(
      [&](const auto& x) -> Action {
        using T = std::decay_t<decltype(x)>;
        T y = x;
        y.policy = *to.p[index(x.policy.get())];
        if constexpr (std::is_same_v<T, Forward>) {
          if (x.black_holes != nullptr) {
            to.keep.push_back(std::make_shared<const Eligibility>(*x.black_holes));
            y.black_holes = to.keep.back().get();
          }
        }
        return y;
      },
      a);
}

void expect_same_worlds(const World& a, const World& b, const std::string& where) {
  expect_same(ledger_state(*a.ledger), ledger_state(*b.ledger), where);
  ASSERT_EQ(a.p.size(), b.p.size()) << where;
  for (std::size_t i = 0; i < a.p.size(); ++i) expect_same(epoch(*a.p[i]), epoch(*b.p[i]), where + ", policy " + std::to_string(i));
}

// An ORDERED list gives exactly the results, records (FlowIds included), ledger and epochs of
// the same scalar calls in the same order on an identical world.
TEST(ApplyActions, OrderedEqualsScalarCalls) {
  int checked = 0;
  int forwards = 0;
  int holds = 0;
  for (std::uint32_t seed = 1; seed <= 200; ++seed) {
    auto a = generate(seed);
    if (!a.w) continue;
    auto b = generate(seed, a.w->g);
    std::mt19937 rng(seed * 7919U);
    std::vector<std::size_t> order(a.w->p.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::shuffle(order.begin(), order.end(), rng);
    order.resize(2 + rng() % 6);
    Actions list;
    for (auto i : order) list.push_back(ordered_action(rng, *a.w, i));
    const std::string where = "seed " + std::to_string(seed);
    std::vector<PlacementResult> got;
    std::string error;
    try {
      got = apply_actions(*a.w->ledger, list, ORDERED, kAll);
    } catch (const std::invalid_argument& e) {
      error = e.what();
    }
    std::string scalar_error;
    std::vector<PlacementResult> want;
    for (const auto& x : list) {
      try {
        want.push_back(scalar(retarget(x, *a.w, *b.w), kAll));
      } catch (const std::invalid_argument& e) {
        scalar_error = e.what();
        break;
      }
    }
    ASSERT_EQ(error, scalar_error) << where;
    if (!error.empty()) continue;
    ASSERT_EQ(got.size(), want.size()) << where;
    for (std::size_t k = 0; k < got.size(); ++k) expect_same(got[k], want[k], where + ", action " + std::to_string(k));
    expect_same_worlds(*a.w, *b.w, where);
    for (const auto& x : list) forwards += std::holds_alternative<Forward>(x) ? 1 : 0;
    holds += a.w->ledger->reserved() != Amounts(a.w->ledger->reserved().size(), 0) ? 1 : 0;
    ++checked;
    if (::testing::Test::HasFailure()) return;
  }
  RecordProperty("checked", checked);
  EXPECT_GT(checked, 150);
  EXPECT_GT(forwards, 80);
  EXPECT_GT(holds, 30);
}

// One SIMULTANEOUS participant (and a second one on the separate link) equals its scalar
// replace or forward on an identical world: results, records, ledger and epochs.
TEST(ApplyActions, OneParticipantEqualsTheScalarCall) {
  int checked = 0;
  int pooled = 0;
  int black = 0;
  for (std::uint32_t seed = 1; seed <= 300; ++seed) {
    auto a = generate(seed);
    if (!a.w) continue;
    auto b = generate(seed, a.w->g);
    std::mt19937 rng(seed * 104729U);
    const auto lossy_forward = a.forwardable;
    const bool forward = rng() % 2 == 0;
    const std::size_t i = forward ? lossy_forward[rng() % lossy_forward.size()] : a.lossy_replace[rng() % a.lossy_replace.size()];
    const Amount v = static_cast<Amount>(rng() % 13);
    Actions list;
    FlowPolicy& p = *a.w->p[i];
    if (forward) {
      if (p.program() == nullptr) (void)p.replace(1), (void)b.w->p[i]->replace(1);
      const Eligibility* bh = nullptr;
      if (rng() % 2 == 0) {
        a.w->keep.push_back(no_edges(*a.w->g, {static_cast<EdgeId>(rng() % (a.w->g->num_edges() - 1))}));
        bh = a.w->keep.back().get();
        ++black;
      }
      list.push_back(Forward {p, v, bh});
    } else {
      list.push_back(Replace {p, v});
      pooled += i <= 5 ? 1 : 0;
    }
    const bool pair = rng() % 2 == 0;
    if (pair) list.push_back(Replace {*a.w->p[a.separate], static_cast<Amount>(rng() % 13)});
    const std::string where = "seed " + std::to_string(seed) + " policy " + std::to_string(i);
    const auto got = apply_actions(*a.w->ledger, list, SIMULTANEOUS, kAll);
    std::vector<PlacementResult> want;
    for (const auto& x : list) want.push_back(scalar(retarget(x, *a.w, *b.w), kAll));
    // The scalar calls' outputs read the ledger after each; the list's after both: with the pair,
    // each entry's residual is the one after the list (the other participant is disjoint).
    if (pair) {
      want[0].residual = want[1].residual;
      if (want[0].saturated_boundary_edges) {
        const auto& pe = detail::PolicyAccess::engine(*b.w->p[i]);
        want[0].saturated_boundary_edges = detail::saturated_boundary_edges_of(*b.w->g, *want[1].residual, pe.sources, pe.eligibility);
      }
    }
    for (std::size_t k = 0; k < got.size(); ++k) expect_same(got[k], want[k], where + ", entry " + std::to_string(k));
    expect_same_worlds(*a.w, *b.w, where);
    ++checked;
    if (::testing::Test::HasFailure()) return;
  }
  RecordProperty("checked", checked);
  EXPECT_GT(checked, 200);
  EXPECT_GT(pooled, 60);
  EXPECT_GT(black, 40);
}

} // namespace
