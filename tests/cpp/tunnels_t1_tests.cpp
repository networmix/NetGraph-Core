// T1 (CapacityAware(Tunnels(N))): every head sets up its own strict tunnels by CSPF (contract
// section 6), the two ledger projections (10.1, 10.3), T1's programs, replays and stale forwarding
// (8.1, 8.4, 8.6), records, rows, flags and terminations (10.4, 14.1, 14.2), T1's evaluators (12)
// and the numeric ranges of its sums. Every expected value is derived by hand or copied from the
// contract's worked examples, never from a build; the property test compares random small cases
// with t1_reference (test_utils.hpp), an independent restatement of 6.2 to 6.5 that calls none of
// the library's searches or apportionment helpers.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/program.hpp"
#include "state_comparators.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

constexpr Admission LOSSLESS = Admission::Lossless;
constexpr Admission LOSSY = Admission::Lossy;
constexpr TunnelStatus UP = TunnelStatus::Up;
constexpr TunnelStatus NB = TunnelStatus::NoBandwidth;
constexpr TunnelStatus CC = TunnelStatus::CostCeiling;
constexpr TunnelStatus FAILED = TunnelStatus::Failed;
constexpr TunnelStatus BLACKHOLED = TunnelStatus::Blackholed;
const std::optional<Amount> INF;  // +inf

using Amounts = std::vector<Amount>;
using Edges = std::vector<EdgeId>;
using Pairs = std::vector<std::pair<EdgeId, Amount>>;
// A by_tunnel row: (index, head, status, assigned, delivered, dropped).
using Row = std::tuple<std::int32_t, NodeId, TunnelStatus, Amount, Amount, Amount>;
using Rows = std::vector<Row>;
// A by_source row: (source, offered or -1 when Lossless, admitted, delivered, source_loss, dropped).
using SourceRow = std::tuple<NodeId, Amount, Amount, Amount, Amount, Amount>;

// A graph whose rows are given in EdgeId order (the ledger's eK notation).
std::shared_ptr<const Graph> ordered(std::int32_t n, const std::vector<InEdge>& rows) {
  auto g = shared(make_graph(n, rows));
  for (std::size_t i = 0; i < rows.size(); ++i) EXPECT_EQ(g->input_index()[i], static_cast<std::int32_t>(i));
  return g;
}

FlowPolicyConfig t1(std::int64_t n, Admission a = LOSSLESS, std::optional<CostCeiling> ceiling = std::nullopt,
                    Amount minimum = 0, std::optional<PathSelection> sel = std::nullopt) {
  return FlowPolicyConfig {.routing = CapacityAware {Tunnels {n}, ceiling, minimum}, .admission = a, .selection = sel};
}

PathSelection sel(EdgeTieBreak t, EdgePriority p, bool multi_edge = true) { return PathSelection(t, p, false, multi_edge); }
constexpr EdgeTieBreak DET = EdgeTieBreak::Deterministic;
constexpr EdgeTieBreak PHP = EdgeTieBreak::PreferHigherPriority;

Amount lo(const Counter128& c) {
  EXPECT_EQ(c.hi, 0U);
  return static_cast<Amount>(c.lo);
}
Amount lo(const std::optional<Counter128>& c) { return c.has_value() ? lo(*c) : -1; }

Rows rows(const PlacementResult& r) {
  Rows out;
  for (const auto& t : r.by_tunnel) out.emplace_back(t.index, t.head, t.status, lo(t.assigned), t.delivered, lo(t.dropped));
  return out;
}
std::vector<SourceRow> source_rows(const PlacementResult& r) {
  std::vector<SourceRow> out;
  for (const auto& s : r.by_source) {
    out.emplace_back(s.source, lo(s.offered), lo(s.admitted), s.delivered, lo(s.source_loss), lo(s.dropped));
  }
  return out;
}
std::vector<TunnelStatus> statuses(const ForwardingProgram& p) {
  std::vector<TunnelStatus> out;
  for (const auto& t : p.tunnels()) out.push_back(t.status);
  return out;
}
// Each tunnel's path (empty when down), in index order.
std::vector<Edges> routes(const ForwardingProgram& p) {
  std::vector<Edges> out;
  for (const auto& t : p.tunnels()) out.push_back(t.edges.value_or(Edges {}));
  return out;
}
Amounts requests(const ForwardingProgram& p) {
  Amounts out;
  for (const auto& t : p.tunnels()) out.push_back(t.requested);
  return out;
}
std::vector<std::optional<std::int32_t>> entries(const std::vector<FlowRecord>& records) {
  std::vector<std::optional<std::int32_t>> out;
  for (const auto& r : records) out.push_back(r.entry);
  return out;
}
std::vector<Pairs> holds(const std::vector<FlowRecord>& records) {
  std::vector<Pairs> out;
  for (const auto& r : records) out.push_back(r.reserved_edges);
  return out;
}

// One policy on its own ledger.
struct Pol {
  std::shared_ptr<FlowLedger> ledger;
  std::unique_ptr<FlowPolicy> p;
  FlowPolicy* operator->() const { return p.get(); }
};
Pol make(const std::shared_ptr<const Graph>& g, std::vector<NodeId> s, std::vector<NodeId> d, const FlowPolicyConfig& c,
         const Eligibility* elig = nullptr, std::shared_ptr<FlowLedger> ledger = nullptr) {
  Pol out;
  out.ledger = ledger ? std::move(ledger) : std::make_shared<FlowLedger>(g);
  out.p = std::make_unique<FlowPolicy>(*out.ledger, s, d, c, elig);
  return out;
}

// The graphs of the hand derivations.
std::shared_ptr<const Graph> one10() { return ordered(2, {{0, 1, 10, 1}}); }
std::shared_ptr<const Graph> twoh(Amount a, Amount b) { return ordered(3, {{0, 2, a, 1}, {1, 2, b, 1}}); }
std::shared_ptr<const Graph> par612() { return ordered(2, {{0, 1, 6, 1}, {0, 1, 6, 2}}); }
std::shared_ptr<const Graph> diam5() {
  return ordered(4, {{0, 1, 6, 1}, {0, 2, 6, 3}, {1, 2, 4, 1}, {1, 3, 6, 3}, {2, 3, 6, 1}});
}
std::shared_ptr<const Graph> x8() { return ordered(3, {{0, 2, 0, 1}, {1, 2, 8, 3}}); }
std::shared_ptr<const Graph> zl() { return ordered(2, {{0, 1, 0, 1}, {0, 1, 10, 2}}); }

const OutputOptions kHolds {.reservation_edge_flows = true};

} // namespace

// ---------------------------------------------------------------- the contract's examples

// Every T1 example of the contract (6.2 to 6.4, 3.4), reproduced with the example's own inputs
// and expected values: the requests in round order, every tunnel's status and path, the holds per
// edge and in all, the lossless bound (Lossless replace(+inf, D)), and with an offer O the Lossy
// shares per entry, delivery and source loss of replace(O, D).
TEST(T1Setup, ContractExamples) {
  struct Example {
    const char* id;
    std::int32_t n;
    std::vector<InEdge> edges;
    std::vector<NodeId> heads;
    NodeId dst;
    Amount D;
    std::int64_t N;
    std::optional<CostCeiling> ceiling;
    Amounts requests;
    std::vector<Edges> routes;
    std::vector<TunnelStatus> status;
    Amount reserved;
    Amounts reserved_edges;
    Amount bound;
    std::optional<Amount> offer {};
    Amounts shares {};  // Lossy offered per entry
    Amount delivery {0};
    Amount source_loss {0};
  };
  const std::vector<InEdge> diam = {{0, 1, 6, 1}, {0, 2, 6, 3}, {1, 2, 4, 1}, {1, 3, 6, 3}, {2, 3, 6, 1}};
  const std::vector<InEdge> g1 = {{0, 1, 8192, 1}, {0, 2, 40960, 1}, {1, 3, 8192, 1}, {2, 3, 40960, 1}};
  std::vector<Edges> branch_routes(16, Edges {1, 3});
  branch_routes[0] = branch_routes[1] = Edges {0, 2};
  branch_routes[15] = {};
  std::vector<TunnelStatus> branch_status(16, UP);
  branch_status[15] = NB;
  std::vector<Edges> large_routes(16, Edges {});
  large_routes[0] = Edges {1, 3};
  std::vector<TunnelStatus> large_status(16, NB);
  large_status[0] = UP;
  Amounts branch_shares(16, 3277);
  branch_shares[12] = branch_shares[13] = branch_shares[14] = 3276;
  branch_shares[15] = 0;
  const std::vector<Example> examples = {
      {"four tunnels on two parallel links",
       2, {{0, 1, 6, 1}, {0, 1, 6, 2}}, {0}, 1, 12, 4, std::nullopt, {3, 3, 3, 3}, {{0}, {0}, {1}, {1}},
       {UP, UP, UP, UP}, 12, {6, 6}, 12, 13, {4, 3, 3, 3}, 12, 0},
      {"uneven requests 3 and 2",
       2, {{0, 1, 3, 1}, {0, 1, 2, 2}}, {0}, 1, 5, 2, std::nullopt, {3, 2}, {{0}, {1}}, {UP, UP}, 5, {3, 2}, 5, 5,
       {3, 2}, 5, 0},
      {"a share below the tunnel count",
       2, {{0, 1, 6, 1}, {0, 1, 6, 2}}, {0}, 1, 1, 2, std::nullopt, {1, 0}, {{0}, {0}}, {UP, UP}, 1,
       {1, 0}, 6, 1, {1, 0}, 1, 0},
      {"carried above the holds",
       2, {{0, 1, 5, 1}, {0, 1, 5, 1}}, {0}, 1, 8, 2, std::nullopt, {4, 4}, {{0}, {1}}, {UP, UP},
       8, {4, 4}, 10, 10, {5, 5}, 10, 0},
      {"a zero request on a zero-capacity link",
       2, {{0, 1, 0, 1}}, {0}, 1, 0, 1, std::nullopt, {0}, {{0}}, {UP}, 0, {0}, 0, 1, {1}, 0, 0},
      {"the last tunnel without bandwidth",
       2, {{0, 1, 2, 1}, {0, 1, 10, 2}}, {0}, 1, 12, 4, std::nullopt, {3, 3, 3, 3}, {{1}, {1}, {1}, {}},
       {UP, UP, UP, NB}, 9, {0, 9}, 10, 12, {4, 4, 4, 0}, 10, 0},
      {"two heads sharing a link",
       4, {{0, 2, 6, 1}, {1, 2, 6, 1}, {2, 3, 6, 1}}, {0, 1}, 3, 12, 2, std::nullopt, {3, 3, 3, 3},
       {{0, 2}, {1, 2}, {}, {}}, {UP, UP, NB, NB}, 6, {3, 3, 6}, 6, 12, {6, 6, 0, 0}, 6, 0},
      {"three heads, one with a zero request",
       4, {{0, 3, 20, 1}, {1, 3, 20, 1}, {2, 3, 20, 1}}, {0, 1, 2}, 3, 2, 1, std::nullopt, {1, 1, 0},
       {{0}, {1}, {2}}, {UP, UP, UP}, 2, {1, 1, 0}, 60, 3, {1, 1, 1}, 3, 0},
      {"two equal links, one tunnel each",
       2, {{0, 1, 5, 1}, {0, 1, 5, 1}}, {0}, 1, 8, 2, std::nullopt, {4, 4}, {{0}, {1}}, {UP, UP}, 8, {4, 4}, 10,
       8, {4, 4}, 8, 0},
      {"a head behind a zero-capacity link",
       3, {{0, 2, 0, 1}, {1, 2, 8, 3}}, {0, 1}, 2, 4, 1, std::nullopt, {2, 2}, {{}, {1}}, {NB, UP}, 2, {0, 2}, 0,
       4, {0, 2}, 2, 2},
      {"a head behind a zero-capacity link, Absolute(3)",
       3, {{0, 2, 0, 1}, {1, 2, 8, 3}}, {0, 1}, 2, 4, 1, Absolute {3}, {2, 2}, {{}, {1}}, {NB, UP}, 2,
       {0, 2}, 0},
      {"a head behind a zero-capacity link, Stretch(2)",
       3, {{0, 2, 0, 1}, {1, 2, 8, 3}}, {0, 1}, 2, 4, 1, Stretch {2.0}, {2, 2}, {{}, {}}, {NB, CC}, 0,
       {0, 0}, 0},
      {"Stretch(1.2) just below the next path",
       2, {{0, 1, 0, 5}, {0, 1, 1, 6}}, {0}, 1, 1, 1, Stretch {1.2}, {1}, {{}}, {CC}, 0, {0, 0}, 0,
       1, {0}, 0, 1},
      {"the diamond, design 8",
       4, diam, {0}, 3, 8, 2, std::nullopt, {4, 4}, {{0, 2, 4}, {}}, {UP, NB}, 4, {4, 0, 4, 0, 4}, 4},
      {"the diamond, design 10",
       4, diam, {0}, 3, 10, 2, std::nullopt, {5, 5}, {{0, 3}, {1, 4}}, {UP, UP}, 10, {5, 5, 0, 5, 5}, 12},
      {"the diamond, design 14",
       4, diam, {0}, 3, 14, 2, std::nullopt, {7, 7}, {{}, {}}, {NB, NB}, 0, {0, 0, 0, 0, 0}, 0},
      {"two heads, the second with bandwidth",
       3, {{0, 2, 2, 1}, {1, 2, 10, 1}}, {0, 1}, 2, 12, 2, std::nullopt, {3, 3, 3, 3}, {{}, {1}, {}, {1}},
       {NB, UP, NB, UP}, 6, {0, 6}, 0, 12, {0, 3, 0, 3}, 6, 6},
      {"more tunnels than units",
       2, {{0, 1, 10, 1}}, {0}, 1, 3, 4, std::nullopt, {1, 1, 1, 0}, {{0}, {0}, {0}, {0}}, {UP, UP, UP, UP}, 3,
       {3}, 10, 3, {1, 1, 1, 0}, 3, 0},
      {"an idle tunnel on a zero-capacity link",
       2, {{0, 1, 0, 1}, {0, 1, 10, 2}}, {0}, 1, 1, 2, std::nullopt, {1, 0}, {{1}, {0}}, {UP, UP}, 1,
       {0, 1}, 1, 2, {1, 1}, 1, 0},
      {"a zero design",
       2, {{0, 1, 10, 1}}, {0}, 1, 0, 2, std::nullopt, {0, 0}, {{0}, {0}}, {UP, UP}, 0, {0}, 10, 4,
       {2, 2}, 4, 0},
      // The sixteen tunnels' Lossy delivery (not in the example), by hand: tunnels 0 and 1 carry their 3277 + 3277 =
      // 6554 within 8192; tunnels 2 to 14 bring 10 * 3277 + 3 * 3276 = 42598 to the 40960 link: 47514.
      {"sixteen tunnels on two branches",
       4, g1, {0}, 3, 49152, 16, std::nullopt, Amounts(16, 3072), branch_routes, branch_status, 46080,
       {6144, 39936, 6144, 39936}, 47262, 49152, branch_shares, 47514, 0},
      {"sixteen tunnels, a design above the network",
       4, g1, {0}, 3, 409600, 16, std::nullopt, Amounts(16, 25600), large_routes, large_status,
       25600, {0, 25600, 0, 25600}, 40960},
  };
  for (const auto& x : examples) {
    const std::string where = x.id;
    const auto g = ordered(x.n, x.edges);
    // The setup and the lossless bound.
    auto ll = make(g, x.heads, {x.dst}, t1(x.N, LOSSLESS, x.ceiling));
    const auto r = ll->replace(INF, x.D, kHolds);
    ASSERT_TRUE(r.program) << where;
    EXPECT_EQ(requests(*r.program), x.requests) << where;
    EXPECT_EQ(routes(*r.program), x.routes) << where;
    EXPECT_EQ(statuses(*r.program), x.status) << where;
    ASSERT_TRUE(r.plan.has_value()) << where;
    EXPECT_EQ(r.plan->reserved_total, x.reserved) << where;
    EXPECT_EQ(*r.reservation_edge_flows, x.reserved_edges) << where;
    EXPECT_EQ(ll.ledger->reserved(), x.reserved_edges) << where;
    EXPECT_EQ(r.placed, x.bound) << where;
    EXPECT_EQ(r.termination, Termination::Exhausted) << where;  // an unlimited query
    if (!x.offer.has_value()) continue;
    auto ly = make(g, x.heads, {x.dst}, t1(x.N, LOSSY, x.ceiling));
    const auto s = ly->replace(*x.offer, x.D);
    Amounts shares;
    for (const auto& t : s.by_tunnel) shares.push_back(lo(t.assigned));
    EXPECT_EQ(shares, x.shares) << where;
    EXPECT_EQ(s.placed, x.delivery) << where;
    EXPECT_EQ(s.call_source_loss, x.source_loss) << where;
    EXPECT_EQ(s.termination, Termination::Complete) << where;
    // The setup does not depend on the offer (6.5): the same program as the Lossless call's.
    EXPECT_EQ(routes(*s.program), x.routes) << where;
    EXPECT_EQ(statuses(*s.program), x.status) << where;
  }
}

// The split of 6.2: D = 5 over two heads and N = 2: shares RR(5, 2) = (3, 2); requests in round
// order 2, 1, 1, 1. And TWOH(10, 10), D = 10: shares (5, 5), requests 3, 3, 2, 2.
TEST(T1Setup, RoundOrder) {
  auto p = make(twoh(10, 10), {0, 1}, {2}, t1(2));
  const auto split = p->replace(5, 5);
  EXPECT_EQ(requests(*split.program), (Amounts {2, 1, 1, 1}));
  std::vector<NodeId> heads;
  std::vector<std::int32_t> indices, locals;
  for (const auto& t : split.program->tunnels()) {
    heads.push_back(t.head);
    indices.push_back(t.index);
    locals.push_back(t.local_index);
  }
  EXPECT_EQ(heads, v({0, 1, 0, 1}));
  EXPECT_EQ(indices, (std::vector<std::int32_t> {0, 1, 2, 3}));
  EXPECT_EQ(locals, (std::vector<std::int32_t> {0, 0, 1, 1}));
  EXPECT_EQ(split.plan->head_design, (std::vector<std::pair<NodeId, Amount>> {{0, 3}, {1, 2}}));
  EXPECT_EQ(requests(*p->replace(10, 10).program), (Amounts {3, 3, 2, 2}));
}

// ---------------------------------------------------------------- setup (6.4)

// PARX (e0 cap 4, e1 cap 10, both cost 1), N = 2, D = 8 (requests 4 and 4). DETERMINISTIC
// ignores priorities: e0 (lowest EdgeId), then e1 (e0 is filtered at U = 0). PREFER_HIGHER_PRIORITY:
// CAPACITY reads C (10 > 4): e1 twice; RESIDUAL reads U (10 > 4, then 6 > 4): e1 twice; HEADROOM
// reads U/C (4/4 = 10/10: a tie, then the lower EdgeId e0; then e0 is filtered): e0, e1.
TEST(T1Setup, PriorityTies) {
  const auto g = ordered(2, {{0, 1, 4, 1}, {0, 1, 10, 1}});
  struct Case {
    PathSelection s;
    std::vector<Edges> routes;
    Amounts holds;
  };
  for (const auto& c : {Case {sel(DET, EdgePriority::Capacity), {{0}, {1}}, {4, 4}},
                        Case {sel(PHP, EdgePriority::Capacity), {{1}, {1}}, {0, 8}},
                        Case {sel(PHP, EdgePriority::Residual), {{1}, {1}}, {0, 8}},
                        Case {sel(PHP, EdgePriority::Headroom), {{0}, {1}}, {4, 4}}}) {
    auto p = make(g, {0}, {1}, t1(2, LOSSLESS, std::nullopt, 0, c.s));
    const auto r = p->replace(8, 8);
    EXPECT_EQ(routes(*r.program), c.routes);
    EXPECT_EQ(p.ledger->reserved(), c.holds);
  }
}

// PVIEW (e0 cap 10, e1 cap 8) with FlowLedger(g, unreserved = [5, 8]), N = 1, D = 4:
// DETERMINISTIC e0; PHP CAPACITY e0 (C 10 > 8); PHP RESIDUAL e1 (U 8 > 5); PHP HEADROOM e1 (8/8 >
// 5/10). RESIDUAL and HEADROOM read U, not the service residual (both 10 and 8 here).
TEST(T1Setup, PriorityReadsUnreserved) {
  const auto g = ordered(2, {{0, 1, 10, 1}, {0, 1, 8, 1}});
  const Amounts unreserved = {5, 8};
  for (const auto& [s, want] : {std::pair {sel(DET, EdgePriority::Capacity), Edges {0}},
                                std::pair {sel(PHP, EdgePriority::Capacity), Edges {0}},
                                std::pair {sel(PHP, EdgePriority::Residual), Edges {1}},
                                std::pair {sel(PHP, EdgePriority::Headroom), Edges {1}}}) {
    auto ledger = std::make_shared<FlowLedger>(g, g->capacity_units(), unreserved);
    auto p = make(g, {0}, {1}, t1(1, LOSSLESS, std::nullopt, 0, s), nullptr, ledger);
    EXPECT_EQ(routes(*p->replace(4, 4).program), (std::vector<Edges> {want}));
  }
}

// FILT (e0 cap 2, e1 cap 10), N = 1, D = 6: the filter (e0: 2 < 6) runs before the
// collapse of parallel members, which would otherwise keep e0 (the lowest EdgeId at equal cost).
TEST(T1Setup, FilterBeforeCollapse) {
  const auto g = ordered(2, {{0, 1, 2, 1}, {0, 1, 10, 1}});
  for (const bool multi_edge : {false, true}) {
    auto p = make(g, {0}, {1}, t1(1, LOSSLESS, std::nullopt, 0, sel(DET, EdgePriority::Capacity, multi_edge)));
    const auto r = p->replace(6, 6);
    EXPECT_EQ(routes(*r.program), (std::vector<Edges> {{1}})) << multi_edge;
    EXPECT_EQ(statuses(*r.program), (std::vector<TunnelStatus> {UP})) << multi_edge;
  }
  // The last-tunnel-without-bandwidth example with multi_edge = false (the filter comes before the
  // collapse): UP, UP, UP, NB as with true.
  const auto ax3 = ordered(2, {{0, 1, 2, 1}, {0, 1, 10, 2}});
  auto p = make(ax3, {0}, {1}, t1(4, LOSSLESS, std::nullopt, 0, sel(DET, EdgePriority::Capacity, false)));
  EXPECT_EQ(statuses(*p->replace(12, 12).program), (std::vector<TunnelStatus> {UP, UP, UP, NB}));
}

// min_unreserved_bandwidth filters every attempt at max(b, minimum) and is never a
// reservation. PARC (e0 cap 4 cost 1, e1 cap 10 cost 2), N = 2, D = 2: minimum 5 filters e0 (4 <
// 5): both on e1, holds [0, 2], bound 10; minimum 0: both on e0, holds [2, 0], bound 4. ZL with
// minimum 1, D = 1, Lossy O = 2: requests 1 and 0, the request-0 tunnel now filtered at 1 too:
// both on e1, holds [0, 1], delivered 2 (without the minimum, 1 is delivered).
TEST(T1Setup, MinUnreserved) {
  const auto parc = ordered(2, {{0, 1, 4, 1}, {0, 1, 10, 2}});
  for (const auto& [minimum, path, holding, bound] :
       {std::tuple {Amount {5}, Edges {1}, Amounts {0, 2}, Amount {10}}, std::tuple {Amount {0}, Edges {0}, Amounts {2, 0}, Amount {4}}}) {
    auto p = make(parc, {0}, {1}, t1(2, LOSSLESS, std::nullopt, minimum));
    const auto r = p->replace(INF, 2);
    EXPECT_EQ(routes(*r.program), (std::vector<Edges> {path, path})) << minimum;
    EXPECT_EQ(p.ledger->reserved(), holding) << minimum;
    EXPECT_EQ(r.placed, bound) << minimum;
    EXPECT_EQ(r.plan->reserved_total, 2) << minimum;  // the requests, not the minimum
  }
  auto z = make(zl(), {0}, {1}, t1(2, LOSSY, std::nullopt, 1));
  const auto r = z->replace(2, 1);
  EXPECT_EQ(requests(*r.program), (Amounts {1, 0}));
  EXPECT_EQ(routes(*r.program), (std::vector<Edges> {{1}, {1}}));
  EXPECT_EQ(z.ledger->reserved(), (Amounts {0, 1}));
  EXPECT_EQ(r.placed, 2);
}

// Other owners' holds reduce U (10.1). ONE10: policy A (N = 1) replace(6, 6) holds 6
// (reserved [6], unreserved [4]); policy B (N = 1, Lossy) replace(6, 6) sees U = 4 < 6: NB, its
// share lost at the head; after A.remove() B's replace(6, 6) is up with the hold [6].
TEST(T1Setup, OtherOwnersHolds) {
  const auto g = one10();
  auto ledger = std::make_shared<FlowLedger>(g);
  auto a = make(g, {0}, {1}, t1(1), nullptr, ledger);
  auto b = make(g, {0}, {1}, t1(1, LOSSY), nullptr, ledger);
  EXPECT_EQ(statuses(*a->replace(6, 6).program), (std::vector<TunnelStatus> {UP}));
  EXPECT_EQ(ledger->reserved(), (Amounts {6}));
  EXPECT_EQ(ledger->unreserved(), (Amounts {4}));
  const auto down = b->replace(6, 6);
  EXPECT_EQ(statuses(*down.program), (std::vector<TunnelStatus> {NB}));
  EXPECT_EQ(down.call_source_loss, 6);
  a->remove();
  const auto up = b->replace(6, 6);
  EXPECT_EQ(statuses(*up.program), (std::vector<TunnelStatus> {UP}));
  EXPECT_EQ(ledger->reserved(), (Amounts {6}));
  expect_balanced(*ledger, "after B's setup");
}

// TWOH(10, 10), N = 2, D = 10 with node 0 masked: the roster is [1] (requests 5, 5,
// head_design ((1, 10))); source 0 keeps a zero by_source row (Lossy: offered 0).
TEST(T1Setup, MaskedSource) {
  const auto g = twoh(10, 10);
  const auto mask = g->eligibility(std::vector<NodeId> {0}, {});
  auto p = make(g, {0, 1}, {2}, t1(2, LOSSY), &mask);
  const auto r = p->replace(10, 10);
  EXPECT_EQ(requests(*r.program), (Amounts {5, 5}));
  EXPECT_EQ(r.plan->head_design, (std::vector<std::pair<NodeId, Amount>> {{1, 10}}));
  EXPECT_EQ(source_rows(r), (std::vector<SourceRow> {{0, 0, 0, 0, 0, 0}, {1, 10, 10, 10, 0, 0}}));
  EXPECT_EQ(vec(r.program->sources()), v({0, 1}));
}

// Stretch(2) on X8: c0 = 1 is the pool minimum, through head 0's
// zero-capacity link (before any filter): head 1's only path costs 3 > 2: COST_CEILING (T0 took a
// per-head baseline of 3 and admitted it). Stretch(1e300) does not overflow: NB, UP.
TEST(T1Setup, StretchPoolBaseline) {
  for (const auto& [f, second] : {std::pair {2.0, CC}, std::pair {1e300, UP}}) {
    auto p = make(x8(), {0, 1}, {2}, t1(1, LOSSY, Stretch {f}));
    EXPECT_EQ(statuses(*p->replace(4, 4).program), (std::vector<TunnelStatus> {NB, second})) << f;
  }
  // Absolute is inclusive: Absolute(3) admits the cost-3 path, Absolute(2) not.
  auto in = make(x8(), {0, 1}, {2}, t1(1, LOSSY, Absolute {3}));
  EXPECT_EQ(statuses(*in->replace(4, 4).program), (std::vector<TunnelStatus> {NB, UP}));
  auto out = make(x8(), {0, 1}, {2}, t1(1, LOSSY, Absolute {2}));
  EXPECT_EQ(statuses(*out->replace(4, 4).program), (std::vector<TunnelStatus> {NB, CC}));
  // Large costs compare exactly: two parallel paths of 2^60 and 2^60 + 1 (two units each), D = 4,
  // N = 2: Stretch(1) admits only the cheaper one for the first tunnel, the second tunnel's only
  // remaining path costs 2^60 + 1: CostCeiling; Stretch(1 + 2^-52) admits it.
  const Cost half = Cost {1} << 59;
  const auto big = ordered(5, {{0, 1, 2, half}, {0, 2, 2, half}, {1, 4, 2, half}, {2, 4, 2, half + 1}});
  for (const auto& [f, want] : {std::pair {1.0, std::vector<TunnelStatus> {UP, CC}},
                                std::pair {std::nextafter(1.0, 2.0), std::vector<TunnelStatus> {UP, UP}},
                                std::pair {std::ldexp(1.0, 70), std::vector<TunnelStatus> {UP, UP}}}) {
    auto p = make(big, {0}, {4}, t1(2, LOSSLESS, Stretch {f}));
    EXPECT_EQ(statuses(*p->replace(4, 4).program), want) << f;
  }
  // The baseline ignores the setup view: a third party's holds on the cheap path leave c0 = 2^60.
  auto ledger = std::make_shared<FlowLedger>(big);
  auto third = make(big, {0}, {4}, t1(1), nullptr, ledger);
  EXPECT_EQ(statuses(*third->replace(2, 2).program), (std::vector<TunnelStatus> {UP}));
  auto p = make(big, {0}, {4}, t1(1, LOSSLESS, Stretch {1.0}), nullptr, ledger);
  EXPECT_EQ(statuses(*p->replace(2, 2).program), (std::vector<TunnelStatus> {CC}));
  // A fractional factor is inclusive at equality: Stretch(1.5) on c0 = 2 admits a path of cost 3
  // exactly and rejects one of cost 4 (the cheaper links are too small for the request of 5).
  auto at = make(ordered(2, {{0, 1, 1, 2}, {0, 1, 10, 3}}), {0}, {1}, t1(1, LOSSLESS, Stretch {1.5}));
  const auto ra = at->replace(5, 5);
  EXPECT_EQ(statuses(*ra.program), (std::vector<TunnelStatus> {UP}));
  EXPECT_EQ(ra.program->tunnels()[0].route_cost, std::optional<Cost> {3});
  auto above = make(ordered(2, {{0, 1, 1, 2}, {0, 1, 1, 3}, {0, 1, 10, 4}}), {0}, {1}, t1(1, LOSSLESS, Stretch {1.5}));
  EXPECT_EQ(statuses(*above->replace(5, 5).program), (std::vector<TunnelStatus> {CC}));
}

// ---------------------------------------------------------------- forwarding (6.5)

// TWOH(2, 10) Lossless replace(12, 12): head 0 has no up tunnel, so only the total
// 0 gives it nothing: placed 0, unadmitted 12, Exhausted; the holds stand ([0, 6]).
TEST(T1Forwarding, LosslessDeadHead) {
  auto p = make(twoh(2, 10), {0, 1}, {2}, t1(2));
  const auto r = p->replace(12, 12);
  EXPECT_EQ(r.placed, 0);
  EXPECT_EQ(lo(r.unadmitted), 12);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  EXPECT_EQ(p.ledger->reserved(), (Amounts {0, 6}));
}

// Holds do not police and have no priority. ONE10: a CostOnly Lossy policy
// carries 10 (residual 0); T1 (N = 1, Lossy) replace(4, 4) sets up on unreserved 10 (not the
// residual): up, holding 4; its forwarding meets the residual 0: 4 dropped on e0, nothing delivered.
TEST(T1Forwarding, HoldsDoNotQueue) {
  const auto g = one10();
  auto ledger = std::make_shared<FlowLedger>(g);
  auto c = make(g, {0}, {1}, FlowPolicyConfig {CostOnly {}, FlowPlacement::EqualBalanced, LOSSY}, nullptr, ledger);
  EXPECT_EQ(c->replace(10).placed, 10);
  auto p = make(g, {0}, {1}, t1(1, LOSSY), nullptr, ledger);
  const auto r = p->replace(4, 4);
  EXPECT_EQ(statuses(*r.program), (std::vector<TunnelStatus> {UP}));
  EXPECT_EQ(ledger->reserved(), (Amounts {4}));
  EXPECT_EQ(r.call_drops, (Pairs {{0, 4}}));
  EXPECT_EQ(r.placed, 0);
  EXPECT_EQ(ledger->residual(), (Amounts {0}));
  // Carried above the holds (2.5): each tunnel holds 4 and carries 5.
  auto o = make(ordered(2, {{0, 1, 5, 1}, {0, 1, 5, 1}}), {0}, {1}, t1(2, LOSSY));
  const auto over = o->replace(10, 8);
  EXPECT_EQ(over.placed, 10);
  for (const auto& rec : o->records()) {
    EXPECT_EQ(rec.delivered, 5);
    ASSERT_EQ(rec.reserved_edges.size(), 1U);
    EXPECT_EQ(rec.reserved_edges[0].second, 4);
  }
}

// DIAM5 with D = 14: both tunnels down; Lossy replace(5, 14) loses
// the 5 at the head (Complete: 14.2 rule (3) needs an entry with a path); Lossless places 0
// (Exhausted). The contract's offer 14 loses 14.
TEST(T1Forwarding, AllDown) {
  for (const Amount offer : {5, 14}) {
    auto ly = make(diam5(), {0}, {3}, t1(2, LOSSY));
    const auto r = ly->replace(offer, 14);
    EXPECT_EQ(statuses(*r.program), (std::vector<TunnelStatus> {NB, NB}));
    EXPECT_EQ(r.call_source_loss, offer);
    EXPECT_EQ(r.placed, 0);
    EXPECT_EQ(r.termination, Termination::Complete);
    auto ll = make(diam5(), {0}, {3}, t1(2));
    const auto s = ll->replace(offer, 14);
    EXPECT_EQ(s.placed, 0);
    EXPECT_EQ(s.termination, Termination::Exhausted);
  }
}

// PAR612, N = 4, Lossy replace(13, 12): shares RR(13, 4) = (4, 3, 3, 3); e0
// receives 7 > 6 and carries H(6, [4, 3]) = (4, 2): tunnel 1 drops 1.
TEST(T1Forwarding, ByTunnel) {
  auto p = make(par612(), {0}, {1}, t1(4, LOSSY));
  const auto r = p->replace(13, 12);
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 4, 4, 0}, {1, 0, UP, 3, 2, 1}, {2, 0, UP, 3, 3, 0}, {3, 0, UP, 3, 3, 0}}));
  EXPECT_EQ(r.call_drops, (Pairs {{0, 1}}));
  EXPECT_EQ(r.placed, 12);
}

// ---------------------------------------------------------------- operations (6.6, 10.2, 10.3)

// The contract's call of 6.6: TWOH(2, 10), N = 2, Lossy replace(12, 12),
// then forward(8).
TEST(T1Operations, ContractCall) {
  auto p = make(twoh(2, 10), {0, 1}, {2}, t1(2, LOSSY));
  const auto r = p->replace(12, 12, kHolds);
  EXPECT_EQ(statuses(*r.program), (std::vector<TunnelStatus> {NB, UP, NB, UP}));
  EXPECT_EQ(r.placed, 6);
  EXPECT_EQ(r.call_source_loss, 6);
  EXPECT_EQ(r.plan->reserved_total, 6);
  EXPECT_EQ(r.plan->count, 2);
  EXPECT_EQ(r.plan->design_demand, std::optional<Amount> {12});
  EXPECT_EQ(r.plan->head_design, (std::vector<std::pair<NodeId, Amount>> {{0, 6}, {1, 6}}));
  EXPECT_EQ(*r.reservation_edge_flows, (Amounts {0, 6}));
  EXPECT_EQ(rows(r), (Rows {{0, 0, NB, 0, 0, 0}, {1, 1, UP, 3, 3, 0}, {2, 0, NB, 0, 0, 0}, {3, 1, UP, 3, 3, 0}}));
  EXPECT_EQ(source_rows(r), (std::vector<SourceRow> {{0, 6, 6, 0, 6, 0}, {1, 6, 6, 6, 0, 0}}));
  const auto recs = p->records();
  EXPECT_EQ(entries(recs), (std::vector<std::optional<std::int32_t>> {0, 1, 2, 3}));
  EXPECT_EQ(holds(recs), (std::vector<Pairs> {{}, {{1, 3}}, {}, {{1, 3}}}));
  for (const std::size_t down : {0U, 2U}) {
    EXPECT_TRUE(ascending_edges(*recs[down].bundle).empty());
    EXPECT_EQ(recs[down].route_cost, kUnreachable);
  }
  EXPECT_EQ(ascending_edges(*recs[1].bundle), (Edges {1}));
  EXPECT_EQ(recs[1].route_cost, 1);
  std::vector<std::int64_t> before;
  for (const auto& rec : recs) before.push_back(rec.flow_id.second);
  // forward(8): the same tunnels and holds; head 0 still loses its share (4).
  const auto f = p->forward(8);
  EXPECT_EQ(f.placed, 4);
  EXPECT_EQ(f.call_source_loss, 4);
  EXPECT_EQ(f.delivery_delta, -2);
  EXPECT_FALSE(f.routes_changed || f.weights_changed || f.reservations_changed);
  EXPECT_EQ(f.program, r.program);
  EXPECT_TRUE(f.plan == r.plan);
  EXPECT_EQ(p.ledger->reserved(), (Amounts {0, 6}));
  Amounts shares;
  for (const auto& t : f.by_tunnel) shares.push_back(lo(t.assigned));
  EXPECT_EQ(shares, (Amounts {0, 2, 0, 2}));
  const auto after = p->records();
  EXPECT_EQ(holds(after), holds(recs));
  for (std::size_t i = 0; i < after.size(); ++i) EXPECT_NE(after[i].flow_id.second, before[i]);
  expect_balanced(*p.ledger, "after forward");
}

// ONE10, N = 1, Lossless replace(6, 6) twice: the second sees U = 10 (its own hold
// released): up again, reserved [6] (not 12); all three flags false; the record id is new.
TEST(T1Operations, ReplaceReleasesOwnHolds) {
  auto p = make(one10(), {0}, {1}, t1(1));
  (void)p->replace(6, 6);
  const auto id = p->records().front().flow_id;
  const auto r = p->replace(6, 6);
  EXPECT_EQ(statuses(*r.program), (std::vector<TunnelStatus> {UP}));
  EXPECT_EQ(p.ledger->reserved(), (Amounts {6}));
  EXPECT_FALSE(r.routes_changed || r.weights_changed || r.reservations_changed);
  EXPECT_NE(p->records().front().flow_id, id);
}

// After replace(6, 6) on ONE10, evaluate(3, 8) sets up for 8 (its result holds 8) and publishes
// nothing; replace(0) with N = 2 (D = 0): both up without a hold, Complete; replace(0, 12): requests
// 6, 6, tunnel 0 up (U = 10), tunnel 1 NB (U = 4 < 6), holds [6]: a zero offer still books holds
// (10.3).
TEST(T1Operations, EvaluateAndZero) {
  auto p = make(one10(), {0}, {1}, t1(1));
  (void)p->replace(6, 6);
  const auto ledger = ledger_state(*p.ledger);
  const auto ep = epoch(*p.p);
  const auto e = p->evaluate(3, 8, kHolds);
  EXPECT_EQ(statuses(*e.program), (std::vector<TunnelStatus> {UP}));
  EXPECT_EQ(requests(*e.program), (Amounts {8}));
  EXPECT_EQ(e.program->tunnels()[0].reserved, 8);
  EXPECT_EQ(e.plan->reserved_total, 8);
  EXPECT_EQ(*e.reservation_edge_flows, (Amounts {8}));
  EXPECT_EQ(e.placed, 3);
  expect_same(ledger_state(*p.ledger), ledger, "evaluate");
  expect_unchanged(epoch(*p.p), ep, "evaluate");
  auto q = make(one10(), {0}, {1}, t1(2));
  const auto z = q->replace(0);
  EXPECT_EQ(statuses(*z.program), (std::vector<TunnelStatus> {UP, UP}));
  EXPECT_EQ(requests(*z.program), (Amounts {0, 0}));
  EXPECT_EQ(q.ledger->reserved(), (Amounts {0}));
  EXPECT_EQ(z.termination, Termination::Complete);
  EXPECT_EQ(z.placed, 0);
  const auto z12 = q->replace(0, 12);
  EXPECT_EQ(requests(*z12.program), (Amounts {6, 6}));
  EXPECT_EQ(statuses(*z12.program), (std::vector<TunnelStatus> {UP, NB}));
  EXPECT_EQ(q.ledger->reserved(), (Amounts {6}));
  EXPECT_EQ(z12.termination, Termination::Complete);
}

// PAR612, N = 4, Lossless replace(+inf, 12): placed 12, remaining +inf, requested +inf,
// unadmitted absent, Exhausted.
TEST(T1Operations, InfiniteOffer) {
  auto p = make(par612(), {0}, {1}, t1(4));
  const auto r = p->replace(INF, 12);
  EXPECT_EQ(r.placed, 12);
  EXPECT_FALSE(r.remaining.has_value());
  EXPECT_FALSE(r.requested.has_value());
  EXPECT_FALSE(r.unadmitted.has_value());
  EXPECT_EQ(r.termination, Termination::Exhausted);
}

// EMPTY (node 0 has no edge), N = 2: replace(5, 5) Lossy is NO_ROUTE with the empty
// program (sources (0,), destinations (2,), no tunnel), the plan (count 2, design 5, no head, no
// hold), the 5 lost outside every row; replace(0, 5), evaluate(0, 5) and forward(0) are Complete.
TEST(T1Operations, EmptyRoster) {
  const auto g = ordered(3, {{1, 2, 5, 1}});
  auto p = make(g, {0}, {2}, t1(2, LOSSY));
  const auto r = p->replace(5, 5);
  EXPECT_EQ(r.termination, Termination::NoRoute);
  ASSERT_TRUE(r.program);
  EXPECT_EQ(vec(r.program->sources()), v({0}));
  EXPECT_EQ(vec(r.program->destinations()), v({2}));
  EXPECT_TRUE(r.program->tunnels().empty());
  EXPECT_TRUE(r.plan == (TunnelPlan {2, 5, {}, 0}));
  EXPECT_EQ(r.call_source_loss, 5);
  EXPECT_EQ(source_rows(r), (std::vector<SourceRow> {{0, 0, 0, 0, 0, 0}}));
  EXPECT_TRUE(p->records().empty());
  for (const auto& z : {p->replace(0, 5), p->evaluate(0, 5), p->forward(0)}) {
    EXPECT_EQ(z.termination, Termination::Complete);
    EXPECT_TRUE(z.program && z.program->tunnels().empty());
    EXPECT_TRUE(z.plan.has_value());
  }
  EXPECT_TRUE(p->records().empty());
}

// DIAM5, N = 2: replace(8, 8) (one tunnel up), replace(10, 10) (both up: routes,
// weights 1, 0 to 1, 1, and reservations all change), forward(10) (none). The first program of
// the diamond with design 14 (all down) changes routes and reservations (requests 7, 7) but no weight.
TEST(T1Operations, Flags) {
  auto p = make(diam5(), {0}, {3}, t1(2));
  const auto first = p->replace(8, 8);
  EXPECT_TRUE(first.routes_changed && first.weights_changed && first.reservations_changed);
  const auto second = p->replace(10, 10);
  EXPECT_TRUE(second.routes_changed && second.weights_changed && second.reservations_changed);
  const auto f = p->forward(10);
  EXPECT_FALSE(f.routes_changed || f.weights_changed || f.reservations_changed);
  auto down = make(diam5(), {0}, {3}, t1(2, LOSSY));
  const auto d = down->replace(5, 14);
  EXPECT_TRUE(d.routes_changed);
  EXPECT_FALSE(d.weights_changed);
  EXPECT_TRUE(d.reservations_changed);
  // Only the requests change (D = 15: requests 8 and 7, still all down): reservations alone.
  const auto d2 = down->replace(5, 15);
  EXPECT_FALSE(d2.routes_changed || d2.weights_changed);
  EXPECT_TRUE(d2.reservations_changed);
  // A controller's tunnels hold nothing: a change in their number is a route change, never a
  // reservation change (a tunnel absent from one program counts as (0, 0) there).
  auto gl = make(diam5(), {0}, {3}, global_config());
  const auto many = gl->replace(INF);
  ASSERT_GT(many.program->tunnels().size(), 1u);
  EXPECT_FALSE(many.reservations_changed);
  const auto one = gl->replace(1);
  ASSERT_EQ(one.program->tunnels().size(), 1u);
  EXPECT_TRUE(one.routes_changed);
  EXPECT_FALSE(one.reservations_changed);
  const auto again = gl->replace(INF);
  EXPECT_TRUE(again.routes_changed);
  EXPECT_FALSE(again.reservations_changed);
  // Equal requests, changed holds: two links of 4, N = 2, D = 8 (requests 4, 4). Another owner's
  // hold of 1 on e0 leaves tunnel 0 on e1 and tunnel 1 down (reserved 4, 0); once it leaves, the
  // same design brings both up (reserved 4, 4): the reservations changed though no request did.
  const auto two = ordered(2, {{0, 1, 4, 1}, {0, 1, 4, 1}});
  auto ledger = std::make_shared<FlowLedger>(two);
  auto other = make(two, {0}, {1}, t1(1), nullptr, ledger);
  (void)other->replace(1, 1);
  auto mine = make(two, {0}, {1}, t1(2), nullptr, ledger);
  const auto held = mine->replace(8, 8);
  EXPECT_EQ(statuses(*held.program), (std::vector<TunnelStatus> {UP, NB}));
  other->remove();
  const auto freed = mine->replace(8, 8);
  EXPECT_EQ(statuses(*freed.program), (std::vector<TunnelStatus> {UP, UP}));
  EXPECT_EQ(requests(*held.program), requests(*freed.program));
  EXPECT_TRUE(freed.reservations_changed);
}

// The comparison behind the flags (10.4): a tunnel absent from one program counts as (0, 0) there,
// so programs of zero requests and holds with different tunnel counts (or none) are no reservation
// change, and a program that holds something against none is one. No T1 call changes its tunnel
// count, so only this direct check reaches the absent-tunnel branches.
TEST(T1Operations, AbsentTunnelsCountAsZero) {
  auto two = make(one10(), {0}, {1}, t1(2));
  const auto zero2 = two->replace(0, 0).program;  // two tunnels, requests 0
  auto one = make(one10(), {0}, {1}, t1(1));
  const auto zero1 = one->replace(0, 0).program;  // one tunnel, request 0
  const auto held1 = one->replace(0, 6).program;  // one tunnel holding 6
  ASSERT_EQ(zero2->tunnels().size(), 2U);
  ASSERT_EQ(held1->tunnels()[0].reserved, 6);
  EXPECT_FALSE(detail::compare_programs(zero2.get(), zero1.get()).reservations);  // absent after
  EXPECT_FALSE(detail::compare_programs(zero1.get(), zero2.get()).reservations);  // absent before
  EXPECT_FALSE(detail::compare_programs(zero2.get(), nullptr).reservations);
  EXPECT_FALSE(detail::compare_programs(nullptr, zero2.get()).reservations);
  EXPECT_TRUE(detail::compare_programs(held1.get(), nullptr).reservations);
  EXPECT_TRUE(detail::compare_programs(nullptr, held1.get()).reservations);
  EXPECT_TRUE(detail::compare_programs(held1.get(), zero2.get()).reservations);
}

// T1 has no incremental place (also 0), and forward needs an installed program.
TEST(T1Operations, NoIncrementsAndForwardNeedsAProgram) {
  auto p = make(one10(), {0}, {1}, t1(1));
  EXPECT_EQ(error_of([&] { (void)p->place(5); }), "T1 does not support incremental placement");
  EXPECT_EQ(error_of([&] { (void)p->place(0); }), "T1 does not support incremental placement");
  EXPECT_EQ(error_of([&] { (void)p->forward(5); }), "forward requires an installed program");
  (void)p->replace(1, 1);
  p->remove();
  EXPECT_EQ(error_of([&] { (void)p->forward(5); }), "forward requires an installed program");
}

// The call's checks in their order (the offer, LOSSY's finite offer, then the
// design: the model, its domain, finiteness) and the terminations.
TEST(T1Operations, DesignValidationAndTerminations) {
  const auto g = one10();
  for (const auto& c : {cost_only_config(), global_config(),
                        FlowPolicyConfig {.routing = pinned_paths(*g, {0}, {1}, {{0, {0}, 1}})}}) {
    auto p = make(g, {0}, {1}, c);
    const std::string model = detail::model_name(c.routing);
    EXPECT_EQ(error_of([&] { (void)p->replace(5, 5); }), model + " does not accept design_demand");
    EXPECT_EQ(error_of([&] { (void)p->evaluate(5, 5); }), model + " does not accept design_demand");
    EXPECT_EQ(error_of([&] { (void)p->replace(5, std::optional<Amount> {}); }), model + " does not accept design_demand");
  }
  auto ll = make(g, {0}, {1}, t1(1));
  const std::string domain = "design_demand must be an integer in [0, 2^62) or permitted +inf";
  EXPECT_EQ(error_of([&] { (void)ll->replace(5, -1); }), domain);
  EXPECT_EQ(error_of([&] { (void)ll->replace(5, kL62); }), domain);
  EXPECT_EQ(error_of([&] { (void)ll->replace(-1, -1); }), "offer must be an integer in [0, 2^62) or permitted +inf");
  EXPECT_EQ(error_of([&] { (void)ll->replace(5, std::optional<Amount> {}); }), "T1 requires a finite design demand");
  EXPECT_EQ(error_of([&] { (void)ll->replace(INF); }), "T1 requires a finite design demand");
  EXPECT_EQ(error_of([&] { (void)ll->evaluate(INF); }), "T1 requires a finite design demand");
  auto ly = make(g, {0}, {1}, t1(1, LOSSY));
  EXPECT_EQ(error_of([&] { (void)ly->replace(INF, 5); }), "LOSSY forwarding requires a finite offer");
  EXPECT_TRUE(ll.ledger->records().empty() && ly.ledger->records().empty());
  // Terminations: zero offer, Lossy finite, Lossless full and partial, unlimited.
  EXPECT_EQ(ll->replace(0, 5).termination, Termination::Complete);
  EXPECT_EQ(ll->replace(10, 10).termination, Termination::Complete);
  EXPECT_EQ(ll->replace(11, 10).termination, Termination::Exhausted);
  EXPECT_EQ(ll->replace(INF, 10).termination, Termination::Exhausted);
  EXPECT_EQ(ly->replace(11, 10).termination, Termination::Complete);
  // A zero-weight route, the case rule (3) keeps: one supplied path of weight 0, both admissions.
  for (const auto a : {LOSSLESS, LOSSY}) {
    auto z = make(g, {0}, {1}, FlowPolicyConfig {.routing = pinned_paths(*g, {0}, {1}, {{0, {0}, 0}}), .admission = a});
    const auto r = z->replace(4);
    EXPECT_EQ(r.termination, Termination::NoActiveWeight);
    if (a == LOSSY) EXPECT_EQ(r.call_source_loss, 4);
  }
}

// No certificate for T1 (total_cost, min_cut, cost_distribution absent), and plan
// only for T1 (absent for the other models, Pinned over a T1 program too).
TEST(T1Operations, NoCertificatesAndPlanOnlyForT1) {
  const auto g = twoh(2, 10);
  auto p = make(g, {0, 1}, {2}, t1(2));
  const auto r = p->replace(INF, 12, OutputOptions {.min_cut = true});
  EXPECT_FALSE(r.total_cost.has_value() || r.min_cut.has_value() || r.cost_distribution.has_value());
  ASSERT_TRUE(r.plan.has_value());
  auto pin = make(g, {0, 1}, {2}, FlowPolicyConfig {.routing = Pinned {r.program, nullptr}});
  EXPECT_FALSE(pin->replace(12).plan.has_value());
  for (const auto& c : {cost_only_config(), global_config()}) {
    auto q = make(g, {0, 1}, {2}, c);
    EXPECT_FALSE(q->place(5).plan.has_value());
    EXPECT_FALSE(q->replace(5).plan.has_value());
  }
}

// ---------------------------------------------------------------- the two projections (10.1, 10.3)

// The two projections of one link: ONE10, N = 1, Lossless replace(4, 6): carried 4, held 6 (allocated [4],
// residual [6], reserved [6], unreserved [4]); then replace(2, 6): placed 2, delivery_delta -2. A
// CostOnly policy on the same ledger then places 6 (the residual): holds are never added to
// carried load.
TEST(LedgerProjections, ELedger) {
  const auto g = one10();
  auto ledger = std::make_shared<FlowLedger>(g);
  auto p = make(g, {0}, {1}, t1(1), nullptr, ledger);
  (void)p->replace(4, 6);
  EXPECT_EQ(ledger->allocated(), (Amounts {4}));
  EXPECT_EQ(ledger->residual(), (Amounts {6}));
  EXPECT_EQ(ledger->reserved(), (Amounts {6}));
  EXPECT_EQ(ledger->unreserved(), (Amounts {4}));
  auto c = make(g, {0}, {1}, cost_only_config(), nullptr, ledger);
  EXPECT_EQ(c->place(10).placed, 6);
  c->remove();
  const auto r = p->replace(2, 6);
  EXPECT_EQ(r.placed, 2);
  EXPECT_EQ(r.delivery_delta, -2);
  expect_balanced(*ledger, "the two projections of one link");
}

// The reservation baseline is the setup's U: unreserved [3] puts
// a request of 4 down while the residual stays 10; its texts are the residual baseline's with the
// name unreserved, checked after the residual's.
TEST(LedgerProjections, Baselines) {
  const auto g = one10();
  auto ledger = std::make_shared<FlowLedger>(g, g->capacity_units(), Amounts {3});
  auto p = make(g, {0}, {1}, t1(1), nullptr, ledger);
  EXPECT_EQ(statuses(*p->replace(4, 4).program), (std::vector<TunnelStatus> {NB}));
  EXPECT_EQ(ledger->residual(), (Amounts {10}));
  EXPECT_EQ(ledger->unreserved(), (Amounts {3}));
  EXPECT_EQ(error_of([&] { FlowLedger x(g, g->capacity_units(), Amounts {3, 1}); }),
            "unreserved has length 2, expected 1 (one entry per edge)");
  EXPECT_EQ(error_of([&] { FlowLedger x(g, g->capacity_units(), Amounts {11}); }),
            "unreserved[0] = 11 is outside [0, capacity_units[0] = 10]");
  EXPECT_EQ(error_of([&] { FlowLedger x(g, Amounts {11}, Amounts {11}); }),
            "residual[0] = 11 is outside [0, capacity_units[0] = 10]");
  // The contract's two-edge cases of the unreserved texts.
  const auto two = ordered(2, {{0, 1, 10, 1}, {0, 1, 10, 1}});
  EXPECT_EQ(error_of([&] { FlowLedger x(two, two->capacity_units(), Amounts {1, 2, 3}); }),
            "unreserved has length 3, expected 2 (one entry per edge)");
  EXPECT_EQ(error_of([&] { FlowLedger x(two, two->capacity_units(), Amounts {5, 11}); }),
            "unreserved[1] = 11 is outside [0, capacity_units[1] = 10]");
  EXPECT_EQ(error_of([&] { FlowLedger x(two, two->capacity_units(), Amounts {10, 0}); }), "");
}

// remove() releases both projections and clears the program.
TEST(LedgerProjections, RemoveReleasesBoth) {
  auto p = make(one10(), {0}, {1}, t1(1));
  (void)p->replace(4, 6);
  p->remove();
  EXPECT_EQ(p.ledger->reserved(), (Amounts {0}));
  EXPECT_EQ(p.ledger->allocated(), (Amounts {0}));
  EXPECT_EQ(p->program(), nullptr);
  EXPECT_TRUE(p->records().empty());
  p->remove();  // a repeated remove does nothing
  EXPECT_EQ(p.ledger->reserved(), (Amounts {0}));
}

// A policy destroyed without remove() leaves its holds as orphans (records readable)
// until reset(); reset(residual, unreserved) sets both baselines; a copy of an unbound ledger keeps
// the holds and records under a new uid; a bound ledger cannot be copied or reset.
TEST(LedgerProjections, OrphansResetAndCopy) {
  const auto g = one10();
  auto ledger = std::make_shared<FlowLedger>(g);
  {
    auto p = std::make_unique<FlowPolicy>(*ledger, v({0}), v({1}), t1(1));
    (void)p->replace(4, 6);
    EXPECT_EQ(error_of([&] { FlowLedger copy(*ledger); }), "FlowLedger: a ledger with bound policies cannot be copied");
    EXPECT_EQ(error_of([&] { ledger->reset(); }), "FlowLedger.reset is rejected while a policy is bound (1 bound)");
  }
  EXPECT_EQ(ledger->reserved(), (Amounts {6}));
  ASSERT_EQ(ledger->records().size(), 1U);
  EXPECT_EQ(ledger->records().front().reserved_edges, (Pairs {{0, 6}}));
  const FlowLedger copy(*ledger);
  EXPECT_NE(copy.uid(), ledger->uid());
  EXPECT_EQ(copy.reserved(), (Amounts {6}));
  EXPECT_EQ(copy.unreserved(), (Amounts {4}));
  ASSERT_EQ(copy.records().size(), 1U);
  EXPECT_EQ(copy.records().front().reserved_edges, (Pairs {{0, 6}}));
  EXPECT_EQ(copy.records().front().flow_id.first, copy.uid());
  expect_balanced(copy, "the copy");
  ledger->reset();
  EXPECT_EQ(ledger->reserved(), (Amounts {0}));
  EXPECT_TRUE(ledger->records().empty());
  EXPECT_EQ(copy.reserved(), (Amounts {6}));  // the copy owns its holds
  ledger->reset(g->capacity_units(), Amounts {3});
  EXPECT_EQ(ledger->unreserved(), (Amounts {3}));
  EXPECT_EQ(ledger->residual(), (Amounts {10}));
  ledger->reset(Amounts {7});
  EXPECT_EQ(ledger->unreserved(), (Amounts {10}));  // an absent reservation baseline: the capacities
  EXPECT_EQ(ledger->residual(), (Amounts {7}));
}

// A ledger stores no holds until something holds (a policy without holds pays nothing for the
// reservation projection): it reads holds 0 and the capacities, or the given reservation
// baseline, as unreserved; the first hold stores them; a reset empties them again. A copy of an
// unheld ledger owns its own storage once either side holds.
TEST(LedgerProjections, HoldsAreStoredFromTheFirstHold) {
  const auto g = twoh(2, 10);
  FlowLedger fg(g);
  EXPECT_TRUE(detail::LedgerAccess::reserved(fg).empty());
  EXPECT_EQ(fg.reserved(), (Amounts {0, 0}));
  EXPECT_EQ(fg.unreserved(), (Amounts {2, 10}));
  EXPECT_EQ(vec(detail::LedgerAccess::reservation_baseline(fg)), (Amounts {2, 10}));
  {
    FlowPolicy cost(fg, v({0, 1}), v({2}), cost_only_config());
    (void)cost.place(3);
    EXPECT_TRUE(detail::LedgerAccess::reserved(fg).empty());  // carried traffic holds nothing
    cost.remove();
  }
  const FlowLedger unheld(fg);
  {
    FlowPolicy p(fg, v({0, 1}), v({2}), t1(2, LOSSY));
    (void)p.replace(12, 12);
    EXPECT_EQ(fg.reserved(), (Amounts {0, 6}));
    EXPECT_EQ(fg.unreserved(), (Amounts {2, 4}));
    EXPECT_EQ(unheld.reserved(), (Amounts {0, 0}));  // the copy made before the hold has none
    EXPECT_TRUE(detail::LedgerAccess::reserved(unheld).empty());
    p.remove();
    EXPECT_EQ(fg.reserved(), (Amounts {0, 0}));
  }
  fg.reset();
  EXPECT_TRUE(detail::LedgerAccess::reserved(fg).empty());
  EXPECT_EQ(fg.unreserved(), (Amounts {2, 10}));
  FlowLedger given(g, g->capacity_units(), Amounts {1, 7});
  EXPECT_TRUE(detail::LedgerAccess::reserved(given).empty());
  EXPECT_EQ(given.reserved(), (Amounts {0, 0}));
  EXPECT_EQ(given.unreserved(), (Amounts {1, 7}));
  FlowPolicy q(given, v({1}), v({2}), t1(1));
  (void)q.replace(0, 5);
  EXPECT_EQ(given.reserved(), (Amounts {0, 5}));
  EXPECT_EQ(given.unreserved(), (Amounts {1, 2}));
  expect_balanced(given, "the given baseline");
}

// T1's setup reads an absent reservation baseline as the capacities, also when the service
// baseline is smaller: TWOH(2, 10) with residual (1, 4), head 1, N = 1, Lossy replace(6, 6) holds 6
// (U = 10 - 0), and the same after reset(residual) of a ledger that had a given baseline.
TEST(LedgerProjections, AbsentReservationBaselineIsTheCapacities) {
  const auto g = twoh(2, 10);
  for (const bool after_reset : {false, true}) {
    FlowLedger fg = after_reset ? FlowLedger(g, Amounts {1, 4}, Amounts {1, 3}) : FlowLedger(g, Amounts {1, 4});
    if (after_reset) fg.reset(Amounts {1, 4});
    FlowPolicy p(fg, v({1}), v({2}), t1(1, LOSSY));
    const auto r = p.replace(6, 6);
    ASSERT_EQ(r.program->tunnels().size(), 1U);
    EXPECT_EQ(r.program->tunnels()[0].status, UP) << after_reset;
    EXPECT_EQ(fg.reserved(), (Amounts {0, 6})) << after_reset;
    EXPECT_EQ(fg.unreserved(), (Amounts {2, 4})) << after_reset;
    expect_balanced(fg, after_reset ? "after reset" : "two-argument ledger");
    p.remove();
  }
}

// reservation_edge_flows is zeros for a model without holds; an imported T1
// program books no hold; the bottleneck report reads the service view only.
TEST(LedgerProjections, OutputsImportsAndBottleneck) {
  const auto g = twoh(2, 10);
  auto c = make(g, {0, 1}, {2}, cost_only_config());
  EXPECT_EQ(*c->place(4, kHolds).reservation_edge_flows, (Amounts {0, 0}));
  auto p = make(g, {0, 1}, {2}, t1(2, LOSSY));
  const auto r = p->replace(12, 12);
  auto pin = make(g, {0, 1}, {2}, FlowPolicyConfig {.routing = Pinned {r.program, nullptr}, .admission = LOSSY});
  const auto q = pin->replace(12, std::nullopt, kHolds);
  EXPECT_EQ(*q.reservation_edge_flows, (Amounts {0, 0}));
  EXPECT_EQ(pin.ledger->reserved(), (Amounts {0, 0}));
  for (const auto& rec : pin->records()) EXPECT_TRUE(rec.reserved_edges.empty());
  // ONE10 held whole (replace(0, 10): unreserved 0) is not a bottleneck: the residual is 10.
  auto h = make(one10(), {0}, {1}, t1(1));
  (void)h->replace(0, 10);
  EXPECT_EQ(h.ledger->unreserved(), (Amounts {0}));
  EXPECT_TRUE(h.ledger->bottleneck(v({0})).empty());
  EXPECT_EQ(*h->replace(0, 10, OutputOptions {.bottleneck = true}).bottleneck, Edges {});
}

// ---------------------------------------------------------------- programs (8.1, 8.4, 8.6)

// T1's tunnels (the 6.6 call): a down tunnel has no destination, path or cost, weight 0
// and no hold; an up one holds its request (8.1); requested 3 each. The program's roster is the
// policy's (every head), its sources and destinations the policy's, equal sources.
TEST(TunnelPrograms, T1Fields) {
  auto p = make(twoh(2, 10), {0, 1}, {2}, t1(2, LOSSY));
  const auto prog = p->replace(12, 12).program;
  ASSERT_EQ(prog->tunnels().size(), 4U);
  for (std::size_t i = 0; i < 4; ++i) {
    const auto& t = prog->tunnels()[i];
    const bool up = i % 2 == 1;
    EXPECT_EQ(t.index, static_cast<std::int32_t>(i));
    EXPECT_EQ(t.head, static_cast<NodeId>(i % 2));
    EXPECT_EQ(t.local_index, static_cast<std::int32_t>(i / 2));
    EXPECT_EQ(t.status, up ? UP : NB);
    EXPECT_EQ(t.requested, 3);
    EXPECT_EQ(t.reserved, up ? 3 : 0);
    EXPECT_EQ(t.weight, up ? 1 : 0);
    EXPECT_EQ(t.destination, up ? std::optional<NodeId> {2} : std::nullopt);
    EXPECT_EQ(t.edges, up ? std::optional<Edges> {{1}} : std::nullopt);
    EXPECT_EQ(t.route_cost, up ? std::optional<Cost> {1} : std::nullopt);
  }
  EXPECT_EQ(vec(prog->sources()), v({0, 1}));
  EXPECT_EQ(vec(prog->destinations()), v({2}));
  EXPECT_FALSE(prog->source_weights().has_value());
  EXPECT_EQ(vec(detail::ProgramAccess::data(*prog).roster), v({0, 1}));
  // Two-hop tunnels end at the destination, not at the first link's end: DIAM5, N = 2, D = 10
  // (requests 5, 5; e2 of 4 is filtered out): tunnel 0 takes e0 e3, tunnel 1 e1 e4, both cost 4.
  auto d = make(diam5(), {0}, {3}, t1(2));
  const auto two = d->replace(10, 10).program;
  EXPECT_EQ(routes(*two), (std::vector<Edges> {{0, 3}, {1, 4}}));
  for (const auto& t : two->tunnels()) {
    EXPECT_EQ(t.destination, std::optional<NodeId> {3});
    EXPECT_EQ(t.route_cost, std::optional<Cost> {4});
  }
}

// Pinned over an imported T1 program books no holds, so it never reports a reservation change:
// the first evaluate, replace and place over an up program (requested and reserved 6) and over an
// all-down one with a positive request (20 on one link of 10) leave the ledger's and the records'
// holds empty, and the installed program keeps the imported amounts as data.
TEST(TunnelPrograms, PinnedOverT1ReservesNothing) {
  const auto g = one10();
  auto u = make(g, {0}, {1}, t1(1));
  const auto up = u->replace(0, 6).program;
  auto d = make(g, {0}, {1}, t1(1));
  const auto down = d->replace(0, 20).program;
  ASSERT_EQ(statuses(*up), (std::vector<TunnelStatus> {UP}));
  ASSERT_EQ(statuses(*down), (std::vector<TunnelStatus> {NB}));
  EXPECT_EQ(up->tunnels()[0].reserved, 6);
  EXPECT_EQ(down->tunnels()[0].requested, 20);
  for (const auto& prog : {up, down}) {
    for (int op = 0; op < 3; ++op) {
      auto p = make(g, {0}, {1}, FlowPolicyConfig {.routing = Pinned {prog, nullptr}});
      const auto r = op == 0 ? p->evaluate(1) : op == 1 ? p->replace(1) : p->place(1);
      EXPECT_FALSE(r.reservations_changed) << op;
      EXPECT_FALSE(r.plan.has_value()) << op;
      EXPECT_EQ(r.program->tunnels()[0].requested, prog->tunnels()[0].requested) << op;
      EXPECT_EQ(r.program->tunnels()[0].reserved, prog->tunnels()[0].reserved) << op;
      EXPECT_EQ(p.ledger->reserved(), (Amounts {0})) << op;
      for (const auto& rec : p->records()) EXPECT_TRUE(rec.reserved_edges.empty()) << op;
    }
  }
}

// The builder alone (no forwarding runs): an up tunnel has weight 1, a down one weight 0 and no
// path, so a call never forwards a share onto a tunnel without a path.
TEST(TunnelPrograms, BuilderWeighsOnlyUpTunnels) {
  const auto g = twoh(2, 10);
  const std::vector<NodeId> heads {0, 1}, dests {2};
  std::vector<detail::SetupTunnel> setup(2);
  setup[0] = {.slot = 0, .round = 0, .requested = 3, .status = NB};
  setup[1] = {.slot = 1, .round = 0, .requested = 3, .status = UP, .edges = {1}, .cost = 1};
  const auto prog = detail::tunnel_program(*g, heads, dests, heads, std::move(setup), {}, {});
  ASSERT_EQ(prog->tunnels().size(), 2U);
  EXPECT_EQ(prog->tunnels()[0].weight, 0);
  EXPECT_EQ(prog->tunnels()[0].edges, std::nullopt);
  EXPECT_EQ(prog->tunnels()[1].weight, 1);
  EXPECT_EQ(prog->tunnels()[1].edges, std::optional<Edges> {{1}});
}

// Pinned over T1's program on PAR612 (N = 4, replace(13, 12)) with e0 masked: tunnels 0 and 1
// FAILED; the head splits RR(13, 2) = (7, 6) over tunnels 2 and 3 on e1 (capacity 6): H(6, [7, 6])
// = (3, 3) carried, 6 delivered, 7 dropped; Lossless admits 6.
TEST(TunnelPrograms, T1Replay) {
  const auto g = par612();
  auto p = make(g, {0}, {1}, t1(4, LOSSY));
  const auto prog = p->replace(13, 12).program;
  const auto mask = g->eligibility({}, std::vector<EdgeId> {0});
  auto ly = make(g, {0}, {1}, FlowPolicyConfig {.routing = Pinned {prog, nullptr}, .admission = LOSSY}, &mask);
  const auto r = ly->replace(13);
  EXPECT_EQ(rows(r), (Rows {{0, 0, FAILED, 0, 0, 0}, {1, 0, FAILED, 0, 0, 0}, {2, 0, UP, 7, 3, 4}, {3, 0, UP, 6, 3, 3}}));
  EXPECT_EQ(r.placed, 6);
  EXPECT_EQ(r.call_drops, (Pairs {{1, 7}}));
  auto ll = make(g, {0}, {1}, FlowPolicyConfig {.routing = Pinned {prog, nullptr}}, &mask);
  EXPECT_EQ(ll->replace(13).placed, 6);
  // The replay roster comes from the topology basis, not from the surviving tunnels: REPLAY, N = 1,
  // D = 10: head 0's tunnel e0 fails, but head 0
  // still reaches node 2 over the topology (e1, e3): it keeps its share 5 and loses it; head 1
  // delivers its 5.
  const auto rg = ordered(4, {{0, 2, 10, 1}, {0, 3, 10, 5}, {1, 2, 10, 1}, {3, 2, 10, 5}});
  auto t = make(rg, {0, 1}, {2}, t1(1, LOSSY));
  const auto rp = t->replace(10, 10).program;
  EXPECT_EQ(routes(*rp), (std::vector<Edges> {{0}, {2}}));
  const auto m0 = rg->eligibility({}, std::vector<EdgeId> {0});
  auto replay = make(rg, {0, 1}, {2}, FlowPolicyConfig {.routing = Pinned {rp, nullptr}, .admission = LOSSY}, &m0);
  const auto x = replay->replace(10);
  EXPECT_EQ(x.placed, 5);
  EXPECT_EQ(x.call_source_loss, 5);
  EXPECT_EQ(source_rows(x), (std::vector<SourceRow> {{0, 5, 5, 0, 5, 0}, {1, 5, 5, 5, 0, 0}}));
  EXPECT_EQ(x.termination, Termination::Complete);
}

// DIAM5's design-10 program (tunnels 0 -> 1 -> 3 and 0 -> 2 -> 3) with node 2 as a black
// hole: Lossy replace(10): 5 and 5 assigned, the second tunnel BLACKHOLED with its 5 dropped on e1
// (the link entering node 2), 5 delivered; Lossless admits 1. The live T1's forward(10,
// black_holes) gives the same and keeps its holds and its program's statuses.
TEST(TunnelPrograms, T1Stale) {
  const auto g = diam5();
  auto p = make(g, {0}, {3}, t1(2, LOSSY));
  const auto prog = p->replace(10, 10).program;
  const auto dark = std::make_shared<const Eligibility>(g->eligibility(std::vector<NodeId> {2}, {}));
  auto ly = make(g, {0}, {3}, FlowPolicyConfig {.routing = Pinned {prog, dark}, .admission = LOSSY});
  const auto r = ly->replace(10);
  EXPECT_EQ(rows(r), (Rows {{0, 0, UP, 5, 5, 0}, {1, 0, BLACKHOLED, 5, 0, 5}}));
  EXPECT_EQ(r.call_drops, (Pairs {{1, 5}}));
  EXPECT_EQ(r.placed, 5);
  auto ll = make(g, {0}, {3}, FlowPolicyConfig {.routing = Pinned {prog, dark}});
  EXPECT_EQ(ll->replace(10).placed, 1);
  const auto holds_before = p.ledger->reserved();
  const auto f = p->forward(10, dark.get());
  EXPECT_EQ(rows(f), rows(r));
  EXPECT_EQ(f.placed, 5);
  EXPECT_EQ(p.ledger->reserved(), holds_before);
  EXPECT_EQ(statuses(*p->program()), (std::vector<TunnelStatus> {UP, UP}));
}

// Pinned over the 6.6 program forwards like the live T1 (3.4): 6
// delivered, head 0's 6 lost at the head; the setup-down tunnels keep their status in by_tunnel;
// no plan, no hold.
TEST(TunnelPrograms, T1DownTunnels) {
  const auto g = twoh(2, 10);
  auto p = make(g, {0, 1}, {2}, t1(2, LOSSY));
  const auto live = p->replace(12, 12);
  auto pin = make(g, {0, 1}, {2}, FlowPolicyConfig {.routing = Pinned {live.program, nullptr}, .admission = LOSSY});
  const auto r = pin->replace(12);
  EXPECT_EQ(r.placed, 6);
  EXPECT_EQ(r.call_source_loss, 6);
  EXPECT_EQ(rows(r), rows(live));
  EXPECT_FALSE(r.plan.has_value());
  EXPECT_EQ(pin.ledger->reserved(), (Amounts {0, 0}));
  EXPECT_EQ(r.termination, Termination::Complete);
  // A replay never brings a down tunnel up: Pinned over the diamond's all-down program (design 14) with a mask
  // that leaves every path: still NB, NB; Lossy loses the offer at the head (not rule (3)).
  const auto dg = diam5();
  auto d = make(dg, {0}, {3}, t1(2, LOSSY));
  const auto down = d->replace(5, 14).program;
  const auto mask = dg->eligibility({}, std::vector<EdgeId> {2});
  auto replay = make(dg, {0}, {3}, FlowPolicyConfig {.routing = Pinned {down, nullptr}, .admission = LOSSY}, &mask);
  const auto x = replay->replace(5);
  EXPECT_EQ(rows(x), (Rows {{0, 0, NB, 0, 0, 0}, {1, 0, NB, 0, 0, 0}}));
  EXPECT_EQ(x.termination, Termination::Complete);
  EXPECT_EQ(x.call_source_loss, 5);
}

// The other producers' tunnels carry the new fields' plain values: status UP, requested
// and reserved 0, a destination, a path and a cost.
TEST(TunnelPrograms, OtherProducersFields) {
  const auto g = ordered(3, {{0, 2, 2, 1}, {1, 2, 10, 1}});
  auto gl = make(g, {0, 1}, {2}, global_config());
  for (const auto& t : gl->place(12).program->tunnels()) {
    EXPECT_EQ(t.status, UP);
    EXPECT_EQ(t.requested, 0);
    EXPECT_EQ(t.reserved, 0);
    EXPECT_TRUE(t.destination.has_value() && t.edges.has_value() && t.route_cost.has_value());
  }
  const auto supplied = pinned_paths(*g, {0, 1}, {2}, {{0, {0}, 1}, {1, {1}, 0}});
  for (const auto& t : supplied.program->tunnels()) {
    EXPECT_EQ(t.status, UP);
    EXPECT_EQ(t.requested, 0);
    EXPECT_EQ(t.reserved, 0);
    EXPECT_EQ(t.destination, std::optional<NodeId> {2});
  }
}

// ---------------------------------------------------------------- zero offers and everything masked (E)

// The 6.6 policy with zero offers: replace(0, 12) programs the tunnels and books the
// holds, four zero records; evaluate(0, 12) the same without publishing; forward(0) keeps program
// and holds, recreates the records (new ids), sets no flag. Both admissions; every call Complete.
TEST(TunnelZero, T1) {
  for (const auto a : {LOSSLESS, LOSSY}) {
    auto p = make(twoh(2, 10), {0, 1}, {2}, t1(2, a));
    const auto r = p->replace(0, 12);
    EXPECT_EQ(r.termination, Termination::Complete);
    EXPECT_EQ(r.placed, 0);
    EXPECT_EQ(statuses(*r.program), (std::vector<TunnelStatus> {NB, UP, NB, UP}));
    EXPECT_EQ(p.ledger->reserved(), (Amounts {0, 6}));
    EXPECT_EQ(entries(p->records()), (std::vector<std::optional<std::int32_t>> {0, 1, 2, 3}));
    for (const auto& rec : p->records()) EXPECT_EQ(rec.delivered, 0);
    const auto state = ledger_state(*p.ledger);
    const auto e = p->evaluate(0, 12);
    EXPECT_EQ(e.termination, Termination::Complete);
    EXPECT_EQ(statuses(*e.program), (std::vector<TunnelStatus> {NB, UP, NB, UP}));
    expect_same(ledger_state(*p.ledger), state, "evaluate(0, 12)");
    const auto before = p->records();
    const auto f = p->forward(0);
    EXPECT_EQ(f.termination, Termination::Complete);
    EXPECT_FALSE(f.routes_changed || f.weights_changed || f.reservations_changed);
    EXPECT_EQ(f.program, r.program);
    EXPECT_EQ(p.ledger->reserved(), (Amounts {0, 6}));
    // The records are retired and recreated: the same entries and holds under new ids.
    const auto after = p->records();
    EXPECT_EQ(entries(after), entries(before));
    EXPECT_EQ(holds(after), holds(before));
    for (const auto& rec : after) {
      EXPECT_TRUE(std::none_of(before.begin(), before.end(), [&](const FlowRecord& b) { return b.flow_id == rec.flow_id; }));
    }
  }
}

// DIAM5 all down at D = 14: replace(0, 14) Complete; replace(5, 14) Lossy Complete, Lossless
// Exhausted; two zero records (entries 0 and 1) every time; forward(0) Complete.
TEST(TunnelZero, T1AllDown) {
  for (const auto a : {LOSSLESS, LOSSY}) {
    auto p = make(diam5(), {0}, {3}, t1(2, a));
    EXPECT_EQ(p->replace(0, 14).termination, Termination::Complete);
    EXPECT_EQ(entries(p->records()), (std::vector<std::optional<std::int32_t>> {0, 1}));
    const auto r = p->replace(5, 14);
    EXPECT_EQ(r.termination, a == LOSSY ? Termination::Complete : Termination::Exhausted);
    EXPECT_EQ(entries(p->records()), (std::vector<std::optional<std::int32_t>> {0, 1}));
    EXPECT_EQ(p->forward(0).termination, Termination::Complete);
    EXPECT_EQ(p->records().size(), 2U);
  }
}

// TWOH(10, 10) with every edge masked: an empty roster: replace(5, 5) NO_ROUTE, then
// replace(0, 5) Complete; the empty program; no record.
TEST(TunnelZero, T1AllMasked) {
  const auto g = twoh(10, 10);
  const auto mask = g->eligibility({}, std::vector<EdgeId> {0, 1});
  auto p = make(g, {0, 1}, {2}, t1(2), &mask);
  const auto r = p->replace(5, 5);
  EXPECT_EQ(r.termination, Termination::NoRoute);
  EXPECT_TRUE(r.program->tunnels().empty());
  EXPECT_EQ(p->replace(0, 5).termination, Termination::Complete);
  EXPECT_TRUE(p->records().empty());
}

// Pinned over the 6.6 program. With e1 masked, tunnels 1 and 3 fail and head
// 1 cannot reach node 2 (it leaves the replay roster); head 0 still can, over e0, but has no up
// tunnel: it keeps the whole 12 and loses it. The program has up tunnels and no head can use one:
// rule (3), NO_ACTIVE_WEIGHT in both admissions; replace(0) is Complete; four zero records. With
// nodes 0 and 1 masked the replay roster is empty: NO_ROUTE, then Complete, four records still.
TEST(TunnelZero, PinnedT1AllFailedOrNoRoster) {
  const auto g = twoh(2, 10);
  auto p = make(g, {0, 1}, {2}, t1(2, LOSSY));
  const auto prog = p->replace(12, 12).program;
  const auto no_e1 = g->eligibility({}, std::vector<EdgeId> {1});
  const auto no_heads = g->eligibility(std::vector<NodeId> {0, 1}, {});
  for (const auto a : {LOSSLESS, LOSSY}) {
    auto pin = make(g, {0, 1}, {2}, FlowPolicyConfig {.routing = Pinned {prog, nullptr}, .admission = a}, &no_e1);
    const auto r = pin->replace(12);
    EXPECT_EQ(r.termination, Termination::NoActiveWeight);
    if (a == LOSSY) {
      EXPECT_EQ(r.call_source_loss, 12);
      EXPECT_EQ(source_rows(r), (std::vector<SourceRow> {{0, 12, 12, 0, 12, 0}, {1, 0, 0, 0, 0, 0}}));
    }
    EXPECT_EQ(rows(r), (Rows {{0, 0, NB, 0, 0, 0}, {1, 1, FAILED, 0, 0, 0}, {2, 0, NB, 0, 0, 0}, {3, 1, FAILED, 0, 0, 0}}));
    EXPECT_EQ(pin->replace(0).termination, Termination::Complete);
    EXPECT_EQ(entries(pin->records()), (std::vector<std::optional<std::int32_t>> {0, 1, 2, 3}));
    auto none = make(g, {0, 1}, {2}, FlowPolicyConfig {.routing = Pinned {prog, nullptr}, .admission = a}, &no_heads);
    EXPECT_EQ(none->replace(12).termination, Termination::NoRoute);
    EXPECT_EQ(none->replace(0).termination, Termination::Complete);
    EXPECT_EQ(none->records().size(), 4U);
  }
}

// ---------------------------------------------------------------- evaluators (12)

// Capacity of T1 is the lossless bound of the tunnels set up for the design
// (required and finite); Lossy forwards `offered` (required and finite); the other models reject a
// design; batch_capacity shares nothing between pairs.
TEST(TunnelEvaluators, T1Capacity) {
  const auto g = par612();
  QueryOptions q {t1(4)};
  q.design_demand = Amount {12};
  EXPECT_EQ(capacity(*g, v({0}), v({1}), q).placed, 12);
  QueryOptions ly {t1(4, LOSSY)};
  ly.design_demand = Amount {12};
  ly.offered = Amount {13};
  const auto r = capacity(*g, v({0}), v({1}), ly);
  EXPECT_EQ(r.placed, 12);
  EXPECT_EQ(r.call_drops, (Pairs {{0, 1}}));
  ASSERT_TRUE(r.plan.has_value());
  EXPECT_EQ(r.plan->reserved_total, 12);
  auto err = [&](const QueryOptions& o) { return error_of([&] { (void)capacity(*g, v({0}), v({1}), o); }); };
  for (const auto a : {LOSSLESS, LOSSY}) EXPECT_EQ(err(QueryOptions {t1(4, a)}), "T1 capacity requires design_demand");
  QueryOptions inf {t1(4)};
  inf.design_demand = std::optional<Amount> {};
  EXPECT_EQ(err(inf), "T1 requires a finite design demand");
  QueryOptions no_offer {t1(4, LOSSY)};
  no_offer.design_demand = Amount {12};
  EXPECT_EQ(err(no_offer), "T1 LOSSY capacity requires finite offered");
  no_offer.offered = std::optional<Amount> {};
  EXPECT_EQ(err(no_offer), "T1 LOSSY capacity requires finite offered");
  QueryOptions with_offer {t1(4)};
  with_offer.design_demand = Amount {12};
  with_offer.offered = Amount {5};
  EXPECT_EQ(err(with_offer), "T1 LOSSLESS capacity does not accept offered");
  QueryOptions bad {t1(4)};
  bad.design_demand = Amount {-1};
  EXPECT_EQ(err(bad), "design_demand must be an integer in [0, 2^62) or permitted +inf");
  // The domain of offered comes first.
  bad.offered = kL62;
  EXPECT_EQ(err(bad), "offered must be an integer in [0, 2^62) or permitted +inf");
  for (const auto& c : {cost_only_config(), global_config()}) {
    QueryOptions o {c};
    o.design_demand = Amount {5};
    EXPECT_EQ(err(o), std::string(detail::model_name(c.routing)) + " does not accept design_demand");
    const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{v({0}), v({1})}};
    EXPECT_EQ(error_of([&] { (void)batch_capacity(*g, pairs, o, {}, 1); }),
              std::string(detail::model_name(c.routing)) + " does not accept design_demand");
  }
  QueryOptions pin {FlowPolicyConfig {.routing = pinned_paths(*g, {0}, {1}, {{0, {0}, 1}})}};
  pin.design_demand = Amount {5};
  EXPECT_EQ(err(pin), "Pinned does not accept design_demand");
  // Each pair sets up on a fresh ledger: ONE10, N = 1, design 6, twice the same pair: 10 and 10.
  const auto one = one10();
  QueryOptions b {t1(1)};
  b.design_demand = Amount {6};
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{v({0}), v({1})}, {v({0}), v({1})}};
  for (const std::int32_t workers : {1, 2}) {
    const auto rs = batch_capacity(*one, pairs, b, {}, workers);
    ASSERT_EQ(rs.size(), 2U);
    EXPECT_EQ(rs[0].placed, 10);
    EXPECT_EQ(rs[1].placed, 10);
  }
}

// sensitivity(X8, [0, 1], 2, T1 N = 1, design 4): baseline 0
// (head 0 is down and first); deleting e0 takes head 0 out of the roster and head 1's tunnel carries
// 8: loss -8; deleting e1: 0. Every allowed link is a candidate (no shortcut), here an edge outside
// every route too.
TEST(TunnelEvaluators, T1Sensitivity) {
  QueryOptions q {t1(1)};
  q.design_demand = Amount {4};
  for (const std::int32_t workers : {1, 2}) {
    const auto s = sensitivity(*x8(), v({0, 1}), v({2}), q, workers);
    EXPECT_EQ(s.baseline.placed, 0);
    std::vector<std::pair<EdgeId, Amount>> got;
    for (const auto& e : s.entries) got.emplace_back(e.edge, e.loss);
    EXPECT_EQ(got, (std::vector<std::pair<EdgeId, Amount>> {{0, -8}, {1, 0}}));
  }
  const auto g = ordered(3, {{0, 1, 6, 1}, {0, 1, 6, 2}, {2, 0, 5, 1}});
  QueryOptions all {t1(4)};
  all.design_demand = Amount {12};
  const auto s = sensitivity(*g, v({0}), v({1}), all, 1);
  std::vector<std::pair<EdgeId, Amount>> got;
  for (const auto& e : s.entries) got.emplace_back(e.edge, e.loss);
  // Deleting e0: the four tunnels set up on e1 alone (two fit: 6 of 3 + 3), the bound 6.
  EXPECT_EQ(got, (std::vector<std::pair<EdgeId, Amount>> {{0, 6}, {1, 6}, {2, 0}}));
  // A zero offered (Lossy) is Complete and every loss is 0.
  QueryOptions zero {t1(1, LOSSY)};
  zero.design_demand = Amount {4};
  zero.offered = Amount {0};
  const auto c0 = capacity(*x8(), v({0, 1}), v({2}), zero);
  EXPECT_EQ(c0.termination, Termination::Complete);
  EXPECT_EQ(c0.placed, 0);
  for (const auto& e : sensitivity(*x8(), v({0, 1}), v({2}), zero, 1).entries) EXPECT_EQ(e.loss, 0);
  // An Eligibility that excludes every edge: an empty roster, NO_ROUTE, placed 0.
  const auto gm = x8();
  const auto mask = gm->eligibility({}, std::vector<EdgeId> {0, 1});
  const auto none = capacity(*gm, v({0, 1}), v({2}), QueryOptions {t1(1), &mask, nullptr, {}, Amount {4}});
  EXPECT_EQ(none.termination, Termination::NoRoute);
  EXPECT_EQ(none.placed, 0);
}

// ---------------------------------------------------------------- validation (V)

// T1's row of the matrix and its own parameters.
TEST(TunnelValidation, T1Matrix) {
  const auto g = one10();
  FlowLedger fg(g);
  auto rejects = [&](const FlowPolicyConfig& c, const std::string& text, DagCache* cache = nullptr) {
    EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), c, nullptr, cache); }), text);
    EXPECT_TRUE(detail::LedgerAccess::bound_owners(fg).empty());
  };
  rejects(t1(2, LOSSLESS, std::nullopt, 0, PathSelection(DET, EdgePriority::Capacity, true)), "T1 requires multipath=False");
  for (const std::optional<bool> mp : {std::optional<bool> {}, std::optional<bool> {false}}) {
    EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), t1(2, LOSSLESS, std::nullopt, 0, PathSelection(DET, EdgePriority::Capacity, mp))); }), "");
  }
  {
    FlowPolicy p(fg, v({0}), v({1}), t1(2));
    EXPECT_TRUE(p.config().selection == PathSelection());
    EXPECT_FALSE(p.config().selection->multipath.has_value());
  }
  rejects({.routing = CapacityAware {Tunnels {2}}, .weights = WeightBasis::Capacity}, "T1 does not accept weights");
  rejects({.routing = CapacityAware {Tunnels {2}}, .edge_weights = std::make_shared<const std::vector<Amount>>(Amounts {1})},
          "T1 does not accept edge_weights");
  DagCache cache(*g);
  rejects(t1(2), "T1 does not use DagCache", &cache);
  rejects({.routing = CapacityAware {Tunnels {2, true}}}, "T1 requires an exact tunnel count");
  for (const std::int64_t n : {std::int64_t {0}, std::int64_t {10001}, std::int64_t {-1}}) rejects(t1(n), "tunnel count must be in [1, 10000]");
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), t1(10000)); }), "");
  rejects({.routing = CapacityAware {Tunnels {1}}, .origin = Origin::Optimized}, "T1 requires origin=EQUAL");
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({1}), t1(2, LOSSY)); }), "");
}

// The tunnel limit counts the roster's heads (after the mask) times the count:
// FAN(107375) with N = 10000 is above 2^30 - 1; FAN(107374), FAN(112871) with N = 9513 (exactly
// 2^30 - 1: the limit is inclusive) and FAN(107375) with head 0 masked construct (no tunnel is
// made at construction).
TEST(TunnelValidation, TunnelLimit) {
  const auto fan = [](std::int32_t h) {
    std::vector<InEdge> rows;
    rows.reserve(static_cast<std::size_t>(h));
    for (std::int32_t i = 0; i < h; ++i) rows.push_back({i, h, 1, 1});
    return shared(make_graph(h + 1, rows));
  };
  const auto heads = [](std::int32_t h) {
    std::vector<NodeId> out(static_cast<std::size_t>(h));
    for (std::int32_t i = 0; i < h; ++i) out[static_cast<std::size_t>(i)] = i;
    return out;
  };
  const auto big = fan(107375);
  FlowLedger fg(big);
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, heads(107375), v({107375}), t1(10000)); }),
            "tunnel limit exceeded: 107375 heads x 10000 tunnels > 1073741823");
  const auto mask = big->eligibility(std::vector<NodeId> {0}, {});
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, heads(107375), v({107375}), t1(10000), &mask); }), "");
  const auto ok = fan(107374);
  FlowLedger fo(ok);
  EXPECT_EQ(error_of([&] { FlowPolicy p(fo, heads(107374), v({107374}), t1(10000)); }), "");
  // Exactly at the limit (inclusive): FAN(112871) with N = 9513 is 112871 x 9513 = 2^30 - 1 tunnels.
  const auto edge = fan(112871);
  FlowLedger fe(edge);
  EXPECT_EQ(error_of([&] { FlowPolicy p(fe, heads(112871), v({112871}), t1(9513)); }), "");
  // The limit at an evaluator, per pair.
  QueryOptions q {t1(10000)};
  q.design_demand = Amount {1};
  EXPECT_EQ(error_of([&] { (void)capacity(*big, heads(107375), v({107375}), q); }),
            "tunnel limit exceeded: 107375 heads x 10000 tunnels > 1073741823");
}

// ---------------------------------------------------------------- numeric ranges (W)

// ONEMAX (one link of 2^62 - 1): one tunnel holds 2^62 - 1 and carries it; two
// owners near the baseline: A holds 2^61, B's request of 2^61 sees U = 2^61 - 1 (NB); the commit's
// sum of holds stays below 2^63.
TEST(TunnelRanges, OneMaxAndHoldsNearTheBaseline) {
  const Amount max = kL62 - 1;
  const auto g = ordered(2, {{0, 1, max, 1}});
  auto p = make(g, {0}, {1}, t1(1, LOSSY));
  const auto r = p->replace(max, max, kHolds);
  EXPECT_EQ(statuses(*r.program), (std::vector<TunnelStatus> {UP}));
  EXPECT_EQ(r.plan->reserved_total, max);
  EXPECT_EQ(p.ledger->reserved(), (Amounts {max}));
  EXPECT_EQ(p.ledger->unreserved(), (Amounts {0}));
  EXPECT_EQ(r.placed, max);
  auto ll = make(g, {0}, {1}, t1(1));
  EXPECT_EQ(ll->replace(INF, max).placed, max);
  auto ledger = std::make_shared<FlowLedger>(g);
  auto a = make(g, {0}, {1}, t1(1), nullptr, ledger);
  auto b = make(g, {0}, {1}, t1(1), nullptr, ledger);
  const Amount half = Amount {1} << 61;
  EXPECT_EQ(statuses(*a->replace(half, half).program), (std::vector<TunnelStatus> {UP}));
  EXPECT_EQ(statuses(*b->replace(half, half).program), (std::vector<TunnelStatus> {NB}));
  EXPECT_EQ(ledger->unreserved(), (Amounts {half - 1}));
  expect_balanced(*ledger, "near the baseline");
}

// FAN4MAX (four heads, private links of 2^60 - 1), N = 10000, Lossless replace(2^62 - 4,
// 2^62 - 4): 40000 tunnels, all up; each head's share 2^60 - 1 splits as 6975 requests of
// 115292150460685 and 3025 of 115292150460684, which fill each link exactly; reserved_total 2^62 -
// 4; the whole design is carried (6.5: every tunnel up and U <= V).
TEST(TunnelRanges, ManyTunnels) {
  const Amount cap = (Amount {1} << 60) - 1;
  const auto g = ordered(5, {{0, 4, cap, 1}, {1, 4, cap, 1}, {2, 4, cap, 1}, {3, 4, cap, 1}});
  auto p = make(g, {0, 1, 2, 3}, {4}, t1(10000));
  const Amount D = kL62 - 4;
  const auto r = p->replace(D, D);
  ASSERT_EQ(r.program->tunnels().size(), 40000U);
  std::map<Amount, std::int64_t> counts;
  for (const auto& t : r.program->tunnels()) {
    EXPECT_EQ(t.status, UP);
    ++counts[t.requested];
  }
  EXPECT_EQ(counts, (std::map<Amount, std::int64_t> {{115292150460684, 4 * 3025}, {115292150460685, 4 * 6975}}));
  EXPECT_EQ(r.plan->reserved_total, D);
  EXPECT_EQ(p.ledger->reserved(), (Amounts {cap, cap, cap, cap}));
  EXPECT_EQ(r.placed, D);
  EXPECT_EQ(r.termination, Termination::Complete);
}

// ---------------------------------------------------------------- property

// The reference's exact Stretch comparison where a careless product overflows: factors far above
// 2^63 with a large or a zero baseline, the shortcut's boundary (f = 2^62 is computed, f = 2^63
// is above every cost), a factor just above 1 on a baseline near 2^62 and a fractional factor at
// equality. The library's comparison (baseline >= 1, as costs are positive) gives the same.
TEST(T1Reference, StretchComparisonIsExact) {
  using t1_reference::stretch_admits;
  constexpr Cost kMax = std::numeric_limits<Cost>::max();
  constexpr Cost p62 = Cost {1} << 62;
  struct Case {
    Cost cost, base;
    double f;
    bool admits;
  };
  const std::vector<Case> cases {
      {Cost {1} << 60, Cost {1} << 60, std::ldexp(1.0, 80), true},  // 2^60 <= 2^140
      {kMax, Cost {1} << 60, std::ldexp(1.0, 120), true},
      {kMax, 1, std::ldexp(1.0, 63), true},  // 2^63 - 1 < 2^63
      {p62, 1, std::ldexp(1.0, 62), true},  // equality
      {p62 + 1, 1, std::ldexp(1.0, 62), false},
      {kMax, p62, std::ldexp(1.0, 40), true},
      {kMax, Cost {1} << 22, std::ldexp(1.0, 40), false},  // 2^62 < 2^63 - 1
      {p62 + 1022, p62 - 1, 1.0 + std::ldexp(1.0, -52), true},  // f * base = 2^62 + 1023 - 2^-52
      {p62 + 1023, p62 - 1, 1.0 + std::ldexp(1.0, -52), false},
      {3, 2, 1.5, true},
      {4, 2, 1.5, false},
  };
  for (const auto& c : cases) {
    EXPECT_EQ(stretch_admits(c.cost, c.base, c.f), c.admits) << c.cost << " " << c.base << " " << c.f;
    EXPECT_EQ(detail::stretch_le(c.cost, c.base, c.f), c.admits) << c.cost << " " << c.base << " " << c.f;
  }
  // A zero baseline admits only a zero cost, whatever the factor (the reference only).
  EXPECT_TRUE(stretch_admits(0, 0, std::ldexp(1.0, 120)));
  EXPECT_FALSE(stretch_admits(1, 0, std::ldexp(1.0, 120)));
  EXPECT_FALSE(stretch_admits(1, 0, std::ldexp(1.0, 30)));
}

// Random small graphs (costs 1 to 4, capacities 0 to 6, one or two heads and
// destinations, a node or edge mask, every selection field T1 reads, the three ceilings, a search
// minimum, preloaded service and reservation baselines, another owner's holds and another owner's
// carried load) against t1_reference (test_utils.hpp): the roster, every tunnel's request, status
// and path, the holds per edge (the up requests, within the reservation baseline), the lossless
// bound on the service view; the setup does not depend on the offer; when every tunnel is up and U
// <= V on every link, D itself is admitted (6.5); Lossy accounting and per-record conservation.
TEST(T1Property, SetupAndBoundAgainstTheReference) {
  std::mt19937 rng(20261018);
  auto below = [&](int k) { return static_cast<int>(rng() % static_cast<unsigned>(k)); };
  int up = 0, nb = 0, cc = 0, bounded = 0, design_carried = 0, dropped = 0, preferred = 0, collapsed = 0, held = 0;
  for (int trial = 0; trial < 700; ++trial) {
    const int n = 4 + below(3);
    const int m = 3 + below(10);
    std::vector<InEdge> rows;
    for (int i = 0; i < m; ++i) {
      const NodeId u = below(n);
      NodeId w = below(n - 1);
      if (w >= u) ++w;
      rows.push_back({u, w, static_cast<Amount>(below(7)), static_cast<Cost>(1 + below(4))});
    }
    const auto g = shared(make_graph(n, rows));
    std::vector<NodeId> src = {0};
    std::vector<NodeId> dst = {n - 1};
    if (below(2) == 0) src.push_back(1);
    if (below(2) == 0) dst.insert(dst.begin(), n - 2);
    std::optional<Eligibility> mask;
    switch (below(4)) {
      case 0: mask = g->eligibility(std::vector<NodeId> {static_cast<NodeId>(below(n))}, {}); break;
      case 1: mask = g->eligibility({}, std::vector<EdgeId> {static_cast<EdgeId>(below(m))}); break;
      default: break;
    }
    const Eligibility* elig = mask ? &*mask : nullptr;
    const t1_reference::Selection rs {below(2) == 0, std::vector<EdgePriority> {EdgePriority::Residual, EdgePriority::Capacity,
                                                                                EdgePriority::Headroom}[static_cast<std::size_t>(below(3))],
                                      below(3) != 0};
    std::optional<CostCeiling> ceiling;
    switch (below(5)) {
      case 0: ceiling = Absolute {1 + below(8)}; break;
      case 1: ceiling = Stretch {std::vector<double> {1.0, 1.2, 1.5, 2.0}[static_cast<std::size_t>(below(4))]}; break;
      default: break;
    }
    const Amount minimum = below(4) == 0 ? below(3) : 0;
    const std::int64_t N = 1 + below(3);
    const Amount D = below(13);
    const auto selection = PathSelection(rs.prefer ? PHP : DET, rs.priority, std::nullopt, rs.multi_edge);
    std::vector<Amount> service(g->capacity_units().begin(), g->capacity_units().end());
    std::vector<Amount> reservation = service;
    if (below(2) == 0) {
      for (auto& x : service) x = x == 0 ? 0 : below(static_cast<int>(x) + 1);
    }
    if (below(2) == 0) {
      for (auto& x : reservation) x = x == 0 ? 0 : below(static_cast<int>(x) + 1);
    }
    const bool other_holds = below(3) == 0;
    const bool other_load = below(3) == 0;
    const Amount other_d = below(6);
    // One scenario: the baselines and the other owners, then the policy under test.
    struct Scenario {
      std::shared_ptr<FlowLedger> ledger;
      std::unique_ptr<FlowPolicy> holds, load, p;
    };
    auto build = [&](Admission a) {
      Scenario s;
      s.ledger = std::make_shared<FlowLedger>(g, service, reservation);
      if (other_holds) {
        s.holds = std::make_unique<FlowPolicy>(*s.ledger, src, dst, t1(1));
        (void)s.holds->replace(0, other_d);
      }
      if (other_load) {
        s.load = std::make_unique<FlowPolicy>(*s.ledger, src, dst, cost_only_config());
        (void)s.load->place(other_d);
      }
      s.p = std::make_unique<FlowPolicy>(*s.ledger, src, dst, t1(N, a, ceiling, minimum, selection), elig);
      return s;
    };
    const std::string where = "trial " + std::to_string(trial);
    auto s = build(LOSSLESS);
    const auto U = s.ledger->unreserved();
    const auto V = s.ledger->residual();
    const auto R0 = detail::LedgerAccess::reservation_baseline(*s.ledger);
    const auto held_before = s.ledger->reserved();
    std::vector<NodeId> dests_ok;
    for (NodeId d : dst) {
      if (elig == nullptr || elig->node_ok(d)) dests_ok.push_back(d);
    }
    const auto heads = t1_reference::roster(*g, elig, src, dests_ok);
    if (heads.empty()) continue;
    const auto want = t1_reference::setup(*g, elig, heads, dests_ok, rs, ceiling, minimum, D, N, U);
    const auto r = s.p->replace(INF, D);
    ASSERT_EQ(vec(detail::ProgramAccess::data(*r.program).roster), heads) << where;
    ASSERT_EQ(r.program->tunnels().size(), want.size()) << where;
    std::vector<Amount> sum(U.size(), 0);
    bool all_up = true;
    for (std::size_t i = 0; i < want.size(); ++i) {
      const auto& t = r.program->tunnels()[i];
      EXPECT_EQ(t.head, want[i].head) << where << " tunnel " << i;
      EXPECT_EQ(t.requested, want[i].request) << where << " tunnel " << i;
      EXPECT_EQ(t.status, want[i].status) << where << " tunnel " << i;
      EXPECT_EQ(t.edges.value_or(Edges {}), want[i].path) << where << " tunnel " << i;
      if (t.status == UP) {
        ++up;
        for (EdgeId e : *t.edges) {
          sum[static_cast<std::size_t>(e)] += t.requested;
          EXPECT_GE(U[static_cast<std::size_t>(e)], std::max(t.requested, minimum)) << where;  // it passed its filter
        }
      } else {
        all_up = false;
        ++(t.status == NB ? nb : cc);
      }
    }
    const auto held_after = s.ledger->reserved();
    for (std::size_t e = 0; e < sum.size(); ++e) {
      EXPECT_EQ(held_after[e], held_before[e] + sum[e]) << where << " edge " << e;
      EXPECT_LE(held_after[e], R0[e]) << where << " edge " << e;
      held += sum[e] > 0 ? 1 : 0;
    }
    expect_balanced(*s.ledger, where);
    const Amount bound = t1_reference::lossless_bound(heads.size(), want, V);
    EXPECT_EQ(r.placed, bound) << where;
    bounded += bound > 0 ? 1 : 0;
    // The setup does not depend on the offer (evaluate releases the policy's own holds).
    const auto e = s.p->evaluate(below(9), D);
    EXPECT_EQ(routes(*e.program), routes(*r.program)) << where;
    EXPECT_EQ(statuses(*e.program), statuses(*r.program)) << where;
    // 6.5: every tunnel up and U <= V on every link: D is carried without loss.
    bool covered = all_up;
    for (std::size_t x = 0; covered && x < U.size(); ++x) covered = U[x] <= V[x];
    if (covered) {
      EXPECT_EQ(s.p->evaluate(D, D).placed, D) << where;
      ++design_carried;
    }
    preferred += rs.prefer ? 1 : 0;
    collapsed += rs.multi_edge ? 0 : 1;
    // Lossy: the same setup; the accounting and every record's conservation.
    auto ly = build(LOSSY);
    const Amount offer = below(15);
    const auto lr = ly.p->replace(offer, D);
    EXPECT_EQ(routes(*lr.program), routes(*r.program)) << where;
    Amount drops = 0;
    for (const auto& [x, d] : lr.drops) drops += lo(d);
    EXPECT_EQ(lo(lr.offered), offer) << where;
    EXPECT_EQ(lo(lr.source_loss) + drops + lr.delivered, offer) << where;
    dropped += drops > 0 ? 1 : 0;
    for (const auto& rec : ly.p->records()) (void)expect_conserved(*g, rec, where);
    expect_balanced(*ly.ledger, where);
    ASSERT_FALSE(::testing::Test::HasFailure()) << where;
  }
  // Not vacuous.
  EXPECT_GT(up, 300);
  EXPECT_GT(nb, 100);
  EXPECT_GT(cc, 20);
  EXPECT_GT(bounded, 100);
  EXPECT_GT(design_carried, 30);
  EXPECT_GT(dropped, 30);
  EXPECT_GT(preferred, 100);
  EXPECT_GT(collapsed, 100);
  EXPECT_GT(held, 100);
}
