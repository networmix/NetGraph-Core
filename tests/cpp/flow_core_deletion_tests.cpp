// The deletion re-solve of the flow core and the exact sensitivity model that uses it. Every check
// is independent of the repair: a solve from zero on the deleted view, the Bellman-Ford curve of
// test_utils.hpp, residuals rebuilt from flows, divergences and reduced costs computed here in
// the exact test integer of test_int128.hpp. Without the private hooks, check_deletion checks
// each re-solve's final result only: the probe that observes every intermediate state of the
// repairs is a hook. Three tests check the test tools themselves: the portable 128-bit classes of
// test_int128.hpp against the built-in type, the equality helpers of test_utils.hpp field by
// field, and this file's check_residual against a broken state per clause.
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "netgraph/core/detail/flow_core.hpp"
#include "netgraph/core/evaluators.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

using detail::DeletionCapture;
using detail::DeletionResolver;
using detail::DeletionResult;
using detail::DeletionState;
using detail::FlowCoreCounters;
using detail::I128;
using detail::Network;

constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();

Int128 wide(I128 x) { return int128_of(x.hi, x.lo); }

// A solve to capture: the graph, the view, the endpoint sets, the eligibility and the cap.
struct Case {
  std::shared_ptr<const Graph> g;
  std::vector<Amount> view;
  std::vector<NodeId> sources;
  std::vector<NodeId> destinations;
  std::optional<Eligibility> elig {};
  std::optional<Amount> cap {};
  [[nodiscard]] const Eligibility* el() const { return elig ? &*elig : nullptr; }
};

Case make_case(std::int32_t n, const std::vector<InEdge>& edges, std::vector<NodeId> s, std::vector<NodeId> d) {
  Case c;
  c.g = std::make_shared<const Graph>(make_graph(n, edges));
  c.view.assign(c.g->capacity_units().begin(), c.g->capacity_units().end());
  c.sources = std::move(s);
  c.destinations = std::move(d);
  return c;
}

DeletionCapture capture_of(const Case& c, FlowCoreResult* payload = nullptr) {
  DeletionCapture req;
  req.graph = c.g;
  auto r = place_global_units(*c.g, c.view, c.sources, c.destinations, c.el(), c.cap, &req);
  if (payload != nullptr) *payload = std::move(r);
  return req;
}

Int128 cost_sum(const Graph& g) {
  Int128 c = 0;
  for (Cost w : g.cost()) c += w;
  return c;
}

// The local index of edge e in the state's network.
std::size_t local_of(const Network& net, EdgeId e) {
  const auto it = std::lower_bound(net.edges.begin(), net.edges.end(), e);
  return it != net.edges.end() && *it == e ? static_cast<std::size_t>(it - net.edges.begin()) : kNone;
}

// The tail of every slot.
std::vector<std::int32_t> slot_tails(const Network& net) {
  std::vector<std::int32_t> tail(net.arc_head.size());
  for (std::int32_t x = 0; x < net.num_nodes + 2; ++x) {
    for (auto p = net.off[static_cast<std::size_t>(x)]; p < net.off[static_cast<std::size_t>(x) + 1]; ++p) {
      tail[static_cast<std::size_t>(p)] = x;
    }
  }
  return tail;
}

// The residual of a flow (per EdgeId) on the network: real arcs from the flow and the view
// capacity, virtual arcs from the endpoints' divergences; the deleted edge's arcs zero.
std::vector<Amount> residual_of(const Network& net, const std::vector<Amount>& flows, std::size_t deleted) {
  std::vector<Amount> res(net.arc_head.size(), 0);
  std::vector<Amount> div(static_cast<std::size_t>(net.num_nodes) + 2, 0);
  for (std::size_t k = 0; k < net.edges.size(); ++k) {
    const Amount f = flows[static_cast<std::size_t>(net.edges[k])];
    const auto fwd = static_cast<std::size_t>(net.edge_slot[k]);
    const auto bwd = static_cast<std::size_t>(net.rev[fwd]);
    if (k == deleted) continue;
    res[fwd] = net.cap[k] - f;
    res[bwd] = f;
    div[static_cast<std::size_t>(net.tail[k])] += f;
    div[static_cast<std::size_t>(net.head[k])] -= f;
  }
  for (std::size_t i = 0; i < net.sources.size(); ++i) {
    const Amount x = div[static_cast<std::size_t>(net.sources[i])];
    const auto a = static_cast<std::size_t>(net.src_rev[i]);
    res[a] = x;
    res[static_cast<std::size_t>(net.rev[a])] = kMaxFlowCap - x;
  }
  for (std::size_t i = 0; i < net.destinations.size(); ++i) {
    const Amount y = -div[static_cast<std::size_t>(net.destinations[i])];
    const auto a = static_cast<std::size_t>(net.dst_fwd[i]);
    res[a] = kMaxFlowCap - y;
    res[static_cast<std::size_t>(net.rev[a])] = y;
  }
  return res;
}

// The divergences (outflow minus inflow) every vertex should have: `at` lists (vertex, value).
std::vector<Int128> divergences(const Network& net, std::initializer_list<std::pair<std::int32_t, Int128>> at) {
  std::vector<Int128> d(static_cast<std::size_t>(net.num_nodes) + 2, 0);
  for (auto [v, x] : at) d[static_cast<std::size_t>(v)] += x;
  return d;
}

// The failed checks of a residual with potentials on the network: capacity
// bounds (both arcs of a pair non-negative, summing to the pair's capacity; the deleted edge's
// pair zero), the divergence vector `want`, non-negative reduced cost on every positive
// residual arc (real and virtual), and the bound `bound` on every potential (non-negative).
std::vector<std::string> check_residual(const Network& net, std::span<const Amount> res, std::span<const I128> pi,
                                        std::size_t deleted, const std::vector<Int128>& want, Int128 bound) {
  std::vector<std::string> fails;
  const auto n = static_cast<std::size_t>(net.num_nodes) + 2;
  if (res.size() != net.arc_head.size() || pi.size() != n) return {"shape"};
  std::vector<Int128> div(n, 0);
  for (std::size_t k = 0; k < net.edges.size(); ++k) {
    const auto fwd = static_cast<std::size_t>(net.edge_slot[k]);
    const auto bwd = static_cast<std::size_t>(net.rev[fwd]);
    if (k == deleted) {
      if (res[fwd] != 0 || res[bwd] != 0) fails.push_back("the deleted edge keeps a residual arc");
      continue;
    }
    if (res[fwd] < 0 || res[bwd] < 0 || res[fwd] + res[bwd] != net.cap[k]) {
      fails.push_back("capacity bounds at edge " + std::to_string(net.edges[k]));
    }
    div[static_cast<std::size_t>(net.tail[k])] += res[bwd];
    div[static_cast<std::size_t>(net.head[k])] -= res[bwd];
  }
  const auto S = static_cast<std::size_t>(net.num_nodes);
  const std::size_t T = S + 1;
  auto virtual_pair = [&](std::size_t a, std::size_t from, std::size_t to, Amount flow) {
    if (res[a] < 0 || res[static_cast<std::size_t>(net.rev[a])] < 0 ||
        static_cast<Int128>(res[a]) + res[static_cast<std::size_t>(net.rev[a])] != kMaxFlowCap) {
      fails.push_back("capacity bounds at a virtual arc");
    }
    div[from] += flow;
    div[to] -= flow;
  };
  for (std::size_t i = 0; i < net.sources.size(); ++i) {
    const auto a = static_cast<std::size_t>(net.src_rev[i]);
    virtual_pair(a, S, static_cast<std::size_t>(net.sources[i]), res[a]);
  }
  for (std::size_t i = 0; i < net.destinations.size(); ++i) {
    const auto a = static_cast<std::size_t>(net.dst_fwd[i]);
    virtual_pair(a, static_cast<std::size_t>(net.destinations[i]), T, res[static_cast<std::size_t>(net.rev[a])]);
  }
  if (div != want) fails.push_back("divergence vector");
  const auto tail = slot_tails(net);
  for (std::size_t a = 0; a < res.size(); ++a) {
    if (res[a] <= 0) continue;
    const auto u = static_cast<std::size_t>(tail[a]);
    const auto v = static_cast<std::size_t>(net.arc_head[a]);
    if (static_cast<Int128>(net.arc_cost[a]) + wide(pi[u]) - wide(pi[v]) < 0) {
      fails.push_back("negative reduced cost on slot " + std::to_string(a) + (u >= S || v >= S ? " (virtual)" : ""));
    }
  }
  for (std::size_t v = 0; v < n; ++v) {
    if (wide(pi[v]) < 0 || wide(pi[v]) > bound) fails.push_back("potential outside [0, bound] at " + std::to_string(v));
  }
  return fails;
}

// The vertices the real sources reach over positive residual arcs (their residual closure).
std::vector<char> closure_of(const Network& net, std::span<const Amount> res) {
  const auto n = static_cast<std::size_t>(net.num_nodes) + 2;
  std::vector<char> seen(n, 0);
  std::vector<std::size_t> queue;
  for (std::int32_t s : net.sources) {
    seen[static_cast<std::size_t>(s)] = 1;
    queue.push_back(static_cast<std::size_t>(s));
  }
  for (std::size_t qi = 0; qi < queue.size(); ++qi) {
    const auto u = queue[qi];
    for (auto p = net.off[u]; p < net.off[u + 1]; ++p) {
      const auto v = static_cast<std::size_t>(net.arc_head[static_cast<std::size_t>(p)]);
      if (res[static_cast<std::size_t>(p)] > 0 && seen[v] == 0) {
        seen[v] = 1;
        queue.push_back(v);
      }
    }
  }
  return seen;
}

// Whether the real sources reach t* over positive residual arcs.
bool connected(const Network& net, std::span<const Amount> res) { return closure_of(net, res).back() != 0; }

// A solve from zero with e deleted as sensitivity deletes it: excluded by an Eligibility
// (intersected with the case's), its view 0. An excluded edge is no edge of the network, so it is
// never a cut edge; an allowed edge of view 0 leaving the closure is one.
FlowCoreResult fresh_without(const Case& c, EdgeId e) {
  auto view = c.view;
  view[static_cast<std::size_t>(e)] = 0;
  const EdgeId one[1] = {e};
  Eligibility deleted = c.g->eligibility({}, one);
  if (c.el() != nullptr) deleted = Eligibility::intersect(*c.el(), deleted);
  return place_global_units(*c.g, view, c.sources, c.destinations, &deleted, std::nullopt);
}

// The cut of a residual: the network's edges, the deleted one aside, from the closure of the real
// sources to its complement, whatever their view capacity, ascending EdgeId.
std::vector<EdgeId> closure_cut(const Network& net, std::span<const Amount> res, std::size_t deleted) {
  const auto seen = closure_of(net, res);
  std::vector<EdgeId> cut;
  for (std::size_t k = 0; k < net.edges.size(); ++k) {
    if (k == deleted) continue;
    if (seen[static_cast<std::size_t>(net.tail[k])] != 0 && seen[static_cast<std::size_t>(net.head[k])] == 0) {
      cut.push_back(net.edges[k]);
    }
  }
  return cut;
}

// The Dijkstra-run bound of a re-solve of an edge carrying f units: at most m + 2 r runs, where r
// is the remainder (the units lost) and m = f - r the units repair 1 reroutes. The failed clause,
// "" when the bound holds.
std::string run_bound_failure(const DeletionState& st, Amount f, const DeletionResult& r) {
  const Amount rem = st.value() - r.value;
  if (rem < 0) return "the re-solve gained value";
  if (r.counters.dijkstra_runs > (f - rem) + 2 * rem) return "more Dijkstra runs than m + 2 r";
  return "";
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
// What the probe reported, per repair (0 after the deletion, 1 to 3 the repairs, 4 at
// completion): the events of each kind and the amount the repair had moved at its last event.
struct Observation {
  int updates[5] = {0, 0, 0, 0, 0};
  int pushes[5] = {0, 0, 0, 0, 0};
  Amount moved[5] = {0, 0, 0, 0, 0};
};

// The re-solve of e observed by the probe: every observed state (after the deletion, after every
// potential update and push of the three repairs, at completion) is checked against the
// capacity bounds, the exact divergence vector for the amounts moved so far and dual
// feasibility on every positive residual arc. `updates` collects the counters after every
// potential update with its repair; `seen` counts the events.
DeletionResult observed(DeletionResolver& rs, const DeletionState& st, EdgeId e, std::vector<std::string>& fails,
                        bool rerun = false, std::vector<std::pair<int, FlowCoreCounters>>* updates = nullptr,
                        Observation* seen = nullptr) {
  const Network& net = st.network();
  const std::size_t k = local_of(net, e);
  const std::int32_t u = net.tail[k];
  const std::int32_t v = net.head[k];
  const Amount F = st.value();
  const Amount f = st.edge_flows()[static_cast<std::size_t>(e)];
  const Int128 n = net.num_nodes + 2;
  const Int128 bound = (22 * n + 15) * cost_sum(st.graph());
  Amount moved[5] = {0, 0, 0, 0, 0};
  detail::DeletionProbe probe;
  probe.rerun_failed_search = rerun;
  probe.observe = [&](const detail::DeletionEvent& ev) {
    if (ev.repair >= 1 && ev.repair <= 3) moved[ev.repair] = ev.moved;
    const auto want = divergences(net, {{net.num_nodes, F - moved[2]},
                                        {net.num_nodes + 1, -static_cast<Int128>(F) + moved[3]},
                                        {u, -static_cast<Int128>(f) + moved[1] + moved[2]},
                                        {v, static_cast<Int128>(f) - moved[1] - moved[3]}});
    for (const auto& s : check_residual(net, ev.residual, ev.potentials, k, want, bound)) {
      fails.push_back("repair " + std::to_string(ev.repair) + (ev.push ? " push" : " update") + ": " + s);
    }
    if (updates != nullptr && !ev.push && ev.repair >= 1 && ev.repair <= 3) updates->emplace_back(ev.repair, *ev.counters);
    if (seen != nullptr && ev.repair >= 0 && ev.repair <= 4) {
      ++(ev.push ? seen->pushes : seen->updates)[ev.repair];
      seen->moved[ev.repair] = ev.moved;
    }
  };
  return rs.resolve(e, probe);
}
#endif

// A completed deletion's certificate, separate from the solve so broken answers can be fed to
// every clause. The view, a fresh solve and an independent cost curve supply the expected values.
void expect_deletion(const Case& c, const DeletionState& st, EdgeId e, const DeletionResult& r,
                     const std::string& where) {
  const Network& net = st.network();
  const std::size_t k = local_of(net, e);
  // Feasible on V \ {e} read per EdgeId: within the view, zero on the deleted edge and on every
  // edge the Eligibility excludes.
  auto view = c.view;
  view[static_cast<std::size_t>(e)] = 0;
  ASSERT_EQ(r.edge_flows.size(), view.size()) << where;
  for (std::size_t x = 0; x < view.size(); ++x) {
    EXPECT_TRUE(r.edge_flows[x] >= 0 && r.edge_flows[x] <= view[x]) << where << " edge " << x;
    if (c.el() != nullptr && !c.el()->edge_ok(static_cast<EdgeId>(x))) EXPECT_EQ(r.edge_flows[x], 0) << where;
  }
  // From zero on V \ {e} with e excluded, and the independent Bellman-Ford curve.
  const auto fresh = fresh_without(c, e);
  const auto ref = independent_curve(*c.g, view, c.sources, c.destinations, c.el());
  EXPECT_EQ(r.value, fresh.placed) << where;
  EXPECT_EQ(r.value, ref.back().first) << where;
  EXPECT_TRUE(flow_cost(*c.g, r.edge_flows) == flow_cost(*c.g, fresh.edge_flows)) << where;
  EXPECT_TRUE(flow_cost(*c.g, r.edge_flows) == ref.back().second) << where;
  ASSERT_TRUE(fresh.min_cut.has_value()) << where;
  EXPECT_EQ(r.cut, *fresh.min_cut) << where;
  Amount cut_value = 0;
  for (EdgeId x : r.cut) cut_value += view[static_cast<std::size_t>(x)];
  EXPECT_EQ(cut_value, r.value) << where;
  EXPECT_EQ(r.edge_flows[static_cast<std::size_t>(e)], 0) << where;
  // The final residual from the flow: bounds, divergences (value at s*), working potentials
  // valid on every positive residual arc, the destinations disconnected, and the cut its closure's.
  const auto res = residual_of(net, r.edge_flows, k);
  const Int128 n = net.num_nodes + 2;
  const auto want = divergences(net, {{net.num_nodes, r.value}, {net.num_nodes + 1, -static_cast<Int128>(r.value)}});
  for (const auto& s : check_residual(net, res, r.potentials, k, want, (22 * n + 15) * cost_sum(*c.g))) {
    ADD_FAILURE() << where << ": final: " << s;
  }
  EXPECT_FALSE(connected(net, res)) << where;
  EXPECT_EQ(r.cut, closure_cut(net, res, k)) << where;  // the cut of that residual
  EXPECT_EQ(run_bound_failure(st, st.edge_flows()[static_cast<std::size_t>(e)], r), "") << where;
}

// Every certificate check for the deletion of positive-flow edge e from the captured state of c.
void check_deletion(const Case& c, const DeletionState& st, EdgeId e, const std::string& where) {
  const Network& net = st.network();
  const std::size_t k = local_of(net, e);
  ASSERT_NE(k, kNone) << where;
  const Amount f = st.edge_flows()[static_cast<std::size_t>(e)];
  ASSERT_GT(f, 0) << where;
  DeletionResolver rs(st);
  const auto r = rs.resolve(e);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  // rem is the remainder (the units lost) and m the units repair 1 reroutes.
  const Amount rem = st.value() - r.value;
  const Amount m = f - rem;
  std::vector<std::string> fails;
  Observation seen;
  const auto ro = observed(rs, st, e, fails, false, nullptr, &seen);
  for (const auto& s : fails) ADD_FAILURE() << where << ": " << s;
  EXPECT_TRUE(same_result(ro, r)) << where;
  // The probe fired at every stage: one event after the deletion and one at completion; repair 1
  // searched at least once and pushed exactly when it moved units, m in all; repairs 2 and 3
  // pushed the remainder when there is one and are silent otherwise (repair 2 may reuse repair
  // 1's failed search, so it need not search).
  EXPECT_EQ(seen.updates[0], 1) << where;
  EXPECT_EQ(seen.updates[4], 1) << where;
  EXPECT_EQ(seen.pushes[0] + seen.pushes[4], 0) << where;
  EXPECT_GE(seen.updates[1], 1) << where;
  EXPECT_EQ(seen.pushes[1] > 0, m > 0) << where;
  EXPECT_EQ(seen.moved[1], m) << where;
  for (const int i : {2, 3}) {
    EXPECT_EQ(seen.pushes[i] > 0, rem > 0) << where << " repair " << i;
    EXPECT_EQ(seen.moved[i], rem) << where << " repair " << i;
  }
  EXPECT_EQ(seen.updates[3] > 0, rem > 0) << where;
  EXPECT_TRUE(rem > 0 || seen.updates[2] == 0) << where;
#endif
  expect_deletion(c, st, e, r, where);
}

// The captured state's own checks: bounds, divergences, valid working potentials within
// (n + 1) C, and the certificate.
void check_state(const Case& c, const DeletionState& st, const FlowCoreResult& payload, const std::string& where) {
  const Network& net = st.network();
  EXPECT_EQ(st.value(), payload.placed) << where;
  EXPECT_EQ(st.edge_flows(), payload.edge_flows) << where;
  const Int128 n = net.num_nodes + 2;
  const auto want =
      divergences(net, {{net.num_nodes, st.value()}, {net.num_nodes + 1, -static_cast<Int128>(st.value())}});
  for (const auto& s : check_residual(net, st.residual(), st.potentials(), kNone, want, (n + 1) * cost_sum(*c.g))) {
    ADD_FAILURE() << where << ": state: " << s;
  }
  EXPECT_EQ(st.maximum(), !connected(net, st.residual())) << where;
  EXPECT_TRUE(same_counters(st.counters(), payload.counters)) << where;
}

// Every allowed edge of the captured state: positive flow by the re-solve, zero flow by a
// solve from zero on V \ {e} and the independent Bellman-Ford curve there (the baseline's value
// and cost; no re-solve, which rejects it).
// `resolved`, when given, counts the edges checked by a re-solve.
void check_every_edge(const Case& c, const std::string& where, int* resolved = nullptr) {
  FlowCoreResult payload;
  const auto req = capture_of(c, &payload);
  ASSERT_TRUE(req.state.has_value()) << where;
  const DeletionState& st = *req.state;
  check_state(c, st, payload, where);
  ASSERT_TRUE(st.maximum()) << where;
  const auto flows = st.edge_flows();
  DeletionResolver reused(st);
  const auto positive = std::find_if(st.network().edges.begin(), st.network().edges.end(),
                                    [&](EdgeId e) { return flows[static_cast<std::size_t>(e)] > 0; });
  for (EdgeId e : st.network().edges) {
    const std::string at = where + " edge " + std::to_string(e);
    if (flows[static_cast<std::size_t>(e)] > 0) {
      check_deletion(c, st, e, at);
      if (resolved != nullptr) ++*resolved;
    } else {
      auto view = c.view;
      view[static_cast<std::size_t>(e)] = 0;
      const auto fresh = place_global_units(*c.g, view, c.sources, c.destinations, c.el(), std::nullopt);
      EXPECT_EQ(fresh.placed, st.value()) << at;
      EXPECT_TRUE(flow_cost(*c.g, fresh.edge_flows) == flow_cost(*c.g, flows)) << at;
      const auto ref = independent_curve(*c.g, view, c.sources, c.destinations, c.el());
      EXPECT_EQ(ref.back().first, st.value()) << at;
      EXPECT_TRUE(ref.back().second == flow_cost(*c.g, flows)) << at;
      EXPECT_EQ(error_of([&] { (void)reused.resolve(e); }), "deletion re-solve: edge " + std::to_string(e) + " carries no flow")
          << at;
      if (positive != st.network().edges.end()) {
        EXPECT_TRUE(same_result(reused.resolve(*positive), DeletionResolver(st).resolve(*positive))) << at;
      }
    }
    if (::testing::Test::HasFailure()) return;
  }
}

Case random_case(std::mt19937_64& rng) {
  auto pick = [&](std::int64_t lo, std::int64_t hi) { return std::uniform_int_distribution<std::int64_t>(lo, hi)(rng); };
  const auto n = static_cast<std::int32_t>(pick(3, 8));
  std::vector<NodeId> nodes(static_cast<std::size_t>(n));
  for (std::int32_t v = 0; v < n; ++v) nodes[static_cast<std::size_t>(v)] = v;
  std::shuffle(nodes.begin(), nodes.end(), rng);
  const auto ns = static_cast<std::size_t>(pick(1, std::min<std::int64_t>(3, n - 1)));
  const auto nd = static_cast<std::size_t>(pick(1, static_cast<std::int64_t>(nodes.size() - ns)));
  std::vector<NodeId> s(nodes.begin(), nodes.begin() + static_cast<std::ptrdiff_t>(ns));
  std::vector<NodeId> d(nodes.begin() + static_cast<std::ptrdiff_t>(ns),
                        nodes.begin() + static_cast<std::ptrdiff_t>(ns + nd));
  std::vector<InEdge> edges;
  const auto m = pick(1, 24);
  for (std::int64_t i = 0; i < m; ++i) {
    edges.push_back({static_cast<NodeId>(pick(0, n - 1)), static_cast<NodeId>(pick(0, n - 1)), pick(0, 4), pick(1, 13)});
  }
  const Cost factor = m % 3 == 0 ? 1000003 : 1;
  for (auto& edge : edges) edge.cost *= factor;
  Case c = make_case(n, edges, s, d);
  for (auto& x : c.view) {
    if (pick(0, 1) == 1) x = pick(0, x);  // a reduced view
  }
  if (pick(0, 1) == 1) {
    std::vector<NodeId> xn;
    std::vector<EdgeId> xe;
    for (NodeId v = 0; v < n; ++v) {
      if (pick(0, 10) == 0) xn.push_back(v);
    }
    for (EdgeId e = 0; e < c.g->num_edges(); ++e) {
      if (pick(0, 9) == 0) xe.push_back(e);
    }
    c.elig = c.g->eligibility(xn, xe);
  }
  return c;
}

// The positive-flow edges of a state.
std::vector<EdgeId> positive_edges(const DeletionState& st) {
  std::vector<EdgeId> out;
  const auto flows = st.edge_flows();
  for (EdgeId e : st.network().edges) {
    if (flows[static_cast<std::size_t>(e)] > 0) out.push_back(e);
  }
  return out;
}

} // namespace

// ------------------------------------------------------------------ the arithmetic

TEST(FlowCoreDeletion, SignedWideArithmeticMatchesInt128) {
  std::mt19937_64 rng(11);
  auto any = [&] {
    I128 x;
    x.hi = static_cast<std::int64_t>(rng()) >> (rng() % 64);
    x.lo = rng();
    return x;
  };
  int equal = 0;
  for (int i = 0; i < 20000; ++i) {
    const I128 a = any();
    const I128 b = i % 8 == 0 ? a : any();  // equal operands too (two random ones never are)
    const Int128 x = wide(a);
    const Int128 y = wide(b);
    // Modulo 2^128, through the unsigned type (a signed built-in overflow is undefined).
    EXPECT_TRUE(wide(a + b) == static_cast<Int128>(static_cast<UInt128>(x) + static_cast<UInt128>(y)));
    EXPECT_TRUE(wide(a - b) == static_cast<Int128>(static_cast<UInt128>(x) - static_cast<UInt128>(y)));
    EXPECT_EQ(a < b, x < y);
    EXPECT_EQ(a == b, x == y);
    equal += a == b ? 1 : 0;
    const auto c = static_cast<std::int64_t>(rng());
    EXPECT_TRUE(wide(I128 {c}) == c);
  }
  EXPECT_GT(equal, 2000);
}

#if defined(__SIZEOF_INT128__)
// The portable classes of test_int128.hpp (the test integer where the compiler has no built-in
// 128-bit type) against the built-in type: every operation, the compound assignments included,
// on operands of every magnitude, the signed and unsigned orders and the conversions. The signed
// built-in sums and products are formed in the unsigned type (a signed overflow is undefined; the
// classes wrap).
TEST(PortableInt128, PortableClassesMatchTheBuiltInType) {
  using B = unsigned __int128;
  using BS = __int128;
  using P = PortableUInt128;
  using PS = PortableInt128;
  auto to_p = [](B x) { return P::of_words(static_cast<std::uint64_t>(x >> 64), static_cast<std::uint64_t>(x)); };
  auto same = [](P p, B x) { return p.high() == static_cast<std::uint64_t>(x >> 64) && p.low() == static_cast<std::uint64_t>(x); };
  auto same_s = [&](PS p, BS x) { return same(static_cast<P>(p), static_cast<B>(x)); };
  std::mt19937_64 rng(128);
  // Words of every bit length, so carries, borrows, one-word and two-word operands all occur.
  auto word = [&] { return rng() >> (rng() % 64); };
  auto any = [&] {
    const int shape = static_cast<int>(rng() % 4);
    const std::uint64_t hi = shape == 0 ? 0 : (shape == 1 ? ~std::uint64_t {0} : word());
    return (static_cast<B>(hi) << 64) | word();
  };
  int wide_divisors = 0;
  for (int i = 0; i < 20000; ++i) {
    const B x = any();
    const B y = i % 5 == 0 ? static_cast<B>(1 + rng() % 16) : any();
    const P a = to_p(x);
    const P b = to_p(y);
    ASSERT_TRUE(same(a + b, x + y) && same(a - b, x - y) && same(a * b, x * y)) << i;
    ASSERT_TRUE(same(a | b, x | y) && same(a & b, x & y)) << i;
    P t = a;  // the compound assignments, in sequence on one value
    ASSERT_TRUE(same(t += b, x + y) && same(t -= b, x) && same(t *= b, x * y)) << i;
    ASSERT_EQ(a == b, x == y);
    ASSERT_EQ(a < b, x < y);
    ASSERT_EQ(a <= b, x <= y);
    const int s = static_cast<int>(rng() % 128);
    ASSERT_TRUE(same(a << s, x << s) && same(a >> s, x >> s)) << i << " shift " << s;
    if (y != 0) {
      ASSERT_TRUE(same(a / b, x / y) && same(a % b, x % y)) << i;
      P q = a;
      ASSERT_TRUE(same(q /= b, x / y)) << i;
      wide_divisors += (y >> 64) != 0 ? 1 : 0;  // the two-word long division
    }
    ASSERT_EQ(static_cast<std::uint64_t>(a), static_cast<std::uint64_t>(x));
    ASSERT_EQ(static_cast<std::int32_t>(a), static_cast<std::int32_t>(x));
    const auto bx = static_cast<BS>(x);
    const auto by = static_cast<BS>(y);
    const PS sa(a);
    const PS sb(b);
    ASSERT_TRUE(same_s(sa + sb, static_cast<BS>(x + y)) && same_s(sa - sb, static_cast<BS>(x - y)) &&
                same_s(sa * sb, static_cast<BS>(x * y)) && same_s(-sa, static_cast<BS>(B {0} - x)))
        << i;
    PS st = sa;
    ASSERT_TRUE(same_s(st += sb, static_cast<BS>(x + y)) && same_s(st -= sb, bx) &&
                same_s(st *= sb, static_cast<BS>(x * y)))
        << i;
    ASSERT_EQ(sa == sb, bx == by);
    ASSERT_EQ(sa < sb, bx < by);
    ASSERT_EQ(sa > sb, bx > by);
    const BS min = static_cast<BS>(B {1} << 127);
    if (by != 0 && !(bx == min && by == -1)) {
      ASSERT_TRUE(same_s(sa / sb, bx / by) && same_s(sa % sb, bx % by)) << i;
    }
    ASSERT_EQ(static_cast<std::int64_t>(sa), static_cast<std::int64_t>(bx));
    const auto c = static_cast<std::int64_t>(rng()) >> (rng() % 64);
    const auto u = rng() >> (rng() % 64);
    ASSERT_TRUE(same_s(PS {c}, BS {c}) && same(P {c}, static_cast<B>(c)) && same(P {u}, B {u})) << i;
    ASSERT_TRUE(same_s(PS {u}, BS {u}) && same_s(PS {static_cast<std::int32_t>(c)}, BS {static_cast<std::int32_t>(c)}));
  }
  EXPECT_GT(wide_divisors, 5000);
}
#endif

// ------------------------------------------------------------------ the capture

// The capture changes no field of the ordinary payload (every field, all six counters): cap 0
// (no search), an infeasible positive cap, caps below, at and above the maximum, uncapped;
// never-reached sentinels and the final failed tier included. Every captured state also passes
// the state checks, whose reachability test must then say "reachable" for each state stopped
// below the maximum.
TEST(FlowCoreDeletion, CaptureLeavesTheFromZeroPayloadUnchanged) {
  std::mt19937_64 rng(5);
  int below_maximum = 0;
  for (int trial = 0; trial < 250; ++trial) {
    Case c = random_case(rng);
    const auto ref = independent_curve(*c.g, c.view, c.sources, c.destinations, c.el());
    const Amount fmax = ref.back().first;
    for (std::optional<Amount> cap : {std::optional<Amount> {}, std::optional<Amount> {0},
                                      std::optional<Amount> {fmax / 2}, std::optional<Amount> {fmax},
                                      std::optional<Amount> {fmax + 3}}) {
      c.cap = cap;
      const auto plain = place_global_units(*c.g, c.view, c.sources, c.destinations, c.el(), cap);
      FlowCoreResult with;
      const auto req = capture_of(c, &with);
      ASSERT_TRUE(same_payload(plain, with)) << "trial " << trial;
      ASSERT_TRUE(req.state.has_value());
      EXPECT_EQ(req.state->maximum(), !cap.has_value() || *cap >= fmax || fmax == 0) << "trial " << trial;
      check_state(c, *req.state, with, "trial " + std::to_string(trial));
      if (::testing::Test::HasFailure()) return;
      below_maximum += req.state->maximum() ? 0 : 1;
    }
  }
  // An infeasible positive cap: no path, the failed first tier.
  Case none = make_case(3, {{0, 1, 2, 1}}, {0}, {2});
  none.cap = 5;
  const auto plain = place_global_units(*none.g, none.view, none.sources, none.destinations, nullptr, none.cap);
  FlowCoreResult with;
  const auto req = capture_of(none, &with);
  EXPECT_TRUE(same_payload(plain, with));
  EXPECT_EQ(with.placed, 0);
  EXPECT_EQ(with.counters.dijkstra_runs, 1);
  EXPECT_TRUE(req.state->maximum());
  EXPECT_GT(below_maximum, 100);
}

// The shared equality helpers of test_utils.hpp compare every field: a copy with any one field
// changed is unequal, a vector field both by one more entry and by its last entry's value (so a
// helper comparing sizes only, or a prefix only, fails).
TEST(TestHelpers, EqualityHelpersCompareEveryField) {
  Case c = make_case(5, {{0, 1, 3, 3}, {0, 2, 1, 3}, {1, 2, 2, 7}, {2, 3, 1, 3}, {3, 1, 1, 2}}, {0}, {3});
  FlowCoreResult payload;
  const auto req = capture_of(c, &payload);
  ASSERT_TRUE(req.state.has_value());
  const auto& st = *req.state;
  ASSERT_GT(st.value(), 0);
  const std::vector<std::function<void(FlowCoreCounters&)>> counters = {
      [](auto& x) { ++x.dijkstra_runs; }, [](auto& x) { ++x.dijkstra_pops; }, [](auto& x) { ++x.settled; },
      [](auto& x) { ++x.arc_scans; },     [](auto& x) { ++x.dinic_phases; },  [](auto& x) { ++x.augmentations; },
  };
  for (std::size_t k = 0; k < counters.size(); ++k) {
    auto x = payload.counters;
    counters[k](x);
    EXPECT_FALSE(same_counters(x, payload.counters)) << "counter " << k;
  }
  const std::vector<std::function<void(FlowCoreResult&)>> payloads = {
      [](auto& x) { ++x.placed; },
      [](auto& x) { x.cost_distribution.emplace_back(1, 1); },
      [](auto& x) { ++x.cost_distribution.back().first; },
      [](auto& x) { ++x.cost_distribution.back().second; },
      [](auto& x) { x.edge_flows.push_back(0); },
      [](auto& x) { ++x.edge_flows.back(); },
      [](auto& x) { x.min_cut.reset(); },
      [](auto& x) { x.min_cut->push_back(0); },
      [](auto& x) { ++x.min_cut->back(); },
      [](auto& x) { x.potentials.push_back(0); },
      [](auto& x) { ++x.potentials.back(); },
      [](auto& x) { ++x.counters.settled; },
  };
  ASSERT_TRUE(payload.min_cut.has_value() && !payload.min_cut->empty());
  ASSERT_FALSE(payload.cost_distribution.empty());
  for (std::size_t k = 0; k < payloads.size(); ++k) {
    auto x = payload;
    payloads[k](x);
    EXPECT_FALSE(same_payload(x, payload)) << "payload field " << k;
  }
  EdgeId carrying = 0;
  while (st.edge_flows()[static_cast<std::size_t>(carrying)] == 0) ++carrying;
  const auto result = DeletionResolver(st).resolve(carrying);
  const std::vector<std::function<void(DeletionResult&)>> results = {
      [](auto& x) { ++x.value; },
      [](auto& x) { x.potentials.push_back(I128 {0}); },
      [](auto& x) { x.potentials.back() += I128 {1}; },
      [](auto& x) { x.cut.push_back(0); },
      [](auto& x) { ++x.cut.back(); },
      [](auto& x) { x.edge_flows.push_back(0); },
      [](auto& x) { ++x.edge_flows.back(); },
      [](auto& x) { ++x.counters.augmentations; },
  };
  ASSERT_FALSE(result.cut.empty());
  for (std::size_t k = 0; k < results.size(); ++k) {
    auto x = result;
    results[k](x);
    EXPECT_FALSE(same_result(x, result)) << "result field " << k;
  }
  const Network& net = st.network();
  const std::vector<std::function<void(Network&)>> networks = {
      [](auto& x) { ++x.num_nodes; },          [](auto& x) { x.edges.push_back(0); },
      [](auto& x) { x.tail.push_back(0); },    [](auto& x) { x.head.push_back(0); },
      [](auto& x) { x.cost.push_back(0); },    [](auto& x) { x.cap.push_back(0); },
      [](auto& x) { x.sources.push_back(0); }, [](auto& x) { x.destinations.push_back(0); },
      [](auto& x) { x.off.push_back(0); },     [](auto& x) { x.arc_head.push_back(0); },
      [](auto& x) { x.rev.push_back(0); },     [](auto& x) { x.arc_cost.push_back(0); },
      [](auto& x) { x.edge_slot.push_back(0); }, [](auto& x) { x.src_rev.push_back(0); },
      [](auto& x) { x.dst_fwd.push_back(0); },
      [](auto& x) { ++x.edges.back(); },        [](auto& x) { ++x.tail.back(); },
      [](auto& x) { ++x.head.back(); },         [](auto& x) { ++x.cost.back(); },
      [](auto& x) { ++x.cap.back(); },          [](auto& x) { ++x.sources.back(); },
      [](auto& x) { ++x.destinations.back(); }, [](auto& x) { ++x.off.back(); },
      [](auto& x) { ++x.arc_head.back(); },     [](auto& x) { ++x.rev.back(); },
      [](auto& x) { ++x.arc_cost.back(); },     [](auto& x) { ++x.edge_slot.back(); },
      [](auto& x) { ++x.src_rev.back(); },      [](auto& x) { ++x.dst_fwd.back(); },
  };
  for (std::size_t k = 0; k < networks.size(); ++k) {
    Network x = net;
    networks[k](x);
    EXPECT_FALSE(same_network(x, net)) << "network field " << k;
  }
  EXPECT_TRUE(same_network(net, net));
}

// check_residual is not vacuous: a captured state passes it, and each mutation below breaks the
// clause it names.
TEST(FlowCoreDeletion, ResidualCheckRejectsEachBrokenClause) {
  Case c = make_case(5, {{0, 1, 3, 3}, {0, 2, 1, 3}, {1, 2, 2, 7}, {2, 3, 1, 3}, {3, 1, 1, 2}}, {0}, {3});
  const auto req = capture_of(c);
  ASSERT_TRUE(req.state.has_value());
  const auto& st = *req.state;
  const Network& net = st.network();
  const std::vector<Amount> res(st.residual().begin(), st.residual().end());
  const std::vector<I128> pi(st.potentials().begin(), st.potentials().end());
  const auto want =
      divergences(net, {{net.num_nodes, st.value()}, {net.num_nodes + 1, -static_cast<Int128>(st.value())}});
  const Int128 bound = (static_cast<Int128>(net.num_nodes) + 3) * cost_sum(*c.g);  // (n + 1) C, n = num_nodes + 2
  auto check = [&](const std::vector<Amount>& r, const std::vector<I128>& p, std::size_t deleted = kNone) {
    return check_residual(net, r, p, deleted, want, bound);
  };
  ASSERT_TRUE(check(res, pi).empty());
  ASSERT_GT(st.value(), 0);
  // An edge that carries flow: its backward arc is positive.
  std::size_t k = 0;
  while (k < net.edges.size() && res[static_cast<std::size_t>(net.rev[static_cast<std::size_t>(net.edge_slot[k])])] == 0) {
    ++k;
  }
  ASSERT_LT(k, net.edges.size());
  const auto fwd = static_cast<std::size_t>(net.edge_slot[k]);
  const auto bwd = static_cast<std::size_t>(net.rev[fwd]);
  // One more unit on the forward arc: the pair exceeds the edge's capacity (the divergences read
  // the backward arc only).
  auto over = res;
  ++over[fwd];
  EXPECT_TRUE(has_failure(check(over, pi), "capacity bounds at edge " + std::to_string(net.edges[k])));
  // One unit from the backward to the forward arc: the pair's sum holds, the edge carries one unit
  // less, so its endpoints' divergences are off.
  auto less = res;
  --less[bwd];
  ++less[fwd];
  EXPECT_TRUE(has_failure(check(less, pi), "divergence vector"));
  // A source's virtual pair no longer sums to kMaxFlowCap.
  auto virt = res;
  ++virt[static_cast<std::size_t>(net.src_rev[0])];
  EXPECT_TRUE(has_failure(check(virt, pi), "capacity bounds at a virtual arc"));
  // The same residual read as if edge k were deleted: its arcs should be zero.
  EXPECT_TRUE(has_failure(check(res, pi, k), "the deleted edge keeps a residual arc"));
  // A positive real arc u -> v with pi(v) raised to pi(u) + cost + 1: reduced cost -1.
  const auto tail = slot_tails(net);
  std::size_t a = 0;
  while (a < res.size() && (res[a] <= 0 || tail[a] >= net.num_nodes || net.arc_head[a] >= net.num_nodes)) ++a;
  ASSERT_LT(a, res.size());
  auto raised = pi;
  raised[static_cast<std::size_t>(net.arc_head[a])] =
      raised[static_cast<std::size_t>(tail[a])] + I128 {net.arc_cost[a] + 1};
  EXPECT_TRUE(has_failure(check(res, raised), "negative reduced cost on slot " + std::to_string(a)));
  // A potential below 0 and one above the bound.
  auto low = pi;
  low[1] = I128 {-1};
  EXPECT_TRUE(has_failure(check(res, low), "potential outside [0, bound] at 1"));
  auto high = pi;
  high[1] = I128 {static_cast<std::int64_t>(bound) + 1};
  EXPECT_TRUE(has_failure(check(res, high), "potential outside [0, bound] at 1"));
  // Either array one entry too long (a zero residual entry, read by no other clause).
  auto longer_res = res;
  longer_res.push_back(0);
  EXPECT_EQ(check(longer_res, pi), std::vector<std::string> {"shape"});
  auto longer_pi = pi;
  longer_pi.push_back(I128 {0});
  EXPECT_EQ(check(res, longer_pi), std::vector<std::string> {"shape"});
  // Each arc of a pair negative while the pair's sum holds: edge k's, then the source's virtual
  // pair (the moved units also shift the divergences).
  auto negative = [&](std::size_t x) {
    auto r = res;
    const auto y = static_cast<std::size_t>(net.rev[x]);
    r[y] += r[x] + 1;
    r[x] = -1;
    return r;
  };
  EXPECT_TRUE(has_failure(check(negative(fwd), pi), "capacity bounds at edge " + std::to_string(net.edges[k])));
  EXPECT_TRUE(has_failure(check(negative(bwd), pi), "capacity bounds at edge " + std::to_string(net.edges[k])));
  const auto src_arc = static_cast<std::size_t>(net.src_rev[0]);
  EXPECT_TRUE(has_failure(check(negative(src_arc), pi), "capacity bounds at a virtual arc"));
  EXPECT_TRUE(has_failure(check(negative(static_cast<std::size_t>(net.rev[src_arc])), pi),
                          "capacity bounds at a virtual arc"));
  // Read as deleted, edge k with only its forward arc, then only its backward arc, positive.
  for (const bool forward : {true, false}) {
    auto kept = res;
    kept[fwd] = forward ? 1 : 0;
    kept[bwd] = forward ? 0 : 1;
    EXPECT_TRUE(has_failure(check(kept, pi, k), "the deleted edge keeps a residual arc")) << forward;
  }
  // A positive virtual arc with its head's potential raised by its cost + 1 over its tail's.
  std::size_t b = 0;
  while (b < res.size() && (res[b] <= 0 || (tail[b] < net.num_nodes && net.arc_head[b] < net.num_nodes))) ++b;
  ASSERT_LT(b, res.size());
  auto raised_virtual = pi;
  raised_virtual[static_cast<std::size_t>(net.arc_head[b])] =
      raised_virtual[static_cast<std::size_t>(tail[b])] + I128 {net.arc_cost[b] + 1};
  EXPECT_TRUE(
      has_failure(check(res, raised_virtual), "negative reduced cost on slot " + std::to_string(b) + " (virtual)"));
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
namespace {

// Captures c with the independent iterative all-vertex update beside: after every search, each
// vertex by its distance or, unreached, by the search's largest finite distance, from 0.
DeletionCapture capture_iterative(const Case& c, std::vector<Int128>& pi, int& searches) {
  pi.assign(static_cast<std::size_t>(c.g->num_nodes()) + 2, 0);
  searches = 0;
  DeletionCapture req;
  req.graph = c.g;
  req.on_search = [&](std::span<const Cost> d) {
    Cost top = 0;
    for (Cost x : d) {
      if (x != std::numeric_limits<Cost>::max()) top = std::max(top, x);
    }
    for (std::size_t v = 0; v < d.size(); ++v) pi[v] += d[v] != std::numeric_limits<Cost>::max() ? d[v] : top;
    ++searches;
  };
  (void)place_global_units(*c.g, c.view, c.sources, c.destinations, c.el(), c.cap, &req);
  return req;
}

} // namespace
#endif

// The running-sum checkpoints equal the iterative all-vertex update exactly, real and virtual
// vertices, a final failed search, never-settled vertices and the no-search cap-zero case
// included; a common offset keeps every reduced cost valid and fails this convention check.
TEST(FlowCoreDeletion, CheckpointsEqualTheIterativeUpdate) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  std::mt19937_64 rng(17);
  int several = 0;
  for (int trial = 0; trial < 300; ++trial) {
    Case c = random_case(rng);
    const Amount fmax = independent_curve(*c.g, c.view, c.sources, c.destinations, c.el()).back().first;
    for (std::optional<Amount> cap : {std::optional<Amount> {}, std::optional<Amount> {0},
                                      std::optional<Amount> {fmax / 2}, std::optional<Amount> {fmax + 1}}) {
      c.cap = cap;
      std::vector<Int128> pi;
      int searches = 0;
      const auto req = capture_iterative(c, pi, searches);
      const auto& st = *req.state;
      ASSERT_EQ(st.potentials().size(), pi.size());
      for (std::size_t v = 0; v < pi.size(); ++v) ASSERT_TRUE(wide(st.potentials()[v]) == pi[v]) << trial << " " << v;
      EXPECT_EQ(searches, st.counters().dijkstra_runs);
      if (cap == 0) EXPECT_EQ(searches, 0);
      several += searches >= 2 ? 1 : 0;  // a running sum of more than one search
    }
  }
  EXPECT_GT(several, 150);
  // The common-offset mutant: valid reduced costs, the wrong convention.
  Case c = make_case(5, {{0, 1, 3, 3}, {0, 2, 1, 3}, {1, 2, 2, 7}, {2, 3, 1, 3}, {3, 1, 1, 2}}, {0}, {3});
  std::vector<Int128> pi;
  int searches = 0;
  const auto req = capture_iterative(c, pi, searches);
  const auto& st = *req.state;
  std::vector<I128> shifted(st.potentials().begin(), st.potentials().end());
  for (auto& p : shifted) p += 1000;
  const Network& net = st.network();
  const auto want =
      divergences(net, {{net.num_nodes, st.value()}, {net.num_nodes + 1, -static_cast<Int128>(st.value())}});
  EXPECT_TRUE(check_residual(net, st.residual(), shifted, kNone, want, cost_sum(*c.g) * 1000).empty());
  bool same = true;
  for (std::size_t v = 0; v < pi.size(); ++v) same = same && wide(shifted[v]) == pi[v];
  EXPECT_FALSE(same);
#else
  GTEST_SKIP() << "needs the private hooks (NETGRAPH_CORE_PRIVATE_HOOKS)";
#endif
}

// ------------------------------------------------------------------ random views

TEST(FlowCoreDeletion, EveryEdgeOfRandomViewsAgreesWithAFreshSolve) {
  std::mt19937_64 rng(3005);
  int deletions = 0, masked = 0, several = 0, scaled = 0, rerouted = 0, partial = 0, cut = 0, capped = 0;
  for (int trial = 0; trial < 700; ++trial) {
    Case c = random_case(rng);
    const Amount fmax = independent_curve(*c.g, c.view, c.sources, c.destinations, c.el()).back().first;
    // Uncapped, and positive caps at and above the independently established maximum.
    for (std::optional<Amount> cap : {std::optional<Amount> {}, std::optional<Amount> {std::max<Amount>(fmax, 1)},
                                      std::optional<Amount> {fmax + 2}}) {
      c.cap = cap;
      capped += cap.has_value() && *cap == fmax && fmax > 0 ? 1 : 0;
      check_every_edge(c, "trial " + std::to_string(trial));
      if (::testing::Test::HasFailure()) return;
    }
    const auto req = capture_of(c);
    const auto& st = *req.state;
    for (const EdgeId e : positive_edges(st)) {
      ++deletions;
      const auto result = DeletionResolver(st).resolve(e);
      const Amount flow = st.edge_flows()[static_cast<std::size_t>(e)];
      const Amount lost = st.value() - result.value;
      masked += c.elig.has_value() && c.elig->uid() != 0 ? 1 : 0;
      several += c.sources.size() > 1 || c.destinations.size() > 1 ? 1 : 0;
      scaled += c.g->max_cost() > 1000000 ? 1 : 0;
      rerouted += lost == 0 ? 1 : 0;
      partial += lost > 0 && lost < flow ? 1 : 0;
      cut += lost == flow ? 1 : 0;
    }
  }
  EXPECT_GT(deletions, 600);
  EXPECT_GT(masked, 30);
  EXPECT_GT(several, 100);
  EXPECT_GT(scaled, 20);
  EXPECT_GT(rerouted, 10);
  EXPECT_GT(partial, 2);
  EXPECT_GT(cut, 100);
  EXPECT_GT(capped, 100);
}

// Deliberately damaged answers must trip each of the final-result certificate's clauses.
TEST(FlowCoreDeletion, CertificateRejectsDamagedAnswers) {
  Case c = make_case(4, {{0, 1, 3, 1}, {1, 3, 3, 1}, {0, 2, 3, 3}, {2, 3, 3, 1},
                          {0, 3, 1, 9}, {2, 2, 1, 7}}, {0}, {3});
  const auto ids = id_of_input(*c.g);
  c.elig = c.g->eligibility({}, std::vector<EdgeId> {ids[5]});
  const auto captured = capture_of(c);
  const auto& st = *captured.state;
  const EdgeId deleted = ids[4];
  const auto good = DeletionResolver(st).resolve(deleted);
  ASSERT_EQ(good.value, 6);
  EXPECT_EQ(failures_of([&] { expect_deletion(c, st, deleted, good, "good"); }), 0);
  auto rejects = [&](const DeletionResult& bad, const std::vector<std::string>& clauses) {
    ::testing::TestPartResultArray results;
    {
      ::testing::ScopedFakeTestPartResultReporter capture(
          ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD, &results);
      expect_deletion(c, st, deleted, bad, "damaged");
    }
    for (const auto& clause : clauses) {
      bool found = false;
      for (int i = 0; i < results.size(); ++i) {
        found = found || std::string(results.GetTestPartResult(i).message()).find(clause) != std::string::npos;
      }
      EXPECT_TRUE(found) << "missing certificate failure: " << clause;
    }
  };
  auto bad = good;
  bad.edge_flows.pop_back();
  rejects(bad, {"r.edge_flows.size()"});
  bad = good;
  bad.edge_flows[static_cast<std::size_t>(deleted)] = 1;
  rejects(bad, {"r.edge_flows[x] >= 0", "r.edge_flows[static_cast<std::size_t>(e)]",
                "flow_cost(*c.g, r.edge_flows) == flow_cost(*c.g, fresh.edge_flows)",
                "flow_cost(*c.g, r.edge_flows) == ref.back().second"});
  bad = good;
  bad.edge_flows[static_cast<std::size_t>(ids[5])] = 1;
  rejects(bad, {"r.edge_flows[x]"});
  bad = good;
  ++bad.value;
  rejects(bad, {"fresh.placed", "ref.back().first", "cut_value"});
  bad = good;
  bad.cut.clear();
  rejects(bad, {"*fresh.min_cut", "closure_cut(net, res, k)"});
  bad = good;
  std::fill(bad.edge_flows.begin(), bad.edge_flows.end(), 0);
  rejects(bad, {"connected(net, res)", "final:"});
}

// ------------------------------------------------------------------ hand-built examples

// Repair 1 reroutes everything (m = f, r = 0): u -> w -> v replaces u -> v.
TEST(FlowCoreDeletion, FullRepairOneReroute) {
  const Case c = make_case(5, {{0, 1, 1, 1}, {1, 2, 1, 1}, {1, 3, 1, 5}, {3, 2, 1, 5}, {2, 4, 1, 1}}, {0}, {4});
  check_every_edge(c, "full reroute");
  const auto req = capture_of(c);
  const EdgeId uv = id_of_input(*c.g)[1];
  const auto r = DeletionResolver(*req.state).resolve(uv);
  EXPECT_EQ(r.value, 1);
  EXPECT_TRUE(flow_cost(*c.g, r.edge_flows) == 12);
  EXPECT_EQ(r.counters.dijkstra_runs, 1);
}

// No reroute (m = 0, r = 1): deleting the sole unit-capacity source-to-destination edge takes
// exactly two Dijkstra runs.
TEST(FlowCoreDeletion, SoleEdgeTakesTwoSearches) {
  const Case c = make_case(2, {{0, 1, 1, 1}}, {0}, {1});
  check_every_edge(c, "sole edge");
  const auto req = capture_of(c);
  DeletionResolver rs(*req.state);
  const auto r = rs.resolve(0);
  EXPECT_EQ(r.value, 0);
  EXPECT_EQ(r.counters.dijkstra_runs, 2);
  // The run bound of check_deletion (m + 2 r = 2 here) holds, and each of its two clauses reports
  // a result that breaks it.
  EXPECT_EQ(run_bound_failure(*req.state, 1, r), "");
  DeletionResult gained = r;
  gained.value = req.state->value() + 1;
  EXPECT_EQ(run_bound_failure(*req.state, 1, gained), "the re-solve gained value");
  DeletionResult slower = r;
  ++slower.counters.dijkstra_runs;
  EXPECT_EQ(run_bound_failure(*req.state, 1, slower), "more Dijkstra runs than m + 2 r");
}

// The same deletion through the probe, with the control that reruns repair 1's failed search
// instead of reusing it: three searches, failing the m + 2 r = 2 bound, with the same flow,
// potentials, phases and paths; the reuse skips exactly that search's work.
TEST(FlowCoreDeletion, RerunningTheFailedSearchAddsExactlyOneSearch) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  const Case c = make_case(2, {{0, 1, 1, 1}}, {0}, {1});
  const auto req = capture_of(c);
  DeletionResolver rs(*req.state);
  std::vector<std::string> fails;
  std::vector<std::pair<int, FlowCoreCounters>> reuse_updates;
  std::vector<std::pair<int, FlowCoreCounters>> rerun_updates;
  const auto a = observed(rs, *req.state, 0, fails, false, &reuse_updates);
  const auto b = observed(rs, *req.state, 0, fails, true, &rerun_updates);
  for (const auto& s : fails) ADD_FAILURE() << s;
  EXPECT_EQ(a.counters.dijkstra_runs, 2);
  EXPECT_EQ(b.counters.dijkstra_runs, 3);
  // m + 2 r = 2: the bound check of check_deletion accepts the reuse and rejects the control.
  EXPECT_EQ(run_bound_failure(*req.state, 1, a), "");
  EXPECT_EQ(run_bound_failure(*req.state, 1, b), "more Dijkstra runs than m + 2 r");
  EXPECT_EQ(a.edge_flows, b.edge_flows);
  EXPECT_EQ(a.potentials, b.potentials);
  EXPECT_EQ(a.cut, b.cut);
  EXPECT_EQ(a.counters.dinic_phases, b.counters.dinic_phases);
  EXPECT_EQ(a.counters.augmentations, b.counters.augmentations);
  // The control's extra search: the counters between repair 1's last update and repair 2's
  // first (no Dinic work between them).
  ASSERT_EQ(rerun_updates.size(), reuse_updates.size() + 1);
  const auto first_of_2 =
      std::find_if(rerun_updates.begin(), rerun_updates.end(), [](const auto& u) { return u.first == 2; });
  ASSERT_NE(first_of_2, rerun_updates.end());
  ASSERT_NE(first_of_2, rerun_updates.begin());
  const auto& before = std::prev(first_of_2)->second;
  const auto& after = first_of_2->second;
  EXPECT_EQ(after.dinic_phases, before.dinic_phases);
  EXPECT_EQ(after.augmentations, before.augmentations);
  EXPECT_EQ(b.counters.dijkstra_runs - (after.dijkstra_runs - before.dijkstra_runs), a.counters.dijkstra_runs);
  EXPECT_EQ(b.counters.dijkstra_pops - (after.dijkstra_pops - before.dijkstra_pops), a.counters.dijkstra_pops);
  EXPECT_EQ(b.counters.settled - (after.settled - before.settled), a.counters.settled);
  EXPECT_EQ(b.counters.arc_scans - (after.arc_scans - before.arc_scans), a.counters.arc_scans);
  EXPECT_GT(after.settled - before.settled, 0);
#else
  GTEST_SKIP() << "needs the private hooks (NETGRAPH_CORE_PRIVATE_HOOKS)";
#endif
}

// Repair 1 moves part (m = 1 of f = 2); the remainder 1 goes back to s* and then from t* to v.
TEST(FlowCoreDeletion, PartialRerouteThenBothRemainderPushes) {
  const Case c = make_case(5, {{0, 1, 2, 1}, {1, 2, 2, 1}, {1, 3, 1, 2}, {3, 2, 1, 2}, {2, 4, 2, 1}}, {0}, {4});
  check_every_edge(c, "partial reroute");
  const auto req = capture_of(c);
  const EdgeId uv = id_of_input(*c.g)[1];
  ASSERT_EQ(req.state->edge_flows()[static_cast<std::size_t>(uv)], 2);
  const auto r = DeletionResolver(*req.state).resolve(uv);
  EXPECT_EQ(r.value, 1);  // r = 1, m = 1
  EXPECT_TRUE(flow_cost(*c.g, r.edge_flows) == 6);
  EXPECT_LE(r.counters.dijkstra_runs, 1 + 2 * 1);
}

// Every cost times a factor: each re-solve of the scaled graph passes the independent checks, and
// against the unscaled one its working potentials scale while its value, flows, cut and operation
// counters do not (every search meets the same keys in the same order).
TEST(FlowCoreDeletion, ScalingEveryCostScalesOnlyThePotentials) {
  const std::vector<InEdge> edges = {{0, 1, 2, 1}, {1, 2, 2, 1}, {1, 3, 1, 2}, {3, 2, 1, 2}, {2, 4, 2, 1}, {0, 3, 1, 5}};
  const Case base = make_case(5, edges, {0}, {4});
  const auto base_req = capture_of(base);
  int compared = 0;
  for (Cost factor : {3, 1000003}) {
    std::vector<InEdge> scaled_edges = edges;
    for (auto& e : scaled_edges) e.cost *= factor;
    const Case scaled = make_case(5, scaled_edges, {0}, {4});
    check_every_edge(scaled, "factor " + std::to_string(factor));
    const auto req = capture_of(scaled);
    for (EdgeId e : base_req.state->network().edges) {
      if (base_req.state->edge_flows()[static_cast<std::size_t>(e)] == 0) continue;
      const auto want = DeletionResolver(*base_req.state).resolve(e);
      const auto got = DeletionResolver(*req.state).resolve(e);
      const std::string where = "factor " + std::to_string(factor) + " edge " + std::to_string(e);
      EXPECT_EQ(got.value, want.value) << where;
      EXPECT_EQ(got.edge_flows, want.edge_flows) << where;
      EXPECT_EQ(got.cut, want.cut) << where;
      ASSERT_EQ(got.potentials.size(), want.potentials.size()) << where;
      for (std::size_t v = 0; v < want.potentials.size(); ++v) {
        EXPECT_TRUE(wide(got.potentials[v]) == wide(want.potentials[v]) * factor) << where << " vertex " << v;
      }
      EXPECT_TRUE(same_counters(got.counters, want.counters)) << where;
      ++compared;
    }
  }
  // Both units take 0 -> 1 -> 2 -> 4 (cost 3; 2 -> 4 caps the value at 2): three positive-flow
  // edges, for each factor.
  EXPECT_EQ(compared, 6);
}

// Total loss: a bridge carrying every unit.
TEST(FlowCoreDeletion, TotalLoss) {
  const Case c = make_case(3, {{0, 1, 3, 1}, {1, 2, 2, 1}, {1, 2, 1, 3}}, {0}, {2});
  check_every_edge(c, "total loss");
  const auto req = capture_of(c);
  const auto r = DeletionResolver(*req.state).resolve(id_of_input(*c.g)[0]);
  EXPECT_EQ(r.value, 0);
  EXPECT_TRUE(r.cut.empty());
}

// Zero-flow deletions with positive view capacity (the deletion changes the canonical cut, not
// the value) and with zero view capacity (an arc crossing the closure without capacity).
TEST(FlowCoreDeletion, ZeroFlowDeletions) {
  const Case positive = make_case(4, {{0, 1, 1, 1}, {1, 2, 1, 1}, {2, 3, 1, 1}, {0, 2, 1, 3}}, {0}, {3});
  check_every_edge(positive, "zero flow, positive capacity");
  const EdgeId spare = id_of_input(*positive.g)[3];
  const auto req = capture_of(positive);
  EXPECT_EQ(req.state->edge_flows()[static_cast<std::size_t>(spare)], 0);
  FlowCoreResult base;
  (void)capture_of(positive, &base);
  const auto fresh = fresh_without(positive, spare);
  ASSERT_TRUE(base.min_cut.has_value() && fresh.min_cut.has_value());
  EXPECT_NE(*base.min_cut, *fresh.min_cut);
  const Case zero = make_case(2, {{0, 1, 0, 1}}, {0}, {1});
  check_every_edge(zero, "zero view capacity");
  const auto zreq = capture_of(zero);
  EXPECT_TRUE(zreq.state->maximum());
  EXPECT_EQ(zreq.state->value(), 0);
}

// Parallel equal-cost ties: the re-solve's flow differs from the fresh solve's, its value, cost
// and cut do not.
TEST(FlowCoreDeletion, TiedFlowsDifferValueCostAndCutAgree) {
  const Case c = make_case(3, {{0, 1, 1, 2}, {0, 1, 2, 3}, {1, 2, 1, 1}, {1, 2, 1, 1}}, {0}, {2});
  check_every_edge(c, "ties");
  const auto req = capture_of(c);
  const EdgeId dear = id_of_input(*c.g)[1];
  const auto r = DeletionResolver(*req.state).resolve(dear);
  const auto fresh = fresh_without(c, dear);
  EXPECT_EQ(r.edge_flows, (std::vector<Amount> {1, 0, 0, 1}));
  EXPECT_EQ(fresh.edge_flows, (std::vector<Amount> {1, 0, 1, 0}));
  EXPECT_EQ(r.value, 1);
  EXPECT_TRUE(flow_cost(*c.g, r.edge_flows) == 3);
  ASSERT_TRUE(fresh.min_cut.has_value());
  EXPECT_EQ(r.cut, *fresh.min_cut);
  EXPECT_EQ(r.cut, (std::vector<EdgeId> {0}));
}

// Several sources and destinations; edges out of a source, into a destination and between
// endpoints; masked endpoints; vertices no search reaches.
TEST(FlowCoreDeletion, EndpointsMasksAndUnreachedVertices) {
  Case c = make_case(8,
                     {{0, 2, 2, 1}, {1, 2, 1, 2}, {2, 4, 2, 1}, {2, 5, 2, 3}, {1, 5, 1, 7}, {0, 4, 1, 9},
                      {4, 5, 1, 1}, {3, 6, 4, 1}, {6, 7, 4, 1}, {5, 2, 1, 1}},
                     {0, 1, 3}, {4, 5});
  // Each configuration re-solves at least one deletion (its flow is positive).
  int resolved = 0;
  check_every_edge(c, "endpoints", &resolved);
  EXPECT_GT(resolved, 0);
  c.elig = c.g->eligibility(std::vector<NodeId> {1}, std::vector<EdgeId> {});
  resolved = 0;
  check_every_edge(c, "masked source", &resolved);
  EXPECT_GT(resolved, 0);
  c.elig = c.g->eligibility(std::vector<NodeId> {4}, std::vector<EdgeId> {id_of_input(*c.g)[3]});
  resolved = 0;
  check_every_edge(c, "masked destination and edge", &resolved);
  EXPECT_GT(resolved, 0);
}

// A near-domain cost witness: H = 2^61 and a cost sum H + 11
// below 2^62; the working potentials exceed int64 after the deletion, vertex 4 is never
// reached (the published sentinel, a numeric working potential).
TEST(FlowCoreDeletion, NearDomainCostsNeedWideWorkingPotentials) {
  const Cost H = Cost {1} << 61;
  const Case c = make_case(5, {{0, 1, 3, 3}, {0, 2, 1, 3}, {1, 2, 2, H}, {2, 3, 1, 3}, {3, 1, 1, 2}}, {0}, {3});
  ASSERT_TRUE(cost_sum(*c.g) == static_cast<Int128>(H) + 11);
  FlowCoreResult payload;
  const auto req = capture_of(c, &payload);
  const auto& st = *req.state;
  EXPECT_EQ(payload.potentials[4], kNeverReached);
  const std::vector<Int128> base {0, 3, static_cast<Int128>(H) + 3, static_cast<Int128>(H) + 6,
                                    static_cast<Int128>(H) + 6, 0, static_cast<Int128>(H) + 6};
  for (std::size_t v = 0; v < base.size(); ++v) EXPECT_TRUE(wide(st.potentials()[v]) == base[v]) << v;
  const EdgeId e = id_of_input(*c.g)[3];
  check_deletion(c, st, e, "near domain");
  const auto r = DeletionResolver(st).resolve(e);
  const Int128 h = H;
  const std::vector<Int128> after {3 * h + 5, 2 * h + 8, 3 * h + 8, 2 * h + 6, 4 * h + 11, 3 * h + 5, 2 * h + 6};
  for (std::size_t v = 0; v < after.size(); ++v) EXPECT_TRUE(wide(r.potentials[v]) == after[v]) << v;
  EXPECT_TRUE(wide(r.potentials[4]) > std::numeric_limits<std::int64_t>::max());
}

// The precondition: a state stopped by a cap below the maximum with a residual path left is
// rejected before any change; maximum states (uncapped, cap at the maximum) are accepted.
TEST(FlowCoreDeletion, OnlyMaximumStatesAreAccepted) {
  Case c = make_case(2, {{0, 1, 1, 1}, {0, 1, 3, 2}}, {0}, {1});
  c.cap = 1;
  const auto req = capture_of(c);
  ASSERT_FALSE(req.state->maximum());
  DeletionResolver rs(*req.state);
  const EdgeId cheap = id_of_input(*c.g)[0];
  // Rejected before any repair (the resolver works on its own copy of the state).
  EXPECT_EQ(error_of([&] { (void)rs.resolve(cheap); }),
            "deletion re-solve: the state is not certified maximum (a residual path from the sources to the "
            "destinations remains)");
  for (std::optional<Amount> cap : {std::optional<Amount> {}, std::optional<Amount> {4}}) {
    c.cap = cap;
    const auto ok = capture_of(c);
    ASSERT_TRUE(ok.state->maximum());
    EXPECT_EQ(DeletionResolver(*ok.state).resolve(cheap).value, 3);
  }
  // Not an edge of the view (excluded by the Eligibility), and an empty state (an edge of the
  // view that carries no flow: check_every_edge).
  c.cap = std::nullopt;
  c.elig = c.g->eligibility(std::vector<NodeId> {}, std::vector<EdgeId> {cheap});
  const auto masked = capture_of(c);
  DeletionResolver reused(*masked.state);
  EXPECT_EQ(error_of([&] { (void)reused.resolve(cheap); }),
            "deletion re-solve: edge 0 is not in the solve's view");
  const EdgeId allowed = id_of_input(*c.g)[1];
  EXPECT_TRUE(same_result(reused.resolve(allowed), DeletionResolver(*masked.state).resolve(allowed)));
  EXPECT_EQ(error_of([] { DeletionResolver {DeletionState {}}; }),
            "deletion re-solve: an empty state (no solve captured it)");
}

// ------------------------------------------------------------------ the saved state

// A saved state survives another solve on the same thread and serves candidates A, B, A: it
// never changes and both A results (flow, certificates, counters) agree, on 1, 2 and 4 threads
// each with its own resolver over one shared state.
TEST(FlowCoreDeletion, ASavedStateServesCandidatesIndependently) {
  std::mt19937_64 rng(29);
  int checked = 0;
  for (int trial = 0; trial < 400 && checked < 20; ++trial) {
    const Case c = random_case(rng);
    const auto req = capture_of(c);
    const DeletionState& st = *req.state;
    const std::vector<EdgeId> positive = positive_edges(st);
    if (positive.size() < 2) continue;
    ++checked;
    const std::vector<Amount> res0(st.residual().begin(), st.residual().end());
    const std::vector<I128> pi0(st.potentials().begin(), st.potentials().end());
    const DeletionState copy = st;  // a copy shares the network only
    const Case other = random_case(rng);
    (void)place_global_units(*other.g, other.view, other.sources, other.destinations, other.el(), std::nullopt);
    const EdgeId A = positive.front();
    const EdgeId B = positive.back();
    DeletionResolver one(copy);
    const auto a1 = one.resolve(A);
    const auto b1 = one.resolve(B);
    const auto a2 = one.resolve(A);
    EXPECT_TRUE(same_result(a1, a2)) << trial;
    EXPECT_EQ(b1.value, DeletionResolver(st).resolve(B).value) << trial;
    for (int threads : {2, 4}) {
      std::vector<DeletionResolver> resolvers;
      for (int t = 0; t < threads; ++t) resolvers.emplace_back(st);
      std::vector<std::vector<DeletionResult>> got(static_cast<std::size_t>(threads));
      std::vector<std::thread> pool;
      for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t] {
          auto& rs = resolvers[static_cast<std::size_t>(t)];
          auto& out = got[static_cast<std::size_t>(t)];
          out.push_back(rs.resolve(A));
          out.push_back(rs.resolve(B));
          out.push_back(rs.resolve(A));
        });
      }
      for (auto& th : pool) th.join();
      for (const auto& out : got) {
        EXPECT_TRUE(same_result(out[0], a1) && same_result(out[2], a1) && same_result(out[1], b1))
            << trial << " threads " << threads;
      }
    }
    EXPECT_EQ(std::vector<Amount>(st.residual().begin(), st.residual().end()), res0) << trial;
    EXPECT_EQ(std::vector<I128>(st.potentials().begin(), st.potentials().end()), pi0) << trial;
    EXPECT_EQ(std::vector<Amount>(copy.residual().begin(), copy.residual().end()), res0) << trial;
    EXPECT_EQ(&copy.network(), &st.network());
  }
  EXPECT_EQ(checked, 20);
}

namespace {

// Whether a re-solve of e agrees with a solve from zero on the view without e: the value, the
// cost and the minimum cut (independent of the resolver).
bool matches_fresh_solve(const Case& c, EdgeId e, const DeletionResult& r) {
  const auto fresh = fresh_without(c, e);
  return r.value == fresh.placed && flow_cost(*c.g, r.edge_flows) == flow_cost(*c.g, fresh.edge_flows) &&
         fresh.min_cut.has_value() && r.cut == *fresh.min_cut;
}

} // namespace

// A resolver reads its state's network in place, so the network must live as
// long as the resolver by the resolver's own copy of the state. Built from a copy that then goes
// away with the capture request, moved by construction, by the reallocations of a growing vector
// and by an assignment over a resolver of another state (which drops that resolver's share of the
// other network; the network goes with the last of grow's resolvers), and after a solve of another
// network on the same thread, it serves every candidate twice, in reverse and in order, exactly as
// a resolver that never moved, whose results agree with solves from zero. AddressSanitizer
// reports a read of a released network.
TEST(FlowCoreDeletion, AResolverKeepsItsNetworkAliveAcrossMoves) {
  std::mt19937_64 rng(31);
  int checked = 0;
  for (int trial = 0; trial < 400 && checked < 20; ++trial) {
    const Case c = random_case(rng);
    std::vector<EdgeId> positive;
    std::vector<DeletionResult> want;
    std::optional<DeletionResolver> kept;
    {
      auto req = capture_of(c);
      ASSERT_TRUE(req.state->maximum()) << "trial " << trial;  // random_case sets no cap
      positive = positive_edges(*req.state);
      if (positive.empty()) continue;
      {
        DeletionResolver still(*req.state);
        for (EdgeId e : positive) want.push_back(still.resolve(e));
      }
      for (std::size_t i = 0; i < positive.size(); ++i) {
        EXPECT_TRUE(matches_fresh_solve(c, positive[i], want[i])) << "trial " << trial << " edge " << positive[i];
      }
      auto copy = std::make_unique<DeletionState>(*req.state);
      req.state = std::nullopt;
      kept.emplace(*copy);
    }  // the copy and the request are gone: the resolver's state is the network's only owner
    ++checked;
    const Case other = random_case(rng);
    auto other_req = capture_of(other);  // rebuilds this thread's solver for another network
    std::vector<DeletionResolver> grow;
    grow.emplace_back(std::move(*kept));
    kept = std::nullopt;
    for (int i = 0; i < 6; ++i) grow.emplace_back(*other_req.state);  // reallocations move grow[0]
    DeletionResolver target(*other_req.state);
    other_req.state = std::nullopt;
    target = std::move(grow.front());
    grow.clear();  // the other network's last owners go
    for (std::size_t i = positive.size(); i-- > 0;) {
      EXPECT_TRUE(same_result(target.resolve(positive[i]), want[i])) << "trial " << trial << " edge " << positive[i];
    }
    for (std::size_t i = 0; i < positive.size(); ++i) {
      EXPECT_TRUE(same_result(target.resolve(positive[i]), want[i])) << "trial " << trial << " edge " << positive[i];
    }
  }
  EXPECT_EQ(checked, 20);
}

// Resolvers of one state on 1, 2 and 4 threads read the one network the state
// copies share and never write it; each keeps its own residual, potentials and scratch. Each
// thread runs every candidate in its own rotation of the order, twice; every result equals the
// single-thread one (which agrees with a solve from zero), and the network equals its deep copy
// taken before, field by field.
// ThreadSanitizer checks the concurrent reads.
TEST(FlowCoreDeletion, WorkersReadOneNetworkAndNeverWriteIt) {
  std::mt19937_64 rng(37);
  int checked = 0;
  for (int trial = 0; trial < 400 && checked < 15; ++trial) {
    const Case c = random_case(rng);
    const auto req = capture_of(c);
    const DeletionState& st = *req.state;
    ASSERT_TRUE(st.maximum()) << "trial " << trial;  // random_case sets no cap
    const auto positive = positive_edges(st);
    if (positive.size() < 2) continue;
    ++checked;
    const Network before = st.network();
    std::vector<DeletionResult> want;
    for (EdgeId e : positive) {
      want.push_back(DeletionResolver(st).resolve(e));
      EXPECT_TRUE(matches_fresh_solve(c, e, want.back())) << trial << " edge " << e;
    }
    for (int threads : {1, 2, 4}) {
      std::vector<DeletionResolver> resolvers;
      for (int t = 0; t < threads; ++t) resolvers.emplace_back(st);
      std::vector<int> wrong(static_cast<std::size_t>(threads), 0);
      std::vector<std::thread> pool;
      for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t] {
          auto& rs = resolvers[static_cast<std::size_t>(t)];
          for (int round = 0; round < 2; ++round) {
            for (std::size_t j = 0; j < positive.size(); ++j) {
              const std::size_t i = (j + static_cast<std::size_t>(t)) % positive.size();
              if (!same_result(rs.resolve(positive[i]), want[i])) ++wrong[static_cast<std::size_t>(t)];
            }
          }
        });
      }
      for (auto& th : pool) th.join();
      for (int t = 0; t < threads; ++t) EXPECT_EQ(wrong[static_cast<std::size_t>(t)], 0) << trial << " thread " << t;
      EXPECT_TRUE(same_network(st.network(), before)) << trial << " threads " << threads;
    }
  }
  EXPECT_EQ(checked, 15);
}

// ------------------------------------------------------------------ sensitivity

// The exact model's loss is placed minus capacity() with the edge excluded, for every
// candidate (the baseline's positive-flow edges), on 1, 2 and 4 workers; zero-flow edges are
// not candidates. Both maxima also equal the independent Bellman-Ford curve's.
TEST(Sensitivity, GlobalLossesEqualCapacityWithTheEdgeExcluded) {
  std::mt19937_64 rng(41);
  int entries = 0, masked = 0, several = 0, scaled = 0, rerouted = 0, partial = 0, cut = 0;
  for (int trial = 0; trial < 300; ++trial) {
    const Case c = random_case(rng);
    const QueryOptions opts {global_config(), c.el()};
    const auto base = capacity(*c.g, c.sources, c.destinations, QueryOptions {global_config(), c.el(), nullptr,
                                                                                 OutputOptions {true, false, false}});
    const auto s1 = sensitivity(*c.g, c.sources, c.destinations, opts, 1);
    EXPECT_EQ(s1.baseline.placed, base.placed) << trial;
    const std::vector<Amount> caps(c.g->capacity_units().begin(), c.g->capacity_units().end());
    EXPECT_EQ(base.placed, independent_curve(*c.g, caps, c.sources, c.destinations, c.el()).back().first) << trial;
    std::vector<EdgeId> want;
    for (EdgeId e = 0; e < c.g->num_edges(); ++e) {
      if ((*base.edge_flows)[static_cast<std::size_t>(e)] > 0) want.push_back(e);
    }
    ASSERT_EQ(s1.entries.size(), want.size()) << trial;
    for (std::size_t i = 0; i < want.size(); ++i) {
      const auto& entry = s1.entries[i];
      EXPECT_EQ(entry.edge, want[i]) << trial;
      const EdgeId one[1] = {entry.edge};
      Eligibility without = c.g->eligibility({}, one);
      if (c.elig) without = Eligibility::intersect(*c.elig, without);
      const auto fresh = capacity(*c.g, c.sources, c.destinations, QueryOptions {global_config(), &without});
      EXPECT_EQ(entry.loss, base.placed - fresh.placed) << trial << " edge " << entry.edge;
      auto without_view = caps;
      without_view[static_cast<std::size_t>(entry.edge)] = 0;
      EXPECT_EQ(fresh.placed, independent_curve(*c.g, without_view, c.sources, c.destinations, c.el()).back().first)
          << trial << " edge " << entry.edge;
      ++entries;
      const Amount flow = (*base.edge_flows)[static_cast<std::size_t>(entry.edge)];
      masked += c.elig.has_value() && c.elig->uid() != 0 ? 1 : 0;
      several += c.sources.size() > 1 || c.destinations.size() > 1 ? 1 : 0;
      scaled += c.g->max_cost() > 1000000 ? 1 : 0;
      rerouted += entry.loss == 0 ? 1 : 0;
      partial += entry.loss > 0 && entry.loss < flow ? 1 : 0;
      cut += entry.loss == flow ? 1 : 0;
    }
    for (std::int32_t workers : {2, 4}) {
      const auto sw = sensitivity(*c.g, c.sources, c.destinations, opts, workers);
      ASSERT_EQ(sw.entries.size(), s1.entries.size());
      for (std::size_t i = 0; i < sw.entries.size(); ++i) {
        EXPECT_EQ(sw.entries[i].edge, s1.entries[i].edge);
        EXPECT_EQ(sw.entries[i].loss, s1.entries[i].loss);
      }
    }
  }
  EXPECT_GT(entries, 300);
  EXPECT_GT(masked, 20);
  EXPECT_GT(several, 70);
  EXPECT_GT(scaled, 10);
  EXPECT_GT(rerouted, 10);
  EXPECT_GT(partial, 1);
  EXPECT_GT(cut, 100);
}

// OutputOptions controls only the published baseline: every combination captures the same state
// (the same entries) and the baseline's optional fields are capacity()'s.
TEST(Sensitivity, GlobalBaselineKeepsEveryOutputsCombination) {
  const Case c = make_case(5, {{0, 1, 2, 1}, {1, 2, 2, 1}, {1, 3, 1, 2}, {3, 2, 1, 2}, {2, 4, 2, 1}, {0, 3, 1, 5}},
                           {0}, {4});
  const auto reference = sensitivity(*c.g, c.sources, c.destinations, QueryOptions {global_config()}, 1);
  ASSERT_FALSE(reference.entries.empty());
  for (int bits = 0; bits < 16; ++bits) {
    const OutputOptions outputs {(bits & 1) != 0, (bits & 2) != 0, (bits & 4) != 0, (bits & 8) != 0};
    QueryOptions opts {global_config()};
    opts.outputs = outputs;
    const auto s = sensitivity(*c.g, c.sources, c.destinations, opts, 1);
    const auto cap = capacity(*c.g, c.sources, c.destinations, QueryOptions {global_config(), nullptr, nullptr, outputs});
    EXPECT_EQ(s.baseline.placed, cap.placed) << bits;
    EXPECT_EQ(s.baseline.cost_distribution, cap.cost_distribution) << bits;
    EXPECT_EQ(s.baseline.edge_flows, cap.edge_flows) << bits;
    EXPECT_EQ(s.baseline.residual, cap.residual) << bits;
    EXPECT_EQ(s.baseline.bottleneck, cap.bottleneck) << bits;
    EXPECT_EQ(s.baseline.min_cut, cap.min_cut) << bits;
    // Each requested output is present, each other one absent; Global's curve always, and its
    // min_cut when requested (the uncapped baseline is maximal).
    EXPECT_EQ(s.baseline.edge_flows.has_value(), outputs.edge_flows) << bits;
    EXPECT_EQ(s.baseline.residual.has_value(), outputs.residual) << bits;
    EXPECT_EQ(s.baseline.bottleneck.has_value(), outputs.bottleneck) << bits;
    EXPECT_EQ(s.baseline.min_cut.has_value(), outputs.min_cut) << bits;
    EXPECT_TRUE(s.baseline.cost_distribution.has_value()) << bits;
    ASSERT_EQ(s.entries.size(), reference.entries.size()) << bits;
    for (std::size_t i = 0; i < s.entries.size(); ++i) {
      EXPECT_EQ(s.entries[i].edge, reference.entries[i].edge) << bits;
      EXPECT_EQ(s.entries[i].loss, reference.entries[i].loss) << bits;
    }
  }
}
