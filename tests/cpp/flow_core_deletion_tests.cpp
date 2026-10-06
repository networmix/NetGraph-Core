// The deletion re-solve of the flow core (APPROACH 2.2) and the exact sensitivity model that uses
// it (APPROACH 11). Every check is independent of the repair: a solve from zero on the deleted
// view, the Bellman-Ford curve of test_utils.hpp, residuals rebuilt from flows, divergences and
// reduced costs computed here in __int128.
#include <gtest/gtest.h>

#include <cstdint>
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

__int128 wide(I128 x) {
  return static_cast<__int128>((static_cast<unsigned __int128>(static_cast<std::uint64_t>(x.hi)) << 64) | x.lo);
}

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

__int128 cost_of(const Graph& g, const std::vector<Amount>& flows) {
  __int128 total = 0;
  for (std::size_t e = 0; e < flows.size(); ++e) total += static_cast<__int128>(g.cost()[e]) * flows[e];
  return total;
}

__int128 cost_sum(const Graph& g) {
  __int128 c = 0;
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
std::vector<__int128> divergences(const Network& net, std::initializer_list<std::pair<std::int32_t, __int128>> at) {
  std::vector<__int128> d(static_cast<std::size_t>(net.num_nodes) + 2, 0);
  for (auto [v, x] : at) d[static_cast<std::size_t>(v)] += x;
  return d;
}

// The failed checks of a residual with potentials on the network (APPROACH 2.2): capacity
// bounds (both arcs of a pair non-negative, summing to the pair's capacity; the deleted edge's
// pair zero), the divergence vector `want`, non-negative reduced cost on every positive
// residual arc (real and virtual), and the bound `bound` on every potential (non-negative).
std::vector<std::string> check_residual(const Network& net, std::span<const Amount> res, std::span<const I128> pi,
                                        std::size_t deleted, const std::vector<__int128>& want, __int128 bound) {
  std::vector<std::string> fails;
  const auto n = static_cast<std::size_t>(net.num_nodes) + 2;
  if (res.size() != net.arc_head.size() || pi.size() != n) return {"shape"};
  std::vector<__int128> div(n, 0);
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
        static_cast<__int128>(res[a]) + res[static_cast<std::size_t>(net.rev[a])] != kMaxFlowCap) {
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
    if (static_cast<__int128>(net.arc_cost[a]) + wide(pi[u]) - wide(pi[v]) < 0) {
      fails.push_back("negative reduced cost on slot " + std::to_string(a) + (u >= S || v >= S ? " (virtual)" : ""));
    }
  }
  for (std::size_t v = 0; v < n; ++v) {
    if (wide(pi[v]) < 0 || wide(pi[v]) > bound) fails.push_back("potential outside [0, bound] at " + std::to_string(v));
  }
  return fails;
}

// Whether the real sources reach t* over positive residual arcs.
bool connected(const Network& net, std::span<const Amount> res) {
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
  return seen[n - 1] != 0;
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
// The re-solve of e observed by the probe: every observed state (after the deletion, after every
// potential update and push of the three repairs, at completion) is checked against the
// capacity bounds, the exact divergence vector of 2.2 for the amounts moved so far and dual
// feasibility on every positive residual arc. `updates` collects the counters after every
// potential update with its repair.
DeletionResult observed(DeletionResolver& rs, const DeletionState& st, EdgeId e, std::vector<std::string>& fails,
                        bool rerun = false, std::vector<std::pair<int, FlowCoreCounters>>* updates = nullptr) {
  const Network& net = st.network();
  const std::size_t k = local_of(net, e);
  const std::int32_t u = net.tail[k];
  const std::int32_t v = net.head[k];
  const Amount F = st.value();
  const Amount f = st.edge_flows()[static_cast<std::size_t>(e)];
  const __int128 n = net.num_nodes + 2;
  const __int128 bound = (22 * n + 15) * cost_sum(st.graph());
  Amount moved[5] = {0, 0, 0, 0, 0};
  detail::DeletionProbe probe;
  probe.rerun_failed_search = rerun;
  probe.observe = [&](const detail::DeletionEvent& ev) {
    if (ev.repair >= 1 && ev.repair <= 3) moved[ev.repair] = ev.moved;
    const auto want = divergences(net, {{net.num_nodes, F - moved[2]},
                                        {net.num_nodes + 1, -static_cast<__int128>(F) + moved[3]},
                                        {u, -static_cast<__int128>(f) + moved[1] + moved[2]},
                                        {v, static_cast<__int128>(f) - moved[1] - moved[3]}});
    for (const auto& s : check_residual(net, ev.residual, ev.potentials, k, want, bound)) {
      fails.push_back("repair " + std::to_string(ev.repair) + (ev.push ? " push" : " update") + ": " + s);
    }
    if (updates != nullptr && !ev.push && ev.repair >= 1 && ev.repair <= 3) updates->emplace_back(ev.repair, *ev.counters);
  };
  return rs.resolve(e, probe);
}
#endif

// Every check of 2.2 for the deletion of positive-flow edge e from the captured state of c.
void check_deletion(const Case& c, const DeletionState& st, EdgeId e, const std::string& where) {
  const Network& net = st.network();
  const std::size_t k = local_of(net, e);
  ASSERT_NE(k, kNone) << where;
  const Amount f = st.edge_flows()[static_cast<std::size_t>(e)];
  ASSERT_GT(f, 0) << where;
  DeletionResolver rs(st);
  const auto r = rs.resolve(e);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  std::vector<std::string> fails;
  const auto ro = observed(rs, st, e, fails);
  for (const auto& s : fails) ADD_FAILURE() << where << ": " << s;
  EXPECT_EQ(ro.value, r.value) << where;
  EXPECT_EQ(ro.edge_flows, r.edge_flows) << where;
#endif
  // From zero on V \ {e}, and the independent Bellman-Ford curve.
  auto view = c.view;
  view[static_cast<std::size_t>(e)] = 0;
  const auto fresh = place_global_units(*c.g, view, c.sources, c.destinations, c.el(), std::nullopt);
  const auto ref = independent_curve(*c.g, view, c.sources, c.destinations, c.el());
  EXPECT_EQ(r.value, fresh.placed) << where;
  EXPECT_EQ(r.value, ref.back().first) << where;
  EXPECT_TRUE(cost_of(*c.g, r.edge_flows) == cost_of(*c.g, fresh.edge_flows)) << where;
  EXPECT_TRUE(cost_of(*c.g, r.edge_flows) == ref.back().second) << where;
  ASSERT_TRUE(fresh.min_cut.has_value()) << where;
  EXPECT_EQ(r.cut, *fresh.min_cut) << where;
  Amount cut_value = 0;
  for (EdgeId x : r.cut) cut_value += view[static_cast<std::size_t>(x)];
  EXPECT_EQ(cut_value, r.value) << where;
  EXPECT_EQ(r.edge_flows[static_cast<std::size_t>(e)], 0) << where;
  // The final residual from the flow: bounds, divergences (value at s*), working potentials
  // valid on every positive residual arc, the destinations disconnected.
  const auto res = residual_of(net, r.edge_flows, k);
  const __int128 n = net.num_nodes + 2;
  const auto want = divergences(net, {{net.num_nodes, r.value}, {net.num_nodes + 1, -static_cast<__int128>(r.value)}});
  for (const auto& s : check_residual(net, res, r.potentials, k, want, (22 * n + 15) * cost_sum(*c.g))) {
    ADD_FAILURE() << where << ": final: " << s;
  }
  EXPECT_FALSE(connected(net, res)) << where;
  // At most m + 2 r Dijkstra runs.
  const Amount rem = st.value() - r.value;
  const Amount m = f - rem;
  EXPECT_GE(rem, 0) << where;
  EXPECT_LE(r.counters.dijkstra_runs, m + 2 * rem) << where;
}

// The captured state's own checks: bounds, divergences, valid working potentials within
// (n + 1) C, and the certificate.
void check_state(const Case& c, const DeletionState& st, const FlowCoreResult& payload, const std::string& where) {
  const Network& net = st.network();
  EXPECT_EQ(st.value(), payload.placed) << where;
  EXPECT_EQ(st.edge_flows(), payload.edge_flows) << where;
  const __int128 n = net.num_nodes + 2;
  const auto want =
      divergences(net, {{net.num_nodes, st.value()}, {net.num_nodes + 1, -static_cast<__int128>(st.value())}});
  for (const auto& s : check_residual(net, st.residual(), st.potentials(), kNone, want, (n + 1) * cost_sum(*c.g))) {
    ADD_FAILURE() << where << ": state: " << s;
  }
  EXPECT_EQ(st.maximum(), !connected(net, st.residual())) << where;
  const auto& x = st.counters();
  const auto& y = payload.counters;
  EXPECT_TRUE(x.dijkstra_runs == y.dijkstra_runs && x.dijkstra_pops == y.dijkstra_pops && x.settled == y.settled &&
              x.arc_scans == y.arc_scans && x.dinic_phases == y.dinic_phases && x.augmentations == y.augmentations)
      << where;
}

// Every allowed edge of the captured state: positive flow by the re-solve, zero flow by a
// solve from zero on V \ {e} (the baseline's value and cost; no re-solve, which rejects it).
void check_every_edge(const Case& c, const std::string& where) {
  FlowCoreResult payload;
  const auto req = capture_of(c, &payload);
  ASSERT_TRUE(req.state.has_value()) << where;
  const DeletionState& st = *req.state;
  check_state(c, st, payload, where);
  ASSERT_TRUE(st.maximum()) << where;
  const auto flows = st.edge_flows();
  for (EdgeId e : st.network().edges) {
    const std::string at = where + " edge " + std::to_string(e);
    if (flows[static_cast<std::size_t>(e)] > 0) {
      check_deletion(c, st, e, at);
    } else {
      auto view = c.view;
      view[static_cast<std::size_t>(e)] = 0;
      const auto fresh = place_global_units(*c.g, view, c.sources, c.destinations, c.el(), std::nullopt);
      EXPECT_EQ(fresh.placed, st.value()) << at;
      EXPECT_TRUE(cost_of(*c.g, fresh.edge_flows) == cost_of(*c.g, flows)) << at;
      DeletionResolver rs(st);
      EXPECT_THROW((void)rs.resolve(e), std::invalid_argument) << at;
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
  for (int i = 0; i < 20000; ++i) {
    const I128 a = any();
    const I128 b = any();
    const __int128 x = wide(a);
    const __int128 y = wide(b);
    const auto wrap = [](__int128 z) { return static_cast<__int128>(static_cast<unsigned __int128>(z)); };
    EXPECT_TRUE(wide(a + b) == wrap(static_cast<__int128>(static_cast<unsigned __int128>(x) + static_cast<unsigned __int128>(y))));
    EXPECT_TRUE(wide(a - b) == wrap(static_cast<__int128>(static_cast<unsigned __int128>(x) - static_cast<unsigned __int128>(y))));
    EXPECT_EQ(a < b, x < y);
    EXPECT_EQ(a == b, x == y);
    const auto c = static_cast<std::int64_t>(rng());
    EXPECT_TRUE(wide(I128 {c}) == c);
  }
}

// ------------------------------------------------------------------ the capture

// The capture changes no field of the ordinary payload (every field, all six counters): cap 0
// (no search), an infeasible positive cap, caps below, at and above the maximum, uncapped;
// never-reached sentinels and the final failed tier included.
TEST(FlowCoreDeletion, CaptureLeavesTheFromZeroPayloadUnchanged) {
  std::mt19937_64 rng(5);
  int compared = 0;
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
      ++compared;
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
  EXPECT_GT(compared, 1000);
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
namespace {

// Captures c with the independent iterative all-vertex update beside: after every search, each
// vertex by its distance or, unreached, by the search's largest finite distance, from 0.
DeletionCapture capture_iterative(const Case& c, std::vector<__int128>& pi, int& searches) {
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
  int searched = 0;
  for (int trial = 0; trial < 300; ++trial) {
    Case c = random_case(rng);
    const Amount fmax = independent_curve(*c.g, c.view, c.sources, c.destinations, c.el()).back().first;
    for (std::optional<Amount> cap : {std::optional<Amount> {}, std::optional<Amount> {0},
                                      std::optional<Amount> {fmax / 2}, std::optional<Amount> {fmax + 1}}) {
      c.cap = cap;
      std::vector<__int128> pi;
      int searches = 0;
      const auto req = capture_iterative(c, pi, searches);
      const auto& st = *req.state;
      ASSERT_EQ(st.potentials().size(), pi.size());
      for (std::size_t v = 0; v < pi.size(); ++v) ASSERT_TRUE(wide(st.potentials()[v]) == pi[v]) << trial << " " << v;
      EXPECT_EQ(searches, st.counters().dijkstra_runs);
      if (cap == 0) EXPECT_EQ(searches, 0);
      searched += searches;
    }
  }
  EXPECT_GT(searched, 1000);
  // The common-offset mutant: valid reduced costs, the wrong convention.
  Case c = make_case(5, {{0, 1, 3, 3}, {0, 2, 1, 3}, {1, 2, 2, 7}, {2, 3, 1, 3}, {3, 1, 1, 2}}, {0}, {3});
  std::vector<__int128> pi;
  int searches = 0;
  const auto req = capture_iterative(c, pi, searches);
  const auto& st = *req.state;
  std::vector<I128> shifted(st.potentials().begin(), st.potentials().end());
  for (auto& p : shifted) p += 1000;
  const Network& net = st.network();
  const auto want =
      divergences(net, {{net.num_nodes, st.value()}, {net.num_nodes + 1, -static_cast<__int128>(st.value())}});
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
  int deletions = 0;
  for (int trial = 0; trial < 700; ++trial) {
    Case c = random_case(rng);
    const Amount fmax = independent_curve(*c.g, c.view, c.sources, c.destinations, c.el()).back().first;
    // Uncapped, and positive caps at and above the independently established maximum.
    for (std::optional<Amount> cap : {std::optional<Amount> {}, std::optional<Amount> {std::max<Amount>(fmax, 1)},
                                      std::optional<Amount> {fmax + 2}}) {
      c.cap = cap;
      check_every_edge(c, "trial " + std::to_string(trial));
      if (::testing::Test::HasFailure()) return;
    }
    const auto req = capture_of(c);
    for (Amount x : req.state->edge_flows()) deletions += x > 0 ? 1 : 0;
  }
  EXPECT_GT(deletions, 600);
}

// ------------------------------------------------------------------ witnesses

// Repair 1 reroutes everything (m = f, r = 0): u -> w -> v replaces u -> v.
TEST(FlowCoreDeletion, FullRepairOneReroute) {
  const Case c = make_case(5, {{0, 1, 1, 1}, {1, 2, 1, 1}, {1, 3, 1, 5}, {3, 2, 1, 5}, {2, 4, 1, 1}}, {0}, {4});
  check_every_edge(c, "full reroute");
  const auto req = capture_of(c);
  const EdgeId uv = id_of_input(*c.g)[1];
  const auto r = DeletionResolver(*req.state).resolve(uv);
  EXPECT_EQ(r.value, 1);
  EXPECT_TRUE(cost_of(*c.g, r.edge_flows) == 12);
  EXPECT_EQ(r.counters.dijkstra_runs, 1);
}

// No reroute (m = 0, r = 1): deleting the sole unit-capacity source-to-destination edge takes
// exactly two Dijkstra runs; the control that reruns repair 1's failed search takes three and
// fails the m + 2 r bound, with the same flow, potentials, phases and paths; the reuse skips
// exactly that search's work.
TEST(FlowCoreDeletion, SoleEdgeTakesTwoSearches) {
  const Case c = make_case(2, {{0, 1, 1, 1}}, {0}, {1});
  check_every_edge(c, "sole edge");
  const auto req = capture_of(c);
  DeletionResolver rs(*req.state);
  const auto r = rs.resolve(0);
  EXPECT_EQ(r.value, 0);
  EXPECT_EQ(r.counters.dijkstra_runs, 2);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  std::vector<std::string> fails;
  std::vector<std::pair<int, FlowCoreCounters>> reuse_updates;
  std::vector<std::pair<int, FlowCoreCounters>> rerun_updates;
  const auto a = observed(rs, *req.state, 0, fails, false, &reuse_updates);
  const auto b = observed(rs, *req.state, 0, fails, true, &rerun_updates);
  for (const auto& s : fails) ADD_FAILURE() << s;
  EXPECT_EQ(a.counters.dijkstra_runs, 2);
  EXPECT_EQ(b.counters.dijkstra_runs, 3);  // > m + 2 r = 2: the control fails the bound
  EXPECT_EQ(a.edge_flows, b.edge_flows);
  EXPECT_EQ(a.potentials, b.potentials);
  EXPECT_EQ(a.cut, b.cut);
  EXPECT_EQ(a.counters.dinic_phases, b.counters.dinic_phases);
  EXPECT_EQ(a.counters.augmentations, b.counters.augmentations);
  // The control's extra search: the counters between repair 1's last update and repair 2's
  // first (no Dinic work between them).
  ASSERT_EQ(rerun_updates.size(), reuse_updates.size() + 1);
  std::size_t i = 0;
  while (rerun_updates[i].first != 2) ++i;
  const auto& before = rerun_updates[i - 1].second;
  const auto& after = rerun_updates[i].second;
  EXPECT_EQ(after.dinic_phases, before.dinic_phases);
  EXPECT_EQ(after.augmentations, before.augmentations);
  EXPECT_EQ(b.counters.dijkstra_runs - (after.dijkstra_runs - before.dijkstra_runs), a.counters.dijkstra_runs);
  EXPECT_EQ(b.counters.dijkstra_pops - (after.dijkstra_pops - before.dijkstra_pops), a.counters.dijkstra_pops);
  EXPECT_EQ(b.counters.settled - (after.settled - before.settled), a.counters.settled);
  EXPECT_EQ(b.counters.arc_scans - (after.arc_scans - before.arc_scans), a.counters.arc_scans);
  EXPECT_GT(after.settled - before.settled, 0);
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
  EXPECT_TRUE(cost_of(*c.g, r.edge_flows) == 6);
  EXPECT_LE(r.counters.dijkstra_runs, 1 + 2 * 1);
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
  auto view = positive.view;
  view[static_cast<std::size_t>(spare)] = 0;
  const auto fresh = place_global_units(*positive.g, view, positive.sources, positive.destinations, nullptr, std::nullopt);
  EXPECT_NE(*base.min_cut, *fresh.min_cut);
  const Case zero = make_case(2, {{0, 1, 0, 1}}, {0}, {1});
  check_every_edge(zero, "zero view capacity");
  const auto zreq = capture_of(zero);
  EXPECT_TRUE(zreq.state->maximum());
  EXPECT_EQ(zreq.state->value(), 0);
}

// Parallel equal-cost ties: the re-solve's flow differs from the fresh solve's, its value, cost
// and cut do not (review-resolve/focused.py).
TEST(FlowCoreDeletion, TiedFlowsDifferValueCostAndCutAgree) {
  const Case c = make_case(3, {{0, 1, 1, 2}, {0, 1, 2, 3}, {1, 2, 1, 1}, {1, 2, 1, 1}}, {0}, {2});
  check_every_edge(c, "ties");
  const auto req = capture_of(c);
  const EdgeId dear = id_of_input(*c.g)[1];
  const auto r = DeletionResolver(*req.state).resolve(dear);
  auto view = c.view;
  view[static_cast<std::size_t>(dear)] = 0;
  const auto fresh = place_global_units(*c.g, view, c.sources, c.destinations, nullptr, std::nullopt);
  EXPECT_EQ(r.edge_flows, (std::vector<Amount> {1, 0, 0, 1}));
  EXPECT_EQ(fresh.edge_flows, (std::vector<Amount> {1, 0, 1, 0}));
  EXPECT_EQ(r.value, 1);
  EXPECT_TRUE(cost_of(*c.g, r.edge_flows) == 3);
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
  check_every_edge(c, "endpoints");
  c.elig = c.g->eligibility(std::vector<NodeId> {1}, std::vector<EdgeId> {});
  check_every_edge(c, "masked source");
  c.elig = c.g->eligibility(std::vector<NodeId> {4}, std::vector<EdgeId> {id_of_input(*c.g)[3]});
  check_every_edge(c, "masked destination and edge");
}

// The near-domain cost witness of review-resolve/REPORT.md C2: H = 2^61 and a cost sum H + 11
// below 2^62; the working potentials exceed int64 after the deletion, vertex 4 is never
// reached (the published sentinel, a numeric working potential).
TEST(FlowCoreDeletion, NearDomainCostsNeedWideWorkingPotentials) {
  const Cost H = Cost {1} << 61;
  const Case c = make_case(5, {{0, 1, 3, 3}, {0, 2, 1, 3}, {1, 2, 2, H}, {2, 3, 1, 3}, {3, 1, 1, 2}}, {0}, {3});
  ASSERT_TRUE(cost_sum(*c.g) == static_cast<__int128>(H) + 11);
  FlowCoreResult payload;
  const auto req = capture_of(c, &payload);
  const auto& st = *req.state;
  EXPECT_EQ(payload.potentials[4], kNeverReached);
  const std::vector<__int128> base {0, 3, static_cast<__int128>(H) + 3, static_cast<__int128>(H) + 6,
                                    static_cast<__int128>(H) + 6, 0, static_cast<__int128>(H) + 6};
  for (std::size_t v = 0; v < base.size(); ++v) EXPECT_TRUE(wide(st.potentials()[v]) == base[v]) << v;
  const EdgeId e = id_of_input(*c.g)[3];
  check_deletion(c, st, e, "near domain");
  const auto r = DeletionResolver(st).resolve(e);
  const __int128 h = H;
  const std::vector<__int128> after {3 * h + 5, 2 * h + 8, 3 * h + 8, 2 * h + 6, 4 * h + 11, 3 * h + 5, 2 * h + 6};
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
  const DeletionState saved = *req.state;
  DeletionResolver rs(*req.state);
  const EdgeId cheap = id_of_input(*c.g)[0];
  EXPECT_THROW((void)rs.resolve(cheap), std::invalid_argument);
  EXPECT_EQ(std::vector<Amount>(saved.residual().begin(), saved.residual().end()),
            std::vector<Amount>(req.state->residual().begin(), req.state->residual().end()));
  for (std::optional<Amount> cap : {std::optional<Amount> {}, std::optional<Amount> {4}}) {
    c.cap = cap;
    const auto ok = capture_of(c);
    ASSERT_TRUE(ok.state->maximum());
    EXPECT_EQ(DeletionResolver(*ok.state).resolve(cheap).value, 3);
  }
  // Not an edge of the view, no flow, an empty state.
  c.cap = std::nullopt;
  c.elig = c.g->eligibility(std::vector<NodeId> {}, std::vector<EdgeId> {cheap});
  const auto masked = capture_of(c);
  EXPECT_THROW((void)DeletionResolver(*masked.state).resolve(cheap), std::invalid_argument);
  EXPECT_THROW(DeletionResolver {DeletionState {}}, std::invalid_argument);
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
    std::vector<EdgeId> positive;
    for (EdgeId e : st.network().edges) {
      if (st.edge_flows()[static_cast<std::size_t>(e)] > 0) positive.push_back(e);
    }
    if (positive.size() < 2) continue;
    ++checked;
    const std::vector<Amount> res0(st.residual().begin(), st.residual().end());
    const std::vector<I128> pi0(st.potentials().begin(), st.potentials().end());
    const DeletionState copy = st;  // a copy shares the network only
    const Case other = random_case(rng);
    (void)place_global_units(*other.g, other.view, other.sources, other.destinations, other.el(), std::nullopt);
    const EdgeId A = positive.front();
    const EdgeId B = positive.back();
    auto same = [](const DeletionResult& x, const DeletionResult& y) {
      return x.value == y.value && x.potentials == y.potentials && x.cut == y.cut && x.edge_flows == y.edge_flows &&
             x.counters.dijkstra_runs == y.counters.dijkstra_runs && x.counters.dijkstra_pops == y.counters.dijkstra_pops &&
             x.counters.settled == y.counters.settled && x.counters.arc_scans == y.counters.arc_scans &&
             x.counters.dinic_phases == y.counters.dinic_phases && x.counters.augmentations == y.counters.augmentations;
    };
    DeletionResolver one(copy);
    const auto a1 = one.resolve(A);
    const auto b1 = one.resolve(B);
    const auto a2 = one.resolve(A);
    EXPECT_TRUE(same(a1, a2)) << trial;
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
        EXPECT_TRUE(same(out[0], a1) && same(out[2], a1) && same(out[1], b1)) << trial << " threads " << threads;
      }
    }
    EXPECT_EQ(std::vector<Amount>(st.residual().begin(), st.residual().end()), res0) << trial;
    EXPECT_EQ(std::vector<I128>(st.potentials().begin(), st.potentials().end()), pi0) << trial;
    EXPECT_EQ(std::vector<Amount>(copy.residual().begin(), copy.residual().end()), res0) << trial;
    EXPECT_EQ(&copy.network(), &st.network());
  }
  EXPECT_EQ(checked, 20);
}

// ------------------------------------------------------------------ sensitivity

namespace {

FlowPolicyConfig global_config() {
  return FlowPolicyConfig {Global {}, FlowPlacement::MaxFlow, Admission::Lossless,
                           PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity)};
}

} // namespace

// The exact model's loss is placed minus capacity() with the edge excluded, for every
// candidate (the baseline's positive-flow edges), on 1, 2 and 4 workers; zero-flow edges are
// not candidates.
TEST(Sensitivity, GlobalLossesEqualCapacityWithTheEdgeExcluded) {
  std::mt19937_64 rng(41);
  int entries = 0;
  for (int trial = 0; trial < 300; ++trial) {
    const Case c = random_case(rng);
    const SensitivityOptions opts {global_config(), c.el()};
    const auto base = capacity(*c.g, c.sources, c.destinations, CapacityOptions {global_config(), c.el(), nullptr,
                                                                                 OutputOptions {true, false, false}});
    const auto s1 = sensitivity(*c.g, c.sources, c.destinations, opts, 1);
    EXPECT_EQ(s1.baseline.placed, base.placed) << trial;
    std::vector<EdgeId> want;
    for (EdgeId e = 0; e < c.g->num_edges(); ++e) {
      if ((*base.edge_flows)[static_cast<std::size_t>(e)] > 0) want.push_back(e);
    }
    ASSERT_EQ(s1.entries.size(), want.size()) << trial;
    for (std::size_t i = 0; i < want.size(); ++i) {
      const auto& entry = s1.entries[i];
      EXPECT_EQ(entry.edge, want[i]) << trial;
      EXPECT_EQ(entry.status, SensitivityStatus::Evaluated) << trial;
      const EdgeId one[1] = {entry.edge};
      Eligibility without = c.g->eligibility({}, one);
      if (c.elig) without = Eligibility::intersect(*c.elig, without);
      const auto fresh = capacity(*c.g, c.sources, c.destinations, CapacityOptions {global_config(), &without});
      EXPECT_EQ(entry.loss, base.placed - fresh.placed) << trial << " edge " << entry.edge;
      ++entries;
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
}

// OutputOptions controls only the published baseline: every combination captures the same state
// (the same entries) and the baseline's optional fields are capacity()'s.
TEST(Sensitivity, GlobalBaselineKeepsEveryOutputsCombination) {
  const Case c = make_case(5, {{0, 1, 2, 1}, {1, 2, 2, 1}, {1, 3, 1, 2}, {3, 2, 1, 2}, {2, 4, 2, 1}, {0, 3, 1, 5}},
                           {0}, {4});
  const auto reference = sensitivity(*c.g, c.sources, c.destinations, SensitivityOptions {global_config()}, 1);
  ASSERT_FALSE(reference.entries.empty());
  for (int bits = 0; bits < 8; ++bits) {
    const OutputOptions outputs {(bits & 1) != 0, (bits & 2) != 0, (bits & 4) != 0};
    SensitivityOptions opts {global_config()};
    opts.outputs = outputs;
    const auto s = sensitivity(*c.g, c.sources, c.destinations, opts, 1);
    const auto cap = capacity(*c.g, c.sources, c.destinations, CapacityOptions {global_config(), nullptr, nullptr, outputs});
    EXPECT_EQ(s.baseline.placed, cap.placed) << bits;
    EXPECT_EQ(s.baseline.cost_distribution, cap.cost_distribution) << bits;
    EXPECT_EQ(s.baseline.edge_flows, cap.edge_flows) << bits;
    EXPECT_EQ(s.baseline.residual, cap.residual) << bits;
    EXPECT_EQ(s.baseline.bottleneck, cap.bottleneck) << bits;
    EXPECT_EQ(s.baseline.min_cut, cap.min_cut) << bits;
    EXPECT_EQ(s.baseline.edge_flows.has_value(), outputs.edge_flows) << bits;
    ASSERT_EQ(s.entries.size(), reference.entries.size()) << bits;
    for (std::size_t i = 0; i < s.entries.size(); ++i) {
      EXPECT_EQ(s.entries[i].edge, reference.entries[i].edge) << bits;
      EXPECT_EQ(s.entries[i].loss, reference.entries[i].loss) << bits;
    }
  }
}
