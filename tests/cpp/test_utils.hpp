#pragma once

#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <limits>
#include <map>
#include <numeric>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/detail/flow_core.hpp"
#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/program.hpp"
#include "test_int128.hpp"

namespace netgraph::core::test {

using detail::FlowCoreResult;
using detail::kMaxFlowCap;
using detail::kNeverReached;
using detail::place_global_units;

// An edge in input coordinates.
struct InEdge {
  NodeId u;
  NodeId v;
  Amount cap;
  Cost cost;
};

// 2^62: the bound of capacities and of the capacity and cost sums, in units.
inline constexpr Amount kL62 = Amount {1} << 62;

// The two policy configurations most tests use, both Lossless with deterministic ties and the
// capacity priority: the exact maximum flow (Global routing, the optimized origin) and CostOnly
// EqualBalanced.
inline FlowPolicyConfig global_config() { return FlowPolicyConfig {.routing = Global {}}; }
inline FlowPolicyConfig cost_only_config() {
  return FlowPolicyConfig {CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless,
                           PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity)};
}

// The five statistics of a DagCache, comparable as one value.
inline std::tuple<std::int64_t, std::int64_t, std::int64_t, std::int64_t, std::int64_t> stats_of(
    const DagCache::Stats& s) {
  return {s.hits, s.misses, s.entries, s.bytes, s.live_values};
}

// The kernel-combination errors of the configuration check.
inline const std::string kLossyKernelError = "MAX_FLOW requires admission=LOSSLESS";
inline const std::string kHeadroomKernelError = "the Weighted kernel is rejected with priority Headroom (weights are units)";

inline Graph make_graph(std::int32_t n, const std::vector<InEdge>& edges, double q = 1.0) {
  std::vector<NodeId> src, dst;
  std::vector<Amount> cap;
  std::vector<Cost> cost;
  for (const auto& e : edges) {
    src.push_back(e.u);
    dst.push_back(e.v);
    cap.push_back(e.cap);
    cost.push_back(e.cost);
  }
  return Graph(n, src, dst, cap, cost, q);
}

// m random edges on n nodes (self-loops and parallel edges allowed): capacities min_capacity to 4
// units, costs 1 to max_cost times scale. The draws are in a braced list, so their order is fixed.
inline Graph random_graph(std::mt19937& rng, int n, int m, Cost max_cost, Cost scale = 1, Amount min_capacity = 0) {
  std::vector<InEdge> edges;
  for (int i = 0; i < m; ++i) {
    edges.push_back({static_cast<NodeId>(rng() % static_cast<unsigned>(n)), static_cast<NodeId>(rng() % static_cast<unsigned>(n)),
                     min_capacity + static_cast<Amount>(rng() % static_cast<unsigned>(5 - min_capacity)),
                     scale * static_cast<Cost>(1 + rng() % static_cast<unsigned>(max_cost))});
  }
  return make_graph(n, edges);
}

// The message of the std::invalid_argument fn throws; "" when it throws none.
inline std::string error_of(const std::function<void()>& fn) {
  try {
    fn();
  } catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
}

// The failures fn records, intercepted (they do not fail the calling test): the negative
// controls of the tests' own checks and comparators.
inline int failures_of(const std::function<void()>& fn) {
  ::testing::TestPartResultArray results;
  {
    const ::testing::ScopedFakeTestPartResultReporter reporter(
        ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD, &results);
    fn();
  }
  return results.size();
}

inline bool contains(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

// Whether a checker's failed clauses (a list of messages) include one that starts with `clause`.
inline bool has_failure(const std::vector<std::string>& fails, const std::string& clause) {
  return std::any_of(fails.begin(), fails.end(), [&](const std::string& f) { return f.rfind(clause, 0) == 0; });
}

inline std::shared_ptr<const Graph> shared(Graph g) { return std::make_shared<const Graph>(std::move(g)); }

// A node list written inline: v({0, 1}).
inline std::vector<NodeId> v(std::initializer_list<NodeId> x) { return x; }

// Per-edge amount arrays with one or several defects (one error; detail::check_edge_amounts) for
// edge_amount_graph(), three edges of capacity 2, so the texts hold whatever the EdgeId order:
// the length is reported before any entry, then the first entry outside [0, 2] in EdgeId order.
// `error` follows the caller's label ("residual", "view").
struct EdgeAmountDefect {
  std::vector<Amount> amounts;
  std::string error;
};

inline Graph edge_amount_graph() { return make_graph(4, {{0, 1, 2, 1}, {1, 2, 2, 1}, {2, 3, 2, 1}}); }

inline std::vector<EdgeAmountDefect> edge_amount_defects() {
  return {
      {{9, -1, 0, 0}, " has length 4, expected 3 (one entry per edge)"},
      {{-5, 9}, " has length 2, expected 3 (one entry per edge)"},
      {{}, " has length 0, expected 3 (one entry per edge)"},
      {{0, -1, 3}, "[1] = -1 is outside [0, capacity_units[1] = 2]"},
      {{0, 3, -1}, "[1] = 3 is outside [0, capacity_units[1] = 2]"},
      {{2, 0, 3}, "[2] = 3 is outside [0, capacity_units[2] = 2]"},
      {{-1, -1, -1}, "[0] = -1 is outside [0, capacity_units[0] = 2]"},
  };
}

// Both bounds of [0, capacity] are inside.
inline std::vector<Amount> edge_amounts_at_the_bounds() { return {0, 2, 1}; }

// A span (a public array, a plan's structure field) as a vector.
template <class T>
std::vector<std::remove_const_t<T>> vec(std::span<T> s) {
  return std::vector<std::remove_const_t<T>>(s.begin(), s.end());
}

// Every EdgeId of g, ascending.
inline std::vector<EdgeId> all_ids(const Graph& g) {
  std::vector<EdgeId> out(static_cast<std::size_t>(g.num_edges()));
  for (std::size_t i = 0; i < out.size(); ++i) out[i] = static_cast<EdgeId>(i);
  return out;
}

// The PredDAG of a set of EdgeIds (the builder's of_edges).
inline std::shared_ptr<const PredDAG> dag_of(const Graph& g, std::vector<EdgeId> edges) {
  return std::make_shared<const PredDAG>(detail::PredDAGBuilder::of_edges(g, edges));
}

// A DAG's edges in ascending EdgeId (the builder's view of the storage).
inline std::vector<EdgeId> ascending_edges(const PredDAG& dag) {
  const auto e = detail::PredDAGBuilder::ascending(dag);
  return {e.begin(), e.end()};
}

// The sorted, distinct nodes of `set` that `elig` allows (all of them without an Eligibility).
inline std::vector<NodeId> allowed_nodes(std::vector<NodeId> set, const Eligibility* elig) {
  std::sort(set.begin(), set.end());
  set.erase(std::unique(set.begin(), set.end()), set.end());
  std::vector<NodeId> out;
  for (NodeId v : set) {
    if (elig == nullptr || elig->node_ok(v)) out.push_back(v);
  }
  return out;
}

// The bundle of every edge of g.
inline std::shared_ptr<const PredDAG> all_edges(const Graph& g) { return dag_of(g, all_ids(g)); }

// Pinned routing over a supplied program of split routes of one head, in order (weight 1 each;
// placement absent: EqualBalanced).
inline Pinned pinned_dags(const Graph& g, NodeId src, NodeId dst, const std::vector<std::shared_ptr<const PredDAG>>& dags,
                          std::optional<FlowPlacement> placement = std::nullopt, Admission admission = Admission::Lossless) {
  (void)admission;
  std::vector<Route> routes;
  for (const auto& d : dags) routes.push_back(Route {src, d, 1});
  RouteOptions options;
  options.placement = placement;
  const NodeId s[1] = {src};
  const NodeId t[1] = {dst};
  return Pinned {ForwardingProgram::from_routes(g, s, t, routes, options), nullptr};
}

// Pinned routing over a supplied program of tunnels: (head, EdgeIds in path order, weight).
inline Pinned pinned_paths(const Graph& g, const std::vector<NodeId>& sources, const std::vector<NodeId>& destinations,
                           const std::vector<std::tuple<NodeId, std::vector<EdgeId>, Amount>>& tunnels) {
  std::vector<Route> routes;
  for (const auto& [head, edges, weight] : tunnels) {
    std::vector<std::vector<EdgeId>> hops;
    for (EdgeId e : edges) hops.push_back({e});
    routes.push_back(Route {head, std::make_shared<const Path>(g, hops), weight});
  }
  return Pinned {ForwardingProgram::from_routes(g, sources, destinations, routes), nullptr};
}

// Every operation counter of two solves (the counter struct has no ==).
inline bool same_counters(const detail::FlowCoreCounters& x, const detail::FlowCoreCounters& y) {
  return x.dijkstra_runs == y.dijkstra_runs && x.dijkstra_pops == y.dijkstra_pops && x.settled == y.settled &&
         x.arc_scans == y.arc_scans && x.dinic_phases == y.dinic_phases && x.augmentations == y.augmentations;
}

// Every field of a network against its copy `in`.
inline bool same_network(const detail::Network& net, const detail::Network& in) {
  return net.num_nodes == in.num_nodes && net.edges == in.edges && net.tail == in.tail && net.head == in.head &&
         net.cost == in.cost && net.cap == in.cap && net.sources == in.sources && net.destinations == in.destinations &&
         net.off == in.off && net.arc_head == in.arc_head && net.rev == in.rev && net.arc_cost == in.arc_cost &&
         net.edge_slot == in.edge_slot && net.src_rev == in.src_rev && net.dst_fwd == in.dst_fwd;
}

// Every field of a deletion re-solve's result, its operation counters included.
inline bool same_result(const detail::DeletionResult& a, const detail::DeletionResult& b) {
  return a.value == b.value && a.potentials == b.potentials && a.cut == b.cut && a.edge_flows == b.edge_flows &&
         same_counters(a.counters, b.counters);
}

// The caller's inverse of input_index(): input position -> EdgeId.
inline std::vector<EdgeId> id_of_input(const Graph& g) {
  auto ii = g.input_index();
  std::vector<EdgeId> inv(ii.size());
  for (std::size_t e = 0; e < ii.size(); ++e) inv[static_cast<std::size_t>(ii[e])] = static_cast<EdgeId>(e);
  return inv;
}

// The path 0 -> 1 -> ... -> n - 1 (cost 1 per edge) at capacity `cap`.
inline std::shared_ptr<const Graph> chain_graph(std::int32_t n, Amount cap) {
  std::vector<InEdge> edges;
  for (NodeId v = 0; v + 1 < n; ++v) edges.push_back({v, v + 1, cap, 1});
  return std::make_shared<const Graph>(make_graph(n, edges));
}

// A per-edge array in input order, permuted to id order.
template <typename T>
std::vector<T> to_ids(const Graph& g, const std::vector<T>& by_input) {
  auto ii = g.input_index();
  std::vector<T> out(ii.size());
  for (std::size_t e = 0; e < ii.size(); ++e) out[e] = by_input[static_cast<std::size_t>(ii[e])];
  return out;
}

// A per-edge array in id order, permuted to input order.
template <typename T>
std::vector<T> to_input(const Graph& g, std::span<const T> by_id) {
  auto ii = g.input_index();
  std::vector<T> out(ii.size());
  for (std::size_t e = 0; e < ii.size(); ++e) out[static_cast<std::size_t>(ii[e])] = by_id[e];
  return out;
}

template <typename T>
std::vector<T> to_input(const Graph& g, const std::vector<T>& by_id) {
  return to_input(g, std::span<const T>(by_id));
}

// EdgeIds translated to sorted input positions.
inline std::vector<std::int32_t> ids_to_input(const Graph& g, const std::vector<EdgeId>& ids) {
  std::vector<std::int32_t> out;
  for (EdgeId e : ids) out.push_back(g.input_index()[static_cast<std::size_t>(e)]);
  std::sort(out.begin(), out.end());
  return out;
}

inline std::vector<std::pair<Cost, Amount>> curve(std::initializer_list<std::pair<Cost, Amount>> c) { return c; }

inline FlowCoreResult solve(const Graph& g, std::vector<NodeId> s, std::vector<NodeId> d,
                            std::optional<Amount> cap, const std::vector<Amount>* view = nullptr,
                            const Eligibility* elig = nullptr) {
  std::vector<Amount> full(g.capacity_units().begin(), g.capacity_units().end());
  return place_global_units(g, view ? std::span<const Amount>(*view) : std::span<const Amount>(full), s, d, elig, cap);
}

inline bool same_payload(const FlowCoreResult& a, const FlowCoreResult& b) {
  return a.placed == b.placed && a.cost_distribution == b.cost_distribution && a.edge_flows == b.edge_flows &&
         a.min_cut == b.min_cut && a.potentials == b.potentials && same_counters(a.counters, b.counters);
}

// Value -> minimum cost of every integer flow value, by enumerating every integer flow on the
// eligible edges (small graphs only). Independent of the solver.
inline std::map<Amount, Int128> exhaustive_costs(const Graph& g, std::span<const Amount> view,
                                                   const std::vector<NodeId>& srcs, const std::vector<NodeId>& dsts,
                                                   const Eligibility* elig) {
  const auto m = static_cast<std::size_t>(g.num_edges());
  const auto n = static_cast<std::size_t>(g.num_nodes());
  std::vector<std::size_t> free;
  for (std::size_t e = 0; e < m; ++e) {
    if (elig == nullptr || elig->edge_ok(static_cast<EdgeId>(e))) free.push_back(e);
  }
  std::set<NodeId> S, D;
  for (NodeId s : srcs) if (elig == nullptr || elig->node_ok(s)) S.insert(s);
  for (NodeId d : dsts) if (elig == nullptr || elig->node_ok(d)) D.insert(d);
  std::map<Amount, Int128> best;
  std::vector<Amount> f(free.size(), 0);
  while (true) {
    std::vector<Amount> bal(n, 0);
    Int128 cost = 0;
    for (std::size_t k = 0; k < free.size(); ++k) {
      const auto e = free[k];
      bal[static_cast<std::size_t>(g.edge_src()[e])] -= f[k];
      bal[static_cast<std::size_t>(g.edge_dst()[e])] += f[k];
      cost += static_cast<Int128>(g.cost()[e]) * f[k];
    }
    bool ok = true;
    Amount value = 0;
    for (std::size_t v = 0; v < n && ok; ++v) {
      const auto node = static_cast<NodeId>(v);
      if (S.count(node)) { ok = bal[v] <= 0; value -= bal[v]; }
      else if (D.count(node)) ok = bal[v] >= 0;
      else ok = bal[v] == 0;
    }
    if (ok) {
      auto it = best.find(value);
      if (it == best.end() || cost < it->second) best[value] = cost;
    }
    std::size_t k = 0;
    while (k < free.size() && f[k] == view[free[k]]) f[k++] = 0;
    if (k == free.size()) break;
    ++f[k];
  }
  return best;
}

// Breakpoints [(value, cost)] of the min-cost curve by successive shortest paths with
// Bellman-Ford on the residual (no potentials, no Dinic, its own arc order); run to
// disconnection. Independent of the solver.
inline std::vector<std::pair<Amount, Int128>> independent_curve(const Graph& g,
                                                                  std::span<const Amount> view,
                                                                  const std::vector<NodeId>& srcs,
                                                                  const std::vector<NodeId>& dsts,
                                                                  const Eligibility* elig) {
  struct Arc { int u, v; Amount cap; Cost w; Amount flow; };
  const int n = g.num_nodes();
  const int S = n, T = n + 1;
  std::vector<Arc> arcs;
  for (std::size_t e = 0; e < static_cast<std::size_t>(g.num_edges()); ++e) {
    if (elig != nullptr && !elig->edge_ok(static_cast<EdgeId>(e))) continue;
    arcs.push_back({g.edge_src()[e], g.edge_dst()[e], view[e], g.cost()[e], 0});
  }
  for (NodeId s : srcs) if (elig == nullptr || elig->node_ok(s)) arcs.push_back({S, s, kMaxFlowCap, 0, 0});
  for (NodeId d : dsts) if (elig == nullptr || elig->node_ok(d)) arcs.push_back({d, T, kMaxFlowCap, 0, 0});
  std::vector<std::pair<Amount, Int128>> pts{{0, 0}};
  while (true) {
    const Int128 inf = int128_max();
    std::vector<Int128> dist(static_cast<std::size_t>(n + 2), inf);
    std::vector<std::pair<int, int>> parent(static_cast<std::size_t>(n + 2), {-1, 0});
    dist[static_cast<std::size_t>(S)] = 0;
    for (int it = 0; it < n + 2; ++it) {
      bool changed = false;
      for (int i = 0; i < static_cast<int>(arcs.size()); ++i) {
        const auto& a = arcs[static_cast<std::size_t>(i)];
        const auto du = dist[static_cast<std::size_t>(a.u)];
        const auto dv = dist[static_cast<std::size_t>(a.v)];
        if (a.flow < a.cap && du != inf && du + a.w < dv) {
          dist[static_cast<std::size_t>(a.v)] = du + a.w;
          parent[static_cast<std::size_t>(a.v)] = {i, 1};
          changed = true;
        }
        if (a.flow > 0 && dist[static_cast<std::size_t>(a.v)] != inf &&
            dist[static_cast<std::size_t>(a.v)] - a.w < dist[static_cast<std::size_t>(a.u)]) {
          dist[static_cast<std::size_t>(a.u)] = dist[static_cast<std::size_t>(a.v)] - a.w;
          parent[static_cast<std::size_t>(a.u)] = {i, -1};
          changed = true;
        }
      }
      if (!changed) break;
    }
    if (dist[static_cast<std::size_t>(T)] == inf) return pts;
    std::vector<std::pair<int, int>> path;
    for (int v = T; v != S;) {
      auto [i, sign] = parent[static_cast<std::size_t>(v)];
      path.emplace_back(i, sign);
      const auto& a = arcs[static_cast<std::size_t>(i)];
      v = sign == 1 ? a.u : a.v;
    }
    Amount amount = std::numeric_limits<Amount>::max();
    for (auto [i, sign] : path) {
      const auto& a = arcs[static_cast<std::size_t>(i)];
      amount = std::min(amount, sign == 1 ? a.cap - a.flow : a.flow);
    }
    for (auto [i, sign] : path) arcs[static_cast<std::size_t>(i)].flow += sign * amount;
    pts.emplace_back(pts.back().first + amount, pts.back().second + dist[static_cast<std::size_t>(T)] * amount);
  }
}

// Cost of a piecewise-linear curve at an integer value inside its domain. Every segment of both
// curves has an integer slope (the cost per unit of one tier or one augmenting path), so the
// value is exact in 128 bits: costs stay below 2^62 per unit and amounts below 2^62. An empty or
// backward segment (a payload entry of volume 0 or less, which the curve-shape clause rejects) is
// skipped rather than divided by its length.
inline Int128 curve_at(const std::vector<std::pair<Amount, Int128>>& pts, Amount x) {
  for (std::size_t i = 1; i < pts.size(); ++i) {
    const auto [a, ca] = pts[i - 1];
    const auto [b, cb] = pts[i];
    if (a < b && a <= x && x <= b) {
      const Int128 slope = (cb - ca) / (b - a);
      return ca + slope * (x - a);
    }
  }
  return pts.back().second;
}

// The curve clauses shared by the flow certificates (certify below, the MaxFlow kernel step's in
// kernels_tests.cpp and the Global policy's cumulative curve in flow_policy_tests.cpp). The shape: every volume positive and the costs strictly increasing; the
// curve's total equals the delivered value and its moment (the sum of cost x volume) the flow's
// cost. Returns the curve's breakpoints [(value, cost)] from (0, 0).
inline std::vector<std::pair<Amount, Int128>> curve_shape_clauses(const std::vector<std::pair<Cost, Amount>>& curve,
                                                                  Amount placed, Int128 moment,
                                                                  std::vector<std::string>& fails) {
  std::vector<std::pair<Amount, Int128>> pts {{0, 0}};
  for (std::size_t i = 0; i < curve.size(); ++i) {
    const auto [c, u] = curve[i];
    if (u <= 0) fails.push_back("curve volume not positive");
    if (i > 0 && c <= curve[i - 1].first) fails.push_back("curve costs not strictly increasing");
    pts.emplace_back(pts.back().first + u, pts.back().second + static_cast<Int128>(c) * u);
  }
  if (pts.back().first != placed) fails.push_back("curve total differs from placed");
  if (pts.back().second != moment) fails.push_back("curve moment differs from the flow cost");
  return pts;
}

// The curve against an optimum's breakpoints: equal at every breakpoint of either up to `upto`
// (both are piecewise linear, so equal everywhere there); the first value that differs is
// reported, naming the optimum.
inline void curve_optimum_clause(const std::vector<std::pair<Amount, Int128>>& pts,
                                 const std::vector<std::pair<Amount, Int128>>& optimum, Amount upto,
                                 const std::string& name, std::vector<std::string>& fails) {
  std::set<Amount> xs;
  for (const auto& [x, c] : pts) xs.insert(x);
  for (const auto& [x, c] : optimum) xs.insert(x);
  for (Amount x : xs) {
    if (x > upto) continue;
    if (curve_at(pts, x) != curve_at(optimum, x)) {
      fails.push_back("curve differs from " + name + " at " + std::to_string(x));
      return;
    }
  }
}

// The cost of a flow given per EdgeId: the sum of cost x flow, exact in 128 bits.
inline Int128 flow_cost(const Graph& g, std::span<const Amount> flows) {
  Int128 total = 0;
  for (std::size_t e = 0; e < flows.size(); ++e) total += static_cast<Int128>(g.cost()[e]) * flows[e];
  return total;
}

// Per node, what a flow's sources send out and what its destinations receive (net amounts; 0 at
// every other node).
struct FlowBalance {
  std::vector<Int128> supply;
  std::vector<Int128> inflow;
};

// The flow clauses of the flow certificates (certify below, the deletion re-solve's in
// flow_core_deletion_tests.cpp and the Global policy's in flow_policy_tests.cpp), on a flow given
// per EdgeId: every edge flow within [0, view] and zero on an edge the Eligibility excludes; a
// source absorbs nothing and a destination emits nothing (a node in both sets counts as a
// source), an excluded one carries nothing, and every other node passes on what it receives;
// what the sources send out and what the destinations receive both equal `placed`. Appends the
// failed clauses to `fails`.
inline FlowBalance flow_clauses(const Graph& g, std::span<const Amount> view, std::span<const Amount> flows,
                                const std::vector<NodeId>& srcs, const std::vector<NodeId>& dsts,
                                const Eligibility* elig, Amount placed, std::vector<std::string>& fails) {
  const auto m = static_cast<std::size_t>(g.num_edges());
  const auto n = static_cast<std::size_t>(g.num_nodes());
  auto edge_ok = [&](std::size_t e) { return elig == nullptr || elig->edge_ok(static_cast<EdgeId>(e)); };
  auto node_ok = [&](NodeId v) { return elig == nullptr || elig->node_ok(v); };
  std::vector<Int128> bal(n, 0);  // inflow minus outflow
  for (std::size_t e = 0; e < m; ++e) {
    const Amount f = flows[e];
    if (f < 0 || f > view[e]) fails.push_back("feasibility at edge " + std::to_string(e));
    if (!edge_ok(e) && f != 0) fails.push_back("flow on an ineligible edge " + std::to_string(e));
    bal[static_cast<std::size_t>(g.edge_src()[e])] -= f;
    bal[static_cast<std::size_t>(g.edge_dst()[e])] += f;
  }
  const std::set<NodeId> S(srcs.begin(), srcs.end()), D(dsts.begin(), dsts.end());
  FlowBalance out {std::vector<Int128>(n, 0), std::vector<Int128>(n, 0)};
  Int128 supply = 0, inflow = 0;
  for (std::size_t v = 0; v < n; ++v) {
    const auto node = static_cast<NodeId>(v);
    if (S.count(node)) {
      if (bal[v] > 0) fails.push_back("source absorbs flow");
      if (!node_ok(node) && bal[v] != 0) fails.push_back("excluded source supplies");
      out.supply[v] = -bal[v];
      supply += -bal[v];
    } else if (D.count(node)) {
      if (bal[v] < 0) fails.push_back("destination emits flow");
      if (!node_ok(node) && bal[v] != 0) fails.push_back("excluded destination receives");
      out.inflow[v] = bal[v];
      inflow += bal[v];
    } else if (bal[v] != 0) {
      fails.push_back("conservation at node " + std::to_string(v));
    }
  }
  if (supply != placed || inflow != placed) fails.push_back("value differs from supply or inflow");
  return out;
}

// The residual closure of a flow given per EdgeId: the nodes that the allowed nodes of `seeds`
// reach over allowed edges, forward along an edge whose flow leaves room in the view and backward
// along an edge that carries flow (a unit there can be undone). Marks them with 1.
inline std::vector<char> residual_closure(const Graph& g, std::span<const Amount> view,
                                          std::span<const Amount> flows, const std::vector<NodeId>& seeds,
                                          const Eligibility* elig) {
  const auto m = static_cast<std::size_t>(g.num_edges());
  const auto n = static_cast<std::size_t>(g.num_nodes());
  std::vector<std::vector<std::size_t>> next(n);
  for (std::size_t e = 0; e < m; ++e) {
    if (elig != nullptr && !elig->edge_ok(static_cast<EdgeId>(e))) continue;
    const auto u = static_cast<std::size_t>(g.edge_src()[e]);
    const auto v = static_cast<std::size_t>(g.edge_dst()[e]);
    if (flows[e] < view[e]) next[u].push_back(v);
    if (flows[e] > 0) next[v].push_back(u);
  }
  std::vector<char> seen(n, 0);
  std::vector<std::size_t> queue;
  for (NodeId s : seeds) {
    const auto x = static_cast<std::size_t>(s);
    if ((elig == nullptr || elig->node_ok(s)) && seen[x] == 0) {
      seen[x] = 1;
      queue.push_back(x);
    }
  }
  for (std::size_t qi = 0; qi < queue.size(); ++qi) {
    for (const std::size_t v : next[queue[qi]]) {
      if (seen[v] == 0) {
        seen[v] = 1;
        queue.push_back(v);
      }
    }
  }
  return seen;
}

// The cut of a closure: every allowed edge from a node of the closure to a node outside it,
// whatever its view, ascending EdgeId (an edge of view 0 adds 0 to the value; it is listed
// because every edge leaving the closure is, also where it leads to a dead end). On a maximal
// flow, the closure of the sources holds no destination and this is the minimum cut the flow
// certifies.
inline std::vector<EdgeId> closure_cut(const Graph& g, const std::vector<char>& closure, const Eligibility* elig) {
  std::vector<EdgeId> cut;
  for (std::size_t e = 0; e < static_cast<std::size_t>(g.num_edges()); ++e) {
    if (elig != nullptr && !elig->edge_ok(static_cast<EdgeId>(e))) continue;
    if (closure[static_cast<std::size_t>(g.edge_src()[e])] != 0 &&
        closure[static_cast<std::size_t>(g.edge_dst()[e])] == 0) {
      cut.push_back(static_cast<EdgeId>(e));
    }
  }
  return cut;
}

// The flow core's optimality certificate in two ordered branches. Returns the failed clauses. With
// `tier_a`, the complete curve is compared with the exhaustive optimum of every value; else with
// the independent Bellman-Ford curve at the union of both breakpoint sets.
inline std::vector<std::string> certify(const Graph& g, std::span<const Amount> view,
                                        std::vector<NodeId> srcs, std::vector<NodeId> dsts,
                                        std::optional<Amount> cap, const Eligibility* elig,
                                        const FlowCoreResult& r, bool tier_a) {
  std::vector<std::string> fails;
  const auto m = static_cast<std::size_t>(g.num_edges());
  const auto n = static_cast<std::size_t>(g.num_nodes());
  std::sort(srcs.begin(), srcs.end());
  srcs.erase(std::unique(srcs.begin(), srcs.end()), srcs.end());
  std::sort(dsts.begin(), dsts.end());
  dsts.erase(std::unique(dsts.begin(), dsts.end()), dsts.end());
  auto edge_ok = [&](std::size_t e) { return elig == nullptr || elig->edge_ok(static_cast<EdgeId>(e)); };
  auto node_ok = [&](NodeId v) { return elig == nullptr || elig->node_ok(v); };
  if (r.edge_flows.size() != m || r.potentials.size() != n) {
    fails.push_back("shape");
    return fails;
  }
  if (cap.has_value() && *cap == 0) {
    if (r.placed != 0 || !r.cost_distribution.empty() || r.min_cut.has_value() ||
        std::any_of(r.edge_flows.begin(), r.edge_flows.end(), [](Amount x) { return x != 0; }) ||
        std::any_of(r.potentials.begin(), r.potentials.end(), [](Cost p) { return p != kNeverReached; }) ||
        !same_counters(r.counters, detail::FlowCoreCounters {})) {
      fails.push_back("zero cap: the prescribed zero payload");
    }
    return fails;
  }
  const FlowBalance balance = flow_clauses(g, view, r.edge_flows, srcs, dsts, elig, r.placed, fails);
  const Int128 moment = flow_cost(g, r.edge_flows);
  if (r.placed < 0 || (cap.has_value() && r.placed > *cap)) fails.push_back("placed outside [0, cap]");
  // The full residual closure, with a virtual source s* (an arc to every allowed source, and one
  // back from a source that sends) and a virtual sink t* (an arc from every allowed destination,
  // and one back to a destination that receives): t* is in it when an allowed destination is
  // (the destinations are still connected to the sources), and then so is everything a receiving
  // destination reaches; s* adds no node, every allowed source being a seed.
  std::vector<char> seen = residual_closure(g, view, r.edge_flows, srcs, elig);
  const bool connected = std::any_of(dsts.begin(), dsts.end(), [&](NodeId d) {
    return node_ok(d) && seen[static_cast<std::size_t>(d)] != 0;
  });
  if (connected) {
    std::vector<NodeId> seeds = srcs;
    for (NodeId d : dsts) {
      if (balance.inflow[static_cast<std::size_t>(d)] > 0) seeds.push_back(d);
    }
    seen = residual_closure(g, view, r.edge_flows, seeds, elig);
  }
  if (!cap.has_value() && connected) fails.push_back("uncapped solve left the destinations connected");
  if (cap.has_value() && r.placed != *cap && connected) fails.push_back("attained value: early stop");
  if (connected) {
    if (r.min_cut.has_value()) fails.push_back("min_cut present on a connected solve");
  } else {
    const std::vector<EdgeId> expect = closure_cut(g, seen, elig);
    Int128 cut_value = 0;
    for (const EdgeId e : expect) cut_value += view[static_cast<std::size_t>(e)];
    if (!r.min_cut.has_value()) fails.push_back("min_cut missing on a maximal solve");
    else if (*r.min_cut != expect) fails.push_back("min_cut differs from the closure cut");
    if (cut_value != r.placed) fails.push_back("cut value differs from placed");
  }
  // The vertices some search settled: the first search runs on the zero flow and settles every
  // vertex an allowed source reaches over eligible edges of positive view; a later search adds
  // only the reverses of carrying edges, whose ends that one already settled.
  std::vector<char> settled(n, 0);
  std::vector<std::size_t> from_zero;
  for (NodeId s : srcs) {
    const auto x = static_cast<std::size_t>(s);
    if (node_ok(s) && !settled[x]) { settled[x] = 1; from_zero.push_back(x); }
  }
  for (std::size_t qi = 0; qi < from_zero.size(); ++qi) {
    for (std::size_t e = 0; e < m; ++e) {
      const auto u = static_cast<std::size_t>(g.edge_src()[e]);
      const auto v = static_cast<std::size_t>(g.edge_dst()[e]);
      if (u == from_zero[qi] && edge_ok(e) && view[e] > 0 && !settled[v]) { settled[v] = 1; from_zero.push_back(v); }
    }
  }
  for (std::size_t v = 0; v < n; ++v) {
    if (r.potentials[v] < kNeverReached) fails.push_back("potential below the sentinel");
    if (seen[v] && r.potentials[v] == kNeverReached) fails.push_back("closure vertex holds NEVER_REACHED");
    if ((r.potentials[v] == kNeverReached) == (settled[v] != 0)) {
      fails.push_back("NEVER_REACHED differs from reachability at zero flow at node " + std::to_string(v));
    }
  }
  for (NodeId s : srcs) {
    const auto p = r.potentials[static_cast<std::size_t>(s)];
    if (node_ok(s) && p != 0) fails.push_back("allowed root with potential other than 0");
  }
  // Every residual arc between two closure vertices with potentials: the forward arc of an
  // allowed edge with room (its cost) and the backward arc of one that carries flow (minus it).
  auto reduced_cost_clause = [&](std::size_t u, std::size_t v, Int128 w) {
    if (!seen[u] || !seen[v] || r.potentials[u] == kNeverReached || r.potentials[v] == kNeverReached) return;
    if (w + r.potentials[u] - r.potentials[v] < 0) fails.push_back("negative reduced cost inside the closure");
  };
  for (std::size_t e = 0; e < m; ++e) {
    if (!edge_ok(e)) continue;
    const auto u = static_cast<std::size_t>(g.edge_src()[e]);
    const auto v = static_cast<std::size_t>(g.edge_dst()[e]);
    if (r.edge_flows[e] < view[e]) reduced_cost_clause(u, v, g.cost()[e]);
    if (r.edge_flows[e] > 0) reduced_cost_clause(v, u, -static_cast<Int128>(g.cost()[e]));
  }
  std::optional<Cost> pt;
  for (NodeId d : dsts) {
    const auto p = r.potentials[static_cast<std::size_t>(d)];
    if (p != kNeverReached) pt = pt ? std::min(*pt, p) : p;
  }
  if (pt) {
    for (NodeId d : dsts) {
      const auto p = r.potentials[static_cast<std::size_t>(d)];
      if (p == kNeverReached) continue;
      if (balance.inflow[static_cast<std::size_t>(d)] > 0 && p != *pt) {
        fails.push_back("inflow destination potential differs from pi(t*)");
      }
    }
    if (!r.cost_distribution.empty() && r.cost_distribution.back().first != *pt) {
      fails.push_back("pi(t*) differs from the last curve cost");
    }
  } else if (r.placed != 0) {
    fails.push_back("flow delivered but no destination reached");
  }
  const auto pts = curve_shape_clauses(r.cost_distribution, r.placed, moment, fails);
  std::vector<std::pair<Amount, Int128>> optimum;
  if (tier_a) {
    for (const auto& [x, c] : exhaustive_costs(g, view, srcs, dsts, elig)) optimum.emplace_back(x, c);
  } else {
    optimum = independent_curve(g, view, srcs, dsts, elig);
  }
  const Amount fmax = optimum.back().first;
  const Amount top = cap.has_value() ? std::min(*cap, fmax) : fmax;
  if (r.placed != top) fails.push_back("placed differs from min(cap, F_max)");
  curve_optimum_clause(pts, optimum, std::min(top, r.placed),
                       tier_a ? "the exhaustive optimum" : "the independent solver", fails);
  return fails;
}

// The fields in which two kernel plans differ, by name. `with_route` false leaves the route out: a
// plan built on a larger DAG keeps that DAG as its route but splits over the same members.
inline std::vector<std::string> plan_differences(const detail::KernelPlan& got, const detail::KernelPlan& want,
                                                 bool with_route = true) {
  std::vector<std::string> out;
  auto differ = [&](bool same, const char* name) {
    if (!same) out.emplace_back(name);
  };
  differ(got.placement == want.placement, "placement");
  differ(got.admission == want.admission, "admission");
  differ(vec(got.roots) == vec(want.roots), "roots");
  differ(vec(got.root_pos) == vec(want.root_pos), "root_pos");
  differ(vec(got.dsts) == vec(want.dsts), "dsts");
  if (with_route) differ(got.route == want.route, "route");
  differ(got.members == want.members, "members");
  differ(got.weights == want.weights, "weights");
  differ(vec(got.topo) == vec(want.topo), "topo");
  differ(vec(got.out_begin) == vec(want.out_begin), "out_begin");
  differ(vec(got.out_end) == vec(want.out_end), "out_end");
  differ(vec(got.head_pos) == vec(want.head_pos), "head_pos");
  differ(vec(got.terminal) == vec(want.terminal), "terminal");
  return out;
}

// One record on its own, from its fields alone: `edges` and `drops` hold positive amounts in
// strictly ascending EdgeId, only a Lossy record has `offered` and drops, and the edges carry
// the flow. Lossless: every node other than a source or a destination passes on what it
// receives, and `delivered` is the net inflow of the destinations (the net outflow of the
// sources). Lossy: a transit node passes on what it receives less what is dropped on the edges
// it forwards over, `delivered` is still the destinations' net inflow, and what was offered
// covers what was delivered and dropped (the rest is the source drops). Returns whether the
// record carries something (delivered or dropped).
inline bool expect_conserved(const Graph& g, const FlowRecord& r, const std::string& where) {
  const auto n = static_cast<std::size_t>(g.num_nodes());
  std::vector<Amount> in(n, 0), out(n, 0), dropped_out(n, 0);
  Amount total_drops = 0;
  for (std::size_t i = 0; i < r.edges.size(); ++i) {
    const auto [e, x] = r.edges[i];
    EXPECT_GT(x, 0) << where << " edge " << e;
    if (i > 0) EXPECT_LT(r.edges[i - 1].first, e) << where;
    out[static_cast<std::size_t>(g.edge_src()[static_cast<std::size_t>(e)])] += x;
    in[static_cast<std::size_t>(g.edge_dst()[static_cast<std::size_t>(e)])] += x;
  }
  const bool lossy = r.admission == Admission::Lossy;
  EXPECT_EQ(r.offered.has_value(), lossy) << where;
  if (!lossy) EXPECT_TRUE(r.drops.empty()) << where;
  for (std::size_t i = 0; i < r.drops.size(); ++i) {
    const auto& [e, c] = r.drops[i];
    EXPECT_EQ(c.hi, 0u) << where;  // the trials' drops stay below 2^64
    EXPECT_GT(c.lo, 0u) << where << " drop " << e;
    if (i > 0) EXPECT_LT(r.drops[i - 1].first, e) << where;
    const auto x = static_cast<Amount>(c.lo);
    dropped_out[static_cast<std::size_t>(g.edge_src()[static_cast<std::size_t>(e)])] += x;
    total_drops += x;
  }
  auto in_set = [](const std::vector<NodeId>& set, std::size_t v) {
    return std::find(set.begin(), set.end(), static_cast<NodeId>(v)) != set.end();
  };
  Amount into_destinations = 0, out_of_sources = 0;
  for (std::size_t v = 0; v < n; ++v) {
    if (in_set(r.destinations, v)) {
      into_destinations += in[v] - out[v];
    } else if (in_set(r.sources, v)) {
      out_of_sources += out[v] - in[v];
    } else {
      EXPECT_EQ(in[v], out[v] + dropped_out[v]) << where << " node " << v;
    }
  }
  EXPECT_EQ(into_destinations, r.delivered) << where;
  if (lossy) {
    if (r.offered.has_value() && r.offered->hi == 0) {
      EXPECT_GE(static_cast<Amount>(r.offered->lo), r.delivered + total_drops) << where;
    } else {
      ADD_FAILURE() << where << ": a Lossy record without offered, or offered beyond 2^64";
    }
  } else {
    EXPECT_EQ(out_of_sources, r.delivered) << where;
  }
  return r.delivered > 0 || total_drops > 0;
}

// An independent reference of a tunnel program's forwarding (contract 8.1 and 8.3): equal
// sources inject a total by round robin over the roster (the heads with a tunnel, ascending), each
// head splits its share by highest averages over its tunnels' weights (ascending index), and the
// lossless bound is the largest total whose loads fit the capacities with no share lost at a
// source. Its own round robin and unit-by-unit highest averages: it calls none of the library's
// apportionment helpers.
namespace program_reference {

struct RefTunnel {
  NodeId head {0};
  std::vector<EdgeId> edges {};
  Amount weight {0};
};

// x over k slots: x / k each, the first x mod k one more.
inline std::vector<Amount> round_robin(Amount x, std::size_t k) {
  std::vector<Amount> out(k, 0);
  for (Amount u = 0; u < x; ++u) ++out[static_cast<std::size_t>(u) % k];
  return out;
}

// x units one at a time, each to the slot with the largest w / (a + 1) (ties: the lower slot);
// slots of weight 0 get nothing (x must be 0 when every weight is 0).
inline std::vector<Amount> highest_averages(Amount x, const std::vector<Amount>& w) {
  std::vector<Amount> a(w.size(), 0);
  for (Amount u = 0; u < x; ++u) {
    std::size_t best = w.size();
    for (std::size_t i = 0; i < w.size(); ++i) {
      if (w[i] <= 0) continue;
      if (best == w.size() || Int128(w[i]) * (a[best] + 1) > Int128(w[best]) * (a[i] + 1)) best = i;
    }
    if (best == w.size()) return {};
    ++a[best];
  }
  return a;
}

// The roster: the heads with a tunnel, ascending.
inline std::vector<NodeId> roster(const std::vector<RefTunnel>& t) {
  std::set<NodeId> heads;
  for (const auto& x : t) heads.insert(x.head);
  return {heads.begin(), heads.end()};
}

// Per tunnel its amount of the total x; nullopt when a positive head share has no tunnel of
// positive weight (lost at the source).
inline std::optional<std::vector<Amount>> tunnel_amounts(Amount x, const std::vector<RefTunnel>& t) {
  const auto heads = roster(t);
  const auto shares = round_robin(x, heads.size());
  std::vector<Amount> out(t.size(), 0);
  for (std::size_t h = 0; h < heads.size(); ++h) {
    std::vector<std::size_t> ids;
    std::vector<Amount> w;
    for (std::size_t i = 0; i < t.size(); ++i) {
      if (t[i].head == heads[h]) {
        ids.push_back(i);
        w.push_back(t[i].weight);
      }
    }
    const auto split = highest_averages(shares[h], w);
    if (split.empty() && shares[h] > 0) return std::nullopt;
    for (std::size_t j = 0; j < split.size(); ++j) out[ids[j]] = split[j];
  }
  return out;
}

// The per-edge loads of per-tunnel amounts.
inline std::vector<Amount> loads(const Graph& g, const std::vector<RefTunnel>& t, const std::vector<Amount>& amounts) {
  std::vector<Amount> out(static_cast<std::size_t>(g.num_edges()), 0);
  for (std::size_t i = 0; i < t.size(); ++i) {
    for (EdgeId e : t[i].edges) out[static_cast<std::size_t>(e)] += amounts[i];
  }
  return out;
}

inline bool fits(const Graph& g, const std::vector<RefTunnel>& t, Amount x) {
  const auto a = tunnel_amounts(x, t);
  if (!a.has_value()) return false;
  const auto l = loads(g, t, *a);
  for (std::size_t e = 0; e < l.size(); ++e) {
    if (l[e] > g.capacity_units()[e]) return false;
  }
  return true;
}

// The largest total that fits, scanning up from 0 (the loads only grow with the total); a total
// above the sum of the capacities never fits.
inline Amount lossless_bound(const Graph& g, const std::vector<RefTunnel>& t) {
  Amount limit = 1;
  for (Amount c : g.capacity_units()) limit += c;
  Amount x = 0;
  while (x < limit && fits(g, t, x + 1)) ++x;
  return x;
}

} // namespace program_reference

// An independent reference of T1 (contract 6.2 to 6.5, 8.3): the setup by brute force over the
// simple paths of small graphs and the lossless bound with the policy's whole roster. It calls
// none of the library's searches or apportionment helpers (its own round robin, Bellman-Ford
// distances, exact comparisons in Int128).
namespace t1_reference {

struct Selection {
  bool prefer {false};  // PreferHigherPriority
  EdgePriority priority {EdgePriority::Capacity};
  bool multi_edge {true};
};

struct RefTunnel {
  NodeId head {0};
  std::size_t slot {0};
  Amount request {0};
  TunnelStatus status {TunnelStatus::Up};
  std::vector<EdgeId> path {};  // walk order; empty when down
  Cost cost {0};
};

// RR(x, k, i) = x // k + [i < x % k].
inline Amount rr(Amount x, std::int64_t k, std::int64_t i) { return x / k + (i < x % k ? 1 : 0); }

// The sources that reach an allowed destination over allowed edges (capacity ignored), ascending.
inline std::vector<NodeId> roster(const Graph& g, const Eligibility* elig, const std::vector<NodeId>& sources,
                                  const std::vector<NodeId>& dests) {
  const auto n = static_cast<std::size_t>(g.num_nodes());
  std::vector<char> is_dst(n, 0);
  for (NodeId d : dests) {
    if (elig == nullptr || elig->node_ok(d)) is_dst[static_cast<std::size_t>(d)] = 1;
  }
  std::vector<NodeId> out;
  for (NodeId s : sources) {
    if (elig != nullptr && !elig->node_ok(s)) continue;
    std::vector<char> seen(n, 0);
    std::vector<NodeId> stack {s};
    seen[static_cast<std::size_t>(s)] = 1;
    bool reached = false;
    while (!stack.empty() && !reached) {
      const NodeId u = stack.back();
      stack.pop_back();
      for (EdgeId e = 0; e < g.num_edges(); ++e) {
        if (g.edge_src()[static_cast<std::size_t>(e)] != u || (elig != nullptr && !elig->edge_ok(e))) continue;
        const NodeId w = g.edge_dst()[static_cast<std::size_t>(e)];
        if (is_dst[static_cast<std::size_t>(w)] != 0) reached = true;
        if (seen[static_cast<std::size_t>(w)] == 0) {
          seen[static_cast<std::size_t>(w)] = 1;
          stack.push_back(w);
        }
      }
    }
    if (reached) out.push_back(s);
  }
  return out;
}

// The pool baseline c0: the smallest static distance from a roster head to an allowed destination
// over allowed edges, every capacity (Bellman-Ford).
inline Cost pool_baseline(const Graph& g, const Eligibility* elig, const std::vector<NodeId>& heads,
                          const std::vector<NodeId>& dests) {
  const auto n = static_cast<std::size_t>(g.num_nodes());
  constexpr Cost kInf = std::numeric_limits<Cost>::max();
  Cost best = kInf;
  for (NodeId h : heads) {
    std::vector<Cost> d(n, kInf);
    d[static_cast<std::size_t>(h)] = 0;
    for (std::size_t round = 0; round < n; ++round) {
      for (EdgeId e = 0; e < g.num_edges(); ++e) {
        if (elig != nullptr && !elig->edge_ok(e)) continue;
        const auto u = static_cast<std::size_t>(g.edge_src()[static_cast<std::size_t>(e)]);
        const auto w = static_cast<std::size_t>(g.edge_dst()[static_cast<std::size_t>(e)]);
        if (d[u] == kInf) continue;
        d[w] = std::min(d[w], d[u] + g.cost()[static_cast<std::size_t>(e)]);
      }
    }
    for (NodeId t : dests) {
      if (elig == nullptr || elig->node_ok(t)) best = std::min(best, d[static_cast<std::size_t>(t)]);
    }
  }
  return best;
}

// cost <= f * base exactly, for costs in [0, 2^63) and f finite and >= 1: f = m * 2^x with an
// integer m in [2^52, 2^53) (frexp). The magnitudes are compared before any product is formed, so
// every product stays below 2^127: base * m < 2^116; for x >= 0, a positive base gives
// f * base >= 2^(52 + x), above every cost once x >= 11, and below that rhs * 2^x < 2^126; for
// x < 0 (x >= -52, as f >= 1), cost * 2^-x < 2^115.
inline bool stretch_admits(Cost cost, Cost base, double f) {
  int exp = 0;
  const double frac = std::frexp(f, &exp);  // f = frac * 2^exp, frac in [0.5, 1)
  const auto m = static_cast<std::int64_t>(std::ldexp(frac, 53));  // f = m * 2^(exp - 53)
  const int x = exp - 53;
  if (base == 0) return cost == 0;
  if (x >= 11) return true;
  const Int128 rhs = Int128(base) * m;
  auto pow2 = [](int k) {  // the portable Int128 has no signed shift
    Int128 v = 1;
    for (int i = 0; i < k; ++i) v = v * 2;
    return v;
  };
  if (x >= 0) return Int128(cost) <= rhs * pow2(x);
  return Int128(cost) * pow2(-x) <= rhs;
}

// A priority as the exact rational num / den (Headroom: U / C, 0 on an edge of capacity 0).
inline std::pair<Amount, Amount> priority(const Graph& g, EdgePriority p, const std::vector<Amount>& U, EdgeId e) {
  const auto ei = static_cast<std::size_t>(e);
  const Amount c = g.capacity_units()[ei];
  if (p == EdgePriority::Residual) return {U[ei], 1};
  if (p == EdgePriority::Capacity) return {c, 1};
  return c == 0 ? std::pair<Amount, Amount> {0, 1} : std::pair<Amount, Amount> {U[ei], c};
}
inline int compare(std::pair<Amount, Amount> a, std::pair<Amount, Amount> b) {
  const Int128 l = Int128(a.first) * b.second;
  const Int128 r = Int128(b.first) * a.second;
  return l < r ? -1 : l > r ? 1 : 0;
}

// The setup of 6.2 to 6.4 for the roster `heads` (ascending) on the view U (consumed).
inline std::vector<RefTunnel> setup(const Graph& g, const Eligibility* elig, const std::vector<NodeId>& heads,
                                    const std::vector<NodeId>& dests, const Selection& sel,
                                    const std::optional<CostCeiling>& ceiling, Amount minimum, Amount D, std::int64_t N,
                                    std::vector<Amount> U) {
  const auto n = static_cast<std::size_t>(g.num_nodes());
  std::vector<char> is_dst(n, 0);
  for (NodeId d : dests) {
    if (elig == nullptr || elig->node_ok(d)) is_dst[static_cast<std::size_t>(d)] = 1;
  }
  const Cost c0 = pool_baseline(g, elig, heads, dests);
  const auto H = static_cast<std::int64_t>(heads.size());
  std::vector<RefTunnel> out;
  for (std::int64_t j = 0; j < N; ++j) {
    for (std::int64_t i = 0; i < H; ++i) {
      RefTunnel t;
      t.head = heads[static_cast<std::size_t>(i)];
      t.slot = static_cast<std::size_t>(i);
      t.request = rr(rr(D, H, i), N, j);
      const Amount threshold = std::max(t.request, minimum);
      // The filtered member graph, then one member per (u, v) without multi_edge.
      std::vector<EdgeId> usable;
      for (EdgeId e = 0; e < g.num_edges(); ++e) {
        if (elig != nullptr && !elig->edge_ok(e)) continue;
        if (U[static_cast<std::size_t>(e)] < threshold) continue;
        usable.push_back(e);
      }
      if (!sel.multi_edge) {
        std::map<std::pair<NodeId, NodeId>, EdgeId> pick;
        for (EdgeId e : usable) {
          const auto key = std::pair {g.edge_src()[static_cast<std::size_t>(e)], g.edge_dst()[static_cast<std::size_t>(e)]};
          const auto it = pick.find(key);
          if (it == pick.end()) {
            pick.emplace(key, e);
            continue;
          }
          const EdgeId cur = it->second;
          const Cost ce = g.cost()[static_cast<std::size_t>(e)];
          const Cost cc = g.cost()[static_cast<std::size_t>(cur)];
          if (ce < cc || (ce == cc && sel.prefer && compare(priority(g, sel.priority, U, e), priority(g, sel.priority, U, cur)) > 0)) {
            it->second = e;
          }
        }
        usable.clear();
        for (const auto& [key, e] : pick) usable.push_back(e);
        std::sort(usable.begin(), usable.end());
      }
      // Every simple path from the head that ends at the first destination it meets.
      std::vector<std::vector<EdgeId>> found;
      std::vector<EdgeId> walk;
      std::vector<char> on(n, 0);
      std::function<void(NodeId)> dfs = [&](NodeId u) {
        on[static_cast<std::size_t>(u)] = 1;
        for (EdgeId e : usable) {
          if (g.edge_src()[static_cast<std::size_t>(e)] != u) continue;
          const NodeId w = g.edge_dst()[static_cast<std::size_t>(e)];
          if (on[static_cast<std::size_t>(w)] != 0) continue;
          walk.push_back(e);
          if (is_dst[static_cast<std::size_t>(w)] != 0) {
            found.push_back(walk);
          } else {
            dfs(w);
          }
          walk.pop_back();
        }
        on[static_cast<std::size_t>(u)] = 0;
      };
      dfs(t.head);
      // The single-path order: cost, the bottleneck priority (PreferHigherPriority, larger first),
      // hops, the EdgeId sequence.
      const std::vector<EdgeId>* best = nullptr;
      Cost best_cost = 0;
      std::pair<Amount, Amount> best_bottleneck {0, 1};
      for (const auto& p : found) {
        Cost cost = 0;
        for (EdgeId e : p) cost += g.cost()[static_cast<std::size_t>(e)];
        std::pair<Amount, Amount> bottleneck {std::numeric_limits<Amount>::max(), 1};
        for (EdgeId e : p) {
          const auto pr = priority(g, sel.priority, U, e);
          if (compare(pr, bottleneck) < 0) bottleneck = pr;
        }
        bool better = best == nullptr;
        if (!better && cost != best_cost) better = cost < best_cost;
        else if (!better && sel.prefer && compare(bottleneck, best_bottleneck) != 0) better = compare(bottleneck, best_bottleneck) > 0;
        else if (!better && p.size() != best->size()) better = p.size() < best->size();
        else if (!better) better = p < *best;
        if (better) {
          best = &p;
          best_cost = cost;
          best_bottleneck = bottleneck;
        }
      }
      if (best == nullptr) {
        t.status = TunnelStatus::NoBandwidth;
      } else if (ceiling.has_value() && (std::holds_alternative<Absolute>(*ceiling)
                                             ? best_cost > std::get<Absolute>(*ceiling).cost
                                             : !stretch_admits(best_cost, c0, std::get<Stretch>(*ceiling).factor))) {
        t.status = TunnelStatus::CostCeiling;
      } else {
        t.path = *best;
        t.cost = best_cost;
        for (EdgeId e : t.path) U[static_cast<std::size_t>(e)] -= t.request;
      }
      out.push_back(std::move(t));
    }
  }
  return out;
}

// The lossless bound (8.3) of the tunnels on the free capacity V: the largest total whose round
// robin over the H heads, then over each head's up tunnels in index order, loads no edge above V
// and gives no positive share to a head without an up tunnel (scanning up from 0).
inline Amount lossless_bound(std::size_t H, const std::vector<RefTunnel>& tunnels, const std::vector<Amount>& V) {
  Amount limit = 1;
  for (Amount c : V) limit += c;
  auto fits = [&](Amount x) {
    std::vector<Amount> load(V.size(), 0);
    for (std::size_t h = 0; h < H; ++h) {
      const Amount share = rr(x, static_cast<std::int64_t>(H), static_cast<std::int64_t>(h));
      std::vector<const RefTunnel*> up;
      for (const auto& t : tunnels) {
        if (t.slot == h && t.status == TunnelStatus::Up) up.push_back(&t);
      }
      if (share > 0 && up.empty()) return false;
      for (std::size_t k = 0; k < up.size(); ++k) {
        const Amount part = rr(share, static_cast<std::int64_t>(up.size()), static_cast<std::int64_t>(k));
        for (EdgeId e : up[k]->path) load[static_cast<std::size_t>(e)] += part;
      }
    }
    for (std::size_t e = 0; e < V.size(); ++e) {
      if (load[e] > V[e]) return false;
    }
    return true;
  };
  Amount x = 0;
  while (x < limit && fits(x + 1)) ++x;
  return x;
}

} // namespace t1_reference

// An independent reference of T2 (contract 7.2, 7.3, 7.7) for small graphs: every multiset of k
// simple paths per live head (a count vector x) by brute force; U = the largest
// min over used edges of floor(V / x), s* = the largest min of V / x (an exact fraction), and the
// minimum (J, Phi) over the vectors with x u <= V. It calls none of the library's flow, search or
// apportionment code (its own path enumeration, Int128 arithmetic).
namespace t2_reference {

// Every simple path from `src` that ends at the first destination it meets, over the allowed edges
// that do not leave a destination.
inline std::vector<std::vector<EdgeId>> paths(const Graph& g, const Eligibility* elig, NodeId src,
                                              const std::vector<char>& is_dst) {
  const auto n = static_cast<std::size_t>(g.num_nodes());
  std::vector<std::vector<EdgeId>> found;
  std::vector<EdgeId> walk;
  std::vector<char> on(n, 0);
  std::function<void(NodeId)> dfs = [&](NodeId u) {
    on[static_cast<std::size_t>(u)] = 1;
    for (EdgeId e = 0; e < g.num_edges(); ++e) {
      const auto ei = static_cast<std::size_t>(e);
      if (g.edge_src()[ei] != u) continue;
      const NodeId w = g.edge_dst()[ei];
      if (elig != nullptr && (!elig->edge_ok(e) || !elig->node_ok(u) || !elig->node_ok(w))) continue;
      if (on[static_cast<std::size_t>(w)] != 0) continue;
      walk.push_back(e);
      if (is_dst[static_cast<std::size_t>(w)] != 0) {
        found.push_back(walk);
      } else {
        dfs(w);
      }
      walk.pop_back();
    }
    on[static_cast<std::size_t>(u)] = 0;
  };
  if (is_dst[static_cast<std::size_t>(src)] == 0) dfs(src);
  return found;
}

struct Result {
  std::vector<NodeId> live {};
  Amount U {0};
  std::pair<Amount, Amount> fluid {0, 1};  // reduced
  std::optional<std::pair<Int128, Int128>> best {};  // the minimum (J, Phi) at u
};

// floor(2^32 (B - V + j u) / B).
inline Int128 slot_price(Amount B, Amount V, Amount j, Amount u) {
  return (Int128(B - V + j * u) * Int128(std::int64_t {1} << 32)) / Int128(B);
}

// T2 at count k for the roster `heads` (ascending) on V; `u_of(U)` gives the working size at U.
inline Result solve(const Graph& g, const Eligibility* elig, const std::vector<NodeId>& heads,
                    const std::vector<NodeId>& dests, std::int64_t k, const std::vector<Amount>& V,
                    const std::function<Amount(Amount)>& u_of, bool headroom) {
  std::vector<char> is_dst(static_cast<std::size_t>(g.num_nodes()), 0);
  for (NodeId d : dests) {
    if (elig == nullptr || elig->node_ok(d)) is_dst[static_cast<std::size_t>(d)] = 1;
  }
  Result out;
  std::vector<std::vector<std::vector<EdgeId>>> per;
  for (NodeId h : heads) {
    auto ps = paths(g, elig, h, is_dst);
    const bool live = std::any_of(ps.begin(), ps.end(), [&](const auto& p) {
      return std::all_of(p.begin(), p.end(), [&](EdgeId e) { return V[static_cast<std::size_t>(e)] >= 1; });
    });
    if (!live) continue;
    out.live.push_back(h);
    per.push_back(std::move(ps));
  }
  if (out.live.empty()) return out;
  // Every count vector: per live head a nondecreasing choice of k path indices.
  std::vector<std::vector<Amount>> xs;
  std::vector<Amount> x(static_cast<std::size_t>(g.num_edges()), 0);
  std::function<void(std::size_t, std::int64_t, std::size_t)> pick = [&](std::size_t h, std::int64_t left,
                                                                         std::size_t from) {
    if (h == per.size()) {
      xs.push_back(x);
      return;
    }
    if (left == 0) {
      pick(h + 1, k, 0);
      return;
    }
    for (std::size_t i = from; i < per[h].size(); ++i) {
      for (EdgeId e : per[h][i]) ++x[static_cast<std::size_t>(e)];
      pick(h, left - 1, i);
      for (EdgeId e : per[h][i]) --x[static_cast<std::size_t>(e)];
    }
  };
  pick(0, k, 0);
  // U and s*: the largest common size over the vectors (fractions compared exactly).
  bool have = false;
  for (const auto& c : xs) {
    Amount size = std::numeric_limits<Amount>::max();
    std::pair<Amount, Amount> frac {1, 0};  // +inf
    for (std::size_t e = 0; e < c.size(); ++e) {
      if (c[e] == 0) continue;
      size = std::min(size, V[e] / c[e]);
      if (frac.second == 0 || Int128(V[e]) * frac.second < Int128(frac.first) * c[e]) frac = {V[e], c[e]};
    }
    out.U = std::max(out.U, size);
    if (!have || Int128(frac.first) * out.fluid.second > Int128(out.fluid.first) * frac.second) {
      out.fluid = frac;
      have = true;
    }
  }
  const Amount d = std::gcd(out.fluid.first, out.fluid.second);
  out.fluid = {out.fluid.first / d, out.fluid.second / d};
  if (out.U == 0) return out;
  const Amount u = u_of(out.U);
  for (const auto& c : xs) {
    bool fits = true;
    Int128 J = 0;
    Int128 phi = 0;
    for (std::size_t e = 0; e < c.size(); ++e) {
      if (c[e] == 0) continue;
      fits = fits && c[e] * u <= V[e];
      J = J + Int128(g.cost()[e]) * c[e];
      const Amount B = g.capacity_units()[e];
      for (Amount j = 1; headroom && B > 0 && j <= c[e]; ++j) phi = phi + slot_price(B, V[e], j, u);
    }
    if (!fits) continue;
    if (!out.best.has_value() || std::pair {J, phi} < *out.best) out.best = std::pair {J, phi};
  }
  return out;
}

} // namespace t2_reference

} // namespace netgraph::core::test
