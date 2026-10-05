#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <vector>

#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/strict_multidigraph.hpp"

namespace netgraph::core::test {

// An edge in input coordinates.
struct InEdge {
  NodeId u;
  NodeId v;
  Amount cap;
  Cost cost;
};

inline StrictMultiDiGraph make_graph(std::int32_t n, const std::vector<InEdge>& edges, double q = 1.0) {
  std::vector<NodeId> src, dst;
  std::vector<Amount> cap;
  std::vector<Cost> cost;
  for (const auto& e : edges) {
    src.push_back(e.u);
    dst.push_back(e.v);
    cap.push_back(e.cap);
    cost.push_back(e.cost);
  }
  return StrictMultiDiGraph(n, src, dst, cap, cost, q);
}

// The caller's inverse of input_index(): input position -> EdgeId.
inline std::vector<EdgeId> id_of_input(const StrictMultiDiGraph& g) {
  auto ii = g.input_index();
  std::vector<EdgeId> inv(ii.size());
  for (std::size_t e = 0; e < ii.size(); ++e) inv[static_cast<std::size_t>(ii[e])] = static_cast<EdgeId>(e);
  return inv;
}

// A per-edge array in input order, permuted to id order.
template <typename T>
std::vector<T> to_ids(const StrictMultiDiGraph& g, const std::vector<T>& by_input) {
  auto ii = g.input_index();
  std::vector<T> out(ii.size());
  for (std::size_t e = 0; e < ii.size(); ++e) out[e] = by_input[static_cast<std::size_t>(ii[e])];
  return out;
}

// A per-edge array in id order, permuted to input order.
template <typename T>
std::vector<T> to_input(const StrictMultiDiGraph& g, std::span<const T> by_id) {
  auto ii = g.input_index();
  std::vector<T> out(ii.size());
  for (std::size_t e = 0; e < ii.size(); ++e) out[static_cast<std::size_t>(ii[e])] = by_id[e];
  return out;
}

template <typename T>
std::vector<T> to_input(const StrictMultiDiGraph& g, const std::vector<T>& by_id) {
  return to_input(g, std::span<const T>(by_id));
}

// EdgeIds translated to sorted input positions.
inline std::vector<std::int32_t> ids_to_input(const StrictMultiDiGraph& g, const std::vector<EdgeId>& ids) {
  std::vector<std::int32_t> out;
  for (EdgeId e : ids) out.push_back(g.input_index()[static_cast<std::size_t>(e)]);
  std::sort(out.begin(), out.end());
  return out;
}

inline std::vector<std::pair<Cost, Amount>> curve(std::initializer_list<std::pair<Cost, Amount>> c) { return c; }

inline FlowCoreResult solve(const StrictMultiDiGraph& g, std::vector<NodeId> s, std::vector<NodeId> d,
                            std::optional<Amount> cap, const std::vector<Amount>* view = nullptr,
                            const Eligibility* elig = nullptr) {
  std::vector<Amount> full(g.capacity_units().begin(), g.capacity_units().end());
  return place_global_units(g, view ? std::span<const Amount>(*view) : std::span<const Amount>(full), s, d, elig, cap);
}

inline bool same_payload(const FlowCoreResult& a, const FlowCoreResult& b) {
  const auto& x = a.counters;
  const auto& y = b.counters;
  return a.placed == b.placed && a.cost_distribution == b.cost_distribution && a.edge_flows == b.edge_flows &&
         a.min_cut == b.min_cut && a.potentials == b.potentials && x.dijkstra_runs == y.dijkstra_runs &&
         x.dijkstra_pops == y.dijkstra_pops && x.settled == y.settled && x.arc_scans == y.arc_scans &&
         x.dinic_phases == y.dinic_phases && x.augmentations == y.augmentations;
}

// Value -> minimum cost of every integer flow value, by enumerating every integer flow on the
// eligible edges (small graphs only). Independent of the solver.
inline std::map<Amount, __int128> exhaustive_costs(const StrictMultiDiGraph& g, std::span<const Amount> view,
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
  std::map<Amount, __int128> best;
  std::vector<Amount> f(free.size(), 0);
  while (true) {
    std::vector<Amount> bal(n, 0);
    __int128 cost = 0;
    for (std::size_t k = 0; k < free.size(); ++k) {
      const auto e = free[k];
      bal[static_cast<std::size_t>(g.edge_src()[e])] -= f[k];
      bal[static_cast<std::size_t>(g.edge_dst()[e])] += f[k];
      cost += static_cast<__int128>(g.cost()[e]) * f[k];
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
inline std::vector<std::pair<Amount, __int128>> independent_curve(const StrictMultiDiGraph& g,
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
  std::vector<std::pair<Amount, __int128>> pts{{0, 0}};
  while (true) {
    const __int128 inf = std::numeric_limits<__int128>::max();
    std::vector<__int128> dist(static_cast<std::size_t>(n + 2), inf);
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
// value is exact in 128 bits: costs stay below 2^62 per unit and amounts below 2^62.
inline __int128 curve_at(const std::vector<std::pair<Amount, __int128>>& pts, Amount x) {
  for (std::size_t i = 1; i < pts.size(); ++i) {
    const auto [a, ca] = pts[i - 1];
    const auto [b, cb] = pts[i];
    if (a <= x && x <= b) {
      const __int128 slope = (cb - ca) / (b - a);
      return ca + slope * (x - a);
    }
  }
  return pts.back().second;
}

inline std::vector<std::pair<Amount, __int128>> payload_points(const FlowCoreResult& r) {
  std::vector<std::pair<Amount, __int128>> pts{{0, 0}};
  for (const auto& [c, u] : r.cost_distribution) {
    pts.emplace_back(pts.back().first + u, pts.back().second + static_cast<__int128>(c) * u);
  }
  return pts;
}

// The certificate of APPROACH 2.2 in two ordered branches. Returns the failed clauses. With
// `tier_a`, the complete curve is compared with the exhaustive optimum of every value; else with
// the independent Bellman-Ford curve at the union of both breakpoint sets.
inline std::vector<std::string> certify(const StrictMultiDiGraph& g, std::span<const Amount> view,
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
    const auto& c = r.counters;
    if (r.placed != 0 || !r.cost_distribution.empty() || r.min_cut.has_value() ||
        std::any_of(r.edge_flows.begin(), r.edge_flows.end(), [](Amount x) { return x != 0; }) ||
        std::any_of(r.potentials.begin(), r.potentials.end(), [](Cost p) { return p != kNeverReached; }) ||
        c.dijkstra_runs || c.dijkstra_pops || c.settled || c.arc_scans || c.dinic_phases || c.augmentations) {
      fails.push_back("zero cap: the prescribed zero payload");
    }
    return fails;
  }
  std::vector<__int128> bal(n, 0);
  __int128 moment = 0;
  for (std::size_t e = 0; e < m; ++e) {
    const Amount f = r.edge_flows[e];
    if (f < 0 || f > view[e]) fails.push_back("feasibility at edge " + std::to_string(e));
    if (!edge_ok(e) && f != 0) fails.push_back("flow on an ineligible edge " + std::to_string(e));
    bal[static_cast<std::size_t>(g.edge_src()[e])] -= f;
    bal[static_cast<std::size_t>(g.edge_dst()[e])] += f;
    moment += static_cast<__int128>(g.cost()[e]) * f;
  }
  std::set<NodeId> S(srcs.begin(), srcs.end()), D(dsts.begin(), dsts.end());
  __int128 supply = 0, inflow = 0;
  std::vector<__int128> sup(n, 0), inf_d(n, 0);
  for (std::size_t v = 0; v < n; ++v) {
    const auto node = static_cast<NodeId>(v);
    if (S.count(node)) {
      if (bal[v] > 0) fails.push_back("source absorbs flow");
      if (!node_ok(node) && bal[v] != 0) fails.push_back("excluded source supplies");
      sup[v] = -bal[v];
      supply += -bal[v];
    } else if (D.count(node)) {
      if (bal[v] < 0) fails.push_back("destination emits flow");
      if (!node_ok(node) && bal[v] != 0) fails.push_back("excluded destination receives");
      inf_d[v] = bal[v];
      inflow += bal[v];
    } else if (bal[v] != 0) {
      fails.push_back("conservation at node " + std::to_string(v));
    }
  }
  if (supply != r.placed || inflow != r.placed) fails.push_back("value differs from supply or inflow");
  if (r.placed < 0 || (cap.has_value() && r.placed > *cap)) fails.push_back("placed outside [0, cap]");
  // Full residual closure seeded at the allowed real sources.
  const std::size_t s_star = n, t_star = n + 1;
  std::vector<std::vector<std::pair<std::size_t, Cost>>> out(n + 2);
  for (std::size_t e = 0; e < m; ++e) {
    if (!edge_ok(e)) continue;
    const auto u = static_cast<std::size_t>(g.edge_src()[e]);
    const auto v = static_cast<std::size_t>(g.edge_dst()[e]);
    if (r.edge_flows[e] < view[e]) out[u].emplace_back(v, g.cost()[e]);
    if (r.edge_flows[e] > 0) out[v].emplace_back(u, -g.cost()[e]);
  }
  for (NodeId s : srcs) {
    if (!node_ok(s)) continue;
    out[s_star].emplace_back(static_cast<std::size_t>(s), 0);
    if (sup[static_cast<std::size_t>(s)] > 0) out[static_cast<std::size_t>(s)].emplace_back(s_star, 0);
  }
  for (NodeId d : dsts) {
    if (!node_ok(d)) continue;
    out[static_cast<std::size_t>(d)].emplace_back(t_star, 0);
    if (inf_d[static_cast<std::size_t>(d)] > 0) out[t_star].emplace_back(static_cast<std::size_t>(d), 0);
  }
  std::vector<char> seen(n + 2, 0);
  std::vector<std::size_t> queue;
  for (NodeId s : srcs) {
    if (node_ok(s)) { seen[static_cast<std::size_t>(s)] = 1; queue.push_back(static_cast<std::size_t>(s)); }
  }
  for (std::size_t qi = 0; qi < queue.size(); ++qi) {
    for (auto [v, w] : out[queue[qi]]) {
      (void)w;
      if (!seen[v]) { seen[v] = 1; queue.push_back(v); }
    }
  }
  const bool connected = seen[t_star] != 0;
  if (!cap.has_value() && connected) fails.push_back("uncapped solve left the destinations connected");
  if (cap.has_value() && r.placed != *cap && connected) fails.push_back("attained value: early stop");
  if (connected) {
    if (r.min_cut.has_value()) fails.push_back("min_cut present on a connected solve");
  } else {
    std::vector<EdgeId> expect;
    __int128 cut_value = 0;
    for (std::size_t e = 0; e < m; ++e) {
      if (!edge_ok(e)) continue;
      if (seen[static_cast<std::size_t>(g.edge_src()[e])] && !seen[static_cast<std::size_t>(g.edge_dst()[e])] &&
          view[e] > 0) {
        expect.push_back(static_cast<EdgeId>(e));
        cut_value += view[e];
      }
    }
    if (!r.min_cut.has_value()) fails.push_back("min_cut missing on a maximal solve");
    else if (*r.min_cut != expect) fails.push_back("min_cut differs from the closure cut");
    if (cut_value != r.placed) fails.push_back("cut value differs from placed");
  }
  for (std::size_t v = 0; v < n; ++v) {
    if (r.potentials[v] < kNeverReached) fails.push_back("potential below the sentinel");
    if (seen[v] && r.potentials[v] == kNeverReached) fails.push_back("closure vertex holds NEVER_REACHED");
  }
  for (NodeId s : srcs) {
    const auto p = r.potentials[static_cast<std::size_t>(s)];
    if (node_ok(s) && p != 0) fails.push_back("allowed root with potential other than 0");
  }
  for (std::size_t u = 0; u < n; ++u) {
    if (!seen[u] || r.potentials[u] == kNeverReached) continue;
    for (auto [v, w] : out[u]) {
      if (v >= n || !seen[v] || r.potentials[v] == kNeverReached) continue;
      if (static_cast<__int128>(w) + r.potentials[u] - r.potentials[v] < 0) {
        fails.push_back("negative reduced cost inside the closure");
      }
    }
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
      if (inf_d[static_cast<std::size_t>(d)] > 0 && p != *pt) fails.push_back("inflow destination potential differs from pi(t*)");
    }
    if (!r.cost_distribution.empty() && r.cost_distribution.back().first != *pt) {
      fails.push_back("pi(t*) differs from the last curve cost");
    }
  } else if (r.placed != 0) {
    fails.push_back("flow delivered but no destination reached");
  }
  __int128 total = 0, curve_moment = 0;
  for (std::size_t i = 0; i < r.cost_distribution.size(); ++i) {
    const auto [c, u] = r.cost_distribution[i];
    if (u <= 0) fails.push_back("curve volume not positive");
    if (i > 0 && c <= r.cost_distribution[i - 1].first) fails.push_back("curve costs not strictly increasing");
    total += u;
    curve_moment += static_cast<__int128>(c) * u;
  }
  if (total != r.placed) fails.push_back("curve total differs from placed");
  if (curve_moment != moment) fails.push_back("curve moment differs from the flow cost");
  const auto pts = payload_points(r);
  if (tier_a) {
    const auto best = exhaustive_costs(g, view, srcs, dsts, elig);
    const Amount fmax = best.rbegin()->first;
    const Amount top = cap.has_value() ? std::min(*cap, fmax) : fmax;
    if (r.placed != top) fails.push_back("placed differs from min(cap, F_max)");
    for (Amount x = 0; x <= std::min(top, r.placed); ++x) {
      auto it = best.find(x);
      if (it == best.end() || curve_at(pts, x) != it->second) {
        fails.push_back("curve differs from the exhaustive optimum at " + std::to_string(x));
        break;
      }
    }
  } else {
    const auto ref = independent_curve(g, view, srcs, dsts, elig);
    const Amount fmax = ref.back().first;
    const Amount top = cap.has_value() ? std::min(*cap, fmax) : fmax;
    if (r.placed != top) fails.push_back("placed differs from min(cap, F_max)");
    std::set<Amount> xs;
    for (auto [x, c] : pts) { (void)c; xs.insert(x); }
    for (auto [x, c] : ref) { (void)c; xs.insert(x); }
    for (Amount x : xs) {
      if (x > std::min(top, r.placed)) continue;
      if (curve_at(pts, x) != curve_at(ref, x)) {
        fails.push_back("curve differs from the independent solver at " + std::to_string(x));
        break;
      }
    }
  }
  return fails;
}

} // namespace netgraph::core::test
