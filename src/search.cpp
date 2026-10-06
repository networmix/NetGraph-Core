/*
  The search semantics of APPROACH 6 as internal routines, and the constructor of their
  PathSelection.

  Distances are int64 path costs: every cost is positive and the sum of all costs is below
  2^62, so no simple path cost overflows and kUnreachable (INT64_MAX) never collides. Dijkstra
  uses a binary heap of (distance, node); distances are unique, so the heap's order among
  equal keys changes no output.
*/
#include "netgraph/core/detail/search.hpp"

#include <algorithm>
#include <functional>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>

#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"

namespace netgraph::core::detail {

int compare(const Priority& a, const Priority& b) noexcept {
  if (a.infinite || b.infinite) return a.infinite == b.infinite ? 0 : (a.infinite ? 1 : -1);
  const U128 left = mul_u64(static_cast<std::uint64_t>(a.num), static_cast<std::uint64_t>(b.den));
  const U128 right = mul_u64(static_cast<std::uint64_t>(b.num), static_cast<std::uint64_t>(a.den));
  if (left == right) return 0;
  return left < right ? -1 : 1;
}

Priority priority_value(EdgePriority priority, Amount residual, Amount capacity) noexcept {
  if (priority == EdgePriority::Capacity) return {capacity, 1, false};
  if (priority == EdgePriority::Residual) return {residual, 1, false};
  return {residual, capacity, false};
}

Priority PriorityValues::of(EdgeId e) const noexcept {
  const auto ei = static_cast<std::size_t>(e);
  const Amount cap = g->capacity_units()[ei];
  if (residual.empty()) return {cap, 1, false};
  return priority_value(priority, residual[ei], cap);
}

std::vector<std::uint8_t> member_mask(const Graph& g, const Eligibility* elig,
                                      const PathSelection& selection, std::span<const Amount> residual,
                                      Amount min_residual) {
  const auto m = static_cast<std::size_t>(g.num_edges());
  const auto capacity = g.capacity_units();
  std::vector<std::uint8_t> member(m, 0);
  for (std::size_t e = 0; e < m; ++e) {
    if (elig != nullptr && !elig->edge_ok(static_cast<EdgeId>(e))) continue;
    if (capacity[e] < 1) continue;
    if (!residual.empty() && residual[e] < min_residual) continue;
    member[e] = 1;
  }
  if (selection.multi_edge) return member;
  // Collapse: one member per (u, v) by select_adjacency.
  const PriorityValues prio {&g, selection.priority, residual};
  const auto value = [&prio](EdgeId e) { return prio.of(e); };
  const bool prefer = selection.tie_break == EdgeTieBreak::PreferHigherPriority;
  const auto edge_dst = g.edge_dst();
  const auto row = g.row_offsets();
  std::vector<EdgeId> best(static_cast<std::size_t>(g.num_nodes()), -1);
  std::vector<NodeId> touched;
  for (std::int32_t u = 0; u < g.num_nodes(); ++u) {
    touched.clear();
    for (auto e = row[static_cast<std::size_t>(u)]; e < row[static_cast<std::size_t>(u) + 1]; ++e) {
      const auto ei = static_cast<std::size_t>(e);
      if (member[ei] == 0) continue;
      const auto v = static_cast<std::size_t>(edge_dst[ei]);
      const EdgeId cur = best[v];
      if (cur < 0) {
        best[v] = e;
        touched.push_back(static_cast<NodeId>(v));
        continue;
      }
      if (replaces(g, e, cur, prefer, value)) best[v] = e;
    }
    for (std::size_t e = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
         e < static_cast<std::size_t>(row[static_cast<std::size_t>(u) + 1]); ++e) {
      if (member[e] != 0 && best[static_cast<std::size_t>(edge_dst[e])] != static_cast<EdgeId>(e)) member[e] = 0;
    }
    for (NodeId v : touched) best[static_cast<std::size_t>(v)] = -1;
  }
  return member;
}

namespace {

using HeapEntry = std::pair<Cost, NodeId>;
using MinHeap = std::priority_queue<HeapEntry, std::vector<HeapEntry>, std::greater<>>;

} // namespace

std::vector<Cost> distances(const Graph& g, std::span<const std::uint8_t> member,
                            std::span<const NodeId> roots, Direction direction) {
  const auto n = static_cast<std::size_t>(g.num_nodes());
  std::vector<Cost> dist(n, kUnreachable);
  const auto cost = g.cost();
  const auto edge_dst = g.edge_dst();
  const auto edge_src = g.edge_src();
  const auto row = g.row_offsets();
  const auto in_row = g.in_row_offsets();
  const auto in_adj = g.in_adj_edge_index();
  MinHeap heap;
  for (NodeId r : roots) {
    dist[static_cast<std::size_t>(r)] = 0;
    heap.emplace(0, r);
  }
  while (!heap.empty()) {
    const auto [d, u] = heap.top();
    heap.pop();
    const auto ui = static_cast<std::size_t>(u);
    if (d != dist[ui]) continue;
    if (direction == Direction::Forward) {
      for (auto e = row[ui]; e < row[ui + 1]; ++e) {
        const auto ei = static_cast<std::size_t>(e);
        if (member[ei] == 0) continue;
        const auto v = static_cast<std::size_t>(edge_dst[ei]);
        const Cost nd = d + cost[ei];
        if (nd < dist[v]) {
          dist[v] = nd;
          heap.emplace(nd, static_cast<NodeId>(v));
        }
      }
    } else {
      for (auto s = in_row[ui]; s < in_row[ui + 1]; ++s) {
        const auto ei = static_cast<std::size_t>(in_adj[static_cast<std::size_t>(s)]);
        if (member[ei] == 0) continue;
        const auto v = static_cast<std::size_t>(edge_src[ei]);
        const Cost nd = d + cost[ei];
        if (nd < dist[v]) {
          dist[v] = nd;
          heap.emplace(nd, static_cast<NodeId>(v));
        }
      }
    }
  }
  return dist;
}

std::vector<EdgeId> full_dag(const Graph& g, std::span<const std::uint8_t> member,
                             std::span<const Cost> dist, Direction direction) {
  const auto cost = g.cost();
  const auto edge_src = g.edge_src();
  const auto edge_dst = g.edge_dst();
  std::vector<EdgeId> out;
  for (std::size_t e = 0; e < member.size(); ++e) {
    if (member[e] == 0) continue;
    const Cost du = dist[static_cast<std::size_t>(edge_src[e])];
    const Cost dv = dist[static_cast<std::size_t>(edge_dst[e])];
    if (du == kUnreachable || dv == kUnreachable) continue;
    const bool tight = direction == Direction::Forward ? du + cost[e] == dv : dv + cost[e] == du;
    if (tight) out.push_back(static_cast<EdgeId>(e));
  }
  return out;
}

TightDag tight_dag(const Graph& g, std::span<const std::uint8_t> member, std::span<const NodeId> sources,
                   std::span<const NodeId> destinations) {
  const auto n = static_cast<std::size_t>(g.num_nodes());
  const auto cost = g.cost();
  const auto edge_dst = g.edge_dst();
  const auto edge_src = g.edge_src();
  const auto row = g.row_offsets();
  const auto in_row = g.in_row_offsets();
  const auto in_adj = g.in_adj_edge_index();
  std::vector<std::uint8_t> is_dst(n, 0);
  for (NodeId t : destinations) is_dst[static_cast<std::size_t>(t)] = 1;
  // Forward Dijkstra from the sources until every node at distance <= d* is final.
  std::vector<Cost> dist(n, kUnreachable);
  std::vector<std::uint8_t> settled(n, 0);
  MinHeap heap;
  for (NodeId s : sources) {
    dist[static_cast<std::size_t>(s)] = 0;
    heap.emplace(0, s);
  }
  Cost dstar = kUnreachable;
  while (!heap.empty()) {
    const auto [d, u] = heap.top();
    if (d > dstar) break;
    heap.pop();
    const auto ui = static_cast<std::size_t>(u);
    if (d != dist[ui] || settled[ui] != 0) continue;
    settled[ui] = 1;
    if (is_dst[ui] != 0) {
      dstar = std::min(dstar, d);
      continue;  // destinations are terminal
    }
    for (auto e = row[ui]; e < row[ui + 1]; ++e) {
      const auto ei = static_cast<std::size_t>(e);
      if (member[ei] == 0) continue;
      const auto v = static_cast<std::size_t>(edge_dst[ei]);
      const Cost nd = d + cost[ei];
      if (nd < dist[v]) {
        dist[v] = nd;
        heap.emplace(nd, static_cast<NodeId>(v));
      }
    }
  }
  TightDag out;
  if (dstar == kUnreachable) {
    out.cost = kUnreachable;
    return out;
  }
  out.cost = dstar;
  // Back-prune from the destinations at d* over tight in-edges whose tail is settled below d*.
  std::vector<std::uint8_t> keep(n, 0);
  std::vector<NodeId> stack;
  for (NodeId t : destinations) {
    const auto ti = static_cast<std::size_t>(t);
    if (dist[ti] == dstar && keep[ti] == 0) {
      keep[ti] = 1;
      stack.push_back(t);
    }
  }
  while (!stack.empty()) {
    const auto v = static_cast<std::size_t>(stack.back());
    stack.pop_back();
    for (auto s = in_row[v]; s < in_row[v + 1]; ++s) {
      const auto ei = static_cast<std::size_t>(in_adj[static_cast<std::size_t>(s)]);
      if (member[ei] == 0) continue;
      const auto u = static_cast<std::size_t>(edge_src[ei]);
      if (settled[u] == 0 || is_dst[u] != 0 || dist[u] >= dstar || dist[u] + cost[ei] != dist[v]) continue;
      out.edges.push_back(static_cast<EdgeId>(ei));
      if (keep[u] == 0) {
        keep[u] = 1;
        stack.push_back(static_cast<NodeId>(u));
      }
    }
  }
  std::sort(out.edges.begin(), out.edges.end());
  return out;
}

std::optional<SinglePath> canonical_path(const Graph& g, std::span<const std::uint8_t> member,
                                         std::span<const NodeId> sources, std::span<const NodeId> destinations,
                                         const PriorityValues* priority) {
  TightDag tight = tight_dag(g, member, sources, destinations);
  if (tight.cost == kUnreachable) return std::nullopt;
  const auto n = static_cast<std::size_t>(g.num_nodes());
  const auto cost = g.cost();
  const auto edge_src = g.edge_src();
  const auto edge_dst = g.edge_dst();
  std::vector<std::uint8_t> is_dst(n, 0);
  for (NodeId t : destinations) is_dst[static_cast<std::size_t>(t)] = 1;
  std::vector<EdgeId> dag = std::move(tight.edges);
  SinglePath out;
  out.cost = tight.cost;
  if (priority != nullptr) {
    // Stage 2: B*, the max-min bottleneck over the walks of the tight DAG. Distances from the
    // sources inside the DAG order its nodes topologically.
    std::vector<std::uint8_t> dag_member(member.size(), 0);
    for (EdgeId e : dag) dag_member[static_cast<std::size_t>(e)] = 1;
    const auto dist = distances(g, dag_member, sources, Direction::Forward);
    std::vector<EdgeId> by_tail(dag);
    std::sort(by_tail.begin(), by_tail.end(), [&](EdgeId a, EdgeId b) {
      const Cost da = dist[static_cast<std::size_t>(edge_src[static_cast<std::size_t>(a)])];
      const Cost db = dist[static_cast<std::size_t>(edge_src[static_cast<std::size_t>(b)])];
      return da != db ? da < db : a < b;
    });
    std::vector<std::optional<Priority>> best(n);
    for (NodeId s : sources) best[static_cast<std::size_t>(s)] = Priority {0, 1, true};
    for (EdgeId e : by_tail) {
      const auto u = static_cast<std::size_t>(edge_src[static_cast<std::size_t>(e)]);
      const auto v = static_cast<std::size_t>(edge_dst[static_cast<std::size_t>(e)]);
      if (!best[u].has_value()) continue;
      Priority cand = priority->of(e);
      if (compare(*best[u], cand) < 0) cand = *best[u];
      if (!best[v].has_value() || compare(cand, *best[v]) > 0) best[v] = cand;
    }
    std::optional<Priority> b_star;
    for (NodeId t : destinations) {
      const auto ti = static_cast<std::size_t>(t);
      if (dist[ti] != tight.cost || !best[ti].has_value()) continue;
      if (!b_star.has_value() || compare(*best[ti], *b_star) > 0) b_star = best[ti];
    }
    if (!b_star.has_value()) return std::nullopt;
    std::vector<EdgeId> kept;
    for (EdgeId e : dag) {
      if (compare(priority->of(e), *b_star) >= 0) kept.push_back(e);
    }
    dag = std::move(kept);
  }
  // Stage 3: hop counts toward the destinations inside the (restricted) DAG, then the greedy
  // walk by the smallest tight edge whose head is one hop closer.
  constexpr std::int32_t kNone = -1;
  std::vector<std::int32_t> h(n, kNone);
  std::vector<std::vector<EdgeId>> in_dag(n), out_dag(n);
  for (EdgeId e : dag) {  // ascending: out_dag lists stay ascending
    in_dag[static_cast<std::size_t>(edge_dst[static_cast<std::size_t>(e)])].push_back(e);
    out_dag[static_cast<std::size_t>(edge_src[static_cast<std::size_t>(e)])].push_back(e);
  }
  std::vector<NodeId> queue;
  for (NodeId t : destinations) {
    h[static_cast<std::size_t>(t)] = 0;
    queue.push_back(t);
  }
  for (std::size_t qi = 0; qi < queue.size(); ++qi) {
    const auto v = static_cast<std::size_t>(queue[qi]);
    for (EdgeId e : in_dag[v]) {
      const auto u = static_cast<std::size_t>(edge_src[static_cast<std::size_t>(e)]);
      if (h[u] == kNone) {
        h[u] = h[v] + 1;
        queue.push_back(static_cast<NodeId>(u));
      }
    }
  }
  std::int32_t hmin = kNone;
  for (NodeId s : sources) {
    const std::int32_t hs = h[static_cast<std::size_t>(s)];
    if (hs != kNone && (hmin == kNone || hs < hmin)) hmin = hs;
  }
  if (hmin == kNone) return std::nullopt;
  std::vector<NodeId> frontier;
  for (NodeId s : sources) {
    if (h[static_cast<std::size_t>(s)] == hmin) frontier.push_back(s);
  }
  while (true) {
    EdgeId pick = -1;
    for (NodeId x : frontier) {
      const auto xi = static_cast<std::size_t>(x);
      for (EdgeId e : out_dag[xi]) {
        const auto v = static_cast<std::size_t>(edge_dst[static_cast<std::size_t>(e)]);
        if (h[v] != kNone && h[v] == h[xi] - 1) {
          if (pick < 0 || e < pick) pick = e;
          break;  // ascending: the first qualifying edge of x is its smallest
        }
      }
    }
    if (pick < 0) break;
    out.walk.push_back(pick);
    const NodeId u = edge_dst[static_cast<std::size_t>(pick)];
    if (is_dst[static_cast<std::size_t>(u)] != 0) break;
    frontier.assign(1, u);
  }
  if (out.walk.empty()) return std::nullopt;
  Cost walk_cost = 0;
  for (EdgeId e : out.walk) walk_cost += cost[static_cast<std::size_t>(e)];
  if (walk_cost != out.cost) throw std::logic_error("canonical path: the walk is not tight");
  return out;
}

SpfOut spf(const Graph& g, std::span<const NodeId> roots, Direction direction,
           std::optional<std::span<const NodeId>> targets, const PathSelection& selection,
           std::optional<std::span<const Amount>> residual_arg, Amount min_residual, const Eligibility* elig) {
  check_member(direction, "direction");
  check_selection(selection);
  const auto root_set = node_set(roots, g.num_nodes(), "roots");
  std::vector<NodeId> target_set;
  if (targets.has_value()) {
    target_set = node_set(*targets, g.num_nodes(), "targets");
    check_disjoint(root_set, target_set, "roots", "targets");
  } else if (!selection.multipath) {
    throw std::invalid_argument("multipath = false requires targets (the single path joins a root to a target)");
  }
  const auto m = static_cast<std::size_t>(g.num_edges());
  const std::span<const Amount> residual = residual_arg.value_or(std::span<const Amount> {});
  if (residual_arg.has_value() && residual.size() != m) {
    throw std::invalid_argument("residual has length " + std::to_string(residual.size()) + ", expected " +
                                std::to_string(m) + " (one entry per edge)");
  }
  const auto capacity = g.capacity_units();
  for (std::size_t e = 0; e < residual.size(); ++e) {
    if (residual[e] < 0 || residual[e] > capacity[e]) {
      throw std::invalid_argument("residual[" + std::to_string(e) + "] = " + std::to_string(residual[e]) +
                                  " is outside [0, capacity_units[" + std::to_string(e) + "] = " +
                                  std::to_string(capacity[e]) + "]");
    }
  }
  if (min_residual < 1) {
    throw std::invalid_argument("min_residual must be at least 1 unit (got " + std::to_string(min_residual) + ")");
  }
  if (elig != nullptr && elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  const auto member = member_mask(g, elig, selection, residual, min_residual);
  std::vector<NodeId> live_roots;
  for (NodeId r : root_set) {
    if (elig == nullptr || elig->node_ok(r)) live_roots.push_back(r);
  }
  SpfOut out;
  if (!targets.has_value()) {
    auto dist = distances(g, member, live_roots, direction);
    out.edges = full_dag(g, member, dist, direction);
    out.costs = std::move(dist);
    return out;
  }
  std::vector<NodeId> live_targets;
  for (NodeId t : target_set) {
    if (elig == nullptr || elig->node_ok(t)) live_targets.push_back(t);
  }
  out.target_cost = std::nullopt;
  if (live_roots.empty() || live_targets.empty()) return out;
  const auto& src = direction == Direction::Forward ? live_roots : live_targets;
  const auto& dst = direction == Direction::Forward ? live_targets : live_roots;
  if (selection.multipath) {
    auto tight = tight_dag(g, member, src, dst);
    if (tight.cost != kUnreachable) {
      out.target_cost = tight.cost;
      out.edges = std::move(tight.edges);
    }
    return out;
  }
  const PriorityValues prio {&g, selection.priority, residual};
  const bool prefer = selection.tie_break == EdgeTieBreak::PreferHigherPriority;
  auto path = canonical_path(g, member, src, dst, prefer ? &prio : nullptr);
  if (path.has_value()) {
    out.target_cost = path->cost;
    out.edges = path->walk;
    std::sort(out.edges.begin(), out.edges.end());
    out.path = std::move(path);
  }
  return out;
}

PredDAG static_support(const Graph& g, std::span<const NodeId> roots, Direction direction,
                       const Eligibility* elig) {
  const PathSelection full(EdgeTieBreak::Deterministic, EdgePriority::Capacity, true, true);
  const auto out = spf(g, roots, direction, std::nullopt, full, std::nullopt, 1, elig);
  return PredDAGBuilder::of_edges(g, out.edges);
}

} // namespace netgraph::core::detail

namespace netgraph::core {

// The selection of a search (options.hpp): its two enum-typed fields inside their enumerations.
PathSelection::PathSelection(EdgeTieBreak tie_break_, EdgePriority priority_, bool multipath_, bool multi_edge_)
    : tie_break(tie_break_), priority(priority_), multipath(multipath_), multi_edge(multi_edge_) {
  detail::check_member(tie_break, "tie_break");
  detail::check_member(priority, "priority");
}

} // namespace netgraph::core
