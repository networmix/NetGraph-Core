// An independent reference of the search semantics (docs/semantics.md#4-search) for the tests: the member graph
// by an explicit per-pair selection, Bellman-Ford distances, the tight-edge DAG and the pruned
// DAG by two-sided distances, and the canonical single path by enumerating every simple walk.
// It shares no code with the traversal core (no queue, no lists, no extractor).
#pragma once

#include <algorithm>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/types.hpp"
#include "test_int128.hpp"

namespace netgraph::core::test::ref {

// An exact priority: residual / capacity as a rational, or an integer (den 1); `inf` above all.
struct Prio {
  Int128 num {0};
  Int128 den {1};
  bool inf {false};
};
inline int cmp(const Prio& a, const Prio& b) {
  if (a.inf || b.inf) return a.inf == b.inf ? 0 : (a.inf ? 1 : -1);
  const Int128 l = a.num * b.den;
  const Int128 r = b.num * a.den;
  return l == r ? 0 : (l < r ? -1 : 1);
}

struct Query {
  const Graph* g {nullptr};
  const Eligibility* elig {nullptr};
  const std::vector<Amount>* residual {nullptr};
  Amount min_residual {1};
  bool multi_edge {true};
  bool prefer {false};
  EdgePriority priority {EdgePriority::Capacity};
};

inline Prio prio_of(const Query& q, EdgeId e) {
  const auto ei = static_cast<std::size_t>(e);
  const Amount cap = q.g->capacity_units()[ei];
  if (q.residual == nullptr || q.priority == EdgePriority::Capacity) return {cap, 1, false};
  const Amount r = (*q.residual)[ei];
  if (q.priority == EdgePriority::Residual) return {r, 1, false};
  return {r, cap, false};
}

// The member flags per edge.
inline std::vector<char> members(const Query& q) {
  const Graph& g = *q.g;
  const auto m = static_cast<std::size_t>(g.num_edges());
  std::vector<char> cand(m, 0);
  for (std::size_t e = 0; e < m; ++e) {
    if (q.elig != nullptr && !q.elig->edge_ok(static_cast<EdgeId>(e))) continue;
    if (g.capacity_units()[e] < 1) continue;
    if (q.residual != nullptr && (*q.residual)[e] < q.min_residual) continue;
    cand[e] = 1;
  }
  if (q.multi_edge) return cand;
  // select_adjacency: the minimum cost, then the higher priority under prefer, then the lower id.
  auto better = [&](EdgeId a, EdgeId b) {
    const Cost ca = g.cost()[static_cast<std::size_t>(a)];
    const Cost cb = g.cost()[static_cast<std::size_t>(b)];
    if (ca != cb) return ca < cb;
    if (q.prefer) {
      const int c = cmp(prio_of(q, a), prio_of(q, b));
      if (c != 0) return c > 0;
    }
    return a < b;
  };
  std::map<std::pair<NodeId, NodeId>, EdgeId> best;
  for (std::size_t e = 0; e < m; ++e) {
    if (cand[e] == 0) continue;
    const auto key = std::make_pair(g.edge_src()[e], g.edge_dst()[e]);
    auto it = best.find(key);
    if (it == best.end()) {
      best[key] = static_cast<EdgeId>(e);
    } else if (better(static_cast<EdgeId>(e), it->second)) {
      it->second = static_cast<EdgeId>(e);
    }
  }
  std::vector<char> out(m, 0);
  for (const auto& [k, e] : best) out[static_cast<std::size_t>(e)] = 1;
  return out;
}

constexpr Cost kNone = kUnreachable;

// Bellman-Ford over the members: from the set (forward) or toward it (reverse).
inline std::vector<Cost> distances(const Graph& g, const std::vector<char>& mem, const std::vector<NodeId>& set,
                                   bool reverse) {
  std::vector<Cost> d(static_cast<std::size_t>(g.num_nodes()), kNone);
  for (NodeId v : set) d[static_cast<std::size_t>(v)] = 0;
  for (bool changed = true; changed;) {
    changed = false;
    for (std::size_t e = 0; e < mem.size(); ++e) {
      if (mem[e] == 0) continue;
      const auto a = static_cast<std::size_t>(reverse ? g.edge_dst()[e] : g.edge_src()[e]);
      const auto b = static_cast<std::size_t>(reverse ? g.edge_src()[e] : g.edge_dst()[e]);
      if (d[a] == kNone) continue;
      const Cost nd = d[a] + g.cost()[e];
      if (nd < d[b]) {
        d[b] = nd;
        changed = true;
      }
    }
  }
  return d;
}

// The full DAG of a search from (forward) or toward (reverse) the roots: tight member edges.
inline std::set<EdgeId> full_dag(const Graph& g, const std::vector<char>& mem, const std::vector<Cost>& d, bool reverse) {
  std::set<EdgeId> out;
  for (std::size_t e = 0; e < mem.size(); ++e) {
    if (mem[e] == 0) continue;
    const Cost du = d[static_cast<std::size_t>(g.edge_src()[e])];
    const Cost dv = d[static_cast<std::size_t>(g.edge_dst()[e])];
    if (du == kNone || dv == kNone) continue;
    if (reverse ? dv + g.cost()[e] == du : du + g.cost()[e] == dv) out.insert(static_cast<EdgeId>(e));
  }
  return out;
}

// The walks from the start set to the end set of minimum cost d*: (d*, their edges); d* kNone
// when no end is reachable.
inline std::pair<Cost, std::set<EdgeId>> pruned(const Graph& g, const std::vector<char>& mem,
                                                const std::vector<NodeId>& starts, const std::vector<NodeId>& ends) {
  const auto ds = distances(g, mem, starts, false);
  const auto de = distances(g, mem, ends, true);
  Cost dstar = kNone;
  for (NodeId t : ends) dstar = std::min(dstar, ds[static_cast<std::size_t>(t)]);
  std::set<EdgeId> out;
  if (dstar == kNone) return {dstar, out};
  for (std::size_t e = 0; e < mem.size(); ++e) {
    if (mem[e] == 0) continue;
    const Cost du = ds[static_cast<std::size_t>(g.edge_src()[e])];
    const Cost dv = de[static_cast<std::size_t>(g.edge_dst()[e])];
    if (du != kNone && dv != kNone && du + g.cost()[e] + dv == dstar) out.insert(static_cast<EdgeId>(e));
  }
  return {dstar, out};
}

// Every simple walk from a start to an end (stopping at the first end entered) over the members.
inline void walks(const Graph& g, const std::vector<char>& mem, NodeId u, const std::vector<char>& is_end,
                  std::vector<EdgeId>& walk, std::vector<char>& on, std::vector<std::vector<EdgeId>>& out) {
  if (is_end[static_cast<std::size_t>(u)] != 0) {
    out.push_back(walk);
    return;
  }
  for (std::size_t e = 0; e < mem.size(); ++e) {
    if (mem[e] == 0 || g.edge_src()[e] != u) continue;
    const NodeId v = g.edge_dst()[e];
    if (on[static_cast<std::size_t>(v)] != 0) continue;
    on[static_cast<std::size_t>(v)] = 1;
    walk.push_back(static_cast<EdgeId>(e));
    walks(g, mem, v, is_end, walk, on, out);
    walk.pop_back();
    on[static_cast<std::size_t>(v)] = 0;
  }
}

// Every simple walk from any start to an end over the members (walks above, from each start).
inline std::vector<std::vector<EdgeId>> all_walks(const Graph& g, const std::vector<char>& mem,
                                                  const std::vector<NodeId>& starts, const std::vector<NodeId>& ends) {
  std::vector<char> is_end(static_cast<std::size_t>(g.num_nodes()), 0);
  for (NodeId t : ends) is_end[static_cast<std::size_t>(t)] = 1;
  std::vector<std::vector<EdgeId>> all;
  for (NodeId s : starts) {
    std::vector<EdgeId> walk;
    std::vector<char> on(static_cast<std::size_t>(g.num_nodes()), 0);
    on[static_cast<std::size_t>(s)] = 1;
    walks(g, mem, s, is_end, walk, on, all);
  }
  return all;
}

inline Cost walk_cost(const Graph& g, const std::vector<EdgeId>& walk) {
  Cost x = 0;
  for (EdgeId e : walk) x += g.cost()[static_cast<std::size_t>(e)];
  return x;
}

// A walk's bottleneck: its lowest priority value (cmp order).
inline Prio bottleneck(const Query& q, const std::vector<EdgeId>& walk) {
  Prio b {0, 1, true};
  for (EdgeId e : walk) {
    if (cmp(prio_of(q, e), b) < 0) b = prio_of(q, e);
  }
  return b;
}

// The canonical key (cost, -bottleneck under prefer, hops, edge-id sequence) as a strict order.
inline bool key_less(const Query& q, const std::vector<EdgeId>& a, const std::vector<EdgeId>& b) {
  const Cost ca = walk_cost(*q.g, a);
  const Cost cb = walk_cost(*q.g, b);
  if (ca != cb) return ca < cb;
  if (q.prefer) {
    const int c = cmp(bottleneck(q, a), bottleneck(q, b));
    if (c != 0) return c > 0;
  }
  if (a.size() != b.size()) return a.size() < b.size();
  return a < b;
}

// The canonical single path from the start set to the end set, or nullopt.
inline std::optional<std::vector<EdgeId>> canonical(const Query& q, const std::vector<char>& mem,
                                                    const std::vector<NodeId>& starts, const std::vector<NodeId>& ends) {
  const auto all = all_walks(*q.g, mem, starts, ends);
  if (all.empty()) return std::nullopt;
  return *std::min_element(all.begin(), all.end(), [&](const auto& a, const auto& b) { return key_less(q, a, b); });
}

} // namespace netgraph::core::test::ref
