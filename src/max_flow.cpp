/*
  The one integer minimum-cost flow core (APPROACH 2.2).

  One objective on every view: minimum cost for the maximum value attainable under the value
  cap. Successive shortest paths with potentials and blocking flows. Each tier runs Dijkstra on
  reduced costs cost(e) + pi(u) - pi(v) to completion over the residual network reachable from
  the virtual source s* (through the virtual sink t* and onward), adds the distances to the
  potentials of the vertices it settled, then drains the zero-reduced-cost admissible subgraph
  by Dinic phases (BFS levels, DFS with current-arc pointers) until no admissible augmenting
  path remains or the cap is reached; then the next tier. Positive original costs do not
  exclude zero-reduced-cost residual cycles, hence the repeated phases. The solve ends when a
  Dijkstra does not reach t* or the cap is reached; there is no completion phase.

  Arc identity: every eligible edge is its own arc pair (forward with the view capacity,
  reverse with the flow). At a node the forward arcs of its out-edges and the reverse arcs of
  its in-edges are visited in one sequence by ascending EdgeId (a self-loop's forward arc before
  its reverse arc); at a real endpoint its virtual arc comes first; s* visits its arcs by
  ascending source NodeId and t* by ascending destination NodeId. Virtual arcs cost 0 and carry
  2^63 - 2 units, above every attainable value. Every solve starts from zero flow, so the zero
  potentials are valid (all costs are positive).

  Dijkstra's queue is a radix heap: keys pop in increasing order and entries with equal keys pop
  in a fixed order (last pushed first within the lowest bucket), so the settle order, and with it
  every counter, is a function of the input alone. Equal-key order affects only the counters,
  never a distance, a potential or an augmenting path. The solve runs to the cap or to
  disconnection; there is no clamp of the cap (APPROACH 2.2, revision 25).

  Arithmetic stays in int64: every true residual distance is the cost of a simple residual path,
  below the cost sum bound 2^62; potentials are such distances; a path's bottleneck contains a
  real arc (sources and destinations are disjoint), so every amount is below 2^62.
*/
#include "netgraph/core/max_flow.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <limits>
#include <stdexcept>
#include <string>

namespace netgraph::core {

namespace {

constexpr Cost kInfDist = std::numeric_limits<Cost>::max();
constexpr Amount kUncapped = std::numeric_limits<Amount>::max();

// Sorted, duplicate-free endpoint set; throws on an empty set or an id out of range.
std::vector<NodeId> endpoint_set(std::span<const NodeId> nodes, std::int32_t n, const char* name) {
  if (nodes.empty()) {
    throw std::invalid_argument(std::string(name) + ": the endpoint set is empty");
  }
  std::vector<NodeId> out(nodes.begin(), nodes.end());
  for (NodeId v : out) {
    if (v < 0 || v >= n) {
      throw std::invalid_argument(std::string(name) + ": node " + std::to_string(v) +
                                  " out of node range [0, " + std::to_string(n) + ")");
    }
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// Monotone priority queue for Dijkstra on non-negative integer keys (a radix heap): a key goes
// to the bucket of the highest bit in which it differs from the last popped key, and a pop
// redistributes the first non-empty bucket below its minimum. Entries are (key, vertex); equal
// keys pop in a fixed order, so every run is a function of its input.
class RadixHeap {
public:
  void clear() {
    for (auto& b : buckets_) b.clear();
    last_ = 0;
    size_ = 0;
  }
  [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
  void push(Cost key, std::int32_t v) {
    buckets_[bucket(key)].emplace_back(key, v);
    ++size_;
  }
  std::pair<Cost, std::int32_t> pop() {
    if (buckets_[0].empty()) {
      std::size_t i = 1;
      while (buckets_[i].empty()) ++i;
      auto& from = buckets_[i];
      Cost lowest = from.front().first;
      for (const auto& entry : from) lowest = std::min(lowest, entry.first);
      last_ = lowest;
      for (const auto& entry : from) buckets_[bucket(entry.first)].push_back(entry);
      from.clear();
    }
    const auto top = buckets_[0].back();
    buckets_[0].pop_back();
    --size_;
    return top;
  }

private:
  [[nodiscard]] std::size_t bucket(Cost key) const noexcept {
    return static_cast<std::size_t>(
        std::bit_width(static_cast<std::uint64_t>(key) ^ static_cast<std::uint64_t>(last_)));
  }
  std::array<std::vector<std::pair<Cost, std::int32_t>>, 65> buckets_ {};
  Cost last_ {0};
  std::size_t size_ {0};
};

// The residual network of one solve: real vertices 0..N-1, s* = N, t* = N + 1. Arcs are stored
// as slots in adjacency order (the slots of u are [off_[u], off_[u + 1])), each with its head,
// its cost, its residual and the slot of its partner arc.
class Solver {
public:
  Solver(const StrictMultiDiGraph& g, std::span<const Amount> view, const std::vector<NodeId>& srcs,
         const std::vector<NodeId>& dsts, const Eligibility* elig)
      : g_(g), n_(g.num_nodes() + 2), s_(g.num_nodes()), t_(g.num_nodes() + 1) {
    const auto num_edges = static_cast<std::size_t>(g.num_edges());
    const auto edge_src = g.edge_src();
    const auto edge_dst = g.edge_dst();
    const auto cost = g.cost();
    active_.reserve(num_edges);
    for (std::size_t e = 0; e < num_edges; ++e) {
      if (elig == nullptr || elig->edge_ok(static_cast<EdgeId>(e))) active_.push_back(static_cast<EdgeId>(e));
    }
    for (NodeId s : srcs) {
      if (elig == nullptr || elig->node_ok(s)) srcs_.push_back(s);
    }
    for (NodeId d : dsts) {
      if (elig == nullptr || elig->node_ok(d)) dsts_.push_back(d);
    }
    const auto n = static_cast<std::size_t>(n_);
    off_.assign(n + 1, 0);
    auto deg = [&](std::int32_t v) { ++off_[static_cast<std::size_t>(v) + 1]; };
    for (EdgeId e : active_) {
      deg(edge_src[static_cast<std::size_t>(e)]);
      deg(edge_dst[static_cast<std::size_t>(e)]);
    }
    for (NodeId s : srcs_) { deg(s); deg(s_); }
    for (NodeId d : dsts_) { deg(d); deg(t_); }
    for (std::size_t v = 0; v < n; ++v) off_[v + 1] += off_[v];
    const auto num_slots = static_cast<std::size_t>(off_[n]);
    head_.resize(num_slots);
    rev_.resize(num_slots);
    cost_.resize(num_slots);
    res_.resize(num_slots);
    std::vector<std::int32_t> pos(off_.begin(), off_.end() - 1);
    auto take = [&](std::int32_t v) { return pos[static_cast<std::size_t>(v)]++; };
    auto link = [&](std::int32_t a, std::int32_t to_a, Amount cap, Cost w, std::int32_t b, std::int32_t to_b) {
      const auto ai = static_cast<std::size_t>(a);
      const auto bi = static_cast<std::size_t>(b);
      head_[ai] = to_a;
      rev_[ai] = b;
      cost_[ai] = w;
      res_[ai] = cap;
      head_[bi] = to_b;
      rev_[bi] = a;
      cost_[bi] = -w;
      res_[bi] = 0;
    };
    // Virtual arcs first at the real endpoints, then the real arcs by ascending EdgeId (forward
    // slot at the tail, reverse slot at the head), then the slots of s* and t*.
    std::vector<std::int32_t> src_rev(srcs_.size()), dst_fwd(dsts_.size());
    for (std::size_t i = 0; i < srcs_.size(); ++i) src_rev[i] = take(srcs_[i]);
    for (std::size_t i = 0; i < dsts_.size(); ++i) dst_fwd[i] = take(dsts_[i]);
    edge_slot_.resize(active_.size());
    for (std::size_t k = 0; k < active_.size(); ++k) {
      const auto ei = static_cast<std::size_t>(active_[k]);
      const NodeId u = edge_src[ei];
      const NodeId v = edge_dst[ei];
      const std::int32_t fwd = take(u);
      const std::int32_t bwd = take(v);
      link(fwd, v, view[ei], cost[ei], bwd, u);
      edge_slot_[k] = fwd;
    }
    for (std::size_t i = 0; i < srcs_.size(); ++i) link(take(s_), srcs_[i], kMaxFlowCap, 0, src_rev[i], s_);
    for (std::size_t i = 0; i < dsts_.size(); ++i) link(dst_fwd[i], t_, kMaxFlowCap, 0, take(t_), dsts_[i]);
    pi_.assign(n, 0);
    ever_.assign(n, 0);
    dist_.resize(n);
    level_.resize(n);
    reach_.assign(n, 0);
    cur_.resize(n);
  }

  FlowCoreResult run(std::span<const Amount> view, Amount limit) {
    FlowCoreResult out;
    auto& c = out.counters;
    Amount total = 0;
    while (total < limit) {
      if (!dijkstra(c)) break;
      const Cost tier_cost = pi_[static_cast<std::size_t>(t_)] - pi_[static_cast<std::size_t>(s_)];
      const Amount before = total;
      while (total < limit) {
        if (!build_levels(c)) break;
        std::copy(off_.begin(), off_.end() - 1, cur_.begin());
        while (total < limit) {
          if (!find_path(c)) break;
          Amount amount = limit - total;
          for (std::int32_t a : path_) amount = std::min(amount, res_[static_cast<std::size_t>(a)]);
          for (std::int32_t a : path_) {
            res_[static_cast<std::size_t>(a)] -= amount;
            res_[static_cast<std::size_t>(rev_[static_cast<std::size_t>(a)])] += amount;
          }
          total += amount;
          ++c.augmentations;
        }
      }
      // A successful Dijkstra reaches t* over arcs with positive residual and zero reduced cost
      // along its shortest paths, so the drain moves at least one unit.
      if (total == before) {
        throw std::logic_error("flow core: a successful Dijkstra drained nothing");
      }
      out.cost_distribution.emplace_back(tier_cost, total - before);
    }
    out.placed = total;

    const auto num_edges = static_cast<std::size_t>(g_.num_edges());
    out.edge_flows.assign(num_edges, 0);
    for (std::size_t k = 0; k < active_.size(); ++k) {
      const auto fwd = static_cast<std::size_t>(edge_slot_[k]);
      out.edge_flows[static_cast<std::size_t>(active_[k])] = res_[static_cast<std::size_t>(rev_[fwd])];
    }
    // Full residual closure seeded at the real sources (reverse and virtual arcs included).
    std::vector<std::uint8_t> seen(static_cast<std::size_t>(n_), 0);
    std::vector<std::int32_t> queue;
    queue.reserve(static_cast<std::size_t>(n_));
    for (NodeId s : srcs_) {
      seen[static_cast<std::size_t>(s)] = 1;
      queue.push_back(s);
    }
    for (std::size_t qi = 0; qi < queue.size(); ++qi) {
      const auto u = static_cast<std::size_t>(queue[qi]);
      for (auto p = off_[u]; p < off_[u + 1]; ++p) {
        const auto pi = static_cast<std::size_t>(p);
        const auto v = static_cast<std::size_t>(head_[pi]);
        if (res_[pi] > 0 && seen[v] == 0) {
          seen[v] = 1;
          queue.push_back(static_cast<std::int32_t>(v));
        }
      }
    }
    if (seen[static_cast<std::size_t>(t_)] == 0) {
      const auto edge_src = g_.edge_src();
      const auto edge_dst = g_.edge_dst();
      std::vector<EdgeId> cut;
      for (EdgeId e : active_) {
        const auto ei = static_cast<std::size_t>(e);
        if (seen[static_cast<std::size_t>(edge_src[ei])] != 0 &&
            seen[static_cast<std::size_t>(edge_dst[ei])] == 0 && view[ei] > 0) {
          cut.push_back(e);
        }
      }
      out.min_cut = std::move(cut);
    }
    out.potentials.resize(static_cast<std::size_t>(g_.num_nodes()));
    for (std::size_t v = 0; v < out.potentials.size(); ++v) {
      out.potentials[v] = ever_[v] != 0 ? pi_[v] : kNeverReached;
    }
    return out;
  }

private:
  // Dijkstra on reduced costs to completion; potentials of the settled vertices updated.
  // Returns whether t* was reached.
  bool dijkstra(FlowCoreCounters& c) {
    ++c.dijkstra_runs;
    std::fill(dist_.begin(), dist_.end(), kInfDist);
    settled_.clear();
    heap_.clear();
    dist_[static_cast<std::size_t>(s_)] = 0;
    heap_.push(0, s_);
    while (!heap_.empty()) {
      const auto [du, ui] = heap_.pop();
      ++c.dijkstra_pops;
      const auto u = static_cast<std::size_t>(ui);
      if (du != dist_[u]) continue;
      settled_.push_back(ui);
      ++c.settled;
      const Cost base = du + pi_[u];
      const auto end = off_[u + 1];
      c.arc_scans += end - off_[u];
      for (auto p = off_[u]; p < end; ++p) {
        const auto pi = static_cast<std::size_t>(p);
        if (res_[pi] <= 0) continue;
        const auto v = static_cast<std::size_t>(head_[pi]);
        const Cost nd = base + cost_[pi] - pi_[v];
        if (nd < dist_[v]) {
          dist_[v] = nd;
          heap_.push(nd, static_cast<std::int32_t>(v));
        }
      }
    }
    for (std::int32_t v : settled_) {
      const auto vi = static_cast<std::size_t>(v);
      pi_[vi] += dist_[vi];
      ever_[vi] = 1;
    }
    return dist_[static_cast<std::size_t>(t_)] != kInfDist;
  }

  // BFS levels over admissible arcs (positive residual, zero reduced cost). Counts one phase.
  // Only vertices that can reach t* over admissible arcs are levelled: every vertex on a
  // shortest admissible path from s* to such a vertex reaches t* too, so their levels equal the
  // levels of a full BFS, and the other vertices are dead ends that the DFS would only enter and
  // leave. The forward search also stops once t*'s level is complete. The augmenting paths are
  // therefore the ones a full BFS gives.
  bool build_levels(FlowCoreCounters& c) {
    const auto t = static_cast<std::size_t>(t_);
    // Reverse search from t*: the arc into v from head_[p] is the partner of v's slot p.
    ++mark_epoch_;
    queue_.clear();
    reach_[t] = mark_epoch_;
    queue_.push_back(t_);
    for (std::size_t qi = 0; qi < queue_.size(); ++qi) {
      const auto v = static_cast<std::size_t>(queue_[qi]);
      const Cost pv = pi_[v];
      const auto end = off_[v + 1];
      c.arc_scans += end - off_[v];
      for (auto p = off_[v]; p < end; ++p) {
        const auto pi = static_cast<std::size_t>(p);
        const auto u = static_cast<std::size_t>(head_[pi]);
        if (reach_[u] == mark_epoch_) continue;
        const auto in = static_cast<std::size_t>(rev_[pi]);
        if (res_[in] > 0 && cost_[in] + pi_[u] == pv) {
          reach_[u] = mark_epoch_;
          queue_.push_back(static_cast<std::int32_t>(u));
        }
      }
    }
    std::fill(level_.begin(), level_.end(), -1);
    ++c.dinic_phases;
    const auto s = static_cast<std::size_t>(s_);
    if (reach_[s] != mark_epoch_) return false;
    queue_.clear();
    level_[s] = 0;
    queue_.push_back(s_);
    for (std::size_t qi = 0; qi < queue_.size(); ++qi) {
      const auto u = static_cast<std::size_t>(queue_[qi]);
      if (level_[t] >= 0 && level_[u] >= level_[t]) break;
      const Cost pu = pi_[u];
      const auto next = level_[u] + 1;
      const auto end = off_[u + 1];
      c.arc_scans += end - off_[u];
      for (auto p = off_[u]; p < end; ++p) {
        const auto pi = static_cast<std::size_t>(p);
        const auto v = static_cast<std::size_t>(head_[pi]);
        if (level_[v] < 0 && reach_[v] == mark_epoch_ && res_[pi] > 0 && cost_[pi] + pu == pi_[v]) {
          level_[v] = next;
          queue_.push_back(static_cast<std::int32_t>(v));
        }
      }
    }
    return level_[t] >= 0;
  }

  // One augmenting path s* -> t* in the level graph into path_ (current-arc pointers kept: an
  // arc that is pushed through stays current; a dead end advances its parent's pointer).
  bool find_path(FlowCoreCounters& c) {
    stack_.clear();
    path_.clear();
    stack_.push_back(s_);
    while (!stack_.empty()) {
      const auto u = static_cast<std::size_t>(stack_.back());
      if (static_cast<std::int32_t>(u) == t_) return true;
      bool advanced = false;
      const auto end = off_[u + 1];
      const auto next = level_[u] + 1;
      const Cost pu = pi_[u];
      while (cur_[u] < end) {
        ++c.arc_scans;
        const auto pi = static_cast<std::size_t>(cur_[u]);
        const auto v = static_cast<std::size_t>(head_[pi]);
        if (res_[pi] > 0 && level_[v] == next && cost_[pi] + pu == pi_[v]) {
          stack_.push_back(static_cast<std::int32_t>(v));
          path_.push_back(static_cast<std::int32_t>(pi));
          advanced = true;
          break;
        }
        ++cur_[u];
      }
      if (!advanced) {
        stack_.pop_back();
        if (!path_.empty()) {
          path_.pop_back();
          ++cur_[static_cast<std::size_t>(stack_.back())];
        }
      }
    }
    return false;
  }

  const StrictMultiDiGraph& g_;
  std::int32_t n_;
  std::int32_t s_;
  std::int32_t t_;
  std::vector<EdgeId> active_ {};
  std::vector<NodeId> srcs_ {};
  std::vector<NodeId> dsts_ {};
  std::vector<std::int32_t> off_ {};
  std::vector<std::int32_t> head_ {};
  std::vector<std::int32_t> rev_ {};
  std::vector<Cost> cost_ {};
  std::vector<Amount> res_ {};
  std::vector<std::int32_t> edge_slot_ {};
  std::vector<Cost> pi_ {};
  std::vector<std::uint8_t> ever_ {};
  std::vector<Cost> dist_ {};
  std::vector<std::int32_t> settled_ {};
  RadixHeap heap_ {};
  std::vector<std::int32_t> level_ {};
  std::vector<std::uint32_t> reach_ {};
  std::uint32_t mark_epoch_ {0};
  std::vector<std::int32_t> queue_ {};
  std::vector<std::int32_t> cur_ {};
  std::vector<std::int32_t> stack_ {};
  std::vector<std::int32_t> path_ {};
};

} // namespace

FlowCoreResult place_global_units(const StrictMultiDiGraph& g,
                                  std::span<const Amount> view,
                                  std::span<const NodeId> sources,
                                  std::span<const NodeId> destinations,
                                  const Eligibility* elig,
                                  std::optional<Amount> cap) {
  const auto num_edges = static_cast<std::size_t>(g.num_edges());
  if (view.size() != num_edges) {
    throw std::invalid_argument("view has length " + std::to_string(view.size()) + ", expected " +
                                std::to_string(num_edges) + " (one entry per edge)");
  }
  const auto capacity = g.capacity_units();
  for (std::size_t e = 0; e < num_edges; ++e) {
    if (view[e] < 0 || view[e] > capacity[e]) {
      throw std::invalid_argument("view[" + std::to_string(e) + "] = " + std::to_string(view[e]) +
                                  " is outside [0, capacity_units[" + std::to_string(e) + "] = " +
                                  std::to_string(capacity[e]) + "]");
    }
  }
  const auto srcs = endpoint_set(sources, g.num_nodes(), "sources");
  const auto dsts = endpoint_set(destinations, g.num_nodes(), "destinations");
  for (std::size_t i = 0, k = 0; i < srcs.size() && k < dsts.size();) {
    if (srcs[i] == dsts[k]) {
      throw std::invalid_argument("node " + std::to_string(srcs[i]) +
                                  " is in both endpoint sets: sources and destinations must not overlap");
    }
    if (srcs[i] < dsts[k]) ++i; else ++k;
  }
  if (cap.has_value() && (*cap < 0 || *cap > kMaxFlowCap)) {
    throw std::invalid_argument("cap " + std::to_string(*cap) + " is outside [0, 2^63 - 2]");
  }
  if (elig != nullptr && elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  if (cap.has_value() && *cap == 0) {
    // The input cap 0 runs nothing: the prescribed zero payload.
    FlowCoreResult out;
    out.edge_flows.assign(num_edges, 0);
    out.potentials.assign(static_cast<std::size_t>(g.num_nodes()), kNeverReached);
    return out;
  }
  Solver solver(g, view, srcs, dsts, elig);
  return solver.run(view, cap.value_or(kUncapped));
}

} // namespace netgraph::core
