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

  The deletion re-solve (the exact sensitivity model, APPROACH 11) continues from a solve's
  final state instead. A solve asked to capture its state also records, as a shadow that
  changes none of its searches, the working potentials: the running sum of the tier maxima D_k
  and, per settled vertex, that sum at its last settlement; the working potential is the
  published one plus the difference (a never-settled vertex: the whole sum, from numeric zero).
  They equal the update of every vertex after every Dijkstra, the reached ones by their distance
  and the others by D_k, so the reduced cost is non-negative on every positive residual arc,
  real or virtual. The re-solve deletes a real arc e = (u, v) with flow f_e > 0 from a copy of
  the state and runs three fixed-origin repairs, each a sequence of tiers as above from its
  origin to its target, with that all-vertex update before the Dinic phases: u -> v for f_e,
  moving m; when r = f_e - m > 0, u -> s* for r, whose first tier reuses repair 1's failed
  search (nothing changed since), then t* -> v for r. The result is the minimum-cost flow of the
  maximum value placed - r on the view without e, certified by the working potentials and the
  real sources' residual closure; at most m + 2 r Dijkstra runs.

  Dijkstra's queue is a radix heap: keys pop in increasing order and entries with equal keys pop
  in a fixed order (last pushed first within the lowest bucket), so the settle order, and with it
  every counter, is a function of the input alone. Equal-key order affects only the counters,
  never a distance, a potential or an augmenting path. The solve runs to the cap or to
  disconnection; there is no clamp of the cap (APPROACH 2.2, revision 25).

  Arithmetic. A solve from zero stays in int64: every true residual distance is the cost of a
  simple residual path, below the cost sum bound 2^62; potentials are such distances; a path's
  bottleneck contains a real arc (sources and destinations are disjoint), so every amount is
  below 2^62. The deletion mode's working potentials, distance keys, tentative distances and
  sums of tier maxima are signed 128-bit (I128; reachability kept apart from the value by a
  stamp): they stay within (n + 1) C at capture and (22 n + 15) C across the repairs for n
  vertices and the cost sum C < 2^62, and its heap takes 128-bit keys. Amounts stay int64.
*/
#include "netgraph/core/detail/flow_core.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/node_map.hpp"

namespace netgraph::core::detail {

namespace {

constexpr Cost kInfDist = std::numeric_limits<Cost>::max();
constexpr Amount kUncapped = std::numeric_limits<Amount>::max();

// The bucket of a key: the number of the highest bit in which it differs from the last popped
// key (keys are non-negative).
std::size_t bucket_of(Cost key, Cost last) noexcept {
  return static_cast<std::size_t>(std::bit_width(static_cast<std::uint64_t>(key) ^ static_cast<std::uint64_t>(last)));
}
std::size_t bucket_of(I128 key, I128 last) noexcept {
  return static_cast<std::size_t>(bits(U128 {static_cast<std::uint64_t>(key.hi ^ last.hi), key.lo ^ last.lo}));
}

// Monotone priority queue for Dijkstra on non-negative integer keys (a radix heap): a key goes
// to the bucket of the highest bit in which it differs from the last popped key, and a pop
// redistributes the first non-empty bucket below its minimum. Entries are (key, vertex); equal
// keys pop in a fixed order, so every run is a function of its input. Keys are int64 (a solve
// from zero) or I128 (the deletion mode).
template <class Key>
class RadixHeap {
public:
  void clear() {
    for (auto& b : buckets_) b.clear();
    last_ = 0;
    size_ = 0;
  }
  [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
  void push(Key key, std::int32_t v) {
    buckets_[bucket_of(key, last_)].emplace_back(key, v);
    ++size_;
  }
  std::pair<Key, std::int32_t> pop() {
    if (buckets_[0].empty()) {
      std::size_t i = 1;
      while (buckets_[i].empty()) ++i;
      auto& from = buckets_[i];
      Key lowest = from.front().first;
      for (const auto& entry : from) lowest = std::min(lowest, entry.first);
      last_ = lowest;
      for (const auto& entry : from) buckets_[bucket_of(entry.first, last_)].push_back(entry);
      from.clear();
    }
    const auto top = buckets_[0].back();
    buckets_[0].pop_back();
    --size_;
    return top;
  }

private:
  std::array<std::vector<std::pair<Key, std::int32_t>>, 8 * sizeof(Key) + 1> buckets_ {};
  Key last_ {0};
  std::size_t size_ {0};
};

template <class P>
inline constexpr bool kWide = std::is_same_v<P, I128>;

} // namespace

// The residual network of one solve (Network: real vertices 0..N-1, s* = N, t* = N + 1) with
// its flow, potentials and scratch; also the working storage of a deletion re-solve.
class FlowSolver {
public:
  // One reusable instance per thread (`local()`) for the solves from zero: build() resizes every
  // buffer for the input, keeping capacity across calls, so repeated solves do not allocate
  // once warm. Not reentrant (a solve never starts another). A deletion resolver owns its own.
  static FlowSolver& local() {
    thread_local FlowSolver solver;
    return solver;
  }

  // The input part of the network (num_nodes, the edges, the endpoint sets) is the caller's to
  // fill; build() adds the arcs.
  Network& network() noexcept { return net_; }

  // The arcs of the input, zero flow and zero potentials. `shadow`: record the working
  // potentials for a capture.
  void build(bool shadow) {
    auto& in = net_;
    size(in.num_nodes);
    const auto n = static_cast<std::size_t>(n_);
    in.off.assign(n + 1, 0);
    auto deg = [&](std::int32_t v) { ++in.off[static_cast<std::size_t>(v) + 1]; };
    for (std::size_t k = 0; k < in.edges.size(); ++k) {
      deg(in.tail[k]);
      deg(in.head[k]);
    }
    for (std::int32_t s : in.sources) { deg(s); deg(s_); }
    for (std::int32_t d : in.destinations) { deg(d); deg(t_); }
    for (std::size_t v = 0; v < n; ++v) in.off[v + 1] += in.off[v];
    const auto num_slots = static_cast<std::size_t>(in.off[n]);
    in.arc_head.resize(num_slots);
    in.rev.resize(num_slots);
    in.arc_cost.resize(num_slots);
    res_.resize(num_slots);
    pos_.assign(in.off.begin(), in.off.end() - 1);
    auto take = [&](std::int32_t v) { return pos_[static_cast<std::size_t>(v)]++; };
    auto link = [&](std::int32_t a, std::int32_t to_a, Amount cap, Cost w, std::int32_t b, std::int32_t to_b) {
      const auto ai = static_cast<std::size_t>(a);
      const auto bi = static_cast<std::size_t>(b);
      in.arc_head[ai] = to_a;
      in.rev[ai] = b;
      in.arc_cost[ai] = w;
      res_[ai] = cap;
      in.arc_head[bi] = to_b;
      in.rev[bi] = a;
      in.arc_cost[bi] = -w;
      res_[bi] = 0;
    };
    // Virtual arcs first at the real endpoints, then the real arcs by ascending EdgeId (forward
    // slot at the tail, reverse slot at the head), then the slots of s* and t*.
    in.src_rev.resize(in.sources.size());
    in.dst_fwd.resize(in.destinations.size());
    for (std::size_t i = 0; i < in.sources.size(); ++i) in.src_rev[i] = take(in.sources[i]);
    for (std::size_t i = 0; i < in.destinations.size(); ++i) in.dst_fwd[i] = take(in.destinations[i]);
    in.edge_slot.resize(in.edges.size());
    for (std::size_t k = 0; k < in.edges.size(); ++k) {
      const std::int32_t u = in.tail[k];
      const std::int32_t v = in.head[k];
      const std::int32_t fwd = take(u);
      const std::int32_t bwd = take(v);
      link(fwd, v, in.cap[k], in.cost[k], bwd, u);
      in.edge_slot[k] = fwd;
    }
    for (std::size_t i = 0; i < in.sources.size(); ++i) {
      link(take(s_), in.sources[i], kMaxFlowCap, 0, in.src_rev[i], s_);
    }
    for (std::size_t i = 0; i < in.destinations.size(); ++i) {
      link(in.dst_fwd[i], t_, kMaxFlowCap, 0, take(t_), in.destinations[i]);
    }
    pi_.assign(n, 0);
    ever_.assign(n, 0);
    dist_.resize(n);
    shadow_ = shadow;
    if (shadow) {
      tier_sum_ = 0;
      checkpoint_.assign(n, I128 {});
    }
  }

  // A captured state's network, as the working storage of its deletion re-solves.
  void load(const Network& net) {
    net_ = net;
    size(net_.num_nodes);
    res_.resize(net_.arc_head.size());
    wpi_.resize(static_cast<std::size_t>(n_));
    wdist_.resize(static_cast<std::size_t>(n_));
    wseen_.assign(static_cast<std::size_t>(n_), 0);
    wepoch_ = 0;
  }

  // The solve to the cap or to disconnection. Output in local coordinates: `flows` aligned with
  // the input edges, `cut` as input edge indices (positive capacity, from the full residual
  // closure of the real sources to its complement, ascending EdgeId) when the destinations are
  // outside the closure, potentials per local vertex (kNeverReached when never settled).
  struct Output {
    Amount placed {0};
    std::vector<std::pair<Cost, Amount>> cost_distribution {};
    std::vector<Amount> flows {};
    std::optional<std::vector<std::int32_t>> cut {};
    std::vector<Cost> potentials {};
    FlowCoreCounters counters {};
  };

  // `closure` adds the cut; `potentials` the per-vertex potentials (the edge-subset entry
  // needs neither).
  Output run(Amount limit, bool closure, bool potentials = true) {
    Output out;
    auto& c = out.counters;
    Amount total = 0;
    origin_ = s_;
    target_ = t_;
    while (total < limit) {
      const bool reached = dijkstra<Cost>(c);
      settle();
      if (!reached) break;
      const Cost tier_cost = pi_[static_cast<std::size_t>(t_)] - pi_[static_cast<std::size_t>(s_)];
      const Amount before = total;
      total = drain<Cost>(c, limit, total, [](Amount) {});
      // A successful Dijkstra reaches t* over arcs with positive residual and zero reduced cost
      // along its shortest paths, so the drain moves at least one unit.
      if (total == before) {
        throw std::logic_error("flow core: a successful Dijkstra drained nothing");
      }
      out.cost_distribution.emplace_back(tier_cost, total - before);
    }
    out.placed = total;

    out.flows.resize(net_.edges.size());
    for (std::size_t k = 0; k < net_.edges.size(); ++k) out.flows[k] = flow_of(k);
    if (closure && close_sources()) out.cut = cut_edges(-1);
    if (potentials) {
      out.potentials.resize(static_cast<std::size_t>(net_.num_nodes));
      for (std::size_t v = 0; v < out.potentials.size(); ++v) {
        out.potentials[v] = ever_[v] != 0 ? pi_[v] : kNeverReached;
      }
    }
    return out;
  }

  // The final state of the last solve (its value and counters given) into `req`.
  void capture(DeletionCapture& req, Amount placed, const FlowCoreCounters& counters) {
    DeletionState state;
    state.graph_ = req.graph;
    state.net_ = std::make_shared<const Network>(net_);
    state.res_ = res_;
    state.pi_.resize(static_cast<std::size_t>(n_));
    for (std::size_t v = 0; v < state.pi_.size(); ++v) state.pi_[v] = I128 {pi_[v]} + (tier_sum_ - checkpoint_[v]);
    state.value_ = placed;
    state.maximum_ = close_sources();
    state.counters_ = counters;
    req.state = std::move(state);
  }

#if NETGRAPH_CORE_PRIVATE_HOOKS
  void observe(const std::function<void(std::span<const Cost>)>* on_search) noexcept { on_search_ = on_search; }
  // The tests' probe of the re-solves until reset (nullptr: none).
  void probe(const DeletionProbe* p) noexcept { probe_ = p; }
#endif

  // The deletion re-solve of input edge k (a real arc with positive flow) from `state`, whose
  // network this solver holds (load()).
  DeletionResult resolve(const DeletionState& state, std::size_t k) {
    std::copy(state.residual().begin(), state.residual().end(), res_.begin());
    std::copy(state.potentials().begin(), state.potentials().end(), wpi_.begin());
    DeletionResult out;
    auto& c = out.counters;
    const auto fwd = static_cast<std::size_t>(net_.edge_slot[k]);
    const auto bwd = static_cast<std::size_t>(net_.rev[fwd]);
    const Amount f = res_[bwd];
    const std::int32_t u = net_.tail[k];
    const std::int32_t v = net_.head[k];
    // The deletion: a pseudo-flow with an excess of f at u and a deficit of f at v, its residual
    // the state's without e's two arcs.
    res_[fwd] = 0;
    res_[bwd] = 0;
    emit(0, false, 0, c);
    const Amount m = repair(1, u, v, f, false, c);
    const Amount r = f - m;
    if (r > 0) {
      bool reuse = true;
#if NETGRAPH_CORE_PRIVATE_HOOKS
      if (probe_ != nullptr && probe_->rerun_failed_search) reuse = false;
#endif
      if (repair(2, u, s_, r, reuse, c) != r) {
        throw std::logic_error("deletion re-solve: the remainder has no residual path back to the sources");
      }
      if (repair(3, t_, v, r, false, c) != r) {
        throw std::logic_error("deletion re-solve: the deficit has no residual path from the destinations");
      }
    }
    emit(4, false, 0, c);
    if (!close_sources()) {
      throw std::logic_error("deletion re-solve: the destinations stay connected to the sources");
    }
    out.value = state.value() - r;
    out.potentials = wpi_;
    for (std::int32_t j : cut_edges(static_cast<std::int64_t>(k))) out.cut.push_back(net_.edges[static_cast<std::size_t>(j)]);
    out.edge_flows.assign(static_cast<std::size_t>(state.graph().num_edges()), 0);
    for (std::size_t j = 0; j < net_.edges.size(); ++j) {
      out.edge_flows[static_cast<std::size_t>(net_.edges[j])] = j == k ? 0 : flow_of(j);
    }
    return out;
  }

private:
  void size(std::int32_t num_nodes) {
    n_ = num_nodes + 2;
    s_ = num_nodes;
    t_ = num_nodes + 1;
    const auto n = static_cast<std::size_t>(n_);
    level_.resize(n);
    reach_.assign(n, 0);
    mark_epoch_ = 0;
    cur_.resize(n);
  }

  // The flow of input edge k: the residual of its reverse slot.
  [[nodiscard]] Amount flow_of(std::size_t k) const {
    return res_[static_cast<std::size_t>(net_.rev[static_cast<std::size_t>(net_.edge_slot[k])])];
  }

  template <class P>
  std::vector<P>& potentials() {
    if constexpr (kWide<P>) {
      return wpi_;
    } else {
      return pi_;
    }
  }
  template <class P>
  std::vector<P>& distances() {
    if constexpr (kWide<P>) {
      return wdist_;
    } else {
      return dist_;
    }
  }
  template <class P>
  RadixHeap<P>& heap() {
    if constexpr (kWide<P>) {
      return wheap_;
    } else {
      return heap_;
    }
  }
  // Whether the last search of the mode reached x (int64: a finite distance; I128: the stamp).
  template <class P>
  [[nodiscard]] bool reached(std::int32_t x) const {
    const auto xi = static_cast<std::size_t>(x);
    if constexpr (kWide<P>) {
      return wseen_[xi] == wepoch_;
    } else {
      return dist_[xi] != kInfDist;
    }
  }

  // Dijkstra on reduced costs from origin_ to completion over the residual network. Returns
  // whether target_ was reached; settled_ lists the vertices it settled.
  template <class P>
  bool dijkstra(FlowCoreCounters& c) {
    auto& pi = potentials<P>();
    auto& dist = distances<P>();
    auto& queue = heap<P>();
    ++c.dijkstra_runs;
    if constexpr (kWide<P>) {
      if (++wepoch_ == 0) {
        std::fill(wseen_.begin(), wseen_.end(), 0);
        wepoch_ = 1;
      }
    } else {
      std::fill(dist.begin(), dist.end(), kInfDist);
    }
    settled_.clear();
    queue.clear();
    const auto o = static_cast<std::size_t>(origin_);
    dist[o] = P {0};
    if constexpr (kWide<P>) wseen_[o] = wepoch_;
    queue.push(P {0}, origin_);
    while (!queue.empty()) {
      const auto [du, ui] = queue.pop();
      ++c.dijkstra_pops;
      const auto u = static_cast<std::size_t>(ui);
      if (du != dist[u]) continue;
      settled_.push_back(ui);
      ++c.settled;
      const P base = du + pi[u];
      const auto end = net_.off[u + 1];
      c.arc_scans += end - net_.off[u];
      for (auto p = net_.off[u]; p < end; ++p) {
        const auto a = static_cast<std::size_t>(p);
        if (res_[a] <= 0) continue;
        const auto v = static_cast<std::size_t>(net_.arc_head[a]);
        const P nd = base + net_.arc_cost[a] - pi[v];
        if constexpr (kWide<P>) {
          if (wseen_[v] != wepoch_ || nd < dist[v]) {
            wseen_[v] = wepoch_;
            dist[v] = nd;
            queue.push(nd, static_cast<std::int32_t>(v));
          }
        } else if (nd < dist[v]) {
          dist[v] = nd;
          queue.push(nd, static_cast<std::int32_t>(v));
        }
      }
    }
    return reached<P>(target_);
  }

  // The update of a solve from zero: the settled vertices' potentials by their distances; with
  // a capture, the running sum of the tier maxima and the settled vertices' checkpoints.
  void settle() {
    for (std::int32_t v : settled_) {
      const auto vi = static_cast<std::size_t>(v);
      pi_[vi] += dist_[vi];
      ever_[vi] = 1;
    }
    if (shadow_) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
      if (on_search_ != nullptr && *on_search_) (*on_search_)(dist_);
#endif
      Cost top = 0;
      for (std::int32_t v : settled_) top = std::max(top, dist_[static_cast<std::size_t>(v)]);
      tier_sum_ += top;
      for (std::int32_t v : settled_) checkpoint_[static_cast<std::size_t>(v)] = tier_sum_;
    }
  }

  // The update of the deletion mode: every vertex, the reached ones by their distance and the
  // others by the tier's largest finite distance.
  void update_all() {
    I128 top = 0;
    for (std::int32_t v : settled_) top = std::max(top, wdist_[static_cast<std::size_t>(v)]);
    for (std::size_t v = 0; v < wpi_.size(); ++v) wpi_[v] += wseen_[v] == wepoch_ ? wdist_[v] : top;
  }

  // The Dinic phases from origin_ to target_ over the admissible subgraph (positive residual,
  // zero reduced cost) until no admissible path remains or `limit` is reached; on_push(total)
  // after every push.
  template <class P, class OnPush>
  Amount drain(FlowCoreCounters& c, Amount limit, Amount total, OnPush on_push) {
    while (total < limit) {
      if (!build_levels<P>(c)) break;
      std::copy(net_.off.begin(), net_.off.end() - 1, cur_.begin());
      while (total < limit) {
        if (!find_path<P>(c)) break;
        Amount amount = limit - total;
        for (std::int32_t a : path_) amount = std::min(amount, res_[static_cast<std::size_t>(a)]);
        for (std::int32_t a : path_) {
          res_[static_cast<std::size_t>(a)] -= amount;
          res_[static_cast<std::size_t>(net_.rev[static_cast<std::size_t>(a)])] += amount;
        }
        total += amount;
        ++c.augmentations;
        on_push(total);
      }
    }
    return total;
  }

  // BFS levels over admissible arcs (positive residual, zero reduced cost). Counts one phase.
  // Only vertices that can reach target_ over admissible arcs are levelled: every vertex on a
  // shortest admissible path from origin_ to such a vertex reaches the target too, so their
  // levels equal the levels of a full BFS, and the other vertices are dead ends that the DFS
  // would only enter and leave. The forward search also stops once the target's level is
  // complete. The augmenting paths are therefore the ones a full BFS gives.
  template <class P>
  bool build_levels(FlowCoreCounters& c) {
    const auto& pi = potentials<P>();
    const auto t = static_cast<std::size_t>(target_);
    // Reverse search from the target: the arc into v from arc_head[p] is the partner of v's slot p.
    ++mark_epoch_;
    queue_.clear();
    reach_[t] = mark_epoch_;
    queue_.push_back(target_);
    for (std::size_t qi = 0; qi < queue_.size(); ++qi) {
      const auto v = static_cast<std::size_t>(queue_[qi]);
      const P pv = pi[v];
      const auto end = net_.off[v + 1];
      c.arc_scans += end - net_.off[v];
      for (auto p = net_.off[v]; p < end; ++p) {
        const auto a = static_cast<std::size_t>(p);
        const auto u = static_cast<std::size_t>(net_.arc_head[a]);
        if (reach_[u] == mark_epoch_) continue;
        const auto in = static_cast<std::size_t>(net_.rev[a]);
        if (res_[in] > 0 && net_.arc_cost[in] + pi[u] == pv) {
          reach_[u] = mark_epoch_;
          queue_.push_back(static_cast<std::int32_t>(u));
        }
      }
    }
    std::fill(level_.begin(), level_.end(), -1);
    ++c.dinic_phases;
    const auto s = static_cast<std::size_t>(origin_);
    if (reach_[s] != mark_epoch_) return false;
    queue_.clear();
    level_[s] = 0;
    queue_.push_back(origin_);
    for (std::size_t qi = 0; qi < queue_.size(); ++qi) {
      const auto u = static_cast<std::size_t>(queue_[qi]);
      if (level_[t] >= 0 && level_[u] >= level_[t]) break;
      const P pu = pi[u];
      const auto next = level_[u] + 1;
      const auto end = net_.off[u + 1];
      c.arc_scans += end - net_.off[u];
      for (auto p = net_.off[u]; p < end; ++p) {
        const auto a = static_cast<std::size_t>(p);
        const auto v = static_cast<std::size_t>(net_.arc_head[a]);
        if (level_[v] < 0 && reach_[v] == mark_epoch_ && res_[a] > 0 && net_.arc_cost[a] + pu == pi[v]) {
          level_[v] = next;
          queue_.push_back(static_cast<std::int32_t>(v));
        }
      }
    }
    return level_[t] >= 0;
  }

  // One augmenting path origin_ -> target_ in the level graph into path_ (current-arc pointers
  // kept: an arc that is pushed through stays current; a dead end advances its parent's pointer).
  template <class P>
  bool find_path(FlowCoreCounters& c) {
    const auto& pi = potentials<P>();
    stack_.clear();
    path_.clear();
    stack_.push_back(origin_);
    while (!stack_.empty()) {
      const auto u = static_cast<std::size_t>(stack_.back());
      if (static_cast<std::int32_t>(u) == target_) return true;
      bool advanced = false;
      const auto end = net_.off[u + 1];
      const auto next = level_[u] + 1;
      const P pu = pi[u];
      while (cur_[u] < end) {
        ++c.arc_scans;
        const auto a = static_cast<std::size_t>(cur_[u]);
        const auto v = static_cast<std::size_t>(net_.arc_head[a]);
        if (res_[a] > 0 && level_[v] == next && net_.arc_cost[a] + pu == pi[v]) {
          stack_.push_back(static_cast<std::int32_t>(v));
          path_.push_back(static_cast<std::int32_t>(a));
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

  // One fixed-origin repair of the deletion re-solve: tiers from `origin` to `target` until
  // `amount` is moved or the target is unreachable; with `reuse`, the first tier reuses the
  // last search (from the same origin, nothing changed since). Returns the amount moved.
  Amount repair(int index, std::int32_t origin, std::int32_t target, Amount amount, bool reuse, FlowCoreCounters& c) {
    origin_ = origin;
    target_ = target;
    Amount moved = 0;
    bool search = !reuse;
    while (moved < amount) {
      if (search) {
        const bool hit = dijkstra<I128>(c);
        update_all();
        emit(index, false, moved, c);
        if (!hit) break;
      } else if (!reached<I128>(target_)) {
        throw std::logic_error("deletion re-solve: the reused search did not reach the repair's target");
      }
      search = true;
      const Amount before = moved;
      moved = drain<I128>(c, amount, moved, [&](Amount total) { emit(index, true, total, c); });
      if (moved == before) throw std::logic_error("flow core: a successful Dijkstra drained nothing");
    }
    return moved;
  }

  // The full residual closure of the real sources (reverse and virtual arcs included) into
  // seen_. Returns whether t* is outside it.
  bool close_sources() {
    seen_.assign(static_cast<std::size_t>(n_), 0);
    queue_.clear();
    for (std::int32_t s : net_.sources) {
      seen_[static_cast<std::size_t>(s)] = 1;
      queue_.push_back(s);
    }
    for (std::size_t qi = 0; qi < queue_.size(); ++qi) {
      const auto u = static_cast<std::size_t>(queue_[qi]);
      for (auto p = net_.off[u]; p < net_.off[u + 1]; ++p) {
        const auto a = static_cast<std::size_t>(p);
        const auto v = static_cast<std::size_t>(net_.arc_head[a]);
        if (res_[a] > 0 && seen_[v] == 0) {
          seen_[v] = 1;
          queue_.push_back(static_cast<std::int32_t>(v));
        }
      }
    }
    return seen_[static_cast<std::size_t>(t_)] == 0;
  }

  // The input edges of positive view capacity from the closure to its complement (ascending),
  // the input edge `skip` (a deleted edge) left out.
  [[nodiscard]] std::vector<std::int32_t> cut_edges(std::int64_t skip) const {
    std::vector<std::int32_t> cut;
    for (std::size_t k = 0; k < net_.edges.size(); ++k) {
      if (static_cast<std::int64_t>(k) != skip && seen_[static_cast<std::size_t>(net_.tail[k])] != 0 &&
          seen_[static_cast<std::size_t>(net_.head[k])] == 0 && net_.cap[k] > 0) {
        cut.push_back(static_cast<std::int32_t>(k));
      }
    }
    return cut;
  }

#if NETGRAPH_CORE_PRIVATE_HOOKS
  void emit(int repair, bool push, Amount moved, const FlowCoreCounters& c) const {
    if (probe_ != nullptr && probe_->observe) probe_->observe(DeletionEvent {repair, push, moved, res_, wpi_, &c});
  }
#else
  void emit(int /*repair*/, bool /*push*/, Amount /*moved*/, const FlowCoreCounters& /*c*/) const noexcept {}
#endif

  Network net_ {};
  std::int32_t n_ {0};
  std::int32_t s_ {0};
  std::int32_t t_ {0};
  std::int32_t origin_ {0};
  std::int32_t target_ {0};
  std::vector<std::int32_t> pos_ {};
  std::vector<Amount> res_ {};
  // A solve from zero.
  std::vector<Cost> pi_ {};
  std::vector<std::uint8_t> ever_ {};
  std::vector<Cost> dist_ {};
  RadixHeap<Cost> heap_ {};
  bool shadow_ {false};
  I128 tier_sum_ {};
  std::vector<I128> checkpoint_ {};
  // The deletion mode.
  std::vector<I128> wpi_ {};
  std::vector<I128> wdist_ {};
  std::vector<std::uint32_t> wseen_ {};
  std::uint32_t wepoch_ {0};
  RadixHeap<I128> wheap_ {};
#if NETGRAPH_CORE_PRIVATE_HOOKS
  const DeletionProbe* probe_ {nullptr};
  const std::function<void(std::span<const Cost>)>* on_search_ {nullptr};
#endif
  // Shared scratch.
  std::vector<std::int32_t> settled_ {};
  std::vector<std::int32_t> level_ {};
  std::vector<std::uint32_t> reach_ {};
  std::uint32_t mark_epoch_ {0};
  std::vector<std::int32_t> queue_ {};
  std::vector<std::int32_t> cur_ {};
  std::vector<std::int32_t> stack_ {};
  std::vector<std::int32_t> path_ {};
  std::vector<std::uint8_t> seen_ {};
};

std::vector<Amount> DeletionState::edge_flows() const {
  std::vector<Amount> flows(static_cast<std::size_t>(graph_->num_edges()), 0);
  for (std::size_t k = 0; k < net_->edges.size(); ++k) {
    const auto fwd = static_cast<std::size_t>(net_->edge_slot[k]);
    flows[static_cast<std::size_t>(net_->edges[k])] = res_[static_cast<std::size_t>(net_->rev[fwd])];
  }
  return flows;
}

DeletionState& DeletionState::operator=(const DeletionState& other) {
  DeletionState copy(other);  // every allocation; *this is untouched until it succeeds
  return *this = std::move(copy);
}

DeletionResolver::DeletionResolver(const DeletionState& state) : state_(state), solver_(std::make_unique<FlowSolver>()) {
  if (state_.net_ == nullptr) throw std::invalid_argument("deletion re-solve: an empty state (no solve captured it)");
  solver_->load(state_.network());
}

DeletionResolver::~DeletionResolver() = default;
DeletionResolver::DeletionResolver(DeletionResolver&&) noexcept = default;
DeletionResolver& DeletionResolver::operator=(DeletionResolver&&) noexcept = default;

DeletionResult DeletionResolver::resolve(EdgeId e) { return solver_->resolve(state_, index_of(e)); }

#if NETGRAPH_CORE_PRIVATE_HOOKS
DeletionResult DeletionResolver::resolve(EdgeId e, const DeletionProbe& probe) {
  const std::size_t k = index_of(e);
  struct Reset {
    FlowSolver& solver;
    ~Reset() { solver.probe(nullptr); }
  } reset {*solver_};
  solver_->probe(&probe);
  return solver_->resolve(state_, k);
}
#endif

std::size_t DeletionResolver::index_of(EdgeId e) const {
  if (!state_.maximum()) {
    throw std::invalid_argument(
        "deletion re-solve: the state is not certified maximum (a residual path from the sources to the "
        "destinations remains)");
  }
  const Network& net = state_.network();
  const auto it = std::lower_bound(net.edges.begin(), net.edges.end(), e);
  if (it == net.edges.end() || *it != e) {
    throw std::invalid_argument("deletion re-solve: edge " + std::to_string(e) + " is not in the solve's view");
  }
  const auto k = static_cast<std::size_t>(it - net.edges.begin());
  const auto fwd = static_cast<std::size_t>(net.edge_slot[k]);
  if (state_.residual()[static_cast<std::size_t>(net.rev[fwd])] <= 0) {
    throw std::invalid_argument("deletion re-solve: edge " + std::to_string(e) + " carries no flow");
  }
  return k;
}

FlowCoreResult place_global_units(const Graph& g,
                                  std::span<const Amount> view,
                                  std::span<const NodeId> sources,
                                  std::span<const NodeId> destinations,
                                  const Eligibility* elig,
                                  std::optional<Amount> cap,
                                  DeletionCapture* capture) {
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
  const auto srcs = node_set(sources, g.num_nodes(), "sources");
  const auto dsts = node_set(destinations, g.num_nodes(), "destinations");
  check_disjoint(srcs, dsts, "sources", "destinations");
  if (cap.has_value() && (*cap < 0 || *cap > kMaxFlowCap)) {
    throw std::invalid_argument("cap " + std::to_string(*cap) + " is outside [0, 2^63 - 2]");
  }
  if (elig != nullptr && elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  if (capture != nullptr && capture->graph.get() != &g) {
    throw std::logic_error("flow core: a capture request bound to another graph");
  }
  FlowCoreResult out;
  const bool zero = cap.has_value() && *cap == 0;
  if (zero) {
    // The input cap 0 runs nothing: the prescribed zero payload.
    out.edge_flows.assign(num_edges, 0);
    out.potentials.assign(static_cast<std::size_t>(g.num_nodes()), kNeverReached);
    if (capture == nullptr) return out;
  }
  // Identity coordinates: every real vertex, the allowed edges and endpoints.
  FlowSolver& solver = FlowSolver::local();
  Network& in = solver.network();
  in.num_nodes = g.num_nodes();
  const auto edge_src = g.edge_src();
  const auto edge_dst = g.edge_dst();
  const auto cost = g.cost();
  in.edges.clear();
  in.tail.clear();
  in.head.clear();
  in.cost.clear();
  in.cap.clear();
  in.sources.clear();
  in.destinations.clear();
  in.edges.reserve(num_edges);
  for (std::size_t e = 0; e < num_edges; ++e) {
    if (elig != nullptr && !elig->edge_ok(static_cast<EdgeId>(e))) continue;
    in.edges.push_back(static_cast<EdgeId>(e));
    in.tail.push_back(edge_src[e]);
    in.head.push_back(edge_dst[e]);
    in.cost.push_back(cost[e]);
    in.cap.push_back(view[e]);
  }
  for (NodeId s : srcs) {
    if (elig == nullptr || elig->node_ok(s)) in.sources.push_back(s);
  }
  for (NodeId d : dsts) {
    if (elig == nullptr || elig->node_ok(d)) in.destinations.push_back(d);
  }
  solver.build(capture != nullptr);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  solver.observe(capture != nullptr ? &capture->on_search : nullptr);
#endif
  if (!zero) {
    auto r = solver.run(cap.value_or(kUncapped), /*closure=*/true);
    out.placed = r.placed;
    out.cost_distribution = std::move(r.cost_distribution);
    out.edge_flows.assign(num_edges, 0);
    for (std::size_t k = 0; k < in.edges.size(); ++k) {
      out.edge_flows[static_cast<std::size_t>(in.edges[k])] = r.flows[k];
    }
    if (r.cut.has_value()) {
      std::vector<EdgeId> cut;
      cut.reserve(r.cut->size());
      for (std::int32_t k : *r.cut) cut.push_back(in.edges[static_cast<std::size_t>(k)]);
      out.min_cut = std::move(cut);
    }
    out.potentials = std::move(r.potentials);
    out.counters = r.counters;
  }
  if (capture != nullptr) solver.capture(*capture, out.placed, out.counters);
  return out;
}

SubsetFlow solve_subset(const Graph& g, std::span<const EdgeId> edges, std::span<const Amount> caps,
                        std::span<const NodeId> sources, std::span<const NodeId> destinations,
                        std::optional<Amount> cap) {
  SubsetFlow out;
  out.flows.assign(edges.size(), 0);
  if (cap.has_value() && *cap == 0) return out;
  const auto edge_src = g.edge_src();
  const auto edge_dst = g.edge_dst();
  const auto cost = g.cost();
  // Compressed vertices in first-touch order. The numbering changes no flow, curve or path:
  // the arc order at every vertex is by EdgeId (virtual arc first) and the virtual arcs of s*
  // and t* follow the endpoint lists, which stay in ascending NodeId.
  NodeMap& map = NodeMap::local();
  map.reset(static_cast<std::size_t>(g.num_nodes()));
  // The input's storage is the thread's solver's (every field is rewritten here).
  FlowSolver& solver = FlowSolver::local();
  Network& in = solver.network();
  in.edges.assign(edges.begin(), edges.end());
  in.tail.clear();
  in.head.clear();
  in.cost.clear();
  in.sources.clear();
  in.destinations.clear();
  for (EdgeId e : edges) {
    const auto ei = static_cast<std::size_t>(e);
    in.tail.push_back(map.add(edge_src[ei]));
    in.head.push_back(map.add(edge_dst[ei]));
    in.cost.push_back(cost[ei]);
  }
  for (NodeId v : sources) in.sources.push_back(map.add(v));
  for (NodeId v : destinations) in.destinations.push_back(map.add(v));
  in.num_nodes = static_cast<std::int32_t>(map.size());
  in.cap.assign(caps.begin(), caps.end());
  solver.build(false);
  auto r = solver.run(cap.value_or(kUncapped), /*closure=*/false, /*potentials=*/false);
  out.placed = r.placed;
  out.flows = std::move(r.flows);
  out.cost_distribution = std::move(r.cost_distribution);
  return out;
}

} // namespace netgraph::core::detail
