/*
  The one integer minimum-cost flow core.

  Find the maximum flow allowed by the cap at minimum cost. The residual network records how
  much flow a forward arc can add and how much its reverse arc can cancel. A potential pi is
  a vertex cost offset; the reduced cost of u -> v is cost(e) + pi(u) - pi(v). These offsets
  make searchable arc costs non-negative without changing which path between two endpoints
  is cheapest. Each tier runs Dijkstra, which finalises the nearest pending vertex first,
  over all residual vertices reachable from virtual source s*, through sink t* and onward.
  Add the final distances to settled vertices' potentials. Arcs on shortest paths now have
  zero reduced cost; those with positive residual capacity form the admissible subgraph.
  Dinic phases push flow through it: breadth-first search (BFS) assigns minimum-hop levels;
  depth-first search (DFS) follows level-increasing paths, keeping each node's next arc to try.
  A blocking flow fills these paths until none can take more. Repeat phases until no admissible
  augmenting path (a residual source-to-sink path) remains, or the cap is reached. Repetition is
  needed because a phase's level graph keeps only the admissible paths with the fewest arcs.
  Start the next tier; stop when Dijkstra cannot reach t* or the cap is reached. No completion
  phase follows.

  Arc identity: each eligible edge has a forward arc with the view capacity and a reverse arc
  with its flow. At each node, visit forward out-arcs and reverse in-arcs together by ascending
  EdgeId; for a self-loop visit forward before reverse. Virtual vertices s* and t* connect the
  real source and destination sets. A real endpoint's virtual arc comes first; s* visits
  sources and t* visits destinations by ascending NodeId. Virtual arcs cost 0 and carry
  2^63 - 2 units, above every attainable value, unless the solve gives the sources' arcs their own
  capacities (the equal source allocation: source i's arc carries at most its round-robin share). Starting
  from zero flow makes zero potentials valid: usable real arcs have positive costs and virtual
  arcs cost 0.

  Deletion re-solves start from a captured final state. During capture, track working
  potentials without changing the original searches. Keep the sum of tier maxima D_k and
  each settled vertex's last such sum; D_k is that search's largest finite distance.
  Working potential = published potential + tier maxima accumulated since its last settlement.
  For a never-settled vertex, start at zero and add the whole sum. This equals updating every
  vertex after every Dijkstra: reached vertices by their distance, the others by D_k.
  Thus every positive residual arc, real or virtual, has non-negative reduced cost.
  Delete e = (u, v), with flow f_e > 0, from a copy and run three repairs, each with a fixed
  origin and target and that all-vertex update before Dinic phases. First move up to f_e from
  u to v; call the amount moved m. If r = f_e - m > 0, move r from u to s*, then from t* to v.
  The u -> s* repair reuses the first repair's failed search because nothing has changed.
  The result has maximum value (the solved value minus r) and minimum cost on the view without e. Working
  potentials certify its cost; the real sources' full residual closure certifies its value.
  There are at most m + 2 r Dijkstra runs.

  Dijkstra's queue is a radix heap: keys pop in increasing order and entries with equal keys pop
  in a fixed order (last pushed first within the lowest bucket), so the settle order, and with it
  every counter, is a function of the input alone. Equal-key order affects only the counters,
  never a distance, a potential or an augmenting path. The solve runs to the cap or to
  disconnection; there is no clamp of the cap.

  Arithmetic. A solve from zero stays in int64: each true residual distance is the cost of a
  simple path (no repeated vertex), below the cost sum bound 2^62; potentials are such distances.
  A path's bottleneck is its smallest residual capacity. It includes a real arc because sources
  and destinations are disjoint, so every amount is below 2^62. The deletion mode's working
  potentials, distance keys, tentative distances and sums of tier maxima are signed 128-bit
  (I128; reachability kept apart from the value by a stamp): they stay within (n + 1) C at
  capture and (22 n + 15) C across the repairs for n vertices and the cost sum C < 2^62, and its
  heap takes 128-bit keys. Amounts stay int64.

  Size. Vertex ids, arc slots and slot offsets are int32. A network of N real vertices has
  N + 2 vertices, and every edge and every virtual arc has two slots (forward and reverse), so
  N <= INT32_MAX - 2 and edges + sources + destinations <= 2^30 - 1 (twice that is the largest
  even int32). network_vertices() refuses a larger network; both ways of sizing the solver (a
  solve from zero, a captured state's re-solves) take their vertex count from it before anything
  else, so every vertex id, slot and offset fits.
*/
#include "netgraph/core/detail/min_cost_flow.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "netgraph/core/detail/edge_amounts.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/node_map.hpp"

namespace netgraph::core::detail {

namespace {

constexpr Cost kInfDist = std::numeric_limits<Cost>::max();
constexpr Amount kUncapped = std::numeric_limits<Amount>::max();

// The bucket of a key: the number of the highest bit in which it differs from the last popped
// key (keys are non-negative). A 128-bit key differs in its high word from 2^64 on, which the
// deletion mode's bound (22 n + 15) C allows (see Arithmetic above).
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

std::int32_t network_vertices(std::int32_t num_nodes, std::size_t num_edges, std::size_t num_sources,
                              std::size_t num_destinations) {
  constexpr std::int32_t kMax = std::numeric_limits<std::int32_t>::max();
  if (num_nodes > kMax - 2) {
    throw std::invalid_argument("flow solver node limit exceeded");
  }
  // Each term is compared with what the earlier ones leave, so no sum of counts can wrap.
  constexpr auto kPairs = static_cast<std::size_t>(kMax / 2);  // 2^30 - 1
  if (num_edges > kPairs || num_sources > kPairs - num_edges || num_destinations > kPairs - num_edges - num_sources) {
    throw std::invalid_argument("flow solver arc limit exceeded");
  }
  return num_nodes + 2;
}

// The residual network of one solve (Network: real vertices 0..N-1, s* = N, t* = N + 1) with
// its flow, potentials and scratch; also the working storage of a deletion re-solve.
class MinCostFlowSolver {
public:
  // One reusable instance per thread (`local()`) for the solves from zero: build() resizes every
  // buffer for the input, keeping capacity across calls, so the solver's own buffers are not
  // reallocated once warm (each solve still allocates its outputs). Not reentrant (a solve never
  // starts another). A deletion resolver owns its own.
  static MinCostFlowSolver& local() {
    thread_local MinCostFlowSolver solver;
    return solver;
  }

  MinCostFlowSolver() = default;
  // Never copied or moved: net_ may point at own_.
  MinCostFlowSolver(const MinCostFlowSolver&) = delete;
  MinCostFlowSolver& operator=(const MinCostFlowSolver&) = delete;

  // The input part of the network (num_nodes, the edges, the endpoint sets) is the caller's to
  // fill; build() adds the arcs.
  ResidualNetwork& network() noexcept { return owned_network_; }

  // The arcs of the input, zero flow and zero potentials. `shadow`: record the working
  // potentials for a capture.
  void build(bool shadow) {
    auto& in = owned_network_;
    size(in);
    network_ = &owned_network_;
    const auto n = static_cast<std::size_t>(num_vertices_);
    in.arc_offsets.assign(n + 1, 0);
    auto deg = [&](std::int32_t v) { ++in.arc_offsets[static_cast<std::size_t>(v) + 1]; };
    for (std::size_t k = 0; k < in.edges.size(); ++k) {
      deg(in.edge_sources[k]);
      deg(in.edge_destinations[k]);
    }
    for (std::int32_t s : in.sources) { deg(s); deg(virtual_source_); }
    for (std::int32_t d : in.destinations) { deg(d); deg(virtual_sink_); }
    for (std::size_t v = 0; v < n; ++v) in.arc_offsets[v + 1] += in.arc_offsets[v];
    const auto num_slots = static_cast<std::size_t>(in.arc_offsets[n]);
    in.arc_destinations.resize(num_slots);
    in.reverse_arcs.resize(num_slots);
    in.arc_cost.resize(num_slots);
    arc_residual_.resize(num_slots);
    arc_fill_positions_.assign(in.arc_offsets.begin(), in.arc_offsets.end() - 1);
    auto take = [&](std::int32_t v) { return arc_fill_positions_[static_cast<std::size_t>(v)]++; };
    auto link = [&](std::int32_t a, std::int32_t to_a, Amount cap, Cost w, std::int32_t b, std::int32_t to_b) {
      const auto ai = static_cast<std::size_t>(a);
      const auto bi = static_cast<std::size_t>(b);
      in.arc_destinations[ai] = to_a;
      in.reverse_arcs[ai] = b;
      in.arc_cost[ai] = w;
      arc_residual_[ai] = cap;
      in.arc_destinations[bi] = to_b;
      in.reverse_arcs[bi] = a;
      in.arc_cost[bi] = -w;
      arc_residual_[bi] = 0;
    };
    // Virtual arcs first at the real endpoints, then the real arcs by ascending EdgeId (forward
    // slot at the tail, reverse slot at the head), then the slots of s* and t*.
    in.source_reverse_arcs.resize(in.sources.size());
    in.destination_forward_arcs.resize(in.destinations.size());
    for (std::size_t i = 0; i < in.sources.size(); ++i) in.source_reverse_arcs[i] = take(in.sources[i]);
    for (std::size_t i = 0; i < in.destinations.size(); ++i) in.destination_forward_arcs[i] = take(in.destinations[i]);
    in.edge_forward_arcs.resize(in.edges.size());
    for (std::size_t k = 0; k < in.edges.size(); ++k) {
      const std::int32_t u = in.edge_sources[k];
      const std::int32_t v = in.edge_destinations[k];
      const std::int32_t fwd = take(u);
      const std::int32_t bwd = take(v);
      link(fwd, v, in.edge_capacities[k], in.cost[k], bwd, u);
      in.edge_forward_arcs[k] = fwd;
    }
    for (std::size_t i = 0; i < in.sources.size(); ++i) {
      link(take(virtual_source_), in.sources[i], in.source_capacities.empty() ? kMaxFlowCap : in.source_capacities[i], 0, in.source_reverse_arcs[i], virtual_source_);
    }
    for (std::size_t i = 0; i < in.destinations.size(); ++i) {
      link(in.destination_forward_arcs[i], virtual_sink_, kMaxFlowCap, 0, take(virtual_sink_), in.destinations[i]);
    }
    potentials_.assign(n, 0);
    ever_settled_.assign(n, 0);
    distances_.resize(n);
    shadow_ = shadow;
    if (shadow) {
      tier_sum_ = 0;
      checkpoint_.assign(n, I128 {});
    }
  }

  // A captured state's network, read in place by its deletion re-solves (never written; the
  // caller keeps it alive while this solver uses it). The residual, the potentials and the
  // scratch stay this solver's own.
  void load(const ResidualNetwork& net) {
    size(net);
    network_ = &net;
    arc_residual_.resize(network_->arc_destinations.size());
    repair_potentials_.resize(static_cast<std::size_t>(num_vertices_));
    repair_distances_.resize(static_cast<std::size_t>(num_vertices_));
    repair_seen_generations_.assign(static_cast<std::size_t>(num_vertices_), 0);
    repair_generation_ = 0;
  }

  // The solve to the cap or to disconnection. Output in local coordinates: `flows` aligned with
  // the input edges, `cut` as input edge indices (every edge from the full residual closure of
  // the real sources to its complement, ascending EdgeId) when the destinations are
  // outside the closure, potentials per local vertex (kNeverReached when never settled).
  struct Output {
    Amount flow_value {0};
    std::vector<std::pair<Cost, Amount>> marginal_cost_curve {};
    std::vector<Amount> edge_flows {};
    std::optional<std::vector<std::int32_t>> min_cut_positions {};
    std::vector<Cost> potentials {};
    MinCostFlowCounters counters {};
  };

  // `closure` adds the cut; `potentials` the per-vertex potentials (the edge-subset entry
  // needs neither).
  Output run(Amount limit, bool closure, bool potentials = true) {
    Output out;
    auto& c = out.counters;
    Amount total = 0;
    search_source_ = virtual_source_;
    search_target_ = virtual_sink_;
    while (total < limit) {
      const bool reached = dijkstra<Cost>(c);
      settle();
      if (!reached) break;
      const Cost tier_cost = potentials_[static_cast<std::size_t>(virtual_sink_)] - potentials_[static_cast<std::size_t>(virtual_source_)];
      const Amount before = total;
      total = drain<Cost>(c, limit, total, [](Amount) {});
      // A successful Dijkstra reaches t* over arcs with positive residual and zero reduced cost
      // along its shortest paths, so the drain moves at least one unit.
      if (total == before) {
        throw std::logic_error("flow core: a successful Dijkstra drained nothing");
      }
      out.marginal_cost_curve.emplace_back(tier_cost, total - before);
    }
    out.flow_value = total;

    out.edge_flows.resize(network_->edges.size());
    for (std::size_t k = 0; k < network_->edges.size(); ++k) out.edge_flows[k] = flow_of(k);
    if (closure && close_sources()) out.min_cut_positions = cut_edges(-1);
    if (potentials) {
      out.potentials.resize(static_cast<std::size_t>(network_->num_nodes));
      for (std::size_t v = 0; v < out.potentials.size(); ++v) {
        out.potentials[v] = ever_settled_[v] != 0 ? potentials_[v] : kNeverReached;
      }
    }
    return out;
  }

  // The final state of the last solve (its value and counters given) into `req`.
  void capture(DeletionCapture& req, Amount value, const MinCostFlowCounters& counters) {
    DeletionState state;
    state.graph_ = req.graph;
    state.network_ = std::make_shared<const ResidualNetwork>(*network_);
    state.arc_residual_ = arc_residual_;
    state.potentials_.resize(static_cast<std::size_t>(num_vertices_));
    for (std::size_t v = 0; v < state.potentials_.size(); ++v) state.potentials_[v] = I128 {potentials_[v]} + (tier_sum_ - checkpoint_[v]);
    state.value_ = value;
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
    std::copy(state.residual().begin(), state.residual().end(), arc_residual_.begin());
    std::copy(state.potentials().begin(), state.potentials().end(), repair_potentials_.begin());
    DeletionResult out;
    auto& c = out.counters;
    const auto fwd = static_cast<std::size_t>(network_->edge_forward_arcs[k]);
    const auto bwd = static_cast<std::size_t>(network_->reverse_arcs[fwd]);
    const Amount f = arc_residual_[bwd];
    const std::int32_t u = network_->edge_sources[k];
    const std::int32_t v = network_->edge_destinations[k];
    // Removing e leaves f extra units at u and f missing units at v; conservation is broken
    // until repairs reconnect or cancel them. Remove both residual arcs of e.
    arc_residual_[fwd] = 0;
    arc_residual_[bwd] = 0;
    emit(0, false, 0, c);
    const Amount m = repair(1, u, v, f, false, c);
    const Amount r = f - m;
    if (r > 0) {
      bool reuse = true;
#if NETGRAPH_CORE_PRIVATE_HOOKS
      if (probe_ != nullptr && probe_->rerun_failed_search) reuse = false;
#endif
      if (repair(2, u, virtual_source_, r, reuse, c) != r) {
        throw std::logic_error("deletion re-solve: the remainder has no residual path back to the sources");
      }
      if (repair(3, virtual_sink_, v, r, false, c) != r) {
        throw std::logic_error("deletion re-solve: the deficit has no residual path from the destinations");
      }
    }
    emit(4, false, 0, c);
    if (!close_sources()) {
      throw std::logic_error("deletion re-solve: the destinations stay connected to the sources");
    }
    out.value = state.value() - r;
    out.potentials = repair_potentials_;
    for (std::int32_t j : cut_edges(static_cast<std::int64_t>(k))) out.min_cut.push_back(network_->edges[static_cast<std::size_t>(j)]);
    out.edge_flows.assign(static_cast<std::size_t>(state.graph().num_edges()), 0);
    for (std::size_t j = 0; j < network_->edges.size(); ++j) {
      out.edge_flows[static_cast<std::size_t>(network_->edges[j])] = j == k ? 0 : flow_of(j);
    }
    return out;
  }

private:
  // The vertex count and the virtual vertices of `net` (refused beyond the limits of
  // network_vertices(), before any change), then the scratch of that many vertices.
  void size(const ResidualNetwork& net) {
    num_vertices_ = network_vertices(net.num_nodes, net.edges.size(), net.sources.size(), net.destinations.size());
    virtual_source_ = num_vertices_ - 2;
    virtual_sink_ = num_vertices_ - 1;
    const auto n = static_cast<std::size_t>(num_vertices_);
    level_.resize(n);
    reach_.assign(n, 0);
    reachability_generation_ = 0;
    current_arcs_.resize(n);
  }

  // The flow of input edge k: the residual of its reverse slot.
  [[nodiscard]] Amount flow_of(std::size_t k) const {
    return arc_residual_[static_cast<std::size_t>(network_->reverse_arcs[static_cast<std::size_t>(network_->edge_forward_arcs[k])])];
  }

  template <class P>
  std::vector<P>& potentials() {
    if constexpr (kWide<P>) {
      return repair_potentials_;
    } else {
      return potentials_;
    }
  }
  template <class P>
  std::vector<P>& distances() {
    if constexpr (kWide<P>) {
      return repair_distances_;
    } else {
      return distances_;
    }
  }
  template <class P>
  RadixHeap<P>& heap() {
    if constexpr (kWide<P>) {
      return repair_heap_;
    } else {
      return heap_;
    }
  }
  // Whether the last search of the mode reached x (int64: a finite distance; I128: the stamp).
  template <class P>
  [[nodiscard]] bool reached(std::int32_t x) const {
    const auto xi = static_cast<std::size_t>(x);
    if constexpr (kWide<P>) {
      return repair_seen_generations_[xi] == repair_generation_;
    } else {
      return distances_[xi] != kInfDist;
    }
  }

  // Dijkstra on reduced costs from origin_ to completion over the residual network. Returns
  // whether target_ was reached; settled_ lists the vertices it settled.
  template <class P>
  bool dijkstra(MinCostFlowCounters& c) {
    auto& pi = potentials<P>();
    auto& dist = distances<P>();
    auto& queue = heap<P>();
    // The network is never resized during a search: raw pointers keep the compiler from
    // reloading net_ and the vector bases around every store into the solver's own buffers.
    const std::int32_t* const off = network_->arc_offsets.data();
    const std::int32_t* const arc_destinations = network_->arc_destinations.data();
    const Cost* const arc_cost = network_->arc_cost.data();
    ++c.dijkstra_runs;
    if constexpr (kWide<P>) {
      // Wraparound after 2^32 wide searches since load(): one resolver's searches over all its
      // deletions, which no stated limit keeps below 2^32.
      if (++repair_generation_ == 0) {
        std::fill(repair_seen_generations_.begin(), repair_seen_generations_.end(), 0);
        repair_generation_ = 1;
      }
    } else {
      std::fill(dist.begin(), dist.end(), kInfDist);
    }
    settled_.clear();
    queue.clear();
    const auto o = static_cast<std::size_t>(search_source_);
    dist[o] = P {0};
    if constexpr (kWide<P>) repair_seen_generations_[o] = repair_generation_;
    queue.push(P {0}, search_source_);
    while (!queue.empty()) {
      const auto [du, ui] = queue.pop();
      ++c.dijkstra_pops;
      const auto u = static_cast<std::size_t>(ui);
      if (du != dist[u]) continue;
      settled_.push_back(ui);
      ++c.settled;
      const P base = du + pi[u];
      const auto end = off[u + 1];
      c.arc_scans += end - off[u];
      for (auto p = off[u]; p < end; ++p) {
        const auto a = static_cast<std::size_t>(p);
        if (arc_residual_[a] <= 0) continue;
        const auto v = static_cast<std::size_t>(arc_destinations[a]);
        const P nd = base + arc_cost[a] - pi[v];
        if constexpr (kWide<P>) {
          if (repair_seen_generations_[v] != repair_generation_ || nd < dist[v]) {
            repair_seen_generations_[v] = repair_generation_;
            dist[v] = nd;
            queue.push(nd, static_cast<std::int32_t>(v));
          }
        } else if (nd < dist[v]) {
          dist[v] = nd;
          queue.push(nd, static_cast<std::int32_t>(v));
        }
      }
    }
    return reached<P>(search_target_);
  }

  // The update of a solve from zero: the settled vertices' potentials by their distances; with
  // a capture, the running sum of the tier maxima and the settled vertices' checkpoints.
  void settle() {
    for (std::int32_t v : settled_) {
      const auto vi = static_cast<std::size_t>(v);
      potentials_[vi] += distances_[vi];
      ever_settled_[vi] = 1;
    }
    if (shadow_) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
      if (on_search_ != nullptr && *on_search_) (*on_search_)(distances_);
#endif
      Cost top = 0;
      for (std::int32_t v : settled_) top = std::max(top, distances_[static_cast<std::size_t>(v)]);
      tier_sum_ += top;
      for (std::int32_t v : settled_) checkpoint_[static_cast<std::size_t>(v)] = tier_sum_;
    }
  }

  // The update of the deletion mode: every vertex, the reached ones by their distance and the
  // others by the tier's largest finite distance.
  void update_all() {
    I128 top = 0;
    for (std::int32_t v : settled_) top = std::max(top, repair_distances_[static_cast<std::size_t>(v)]);
    for (std::size_t v = 0; v < repair_potentials_.size(); ++v) repair_potentials_[v] += repair_seen_generations_[v] == repair_generation_ ? repair_distances_[v] : top;
  }

  // The Dinic phases from origin_ to target_ over the admissible subgraph (positive residual,
  // zero reduced cost) until no admissible path remains or `limit` is reached; on_push(total)
  // after every push.
  template <class P, class OnPush>
  Amount drain(MinCostFlowCounters& c, Amount limit, Amount total, OnPush on_push) {
    const std::int32_t* const rev = network_->reverse_arcs.data();
    while (total < limit) {
      if (!build_levels<P>(c)) break;
      std::copy(network_->arc_offsets.begin(), network_->arc_offsets.end() - 1, current_arcs_.begin());
      while (total < limit) {
        if (!find_path<P>(c)) break;
        Amount amount = limit - total;
        for (std::int32_t a : path_) amount = std::min(amount, arc_residual_[static_cast<std::size_t>(a)]);
        for (std::int32_t a : path_) {
          arc_residual_[static_cast<std::size_t>(a)] -= amount;
          arc_residual_[static_cast<std::size_t>(rev[static_cast<std::size_t>(a)])] += amount;
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
  bool build_levels(MinCostFlowCounters& c) {
    const auto& pi = potentials<P>();
    const std::int32_t* const off = network_->arc_offsets.data();
    const std::int32_t* const arc_destinations = network_->arc_destinations.data();
    const std::int32_t* const rev = network_->reverse_arcs.data();
    const Cost* const arc_cost = network_->arc_cost.data();
    const auto t = static_cast<std::size_t>(search_target_);
    // Wraparound after 2^32 phases since size(): one solve from zero, or one resolver's phases
    // over all its deletions, which no stated limit keeps below 2^32.
    if (++reachability_generation_ == 0) {
      std::fill(reach_.begin(), reach_.end(), 0);
      reachability_generation_ = 1;
    }
    // Reverse search from the target: the arc into v from arc_destinations[p] is the partner of v's slot p.
    queue_.clear();
    reach_[t] = reachability_generation_;
    queue_.push_back(search_target_);
    for (std::size_t qi = 0; qi < queue_.size(); ++qi) {
      const auto v = static_cast<std::size_t>(queue_[qi]);
      const P pv = pi[v];
      const auto end = off[v + 1];
      c.arc_scans += end - off[v];
      for (auto p = off[v]; p < end; ++p) {
        const auto a = static_cast<std::size_t>(p);
        const auto u = static_cast<std::size_t>(arc_destinations[a]);
        if (reach_[u] == reachability_generation_) continue;
        const auto in = static_cast<std::size_t>(rev[a]);
        if (arc_residual_[in] > 0 && arc_cost[in] + pi[u] == pv) {
          reach_[u] = reachability_generation_;
          queue_.push_back(static_cast<std::int32_t>(u));
        }
      }
    }
    std::fill(level_.begin(), level_.end(), -1);
    ++c.dinic_phases;
    const auto s = static_cast<std::size_t>(search_source_);
    if (reach_[s] != reachability_generation_) return false;
    queue_.clear();
    level_[s] = 0;
    queue_.push_back(search_source_);
    for (std::size_t qi = 0; qi < queue_.size(); ++qi) {
      const auto u = static_cast<std::size_t>(queue_[qi]);
      if (level_[t] >= 0 && level_[u] >= level_[t]) break;
      const P pu = pi[u];
      const auto next = level_[u] + 1;
      const auto end = off[u + 1];
      c.arc_scans += end - off[u];
      for (auto p = off[u]; p < end; ++p) {
        const auto a = static_cast<std::size_t>(p);
        const auto v = static_cast<std::size_t>(arc_destinations[a]);
        if (level_[v] < 0 && reach_[v] == reachability_generation_ && arc_residual_[a] > 0 && arc_cost[a] + pu == pi[v]) {
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
  bool find_path(MinCostFlowCounters& c) {
    const auto& pi = potentials<P>();
    const std::int32_t* const off = network_->arc_offsets.data();
    const std::int32_t* const arc_destinations = network_->arc_destinations.data();
    const Cost* const arc_cost = network_->arc_cost.data();
    stack_.clear();
    path_.clear();
    stack_.push_back(search_source_);
    while (!stack_.empty()) {
      const auto u = static_cast<std::size_t>(stack_.back());
      if (static_cast<std::int32_t>(u) == search_target_) return true;
      bool advanced = false;
      const auto end = off[u + 1];
      const auto next = level_[u] + 1;
      const P pu = pi[u];
      while (current_arcs_[u] < end) {
        ++c.arc_scans;
        const auto a = static_cast<std::size_t>(current_arcs_[u]);
        const auto v = static_cast<std::size_t>(arc_destinations[a]);
        if (arc_residual_[a] > 0 && level_[v] == next && arc_cost[a] + pu == pi[v]) {
          stack_.push_back(static_cast<std::int32_t>(v));
          path_.push_back(static_cast<std::int32_t>(a));
          advanced = true;
          break;
        }
        ++current_arcs_[u];
      }
      if (!advanced) {
        stack_.pop_back();
        if (!path_.empty()) {
          path_.pop_back();
          ++current_arcs_[static_cast<std::size_t>(stack_.back())];
        }
      }
    }
    return false;
  }

  // One fixed-origin repair of the deletion re-solve: tiers from `origin` to `target` until
  // `amount` is moved or the target is unreachable; with `reuse`, the first tier reuses the
  // last search (from the same origin, nothing changed since). Returns the amount moved.
  Amount repair(int index, std::int32_t origin, std::int32_t target, Amount amount, bool reuse, MinCostFlowCounters& c) {
    search_source_ = origin;
    search_target_ = target;
    Amount moved = 0;
    bool search = !reuse;
    while (moved < amount) {
      if (search) {
        const bool hit = dijkstra<I128>(c);
        update_all();
        emit(index, false, moved, c);
        if (!hit) break;
      } else if (!reached<I128>(search_target_)) {
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
    const std::int32_t* const off = network_->arc_offsets.data();
    const std::int32_t* const arc_destinations = network_->arc_destinations.data();
    seen_.assign(static_cast<std::size_t>(num_vertices_), 0);
    queue_.clear();
    for (std::int32_t s : network_->sources) {
      seen_[static_cast<std::size_t>(s)] = 1;
      queue_.push_back(s);
    }
    for (std::size_t qi = 0; qi < queue_.size(); ++qi) {
      const auto u = static_cast<std::size_t>(queue_[qi]);
      for (auto p = off[u]; p < off[u + 1]; ++p) {
        const auto a = static_cast<std::size_t>(p);
        const auto v = static_cast<std::size_t>(arc_destinations[a]);
        if (arc_residual_[a] > 0 && seen_[v] == 0) {
          seen_[v] = 1;
          queue_.push_back(static_cast<std::int32_t>(v));
        }
      }
    }
    return seen_[static_cast<std::size_t>(virtual_sink_)] == 0;
  }

  // The input edges from the closure to its complement (ascending), whatever their view
  // capacity: an edge of view 0 there adds nothing to the cut's value and is listed because every
  // edge leaving the closure is; the input edge `skip` (a deleted edge) left out.
  [[nodiscard]] std::vector<std::int32_t> cut_edges(std::int64_t skip) const {
    const std::int32_t* const tail = network_->edge_sources.data();
    const std::int32_t* const head = network_->edge_destinations.data();
    std::vector<std::int32_t> cut;
    for (std::size_t k = 0; k < network_->edges.size(); ++k) {
      if (static_cast<std::int64_t>(k) != skip && seen_[static_cast<std::size_t>(tail[k])] != 0 &&
          seen_[static_cast<std::size_t>(head[k])] == 0) {
        cut.push_back(static_cast<std::int32_t>(k));
      }
    }
    return cut;
  }

#if NETGRAPH_CORE_PRIVATE_HOOKS
  void emit(int repair, bool push, Amount moved, const MinCostFlowCounters& c) const {
    if (probe_ != nullptr && probe_->observe) probe_->observe(DeletionEvent {repair, push, moved, arc_residual_, repair_potentials_, &c});
  }
#else
  void emit(int /*repair*/, bool /*push*/, Amount /*moved*/, const MinCostFlowCounters& /*c*/) const noexcept {}
#endif

  ResidualNetwork owned_network_ {};              // a solve from zero's network (build())
  const ResidualNetwork* network_ {&owned_network_};  // the network the searches read: own_, or a captured one (load())
  std::int32_t num_vertices_ {0};
  std::int32_t virtual_source_ {0};
  std::int32_t virtual_sink_ {0};
  std::int32_t search_source_ {0};
  std::int32_t search_target_ {0};
  std::vector<std::int32_t> arc_fill_positions_ {};
  std::vector<Amount> arc_residual_ {};
  // A solve from zero.
  std::vector<Cost> potentials_ {};
  std::vector<std::uint8_t> ever_settled_ {};
  std::vector<Cost> distances_ {};
  RadixHeap<Cost> heap_ {};
  bool shadow_ {false};
  I128 tier_sum_ {};
  std::vector<I128> checkpoint_ {};
  // The deletion mode.
  std::vector<I128> repair_potentials_ {};
  std::vector<I128> repair_distances_ {};
  std::vector<std::uint32_t> repair_seen_generations_ {};
  std::uint32_t repair_generation_ {0};
  RadixHeap<I128> repair_heap_ {};
#if NETGRAPH_CORE_PRIVATE_HOOKS
  const DeletionProbe* probe_ {nullptr};
  const std::function<void(std::span<const Cost>)>* on_search_ {nullptr};
#endif
  // Shared scratch.
  std::vector<std::int32_t> settled_ {};
  std::vector<std::int32_t> level_ {};
  std::vector<std::uint32_t> reach_ {};
  std::uint32_t reachability_generation_ {0};
  std::vector<std::int32_t> queue_ {};
  std::vector<std::int32_t> current_arcs_ {};
  std::vector<std::int32_t> stack_ {};
  std::vector<std::int32_t> path_ {};
  std::vector<std::uint8_t> seen_ {};
};

std::vector<Amount> DeletionState::edge_flows() const {
  std::vector<Amount> flows(static_cast<std::size_t>(graph_->num_edges()), 0);
  for (std::size_t k = 0; k < network_->edges.size(); ++k) {
    const auto fwd = static_cast<std::size_t>(network_->edge_forward_arcs[k]);
    flows[static_cast<std::size_t>(network_->edges[k])] = arc_residual_[static_cast<std::size_t>(network_->reverse_arcs[fwd])];
  }
  return flows;
}

DeletionState& DeletionState::operator=(const DeletionState& other) {
  DeletionState copy(other);  // every allocation; *this is untouched until it succeeds
  return *this = std::move(copy);
}

// The solver reads the network that state_ shares: it lives as long as state_, a move included
// (a move transfers the owner, the network stays where it is).
DeletionResolver::DeletionResolver(const DeletionState& state) : state_(state), solver_(std::make_unique<MinCostFlowSolver>()) {
  if (state_.network_ == nullptr) throw std::invalid_argument("deletion re-solve: an empty state (no solve captured it)");
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
    MinCostFlowSolver& solver;
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
  const ResidualNetwork& net = state_.network();
  const auto it = std::lower_bound(net.edges.begin(), net.edges.end(), e);
  if (it == net.edges.end() || *it != e) {
    throw std::invalid_argument("deletion re-solve: edge " + std::to_string(e) + " is not in the solve's view");
  }
  const auto k = static_cast<std::size_t>(it - net.edges.begin());
  const auto fwd = static_cast<std::size_t>(net.edge_forward_arcs[k]);
  if (state_.residual()[static_cast<std::size_t>(net.reverse_arcs[fwd])] <= 0) {
    throw std::invalid_argument("deletion re-solve: edge " + std::to_string(e) + " carries no flow");
  }
  return k;
}

MinCostFlowResult solve_min_cost_flow(const Graph& graph,
                                  std::span<const Amount> view,
                                  std::span<const NodeId> sources,
                                  std::span<const NodeId> destinations,
                                  const Eligibility* eligibility,
                                  std::optional<Amount> cap,
                                  DeletionCapture* capture,
                                  std::span<const Amount> source_caps) {
  const auto num_edges = static_cast<std::size_t>(graph.num_edges());
  check_edge_amounts(graph, view, "view");
  const auto srcs = node_set(sources, graph.num_nodes(), "sources");
  const auto dsts = node_set(destinations, graph.num_nodes(), "destinations");
  check_disjoint(srcs, dsts, "sources", "destinations");
  if (cap.has_value() && (*cap < 0 || *cap > kMaxFlowCap)) {
    throw std::invalid_argument("cap " + std::to_string(*cap) + " is outside [0, 2^63 - 2]");
  }
  if (eligibility != nullptr && eligibility->graph_uid() != graph.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  if (capture != nullptr && capture->graph.get() != &graph) {
    throw std::logic_error("flow core: a capture request bound to another graph");
  }
  if (!source_caps.empty() && (source_caps.size() != srcs.size() || capture != nullptr)) {
    throw std::logic_error("flow core: source capacities not aligned with the sources, or captured");
  }
  MinCostFlowResult out;
  const bool zero = cap.has_value() && *cap == 0;
  if (zero) {
    // The input cap 0 runs nothing: the prescribed zero payload.
    out.edge_flows.assign(num_edges, 0);
    out.potentials.assign(static_cast<std::size_t>(graph.num_nodes()), kNeverReached);
    if (capture == nullptr) return out;
  }
  // Identity coordinates: every real vertex, the allowed edges and endpoints.
  MinCostFlowSolver& solver = MinCostFlowSolver::local();
  ResidualNetwork& in = solver.network();
  in.num_nodes = graph.num_nodes();
  const auto edge_sources = graph.edge_sources();
  const auto edge_destinations = graph.edge_destinations();
  const auto cost = graph.edge_costs();
  in.edges.clear();
  in.edge_sources.clear();
  in.edge_destinations.clear();
  in.cost.clear();
  in.edge_capacities.clear();
  in.sources.clear();
  in.source_capacities.clear();
  in.destinations.clear();
  in.edges.reserve(num_edges);
  for (std::size_t e = 0; e < num_edges; ++e) {
    if (eligibility != nullptr && !eligibility->allows_edge(static_cast<EdgeId>(e))) continue;
    in.edges.push_back(static_cast<EdgeId>(e));
    in.edge_sources.push_back(edge_sources[e]);
    in.edge_destinations.push_back(edge_destinations[e]);
    in.cost.push_back(cost[e]);
    in.edge_capacities.push_back(view[e]);
  }
  for (std::size_t i = 0; i < srcs.size(); ++i) {
    if (eligibility != nullptr && !eligibility->allows_node(srcs[i])) continue;
    in.sources.push_back(srcs[i]);
    if (!source_caps.empty()) in.source_capacities.push_back(source_caps[i]);
  }
  for (NodeId d : dsts) {
    if (eligibility == nullptr || eligibility->allows_node(d)) in.destinations.push_back(d);
  }
  solver.build(capture != nullptr);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  solver.observe(capture != nullptr ? &capture->on_search : nullptr);
#endif
  if (!zero) {
    auto r = solver.run(cap.value_or(kUncapped), /*closure=*/true);
    out.flow_value = r.flow_value;
    out.marginal_cost_curve = std::move(r.marginal_cost_curve);
    out.edge_flows.assign(num_edges, 0);
    for (std::size_t k = 0; k < in.edges.size(); ++k) {
      out.edge_flows[static_cast<std::size_t>(in.edges[k])] = r.edge_flows[k];
    }
    if (r.min_cut_positions.has_value()) {
      std::vector<EdgeId> cut;
      cut.reserve(r.min_cut_positions->size());
      for (std::int32_t k : *r.min_cut_positions) cut.push_back(in.edges[static_cast<std::size_t>(k)]);
      out.min_cut = std::move(cut);
    }
    out.potentials = std::move(r.potentials);
    out.counters = r.counters;
  }
  if (capture != nullptr) solver.capture(*capture, out.flow_value, out.counters);
  return out;
}

namespace {

// solve_subset's body; `by_cost` false takes every arc's cost as 0 (max_flow_subset).
SubsetFlow subset_flow(const Graph& g, std::span<const EdgeId> edges, std::span<const Amount> caps,
                       std::span<const NodeId> sources, std::span<const NodeId> destinations, std::optional<Amount> cap,
                       std::span<const Amount> source_caps, bool by_cost) {
  SubsetFlow out;
  out.edge_flows.assign(edges.size(), 0);
  if (cap.has_value() && *cap == 0) return out;
  const auto edge_sources = g.edge_sources();
  const auto edge_destinations = g.edge_destinations();
  const auto cost = g.edge_costs();
  // Compressed vertices in first-touch order. The numbering changes no flow, curve or path:
  // the arc order at every vertex is by EdgeId (virtual arc first) and the virtual arcs of s*
  // and t* follow the endpoint lists, which stay in ascending NodeId.
  NodeMap& map = NodeMap::local();
  map.reset(static_cast<std::size_t>(g.num_nodes()));
  // The input's storage is the thread's solver's (every field is rewritten here).
  MinCostFlowSolver& solver = MinCostFlowSolver::local();
  ResidualNetwork& in = solver.network();
  in.edges.assign(edges.begin(), edges.end());
  in.edge_sources.clear();
  in.edge_destinations.clear();
  in.cost.clear();
  in.sources.clear();
  in.destinations.clear();
  for (EdgeId e : edges) {
    const auto ei = static_cast<std::size_t>(e);
    in.edge_sources.push_back(map.add(edge_sources[ei]));
    in.edge_destinations.push_back(map.add(edge_destinations[ei]));
    in.cost.push_back(by_cost ? cost[ei] : 0);
  }
  for (NodeId v : sources) in.sources.push_back(map.add(v));
  in.source_capacities.assign(source_caps.begin(), source_caps.end());
  for (NodeId v : destinations) in.destinations.push_back(map.add(v));
  in.num_nodes = static_cast<std::int32_t>(map.size());
  in.edge_capacities.assign(caps.begin(), caps.end());
  solver.build(false);
  auto r = solver.run(cap.value_or(kUncapped), /*closure=*/false, /*potentials=*/false);
  out.flow_value = r.flow_value;
  out.edge_flows = std::move(r.edge_flows);
  out.marginal_cost_curve = std::move(r.marginal_cost_curve);
  return out;
}

} // namespace

SubsetFlow solve_subset(const Graph& graph, std::span<const EdgeId> edges, std::span<const Amount> caps,
                        std::span<const NodeId> sources, std::span<const NodeId> destinations,
                        std::optional<Amount> cap, std::span<const Amount> source_caps) {
  return subset_flow(graph, edges, caps, sources, destinations, cap, source_caps, true);
}

Amount max_flow_subset(const Graph& g, std::span<const EdgeId> edges, std::span<const Amount> caps,
                       std::span<const NodeId> sources, std::span<const NodeId> destinations, Amount cap,
                       std::span<const Amount> source_caps) {
  return subset_flow(g, edges, caps, sources, destinations, cap, source_caps, false).flow_value;
}

} // namespace netgraph::core::detail
