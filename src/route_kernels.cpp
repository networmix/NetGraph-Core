/*
  Apportionment and the placement kernels in units. T is a plan's total injection.

  Every forwarded Amount is an int64 number of units below 2^62: T is bounded by the view on the roots'
  members (the graph's capacities sum below 2^62), each node's inflow by T, and each share by that
  inflow. Highest averages compares
  w_i/d_i and w_j/d_j through w_i*d_j and w_j*d_i (operands < 2^63, products < 2^126). The
  lower quota floor(I*w/W) uses a 128-bit product, with inflow I, weight w < 2^62 and the exact
  128-bit weight sum W (explicit weights are each below 2^62, not their sum).
  Cumulative source apportionment instead uses totals below 2^127 and 192-bit products.
*/
#include "netgraph/core/detail/route_kernels.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>

#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/concurrency.hpp"
#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/ledger_entry.hpp"
#include "netgraph/core/detail/route_builders.hpp"

namespace netgraph::core::detail {

// ---------------------------------------------------------------- apportionment

namespace {

// w_i / d_i > w_j / d_j, ties by the lower index: the next unit's recipient order.
bool before(std::uint64_t wi, std::uint64_t di, std::size_t i, std::uint64_t wj, std::uint64_t dj, std::size_t j) {
  const U128 left = mul_u64(wi, dj);
  const U128 right = mul_u64(wj, di);
  if (!(left == right)) return left > right;
  return i < j;
}

Amount quota(Amount total, Amount w, U128 sum) noexcept {
  return static_cast<Amount>(div_u128(mul_u64(static_cast<std::uint64_t>(total), static_cast<std::uint64_t>(w)), sum));
}
U128 quota(U128 total, Amount w, U128 sum) noexcept {
  return div_wide(mul_wide(total, static_cast<std::uint64_t>(w)), sum);
}
Amount plus(Amount a, Amount b) noexcept { return a + b; }
U128 plus(U128 a, U128 b) noexcept { return add(a, b); }
Amount remainder(Amount a, Amount b) noexcept { return a - b; }
Amount remainder(U128 a, U128 b) noexcept { return static_cast<Amount>(sub(a, b).lo); }
void award(Amount& a) noexcept { ++a; }
void award(U128& a) noexcept { a = add(a, U128 {0, 1}); }
bool next_before(Amount wi, Amount ai, std::size_t i, Amount wj, Amount aj, std::size_t j) noexcept {
  return before(static_cast<std::uint64_t>(wi), static_cast<std::uint64_t>(ai) + 1, i,
                static_cast<std::uint64_t>(wj), static_cast<std::uint64_t>(aj) + 1, j);
}
bool next_before(Amount wi, U128 ai, std::size_t i, Amount wj, U128 aj, std::size_t j) noexcept {
  const auto left = mul_wide(add(aj, U128 {0, 1}), static_cast<std::uint64_t>(wi));
  const auto right = mul_wide(add(ai, U128 {0, 1}), static_cast<std::uint64_t>(wj));
  return left != right ? left > right : i < j;
}

// One rule for both widths: only quota and priority arithmetic differ. The
// unassigned remainder is less than k, even for a cumulative 127-bit total.
template<class Share>
void apportion(Share total, std::span<const Amount> weights, std::span<Share> out) {
  const std::size_t k = weights.size();
  U128 wsum {};
  for (Amount w : weights) wsum = add(wsum, U128 {0, static_cast<std::uint64_t>(w)});
  if (wsum == U128 {}) {  // all-zero weights give nothing
    std::fill(out.begin(), out.begin() + static_cast<std::ptrdiff_t>(k), Share {});
    return;
  }
  Share given {};
  for (std::size_t m = 0; m < k; ++m) {
    // Lower quota floor(total * w / W): exactly the quotients w / j >= W / total.
    out[m] = quota(total, weights[m], wsum);
    given = plus(given, out[m]);
  }
  // Fewer than k (the positive weights' count) units are left; a zero weight's quotient 0 never
  // beats a positive one, so it gets none of them.
  Amount left = remainder(total, given);
  if (left == 0) return;
  if (k <= 16) {
    while (left-- > 0) {
      std::size_t best = 0;
      for (std::size_t m = 1; m < k; ++m) {
        if (next_before(weights[m], out[m], m, weights[best], out[best], best)) {
          best = m;
        }
      }
      award(out[best]);
    }
    return;
  }
  std::vector<std::size_t> heap(k);
  std::iota(heap.begin(), heap.end(), std::size_t{0});
  // std heaps put the largest first: "a < b" when b goes first.
  auto later = [&](std::size_t a, std::size_t b) {
    return next_before(weights[b], out[b], b, weights[a], out[a], a);
  };
  std::make_heap(heap.begin(), heap.end(), later);
  while (left-- > 0) {
    std::pop_heap(heap.begin(), heap.end(), later);
    award(out[heap.back()]);
    std::push_heap(heap.begin(), heap.end(), later);
  }
}

} // namespace

void highest_averages(Amount total, std::span<const Amount> weights, std::span<Amount> out) {
  apportion(total, weights, out);
}

void highest_averages(U128 total, std::span<const Amount> weights, std::span<U128> out) {
  apportion(total, weights, out);
}

void round_robin_increments(const Counter128& total, Amount v, std::span<Amount> out) noexcept {
  const auto k = static_cast<std::uint64_t>(out.size());
  if (k == 0) return;
  const std::uint64_t r = mod_u128(U128 {total.hi, total.lo}, k);
  const std::uint64_t t = static_cast<std::uint64_t>(v) + r;
  const auto base = static_cast<Amount>(t / k);
  const std::uint64_t extra = t % k;
  for (std::size_t i = 0; i < out.size(); ++i) out[i] = base + (i < extra ? 1 : 0) - (i < r ? 1 : 0);
}

// ---------------------------------------------------------------- plans

void check_kernel_combination(RouteAllocation allocation, Admission admission) {
  if (admission == Admission::Lossy && allocation == RouteAllocation::Optimized) {
    throw std::invalid_argument("OPTIMIZED requires admission=LOSSLESS");
  }
}

void check_kernel(RouteAllocation allocation, Admission admission) {
  check_member(allocation, "allocation");
  check_member(admission, "admission");
  check_kernel_combination(allocation, admission);
}

PreparedRoute::Layout PreparedRoute::point(std::size_t num_roots, std::size_t num_dsts, std::size_t positions,
                                     std::size_t num_members) noexcept {
  std::int32_t* at = store_.data();
  Layout w;
  sources = std::span<const NodeId>(at, num_roots);
  at += num_roots;
  w.source_positions = std::span<std::int32_t>(at, num_roots);
  source_positions = w.source_positions;
  at += num_roots;
  destinations = std::span<const NodeId>(at, num_dsts);
  at += num_dsts;
  w.topological_order = std::span<NodeId>(at, positions);
  at += positions;
  w.outgoing_begin = std::span<std::int32_t>(at, positions);
  at += positions;
  w.outgoing_end = std::span<std::int32_t>(at, positions);
  at += positions;
  w.edge_destination_positions = std::span<std::int32_t>(at, num_members);
  at += num_members;
  w.is_destination = std::span<std::uint8_t>(reinterpret_cast<std::uint8_t*>(at), positions);
  topological_order = w.topological_order;
  outgoing_begin = w.outgoing_begin;
  outgoing_end = w.outgoing_end;
  edge_destination_positions = w.edge_destination_positions;
  is_destination = w.is_destination;
  return w;
}

PreparedRoute::Layout PreparedRoute::layout(std::span<const NodeId> roots_in, std::span<const NodeId> destinations,
                                      std::size_t positions, std::size_t num_members) {
  const std::size_t words =
      2 * roots_in.size() + destinations.size() + 3 * positions + num_members + (positions + 3) / 4;
  std::vector<std::int32_t> store(words, 0);  // the one allocation
  store_.swap(store);
  const Layout w = point(roots_in.size(), destinations.size(), positions, num_members);
  std::copy(roots_in.begin(), roots_in.end(), store_.begin());
  std::fill(w.source_positions.begin(), w.source_positions.end(), -1);
  std::copy(destinations.begin(), destinations.end(), store_.begin() + static_cast<std::ptrdiff_t>(2 * roots_in.size()));
  return w;
}

PreparedRoute::PreparedRoute(const PreparedRoute& other)
    : allocation(other.allocation), support_edges(other.support_edges), forwarding_edges(other.forwarding_edges),
      weights(other.weights), store_(other.store_) {
  (void)point(other.sources.size(), other.destinations.size(), other.topological_order.size(), other.edge_destination_positions.size());
}

PreparedRoute::PreparedRoute(PreparedRoute&& other) noexcept
    : allocation(other.allocation), support_edges(std::move(other.support_edges)),
      forwarding_edges(std::move(other.forwarding_edges)), weights(std::move(other.weights)), store_(std::move(other.store_)) {
  (void)point(other.sources.size(), other.destinations.size(), other.topological_order.size(), other.edge_destination_positions.size());
  (void)other.point(0, 0, 0, 0);
}

PreparedRoute& PreparedRoute::operator=(const PreparedRoute& other) {
  if (this != &other) {
    PreparedRoute copy(other);
    *this = std::move(copy);
  }
  return *this;
}

PreparedRoute& PreparedRoute::operator=(PreparedRoute&& other) noexcept {
  if (this != &other) {
    allocation = other.allocation;
    support_edges = std::move(other.support_edges);
    forwarding_edges = std::move(other.forwarding_edges);
    weights = std::move(other.weights);
    const std::size_t nr = other.sources.size();
    const std::size_t nd = other.destinations.size();
    const std::size_t np = other.topological_order.size();
    const std::size_t nm = other.edge_destination_positions.size();
    store_ = std::move(other.store_);
    (void)point(nr, nd, np, nm);
    (void)other.point(0, 0, 0, 0);
  }
  return *this;
}

namespace {

// Per-thread build_split state, by NodeId, for its backward depth-first search (DFS) and
// forward pass. A call advances an epoch (a numeric stamp); a slot is valid only if its stamp
// matches. Later calls ignore stale entries, including after exceptions; only stamp wraparound
// (after 2^32 calls on one thread, which a long-running process reaches) clears all stamps. One
// array holds every node's fields, so a failed growth leaves the old array intact (the strong
// guarantee); separate arrays could end at mismatched sizes.
struct DagScratch {
  struct Node {
    std::uint32_t stamp {0};
    std::uint8_t state {0};   // backward DFS: 1 open, 2 done (reaches a destination)
    std::uint8_t fwd {0};     // reached from the source
    std::uint8_t target {0};  // a destination
    std::uint8_t root {0};    // a root of the plan
    std::int32_t position {0};
    std::int32_t first {0};   // the node's first member
    std::int32_t count {0};   // the node's number of members
    std::int32_t entry_begin {0};  // the node's predecessor entries in the DAG's sparse storage
    std::int32_t entry_end {0};
    Cost dist {0};
  };
  std::uint32_t generation {0};
  std::vector<Node> node;
  std::vector<NodeId> order;
  std::vector<std::pair<NodeId, std::int32_t>> stack;
  std::vector<NodeId> member_tail;
  std::vector<NodeId> member_head;

  void begin(std::size_t n) {
    if (node.size() < n) node.resize(n);  // strong guarantee; new slots carry stamp 0
    if (++generation == 0) {
      for (auto& v : node) v.stamp = 0;
      generation = 1;
    }
  }
  [[nodiscard]] bool seen(NodeId v) const noexcept { return node[static_cast<std::size_t>(v)].stamp == generation; }
  Node& at(NodeId v) noexcept { return node[static_cast<std::size_t>(v)]; }
  void touch(NodeId v) noexcept {
    auto& x = at(v);
    x = Node {};
    x.stamp = generation;
  }
};

DagScratch& dag_scratch() {
  thread_local DagScratch s;
  return s;
}

// Build plan.forwarding_edges from the DAG's incoming-edge lists: keep edges accepted by keep(e) that
// lie on a walk from a root to a destination; exclude edges leaving destinations. No node
// renumbering, extra arc arrays or adjacency offsets are needed. A backward DFS from destinations
// lists each node after its parents (postorder, also a topological order here).
// A forward pass then marks nodes reachable from the roots. Members ascend by EdgeId, so one
// tail's members are consecutive. Any topological order gives the same allocation because each
// node splits only after receiving all its inflow.
// With `cost`, the forward pass computes minimum cost to each node from its parents' costs
// (dynamic programming, DP), giving the minimum root-to-destination cost over members, or
// kUnreachable if none exist. Without `out` (min_route_cost), stop there and build no plan.
template <class Keep>
void build_split(const Graph& graph, const PredecessorDAG& dag, std::span<const NodeId> roots, std::span<const NodeId> destinations,
                 Keep keep, PreparedRoute* out, Cost* cost) {
  concurrent_event(ConcurrentStage::RouteWalk, &dag);
  if (cost != nullptr) *cost = kUnreachable;
  const auto parents = dag.parents();
  const auto via = dag.via_edges();
  DagScratch& s = dag_scratch();
  s.begin(static_cast<std::size_t>(graph.num_nodes()));
  // The DAG's sparse storage scattered into the per-node slots (O(member nodes)); a node outside
  // it has no entries when first touched (a destination, source or a parent the search meets).
  {
    const auto nodes = PredecessorDAGBuilder::member_nodes(dag);
    const auto offs = PredecessorDAGBuilder::member_offsets(dag);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      s.touch(nodes[i]);
      s.at(nodes[i]).entry_begin = offs[i];
      s.at(nodes[i]).entry_end = offs[i + 1];
    }
  }
  for (NodeId d : destinations) {
    if (!s.seen(d)) s.touch(d);
    s.at(d).target = 1;
  }
  for (NodeId r : roots) {
    if (!s.seen(r)) s.touch(r);
    s.at(r).root = 1;
  }
  auto is_dst = [&](NodeId u) { return s.seen(u) && s.at(u).target != 0; };
  // Backward DFS from the destinations, postorder (every node after all of its parents), over
  // the entries (u, e) that `keep` accepts and whose parent is not a destination.
  s.order.clear();
  for (NodeId d : destinations) {
    if (s.at(d).state != 0) continue;
    s.at(d).state = 1;
    s.stack.clear();
    s.stack.emplace_back(d, s.at(d).entry_begin);
    while (!s.stack.empty()) {
      const auto [v, k] = s.stack.back();
      if (k < s.at(v).entry_end) {
        s.stack.back().second = k + 1;
        const auto ki = static_cast<std::size_t>(k);
        const NodeId u = parents[ki];
        if (is_dst(u) || !keep(via[ki])) continue;
        if (!s.seen(u)) s.touch(u);
        const auto su = s.at(u).state;
        if (su == 0) {
          s.at(u).state = 1;
          s.stack.emplace_back(u, s.at(u).entry_begin);
        } else if (su == 1) {
          throw std::logic_error("kernel plan: the route is cyclic");
        }
      } else {
        s.at(v).state = 2;
        s.order.push_back(v);
        s.stack.pop_back();
      }
    }
  }
  // Forward over the postorder, with the cost DP when asked.
  const auto edge_cost = graph.edge_costs();
  Cost best = kUnreachable;
  for (NodeId r : roots) {
    if (s.at(r).state == 2) {
      s.at(r).fwd = 1;
      s.at(r).dist = 0;
    }
  }
  for (NodeId v : s.order) {
    auto& nv = s.at(v);
    if (nv.root == 0) {
      Cost dv = kUnreachable;
      bool reached = false;
      for (auto k = nv.entry_begin; k < nv.entry_end; ++k) {
        const auto ki = static_cast<std::size_t>(k);
        const NodeId u = parents[ki];
        if (!s.seen(u)) continue;
        const auto& nu = s.at(u);
        if (nu.fwd == 0 || nu.target != 0 || !keep(via[ki])) continue;
        reached = true;
        if (cost == nullptr) break;
        dv = std::min(dv, nu.dist + edge_cost[static_cast<std::size_t>(via[ki])]);
      }
      if (!reached) continue;
      nv.fwd = 1;
      nv.dist = dv;
    }
    if (cost != nullptr && nv.fwd != 0 && nv.target != 0) best = std::min(best, nv.dist);
  }
  if (cost != nullptr) *cost = best;
  if (out == nullptr) return;
  PreparedRoute& plan = *out;
  plan.forwarding_edges.clear();
  // Members, ascending: the DAG's edges kept, from a reached tail (not a destination) to a head
  // that reaches a destination.
  const auto edge_sources = graph.edge_sources();
  const auto edge_destinations = graph.edge_destinations();
  // Size vectors to known bounds (DAG edges here, reached nodes below), then write by index
  // to avoid push_back overhead in these loops.
  const auto asc = PredecessorDAGBuilder::ascending(dag);
  const std::size_t bound = asc.size();
  plan.forwarding_edges.resize(bound);
  s.member_tail.resize(bound);
  s.member_head.resize(bound);
  std::size_t n_members = 0;
  for (EdgeId e : asc) {
    const auto ei = static_cast<std::size_t>(e);
    const NodeId u = edge_sources[ei];
    const NodeId v = edge_destinations[ei];
    if (!s.seen(u) || s.at(u).fwd == 0 || s.at(u).target != 0) continue;
    if (!s.seen(v) || s.at(v).state != 2 || !keep(e)) continue;
    plan.forwarding_edges[n_members] = e;
    s.member_tail[n_members] = u;
    s.member_head[n_members] = v;
    ++n_members;
  }
  plan.forwarding_edges.resize(n_members);
  s.member_tail.resize(n_members);
  s.member_head.resize(n_members);
  const std::size_t num_members = n_members;
  if (num_members == 0) {
    (void)plan.layout(roots, destinations, 0, 0);  // every root without a member
    return;
  }
  std::int32_t num_pos = 0;
  for (NodeId v : s.order) {
    if (s.at(v).fwd != 0) s.at(v).position = num_pos++;
  }
  for (std::size_t m = 0; m < num_members; ++m) {
    auto& nt = s.at(s.member_tail[m]);
    if (nt.count == 0) {
      nt.first = static_cast<std::int32_t>(m);
    } else if (static_cast<std::size_t>(nt.first + nt.count) != m) {
      throw std::logic_error("kernel plan: the members of a node are not consecutive (member EdgeIds not grouped by source node)");
    }
    ++nt.count;
  }
  const auto P = static_cast<std::size_t>(num_pos);
  const PreparedRoute::Layout w = plan.layout(roots, destinations, P, num_members);
  std::size_t pos = 0;
  for (NodeId v : s.order) {
    const auto& nv = s.at(v);
    if (nv.fwd == 0) continue;
    w.topological_order[pos] = v;
    w.is_destination[pos] = nv.target;
    const std::int32_t b = nv.count == 0 ? 0 : nv.first;
    w.outgoing_begin[pos] = b;
    w.outgoing_end[pos] = b + nv.count;
    ++pos;
  }
  // Positions follow the postorder (every node after its parents; the DFS rejected cycles), so a
  // member's tail precedes its head.
  for (std::size_t m = 0; m < num_members; ++m) w.edge_destination_positions[m] = s.at(s.member_head[m]).position;
  // A root with a member has a position (it reaches a destination and is reached); a root
  // without one keeps -1.
  for (std::size_t i = 0; i < roots.size(); ++i) {
    const auto& nr = s.at(roots[i]);
    if (nr.fwd != 0 && nr.target == 0) w.source_positions[i] = nr.position;
  }
}

// Distances to the first destination, shared by all roots. Outgoing DFS finishes each child
// before its parent: d(u) = min(cost(u,v) + d(v)). With unit costs it gives the fewest hops.
// The same stamped slots serve arbitrary DAG costs and canonical ShortestPaths hop selection;
// no graph-sized buffer is cleared on each call. A failed call's slots expire at begin().
DagScratch& route_distances(const Graph& graph, std::span<const EdgeId> edges, std::span<const NodeId> roots,
                            std::span<const NodeId> destinations, bool unit_cost) {
  DagScratch& s = dag_scratch();
  s.begin(static_cast<std::size_t>(graph.num_nodes()));
  s.stack.clear();
  auto touch = [&](NodeId v) -> DagScratch::Node& {
    if (!s.seen(v)) {
      s.touch(v);
      s.at(v).dist = kUnreachable;
    }
    return s.at(v);
  };
  const auto tails = graph.edge_sources();
  const auto heads = graph.edge_destinations();
  const auto costs = graph.edge_costs();
  for (std::size_t lo = 0; lo < edges.size();) {
    const NodeId u = tails[static_cast<std::size_t>(edges[lo])];
    std::size_t hi = lo + 1;
    while (hi < edges.size() && tails[static_cast<std::size_t>(edges[hi])] == u) ++hi;
    auto& nu = touch(u);
    nu.first = static_cast<std::int32_t>(lo);
    nu.count = static_cast<std::int32_t>(hi - lo);
    lo = hi;
  }
  for (NodeId d : destinations) {
    auto& nd = touch(d);
    nd.target = 1;
    nd.state = 2;
    nd.dist = 0;
  }
  for (NodeId root : roots) {
    auto& nr = touch(root);
    if (nr.state == 2) continue;
    nr.state = 1;
    s.stack.emplace_back(root, nr.first);
    while (!s.stack.empty()) {
      const auto [u, at] = s.stack.back();
      auto& nu = s.at(u);
      if (at == nu.first + nu.count) {
        nu.state = 2;
        s.stack.pop_back();
        continue;
      }
      const EdgeId e = edges[static_cast<std::size_t>(at)];
      const NodeId v = heads[static_cast<std::size_t>(e)];
      auto& nv = touch(v);
      if (nv.state == 0) {
        nv.state = 1;
        s.stack.emplace_back(v, nv.first);
        continue;  // read this edge after the child's distance is complete
      }
      if (nv.state == 1) throw std::logic_error("kernel plan: the route is cyclic");
      if (nv.dist != kUnreachable) {
        nu.dist = std::min(nu.dist, nv.dist + (unit_cost ? Cost {1} : costs[static_cast<std::size_t>(e)]));
      }
      ++s.stack.back().second;
    }
  }
  return s;
}

} // namespace

std::vector<EdgeId> prune(const Graph& graph, const PredecessorDAG& dag, NodeId source, std::span<const NodeId> destinations,
                          const Eligibility* eligibility) {
  PreparedRoute plan;
  build_split(graph, dag, {&source, 1}, destinations, [&](EdgeId e) { return eligibility == nullptr || eligibility->allows_edge(e); }, &plan, nullptr);
  return std::move(plan.forwarding_edges);
}

std::vector<EdgeId> prune(const Graph& graph, const PredecessorDAG& dag, std::span<const NodeId> roots,
                          std::span<const NodeId> destinations, const Eligibility* eligibility) {
  PreparedRoute plan;
  build_split(graph, dag, roots, destinations, [&](EdgeId e) { return eligibility == nullptr || eligibility->allows_edge(e); }, &plan, nullptr);
  return std::move(plan.forwarding_edges);
}

Cost min_route_cost(const Graph& graph, const PredecessorDAG& route, NodeId source, std::span<const NodeId> destinations) {
  Cost cost = kUnreachable;
  build_split(graph, route, {&source, 1}, destinations, [](EdgeId) { return true; }, nullptr, &cost);
  return cost;
}

std::vector<Cost> min_route_costs(const Graph& graph, const PredecessorDAG& route, std::span<const NodeId> roots,
                                std::span<const NodeId> destinations) {
  if (roots.empty()) return {};
  const auto& s = route_distances(graph, PredecessorDAGBuilder::ascending(route), roots, destinations, false);
  std::vector<Cost> out;
  out.reserve(roots.size());
  for (NodeId root : roots) out.push_back(s.node[static_cast<std::size_t>(root)].dist);
  return out;
}

namespace {

PreparedRoute plan_of(const Graph& graph, RouteAllocation allocation, Admission admission, const WeightSource& weights,
                   std::span<const NodeId> roots, std::span<const NodeId> destinations, const PredecessorDAG& route, Cost* cost) {
  check_kernel(allocation, admission);
  const auto edges = PredecessorDAGBuilder::ascending(route);
  const bool weighted = allocation == RouteAllocation::Weighted;
  const bool by_view = weighted && weights.basis == WeightBasis::ResidualSnapshot;
  if (by_view && weights.values.size() != edges.size()) {
    throw std::logic_error("prepare_route: the creation view is not route-aligned");
  }
  if (weighted && weights.basis == WeightBasis::Explicit &&
      weights.values.size() != static_cast<std::size_t>(graph.num_edges())) {
    throw std::logic_error("prepare_route: the explicit weights are not edge-aligned");
  }
  const auto capacity = graph.edge_capacity_units();
  PreparedRoute plan;
  plan.allocation = allocation;
  plan.support_edges.assign(edges.begin(), edges.end());
  // A route edge's Weighted weight: its static capacity, its explicit weight (0 on an edge of
  // capacity 0) or its creation view (the route ascends: a binary search).
  auto weight = [&](EdgeId e) -> Amount {
    const auto ei = static_cast<std::size_t>(e);
    switch (weights.basis) {
      case WeightBasis::Capacity: return capacity[ei];
      case WeightBasis::Explicit: return capacity[ei] >= 1 ? weights.values[ei] : 0;
      case WeightBasis::ResidualSnapshot:
        return weights.values[static_cast<std::size_t>(std::lower_bound(edges.begin(), edges.end(), e) - edges.begin())];
    }
    return 0;
  };
  // Membership: Optimized and Equal every route edge (an edge of capacity 0 is a member
  // that carries nothing), Weighted the edges of positive weight.
  auto member = [&](EdgeId e) { return !weighted || weight(e) >= 1; };
  build_split(graph, route, roots, destinations, member, &plan, cost);
  // The route's cost is the members' only when every route edge is a member.
  if (cost != nullptr && plan.forwarding_edges.size() != plan.support_edges.size()) *cost = min_route_cost(graph, route, roots[0], destinations);
  if (weighted) {
    plan.weights.reserve(plan.forwarding_edges.size());
    for (EdgeId e : plan.forwarding_edges) plan.weights.push_back(weight(e));
  }
  return plan;
}

} // namespace

PreparedRoute prepare_route(const Graph& graph, RouteAllocation allocation, Admission admission, const WeightSource& weights,
                     NodeId source, std::span<const NodeId> destinations, const PredecessorDAG& route, Cost* cost) {
  return plan_of(graph, allocation, admission, weights, {&source, 1}, destinations, route, cost);
}

PreparedRoute prepare_route(const Graph& graph, RouteAllocation allocation, Admission admission, const WeightSource& weights,
                     std::span<const NodeId> roots, std::span<const NodeId> destinations, const PredecessorDAG& route) {
  return plan_of(graph, allocation, admission, weights, roots, destinations, route, nullptr);
}

PreparedRoute prepare_supplied_route(const Graph& graph, RouteAllocation allocation, WeightBasis basis, std::span<const Amount> edge_values,
                           NodeId source, std::span<const NodeId> destinations, const PredecessorDAG& dag, Cost* cost) {
  check_kernel(allocation, Admission::Lossless);
  const auto capacity = graph.edge_capacity_units();
  PreparedRoute plan;
  plan.allocation = allocation;
  build_split(graph, dag, {&source, 1}, destinations, [](EdgeId) { return true; }, &plan, cost);
  plan.support_edges = plan.forwarding_edges;
  if (allocation != RouteAllocation::Weighted) return plan;
  auto weight = [&](EdgeId e) -> Amount {
    const auto ei = static_cast<std::size_t>(e);
    switch (basis) {
      case WeightBasis::Capacity: return capacity[ei];
      case WeightBasis::Explicit: return capacity[ei] >= 1 ? edge_values[ei] : 0;
      case WeightBasis::ResidualSnapshot: return edge_values[ei];
    }
    return 0;
  };
  // When every route edge has a positive weight the members are the route itself; otherwise a
  // second pass keeps the walks over the positive ones (the route and its cost stay).
  if (std::all_of(plan.support_edges.begin(), plan.support_edges.end(), [&](EdgeId e) { return weight(e) >= 1; })) {
    plan.weights.reserve(plan.forwarding_edges.size());
    for (EdgeId e : plan.forwarding_edges) plan.weights.push_back(weight(e));
    return plan;
  }
  PreparedRoute weighted;
  build_split(graph, dag, {&source, 1}, destinations, [&](EdgeId e) { return weight(e) >= 1; }, &weighted, nullptr);
  weighted.allocation = allocation;
  weighted.support_edges = std::move(plan.support_edges);
  weighted.weights.reserve(weighted.forwarding_edges.size());
  for (EdgeId e : weighted.forwarding_edges) weighted.weights.push_back(weight(e));
  return weighted;
}

PreparedRoute restrict_prepared_route(const Graph& graph, const PreparedRoute& plan, std::span<const NodeId> roots,
                         std::span<const NodeId> destinations, const Eligibility* eligibility) {
  PreparedRoute out;
  out.allocation = plan.allocation;
  const PredecessorDAG members = PredecessorDAGBuilder::of_edges(graph, plan.forwarding_edges);
  build_split(graph, members, roots, destinations, [&](EdgeId e) { return eligibility == nullptr || eligibility->allows_edge(e); }, &out, nullptr);
  out.support_edges = plan.support_edges;
  if (plan.allocation == RouteAllocation::Weighted) {
    // The kept members' frozen weights (both member lists ascend).
    out.weights.resize(out.forwarding_edges.size());
    std::size_t k = 0;
    for (std::size_t m = 0; m < out.forwarding_edges.size(); ++m) {
      while (plan.forwarding_edges[k] != out.forwarding_edges[m]) ++k;
      out.weights[m] = plan.weights[k];
    }
  }
  return out;
}

// ---------------------------------------------------------------- ShortestPaths creation

namespace {

// Per-thread state of the roster route walk. Per node (by NodeId) a slot valid only under the
// current epoch (DagScratch's pattern: nothing is cleared between calls except on stamp
// wraparound, a call that throws leaves stale slots behind an old epoch, and the one array of
// slots grows whole or not at all). The reached list is the queue, then sorted tails; each
// touched node's outgoing range is rewritten before use.
struct CreationScratch {
  struct Node {
    std::uint32_t stamp {0};
    std::uint8_t reached {0};
    std::uint8_t target {0};
    std::int32_t lo {0};     // the node's support out-edges: support[lo, hi)
    std::int32_t hi {0};
  };
  std::uint32_t generation {0};
  std::vector<Node> node;
  std::vector<NodeId> reached;  // the walk's nodes (its queue), then ascending

  void begin(std::size_t n) {
    if (node.size() < n) node.resize(n);  // strong guarantee; new slots carry stamp 0
    if (++generation == 0) {
      for (auto& v : node) v.stamp = 0;
      generation = 1;
    }
  }
  Node& at(NodeId v) noexcept { return node[static_cast<std::size_t>(v)]; }
  // The slot of v, reset on its first use in this call.
  Node& touch(NodeId v) noexcept {
    auto& x = at(v);
    if (x.stamp != generation) {
      x = Node {};
      x.stamp = generation;
    }
    return x;
  }
};

CreationScratch& creation_scratch() {
  thread_local CreationScratch s;
  return s;
}

} // namespace

std::size_t lower_bound_from(std::span<const EdgeId> a, std::size_t hint, EdgeId x) noexcept {
  const std::size_t n = a.size();
  std::size_t lo = 0;
  std::size_t hi = 0;
  if (hint < n && a[hint] < x) {
    // The answer lies in (hint, n]: a[lo - 1] < x throughout.
    lo = hint + 1;
    hi = lo;
    for (std::size_t step = 1; hi < n && a[hi] < x; step <<= 1) {
      lo = hi + 1;
      hi = lo + step;
    }
    hi = std::min(hi, n);
  } else {
    // The answer lies in [0, min(hint, n)] (hint >= n, or a[hint] >= x): a[hi] >= x or hi == n
    // throughout.
    hi = std::min(hint, n);
    for (std::size_t step = 1;; step <<= 1) {
      if (hi < step) {
        lo = 0;
        break;
      }
      const std::size_t probe = hi - step;
      if (a[probe] < x) {
        lo = probe + 1;
        break;
      }
      hi = probe;
    }
  }
  return static_cast<std::size_t>(std::lower_bound(a.begin() + static_cast<std::ptrdiff_t>(lo),
                                                   a.begin() + static_cast<std::ptrdiff_t>(hi), x) -
                                  a.begin());
}

namespace {

std::vector<EdgeId> select_routes(const Graph& graph, std::span<const EdgeId> sup,
                                       std::span<const NodeId> roots, std::span<const NodeId> destinations,
                                       const RouteChoice& choice, PriorityReader view, std::span<const std::uint64_t> allowed = {}) {
  if (roots.empty() || destinations.empty()) return {};
  if (choice.single_path() && choice.prefers_score()) {
    if (roots.size() != 1) throw std::logic_error("ShortestPaths bottleneck route requires one root");
    return support_route(graph, sup, roots.front(), destinations, choice, view).edges;
  }
  const auto rows = graph.outgoing_offsets();
  const auto heads = graph.edge_destinations();
  CreationScratch& s = creation_scratch();
  s.begin(static_cast<std::size_t>(graph.num_nodes()));
  s.reached.clear();
  for (NodeId d : destinations) s.touch(d).target = 1;
  for (NodeId root : roots) {
    s.touch(root).reached = 1;
    s.reached.push_back(root);
  }
  std::size_t hint = 0;
  for (std::size_t q = 0; q < s.reached.size(); ++q) {
    const NodeId u = s.reached[q];
    auto& nu = s.at(u);
    if (nu.target != 0) continue;
    const auto ui = static_cast<std::size_t>(u);
    const std::size_t lo = lower_bound_from(sup, hint, rows[ui]);
    std::size_t hi = lo;
    for (; hi < sup.size() && sup[hi] < rows[ui + 1]; ++hi) {
      if (!allowed.empty() && !bit_at(allowed, hi)) continue;
      const NodeId v = heads[static_cast<std::size_t>(sup[hi])];
      auto& nv = s.touch(v);
      if (nv.reached != 0) continue;
      nv.reached = 1;
      s.reached.push_back(v);
    }
    nu.lo = static_cast<std::int32_t>(lo);
    nu.hi = static_cast<std::int32_t>(hi);
    hint = lo;
  }
  std::sort(s.reached.begin(), s.reached.end());
  std::vector<EdgeId> out;
  for (NodeId u : s.reached) {
    const auto& nu = s.at(u);
    for (auto at = nu.lo; at < nu.hi;) {
      if (!allowed.empty() && !bit_at(allowed, static_cast<std::size_t>(at))) { ++at; continue; }
      EdgeId best = sup[static_cast<std::size_t>(at)];
      auto next = at + 1;
      if (!choice.parallel_edges) {
        const NodeId v = heads[static_cast<std::size_t>(best)];
        for (; next < nu.hi && heads[static_cast<std::size_t>(sup[static_cast<std::size_t>(next)])] == v; ++next) {
          if ((allowed.empty() || bit_at(allowed, static_cast<std::size_t>(next))) && replaces(graph, sup[static_cast<std::size_t>(next)], best, choice.prefers_score(), view)) {
            best = sup[static_cast<std::size_t>(next)];
          }
        }
      }
      out.push_back(best);
      at = next;
    }
  }
  if (!choice.single_path()) return out;
  auto& distances = route_distances(graph, out, roots, destinations, true);
  for (NodeId u : s.reached) distances.at(u).position = -1;
  const auto tails = graph.edge_sources();
  for (EdgeId e : out) {
    const NodeId u = tails[static_cast<std::size_t>(e)];
    const NodeId v = heads[static_cast<std::size_t>(e)];
    auto& nu = distances.at(u);
    const auto& nv = distances.at(v);
    if (nu.position < 0 && nv.dist != kUnreachable && nu.dist == nv.dist + 1) nu.position = e;
  }
  // The tree may leave a reached support branch unused: walk next() from the roots again.
  for (NodeId u : s.reached) s.at(u).reached = 0;
  s.reached.clear();
  for (NodeId root : roots) {
    s.at(root).reached = 1;
    s.reached.push_back(root);
  }
  for (std::size_t q = 0; q < s.reached.size(); ++q) {
    const auto next = distances.at(s.reached[q]).position;
    if (next < 0) continue;
    const NodeId v = heads[static_cast<std::size_t>(next)];
    if (s.at(v).reached != 0) continue;
    s.at(v).reached = 1;
    s.reached.push_back(v);
  }
  std::sort(s.reached.begin(), s.reached.end());
  out.clear();
  for (NodeId u : s.reached) {
    const auto next = distances.at(u).position;
    if (next >= 0) out.push_back(next);
  }
  return out;
}

} // namespace

std::vector<EdgeId> shortest_path_routes(const Graph& graph, const PredecessorDAG& support,
                                       std::span<const NodeId> roots, std::span<const NodeId> destinations,
                                       const RouteChoice& choice, PriorityReader view) {
  return select_routes(graph, PredecessorDAGBuilder::ascending(support), roots, destinations, choice, view);
}

std::vector<EdgeId> reachable_route(const Graph& graph, std::span<const EdgeId> basis,
                                 NodeId source, std::span<const NodeId> destinations, std::span<const std::uint64_t> allowed) {
  return select_routes(graph, basis, {&source, 1}, destinations, RouteChoice {}, PriorityReader {}, allowed);
}

// ---------------------------------------------------------------- allocation

namespace {

// The view summed over the roots' out-members, and at most the number of roots before the first
// root without a member: no larger total fits (that root would receive a share with no route).
Amount total_bound(const PreparedRoute& plan, std::span<const Amount> view_m) {
  Amount hi = 0;
  for (std::size_t i = 0; i < plan.sources.size(); ++i) {
    const auto p = plan.source_positions[i];
    if (p < 0) return std::min(hi, static_cast<Amount>(i));
    for (auto m = plan.outgoing_begin[static_cast<std::size_t>(p)]; m < plan.outgoing_end[static_cast<std::size_t>(p)]; ++m) {
      hi += view_m[static_cast<std::size_t>(m)];
    }
  }
  return hi;
}

} // namespace

bool allocate(const PreparedRoute& plan, Amount total, std::span<Amount> out, std::span<const Amount> limit) {
  std::fill(out.begin(), out.end(), 0);
  if (total == 0) return true;
  const bool weighted = plan.allocation == RouteAllocation::Weighted;
  const bool limited = !limit.empty();
  thread_local std::vector<Amount> inflow;
  inflow.assign(plan.topological_order.size(), 0);
  // The roots' round-robin shares (only the positive ones are visited). A positive share of a
  // root without a member has no route; every caller stays at or below total_bound, which no such
  // share reaches.
  bool routed = true;
  round_robin(total, plan.sources.size(), [&](std::size_t i, Amount x) {
    const auto p = plan.source_positions[i];
    if (p < 0) {
      routed = false;
      return false;
    }
    inflow[static_cast<std::size_t>(p)] += x;
    return true;
  });
  if (!routed) throw std::logic_error("allocate: a root without a member receives a share");
  for (std::size_t p = 0; p < plan.topological_order.size(); ++p) {
    const Amount amount = inflow[p];
    if (amount == 0 || plan.is_destination[p] != 0) continue;
    const auto b = static_cast<std::size_t>(plan.outgoing_begin[p]);
    const auto e = static_cast<std::size_t>(plan.outgoing_end[p]);
    const std::size_t k = e - b;
    if (weighted) {
      highest_averages(amount, std::span<const Amount>(plan.weights.data() + b, k), out.subspan(b, k));
      for (std::size_t m = b; m < e; ++m) {
        const Amount x = out[m];
        if (limited && x > limit[m]) return false;
        inflow[static_cast<std::size_t>(plan.edge_destination_positions[m])] += x;
      }
    } else {
      // Only the members with a positive share are written and visited (the rest stay 0).
      bool fits = true;
      round_robin(amount, k, [&](std::size_t i, Amount x) {
        const std::size_t m = b + i;
        if (limited && x > limit[m]) {
          fits = false;
          return false;
        }
        out[m] = x;
        inflow[static_cast<std::size_t>(plan.edge_destination_positions[m])] += x;
        return true;
      });
      if (!fits) return false;
    }
  }
  return true;
}

namespace {

// A floating-point estimate of T_max: the fraction f of the total that each member carries
// under the plan's split (1/H at each of H roots, 1/k per round-robin member, w/W per Weighted
// member), then min over members of view / f. It only chooses the first probe of the exact
// search below (approximate arithmetic may choose probes, nothing else).
double estimate_max_lossless_offer(const PreparedRoute& plan, std::span<const Amount> view_m) {
  thread_local std::vector<double> f;
  f.assign(plan.topological_order.size(), 0.0);
  for (auto p : plan.source_positions) {
    if (p >= 0) f[static_cast<std::size_t>(p)] += 1.0 / static_cast<double>(plan.sources.size());
  }
  const bool weighted = plan.allocation == RouteAllocation::Weighted;
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t p = 0; p < plan.topological_order.size(); ++p) {
    if (f[p] == 0.0 || plan.is_destination[p] != 0) continue;
    const auto b = static_cast<std::size_t>(plan.outgoing_begin[p]);
    const auto e = static_cast<std::size_t>(plan.outgoing_end[p]);
    double wsum = 0.0;
    if (weighted) {
      for (std::size_t m = b; m < e; ++m) wsum += static_cast<double>(plan.weights[m]);
    }
    for (std::size_t m = b; m < e; ++m) {
      const double share = weighted ? f[p] * static_cast<double>(plan.weights[m]) / wsum
                                    : f[p] / static_cast<double>(e - b);
      if (share <= 0.0) continue;
      best = std::min(best, static_cast<double>(view_m[m]) / share);
      f[static_cast<std::size_t>(plan.edge_destination_positions[m])] += share;
    }
  }
  return best;
}

constexpr std::int32_t kProbeBudget = 62;  // a plain bisection over units below 2^62

// Find the largest fitting T in [lo, hi): alloc(lo) fits view_m; alloc(hi) does not.
// Fit checks compare integer amounts per member exactly. Increasing T never reduces a
// member's allocation, so these bounds remain valid regardless of probe choices.
// First try a floating-point estimate, then double step sizes to narrow the bounds, then
// bisect. Only probe choice is approximate; bounds always come from exact checks or the given
// ends. Allow an estimate or doubling probe only if it and the worst-case remaining bisection
// fit within kProbeBudget = 62, counting probes already made. Check the budget on return.
// Swap each fitting allocation into `best` when supplied; on return it holds alloc(result)
// whenever result is greater than the original lo.
Amount largest_fit(const PreparedRoute& plan, std::span<const Amount> view_m, Amount lo, Amount hi,
                   std::vector<Amount>& probe, std::vector<Amount>* best, std::int32_t& probes) {
  auto fits = [&](Amount t) {
    ++probes;
    if (!allocate(plan, t, probe, view_m)) return false;
    if (best != nullptr) std::swap(probe, *best);
    return true;
  };
  auto bisection_steps = [](Amount a, Amount b) {  // probes a plain bisection of (a, b) needs
    return static_cast<std::int32_t>(std::bit_width(static_cast<std::uint64_t>(b - a - 1)));
  };
  const std::int32_t start = probes;
  auto spare = [&](Amount a, Amount b) { return kProbeBudget - (probes - start) - bisection_steps(a, b) - 1; };
  if (hi - lo > 2 && spare(lo, hi) >= 0) {
    const double est = estimate_max_lossless_offer(plan, view_m);
    Amount t = lo + 1;
    if (est > static_cast<double>(lo) + 1.0) {
      t = est >= static_cast<double>(hi) - 1.0 ? hi - 1 : static_cast<Amount>(est);
    }
    if (fits(t)) {
      lo = t;
      for (Amount step = 1; hi - lo > 1 && spare(lo, hi) >= 0; step *= 2) {
        const Amount next = hi - lo > step ? lo + step : hi - 1;
        if (fits(next)) {
          lo = next;
        } else {
          hi = next;
          break;
        }
      }
    } else {
      hi = t;
      for (Amount step = 1; hi - lo > 1 && spare(lo, hi) >= 0; step *= 2) {
        const Amount next = hi - lo > step ? hi - step : lo + 1;
        if (fits(next)) {
          lo = next;
          break;
        }
        hi = next;
      }
    }
  }
  while (hi - lo > 1) {
    const Amount mid = lo + (hi - lo) / 2;
    if (fits(mid)) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  if (probes - start > kProbeBudget) {
    throw std::logic_error("T_max search: " + std::to_string(probes - start) + " probes, above the bisection bound 62");
  }
  return lo;
}

std::vector<Amount>& probe_scratch(std::size_t size) {
  thread_local std::vector<Amount> v;
  v.resize(size);
  return v;
}

} // namespace

LosslessBound max_lossless_offer(const PreparedRoute& plan, std::span<const Amount> view_m) {
  LosslessBound out;
  if (plan.forwarding_edges.empty()) return out;
  const Amount hi = total_bound(plan, view_m);
  auto& probe = probe_scratch(plan.forwarding_edges.size());
  // alloc(hi) fits only when every member can take it; else search below hi (the search's
  // probes are counted, this one is not, like the target probe of a Lossless step).
  if (allocate(plan, hi, probe, view_m)) {
    out.value = hi;
    return out;
  }
  out.value = largest_fit(plan, view_m, 0, hi, probe, nullptr, out.probes);
  return out;
}

LosslessStep lossless_step(const PreparedRoute& plan, Amount total, std::span<const Amount> view_m,
                           std::span<const Amount> own_m, std::optional<Amount> request) {
  LosslessStep out;
  if (plan.forwarding_edges.empty()) {
    out.edge_flows.assign(own_m.begin(), own_m.end());
    out.admitted = total;
    return out;
  }
  const Amount hi = total_bound(plan, view_m);
  const Amount target = request.has_value() ? std::min(total + *request, hi) : hi;
  if (target <= total) {
    out.admitted = total;
    out.edge_flows.assign(own_m.begin(), own_m.end());
    return out;
  }
  out.edge_flows.assign(plan.forwarding_edges.size(), 0);
  if (allocate(plan, target, out.edge_flows, view_m)) {
    out.admitted = target;
    return out;
  }
  // T_max in [total, target - 1]: alloc(total) = own_m fits, alloc(target) does not. Every
  // fitting probe leaves its allocation in out.edge_flows, so a result above total needs no further
  // evaluation.
  auto& probe = probe_scratch(plan.forwarding_edges.size());
  const Amount lo = largest_fit(plan, view_m, total, target, probe, &out.edge_flows, out.probes);
  out.admitted = lo;
  if (lo == total) out.edge_flows.assign(own_m.begin(), own_m.end());
  return out;
}

namespace {

// The per-root bookkeeping of an attributed forwarding (per thread): each position's arrivals as
// a list of (root, amount) entries, and a dense per-root accumulator with its touched roots.
struct Arrivals {
  std::vector<std::int32_t> arrival_list_heads;  // per position: its first entry, -1 none
  std::vector<std::int32_t> source_indices;  // per entry
  std::vector<Amount> amount;
  std::vector<std::int32_t> next;
  std::vector<Amount> source_totals;  // per root
  std::vector<std::int32_t> touched;
  std::vector<Amount> source_weights;  // per touched root
  std::vector<Amount> row;
  std::vector<Amount> arrive;
  std::vector<Amount> keep;

  void begin(std::size_t positions, std::size_t roots) {
    arrival_list_heads.assign(positions, -1);
    source_indices.clear();
    amount.clear();
    next.clear();
    source_totals.assign(roots, 0);
  }
  void push(std::int32_t p, std::int32_t r, Amount x) {
    source_indices.push_back(r);
    amount.push_back(x);
    next.push_back(arrival_list_heads[static_cast<std::size_t>(p)]);
    arrival_list_heads[static_cast<std::size_t>(p)] = static_cast<std::int32_t>(source_indices.size()) - 1;
  }
  // The roots arriving at position p (plus `own`, a root injecting there, when >= 0) with their
  // amounts in `acc`, in ascending root order.
  void gather(std::int32_t p, std::int32_t own) {
    touched.clear();
    for (auto k = arrival_list_heads[static_cast<std::size_t>(p)]; k >= 0; k = next[static_cast<std::size_t>(k)]) {
      const auto r = static_cast<std::size_t>(source_indices[static_cast<std::size_t>(k)]);
      if (source_totals[r] == 0) touched.push_back(static_cast<std::int32_t>(r));
      source_totals[r] += amount[static_cast<std::size_t>(k)];
    }
    if (own >= 0 && source_totals[static_cast<std::size_t>(own)] == 0) touched.push_back(own);
    std::sort(touched.begin(), touched.end());
  }
  void clear_touched() {
    for (auto r : touched) source_totals[static_cast<std::size_t>(r)] = 0;
  }
};

Arrivals& arrivals_scratch() {
  thread_local Arrivals a;
  return a;
}

} // namespace

namespace {

// The forwarding, compiled with and without the attribution (kAttribute: attr is non-null), so the
// one-source and Pinned calls carry none of its code.
template <bool kAttribute>
Forwarding forward_on(const PreparedRoute& plan, const Injection& in, std::span<const Amount> w_m, Attribution* attr) {
  const std::size_t H = plan.sources.size();
  const std::size_t P = plan.topological_order.size();
  const bool capped = !w_m.empty();
  if (!in.first_hop && in.per_root.size() != H) throw std::logic_error("forward: the injection is not root-aligned");
  // The first-hop experiment reads the view of the roots' members (a program without a member
  // injects nothing).
  if (in.first_hop && !capped && !plan.forwarding_edges.empty()) {
    throw std::logic_error("forward: the first-hop injection needs a view");
  }
  Forwarding out;
  out.edge_flows.assign(plan.forwarding_edges.size(), 0);
  out.edge_drops.assign(plan.forwarding_edges.size(), 0);
  if constexpr (kAttribute) {
    attr->offered.assign(H, 0);
    attr->source_drops_by_source.assign(H, 0);
    attr->delivered.assign(H, 0);
    attr->entries.assign(H, {});
  }
  const bool weighted = plan.allocation == RouteAllocation::Weighted;
  thread_local std::vector<Amount> inflow;
  thread_local std::vector<Amount> part;
  thread_local std::vector<std::int32_t> own_at;  // first hop: the root injecting at each position
  inflow.assign(P, 0);
  Arrivals& a = arrivals_scratch();
  if constexpr (kAttribute) a.begin(P, H);
  if (in.first_hop) own_at.assign(P, -1);
  // The roots' injections: a root without a member loses its own (none in the first-hop
  // experiment, where a root injects only on its members).
  for (std::size_t i = 0; i < H; ++i) {
    const auto p = plan.source_positions[i];
    const auto ri = static_cast<std::int32_t>(i);
    if (in.first_hop) {
      if (p < 0) continue;
      own_at[static_cast<std::size_t>(p)] = ri;
      Amount x = 0;
      for (auto m = plan.outgoing_begin[static_cast<std::size_t>(p)]; m < plan.outgoing_end[static_cast<std::size_t>(p)]; ++m) {
        x += w_m[static_cast<std::size_t>(m)];
      }
      out.offered += x;
      out.network_entered += x;
      if constexpr (kAttribute) attr->offered[i] = x;
      continue;
    }
    const Amount x = in.per_root[i];
    out.offered += x;
    if constexpr (kAttribute) attr->offered[i] = x;
    if (x == 0) continue;
    if (p < 0) {
      out.source_drops += x;
      if constexpr (kAttribute) attr->source_drops_by_source[i] = x;
      continue;
    }
    out.network_entered += x;
    inflow[static_cast<std::size_t>(p)] += x;
    if constexpr (kAttribute) a.push(p, ri, x);
  }
  for (std::size_t p = 0; p < P; ++p) {
    const auto pi = static_cast<std::int32_t>(p);
    const std::int32_t own = in.first_hop ? own_at[p] : -1;
    const Amount amount = inflow[p];
    if (plan.is_destination[p] != 0) {
      out.delivered += amount;
      if (kAttribute && amount != 0) {
        a.gather(pi, -1);
        for (auto r : a.touched) attr->delivered[static_cast<std::size_t>(r)] += a.source_totals[static_cast<std::size_t>(r)];
        a.clear_touched();
      }
      continue;
    }
    if (amount == 0 && own < 0) continue;
    const auto b = static_cast<std::size_t>(plan.outgoing_begin[p]);
    const auto e = static_cast<std::size_t>(plan.outgoing_end[p]);
    const std::size_t k = e - b;
    // The members with a share: round robin writes the first n (the rest receive nothing),
    // highest averages all k; a root of the first-hop experiment injects on every member.
    part.resize(k);
    std::size_t n = 0;
    if (amount != 0) {
      if (weighted) {
        highest_averages(amount, std::span<const Amount>(plan.weights.data() + b, k), part);
        n = k;
      } else {
        n = round_robin(amount, part);
      }
    }
    if (own >= 0) {
      std::fill(part.begin() + static_cast<std::ptrdiff_t>(n), part.end(), Amount {0});
      n = k;
    }
    if constexpr (kAttribute) {
      a.gather(pi, own);
      a.source_weights.resize(a.touched.size());
      a.row.resize(a.touched.size());
      a.arrive.resize(a.touched.size());
      a.keep.resize(a.touched.size());
    }
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t m = b + i;
      const Amount mine = own >= 0 ? w_m[m] : 0;
      const Amount x = part[i] + mine;
      if (x == 0) continue;
      const Amount c = capped ? std::min(x, w_m[m]) : x;
      out.edge_flows[m] = c;
      out.edge_drops[m] = x - c;
      const auto head = plan.edge_destination_positions[m];
      inflow[static_cast<std::size_t>(head)] += c;
      if constexpr (!kAttribute) continue;
      // The member's amount by highest averages over the remaining inflows (then the root's own
      // first-hop injection on it), and what it carries by highest averages over those arrivals.
      const std::size_t t = a.touched.size();
      for (std::size_t j = 0; j < t; ++j) a.source_weights[j] = a.source_totals[static_cast<std::size_t>(a.touched[j])];
      highest_averages(part[i], a.source_weights, a.row);
      for (std::size_t j = 0; j < t; ++j) {
        const auto r = a.touched[j];
        a.source_totals[static_cast<std::size_t>(r)] -= a.row[j];
        a.arrive[j] = a.row[j] + (r == own ? mine : 0);
      }
      if (c == x) {
        std::copy(a.arrive.begin(), a.arrive.end(), a.keep.begin());
      } else {
        highest_averages(c, a.arrive, a.keep);
      }
      for (std::size_t j = 0; j < t; ++j) {
        if (a.arrive[j] == 0) continue;
        const auto r = a.touched[j];
        attr->entries[static_cast<std::size_t>(r)].push_back(
            Attribution::Entry {static_cast<std::int32_t>(m), a.keep[j], a.arrive[j] - a.keep[j]});
        if (a.keep[j] != 0) a.push(head, r, a.keep[j]);
      }
    }
    if constexpr (kAttribute) {
      for (auto r : a.touched) {
        if (a.source_totals[static_cast<std::size_t>(r)] != 0) {
          throw std::logic_error("forward: a router's inflow is not fully attributed");
        }
      }
    }
  }
  if (in.first_hop) out.offered = out.network_entered;
  if constexpr (kAttribute) {
    for (auto& rows : attr->entries) {
      std::sort(rows.begin(), rows.end(), [](const auto& x, const auto& y) { return x.edge_position < y.edge_position; });
    }
  }
  return out;
}

} // namespace

Forwarding forward(const PreparedRoute& plan, const Injection& in, std::span<const Amount> w_m, Attribution* attr) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  ++thread_entry_work().forwards;
  if (attr) ++thread_entry_work().attribution;
#endif
  return attr != nullptr ? forward_on<true>(plan, in, w_m, attr) : forward_on<false>(plan, in, w_m, attr);
}

} // namespace netgraph::core::detail
