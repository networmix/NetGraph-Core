/*
  Apportionment and the placement kernels in units. T is a plan's total injection.

  Every amount is an int64 number of units below 2^62: T is bounded by the view on the roots'
  members (the graph's capacities sum below 2^62), each node's inflow by T, and each share by that
  inflow. Highest averages compares
  w_i/d_i and w_j/d_j through w_i*d_j and w_j*d_i (operands < 2^63, products < 2^126). The
  lower quota floor(I*w/W) uses a 128-bit product, with inflow I, weight w < 2^62 and the exact
  128-bit weight sum W (explicit weights are each below 2^62, not their sum).
*/
#include "netgraph/core/detail/kernels.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>

#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/flow_core.hpp"
#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"

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

} // namespace

void highest_averages(Amount total, std::span<const Amount> weights, std::span<Amount> out) {
  const std::size_t k = weights.size();
  U128 wsum {};
  for (Amount w : weights) wsum = add(wsum, U128 {0, static_cast<std::uint64_t>(w)});
  if (wsum == U128 {}) {  // all-zero weights give nothing
    std::fill(out.begin(), out.begin() + static_cast<std::ptrdiff_t>(k), 0);
    return;
  }
  Amount given = 0;
  for (std::size_t m = 0; m < k; ++m) {
    // Lower quota floor(total * w / W): exactly the quotients w / j >= W / total.
    out[m] = static_cast<Amount>(div_u128(mul_u64(static_cast<std::uint64_t>(total),
                                                  static_cast<std::uint64_t>(weights[m])),
                                          wsum));
    given += out[m];
  }
  // Fewer than k (the positive weights' count) units are left; a zero weight's quotient 0 never
  // beats a positive one, so it gets none of them.
  Amount left = total - given;
  if (left == 0) return;
  if (k <= 16) {
    while (left-- > 0) {
      std::size_t best = 0;
      for (std::size_t m = 1; m < k; ++m) {
        if (before(static_cast<std::uint64_t>(weights[m]), static_cast<std::uint64_t>(out[m]) + 1, m,
                   static_cast<std::uint64_t>(weights[best]), static_cast<std::uint64_t>(out[best]) + 1, best)) {
          best = m;
        }
      }
      ++out[best];
    }
    return;
  }
  std::vector<std::size_t> heap(k);
  std::iota(heap.begin(), heap.end(), std::size_t{0});
  // std heaps put the largest first: "a < b" when b goes first.
  auto later = [&](std::size_t a, std::size_t b) {
    return before(static_cast<std::uint64_t>(weights[b]), static_cast<std::uint64_t>(out[b]) + 1, b,
                  static_cast<std::uint64_t>(weights[a]), static_cast<std::uint64_t>(out[a]) + 1, a);
  };
  std::make_heap(heap.begin(), heap.end(), later);
  while (left-- > 0) {
    std::pop_heap(heap.begin(), heap.end(), later);
    ++out[heap.back()];
    std::push_heap(heap.begin(), heap.end(), later);
  }
}

// ---------------------------------------------------------------- plans

void check_kernel_combination(FlowPlacement placement, Admission admission) {
  if (admission == Admission::Lossy && placement == FlowPlacement::MaxFlow) {
    throw std::invalid_argument("MAX_FLOW requires admission=LOSSLESS");
  }
}

void check_kernel(FlowPlacement placement, Admission admission) {
  check_member(placement, "placement");
  check_member(admission, "admission");
  check_kernel_combination(placement, admission);
}

KernelPlan::Layout KernelPlan::point(std::size_t num_roots, std::size_t num_dsts, std::size_t positions,
                                     std::size_t num_members) noexcept {
  std::int32_t* at = store_.data();
  Layout w;
  roots = std::span<const NodeId>(at, num_roots);
  at += num_roots;
  w.root_pos = std::span<std::int32_t>(at, num_roots);
  root_pos = w.root_pos;
  at += num_roots;
  dsts = std::span<const NodeId>(at, num_dsts);
  at += num_dsts;
  w.topo = std::span<NodeId>(at, positions);
  at += positions;
  w.out_begin = std::span<std::int32_t>(at, positions);
  at += positions;
  w.out_end = std::span<std::int32_t>(at, positions);
  at += positions;
  w.head_pos = std::span<std::int32_t>(at, num_members);
  at += num_members;
  w.terminal = std::span<std::uint8_t>(reinterpret_cast<std::uint8_t*>(at), positions);
  topo = w.topo;
  out_begin = w.out_begin;
  out_end = w.out_end;
  head_pos = w.head_pos;
  terminal = w.terminal;
  return w;
}

KernelPlan::Layout KernelPlan::layout(std::span<const NodeId> roots_in, std::span<const NodeId> destinations,
                                      std::size_t positions, std::size_t num_members) {
  const std::size_t words =
      2 * roots_in.size() + destinations.size() + 3 * positions + num_members + (positions + 3) / 4;
  std::vector<std::int32_t> store(words, 0);  // the one allocation
  store_.swap(store);
  const Layout w = point(roots_in.size(), destinations.size(), positions, num_members);
  std::copy(roots_in.begin(), roots_in.end(), store_.begin());
  std::fill(w.root_pos.begin(), w.root_pos.end(), -1);
  std::copy(destinations.begin(), destinations.end(), store_.begin() + static_cast<std::ptrdiff_t>(2 * roots_in.size()));
  return w;
}

KernelPlan::KernelPlan(const KernelPlan& other)
    : placement(other.placement), admission(other.admission), route(other.route), members(other.members),
      weights(other.weights), store_(other.store_) {
  (void)point(other.roots.size(), other.dsts.size(), other.topo.size(), other.head_pos.size());
}

KernelPlan::KernelPlan(KernelPlan&& other) noexcept
    : placement(other.placement), admission(other.admission), route(std::move(other.route)),
      members(std::move(other.members)), weights(std::move(other.weights)), store_(std::move(other.store_)) {
  (void)point(other.roots.size(), other.dsts.size(), other.topo.size(), other.head_pos.size());
  (void)other.point(0, 0, 0, 0);
}

KernelPlan& KernelPlan::operator=(const KernelPlan& other) {
  if (this != &other) {
    KernelPlan copy(other);
    *this = std::move(copy);
  }
  return *this;
}

KernelPlan& KernelPlan::operator=(KernelPlan&& other) noexcept {
  if (this != &other) {
    placement = other.placement;
    admission = other.admission;
    route = std::move(other.route);
    members = std::move(other.members);
    weights = std::move(other.weights);
    const std::size_t nr = other.roots.size();
    const std::size_t nd = other.dsts.size();
    const std::size_t np = other.topo.size();
    const std::size_t nm = other.head_pos.size();
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
  std::uint32_t epoch {0};
  std::vector<Node> node;
  std::vector<NodeId> order;
  std::vector<std::pair<NodeId, std::int32_t>> stack;
  std::vector<NodeId> member_tail;
  std::vector<NodeId> member_head;

  void begin(std::size_t n) {
    if (node.size() < n) node.resize(n);  // strong guarantee; new slots carry stamp 0
    if (++epoch == 0) {
      for (auto& v : node) v.stamp = 0;
      epoch = 1;
    }
  }
  [[nodiscard]] bool seen(NodeId v) const noexcept { return node[static_cast<std::size_t>(v)].stamp == epoch; }
  Node& at(NodeId v) noexcept { return node[static_cast<std::size_t>(v)]; }
  void touch(NodeId v) noexcept {
    auto& x = at(v);
    x = Node {};
    x.stamp = epoch;
  }
};

DagScratch& dag_scratch() {
  thread_local DagScratch s;
  return s;
}

// Build plan.members from the DAG's incoming-edge lists: keep edges accepted by keep(e) that
// lie on a walk from a root to a destination; exclude edges leaving destinations. No node
// renumbering, extra arc arrays or adjacency offsets are needed. A backward DFS from destinations
// lists each node after its parents (postorder, also a topological order here).
// A forward pass then marks nodes reachable from the roots. Members ascend by EdgeId, so one
// tail's members are consecutive. Any topological order gives the same allocation because each
// node splits only after receiving all its inflow.
// With `cost`, the forward pass computes minimum cost to each node from its parents' costs
// (dynamic programming, DP), giving the minimum root-to-destination cost over members, or
// kUnreachable if none exist. Without `out` (route_cost), stop there and build no plan.
template <class Keep>
void build_split(const Graph& g, const PredDAG& dag, std::span<const NodeId> roots, std::span<const NodeId> dsts,
                 Keep keep, KernelPlan* out, Cost* cost) {
  if (cost != nullptr) *cost = kUnreachable;
  const auto parents = dag.parents();
  const auto via = dag.via_edges();
  DagScratch& s = dag_scratch();
  s.begin(static_cast<std::size_t>(g.num_nodes()));
  // The DAG's sparse storage scattered into the per-node slots (O(member nodes)); a node outside
  // it has no entries when first touched (a destination, src or a parent the search meets).
  {
    const auto nodes = PredDAGBuilder::member_nodes(dag);
    const auto offs = PredDAGBuilder::member_offsets(dag);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      s.touch(nodes[i]);
      s.at(nodes[i]).entry_begin = offs[i];
      s.at(nodes[i]).entry_end = offs[i + 1];
    }
  }
  for (NodeId d : dsts) {
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
  for (NodeId d : dsts) {
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
  const auto edge_cost = g.cost();
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
  KernelPlan& plan = *out;
  plan.members.clear();
  // Members, ascending: the DAG's edges kept, from a reached tail (not a destination) to a head
  // that reaches a destination.
  const auto edge_src = g.edge_src();
  const auto edge_dst = g.edge_dst();
  // Size vectors to known bounds (DAG edges here, reached nodes below), then write by index
  // to avoid push_back overhead in these loops.
  const auto asc = PredDAGBuilder::ascending(dag);
  const std::size_t bound = asc.size();
  plan.members.resize(bound);
  s.member_tail.resize(bound);
  s.member_head.resize(bound);
  std::size_t n_members = 0;
  for (EdgeId e : asc) {
    const auto ei = static_cast<std::size_t>(e);
    const NodeId u = edge_src[ei];
    const NodeId v = edge_dst[ei];
    if (!s.seen(u) || s.at(u).fwd == 0 || s.at(u).target != 0) continue;
    if (!s.seen(v) || s.at(v).state != 2 || !keep(e)) continue;
    plan.members[n_members] = e;
    s.member_tail[n_members] = u;
    s.member_head[n_members] = v;
    ++n_members;
  }
  plan.members.resize(n_members);
  s.member_tail.resize(n_members);
  s.member_head.resize(n_members);
  const std::size_t num_members = n_members;
  if (num_members == 0) {
    (void)plan.layout(roots, dsts, 0, 0);  // every root without a member
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
  const KernelPlan::Layout w = plan.layout(roots, dsts, P, num_members);
  std::size_t pos = 0;
  for (NodeId v : s.order) {
    const auto& nv = s.at(v);
    if (nv.fwd == 0) continue;
    w.topo[pos] = v;
    w.terminal[pos] = nv.target;
    const std::int32_t b = nv.count == 0 ? 0 : nv.first;
    w.out_begin[pos] = b;
    w.out_end[pos] = b + nv.count;
    ++pos;
  }
  // Positions follow the postorder (every node after its parents; the DFS rejected cycles), so a
  // member's tail precedes its head.
  for (std::size_t m = 0; m < num_members; ++m) w.head_pos[m] = s.at(s.member_head[m]).position;
  // A root with a member has a position (it reaches a destination and is reached); a root
  // without one keeps -1.
  for (std::size_t i = 0; i < roots.size(); ++i) {
    const auto& nr = s.at(roots[i]);
    if (nr.fwd != 0 && nr.target == 0) w.root_pos[i] = nr.position;
  }
}

} // namespace

std::vector<EdgeId> prune(const Graph& g, const PredDAG& dag, NodeId src, std::span<const NodeId> dsts,
                          const Eligibility* elig) {
  KernelPlan plan;
  build_split(g, dag, {&src, 1}, dsts, [&](EdgeId e) { return elig == nullptr || elig->edge_ok(e); }, &plan, nullptr);
  return std::move(plan.members);
}

Cost route_cost(const Graph& g, const PredDAG& route, NodeId src, std::span<const NodeId> dsts) {
  Cost cost = kUnreachable;
  build_split(g, route, {&src, 1}, dsts, [](EdgeId) { return true; }, nullptr, &cost);
  return cost;
}

namespace {

KernelPlan plan_of(const Graph& g, FlowPlacement placement, Admission admission, const WeightSource& weights,
                   std::span<const NodeId> roots, std::span<const NodeId> dsts, const PredDAG& route, Cost* cost) {
  check_kernel(placement, admission);
  const auto edges = PredDAGBuilder::ascending(route);
  const bool weighted = placement == FlowPlacement::Weighted;
  const bool by_view = weighted && weights.basis == WeightBasis::ResidualSnapshot;
  if (by_view && weights.values.size() != edges.size()) {
    throw std::logic_error("make_plan: the creation view is not route-aligned");
  }
  if (weighted && weights.basis == WeightBasis::Explicit &&
      weights.values.size() != static_cast<std::size_t>(g.num_edges())) {
    throw std::logic_error("make_plan: the explicit weights are not edge-aligned");
  }
  const auto capacity = g.capacity_units();
  KernelPlan plan;
  plan.placement = placement;
  plan.admission = admission;
  plan.route.assign(edges.begin(), edges.end());
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
  // Membership: MaxFlow and EqualBalanced every route edge (an edge of capacity 0 is a member
  // that carries nothing), Weighted the edges of positive weight.
  auto member = [&](EdgeId e) { return !weighted || weight(e) >= 1; };
  build_split(g, route, roots, dsts, member, &plan, cost);
  // The route's cost is the members' only when every route edge is a member.
  if (cost != nullptr && plan.members.size() != plan.route.size()) *cost = route_cost(g, route, roots[0], dsts);
  if (weighted) {
    plan.weights.reserve(plan.members.size());
    for (EdgeId e : plan.members) plan.weights.push_back(weight(e));
  }
  return plan;
}

} // namespace

KernelPlan make_plan(const Graph& g, FlowPlacement placement, Admission admission, const WeightSource& weights,
                     NodeId src, std::span<const NodeId> dsts, const PredDAG& route, Cost* cost) {
  return plan_of(g, placement, admission, weights, {&src, 1}, dsts, route, cost);
}

KernelPlan make_plan(const Graph& g, FlowPlacement placement, Admission admission, const WeightSource& weights,
                     std::span<const NodeId> roots, std::span<const NodeId> dsts, const PredDAG& route) {
  return plan_of(g, placement, admission, weights, roots, dsts, route, nullptr);
}

KernelPlan make_route_plan(const Graph& g, FlowPlacement placement, WeightBasis basis, std::span<const Amount> edge_values,
                           NodeId src, std::span<const NodeId> dsts, const PredDAG& dag, Cost* cost) {
  check_kernel(placement, Admission::Lossless);
  const auto capacity = g.capacity_units();
  KernelPlan plan;
  plan.placement = placement;
  plan.admission = Admission::Lossless;
  build_split(g, dag, {&src, 1}, dsts, [](EdgeId) { return true; }, &plan, cost);
  plan.route = plan.members;
  if (placement != FlowPlacement::Weighted) return plan;
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
  if (std::all_of(plan.route.begin(), plan.route.end(), [&](EdgeId e) { return weight(e) >= 1; })) {
    plan.weights.reserve(plan.members.size());
    for (EdgeId e : plan.members) plan.weights.push_back(weight(e));
    return plan;
  }
  KernelPlan weighted;
  build_split(g, dag, {&src, 1}, dsts, [&](EdgeId e) { return weight(e) >= 1; }, &weighted, nullptr);
  weighted.placement = placement;
  weighted.admission = Admission::Lossless;
  weighted.route = std::move(plan.route);
  weighted.weights.reserve(weighted.members.size());
  for (EdgeId e : weighted.members) weighted.weights.push_back(weight(e));
  return weighted;
}

KernelPlan restrict_plan(const Graph& g, const KernelPlan& plan, std::span<const NodeId> roots,
                         std::span<const NodeId> dsts, const Eligibility* elig) {
  KernelPlan out;
  out.placement = plan.placement;
  out.admission = plan.admission;
  const PredDAG members = PredDAGBuilder::of_edges(g, plan.members);
  build_split(g, members, roots, dsts, [&](EdgeId e) { return elig == nullptr || elig->edge_ok(e); }, &out, nullptr);
  out.route = plan.route;
  if (plan.placement == FlowPlacement::Weighted) {
    // The kept members' frozen weights (both member lists ascend).
    out.weights.resize(out.members.size());
    std::size_t k = 0;
    for (std::size_t m = 0; m < out.members.size(); ++m) {
      while (plan.members[k] != out.members[m]) ++k;
      out.weights[m] = plan.weights[k];
    }
  }
  return out;
}

// ---------------------------------------------------------------- CostOnly creation

namespace {

// Per-thread state of cost_only_creation. Per node (by NodeId) a slot valid only under the
// current epoch (DagScratch's pattern: nothing is cleared between calls except on stamp
// wraparound, a call that throws leaves stale slots behind an old epoch, and the one array of
// slots grows whole or not at all); every other buffer is indexed by the reached nodes' local
// indices (their rank in ascending NodeId) or by route position, and is rewritten before it is
// read.
struct CreationScratch {
  struct Node {
    std::uint32_t stamp {0};
    std::uint8_t reached {0};
    std::uint8_t target {0};
    std::int32_t local {0};  // rank among the reached nodes (ascending NodeId)
    std::int32_t lo {0};     // the node's support out-edges: support[lo, hi)
    std::int32_t hi {0};
  };
  std::uint32_t epoch {0};
  std::vector<Node> node;
  std::vector<NodeId> reached;  // the walk's nodes (its queue), then ascending
  std::vector<EdgeId> route;                               // ascending
  std::vector<std::int32_t> tail;                          // per route edge: its tail's local index
  std::vector<std::int32_t> head;                          // per route edge: its head's local index
  std::vector<std::int32_t> parent;                        // per entry: its parent's local index
  std::vector<std::int32_t> entry;                         // per local node + 1: its entries (by child)
  std::vector<std::int32_t> fill;
  std::vector<EdgeId> via;                                 // the entries: by child, ascending via
  std::vector<std::uint8_t> state;                         // per local node: DFS 1 open, 2 done
  std::vector<std::pair<std::int32_t, std::int32_t>> dfs;  // (local node, next entry)
  std::vector<std::int32_t> order;                         // local nodes, every one after its parents
  std::vector<std::int32_t> position;                      // per local node
  std::vector<std::int32_t> first;                         // per local node: its first member
  std::vector<std::int32_t> count;                         // per local node: its number of members
  std::vector<Cost> dist;                                  // per local node

  void begin(std::size_t n) {
    if (node.size() < n) node.resize(n);  // strong guarantee; new slots carry stamp 0
    if (++epoch == 0) {
      for (auto& v : node) v.stamp = 0;
      epoch = 1;
    }
  }
  Node& at(NodeId v) noexcept { return node[static_cast<std::size_t>(v)]; }
  // The slot of v, reset on its first use in this call.
  Node& touch(NodeId v) noexcept {
    auto& x = at(v);
    if (x.stamp != epoch) {
      x = Node {};
      x.stamp = epoch;
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

std::optional<CreationProducts> cost_only_creation(const Graph& g, const std::shared_ptr<const PredDAG>& support,
                                                   NodeId src, std::span<const NodeId> dsts,
                                                   const PathSelection& selection, PriorityFn view,
                                                   FlowPlacement placement, Admission admission,
                                                   const WeightSource& weights) {
  const bool weighted = placement == FlowPlacement::Weighted;
  if (!selection.multipath || (weighted && weights.basis == WeightBasis::ResidualSnapshot)) return std::nullopt;
  check_kernel(placement, admission);
  // A Weighted route edge's weight (static membership: Capacity or Explicit).
  auto weight = [&](EdgeId e) -> Amount {
    const auto ei = static_cast<std::size_t>(e);
    const Amount cap = g.capacity_units()[ei];
    if (weights.basis == WeightBasis::Explicit) return cap >= 1 ? weights.values[ei] : 0;
    return cap;
  };
  const auto sup = PredDAGBuilder::ascending(*support);
  const auto row = g.row_offsets();
  const auto edst = g.edge_dst();
  const auto edge_cost = g.cost();
  CreationScratch& s = creation_scratch();
  s.begin(static_cast<std::size_t>(g.num_nodes()));
  for (NodeId d : dsts) s.touch(d).target = 1;
  // support_route's walk over the support from src: every reached node keeps its support
  // out-edges as one range (EdgeIds rank edges by source first: one tail's edges are consecutive).
  s.reached.clear();
  s.touch(src).reached = 1;
  s.reached.push_back(src);
  bool any_target = false;
  std::size_t m = 0;
  std::size_t hint = 0;
  for (std::size_t q = 0; q < s.reached.size(); ++q) {  // the reached list is the walk's queue
    const NodeId u = s.reached[q];
    const auto ui = static_cast<std::size_t>(u);
    const std::size_t lo = lower_bound_from(sup, hint, row[ui]);
    std::size_t hi = lo;
    for (; hi < sup.size() && sup[hi] < row[ui + 1]; ++hi) {
      const NodeId v = edst[static_cast<std::size_t>(sup[hi])];
      auto& nv = s.touch(v);
      if (nv.reached != 0) continue;
      nv.reached = 1;
      s.reached.push_back(v);
    }
    auto& nu = s.at(u);
    any_target = any_target || nu.target != 0;
    nu.lo = static_cast<std::int32_t>(lo);
    nu.hi = static_cast<std::int32_t>(hi);
    m += hi - lo;
    hint = lo;
  }
  if (!any_target || m == 0) return std::nullopt;
  // The route in ascending EdgeId: the reached nodes' ranges in ascending NodeId (ascending in the
  // support, checked), each with Collapse when multi_edge is false (support_route's: one member per
  // (u, v) by select_adjacency on the view), every route edge with its tail's local index. A
  // destination with route edges, or a route edge that is no member candidate, means build_split
  // would not keep the whole route: the caller's path then runs.
  std::sort(s.reached.begin(), s.reached.end());
  const std::size_t r = s.reached.size();
  for (std::size_t i = 0; i < r; ++i) s.at(s.reached[i]).local = static_cast<std::int32_t>(i);
  s.route.resize(m);
  s.tail.resize(m);
  s.first.resize(r);
  s.count.resize(r);
  const bool collapse = !selection.multi_edge;
  const bool prefer = selection.tie_break == EdgeTieBreak::PreferHigherPriority;
  std::size_t w = 0;
  std::int32_t prev_hi = 0;
  for (std::size_t i = 0; i < r; ++i) {
    const auto& nu = s.at(s.reached[i]);
    s.first[i] = static_cast<std::int32_t>(w);
    if (nu.hi > nu.lo) {
      if (nu.lo < prev_hi || nu.target != 0) return std::nullopt;
      prev_hi = nu.hi;
    }
    for (auto k = nu.lo; k < nu.hi;) {
      EdgeId best = sup[static_cast<std::size_t>(k)];
      auto j = k + 1;
      if (collapse) {
        const NodeId v = edst[static_cast<std::size_t>(best)];
        for (; j < nu.hi && edst[static_cast<std::size_t>(sup[static_cast<std::size_t>(j)])] == v; ++j) {
          if (replaces(g, sup[static_cast<std::size_t>(j)], best, prefer, view)) best = sup[static_cast<std::size_t>(j)];
        }
      }
      if (weighted && weight(best) < 1) return std::nullopt;
      s.route[w] = best;
      s.tail[w] = static_cast<std::int32_t>(i);
      ++w;
      k = j;
    }
    s.count[i] = static_cast<std::int32_t>(w) - s.first[i];
  }
  const std::span<const EdgeId> route(s.route.data(), w);
  // The entries grouped by child (ascending NodeId), ascending via within a child: a counting pass
  // over the heads' local indices, the route scattered in ascending order with each entry's parent.
  s.head.resize(w);
  s.entry.assign(r + 1, 0);
  for (std::size_t k = 0; k < w; ++k) {
    const auto h = s.at(edst[static_cast<std::size_t>(route[k])]).local;
    s.head[k] = h;
    ++s.entry[static_cast<std::size_t>(h) + 1];
  }
  std::size_t children = 0;
  for (std::size_t i = 0; i < r; ++i) {
    if (s.entry[i + 1] != 0) ++children;
    s.entry[i + 1] += s.entry[i];
  }
  s.fill.assign(s.entry.begin(), s.entry.end() - 1);
  s.via.resize(w);
  s.parent.resize(w);
  for (std::size_t k = 0; k < w; ++k) {
    const auto at = static_cast<std::size_t>(s.fill[static_cast<std::size_t>(s.head[k])]++);
    s.via[at] = route[k];
    s.parent[at] = s.tail[k];
  }
  // build_split's order: a backward DFS from the destinations, in their order, over the entries,
  // listing every node after its parents (postorder). A destination the walk did not reach has no
  // entry and no position there, so it is skipped.
  s.state.assign(r, 0);
  s.order.clear();
  for (NodeId d : dsts) {
    const auto& nd = s.at(d);
    if (nd.reached == 0) continue;
    const auto di = nd.local;
    if (s.state[static_cast<std::size_t>(di)] != 0) continue;
    s.state[static_cast<std::size_t>(di)] = 1;
    s.dfs.clear();
    s.dfs.emplace_back(di, s.entry[static_cast<std::size_t>(di)]);
    while (!s.dfs.empty()) {
      const auto [v, k] = s.dfs.back();
      if (k < s.entry[static_cast<std::size_t>(v) + 1]) {
        s.dfs.back().second = k + 1;
        const auto u = s.parent[static_cast<std::size_t>(k)];
        const auto su = s.state[static_cast<std::size_t>(u)];
        if (su == 0) {
          s.state[static_cast<std::size_t>(u)] = 1;
          s.dfs.emplace_back(u, s.entry[static_cast<std::size_t>(u)]);
        } else if (su == 1) {
          return std::nullopt;  // a cycle (build_split throws)
        }
      } else {
        s.state[static_cast<std::size_t>(v)] = 2;
        s.order.push_back(v);
        s.dfs.pop_back();
      }
    }
  }
  const auto src_local = s.at(src).local;
  if (s.order.size() != r || s.order.front() != src_local) return std::nullopt;
  s.position.resize(r);
  for (std::size_t p = 0; p < r; ++p) s.position[static_cast<std::size_t>(s.order[p])] = static_cast<std::int32_t>(p);
  // make_plan's cost: the minimum walk cost from src to a destination over the members, by the
  // order (every member is a route edge, so route_cost is not called).
  s.dist.resize(r);
  Cost best = kUnreachable;
  for (const auto v : s.order) {
    const auto vi = static_cast<std::size_t>(v);
    Cost dv = 0;
    if (v != src_local) {
      dv = kUnreachable;
      for (auto k = s.entry[vi]; k < s.entry[vi + 1]; ++k) {
        const auto ki = static_cast<std::size_t>(k);
        dv = std::min(dv, s.dist[static_cast<std::size_t>(s.parent[ki])] + edge_cost[static_cast<std::size_t>(s.via[ki])]);
      }
    }
    s.dist[vi] = dv;
    if (s.at(s.reached[vi]).target != 0) best = std::min(best, dv);
  }
  CreationProducts out;
  out.cost = best;
  KernelPlan& plan = out.plan;
  plan.placement = placement;
  plan.admission = admission;
  plan.route.assign(route.begin(), route.end());
  plan.members.assign(route.begin(), route.end());
  if (weighted) {
    plan.weights.resize(route.size());
    for (std::size_t k = 0; k < route.size(); ++k) plan.weights[k] = weight(route[k]);
  }
  const KernelPlan::Layout lay = plan.layout({&src, 1}, dsts, r, route.size());
  lay.root_pos[0] = 0;  // src comes first in the order
  for (std::size_t p = 0; p < r; ++p) {
    const auto v = static_cast<std::size_t>(s.order[p]);
    const NodeId x = s.reached[v];
    lay.topo[p] = x;
    lay.terminal[p] = s.at(x).target;
    const std::int32_t b = s.count[v] == 0 ? 0 : s.first[v];
    lay.out_begin[p] = b;
    lay.out_end[p] = b + s.count[v];
  }
  for (std::size_t k = 0; k < route.size(); ++k) lay.head_pos[k] = s.position[static_cast<std::size_t>(s.head[k])];
  // The record bundle: the support itself when the route is all of it (bundle_of's share).
  if (route.size() == sup.size()) {
    out.bundle = support;
  } else {
    std::vector<NodeId> nodes(children);
    std::vector<std::int32_t> offsets(children + 1);
    offsets[0] = 0;
    for (std::size_t i = 0, c = 0; i < r; ++i) {
      if (s.entry[i + 1] == s.entry[i]) continue;
      nodes[c] = s.reached[i];
      offsets[++c] = s.entry[i + 1];
    }
    std::vector<EdgeId> via(s.via.begin(), s.via.begin() + static_cast<std::ptrdiff_t>(w));
    std::vector<NodeId> parents(w);
    for (std::size_t k = 0; k < w; ++k) parents[k] = s.reached[static_cast<std::size_t>(s.parent[k])];
    std::vector<EdgeId> asc(route.begin(), route.end());
    out.bundle = std::make_shared<const PredDAG>(PredDAGBuilder::of_parts(
        g, std::move(nodes), std::move(offsets), std::move(parents), std::move(via), std::move(asc)));
  }
  return out;
}

// ---------------------------------------------------------------- allocation

namespace {

// The view summed over the roots' out-members, and at most the number of roots before the first
// root without a member: no larger total fits (that root would receive a share with no route).
Amount total_bound(const KernelPlan& plan, std::span<const Amount> view_m) {
  Amount hi = 0;
  for (std::size_t i = 0; i < plan.roots.size(); ++i) {
    const auto p = plan.root_pos[i];
    if (p < 0) return std::min(hi, static_cast<Amount>(i));
    for (auto m = plan.out_begin[static_cast<std::size_t>(p)]; m < plan.out_end[static_cast<std::size_t>(p)]; ++m) {
      hi += view_m[static_cast<std::size_t>(m)];
    }
  }
  return hi;
}

} // namespace

bool allocate(const KernelPlan& plan, Amount total, std::span<Amount> out, std::span<const Amount> limit) {
  std::fill(out.begin(), out.end(), 0);
  if (total == 0) return true;
  const bool weighted = plan.placement == FlowPlacement::Weighted;
  const bool limited = !limit.empty();
  thread_local std::vector<Amount> inflow;
  inflow.assign(plan.topo.size(), 0);
  // The roots' round-robin shares (only the positive ones are visited). A positive share of a
  // root without a member has no route; every caller stays at or below total_bound, which no such
  // share reaches.
  bool routed = true;
  round_robin(total, plan.roots.size(), [&](std::size_t i, Amount x) {
    const auto p = plan.root_pos[i];
    if (p < 0) {
      routed = false;
      return false;
    }
    inflow[static_cast<std::size_t>(p)] += x;
    return true;
  });
  if (!routed) throw std::logic_error("allocate: a root without a member receives a share");
  for (std::size_t p = 0; p < plan.topo.size(); ++p) {
    const Amount amount = inflow[p];
    if (amount == 0 || plan.terminal[p] != 0) continue;
    const auto b = static_cast<std::size_t>(plan.out_begin[p]);
    const auto e = static_cast<std::size_t>(plan.out_end[p]);
    const std::size_t k = e - b;
    if (weighted) {
      highest_averages(amount, std::span<const Amount>(plan.weights.data() + b, k), out.subspan(b, k));
      for (std::size_t m = b; m < e; ++m) {
        const Amount x = out[m];
        if (limited && x > limit[m]) return false;
        inflow[static_cast<std::size_t>(plan.head_pos[m])] += x;
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
        inflow[static_cast<std::size_t>(plan.head_pos[m])] += x;
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
double estimate_t_max(const KernelPlan& plan, std::span<const Amount> view_m) {
  thread_local std::vector<double> f;
  f.assign(plan.topo.size(), 0.0);
  for (auto p : plan.root_pos) {
    if (p >= 0) f[static_cast<std::size_t>(p)] += 1.0 / static_cast<double>(plan.roots.size());
  }
  const bool weighted = plan.placement == FlowPlacement::Weighted;
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t p = 0; p < plan.topo.size(); ++p) {
    if (f[p] == 0.0 || plan.terminal[p] != 0) continue;
    const auto b = static_cast<std::size_t>(plan.out_begin[p]);
    const auto e = static_cast<std::size_t>(plan.out_end[p]);
    double wsum = 0.0;
    if (weighted) {
      for (std::size_t m = b; m < e; ++m) wsum += static_cast<double>(plan.weights[m]);
    }
    for (std::size_t m = b; m < e; ++m) {
      const double share = weighted ? f[p] * static_cast<double>(plan.weights[m]) / wsum
                                    : f[p] / static_cast<double>(e - b);
      if (share <= 0.0) continue;
      best = std::min(best, static_cast<double>(view_m[m]) / share);
      f[static_cast<std::size_t>(plan.head_pos[m])] += share;
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
Amount largest_fit(const KernelPlan& plan, std::span<const Amount> view_m, Amount lo, Amount hi,
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
    const double est = estimate_t_max(plan, view_m);
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

TMax t_max(const KernelPlan& plan, std::span<const Amount> view_m) {
  TMax out;
  if (plan.members.empty()) return out;
  const Amount hi = total_bound(plan, view_m);
  auto& probe = probe_scratch(plan.members.size());
  // alloc(hi) fits only when every member can take it; else search below hi (the search's
  // probes are counted, this one is not, like the target probe of a Lossless step).
  if (allocate(plan, hi, probe, view_m)) {
    out.value = hi;
    return out;
  }
  out.value = largest_fit(plan, view_m, 0, hi, probe, nullptr, out.iterations);
  return out;
}

LosslessStep lossless_step(const Graph& g, const KernelPlan& plan, Amount total,
                           std::span<const Amount> view_m, std::span<const Amount> own_m,
                           std::optional<Amount> request) {
  LosslessStep out;
  if (plan.placement == FlowPlacement::MaxFlow) {
    // A fresh solve on the replacement view (own_m is part of view_m, so the flow's own
    // allocation is cancellable; other owners are not in the view).
    const std::optional<Amount> cap = request.has_value() ? std::optional<Amount>(total + *request) : std::nullopt;
    auto r = solve_subset(g, plan.members, view_m, plan.roots, plan.dsts, cap);
    if (r.placed < total) throw std::logic_error("MaxFlow kernel: the replacement solve lost committed flow");
    out.total = r.placed;
    out.alloc = std::move(r.flows);
    out.cost_distribution = std::move(r.cost_distribution);
    return out;
  }
  if (plan.members.empty()) {
    out.alloc.assign(own_m.begin(), own_m.end());
    out.total = total;
    return out;
  }
  const Amount hi = total_bound(plan, view_m);
  const Amount target = request.has_value() ? std::min(total + *request, hi) : hi;
  if (target <= total) {
    out.total = total;
    out.alloc.assign(own_m.begin(), own_m.end());
    return out;
  }
  out.alloc.assign(plan.members.size(), 0);
  if (allocate(plan, target, out.alloc, view_m)) {
    out.total = target;
    return out;
  }
  // T_max in [total, target - 1]: alloc(total) = own_m fits, alloc(target) does not. Every
  // fitting probe leaves its allocation in out.alloc, so a result above total needs no further
  // evaluation.
  auto& probe = probe_scratch(plan.members.size());
  const Amount lo = largest_fit(plan, view_m, total, target, probe, &out.alloc, out.iterations);
  out.total = lo;
  if (lo == total) out.alloc.assign(own_m.begin(), own_m.end());
  return out;
}

namespace {

// The per-root bookkeeping of an attributed forwarding (per thread): each position's arrivals as
// a list of (root, amount) entries, and a dense per-root accumulator with its touched roots.
struct Arrivals {
  std::vector<std::int32_t> head;  // per position: its first entry, -1 none
  std::vector<std::int32_t> root;  // per entry
  std::vector<Amount> amount;
  std::vector<std::int32_t> next;
  std::vector<Amount> acc;  // per root
  std::vector<std::int32_t> touched;
  std::vector<Amount> weight;  // per touched root
  std::vector<Amount> row;
  std::vector<Amount> arrive;
  std::vector<Amount> keep;

  void begin(std::size_t positions, std::size_t roots) {
    head.assign(positions, -1);
    root.clear();
    amount.clear();
    next.clear();
    acc.assign(roots, 0);
  }
  void push(std::int32_t p, std::int32_t r, Amount x) {
    root.push_back(r);
    amount.push_back(x);
    next.push_back(head[static_cast<std::size_t>(p)]);
    head[static_cast<std::size_t>(p)] = static_cast<std::int32_t>(root.size()) - 1;
  }
  // The roots arriving at position p (plus `own`, a root injecting there, when >= 0) with their
  // amounts in `acc`, in ascending root order.
  void gather(std::int32_t p, std::int32_t own) {
    touched.clear();
    for (auto k = head[static_cast<std::size_t>(p)]; k >= 0; k = next[static_cast<std::size_t>(k)]) {
      const auto r = static_cast<std::size_t>(root[static_cast<std::size_t>(k)]);
      if (acc[r] == 0) touched.push_back(static_cast<std::int32_t>(r));
      acc[r] += amount[static_cast<std::size_t>(k)];
    }
    if (own >= 0 && acc[static_cast<std::size_t>(own)] == 0) touched.push_back(own);
    std::sort(touched.begin(), touched.end());
  }
  void clear_touched() {
    for (auto r : touched) acc[static_cast<std::size_t>(r)] = 0;
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
Forwarding forward_on(const KernelPlan& plan, const Injection& in, std::span<const Amount> w_m, Attribution* attr) {
  const std::size_t H = plan.roots.size();
  const std::size_t P = plan.topo.size();
  const bool capped = !w_m.empty();
  if (!in.first_hop && in.per_root.size() != H) throw std::logic_error("forward: the injection is not root-aligned");
  // The first-hop experiment reads the view of the roots' members (a program without a member
  // injects nothing).
  if (in.first_hop && !capped && !plan.members.empty()) {
    throw std::logic_error("forward: the first-hop injection needs a view");
  }
  Forwarding out;
  out.carried.assign(plan.members.size(), 0);
  out.drops.assign(plan.members.size(), 0);
  if constexpr (kAttribute) {
    attr->offered.assign(H, 0);
    attr->source_loss.assign(H, 0);
    attr->delivered.assign(H, 0);
    attr->entries.assign(H, {});
  }
  const bool weighted = plan.placement == FlowPlacement::Weighted;
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
    const auto p = plan.root_pos[i];
    const auto ri = static_cast<std::int32_t>(i);
    if (in.first_hop) {
      if (p < 0) continue;
      own_at[static_cast<std::size_t>(p)] = ri;
      Amount x = 0;
      for (auto m = plan.out_begin[static_cast<std::size_t>(p)]; m < plan.out_end[static_cast<std::size_t>(p)]; ++m) {
        x += w_m[static_cast<std::size_t>(m)];
      }
      out.offered += x;
      out.injected += x;
      if constexpr (kAttribute) attr->offered[i] = x;
      continue;
    }
    const Amount x = in.per_root[i];
    out.offered += x;
    if constexpr (kAttribute) attr->offered[i] = x;
    if (x == 0) continue;
    if (p < 0) {
      out.source_loss += x;
      if constexpr (kAttribute) attr->source_loss[i] = x;
      continue;
    }
    out.injected += x;
    inflow[static_cast<std::size_t>(p)] += x;
    if constexpr (kAttribute) a.push(p, ri, x);
  }
  for (std::size_t p = 0; p < P; ++p) {
    const auto pi = static_cast<std::int32_t>(p);
    const std::int32_t own = in.first_hop ? own_at[p] : -1;
    const Amount amount = inflow[p];
    if (plan.terminal[p] != 0) {
      out.delivered += amount;
      if (kAttribute && amount != 0) {
        a.gather(pi, -1);
        for (auto r : a.touched) attr->delivered[static_cast<std::size_t>(r)] += a.acc[static_cast<std::size_t>(r)];
        a.clear_touched();
      }
      continue;
    }
    if (amount == 0 && own < 0) continue;
    const auto b = static_cast<std::size_t>(plan.out_begin[p]);
    const auto e = static_cast<std::size_t>(plan.out_end[p]);
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
      a.weight.resize(a.touched.size());
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
      out.carried[m] = c;
      out.drops[m] = x - c;
      const auto head = plan.head_pos[m];
      inflow[static_cast<std::size_t>(head)] += c;
      if constexpr (!kAttribute) continue;
      // The member's amount by highest averages over the remaining inflows (then the root's own
      // first-hop injection on it), and what it carries by highest averages over those arrivals.
      const std::size_t t = a.touched.size();
      for (std::size_t j = 0; j < t; ++j) a.weight[j] = a.acc[static_cast<std::size_t>(a.touched[j])];
      highest_averages(part[i], a.weight, a.row);
      for (std::size_t j = 0; j < t; ++j) {
        const auto r = a.touched[j];
        a.acc[static_cast<std::size_t>(r)] -= a.row[j];
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
        if (a.acc[static_cast<std::size_t>(r)] != 0) {
          throw std::logic_error("forward: a router's inflow is not fully attributed");
        }
      }
    }
  }
  if (in.first_hop) out.offered = out.injected;
  if constexpr (kAttribute) {
    for (auto& rows : attr->entries) {
      std::sort(rows.begin(), rows.end(), [](const auto& x, const auto& y) { return x.member < y.member; });
    }
  }
  return out;
}

} // namespace

Forwarding forward(const KernelPlan& plan, const Injection& in, std::span<const Amount> w_m, Attribution* attr) {
  return attr != nullptr ? forward_on<true>(plan, in, w_m, attr) : forward_on<false>(plan, in, w_m, attr);
}

} // namespace netgraph::core::detail
