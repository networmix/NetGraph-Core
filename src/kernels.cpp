/*
  Apportionment and the placement kernels in units (APPROACH 2.1, 7 and 16).

  Every amount is an int64 number of units below 2^62: a flow's T is bounded by its source
  members' view, every inflow by T, every share by its node's inflow. Highest averages compares
  quotients by exact 128-bit cross products (weights and divisors below 2^63, products below
  2^126); lower quotas floor(I * w / W) are 128-bit products divided by W < 2^62.
*/
#include "netgraph/core/detail/kernels.hpp"

#include <algorithm>
#include <bit>
#include <limits>
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
  std::uint64_t wsum = 0;
  for (Amount w : weights) wsum += static_cast<std::uint64_t>(w);
  Amount given = 0;
  for (std::size_t m = 0; m < k; ++m) {
    // Lower quota floor(total * w / W): exactly the quotients w / j >= W / total.
    out[m] = static_cast<Amount>(div_u128(mul_u64(static_cast<std::uint64_t>(total),
                                                  static_cast<std::uint64_t>(weights[m])),
                                          wsum));
    given += out[m];
  }
  Amount left = total - given;  // fewer than k
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

void check_kernel(FlowPlacement placement, Admission admission, EdgePriority priority) {
  check_member(placement, "placement");
  check_member(admission, "admission");
  check_member(priority, "priority");
  if (admission == Admission::Lossy && placement == FlowPlacement::MaxFlow) {
    throw std::invalid_argument(
        "Lossy admission is rejected with the MaxFlow kernel: a controller never drops (EqualBalanced or "
        "Weighted only)");
  }
  if (placement == FlowPlacement::Weighted && priority == EdgePriority::Headroom) {
    throw std::invalid_argument("the Weighted kernel is rejected with priority Headroom (weights are units)");
  }
}

namespace {

// Per-thread per-node state of build_split, indexed by NodeId and valid for a node only
// when its stamp is the current epoch (no clearing between calls; a call that throws leaves stale
// entries behind an old epoch). One array of slots, so growing it either succeeds whole or leaves
// the previous array in place (review-1b-fixes N1: separately grown arrays could be left at
// different sizes by a failed allocation).
struct DagScratch {
  struct Node {
    std::uint32_t stamp {0};
    std::uint8_t state {0};   // backward DFS: 1 open, 2 done (reaches a destination)
    std::uint8_t fwd {0};     // reached from the source
    std::uint8_t target {0};  // a destination
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

// The split structure of a flow from src to the destination set over a DAG, on the DAG's own
// predecessor lists (no local renumbering, no arc arrays, no CSR to build), restricted to the
// edges `keep(e)` accepts: the kept edges on some walk from src to a destination (destinations
// terminal) become plan.members, with the plan's topological fields. One pass: a backward DFS
// from the destinations over the predecessor entries lists the nodes that reach a destination in
// postorder (every node after its parents, a topological order); a forward pass over that order
// marks the nodes the source reaches. Members ascend by EdgeId and B' ids rank edges by source
// first, so the members of one tail are consecutive. Any topological order gives the same
// allocation (a node splits once its inflow is complete).
//
// With `cost`, the forward pass also computes the minimum walk cost from the source to a
// destination over the members (kUnreachable without members): the DP runs in the same
// topological order. Without `out` (route_cost) the pass ends with the DP: no plan is built.
template <class Keep>
void build_split(const Graph& g, const PredDAG& dag, NodeId src, std::span<const NodeId> dsts, Keep keep,
                 KernelPlan* out, Cost* cost) {
  if (cost != nullptr) *cost = kUnreachable;
  const auto parents = dag.parents();
  const auto via = dag.via_edges();
  DagScratch& s = dag_scratch();
  s.begin(static_cast<std::size_t>(g.num_nodes()));
  // The DAG's sparse storage scattered into the per-node slots (O(member nodes)); every other
  // node is touched with no entries.
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
  if (!s.seen(src)) s.touch(src);
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
  if (s.at(src).state == 2) {
    s.at(src).fwd = 1;
    s.at(src).dist = 0;
  }
  for (NodeId v : s.order) {
    auto& nv = s.at(v);
    if (v != src) {
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
  plan.topo.clear();
  plan.terminal.clear();
  plan.out_begin.clear();
  plan.out_end.clear();
  plan.head_pos.clear();
  plan.src_out_begin = 0;
  plan.src_out_end = 0;
  // Members, ascending: the DAG's edges kept, from a reached tail (not a destination) to a head
  // that reaches a destination.
  const auto edge_src = g.edge_src();
  const auto edge_dst = g.edge_dst();
  // The member and node loops write by index into vectors sized to a known bound (the DAG's edges
  // here, the reached nodes below): libc++ 22's emplace_back (a fast-path and a slow-path lambda)
  // made the push_back form about ten percent slower on the kernel rows under Xcode 27.
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
  if (num_members == 0) return;
  std::int32_t num_pos = 0;
  for (NodeId v : s.order) {
    if (s.at(v).fwd != 0) s.at(v).position = num_pos++;
  }
  if (s.at(src).position != 0) {
    throw std::logic_error("kernel plan: the split graph does not start at its source");
  }
  for (std::size_t m = 0; m < num_members; ++m) {
    auto& nt = s.at(s.member_tail[m]);
    if (nt.count == 0) {
      nt.first = static_cast<std::int32_t>(m);
    } else if (static_cast<std::size_t>(nt.first + nt.count) != m) {
      throw std::logic_error("kernel plan: the members of a node are not consecutive (EdgeIds not in B' order)");
    }
    ++nt.count;
  }
  const auto P = static_cast<std::size_t>(num_pos);
  plan.topo.resize(P);
  plan.terminal.resize(P);
  plan.out_begin.resize(P);
  plan.out_end.resize(P);
  std::size_t pos = 0;
  for (NodeId v : s.order) {
    const auto& nv = s.at(v);
    if (nv.fwd == 0) continue;
    plan.topo[pos] = v;
    plan.terminal[pos] = nv.target;
    const std::int32_t b = nv.count == 0 ? 0 : nv.first;
    plan.out_begin[pos] = b;
    plan.out_end[pos] = b + nv.count;
    ++pos;
  }
  plan.head_pos.resize(num_members);
  for (std::size_t m = 0; m < num_members; ++m) {
    plan.head_pos[m] = s.at(s.member_head[m]).position;
  }
  plan.src_out_begin = plan.out_begin[0];
  plan.src_out_end = plan.out_end[0];
}

} // namespace

std::vector<EdgeId> prune(const Graph& g, const PredDAG& dag, NodeId src, std::span<const NodeId> dsts) {
  KernelPlan plan;
  build_split(g, dag, src, dsts, [](EdgeId) { return true; }, &plan, nullptr);
  return std::move(plan.members);
}

Cost route_cost(const Graph& g, const PredDAG& route, NodeId src, std::span<const NodeId> dsts) {
  Cost cost = kUnreachable;
  build_split(g, route, src, dsts, [](EdgeId) { return true; }, nullptr, &cost);
  return cost;
}

KernelPlan make_plan(const Graph& g, FlowPlacement placement, Admission admission, EdgePriority priority,
                     NodeId src, std::span<const NodeId> dsts, const PredDAG& route,
                     std::span<const Amount> creation_view, Cost* cost) {
  check_kernel(placement, admission, priority);
  const auto edges = PredDAGBuilder::ascending(route);
  const bool by_view = placement == FlowPlacement::Weighted && priority == EdgePriority::Residual;
  if (by_view && creation_view.size() != edges.size()) {
    throw std::logic_error("make_plan: the creation view is not route-aligned");
  }
  const auto capacity = g.capacity_units();
  KernelPlan plan;
  plan.placement = placement;
  plan.admission = admission;
  plan.src = src;
  plan.dsts.assign(dsts.begin(), dsts.end());
  plan.route.assign(edges.begin(), edges.end());
  // The creation view of a route edge (the route ascends: a binary search).
  auto creation = [&](EdgeId e) {
    return creation_view[static_cast<std::size_t>(std::lower_bound(edges.begin(), edges.end(), e) - edges.begin())];
  };
  // Membership: MaxFlow every route edge (its view decides), Weighted Residual creation view
  // >= 1 unit, every other kernel static capacity >= 1 unit.
  auto member = [&](EdgeId e) {
    if (by_view) return creation(e) >= 1;
    return placement == FlowPlacement::MaxFlow || capacity[static_cast<std::size_t>(e)] >= 1;
  };
  build_split(g, route, src, dsts, member, &plan, cost);
  // The route's cost is the members' only when every route edge is a member.
  if (cost != nullptr && plan.members.size() != plan.route.size()) *cost = route_cost(g, route, src, dsts);
  if (placement == FlowPlacement::Weighted) {
    plan.weights.reserve(plan.members.size());
    for (EdgeId e : plan.members) plan.weights.push_back(by_view ? creation(e) : capacity[static_cast<std::size_t>(e)]);
  }
  return plan;
}

KernelPlan make_route_plan(const Graph& g, FlowPlacement placement, Admission admission,
                           EdgePriority priority, NodeId src, std::span<const NodeId> dsts, const PredDAG& bundle,
                           const Eligibility* elig, Cost* cost) {
  check_kernel(placement, admission, priority);
  const auto capacity = g.capacity_units();
  KernelPlan plan;
  plan.placement = placement;
  plan.admission = admission;
  plan.src = src;
  plan.dsts.assign(dsts.begin(), dsts.end());
  auto usable = [&](EdgeId e) {
    return capacity[static_cast<std::size_t>(e)] >= 1 && (elig == nullptr || elig->edge_ok(e));
  };
  build_split(g, bundle, src, dsts, usable, &plan, cost);
  plan.route = plan.members;
  if (placement == FlowPlacement::Weighted && priority == EdgePriority::Capacity) {
    plan.weights.reserve(plan.members.size());
    for (EdgeId e : plan.members) plan.weights.push_back(capacity[static_cast<std::size_t>(e)]);
  }
  return plan;
}

// ---------------------------------------------------------------- allocation

bool allocate(const KernelPlan& plan, Amount total, std::span<Amount> out, std::span<const Amount> limit) {
  std::fill(out.begin(), out.end(), 0);
  if (total == 0 || plan.members.empty()) {
    if (total != 0) throw std::logic_error("allocate: an empty split graph carries nothing");
    return true;
  }
  const bool weighted = plan.placement == FlowPlacement::Weighted;
  const bool limited = !limit.empty();
  thread_local std::vector<Amount> inflow;
  inflow.assign(plan.topo.size(), 0);
  inflow[0] = total;
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

// A floating-point estimate of T_max: the fraction f of the source's inflow that each member
// carries under the plan's split (1/k per round-robin member, w/W per Weighted member), then
// min over members of view / f. It only chooses the first probe of the exact search below
// (APPROACH 7, revision 26: approximate arithmetic may choose probes, nothing else).
double estimate_t_max(const KernelPlan& plan, std::span<const Amount> view_m) {
  thread_local std::vector<double> f;
  f.assign(plan.topo.size(), 0.0);
  f[0] = 1.0;
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

// The largest T in [lo, hi) whose allocation fits view_m, given that alloc(lo) fits and alloc(hi)
// does not (APPROACH 7 and 16.3, revision 26). fits(T) is exact (integer allocation compared
// member by member with the view) and monotone in T (the allocation is non-decreasing member by
// member: house monotonicity), so the result is exact whatever probes are taken. The probes: one
// at the floating-point estimate (the only approximate arithmetic, which chooses that probe and
// nothing else), then a galloping bracket, then bisection; both bracket ends are always exact
// probe results (or the given ends). The number of probes stays within kProbeBudget = 62, the
// bound of a plain bisection: a seed or galloping probe is taken only while the probes made so
// far, that probe and a worst-case bisection of the current bracket fit in the budget; the bound
// is checked on return. With `best`, every fitting probe's allocation is kept there (swapped in),
// so on return it holds alloc(result) whenever result > the given lo.
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
  Amount hi = 0;
  for (auto k = plan.src_out_begin; k < plan.src_out_end; ++k) hi += view_m[static_cast<std::size_t>(k)];
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
  if (plan.admission != Admission::Lossless) throw std::logic_error("lossless_step on a Lossy plan");
  LosslessStep out;
  if (plan.placement == FlowPlacement::MaxFlow) {
    // A fresh solve on the replacement view (own_m is part of view_m, so the flow's own
    // allocation is cancellable; other owners are not in the view).
    const std::optional<Amount> cap = request.has_value() ? std::optional<Amount>(total + *request) : std::nullopt;
    const NodeId src[1] = {plan.src};
    auto r = solve_subset(g, plan.members, view_m, src, plan.dsts, cap);
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
  Amount hi = 0;
  for (auto k = plan.src_out_begin; k < plan.src_out_end; ++k) hi += view_m[static_cast<std::size_t>(k)];
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

LossyStep lossy_step(const KernelPlan& plan, std::optional<Amount> offer, std::span<const Amount> w_m) {
  if (plan.admission != Admission::Lossy) throw std::logic_error("lossy_step on a Lossless plan");
  LossyStep out;
  out.carried.assign(plan.members.size(), 0);
  out.drops.assign(plan.members.size(), 0);
  if (plan.members.empty() || plan.src_out_begin == plan.src_out_end) {
    out.offered = offer.value_or(0);
    return out;
  }
  const bool weighted = plan.placement == FlowPlacement::Weighted;
  thread_local std::vector<Amount> inflow;
  thread_local std::vector<Amount> part;
  inflow.assign(plan.topo.size(), 0);
  for (std::size_t p = 0; p < plan.topo.size(); ++p) {
    if (plan.terminal[p] != 0) continue;
    const auto b = static_cast<std::size_t>(plan.out_begin[p]);
    const auto e = static_cast<std::size_t>(plan.out_end[p]);
    const std::size_t k = e - b;
    part.resize(k);
    Amount amount = p == 0 ? offer.value_or(0) : inflow[p];
    std::size_t n = k;  // the members with a share (round robin leaves the rest at 0, unwritten)
    if (p == 0 && !offer.has_value()) {
      for (std::size_t i = 0; i < k; ++i) part[i] = w_m[b + i];
    } else {
      if (amount == 0) continue;
      if (weighted) {
        highest_averages(amount, std::span<const Amount>(plan.weights.data() + b, k), std::span<Amount>(part.data(), k));
      } else {
        n = round_robin(amount, k, [&](std::size_t i, Amount x) {
          part[i] = x;
          return true;
        });
      }
    }
    for (std::size_t i = 0; i < n; ++i) {
      const Amount x = part[i];
      if (p == 0) out.injected += x;
      if (x == 0) continue;
      const auto m = b + i;
      const Amount c = std::min(x, w_m[m]);
      out.carried[m] = c;
      out.drops[m] = x - c;
      inflow[static_cast<std::size_t>(plan.head_pos[m])] += c;
    }
  }
  for (std::size_t p = 0; p < plan.topo.size(); ++p) {
    if (plan.terminal[p] != 0) out.delivered += inflow[p];
  }
  out.offered = offer.has_value() ? *offer : out.injected;
  return out;
}

} // namespace netgraph::core::detail
