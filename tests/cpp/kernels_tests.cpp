// Apportionment and the placement kernels in units (APPROACH 2.1, 7, 16): the closed-form
// highest averages against the sequential definition, house monotonicity and the Hamilton hole
// control, plans (membership by kernel and priority, split-graph pruning), alloc(T) at every
// split, exact T_max, the Lossless steps of MaxFlow (replacement on V = W + own, the
// cancellation DAG, the equal-optimum and cancel-upstream controls, the release witness, the
// mixed-cost bundle), Weighted and EqualBalanced, and the Lossy steps (drops, +inf injection).
// Expected values are the closed forms and the prototype-1 fixtures named in each test.
#include <gtest/gtest.h>

#include <algorithm>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>

#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using namespace netgraph::core::detail;

namespace {

std::vector<Amount> rr(Amount total, std::size_t k) {
  std::vector<Amount> out(k);
  round_robin(total, out);
  return out;
}

std::vector<Amount> ha(Amount total, std::vector<Amount> w) {
  std::vector<Amount> out(w.size());
  highest_averages(total, w, out);
  return out;
}

// The sequential D'Hondt definition: each unit to the member maximising w / (a + 1), ties to
// the lower index.
std::vector<Amount> ha_slow(Amount total, const std::vector<Amount>& w) {
  std::vector<Amount> a(w.size(), 0);
  for (Amount u = 0; u < total; ++u) {
    std::size_t best = 0;
    for (std::size_t i = 1; i < w.size(); ++i) {
      if (static_cast<__int128>(w[i]) * (a[best] + 1) > static_cast<__int128>(w[best]) * (a[i] + 1)) best = i;
    }
    ++a[best];
  }
  return a;
}

// Hamilton (largest remainders, rejected): a control only.
std::vector<Amount> hamilton(Amount total, const std::vector<Amount>& w) {
  const Amount sum = std::accumulate(w.begin(), w.end(), Amount {0});
  std::vector<Amount> a(w.size());
  std::vector<std::size_t> order(w.size());
  for (std::size_t i = 0; i < w.size(); ++i) {
    a[i] = total * w[i] / sum;
    order[i] = i;
  }
  std::stable_sort(order.begin(), order.end(),
                   [&](std::size_t x, std::size_t y) { return (total * w[x]) % sum > (total * w[y]) % sum; });
  Amount left = total - std::accumulate(a.begin(), a.end(), Amount {0});
  for (std::size_t i = 0; i < order.size() && left > 0; ++i, --left) ++a[order[i]];
  return a;
}

// Member-aligned view of a dense per-edge array.
std::vector<Amount> on_members(const KernelPlan& plan, const std::vector<Amount>& dense) {
  std::vector<Amount> out;
  for (EdgeId e : plan.members) out.push_back(dense[static_cast<std::size_t>(e)]);
  return out;
}

// A member-aligned allocation as a dense per-edge array in input coordinates.
std::vector<Amount> dense_input(const Graph& g, const KernelPlan& plan, const std::vector<Amount>& alloc_m) {
  std::vector<Amount> by_id(static_cast<std::size_t>(g.num_edges()), 0);
  for (std::size_t m = 0; m < plan.members.size(); ++m) by_id[static_cast<std::size_t>(plan.members[m])] = alloc_m[m];
  return to_input(g, by_id);
}

std::vector<EdgeId> all_ids(const Graph& g) {
  std::vector<EdgeId> out(static_cast<std::size_t>(g.num_edges()));
  std::iota(out.begin(), out.end(), 0);
  return out;
}

KernelPlan plan_of(const Graph& g, FlowPlacement p, NodeId src, std::vector<NodeId> dsts,
                   std::vector<EdgeId> route, const std::vector<Amount>& view,
                   EdgePriority prio = EdgePriority::Residual, Admission adm = Admission::Lossless) {
  std::sort(route.begin(), route.end());
  std::vector<Amount> on_route;  // the creation view is route-aligned
  for (EdgeId e : route) on_route.push_back(view[static_cast<std::size_t>(e)]);
  return make_plan(g, p, adm, prio, src, dsts, PredDAGBuilder::of_edges(g, route), on_route);
}

// The route cost of a plan's route (its record's bundle).
Cost cost_of(const Graph& g, const KernelPlan& plan) {
  return route_cost(g, PredDAGBuilder::of_edges(g, plan.route), plan.src, plan.dsts);
}

std::vector<Amount> caps(const Graph& g) {
  return std::vector<Amount>(g.capacity_units().begin(), g.capacity_units().end());
}

// One Lossless step on a fresh flow at T = 0.
LosslessStep first_step(const Graph& g, const KernelPlan& plan, const std::vector<Amount>& view,
                        std::optional<Amount> request) {
  const auto v = on_members(plan, view);
  const std::vector<Amount> own(plan.members.size(), 0);
  return lossless_step(g, plan, 0, v, own, request);
}

} // namespace

TEST(Apportion, RoundRobinAndReviewValues) {
  EXPECT_EQ(rr(7, 3), (std::vector<Amount>{3, 2, 2}));
  EXPECT_EQ(rr(1, 3), (std::vector<Amount>{1, 0, 0}));
  EXPECT_EQ(rr(0, 4), (std::vector<Amount>{0, 0, 0, 0}));
  const Amount big = kDomainLimit - 1;
  EXPECT_EQ(rr(big, 3), (std::vector<Amount>{big / 3, big / 3, big / 3}));
  EXPECT_EQ(rr(big - 1, 3), (std::vector<Amount>{big / 3, big / 3, big / 3 - 1}));
  EXPECT_EQ(ha(10, {3, 9, 9}), (std::vector<Amount>{1, 5, 4}));
  EXPECT_EQ(ha(11, {3, 9, 9}), (std::vector<Amount>{1, 5, 5}));
  EXPECT_EQ(ha(4, {1, 3, 3}), (std::vector<Amount>{0, 2, 2}));
  EXPECT_EQ(ha(7, {1, 3, 3}), (std::vector<Amount>{1, 3, 3}));
}

// The callback form hands out the span form's shares (positive ones only, in index order) and
// stops at the first share its callback refuses.
TEST(Apportion, RoundRobinCallbackFormEqualsSpanFormAndStopsEarly) {
  for (Amount total : {Amount {0}, Amount {1}, Amount {2}, Amount {7}, Amount {23}, kDomainLimit - 1}) {
    for (std::size_t k : {std::size_t {1}, std::size_t {3}, std::size_t {4}, std::size_t {8}}) {
      std::vector<Amount> span_out(k);
      const std::size_t n = round_robin(total, span_out);
      std::vector<Amount> cb_out(k);
      std::vector<std::size_t> visited;
      const std::size_t m = round_robin(total, k, [&](std::size_t i, Amount x) {
        visited.push_back(i);
        cb_out[i] = x;
        return x > 0;
      });
      std::vector<std::size_t> leading(total < static_cast<Amount>(k) ? static_cast<std::size_t>(total) : k);
      std::iota(leading.begin(), leading.end(), std::size_t {0});
      EXPECT_EQ(n, leading.size()) << total << " over " << k;
      EXPECT_EQ(m, n) << total << " over " << k;
      EXPECT_EQ(visited, leading) << total << " over " << k;
      EXPECT_EQ(cb_out, span_out) << total << " over " << k;
    }
  }
  // 7 over 3 is {3, 2, 2}: refusing the share at index 1 returns 1 and visits nothing after it.
  std::vector<std::size_t> seen;
  const std::size_t stop = round_robin(7, 3, [&](std::size_t i, Amount x) {
    seen.push_back(i);
    return x > 2;
  });
  EXPECT_EQ(stop, 1u);
  EXPECT_EQ(seen, (std::vector<std::size_t>{0, 1}));
  EXPECT_EQ(round_robin(7, 3, [](std::size_t, Amount) { return false; }), 0u);
}

TEST(Apportion, ClosedFormEqualsSequentialAndIsHouseMonotone) {
  std::mt19937 rng(7);
  int cases = 0;
  for (int k = 1; k <= 20; ++k) {
    for (int trial = 0; trial < 40; ++trial) {
      std::vector<Amount> w(static_cast<std::size_t>(k));
      for (auto& x : w) x = 1 + static_cast<Amount>(rng() % (trial < 20 ? 5 : 1000));
      std::vector<Amount> prev(w.size(), 0);
      for (Amount total = 0; total <= 60; ++total) {
        const auto a = ha(total, w);
        ASSERT_EQ(a, ha_slow(total, w)) << "k " << k << " total " << total;
        for (std::size_t i = 0; i < a.size(); ++i) ASSERT_GE(a[i], prev[i]);
        prev = a;
        ++cases;
      }
    }
  }
  // Equal weights reduce to round robin.
  for (std::size_t k = 1; k <= 9; ++k) {
    for (Amount total = 0; total <= 40; ++total) EXPECT_EQ(ha(total, std::vector<Amount>(k, 7)), rr(total, k));
  }
  EXPECT_GT(cases, 40000);
}

TEST(Apportion, HighestAveragesNearTheDomainBound) {
  const std::vector<Amount> w = {(Amount {1} << 60) + 1, (Amount {1} << 60) - 1, 7};
  const Amount total = kDomainLimit - 2;
  const auto a = ha(total, w);
  EXPECT_EQ(std::accumulate(a.begin(), a.end(), Amount {0}), total);
  // Boundary property: every selected quotient w_i / a_i outranks every unselected one
  // w_j / (a_j + 1) (ties by the lower index).
  for (std::size_t i = 0; i < w.size(); ++i) {
    if (a[i] == 0) continue;
    for (std::size_t j = 0; j < w.size(); ++j) {
      const __int128 left = static_cast<__int128>(w[i]) * (a[j] + 1);
      const __int128 right = static_cast<__int128>(w[j]) * a[i];
      EXPECT_TRUE(left > right || (left == right && i <= j)) << i << " " << j;
    }
  }
  // Many members at large totals (the heap path) against the boundary property.
  std::vector<Amount> wide(40);
  for (std::size_t i = 0; i < wide.size(); ++i) wide[i] = (Amount {1} << 55) + static_cast<Amount>(i * 7919);
  const auto b = ha((Amount {1} << 61) + 12345, wide);
  EXPECT_EQ(std::accumulate(b.begin(), b.end(), Amount {0}), (Amount {1} << 61) + 12345);
  for (std::size_t i = 0; i < wide.size(); ++i) {
    for (std::size_t j = 0; j < wide.size(); ++j) {
      if (b[i] == 0) continue;
      const __int128 left = static_cast<__int128>(wide[i]) * (b[j] + 1);
      const __int128 right = static_cast<__int128>(wide[j]) * b[i];
      EXPECT_TRUE(left > right || (left == right && i <= j));
    }
  }
}

TEST(Apportion, HamiltonHoleControl) {
  // Weights [1, 3, 3] on members with view [0, 3, 3]: highest averages is feasible on 0..4
  // (T_max = 4 exactly); Hamilton's rule has the hole at 3 ({0, 1, 2, 4}).
  const std::vector<Amount> w = {1, 3, 3};
  const std::vector<Amount> view = {0, 3, 3};
  auto fits = [&](const std::vector<Amount>& a) {
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (a[i] > view[i]) return false;
    }
    return true;
  };
  std::vector<Amount> ha_ok, ham_ok;
  for (Amount t = 0; t < 8; ++t) {
    if (fits(ha(t, w))) ha_ok.push_back(t);
    if (fits(hamilton(t, w))) ham_ok.push_back(t);
  }
  EXPECT_EQ(ha_ok, (std::vector<Amount>{0, 1, 2, 3, 4}));
  EXPECT_EQ(ham_ok, (std::vector<Amount>{0, 1, 2, 4}));
  EXPECT_EQ(hamilton(3, w), (std::vector<Amount>{1, 1, 1}));
  EXPECT_EQ(hamilton(4, w), (std::vector<Amount>{0, 2, 2}));
  // The kernel's T_max on the same plan (parallel links 1, 3, 3 with Capacity weights, link 0
  // held by a third party): exactly 4, allocation [0, 2, 2].
  auto g = make_graph(2, {{0, 1, 1, 1}, {0, 1, 3, 1}, {0, 1, 3, 1}});
  const auto plan = plan_of(g, FlowPlacement::Weighted, 0, {1}, all_ids(g), caps(g), EdgePriority::Capacity);
  const auto view_dense = to_ids(g, std::vector<Amount> {0, 3, 3});
  const auto tm = t_max(plan, on_members(plan, view_dense));
  EXPECT_EQ(tm.value, 4);
  std::vector<Amount> at4(plan.members.size());
  ASSERT_TRUE(allocate(plan, 4, at4));
  EXPECT_EQ(dense_input(g, plan, at4), (std::vector<Amount>{0, 2, 2}));
  std::vector<Amount> at3(plan.members.size());
  ASSERT_TRUE(allocate(plan, 3, at3));
  EXPECT_EQ(dense_input(g, plan, at3), (std::vector<Amount>{0, 2, 1}));
}

TEST(Plans, MembershipByKernelAndPriority) {
  // Parallel links of capacities [4, 4, 0] with residual [0, 4, 0] on the creation view.
  auto g = make_graph(2, {{0, 1, 4, 1}, {0, 1, 4, 1}, {0, 1, 0, 1}});
  const auto view = to_ids(g, std::vector<Amount> {0, 4, 0});
  const auto inv = id_of_input(g);
  auto mf = plan_of(g, FlowPlacement::MaxFlow, 0, {1}, all_ids(g), view);
  EXPECT_EQ(mf.members, (std::vector<EdgeId>{0, 1, 2}));  // MaxFlow: the view decides
  auto eb = plan_of(g, FlowPlacement::EqualBalanced, 0, {1}, all_ids(g), view);
  EXPECT_EQ(eb.members, (std::vector<EdgeId>{inv[0], inv[1]}));  // static capacity >= 1
  auto wr = plan_of(g, FlowPlacement::Weighted, 0, {1}, all_ids(g), view, EdgePriority::Residual);
  EXPECT_EQ(wr.members, (std::vector<EdgeId>{inv[1]}));
  EXPECT_EQ(wr.weights, (std::vector<Amount>{4}));
  auto wc = plan_of(g, FlowPlacement::Weighted, 0, {1}, all_ids(g), view, EdgePriority::Capacity);
  EXPECT_EQ(wc.members, (std::vector<EdgeId>{inv[0], inv[1]}));
  EXPECT_EQ(wc.weights, (std::vector<Amount>{4, 4}));
  EXPECT_THROW((void)plan_of(g, FlowPlacement::Weighted, 0, {1}, all_ids(g), view, EdgePriority::Headroom),
               std::invalid_argument);
  EXPECT_THROW(
      (void)plan_of(g, FlowPlacement::MaxFlow, 0, {1}, all_ids(g), view, EdgePriority::Residual, Admission::Lossy),
      std::invalid_argument);
}

TEST(Plans, SplitGraphPruning) {
  // E1a: EqualBalanced with a zero-capacity tail. 0->1, 0->2, 1->3, 2->4 (cap 0), 4->3: the
  // branch through 2 cannot reach 3 over members, so it is pruned (no stranded unit).
  auto g = make_graph(5, {{0, 1, 2, 1}, {0, 2, 2, 1}, {1, 3, 2, 1}, {2, 4, 0, 1}, {4, 3, 2, 1}});
  const auto inv = id_of_input(g);
  auto eb = plan_of(g, FlowPlacement::EqualBalanced, 0, {3}, all_ids(g), caps(g));
  std::vector<EdgeId> want = {inv[0], inv[2]};
  std::sort(want.begin(), want.end());
  EXPECT_EQ(eb.members, want);
  std::vector<Amount> a(eb.members.size());
  ASSERT_TRUE(allocate(eb, 2, a));
  EXPECT_EQ(dense_input(g, eb, a), (std::vector<Amount>{2, 0, 2, 0, 0}));
  // Destinations are terminal; a destination set keeps both destinations.
  auto g2 = make_graph(4, {{0, 1, 3, 1}, {1, 2, 3, 1}, {0, 2, 3, 2}, {2, 3, 3, 1}});
  auto any = plan_of(g2, FlowPlacement::EqualBalanced, 0, {1, 2}, all_ids(g2), caps(g2));
  const auto inv2 = id_of_input(g2);
  std::vector<EdgeId> want2 = {inv2[0], inv2[2]};  // 1->2 leaves destination 1; 2->3 leaves 2
  std::sort(want2.begin(), want2.end());
  EXPECT_EQ(any.members, want2);
  EXPECT_EQ(cost_of(g2, any), 1);
  // An empty split graph: route cost unreachable, T_max 0.
  auto none = plan_of(g, FlowPlacement::EqualBalanced, 0, {3}, std::vector<EdgeId> {inv[1], inv[3]}, caps(g));
  EXPECT_TRUE(none.members.empty());
  EXPECT_EQ(cost_of(g, none), kUnreachable);
  EXPECT_EQ(t_max(none, {}).value, 0);
}

TEST(Plans, SplitGraphEqualsTheWalkDefinition) {
  // On random DAGs (parallel edges, zero capacities, statically excluded edges and nodes, edges
  // into the source and out of a destination, one or two destinations) the split graph of a
  // bundle is, by definition, the usable edges (capacity >= 1 unit, allowed) on some walk from
  // the source to a destination that enters no destination before its end: the members equal an
  // independent reachability computation, the route cost equals an independent shortest walk,
  // the topological structure orders every member's tail before its head, and an allocation
  // conserves its total. make_plan on the DAG of the usable edges gives the same plan as
  // make_route_plan filtering the bundle (the two entry points of the one builder).
  std::mt19937 rng(2024);
  int compared = 0;
  for (int trial = 0; trial < 2000; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 8);
    std::vector<InEdge> edges;
    for (int i = 0; i < 2 * n + static_cast<int>(rng() % 10); ++i) {
      const NodeId u = static_cast<NodeId>(rng() % static_cast<std::uint32_t>(n - 1));
      const NodeId v = static_cast<NodeId>(u + 1 + rng() % static_cast<std::uint32_t>(n - 1 - u));
      const Amount cap = (rng() % 6) == 0 ? 0 : 1 + static_cast<Amount>(rng() % 9);
      edges.push_back({u, v, cap, 1 + static_cast<Cost>(rng() % 4)});
    }
    auto g = make_graph(n, edges);
    std::vector<EdgeId> all(static_cast<std::size_t>(g.num_edges()));
    std::iota(all.begin(), all.end(), 0);
    std::vector<EdgeId> in_bundle;
    for (EdgeId e : all) {
      if ((rng() % 5) != 0) in_bundle.push_back(e);
    }
    const auto bundle = PredDAGBuilder::of_edges(g, in_bundle);
    const NodeId src = static_cast<NodeId>(rng() % static_cast<std::uint32_t>(n));
    std::vector<NodeId> dsts;
    for (int k = 0; k < 1 + static_cast<int>(rng() % 2); ++k) {
      const NodeId d = static_cast<NodeId>(rng() % static_cast<std::uint32_t>(n));
      if (d != src && std::find(dsts.begin(), dsts.end(), d) == dsts.end()) dsts.push_back(d);
    }
    if (dsts.empty()) continue;
    std::sort(dsts.begin(), dsts.end());
    std::optional<Eligibility> elig;
    if ((rng() % 3) == 0) {
      std::vector<NodeId> xn;
      std::vector<EdgeId> xe;
      if ((rng() % 2) == 0) xn.push_back(static_cast<NodeId>(rng() % static_cast<std::uint32_t>(n)));
      for (EdgeId e : all) {
        if ((rng() % 7) == 0) xe.push_back(e);
      }
      elig = g.eligibility(xn, xe);
    }
    const auto capacity = g.capacity_units();
    const auto edge_src = g.edge_src();
    const auto edge_dst = g.edge_dst();
    const auto cost = g.cost();
    std::vector<EdgeId> usable;
    for (EdgeId e : PredDAGBuilder::ascending(bundle)) {
      if (capacity[static_cast<std::size_t>(e)] >= 1 && (!elig || elig->edge_ok(e))) usable.push_back(e);
    }
    // The definition: forward from src and backward from the destinations over the usable edges
    // that leave no destination, by fixpoint; the shortest walk by Bellman-Ford over the same.
    std::vector<std::uint8_t> is_dst(static_cast<std::size_t>(n), 0);
    for (NodeId d : dsts) is_dst[static_cast<std::size_t>(d)] = 1;
    std::vector<EdgeId> open;
    for (EdgeId e : usable) {
      if (is_dst[static_cast<std::size_t>(edge_src[static_cast<std::size_t>(e)])] == 0) open.push_back(e);
    }
    std::vector<std::uint8_t> fwd(static_cast<std::size_t>(n), 0), bwd(static_cast<std::size_t>(n), 0);
    fwd[static_cast<std::size_t>(src)] = 1;
    for (NodeId d : dsts) bwd[static_cast<std::size_t>(d)] = 1;
    std::vector<Cost> dist(static_cast<std::size_t>(n), kUnreachable);
    dist[static_cast<std::size_t>(src)] = 0;
    for (int round = 0; round < n; ++round) {
      for (EdgeId e : open) {
        const auto u = static_cast<std::size_t>(edge_src[static_cast<std::size_t>(e)]);
        const auto v = static_cast<std::size_t>(edge_dst[static_cast<std::size_t>(e)]);
        if (fwd[u] != 0) fwd[v] = 1;
        if (bwd[v] != 0) bwd[u] = 1;
        if (dist[u] != kUnreachable) dist[v] = std::min(dist[v], dist[u] + cost[static_cast<std::size_t>(e)]);
      }
    }
    std::vector<EdgeId> want;
    for (EdgeId e : open) {
      if (fwd[static_cast<std::size_t>(edge_src[static_cast<std::size_t>(e)])] != 0 &&
          bwd[static_cast<std::size_t>(edge_dst[static_cast<std::size_t>(e)])] != 0) {
        want.push_back(e);
      }
    }
    Cost want_cost = kUnreachable;
    if (!want.empty()) {
      for (NodeId d : dsts) want_cost = std::min(want_cost, dist[static_cast<std::size_t>(d)]);
    }
    const auto usable_dag = PredDAGBuilder::of_edges(g, usable);
    for (const auto kernel : {FlowPlacement::EqualBalanced, FlowPlacement::Weighted}) {
      Cost got_cost = -1;
      const auto on_bundle = make_route_plan(g, kernel, Admission::Lossless, EdgePriority::Capacity, src, dsts, bundle,
                                             elig ? &*elig : nullptr, &got_cost);
      ASSERT_EQ(on_bundle.members, want) << "trial " << trial;
      EXPECT_EQ(on_bundle.route, want);
      EXPECT_EQ(got_cost, want_cost) << "trial " << trial;
      EXPECT_EQ(prune(g, usable_dag, src, dsts), want);
      Cost plan_cost = -1;
      const auto on_usable =
          make_plan(g, kernel, Admission::Lossless, EdgePriority::Capacity, src, dsts, usable_dag, {}, &plan_cost);
      EXPECT_EQ(on_usable.members, on_bundle.members);
      EXPECT_EQ(on_usable.topo, on_bundle.topo);
      EXPECT_EQ(on_usable.out_begin, on_bundle.out_begin);
      EXPECT_EQ(on_usable.out_end, on_bundle.out_end);
      EXPECT_EQ(on_usable.head_pos, on_bundle.head_pos);
      EXPECT_EQ(on_usable.terminal, on_bundle.terminal);
      EXPECT_EQ(on_usable.weights, on_bundle.weights);
      // The route cost is over the route (every usable edge), not only the members.
      dist.assign(static_cast<std::size_t>(n), kUnreachable);
      dist[static_cast<std::size_t>(src)] = 0;
      for (int round = 0; round < n; ++round) {
        for (EdgeId e : usable) {
          const auto u = static_cast<std::size_t>(edge_src[static_cast<std::size_t>(e)]);
          const auto v = static_cast<std::size_t>(edge_dst[static_cast<std::size_t>(e)]);
          if (is_dst[u] == 0 && dist[u] != kUnreachable) dist[v] = std::min(dist[v], dist[u] + cost[static_cast<std::size_t>(e)]);
        }
      }
      Cost route_want = kUnreachable;
      for (NodeId d : dsts) route_want = std::min(route_want, dist[static_cast<std::size_t>(d)]);
      EXPECT_EQ(plan_cost, route_want) << "trial " << trial;
      EXPECT_EQ(route_cost(g, usable_dag, src, dsts), route_want) << "trial " << trial;  // the cost pass alone
      if (on_bundle.members.empty()) continue;
      // Topology: the source first, every member's tail before its head, members of a node
      // consecutive.
      EXPECT_EQ(on_bundle.topo.front(), src);
      for (std::size_t p = 0; p < on_bundle.topo.size(); ++p) {
        for (auto k = on_bundle.out_begin[p]; k < on_bundle.out_end[p]; ++k) {
          EXPECT_EQ(edge_src[static_cast<std::size_t>(on_bundle.members[static_cast<std::size_t>(k)])], on_bundle.topo[p]);
          EXPECT_GT(on_bundle.head_pos[static_cast<std::size_t>(k)], static_cast<std::int32_t>(p));
        }
      }
      // Conservation: what enters the destinations is the total.
      std::vector<Amount> a(on_bundle.members.size());
      for (Amount t : {Amount {1}, Amount {7}, Amount {50}}) {
        ASSERT_TRUE(allocate(on_bundle, t, a));
        Amount delivered = 0;
        for (std::size_t m = 0; m < a.size(); ++m) {
          if (is_dst[static_cast<std::size_t>(edge_dst[static_cast<std::size_t>(on_bundle.members[m])])] != 0) {
            delivered += a[m];
          }
        }
        EXPECT_EQ(delivered, t) << "trial " << trial;
      }
      ++compared;
    }
  }
  EXPECT_GT(compared, 1000);
}

TEST(Allocation, Review14SplitAndReconvergeEveryIntermediateValue) {
  // Prototype 1 fixture apportion_house_monotone_dag: capacities 20 x weights, baseline view.
  const std::vector<std::pair<NodeId, NodeId>> dag = {{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 5}, {4, 5}};
  const std::vector<Amount> w = {2, 7, 1, 4, 3, 2, 1, 1};
  const std::vector<Amount> view_in = {8, 20, 4, 9, 10, 10, 15, 15};
  std::vector<InEdge> edges;
  for (std::size_t i = 0; i < dag.size(); ++i) edges.push_back({dag[i].first, dag[i].second, 20 * w[i], 1});
  auto g = make_graph(6, edges);
  const auto view = to_ids(g, view_in);
  for (auto kernel : {FlowPlacement::Weighted, FlowPlacement::EqualBalanced}) {
    // Independent per-node reference on input positions (node ids are topological).
    auto ref = [&](Amount total) {
      std::vector<Amount> inflow(6, 0), out(dag.size(), 0);
      inflow[0] = total;
      for (NodeId u = 0; u < 5; ++u) {
        std::vector<std::size_t> mine;
        for (std::size_t i = 0; i < dag.size(); ++i) {
          if (dag[i].first == u) mine.push_back(i);
        }
        if (mine.empty() || inflow[static_cast<std::size_t>(u)] == 0) continue;
        std::vector<Amount> ww;
        for (auto i : mine) ww.push_back(w[i]);
        const auto parts = kernel == FlowPlacement::Weighted ? ha_slow(inflow[static_cast<std::size_t>(u)], ww)
                                                              : rr(inflow[static_cast<std::size_t>(u)], mine.size());
        for (std::size_t j = 0; j < mine.size(); ++j) {
          out[mine[j]] = parts[j];
          inflow[static_cast<std::size_t>(dag[mine[j]].second)] += parts[j];
        }
      }
      return out;
    };
    Amount brute = 0;
    for (Amount t = 0; t <= 100; ++t) {
      const auto a = ref(t);
      bool ok = true;
      for (std::size_t i = 0; i < a.size(); ++i) ok = ok && a[i] <= view_in[i];
      if (ok) brute = t;
    }
    const auto plan = plan_of(g, kernel, 0, {5}, all_ids(g), caps(g), EdgePriority::Capacity);
    const auto vm = on_members(plan, view);
    std::vector<Amount> own(plan.members.size(), 0);
    Amount total = 0;
    int steps = 0;
    while (true) {
      const auto step = lossless_step(g, plan, total, vm, own, 1);
      if (step.total == total) break;
      ++steps;
      for (std::size_t m = 0; m < own.size(); ++m) ASSERT_GE(step.alloc[m], own[m]);  // member-monotone
      EXPECT_EQ(dense_input(g, plan, step.alloc), ref(step.total));
      total = step.total;
      own = step.alloc;
    }
    EXPECT_EQ(total, brute);
    EXPECT_EQ(t_max(plan, vm).value, brute);
    const auto one = first_step(g, plan, view, brute);
    EXPECT_EQ(one.alloc, own);  // one call of T_max equals the calls of 1
    EXPECT_GT(steps, 5);
  }
}

TEST(Allocation, EqualBalancedOnePlusTwoEqualsThree) {
  // Review-13 EB reconvergence DAG (prototype 1 apportion_house_monotone_dag, second part).
  const std::vector<std::pair<NodeId, NodeId>> recon = {{0, 1}, {0, 2}, {0, 3}, {1, 4}, {2, 4},
                                                        {3, 4}, {4, 5}, {4, 6}, {5, 7}, {6, 7}};
  std::vector<InEdge> edges;
  for (auto [u, v] : recon) edges.push_back({u, v, 100, 1});
  auto g = make_graph(8, edges);
  const auto plan = plan_of(g, FlowPlacement::EqualBalanced, 0, {7}, all_ids(g), caps(g));
  const auto vm = on_members(plan, caps(g));
  std::vector<Amount> own(plan.members.size(), 0);
  auto s1 = lossless_step(g, plan, 0, vm, own, 1);
  EXPECT_EQ(dense_input(g, plan, s1.alloc), (std::vector<Amount>{1, 0, 0, 1, 0, 0, 1, 0, 1, 0}));
  auto s2 = lossless_step(g, plan, s1.total, vm, s1.alloc, 2);
  EXPECT_EQ(dense_input(g, plan, s2.alloc), (std::vector<Amount>{1, 1, 1, 1, 1, 1, 2, 1, 2, 1}));
  EXPECT_EQ(first_step(g, plan, caps(g), 3).alloc, s2.alloc);
}

TEST(Allocation, TMaxExactAgainstBruteForce) {
  std::mt19937 rng(99);
  int checked = 0;
  for (int trial = 0; trial < 400; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 5);
    std::vector<InEdge> edges;
    for (int i = 0; i < 2 * n + static_cast<int>(rng() % 6); ++i) {
      const NodeId u = static_cast<NodeId>(rng() % (n - 1));
      const NodeId v = static_cast<NodeId>(u + 1 + rng() % static_cast<std::uint32_t>(n - 1 - u));
      edges.push_back({u, v, 1 + static_cast<Amount>(rng() % 9), 1});
    }
    auto g = make_graph(n, edges);
    std::vector<Amount> view(static_cast<std::size_t>(g.num_edges()));
    for (std::size_t e = 0; e < view.size(); ++e) {
      view[e] = static_cast<Amount>(rng() % static_cast<std::uint32_t>(g.capacity_units()[e] + 1));
    }
    const auto kernel = (rng() % 2) == 0 ? FlowPlacement::Weighted : FlowPlacement::EqualBalanced;
    const auto plan = plan_of(g, kernel, 0, {n - 1}, all_ids(g), caps(g), EdgePriority::Capacity);
    if (plan.members.empty()) continue;
    const auto vm = on_members(plan, view);
    const auto tm = t_max(plan, vm);
    EXPECT_LE(tm.iterations, 62);
    Amount brute = 0;
    std::vector<Amount> a(plan.members.size());
    for (Amount t = 0; t <= 200; ++t) {
      if (allocate(plan, t, a, vm)) brute = t;
    }
    ASSERT_EQ(tm.value, brute) << "trial " << trial;
    // At T_max some member is exactly full and the next unit overfills one.
    ASSERT_TRUE(allocate(plan, brute, a, vm));
    std::vector<Amount> next(plan.members.size());
    EXPECT_FALSE(allocate(plan, brute + 1, next, vm));
    ++checked;
  }
  EXPECT_GT(checked, 300);
}

TEST(Allocation, TMaxSearchIsExactAndWithinTheProbeBudgetOnLargeRanges) {
  // Capacities and views up to 2^58 units (search ranges near 2^62 at the source): the seeded
  // search must return the same T as a plain bisection from 0 and probe at most 62 times, also
  // when its floating-point estimate is far off (a member capped near zero deep in the DAG).
  std::mt19937_64 rng(1234);
  int checked = 0;
  for (int trial = 0; trial < 300; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 6);
    std::vector<InEdge> edges;
    for (int i = 0; i < 2 * n + static_cast<int>(rng() % 8); ++i) {
      const NodeId u = static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n - 1));
      const NodeId v = static_cast<NodeId>(u + 1 + rng() % static_cast<std::uint64_t>(n - 1 - u));
      edges.push_back({u, v, 1 + static_cast<Amount>(rng() % (std::uint64_t {1} << 58)), 1});
    }
    auto g = make_graph(n, edges);
    std::vector<Amount> view(static_cast<std::size_t>(g.num_edges()));
    for (std::size_t e = 0; e < view.size(); ++e) {
      const Amount c = g.capacity_units()[e];
      view[e] = (rng() % 4) == 0 ? static_cast<Amount>(rng() % 3) : static_cast<Amount>(rng() % static_cast<std::uint64_t>(c + 1));
    }
    const auto kernel = (rng() % 2) == 0 ? FlowPlacement::Weighted : FlowPlacement::EqualBalanced;
    const auto plan = plan_of(g, kernel, 0, {n - 1}, all_ids(g), caps(g), EdgePriority::Capacity);
    if (plan.members.empty()) continue;
    const auto vm = on_members(plan, view);
    const auto tm = t_max(plan, vm);
    EXPECT_LE(tm.iterations, 62);
    // Reference: plain bisection on [0, hi] with the same exact predicate.
    Amount hi = 0;
    for (auto k = plan.src_out_begin; k < plan.src_out_end; ++k) hi += vm[static_cast<std::size_t>(k)];
    std::vector<Amount> a(plan.members.size());
    Amount lo = 0;
    while (lo < hi) {
      const Amount mid = lo + (hi - lo + 1) / 2;
      if (allocate(plan, mid, a, vm)) lo = mid; else hi = mid - 1;
    }
    ASSERT_EQ(tm.value, lo) << "trial " << trial;
    // The Lossless step from T = 0 with a request beyond the maximum lands on the same T.
    const std::vector<Amount> zero(plan.members.size(), 0);
    const auto step = lossless_step(g, plan, 0, vm, zero, kDomainLimit - 1);
    EXPECT_EQ(step.total, lo);
    EXPECT_LE(step.iterations, 62);
    ++checked;
  }
  EXPECT_GT(checked, 200);
}

TEST(MaxFlowKernel, CancellationDagGivesTwo) {
  // Review-5 V3 (prototype 1 maxflow_cancellation_dag): place(1) takes 0-1-3-5, the next
  // place(1) re-solves on V = W + own and returns the unique value-2 flow (1->3 goes to 0).
  const std::vector<std::pair<NodeId, NodeId>> pairs = {{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 3}, {3, 5}, {4, 5}};
  std::vector<InEdge> edges;
  for (auto [u, v] : pairs) edges.push_back({u, v, 1, 1});
  auto g = make_graph(6, edges);
  const auto plan = plan_of(g, FlowPlacement::MaxFlow, 0, {5}, all_ids(g), caps(g));
  auto s1 = first_step(g, plan, caps(g), 1);
  EXPECT_EQ(dense_input(g, plan, s1.alloc), (std::vector<Amount>{1, 0, 1, 0, 0, 1, 0}));
  // The flow's own allocation is in its view (W + own = capacities here).
  auto s2 = lossless_step(g, plan, s1.total, on_members(plan, caps(g)), s1.alloc, 1);
  EXPECT_EQ(s2.total, 2);
  EXPECT_EQ(dense_input(g, plan, s2.alloc), (std::vector<Amount>{1, 1, 0, 1, 1, 1, 1}));
  auto s3 = lossless_step(g, plan, s2.total, on_members(plan, caps(g)), s2.alloc, 1);
  EXPECT_EQ(s3.total, 2);
  EXPECT_EQ(first_step(g, plan, caps(g), 2).alloc, s2.alloc);
  EXPECT_EQ(first_step(g, plan, caps(g), std::nullopt).total, 2);
}

namespace {

// Review-9 witness edges (prototype 1 CANCEL_EDGES): s=0; A=1, B=2, C=3; X=4, Y=5; t=6.
std::vector<InEdge> cancel_edges(Cost b_to_y = 9) {
  return {{0, 1, 1, 1}, {0, 2, 1, 1}, {0, 3, 1, 1}, {1, 4, 1, 2}, {1, 5, 1, 9},
          {2, 4, 1, 1}, {2, 5, 1, b_to_y}, {3, 4, 1, 3}, {4, 6, 2, 1}, {5, 6, 1, 1}};
}

} // namespace

TEST(MaxFlowKernel, CancelUpstreamMixedCost) {
  auto g = make_graph(7, cancel_edges());
  const auto plan = plan_of(g, FlowPlacement::MaxFlow, 0, {6}, all_ids(g), caps(g));
  const std::vector<Amount> final_alloc = {1, 1, 1, 0, 1, 1, 0, 1, 2, 1};
  auto one = first_step(g, plan, caps(g), 3);
  EXPECT_EQ(dense_input(g, plan, one.alloc), final_alloc);
  // The sequence 1, 1, 1 with the third call cancelling A->X upstream (a signed delta).
  const std::vector<std::vector<Amount>> seq = {
      {0, 1, 0, 0, 0, 1, 0, 0, 1, 0}, {1, 1, 0, 1, 0, 1, 0, 0, 2, 0}, final_alloc};
  std::vector<Amount> own(plan.members.size(), 0);
  Amount total = 0;
  for (const auto& want : seq) {
    auto s = lossless_step(g, plan, total, on_members(plan, caps(g)), own, 1);
    EXPECT_EQ(s.total, total + 1);
    EXPECT_EQ(dense_input(g, plan, s.alloc), want);
    total = s.total;
    own = s.alloc;
  }
  // The minimum cost over the bundle: curve (3, 1), (4, 1), (12, 1) totals 19.
  Cost cost = 0;
  const auto a = dense_input(g, plan, own);
  for (std::size_t i = 0; i < a.size(); ++i) cost += a[i] * cancel_edges()[i].cost;
  EXPECT_EQ(cost, 19);
}

TEST(MaxFlowKernel, EqualOptimumControlFollowsBPrimeArcOrder) {
  // Rule-20 control: with cost(B->Y) = 8 two value-3 optima cost 19; ascending B' arcs pick the
  // flow with Y fed by A, so a later demand on A->X still finds its unit.
  auto g = make_graph(7, cancel_edges(8));
  const auto plan = plan_of(g, FlowPlacement::MaxFlow, 0, {6}, all_ids(g), caps(g));
  auto s = first_step(g, plan, caps(g), 3);
  EXPECT_EQ(dense_input(g, plan, s.alloc), (std::vector<Amount>{1, 1, 1, 0, 1, 1, 0, 1, 2, 1}));
}

TEST(MaxFlowKernel, ReleaseWitnessAndMixedCostBundle) {
  // Release witness: two parallel links of cost 4 and 2; the cheap one is held by another owner
  // at the first call (the flow takes cost 4); after the release the re-solve on the new view
  // moves the unit to the cheap link (placed 0 in that call, allocation moved).
  auto g = make_graph(2, {{0, 1, 1, 4}, {0, 1, 1, 2}});
  const auto plan = plan_of(g, FlowPlacement::MaxFlow, 0, {1}, all_ids(g), caps(g));
  const auto held = to_ids(g, std::vector<Amount> {1, 0});
  auto s1 = first_step(g, plan, held, 1);
  EXPECT_EQ(dense_input(g, plan, s1.alloc), (std::vector<Amount>{1, 0}));
  const auto released = to_ids(g, std::vector<Amount> {0, 1});  // W after the release; own adds 1 on cost 4
  std::vector<Amount> v(plan.members.size());
  for (std::size_t m = 0; m < v.size(); ++m) v[m] = released[static_cast<std::size_t>(plan.members[m])] + s1.alloc[m];
  auto s2 = lossless_step(g, plan, 1, v, s1.alloc, 1);
  EXPECT_EQ(s2.total, 2);
  // Request 0 more is never sent; a capped re-solve at the same T moves the unit (cost 2).
  auto s3 = lossless_step(g, plan, 1, v, s1.alloc, 0);
  EXPECT_EQ(s3.total, 1);
  EXPECT_EQ(dense_input(g, plan, s3.alloc), (std::vector<Amount>{0, 1}));
  // Review-13 D4: on a mixed-cost bundle the 1 + 1 + 1 chain fills before the cost-4 edge.
  auto g2 = make_graph(4, {{0, 1, 1, 1}, {1, 2, 1, 1}, {2, 3, 1, 1}, {0, 3, 1, 4}});
  const auto p2 = plan_of(g2, FlowPlacement::MaxFlow, 0, {3}, all_ids(g2), caps(g2));
  EXPECT_EQ(dense_input(g2, p2, first_step(g2, p2, caps(g2), 1).alloc), (std::vector<Amount>{1, 1, 1, 0}));
  EXPECT_EQ(cost_of(g2, p2), 3);
}

TEST(WeightedKernel, CapacityVersusResidualMembershipAfterOccupancy) {
  // Prototype 1 weighted_capacity_vs_residual_after_occupancy: parallel links [4, 4], another
  // owner holds 3 on link 0 (residual [1, 4]).
  auto g = make_graph(2, {{0, 1, 4, 1}, {0, 1, 4, 1}});
  auto run = [&](EdgePriority prio, std::vector<Amount> residual_in, Amount request) {
    const auto w = to_ids(g, residual_in);
    const auto plan = plan_of(g, FlowPlacement::Weighted, 0, {1}, all_ids(g), w, prio);
    auto s = first_step(g, plan, w, request);
    return std::make_pair(s.total, dense_input(g, plan, s.alloc));
  };
  EXPECT_EQ(run(EdgePriority::Residual, {1, 4}, 5), std::make_pair(Amount {5}, std::vector<Amount>{1, 4}));
  EXPECT_EQ(run(EdgePriority::Capacity, {1, 4}, 5), std::make_pair(Amount {2}, std::vector<Amount>{1, 1}));
  EXPECT_EQ(run(EdgePriority::Residual, {0, 4}, 5), std::make_pair(Amount {4}, std::vector<Amount>{0, 4}));
  EXPECT_EQ(run(EdgePriority::Capacity, {0, 4}, 5), std::make_pair(Amount {0}, std::vector<Amount>{0, 0}));
  // Frozen plans after the owner leaves: Capacity grows to 5 as [3, 2]; Residual stays at 4.
  const auto held = to_ids(g, std::vector<Amount> {0, 4});
  const auto free = to_ids(g, std::vector<Amount> {4, 4});
  const auto pc = plan_of(g, FlowPlacement::Weighted, 0, {1}, all_ids(g), held, EdgePriority::Capacity);
  auto c = lossless_step(g, pc, 0, on_members(pc, free), std::vector<Amount>(pc.members.size(), 0), 5);
  EXPECT_EQ(dense_input(g, pc, c.alloc), (std::vector<Amount>{3, 2}));
  const auto pr = plan_of(g, FlowPlacement::Weighted, 0, {1}, all_ids(g), held, EdgePriority::Residual);
  auto first = first_step(g, pr, held, 5);
  // After the release V = W + own on the one member is its capacity 4.
  auto again = lossless_step(g, pr, first.total, on_members(pr, free), first.alloc, 5);
  EXPECT_EQ(again.total, 4);  // the plan has one member; weights never renormalised
}

TEST(WeightedKernel, OnePlusFiveEqualsSix) {
  // Additivity on a fixed plan (prototype 1 weighted_1_plus_5_vs_6): weights 2 and 6 at the
  // source.
  auto g = make_graph(3, {{0, 1, 2, 1}, {0, 2, 6, 1}, {1, 2, 9, 1}});
  const auto plan = plan_of(g, FlowPlacement::Weighted, 0, {2}, all_ids(g), caps(g), EdgePriority::Capacity);
  const auto vm = on_members(plan, caps(g));
  std::vector<Amount> zero(plan.members.size(), 0);
  auto a = lossless_step(g, plan, 0, vm, zero, 1);
  auto b = lossless_step(g, plan, a.total, vm, a.alloc, 5);
  auto c = lossless_step(g, plan, 0, vm, zero, 6);
  EXPECT_EQ(b.total, c.total);
  EXPECT_EQ(b.alloc, c.alloc);
}

TEST(LossyKernel, DiamondDropsAndInfiniteInjection) {
  // Prototype 1 lossy_diamond: source links 10/10, tails 1/10.
  auto g = make_graph(4, {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 1, 1}, {2, 3, 10, 1}});
  for (auto kernel : {FlowPlacement::EqualBalanced, FlowPlacement::Weighted}) {
    const auto plan = plan_of(g, kernel, 0, {3}, all_ids(g), caps(g), EdgePriority::Residual, Admission::Lossy);
    auto w = on_members(plan, caps(g));
    auto s = lossy_step(plan, 10, w);
    EXPECT_EQ(s.offered, 10);
    EXPECT_EQ(s.injected, 10);
    EXPECT_EQ(s.delivered, 6);
    EXPECT_EQ(dense_input(g, plan, s.carried), (std::vector<Amount>{5, 5, 1, 5}));
    EXPECT_EQ(dense_input(g, plan, s.drops), (std::vector<Amount>{0, 0, 4, 0}));
    // Infinite offer: every source member injected to its residual.
    auto inf = lossy_step(plan, std::nullopt, w);
    EXPECT_EQ(inf.offered, 20);
    EXPECT_EQ(inf.delivered, 11);
    EXPECT_EQ(dense_input(g, plan, inf.drops), (std::vector<Amount>{0, 0, 9, 0}));
  }
  // Highest averages at every node for Lossy Weighted: Capacity weights [1, 3, 3] split 7 as
  // [1, 3, 3] (prototype 1 lossy_weighted_highest_averages).
  auto g2 = make_graph(2, {{0, 1, 1, 1}, {0, 1, 3, 1}, {0, 1, 3, 1}});
  const auto p2 = plan_of(g2, FlowPlacement::Weighted, 0, {1}, all_ids(g2), caps(g2), EdgePriority::Capacity,
                          Admission::Lossy);
  auto s2 = lossy_step(p2, 7, on_members(p2, caps(g2)));
  EXPECT_EQ(dense_input(g2, p2, s2.carried), (std::vector<Amount>{1, 3, 3}));
  EXPECT_EQ(s2.delivered, 7);
  // A source without members injects nothing (the offer is a source drop).
  auto g3 = make_graph(2, {{0, 1, 0, 1}});
  const auto p3 = plan_of(g3, FlowPlacement::EqualBalanced, 0, {1}, all_ids(g3), caps(g3), EdgePriority::Residual,
                          Admission::Lossy);
  auto s3 = lossy_step(p3, 5, {});
  EXPECT_EQ(s3.offered, 5);
  EXPECT_EQ(s3.injected, 0);
  EXPECT_EQ(s3.delivered, 0);
}
