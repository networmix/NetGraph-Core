// Apportionment and the placement kernels in units: the closed-form
// highest averages against the sequential definition, house monotonicity and the Hamilton hole
// control, the round-robin callback form, the lower-quota division on both sides of 2^64; the
// kernel check (the first defect in order, the combination rule alone); plans (membership by
// kernel and priority, split-graph pruning, the walk definition for every kernel), alloc(T) at
// every split, exact T_max against an independent per-node reference, the Lossless
// steps of MaxFlow (replacement on V = W + own, the cancellation DAG, the equal-optimum and
// cancel-upstream controls, the release witness, the mixed-cost bundle; randomized steps
// certified on their member subset against an exhaustive optimum or an independent solver, with
// every cost multiplied, and the step certificate's own controls), Weighted and EqualBalanced,
// and the Lossy steps (drops, +inf injection; on random DAGs against an independent per-node
// reference). Expected values are the closed forms, the hand derivations given in each test, those
// independent optima and the per-node references.
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "netgraph/core/detail/int128.hpp"
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
      if (static_cast<Int128>(w[i]) * (a[best] + 1) > static_cast<Int128>(w[best]) * (a[i] + 1)) best = i;
    }
    ++a[best];
  }
  return a;
}

// Count the exact rank of each quotient w[i] / k in all positive quotients w[j] / h.
// Quotients above it precede it; equal quotients with j <= i also count. The k-th unit belongs
// to i iff that rank is at most total. Binary search each member's last accepted k, using only
// integer comparisons and division, independently of the library's lower-quota allocation.
std::vector<Amount> ha_by_rank(Amount total, const std::vector<Amount>& w) {
  std::vector<Amount> out(w.size(), 0);
  for (std::size_t i = 0; i < w.size(); ++i) {
    Amount lo = 0, hi = total;
    while (lo < hi) {
      const Amount k = lo + (hi - lo + 1) / 2;
      Int128 rank = 0;
      for (std::size_t j = 0; j < w.size(); ++j) {
        const Int128 product = static_cast<Int128>(w[j]) * k;
        rank += (product - 1) / w[i];
        if (product % w[i] == 0 && j <= i) rank += 1;
      }
      if (rank <= total) lo = k;
      else hi = k - 1;
    }
    out[i] = lo;
  }
  return out;
}

// The per-node reference of alloc(T) on a DAG whose node ids are a topological order,
// independent of the library: `total` enters at node `src`; node by node in id order, each node's
// inflow is split over its arcs in the given order by the sequential definition with the arcs'
// weights. Returns each arc's amount.
std::vector<Amount> per_node_reference(int n, const std::vector<std::pair<NodeId, NodeId>>& arcs,
                                       const std::vector<Amount>& weights, Amount total, NodeId src = 0, bool by_rank = false) {
  std::vector<Amount> inflow(static_cast<std::size_t>(n), 0), out(arcs.size(), 0);
  inflow[static_cast<std::size_t>(src)] = total;
  for (NodeId u = 0; u < n; ++u) {
    std::vector<std::size_t> mine;
    std::vector<Amount> w;
    for (std::size_t i = 0; i < arcs.size(); ++i) {
      if (arcs[i].first != u) continue;
      mine.push_back(i);
      w.push_back(weights[i]);
    }
    if (mine.empty()) continue;
    const auto parts = by_rank ? ha_by_rank(inflow[static_cast<std::size_t>(u)], w)
                               : ha_slow(inflow[static_cast<std::size_t>(u)], w);
    for (std::size_t j = 0; j < mine.size(); ++j) {
      out[mine[j]] = parts[j];
      inflow[static_cast<std::size_t>(arcs[mine[j]].second)] += parts[j];
    }
  }
  return out;
}

// The per-node reference of a Lossy step, independent of the library: on a split graph whose node
// ids are a topological order, given as its arcs in ascending EdgeId with each arc's weight and
// residual, the offer enters at `src` (+inf: each arc out of `src` is offered its residual); node
// by node in id order, a node's inflow is split over its arcs one unit at a time (Weighted: each
// unit to the largest weight / (share + 1), ties to the earlier arc; EqualBalanced: the arcs in
// turn), each arc carries the smaller of its share and its residual and drops the rest, and what
// it carries enters its head. A destination splits nothing (no arc of a split graph leaves one).
struct LossyReference {
  std::vector<Amount> carried;
  std::vector<Amount> drops;
  Amount injected {0};
  Amount delivered {0};
};
LossyReference lossy_reference(int n, const std::vector<std::pair<NodeId, NodeId>>& arcs,
                               const std::vector<Amount>& weights, const std::vector<Amount>& residual, bool weighted,
                               std::optional<Amount> offer, NodeId src, const std::vector<NodeId>& dsts) {
  LossyReference out {std::vector<Amount>(arcs.size(), 0), std::vector<Amount>(arcs.size(), 0)};
  std::vector<Amount> inflow(static_cast<std::size_t>(n), 0);
  inflow[static_cast<std::size_t>(src)] = offer.value_or(0);
  for (NodeId u = 0; u < n; ++u) {
    if (std::find(dsts.begin(), dsts.end(), u) != dsts.end()) continue;
    std::vector<std::size_t> mine;
    std::vector<Amount> w;
    for (std::size_t i = 0; i < arcs.size(); ++i) {
      if (arcs[i].first != u) continue;
      mine.push_back(i);
      w.push_back(weights[i]);
    }
    if (mine.empty()) continue;
    std::vector<Amount> share(mine.size(), 0);
    if (u == src && !offer.has_value()) {
      for (std::size_t j = 0; j < mine.size(); ++j) share[j] = residual[mine[j]];
    } else if (weighted) {
      share = ha_slow(inflow[static_cast<std::size_t>(u)], w);
    } else {
      for (Amount unit = 0; unit < inflow[static_cast<std::size_t>(u)]; ++unit) {
        ++share[static_cast<std::size_t>(unit % static_cast<Amount>(mine.size()))];
      }
    }
    for (std::size_t j = 0; j < mine.size(); ++j) {
      const std::size_t i = mine[j];
      if (u == src) out.injected += share[j];
      out.carried[i] = std::min(share[j], residual[i]);
      out.drops[i] = share[j] - out.carried[i];
      inflow[static_cast<std::size_t>(arcs[i].second)] += out.carried[i];
    }
  }
  for (NodeId d : dsts) out.delivered += inflow[static_cast<std::size_t>(d)];
  return out;
}

// Hamilton (largest remainders; not the library's rule): a control only.
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
  return route_cost(g, PredDAGBuilder::of_edges(g, plan.route), plan.roots.front(), plan.dsts);
}

// One Lossy offer of a single-root plan on the member-aligned view `w` (nullopt: the first-hop
// experiment of +inf).
Forwarding lossy_step(const KernelPlan& plan, std::optional<Amount> offer, std::span<const Amount> w) {
  Injection in;
  if (offer.has_value()) {
    in.per_root = std::span<const Amount>(&*offer, 1);
  } else {
    in.first_hop = true;
  }
  return forward(plan, in, w);
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

// `count` random edges of a DAG whose node ids 0 .. n - 1 are a topological order: each u -> v
// with u < v, then its capacity and its cost, drawn in that order.
template <class Rng, class CapacityOf, class CostOf>
std::vector<InEdge> random_dag_edges(Rng& rng, int n, int count, CapacityOf capacity, CostOf cost) {
  std::vector<InEdge> edges;
  for (int i = 0; i < count; ++i) {
    const auto u = static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n - 1));
    const auto v = static_cast<NodeId>(u + 1 + rng() % static_cast<std::uint64_t>(n - 1 - u));
    const Amount c = capacity();
    const Cost w = cost();
    edges.push_back({u, v, c, w});
  }
  return edges;
}

// The walk definition of a split graph, independent of the plan builder: the edges of `edges`
// on some walk from src to the destination set that leaves no destination (destinations are
// terminal), ascending; forward from src and backward from the destinations by fixpoint.
std::vector<EdgeId> walk_members(const Graph& g, const std::vector<EdgeId>& edges, NodeId src,
                                 const std::vector<NodeId>& dsts) {
  const auto n = static_cast<std::size_t>(g.num_nodes());
  std::vector<std::uint8_t> is_dst(n, 0), fwd(n, 0), bwd(n, 0);
  for (NodeId d : dsts) is_dst[static_cast<std::size_t>(d)] = bwd[static_cast<std::size_t>(d)] = 1;
  fwd[static_cast<std::size_t>(src)] = 1;
  std::vector<EdgeId> open;
  for (EdgeId e : edges) {
    if (is_dst[static_cast<std::size_t>(g.edge_src()[static_cast<std::size_t>(e)])] == 0) open.push_back(e);
  }
  for (std::size_t round = 0; round < n; ++round) {
    for (EdgeId e : open) {
      const auto u = static_cast<std::size_t>(g.edge_src()[static_cast<std::size_t>(e)]);
      const auto v = static_cast<std::size_t>(g.edge_dst()[static_cast<std::size_t>(e)]);
      if (fwd[u] != 0) fwd[v] = 1;
      if (bwd[v] != 0) bwd[u] = 1;
    }
  }
  std::vector<EdgeId> out;
  for (EdgeId e : open) {
    if (fwd[static_cast<std::size_t>(g.edge_src()[static_cast<std::size_t>(e)])] != 0 &&
        bwd[static_cast<std::size_t>(g.edge_dst()[static_cast<std::size_t>(e)])] != 0) {
      out.push_back(e);
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

} // namespace

TEST(Apportion, RoundRobinAndHighestAveragesValues) {
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

// Highest averages with zero weights (pooled attribution divides by what each source holds, and a
// source can hold nothing): a zero weight gets nothing and all-zero weights give nothing. The
// contract's examples H-not-one-unit (not "within one unit of proportional"), H-local-bound and
// RR-local-bound, and H-bound-sweep: every total up to 24 over 1 to 4 members with weights 0 to 4,
// not all zero (19,400 cases), equals the sequential definition, and each member's rounding error
// against its proportional share is below k units (k members), below 1 for equal weights.
TEST(Apportion, HighestAveragesExamplesAndZeroWeights) {
  EXPECT_EQ(ha(50, {100, 1, 1, 1}), (std::vector<Amount>{50, 0, 0, 0}));
  std::vector<Amount> local(11, 1);
  local[0] = 1000000;
  std::vector<Amount> want(11, 0);
  want[0] = 1000000;
  EXPECT_EQ(ha(1000000, local), want);
  EXPECT_EQ(rr(10, 3), (std::vector<Amount>{4, 3, 3}));
  EXPECT_EQ(ha(3, {0, 2, 0, 1}), (std::vector<Amount>{0, 2, 0, 1}));
  EXPECT_EQ(ha(5, {0, 0}), (std::vector<Amount>{0, 0}));
  EXPECT_EQ(ha(0, {0, 0, 0}), (std::vector<Amount>{0, 0, 0}));
  int cases = 0;
  for (std::size_t k = 1; k <= 4; ++k) {
    std::vector<Amount> w(k, 0);
    for (;;) {
      // The next weight vector (base 5 counter), all-zero skipped.
      std::size_t i = 0;
      while (i < k && w[i] == 4) w[i++] = 0;
      if (i == k) break;
      ++w[i];
      Amount sum = 0;
      for (Amount x : w) sum += x;
      const bool equal = std::all_of(w.begin(), w.end(), [&](Amount x) { return x == w[0]; });
      for (Amount total = 0; total <= 24; ++total) {
        ++cases;
        const auto got = ha(total, w);
        ASSERT_EQ(got, ha_slow(total, w)) << total;
        for (std::size_t j = 0; j < k; ++j) {
          if (w[j] == 0) EXPECT_EQ(got[j], 0) << total;
          // |got - total * w / sum| < bound, in integers: |got * sum - total * w| < bound * sum.
          const Amount err = got[j] * sum - total * w[j];
          const Amount bound = (equal ? 1 : static_cast<Amount>(k)) * sum;
          EXPECT_LT(err < 0 ? -err : err, bound) << total << " member " << j;
        }
      }
    }
  }
  EXPECT_EQ(cases, 19400);
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
  int weighted = 0;  // allocations that differ from round robin: the weights mattered
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
        weighted += a != rr(total, static_cast<std::size_t>(k)) ? 1 : 0;
      }
    }
  }
  // Equal weights reduce to round robin.
  for (std::size_t k = 1; k <= 9; ++k) {
    for (Amount total = 0; total <= 40; ++total) EXPECT_EQ(ha(total, std::vector<Amount>(k, 7)), rr(total, k));
  }
  EXPECT_GT(weighted, 20000);
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
      const Int128 left = static_cast<Int128>(w[i]) * (a[j] + 1);
      const Int128 right = static_cast<Int128>(w[j]) * a[i];
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
      const Int128 left = static_cast<Int128>(wide[i]) * (b[j] + 1);
      const Int128 right = static_cast<Int128>(wide[j]) * b[i];
      EXPECT_TRUE(left > right || (left == right && i <= j));
    }
  }
}

// The lower-quota division floor(n / d) on both sides of 2^64: dividends with a zero high word
// (the 64-bit path) and with a nonzero one (the 128-bit path), divisors near 1, near the units
// domain bound 2^62 (the call sites divide by a weight sum below it) and near 2^64, against the
// 128-bit quotient of the test integer (test_int128.hpp) and the defining bounds q * d <= n < (q + 1) * d.
TEST(Apportion, QuotaDivisionOnBothSidesOfTheSixtyFourBitBoundary) {
  using u64 = std::uint64_t;
  using u128 = UInt128;
  constexpr u64 kMax = std::numeric_limits<u64>::max();
  const auto domain = static_cast<u64>(kDomainLimit);
  int narrow = 0, wide = 0;
  auto check = [&](u64 hi, u64 lo, u64 d) {
    ASSERT_GT(d, 0u);
    ASSERT_LT(hi, d);  // the quotient fits 64 bits
    const u128 n = (static_cast<u128>(hi) << 64) | lo;
    const u64 q = div_u128(U128 {hi, lo}, d);
    ASSERT_EQ(q, static_cast<u64>(n / d)) << hi << ":" << lo << " / " << d;
    const u128 qd = static_cast<u128>(q) * d;
    ASSERT_TRUE(qd <= n && n - qd < d) << hi << ":" << lo << " / " << d;
    ++(hi == 0 ? narrow : wide);
  };
  const std::vector<u64> divisors = {1,          2,          3,          7,          domain - 2, domain - 1,
                                     domain,     domain + 1, kMax / 2,   kMax - 1,   kMax};
  const std::vector<u64> lows = {0, 1, 2, 3, domain - 1, domain, u64 {1} << 63, kMax - 1, kMax};
  for (u64 d : divisors) {
    for (u64 lo : lows) {
      check(0, lo, d);
      for (u64 hi : {u64 {1}, u64 {2}, d / 2, d - 2, d - 1}) {
        if (hi > 0 && hi < d) check(hi, lo, d);
      }
    }
    // The dividends next to a multiple of d on both sides of 2^64.
    for (u64 q : {u64 {1}, u64 {2}, kMax / d, kMax - 1, kMax}) {
      const u128 at = static_cast<u128>(q) * d;
      for (u128 n : {at - 1, at, at + 1}) {
        if ((n >> 64) < d) check(static_cast<u64>(n >> 64), static_cast<u64>(n), d);
      }
    }
  }
  // The call-site shape: total * w / W with total and the weights below 2^62, W the weight sum.
  for (u64 total : {u64 {1}, u64 {50}, domain - 2, domain - 1}) {
    for (u64 w : {u64 {1}, u64 {3}, domain / 3, domain - 1}) {
      for (u64 sum : {w, w + 1, domain - 1}) {  // every sum at least w
        const U128 n = mul_u64(total, w);
        check(n.hi, n.lo, sum);
      }
    }
  }
  std::mt19937_64 rng(20261006);
  for (int i = 0; i < 20000; ++i) {
    const u64 d = i % 3 == 0 ? 1 + rng() % 4096 : (i % 3 == 1 ? 1 + rng() % (domain - 1) : 1 + rng() % kMax);
    const u64 hi = i % 2 == 0 ? 0 : rng() % d;
    check(hi, rng(), d);
  }
  EXPECT_GT(narrow, 10000);
  EXPECT_GT(wide, 10000);
}

// The exact 64 x 64 -> 128-bit product: the built-in path where the compiler has a 128-bit
// integer, the 32-bit-limb path elsewhere (MSVC), both against the test integer (test_int128.hpp,
// independent of the library) at the word boundaries and on random operands.
TEST(Apportion, MulU64IsTheExactProduct) {
  using u64 = std::uint64_t;
  constexpr u64 kMax = std::numeric_limits<u64>::max();
  const std::vector<u64> edges = {0, 1, 2, 3, 0xffffffffULL, 0x100000000ULL, 0x100000001ULL, (u64 {1} << 62) - 1,
                                  u64 {1} << 62, u64 {1} << 63, kMax - 1, kMax};
  auto check = [](u64 a, u64 b) {
    const U128 p = mul_u64(a, b);
    const UInt128 want = static_cast<UInt128>(a) * static_cast<UInt128>(b);
    ASSERT_EQ(p.hi, static_cast<u64>(want >> 64)) << a << " * " << b;
    ASSERT_EQ(p.lo, static_cast<u64>(want)) << a << " * " << b;
  };
  for (u64 a : edges) {
    for (u64 b : edges) check(a, b);
  }
  std::mt19937_64 rng(1729);
  for (int i = 0; i < 20000; ++i) {
    const int shift_a = static_cast<int>(rng() % 64), shift_b = static_cast<int>(rng() % 64);
    check(rng() >> shift_a, rng() >> shift_b);
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

// The kernel entries keep their complete check: with several defects at once, one
// error, the first in the order placement, admission, priority (named "priority": a kernel entry
// takes it as a field of its own), then the combinations Lossy with MaxFlow and Weighted with
// Headroom; at check_kernel and at the three entries that call it.
TEST(Plans, KernelEntriesReportTheFirstDefectInOrder) {
  const auto bad_placement = static_cast<FlowPlacement>(99);
  const auto bad_admission = static_cast<Admission>(98);
  const auto bad_priority = static_cast<EdgePriority>(97);
  struct Case {
    FlowPlacement k;
    Admission a;
    EdgePriority p;
    std::string error;
  };
  const std::vector<Case> cases = {
      {bad_placement, bad_admission, bad_priority, "placement: 99 is not a value of FlowPlacement"},
      {bad_placement, Admission::Lossy, EdgePriority::Headroom, "placement: 99 is not a value of FlowPlacement"},
      {FlowPlacement::MaxFlow, bad_admission, bad_priority, "admission: 98 is not a value of Admission"},
      {FlowPlacement::MaxFlow, Admission::Lossy, bad_priority, "priority: 97 is not a value of EdgePriority"},
      {FlowPlacement::Weighted, Admission::Lossy, bad_priority, "priority: 97 is not a value of EdgePriority"},
      {FlowPlacement::MaxFlow, Admission::Lossy, EdgePriority::Headroom, kLossyKernelError},
      {FlowPlacement::Weighted, Admission::Lossy, EdgePriority::Headroom, kHeadroomKernelError},
      {FlowPlacement::Weighted, Admission::Lossless, EdgePriority::Headroom, kHeadroomKernelError},
  };
  const auto g = make_graph(3, {{0, 1, 2, 1}, {1, 2, 2, 1}});
  const auto support = all_edges(g);
  const std::vector<NodeId> dsts = {2};
  const PriorityValues values {&g, EdgePriority::Capacity, {}};
  for (const auto& c : cases) {
    EXPECT_EQ(error_of([&] { check_kernel(c.k, c.a, c.p); }), c.error);
    EXPECT_EQ(error_of([&] { (void)make_plan(g, c.k, c.a, c.p, 0, dsts, *support, {}); }), c.error);
    EXPECT_EQ(error_of([&] { (void)make_route_plan(g, c.k, c.a, c.p, 0, dsts, *support, nullptr); }), c.error);
    PathSelection selection(EdgeTieBreak::Deterministic, EdgePriority::Capacity);
    selection.priority = c.p;
    EXPECT_EQ(error_of([&] { (void)cost_only_creation(g, support, 0, dsts, selection, PriorityFn::of(values), c.k, c.a); }),
              c.error);
  }
}

// validate_config's part of the kernel check: the two combinations only,
// never a domain, which the configuration has checked under its own field names just before.
TEST(Plans, TheCombinationRuleAloneChecksNoDomain) {
  const auto bad_placement = static_cast<FlowPlacement>(99);
  const auto bad_admission = static_cast<Admission>(98);
  const auto bad_priority = static_cast<EdgePriority>(97);
  EXPECT_EQ(error_of([&] { check_kernel_combination(FlowPlacement::MaxFlow, Admission::Lossy, bad_priority); }),
            kLossyKernelError);
  EXPECT_EQ(error_of([&] { check_kernel_combination(FlowPlacement::Weighted, bad_admission, EdgePriority::Headroom); }),
            kHeadroomKernelError);
  EXPECT_EQ(error_of([&] { check_kernel_combination(bad_placement, bad_admission, bad_priority); }), "");
  for (auto k : {FlowPlacement::MaxFlow, FlowPlacement::Weighted, FlowPlacement::EqualBalanced}) {
    for (auto a : {Admission::Lossless, Admission::Lossy}) {
      for (auto pr : {EdgePriority::Residual, EdgePriority::Capacity, EdgePriority::Headroom}) {
        const bool rejected = (a == Admission::Lossy && k == FlowPlacement::MaxFlow) ||
                              (k == FlowPlacement::Weighted && pr == EdgePriority::Headroom);
        EXPECT_EQ(error_of([&] { check_kernel_combination(k, a, pr); }).empty(), !rejected);
        EXPECT_EQ(error_of([&] { check_kernel_combination(k, a, pr); }), error_of([&] { check_kernel(k, a, pr); }));
      }
    }
  }
}

TEST(Plans, MembershipByKernelAndPriority) {
  // Parallel links of capacities [4, 4, 0] with residual [0, 4, 0] on the creation view.
  auto g = make_graph(2, {{0, 1, 4, 1}, {0, 1, 4, 1}, {0, 1, 0, 1}});
  const auto view = to_ids(g, std::vector<Amount> {0, 4, 0});
  const auto inv = id_of_input(g);
  auto mf = plan_of(g, FlowPlacement::MaxFlow, 0, {1}, all_ids(g), view);
  EXPECT_EQ(mf.members, (std::vector<EdgeId>{0, 1, 2}));  // MaxFlow: the view decides
  // EqualBalanced: every route edge, the zero-capacity one included (it takes its round-robin
  // share and carries nothing).
  auto eb = plan_of(g, FlowPlacement::EqualBalanced, 0, {1}, all_ids(g), view);
  EXPECT_EQ(eb.members, (std::vector<EdgeId>{0, 1, 2}));
  // Weighted: the edges of positive weight (residual weights on the creation view, capacity
  // weights on the capacities).
  auto wr = plan_of(g, FlowPlacement::Weighted, 0, {1}, all_ids(g), view, EdgePriority::Residual);
  EXPECT_EQ(wr.members, (std::vector<EdgeId>{inv[1]}));
  EXPECT_EQ(wr.weights, (std::vector<Amount>{4}));
  auto wc = plan_of(g, FlowPlacement::Weighted, 0, {1}, all_ids(g), view, EdgePriority::Capacity);
  EXPECT_EQ(wc.members, (std::vector<EdgeId>{inv[0], inv[1]}));
  EXPECT_EQ(wc.weights, (std::vector<Amount>{4, 4}));
  // Lossy keeps its kernel's membership and weights.
  for (const auto& [kernel, prio, lossless] : {std::tuple {FlowPlacement::EqualBalanced, EdgePriority::Residual, &eb},
                                               std::tuple {FlowPlacement::Weighted, EdgePriority::Residual, &wr},
                                               std::tuple {FlowPlacement::Weighted, EdgePriority::Capacity, &wc}}) {
    const auto lossy = plan_of(g, kernel, 0, {1}, all_ids(g), view, prio, Admission::Lossy);
    EXPECT_EQ(lossy.members, lossless->members) << static_cast<int>(kernel) << " " << static_cast<int>(prio);
    EXPECT_EQ(lossy.weights, lossless->weights) << static_cast<int>(kernel) << " " << static_cast<int>(prio);
  }
}

TEST(Plans, SplitGraphPruning) {
  // A zero-capacity tail. 0->1, 0->2, 1->3, 2->4 (cap 0), 4->3. EqualBalanced keeps the branch
  // through 2: its share is stranded on 2->4, so T = 2 puts a unit there and T_max is 1.
  // Weighted gives 2->4 weight 0: node 2 is left without a member, so the branch is pruned and
  // the whole total takes the other branch.
  auto g = make_graph(5, {{0, 1, 2, 1}, {0, 2, 2, 1}, {1, 3, 2, 1}, {2, 4, 0, 1}, {4, 3, 2, 1}});
  const auto inv = id_of_input(g);
  auto eb = plan_of(g, FlowPlacement::EqualBalanced, 0, {3}, all_ids(g), caps(g));
  EXPECT_EQ(eb.members, all_ids(g));
  std::vector<Amount> a(eb.members.size());
  ASSERT_TRUE(allocate(eb, 2, a));
  EXPECT_EQ(dense_input(g, eb, a), (std::vector<Amount>{1, 1, 1, 1, 1}));
  EXPECT_EQ(t_max(eb, caps(g)).value, 1);  // every edge is a member: the capacities are member-aligned
  auto w = plan_of(g, FlowPlacement::Weighted, 0, {3}, all_ids(g), caps(g), EdgePriority::Capacity);
  std::vector<EdgeId> want = {inv[0], inv[2]};
  std::sort(want.begin(), want.end());
  EXPECT_EQ(w.members, want);
  std::vector<Amount> b(w.members.size());
  ASSERT_TRUE(allocate(w, 2, b));
  EXPECT_EQ(dense_input(g, w, b), (std::vector<Amount>{2, 0, 2, 0, 0}));
  // Destinations are terminal; a destination set keeps both destinations.
  auto g2 = make_graph(4, {{0, 1, 3, 1}, {1, 2, 3, 1}, {0, 2, 3, 2}, {2, 3, 3, 1}});
  auto any = plan_of(g2, FlowPlacement::EqualBalanced, 0, {1, 2}, all_ids(g2), caps(g2));
  const auto inv2 = id_of_input(g2);
  std::vector<EdgeId> want2 = {inv2[0], inv2[2]};  // 1->2 leaves destination 1; 2->3 leaves 2
  std::sort(want2.begin(), want2.end());
  EXPECT_EQ(any.members, want2);
  EXPECT_EQ(cost_of(g2, any), 1);
  // An empty split graph (a route that reaches no destination): route cost unreachable, T_max 0.
  auto none = plan_of(g, FlowPlacement::EqualBalanced, 0, {3}, std::vector<EdgeId> {inv[1], inv[3]}, caps(g));
  EXPECT_TRUE(none.members.empty());
  EXPECT_EQ(cost_of(g, none), kUnreachable);
  EXPECT_EQ(t_max(none, {}).value, 0);
}

TEST(Plans, SplitGraphEqualsTheWalkDefinition) {
  // On random DAGs (parallel edges, zero capacities, statically excluded edges and nodes, edges
  // into the source and out of a destination, one or two destinations) the route of a bundle is,
  // by definition, its allowed edges, whatever their capacity, on some walk from the source to a
  // destination that enters no destination before its end, and the split graph is the route for
  // MaxFlow and EqualBalanced, and for Weighted (capacity weights) the same definition over the
  // allowed edges of at least 1 unit: the route and the members equal an independent
  // reachability computation, the route cost equals an independent shortest walk over the
  // allowed edges, the topological structure orders every member's tail before its head, and an
  // allocation conserves its total and equals the per-node reference (each node's inflow split
  // over its members by the kernel's rule). make_plan on the DAG of the allowed edges gives the
  // same plan as make_route_plan filtering the bundle (the two entry points of the one builder),
  // for every kernel, but for the route field. On the bundle itself (no Eligibility) make_plan's
  // MaxFlow and EqualBalanced members are the walk definition over all the bundle's edges.
  std::mt19937 rng(2024);
  int compared = 0;
  int masked_differs = 0;    // trials where the Eligibility changes the route
  int weighted_differs = 0;  // trials where Weighted's members omit a zero-capacity route edge
  int splitting = 0;        // such plans with a split, and two destinations or a source other than 0
  for (int trial = 0; trial < 2000; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 8);
    auto g = make_graph(n, random_dag_edges(
                               rng, n, 2 * n + static_cast<int>(rng() % 10),
                               [&] { return (rng() % 6) == 0 ? Amount {0} : 1 + static_cast<Amount>(rng() % 9); },
                               [&] { return 1 + static_cast<Cost>(rng() % 4); }));
    std::vector<EdgeId> in_bundle;
    for (EdgeId e : all_ids(g)) {
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
      for (EdgeId e : all_ids(g)) {
        if ((rng() % 7) == 0) xe.push_back(e);
      }
      elig = g.eligibility(xn, xe);
    }
    const auto capacity = g.capacity_units();
    const auto edge_src = g.edge_src();
    const auto edge_dst = g.edge_dst();
    const auto cost = g.cost();
    std::vector<EdgeId> usable;    // allowed
    std::vector<EdgeId> positive;  // allowed, at least 1 unit
    for (EdgeId e : PredDAGBuilder::ascending(bundle)) {
      if (elig && !elig->edge_ok(e)) continue;
      usable.push_back(e);
      if (capacity[static_cast<std::size_t>(e)] >= 1) positive.push_back(e);
    }
    // The definition over the usable edges; the shortest walk by Bellman-Ford over those that
    // leave no destination.
    std::vector<std::uint8_t> is_dst(static_cast<std::size_t>(n), 0);
    for (NodeId d : dsts) is_dst[static_cast<std::size_t>(d)] = 1;
    const auto want = walk_members(g, usable, src, dsts);
    const auto weighted_want = walk_members(g, positive, src, dsts);
    weighted_differs += weighted_want != want ? 1 : 0;
    std::vector<Cost> dist(static_cast<std::size_t>(n), kUnreachable);
    dist[static_cast<std::size_t>(src)] = 0;
    for (int round = 0; round < n; ++round) {
      for (EdgeId e : usable) {
        const auto u = static_cast<std::size_t>(edge_src[static_cast<std::size_t>(e)]);
        const auto v = static_cast<std::size_t>(edge_dst[static_cast<std::size_t>(e)]);
        if (is_dst[u] == 0 && dist[u] != kUnreachable) dist[v] = std::min(dist[v], dist[u] + cost[static_cast<std::size_t>(e)]);
      }
    }
    Cost want_cost = kUnreachable;
    if (!want.empty()) {
      for (NodeId d : dsts) want_cost = std::min(want_cost, dist[static_cast<std::size_t>(d)]);
    }
    const auto usable_dag = PredDAGBuilder::of_edges(g, usable);
    const auto bundle_edges = PredDAGBuilder::ascending(bundle);
    const auto unmasked_want = walk_members(g, std::vector<EdgeId>(bundle_edges.begin(), bundle_edges.end()), src, dsts);
    for (const auto kernel : {FlowPlacement::MaxFlow, FlowPlacement::EqualBalanced}) {
      EXPECT_EQ(make_plan(g, kernel, Admission::Lossless, EdgePriority::Capacity, src, dsts, bundle, {}).members,
                unmasked_want)
          << "trial " << trial;
    }
    masked_differs += unmasked_want != want ? 1 : 0;
    for (const auto kernel : {FlowPlacement::EqualBalanced, FlowPlacement::Weighted, FlowPlacement::MaxFlow}) {
      Cost got_cost = -1;
      const auto on_bundle = make_route_plan(g, kernel, Admission::Lossless, EdgePriority::Capacity, src, dsts, bundle,
                                             elig ? &*elig : nullptr, &got_cost);
      ASSERT_EQ(on_bundle.members, kernel == FlowPlacement::Weighted ? weighted_want : want) << "trial " << trial;
      EXPECT_EQ(on_bundle.route, want);
      EXPECT_EQ(got_cost, want_cost) << "trial " << trial;
      EXPECT_EQ(prune(g, usable_dag, src, dsts), want);
      Cost plan_cost = -1;
      const auto on_usable =
          make_plan(g, kernel, Admission::Lossless, EdgePriority::Capacity, src, dsts, usable_dag, {}, &plan_cost);
      // The same plan but for its route, which is the whole allowed DAG (plan_differences,
      // test_utils.hpp).
      EXPECT_EQ(plan_differences(on_usable, on_bundle, false), std::vector<std::string> {}) << "trial " << trial;
      // The plan built on the usable edges reports the same shortest-walk cost, and so does the
      // cost pass alone.
      EXPECT_EQ(plan_cost, want_cost) << "trial " << trial;
      EXPECT_EQ(route_cost(g, usable_dag, src, dsts), want_cost) << "trial " << trial;
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
      // Conservation: what enters the destinations is the total; and the per-node reference (node
      // ids are a topological order; the members in ascending EdgeId, Weighted's weights the
      // capacities). MaxFlow places by the flow core, not by alloc(T).
      if (kernel == FlowPlacement::MaxFlow) continue;
      std::vector<std::pair<NodeId, NodeId>> arcs;
      std::vector<Amount> weights;
      for (EdgeId e : on_bundle.members) {
        arcs.emplace_back(edge_src[static_cast<std::size_t>(e)], edge_dst[static_cast<std::size_t>(e)]);
        weights.push_back(kernel == FlowPlacement::Weighted ? capacity[static_cast<std::size_t>(e)] : 1);
      }
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
        EXPECT_EQ(a, per_node_reference(n, arcs, weights, t, src)) << "trial " << trial << " T " << t;
      }
      // What the single-destination tests below do not reach: a split at some node, with two
      // destinations or a source other than node 0.
      bool splits = false;
      for (std::size_t p = 0; p < on_bundle.topo.size(); ++p) {
        splits = splits || (on_bundle.terminal[p] == 0 && on_bundle.out_end[p] - on_bundle.out_begin[p] >= 2);
      }
      splitting += splits && (dsts.size() > 1 || src != 0) ? 1 : 0;
      ++compared;
    }
  }
  EXPECT_GT(compared, 1000);
  EXPECT_GT(masked_differs, 100);
  EXPECT_GT(weighted_differs, 100);
  EXPECT_GT(splitting, 300);
}

TEST(Allocation, SplitAndReconvergeMatchesThePerNodeReferenceAtEveryValue) {
  // House monotonicity of the apportionment on a DAG: capacities 20 x weights, baseline view.
  const std::vector<std::pair<NodeId, NodeId>> dag = {{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 5}, {4, 5}};
  const std::vector<Amount> w = {2, 7, 1, 4, 3, 2, 1, 1};
  const std::vector<Amount> view_in = {8, 20, 4, 9, 10, 10, 15, 15};
  std::vector<InEdge> edges;
  for (std::size_t i = 0; i < dag.size(); ++i) edges.push_back({dag[i].first, dag[i].second, 20 * w[i], 1});
  auto g = make_graph(6, edges);
  const auto view = to_ids(g, view_in);
  for (auto kernel : {FlowPlacement::Weighted, FlowPlacement::EqualBalanced}) {
    // The per-node reference on input positions (node ids are topological); EqualBalanced is the
    // sequential definition with equal weights (ties to the lower index).
    const auto weights = kernel == FlowPlacement::EqualBalanced ? std::vector<Amount>(w.size(), 1) : w;
    auto ref = [&](Amount total) { return per_node_reference(6, dag, weights, total); };
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
  // EqualBalanced on a DAG that splits three ways at the source and reconverges at node 4.
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
  // Random DAGs from node 0 to node n - 1 (node ids topological). The reference is independent of
  // the library: the members by the walk definition, alloc(T) node by node in id order, each
  // node's inflow split over its members in ascending EdgeId by the sequential definition
  // (capacities as Weighted's weights, equal weights for EqualBalanced), and T_max the last T
  // whose allocation fits the view (the allocation is member-monotone, so the fitting T are a
  // prefix). allocate() equals the reference at every T up to one past T_max, its limit answer
  // agrees, and t_max returns that T.
  std::mt19937 rng(99);
  int several = 0;  // plans whose T_max is at least 2 units
  for (int trial = 0; trial < 400; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 5);
    auto g = make_graph(n, random_dag_edges(
                               rng, n, 2 * n + static_cast<int>(rng() % 6),
                               [&] { return 1 + static_cast<Amount>(rng() % 9); }, [] { return Cost {1}; }));
    std::vector<Amount> view(static_cast<std::size_t>(g.num_edges()));
    for (std::size_t e = 0; e < view.size(); ++e) {
      view[e] = static_cast<Amount>(rng() % static_cast<std::uint32_t>(g.capacity_units()[e] + 1));
    }
    const auto kernel = (rng() % 2) == 0 ? FlowPlacement::Weighted : FlowPlacement::EqualBalanced;
    const auto plan = plan_of(g, kernel, 0, {n - 1}, all_ids(g), caps(g), EdgePriority::Capacity);
    ASSERT_EQ(plan.members, walk_members(g, all_ids(g), 0, {n - 1})) << "trial " << trial;
    if (plan.members.empty()) continue;
    // The members as arcs in ascending EdgeId, with their weights.
    std::vector<std::pair<NodeId, NodeId>> arcs;
    std::vector<Amount> weights;
    for (EdgeId m : plan.members) {
      const auto e = static_cast<std::size_t>(m);
      arcs.emplace_back(g.edge_src()[e], g.edge_dst()[e]);
      weights.push_back(kernel == FlowPlacement::Weighted ? g.capacity_units()[e] : 1);
    }
    auto reference = [&](Amount total) { return per_node_reference(n, arcs, weights, total); };
    const auto vm = on_members(plan, view);
    const auto tm = t_max(plan, vm);
    EXPECT_LE(tm.iterations, 62);
    Amount brute = -1;
    std::vector<Amount> a(plan.members.size()), limited(plan.members.size());
    for (Amount t = 0; t <= 200; ++t) {
      const auto want = reference(t);
      bool fits = true;
      for (std::size_t k = 0; k < want.size(); ++k) fits = fits && want[k] <= vm[k];
      ASSERT_TRUE(allocate(plan, t, a)) << "trial " << trial << " T " << t;
      ASSERT_EQ(a, want) << "trial " << trial << " T " << t;
      ASSERT_EQ(allocate(plan, t, limited, vm), fits) << "trial " << trial << " T " << t;
      if (!fits) break;
      brute = t;
    }
    ASSERT_GE(brute, 0);
    ASSERT_EQ(tm.value, brute) << "trial " << trial;
    several += tm.value >= 2 ? 1 : 0;
  }
  EXPECT_GT(several, 60);
}

TEST(Allocation, TMaxSearchIsExactAndWithinTheProbeBudgetOnLargeRanges) {
  // Capacities and views up to 2^58 units (search ranges near 2^62 at the source): the seeded
  // search must return the same T as a plain bisection from 0 and probe at most 62 times, also
  // when its floating-point estimate is far off (a member capped near zero deep in the DAG).
  // First cross-check the quotient-rank reference against the unit-by-unit definition,
  // including ties. The large trials then use that reference at each node.
  for (Amount total = 0; total < 50; ++total) {
    for (const std::vector<Amount>& weights : {std::vector<Amount> {1, 1, 1}, {1, 3, 7}, {9, 2, 9, 1}}) {
      ASSERT_EQ(ha_by_rank(total, weights), ha_slow(total, weights));
    }
  }
  std::mt19937_64 rng(1234);
  int large = 0;  // T at least 2^40
  int far = 0;    // T at most 2^-20 of the source side's capacity (the estimate far off)
  for (int trial = 0; trial < 300; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 6);
    auto g = make_graph(n, random_dag_edges(
                               rng, n, 2 * n + static_cast<int>(rng() % 8),
                               [&] { return 1 + static_cast<Amount>(rng() % (std::uint64_t {1} << 58)); },
                               [] { return Cost {1}; }));
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
    // Reference: each node splits by exact quotient rank, then plain bisection over feasible
    // totals. No production allocation or search is used to compute an expected amount.
    Amount hi = 0;
    const auto root = static_cast<std::size_t>(plan.root_pos.front());
    for (auto k = plan.out_begin[root]; k < plan.out_end[root]; ++k) hi += vm[static_cast<std::size_t>(k)];
    const Amount source_side = hi;
    ASSERT_EQ(plan.members, walk_members(g, all_ids(g), 0, {n - 1}));
    std::vector<std::pair<NodeId, NodeId>> arcs;
    std::vector<Amount> weights;
    for (EdgeId edge : plan.members) {
      const auto e = static_cast<std::size_t>(edge);
      arcs.emplace_back(g.edge_src()[e], g.edge_dst()[e]);
      weights.push_back(kernel == FlowPlacement::Weighted ? g.capacity_units()[e] : 1);
    }
    Amount lo = 0;
    while (lo < hi) {
      const Amount mid = lo + (hi - lo + 1) / 2;
      const auto a = per_node_reference(n, arcs, weights, mid, 0, true);
      bool fits = true;
      for (std::size_t j = 0; j < a.size(); ++j) fits = fits && a[j] <= vm[j];
      if (fits) lo = mid; else hi = mid - 1;
    }
    ASSERT_EQ(tm.value, lo) << "trial " << trial;
    // The Lossless step from T = 0 with a request beyond the maximum lands on the same T.
    const std::vector<Amount> zero(plan.members.size(), 0);
    const auto step = lossless_step(g, plan, 0, vm, zero, kDomainLimit - 1);
    EXPECT_EQ(step.total, lo);
    EXPECT_LE(step.iterations, 62);
    large += lo >= (Amount {1} << 40) ? 1 : 0;
    far += source_side >= (Amount {1} << 40) && lo <= (source_side >> 20) ? 1 : 0;
  }
  EXPECT_GT(large, 10);
  EXPECT_GT(far, 100);
}

TEST(MaxFlowKernel, CancellationDagGivesTwo) {
  // A DAG where the second unit needs the first one's edge 1->3 cancelled: place(1) takes
  // 0-1-3-5, the next place(1) re-solves on V = W + own and returns the unique value-2 flow (1->3
  // goes to 0).
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

// A graph whose minimum-cost flow cancels earlier flow: s=0; A=1, B=2, C=3; X=4, Y=5; t=6.
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

TEST(MaxFlowKernel, EqualOptimaFollowTheAscendingEdgeIdOrder) {
  // Ties between equal optima follow the arc order (ascending EdgeId):
  // with cost(B->Y) = 8 two value-3 optima cost 19; ascending arcs pick the
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
  // On a mixed-cost bundle the 1 + 1 + 1 chain fills before the cost-4 edge.
  auto g2 = make_graph(4, {{0, 1, 1, 1}, {1, 2, 1, 1}, {2, 3, 1, 1}, {0, 3, 1, 4}});
  const auto p2 = plan_of(g2, FlowPlacement::MaxFlow, 0, {3}, all_ids(g2), caps(g2));
  EXPECT_EQ(dense_input(g2, p2, first_step(g2, p2, caps(g2), 1).alloc), (std::vector<Amount>{1, 1, 1, 0}));
  EXPECT_EQ(cost_of(g2, p2), 3);
}

namespace {

// The certificate of one MaxFlow-kernel step, independent of the solver: the step's member
// allocation, as a flow on the graph with the members' view and every other edge masked, is
// feasible and conserved from the plan's source to its destinations, has the value
// min(total + request, F_max) (the maximum for a request of +inf) and the minimum cost at that
// value, and its curve has the flow's total and moment and the optimum's cost at every value up
// to it (the curve clauses shared with certify, test_utils.hpp). The optimum is the exhaustive
// one (`tier_a`, small graphs) or the independent Bellman-Ford curve (test_utils.hpp). The
// checker relies on the members lying on walks from the source to its destinations, which end
// there, so that no member enters the source or leaves a destination: the random test below
// checks each plan's members against that walk definition first. Returns the failed clauses.
std::vector<std::string> certify_step(const Graph& g, const KernelPlan& plan, const std::vector<Amount>& view_m,
                                      Amount total, std::optional<Amount> request, const LosslessStep& step,
                                      bool tier_a) {
  std::vector<std::string> fails;
  if (step.alloc.size() != plan.members.size()) {
    fails.push_back("shape");
    return fails;
  }
  const auto m = static_cast<std::size_t>(g.num_edges());
  const auto n = static_cast<std::size_t>(g.num_nodes());
  std::vector<Amount> view(m, 0);
  std::vector<Amount> flow(m, 0);
  std::vector<char> member(m, 0);
  for (std::size_t k = 0; k < plan.members.size(); ++k) {
    const auto e = static_cast<std::size_t>(plan.members[k]);
    member[e] = 1;
    view[e] = view_m[k];
    flow[e] = step.alloc[k];
    if (step.alloc[k] < 0 || step.alloc[k] > view_m[k]) fails.push_back("feasibility at member " + std::to_string(k));
  }
  std::vector<EdgeId> outside;
  for (std::size_t e = 0; e < m; ++e) {
    if (member[e] == 0) outside.push_back(static_cast<EdgeId>(e));
  }
  const auto mask = g.eligibility({}, outside);
  std::vector<Int128> bal(n, 0);
  for (std::size_t e = 0; e < m; ++e) {
    bal[static_cast<std::size_t>(g.edge_src()[e])] -= flow[e];
    bal[static_cast<std::size_t>(g.edge_dst()[e])] += flow[e];
  }
  const Int128 moment = flow_cost(g, flow);
  const std::vector<NodeId> src(plan.roots.begin(), plan.roots.end());
  const std::vector<NodeId> dsts(plan.dsts.begin(), plan.dsts.end());
  Int128 inflow = 0;
  for (std::size_t v = 0; v < n; ++v) {
    if (static_cast<NodeId>(v) == plan.roots.front()) continue;
    if (std::find(dsts.begin(), dsts.end(), static_cast<NodeId>(v)) != dsts.end()) {
      inflow += bal[v];
    } else if (bal[v] != 0) {
      fails.push_back("conservation at node " + std::to_string(v));
    }
  }
  if (-bal[static_cast<std::size_t>(plan.roots.front())] != step.total || inflow != step.total) {
    fails.push_back("value differs from supply or inflow");
  }
  // The optimum's cost at every value up to F_max, as breakpoints.
  std::vector<std::pair<Amount, Int128>> optimum;
  if (tier_a) {
    for (const auto& [x, c] : exhaustive_costs(g, view, src, dsts, &mask)) optimum.emplace_back(x, c);
  } else {
    optimum = independent_curve(g, view, src, dsts, &mask);
  }
  const Amount fmax = optimum.back().first;
  const Amount top = request.has_value() ? std::min(total + *request, fmax) : fmax;
  if (step.total != top) fails.push_back("placed differs from min(cap, F_max)");
  if (step.total >= 0 && step.total <= fmax && moment != curve_at(optimum, step.total)) {
    fails.push_back("cost differs from the optimum");
  }
  const auto pts = curve_shape_clauses(step.cost_distribution, step.total, moment, fails);
  curve_optimum_clause(pts, optimum, std::min(step.total, fmax), "the optimum", fails);
  return fails;
}

} // namespace

TEST(MaxFlowKernel, StepCertificateRejectsMutatedSteps) {
  // The step certificate is not vacuous: each clean step certifies, then each mutation is rejected
  // by the clause it breaks.
  using Fails = std::vector<std::string>;
  auto mutated = [](const Graph& g, const KernelPlan& plan, std::optional<Amount> request,
                    const std::function<void(LosslessStep&)>& mutate) {
    const auto vm = on_members(plan, caps(g));
    auto step = first_step(g, plan, caps(g), request);
    EXPECT_TRUE(certify_step(g, plan, vm, 0, request, step, true).empty());
    mutate(step);
    return certify_step(g, plan, vm, 0, request, step, true);
  };
  // The chain 0 -> 1 -> 2 of capacity 2: two units, members in chain order.
  const auto chain = make_graph(3, {{0, 1, 2, 1}, {1, 2, 2, 1}});
  const auto cp = plan_of(chain, FlowPlacement::MaxFlow, 0, {2}, all_ids(chain), caps(chain));
  ASSERT_EQ(cp.members.size(), 2U);
  EXPECT_EQ(mutated(chain, cp, std::nullopt, [](auto& s) { s.alloc.push_back(0); }), Fails {"shape"});
  EXPECT_TRUE(has_failure(mutated(chain, cp, std::nullopt, [](auto& s) { s.alloc = {3, 3}; }), "feasibility at member 0"));
  EXPECT_TRUE(has_failure(mutated(chain, cp, std::nullopt, [](auto& s) { s.alloc = {-1, -1}; }), "feasibility at member 0"));
  // One unit less into node 1 lowers the supply only, one unit less out of it the inflow only.
  for (const std::size_t k : {std::size_t {0}, std::size_t {1}}) {
    const auto f = mutated(chain, cp, std::nullopt, [&](auto& s) { s.alloc[k] = 1; });
    EXPECT_TRUE(has_failure(f, "conservation at node 1")) << k;
    EXPECT_TRUE(has_failure(f, "value differs from supply or inflow")) << k;
  }
  // A capped step (request 1) certified as if the request were 2: the value alone differs.
  {
    const auto vm = on_members(cp, caps(chain));
    const auto step = first_step(chain, cp, caps(chain), 1);
    EXPECT_EQ(certify_step(chain, cp, vm, 0, 2, step, true), Fails {"placed differs from min(cap, F_max)"});
  }
  // Parallel members of cost 1 and 3 and capacity 1: the first unit on the cheap one, the curve
  // [(1, 1), (3, 1)] for both. The unit moved to the costly member keeps the value; the curve's
  // last volume raised breaks its total; its first cost raised breaks its moment; the one tier
  // [(2, 2)] keeps the total and the moment but not the optimum at 1.
  const auto par = make_graph(2, {{0, 1, 1, 1}, {0, 1, 1, 3}});
  const auto pp = plan_of(par, FlowPlacement::MaxFlow, 0, {1}, all_ids(par), caps(par));
  ASSERT_EQ(pp.members.size(), 2U);
  EXPECT_TRUE(has_failure(mutated(par, pp, 1, [](auto& s) { s.alloc = {0, 1}; }), "cost differs from the optimum"));
  EXPECT_TRUE(has_failure(mutated(par, pp, std::nullopt, [](auto& s) { s.cost_distribution = {{1, 1}, {3, 2}}; }),
                          "curve total differs from placed"));
  EXPECT_TRUE(has_failure(mutated(par, pp, std::nullopt, [](auto& s) { s.cost_distribution = {{2, 1}, {3, 1}}; }),
                          "curve moment differs from the flow cost"));
  EXPECT_EQ(mutated(par, pp, std::nullopt, [](auto& s) { s.cost_distribution = {{2, 2}}; }),
            Fails {"curve differs from the optimum at 1"});
  // The curve's shape, each the only failure: an entry of volume 0 between the two (total, moment
  // and the curve's values unchanged), and the chain's one tier [(2, 2)] split in two at its cost.
  EXPECT_EQ(mutated(par, pp, std::nullopt, [](auto& s) { s.cost_distribution = {{1, 1}, {2, 0}, {3, 1}}; }),
            Fails {"curve volume not positive"});
  EXPECT_EQ(mutated(chain, cp, std::nullopt, [](auto& s) { s.cost_distribution = {{2, 1}, {2, 1}}; }),
            Fails {"curve costs not strictly increasing"});
  // The guard of the replacement solve: a committed total beyond what the view can carry.
  const auto vm = on_members(pp, caps(par));
  const std::vector<Amount> own(pp.members.size(), 0);
  try {
    (void)lossless_step(par, pp, 3, vm, own, 1);
    ADD_FAILURE() << "no exception";
  } catch (const std::logic_error& e) {
    EXPECT_STREQ(e.what(), "MaxFlow kernel: the replacement solve lost committed flow");
  }
}

TEST(MaxFlowKernel, RandomStepsAreCertifiedOnTheirSubset) {
  // Random DAG plans (nonuniform costs, parallel edges, a route that is a random subset of the
  // edges, one or two destinations), each placed in three steps on a changing view V = W + own (W
  // drawn below the capacity left beside the flow's own allocation) with a request of +inf, 0 or
  // a few units: every step certified on its member subset (the other edges masked), the small
  // graphs against the exhaustive optimum, the others against the independent solver. The same
  // steps with every cost multiplied give the same allocations and the curve's costs multiplied.
  std::mt19937_64 rng(5005);
  int positive = 0;   // steps that place
  int capped = 0;     // steps stopped by their request below the maximum
  int cancelled = 0;  // steps that take back part of the flow's own allocation
  for (int trial = 0; trial < 400; ++trial) {
    const bool tier_a = trial % 2 == 0;
    const int n = tier_a ? 3 + static_cast<int>(rng() % 4) : 6 + static_cast<int>(rng() % 9);
    const int m = tier_a ? 3 + static_cast<int>(rng() % 5) : n + static_cast<int>(rng() % (2 * n));
    const Amount max_cap = tier_a ? 2 : 40;
    const Cost max_cost = tier_a ? 6 : 40;
    const auto edges = random_dag_edges(
        rng, n, m, [&] { return static_cast<Amount>(rng() % static_cast<std::uint64_t>(max_cap + 1)); },
        [&] { return static_cast<Cost>(1 + rng() % static_cast<std::uint64_t>(max_cost)); });
    const auto g = make_graph(n, edges);
    auto scaled_edges = edges;
    for (auto& e : scaled_edges) e.cost *= 1000003;
    const auto scaled = make_graph(n, scaled_edges);
    std::vector<EdgeId> route;
    for (EdgeId e = 0; e < g.num_edges(); ++e) {
      if (rng() % 5 != 0) route.push_back(e);
    }
    std::vector<NodeId> dsts = {n - 1};
    if (rng() % 2 == 0) dsts.insert(dsts.begin(), n - 2);
    const auto plan = plan_of(g, FlowPlacement::MaxFlow, 0, dsts, route, caps(g));
    ASSERT_EQ(plan.members, walk_members(g, route, 0, dsts)) << "trial " << trial;  // every route edge on a walk
    const auto scaled_plan = plan_of(scaled, FlowPlacement::MaxFlow, 0, dsts, route, caps(scaled));
    ASSERT_EQ(scaled_plan.members, plan.members) << "trial " << trial;  // scaling keeps the EdgeId order
    std::vector<Amount> own(plan.members.size(), 0);
    Amount total = 0;
    for (int k = 0; k < 3; ++k) {
      std::vector<Amount> vm(plan.members.size());
      for (std::size_t j = 0; j < vm.size(); ++j) {
        const Amount cap = g.capacity_units()[static_cast<std::size_t>(plan.members[j])];
        vm[j] = own[j] + static_cast<Amount>(rng() % static_cast<std::uint64_t>(cap - own[j] + 1));
      }
      const std::uint64_t kind = rng() % 4;
      const std::optional<Amount> request =
          kind == 0 ? std::nullopt : std::optional<Amount>(static_cast<Amount>(kind == 1 ? 0 : rng() % 4));
      const auto step = lossless_step(g, plan, total, vm, own, request);
      const std::string where = "trial " + std::to_string(trial) + " step " + std::to_string(k);
      const auto fails = certify_step(g, plan, vm, total, request, step, tier_a);
      ASSERT_TRUE(fails.empty()) << where << ": " << fails.front();
      const auto other = lossless_step(scaled, scaled_plan, total, vm, own, request);
      ASSERT_EQ(other.total, step.total) << where;
      ASSERT_EQ(other.alloc, step.alloc) << where;
      ASSERT_EQ(other.cost_distribution.size(), step.cost_distribution.size()) << where;
      for (std::size_t i = 0; i < step.cost_distribution.size(); ++i) {
        ASSERT_EQ(other.cost_distribution[i].first, step.cost_distribution[i].first * 1000003) << where;
        ASSERT_EQ(other.cost_distribution[i].second, step.cost_distribution[i].second) << where;
      }
      positive += step.total > 0 ? 1 : 0;
      const Amount fmax = independent_curve(g, [&] {
        std::vector<Amount> view(static_cast<std::size_t>(g.num_edges()), 0);
        for (std::size_t j = 0; j < vm.size(); ++j) view[static_cast<std::size_t>(plan.members[j])] = vm[j];
        return view;
      }(), std::vector<NodeId> {0}, dsts, nullptr).back().first;
      capped += request.has_value() && step.total < fmax ? 1 : 0;
      bool took_back = false;
      for (std::size_t j = 0; j < own.size(); ++j) took_back = took_back || step.alloc[j] < own[j];
      cancelled += took_back ? 1 : 0;
      own = step.alloc;
      total = step.total;
    }
  }
  EXPECT_GT(positive, 240);
  EXPECT_GT(capped, 150);
  EXPECT_GT(cancelled, 30);
}

TEST(WeightedKernel, CapacityVersusResidualMembershipAfterOccupancy) {
  // Weighted membership by capacity versus residual after occupancy: parallel links [4, 4], another
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
  // Additivity on a fixed plan (place(1) then place(5) equals place(6)): weights 2 and 6 at the
  // source.
  auto g = make_graph(3, {{0, 1, 2, 1}, {0, 2, 6, 1}, {1, 2, 9, 1}});
  const auto plan = plan_of(g, FlowPlacement::Weighted, 0, {2}, all_ids(g), caps(g), EdgePriority::Capacity);
  const auto vm = on_members(plan, caps(g));
  std::vector<Amount> zero(plan.members.size(), 0);
  auto a = lossless_step(g, plan, 0, vm, zero, 1);
  auto b = lossless_step(g, plan, a.total, vm, a.alloc, 5);
  auto c = lossless_step(g, plan, 0, vm, zero, 6);
  // By hand, D'Hondt over the source weights [2, 6] (each unit to the larger w / (a + 1), ties to
  // the lower index): the first unit to 0 -> 2; 6 units give [1, 5] (the third unit is the tie
  // 2 / 1 = 6 / 3), and node 1 passes its unit on over 1 -> 2. Every amount fits its capacity,
  // so all 6 are placed.
  EXPECT_EQ(dense_input(g, plan, a.alloc), (std::vector<Amount>{0, 1, 0}));
  EXPECT_EQ(c.total, 6);
  EXPECT_EQ(dense_input(g, plan, c.alloc), (std::vector<Amount>{1, 5, 1}));
  EXPECT_EQ(b.total, c.total);
  EXPECT_EQ(b.alloc, c.alloc);
}

TEST(LossyKernel, DiamondDropsAndInfiniteInjection) {
  // A Lossy diamond: source links 10/10, tails 1/10.
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
  // [1, 3, 3].
  auto g2 = make_graph(2, {{0, 1, 1, 1}, {0, 1, 3, 1}, {0, 1, 3, 1}});
  const auto p2 = plan_of(g2, FlowPlacement::Weighted, 0, {1}, all_ids(g2), caps(g2), EdgePriority::Capacity,
                          Admission::Lossy);
  auto s2 = lossy_step(p2, 7, on_members(p2, caps(g2)));
  EXPECT_EQ(dense_input(g2, p2, s2.carried), (std::vector<Amount>{1, 3, 3}));
  EXPECT_EQ(s2.delivered, 7);
  // One link of capacity 0. EqualBalanced keeps it as a member: the offer enters it and is
  // dropped on it. Weighted gives it weight 0: the source has no member and injects nothing (the
  // offer is a source drop).
  auto g3 = make_graph(2, {{0, 1, 0, 1}});
  const auto p3 = plan_of(g3, FlowPlacement::EqualBalanced, 0, {1}, all_ids(g3), caps(g3), EdgePriority::Residual,
                          Admission::Lossy);
  ASSERT_EQ(p3.members, (std::vector<EdgeId>{0}));
  auto s3 = lossy_step(p3, 5, on_members(p3, caps(g3)));
  EXPECT_EQ(s3.offered, 5);
  EXPECT_EQ(s3.injected, 5);
  EXPECT_EQ(s3.delivered, 0);
  EXPECT_EQ(s3.carried, (std::vector<Amount>{0}));
  EXPECT_EQ(s3.drops, (std::vector<Amount>{5}));
  const auto w3 = plan_of(g3, FlowPlacement::Weighted, 0, {1}, all_ids(g3), caps(g3), EdgePriority::Capacity,
                          Admission::Lossy);
  ASSERT_TRUE(w3.members.empty());
  auto s4 = lossy_step(w3, 5, {});
  EXPECT_EQ(s4.offered, 5);
  EXPECT_EQ(s4.injected, 0);
  EXPECT_EQ(s4.delivered, 0);
}

// Lossy steps on random DAGs (parallel arcs, zero capacities, one or two destinations, a
// destination reached through another), both kernels, finite offers and +inf, against the
// per-node reference: every member's carried amount and drop, the injection, the delivery and
// the offer. The members are first checked against the walk definition, so the reference reads
// the same split graph.
TEST(LossyKernel, RandomStepsEqualThePerNodeReference) {
  std::mt19937 rng(1915);
  int weighted_transit_drops = 0;  // Weighted steps that drop at a node other than the source
  int balanced_transit_drops = 0;  // the same for EqualBalanced
  int zero_member_drops = 0;       // EqualBalanced members of capacity 0 that dropped their share
  int infinite = 0;                // +inf offers that injected something
  int delivering = 0;              // steps that delivered something
  for (int trial = 0; trial < 1000; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 6);
    const auto edges = random_dag_edges(
        rng, n, n + static_cast<int>(rng() % 10),
        [&] { return (rng() % 6) == 0 ? Amount {0} : 1 + static_cast<Amount>(rng() % 6); },
        [&] { return 1 + static_cast<Cost>(rng() % 3); });
    auto g = make_graph(n, edges);
    std::vector<NodeId> dsts = {static_cast<NodeId>(n - 1)};
    if (rng() % 3 == 0) dsts.insert(dsts.begin(), static_cast<NodeId>(1 + rng() % static_cast<std::uint32_t>(n - 2)));
    const bool weighted = trial % 2 == 1;
    const auto kernel = weighted ? FlowPlacement::Weighted : FlowPlacement::EqualBalanced;
    const auto plan = plan_of(g, kernel, 0, dsts, all_ids(g), caps(g), EdgePriority::Capacity, Admission::Lossy);
    // The split graph: every edge for EqualBalanced (a zero-capacity member takes its share and
    // drops it), the edges of positive weight (capacity) for Weighted.
    std::vector<EdgeId> usable;
    for (EdgeId e : all_ids(g)) {
      if (!weighted || g.capacity_units()[static_cast<std::size_t>(e)] >= 1) usable.push_back(e);
    }
    ASSERT_EQ(plan.members, walk_members(g, usable, 0, dsts)) << "trial " << trial;
    std::vector<std::pair<NodeId, NodeId>> arcs;
    std::vector<Amount> weights, residual;
    for (EdgeId e : plan.members) {
      const auto k = static_cast<std::size_t>(e);
      arcs.emplace_back(g.edge_src()[k], g.edge_dst()[k]);
      weights.push_back(g.capacity_units()[k]);
      residual.push_back(static_cast<Amount>(rng() % static_cast<std::uint64_t>(g.capacity_units()[k] + 1)));
    }
    const std::optional<Amount> offer =
        rng() % 4 == 0 ? std::nullopt : std::optional<Amount>(static_cast<Amount>(rng() % 25));
    const auto s = lossy_step(plan, offer, residual);
    const auto want = lossy_reference(n, arcs, weights, residual, weighted, offer, 0, dsts);
    const std::string where = "trial " + std::to_string(trial);
    EXPECT_EQ(s.carried, want.carried) << where;
    EXPECT_EQ(s.drops, want.drops) << where;
    EXPECT_EQ(s.injected, want.injected) << where;
    EXPECT_EQ(s.delivered, want.delivered) << where;
    EXPECT_EQ(s.offered, offer.value_or(want.injected)) << where;
    bool transit_drop = false;
    for (std::size_t i = 0; i < arcs.size(); ++i) transit_drop = transit_drop || (arcs[i].first != 0 && want.drops[i] > 0);
    (weighted ? weighted_transit_drops : balanced_transit_drops) += transit_drop ? 1 : 0;
    for (std::size_t i = 0; i < arcs.size(); ++i) {
      zero_member_drops += g.capacity_units()[static_cast<std::size_t>(plan.members[i])] == 0 && want.drops[i] > 0 ? 1 : 0;
    }
    infinite += !offer.has_value() && want.injected > 0 ? 1 : 0;
    delivering += want.delivered > 0 ? 1 : 0;
  }
  EXPECT_GT(weighted_transit_drops, 0);
  EXPECT_GT(balanced_transit_drops, 0);
  EXPECT_GT(zero_member_drops, 0);
  EXPECT_GT(infinite, 0);
  EXPECT_GT(delivering, 0);
}
