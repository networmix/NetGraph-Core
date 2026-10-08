// The CostOnly creation products in one pass (APPROACH 9.1, 10; task 1d-msd): cost_only_creation
// against the functions it replaces on the policy's path, kept as the reference (support_route's
// route, the record bundle of bundle_of, make_plan's plan and cost on that bundle), on random
// graphs for every source, selection and kernel: it answers exactly where that path's members are
// the whole route, with the same bundle (the support itself when the route is all of it), plan and
// cost. On crafted supports that break its preconditions (edges leaving a destination, dead ends,
// members below one unit, a cycle) it answers nothing or the reference's own result. Also
// lower_bound_from against std::lower_bound from every hint, and the plan's one block under copy
// and move.
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/traversal.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using namespace netgraph::core::detail;

namespace {

template <class T>
std::vector<std::remove_const_t<T>> vec(std::span<T> v) {
  return std::vector<std::remove_const_t<T>>(v.begin(), v.end());
}

void expect_same_plan(const KernelPlan& got, const KernelPlan& want, const std::string& where) {
  EXPECT_EQ(got.placement, want.placement) << where;
  EXPECT_EQ(got.admission, want.admission) << where;
  EXPECT_EQ(got.src, want.src) << where;
  EXPECT_EQ(vec(got.dsts), vec(want.dsts)) << where;
  EXPECT_EQ(got.route, want.route) << where;
  EXPECT_EQ(got.members, want.members) << where;
  EXPECT_EQ(got.weights, want.weights) << where;
  EXPECT_EQ(vec(got.topo), vec(want.topo)) << where;
  EXPECT_EQ(vec(got.out_begin), vec(want.out_begin)) << where;
  EXPECT_EQ(vec(got.out_end), vec(want.out_end)) << where;
  EXPECT_EQ(vec(got.head_pos), vec(want.head_pos)) << where;
  EXPECT_EQ(vec(got.terminal), vec(want.terminal)) << where;
  EXPECT_EQ(got.src_out_begin, want.src_out_begin) << where;
  EXPECT_EQ(got.src_out_end, want.src_out_end) << where;
}

// The policy's CostOnly creation before task 1d-msd (flow_policy.cpp: cost_only_route, bundle_of,
// new_flow's plan_on): the reference the one-pass products must equal.
struct Reference {
  std::vector<EdgeId> route;
  std::shared_ptr<const PredDAG> bundle;
  KernelPlan plan;
  Cost cost {kUnreachable};
};

Reference reference(const Graph& g, const std::shared_ptr<const PredDAG>& support, NodeId src,
                    const std::vector<NodeId>& dsts, const PathSelection& sel, PriorityFn view,
                    const std::vector<Amount>& W, FlowPlacement placement, Admission admission) {
  Reference out;
  auto r = support_route(g, PredDAGBuilder::ascending(*support), src, dsts, sel, view);
  out.route = r.edges;
  Cost cost = !sel.multipath && !r.edges.empty() ? r.cost : kUnreachable;
  out.bundle = PredDAGBuilder::ascending(*support).size() == r.edges.size()
                   ? support
                   : std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, r.edges));
  std::vector<Amount> creation;
  if (placement == FlowPlacement::Weighted && sel.priority == EdgePriority::Residual) {
    for (EdgeId e : PredDAGBuilder::ascending(*out.bundle)) creation.push_back(W[static_cast<std::size_t>(e)]);
  }
  const bool computed = cost == kUnreachable;
  out.plan = make_plan(g, placement, admission, sel.priority, src, dsts, *out.bundle, creation,
                       computed ? &cost : nullptr);
  out.cost = cost;
  return out;
}

struct Kernel {
  FlowPlacement placement;
  Admission admission;
};

const Kernel kKernels[] = {{FlowPlacement::MaxFlow, Admission::Lossless},
                           {FlowPlacement::EqualBalanced, Admission::Lossless},
                           {FlowPlacement::EqualBalanced, Admission::Lossy},
                           {FlowPlacement::Weighted, Admission::Lossless},
                           {FlowPlacement::Weighted, Admission::Lossy}};

std::vector<PathSelection> selections() {
  std::vector<PathSelection> out;
  for (auto tie : {EdgeTieBreak::Deterministic, EdgeTieBreak::PreferHigherPriority}) {
    for (auto prio : {EdgePriority::Capacity, EdgePriority::Residual, EdgePriority::Headroom}) {
      for (bool mp : {true, false}) {
        for (bool me : {true, false}) out.emplace_back(tie, prio, mp, me);
      }
    }
  }
  return out;
}

// Random edges with many ties (costs 1 to 3, times `factor`), parallel bursts, capacities 0 to 4
// units.
Graph random_graph(std::mt19937& rng, int n, int m, Cost factor = 1) {
  std::vector<InEdge> edges;
  while (static_cast<int>(edges.size()) < m) {
    const auto u = static_cast<NodeId>(rng() % static_cast<unsigned>(n));
    const auto v = static_cast<NodeId>(rng() % static_cast<unsigned>(n));
    const int copies = rng() % 4 == 0 ? 1 + static_cast<int>(rng() % 3) : 1;
    for (int c = 0; c < copies; ++c) {
      edges.push_back({u, v, static_cast<Amount>(rng() % 5), static_cast<Cost>(1 + rng() % 3) * factor});
    }
  }
  return make_graph(n, edges);
}

std::vector<NodeId> live(std::vector<NodeId> set, const Eligibility* elig) {
  std::sort(set.begin(), set.end());
  set.erase(std::unique(set.begin(), set.end()), set.end());
  std::vector<NodeId> out;
  for (NodeId v : set) {
    if (elig == nullptr || elig->node_ok(v)) out.push_back(v);
  }
  return out;
}

// The comparison on `trials` random graphs whose edge costs are 1 to 3 times `factor`: for every
// source, every one of the 24 selections and every legal kernel, the one pass answers exactly where
// the reference's members are its whole multipath route, with the same bundle (the support itself
// when the route is all of it), every plan array and the exact cost.
struct Tally {
  long made {0};
  long shared {0};
  long collapsed {0};
  long weighted {0};
  long anycast {0};
  long declined {0};
  long beyond_int32 {0};  // answers whose route cost exceeds INT32_MAX
};

Tally compare_with_reference(std::mt19937& rng, int trials, Cost factor) {
  const auto sels = selections();
  Tally t;
  for (int trial = 0; trial < trials; ++trial) {
    const int n = 2 + static_cast<int>(rng() % 10);
    const Graph g = random_graph(rng, n, static_cast<int>(rng() % 40), factor);
    std::optional<Eligibility> elig;
    if (rng() % 3 == 0) {
      std::vector<NodeId> xn;
      std::vector<EdgeId> xe;
      if (rng() % 2 == 0) xn.push_back(static_cast<NodeId>(rng() % static_cast<unsigned>(n)));
      for (int k = 0; k < 2 && g.num_edges() > 0; ++k) xe.push_back(static_cast<EdgeId>(rng() % static_cast<unsigned>(g.num_edges())));
      elig = g.eligibility(xn, xe);
    }
    const Eligibility* el = elig ? &*elig : nullptr;
    std::vector<NodeId> D;
    for (std::size_t k = 1 + rng() % 3; k > 0; --k) D.push_back(static_cast<NodeId>(rng() % static_cast<unsigned>(n)));
    const auto dset = live(D, nullptr);
    const auto dsts = live(D, el);
    if (dsts.empty()) continue;
    const auto support = std::make_shared<const PredDAG>(static_support(g, dset, Direction::Reverse, el));
    std::vector<Amount> W(static_cast<std::size_t>(g.num_edges()));
    for (std::size_t e = 0; e < W.size(); ++e) {
      const Amount cap = g.capacity_units()[e];
      W[e] = cap == 0 ? 0 : static_cast<Amount>(rng() % static_cast<std::uint32_t>(cap + 1));
    }
    for (NodeId src = 0; src < n; ++src) {
      if (std::binary_search(dset.begin(), dset.end(), src)) continue;
      for (const auto& sel : sels) {
        const PriorityValues values {&g, sel.priority, W};
        const PriorityFn view = PriorityFn::of(values);
        for (const auto& k : kKernels) {
          if (k.placement == FlowPlacement::Weighted && sel.priority == EdgePriority::Headroom) continue;
          const std::string where = "factor " + std::to_string(factor) + " trial " + std::to_string(trial) + " src " + std::to_string(src) + " kernel " +
                                    std::to_string(static_cast<int>(k.placement)) + "/" +
                                    std::to_string(static_cast<int>(k.admission)) + " sel " +
                                    std::to_string(static_cast<int>(sel.tie_break)) +
                                    std::to_string(static_cast<int>(sel.priority)) + std::to_string(sel.multipath) +
                                    std::to_string(sel.multi_edge);
          const auto want = reference(g, support, src, dsts, sel, view, W, k.placement, k.admission);
          const auto got = cost_only_creation(g, support, src, dsts, sel, view, k.placement, k.admission);
          // It answers exactly where the reference's members are its whole multipath route under a
          // kernel of static membership (on a valid support that is every route that exists).
          const bool applies = sel.multipath &&
                               !(k.placement == FlowPlacement::Weighted && sel.priority == EdgePriority::Residual) &&
                               !want.route.empty();
          EXPECT_EQ(got.has_value(), applies) << where;
          if (got.has_value() != applies) return t;
          if (!got.has_value()) {
            ++t.declined;
            continue;
          }
          EXPECT_EQ(want.plan.members, want.route) << where;  // the precondition the shortcut checks
          EXPECT_TRUE(got->bundle) << where;
          if (!got->bundle) return t;
          EXPECT_TRUE(*got->bundle == *want.bundle) << where;
          EXPECT_EQ(vec(PredDAGBuilder::ascending(*got->bundle)), want.route) << where;
          EXPECT_EQ(got->bundle == support, want.bundle == support) << where;  // the same shared value
          expect_same_plan(got->plan, want.plan, where);
          EXPECT_EQ(got->cost, want.cost) << where;
          ++t.made;
          if (want.bundle == support) ++t.shared;
          if (!sel.multi_edge && want.route.size() < PredDAGBuilder::ascending(*support).size()) ++t.collapsed;
          if (k.placement == FlowPlacement::Weighted) ++t.weighted;
          if (dsts.size() > 1) ++t.anycast;
          if (want.cost > std::numeric_limits<std::int32_t>::max()) ++t.beyond_int32;
          if (::testing::Test::HasFailure()) return t;
        }
      }
    }
  }
  return t;
}

} // namespace

TEST(CostOnlyCreation, LowerBoundFromEveryHintIsStdLowerBound) {
  std::mt19937 rng(20261007);
  int checked = 0;
  for (int trial = 0; trial < 400; ++trial) {
    const auto n = static_cast<std::size_t>(rng() % 140);
    std::vector<EdgeId> a(n);
    EdgeId x = static_cast<EdgeId>(rng() % 3);
    for (auto& v : a) {
      v = x;
      x += 1 + static_cast<EdgeId>(rng() % 4);
    }
    const EdgeId top = n == 0 ? 3 : a.back() + 3;
    for (EdgeId key = -2; key <= top; ++key) {
      const auto want = static_cast<std::size_t>(std::lower_bound(a.begin(), a.end(), key) - a.begin());
      for (std::size_t hint = 0; hint <= n + 3; ++hint) {
        ASSERT_EQ(lower_bound_from(a, hint, key), want) << "n " << n << " key " << key << " hint " << hint;
        ++checked;
      }
    }
  }
  EXPECT_GT(checked, 100000);
}

TEST(CostOnlyCreation, EqualsTheReferenceForEverySourceSelectionAndKernel) {
  std::mt19937 rng(1007);
  const Tally t = compare_with_reference(rng, 300, 1);
  // Not vacuous: every branch of the one pass was compared, and the declined inputs too.
  EXPECT_GT(t.made, 20000);
  EXPECT_GT(t.shared, 500);
  EXPECT_GT(t.collapsed, 2000);
  EXPECT_GT(t.weighted, 4000);
  EXPECT_GT(t.anycast, 4000);
  EXPECT_GT(t.declined, 20000);
}

TEST(CostOnlyCreation, WideCostsEqualTheReference) {
  // Repeat the reference comparison with factors 3, 64 and 2^35. These cases catch the
  // int32 output truncation; with 2^35 every answered cost exceeds INT32_MAX.
  // They do not cover every narrower representation: a 48-bit output mask survives.
  for (const Cost factor : {Cost {3}, Cost {64}, Cost {1} << 35}) {
    std::mt19937 rng(static_cast<std::uint32_t>(2000 + factor % 997));
    const Tally t = compare_with_reference(rng, 120, factor);
    if (::testing::Test::HasFailure()) return;
    EXPECT_GT(t.made, 7000) << "factor " << factor;
    EXPECT_GT(t.shared, 150) << "factor " << factor;
    EXPECT_GT(t.collapsed, 700) << "factor " << factor;
    EXPECT_GT(t.declined, 7000) << "factor " << factor;
    if (factor == (Cost {1} << 35)) EXPECT_EQ(t.beyond_int32, t.made);
  }
}

TEST(CostOnlyCreation, ChainCostsBeyondInt32AreExact) {
  // The review's witness (0 -> 1 -> 2, two edges of 2^35: the cost 2^36 = 68,719,476,736), a chain
  // whose total 2,500,000,000 exceeds INT32_MAX without being a power of two, a diamond of 2^35
  // edges with a parallel pair, and a full-width chain of two edges of 2^60 + 1 (the cost
  // 2^61 + 2, below the graph's 2^62 bound; a 48-bit output mask gives 2): the exact cost and every
  // array, for every legal kernel, both multi_edge settings.
  struct Chain {
    std::int32_t n;
    std::vector<InEdge> edges;
    NodeId dst;
    Cost cost;
  };
  const Cost big = Cost {1} << 35;
  const std::vector<Chain> chains = {
      {3, {{0, 1, 3, big}, {1, 2, 3, big}}, 2, Cost {68719476736}},
      {3, {{0, 1, 5, 1500000000}, {1, 2, 5, 1000000000}}, 2, Cost {2500000000}},
      {4, {{0, 1, 4, big}, {0, 2, 3, big}, {1, 3, 2, big}, {2, 3, 1, big}, {2, 3, 2, big}}, 3, 2 * big},
      {3, {{0, 1, 2, (Cost {1} << 60) + 1}, {1, 2, 2, (Cost {1} << 60) + 1}}, 2, Cost {2305843009213693954}},
  };
  int checked = 0;
  for (const auto& c : chains) {
    ASSERT_GT(c.cost, std::numeric_limits<std::int32_t>::max());
    const Graph g = make_graph(c.n, c.edges);
    const std::vector<NodeId> dsts = {c.dst};
    const auto support = std::make_shared<const PredDAG>(static_support(g, dsts, Direction::Reverse, nullptr));
    const std::vector<Amount> W(g.capacity_units().begin(), g.capacity_units().end());
    for (bool me : {true, false}) {
      const PathSelection sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Capacity, true, me);
      const PriorityValues values {&g, sel.priority, W};
      for (const auto& k : kKernels) {
        const std::string where = "chain cost " + std::to_string(c.cost) + " kernel " +
                                  std::to_string(static_cast<int>(k.placement)) + "/" +
                                  std::to_string(static_cast<int>(k.admission)) + " multi_edge " + std::to_string(me);
        const auto want = reference(g, support, 0, dsts, sel, PriorityFn::of(values), W, k.placement, k.admission);
        const auto got =
            cost_only_creation(g, support, 0, dsts, sel, PriorityFn::of(values), k.placement, k.admission);
        ASSERT_TRUE(got.has_value()) << where;
        EXPECT_EQ(want.cost, c.cost) << where;
        EXPECT_EQ(got->cost, c.cost) << where;
        EXPECT_TRUE(*got->bundle == *want.bundle) << where;
        EXPECT_EQ(vec(PredDAGBuilder::ascending(*got->bundle)), want.route) << where;
        EXPECT_EQ(got->bundle == support, want.bundle == support) << where;
        expect_same_plan(got->plan, want.plan, where);
        ++checked;
      }
    }
  }
  EXPECT_EQ(checked, 4 * 2 * 5);
}

TEST(CostOnlyCreation, LargerGraphsAndTheSupportItself) {
  // Larger fabrics (Clos and random) where the walk crosses many ranges: every source to one and
  // to two destinations, the EB and MaxFlow kernels, with and without Collapse.
  std::mt19937 rng(7);
  long made = 0;
  for (int shape = 0; shape < 4; ++shape) {
    Graph g = [&] {
      if (shape < 2) {
        std::vector<InEdge> edges;
        const int leaves = shape == 0 ? 12 : 20;
        const int spines = shape == 0 ? 6 : 8;
        for (int u = 0; u < leaves; ++u) {
          for (int s = 0; s < spines; ++s) {
            const Amount cap = static_cast<Amount>(1 + rng() % 9);
            for (int copy = 0; copy < (rng() % 5 == 0 ? 2 : 1); ++copy) {
              edges.push_back({u, leaves + s, cap, 1});
              edges.push_back({leaves + s, u, cap, 1});
            }
          }
        }
        return make_graph(leaves + spines, edges);
      }
      return random_graph(rng, 200, 2400);
    }();
    const int n = g.num_nodes();
    std::vector<Amount> W(g.capacity_units().begin(), g.capacity_units().end());
    for (int q = 0; q < 6; ++q) {
      std::vector<NodeId> D = {static_cast<NodeId>(rng() % static_cast<unsigned>(n))};
      if (q % 2 == 1) D.push_back(static_cast<NodeId>(rng() % static_cast<unsigned>(n)));
      const auto dsts = live(D, nullptr);
      const auto support = std::make_shared<const PredDAG>(static_support(g, dsts, Direction::Reverse, nullptr));
      for (NodeId src = 0; src < n; src += (n > 50 ? 7 : 1)) {
        if (std::binary_search(dsts.begin(), dsts.end(), src)) continue;
        for (bool me : {true, false}) {
          const PathSelection sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Capacity, true, me);
          const PriorityValues values {&g, sel.priority, W};
          for (auto placement : {FlowPlacement::EqualBalanced, FlowPlacement::MaxFlow}) {
            const auto want = reference(g, support, src, dsts, sel, PriorityFn::of(values), W, placement,
                                        Admission::Lossless);
            const auto got = cost_only_creation(g, support, src, dsts, sel, PriorityFn::of(values), placement,
                                                Admission::Lossless);
            ASSERT_EQ(got.has_value(), !want.route.empty());
            if (!got) continue;
            EXPECT_TRUE(*got->bundle == *want.bundle);
            EXPECT_EQ(got->bundle == support, want.bundle == support);
            expect_same_plan(got->plan, want.plan, "shape " + std::to_string(shape) + " src " + std::to_string(src));
            EXPECT_EQ(got->cost, want.cost);
            ++made;
          }
        }
      }
    }
  }
  EXPECT_GT(made, 300);
}

TEST(CostOnlyCreation, CraftedSupportsGiveNothingOrTheReference) {
  // Acyclic edge sets that are no static support: edges leaving a destination, nodes that reach no
  // destination, members below one unit. The one pass declines them or equals the reference.
  std::mt19937 rng(99);
  long made = 0;
  long declined = 0;
  for (int trial = 0; trial < 1500; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 8);
    const Graph g = random_graph(rng, n, 4 + static_cast<int>(rng() % 30));
    std::vector<NodeId> perm(static_cast<std::size_t>(n));
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    std::vector<int> rank(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) rank[static_cast<std::size_t>(perm[static_cast<std::size_t>(i)])] = i;
    std::vector<EdgeId> chosen;
    for (EdgeId e = 0; e < g.num_edges(); ++e) {
      const auto u = g.edge_src()[static_cast<std::size_t>(e)];
      const auto v = g.edge_dst()[static_cast<std::size_t>(e)];
      if (rank[static_cast<std::size_t>(u)] < rank[static_cast<std::size_t>(v)] && rng() % 4 != 0) chosen.push_back(e);
    }
    const auto support = std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, chosen));
    const NodeId src = perm.front();
    std::vector<NodeId> dsts = {perm.back()};
    if (n > 4 && rng() % 2 == 0) dsts.push_back(perm[static_cast<std::size_t>(n / 2)]);  // may have out-edges
    std::sort(dsts.begin(), dsts.end());
    std::vector<Amount> W(g.capacity_units().begin(), g.capacity_units().end());
    for (bool me : {true, false}) {
      const PathSelection sel(EdgeTieBreak::Deterministic, EdgePriority::Capacity, true, me);
      const PriorityValues values {&g, sel.priority, W};
      for (const auto& k : kKernels) {
        const auto want = reference(g, support, src, dsts, sel, PriorityFn::of(values), W, k.placement, k.admission);
        const auto got =
            cost_only_creation(g, support, src, dsts, sel, PriorityFn::of(values), k.placement, k.admission);
        const std::string where = "trial " + std::to_string(trial);
        if (!got.has_value()) {
          ++declined;
          continue;
        }
        EXPECT_TRUE(*got->bundle == *want.bundle) << where;
        EXPECT_EQ(got->bundle == support, want.bundle == support) << where;
        expect_same_plan(got->plan, want.plan, where);
        EXPECT_EQ(got->cost, want.cost) << where;
        ++made;
        if (::testing::Test::HasFailure()) return;
      }
    }
  }
  EXPECT_GT(made, 1000);
  EXPECT_GT(declined, 1000);
}

TEST(CostOnlyCreation, DeclinesEachBrokenPrecondition) {
  // 0 -> 1 -> 3 and 0 -> 2 -> 3, destination 3: the reference keeps the whole route, and so does
  // the one pass. Each variant breaks one precondition; the one pass then answers nothing (the
  // reference's members differ from its route, or it throws).
  const Graph g = make_graph(6, {{0, 1, 2, 1}, {0, 2, 2, 1}, {1, 3, 2, 1}, {2, 3, 2, 1}, {3, 4, 2, 1},
                                 {1, 5, 2, 1}, {2, 1, 2, 1}, {1, 2, 0, 1}});
  const auto inv = id_of_input(g);
  const std::vector<NodeId> dsts = {3};
  const std::vector<Amount> W(g.capacity_units().begin(), g.capacity_units().end());
  const PathSelection sel(EdgeTieBreak::Deterministic, EdgePriority::Capacity, true, true);
  const PriorityValues values {&g, sel.priority, W};
  auto run = [&](std::vector<EdgeId> edges, FlowPlacement placement) {
    for (auto& e : edges) e = inv[static_cast<std::size_t>(e)];
    std::sort(edges.begin(), edges.end());
    const auto support = std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, edges));
    return std::make_pair(support, cost_only_creation(g, support, 0, dsts, sel, PriorityFn::of(values), placement,
                                                      Admission::Lossless));
  };
  {
    const auto [support, got] = run({0, 1, 2, 3}, FlowPlacement::EqualBalanced);
    ASSERT_TRUE(got.has_value());
    EXPECT_EQ(got->bundle, support);  // the route is the whole support
    EXPECT_EQ(got->cost, 2);
    EXPECT_EQ(vec(got->plan.topo), (std::vector<NodeId> {0, 1, 2, 3}));
  }
  // An edge leaving the destination (3 -> 4).
  EXPECT_FALSE(run({0, 1, 2, 3, 4}, FlowPlacement::EqualBalanced).second.has_value());
  // A node that reaches no destination (1 -> 5).
  EXPECT_FALSE(run({0, 1, 2, 3, 5}, FlowPlacement::EqualBalanced).second.has_value());
  // A member below one unit (1 -> 2 has no capacity): declined for EB, kept by MaxFlow.
  EXPECT_FALSE(run({0, 1, 2, 3, 7}, FlowPlacement::EqualBalanced).second.has_value());
  EXPECT_TRUE(run({0, 1, 2, 3, 7}, FlowPlacement::MaxFlow).second.has_value());
  // A cycle (1 -> 2 -> 1): declined, and the reference throws.
  {
    const auto [support, got] = run({0, 1, 2, 3, 6, 7}, FlowPlacement::MaxFlow);
    EXPECT_FALSE(got.has_value());
    EXPECT_THROW((void)reference(g, support, 0, dsts, sel, PriorityFn::of(values), W, FlowPlacement::MaxFlow,
                                 Admission::Lossless),
                 std::logic_error);
  }
  // Not applicable: multipath = false, Weighted with Residual priority, no destination reached.
  {
    std::vector<EdgeId> edges = {inv[0], inv[1], inv[2], inv[3]};
    std::sort(edges.begin(), edges.end());
    const auto support = std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, edges));
    const PathSelection single(EdgeTieBreak::Deterministic, EdgePriority::Capacity, false, true);
    EXPECT_FALSE(cost_only_creation(g, support, 0, dsts, single, PriorityFn::of(values), FlowPlacement::EqualBalanced,
                                    Admission::Lossless)
                     .has_value());
    const PathSelection residual(EdgeTieBreak::Deterministic, EdgePriority::Residual, true, true);
    EXPECT_FALSE(cost_only_creation(g, support, 0, dsts, residual, PriorityFn::of(values), FlowPlacement::Weighted,
                                    Admission::Lossless)
                     .has_value());
    EXPECT_FALSE(cost_only_creation(g, support, 0, std::vector<NodeId> {5}, sel, PriorityFn::of(values),
                                    FlowPlacement::EqualBalanced, Admission::Lossless)
                     .has_value());
  }
}

TEST(KernelPlanBlock, CopiesAndMovesKeepTheirOwnViews) {
  const Graph g = make_graph(4, {{0, 1, 3, 1}, {0, 2, 2, 1}, {1, 3, 3, 1}, {2, 3, 2, 1}});
  std::vector<EdgeId> all(4);
  std::iota(all.begin(), all.end(), 0);
  const auto dag = PredDAGBuilder::of_edges(g, all);
  const std::vector<NodeId> dsts = {3};
  const KernelPlan original =
      make_plan(g, FlowPlacement::Weighted, Admission::Lossless, EdgePriority::Capacity, 0, dsts, dag, {});
  ASSERT_EQ(original.topo.size(), 4U);
  ASSERT_EQ(original.head_pos.size(), 4U);
  ASSERT_EQ(vec(original.dsts), dsts);
  // A copy's views point into its own block (none into the source's), and stay valid once the
  // source is gone.
  auto own_block = [](const KernelPlan& a, const KernelPlan& b) {
    return a.dsts.data() != b.dsts.data() && a.topo.data() != b.topo.data() &&
           a.out_begin.data() != b.out_begin.data() && a.out_end.data() != b.out_end.data() &&
           a.head_pos.data() != b.head_pos.data() && a.terminal.data() != b.terminal.data();
  };
  std::optional<KernelPlan> source(original);
  const KernelPlan copy(*source);
  EXPECT_TRUE(own_block(copy, *source));
  source.reset();
  expect_same_plan(copy, original, "copy");
  KernelPlan donor(copy);
  KernelPlan moved(std::move(donor));
  expect_same_plan(moved, original, "move");
  EXPECT_TRUE(donor.topo.empty() && donor.dsts.empty() && donor.head_pos.empty() && donor.terminal.empty());
  KernelPlan assigned;
  assigned = copy;
  expect_same_plan(assigned, original, "copy assignment");
  EXPECT_TRUE(own_block(assigned, copy));
  KernelPlan target = make_plan(g, FlowPlacement::EqualBalanced, Admission::Lossless, EdgePriority::Capacity, 1, dsts,
                                dag, {});
  target = std::move(moved);
  expect_same_plan(target, original, "move assignment");
  EXPECT_TRUE(moved.topo.empty() && moved.dsts.empty() && moved.head_pos.empty() && moved.terminal.empty());
  const KernelPlan& self = target;
  target = self;
  expect_same_plan(target, original, "self assignment");
  const KernelPlan empty;
  const KernelPlan empty_copy(empty);
  EXPECT_TRUE(empty_copy.dsts.empty() && empty_copy.topo.empty() && empty_copy.members.empty());
}
