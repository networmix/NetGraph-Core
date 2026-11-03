// Roster route selection and multi-root DAG costs against independent walks; plan ownership.
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "netgraph/core/detail/route_kernels.hpp"
#include "netgraph/core/detail/route_builders.hpp"
#include "netgraph/core/detail/search_detail.hpp"
#include "search_reference.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using namespace netgraph::core::detail;

namespace {

// Every field of two plans equal (plan_differences, test_utils.hpp).
void expect_same_plan(const PreparedRoute& got, const PreparedRoute& want, const std::string& where) {
  EXPECT_EQ(plan_differences(got, want), std::vector<std::string> {}) << where;
}

// The 13 route choices a ShortestPaths configuration accepts: every tie_break with multipath and
// parallel_edges, but no HIGHER_* tie_break with both true.
std::vector<RouteChoice> choices() {
  std::vector<RouteChoice> out;
  for (auto tie : {RouteTieBreak::FewerHops, RouteTieBreak::HigherCapacity, RouteTieBreak::HigherAvailableCapacity,
                   RouteTieBreak::HigherHeadroom}) {
    for (bool mp : {true, false}) {
      for (bool me : {true, false}) {
        if (tie != RouteTieBreak::FewerHops && mp && me) continue;
        out.push_back(RouteChoice {tie, mp, me});
      }
    }
  }
  return out;
}

// Random edges with many ties (costs 1 to 3, times `factor`), parallel bursts, capacities 0 to 4
// units.
Graph tied_graph(std::mt19937& rng, int n, int m, Cost factor = 1) {
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

} // namespace

TEST(ShortestPathsCreation, LowerBoundFromEveryHintIsStdLowerBound) {
  std::mt19937 rng(20261007);
  int forward = 0;   // the answer above the hint
  int backward = 0;  // the answer below the hint
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
        forward += want > hint ? 1 : 0;
        backward += want < hint ? 1 : 0;
      }
    }
  }
  EXPECT_GT(forward, 1000000);
  EXPECT_GT(backward, 1000000);
}

TEST(KernelPlanBlock, CopiesAndMovesKeepTheirOwnViews) {
  const Graph g = make_graph(4, {{0, 1, 3, 1}, {0, 2, 2, 1}, {1, 3, 3, 1}, {2, 3, 2, 1}});
  const auto dag = PredecessorDAGBuilder::of_edges(g, all_ids(g));
  const std::vector<NodeId> dsts = {3};
  const PreparedRoute original =
      prepare_route(g, RouteAllocation::Weighted, Admission::Lossless, WeightSource {}, 0, dsts, dag);
  ASSERT_EQ(original.topological_order.size(), 4U);
  ASSERT_EQ(original.edge_destination_positions.size(), 4U);
  ASSERT_EQ(vec(original.destinations), dsts);
  // A copy's views point into its own block (none into the source's), and stay valid once the
  // source is gone.
  auto own_block = [](const PreparedRoute& a, const PreparedRoute& b) {
    return a.sources.data() != b.sources.data() && a.source_positions.data() != b.source_positions.data() &&
           a.destinations.data() != b.destinations.data() && a.topological_order.data() != b.topological_order.data() &&
           a.outgoing_begin.data() != b.outgoing_begin.data() && a.outgoing_end.data() != b.outgoing_end.data() &&
           a.edge_destination_positions.data() != b.edge_destination_positions.data() && a.is_destination.data() != b.is_destination.data();
  };
  std::optional<PreparedRoute> source(original);
  const PreparedRoute copy(*source);
  EXPECT_TRUE(own_block(copy, *source));
  source.reset();
  expect_same_plan(copy, original, "copy");
  PreparedRoute donor(copy);
  PreparedRoute moved(std::move(donor));
  expect_same_plan(moved, original, "move");
  EXPECT_TRUE(donor.topological_order.empty() && donor.destinations.empty() && donor.edge_destination_positions.empty() && donor.is_destination.empty() &&
              donor.sources.empty() && donor.source_positions.empty());
  PreparedRoute assigned;
  assigned = copy;
  expect_same_plan(assigned, original, "copy assignment");
  EXPECT_TRUE(own_block(assigned, copy));
  PreparedRoute target = prepare_route(g, RouteAllocation::Equal, Admission::Lossless, WeightSource {}, 1, dsts, dag);
  target = std::move(moved);
  expect_same_plan(target, original, "move assignment");
  EXPECT_TRUE(moved.topological_order.empty() && moved.destinations.empty() && moved.edge_destination_positions.empty() && moved.is_destination.empty() &&
              moved.sources.empty() && moved.source_positions.empty());
  const PreparedRoute& self = target;
  target = self;
  expect_same_plan(target, original, "self assignment");
  const PreparedRoute empty;
  const PreparedRoute empty_copy(empty);
  EXPECT_TRUE(empty_copy.destinations.empty() && empty_copy.topological_order.empty() && empty_copy.forwarding_edges.empty() &&
              empty_copy.sources.empty() && empty_copy.source_positions.empty());
}

// A lifetime witness: two identical chains, 0 -> 1 -> 2 and 3 -> 4 -> 5, give plans of one shape
// (one root, one destination, three positions, two members) and other values. A copy (and a copy
// assignment) of the first chain's plan keeps its roots, root positions and destinations after the
// plan it was copied from is destroyed and plans of the same shape for the second chain have been
// built (they reuse freed blocks of that size): a copy whose views still pointed into the
// destroyed block would read the second chain's root 3 and destination 5.
TEST(KernelPlanBlock, ACopyOutlivesItsSource) {
  const Graph g = make_graph(6, {{0, 1, 2, 1}, {1, 2, 2, 1}, {3, 4, 2, 1}, {4, 5, 2, 1}});
  auto chain = [&](NodeId src, NodeId dst) {
    const std::vector<EdgeId> edges = src == 0 ? std::vector<EdgeId> {0, 1} : std::vector<EdgeId> {2, 3};
    const auto dag = PredecessorDAGBuilder::of_edges(g, edges);
    return prepare_route(g, RouteAllocation::Equal, Admission::Lossless, WeightSource {}, src,
                     std::vector<NodeId> {dst}, dag);
  };
  for (const bool assign : {false, true}) {
    std::optional<PreparedRoute> owner(chain(0, 2));
    ASSERT_EQ(vec(owner->sources), (std::vector<NodeId> {0}));
    ASSERT_EQ(vec(owner->source_positions), (std::vector<std::int32_t> {0}));
    std::optional<PreparedRoute> copy;
    if (assign) {
      copy.emplace(chain(3, 5));
      *copy = *owner;
    } else {
      copy.emplace(*owner);
    }
    owner.reset();
    std::vector<PreparedRoute> later;
    for (int i = 0; i < 64; ++i) later.push_back(chain(3, 5));
    const std::string where = assign ? "copy assignment" : "copy construction";
    EXPECT_EQ(vec(copy->sources), (std::vector<NodeId> {0})) << where;
    EXPECT_EQ(vec(copy->source_positions), (std::vector<std::int32_t> {0})) << where;
    EXPECT_EQ(vec(copy->destinations), (std::vector<NodeId> {2})) << where;
    EXPECT_EQ(vec(copy->topological_order), (std::vector<NodeId> {0, 1, 2})) << where;
  }
}

// The plan comparison names every field that differs: a copy of a real plan with one field changed
// differs in exactly that field (a vector or view field both by its last entry and by one more
// entry), and without the route only the route may differ.
TEST(TestHelpers, PlanDifferencesNameEveryField) {
  const Graph g = make_graph(4, {{0, 1, 3, 1}, {0, 2, 2, 1}, {1, 3, 3, 1}, {2, 3, 2, 1}});
  const auto dag = PredecessorDAGBuilder::of_edges(g, all_ids(g));
  const std::vector<NodeId> dsts = {3};
  const PreparedRoute plan =
      prepare_route(g, RouteAllocation::Weighted, Admission::Lossless, WeightSource {}, 0, dsts, dag);
  ASSERT_EQ(plan.weights.size(), 4U);
  EXPECT_EQ(plan_differences(plan, plan), std::vector<std::string> {});
  // The views' changed copies live here, as the plan's own views live in its block.
  std::vector<NodeId> nodes;
  std::vector<std::int32_t> positions;
  std::vector<std::uint8_t> flags;
  using Change = std::function<void(PreparedRoute&)>;
  auto nodes_of = [&](std::span<const NodeId> x, bool longer) {
    nodes.assign(x.begin(), x.end());
    if (longer) nodes.push_back(0);
    else ++nodes.back();
    return std::span<const NodeId>(nodes);
  };
  auto positions_of = [&](std::span<const std::int32_t> x, bool longer) {
    positions.assign(x.begin(), x.end());
    if (longer) positions.push_back(0);
    else ++positions.back();
    return std::span<const std::int32_t>(positions);
  };
  const std::vector<std::pair<std::string, Change>> changes = {
      {"allocation", [](auto& x) { x.allocation = RouteAllocation::Equal; }},
      {"sources", [&](auto& x) { x.sources = nodes_of(x.sources, false); }},
      {"sources", [&](auto& x) { x.sources = nodes_of(x.sources, true); }},
      {"source_positions", [&](auto& x) { x.source_positions = positions_of(x.source_positions, false); }},
      {"source_positions", [&](auto& x) { x.source_positions = positions_of(x.source_positions, true); }},
      {"destinations", [&](auto& x) { x.destinations = nodes_of(x.destinations, false); }},
      {"destinations", [&](auto& x) { x.destinations = nodes_of(x.destinations, true); }},
      {"support_edges", [](auto& x) { x.support_edges.push_back(0); }},
      {"support_edges", [](auto& x) { ++x.support_edges.back(); }},
      {"forwarding_edges", [](auto& x) { x.forwarding_edges.push_back(0); }},
      {"forwarding_edges", [](auto& x) { ++x.forwarding_edges.back(); }},
      {"weights", [](auto& x) { x.weights.push_back(1); }},
      {"weights", [](auto& x) { ++x.weights.back(); }},
      {"topological_order", [&](auto& x) { x.topological_order = nodes_of(x.topological_order, false); }},
      {"topological_order", [&](auto& x) { x.topological_order = nodes_of(x.topological_order, true); }},
      {"outgoing_begin", [&](auto& x) { x.outgoing_begin = positions_of(x.outgoing_begin, false); }},
      {"outgoing_begin", [&](auto& x) { x.outgoing_begin = positions_of(x.outgoing_begin, true); }},
      {"outgoing_end", [&](auto& x) { x.outgoing_end = positions_of(x.outgoing_end, false); }},
      {"outgoing_end", [&](auto& x) { x.outgoing_end = positions_of(x.outgoing_end, true); }},
      {"edge_destination_positions", [&](auto& x) { x.edge_destination_positions = positions_of(x.edge_destination_positions, false); }},
      {"edge_destination_positions", [&](auto& x) { x.edge_destination_positions = positions_of(x.edge_destination_positions, true); }},
      {"is_destination", [&](auto& x) {
         flags.assign(x.is_destination.begin(), x.is_destination.end());
         flags.back() = flags.back() == 0 ? 1 : 0;
         x.is_destination = flags;
       }},
      {"is_destination", [&](auto& x) {
         flags.assign(x.is_destination.begin(), x.is_destination.end());
         flags.push_back(0);
         x.is_destination = flags;
       }},
  };
  for (std::size_t k = 0; k < changes.size(); ++k) {
    PreparedRoute x = plan;
    changes[k].second(x);
    EXPECT_EQ(plan_differences(x, plan), std::vector<std::string> {changes[k].first}) << k;
    const auto without_route = changes[k].first == "support_edges" ? std::vector<std::string> {}
                                                           : std::vector<std::string> {changes[k].first};
    EXPECT_EQ(plan_differences(x, plan, false), without_route) << k;
    // expect_same_plan reports it.
    EXPECT_EQ(failures_of([&] { expect_same_plan(x, plan, "changed"); }), 1) << k;
  }
  EXPECT_EQ(failures_of([&] { expect_same_plan(plan, plan, "same"); }), 0);
}

namespace {

// Independent walk enumeration: no production prune, cost pass, route selection or plan builder.
std::vector<EdgeId> walk_union(const Graph& g, const std::vector<char>& members,
                               const std::vector<NodeId>& roots, const std::vector<NodeId>& dsts) {
  std::set<EdgeId> edges;
  for (const auto& walk : ref::all_walks(g, members, roots, dsts)) edges.insert(walk.begin(), walk.end());
  return {edges.begin(), edges.end()};
}

std::vector<Cost> walk_costs(const Graph& g, const std::vector<char>& members,
                             const std::vector<NodeId>& roots, const std::vector<NodeId>& dsts) {
  std::vector<Cost> out;
  for (NodeId root : roots) {
    Cost cost = kUnreachable;
    for (const auto& walk : ref::all_walks(g, members, {root}, dsts)) cost = std::min(cost, ref::walk_cost(g, walk));
    out.push_back(cost);
  }
  return out;
}

std::vector<char> flags_of(const Graph& g, const std::vector<EdgeId>& edges) {
  std::vector<char> out(static_cast<std::size_t>(g.num_edges()), 0);
  for (EdgeId e : edges) out[static_cast<std::size_t>(e)] = 1;
  return out;
}

// Collapse by the reference's per-pair rule, then enumerate each source's paths on the tight DAG.
// Every walk there is shortest to the destination set. Canonical path order is derived by the
// reference's (cost, bottleneck, hops, EdgeId sequence) comparison of whole walks.
std::vector<EdgeId> shortest_union_reference(const Graph& g, const std::vector<EdgeId>& support,
    const std::vector<NodeId>& roots, const std::vector<NodeId>& dsts,
    const RouteChoice& choice, const std::vector<Amount>& view) {
  ref::Query q {&g, nullptr, &view, 0, choice.parallel_edges, choice.tie_break};
  auto members = flags_of(g, support);
  if (!choice.parallel_edges) {
    // Compare only allowed tight edges: a wider parallel edge excluded by the creation mask
    // must not win the choice and then erase the surviving allowed edge.
    for (EdgeId e : support) {
      for (EdgeId other : support) {
        const auto a = static_cast<std::size_t>(e), b = static_cast<std::size_t>(other);
        if (g.edge_sources()[a] != g.edge_sources()[b] || g.edge_destinations()[a] != g.edge_destinations()[b]) continue;
        const int score = choice.prefers_score() ? ref::cmp(ref::prio_of(q, other), ref::prio_of(q, e)) : 0;
        if (score > 0 || (score == 0 && other < e)) members[a] = 0;
      }
    }
  }
  if (!choice.single_path()) return walk_union(g, members, roots, dsts);
  std::set<EdgeId> edges;
  for (NodeId root : roots) {
    const auto path = ref::canonical(q, members, {root}, dsts);
    if (path) edges.insert(path->begin(), path->end());
  }
  return {edges.begin(), edges.end()};
}

std::vector<Amount> member_capacities(const Graph& g, const PreparedRoute& plan) {
  std::vector<Amount> out;
  for (EdgeId e : plan.forwarding_edges) out.push_back(g.edge_capacity_units()[static_cast<std::size_t>(e)]);
  return out;
}

} // namespace

TEST(RosterRoutes, ChoicesEqualIndependentWalksOnMaskedLoadedGraphs) {
  std::mt19937 rng(20261102);
  int compared = 0, collapsed = 0, trees = 0, bottlenecks = 0, empty = 0, anycast = 0;
  for (int trial = 0; trial < 350; ++trial) {
    const auto n = 3 + static_cast<int>(rng() % 5);
    const Graph g = tied_graph(rng, n, 2 * n + static_cast<int>(rng() % 5), 1);
    std::optional<Eligibility> mask;
    if (trial % 3 == 0) mask = g.eligibility(std::vector<NodeId> {static_cast<NodeId>(rng() % n)},
        std::vector<EdgeId> {static_cast<EdgeId>(rng() % g.num_edges())});
    const Eligibility* elig = mask ? &*mask : nullptr;
    std::vector<NodeId> roots, dsts;
    for (NodeId v = 0; v < n; ++v) {
      if (v == n - 1 || (rng() % 5) == 0) {
        if (elig == nullptr || elig->allows_node(v)) dsts.push_back(v);
      } else if (rng() % 3 != 0 && (elig == nullptr || elig->allows_node(v))) {
        roots.push_back(v);
      }
    }
    // Bellman-Ford and the tight-edge definition construct the supplied support independently.
    const auto allowed = ref::members(ref::Query {&g, elig});
    const auto dist = ref::distances(g, allowed, dsts, true);
    const auto tight = ref::full_dag(g, allowed, dist, true);
    const std::vector<EdgeId> support(tight.begin(), tight.end());
    const auto dag = PredecessorDAGBuilder::of_edges(g, support);
    std::vector<Amount> W(static_cast<std::size_t>(g.num_edges()));
    for (std::size_t e = 0; e < W.size(); ++e) W[e] = static_cast<Amount>(rng() % (g.edge_capacity_units()[e] + 1));
    SCOPED_TRACE("tails " + ::testing::PrintToString(std::vector(g.edge_sources().begin(), g.edge_sources().end()))
        + " heads " + ::testing::PrintToString(std::vector(g.edge_destinations().begin(), g.edge_destinations().end()))
        + " costs " + ::testing::PrintToString(std::vector(g.edge_costs().begin(), g.edge_costs().end()))
        + " caps " + ::testing::PrintToString(std::vector(g.edge_capacity_units().begin(), g.edge_capacity_units().end()))
        + " view " + ::testing::PrintToString(W) + " roots " + ::testing::PrintToString(roots)
        + " dsts " + ::testing::PrintToString(dsts) + " support " + ::testing::PrintToString(support));
    for (const auto& choice : choices()) {
      auto selected = roots;
      if (choice.single_path() && choice.prefers_score() && selected.size() > 1) selected.resize(1);
      const PriorityValues priorities {&g, choice.tie_break, W};
      const auto got = shortest_path_routes(g, dag, selected, dsts, choice, PriorityReader::of(priorities));
      const auto want = shortest_union_reference(g, support, selected, dsts, choice, W);
      ASSERT_EQ(got, want) << "trial " << trial << " tie " << static_cast<int>(choice.tie_break)
                           << " mp " << choice.multipath << " parallel " << choice.parallel_edges;
      const auto result_dag = PredecessorDAGBuilder::of_edges(g, got);
      ASSERT_EQ(min_route_costs(g, result_dag, selected, dsts), walk_costs(g, flags_of(g, want), selected, dsts));
      ++compared;
      empty += got.empty();
      anycast += !got.empty() && dsts.size() > 1;
      if (!choice.single_path() && !choice.parallel_edges) {
        collapsed += got != walk_union(g, flags_of(g, support), selected, dsts);
      }
      if (choice.single_path() && !choice.prefers_score()) {
        trees += got != walk_union(g, flags_of(g, support), selected, dsts);
      }
      bottlenecks += !got.empty() && choice.single_path() && choice.prefers_score();
    }
  }
  EXPECT_GT(compared, 4000);
  EXPECT_GT(collapsed, 40);
  EXPECT_GT(trees, 100);
  EXPECT_GT(bottlenecks, 100);
  EXPECT_GT(empty, 500);
  EXPECT_GT(anycast, 100);
}

TEST(RosterRoutes, MultiRootPruneAndCostsEqualArbitraryDagWalks) {
  std::mt19937 rng(11022026);
  int compared = 0, removed = 0, changed_by_mask = 0, distinct_costs = 0;
  for (int trial = 0; trial < 500; ++trial) {
    const auto n = 4 + static_cast<int>(rng() % 4);
    std::vector<NodeId> rank(static_cast<std::size_t>(n));
    std::iota(rank.begin(), rank.end(), NodeId {0});
    std::shuffle(rank.begin(), rank.end(), rng);
    std::vector<InEdge> arcs;
    for (int i = 0; i < 3 * n; ++i) {
      NodeId u = static_cast<NodeId>(rng() % n);
      NodeId v = static_cast<NodeId>(rng() % n);
      if (u == v) continue;
      if (rank[static_cast<std::size_t>(u)] > rank[static_cast<std::size_t>(v)]) std::swap(u, v);
      const auto cap = static_cast<Amount>(rng() % 8);
      const auto cost = Cost {1} + static_cast<Cost>(rng() % 7);
      arcs.push_back({u, v, cap, cost});
    }
    const auto g = make_graph(n, arcs);
    std::vector<EdgeId> input;
    for (EdgeId e : all_ids(g)) if (rng() % 4 != 0) input.push_back(e);
    const auto dag = PredecessorDAGBuilder::of_edges(g, input);
    std::vector<NodeId> roots, dsts;
    for (NodeId v = 0; v < n; ++v) {
      if (v == n - 1 || rng() % 4 == 0) dsts.push_back(v);
      else if (rng() % 4 != 0) roots.push_back(v);
    }
    auto xn = std::vector<NodeId> {};
    auto xe = std::vector<EdgeId> {};
    if (trial % 3 == 0) xn.push_back(static_cast<NodeId>(rng() % n));
    for (EdgeId e : all_ids(g)) if (rng() % 8 == 0) xe.push_back(e);
    const auto mask = g.eligibility(xn, xe);
    for (const Eligibility* elig : {static_cast<const Eligibility*>(nullptr), &mask}) {
      auto allowed = flags_of(g, input);
      for (EdgeId e : input) if (elig != nullptr && !elig->allows_edge(e)) allowed[static_cast<std::size_t>(e)] = 0;
      const auto want = walk_union(g, allowed, roots, dsts);
      const auto got = prune(g, dag, roots, dsts, elig);
      ASSERT_EQ(got, want) << trial;
      const auto pruned = PredecessorDAGBuilder::of_edges(g, got);
      const auto costs = min_route_costs(g, pruned, roots, dsts);
      ASSERT_EQ(costs, walk_costs(g, allowed, roots, dsts)) << trial;
      // The raw DAG's dead ends and destination out-edges must not change a source's cost.
      if (elig == nullptr) ASSERT_EQ(min_route_costs(g, dag, roots, dsts), costs);
      ++compared;
      removed += got != input;
      changed_by_mask += elig != nullptr && got != walk_union(g, flags_of(g, input), roots, dsts);
      if (costs.size() > 1 && costs[0] != costs[1] && costs[0] != kUnreachable && costs[1] != kUnreachable) ++distinct_costs;
    }
  }
  EXPECT_EQ(compared, 1000);
  EXPECT_GT(removed, 700);
  EXPECT_GT(changed_by_mask, 40);
  EXPECT_GT(distinct_costs, 20);
}

TEST(RosterRoutes, CollapseUsesCreationScoresAndLowestEdgeIdOnTies) {
  // M->T capacities100/10, free10/5: capacity and available choose the first, headroom
  // chooses the second (1/2 > 1/10). Equal scores must choose the first EdgeId.
  const auto g = make_graph(4, {{0, 2, 100, 1}, {1, 2, 100, 1}, {2, 3, 100, 1}, {2, 3, 10, 1}});
  const auto support = PredecessorDAGBuilder::of_edges(g, all_ids(g));
  const auto W = to_ids(g, std::vector<Amount> {100, 100, 10, 5});
  const std::vector<NodeId> roots {0, 1}, dsts {3};
  for (auto tie : {RouteTieBreak::HigherCapacity, RouteTieBreak::HigherAvailableCapacity, RouteTieBreak::HigherHeadroom}) {
    const PriorityValues values {&g, tie, W};
    const auto got = shortest_path_routes(g, support, roots, dsts, {tie, true, false}, PriorityReader::of(values));
    EXPECT_EQ(ids_to_input(g, got), (std::vector<std::int32_t> {0, 1, tie == RouteTieBreak::HigherHeadroom ? 3 : 2}));
  }
  const auto equal = make_graph(4, {{0, 2, 100, 1}, {1, 2, 100, 1}, {2, 3, 10, 1}, {2, 3, 10, 1}});
  const PriorityValues same {&equal, RouteTieBreak::HigherAvailableCapacity, {}};
  EXPECT_EQ(ids_to_input(equal, shortest_path_routes(equal, *all_edges(equal), roots, dsts,
      {same.tie_break, true, false}, PriorityReader::of(same))), (std::vector<std::int32_t> {0, 1, 2}));
  // An excluded wider link cannot remove the allowed narrower one from the collapse oracle.
  const auto inv = id_of_input(g);
  std::vector<EdgeId> allowed {inv[0], inv[1], inv[3]};
  std::sort(allowed.begin(), allowed.end());
  const RouteChoice choice {RouteTieBreak::HigherCapacity, true, false};
  const PriorityValues values {&g, choice.tie_break, W};
  EXPECT_EQ(shortest_union_reference(g, allowed, roots, dsts, choice, W), allowed);
  EXPECT_EQ(shortest_path_routes(g, PredecessorDAGBuilder::of_edges(g, allowed), roots, dsts,
      choice, PriorityReader::of(values)), allowed);
  // Two sources merge at M. Keeping only its capacity10 link admits10; keeping both
  // parallel links would give the capacity1 link half the traffic and admit only2.
  const auto collapse = make_graph(4, {{0, 2, 100, 1}, {1, 2, 100, 1}, {2, 3, 1, 1}, {2, 3, 10, 1}});
  const PriorityValues cv {&collapse, choice.tie_break, {}};
  const auto U = shortest_path_routes(collapse, *all_edges(collapse), roots, dsts, choice, PriorityReader::of(cv));
  const auto cp = prepare_route(collapse, RouteAllocation::Equal, Admission::Lossless, {}, roots, dsts,
      PredecessorDAGBuilder::of_edges(collapse, U));
  EXPECT_EQ(max_lossless_offer(cp, member_capacities(collapse, cp)).value, 10);
}

TEST(RosterRoutes, CanonicalTreeDropsUnusedBranchesAndChoosesFewestHops) {
  // Both two-hop routes from0 cost2; lowest EdgeId takes 0->3->2, limited to3.
  // Root1 takes1->4->2: RR(T,2) gives ceil(T/2)<=3, hence T=6, not all-route10.
  const auto g = make_graph(5, {{0, 3, 100, 1}, {0, 4, 100, 1}, {3, 2, 3, 1}, {4, 2, 100, 1}, {1, 4, 100, 1}});
  const std::vector<NodeId> roots {0, 1}, dsts {2};
  const auto support = PredecessorDAGBuilder::of_edges(g, all_ids(g));
  const auto got = shortest_path_routes(g, support, roots, dsts, {RouteTieBreak::FewerHops, false, true}, {});
  EXPECT_EQ(ids_to_input(g, got), (std::vector<std::int32_t> {0, 2, 3, 4}));
  const auto route = PredecessorDAGBuilder::of_edges(g, got);
  const auto plan = prepare_route(g, RouteAllocation::Equal, Admission::Lossless, {}, roots, dsts, route);
  EXPECT_EQ(max_lossless_offer(plan, member_capacities(g, plan)).value, 6);
  const auto one = shortest_path_routes(g, support, std::vector<NodeId> {0}, dsts,
      {RouteTieBreak::FewerHops, false, true}, {});
  EXPECT_EQ(ids_to_input(g, one), (std::vector<std::int32_t> {0, 2}));
  // Direct cost2 beats two hops of cost1 despite the two-hop path's lower initial EdgeId.
  const auto hops = make_graph(3, {{0, 1, 10, 1}, {1, 2, 10, 1}, {0, 2, 10, 2}});
  EXPECT_EQ(ids_to_input(hops, shortest_path_routes(hops, *all_edges(hops), std::vector<NodeId> {0},
      std::vector<NodeId> {2}, {RouteTieBreak::FewerHops, false, true}, {})), (std::vector<std::int32_t> {2}));
}

TEST(RosterRoutes, SingleRootBottleneckKeepsItsPrefixSensitiveContinuation) {
  // The prefix from0 caps both paths at2, so their hop tie takes the lower-EdgeId branch
  // through3. From1 the prefix100 exposes bottlenecks5/3, so the branch through4 wins.
  const auto g = make_graph(6, {{0, 2, 2, 1}, {1, 2, 100, 1}, {2, 4, 5, 1},
      {2, 3, 3, 1}, {4, 5, 5, 1}, {3, 5, 3, 1}});
  const auto dag = all_edges(g);
  const std::vector<NodeId> dsts {5};
  const PriorityValues values {&g, RouteTieBreak::HigherCapacity, {}};
  const RouteChoice choice {values.tie_break, false, true};
  EXPECT_EQ(ids_to_input(g, shortest_path_routes(g, *dag, std::vector<NodeId> {0}, dsts,
      choice, PriorityReader::of(values))), (std::vector<std::int32_t> {0, 3, 5}));
  EXPECT_EQ(ids_to_input(g, shortest_path_routes(g, *dag, std::vector<NodeId> {1}, dsts,
      choice, PriorityReader::of(values))), (std::vector<std::int32_t> {1, 2, 4}));
  EXPECT_THROW((void)shortest_path_routes(g, *dag, std::vector<NodeId> {0, 1}, dsts,
      choice, PriorityReader::of(values)), std::logic_error);
}

TEST(RosterRoutes, NestedRootsKeepDistinctFullWidthCostsAndAddTheirInjection) {
  // The first source traverses both links; the second begins midway. Their costs are 5/3,
  // not the minimum distance3 for both; RR(7,2) puts4 on the first link and7 on the last.
  const auto g = make_graph(3, {{0, 1, 10, 2}, {1, 2, 10, 3}});
  const std::vector<NodeId> roots {0, 1}, dsts {2};
  const auto route = all_edges(g);
  EXPECT_EQ(min_route_costs(g, *route, roots, dsts), (std::vector<Cost> {5, 3}));
  const auto edges = shortest_path_routes(g, *route, roots, dsts, {}, {});
  EXPECT_EQ(edges, all_ids(g));
  const auto plan = prepare_route(g, RouteAllocation::Equal, Admission::Lossless, {}, roots, dsts, *route);
  std::vector<Amount> amounts(plan.forwarding_edges.size());
  ASSERT_TRUE(allocate(plan, 7, amounts));
  EXPECT_EQ(amounts, (std::vector<Amount> {4, 7}));
  // An all-edge cost sum2^62-1 is legal, including path-cost carries beyond32 and48 bits.
  const auto wide = make_graph(4, {{0, 2, 1, (Cost {1} << 61) - 1}, {1, 2, 1, 1}, {2, 3, 1, (Cost {1} << 61) - 1}});
  EXPECT_EQ(min_route_costs(wide, *all_edges(wide), roots, std::vector<NodeId> {3}),
      (std::vector<Cost> {(Cost {1} << 62) - 2, Cost {1} << 61}));
}

TEST(RosterRoutes, SuppliedDagCostsIgnoreDeadEndsAndStopAtFirstDestination) {
  // Costs0->5: min(1+1,1+1+5)=2; 1->5 costs4. Edge2->4 is a dead end,
  // and5->6 leaves the first destination. Neither may enter the pruned route.
  const auto g = make_graph(7, {{0, 2, 10, 1}, {1, 2, 10, 3}, {2, 5, 10, 1},
      {2, 3, 10, 1}, {3, 5, 10, 5}, {2, 4, 10, 1}, {5, 6, 10, 1}});
  const auto dag = all_edges(g);
  const std::vector<NodeId> roots {0, 1}, dsts {5, 6};
  const auto got = prune(g, *dag, roots, dsts);
  EXPECT_EQ(ids_to_input(g, got), (std::vector<std::int32_t> {0, 1, 2, 3, 4}));
  EXPECT_EQ(min_route_costs(g, *dag, roots, dsts), (std::vector<Cost> {2, 4}));
  const auto inv = id_of_input(g);
  const auto mask = g.eligibility({}, std::vector<EdgeId> {inv[2]});
  const auto repaired = prune(g, *dag, roots, dsts, &mask);
  EXPECT_EQ(ids_to_input(g, repaired), (std::vector<std::int32_t> {0, 1, 3, 4}));
  EXPECT_EQ(min_route_costs(g, PredecessorDAGBuilder::of_edges(g, repaired), roots, dsts), (std::vector<Cost> {7, 9}));
  const auto blocked = g.eligibility(std::vector<NodeId> {2}, {});
  EXPECT_TRUE(prune(g, *dag, roots, dsts, &blocked).empty());
}

TEST(RosterRoutes, WeightedMembersDoNotChangeUnprunedRoutesAndCosts) {
  // All four routes cost2. A zero weight on3->1 removes the whole0->3 branch from
  // members, leaving capacity10, while U and the source's record route keep all four edges.
  const auto g = make_graph(4, {{0, 2, 10, 1}, {0, 3, 10, 1}, {2, 1, 10, 1}, {3, 1, 10, 1}});
  const auto U = shortest_path_routes(g, *all_edges(g), std::vector<NodeId> {0}, std::vector<NodeId> {1}, {}, {});
  EXPECT_EQ(U, all_ids(g));
  const auto weights = to_ids(g, std::vector<Amount> {1, 1, 1, 0});
  const auto dag = PredecessorDAGBuilder::of_edges(g, U);
  const auto plan = prepare_route(g, RouteAllocation::Weighted, Admission::Lossless,
      {WeightBasis::Explicit, weights}, std::vector<NodeId> {0}, std::vector<NodeId> {1}, dag);
  EXPECT_EQ(ids_to_input(g, plan.forwarding_edges), (std::vector<std::int32_t> {0, 2}));
  EXPECT_EQ(plan.support_edges, U);
  EXPECT_EQ(min_route_costs(g, dag, std::vector<NodeId> {0}, std::vector<NodeId> {1}), (std::vector<Cost> {2}));
  EXPECT_EQ(max_lossless_offer(plan, member_capacities(g, plan)).value, 10);

  // The second root's zero weight leaves its roster position -1; RR(T,2) first gives it a
  // unit at T=2, so Lossless admits only1. The route and cost remain present for that root.
  const auto block = make_graph(3, {{0, 2, 10, 1}, {1, 2, 10, 1}});
  const std::vector<Amount> bw {1, 0};
  const auto bp = prepare_route(block, RouteAllocation::Weighted, Admission::Lossless,
      {WeightBasis::Explicit, bw}, std::vector<NodeId> {0, 1}, std::vector<NodeId> {2}, *all_edges(block));
  EXPECT_EQ(std::vector<std::int32_t>(bp.source_positions.begin(), bp.source_positions.end()), (std::vector<std::int32_t> {0, -1}));
  EXPECT_EQ(max_lossless_offer(bp, member_capacities(block, bp)).value, 1);
  EXPECT_EQ(min_route_costs(block, *all_edges(block), std::vector<NodeId> {0, 1}, std::vector<NodeId> {2}), (std::vector<Cost> {1, 1}));

  // Snapshot weights3/9 are aligned with ascending U, not with input position or capacity10/10.
  // H(12,[3,9])=(3,9) fits exactly; equal capacity weights would admit only6.
  const auto snapshot = make_graph(2, {{0, 1, 10, 1}, {0, 1, 10, 1}});
  const std::vector<Amount> free {3, 9};
  const auto frozen = prepare_route(snapshot, RouteAllocation::Weighted, Admission::Lossless,
      {WeightBasis::ResidualSnapshot, free}, std::vector<NodeId> {0}, std::vector<NodeId> {1}, *all_edges(snapshot));
  EXPECT_EQ(frozen.weights, free);
  EXPECT_EQ(max_lossless_offer(frozen, free).value, 12);
}

TEST(RosterRoutes, EmptyUnreachableAndTerminalInputs) {
  const auto g = make_graph(4, {{0, 1, 0, 1}, {1, 2, 0, 1}});
  const auto dag = all_edges(g);
  const std::vector<NodeId> roots {0, 3}, dsts {1, 2};
  EXPECT_EQ(ids_to_input(g, prune(g, *dag, roots, dsts)), (std::vector<std::int32_t> {0}));
  EXPECT_EQ(min_route_costs(g, *dag, roots, dsts), (std::vector<Cost> {1, kUnreachable}));
  EXPECT_EQ(min_route_costs(g, *dag, roots, std::vector<NodeId> {}), (std::vector<Cost> {kUnreachable, kUnreachable}));
  EXPECT_TRUE(min_route_costs(g, *dag, std::vector<NodeId> {}, dsts).empty());
  EXPECT_TRUE(prune(g, *dag, std::vector<NodeId> {}, dsts).empty());
  EXPECT_TRUE(prune(g, *dag, roots, std::vector<NodeId> {}).empty());
  EXPECT_TRUE(shortest_path_routes(g, *dag, std::vector<NodeId> {}, dsts, {}, {}).empty());
  EXPECT_TRUE(shortest_path_routes(g, *dag, roots, std::vector<NodeId> {}, {}, {}).empty());
  // Empty graph and unused nodes cannot leave stale slots from the preceding nonempty calls.
  const auto empty = make_graph(4, {});
  EXPECT_EQ(min_route_costs(empty, *all_edges(empty), roots, dsts), (std::vector<Cost> {kUnreachable, kUnreachable}));
  EXPECT_TRUE(shortest_path_routes(empty, *all_edges(empty), roots, dsts, {}, {}).empty());
}

TEST(RosterRoutes, ConcurrentCallsKeepTheirScratchIndependent) {
  const auto first = make_graph(4, {{0, 2, 10, 2}, {1, 2, 10, 4}, {2, 3, 10, 3}});
  const auto second = make_graph(5, {{0, 1, 10, 1}, {1, 3, 10, 2}, {0, 2, 10, 1}, {2, 3, 10, 2}});
  const auto run = [](const Graph& g, const std::vector<NodeId>& roots, const std::vector<Cost>& costs) {
    const auto dag = all_edges(g);
    const std::vector<NodeId> dsts {3};
    for (int i = 0; i < 100; ++i) {
      EXPECT_EQ(shortest_path_routes(g, *dag, roots, dsts, {}, {}), all_ids(g));
      EXPECT_EQ(prune(g, *dag, roots, dsts), all_ids(g));
      EXPECT_EQ(min_route_costs(g, *dag, roots, dsts), costs);
    }
  };
  // The two roots cost5/7 in the first graph; the second graph's one source has two paths
  // costing3. A shared unstamped/global scratch would mix their different rows and costs.
  std::thread a([&] { run(first, {0, 1}, {5, 7}); });
  std::thread b([&] { run(second, {0}, {3}); });
  a.join();
  b.join();
}

namespace {

// Enumerated walks supply U, per-root costs and positive-weight walks independently. Check every
// plan array by NodeId/EdgeId: internal topological positions are not part of the routing rule.
void expect_preparation(const Graph& g, const ShortestPathPreparation& got,
    const std::vector<EdgeId>& route, const std::vector<NodeId>& roots, const std::vector<NodeId>& dsts,
    RouteAllocation allocation, const std::vector<Amount>& edge_weights) {
  const auto& p = got.plan;
  ASSERT_EQ(p.allocation, allocation);
  ASSERT_EQ(p.support_edges, route);
  ASSERT_EQ(got.costs, walk_costs(g, flags_of(g, route), roots, dsts));
  if (allocation == RouteAllocation::Optimized) {
    EXPECT_TRUE(p.forwarding_edges.empty() && p.weights.empty() && p.sources.empty() &&
        p.source_positions.empty() && p.destinations.empty() && p.topological_order.empty() &&
        p.outgoing_begin.empty() && p.outgoing_end.empty() && p.edge_destination_positions.empty() && p.is_destination.empty());
    return;
  }
  auto positive = flags_of(g, route);
  if (allocation == RouteAllocation::Weighted) {
    for (EdgeId e : route) if (edge_weights[static_cast<std::size_t>(e)] == 0) positive[static_cast<std::size_t>(e)] = 0;
  }
  const auto members = walk_union(g, positive, roots, dsts);
  ASSERT_EQ(p.forwarding_edges, members);
  std::vector<Amount> weights;
  std::set<NodeId> nodes;
  for (EdgeId e : members) {
    nodes.insert(g.edge_sources()[static_cast<std::size_t>(e)]);
    nodes.insert(g.edge_destinations()[static_cast<std::size_t>(e)]);
    if (allocation == RouteAllocation::Weighted) weights.push_back(edge_weights[static_cast<std::size_t>(e)]);
  }
  ASSERT_EQ(p.weights, weights);
  ASSERT_EQ(vec(p.sources), roots);
  ASSERT_EQ(vec(p.destinations), dsts);
  ASSERT_EQ(p.source_positions.size(), roots.size());
  ASSERT_EQ(p.topological_order.size(), nodes.size());
  ASSERT_EQ(p.outgoing_begin.size(), nodes.size());
  ASSERT_EQ(p.outgoing_end.size(), nodes.size());
  ASSERT_EQ(p.is_destination.size(), nodes.size());
  ASSERT_EQ(p.edge_destination_positions.size(), members.size());
  ASSERT_EQ((std::set<NodeId>(p.topological_order.begin(), p.topological_order.end())), nodes);
  std::map<NodeId, std::int32_t> position;
  for (std::size_t i = 0; i < p.topological_order.size(); ++i) {
    const NodeId u = p.topological_order[i];
    position[u] = static_cast<std::int32_t>(i);
    EXPECT_EQ(p.is_destination[i], std::binary_search(dsts.begin(), dsts.end(), u));
    std::vector<EdgeId> outgoing;
    for (EdgeId e : members) if (g.edge_sources()[static_cast<std::size_t>(e)] == u) outgoing.push_back(e);
    const auto lo = p.outgoing_begin[i], hi = p.outgoing_end[i];
    ASSERT_GE(lo, 0);
    ASSERT_GE(hi, lo);
    ASSERT_LE(static_cast<std::size_t>(hi), members.size());
    EXPECT_EQ((std::vector<EdgeId>(members.begin() + lo, members.begin() + hi)), outgoing);
    if (outgoing.empty()) EXPECT_EQ(lo, 0);
  }
  for (std::size_t i = 0; i < roots.size(); ++i) {
    const auto it = position.find(roots[i]);
    EXPECT_EQ(p.source_positions[i], it == position.end() ? -1 : it->second);
  }
  for (std::size_t k = 0; k < members.size(); ++k) {
    const auto e = static_cast<std::size_t>(members[k]);
    EXPECT_EQ(p.edge_destination_positions[k], position.at(g.edge_destinations()[e]));
    EXPECT_LT(position.at(g.edge_sources()[e]), p.edge_destination_positions[k]);
  }
}

} // namespace

TEST(ShortestPathsPreparation, IndependentWalksCheckEveryChoiceWeightAndRootShape) {
  std::mt19937 rng(20261103);
  int compared = 0, pruned = 0, nested = 0, memberless = 0, snapshot_calls = 0;
  for (int trial = 0; trial < 160; ++trial) {
    const int n = 3 + static_cast<int>(rng() % 5);
    const auto g = tied_graph(rng, n, 3 * n);
    const auto mask = g.eligibility(trial % 3 == 0 ? std::vector<NodeId> {static_cast<NodeId>(rng() % n)} : std::vector<NodeId> {},
        std::vector<EdgeId> {static_cast<EdgeId>(rng() % g.num_edges())});
    std::vector<NodeId> roots, dsts;
    for (NodeId v = 0; v < n; ++v) {
      if (!mask.allows_node(v)) continue;
      if (v == n - 1 || rng() % 4 == 0) dsts.push_back(v);
      else if (trial % 4 != 0 || roots.empty()) roots.push_back(v);
    }
    const auto allowed = ref::members(ref::Query {&g, &mask});
    const auto dist = ref::distances(g, allowed, dsts, true);
    const auto tight = ref::full_dag(g, allowed, dist, true);
    const std::vector<EdgeId> support(tight.begin(), tight.end());
    const auto dag = PredecessorDAGBuilder::of_edges(g, support);
    std::vector<Amount> W(static_cast<std::size_t>(g.num_edges())), explicit_weights(W.size());
    for (std::size_t e = 0; e < W.size(); ++e) {
      W[e] = static_cast<Amount>(rng() % (g.edge_capacity_units()[e] + 1));
      explicit_weights[e] = static_cast<Amount>(rng() % 7);
    }
    for (const auto& choice : choices()) {
      auto selected = roots;
      if (choice.single_path() && choice.prefers_score() && selected.size() > 1) selected.resize(1);
      const auto want = shortest_union_reference(g, support, selected, dsts, choice, W);
      const PriorityValues priorities {&g, choice.tie_break, W};
      for (int mode = 0; mode < 5; ++mode) {
        const auto allocation = mode == 0 ? RouteAllocation::Equal : mode == 1 ? RouteAllocation::Optimized : RouteAllocation::Weighted;
        const auto basis = mode == 3 ? WeightBasis::Explicit : mode == 4 ? WeightBasis::ResidualSnapshot : WeightBasis::Capacity;
        const WeightSource weights {basis, mode == 3 ? std::span<const Amount>(explicit_weights) : std::span<const Amount> {}};
        int gathered = 0;
        const auto got = prepare_shortest_path_routes(g, dag, selected, dsts, choice, PriorityReader::of(priorities), allocation, weights,
            [&](std::span<const EdgeId> edges, std::span<Amount> values) {
              ++gathered;
              EXPECT_EQ(vec(edges), want);
              for (std::size_t k = 0; k < edges.size(); ++k) values[k] = W[static_cast<std::size_t>(edges[k])];
            });
        std::vector<Amount> expected_weights(g.edge_capacity_units().begin(), g.edge_capacity_units().end());
        if (mode == 3) for (std::size_t e = 0; e < W.size(); ++e) expected_weights[e] = g.edge_capacity_units()[e] == 0 ? 0 : explicit_weights[e];
        if (mode == 4) expected_weights = W;
        SCOPED_TRACE("trial " + std::to_string(trial) + " choice " + std::to_string(static_cast<int>(choice.tie_break)) +
            " multipath " + std::to_string(choice.multipath) + " parallel " + std::to_string(choice.parallel_edges) + " mode " + std::to_string(mode));
        expect_preparation(g, got, want, selected, dsts, allocation, expected_weights);
        ASSERT_FALSE(::testing::Test::HasFailure());
        EXPECT_EQ(gathered, mode == 4 && !want.empty() ? 1 : 0);
        snapshot_calls += gathered;
        ++compared;
        pruned += allocation == RouteAllocation::Weighted && got.plan.forwarding_edges != want;
        if (allocation == RouteAllocation::Optimized) continue;
        for (auto p : got.plan.source_positions) memberless += p < 0;
        for (EdgeId e : want) nested += std::binary_search(selected.begin(), selected.end(), g.edge_destinations()[static_cast<std::size_t>(e)]);
      }
    }
  }
  EXPECT_EQ(compared, 10400);
  EXPECT_GT(pruned, 100);
  EXPECT_GT(memberless, 100);
  EXPECT_GT(nested, 100);
  EXPECT_GT(snapshot_calls, 100);
}

TEST(ShortestPathsPreparation, HandCostsPoolingWidthsAndOwnership) {
  // Prefixes1/3 + suffix2 give3/5 before weighting; zero-weight source1 stays in the roster.
  const auto g = make_graph(5, {{0, 2, 30, 1}, {1, 2, 30, 3}, {2, 3, 30, 2}});
  const auto dag = all_edges(g);
  const std::vector<NodeId> roots {0, 1, 2}, dsts {3, 4};
  std::vector<Amount> weights {1, 0, 1};
  auto got = prepare_shortest_path_routes(g, *dag, roots, dsts, {}, {}, RouteAllocation::Weighted,
      {WeightBasis::Explicit, weights}, {});
  EXPECT_EQ(got.costs, (std::vector<Cost> {3, 5, 2}));
  expect_preparation(g, got, all_ids(g), roots, dsts, RouteAllocation::Weighted, weights);
  EXPECT_EQ(got.plan.source_positions[1], -1);
  PreparedRoute copy = got.plan;
  const auto before = PreparedRoute(copy);
  got = {};
  for (int i = 0; i < 20; ++i) {
    const auto other = prepare_shortest_path_routes(g, *dag, std::vector<NodeId> {1}, dsts, {}, {}, RouteAllocation::Equal, {}, {});
    EXPECT_EQ(other.costs, (std::vector<Cost> {5}));
  }
  expect_same_plan(copy, before, "copy owns all arrays after donor and scratch reuse");
  auto equal = prepare_shortest_path_routes(g, *dag, std::vector<NodeId> {0, 2}, dsts, {}, {}, RouteAllocation::Equal, {}, {});
  std::vector<Amount> load(equal.plan.forwarding_edges.size());
  ASSERT_TRUE(allocate(equal.plan, 11, load));
  EXPECT_EQ(load, (std::vector<Amount> {6, 11}));  // root0 injects6, nested root2 adds5
  std::vector<NodeId> input_roots {0}, input_dsts {3, 4};
  auto owned = prepare_shortest_path_routes(g, *dag, input_roots, input_dsts, {}, {}, RouteAllocation::Equal, {}, {});
  input_roots[0] = 1;
  input_dsts[0] = 2;
  EXPECT_EQ(vec(owned.plan.sources), (std::vector<NodeId> {0}));
  EXPECT_EQ(vec(owned.plan.destinations), dsts);
  const auto wide = make_graph(3, {{0, 1, 1, Cost {1} << 61}, {1, 2, 1, (Cost {1} << 61) - 1}});
  const auto wide_got = prepare_shortest_path_routes(wide, *all_edges(wide), std::vector<NodeId> {0, 1},
      std::vector<NodeId> {2}, {}, {}, RouteAllocation::Equal, {}, {});
  EXPECT_EQ(wide_got.costs, (std::vector<Cost> {(Cost {1} << 62) - 1, (Cost {1} << 61) - 1}));
}

TEST(ShortestPathsPreparation, ManyRootsKeepIndividualCostsAndOrderedShares) {
  std::vector<InEdge> edges;
  std::vector<NodeId> roots;
  std::vector<Cost> costs;
  for (NodeId i = 0; i < 65; ++i) {
    edges.push_back({i, 65, 3, i + 1});
    roots.push_back(i);
    costs.push_back(i + 2);
  }
  edges.push_back({65, 66, 200, 1});
  const auto g = make_graph(67, edges);
  const auto got = prepare_shortest_path_routes(g, *all_edges(g), roots, std::vector<NodeId> {66}, {}, {}, RouteAllocation::Equal, {}, {});
  EXPECT_EQ(got.costs, costs);
  expect_preparation(g, got, all_ids(g), roots, {66}, RouteAllocation::Equal, {});
  // RR(67,65) gives2 to the first two heads,1 to every other head, then pools67 at65.
  std::vector<Amount> load(got.plan.forwarding_edges.size());
  ASSERT_TRUE(allocate(got.plan, 67, load));
  ASSERT_EQ(load.size(), 66U);
  EXPECT_EQ(load[0], 2);
  EXPECT_EQ(load[1], 2);
  EXPECT_EQ(load[64], 1);
  EXPECT_EQ(load[65], 67);
}

TEST(ShortestPathsPreparation, EmptyAndConcurrentCallsKeepReturnedArraysIndependent) {
  const auto g = make_graph(4, {{0, 2, 8, 1}, {1, 2, 8, 3}, {2, 3, 8, 2}});
  const auto dag = all_edges(g);
  for (const auto& roots : {std::vector<NodeId> {}, std::vector<NodeId> {0, 1}}) {
    for (const auto& dsts : {std::vector<NodeId> {}, std::vector<NodeId> {3}}) {
      if (!roots.empty() && !dsts.empty()) continue;
      const auto got = prepare_shortest_path_routes(g, *dag, roots, dsts, {}, {}, RouteAllocation::Equal, {}, {});
      expect_preparation(g, got, {}, roots, dsts, RouteAllocation::Equal, {});
    }
  }
  auto run = [&](const std::vector<NodeId>& roots, const std::vector<Cost>& costs) {
    for (int i = 0; i < 100; ++i) {
      const auto got = prepare_shortest_path_routes(g, *dag, roots, std::vector<NodeId> {3}, {}, {}, RouteAllocation::Equal, {}, {});
      EXPECT_EQ(got.costs, costs);
      expect_preparation(g, got, roots.size() == 1 ? std::vector<EdgeId> {1, 2} : all_ids(g), roots,
          std::vector<NodeId> {3}, RouteAllocation::Equal, {});
    }
  };
  std::thread a([&] { run({0, 1}, {3, 5}); }), b([&] { run({1}, {5}); });
  a.join();
  b.join();
}
