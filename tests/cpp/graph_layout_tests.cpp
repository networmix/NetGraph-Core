// The graph's stored layout and construction: the independent layout checker on every graph, its
// corruption controls, the differential test against the pinned earlier constructor
// (graph_reference.hpp; every array, every error), the hash/fallback boundary, copies and
// moves, and the retained array payload.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <numeric>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <typeinfo>
#include <vector>

#include "graph_inputs.hpp"
#include "graph_reference.hpp"
#include "layout_checker.hpp"
#include "netgraph/core/shortest_path_cache.hpp"
#include "netgraph/core/detail/graph_access.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/search.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using detail::GraphAccess;

namespace {

// Every public array, its element width, the counts, the resolution and the scans equal the
// reference's; the internal copy equals the reference's edge_sources()[incoming_edges()[slot]].
void expect_same_as_reference(const Graph& g, const Input& in) {
  const auto r = reference_graph(in.n, in.src, in.dst, in.cap, in.cost, in.q);
  static_assert(sizeof(decltype(g.edge_capacity_units())::element_type) == 8);
  static_assert(sizeof(decltype(g.edge_costs())::element_type) == 8);
  static_assert(sizeof(decltype(g.edge_sources())::element_type) == 4);
  static_assert(sizeof(decltype(g.input_positions())::element_type) == 4);
  static_assert(sizeof(decltype(g.incoming_edges())::element_type) == 4);
  EXPECT_EQ(g.num_nodes(), r.num_nodes) << in.name;
  EXPECT_EQ(g.num_edges(), static_cast<std::int32_t>(r.capacity.size())) << in.name;
  EXPECT_EQ(g.resolution(), r.resolution) << in.name;
  EXPECT_EQ(g.max_edge_cost(), r.max_edge_cost) << in.name;
  EXPECT_EQ(g.edge_cost_gcd(), r.edge_cost_gcd) << in.name;
  EXPECT_EQ(vec(g.edge_capacity_units()), r.capacity) << in.name;
  EXPECT_EQ(vec(g.edge_costs()), r.cost) << in.name;
  EXPECT_EQ(vec(g.edge_sources()), r.src) << in.name;
  EXPECT_EQ(vec(g.edge_destinations()), r.dst) << in.name;
  EXPECT_EQ(vec(g.input_positions()), r.input_positions) << in.name;
  EXPECT_EQ(vec(g.outgoing_offsets()), r.outgoing_offsets) << in.name;
  EXPECT_EQ(vec(g.incoming_offsets()), r.incoming_offsets) << in.name;
  EXPECT_EQ(vec(g.incoming_edges()), r.incoming_edges) << in.name;
  EXPECT_EQ(vec(GraphAccess::incoming_edge_sources(g)), reference_in_col_indices(r)) << in.name;
  const auto bad = check_layout(spans_of(g), view(in));
  EXPECT_TRUE(bad.empty()) << in.name << ": " << (bad.empty() ? "" : bad.front());
}

// --------------------------------------------------------------------------- generators

Input permuted(const Input& in, const std::vector<std::size_t>& p, const std::string& suffix) {
  Input out {in.name + suffix, in.n, {}, {}, {}, {}, in.q};
  for (std::size_t i : p) {
    out.src.push_back(in.src[i]);
    out.dst.push_back(in.dst[i]);
    out.cap.push_back(in.cap[i]);
    out.cost.push_back(in.cost[i]);
  }
  return out;
}

// The same edges in EdgeId order, reversed and shuffled (all property arrays permuted together).
std::vector<Input> orders(const Input& in, std::uint64_t seed) {
  const auto r = reference_graph(in.n, in.src, in.dst, in.cap, in.cost, in.q);
  std::vector<std::size_t> sorted(r.input_positions.begin(), r.input_positions.end());
  std::vector<std::size_t> rev(in.src.size());
  std::iota(rev.rbegin(), rev.rend(), 0);
  std::vector<std::size_t> sh(in.src.size());
  std::iota(sh.begin(), sh.end(), 0);
  std::shuffle(sh.begin(), sh.end(), std::mt19937_64(seed));
  return {in, permuted(in, sorted, " (EdgeId order)"), permuted(in, rev, " (reversed)"), permuted(in, sh, " (shuffled)")};
}

} // namespace

// --------------------------------------------------------------------------- differential

TEST(GraphLayout, DifferentialSmallAndDegenerateShapes) {
  std::vector<Input> cases = {
      {"N=0", 0},
      {"N=1, E=0", 1},
      {"N=1 self-loops", 1, {0, 0, 0}, {0, 0, 0}, {3, 1, 2}, {2, 2, 1}},
      {"E=0, N=1M", 1 << 20},
      random_input("N >> E", 1 << 18, 64, 9, 3),
      {"parallel edges, distinct capacities", 2, {0, 0, 0, 0}, {1, 1, 1, 1}, {4, 3, 2, 1}, {5, 5, 5, 5}},
  };
  for (const auto& in : cases) expect_same_as_reference(build(in), in);
}

TEST(GraphLayout, DifferentialRandomGraphsInFourOrders) {
  for (std::uint64_t seed = 1; seed <= 200; ++seed) {
    const auto base = random_input("seed " + std::to_string(seed), 1 + static_cast<std::int32_t>(seed % 37),
                                   1 + seed % 211, 1 + static_cast<Cost>(seed % 13), seed);
    for (const auto& in : orders(base, seed)) expect_same_as_reference(build(in), in);
  }
}

TEST(GraphLayout, DifferentialOneLongRowEachDirection) {
  Input out {"long out-row", 50};
  for (int i = 0; i < 20000; ++i) {
    out.src.push_back(7);
    out.dst.push_back(i % 50);
    out.cap.push_back(i);
    out.cost.push_back(1 + (i * 7919) % 31);
  }
  for (const auto& in : orders(out, 11)) expect_same_as_reference(build(in), in);
  Input in = out;
  in.name = "long in-row";
  std::swap(in.src, in.dst);
  for (const auto& x : orders(in, 12)) expect_same_as_reference(build(x), x);
}

TEST(GraphLayout, DifferentialDistinctCostCounts) {
  for (std::int64_t k : {1, 2, 31, 4096, 65535, 65536, 65537}) {
    const auto in = exact_k(k, 4096, 4096, static_cast<std::uint64_t>(k));
    ASSERT_EQ(distinct(in.cost), k);
    expect_same_as_reference(build(in), in);
  }
  // K about E: every cost distinct.
  const auto all = exact_k(200000, 0, 4096, 7);
  expect_same_as_reference(build(all), all);
}

// Exact units above 2^53 (capacities and costs kept as int64, not through double).
TEST(GraphLayout, DifferentialWideIntegers) {
  // Costs and units above 2^53, sums just below 2^62; separate from the large-K cases (the cost
  // sum bound forbids combining them).
  const Cost big = Cost{1} << 59;  // four of them sum to 2^61: below the bound
  Input wide {"wide costs", 3, {0, 0, 1, 2, 0}, {1, 1, 2, 0, 1}, {(Amount{1} << 53) + 1, Amount{1} << 61, 0, 5, 3},
              {big + 3, big + 1, big + 2, 1, big + 1}};
  expect_same_as_reference(build(wide), wide);
  Input sums {"sums just below 2^62", 2, {0, 0}, {1, 1}, {kL62 / 2, kL62 / 2 - 1}, {kL62 / 2, kL62 / 2 - 1}};
  expect_same_as_reference(build(sums), sums);
  Input units {"units above 2^53", 2, {0, 1}, {1, 0}, {(Amount{1} << 53) + 1, (Amount{1} << 53) + 3}, {1, 1}};
  expect_same_as_reference(build(units), units);
  // Just above 2^53 (no double holds 2^53 + 1), both arrays; EdgeId 0 is the cheaper 0 -> 1.
  const Amount above = (Amount{1} << 53) + 1;
  const Input both {"2^53 + 1", 2, {0, 0, 1}, {1, 1, 0}, {above, above + 2, 0}, {(Cost{1} << 53) + 1, (Cost{1} << 53) + 3, 1}};
  const Graph g = build(both);
  expect_same_as_reference(g, both);
  EXPECT_EQ(g.edge_capacity_units()[0], above);
  EXPECT_EQ(g.edge_costs()[0], (Cost{1} << 53) + 1);
}

TEST(GraphLayout, DifferentialHashCollisions) {
  // 40 costs colliding in the initial 64-slot table (their low 6 hash bits equal), and 40 sharing
  // their low 18 hash bits, so they collide in every table these inputs reach.
  for (std::uint64_t mask : {std::uint64_t {63}, (std::uint64_t {1} << 18) - 1}) {
    Input in {"collisions mask " + std::to_string(mask), 16};
    in.cost = colliding_costs(mask, 40);
    for (std::size_t i = 0; i < in.cost.size(); ++i) {
      in.src.push_back(static_cast<NodeId>(i % 16));
      in.dst.push_back(static_cast<NodeId>((i * 5) % 16));
      in.cap.push_back(static_cast<Amount>(i));
    }
    expect_same_as_reference(build(in), in);
  }
}

// --------------------------------------------------------------------------- errors

namespace {

// The exception a construction raises, as (type name, message); ("", "") when none.
template <typename F>
std::pair<std::string, std::string> outcome(F&& fn) {
  try {
    fn();
  } catch (const std::exception& e) {
    return {typeid(e).name(), e.what()};
  }
  return {"", ""};
}

// The defects of the constructor's validation steps, each applied to a copy of a valid input.
using Defect = std::function<void(Input&)>;

std::vector<std::pair<std::string, Defect>> defects() {
  return {
      {"resolution", [](Input& x) { x.q = 3.0; }},
      {"capacity negative at 2", [](Input& x) { x.cap[2] = -1; }},
      {"capacity 2^62 at 1", [](Input& x) { x.cap[1] = kL62; }},
      {"num_nodes negative", [](Input& x) { x.n = -1; }},
      {"dst shorter", [](Input& x) { x.dst.pop_back(); }},
      {"cost longer", [](Input& x) { x.cost.push_back(1); }},
      {"cost zero at 3", [](Input& x) { x.cost[3] = 0; }},
      {"cost negative at 1", [](Input& x) { x.cost[1] = -4; }},
      {"src out of range at 2", [](Input& x) { x.src[2] = x.n; }},
      {"dst negative at 1", [](Input& x) { x.dst[1] = -1; }},
      {"capacity sum", [](Input& x) {
         x.cap[0] = kL62 / 2;
         x.cap[3] = kL62 / 2;
       }},
      {"cost sum", [](Input& x) {
         x.cost[1] = kL62 / 2;
         x.cost[2] = kL62 / 2;
       }},
  };
}

} // namespace

TEST(GraphLayout, DifferentialErrorsSingleAndPairwise) {
  const Input base {"base", 4, {0, 1, 2, 3, 0}, {1, 2, 3, 0, 2}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}};
  const auto ds = defects();
  for (std::size_t a = 0; a < ds.size(); ++a) {
    for (std::size_t b = a; b < ds.size(); ++b) {
      Input in = base;
      ds[a].second(in);
      if (b != a) ds[b].second(in);
      const auto got = outcome([&] { (void)build(in); });
      const auto want = outcome([&] { (void)reference_graph(in.n, in.src, in.dst, in.cap, in.cost, in.q); });
      ASSERT_FALSE(want.first.empty()) << ds[a].first << " + " << ds[b].first;
      EXPECT_EQ(got, want) << ds[a].first << " + " << ds[b].first;
    }
  }
}

TEST(GraphLayout, DifferentialErrorsFirstInvalidIndexAndSumOrder) {
  // Two defects of one kind: the lower input index is named.
  Input in {"two costs", 3, {0, 1, 2, 0}, {1, 2, 0, 2}, {1, 1, 1, 1}, {1, -2, 0, 1}};
  auto got = outcome([&] { (void)build(in); });
  EXPECT_EQ(got, outcome([&] { (void)reference_graph(in.n, in.src, in.dst, in.cap, in.cost, in.q); }));
  EXPECT_NE(got.second.find("edge_costs[1] = -2"), std::string::npos) << got.second;
  // Capacity and cost sums both over: the one whose running sum crosses first, in input order.
  Input both {"sums", 2, {0, 0, 0}, {1, 1, 1}, {1, kL62 / 2, kL62 / 2}, {kL62 / 2, kL62 / 2, 1}};
  got = outcome([&] { (void)build(both); });
  EXPECT_EQ(got, outcome([&] { (void)reference_graph(both.n, both.src, both.dst, both.cap, both.cost, both.q); }));
  EXPECT_NE(got.second.find("edge costs"), std::string::npos) << got.second;
  // Both running sums reach 2^62 at the same edge: the capacity sum is reported.
  Input tie {"sums at one edge", 2, {0, 0}, {1, 1}, {kL62 / 2, kL62 / 2}, {kL62 / 2, kL62 / 2}};
  got = outcome([&] { (void)build(tie); });
  EXPECT_EQ(got, outcome([&] { (void)reference_graph(tie.n, tie.src, tie.dst, tie.cap, tie.cost, tie.q); }));
  EXPECT_NE(got.second.find("sum of all capacities"), std::string::npos) << got.second;
}

// --------------------------------------------------------------------------- the layout checker

// The corruption controls of the checker, on a graph built by the constructor.
// The fixture repeats (source, cost, destination) keys with distinct capacities and input positions, and its
// reverse slot order differs from its id order, with different source values.
TEST(GraphLayout, CheckerRejectsEveryCorruption) {
  const Input in {"controls", 4, {0, 2, 0, 1, 0, 3, 2, 0}, {1, 1, 1, 1, 1, 1, 3, 3}, {5, 6, 7, 8, 9, 1, 2, 3},
                  {2, 1, 2, 1, 2, 1, 4, 4}};
  const Graph g = build(in);
  ASSERT_TRUE(check_layout(spans_of(g), view(in)).empty());
  using Copy = GraphArrays;
  const Copy clean = arrays_of(g);
  bool reverse_differs = false;
  for (std::size_t s = 0; s < clean.in_adj.size(); ++s) reverse_differs |= clean.in_adj[s] != static_cast<EdgeId>(s);
  ASSERT_TRUE(reverse_differs);
  auto check = [&](const Copy& c) {
    return check_layout({in.n, static_cast<std::int64_t>(in.src.size()), c.capacity, c.cost, c.src, c.dst, c.input_positions,
                         c.outgoing_offsets, c.incoming_offsets, c.in_adj, std::span<const NodeId>(c.in_col), c.max_edge_cost,
                         c.edge_cost_gcd},
                        view(in));
  };
  const std::vector<std::tuple<std::string, std::string, std::function<void(Copy&)>>> controls = {
      // Wrong source or row boundary
      {"wrong source", "forward:membership", [](Copy& c) { c.src[0] = 1; }},
      {"row boundary moved", "forward:membership", [](Copy& c) { ++c.outgoing_offsets[1]; }},
      // All-zero offsets (a graph with edges)
      {"all-zero offsets", "row_offsets:end", [](Copy& c) { std::fill(c.outgoing_offsets.begin(), c.outgoing_offsets.end(), 0); }},
      // Missing or duplicated reverse coverage
      {"duplicated reverse slot", "slot:permutation", [](Copy& c) { c.in_adj[1] = c.in_adj[0]; }},
      {"wrong reverse destination", "reverse:membership", [](Copy& c) { ++c.incoming_offsets[1]; }},
      // (source, cost) substituted for (cost, source)
      {"(src, cost) reverse order", "reverse:order",
       [](Copy& c) {
         for (std::size_t v = 0; v + 1 < c.incoming_offsets.size(); ++v) {
           std::sort(c.in_adj.begin() + c.incoming_offsets[v], c.in_adj.begin() + c.incoming_offsets[v + 1],
                     [&](EdgeId p, EdgeId q) {
                       return std::tuple(c.src[static_cast<std::size_t>(p)], c.cost[static_cast<std::size_t>(p)], p) <
                              std::tuple(c.src[static_cast<std::size_t>(q)], c.cost[static_cast<std::size_t>(q)], q);
                     });
         }
         for (std::size_t s = 0; s < c.in_col.size(); ++s) c.in_col[s] = c.src[static_cast<std::size_t>(c.in_adj[s])];
       }},
      // Unstable equal-key members
      {"unstable equal keys", "forward:order",
       [](Copy& c) {
         for (std::size_t e = 0; e + 1 < c.cost.size(); ++e) {
           if (c.src[e] == c.src[e + 1] && c.cost[e] == c.cost[e + 1] && c.dst[e] == c.dst[e + 1]) {
             std::swap(c.input_positions[e], c.input_positions[e + 1]);
             std::swap(c.capacity[e], c.capacity[e + 1]);
             return;
           }
         }
         FAIL() << "the fixture has no equal-key pair";
       }},
      // Corrupted map
      {"map entries swapped", "property:capacity", [](Copy& c) { std::swap(c.input_positions[0], c.input_positions[1]); }},
      // Property-only mispermutation
      {"capacities swapped", "property:capacity", [](Copy& c) { std::swap(c.capacity[1], c.capacity[2]); }},
      // Out-of-range map or slot before access
      {"negative map entry", "map:permutation", [](Copy& c) { c.input_positions[1] = -3; }},
      {"duplicated map entry", "map:permutation", [](Copy& c) { c.input_positions[1] = c.input_positions[0]; }},
      {"map entry E", "map:permutation", [](Copy& c) { c.input_positions[1] = static_cast<std::int32_t>(c.input_positions.size()); }},
      {"slot out of range", "slot:permutation", [](Copy& c) { c.in_adj[2] = 1 << 30; }},
      {"negative slot", "slot:permutation", [](Copy& c) { c.in_adj[2] = -1; }},
      // Additional size, offset and copied-source corruption controls.
      {"truncated capacity", "size:capacity", [](Copy& c) { c.capacity.pop_back(); }},
      {"oversized slots", "size:in_adj_edge_index", [](Copy& c) { c.in_adj.push_back(0); }},
      {"truncated offsets", "size:in_row_offsets", [](Copy& c) { c.incoming_offsets.pop_back(); }},
      {"truncated copy", "size:in_col_indices", [](Copy& c) { c.in_col.pop_back(); }},
      {"negative offset", "row_offsets:range", [](Copy& c) { c.outgoing_offsets[1] = -1; }},
      {"excessive offset", "in_row_offsets:range", [](Copy& c) { c.incoming_offsets[2] = 1 << 30; }},
      {"wrong copied source", "copy:source", [](Copy& c) { c.in_col[0] = (c.in_col[0] + 1) % 4; }},
      {"copy indexed by id", "copy:source",
       [](Copy& c) {
         for (std::size_t s = 0; s < c.in_col.size(); ++s) c.in_col[s] = c.src[s];
       }},
      {"copied source out of range", "range:in_col_indices", [](Copy& c) { c.in_col[0] = 4; }},
      {"changed maximum", "scan:max_cost", [](Copy& c) { ++c.max_edge_cost; }},
      {"changed GCD", "scan:cost_gcd", [](Copy& c) { c.edge_cost_gcd = 2; }},
      // Every other code: each array's size, each offset rule of both offset arrays, the endpoint
      // ranges and each permuted property.
      {"truncated costs", "size:cost", [](Copy& c) { c.cost.pop_back(); }},
      {"truncated sources", "size:src", [](Copy& c) { c.src.pop_back(); }},
      {"truncated destinations", "size:dst", [](Copy& c) { c.dst.pop_back(); }},
      {"truncated map", "size:input_index", [](Copy& c) { c.input_positions.pop_back(); }},
      {"truncated row offsets", "size:row_offsets", [](Copy& c) { c.outgoing_offsets.pop_back(); }},
      {"row offsets start at 1", "row_offsets:start", [](Copy& c) { c.outgoing_offsets[0] = 1; }},
      {"reverse offsets start at 1", "in_row_offsets:start", [](Copy& c) { c.incoming_offsets[0] = 1; }},
      {"reverse offsets end early", "in_row_offsets:end", [](Copy& c) { --c.incoming_offsets.back(); }},
      {"row offsets decrease", "row_offsets:monotone", [](Copy& c) { c.outgoing_offsets[1] = c.outgoing_offsets[2] + 1; }},
      {"reverse offsets decrease", "in_row_offsets:monotone",
       [](Copy& c) { c.incoming_offsets[2] = c.incoming_offsets[1] - 1; }},
      {"source out of range", "range:src", [](Copy& c) { c.src[0] = 4; }},
      {"negative destination", "range:dst", [](Copy& c) { c.dst[0] = -1; }},
      // Two different values exchanged (the first one and the first that differs from it).
      {"costs swapped", "property:cost",
       [](Copy& c) { std::swap(c.cost[0], *std::find_if(c.cost.begin(), c.cost.end(), [&](Cost x) { return x != c.cost[0]; })); }},
      {"sources swapped", "property:src",
       [](Copy& c) { std::swap(c.src[0], *std::find_if(c.src.begin(), c.src.end(), [&](NodeId x) { return x != c.src[0]; })); }},
      {"destinations swapped", "property:dst",
       [](Copy& c) { std::swap(c.dst[0], *std::find_if(c.dst.begin(), c.dst.end(), [&](NodeId x) { return x != c.dst[0]; })); }},
  };
  for (const auto& [name, code, corrupt] : controls) {
    Copy c = clean;
    corrupt(c);
    const auto bad = check(c);
    EXPECT_NE(std::find(bad.begin(), bad.end(), code), bad.end()) << name << ": expected " << code;
  }
  // The counts and the caller's arrays: the node and edge counts, and each input array's size.
  auto spans = spans_of(g);
  ++spans.num_nodes;
  EXPECT_EQ(check_layout(spans, view(in)), std::vector<std::string> {"size:num_nodes"});
  spans = spans_of(g);
  ++spans.num_edges;
  EXPECT_EQ(check_layout(spans, view(in)), std::vector<std::string> {"size:num_edges"});
  // A negative node count fails even when the caller claims the same count.
  spans = spans_of(g);
  spans.num_nodes = -1;
  auto negative_input = view(in);
  negative_input.num_nodes = -1;
  EXPECT_EQ(check_layout(spans, negative_input), std::vector<std::string> {"size:num_nodes"});
  for (const std::string name : {"input.dst", "input.capacity", "input.cost"}) {
    auto short_input = view(in);
    if (name == "input.dst") short_input.dst = short_input.dst.first(short_input.dst.size() - 1);
    if (name == "input.capacity") short_input.capacity = short_input.capacity.first(short_input.capacity.size() - 1);
    if (name == "input.cost") short_input.cost = short_input.cost.first(short_input.cost.size() - 1);
    EXPECT_EQ(check_layout(spans_of(g), short_input), std::vector<std::string> {"size:" + name});
  }
  // Empty graphs are positive controls: all-zero offsets are correct without edges.
  for (std::int32_t n : {0, 1, 1000}) {
    const Input e {"empty", n};
    EXPECT_TRUE(check_layout(spans_of(build(e)), view(e)).empty()) << n;
  }
}

// The hash and fallback boundary at 65,536 distinct costs (the branch observation), on
// inputs that do not exhaust the probe budget. In this test and the five after it, the arrays are
// compared with the reference in every build; only the branch counters need the private hooks, and
// without them the test reports itself skipped after the comparisons have passed.
TEST(GraphLayout, HashAndFallbackBoundaryAt65536DistinctCosts) {
  for (const auto& [k, hashed] : std::vector<std::pair<std::int64_t, bool>> {{65535, true}, {65536, true}, {65537, false}}) {
    const auto in = exact_k(k, 1000, 4096, static_cast<std::uint64_t>(k));
    ASSERT_EQ(distinct(in.cost), k);
#if NETGRAPH_CORE_PRIVATE_HOOKS
    auto& counts = detail::thread_construction_counts();
    const auto before = counts;
    const Graph g = build(in);
    EXPECT_EQ(counts.hashed - before.hashed, hashed ? 1 : 0) << k;
    EXPECT_EQ(counts.fallback - before.fallback, hashed ? 0 : 1) << k;
    EXPECT_EQ(counts.distinct_limit - before.distinct_limit, hashed ? 0 : 1) << k;
    EXPECT_EQ(counts.exhausted_collection - before.exhausted_collection, 0) << k;
    EXPECT_EQ(counts.exhausted_growth - before.exhausted_growth, 0) << k;
    EXPECT_LE(counts.last_probes, counts.last_budget) << k;
    if (hashed) {
      EXPECT_EQ(counts.last_distinct, k);
      EXPECT_LE(counts.last_table_slots, std::int64_t {1} << 18);  // load <= 1/4 at K = 65,536
      EXPECT_GE(counts.last_table_slots, 4 * k);
    }
#else
    const Graph g = build(in);
#endif
    expect_same_as_reference(g, in);
  }
#if !NETGRAPH_CORE_PRIVATE_HOOKS
  GTEST_SKIP() << "branch observations need the private hooks";
#endif
}

// The combined-key boundary: with N K < E the passes run on combined
// keys; N K = E must take the five single-key passes. Both give the reference arrays, in every
// input order, and the random differential inputs reach both branches.
TEST(GraphLayout, CombinedKeyPassesOnlyBelowNTimesKEqualsE) {
  auto shaped = [](std::int32_t n, std::int64_t k, std::size_t e, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    Input in {"N=" + std::to_string(n) + " K=" + std::to_string(k) + " E=" + std::to_string(e), n};
    for (std::size_t i = 0; i < e; ++i) {
      in.src.push_back(static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n)));
      in.dst.push_back(static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n)));
      in.cap.push_back(static_cast<Amount>(1 + i));
      const auto j = i < static_cast<std::size_t>(k) ? static_cast<std::int64_t>(i) : static_cast<std::int64_t>(rng() % static_cast<std::uint64_t>(k));
      in.cost.push_back(3 * j + 1);
    }
    return in;
  };
  for (const auto& [n, k] : std::vector<std::pair<std::int32_t, std::int64_t>> {{5, 1}, {10, 3}, {400, 2}, {64, 31}}) {
    const auto nk = static_cast<std::size_t>(n) * static_cast<std::size_t>(k);
    for (const auto& [e, fused] : std::vector<std::pair<std::size_t, bool>> {{nk + 1, true}, {nk, false}}) {
      for (const auto& in : orders(shaped(n, k, e, nk + e), e)) {
        ASSERT_EQ(distinct(in.cost), k) << in.name;
#if NETGRAPH_CORE_PRIVATE_HOOKS
        auto& counts = detail::thread_construction_counts();
        const auto before = counts.fused;
        const Graph g = build(in);
        EXPECT_EQ(counts.fused - before, fused ? 1 : 0) << in.name;
#else
        const Graph g = build(in);
#endif
        expect_same_as_reference(g, in);
      }
    }
  }
#if NETGRAPH_CORE_PRIVATE_HOOKS
  std::int64_t fused = 0, plain = 0;
  for (std::uint64_t seed = 1; seed <= 200; ++seed) {
    const auto in = random_input("seed " + std::to_string(seed), 1 + static_cast<std::int32_t>(seed % 37),
                                 1 + seed % 211, 1 + static_cast<Cost>(seed % 13), seed);
    const auto before = detail::thread_construction_counts().fused;
    (void)build(in);
    (detail::thread_construction_counts().fused > before ? fused : plain) += 1;
  }
  EXPECT_GT(fused, 0);
  EXPECT_GT(plain, 0);
#else
  GTEST_SKIP() << "branch observations need the private hooks";
#endif
}

// The fallback's 16-bit digits: N = 2, two copies of [1, ..., 65,536,
// 2^59] (K = 65,537, E = 131,074 = N K: the single-key passes; the cost sum is below 2^62), then
// one more cost of 1 (E = N K + 1: the combined keys). Both take the distinct-limit fallback (observed); by the
// costs' values its radix sorts on the digits at 0, 16 and 48 bits and skips the one at 32 (one
// digit value), with arrays equal to the reference in four input orders.
TEST(GraphLayout, RadixFallbackDigitsOnBothSidesOfTheCombinedKeyBoundary) {
  Input base {"two copies of 1..65,536 and 2^59", 2};
  for (int copy = 0; copy < 2; ++copy) {
    for (Cost c = 1; c <= 65536; ++c) base.cost.push_back(c);
    base.cost.push_back(Cost {1} << 59);
  }
  Input plus = base;
  plus.name += " + one cost 1";
  plus.cost.push_back(1);
  for (Input* in : {&base, &plus}) {
    for (std::size_t i = 0; i < in->cost.size(); ++i) {
      in->src.push_back(static_cast<NodeId>(i % 2));
      in->dst.push_back(static_cast<NodeId>((i / 2) % 2));
      in->cap.push_back(static_cast<Amount>(1 + i));
    }
  }
  ASSERT_EQ(base.cost.size(), 131074u);
  for (const auto& [x, fused] : std::vector<std::pair<const Input*, bool>> {{&base, false}, {&plus, true}}) {
    ASSERT_EQ(distinct(x->cost), 65537);
    for (const auto& in : orders(*x, x->cost.size())) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
      auto& counts = detail::thread_construction_counts();
      const auto before = counts;
      const Graph g = build(in);
      EXPECT_EQ(counts.distinct_limit - before.distinct_limit, 1) << in.name;
      EXPECT_EQ(counts.fallback - before.fallback, 1) << in.name;
      EXPECT_EQ(counts.fused - before.fused, fused ? 1 : 0) << in.name;
      EXPECT_EQ(counts.last_distinct, 65537) << in.name;
#else
      const Graph g = build(in);
#endif
      expect_same_as_reference(g, in);
    }
  }
#if !NETGRAPH_CORE_PRIVATE_HOOKS
  GTEST_SKIP() << "branch observations need the private hooks";
#endif
}

// A run of equal costs reuses the previous ordinal without a lookup: one cost repeated costs
// exactly one probe; two alternating costs in different slots cost one probe per edge.
TEST(GraphLayout, RunsOfEqualCostsSkipTheLookup) {
  Cost a = 1, b = 2;
  while ((detail::cost_hash(a) & 63) == (detail::cost_hash(b) & 63)) ++b;
  for (const bool alternate : {false, true}) {
    Input in {alternate ? "alternating" : "one run", 9};
    for (int i = 0; i < 1000; ++i) {
      in.src.push_back(i % 9);
      in.dst.push_back((i * 4) % 9);
      in.cap.push_back(i);
      in.cost.push_back(alternate && i % 2 == 1 ? b : a);
    }
    const Graph g = build(in);
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(detail::thread_construction_counts().last_probes, alternate ? 1000 : 1) << in.name;
#endif
    expect_same_as_reference(g, in);
  }
#if !NETGRAPH_CORE_PRIVATE_HOOKS
  GTEST_SKIP() << "branch observations need the private hooks";
#endif
}

namespace {

enum class BudgetOutcome { Hashed, Collection, Growth };

struct ProbeCount {
  BudgetOutcome outcome {BudgetOutcome::Hashed};
  std::int64_t probes {0};          // every examined slot: exactly the budget when it ran out
  std::vector<std::int64_t> after;  // the probes after each hashed edge, in input order
  std::size_t distinct {0};         // the distinct costs inserted
  std::size_t slots {64};           // the table's size
};

// The probe accounting of the cost hashing, modelled from its definition rather than from the
// constructor: linear probing from 64 slots, every examined slot counted against the budget; a
// cost equal to the previous one reuses its ordinal without a lookup; before a new cost would
// raise the load above 1/4, the table doubles, the distinct costs are reinserted in
// first-appearance order and the new cost is looked up again. A probe beyond the budget stops the
// hashing, inside an input lookup (Collection) or inside a growth, its reinsertion or the repeated
// lookup (Growth). Without a budget every cost is hashed.
ProbeCount count_probes(std::span<const Cost> cost, std::optional<std::int64_t> budget) {
  ProbeCount out;
  std::vector<Cost> table(64, 0);
  std::vector<Cost> order;
  auto locate = [&](Cost c) -> std::optional<std::size_t> {
    for (std::size_t at = detail::cost_hash(c) & (table.size() - 1);; at = (at + 1) & (table.size() - 1)) {
      if (budget.has_value() && out.probes == *budget) return std::nullopt;
      ++out.probes;
      if (table[at] == 0 || table[at] == c) return at;
    }
  };
  Cost prev = 0;
  for (const Cost c : cost) {
    if (c != prev) {
      auto at = locate(c);
      if (!at) {
        out.outcome = BudgetOutcome::Collection;
        return out;
      }
      if (table[*at] == 0) {
        if (4 * (order.size() + 1) > table.size()) {
          table.assign(2 * table.size(), 0);
          for (const Cost d : order) {
            const auto slot = locate(d);
            if (!slot) {
              out.outcome = BudgetOutcome::Growth;
              return out;
            }
            table[*slot] = d;
          }
          at = locate(c);
          if (!at) {
            out.outcome = BudgetOutcome::Growth;
            return out;
          }
        }
        table[*at] = c;
        order.push_back(c);
      }
      prev = c;
    }
    out.after.push_back(out.probes);
  }
  out.distinct = order.size();
  out.slots = table.size();
  return out;
}

// The construction's budget B = 8E + 8 x 65,536 for E edges.
std::int64_t probe_budget(std::size_t edges) { return 8 * static_cast<std::int64_t>(edges) + 8 * 65536; }

BudgetOutcome budget_outcome(std::span<const Cost> cost) { return count_probes(cost, probe_budget(cost.size())).outcome; }

} // namespace

// The cumulative probe budget B = 8E + 8 x 65,536 over the whole construction. Costs whose hash has
// 11 low zero bits start at slot 0 in every table up to 2,048 slots, so they form one cluster:
// the j-th of them (from 0) is found after j + 1 probes. After the 64 cluster costs come D
// duplicates alternating z[62], z[63] (a run of one cost reuses its ordinal without a lookup, so
// the duplicates must alternate), each costing 63 or 64 probes while B grows by 8 per edge. The
// boundaries are those of budget_outcome above (a model of the probe accounting, independent
// of the constructor), which every case checks first:
//  * collection: exhausted at an input lookup exactly when D >= 9,406;
//  * growth: with a 65th cluster cost at the end (its lookup, then the reinsertion into 512
//    slots, then the lookup repeated), exhausted inside the growth for 9,366 <= D <= 9,404, not at
//    all for D <= 9,365, and at the 65th cost's first lookup (collection) for D = 9,405.
// Exhausted constructions use exactly B probes and none more; the arrays equal the reference.
TEST(GraphLayout, ProbeBudgetFallbackInCollectionAndGrowth) {
  auto shape = [](std::size_t dups, bool tail) {
    return budget_shape("cluster D=" + std::to_string(dups) + (tail ? " + 65th" : ""), dups, tail);
  };
  enum Where { None, Collection, Growth };
  const std::vector<std::tuple<std::size_t, bool, Where>> cases = {
      {9000, false, None},  {9405, false, None},  {9406, false, Collection}, {10000, false, Collection},
      {9300, true, None},   {9365, true, None},   {9366, true, Growth},       {9380, true, Growth},
      {9404, true, Growth}, {9405, true, Collection},
  };

  for (const auto& [dups, tail, where] : cases) {
    const auto in = shape(dups, tail);
    const BudgetOutcome modelled[] = {BudgetOutcome::Hashed, BudgetOutcome::Collection, BudgetOutcome::Growth};
    EXPECT_EQ(budget_outcome(in.cost), modelled[where]) << in.name;
#if NETGRAPH_CORE_PRIVATE_HOOKS
    auto& counts = detail::thread_construction_counts();
    const auto before = counts;
    const Graph g = build(in);
    const std::int64_t budget = probe_budget(in.src.size());
    EXPECT_EQ(counts.last_budget, budget) << in.name;
    EXPECT_LE(counts.last_probes, budget) << in.name;  // no lookup examines more than B in total
    EXPECT_EQ(counts.exhausted_collection - before.exhausted_collection, where == Collection ? 1 : 0) << in.name;
    EXPECT_EQ(counts.exhausted_growth - before.exhausted_growth, where == Growth ? 1 : 0) << in.name;
    EXPECT_EQ(counts.fallback - before.fallback, where == None ? 0 : 1) << in.name;
    EXPECT_EQ(counts.hashed - before.hashed, where == None ? 1 : 0) << in.name;
    if (where != None) EXPECT_EQ(counts.last_probes, budget) << in.name;
#else
    const Graph g = build(in);
#endif
    expect_same_as_reference(g, in);
  }
#if !NETGRAPH_CORE_PRIVATE_HOOKS
  GTEST_SKIP() << "branch observations need the private hooks";
#endif
}

// The budget stays cumulative across the large table's growths: 64 costs share their
// low 18 hash bits, so they collide in every table up to 262,144 slots; further costs grow the
// table to K = 32,769; then the last two colliding costs alternate. The probe model above
// (count_probes, independent of the constructor) finds the first edge count whose probes exceed
// B: the collection exhausts the budget there, and one edge fewer stays hashed with exactly the
// counted probes.
TEST(GraphLayout, ProbeBudgetIsCumulativeAcrossTheLargeTableGrowths) {
  const std::vector<Cost> z = colliding_costs(262143, 64, 0);
  std::vector<Cost> seq = z;
  for (Cost c = 1; seq.size() < 32769; ++c) {
    if (!std::binary_search(z.begin(), z.end(), c)) seq.push_back(c);
  }
  const std::size_t distinct_end = seq.size();
  for (std::size_t d = 0; d < 20000; ++d) seq.push_back(z[d % 2 == 0 ? 62 : 63]);
  const ProbeCount unbounded = count_probes(seq, std::nullopt);
  ASSERT_EQ(unbounded.distinct, 32769U);
  ASSERT_EQ(unbounded.slots, 262144U);
  std::size_t first = 0;  // the first edge count whose probes exceed its own budget
  for (std::size_t e = 1; e <= seq.size() && first == 0; ++e) {
    if (unbounded.after[e - 1] > probe_budget(e)) first = e;
  }
  ASSERT_GT(first, distinct_end);  // inside the alternation
  seq.resize(first);
  const std::int64_t prefix_probes = unbounded.after[first - 2];  // all probes of one edge fewer
  // The same model run with each length's budget stops there and not one edge earlier.
  EXPECT_EQ(count_probes(seq, probe_budget(first)).outcome, BudgetOutcome::Collection);
  const ProbeCount shorter = count_probes(std::span<const Cost>(seq).first(first - 1), probe_budget(first - 1));
  EXPECT_EQ(shorter.outcome, BudgetOutcome::Hashed);
  EXPECT_EQ(shorter.probes, prefix_probes);
  EXPECT_EQ(seq.size(), 44436U);  // the count of the probe model

  auto input = [&](std::size_t e) {
    Input in {"late growth E=" + std::to_string(e), 31};
    for (std::size_t i = 0; i < e; ++i) {
      in.src.push_back(static_cast<NodeId>(i % 31));
      in.dst.push_back(static_cast<NodeId>((i * 7) % 31));
      in.cap.push_back(static_cast<Amount>(i % 13));
      in.cost.push_back(seq[i]);
    }
    return in;
  };
  for (const std::size_t e : {seq.size() - 1, seq.size()}) {
    const auto in = input(e);
#if NETGRAPH_CORE_PRIVATE_HOOKS
    const bool exhausts = e == seq.size();
    auto& counts = detail::thread_construction_counts();
    const auto before = counts;
    const Graph g = build(in);
    const std::int64_t budget = probe_budget(e);
    EXPECT_EQ(counts.last_budget, budget) << in.name;
    EXPECT_EQ(counts.last_probes, exhausts ? budget : prefix_probes) << in.name;
    EXPECT_EQ(counts.exhausted_collection - before.exhausted_collection, exhausts ? 1 : 0) << in.name;
    EXPECT_EQ(counts.exhausted_growth - before.exhausted_growth, 0) << in.name;
    EXPECT_EQ(counts.hashed - before.hashed, exhausts ? 0 : 1) << in.name;
    if (!exhausts) EXPECT_EQ(counts.last_distinct, 32769) << in.name;
#else
    (void)prefix_probes;
    const Graph g = build(in);
#endif
    expect_same_as_reference(g, in);
  }
#if !NETGRAPH_CORE_PRIVATE_HOOKS
  GTEST_SKIP() << "branch observations need the private hooks";
#endif
}

// A copy keeps every array, the map and the resolution (one input at a resolution other than 1)
// and draws a new uid.
TEST(GraphLayout, CopyKeepsEveryArrayAndDrawsANewUid) {
  Input fine = exact_k(31, 500, 64, 6);
  fine.name += " at 2^-12";
  fine.q = std::ldexp(1.0, -12);
  for (const auto& in : {exact_k(31, 500, 64, 1), exact_k(65537, 10, 512, 2), Input {"empty", 5}, fine}) {
    const Graph g = build(in);
    EXPECT_GT(g.uid(), 0U);  // a constructed graph's uid is positive (0 marks a moved-from one)
    const Graph c(g);
    EXPECT_NE(c.uid(), g.uid());
    EXPECT_GT(c.uid(), g.uid());  // a relation, not an absolute counter value
    expect_same_as_reference(c, in);
    EXPECT_EQ(GraphAccess::retained(c).size_bytes, GraphAccess::retained(g).size_bytes);
    // A second construction of the same input is a new graph too.
    EXPECT_NE(build(in).uid(), g.uid());
  }
}

// The move: every value and the uid transferred without allocation
// (GraphFaults.MoveIsAllocationFreeUnderAnArmedFailure checks that none is made); the source keeps its resolution, has uid, N, E, maximum and
// GCD zero, empty edge arrays and one zero in each offset view. Copying and moving an already
// moved-from graph give the same state.
TEST(GraphLayout, MoveTransfersEverythingAndLeavesAnEmptyGraph) {
  const auto in = exact_k(31, 500, 64, 3);
  const double q = std::ldexp(1.0, -12);
  Graph g(in.n, in.src, in.dst, in.cap, in.cost, q);
  const auto uid = g.uid();
  const Graph m(std::move(g));
  EXPECT_EQ(m.uid(), uid);
  EXPECT_EQ(m.resolution(), q);
  Input with_q = in;
  with_q.q = q;
  expect_same_as_reference(m, with_q);
  auto expect_moved_from = [&](const Graph& x, const char* what) {
    EXPECT_EQ(x.resolution(), q) << what;
    EXPECT_EQ(x.num_nodes(), 0) << what;
    EXPECT_EQ(x.num_edges(), 0) << what;
    EXPECT_EQ(x.max_edge_cost(), 0) << what;
    EXPECT_EQ(x.edge_cost_gcd(), 0) << what;
    EXPECT_TRUE(x.edge_capacity_units().empty() && x.edge_costs().empty() && x.edge_sources().empty() && x.edge_destinations().empty() &&
                x.input_positions().empty() && x.incoming_edges().empty() && GraphAccess::incoming_edge_sources(x).empty())
        << what;
    EXPECT_EQ(vec(x.outgoing_offsets()), (std::vector<std::int32_t> {0})) << what;
    EXPECT_EQ(vec(x.incoming_offsets()), (std::vector<std::int32_t> {0})) << what;
    EXPECT_TRUE(check_layout(spans_of(x), view(Input {"moved-from", 0})).empty()) << what;
  };
  // NOLINTBEGIN(bugprone-use-after-move): the documented moved-from state.
  EXPECT_EQ(g.uid(), 0u);
  expect_moved_from(g, "source");
  const Graph copied(g);
  EXPECT_GT(copied.uid(), uid);  // a copy is a new graph, even of an empty one
  expect_moved_from(copied, "copy of the source");  // except its uid
  Graph again(std::move(g));
  EXPECT_EQ(again.uid(), 0u);
  expect_moved_from(again, "moved again");
  expect_moved_from(g, "source moved twice");
  // NOLINTEND(bugprone-use-after-move)
}

// The retained array payload: 36 E + 8 (N + 1) bytes, capacities equal to sizes (no slack, no
// retained scratch, no inverse map, no physical-capacity array). The sizes come from the graph's
// own list of its arrays (GraphAccess::retained), whose length is fixed by that list's code and
// so not compared here; that the list is complete, nine arrays and no other allocation held, is
// GraphFaults' count of the allocations a graph holds.
TEST(GraphLayout, RetainedPayloadIs36BytesPerEdgePlusOffsets) {
  std::vector<Input> cases = {Input {"N=0", 0}, Input {"isolated", 1000}, exact_k(31, 4000, 512, 4),
                              exact_k(65537, 100, 4096, 5)};
  for (const auto& in : cases) {
    const Graph g = build(in);
    const Graph c(g);
    Graph moved_from = build(in);
    const Graph m(std::move(moved_from));
    const std::size_t e = in.src.size();
    const std::size_t want = 36 * e + 8 * (static_cast<std::size_t>(in.n) + 1);
    for (const Graph* x : {&g, &c, &m}) {
      const auto r = GraphAccess::retained(*x);
      EXPECT_EQ(r.size_bytes, want) << in.name;
      EXPECT_EQ(r.capacity_bytes, want) << in.name;
    }
  }
}

// Graph-bound objects from another graph are rejected at their existing use sites, on a twin
// graph with identical arrays (only the uid differs, so the check is by identity, not by shape):
// the evaluators' Eligibility and ShortestPathCache checks, and a Path's DAG (to_dag keeps the Path's
// graph) in resolve_to_paths. The flow core, ledger, policy, KSP, search and cache sites have
// their own tests.
TEST(GraphLayout, ForeignGraphBoundObjectsAreRejectedAtTheirUseSites) {
  const std::vector<NodeId> src {0, 0, 1, 2}, dst {1, 2, 3, 3};
  const std::vector<Amount> cap {5, 5, 5, 5};
  const std::vector<Cost> cost {1, 1, 1, 1};
  const Graph g(4, src, dst, cap, cost, 1.0), twin(4, src, dst, cap, cost, 1.0);
  ASSERT_NE(g.uid(), twin.uid());
  const std::vector<NodeId> s {0}, d {3};
  const FlowPolicyConfig cfg {shortest_paths(RouteAllocation::Optimized), Admission::Lossless};
  const Eligibility el = twin.eligibility({}, std::vector<EdgeId> {0});
  ShortestPathCache cache(twin);
  const std::string elig_msg = "foreign Eligibility: built for another graph";
  const std::string cache_msg = "foreign ShortestPathCache: bound to another graph";
  EXPECT_EQ(error_of([&] { (void)evaluate_demand(g, s, d, EvaluationOptions {cfg, &el}); }), elig_msg);
  EXPECT_EQ(error_of([&] { (void)evaluate_demand(g, s, d, EvaluationOptions {cfg, nullptr, &cache}); }), cache_msg);
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs {{s, d}};
  EXPECT_EQ(error_of([&] { (void)evaluate_pairs(g, pairs, EvaluationOptions {cfg, nullptr, &cache}, {}, 1); }), cache_msg);
  EXPECT_EQ(error_of([&] { (void)evaluate_edge_failures(g, s, d, EvaluationOptions {cfg, &el}, 1); }), elig_msg);
  EXPECT_EQ(error_of([&] { (void)evaluate_edge_failures(g, s, d, EvaluationOptions {cfg, nullptr, &cache}, 1); }), cache_msg);
  const Path p(twin, {{0}, {2}});
  EXPECT_EQ(error_of([&] { (void)resolve_to_paths(g, p.to_dag(), 0, 3); }), "foreign PredecessorDAG: built for another graph");
  // The same calls on the owning graph succeed (the controls).
  EXPECT_NO_THROW((void)evaluate_demand(twin, s, d, EvaluationOptions {cfg, &el, &cache}));
  EXPECT_EQ(resolve_to_paths(twin, p.to_dag(), 0, 3).size(), 1U);
}

TEST(GraphLayout, ReferenceComparisonRejectsDifferentInputs) {
  const Input in {"comparison", 4, {0, 1, 0, 2}, {1, 3, 2, 3}, {2, 3, 4, 5}, {2, 4, 6, 8}};
  const auto g = build(in);
  EXPECT_EQ(failures_of([&] { expect_same_as_reference(g, in); }), 0);
  const std::vector<std::function<void(Input&)>> changes = {
      [](auto& x) { ++x.n; }, [](auto& x) { x.q = 2; },
      [](auto& x) { ++x.cap.back(); }, [](auto& x) { x.cost.back() = 9; },
      [](auto& x) { x.src.back() = 1; }, [](auto& x) { x.dst.back() = 0; },
      [](auto& x) { x.src.push_back(3); x.dst.push_back(0); x.cap.push_back(1); x.cost.push_back(1); },
      [](auto& x) { std::swap(x.cost[0], x.cost[2]); },
  };
  for (std::size_t i = 0; i < changes.size(); ++i) {
    auto different = in;
    changes[i](different);
    EXPECT_GT(failures_of([&] { expect_same_as_reference(g, different); }), 0) << i;
  }
}
