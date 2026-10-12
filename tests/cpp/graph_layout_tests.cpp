// Layout B' storage and construction (APPROACH 3; task 1e): the independent layout checker on
// every graph, its corruption controls, the differential test against the pinned 1d constructor
// (every array, every error), the hash/fallback boundary, copies and moves, and the retained
// array payload. The manifest instance 22@1e.storage names 13 obligations; each test below names
// the one it covers ("manifest N").
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <typeinfo>
#include <vector>

#include "graph_reference_1d.hpp"
#include "layout_checker.hpp"
#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/detail/graph_access.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/search.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;
using detail::GraphAccess;

namespace {

constexpr Amount L62 = Amount{1} << 62;

struct Input {
  std::string name;
  std::int32_t n {0};
  std::vector<NodeId> src {}, dst {};
  std::vector<Amount> cap {};
  std::vector<Cost> cost {};
  double q {1.0};
};

LayoutInput view(const Input& in) { return {in.n, in.src, in.dst, in.cap, in.cost}; }

LayoutSpans spans_of(const Graph& g) {
  return {g.num_nodes(), g.num_edges(), g.capacity_units(), g.cost(), g.edge_src(), g.edge_dst(), g.input_index(),
          g.row_offsets(), g.in_row_offsets(), g.in_adj_edge_index(), GraphAccess::in_col_indices(g), g.max_cost(),
          g.cost_gcd()};
}

template <typename T>
std::vector<T> vec(std::span<const T> s) {
  return {s.begin(), s.end()};
}

// Every public array, its element width, the counts, the resolution and the scans equal the
// reference's; the internal copy equals the reference's edge_src[in_adj_edge_index[slot]].
void expect_same_as_reference(const Graph& g, const Input& in) {
  const auto r = reference_1d(in.n, in.src, in.dst, in.cap, in.cost, in.q);
  static_assert(sizeof(decltype(g.capacity_units())::element_type) == 8);
  static_assert(sizeof(decltype(g.cost())::element_type) == 8);
  static_assert(sizeof(decltype(g.edge_src())::element_type) == 4);
  static_assert(sizeof(decltype(g.input_index())::element_type) == 4);
  static_assert(sizeof(decltype(g.in_adj_edge_index())::element_type) == 4);
  EXPECT_EQ(g.num_nodes(), r.num_nodes) << in.name;
  EXPECT_EQ(g.num_edges(), static_cast<std::int32_t>(r.capacity.size())) << in.name;
  EXPECT_EQ(g.resolution(), r.resolution) << in.name;
  EXPECT_EQ(g.max_cost(), r.max_cost) << in.name;
  EXPECT_EQ(g.cost_gcd(), r.cost_gcd) << in.name;
  EXPECT_EQ(vec(g.capacity_units()), r.capacity) << in.name;
  EXPECT_EQ(vec(g.cost()), r.cost) << in.name;
  EXPECT_EQ(vec(g.edge_src()), r.src) << in.name;
  EXPECT_EQ(vec(g.edge_dst()), r.dst) << in.name;
  EXPECT_EQ(vec(g.input_index()), r.input_index) << in.name;
  EXPECT_EQ(vec(g.row_offsets()), r.row_offsets) << in.name;
  EXPECT_EQ(vec(g.in_row_offsets()), r.in_row_offsets) << in.name;
  EXPECT_EQ(vec(g.in_adj_edge_index()), r.in_adj_edge_index) << in.name;
  EXPECT_EQ(vec(GraphAccess::in_col_indices(g)), reference_in_col_indices(r)) << in.name;
  const auto bad = check_layout(spans_of(g), view(in));
  EXPECT_TRUE(bad.empty()) << in.name << ": " << (bad.empty() ? "" : bad.front());
}

Graph build(const Input& in) { return Graph(in.n, in.src, in.dst, in.cap, in.cost, in.q); }

std::int64_t distinct(std::vector<Cost> c) {
  std::sort(c.begin(), c.end());
  return std::unique(c.begin(), c.end()) - c.begin();
}

// --------------------------------------------------------------------------- generators

Input random_input(std::string name, std::int32_t n, std::size_t m, Cost max_cost, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  Input in {std::move(name), n};
  for (std::size_t i = 0; i < m; ++i) {
    in.src.push_back(static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n)));
    in.dst.push_back(static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n)));
    in.cap.push_back(static_cast<Amount>(rng() % 5));  // zero capacities included
    in.cost.push_back(static_cast<Cost>(1 + rng() % static_cast<std::uint64_t>(max_cost)));
  }
  return in;
}

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

// The same edges in B' order, reversed and shuffled (all property arrays permuted together).
std::vector<Input> orders(const Input& in, std::uint64_t seed) {
  const auto r = reference_1d(in.n, in.src, in.dst, in.cap, in.cost, in.q);
  std::vector<std::size_t> sorted(r.input_index.begin(), r.input_index.end());
  std::vector<std::size_t> rev(in.src.size());
  std::iota(rev.rbegin(), rev.rend(), 0);
  std::vector<std::size_t> sh(in.src.size());
  std::iota(sh.begin(), sh.end(), 0);
  std::shuffle(sh.begin(), sh.end(), std::mt19937_64(seed));
  return {in, permuted(in, sorted, " (B' order)"), permuted(in, rev, " (reversed)"), permuted(in, sh, " (shuffled)")};
}

// Exactly k distinct costs, each present, plus `extra` repeats; distinct capacities.
Input exact_k(std::int64_t k, std::size_t extra, std::int32_t n, std::uint64_t seed, Cost stride = 3) {
  std::mt19937_64 rng(seed);
  Input in {"K=" + std::to_string(k), n};
  for (std::int64_t c = 1; c <= k; ++c) in.cost.push_back(c * stride + 1);
  for (std::size_t j = 0; j < extra; ++j) in.cost.push_back(in.cost[rng() % in.cost.size()]);
  std::shuffle(in.cost.begin(), in.cost.end(), rng);
  for (std::size_t i = 0; i < in.cost.size(); ++i) {
    in.src.push_back(static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n)));
    in.dst.push_back(static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n)));
    in.cap.push_back(static_cast<Amount>(1 + i));
  }
  return in;
}

// Costs that share a table slot under the given mask of the construction's hash.
std::vector<Cost> colliding_costs(std::uint64_t mask, std::size_t count) {
  std::vector<Cost> out;
  for (Cost c = 1; out.size() < count; ++c) {
    if ((detail::cost_hash(c) & mask) == 5) out.push_back(c);
  }
  return out;
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

TEST(GraphLayout, DifferentialWideIntegers) {
  // Costs and units above 2^53, sums just below 2^62; separate from the large-K cases (the cost
  // sum bound forbids combining them).
  const Cost big = Cost{1} << 59;  // four of them sum to 2^61: below the bound
  Input wide {"wide costs", 3, {0, 0, 1, 2, 0}, {1, 1, 2, 0, 1}, {(Amount{1} << 53) + 1, Amount{1} << 61, 0, 5, 3},
              {big + 3, big + 1, big + 2, 1, big + 1}};
  expect_same_as_reference(build(wide), wide);
  Input sums {"sums just below 2^62", 2, {0, 0}, {1, 1}, {L62 / 2, L62 / 2 - 1}, {L62 / 2, L62 / 2 - 1}};
  expect_same_as_reference(build(sums), sums);
  Input units {"units above 2^53", 2, {0, 1}, {1, 0}, {(Amount{1} << 53) + 1, (Amount{1} << 53) + 3}, {1, 1}};
  expect_same_as_reference(build(units), units);
}

TEST(GraphLayout, DifferentialHashCollisions) {
  // Costs colliding in the initial table (64 slots) and in the full one (2^18 slots at load 1/4).
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

// The defects of APPROACH 3's validation steps, each applied to a copy of a valid input.
using Defect = std::function<void(Input&)>;

std::vector<std::pair<std::string, Defect>> defects() {
  return {
      {"resolution", [](Input& x) { x.q = 3.0; }},
      {"capacity negative at 2", [](Input& x) { x.cap[2] = -1; }},
      {"capacity 2^62 at 1", [](Input& x) { x.cap[1] = L62; }},
      {"num_nodes negative", [](Input& x) { x.n = -1; }},
      {"dst shorter", [](Input& x) { x.dst.pop_back(); }},
      {"cost longer", [](Input& x) { x.cost.push_back(1); }},
      {"cost zero at 3", [](Input& x) { x.cost[3] = 0; }},
      {"cost negative at 1", [](Input& x) { x.cost[1] = -4; }},
      {"src out of range at 2", [](Input& x) { x.src[2] = x.n; }},
      {"dst negative at 1", [](Input& x) { x.dst[1] = -1; }},
      {"capacity sum", [](Input& x) { x.cap[0] = L62 / 2, x.cap[3] = L62 / 2; }},
      {"cost sum", [](Input& x) { x.cost[1] = L62 / 2, x.cost[2] = L62 / 2; }},
  };
}

} // namespace

TEST(GraphLayout, DifferentialErrorsSingleAndPairwise) {
  const Input base {"base", 4, {0, 1, 2, 3, 0}, {1, 2, 3, 0, 2}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}};
  const auto ds = defects();
  std::size_t checked = 0;
  for (std::size_t a = 0; a < ds.size(); ++a) {
    for (std::size_t b = a; b < ds.size(); ++b) {
      Input in = base;
      ds[a].second(in);
      if (b != a) ds[b].second(in);
      const auto got = outcome([&] { (void)build(in); });
      const auto want = outcome([&] { (void)reference_1d(in.n, in.src, in.dst, in.cap, in.cost, in.q); });
      ASSERT_FALSE(want.first.empty()) << ds[a].first << " + " << ds[b].first;
      EXPECT_EQ(got, want) << ds[a].first << " + " << ds[b].first;
      ++checked;
    }
  }
  EXPECT_EQ(checked, ds.size() * (ds.size() + 1) / 2);
}

TEST(GraphLayout, DifferentialErrorsFirstInvalidIndexAndSumOrder) {
  // Two defects of one kind: the lower input index is named.
  Input in {"two costs", 3, {0, 1, 2, 0}, {1, 2, 0, 2}, {1, 1, 1, 1}, {1, -2, 0, 1}};
  auto got = outcome([&] { (void)build(in); });
  EXPECT_EQ(got, outcome([&] { (void)reference_1d(in.n, in.src, in.dst, in.cap, in.cost, in.q); }));
  EXPECT_NE(got.second.find("cost[1] = -2"), std::string::npos) << got.second;
  // Capacity and cost sums both over: the one whose running sum crosses first, in input order.
  Input both {"sums", 2, {0, 0, 0}, {1, 1, 1}, {1, L62 / 2, L62 / 2}, {L62 / 2, L62 / 2, 1}};
  got = outcome([&] { (void)build(both); });
  EXPECT_EQ(got, outcome([&] { (void)reference_1d(both.n, both.src, both.dst, both.cap, both.cost, both.q); }));
  EXPECT_NE(got.second.find("edge costs"), std::string::npos) << got.second;
}

// --------------------------------------------------------------------------- manifest obligations

// manifest 1-8: the corruption controls of the checker, on a graph built by the constructor.
// The fixture repeats (src, cost, dst) keys with distinct capacities and input positions, and its
// reverse slot order differs from its id order, with different source values.
TEST(GraphLayout, CheckerRejectsEveryCorruption) {
  const Input in {"controls", 4, {0, 2, 0, 1, 0, 3, 2, 0}, {1, 1, 1, 1, 1, 1, 3, 3}, {5, 6, 7, 8, 9, 1, 2, 3},
                  {2, 1, 2, 1, 2, 1, 4, 4}};
  const Graph g = build(in);
  ASSERT_TRUE(check_layout(spans_of(g), view(in)).empty());
  struct Copy {
    std::vector<Amount> capacity;
    std::vector<Cost> cost;
    std::vector<NodeId> src, dst;
    std::vector<std::int32_t> input_index, row_offsets, in_row_offsets;
    std::vector<EdgeId> in_adj;
    std::vector<NodeId> in_col;
    Cost max_cost, cost_gcd;
  };
  const Copy clean {vec(g.capacity_units()), vec(g.cost()), vec(g.edge_src()), vec(g.edge_dst()), vec(g.input_index()),
                    vec(g.row_offsets()), vec(g.in_row_offsets()), vec(g.in_adj_edge_index()),
                    vec(GraphAccess::in_col_indices(g)), g.max_cost(), g.cost_gcd()};
  bool reverse_differs = false;
  for (std::size_t s = 0; s < clean.in_adj.size(); ++s) reverse_differs |= clean.in_adj[s] != static_cast<EdgeId>(s);
  ASSERT_TRUE(reverse_differs);
  auto check = [&](const Copy& c) {
    return check_layout({in.n, static_cast<std::int64_t>(in.src.size()), c.capacity, c.cost, c.src, c.dst, c.input_index,
                         c.row_offsets, c.in_row_offsets, c.in_adj, std::span<const NodeId>(c.in_col), c.max_cost,
                         c.cost_gcd},
                        view(in));
  };
  const std::vector<std::tuple<std::string, std::string, std::function<void(Copy&)>>> controls = {
      // manifest 1: wrong source or row boundary
      {"wrong source", "forward:membership", [](Copy& c) { c.src[0] = 1, c.input_index[0] = c.input_index[0]; }},
      {"row boundary moved", "forward:membership", [](Copy& c) { ++c.row_offsets[1]; }},
      // manifest 2: all-zero offsets (a graph with edges)
      {"all-zero offsets", "row_offsets:end", [](Copy& c) { std::fill(c.row_offsets.begin(), c.row_offsets.end(), 0); }},
      // manifest 3: missing or duplicated reverse coverage
      {"duplicated reverse slot", "slot:permutation", [](Copy& c) { c.in_adj[1] = c.in_adj[0]; }},
      {"wrong reverse destination", "reverse:membership", [](Copy& c) { ++c.in_row_offsets[1]; }},
      // manifest 4: (src, cost) substituted for (cost, src)
      {"(src, cost) reverse order", "reverse:order",
       [](Copy& c) {
         for (std::size_t v = 0; v + 1 < c.in_row_offsets.size(); ++v) {
           std::sort(c.in_adj.begin() + c.in_row_offsets[v], c.in_adj.begin() + c.in_row_offsets[v + 1],
                     [&](EdgeId p, EdgeId q) {
                       return std::tuple(c.src[static_cast<std::size_t>(p)], c.cost[static_cast<std::size_t>(p)], p) <
                              std::tuple(c.src[static_cast<std::size_t>(q)], c.cost[static_cast<std::size_t>(q)], q);
                     });
         }
         for (std::size_t s = 0; s < c.in_col.size(); ++s) c.in_col[s] = c.src[static_cast<std::size_t>(c.in_adj[s])];
       }},
      // manifest 5: unstable equal-key members
      {"unstable equal keys", "forward:order",
       [](Copy& c) {
         for (std::size_t e = 0; e + 1 < c.cost.size(); ++e) {
           if (c.src[e] == c.src[e + 1] && c.cost[e] == c.cost[e + 1] && c.dst[e] == c.dst[e + 1]) {
             std::swap(c.input_index[e], c.input_index[e + 1]);
             std::swap(c.capacity[e], c.capacity[e + 1]);
             return;
           }
         }
         FAIL() << "the fixture has no equal-key pair";
       }},
      // manifest 6: corrupted map
      {"map entries swapped", "property:capacity", [](Copy& c) { std::swap(c.input_index[0], c.input_index[1]); }},
      // manifest 7: property-only mispermutation
      {"capacities swapped", "property:capacity", [](Copy& c) { std::swap(c.capacity[1], c.capacity[2]); }},
      // manifest 8: out-of-range map or slot before access
      {"negative map entry", "map:permutation", [](Copy& c) { c.input_index[1] = -3; }},
      {"map entry E", "map:permutation", [](Copy& c) { c.input_index[1] = static_cast<std::int32_t>(c.input_index.size()); }},
      {"slot out of range", "slot:permutation", [](Copy& c) { c.in_adj[2] = 1 << 30; }},
      {"negative slot", "slot:permutation", [](Copy& c) { c.in_adj[2] = -1; }},
      // the brief's additional controls
      {"truncated capacity", "size:capacity", [](Copy& c) { c.capacity.pop_back(); }},
      {"oversized slots", "size:in_adj_edge_index", [](Copy& c) { c.in_adj.push_back(0); }},
      {"truncated offsets", "size:in_row_offsets", [](Copy& c) { c.in_row_offsets.pop_back(); }},
      {"truncated copy", "size:in_col_indices", [](Copy& c) { c.in_col.pop_back(); }},
      {"negative offset", "row_offsets:range", [](Copy& c) { c.row_offsets[1] = -1; }},
      {"excessive offset", "in_row_offsets:range", [](Copy& c) { c.in_row_offsets[2] = 1 << 30; }},
      {"wrong copied source", "copy:source", [](Copy& c) { c.in_col[0] = (c.in_col[0] + 1) % 4; }},
      {"copy indexed by id", "copy:source",
       [](Copy& c) {
         for (std::size_t s = 0; s < c.in_col.size(); ++s) c.in_col[s] = c.src[s];
       }},
      {"copied source out of range", "range:in_col_indices", [](Copy& c) { c.in_col[0] = 4; }},
      {"changed maximum", "scan:max_cost", [](Copy& c) { ++c.max_cost; }},
      {"changed GCD", "scan:cost_gcd", [](Copy& c) { c.cost_gcd = 2; }},
  };
  for (const auto& [name, code, corrupt] : controls) {
    Copy c = clean;
    corrupt(c);
    const auto bad = check(c);
    EXPECT_NE(std::find(bad.begin(), bad.end(), code), bad.end()) << name << ": expected " << code;
  }
  // Empty graphs are positive controls: all-zero offsets are correct without edges.
  for (std::int32_t n : {0, 1, 1000}) {
    const Input e {"empty", n};
    EXPECT_TRUE(check_layout(spans_of(build(e)), view(e)).empty()) << n;
  }
}

// manifest 9: hash and fallback boundary at 65,536 distinct costs (the branch observation), on
// inputs that do not exhaust the probe budget.
TEST(GraphLayout, HashAndFallbackBoundaryAt65536DistinctCosts) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  for (const auto& [k, hashed] : std::vector<std::pair<std::int64_t, bool>> {{65535, true}, {65536, true}, {65537, false}}) {
    const auto in = exact_k(k, 1000, 4096, static_cast<std::uint64_t>(k));
    ASSERT_EQ(distinct(in.cost), k);
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
    expect_same_as_reference(g, in);
  }
#else
  GTEST_SKIP() << "branch observations need the private hooks";
#endif
}

// Proposal P4: with N K < E the passes run on combined keys; N K = E must take the five single-key
// passes. Both give the reference arrays, in every input order, and the random differential
// inputs reach both branches.
TEST(GraphLayout, CombinedKeyPassesOnlyBelowNTimesKEqualsE) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
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
        auto& counts = detail::thread_construction_counts();
        const auto before = counts.fused;
        const Graph g = build(in);
        EXPECT_EQ(counts.fused - before, fused ? 1 : 0) << in.name;
        expect_same_as_reference(g, in);
      }
    }
  }
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

// The fallback's 16-bit digits (review-1e-consult/REPORT-2): N = 2, two copies of [1, ..., 65,536,
// 2^59] (K = 65,537, E = 131,074 = N K: the single-key passes; the cost sum is below 2^62), then
// one more cost of 1 (E = N K + 1: P4). Both take the distinct-limit fallback (observed); by the
// costs' values its radix sorts on the digits at 0, 16 and 48 bits and skips the one at 32 (one
// digit value), with arrays equal to the reference in four input orders.
TEST(GraphLayout, RadixFallbackDigitsOnBothSidesOfP4) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
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
      auto& counts = detail::thread_construction_counts();
      const auto before = counts;
      const Graph g = build(in);
      EXPECT_EQ(counts.distinct_limit - before.distinct_limit, 1) << in.name;
      EXPECT_EQ(counts.fallback - before.fallback, 1) << in.name;
      EXPECT_EQ(counts.fused - before.fused, fused ? 1 : 0) << in.name;
      EXPECT_EQ(counts.last_distinct, 65537) << in.name;
      expect_same_as_reference(g, in);
    }
  }
#else
  GTEST_SKIP() << "branch observations need the private hooks";
#endif
}

// A run of equal costs reuses the previous ordinal without a lookup: one cost repeated costs
// exactly one probe; two alternating costs in different slots cost one probe per edge.
TEST(GraphLayout, RunsOfEqualCostsSkipTheLookup) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
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
    EXPECT_EQ(detail::thread_construction_counts().last_probes, alternate ? 1000 : 1) << in.name;
    expect_same_as_reference(g, in);
  }
#else
  GTEST_SKIP() << "branch observations need the private hooks";
#endif
}

// The cumulative probe budget B = 8E + 8 x 65,536 (APPROACH 3 revision 43). Costs whose hash has
// 11 low zero bits start at slot 0 in every table up to 2,048 slots, so they form one cluster:
// the j-th of them (from 0) is found after j + 1 probes. After the 64 cluster costs come D
// duplicates alternating z[62], z[63] (a run of one cost reuses its ordinal without a lookup, so
// the duplicates must alternate), each costing 63 or 64 probes while B grows by 8 per edge. The
// boundaries come from an exact model of the probe accounting (build/1e/prep/budget_model.py):
//  * collection: exhausted at an input lookup exactly when D >= 9,406;
//  * growth: with a 65th cluster cost at the end (its lookup, then the reinsertion into 512
//    slots, then the lookup repeated), exhausted inside the growth for 9,366 <= D <= 9,404, not at
//    all for D <= 9,365, and at the 65th cost's first lookup (collection) for D = 9,405.
// Exhausted constructions use exactly B probes and none more; the arrays equal the reference.
TEST(GraphLayout, ProbeBudgetFallbackInCollectionAndGrowth) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  std::vector<Cost> z;
  for (Cost c = 1; z.size() < 65; ++c) {
    if ((detail::cost_hash(c) & 2047) == 0) z.push_back(c);
  }
  auto shape = [&](std::size_t dups, bool tail) {
    Input in {"cluster D=" + std::to_string(dups) + (tail ? " + 65th" : ""), 7};
    std::vector<Cost> c(z.begin(), z.begin() + 64);
    for (std::size_t i = 0; i < dups; ++i) c.push_back(z[i % 2 == 0 ? 62 : 63]);
    if (tail) c.push_back(z[64]);
    for (std::size_t i = 0; i < c.size(); ++i) {
      in.src.push_back(static_cast<NodeId>(i % 7));
      in.dst.push_back(static_cast<NodeId>((i * 3) % 7));
      in.cap.push_back(static_cast<Amount>(i % 11));
      in.cost.push_back(c[i]);
    }
    return in;
  };
  enum Where { None, Collection, Growth };
  const std::vector<std::tuple<std::size_t, bool, Where>> cases = {
      {9000, false, None},  {9405, false, None},  {9406, false, Collection}, {10000, false, Collection},
      {9300, true, None},   {9365, true, None},   {9366, true, Growth},       {9380, true, Growth},
      {9404, true, Growth}, {9405, true, Collection},
  };

  for (const auto& [dups, tail, where] : cases) {
    const auto in = shape(dups, tail);
    auto& counts = detail::thread_construction_counts();
    const auto before = counts;
    const Graph g = build(in);
    const std::int64_t budget = 8 * static_cast<std::int64_t>(in.src.size()) + 8 * 65536;
    EXPECT_EQ(counts.last_budget, budget) << in.name;
    EXPECT_LE(counts.last_probes, budget) << in.name;  // no lookup examines more than B in total
    EXPECT_EQ(counts.exhausted_collection - before.exhausted_collection, where == Collection ? 1 : 0) << in.name;
    EXPECT_EQ(counts.exhausted_growth - before.exhausted_growth, where == Growth ? 1 : 0) << in.name;
    EXPECT_EQ(counts.fallback - before.fallback, where == None ? 0 : 1) << in.name;
    EXPECT_EQ(counts.hashed - before.hashed, where == None ? 1 : 0) << in.name;
    if (where != None) EXPECT_EQ(counts.last_probes, budget) << in.name;
    expect_same_as_reference(g, in);
  }
#else
  GTEST_SKIP() << "branch observations need the private hooks";
#endif
}

// The budget stays cumulative across the large table's growths (review-1e): 64 costs share their
// low 18 hash bits, so they collide in every table up to 262,144 slots; further costs grow the
// table to K = 32,769; then the last two colliding costs alternate. An independent count (linear
// probing, every examined slot counted, a growth reinserting the distinct costs in first-appearance
// order) finds the first edge count whose probes exceed B: the collection exhausts the budget
// there, and one edge fewer stays hashed with exactly the counted probes.
TEST(GraphLayout, ProbeBudgetIsCumulativeAcrossTheLargeTableGrowths) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  std::vector<Cost> z;
  for (Cost c = 1; z.size() < 64; ++c) {
    if ((detail::cost_hash(c) & 262143) == 0) z.push_back(c);
  }
  std::vector<Cost> seq = z;
  for (Cost c = 1; seq.size() < 32769; ++c) {
    if (!std::binary_search(z.begin(), z.end(), c)) seq.push_back(c);
  }
  std::vector<Cost> table(64, 0);
  std::vector<Cost> order;
  std::int64_t probes = 0;
  auto locate = [&](Cost c) {
    for (std::size_t at = detail::cost_hash(c) & (table.size() - 1);; at = (at + 1) & (table.size() - 1)) {
      ++probes;
      if (table[at] == 0 || table[at] == c) return at;
    }
  };
  for (const Cost c : seq) {
    std::size_t at = locate(c);
    if (table[at] != 0) continue;
    if (4 * (order.size() + 1) > table.size()) {
      table.assign(2 * table.size(), 0);
      for (const Cost d : order) table[locate(d)] = d;
      at = locate(c);
    }
    table[at] = c;
    order.push_back(c);
  }
  ASSERT_EQ(order.size(), 32769U);
  ASSERT_EQ(table.size(), 262144U);
  std::int64_t prefix_probes = 0;  // all probes of the longest prefix within its budget
  for (std::size_t d = 0;; ++d) {
    ASSERT_LT(d, 100000U);
    const std::int64_t before = probes;
    seq.push_back(z[d % 2 == 0 ? 62 : 63]);
    locate(seq.back());
    if (probes > 8 * static_cast<std::int64_t>(seq.size()) + 8 * 65536) {
      prefix_probes = before;
      break;
    }
  }
  EXPECT_EQ(seq.size(), 44436U);  // the review's independent count (review-1e/probe.log)

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
  auto& counts = detail::thread_construction_counts();
  for (const std::size_t e : {seq.size() - 1, seq.size()}) {
    const bool exhausts = e == seq.size();
    const auto in = input(e);
    const auto before = counts;
    const Graph g = build(in);
    const std::int64_t budget = 8 * static_cast<std::int64_t>(e) + 8 * 65536;
    EXPECT_EQ(counts.last_budget, budget) << in.name;
    EXPECT_EQ(counts.last_probes, exhausts ? budget : prefix_probes) << in.name;
    EXPECT_EQ(counts.exhausted_collection - before.exhausted_collection, exhausts ? 1 : 0) << in.name;
    EXPECT_EQ(counts.exhausted_growth - before.exhausted_growth, 0) << in.name;
    EXPECT_EQ(counts.hashed - before.hashed, exhausts ? 0 : 1) << in.name;
    if (!exhausts) EXPECT_EQ(counts.last_distinct, 32769) << in.name;
    expect_same_as_reference(g, in);
  }
#else
  GTEST_SKIP() << "branch observations need the private hooks";
#endif
}

// manifest 10: exact units above 2^53 (capacities and costs kept as int64, not through double).
TEST(GraphLayout, ExactUnitsAboveTwo53) {
  const Amount big = (Amount{1} << 53) + 1;
  const Input in {"2^53 + 1", 2, {0, 0, 1}, {1, 1, 0}, {big, big + 2, 0}, {(Cost{1} << 53) + 1, (Cost{1} << 53) + 3, 1}};
  const Graph g = build(in);
  expect_same_as_reference(g, in);
  EXPECT_EQ(g.capacity_units()[0], big);
  EXPECT_EQ(g.cost()[0], (Cost{1} << 53) + 1);
}

// manifest 11: a copy keeps every array and the map and draws a new uid.
TEST(GraphLayout, CopyKeepsEveryArrayAndDrawsANewUid) {
  for (const auto& in : {exact_k(31, 500, 64, 1), exact_k(65537, 10, 512, 2), Input {"empty", 5}}) {
    const Graph g = build(in);
    const Graph c(g);
    EXPECT_NE(c.uid(), g.uid());
    EXPECT_GT(c.uid(), g.uid());  // a relation, not an absolute counter value
    expect_same_as_reference(c, in);
    EXPECT_EQ(GraphAccess::retained(c).size_bytes, GraphAccess::retained(g).size_bytes);
  }
}

// The move (P1, review-1e-consult): every value and the uid transferred without allocation
// (the fault executable's witness); the source keeps its resolution, has uid, N, E, maximum and
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
    EXPECT_EQ(x.max_cost(), 0) << what;
    EXPECT_EQ(x.cost_gcd(), 0) << what;
    EXPECT_TRUE(x.capacity_units().empty() && x.cost().empty() && x.edge_src().empty() && x.edge_dst().empty() &&
                x.input_index().empty() && x.in_adj_edge_index().empty() && GraphAccess::in_col_indices(x).empty())
        << what;
    EXPECT_EQ(vec(x.row_offsets()), (std::vector<std::int32_t> {0})) << what;
    EXPECT_EQ(vec(x.in_row_offsets()), (std::vector<std::int32_t> {0})) << what;
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

// The retained array payload: 36 E + 8 (N + 1) bytes in exactly nine arrays, capacities equal to
// sizes (no slack, no retained scratch, no inverse map, no physical-capacity array).
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
      EXPECT_EQ(r.arrays, 9u) << in.name;
      EXPECT_EQ(r.size_bytes, want) << in.name;
      EXPECT_EQ(r.capacity_bytes, want) << in.name;
    }
  }
}

// Manifest 12: graph-bound objects from another graph are rejected at their existing use sites,
// on a twin graph with identical arrays (only the uid differs, so the check is by identity, not by
// shape): the evaluators' Eligibility and DagCache checks, and a Path's DAG (to_dag keeps the
// Path's graph) in resolve_to_paths. The flow core, ledger, policy, KSP, search and cache sites
// have their own tests (NOTES section 7, item 3).
TEST(GraphLayout, ForeignGraphBoundObjectsAreRejectedAtTheirUseSites) {
  const std::vector<NodeId> src {0, 0, 1, 2}, dst {1, 2, 3, 3};
  const std::vector<Amount> cap {5, 5, 5, 5};
  const std::vector<Cost> cost {1, 1, 1, 1};
  const Graph g(4, src, dst, cap, cost, 1.0), twin(4, src, dst, cap, cost, 1.0);
  ASSERT_NE(g.uid(), twin.uid());
  const auto message = [](const std::function<void()>& f) {
    try {
      f();
    } catch (const std::invalid_argument& e) {
      return std::string(e.what());
    }
    return std::string("no std::invalid_argument");
  };
  const std::vector<NodeId> s {0}, d {3};
  const FlowPolicyConfig cfg {CostOnly {}, FlowPlacement::MaxFlow, Admission::Lossless,
                              PathSelection(EdgeTieBreak::Deterministic, EdgePriority::Capacity)};
  const Eligibility el = twin.eligibility({}, std::vector<EdgeId> {0});
  DagCache cache(twin);
  const std::string elig_msg = "foreign Eligibility: built for another graph";
  const std::string cache_msg = "foreign DagCache: bound to another graph";
  EXPECT_EQ(message([&] { (void)capacity(g, s, d, CapacityOptions {cfg, &el}); }), elig_msg);
  EXPECT_EQ(message([&] { (void)capacity(g, s, d, CapacityOptions {cfg, nullptr, &cache}); }), cache_msg);
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs {{s, d}};
  EXPECT_EQ(message([&] { (void)batch_capacity(g, pairs, CapacityOptions {cfg, nullptr, &cache}, {}, 1); }), cache_msg);
  EXPECT_EQ(message([&] { (void)sensitivity(g, s, d, SensitivityOptions {cfg, &el}, 1); }), elig_msg);
  EXPECT_EQ(message([&] { (void)sensitivity(g, s, d, SensitivityOptions {cfg, nullptr, &cache}, 1); }), cache_msg);
  const Path p(twin, {{0}, {2}});
  EXPECT_EQ(message([&] { (void)resolve_to_paths(g, p.to_dag(), 0, 3); }), "foreign PredDAG: built for another graph");
  // The same calls on the owning graph succeed (the controls).
  EXPECT_NO_THROW((void)capacity(twin, s, d, CapacityOptions {cfg, &el, &cache}));
  EXPECT_EQ(resolve_to_paths(twin, p.to_dag(), 0, 3).size(), 1U);
}
