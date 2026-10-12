// The graph inputs of the three graph test files (graph_layout_tests.cpp, graph_fault_tests.cpp,
// graph_concurrency_tests.cpp): one input in input coordinates, its graph and its view for the
// layout checker, its count of distinct costs, and its generators. Each generator
// is deterministic (seeded or closed-form), so a test file and its readers see the same inputs on
// every run. random_input draws general inputs; each other generator is built for the construction
// path named at it (hash table, growths, collisions, the K-limit and probe-budget fallbacks). The
// branch counters of graph_layout_tests.cpp (hooks builds) are read on exact_k's boundary inputs,
// budget_shape's clusters, a colliding_costs sequence that grows the table, and random_input's
// inputs (whether the fused path ran).
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "layout_checker.hpp"
#include "netgraph/core/detail/graph_access.hpp"
#include "netgraph/core/graph.hpp"

namespace netgraph::core::test {

// A graph input as the constructor takes it.
struct Input {
  std::string name;
  std::int32_t n {0};
  std::vector<NodeId> src {}, dst {};
  std::vector<Amount> cap {};
  std::vector<Cost> cost {};
  double q {1.0};
};

// The graph of an input, at its resolution.
inline Graph build(const Input& in) { return Graph(in.n, in.src, in.dst, in.cap, in.cost, in.q); }

// An input and a graph's stored layout as the layout checker reads them.
inline LayoutInput view(const Input& in) { return {in.n, in.src, in.dst, in.cap, in.cost}; }
inline LayoutSpans spans_of(const Graph& g) {
  return {g.num_nodes(), g.num_edges(), g.capacity_units(), g.cost(), g.edge_src(), g.edge_dst(), g.input_index(),
          g.row_offsets(), g.in_row_offsets(), g.in_adj_edge_index(), detail::GraphAccess::in_col_indices(g),
          g.max_cost(), g.cost_gcd()};
}

// The number of distinct costs of an input (K), which decides its construction path.
inline std::int64_t distinct(std::vector<Cost> c) {
  std::sort(c.begin(), c.end());
  return std::unique(c.begin(), c.end()) - c.begin();
}

// A copy of every stored array of a graph and of its two cost values: the before/after comparison
// of graph_concurrency_tests.cpp, and the copy that the layout checker's corruption controls edit
// (graph_layout_tests.cpp).
struct GraphArrays {
  std::vector<Amount> capacity;
  std::vector<Cost> cost;
  std::vector<NodeId> src, dst;
  std::vector<std::int32_t> input_index, row_offsets, in_row_offsets;
  std::vector<EdgeId> in_adj;
  std::vector<NodeId> in_col;
  Cost max_cost {0}, cost_gcd {0};
  bool operator==(const GraphArrays&) const = default;
};

inline GraphArrays arrays_of(const Graph& g) {
  const auto col = detail::GraphAccess::in_col_indices(g);
  return {{g.capacity_units().begin(), g.capacity_units().end()},
          {g.cost().begin(), g.cost().end()},
          {g.edge_src().begin(), g.edge_src().end()},
          {g.edge_dst().begin(), g.edge_dst().end()},
          {g.input_index().begin(), g.input_index().end()},
          {g.row_offsets().begin(), g.row_offsets().end()},
          {g.in_row_offsets().begin(), g.in_row_offsets().end()},
          {g.in_adj_edge_index().begin(), g.in_adj_edge_index().end()},
          {col.begin(), col.end()},
          g.max_cost(),
          g.cost_gcd()};
}

// m random edges over n nodes, costs 1..max_cost, capacities cap_low .. cap_low + cap_span - 1.
inline Input random_input(std::string name, std::int32_t n, std::size_t m, Cost max_cost, std::uint64_t seed,
                          Amount cap_low = 0, Amount cap_span = 5) {
  std::mt19937_64 rng(seed);
  Input in {std::move(name), n};
  for (std::size_t i = 0; i < m; ++i) {
    in.src.push_back(static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n)));
    in.dst.push_back(static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n)));
    in.cap.push_back(cap_low + static_cast<Amount>(rng() % static_cast<std::uint64_t>(cap_span)));
    in.cost.push_back(static_cast<Cost>(1 + rng() % static_cast<std::uint64_t>(max_cost)));
  }
  return in;
}

// Exactly k distinct costs (3 c + 1), each present, plus `extra` repeats, shuffled by the seed;
// distinct capacities 1, 2, ...
inline Input exact_k(std::int64_t k, std::size_t extra, std::int32_t n, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  Input in {"K=" + std::to_string(k), n};
  for (std::int64_t c = 1; c <= k; ++c) in.cost.push_back(3 * c + 1);
  for (std::size_t j = 0; j < extra; ++j) in.cost.push_back(in.cost[rng() % in.cost.size()]);
  std::shuffle(in.cost.begin(), in.cost.end(), rng);
  for (std::size_t i = 0; i < in.cost.size(); ++i) {
    in.src.push_back(static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n)));
    in.dst.push_back(static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n)));
    in.cap.push_back(static_cast<Amount>(1 + i));
  }
  return in;
}

// Exactly k distinct costs (c * 5 + 2) in ascending order of first appearance, plus `extra`
// repeats, by closed-form strides without a generator; capacities 0 to 8 (zero included).
inline Input exact_k_strided(std::int64_t k, std::size_t extra, std::int32_t n) {
  Input in {"K=" + std::to_string(k), n};
  for (std::int64_t c = 1; c <= k; ++c) in.cost.push_back(c * 5 + 2);
  for (std::size_t j = 0; j < extra; ++j) in.cost.push_back(in.cost[(j * 7919) % in.cost.size()]);
  for (std::size_t i = 0; i < in.cost.size(); ++i) {
    in.src.push_back(static_cast<NodeId>((i * 31) % static_cast<std::size_t>(n)));
    in.dst.push_back(static_cast<NodeId>((i * 17 + 3) % static_cast<std::size_t>(n)));
    in.cap.push_back(static_cast<Amount>(i % 9));
  }
  return in;
}

// The first `count` costs whose construction hash has `residue` under `mask` (they share a table
// slot under that mask).
inline std::vector<Cost> colliding_costs(std::uint64_t mask, std::size_t count, std::uint64_t residue = 5) {
  std::vector<Cost> out;
  for (Cost c = 1; out.size() < count; ++c) {
    if ((detail::cost_hash(c) & mask) == residue) out.push_back(c);
  }
  return out;
}

// `count` edges over 8 nodes whose costs collide in the initial 64-slot table.
inline Input colliding(std::size_t count) {
  Input in {"collisions", 8};
  in.cost = colliding_costs(63, count, 9);
  for (std::size_t i = 0; i < in.cost.size(); ++i) {
    in.src.push_back(static_cast<NodeId>(i % 8));
    in.dst.push_back(static_cast<NodeId>((i + 3) % 8));
    in.cap.push_back(1);
  }
  return in;
}

// The probe-budget shapes (GraphLayout.ProbeBudgetFallbackInCollectionAndGrowth): the 64 first
// costs whose hash has 11 low zero bits (one cluster at slot 0 in every table up to 2,048 slots),
// then `dups` duplicates alternating the 63rd and 64th, optionally a 65th cluster cost at the end.
inline Input budget_shape(const std::string& name, std::size_t dups, bool tail) {
  const std::vector<Cost> z = colliding_costs(2047, 65, 0);
  Input in {name, 7};
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
}

} // namespace netgraph::core::test
