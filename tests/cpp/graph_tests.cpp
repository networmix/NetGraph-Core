// The graph on integer units (APPROACH 2.1, 3): construction, validation order, layout B'
// ids, the reverse CSR, input_index, to_units and as_float.
#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

constexpr Amount L62 = Amount{1} << 62;

// The message of the std::invalid_argument thrown by `fn`, or "" when nothing is thrown.
template <typename F>
std::string error_of(F&& fn) {
  try {
    fn();
  } catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
}

Graph build(std::int32_t n, std::vector<NodeId> src, std::vector<NodeId> dst, std::vector<Amount> cap,
            std::vector<Cost> cost, double q = 1.0) {
  return Graph(n, src, dst, cap, cost, q);
}

bool contains(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

} // namespace

TEST(Graph, BPrimeOrderInputIndexAndReverseCsr) {
  // (u, v, cap, cost) in input order; prototype fixture graph_bprime_order.
  const std::vector<InEdge> edges = {{2, 0, 5, 3}, {0, 2, 4, 2}, {0, 1, 3, 2}, {1, 2, 2, 1},
                                     {0, 1, 1, 2}, {0, 2, 6, 1}, {2, 1, 7, 1}};
  auto g = make_graph(3, edges);
  std::vector<std::int32_t> expect(edges.size());
  std::iota(expect.begin(), expect.end(), 0);
  std::stable_sort(expect.begin(), expect.end(), [&](int a, int b) {
    const auto& x = edges[static_cast<std::size_t>(a)];
    const auto& y = edges[static_cast<std::size_t>(b)];
    return std::tie(x.u, x.cost, x.v) < std::tie(y.u, y.cost, y.v);
  });
  EXPECT_EQ(std::vector<std::int32_t>(g.input_index().begin(), g.input_index().end()), expect);
  EXPECT_EQ(expect, (std::vector<std::int32_t>{5, 2, 4, 1, 3, 6, 0}));
  EXPECT_EQ(std::vector<std::int32_t>(g.row_offsets().begin(), g.row_offsets().end()),
            (std::vector<std::int32_t>{0, 4, 5, 7}));
  // Every per-edge array is the input array permuted by input_index.
  for (std::size_t e = 0; e < edges.size(); ++e) {
    const auto& in = edges[static_cast<std::size_t>(g.input_index()[e])];
    EXPECT_EQ(g.edge_src()[e], in.u);
    EXPECT_EQ(g.edge_dst()[e], in.v);
    EXPECT_EQ(g.capacity_units()[e], in.cap);
    EXPECT_EQ(g.cost()[e], in.cost);
  }
  // The out-edges of u are the id range of its row, in (cost, dst, input) order.
  for (NodeId u = 0; u < 3; ++u) {
    for (auto e = g.row_offsets()[u]; e < g.row_offsets()[u + 1]; ++e) {
      EXPECT_EQ(g.edge_src()[static_cast<std::size_t>(e)], u);
    }
  }
  // input_index round trip through the caller's inverse.
  const auto inv = id_of_input(g);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    EXPECT_EQ(g.input_index()[static_cast<std::size_t>(inv[i])], static_cast<std::int32_t>(i));
  }
  std::vector<Amount> caps_in;
  for (const auto& e : edges) caps_in.push_back(e.cap);
  EXPECT_EQ(to_input(g, g.capacity_units()), caps_in);
  // Reverse CSR in (cost, src, input) order per destination, as input positions.
  const std::vector<std::vector<std::int32_t>> expect_rev = {{0}, {6, 2, 4}, {5, 3, 1}};
  for (NodeId v = 0; v < 3; ++v) {
    std::vector<std::int32_t> got;
    for (auto s = g.in_row_offsets()[v]; s < g.in_row_offsets()[v + 1]; ++s) {
      const EdgeId e = g.in_adj_edge_index()[static_cast<std::size_t>(s)];
      EXPECT_EQ(g.edge_dst()[static_cast<std::size_t>(e)], v);
      got.push_back(g.input_index()[static_cast<std::size_t>(e)]);
    }
    EXPECT_EQ(got, expect_rev[static_cast<std::size_t>(v)]) << "node " << v;
  }
}

TEST(Graph, EqualKeysKeepInputOrder) {
  // Parallel edges with equal (src, cost, dst) are ranked by input position, in both CSRs.
  auto g = make_graph(3, {{0, 1, 4, 2}, {1, 2, 1, 1}, {0, 1, 7, 2}, {0, 1, 5, 1}, {2, 1, 3, 2}, {0, 1, 9, 2}});
  EXPECT_EQ(std::vector<std::int32_t>(g.input_index().begin(), g.input_index().end()),
            (std::vector<std::int32_t>{3, 0, 2, 5, 1, 4}));
  std::vector<std::int32_t> rev;
  for (auto s = g.in_row_offsets()[1]; s < g.in_row_offsets()[2]; ++s) {
    rev.push_back(g.input_index()[static_cast<std::size_t>(g.in_adj_edge_index()[static_cast<std::size_t>(s)])]);
  }
  // (cost, src, input): cost 1 from 0 (input 3), then cost 2 from 0 (inputs 0, 2, 5), then from 2.
  EXPECT_EQ(rev, (std::vector<std::int32_t>{3, 0, 2, 5, 4}));
}

TEST(Graph, ReverseOrderIsCostThenSource) {
  // A scan in id order would give (src, cost, input); the contract is (cost, src, input).
  auto g = make_graph(3, {{0, 2, 1, 5}, {1, 2, 1, 1}});
  ASSERT_EQ(g.in_row_offsets()[2], 0);
  const EdgeId first = g.in_adj_edge_index()[0];
  EXPECT_EQ(g.input_index()[static_cast<std::size_t>(first)], 1);
}

TEST(Graph, ShapesOfEmptyAndIsolated) {
  auto empty = build(0, {}, {}, {}, {});
  EXPECT_EQ(empty.num_nodes(), 0);
  EXPECT_EQ(empty.num_edges(), 0);
  EXPECT_EQ(std::vector<std::int32_t>(empty.row_offsets().begin(), empty.row_offsets().end()),
            (std::vector<std::int32_t>{0}));
  EXPECT_EQ(std::vector<std::int32_t>(empty.in_row_offsets().begin(), empty.in_row_offsets().end()),
            (std::vector<std::int32_t>{0}));
  EXPECT_EQ(empty.max_cost(), 0);
  EXPECT_EQ(empty.cost_gcd(), 0);
  auto isolated = build(4, {}, {}, {}, {});
  EXPECT_EQ(std::vector<std::int32_t>(isolated.row_offsets().begin(), isolated.row_offsets().end()),
            (std::vector<std::int32_t>{0, 0, 0, 0, 0}));
  auto g = make_graph(5, {{3, 1, 1, 6}, {1, 3, 1, 4}, {3, 3, 2, 10}});
  EXPECT_EQ(std::vector<std::int32_t>(g.row_offsets().begin(), g.row_offsets().end()),
            (std::vector<std::int32_t>{0, 0, 1, 1, 3, 3}));
  EXPECT_EQ(std::vector<std::int32_t>(g.in_row_offsets().begin(), g.in_row_offsets().end()),
            (std::vector<std::int32_t>{0, 0, 1, 1, 3, 3}));
  EXPECT_EQ(g.max_cost(), 10);
  EXPECT_EQ(g.cost_gcd(), 2);
}

// The cost loop after the running GCD reaches 1: every later edge is still checked for a
// positive cost (in input order: the first bad index is named) and still counts for the
// maximum; the sums are still validated; the GCD stays 1, and 0 without edges.
TEST(Graph, CostChecksAndMaximumAfterTheGcdReachesOne) {
  // The GCD is 1 from the second edge on; the maximum comes later.
  auto g = build(2, {0, 0, 0, 0}, {1, 1, 1, 1}, {1, 1, 1, 1}, {3, 2, 9, 4});
  EXPECT_EQ(g.cost_gcd(), 1);
  EXPECT_EQ(g.max_cost(), 9);
  auto late = build(2, {0, 0, 0}, {1, 1, 1}, {1, 1, 1}, {4, 6, 9});
  EXPECT_EQ(late.cost_gcd(), 1);  // 1 only at the last edge
  EXPECT_EQ(late.max_cost(), 9);
  auto ones = build(2, {0, 0, 0}, {1, 1, 1}, {1, 1, 1}, {1, 1, L62 / 4});
  EXPECT_EQ(ones.cost_gcd(), 1);
  EXPECT_EQ(ones.max_cost(), L62 / 4);
  EXPECT_EQ(build(2, {0}, {1}, {1}, {6}).cost_gcd(), 6);
  EXPECT_EQ(build(2, {0, 0}, {1, 1}, {1, 1}, {6, 4}).cost_gcd(), 2);
  EXPECT_EQ(build(0, {}, {}, {}, {}).cost_gcd(), 0);
  // A non-positive cost after the GCD reached 1 is still found, at its own index.
  const auto zero = error_of([] { (void)build(2, {0, 0, 0}, {1, 1, 1}, {1, 1, 1}, {2, 3, 0}); });
  EXPECT_TRUE(contains(zero, "cost[2] = 0")) << zero;
  const auto negative = error_of([] { (void)build(2, {0, 0, 0, 0}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 5, -7, 0}); });
  EXPECT_TRUE(contains(negative, "cost[2] = -7")) << negative;
  // The cost sum after the GCD reached 1.
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0, 0, 0}, {1, 1, 1}, {1, 1, 1}, {1, L62 / 2, L62 / 2}); }),
                       "sum of all edge costs"));
}

TEST(Graph, ExactUnitsAboveTwo53) {
  const Amount big = (Amount{1} << 53) + 1;
  auto g = make_graph(2, {{0, 1, big, 1}});
  EXPECT_EQ(g.capacity_units()[0], big);
  EXPECT_EQ(as_float(big, 1.0), 9007199254740992.0);  // the documented export limit
}

TEST(Graph, CopyIsANewGraph) {
  auto g = make_graph(3, {{0, 2, 4, 2}, {0, 1, 3, 1}});
  Graph c(g);
  EXPECT_NE(c.uid(), g.uid());
  EXPECT_EQ(std::vector<std::int32_t>(c.input_index().begin(), c.input_index().end()),
            std::vector<std::int32_t>(g.input_index().begin(), g.input_index().end()));
  EXPECT_EQ(std::vector<Amount>(c.capacity_units().begin(), c.capacity_units().end()),
            std::vector<Amount>(g.capacity_units().begin(), g.capacity_units().end()));
  auto h = make_graph(3, {{0, 2, 4, 2}, {0, 1, 3, 1}});
  EXPECT_NE(h.uid(), g.uid());
  EXPECT_GT(g.uid(), 0u);
}

TEST(Graph, Resolution) {
  for (double q : {std::ldexp(1.0, -32), std::ldexp(1.0, -12), 1.0, std::ldexp(1.0, 32)}) {
    auto g = build(2, {0}, {1}, {1}, {1}, q);
    EXPECT_EQ(g.resolution(), q);
  }
  for (double q : {3.0, std::ldexp(1.0, -33), std::ldexp(1.0, 33), 0.0, -1.0, -0.5,
                   std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), 1e308}) {
    EXPECT_TRUE(contains(error_of([&] { (void)build(2, {0}, {1}, {1}, {1}, q); }), "resolution")) << q;
  }
}

TEST(Graph, EveryValidationError) {
  // The unit values.
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0}, {1}, {-1}, {1}); }), "capacity_units"));
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0}, {1}, {L62}, {1}); }), "units domain"));
  // The shapes.
  EXPECT_TRUE(contains(error_of([] { (void)build(-1, {}, {}, {}, {}); }), "num_nodes"));
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0, 1}, {1}, {1}, {1}); }), "same length"));
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0}, {1}, {1, 2}, {1}); }), "same length"));
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0}, {1}, {1}, {}); }), "same length"));
  // Zero and negative costs (rule 28).
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0}, {1}, {1}, {0}); }), "cost"));
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0}, {1}, {1}, {-5}); }), "costs must be positive"));
  // The remaining bounds.
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {2}, {1}, {1}, {1}); }), "out of node range"));
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0}, {-1}, {1}, {1}); }), "out of node range"));
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0, 0}, {1, 1}, {L62 / 2, L62 / 2}, {1, 1}); }),
                       "sum of all capacities"));
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0, 0}, {1, 1}, {1, 1}, {L62 / 2, L62 / 2}); }),
                       "sum of all edge costs"));
  // Just below the bounds is accepted.
  auto caps = build(2, {0, 0}, {1, 1}, {L62 / 2, L62 / 2 - 1}, {1, 1});
  EXPECT_EQ(caps.capacity_units()[0] + caps.capacity_units()[1], L62 - 1);
  auto costs = build(2, {0, 0}, {1, 1}, {1, 1}, {L62 / 2, L62 / 2 - 1});
  EXPECT_EQ(costs.max_cost(), L62 / 2);
  auto huge_cost = build(2, {0}, {1}, {1}, {L62 - 1});
  EXPECT_EQ(huge_cost.max_cost(), L62 - 1);
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0}, {1}, {1}, {std::numeric_limits<Cost>::max()}); }),
                       "sum of all edge costs"));
}

TEST(Graph, ValidationOrder) {
  // Resolution before the unit values.
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0}, {1}, {-1}, {1}, 3.0); }), "resolution"));
  // The unit values before the shapes.
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0, 1}, {1}, {-1}, {1}); }), "capacity_units"));
  // The shapes before the costs.
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0, 1}, {1}, {1}, {0}); }), "same length"));
  EXPECT_TRUE(contains(error_of([] { (void)build(-1, {0}, {1}, {1}, {0}); }), "num_nodes"));
  // The costs before the remaining bounds (endpoints, sums).
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0}, {5}, {1}, {0}); }), "costs must be positive"));
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0, 0}, {1, 1}, {L62 / 2, L62 / 2}, {1, 0}); }),
                       "costs must be positive"));
  // Endpoints before the sums.
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0, 7}, {1, 1}, {L62 / 2, L62 / 2}, {1, 1}); }),
                       "out of node range"));
}

TEST(Units, HalfUpBoundaries) {
  for (int k : {-12, 0, -32, 32}) {
    const double q = std::ldexp(1.0, k);
    EXPECT_EQ(to_units(2.5 * q, q), 3) << k;
    EXPECT_EQ(to_units(0.5 * q, q), 1) << k;
    EXPECT_EQ(to_units(0.49 * q, q), 0) << k;
    EXPECT_EQ(to_units(7.5 * q, q), 8) << k;
    EXPECT_EQ(to_units(std::nextafter(7.5 * q, 0.0), q), 7) << k;
    EXPECT_EQ(to_units(std::nextafter(0.5 * q, 0.0), q), 0) << k;
    EXPECT_EQ(to_units(0.0, q), 0) << k;
    EXPECT_EQ(to_units(-0.0, q), 0) << k;
    EXPECT_EQ(as_float(12345, q), 12345 * q) << k;
  }
  // Physical 1e15 at q = 2^-12 is 4.096e18 units: below 2^62 and exact.
  EXPECT_EQ(to_units(1e15, std::ldexp(1.0, -12)), Amount{4096000000000000000});
  // Large integral values: the double is the unit count.
  EXPECT_EQ(to_units(9007199254740993.0, 1.0), Amount{9007199254740992});
  EXPECT_EQ(to_units(4503599627370497.5, 1.0), Amount{4503599627370498});  // 2^52 + 1.5 is not a double
  EXPECT_EQ(to_units(4503599627370495.5, 1.0), Amount{4503599627370496});  // 2^52 - 0.5 rounds up
  EXPECT_EQ(to_units(4611686018427387392.0, 1.0), Amount{4611686018427387392});  // 2^62 - 512
}

TEST(Units, RejectedValues) {
  const double q = std::ldexp(1.0, -12);
  EXPECT_TRUE(contains(error_of([&] { (void)to_units(std::numeric_limits<double>::quiet_NaN(), q); }), "NaN"));
  EXPECT_TRUE(contains(error_of([&] { (void)to_units(std::numeric_limits<double>::infinity(), q); }), "finite"));
  EXPECT_TRUE(contains(error_of([&] { (void)to_units(-std::numeric_limits<double>::infinity(), q); }), "finite"));
  EXPECT_TRUE(contains(error_of([&] { (void)to_units(-1e-300, q); }), "negative"));
  EXPECT_TRUE(contains(error_of([&] { (void)to_units(4611686018427387904.0, 1.0); }), "bound"));  // 2^62
  EXPECT_TRUE(contains(error_of([&] { (void)to_units(std::numeric_limits<double>::max(), std::ldexp(1.0, 32)); }),
                       "bound"));
  EXPECT_TRUE(contains(error_of([&] { (void)to_units(1e300, std::ldexp(1.0, -32)); }), "bound"));
  EXPECT_TRUE(contains(error_of([&] { (void)to_units(1.0, 3.0); }), "resolution"));
  EXPECT_TRUE(contains(error_of([&] { (void)as_float(1, 0.0); }), "resolution"));
  // Admission at q = 2^-12: [1e15, 1.25e14] fits, [1e15, 1.26e14] exceeds the unit sum.
  const std::vector<Amount> ok = {to_units(1e15, q), to_units(1.25e14, q)};
  auto g = build(2, {0, 0}, {1, 1}, ok, {1, 1}, q);
  EXPECT_EQ(g.capacity_units()[0] + g.capacity_units()[1], ok[0] + ok[1]);
  const std::vector<Amount> over = {to_units(1e15, q), to_units(1.26e14, q)};
  EXPECT_TRUE(contains(error_of([&] { (void)build(2, {0, 0}, {1, 1}, over, {1, 1}, q); }), "sum of all capacities"));
}
