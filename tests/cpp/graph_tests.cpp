// The graph on integer units: construction, validation order, the EdgeId order
// (source, cost, destination, input), the reverse CSR, input_positions, to_units and as_float.
#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

// A Graph of these arrays (resolution 1 unless given).
Graph build(std::int32_t n, std::vector<NodeId> src, std::vector<NodeId> dst, std::vector<Amount> cap,
            std::vector<Cost> cost, double q = 1.0) {
  return Graph(n, src, dst, cap, cost, q);
}

} // namespace

TEST(Graph, EdgeIdOrderInputIndexAndReverseCsr) {
  // (u, v, cap, cost) in input order. EdgeIds rank edges by source, then cost, then destination,
  // then input position.
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
  EXPECT_EQ(std::vector<std::int32_t>(g.input_positions().begin(), g.input_positions().end()), expect);
  EXPECT_EQ(expect, (std::vector<std::int32_t>{5, 2, 4, 1, 3, 6, 0}));
  EXPECT_EQ(std::vector<std::int32_t>(g.outgoing_offsets().begin(), g.outgoing_offsets().end()),
            (std::vector<std::int32_t>{0, 4, 5, 7}));
  // Every per-edge array is the input array permuted by input_positions.
  for (std::size_t e = 0; e < edges.size(); ++e) {
    const auto& in = edges[static_cast<std::size_t>(g.input_positions()[e])];
    EXPECT_EQ(g.edge_sources()[e], in.u);
    EXPECT_EQ(g.edge_destinations()[e], in.v);
    EXPECT_EQ(g.edge_capacity_units()[e], in.cap);
    EXPECT_EQ(g.edge_costs()[e], in.cost);
  }
  // Row u's id range holds u's out-edges (their (cost, destination, input) order within the row is the
  // input_positions checked above).
  for (NodeId u = 0; u < 3; ++u) {
    for (auto e = g.outgoing_offsets()[u]; e < g.outgoing_offsets()[u + 1]; ++e) {
      EXPECT_EQ(g.edge_sources()[static_cast<std::size_t>(e)], u);
    }
  }
  // The inverse map from input positions to EdgeIds, derived by hand from the order above.
  EXPECT_EQ(id_of_input(g), (std::vector<EdgeId> {6, 3, 1, 4, 2, 0, 5}));
  std::vector<Amount> caps_in;
  for (const auto& e : edges) caps_in.push_back(e.cap);
  EXPECT_EQ(to_input(g, g.edge_capacity_units()), caps_in);
  // Reverse CSR in (cost, source, input) order per destination, as input positions.
  const std::vector<std::vector<std::int32_t>> expect_rev = {{0}, {6, 2, 4}, {5, 3, 1}};
  for (NodeId v = 0; v < 3; ++v) {
    std::vector<std::int32_t> got;
    for (auto s = g.incoming_offsets()[v]; s < g.incoming_offsets()[v + 1]; ++s) {
      const EdgeId e = g.incoming_edges()[static_cast<std::size_t>(s)];
      EXPECT_EQ(g.edge_destinations()[static_cast<std::size_t>(e)], v);
      got.push_back(g.input_positions()[static_cast<std::size_t>(e)]);
    }
    EXPECT_EQ(got, expect_rev[static_cast<std::size_t>(v)]) << "node " << v;
  }
}

TEST(Graph, EqualKeysKeepInputOrder) {
  // Parallel edges with equal (source, cost, destination) are ranked by input position, in both CSRs.
  auto g = make_graph(3, {{0, 1, 4, 2}, {1, 2, 1, 1}, {0, 1, 7, 2}, {0, 1, 5, 1}, {2, 1, 3, 2}, {0, 1, 9, 2}});
  EXPECT_EQ(std::vector<std::int32_t>(g.input_positions().begin(), g.input_positions().end()),
            (std::vector<std::int32_t>{3, 0, 2, 5, 1, 4}));
  std::vector<std::int32_t> rev;
  for (auto s = g.incoming_offsets()[1]; s < g.incoming_offsets()[2]; ++s) {
    rev.push_back(g.input_positions()[static_cast<std::size_t>(g.incoming_edges()[static_cast<std::size_t>(s)])]);
  }
  // (cost, source, input): cost 1 from 0 (input 3), then cost 2 from 0 (inputs 0, 2, 5), then from 2.
  EXPECT_EQ(rev, (std::vector<std::int32_t>{3, 0, 2, 5, 4}));
}

TEST(Graph, ReverseOrderIsCostThenSource) {
  // A scan in id order would give (source, cost, input); the reverse order is (cost, source, input).
  auto g = make_graph(3, {{0, 2, 1, 5}, {1, 2, 1, 1}});
  ASSERT_EQ(g.incoming_offsets()[2], 0);
  const EdgeId first = g.incoming_edges()[0];
  EXPECT_EQ(g.input_positions()[static_cast<std::size_t>(first)], 1);
}

TEST(Graph, ShapesOfEmptyAndIsolated) {
  auto empty = build(0, {}, {}, {}, {});
  EXPECT_EQ(empty.num_nodes(), 0);
  EXPECT_EQ(empty.num_edges(), 0);
  EXPECT_EQ(std::vector<std::int32_t>(empty.outgoing_offsets().begin(), empty.outgoing_offsets().end()),
            (std::vector<std::int32_t>{0}));
  EXPECT_EQ(std::vector<std::int32_t>(empty.incoming_offsets().begin(), empty.incoming_offsets().end()),
            (std::vector<std::int32_t>{0}));
  EXPECT_EQ(empty.max_edge_cost(), 0);
  EXPECT_EQ(empty.edge_cost_gcd(), 0);
  auto isolated = build(4, {}, {}, {}, {});
  EXPECT_EQ(std::vector<std::int32_t>(isolated.outgoing_offsets().begin(), isolated.outgoing_offsets().end()),
            (std::vector<std::int32_t>{0, 0, 0, 0, 0}));
  auto g = make_graph(5, {{3, 1, 1, 6}, {1, 3, 1, 4}, {3, 3, 2, 10}});
  EXPECT_EQ(std::vector<std::int32_t>(g.outgoing_offsets().begin(), g.outgoing_offsets().end()),
            (std::vector<std::int32_t>{0, 0, 1, 1, 3, 3}));
  EXPECT_EQ(std::vector<std::int32_t>(g.incoming_offsets().begin(), g.incoming_offsets().end()),
            (std::vector<std::int32_t>{0, 0, 1, 1, 3, 3}));
  EXPECT_EQ(g.max_edge_cost(), 10);
  EXPECT_EQ(g.edge_cost_gcd(), 2);
}

// The cost loop after the running GCD reaches 1: every later edge is still checked for a
// positive cost (in input order: the first bad index is named) and still counts for the
// maximum; the sums are still validated; the GCD stays 1, and 0 without edges.
TEST(Graph, CostChecksAndMaximumAfterTheGcdReachesOne) {
  // The GCD is 1 from the second edge on; the maximum comes later.
  auto g = build(2, {0, 0, 0, 0}, {1, 1, 1, 1}, {1, 1, 1, 1}, {3, 2, 9, 4});
  EXPECT_EQ(g.edge_cost_gcd(), 1);
  EXPECT_EQ(g.max_edge_cost(), 9);
  auto late = build(2, {0, 0, 0}, {1, 1, 1}, {1, 1, 1}, {4, 6, 9});
  EXPECT_EQ(late.edge_cost_gcd(), 1);  // 1 only at the last edge
  EXPECT_EQ(late.max_edge_cost(), 9);
  auto ones = build(2, {0, 0, 0}, {1, 1, 1}, {1, 1, 1}, {1, 1, kL62 / 4});
  EXPECT_EQ(ones.edge_cost_gcd(), 1);
  EXPECT_EQ(ones.max_edge_cost(), kL62 / 4);
  EXPECT_EQ(build(2, {0}, {1}, {1}, {6}).edge_cost_gcd(), 6);
  EXPECT_EQ(build(2, {0, 0}, {1, 1}, {1, 1}, {6, 4}).edge_cost_gcd(), 2);
  EXPECT_EQ(build(0, {}, {}, {}, {}).edge_cost_gcd(), 0);
  // A non-positive cost after the GCD reached 1 is still found, at its own index.
  const auto zero = error_of([] { (void)build(2, {0, 0, 0}, {1, 1, 1}, {1, 1, 1}, {2, 3, 0}); });
  EXPECT_TRUE(contains(zero, "edge_costs[2] = 0")) << zero;
  const auto negative = error_of([] { (void)build(2, {0, 0, 0, 0}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 5, -7, 0}); });
  EXPECT_TRUE(contains(negative, "edge_costs[2] = -7")) << negative;
  // The cost sum after the GCD reached 1.
  EXPECT_TRUE(contains(error_of([] { (void)build(2, {0, 0, 0}, {1, 1, 1}, {1, 1, 1}, {1, kL62 / 2, kL62 / 2}); }),
                       "sum of all edge costs"));
}

TEST(Graph, Resolution) {
  // The resolution has no default: no constructor takes the arrays without it.
  static_assert(!std::is_constructible_v<Graph, std::int32_t, std::span<const NodeId>, std::span<const NodeId>,
                                         std::span<const Amount>, std::span<const Cost>>);
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
  EXPECT_EQ(error_of([] { (void)build(2, {0}, {1}, {-1}, {1}); }),
            "edge_capacity_units[0] = -1 is outside the units domain [0, 2^62)");
  EXPECT_EQ(error_of([] { (void)build(2, {0}, {1}, {kL62}, {1}); }),
            "edge_capacity_units[0] = 4611686018427387904 is outside the units domain [0, 2^62)");
  // The shapes.
  EXPECT_EQ(error_of([] { (void)build(-1, {}, {}, {}, {}); }), "num_nodes must be >= 0 (got -1)");
  EXPECT_EQ(error_of([] { (void)build(2, {0, 1}, {1}, {1}, {1}); }),
            "edge_sources, edge_destinations, edge_capacity_units and edge_costs must have the same length (got 2, 1, 1, 1)");
  EXPECT_EQ(error_of([] { (void)build(2, {0}, {1}, {1, 2}, {1}); }),
            "edge_sources, edge_destinations, edge_capacity_units and edge_costs must have the same length (got 1, 1, 2, 1)");
  EXPECT_EQ(error_of([] { (void)build(2, {0}, {1}, {1}, {}); }),
            "edge_sources, edge_destinations, edge_capacity_units and edge_costs must have the same length (got 1, 1, 1, 0)");
  // Zero and negative costs.
  EXPECT_EQ(error_of([] { (void)build(2, {0}, {1}, {1}, {0}); }), "edge_costs[0] = 0: costs must be positive (edge_costs >= 1)");
  EXPECT_EQ(error_of([] { (void)build(2, {0}, {1}, {1}, {-5}); }), "edge_costs[0] = -5: costs must be positive (edge_costs >= 1)");
  // The remaining bounds.
  EXPECT_EQ(error_of([] { (void)build(2, {2}, {1}, {1}, {1}); }), "edge 0: endpoint (2, 1) out of node range [0, 2)");
  EXPECT_EQ(error_of([] { (void)build(2, {0}, {-1}, {1}, {1}); }), "edge 0: endpoint (0, -1) out of node range [0, 2)");
  EXPECT_EQ(error_of([] { (void)build(2, {0, 0}, {1, 1}, {kL62 / 2, kL62 / 2}, {1, 1}); }),
            "the sum of all capacities reaches the bound of the units domain (2^62 units)");
  EXPECT_EQ(error_of([] { (void)build(2, {0, 0}, {1, 1}, {1, 1}, {kL62 / 2, kL62 / 2}); }),
            "the sum of all edge costs reaches the bound 2^62: path costs would overflow int64");
  // Just below the bounds is accepted.
  auto caps = build(2, {0, 0}, {1, 1}, {kL62 / 2, kL62 / 2 - 1}, {1, 1});
  EXPECT_EQ(caps.edge_capacity_units()[0] + caps.edge_capacity_units()[1], kL62 - 1);
  auto costs = build(2, {0, 0}, {1, 1}, {1, 1}, {kL62 / 2, kL62 / 2 - 1});
  EXPECT_EQ(costs.max_edge_cost(), kL62 / 2);
  auto huge_cost = build(2, {0}, {1}, {1}, {kL62 - 1});
  EXPECT_EQ(huge_cost.max_edge_cost(), kL62 - 1);
  EXPECT_EQ(error_of([] { (void)build(2, {0}, {1}, {1}, {std::numeric_limits<Cost>::max()}); }),
            "the sum of all edge costs reaches the bound 2^62: path costs would overflow int64");
}

TEST(Graph, ValidationOrder) {
  // Resolution before the unit values.
  EXPECT_EQ(error_of([] { (void)build(2, {0}, {1}, {-1}, {1}, 3.0); }),
            "resolution must be a power of two in [2^-32, 2^32]");
  // The unit values before the shapes (the length message names edge_capacity_units too).
  EXPECT_EQ(error_of([] { (void)build(2, {0, 1}, {1}, {-1}, {1}); }),
            "edge_capacity_units[0] = -1 is outside the units domain [0, 2^62)");
  // The shapes before the costs.
  EXPECT_EQ(error_of([] { (void)build(2, {0, 1}, {1}, {1}, {0}); }),
            "edge_sources, edge_destinations, edge_capacity_units and edge_costs must have the same length (got 2, 1, 1, 1)");
  EXPECT_EQ(error_of([] { (void)build(-1, {0}, {1}, {1}, {0}); }), "num_nodes must be >= 0 (got -1)");
  // The costs before the remaining bounds (endpoints, sums).
  EXPECT_EQ(error_of([] { (void)build(2, {0}, {5}, {1}, {0}); }), "edge_costs[0] = 0: costs must be positive (edge_costs >= 1)");
  EXPECT_EQ(error_of([] { (void)build(2, {0, 0}, {1, 1}, {kL62 / 2, kL62 / 2}, {1, 0}); }),
            "edge_costs[1] = 0: costs must be positive (edge_costs >= 1)");
  // Endpoints before the sums.
  EXPECT_EQ(error_of([] { (void)build(2, {0, 7}, {1, 1}, {kL62 / 2, kL62 / 2}, {1, 1}); }),
            "edge 1: endpoint (7, 1) out of node range [0, 2)");
  // Both running sums reach 2^62 at the same edge: the capacity sum is reported.
  EXPECT_EQ(error_of([] { (void)build(2, {0, 0}, {1, 1}, {kL62 / 2, kL62 / 2}, {kL62 / 2, kL62 / 2}); }),
            "the sum of all capacities reaches the bound of the units domain (2^62 units)");
}

// Uids come from one process-wide counter and are never reused: a graph built after another was
// destroyed (in the same storage) gets a larger uid, and so does a copy.
TEST(Graph, UidsAreNeverReused) {
  std::uint64_t last = 0;
  for (int i = 0; i < 3; ++i) {
    const auto g = build(2, {0}, {1}, {1}, {1});
    EXPECT_GT(g.uid(), last) << i;
    last = g.uid();
  }
  const auto g = build(2, {0}, {1}, {1}, {1});
  const Graph copy(g);
  EXPECT_GT(copy.uid(), g.uid());
  EXPECT_GT(g.uid(), last);
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
  // The export above 2^53 is the nearest double (the documented export limit).
  EXPECT_EQ(as_float((Amount{1} << 53) + 1, 1.0), 9007199254740992.0);
}

TEST(Units, RejectedValues) {
  const double q = std::ldexp(1.0, -12);
  const std::string nan = "to_units: NaN has no unit value";
  const std::string infinite = "to_units: an infinite amount is not finite and has no unit value";
  const std::string bound = "to_units: the amount is at or above the bound of the units domain (2^62 units)";
  const std::string resolution = "resolution must be a power of two in [2^-32, 2^32]";
  EXPECT_EQ(error_of([&] { (void)to_units(std::numeric_limits<double>::quiet_NaN(), q); }), nan);
  EXPECT_EQ(error_of([&] { (void)to_units(std::numeric_limits<double>::infinity(), q); }), infinite);
  EXPECT_EQ(error_of([&] { (void)to_units(-std::numeric_limits<double>::infinity(), q); }), infinite);
  EXPECT_EQ(error_of([&] { (void)to_units(-1e-300, q); }), "to_units: a negative amount has no unit value");
  EXPECT_EQ(error_of([&] { (void)to_units(4611686018427387904.0, 1.0); }), bound);  // 2^62
  EXPECT_EQ(error_of([&] { (void)to_units(std::numeric_limits<double>::max(), std::ldexp(1.0, 32)); }), bound);
  EXPECT_EQ(error_of([&] { (void)to_units(1e300, std::ldexp(1.0, -32)); }), bound);
  // Both resolution checks: a value that is not a power of two, and a power of two out of range
  // (and the values that are not positive and finite).
  for (const double bad : {3.0, std::ldexp(1.0, -33), std::ldexp(1.0, 33), 0.0, -1.0,
                           std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
    EXPECT_EQ(error_of([&] { (void)to_units(1.0, bad); }), resolution) << bad;
    EXPECT_EQ(error_of([&] { (void)as_float(1, bad); }), resolution) << bad;
  }
  // The capacity-sum bound at q = 2^-12: [1e15, 1.25e14] fits, [1e15, 1.26e14] exceeds the unit sum.
  const std::vector<Amount> ok = {to_units(1e15, q), to_units(1.25e14, q)};
  auto g = build(2, {0, 0}, {1, 1}, ok, {1, 1}, q);
  EXPECT_EQ(g.edge_capacity_units()[0] + g.edge_capacity_units()[1], ok[0] + ok[1]);
  const std::vector<Amount> over = {to_units(1e15, q), to_units(1.26e14, q)};
  EXPECT_EQ(error_of([&] { (void)build(2, {0, 0}, {1, 1}, over, {1, 1}, q); }),
            "the sum of all capacities reaches the bound of the units domain (2^62 units)");
}
