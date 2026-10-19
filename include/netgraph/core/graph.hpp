/* An immutable directed graph that allows several edges between the same ordered node pair.
   Capacities use integer units. */
#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

namespace detail {
struct GraphAccess;
} // namespace detail

// Convert a physical value to units of resolution q (the physical size of one unit).
// Round to the nearest integer, with exact halves rounded up. This is the only way a physical
// value enters Core. Throw std::invalid_argument for NaN, infinity, a negative value, a result
// >= 2^62, or a resolution that is not a power of two in [2^-32, 2^32].
[[nodiscard]] Amount to_units(double x, double resolution);
// Units to a physical double (x * q): exact below 2^53 units.
[[nodiscard]] double as_float(Amount x, double resolution);

// Edges arrive in input order. Their public EdgeIds are their positions after sorting by
// (source, cost, destination, input). Node u's out-edges are
// [outgoing_offsets()[u], outgoing_offsets()[u + 1]), ordered by (cost, destination, input);
// input_positions()[e] gives e's input position. All other per-edge arrays use EdgeId order, except incoming_edges() (reverse-CSR slot order). CSR (compressed
// sparse row) stores each node's edges in a flat array and uses offsets to mark its range.
// Reverse CSR lists v's in-edges in (cost, source, input) order.
//
// The constructor validates in this order and reports the first defect (std::invalid_argument):
// the resolution (a power of two in [2^-32, 2^32]); the unit values (every capacity in
// [0, 2^62)); the array shapes (num_nodes >= 0, equal lengths, at most INT32_MAX edges); zero
// or negative costs (costs are positive); the remaining bounds (endpoints in
// [0, num_nodes), then the capacity and cost sums below 2^62, so no path cost overflows int64:
// both are added up edge by edge in input order, and the first edge at which a running sum
// reaches 2^62 reports that sum, the capacity sum before the cost sum at the same edge). The
// order is part of the interface: a language binding validates in the same order before
// narrowing any value, so the reported defect does not depend on the representation of the
// input.
//
// Construction uses N nodes, E edges and K distinct costs. A cost's rank is its position among
// the sorted distinct costs. Stable counting passes group by each key without changing ties:
// destination, cost rank, source give EdgeId order; cost rank, destination give reverse CSR.
// When N K < E, use destination then (source, rank) forward, and (destination, rank) reverse.
// With no edges, only validation and two
// zero-filled offset arrays are needed. Rank costs through a hash table up to 65,536 distinct
// costs and 8E + 8 x 65,536 examined slots; otherwise sort all costs. Both give the same ranks.
// Expected O(E + N) for bounded K; O(N + E log(E + 1)) worst case.
// Retained arrays: 36 bytes per edge plus 8 (N + 1) bytes of offsets.
//
// A copy is a new graph: same ids and input map, a new uid. Uids come from one process-wide
// 64-bit counter and are never reused. A move transfers every array and the uid without
// allocation. The source retains its resolution, has uid, node count, edge count, maximum cost and
// cost GCD zero, exposes empty edge arrays, and exposes one zero in each offset view. Its empty
// offset vectors use an immutable static zero for those views.
class Graph {
public:
  Graph(std::int32_t num_nodes,
        std::span<const NodeId> edge_sources,
        std::span<const NodeId> edge_destinations,
        std::span<const Amount> edge_capacity_units,
        std::span<const Cost> edge_costs,
        double resolution);
  Graph(const Graph& other);
  Graph(Graph&& other) noexcept;
  Graph& operator=(const Graph&) = delete;
  Graph& operator=(Graph&&) = delete;
  ~Graph() noexcept = default;

  [[nodiscard]] std::uint64_t uid() const noexcept { return uid_; }
  [[nodiscard]] double resolution() const noexcept { return resolution_; }
  [[nodiscard]] std::int32_t num_nodes() const noexcept { return num_nodes_; }
  [[nodiscard]] std::int32_t num_edges() const noexcept {
    return static_cast<std::int32_t>(edge_capacity_units_.size());
  }
  [[nodiscard]] std::span<const Amount> edge_capacity_units() const noexcept { return edge_capacity_units_; }
  [[nodiscard]] std::span<const Cost> edge_costs() const noexcept { return edge_costs_; }
  [[nodiscard]] std::span<const NodeId> edge_sources() const noexcept { return edge_sources_; }
  [[nodiscard]] std::span<const NodeId> edge_destinations() const noexcept { return edge_destinations_; }
  // Maximum edge cost and greatest common divisor of all edge costs, computed by full scans
  // at construction; both are 0 on a graph without edges.
  [[nodiscard]] Cost max_edge_cost() const noexcept { return max_edge_cost_; }
  [[nodiscard]] Cost edge_cost_gcd() const noexcept { return edge_cost_gcd_; }
  // id -> input position. The inverse permutation is the caller's (one O(E) scatter).
  [[nodiscard]] std::span<const std::int32_t> input_positions() const noexcept { return input_positions_; }
  // Forward CSR offsets (N + 1 entries, [0] = 0, [N] = E): the id ranges of the out-edges.
  [[nodiscard]] std::span<const std::int32_t> outgoing_offsets() const noexcept { return offsets(outgoing_offsets_); }
  // Reverse CSR: the in-edges of v are incoming_edges()[incoming_offsets()[v] ..
  // incoming_offsets()[v + 1]) (public EdgeIds); the reverse column is
  // edge_sources()[incoming_edges()[slot]].
  [[nodiscard]] std::span<const std::int32_t> incoming_offsets() const noexcept { return offsets(incoming_offsets_); }
  [[nodiscard]] std::span<const EdgeId> incoming_edges() const noexcept { return incoming_edges_; }

  // Static eligibility from excluded ids (duplicates and unsorted lists accepted, out-of-range
  // ids rejected). Empty lists mean all allowed (uid 0). The one public construction of an
  // Eligibility from excluded ids (Eligibility::intersect combines two values).
  [[nodiscard]] Eligibility eligibility(std::span<const NodeId> excluded_nodes = {},
                                        std::span<const EdgeId> excluded_edges = {}) const;

private:
  friend struct detail::GraphAccess;

  // The offsets of a moved-from graph (its vectors are empty): one zero.
  static constexpr std::int32_t kNoOffsets[1] {0};
  [[nodiscard]] static std::span<const std::int32_t> offsets(const std::vector<std::int32_t>& v) noexcept {
    return v.empty() ? std::span<const std::int32_t>(kNoOffsets) : std::span<const std::int32_t>(v);
  }

  std::uint64_t uid_ {0};
  double resolution_ {1.0};
  std::int32_t num_nodes_ {0};
  Cost max_edge_cost_ {0};
  Cost edge_cost_gcd_ {0};
  std::vector<Amount> edge_capacity_units_ {};
  std::vector<Cost> edge_costs_ {};
  std::vector<NodeId> edge_sources_ {};
  std::vector<NodeId> edge_destinations_ {};
  std::vector<std::int32_t> input_positions_ {};
  std::vector<std::int32_t> outgoing_offsets_ {};
  std::vector<std::int32_t> incoming_offsets_ {};
  std::vector<EdgeId> incoming_edges_ {};
  // The source of each reverse slot (edge_sources()[incoming_edges()[slot]]), stored for the
  // reverse search's locality; internal (detail::GraphAccess).
  std::vector<NodeId> incoming_edge_sources_ {};
};

} // namespace netgraph::core
