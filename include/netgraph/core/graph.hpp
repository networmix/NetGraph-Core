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
// (src, cost, dst, input). Node u's out-edges are [row_offsets()[u], row_offsets()[u + 1]),
// ordered by (cost, dst, input); input_index()[e] gives e's input position. All other per-edge
// arrays use EdgeId order, except in_adj_edge_index() (reverse-CSR slot order). CSR (compressed
// sparse row) stores each node's edges in a flat array and uses offsets to mark its range.
// Reverse CSR lists v's in-edges in (cost, src, input) order.
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
// dst, cost rank, src give EdgeId order; cost rank, dst give reverse CSR. When N K < E, use
// dst then (src, rank) forward, and (dst, rank) reverse. With no edges, only validation and two
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
        std::span<const NodeId> src,
        std::span<const NodeId> dst,
        std::span<const Amount> capacity_units,
        std::span<const Cost> cost,
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
    return static_cast<std::int32_t>(capacity_.size());
  }
  [[nodiscard]] std::span<const Amount> capacity_units() const noexcept { return capacity_; }
  [[nodiscard]] std::span<const Cost> cost() const noexcept { return cost_; }
  [[nodiscard]] std::span<const NodeId> edge_src() const noexcept { return src_; }
  [[nodiscard]] std::span<const NodeId> edge_dst() const noexcept { return dst_; }
  // Maximum edge cost and greatest common divisor of all edge costs, computed by full scans
  // at construction; both are 0 on a graph without edges.
  [[nodiscard]] Cost max_cost() const noexcept { return max_cost_; }
  [[nodiscard]] Cost cost_gcd() const noexcept { return cost_gcd_; }
  // id -> input position. The inverse permutation is the caller's (one O(E) scatter).
  [[nodiscard]] std::span<const std::int32_t> input_index() const noexcept { return input_index_; }
  // Forward CSR offsets (N + 1 entries, [0] = 0, [N] = E): the id ranges of the out-edges.
  [[nodiscard]] std::span<const std::int32_t> row_offsets() const noexcept { return offsets(row_offsets_); }
  // Reverse CSR: the in-edges of v are in_adj_edge_index()[in_row_offsets()[v] ..
  // in_row_offsets()[v + 1]) (public EdgeIds); the reverse column is
  // edge_src()[in_adj_edge_index()[slot]].
  [[nodiscard]] std::span<const std::int32_t> in_row_offsets() const noexcept { return offsets(in_row_offsets_); }
  [[nodiscard]] std::span<const EdgeId> in_adj_edge_index() const noexcept { return in_adj_edge_index_; }

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
  Cost max_cost_ {0};
  Cost cost_gcd_ {0};
  std::vector<Amount> capacity_ {};
  std::vector<Cost> cost_ {};
  std::vector<NodeId> src_ {};
  std::vector<NodeId> dst_ {};
  std::vector<std::int32_t> input_index_ {};
  std::vector<std::int32_t> row_offsets_ {};
  std::vector<std::int32_t> in_row_offsets_ {};
  std::vector<EdgeId> in_adj_edge_index_ {};
  // The source of each reverse slot (edge_src()[in_adj_edge_index()[slot]]), stored for the
  // reverse search's locality; internal (detail::GraphAccess).
  std::vector<NodeId> in_col_indices_ {};
};

} // namespace netgraph::core
