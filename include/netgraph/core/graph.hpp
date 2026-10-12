/* The one graph type: an immutable directed multigraph on integer units (APPROACH 3). */
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

// Converts a physical value to units of the resolution q, rounding half up. The only way a
// physical number enters Core. Throws std::invalid_argument for NaN, infinity, a negative value,
// a result at or above the domain bound 2^62, or a resolution that is not a power of two in
// [2^-32, 2^32].
[[nodiscard]] Amount to_units(double x, double resolution);
// Units to a physical double (x * q): exact below 2^53 units.
[[nodiscard]] double as_float(Amount x, double resolution);

// Edges are given in input order. The public EdgeId of an edge is its rank in the stable order
// (src, cost, dst, input): the out-edges of u are exactly the ids [row_offsets()[u],
// row_offsets()[u + 1]) in traversal order (cost, dst, input), and input_index()[e] is the input
// position of e. Every other per-edge array is in id order. The reverse CSR lists the in-edges
// of v in (cost, src, input) order.
//
// The constructor validates in this order and reports the first defect (std::invalid_argument):
// the resolution (a power of two in [2^-32, 2^32]); the unit values (every capacity in
// [0, 2^62)); the array shapes (num_nodes >= 0, equal lengths, at most INT32_MAX edges); zero
// or negative costs (costs are positive, rule 28); the remaining bounds (endpoints in
// [0, num_nodes), then the sum of all capacities and the sum of all costs below 2^62, so no path
// cost overflows int64). This is APPROACH 3's order (revision 25); a language binding validates
// in the same order before narrowing any value, so the reported defect does not depend on the
// representation of the input.
//
// Construction (APPROACH 3): stable counting passes by dst, cost rank and src give the id order,
// two more by cost rank and dst the reverse CSR; when N K < E, combined keys do it in two forward
// passes (dst, then (src, rank)) and one reverse pass ((dst, rank)). An edgeless graph is its
// validation and two zero-filled offset arrays. The cost ranks come from the distinct costs
// collected in a hash table while there are at most 65,536 of them and a cumulative probe budget
// of 8E + 8 x 65,536 table slots lasts, else from a sort of all costs; both give the same ranks.
// Expected linear in E + N for a bounded number of distinct costs, O(N + E log(E + 1)) at worst.
// Retained: 36 bytes per edge plus 8 (N + 1) bytes of offsets.
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
  // Full scans at construction; 0 on a graph without edges.
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
  // Eligibility (eligibility.hpp).
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
