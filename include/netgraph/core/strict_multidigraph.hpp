/* The one graph type: an immutable directed multigraph on integer units (APPROACH 3). */
#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "netgraph/core/types.hpp"

namespace netgraph::core {

class Eligibility;

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
// A copy is a new graph: same ids and input map, a new uid. Uids come from one process-wide
// 64-bit counter and are never reused.
class StrictMultiDiGraph {
public:
  StrictMultiDiGraph(std::int32_t num_nodes,
                     std::span<const NodeId> src,
                     std::span<const NodeId> dst,
                     std::span<const Amount> capacity_units,
                     std::span<const Cost> cost,
                     double resolution);
  StrictMultiDiGraph(const StrictMultiDiGraph& other);
  StrictMultiDiGraph(StrictMultiDiGraph&& other) noexcept;
  StrictMultiDiGraph& operator=(const StrictMultiDiGraph&) = delete;
  StrictMultiDiGraph& operator=(StrictMultiDiGraph&&) = delete;
  ~StrictMultiDiGraph() noexcept = default;

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
  [[nodiscard]] std::span<const std::int32_t> row_offsets() const noexcept { return row_offsets_; }
  // Reverse CSR: the in-edges of v are in_adj_edge_index()[in_row_offsets()[v] ..
  // in_row_offsets()[v + 1]) (public EdgeIds); the reverse column is
  // edge_src()[in_adj_edge_index()[slot]].
  [[nodiscard]] std::span<const std::int32_t> in_row_offsets() const noexcept { return in_row_offsets_; }
  [[nodiscard]] std::span<const EdgeId> in_adj_edge_index() const noexcept { return in_adj_edge_index_; }

  // Static eligibility from excluded ids (duplicates and unsorted lists accepted, out-of-range
  // ids rejected). Empty lists mean all allowed (uid 0).
  [[nodiscard]] Eligibility eligibility(std::span<const NodeId> excluded_nodes = {},
                                        std::span<const EdgeId> excluded_edges = {}) const;

private:
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
};

// Immutable, graph-bound static eligibility: dense node and edge flags. An edge is allowed when
// it is not excluded and both of its endpoints are allowed.
class Eligibility {
public:
  [[nodiscard]] std::uint64_t graph_uid() const noexcept { return graph_uid_; }
  // 0 for the all-allowed value; a fresh uid otherwise.
  [[nodiscard]] std::uint64_t uid() const noexcept { return uid_; }
  [[nodiscard]] std::span<const std::uint8_t> node_allowed() const noexcept { return node_allowed_; }
  [[nodiscard]] std::span<const std::uint8_t> edge_allowed() const noexcept { return edge_allowed_; }
  [[nodiscard]] bool node_ok(NodeId v) const noexcept {
    return node_allowed_[static_cast<std::size_t>(v)] != 0;
  }
  [[nodiscard]] bool edge_ok(EdgeId e) const noexcept {
    return edge_allowed_[static_cast<std::size_t>(e)] != 0;
  }

private:
  friend class StrictMultiDiGraph;
  Eligibility() = default;
  std::uint64_t graph_uid_ {0};
  std::uint64_t uid_ {0};
  std::vector<std::uint8_t> node_allowed_ {};
  std::vector<std::uint8_t> edge_allowed_ {};
};

} // namespace netgraph::core
