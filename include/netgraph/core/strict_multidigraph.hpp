#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "netgraph/core/types.hpp"

namespace netgraph::core {

// Notes on edge identifiers:
// - EdgeId refers to the index of an edge in the graph's internal, compacted
//   representation. Edges are deterministically reordered during construction
//   by (src, dst, cost, link_id) for stable traversal and performance.
// - link_id_of(EdgeId) returns the original user-supplied link identifier
//   (typically the edge "key" in higher-level APIs). Use this to correlate
//   compacted edges back to external expectations and tests.

class StrictMultiDiGraph {
public:
  static StrictMultiDiGraph from_arrays(
      std::int32_t num_nodes,
      std::span<const std::int32_t> src,
      std::span<const std::int32_t> dst,
      std::span<const double> capacity,
      std::span<const double> cost,
      std::span<const std::int64_t> link_ids,
      bool add_reverse);

  std::int32_t num_nodes() const noexcept { return num_nodes_; }
  std::int32_t num_edges() const noexcept { return static_cast<std::int32_t>(edges_); }
  std::int64_t  link_id_of(EdgeId eid) const noexcept;

  std::span<const double> capacity_view() const noexcept { return capacity_; }
  std::span<const double> cost_view() const noexcept { return cost_; }
  std::span<const std::int32_t> row_offsets_view() const noexcept { return row_offsets_; }
  std::span<const std::int32_t> col_indices_view() const noexcept { return col_indices_; }
  std::span<const std::int32_t> adj_edge_index_view() const noexcept { return adj_edge_index_; }
  // Reverse CSR (incoming edges): for each node v, incoming neighbors are stored
  std::span<const std::int32_t> in_row_offsets_view() const noexcept { return in_row_offsets_; }
  std::span<const std::int32_t> in_col_indices_view() const noexcept { return in_col_indices_; }
  std::span<const std::int32_t> in_adj_edge_index_view() const noexcept { return in_adj_edge_index_; }

private:
  // Core storage (edges may be reordered in compact form)
  std::int32_t num_nodes_ {0};
  std::size_t edges_ {0};
  std::vector<std::int64_t> link_ids_ {};
  std::vector<double> capacity_ {};
  std::vector<double> cost_ {};
  std::vector<std::int32_t> src_ {};
  std::vector<std::int32_t> dst_ {};

  // CSR adjacency for deterministic traversal (always built deterministically)
  std::vector<std::int32_t> row_offsets_ {};
  std::vector<std::int32_t> col_indices_ {};
  std::vector<std::int32_t> adj_edge_index_ {}; // map CSR entry -> EdgeId
  // Reverse CSR (incoming adjacency)
  std::vector<std::int32_t> in_row_offsets_ {};
  std::vector<std::int32_t> in_col_indices_ {};
  std::vector<std::int32_t> in_adj_edge_index_ {};
};

} // namespace netgraph::core
