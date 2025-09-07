/* Immutable directed multigraph with CSR and reverse CSR adjacency. */
#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "netgraph/core/types.hpp"

namespace netgraph::core {

// Notes on edge identifiers:
// - EdgeId refers to the index of an edge in the graph's internal, compacted
//   representation. Edges are deterministically reordered during construction
//   by (src, dst, cost) for stable traversal and performance.
// - External link identifiers are managed at higher layers; inside core,
//   the canonical edge identifier is EdgeId (compacted index).

class StrictMultiDiGraph {
public:
  [[nodiscard]] static StrictMultiDiGraph from_arrays(
      std::int32_t num_nodes,
      std::span<const std::int32_t> src,
      std::span<const std::int32_t> dst,
      std::span<const Cap> capacity,
      std::span<const Cost> cost,
      bool add_reverse);
  ~StrictMultiDiGraph() noexcept = default;

  [[nodiscard]] std::int32_t num_nodes() const noexcept { return num_nodes_; }
  [[nodiscard]] std::int32_t num_edges() const noexcept { return static_cast<std::int32_t>(edges_); }

  [[nodiscard]] std::span<const Cap> capacity_view() const noexcept { return capacity_; }
  [[nodiscard]] std::span<const Cost> cost_view() const noexcept { return cost_; }
  [[nodiscard]] std::span<const NodeId> edge_src_view() const noexcept { return src_; }
  [[nodiscard]] std::span<const NodeId> edge_dst_view() const noexcept { return dst_; }
  [[nodiscard]] std::span<const std::int32_t> row_offsets_view() const noexcept { return row_offsets_; }
  [[nodiscard]] std::span<const NodeId> col_indices_view() const noexcept { return col_indices_; }
  [[nodiscard]] std::span<const EdgeId> adj_edge_index_view() const noexcept { return adj_edge_index_; }
  // Reverse CSR (incoming edges): for each node v, incoming neighbors are stored
  [[nodiscard]] std::span<const std::int32_t> in_row_offsets_view() const noexcept { return in_row_offsets_; }
  [[nodiscard]] std::span<const NodeId> in_col_indices_view() const noexcept { return in_col_indices_; }
  [[nodiscard]] std::span<const EdgeId> in_adj_edge_index_view() const noexcept { return in_adj_edge_index_; }

private:
  // Core storage (edges may be reordered in compact form)
  std::int32_t num_nodes_ {0};
  std::size_t edges_ {0};
  // External link identifiers are no longer stored; EdgeId is the canonical id
  std::vector<Cap> capacity_ {};
  std::vector<Cost> cost_ {};
  std::vector<NodeId> src_ {};
  std::vector<NodeId> dst_ {};

  // CSR adjacency for deterministic traversal (always built deterministically)
  std::vector<std::int32_t> row_offsets_ {};
  std::vector<NodeId> col_indices_ {};
  std::vector<EdgeId> adj_edge_index_ {}; // map CSR entry -> EdgeId
  // Reverse CSR (incoming adjacency)
  std::vector<std::int32_t> in_row_offsets_ {};
  std::vector<NodeId> in_col_indices_ {};
  std::vector<EdgeId> in_adj_edge_index_ {};
};

} // namespace netgraph::core
