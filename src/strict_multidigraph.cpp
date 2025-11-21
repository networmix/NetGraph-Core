/*
  StrictMultiDiGraph — immutable directed multigraph with deterministic layout.

  Construction from arrays validates inputs and builds CSR and reverse CSR
  adjacency structures. Edges are reordered by (cost, src, dst) to cluster
  equal-cost tiers for deterministic traversal.
*/
#include "netgraph/core/strict_multidigraph.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace netgraph::core {

StrictMultiDiGraph StrictMultiDiGraph::from_arrays(
    std::int32_t num_nodes,
    std::span<const std::int32_t> src,
    std::span<const std::int32_t> dst,
    std::span<const Cap> capacity,
    std::span<const Cost> cost,
    std::span<const std::int64_t> ext_edge_ids) {

  if (num_nodes < 0) {
    throw std::invalid_argument("num_nodes must be >= 0");
  }
  if (src.size() != dst.size() || src.size() != capacity.size() || src.size() != cost.size()) {
    throw std::invalid_argument("src, dst, capacity, and cost must have the same length");
  }
  if (!ext_edge_ids.empty() && ext_edge_ids.size() != src.size()) {
    throw std::invalid_argument("ext_edge_ids must be empty or have the same length as src/dst/capacity/cost");
  }
  StrictMultiDiGraph g;
  g.num_nodes_ = num_nodes;
  std::size_t m = src.size();
  if (m > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    throw std::invalid_argument("number of edges exceeds INT32_MAX");
  }

  // Invariants: ids within [0, num_nodes), non-negative weights
  for (std::size_t i = 0; i < m; ++i) {
    if (src[i] < 0 || dst[i] < 0 || src[i] >= num_nodes || dst[i] >= num_nodes) {
      throw std::out_of_range("edge index out of range of num_nodes");
    }
    if (capacity[i] < 0.0) {
      throw std::invalid_argument("capacity must be >= 0");
    }
    if (cost[i] < 0) { throw std::invalid_argument("cost must be >= 0"); }
  }
  // Gather initial arrays
  std::vector<NodeId> src_v(src.begin(), src.end());
  std::vector<NodeId> dst_v(dst.begin(), dst.end());
  std::vector<Cap> cap_v(capacity.begin(), capacity.end());
  std::vector<Cost> cost_v(cost.begin(), cost.end());
  std::vector<std::int64_t> ext_v;
  if (!ext_edge_ids.empty()) {
    ext_v.assign(ext_edge_ids.begin(), ext_edge_ids.end());
  }

  // Sort edges by (cost, src, dst) to cluster equal-cost tiers for deterministic traversal.
  // Create an index array [0, 1, 2, ..., m-1] that we'll permute to define the reordering.
  std::vector<std::size_t> idx(m);
  std::iota(idx.begin(), idx.end(), 0);  // iota fills with sequential values

  // Cost-first ordering keeps equal-cost tiers clustered; within the same
  // cost, maintain stable source/destination ordering for reproducibility.
  // stable_sort preserves relative order of equal elements, ensuring determinism.
  // Lambda comparator: [&] captures variables by reference, (a, b) are indices to compare.
  std::stable_sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
    if (cost_v[a] != cost_v[b]) return cost_v[a] < cost_v[b];
    if (src_v[a] != src_v[b]) return src_v[a] < src_v[b];
    return dst_v[a] < dst_v[b];
  });

  // Reorder all edge arrays using the computed permutation.
  // Generic lambda (auto parameters) works with any vector type.
  auto apply_perm = [&](auto& out_vec, const auto& in_vec) {
    out_vec.resize(m);
    for (std::size_t i = 0; i < m; ++i) out_vec[i] = in_vec[idx[i]];
  };
  apply_perm(g.src_, src_v);
  apply_perm(g.dst_, dst_v);
  apply_perm(g.capacity_, cap_v);
  apply_perm(g.cost_, cost_v);
  if (!ext_v.empty()) {
    apply_perm(g.ext_edge_ids_, ext_v);
  }
  g.edges_ = m;

  // Build CSR (Compressed Sparse Row) adjacency structure.
  // For each node u, row_offsets_[u]:row_offsets_[u+1] gives the range in col_indices_
  // that stores u's outgoing neighbors, and adj_edge_index_ stores the corresponding EdgeIds.

  // Step 1: Count outgoing edges per node (histogram).
  g.row_offsets_.assign(static_cast<std::size_t>(num_nodes) + 1, 0);
  for (std::size_t i = 0; i < m; ++i) {
    g.row_offsets_[static_cast<std::size_t>(g.src_[i]) + 1]++;
  }

  // Step 2: Convert counts to cumulative offsets (prefix sum).
  // row_offsets_[u] is the starting index for node u's adjacency list.
  for (std::size_t i = 1; i < g.row_offsets_.size(); ++i) {
    g.row_offsets_[i] += g.row_offsets_[i - 1];
  }

  // Step 3: Fill col_indices_ and adj_edge_index_ arrays.
  g.col_indices_.resize(m);
  g.adj_edge_index_.resize(m);
  // cursor tracks the current write position for each node's adjacency list.
  std::vector<std::int32_t> cursor = g.row_offsets_;  // copy for in-place filling
  for (std::size_t e = 0; e < m; ++e) {
    auto u = g.src_[e];
    auto pos = static_cast<std::size_t>(cursor[static_cast<std::size_t>(u)]++);
    g.col_indices_[pos] = g.dst_[e];        // neighbor node
    g.adj_edge_index_[pos] = static_cast<EdgeId>(e);  // edge index
  }
  // Build reverse CSR (incoming adjacency) for predecessor tracking.
  // Same CSR construction algorithm, but keyed by destination node instead of source.
  // This enables efficient lookup of incoming edges for any node.

  // Step 1: Count incoming edges per node.
  g.in_row_offsets_.assign(static_cast<std::size_t>(num_nodes) + 1, 0);
  for (std::size_t i = 0; i < m; ++i) {
    g.in_row_offsets_[static_cast<std::size_t>(g.dst_[i]) + 1]++;
  }

  // Step 2: Convert counts to cumulative offsets.
  for (std::size_t i = 1; i < g.in_row_offsets_.size(); ++i) {
    g.in_row_offsets_[i] += g.in_row_offsets_[i - 1];
  }

  // Step 3: Fill reverse adjacency arrays.
  g.in_col_indices_.resize(m);
  g.in_adj_edge_index_.resize(m);
  std::vector<std::int32_t> rcursor = g.in_row_offsets_;
  for (std::size_t e = 0; e < m; ++e) {
    auto v = g.dst_[e];
    auto pos = static_cast<std::size_t>(rcursor[static_cast<std::size_t>(v)]++);
    g.in_col_indices_[pos] = g.src_[e];     // predecessor node
    g.in_adj_edge_index_[pos] = static_cast<EdgeId>(e);
  }
  return g;
}

} // namespace netgraph::core
