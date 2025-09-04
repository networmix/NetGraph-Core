/*
  StrictMultiDiGraph — immutable directed multigraph with deterministic layout.

  Construction from arrays validates inputs, optionally duplicates reverse
  edges, and compacts data into CSR adjacency (and reverse CSR) using a stable
  ordering. Cost-first sorting keeps equal-cost tiers clustered for
  deterministic behavior.
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
    bool add_reverse) {

  if (num_nodes < 0) {
    throw std::invalid_argument("num_nodes must be >= 0");
  }
  if (src.size() != dst.size() || src.size() != capacity.size() || src.size() != cost.size()) {
    throw std::invalid_argument("src, dst, capacity, and cost must have the same length");
  }
  StrictMultiDiGraph g;
  g.num_nodes_ = num_nodes;
  std::size_t m = src.size();

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
  // External link identifiers are not used in core.

  // Optionally add reverse edges
  if (add_reverse) {
    std::size_t old_m = m;
    src_v.reserve(2 * m);
    dst_v.reserve(2 * m);
    cap_v.reserve(2 * m);
    cost_v.reserve(2 * m);
    for (std::size_t i = 0; i < old_m; ++i) {
      auto s = src_v[i];
      auto d = dst_v[i];
      src_v.push_back(d);
      dst_v.push_back(s);
      cap_v.push_back(cap_v[i]);
      cost_v.push_back(cost_v[i]);
    }
    m = src_v.size();
  }

  // Compact/sort edges deterministically. Sort primarily by cost, then by
  // (src, dst) to keep equal-cost tiers clustered for reproducibility.
  std::vector<std::size_t> idx(m);
  std::iota(idx.begin(), idx.end(), 0);
  // Cost-first ordering keeps equal-cost tiers clustered; within the same
  // cost, maintain stable source/destination ordering for reproducibility.
  std::stable_sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
    if (cost_v[a] != cost_v[b]) return cost_v[a] < cost_v[b];
    if (src_v[a] != src_v[b]) return src_v[a] < src_v[b];
    return dst_v[a] < dst_v[b];
  });
  auto apply_perm = [&](auto& out_vec, const auto& in_vec) {
    out_vec.resize(m);
    for (std::size_t i = 0; i < m; ++i) out_vec[i] = in_vec[idx[i]];
  };
  apply_perm(g.src_, src_v);
  apply_perm(g.dst_, dst_v);
  apply_perm(g.capacity_, cap_v);
  apply_perm(g.cost_, cost_v);
  g.edges_ = m;

  // Build CSR adjacency
  g.row_offsets_.assign(static_cast<std::size_t>(num_nodes) + 1, 0);
  for (std::size_t i = 0; i < m; ++i) {
    g.row_offsets_[static_cast<std::size_t>(g.src_[i]) + 1]++;
  }
  for (std::size_t i = 1; i < g.row_offsets_.size(); ++i) {
    g.row_offsets_[i] += g.row_offsets_[i - 1];
  }
  g.col_indices_.resize(m);
  g.adj_edge_index_.resize(m);
  // We need a copy of offsets to fill in-place
  std::vector<std::int32_t> cursor = g.row_offsets_;
  for (std::size_t e = 0; e < m; ++e) {
    auto u = g.src_[e];
    auto pos = static_cast<std::size_t>(cursor[static_cast<std::size_t>(u)]++);
    g.col_indices_[pos] = g.dst_[e];
    g.adj_edge_index_[pos] = static_cast<EdgeId>(e);
  }
  // Build reverse CSR (incoming adjacency)
  g.in_row_offsets_.assign(static_cast<std::size_t>(num_nodes) + 1, 0);
  for (std::size_t i = 0; i < m; ++i) {
    g.in_row_offsets_[static_cast<std::size_t>(g.dst_[i]) + 1]++;
  }
  for (std::size_t i = 1; i < g.in_row_offsets_.size(); ++i) {
    g.in_row_offsets_[i] += g.in_row_offsets_[i - 1];
  }
  g.in_col_indices_.resize(m);
  g.in_adj_edge_index_.resize(m);
  std::vector<std::int32_t> rcursor = g.in_row_offsets_;
  for (std::size_t e = 0; e < m; ++e) {
    auto v = g.dst_[e];
    auto pos = static_cast<std::size_t>(rcursor[static_cast<std::size_t>(v)]++);
    g.in_col_indices_[pos] = g.src_[e];
    g.in_adj_edge_index_[pos] = static_cast<EdgeId>(e);
  }
  return g;
}

} // namespace netgraph::core
