/*
  StrictMultiDiGraph: the one graph type on integer units (APPROACH 2.1, 3).

  Construction validates the inputs in the documented order, then orders the edges by the
  stable key (src, cost, dst, input): a counting pass by src places every input position in
  its source's bucket (input order kept), and each bucket is sorted by (cost, dst, input). The
  rank in that order is the public EdgeId, so the out-edges of u are the id range
  [row_offsets[u], row_offsets[u + 1]). The reverse CSR lists the in-edges of v in
  (cost, src, input) order: a counting pass by dst over the ids, then each bucket sorted by
  (cost, src, id). Within one destination, equal (cost, src) means equal (src, cost, dst), whose
  ids are in input order, so (cost, src, id) is (cost, src, input).
*/
#include "netgraph/core/strict_multidigraph.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace netgraph::core {

namespace {

std::atomic<std::uint64_t> uid_counter {0};

std::uint64_t next_uid() noexcept {
  return uid_counter.fetch_add(1, std::memory_order_relaxed) + 1;
}

// The exponent k of a resolution q = 2^k with k in [-32, 32].
int resolution_exponent(double q) {
  if (!std::isfinite(q) || q <= 0.0) {
    throw std::invalid_argument("resolution must be a power of two in [2^-32, 2^32]");
  }
  int exp = 0;
  const double mantissa = std::frexp(q, &exp);  // q = mantissa * 2^exp, mantissa in [0.5, 1)
  const int k = exp - 1;
  if (mantissa != 0.5 || k < -32 || k > 32) {
    throw std::invalid_argument("resolution must be a power of two in [2^-32, 2^32]");
  }
  return k;
}

} // namespace

Amount to_units(double x, double resolution) {
  const int k = resolution_exponent(resolution);
  if (std::isnan(x)) {
    throw std::invalid_argument("to_units: NaN has no unit value");
  }
  if (std::isinf(x)) {
    throw std::invalid_argument("to_units: an infinite amount is not finite and has no unit value");
  }
  if (x < 0.0) {
    throw std::invalid_argument("to_units: a negative amount has no unit value");
  }
  // x / q is exact (q is a power of two) unless it overflows to infinity.
  const double y = std::ldexp(x, -k);
  constexpr double kTwo52 = 4503599627370496.0;  // 2^52: doubles at or above it are integers
  constexpr double kLimit = 4611686018427387904.0;  // 2^62
  double r;
  if (y >= kTwo52) {
    r = y;
  } else {
    const double f = std::floor(y);
    r = (y - f >= 0.5) ? f + 1.0 : f;  // y - f is exact below 2^52: round half up
  }
  if (!(r < kLimit)) {
    throw std::invalid_argument(
        "to_units: the amount is at or above the bound of the units domain (2^62 units)");
  }
  return static_cast<Amount>(r);
}

double as_float(Amount x, double resolution) {
  const int k = resolution_exponent(resolution);
  return std::ldexp(static_cast<double>(x), k);
}

StrictMultiDiGraph::StrictMultiDiGraph(std::int32_t num_nodes,
                                       std::span<const NodeId> src,
                                       std::span<const NodeId> dst,
                                       std::span<const Amount> capacity_units,
                                       std::span<const Cost> cost,
                                       double resolution) {
  (void)resolution_exponent(resolution);
  // The unit values.
  for (std::size_t i = 0; i < capacity_units.size(); ++i) {
    const Amount c = capacity_units[i];
    if (c < 0 || c >= kDomainLimit) {
      throw std::invalid_argument("capacity_units[" + std::to_string(i) + "] = " + std::to_string(c) +
                                  " is outside the units domain [0, 2^62)");
    }
  }
  // The array shapes.
  if (num_nodes < 0) {
    throw std::invalid_argument("num_nodes must be >= 0 (got " + std::to_string(num_nodes) + ")");
  }
  const std::size_t m = src.size();
  if (dst.size() != m || capacity_units.size() != m || cost.size() != m) {
    throw std::invalid_argument("src, dst, capacity_units and cost must have the same length (got " +
                                std::to_string(src.size()) + ", " + std::to_string(dst.size()) + ", " +
                                std::to_string(capacity_units.size()) + ", " +
                                std::to_string(cost.size()) + ")");
  }
  if (m > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    throw std::invalid_argument("the number of edges exceeds INT32_MAX (array length)");
  }
  // Zero or negative costs (rule 28).
  Cost max_cost = 0;
  Cost cost_gcd = 0;
  for (std::size_t i = 0; i < m; ++i) {
    if (cost[i] < 1) {
      throw std::invalid_argument("cost[" + std::to_string(i) + "] = " + std::to_string(cost[i]) +
                                  ": costs must be positive (cost >= 1)");
    }
    max_cost = std::max(max_cost, cost[i]);
    cost_gcd = std::gcd(cost_gcd, cost[i]);
  }
  // The remaining bounds.
  for (std::size_t i = 0; i < m; ++i) {
    if (src[i] < 0 || src[i] >= num_nodes || dst[i] < 0 || dst[i] >= num_nodes) {
      throw std::invalid_argument("edge " + std::to_string(i) + ": endpoint (" + std::to_string(src[i]) +
                                  ", " + std::to_string(dst[i]) + ") out of node range [0, " +
                                  std::to_string(num_nodes) + ")");
    }
  }
  constexpr auto kLimit = static_cast<std::uint64_t>(kDomainLimit);
  std::uint64_t cap_sum = 0;
  std::uint64_t cost_sum = 0;
  for (std::size_t i = 0; i < m; ++i) {
    cap_sum += static_cast<std::uint64_t>(capacity_units[i]);  // each < 2^62: no wrap before the check
    if (cap_sum >= kLimit) {
      throw std::invalid_argument("the sum of all capacities reaches the bound of the units domain (2^62 units)");
    }
    cost_sum += static_cast<std::uint64_t>(cost[i]);
    if (cost_sum >= kLimit) {
      throw std::invalid_argument(
          "the sum of all edge costs reaches the bound 2^62: path costs would overflow int64");
    }
  }

  uid_ = next_uid();
  resolution_ = resolution;
  num_nodes_ = num_nodes;
  max_cost_ = max_cost;
  cost_gcd_ = cost_gcd;
  const auto n = static_cast<std::size_t>(num_nodes);

  // Ids: rank in (src, cost, dst, input).
  row_offsets_.assign(n + 1, 0);
  for (std::size_t i = 0; i < m; ++i) ++row_offsets_[static_cast<std::size_t>(src[i]) + 1];
  for (std::size_t u = 0; u < n; ++u) row_offsets_[u + 1] += row_offsets_[u];
  input_index_.resize(m);
  {
    std::vector<std::int32_t> cursor(row_offsets_.begin(), row_offsets_.end() - 1);
    for (std::size_t i = 0; i < m; ++i) {
      input_index_[static_cast<std::size_t>(cursor[static_cast<std::size_t>(src[i])]++)] =
          static_cast<std::int32_t>(i);
    }
  }
  for (std::size_t u = 0; u < n; ++u) {
    auto first = input_index_.begin() + row_offsets_[u];
    auto last = input_index_.begin() + row_offsets_[u + 1];
    if (last - first < 2) continue;
    std::sort(first, last, [&](std::int32_t a, std::int32_t b) {
      const auto ia = static_cast<std::size_t>(a);
      const auto ib = static_cast<std::size_t>(b);
      if (cost[ia] != cost[ib]) return cost[ia] < cost[ib];
      if (dst[ia] != dst[ib]) return dst[ia] < dst[ib];
      return a < b;
    });
  }
  capacity_.resize(m);
  cost_.resize(m);
  src_.resize(m);
  dst_.resize(m);
  for (std::size_t e = 0; e < m; ++e) {
    const auto i = static_cast<std::size_t>(input_index_[e]);
    capacity_[e] = capacity_units[i];
    cost_[e] = cost[i];
    src_[e] = src[i];
    dst_[e] = dst[i];
  }

  // Reverse CSR in (cost, src, input) order.
  in_row_offsets_.assign(n + 1, 0);
  for (std::size_t e = 0; e < m; ++e) ++in_row_offsets_[static_cast<std::size_t>(dst_[e]) + 1];
  for (std::size_t v = 0; v < n; ++v) in_row_offsets_[v + 1] += in_row_offsets_[v];
  in_adj_edge_index_.resize(m);
  {
    std::vector<std::int32_t> cursor(in_row_offsets_.begin(), in_row_offsets_.end() - 1);
    for (std::size_t e = 0; e < m; ++e) {
      in_adj_edge_index_[static_cast<std::size_t>(cursor[static_cast<std::size_t>(dst_[e])]++)] =
          static_cast<EdgeId>(e);
    }
  }
  for (std::size_t v = 0; v < n; ++v) {
    auto first = in_adj_edge_index_.begin() + in_row_offsets_[v];
    auto last = in_adj_edge_index_.begin() + in_row_offsets_[v + 1];
    if (last - first < 2) continue;
    std::sort(first, last, [&](EdgeId a, EdgeId b) {
      const auto ia = static_cast<std::size_t>(a);
      const auto ib = static_cast<std::size_t>(b);
      if (cost_[ia] != cost_[ib]) return cost_[ia] < cost_[ib];
      if (src_[ia] != src_[ib]) return src_[ia] < src_[ib];
      return a < b;
    });
  }
}

StrictMultiDiGraph::StrictMultiDiGraph(const StrictMultiDiGraph& other)
    : uid_(next_uid()),
      resolution_(other.resolution_),
      num_nodes_(other.num_nodes_),
      max_cost_(other.max_cost_),
      cost_gcd_(other.cost_gcd_),
      capacity_(other.capacity_),
      cost_(other.cost_),
      src_(other.src_),
      dst_(other.dst_),
      input_index_(other.input_index_),
      row_offsets_(other.row_offsets_),
      in_row_offsets_(other.in_row_offsets_),
      in_adj_edge_index_(other.in_adj_edge_index_) {}

StrictMultiDiGraph::StrictMultiDiGraph(StrictMultiDiGraph&& other) noexcept
    : uid_(std::exchange(other.uid_, 0)),
      resolution_(other.resolution_),
      num_nodes_(std::exchange(other.num_nodes_, 0)),
      max_cost_(std::exchange(other.max_cost_, 0)),
      cost_gcd_(std::exchange(other.cost_gcd_, 0)),
      capacity_(std::move(other.capacity_)),
      cost_(std::move(other.cost_)),
      src_(std::move(other.src_)),
      dst_(std::move(other.dst_)),
      input_index_(std::move(other.input_index_)),
      row_offsets_(std::move(other.row_offsets_)),
      in_row_offsets_(std::move(other.in_row_offsets_)),
      in_adj_edge_index_(std::move(other.in_adj_edge_index_)) {
  other.row_offsets_.assign(1, 0);
  other.in_row_offsets_.assign(1, 0);
}

Eligibility StrictMultiDiGraph::eligibility(std::span<const NodeId> excluded_nodes,
                                            std::span<const EdgeId> excluded_edges) const {
  const auto n = static_cast<std::size_t>(num_nodes_);
  const auto m = capacity_.size();
  for (NodeId v : excluded_nodes) {
    if (v < 0 || static_cast<std::size_t>(v) >= n) {
      throw std::invalid_argument("excluded node " + std::to_string(v) + " out of node range [0, " +
                                  std::to_string(n) + ")");
    }
  }
  for (EdgeId e : excluded_edges) {
    if (e < 0 || static_cast<std::size_t>(e) >= m) {
      throw std::invalid_argument("excluded edge " + std::to_string(e) + " out of edge range [0, " +
                                  std::to_string(m) + ")");
    }
  }
  Eligibility el;
  el.graph_uid_ = uid_;
  el.node_allowed_.assign(n, 1);
  el.edge_allowed_.assign(m, 1);
  for (EdgeId e : excluded_edges) el.edge_allowed_[static_cast<std::size_t>(e)] = 0;
  for (NodeId v : excluded_nodes) {
    const auto vi = static_cast<std::size_t>(v);
    if (el.node_allowed_[vi] == 0) continue;  // duplicate
    el.node_allowed_[vi] = 0;
    for (auto e = row_offsets_[vi]; e < row_offsets_[vi + 1]; ++e) {
      el.edge_allowed_[static_cast<std::size_t>(e)] = 0;
    }
    for (auto s = in_row_offsets_[vi]; s < in_row_offsets_[vi + 1]; ++s) {
      el.edge_allowed_[static_cast<std::size_t>(in_adj_edge_index_[static_cast<std::size_t>(s)])] = 0;
    }
  }
  el.uid_ = (excluded_nodes.empty() && excluded_edges.empty()) ? 0 : next_uid();
  return el;
}

} // namespace netgraph::core
