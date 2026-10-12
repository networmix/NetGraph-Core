/*
  Eligibility: the static node and edge flags of one graph, its one public construction from
  excluded ids, Graph::eligibility(excluded_nodes, excluded_edges), and intersect. An excluded
  node excludes its incident edges (its out-edge id range and its reverse-CSR in-edges);
  intersect keeps an operand, and its uid, when the result equals it.
*/
#include "netgraph/core/eligibility.hpp"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

#include "netgraph/core/detail/uid.hpp"
#include "netgraph/core/graph.hpp"

namespace netgraph::core {

namespace {

using detail::next_uid;

} // namespace

Eligibility Graph::eligibility(std::span<const NodeId> excluded_nodes,
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

Eligibility Eligibility::intersect(const Eligibility& a, const Eligibility& b) {
  if (a.graph_uid_ != b.graph_uid_) {
    throw std::invalid_argument("Eligibility.intersect: the two values belong to different graphs");
  }
  Eligibility out;
  out.graph_uid_ = a.graph_uid_;
  out.node_allowed_.resize(a.node_allowed_.size());
  out.edge_allowed_.resize(a.edge_allowed_.size());
  for (std::size_t v = 0; v < out.node_allowed_.size(); ++v) {
    out.node_allowed_[v] = static_cast<std::uint8_t>(a.node_allowed_[v] & b.node_allowed_[v]);
  }
  for (std::size_t e = 0; e < out.edge_allowed_.size(); ++e) {
    out.edge_allowed_[e] = static_cast<std::uint8_t>(a.edge_allowed_[e] & b.edge_allowed_[e]);
  }
  if (out.node_allowed_ == a.node_allowed_ && out.edge_allowed_ == a.edge_allowed_) return a;
  if (out.node_allowed_ == b.node_allowed_ && out.edge_allowed_ == b.edge_allowed_) return b;
  out.uid_ = next_uid();
  return out;
}

} // namespace netgraph::core
