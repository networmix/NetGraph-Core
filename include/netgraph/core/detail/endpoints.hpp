/* Endpoint and root sets (internal): one normalisation, one
   overlap check and one message for each defect, shared by the flow core, the searches, the
   ledger's bottleneck, the policies, the evaluators and the bindings (which check unnarrowed
   Python integers and report the same texts). */
#pragma once

#include <algorithm>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

inline std::string empty_set_error(const char* name) { return std::string(name) + ": the node set is empty"; }

inline std::string node_range_error(const char* name, const std::string& node, std::int32_t num_nodes) {
  return std::string(name) + ": node " + node + " out of node range [0, " + std::to_string(num_nodes) + ")";
}

inline std::string overlap_error(NodeId v, const char* a, const char* b) {
  return "node " + std::to_string(v) + " is in both endpoint sets: " + a + " and " + b + " must not overlap";
}

// node_set's checks without the set: std::invalid_argument for an empty set or an id out of range
// (the first one, in input order; `name` labels the message).
inline void check_range(std::span<const NodeId> nodes, std::int32_t num_nodes, const char* name) {
  if (nodes.empty()) throw std::invalid_argument(empty_set_error(name));
  for (NodeId v : nodes) {
    if (v < 0 || v >= num_nodes) throw std::invalid_argument(node_range_error(name, std::to_string(v), num_nodes));
  }
}

// The sorted, duplicate-free node set; std::invalid_argument for an empty set or an id out of
// range (`name` labels the message).
[[nodiscard]] inline std::vector<NodeId> node_set(std::span<const NodeId> nodes, std::int32_t num_nodes,
                                                  const char* name) {
  check_range(nodes, num_nodes, name);
  std::vector<NodeId> out(nodes.begin(), nodes.end());
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// std::invalid_argument when the two sorted sets share a vertex (the first shared one).
inline void check_disjoint(std::span<const NodeId> a, std::span<const NodeId> b, const char* a_name,
                           const char* b_name) {
  for (std::size_t i = 0, k = 0; i < a.size() && k < b.size();) {
    if (a[i] == b[k]) throw std::invalid_argument(overlap_error(a[i], a_name, b_name));
    if (a[i] < b[k]) {
      ++i;
    } else {
      ++k;
    }
  }
}

} // namespace netgraph::core::detail
