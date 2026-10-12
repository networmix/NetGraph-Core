/* Eligibility: which nodes and edges of one graph may be used.
   Built only by Graph::eligibility(excluded_nodes, excluded_edges) (graph.hpp, which includes
   this header) and by intersect, which combines two values of one graph. */
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "netgraph/core/types.hpp"

namespace netgraph::core {

class Graph;

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
  // The intersection of two values of one graph (std::invalid_argument for values of different
  // graphs): an operand when the result equals it (a first, then b), keeping that operand's
  // uid, so a consumer builds each state once and reuses it; otherwise a new value with a fresh
  // uid (equal intersections built separately get distinct uids).
  [[nodiscard]] static Eligibility intersect(const Eligibility& a, const Eligibility& b);

private:
  friend class Graph;
  Eligibility() = default;
  std::uint64_t graph_uid_ {0};
  std::uint64_t uid_ {0};
  std::vector<std::uint8_t> node_allowed_ {};
  std::vector<std::uint8_t> edge_allowed_ {};
};

} // namespace netgraph::core
