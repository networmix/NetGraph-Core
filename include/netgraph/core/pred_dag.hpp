/* Route structures: the immutable predecessor DAG and the explicit Path (APPROACH 5). */
#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

namespace detail {
struct PredDAGBuilder;
} // namespace detail

// An immutable, graph-bound DAG in predecessor CSR form: the predecessor entries of node v are
// parents()[k] and via_edges()[k] for k in [parent_offsets()[v], parent_offsets()[v + 1]); the
// entry (u, e) means the edge e = u -> v belongs to the DAG. Within a node the entries are in
// strictly ascending via EdgeId, so a DAG value has one canonical form. Equality is structural
// (same graph, same CSR). Every API that takes a route structure takes a PredDAG; a DAG of
// another graph is rejected wherever it is used.
class PredDAG {
public:
  // Validates, in this order, before any traversal (std::invalid_argument on the first defect):
  // offsets of length N + 1, offsets[0] == 0, monotone, offsets[N] == len(parents) ==
  // len(via); then every parent and via in range, every via joining (parent, child), no
  // duplicate (child, via), acyclic; within a node the entries in ascending via EdgeId.
  [[nodiscard]] static PredDAG from_csr(const StrictMultiDiGraph& g, std::span<const std::int32_t> offsets,
                                        std::span<const NodeId> parents, std::span<const EdgeId> via);

  [[nodiscard]] std::uint64_t graph_uid() const noexcept { return graph_uid_; }
  [[nodiscard]] std::int32_t num_nodes() const noexcept { return static_cast<std::int32_t>(offsets_.size()) - 1; }
  [[nodiscard]] std::span<const std::int32_t> parent_offsets() const noexcept { return offsets_; }
  [[nodiscard]] std::span<const NodeId> parents() const noexcept { return parents_; }
  [[nodiscard]] std::span<const EdgeId> via_edges() const noexcept { return via_; }
  // The member edges in ascending EdgeId.
  [[nodiscard]] std::vector<EdgeId> edges() const { return edges_; }

  friend bool operator==(const PredDAG& a, const PredDAG& b) noexcept {
    return a.graph_uid_ == b.graph_uid_ && a.offsets_ == b.offsets_ && a.parents_ == b.parents_ && a.via_ == b.via_;
  }

private:
  friend struct detail::PredDAGBuilder;
  PredDAG() = default;
  std::uint64_t graph_uid_ {0};
  std::vector<std::int32_t> offsets_ {0};
  std::vector<NodeId> parents_ {};
  std::vector<EdgeId> via_ {};
  std::vector<EdgeId> edges_ {};  // via_ in ascending EdgeId (derived once, at construction)
};

// An explicit walk of hops; a hop holds one or more parallel members of one (u, v) pair, in
// the caller's order. src, dst and cost (the sum over hops of the minimum member cost) are
// derived at construction. A single-member walk is [[e1], [e2], ...].
class Path {
public:
  // Validates, in this order (std::invalid_argument on the first defect): every member id in
  // range (before any endpoint is read), no duplicate member, a non-empty hop list, every hop
  // non-empty, the members of a hop sharing (u, v), contiguous hops, no repeated node.
  Path(const StrictMultiDiGraph& g, std::vector<std::vector<EdgeId>> hops);

  [[nodiscard]] std::uint64_t graph_uid() const noexcept { return graph_uid_; }
  [[nodiscard]] const std::vector<std::vector<EdgeId>>& hops() const noexcept { return hops_; }
  [[nodiscard]] NodeId src() const noexcept { return nodes_.front(); }
  [[nodiscard]] NodeId dst() const noexcept { return nodes_.back(); }
  // The visited nodes, src first (hops + 1 entries).
  [[nodiscard]] const std::vector<NodeId>& nodes() const noexcept { return nodes_; }
  [[nodiscard]] Cost cost() const noexcept { return cost_; }
  // The equivalent PredDAG (every member of every hop; entries in ascending via EdgeId).
  [[nodiscard]] PredDAG to_dag() const;

  friend bool operator==(const Path& a, const Path& b) noexcept {
    return a.graph_uid_ == b.graph_uid_ && a.hops_ == b.hops_;
  }

private:
  std::uint64_t graph_uid_ {0};
  std::int32_t num_nodes_ {0};
  std::vector<std::vector<EdgeId>> hops_ {};
  std::vector<NodeId> nodes_ {};
  Cost cost_ {0};
};

} // namespace netgraph::core
