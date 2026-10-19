/* Route structures: the immutable predecessor DAG and the explicit Path. */
#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

namespace detail {
struct PredecessorDAGBuilder;
struct PathBuilder;
} // namespace detail

// An immutable predecessor DAG: a directed graph with no cycles, recording each node's
// incoming edges. It belongs to one graph. CSR (compressed sparse row) groups entries in flat
// arrays: node v uses parents()[k] and via_edges()[k] for
// k in [parent_offsets()[v], parent_offsets()[v + 1]); (u, e) means edge e = u -> v is included.
// Entries per node have strictly ascending via EdgeId, giving one canonical form. Equality
// means the same graph and CSR. Route APIs take PredecessorDAG and reject one from another graph.
//
// Storage is sparse: the entries are kept per member node (the nodes
// with at least one entry, ascending, and their offsets), never as an N + 1 offsets array, so a
// flow's route costs O(its entries) whatever the graph's size. parent_offsets() presents the
// from_csr shape on demand (a fresh N + 1 array); parents() and via_edges() are the stored
// entries (grouped by node in ascending NodeId), which are the same arrays in both forms.
class PredecessorDAG {
public:
  // Validates, in this order, before any traversal (std::invalid_argument on the first defect):
  // parent_offsets of length N + 1, parent_offsets[0] == 0, monotone, parent_offsets[N] ==
  // len(parents) == len(via_edges); then every parent and via edge in range, every via edge joining
  // (parent, child), no duplicate (child, via edge), acyclic; within a node the entries in
  // ascending via EdgeId.
  [[nodiscard]] static PredecessorDAG from_csr(const Graph& graph, std::span<const std::int32_t> parent_offsets,
                                        std::span<const NodeId> parents, std::span<const EdgeId> via);

  [[nodiscard]] std::uint64_t graph_uid() const noexcept { return graph_uid_; }
  // The from_csr parent_offsets (N + 1 entries), built on demand from the sparse storage.
  [[nodiscard]] std::vector<std::int32_t> parent_offsets() const;
  [[nodiscard]] std::span<const NodeId> parents() const noexcept { return parents_; }
  [[nodiscard]] std::span<const EdgeId> via_edges() const noexcept { return via_edges_; }

  friend bool operator==(const PredecessorDAG& a, const PredecessorDAG& b) noexcept {
    return a.graph_uid_ == b.graph_uid_ && a.num_nodes_ == b.num_nodes_ && a.nodes_ == b.nodes_ &&
           a.node_offsets_ == b.node_offsets_ && a.parents_ == b.parents_ && a.via_edges_ == b.via_edges_;
  }

private:
  friend struct detail::PredecessorDAGBuilder;
  PredecessorDAG() = default;
  std::uint64_t graph_uid_ {0};
  std::int32_t num_nodes_ {0};
  std::vector<NodeId> nodes_ {};                // nodes with entries, ascending
  std::vector<std::int32_t> node_offsets_ {0};  // nodes_.size() + 1: the entries of nodes_[i]
  std::vector<NodeId> parents_ {};
  std::vector<EdgeId> via_edges_ {};
  std::vector<EdgeId> edges_ {};  // via_ in ascending EdgeId (derived once, at construction)
};

// An explicit walk of hops; a hop holds one or more parallel members of one (u, v) pair, in
// the caller's order. source(), destination() and cost() (the sum over hops of the minimum member
// cost) are derived at construction. A single-member walk is [[e1], [e2], ...].
class Path {
public:
  // Validates, in this order (std::invalid_argument on the first defect): every member id in
  // range (before any endpoint is read), no duplicate member, a non-empty hop list, every hop
  // non-empty, the members of a hop sharing (u, v), contiguous hops, no repeated node.
  Path(const Graph& graph, std::vector<std::vector<EdgeId>> hops);

  [[nodiscard]] std::uint64_t graph_uid() const noexcept { return graph_uid_; }
  [[nodiscard]] const std::vector<std::vector<EdgeId>>& hops() const noexcept { return hops_; }
  [[nodiscard]] NodeId source() const noexcept { return nodes_.front(); }
  [[nodiscard]] NodeId destination() const noexcept { return nodes_.back(); }
  [[nodiscard]] Cost cost() const noexcept { return cost_; }
  // The equivalent PredecessorDAG (every member of every hop; entries in ascending via EdgeId).
  [[nodiscard]] PredecessorDAG to_dag() const;

  friend bool operator==(const Path& a, const Path& b) noexcept {
    return a.graph_uid_ == b.graph_uid_ && a.hops_ == b.hops_;
  }

private:
  friend struct detail::PathBuilder;
  Path() = default;
  std::uint64_t graph_uid_ {0};
  std::int32_t num_nodes_ {0};
  std::vector<std::vector<EdgeId>> hops_ {};
  std::vector<NodeId> nodes_ {};
  Cost cost_ {0};
};

} // namespace netgraph::core
