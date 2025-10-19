/* Shortest paths (Dijkstra) with multipath predecessor DAG support.
 *
 * For Python developers:
 * - std::pair<A, B>: tuple of two elements (like tuple[A, B])
 * - std::optional<T>: nullable value (like T | None)
 * - std::nullopt: None equivalent for std::optional
 */
#pragma once

#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// PredDAG (Predecessor Directed Acyclic Graph): compact representation of all equal-cost
// shortest paths from a source node. Similar to a defaultdict(list) in Python, but stored
// in CSR format for efficiency.
//
// For each node v, predecessors are stored in parents[parent_offsets[v]:parent_offsets[v+1]]
// with corresponding EdgeIds in via_edges[parent_offsets[v]:parent_offsets[v+1]].
// Multiple parallel edges are represented by multiple entries with the same parent.
// parent_offsets has length N+1 (like CSR row pointers).
struct PredDAG {
  std::vector<std::int32_t> parent_offsets;  // Length N+1 (CSR row pointers)
  std::vector<NodeId> parents;                // Predecessor node IDs
  std::vector<EdgeId> via_edges;              // EdgeId used to reach node from predecessor
};

// Compute shortest paths from src using Dijkstra's algorithm.
// Returns (distances, predecessor_dag) where distances[v] is the shortest cost to reach v
// (or inf if unreachable), and predecessor_dag encodes all equal-cost paths.
//
// Parameters:
// - dst: if provided, algorithm may exit early once destination is reached
// - multipath: if true, keep all equal-cost predecessors; if false, keep only one per node
// - selection: edge selection policy (multi-edge, capacity filtering, tie-breaking)
// - residual: if provided, use these capacities instead of graph's original capacities
// - node_mask: if provided, node_mask[v]==true means node v is allowed (false excludes it)
// - edge_mask: if provided, edge_mask[e]==true means edge e is allowed (false excludes it)
[[nodiscard]] std::pair<std::vector<Cost>, PredDAG>
shortest_paths(const StrictMultiDiGraph& g, NodeId src,
               std::optional<NodeId> dst,
               bool multipath,
               const EdgeSelection& selection,
               std::span<const Cap> residual = {},
               std::span<const bool> node_mask = {},
               std::span<const bool> edge_mask = {});

// Enumerate concrete paths represented by a PredDAG from src to dst.
// Each path is returned as a sequence of (node_id, (edge_ids...)) pairs ending with (dst, ()).
// When split_parallel_edges=false, parallel edges per hop are grouped in the tuple.
// When true, one edge per hop is selected to produce concrete paths; enumeration may be capped with max_paths.
[[nodiscard]] std::vector<std::vector<std::pair<NodeId, std::vector<EdgeId>>>>
resolve_to_paths(const PredDAG& dag, NodeId src, NodeId dst,
                 bool split_parallel_edges = false,
                 std::optional<std::int64_t> max_paths = std::nullopt);

} // namespace netgraph::core
