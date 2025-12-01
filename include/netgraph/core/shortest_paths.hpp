/* Shortest paths (Dijkstra) with multipath predecessor DAG support. */
#pragma once

#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// PredDAG (Predecessor Directed Acyclic Graph): compact representation of all equal-cost
// shortest paths from a source node. Stored in CSR format for efficiency.
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
// (or inf if unreachable), and predecessor_dag encodes the shortest path structure.
//
// Parameters:
// - dst: if provided, algorithm may exit early once destination is reached
// - multipath: if true, keep all equal-cost predecessors; if false, keep one per node
//   (in single-path mode, ties are broken by preferring higher bottleneck capacity)
// - selection: edge selection policy (multi-edge, capacity filtering, tie-breaking)
// - residual: if provided, use these capacities instead of graph's original capacities.
//   Passing a non-empty 'residual' span forces capacity gating:
//   edges with residual < kMinCap are excluded from exploration.
//   For EB semantics this means re-running SPF on updated residuals
//   will remove saturated next-hops and therefore *change* the fixed
//   equal-split set (progressive behavior).
// - node_mask: if provided, node_mask[v]==true means node v is allowed (false excludes it).
//   If the source node is masked (node_mask[src]==false), returns an empty predecessor DAG
//   with all distances at infinity, as no traversal can begin from an excluded source.
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
