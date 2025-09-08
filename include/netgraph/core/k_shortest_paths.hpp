/* Yen-like K-shortest paths with SPF-compatible outputs per-path. */
#pragma once

#include <optional>
#include <vector>

#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// Compute up to k shortest paths from s to t (Yen-like) and return
// SPF-compatible outputs per path: (distances, predecessor DAG).
// Distances are Cost[N] (int64); unreachable entries are set to
// std::numeric_limits<Cost>::max(). PredDAG encodes one concrete path with a
// single parent per node along the path; other nodes have no parents.
// Deterministic tie-breaking across equal-cost edges uses compacted edge order.
[[nodiscard]] std::vector<std::pair<std::vector<Cost>, PredDAG>> k_shortest_paths(
    const StrictMultiDiGraph& g, NodeId src, NodeId dst,
    int k, std::optional<double> max_cost_factor,
    bool unique,
    std::span<const bool> node_mask = {},
    std::span<const bool> edge_mask = {});

} // namespace netgraph::core
