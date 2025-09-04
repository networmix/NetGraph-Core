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
// Distances are float64[N], PredDAG encodes one concrete path with single parent
// per node along the path; other nodes have no parents and dist=inf.
// Deterministic tie-breaking across equal-cost edges uses compacted edge order.
std::vector<std::pair<std::vector<Cost>, PredDAG>> k_shortest_paths(
    const StrictMultiDiGraph& g, NodeId src, NodeId dst,
    int k, std::optional<double> max_cost_factor,
    bool unique,
    const bool* node_mask = nullptr,
    const bool* edge_mask = nullptr);

} // namespace netgraph::core
