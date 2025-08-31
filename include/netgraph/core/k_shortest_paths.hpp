#pragma once

#include <optional>
#include <vector>

#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

struct Path {
  std::vector<NodeId> nodes;
  std::vector<EdgeId> edges;
  double cost {0.0};
};

std::vector<Path> k_shortest_paths(
    const StrictMultiDiGraph& g, NodeId s, NodeId t,
    int k, std::optional<double> max_cost_factor,
    bool unique, double eps);

} // namespace netgraph::core
