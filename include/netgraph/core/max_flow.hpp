#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

struct MinCut { std::vector<EdgeId> edges; };
// CostBucket: 'share' holds the raw flow amount placed at the given total path cost.
// This mirrors NetGraph's cost_distribution semantics (not a normalized ratio).
struct CostBucket { double cost; double share; };
struct CostDistribution { std::vector<CostBucket> buckets; };

struct FlowSummary {
  double total_flow {0.0};
  MinCut min_cut {};
  CostDistribution cost_distribution {};
  std::vector<double> edge_flows; // filled if requested
};

std::pair<double, FlowSummary>
calc_max_flow(const StrictMultiDiGraph& g, NodeId s, NodeId t,
              FlowPlacement placement, bool shortest_path,
              double eps, bool with_edge_flows,
              const bool* node_mask = nullptr,
              const bool* edge_mask = nullptr);

std::vector<FlowSummary>
batch_max_flow(const StrictMultiDiGraph& g,
               const std::vector<std::pair<NodeId,NodeId>>& pairs,
               FlowPlacement placement, bool shortest_path,
               double eps, bool with_edge_flows,
               int threads, std::optional<std::uint64_t> seed,
               const std::vector<const bool*>& node_masks = {},
               const std::vector<const bool*>& edge_masks = {});

} // namespace netgraph::core
