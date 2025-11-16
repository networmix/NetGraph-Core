/* Option structs for algorithm invocations to keep interfaces concise. */
#pragma once

#include <optional>
#include <span>

#include "netgraph/core/types.hpp"

namespace netgraph::core {

struct SpfOptions {
  bool multipath { true };
  EdgeSelection selection {};
  std::optional<NodeId> dst {};
  std::span<const Cap> residual {};
  std::span<const bool> node_mask {};
  std::span<const bool> edge_mask {};
};

struct KspOptions {
  int k { 1 };
  std::optional<double> max_cost_factor {};
  bool unique { true };
  std::span<const bool> node_mask {};
  std::span<const bool> edge_mask {};
};

struct MaxFlowOptions {
  FlowPlacement placement { FlowPlacement::Proportional };
  bool shortest_path { false };
  bool require_capacity { true };  // Require edges to have capacity (SDN/TE). Set false for cost-only routing (IP/IGP)
  bool with_edge_flows { false };
  bool with_reachable { false };
  bool with_residuals { false };
  std::span<const bool> node_mask {};
  std::span<const bool> edge_mask {};
};

} // namespace netgraph::core
