#pragma once

#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"

namespace netgraph::core {

class ExecutionBackend {
public:
  virtual ~ExecutionBackend() = default;
  virtual std::pair<std::vector<double>, PredDAG> shortest_paths(
      const StrictMultiDiGraph& g, NodeId src, std::optional<NodeId> dst,
      EdgeSelect policy, bool multipath, double eps) = 0;

  virtual std::pair<double, FlowSummary> calc_max_flow(
      const StrictMultiDiGraph& g, NodeId s, NodeId t,
      FlowPlacement placement, bool shortest_path,
      double eps, bool with_edge_flows,
      const bool* node_mask = nullptr,
      const bool* edge_mask = nullptr) = 0;
};

std::unique_ptr<ExecutionBackend> make_cpu_backend();

} // namespace netgraph::core
