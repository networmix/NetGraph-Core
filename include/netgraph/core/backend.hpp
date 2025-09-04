/*
  ExecutionBackend interface — abstracts SP/MaxFlow/KSP implementations.

  The default CPU backend delegates to in-process algorithm implementations.
*/
#pragma once

#include <memory>
#include <optional>
#include <utility>
#include <vector>
#include <span>

#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"

namespace netgraph::core {

class ExecutionBackend {
public:
  virtual ~ExecutionBackend() = default;
  virtual std::pair<std::vector<Cost>, PredDAG> shortest_paths(
      const StrictMultiDiGraph& g, NodeId src, std::optional<NodeId> dst,
      const EdgeSelection& selection,
      std::span<const Cap> residual = {},
      const bool* node_mask = nullptr,
      const bool* edge_mask = nullptr) = 0;

  virtual std::pair<Flow, FlowSummary> calc_max_flow(
      const StrictMultiDiGraph& g, NodeId src, NodeId dst,
      FlowPlacement placement, bool shortest_path,
      bool with_edge_flows,
      const bool* node_mask = nullptr,
      const bool* edge_mask = nullptr) = 0;

  virtual std::vector<std::pair<std::vector<Cost>, PredDAG>> k_shortest_paths(
      const StrictMultiDiGraph& g, NodeId s, NodeId t,
      int k, std::optional<double> max_cost_factor,
      bool unique,
      const bool* node_mask = nullptr,
      const bool* edge_mask = nullptr) = 0;
};

std::unique_ptr<ExecutionBackend> make_cpu_backend();

} // namespace netgraph::core
