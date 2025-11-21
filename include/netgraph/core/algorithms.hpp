/* Algorithms façade that forwards to a selected Backend implementation. */
#pragma once

#include <memory>

#include "netgraph/core/backend.hpp"
#include "netgraph/core/options.hpp"

namespace netgraph::core {

class Algorithms {
public:
  explicit Algorithms(BackendPtr backend) : backend_(std::move(backend)) {}

  [[nodiscard]] GraphHandle build_graph(const StrictMultiDiGraph& g) const {
    return backend_->build_graph(g);
  }

  [[nodiscard]] GraphHandle build_graph(std::shared_ptr<const StrictMultiDiGraph> g) const {
    return backend_->build_graph(std::move(g));
  }

  [[nodiscard]] std::pair<std::vector<Cost>, PredDAG>
  spf(const GraphHandle& gh, NodeId src, const SpfOptions& opts) const {
    return backend_->spf(gh, src, opts);
  }

  [[nodiscard]] std::vector<std::pair<std::vector<Cost>, PredDAG>>
  ksp(const GraphHandle& gh, NodeId src, NodeId dst, const KspOptions& opts) const {
    return backend_->ksp(gh, src, dst, opts);
  }

  [[nodiscard]] std::pair<Flow, FlowSummary>
  max_flow(const GraphHandle& gh, NodeId src, NodeId dst, const MaxFlowOptions& opts) const {
    return backend_->max_flow(gh, src, dst, opts);
  }

  [[nodiscard]] std::vector<FlowSummary>
  batch_max_flow(const GraphHandle& gh,
                 const std::vector<std::pair<NodeId,NodeId>>& pairs,
                 const MaxFlowOptions& opts,
                 const std::vector<std::span<const bool>>& node_masks = {},
                 const std::vector<std::span<const bool>>& edge_masks = {}) const {
    return backend_->batch_max_flow(gh, pairs, opts, node_masks, edge_masks);
  }

  [[nodiscard]] std::vector<std::pair<EdgeId, Flow>>
  sensitivity_analysis(const GraphHandle& gh, NodeId src, NodeId dst, const MaxFlowOptions& opts) const {
    return backend_->sensitivity_analysis(gh, src, dst, opts);
  }

private:
  BackendPtr backend_;
};

} // namespace netgraph::core
