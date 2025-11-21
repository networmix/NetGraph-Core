/*
  Backend interface — abstracts SPF/MaxFlow/KSP implementations.

  The default CPU backend delegates to in-process algorithm implementations.
  All execution flows through this interface via an Algorithms façade.

  For Python developers:
  - std::shared_ptr<T>: reference-counted pointer (like Python object references)
  - virtual: method can be overridden in subclasses (like Python's inheritance)
  - = 0: pure virtual (must be implemented by subclass, like @abstractmethod)
*/
#pragma once

#include <memory>
#include <optional>
#include <utility>
#include <vector>
#include <span>

#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/options.hpp"

namespace netgraph::core {

// GraphHandle: opaque handle to a backend-owned graph.
// Uses shared_ptr for automatic lifetime management (like Python's reference counting).
struct GraphHandle {
  std::shared_ptr<const StrictMultiDiGraph> graph {};
};

class Backend {
public:
  virtual ~Backend() noexcept = default;

  // Prepare a backend-specific graph handle from an existing graph reference.
  // CPU backend creates a non-owning shared_ptr with a no-op deleter.
  [[nodiscard]] virtual GraphHandle build_graph(const StrictMultiDiGraph& g) = 0;

  // Prepare a backend-specific graph handle that takes shared ownership of the
  // provided graph instance.
  [[nodiscard]] virtual GraphHandle build_graph(std::shared_ptr<const StrictMultiDiGraph> g) = 0;

  [[nodiscard]] virtual std::pair<std::vector<Cost>, PredDAG> spf(
      const GraphHandle& gh, NodeId src, const SpfOptions& opts) = 0;

  [[nodiscard]] virtual std::pair<Flow, FlowSummary> max_flow(
      const GraphHandle& gh, NodeId src, NodeId dst, const MaxFlowOptions& opts) = 0;

  [[nodiscard]] virtual std::vector<std::pair<std::vector<Cost>, PredDAG>> ksp(
      const GraphHandle& gh, NodeId src, NodeId dst, const KspOptions& opts) = 0;

  [[nodiscard]] virtual std::vector<FlowSummary> batch_max_flow(
      const GraphHandle& gh,
      const std::vector<std::pair<NodeId,NodeId>>& pairs,
      const MaxFlowOptions& opts,
      const std::vector<std::span<const bool>>& node_masks = {},
      const std::vector<std::span<const bool>>& edge_masks = {}) = 0;

  [[nodiscard]] virtual std::vector<std::pair<EdgeId, Flow>> sensitivity_analysis(
      const GraphHandle& gh, NodeId src, NodeId dst, const MaxFlowOptions& opts) = 0;
};

using BackendPtr = std::shared_ptr<Backend>;

[[nodiscard]] BackendPtr make_cpu_backend();

} // namespace netgraph::core
