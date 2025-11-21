/*
  Backend interface — abstracts SPF/MaxFlow/KSP implementations.

  The default CPU backend delegates to in-process algorithm implementations.
  All execution flows through this interface via an Algorithms façade.
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
// Uses shared_ptr for automatic lifetime management.
struct GraphHandle {
  std::shared_ptr<const StrictMultiDiGraph> graph {};
};

class Backend {
public:
  virtual ~Backend() noexcept = default;

  // Prepares a backend-specific graph handle from an existing graph reference.
  //
  // The CPU backend creates a non-owning shared_ptr with a no-op deleter.
  //
  // Arguments:
  //   g: The source graph to wrap.
  //
  // Returns:
  //   A GraphHandle containing the backend-specific graph representation.
  [[nodiscard]] virtual GraphHandle build_graph(const StrictMultiDiGraph& g) = 0;

  // Prepares a backend-specific graph handle that takes shared ownership of the
  // provided graph instance.
  //
  // Arguments:
  //   g: The source graph as a shared_ptr.
  //
  // Returns:
  //   A GraphHandle that shares ownership of the graph.
  [[nodiscard]] virtual GraphHandle build_graph(std::shared_ptr<const StrictMultiDiGraph> g) = 0;

  // Computes shortest paths from a source node.
  //
  // Arguments:
  //   gh: The graph handle.
  //   src: The source node ID.
  //   opts: Configuration options for SPF (e.g., multipath, edge selection).
  //
  // Returns:
  //   A pair containing:
  //     - A vector of costs (distances) from src to all nodes.
  //     - A predecessor DAG encoding the shortest paths.
  [[nodiscard]] virtual std::pair<std::vector<Cost>, PredDAG> spf(
      const GraphHandle& gh, NodeId src, const SpfOptions& opts) = 0;

  // Computes maximum flow between a source and destination node.
  //
  // Arguments:
  //   gh: The graph handle.
  //   src: The source node ID.
  //   dst: The destination node ID.
  //   opts: Configuration options for max flow (e.g., algorithm, capacity constraints).
  //
  // Returns:
  //   A pair containing:
  //     - The total flow amount.
  //     - A FlowSummary struct with detailed results.
  [[nodiscard]] virtual std::pair<Flow, FlowSummary> max_flow(
      const GraphHandle& gh, NodeId src, NodeId dst, const MaxFlowOptions& opts) = 0;

  // Computes K-shortest paths between a source and destination node.
  //
  // Arguments:
  //   gh: The graph handle.
  //   src: The source node ID.
  //   dst: The destination node ID.
  //   opts: Configuration options for KSP (e.g., K, cost constraints).
  //
  // Returns:
  //   A vector of pairs, each representing a path:
  //     - A vector of costs along the path.
  //     - A predecessor DAG representing the path.
  [[nodiscard]] virtual std::vector<std::pair<std::vector<Cost>, PredDAG>> ksp(
      const GraphHandle& gh, NodeId src, NodeId dst, const KspOptions& opts) = 0;

  // Computes maximum flow for a batch of source-destination pairs.
  //
  // Arguments:
  //   gh: The graph handle.
  //   pairs: A vector of (source, destination) pairs.
  //   opts: Common configuration options for all pairs.
  //   node_masks: Optional vector of node masks, one per pair.
  //   edge_masks: Optional vector of edge masks, one per pair.
  //
  // Returns:
  //   A vector of FlowSummary structs, one for each pair.
  [[nodiscard]] virtual std::vector<FlowSummary> batch_max_flow(
      const GraphHandle& gh,
      const std::vector<std::pair<NodeId,NodeId>>& pairs,
      const MaxFlowOptions& opts,
      const std::vector<std::span<const bool>>& node_masks = {},
      const std::vector<std::span<const bool>>& edge_masks = {}) = 0;

  // Performs sensitivity analysis to identify edges that constrain flow.
  //
  // Arguments:
  //   gh: The graph handle.
  //   src: The source node ID.
  //   dst: The destination node ID.
  //   opts: Configuration options.
  //
  // Returns:
  //   A vector of pairs (EdgeId, Flow gain), indicating how much flow would
  //   increase if the edge's capacity were relaxed.
  [[nodiscard]] virtual std::vector<std::pair<EdgeId, Flow>> sensitivity_analysis(
      const GraphHandle& gh, NodeId src, NodeId dst, const MaxFlowOptions& opts) = 0;
};

using BackendPtr = std::shared_ptr<Backend>;

[[nodiscard]] BackendPtr make_cpu_backend();

} // namespace netgraph::core
