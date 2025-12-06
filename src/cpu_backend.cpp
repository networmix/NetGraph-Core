/*
  CPU Backend — thin adapter that delegates to in-process algorithms.
*/
#include "netgraph/core/backend.hpp"
#include "netgraph/core/k_shortest_paths.hpp"
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/shortest_paths.hpp"

namespace netgraph::core {

namespace {
class CpuBackend final : public Backend {
public:
  GraphHandle build_graph(const StrictMultiDiGraph& g) override {
    // Create a non-owning shared_ptr with no-op deleter; lifetime is managed by caller
    return GraphHandle{ std::shared_ptr<const StrictMultiDiGraph>(&g, [](const StrictMultiDiGraph*){}) };
  }

  GraphHandle build_graph(std::shared_ptr<const StrictMultiDiGraph> g) override {
    return GraphHandle{ std::move(g) };
  }

  std::pair<std::vector<Cost>, PredDAG> spf(
      const GraphHandle& gh, NodeId src, const SpfOptions& opts) override {
    const StrictMultiDiGraph& g = *gh.graph;
    // Validate mask lengths strictly; mismatches are user errors.
    // NOTE: Keep this as the single public boundary check; deeper layers
    // (shortest_paths) should rely on this to avoid redundant validation.
    if (!opts.node_mask.empty() && opts.node_mask.size() != static_cast<std::size_t>(g.num_nodes())) {
      throw std::invalid_argument("CpuBackend::spf: node_mask length mismatch");
    }
    if (!opts.edge_mask.empty() && opts.edge_mask.size() != static_cast<std::size_t>(g.num_edges())) {
      throw std::invalid_argument("CpuBackend::spf: edge_mask length mismatch");
    }
    // Forward spans directly (already validated above)
    return netgraph::core::shortest_paths(g, src, opts.dst, opts.multipath, opts.selection,
                                          opts.residual, opts.node_mask, opts.edge_mask);
  }

  std::pair<Flow, FlowSummary> max_flow(
      const GraphHandle& gh, NodeId src, NodeId dst, const MaxFlowOptions& opts) override {
    const StrictMultiDiGraph& g = *gh.graph;
    // Validate mask lengths strictly.
    if (!opts.node_mask.empty() && opts.node_mask.size() != static_cast<std::size_t>(g.num_nodes())) {
      throw std::invalid_argument("CpuBackend::max_flow: node_mask length mismatch");
    }
    if (!opts.edge_mask.empty() && opts.edge_mask.size() != static_cast<std::size_t>(g.num_edges())) {
      throw std::invalid_argument("CpuBackend::max_flow: edge_mask length mismatch");
    }
    // Forward spans directly (already validated above)
    return netgraph::core::calc_max_flow(
        g, src, dst,
        opts.placement, opts.shortest_path,
        opts.require_capacity,
        opts.with_edge_flows,
        opts.with_reachable,
        opts.with_residuals,
        opts.node_mask, opts.edge_mask);
  }

  std::vector<std::pair<std::vector<Cost>, PredDAG>> ksp(
      const GraphHandle& gh, NodeId src, NodeId dst, const KspOptions& opts) override {
    const StrictMultiDiGraph& g = *gh.graph;
    if (opts.k <= 0) { return {}; }
    // Validate mask lengths strictly.
    if (!opts.node_mask.empty() && opts.node_mask.size() != static_cast<std::size_t>(g.num_nodes())) {
      throw std::invalid_argument("CpuBackend::ksp: node_mask length mismatch");
    }
    if (!opts.edge_mask.empty() && opts.edge_mask.size() != static_cast<std::size_t>(g.num_edges())) {
      throw std::invalid_argument("CpuBackend::ksp: edge_mask length mismatch");
    }
    // Forward spans directly (already validated above)
    return netgraph::core::k_shortest_paths(g, src, dst, opts.k, opts.max_cost_factor,
                                            opts.unique, opts.node_mask, opts.edge_mask);
  }

  std::vector<FlowSummary> batch_max_flow(
      const GraphHandle& gh,
      const std::vector<std::pair<NodeId,NodeId>>& pairs,
      const MaxFlowOptions& opts,
      const std::vector<std::span<const bool>>& node_masks,
      const std::vector<std::span<const bool>>& edge_masks) override {
    const StrictMultiDiGraph& g = *gh.graph;
    // Validate batch mask lengths strictly.
    for (const auto& span : node_masks) {
      if (!span.empty() && span.size() != static_cast<std::size_t>(g.num_nodes())) {
        throw std::invalid_argument("CpuBackend::batch_max_flow: node_mask length mismatch");
      }
    }
    for (const auto& span : edge_masks) {
      if (!span.empty() && span.size() != static_cast<std::size_t>(g.num_edges())) {
        throw std::invalid_argument("CpuBackend::batch_max_flow: edge_mask length mismatch");
      }
    }
    // Forward spans directly (already validated above)
    return netgraph::core::batch_max_flow(g, pairs,
                                          opts.placement, opts.shortest_path,
                                          opts.require_capacity,
                                          opts.with_edge_flows, opts.with_reachable, opts.with_residuals,
                                          node_masks, edge_masks);
  }

  std::vector<std::pair<EdgeId, Flow>> sensitivity_analysis(
      const GraphHandle& gh, NodeId src, NodeId dst, const MaxFlowOptions& opts) override {
    const StrictMultiDiGraph& g = *gh.graph;
    // Validate mask lengths strictly.
    if (!opts.node_mask.empty() && opts.node_mask.size() != static_cast<std::size_t>(g.num_nodes())) {
      throw std::invalid_argument("CpuBackend::sensitivity_analysis: node_mask length mismatch");
    }
    if (!opts.edge_mask.empty() && opts.edge_mask.size() != static_cast<std::size_t>(g.num_edges())) {
      throw std::invalid_argument("CpuBackend::sensitivity_analysis: edge_mask length mismatch");
    }
    return netgraph::core::sensitivity_analysis(g, src, dst,
                                                opts.placement, opts.shortest_path,
                                                opts.require_capacity,
                                                opts.node_mask, opts.edge_mask);
  }
};
} // namespace

BackendPtr make_cpu_backend() {
  return std::make_shared<CpuBackend>();
}

} // namespace netgraph::core
