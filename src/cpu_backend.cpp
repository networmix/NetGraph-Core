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
    if (!opts.node_mask.empty() && opts.node_mask.size() != static_cast<std::size_t>(g.num_nodes())) {
      throw std::invalid_argument("CpuBackend::spf: node_mask length mismatch");
    }
    if (!opts.edge_mask.empty() && opts.edge_mask.size() != static_cast<std::size_t>(g.num_edges())) {
      throw std::invalid_argument("CpuBackend::spf: edge_mask length mismatch");
    }
    // Forward spans directly
    std::span<const bool> node_span = (opts.node_mask.size() == static_cast<size_t>(g.num_nodes())) ? opts.node_mask : std::span<const bool>{};
    std::span<const bool> edge_span = (opts.edge_mask.size() == static_cast<size_t>(g.num_edges())) ? opts.edge_mask : std::span<const bool>{};
    return netgraph::core::shortest_paths(g, src, opts.dst, opts.multipath, opts.selection,
                                          opts.residual, node_span, edge_span);
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
    // Forward spans directly
    std::span<const bool> node_span = (opts.node_mask.size() == static_cast<size_t>(g.num_nodes())) ? opts.node_mask : std::span<const bool>{};
    std::span<const bool> edge_span = (opts.edge_mask.size() == static_cast<size_t>(g.num_edges())) ? opts.edge_mask : std::span<const bool>{};
    return netgraph::core::calc_max_flow(
        g, src, dst,
        opts.placement, opts.shortest_path,
        opts.require_capacity,
        opts.with_edge_flows,
        opts.with_reachable,
        opts.with_residuals,
        node_span, edge_span);
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
    // Forward spans directly
    std::span<const bool> node_span = (opts.node_mask.size() == static_cast<size_t>(g.num_nodes())) ? opts.node_mask : std::span<const bool>{};
    std::span<const bool> edge_span = (opts.edge_mask.size() == static_cast<size_t>(g.num_edges())) ? opts.edge_mask : std::span<const bool>{};
    return netgraph::core::k_shortest_paths(g, src, dst, opts.k, opts.max_cost_factor,
                                            opts.unique, node_span, edge_span);
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
    // Forward spans directly
    std::vector<std::span<const bool>> node_ptrs, edge_ptrs;
    node_ptrs.reserve(node_masks.size());
    edge_ptrs.reserve(edge_masks.size());
    for (const auto& span : node_masks) {
      node_ptrs.push_back((span.size() == static_cast<size_t>(g.num_nodes())) ? span : std::span<const bool>{});
    }
    for (const auto& span : edge_masks) {
      edge_ptrs.push_back((span.size() == static_cast<size_t>(g.num_edges())) ? span : std::span<const bool>{});
    }
    return netgraph::core::batch_max_flow(g, pairs,
                                          opts.placement, opts.shortest_path,
                                          opts.require_capacity,
                                          opts.with_edge_flows, opts.with_reachable, opts.with_residuals,
                                          node_ptrs, edge_ptrs);
  }
};
} // namespace

BackendPtr make_cpu_backend() {
  return std::make_shared<CpuBackend>();
}

} // namespace netgraph::core
