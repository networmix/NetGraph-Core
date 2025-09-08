/*
  CPU Backend — thin adapter that delegates to in-process algorithms.
*/
#include "netgraph/core/backend.hpp"
#include "netgraph/core/k_shortest_paths.hpp"
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include <cassert>

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
    // Precondition checks (debug mode only)
    assert(src >= 0 && src < g.num_nodes());
    if (opts.dst.has_value()) {
      assert(*opts.dst >= 0 && *opts.dst < g.num_nodes());
    }
    if (!opts.node_mask.empty()) {
      assert(opts.node_mask.size() == static_cast<size_t>(g.num_nodes()));
    }
    if (!opts.edge_mask.empty()) {
      assert(opts.edge_mask.size() == static_cast<size_t>(g.num_edges()));
    }
    if (!opts.residual.empty()) {
      assert(opts.residual.size() == static_cast<size_t>(g.num_edges()));
    }
    // Convert spans to raw pointers for compatibility with existing free functions
    const bool* node_ptr = opts.node_mask.empty() ? nullptr : opts.node_mask.data();
    const bool* edge_ptr = opts.edge_mask.empty() ? nullptr : opts.edge_mask.data();
    return netgraph::core::shortest_paths(g, src, opts.dst, opts.multipath, opts.selection,
                                          opts.residual, node_ptr, edge_ptr);
  }

  std::pair<Flow, FlowSummary> max_flow(
      const GraphHandle& gh, NodeId src, NodeId dst, const MaxFlowOptions& opts) override {
    const StrictMultiDiGraph& g = *gh.graph;
    // Precondition checks (debug mode only)
    assert(src >= 0 && src < g.num_nodes());
    assert(dst >= 0 && dst < g.num_nodes());
    if (!opts.node_mask.empty()) {
      assert(opts.node_mask.size() == static_cast<size_t>(g.num_nodes()));
    }
    if (!opts.edge_mask.empty()) {
      assert(opts.edge_mask.size() == static_cast<size_t>(g.num_edges()));
    }
    // Convert spans to raw pointers for compatibility with existing free functions
    const bool* node_ptr = opts.node_mask.empty() ? nullptr : opts.node_mask.data();
    const bool* edge_ptr = opts.edge_mask.empty() ? nullptr : opts.edge_mask.data();
    return netgraph::core::calc_max_flow(
        g, src, dst,
        opts.placement, opts.shortest_path,
        opts.with_edge_flows,
        opts.with_reachable,
        opts.with_residuals,
        node_ptr, edge_ptr);
  }

  std::vector<std::pair<std::vector<Cost>, PredDAG>> ksp(
      const GraphHandle& gh, NodeId src, NodeId dst, const KspOptions& opts) override {
    const StrictMultiDiGraph& g = *gh.graph;
    // Precondition checks (debug mode only)
    assert(src >= 0 && src < g.num_nodes());
    assert(dst >= 0 && dst < g.num_nodes());
    assert(opts.k > 0);
    if (!opts.node_mask.empty()) {
      assert(opts.node_mask.size() == static_cast<size_t>(g.num_nodes()));
    }
    if (!opts.edge_mask.empty()) {
      assert(opts.edge_mask.size() == static_cast<size_t>(g.num_edges()));
    }
    // Convert spans to raw pointers for compatibility with existing free functions
    const bool* node_ptr = opts.node_mask.empty() ? nullptr : opts.node_mask.data();
    const bool* edge_ptr = opts.edge_mask.empty() ? nullptr : opts.edge_mask.data();
    return netgraph::core::k_shortest_paths(g, src, dst, opts.k, opts.max_cost_factor,
                                            opts.unique, node_ptr, edge_ptr);
  }

  std::vector<FlowSummary> batch_max_flow(
      const GraphHandle& gh,
      const std::vector<std::pair<NodeId,NodeId>>& pairs,
      const MaxFlowOptions& opts,
      const std::vector<std::span<const bool>>& node_masks,
      const std::vector<std::span<const bool>>& edge_masks) override {
    const StrictMultiDiGraph& g = *gh.graph;
    // Precondition checks (debug mode only)
    for (const auto& [src, dst] : pairs) {
      assert(src >= 0 && src < g.num_nodes());
      assert(dst >= 0 && dst < g.num_nodes());
    }
    assert(node_masks.size() == edge_masks.size() || node_masks.empty() || edge_masks.empty());
    for (const auto& span : node_masks) {
      assert(span.size() == static_cast<size_t>(g.num_nodes()));
    }
    for (const auto& span : edge_masks) {
      assert(span.size() == static_cast<size_t>(g.num_edges()));
    }
    // Convert spans to raw pointers for compatibility with existing free functions
    std::vector<const bool*> node_ptrs, edge_ptrs;
    node_ptrs.reserve(node_masks.size());
    edge_ptrs.reserve(edge_masks.size());
    for (const auto& span : node_masks) node_ptrs.push_back(span.data());
    for (const auto& span : edge_masks) edge_ptrs.push_back(span.data());
    return netgraph::core::batch_max_flow(g, pairs,
                                          opts.placement, opts.shortest_path,
                                          opts.with_edge_flows, opts.with_reachable, opts.with_residuals,
                                          node_ptrs, edge_ptrs);
  }
};
} // namespace

BackendPtr make_cpu_backend() {
  return std::make_shared<CpuBackend>();
}

} // namespace netgraph::core
