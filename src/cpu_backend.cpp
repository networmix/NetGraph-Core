/*
  CPU ExecutionBackend — thin adapter that delegates to in-process algorithms.
*/
#include "netgraph/core/backend.hpp"
#include "netgraph/core/k_shortest_paths.hpp"

namespace netgraph::core {

namespace {
class CpuBackend final : public ExecutionBackend {
public:
  std::pair<std::vector<Cost>, PredDAG> shortest_paths(
      const StrictMultiDiGraph& g, NodeId src, std::optional<NodeId> dst,
      const EdgeSelection& selection,
      std::span<const Cap> residual,
      const bool* node_mask,
      const bool* edge_mask) override {
    return netgraph::core::shortest_paths(g, src, dst, selection, residual, node_mask, edge_mask);
  }

  std::pair<Flow, FlowSummary> calc_max_flow(
      const StrictMultiDiGraph& g, NodeId src, NodeId dst,
      FlowPlacement placement, bool shortest_path,
      bool with_edge_flows,
      bool with_reachable,
      bool with_residuals,
      const bool* node_mask = nullptr,
      const bool* edge_mask = nullptr) override {
    return netgraph::core::calc_max_flow(
        g, src, dst,
        placement, shortest_path,
        with_edge_flows,
        with_reachable,
        with_residuals,
        node_mask, edge_mask);
  }

  std::vector<std::pair<std::vector<Cost>, PredDAG>> k_shortest_paths(
      const StrictMultiDiGraph& g, NodeId src, NodeId dst,
      int k, std::optional<double> max_cost_factor,
      bool unique,
      const bool* node_mask,
      const bool* edge_mask) override {
    return netgraph::core::k_shortest_paths(g, src, dst, k, max_cost_factor, unique, node_mask, edge_mask);
  }
};
} // namespace

std::unique_ptr<ExecutionBackend> make_cpu_backend() {
  return std::make_unique<CpuBackend>();
}

} // namespace netgraph::core
