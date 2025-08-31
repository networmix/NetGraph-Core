#include "netgraph/core/backend.hpp"

namespace netgraph::core {

namespace {
class CpuBackend final : public ExecutionBackend {
public:
  std::pair<std::vector<double>, PredDAG> shortest_paths(
      const StrictMultiDiGraph& g, NodeId src, std::optional<NodeId> dst,
      EdgeSelect policy, bool multipath, double eps) override {
    return netgraph::core::shortest_paths(g, src, dst, policy, multipath, eps);
  }

  std::pair<double, FlowSummary> calc_max_flow(
      const StrictMultiDiGraph& g, NodeId s, NodeId t,
      FlowPlacement placement, bool shortest_path,
      double eps, bool with_edge_flows,
      const bool* node_mask = nullptr,
      const bool* edge_mask = nullptr) override {
    return netgraph::core::calc_max_flow(g, s, t, placement, shortest_path, eps, with_edge_flows, node_mask, edge_mask);
  }
};
} // namespace

std::unique_ptr<ExecutionBackend> make_cpu_backend() {
  return std::make_unique<CpuBackend>();
}

} // namespace netgraph::core
