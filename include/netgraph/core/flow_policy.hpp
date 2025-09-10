/* FlowPolicy manages flows for a single demand. */
#pragma once

#include <cstdint>
#include <limits>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>
#include <stdexcept>

#include "netgraph/core/flow.hpp"
#include "netgraph/core/flow_graph.hpp"
#include "netgraph/core/algorithms.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

enum class PathAlg : std::int32_t { SPF = 1 };

// Execution context with algorithms and graph handle
struct ExecutionContext {
  Algorithms* algorithms;
  GraphHandle graph;

  // Constructor with validation
  ExecutionContext(Algorithms& algs, const GraphHandle& gh) noexcept
      : algorithms(&algs), graph(gh) {}
};

// Configuration struct for FlowPolicy behavior
struct FlowPolicyConfig {
  PathAlg path_alg { PathAlg::SPF };
  FlowPlacement flow_placement { FlowPlacement::Proportional };
  EdgeSelection selection { EdgeSelection{} };
  int min_flow_count { 1 };
  std::optional<int> max_flow_count { std::nullopt };
  std::optional<Cost> max_path_cost { std::nullopt };
  std::optional<double> max_path_cost_factor { std::nullopt };
  bool shortest_path { false };
  bool reoptimize_flows_on_each_placement { false };
  int max_no_progress_iterations { 100 };
  int max_total_iterations { 10000 };
  bool diminishing_returns_enabled { true };
  int diminishing_returns_window { 8 };
  double diminishing_returns_epsilon_frac { 1e-3 };
};

// FlowPolicy manages flow creation, placement, reoptimization for a single demand
class FlowPolicy {
public:
  // New config-based constructor
  FlowPolicy(const ExecutionContext& ctx, const FlowPolicyConfig& cfg)
    : ctx_(ctx),
      path_alg_(cfg.path_alg), flow_placement_(cfg.flow_placement), selection_(cfg.selection),
      shortest_path_(cfg.shortest_path),
      min_flow_count_(cfg.min_flow_count), max_flow_count_(cfg.max_flow_count), max_path_cost_(cfg.max_path_cost),
      max_path_cost_factor_(cfg.max_path_cost_factor), reoptimize_flows_on_each_placement_(cfg.reoptimize_flows_on_each_placement),
      max_no_progress_iterations_(cfg.max_no_progress_iterations), max_total_iterations_(cfg.max_total_iterations),
      diminishing_returns_enabled_(cfg.diminishing_returns_enabled), diminishing_returns_window_(cfg.diminishing_returns_window),
      diminishing_returns_epsilon_frac_(cfg.diminishing_returns_epsilon_frac) {
    if (flow_placement_ == FlowPlacement::EqualBalanced && !max_flow_count_.has_value()) {
      throw std::invalid_argument("max_flow_count must be set for EQUAL_BALANCED placement.");
    }
  }

  FlowPolicy(const ExecutionContext& ctx,
             PathAlg path_alg,
             FlowPlacement flow_placement,
             EdgeSelection selection,
             int min_flow_count = 1,
             std::optional<int> max_flow_count = std::nullopt,
             std::optional<Cost> max_path_cost = std::nullopt,
             std::optional<double> max_path_cost_factor = std::nullopt,
             bool shortest_path = false,
             bool reoptimize_flows_on_each_placement = false,
             int max_no_progress_iterations = 100,
             int max_total_iterations = 10000,
             bool diminishing_returns_enabled = true,
             int diminishing_returns_window = 8,
             double diminishing_returns_epsilon_frac = 1e-3)
    : ctx_(ctx),
      path_alg_(path_alg), flow_placement_(flow_placement), selection_(selection),
      shortest_path_(shortest_path),
      min_flow_count_(min_flow_count), max_flow_count_(max_flow_count), max_path_cost_(max_path_cost),
      max_path_cost_factor_(max_path_cost_factor), reoptimize_flows_on_each_placement_(reoptimize_flows_on_each_placement),
      max_no_progress_iterations_(max_no_progress_iterations), max_total_iterations_(max_total_iterations),
      diminishing_returns_enabled_(diminishing_returns_enabled), diminishing_returns_window_(diminishing_returns_window),
      diminishing_returns_epsilon_frac_(diminishing_returns_epsilon_frac) {
    // Constructor validations enforcing policy invariants
    if (flow_placement_ == FlowPlacement::EqualBalanced && !max_flow_count_.has_value()) {
      throw std::invalid_argument("max_flow_count must be set for EQUAL_BALANCED placement.");
    }
  }
  ~FlowPolicy() noexcept = default;

  [[nodiscard]] int flow_count() const noexcept { return static_cast<int>(flows_.size()); }
  [[nodiscard]] double placed_demand() const noexcept;

  // Core operations
  [[nodiscard]] std::pair<double,double> place_demand(FlowGraph& fg,
                                        NodeId src, NodeId dst,
                                        FlowClass flowClass,
                                        double volume,
                                        std::optional<double> target_per_flow = std::nullopt,
                                        std::optional<double> min_flow = std::nullopt);

  [[nodiscard]] std::pair<double,double> rebalance_demand(FlowGraph& fg,
                                            NodeId src, NodeId dst,
                                            FlowClass flowClass,
                                            double target_per_flow);

  void remove_demand(FlowGraph& fg);

  [[nodiscard]] const std::unordered_map<FlowIndex, FlowRecord, FlowIndexHash>& flows() const noexcept { return flows_; }

// Configure static paths for flow creation. Each entry is (src, dst, dag, cost).
// max_flow_count must equal the number of static paths if set.
  void set_static_paths(std::vector<std::tuple<NodeId, NodeId, PredDAG, Cost>> paths);

private:
  // Helpers
  [[nodiscard]] std::optional<std::pair<PredDAG, Cost>> get_path_bundle(const FlowGraph& fg,
                                                          NodeId src, NodeId dst,
                                                          std::optional<double> min_flow);
  [[nodiscard]] FlowRecord* create_flow(FlowGraph& fg, NodeId src, NodeId dst, FlowClass flowClass,
                    std::optional<double> min_flow);
  [[nodiscard]] FlowRecord* reoptimize_flow(FlowGraph& fg, const FlowIndex& idx, double headroom);

  // Config
  ExecutionContext ctx_;
  PathAlg path_alg_ { PathAlg::SPF };
  FlowPlacement flow_placement_ { FlowPlacement::Proportional };
  EdgeSelection selection_ { EdgeSelection{} };
  bool shortest_path_ { false };
  int min_flow_count_ { 1 };
  std::optional<int> max_flow_count_ {};
  std::optional<Cost> max_path_cost_ {};
  std::optional<double> max_path_cost_factor_ {};
  bool reoptimize_flows_on_each_placement_ { false };
  int max_no_progress_iterations_ { 100 };
  int max_total_iterations_ { 10000 };
  bool diminishing_returns_enabled_ { true };
  int diminishing_returns_window_ { 8 };
  double diminishing_returns_epsilon_frac_ { 1e-3 };

  // State
  std::unordered_map<FlowIndex, FlowRecord, FlowIndexHash> flows_;
  Cost best_path_cost_ { std::numeric_limits<Cost>::max() };
  FlowId next_flow_id_ { 0 };

  // Static paths (optional)
  std::vector<std::tuple<NodeId, NodeId, PredDAG, Cost>> static_paths_;
};

} // namespace netgraph::core
