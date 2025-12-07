/*
  FlowPolicy — manages path selection and flow placement for a single demand.

  Responsibilities:
    - Select shortest-path predecessors (SPF) per policy via `shortest_paths`.
    - Create and track flows with their current DAG and path cost.
    - Place volumes using `FlowPlacement` (Proportional or EqualBalanced).
    - Optionally re-optimize flows and enforce path-cost constraints.

  Notes:
    - "Residual-aware" SPF considers current remaining capacity (or a per-flow
      threshold) when exploring edges.
    - Edge selection (multipath vs single-path) is handled natively by the SPF algorithm.
*/
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/constants.hpp"
#include "netgraph/core/algorithms.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/profiling.hpp"

#include <algorithm>
#include <deque>
#include <limits>
#include <optional>
#include <unordered_set>

namespace netgraph::core {

double FlowPolicy::placed_demand() const noexcept {
  double s = 0.0;
  for (auto const& kv : flows_) s += kv.second.placed_flow;
  return s;
}

/* Compute an SPF predecessor DAG and its destination cost under the current
   policy configuration. Honors static paths, residual capacity, and edge/node masks. */
std::optional<std::pair<PredDAG, Cost>> FlowPolicy::get_path_bundle(const FlowGraph& fg,
                                                                    NodeId src, NodeId dst,
                                                                    std::optional<double> min_flow) {
  // Static path handling: if static paths are configured, use them exclusively.
  if (!static_paths_.empty()) {
    // Search for a static path matching this (src, dst) pair.
    for (auto const& t : static_paths_) {
      if (std::get<0>(t) == src && std::get<1>(t) == dst) {
        return std::make_optional(std::make_pair(std::get<2>(t), std::get<3>(t)));
      }
    }
    return std::nullopt;  // No static path for this pair
  }
  if (path_alg_ != PathAlg::SPF) return std::nullopt;

  // Use configured selection for per-adjacency edge behavior (multi-edge, tie-breaking).
  EdgeSelection sel = selection_;

  // Enforce semantic consistency between multipath and multi_edge:
  // - Tunnel mode (multipath=false): force single edge per hop for true single-path semantics
  // - Hash-ECMP with EqualBalanced: use all equal-cost edges to maximize fanout
  if (!multipath_) {
    sel.multi_edge = false;
  } else if (flow_placement_ == FlowPlacement::EqualBalanced) {
    sel.multi_edge = true;
  }

  // Respect capacity requirements from both config sources
  sel.require_capacity = (selection_.require_capacity || require_capacity_);
  // Decide whether we need residual-aware SPF.
  // Residual awareness is controlled by require_capacity_:
  //   - require_capacity=true: Require edges to have capacity, routes adapt to residuals (SDN/TE behavior)
  //   - require_capacity=false: Routes based on costs only (IP/IGP behavior)
  // Additionally, for EqualBalanced mode with minimum flow threshold, we use residuals.
  const bool require_residual = (require_capacity_ || (flow_placement_ == FlowPlacement::EqualBalanced && min_flow.has_value()));
  const auto residual = fg.residual_view();

  // Edge mask: combine user-provided mask with minimum residual capacity threshold.
  // For Proportional mode with min_flow threshold and max_flow_count limit, we filter edges
  // by minimum capacity to ensure paths can accommodate the required flow. This prevents
  // selecting low-capacity paths when limited to a single flow (max_flow_count=1).
  // For EqualBalanced mode, we skip per-edge thresholds since group-based semantics mean
  // per-edge thresholds can over-prune; capacity gating via residual + kMinCap is sufficient.
  std::unique_ptr<bool[]> combined_edge_mask;
  std::span<const bool> final_edge_mask;

  if (require_residual && min_flow.has_value() && flow_placement_ != FlowPlacement::EqualBalanced) {
    // Need to filter by min_flow threshold for Proportional mode
    combined_edge_mask.reset(new bool[residual.size()]);
    double thr = *min_flow;

    if (!edge_mask_.empty()) {
      // Combine user mask with min_flow mask: both must be true
      for (std::size_t i=0; i<residual.size(); ++i) {
        combined_edge_mask[i] = edge_mask_[i] && (static_cast<double>(residual[i]) >= thr);
      }
    } else {
      // Only min_flow mask
      for (std::size_t i=0; i<residual.size(); ++i) {
        combined_edge_mask[i] = static_cast<double>(residual[i]) >= thr;
      }
    }
    final_edge_mask = std::span<const bool>(combined_edge_mask.get(), residual.size());
  } else if (!edge_mask_.empty()) {
    // Only user-provided mask
    final_edge_mask = edge_mask_;
  }

  SpfOptions opts;
  opts.multipath = multipath_;  // Use configured multipath value (enables/disables flow splitting across equal-cost paths)
  opts.selection = sel;
  opts.dst = dst;
  opts.residual = require_residual ? residual : std::span<const Cap>();
  opts.node_mask = node_mask_;  // Use user-provided node mask
  opts.edge_mask = final_edge_mask;
  auto res = ctx_.algorithms->spf(ctx_.graph, src, opts);
  const auto& dist = res.first;
  PredDAG dag = std::move(res.second);
  if (dst < 0 || static_cast<std::size_t>(dst) >= dist.size()) return std::nullopt;
  Cost dst_cost = dist[static_cast<std::size_t>(dst)];
  if (dst_cost < best_path_cost_) best_path_cost_ = dst_cost;

  // Enforce path cost constraints:
  // 1. In shortest_path mode, only allow paths with cost equal to best discovered cost.
  //    This prevents the policy from stepping up to higher-cost tiers incrementally.
  if (shortest_path_ && dst_cost > best_path_cost_) {
    return std::nullopt;
  }

  // 2. Check absolute and relative path cost limits.
  //    max_path_cost: absolute upper bound on path cost.
  //    max_path_cost_factor: relative multiplier on best path cost (e.g. 1.5 = allow 50% longer).
  if (max_path_cost_.has_value() || max_path_cost_factor_.has_value()) {
    double maxf = max_path_cost_factor_.value_or(1.0);
    Cost absmax = max_path_cost_.value_or(std::numeric_limits<Cost>::max());
    if (dst_cost > std::min<Cost>(absmax, static_cast<Cost>(static_cast<double>(best_path_cost_) * maxf))) return std::nullopt;
  }
  // Ensure there is at least one predecessor for dst
  if (static_cast<std::size_t>(dst) >= dag.parent_offsets.size()-1) return std::nullopt;
  if (dag.parent_offsets[static_cast<std::size_t>(dst)] == dag.parent_offsets[static_cast<std::size_t>(dst)+1]) return std::nullopt;
  // Return DAG and cost as-is; placement logic decides proportional vs equal-balanced behavior.
  return std::make_optional(std::make_pair(dag, dst_cost));
}

/* Create a new flow using the current path bundle. Returns nullptr if no
   admissible path is available given constraints. */
FlowRecord* FlowPolicy::create_flow(FlowGraph& fg, NodeId src, NodeId dst, FlowClass flowClass,
                              std::optional<double> min_flow) {
  // Generate a unique flow index.
  FlowIndex idx{src, dst, flowClass, next_flow_id_++};

  // Request a path DAG from shortest paths algorithm.
  auto pb = get_path_bundle(fg, src, dst, min_flow);
  if (!pb.has_value()) return nullptr;  // No admissible path found

  // Destructure the returned pair (DAG, cost).
  auto [dag, cost] = std::move(pb.value());

  // Create flow record and insert into flows_ map.
  FlowRecord f(idx, src, dst, std::move(dag), cost);
  auto [it, ok] = flows_.emplace(idx, std::move(f));  // emplace returns (iterator, success_flag)
  (void)ok;  // suppress unused variable warning
  return &it->second;
}

/* Re-select a path for an existing flow, requesting at least (current+headroom)
   residual. On failure, restores the flow on its original DAG.

   Reoptimization is useful when a flow's current path becomes suboptimal due to
   network changes or when seeking additional capacity. */
FlowRecord* FlowPolicy::reoptimize_flow(FlowGraph& fg, const FlowIndex& idx, double headroom) {
  auto it = flows_.find(idx);
  if (it == flows_.end()) return nullptr;
  FlowRecord& cur = it->second;
  const double current = cur.placed_flow;
  const double new_min = current + headroom;

  // Temporarily remove this flow's deltas from the graph to compute a path
  // based on available capacity (excluding this flow's own usage).
  fg.remove(idx);
  auto pb = get_path_bundle(fg, cur.src, cur.dst, new_min);
  if (!pb.has_value()) {
    // Reoptimization failed: restore flow on original DAG.
    Flow placed = fg.place(idx, cur.src, cur.dst, cur.dag, current, flow_placement_);
    cur.placed_flow = placed; // may be slightly less if capacity changed; acceptable
    return nullptr;
  }

  // Reoptimization succeeded: update flow to use new DAG.
  auto [dag, cost] = std::move(pb.value());
  cur.dag = std::move(dag);
  cur.cost = cost;
  Flow placed = fg.place(idx, cur.src, cur.dst, cur.dag, current, flow_placement_);
  cur.placed_flow = placed;
  return &cur;
}

/* Place `volume` of demand according to the policy. When `target_per_flow`
   is provided (e.g., during rebalancing), each flow aims for that target.
   Returns (total_placed, leftover). */
std::pair<double,double> FlowPolicy::place_demand(FlowGraph& fg,
                                                  NodeId src, NodeId dst,
                                                  FlowClass flowClass,
                                                  double volume,
                                                  std::optional<double> target_per_flow,
                                                  std::optional<double> min_flow) {
  NGRAPH_PROFILE_SCOPE("place_demand");

  // Compute target flow per flow-record.
  // target: the volume to place per flow (or globally if target_per_flow is unset).
  // per_target: refined target for EqualBalanced mode (considers source capacity).
  double target = target_per_flow.value_or(volume);
  double per_target = target;

  // For EqualBalanced mode with max_flow_count, compute a per-flow target based on
  // available source capacity and the requested volume, divided by the number of flows.
  if (flow_placement_ == FlowPlacement::EqualBalanced && max_flow_count_.has_value()) {
    const auto& g = fg.graph();
    auto row = g.row_offsets_view();
    auto aei = g.adj_edge_index_view();
    auto residual = fg.residual_view();

    // Compute total residual capacity on edges leaving src.
    double src_cap = 0.0;
    if (src >= 0 && src < g.num_nodes()) {
      auto s = static_cast<std::size_t>(row[static_cast<std::size_t>(src)]);
      auto e = static_cast<std::size_t>(row[static_cast<std::size_t>(src)+1]);
      for (std::size_t j = s; j < e; ++j) {
        auto eid = static_cast<std::size_t>(aei[j]);
        src_cap += static_cast<double>(residual[eid]);
      }
    }
    // Compute per-flow target as the minimum of:
    // - requested volume / max_flow_count
    // - source capacity / max_flow_count
    double per_req = target / static_cast<double>(*max_flow_count_);
    double per_src = src_cap / static_cast<double>(*max_flow_count_);
    per_target = std::max(kMinFlow, std::min(per_req, per_src));
  }

  // Initialize flows if none exist yet.
  if (flows_.empty()) {
    if (!static_paths_.empty()) {
      // Static paths: create one flow per static path.
      if (max_flow_count_.has_value() && static_cast<int>(static_paths_.size()) != *max_flow_count_) {
        throw std::invalid_argument("If set, max_flow_count must be equal to the number of static paths.");
      }
      for (auto const& t : static_paths_) {
        if (std::get<0>(t) == src && std::get<1>(t) == dst) {
          [[maybe_unused]] auto* created = create_flow(fg, src, dst, flowClass, std::nullopt);
        } else {
          throw std::invalid_argument("Source and destination nodes of static paths do not match demand.");
        }
      }
    } else {
      // Dynamic paths: seed initial flows.
      int initial = min_flow_count_;
      if (max_flow_count_.has_value()) {
        initial = std::min(initial, *max_flow_count_);
      }
      for (int i=0;i<initial;++i) {
        auto min_req = (flow_placement_ == FlowPlacement::EqualBalanced && max_flow_count_.has_value())
                         ? std::optional<double>(per_target)
                         : min_flow;
        [[maybe_unused]] auto* created = create_flow(fg, src, dst, flowClass, min_req);
      }
    }
  }
  // Round-robin placement: iterate over flows and try to place volume on each.
  std::deque<FlowIndex> q;
  for (auto const& kv : flows_) q.push_back(kv.first);
  double total_placed = 0.0;
  int no_progress = 0;  // counter for consecutive iterations with no progress
  int iters = 0;

  // Diminishing-returns tracking for early exit.
  std::deque<double> recent;
  const double initial_request = volume;

  while (volume >= kMinFlow && !q.empty()) {
    FlowIndex cur_idx = q.front(); q.pop_front();
    auto it_cur = flows_.find(cur_idx);
    if (it_cur == flows_.end()) continue;  // flow removed during iteration
    FlowRecord* f = &it_cur->second;

    // Must have a DAG to place; skip otherwise.
    if (f->dag.parent_offsets.empty()) {
      ++no_progress;
      if (no_progress>=max_no_progress_iterations_) break;
      continue;
    }
    // Refresh DAG based on current residuals for dynamic path selection.
    // This prunes saturated next-hops and updates path selection.
    // For multipath flows, this tracks saturated edges within the DAG.
    // For tunnel flows, this allows different tunnels to discover different paths
    // as residuals change, enabling natural fan-out across equal-cost paths.
    if (flow_placement_ == FlowPlacement::EqualBalanced && static_paths_.empty()) {
      if (auto pb = get_path_bundle(fg, f->src, f->dst, std::optional<double>(per_target))) {
        f->dag = std::move(pb->first);
        f->cost = pb->second;
      }
    }
    double need;
    if (target_per_flow.has_value()) {
      // When a per-flow target is specified (e.g., during rebalancing), cap by remaining per-flow target.
      need = std::max(0.0, target - f->placed_flow);
    } else if (flow_placement_ == FlowPlacement::EqualBalanced && max_flow_count_.has_value()) {
      // For EqualBalanced, request only the remaining deficit toward per-target for this flow.
      need = std::max(0.0, per_target - f->placed_flow);
    } else {
      // Default behavior uses the global target amount.
      need = target;
    }
    const double request = std::min(need, volume);
    Flow placed = fg.place(f->index, f->src, f->dst, f->dag, request, flow_placement_);
    f->placed_flow += placed;
    volume -= placed;
    total_placed += placed;
    ++iters;
    // IP-like mode: perform a single augmentation over the current SPF DAG
    if (shortest_path_) {
      break;
    }
    // track recent placements
    if (diminishing_returns_enabled_) {
      recent.push_back(placed);
      if (static_cast<int>(recent.size()) > diminishing_returns_window_) recent.pop_front();
      if (static_cast<int>(recent.size()) == diminishing_returns_window_) {
        double sum_recent = 0.0; for (double x : recent) sum_recent += x;
        const double threshold = std::max(kMinFlow, diminishing_returns_epsilon_frac_ * initial_request);
        if (sum_recent < threshold) {
          break; // graceful cutoff
        }
      }
    }
    if (placed < kMinFlow) {
      ++no_progress; if (no_progress>=max_no_progress_iterations_) break;
    } else {
      no_progress = 0;
    }
    if (flow_placement_ == FlowPlacement::EqualBalanced) {
      if (max_flow_count_.has_value()) {
        // Bounded EB: add flows up to configured maximum.
        if (static_cast<int>(flows_.size()) < *max_flow_count_) {
          if (auto* nf = create_flow(fg, src, dst, flowClass, std::optional<double>(per_target))) q.push_back(nf->index);
        }
      } else {
        // Unbounded EB: rely on a single flow to equalize over the DAG.
        // Do not create additional flows implicitly.
      }
    } else {
      if (target - f->placed_flow >= kMinFlow) {
        if (!max_flow_count_ || static_cast<int>(flows_.size()) < *max_flow_count_) {
          if (auto* nf = create_flow(fg, src, dst, flowClass, std::nullopt)) q.push_back(nf->index);
        } else {
          if (auto* rf = reoptimize_flow(fg, f->index, kMinFlow)) q.push_back(rf->index);
        }
      }
    }
    if (iters >= max_total_iterations_) break;
  }

  // Reoptimize all flows after placement if enabled
  if (reoptimize_flows_on_each_placement_) {
    for (auto& kv : flows_) {
      (void)reoptimize_flow(fg, kv.first, kMinFlow);
    }
  }

  // For EQUAL_BALANCED placement, rebalance flows to maintain equal volumes.
  if (flow_placement_ == FlowPlacement::EqualBalanced && !flows_.empty()) {
    double target_eq = placed_demand() / static_cast<double>(flows_.size());
    bool unbalanced = false;
    for (auto const& kv : flows_) {
      if (std::abs(target_eq - kv.second.placed_flow) >= kMinFlow) { unbalanced = true; break; }
    }
    if (unbalanced) {
      bool prev_reopt = reoptimize_flows_on_each_placement_;
      reoptimize_flows_on_each_placement_ = false;
      auto pr = rebalance_demand(fg, src, dst, flowClass, target_eq);
      // pr.first = placed in rebalanced pass, pr.second = excess
      volume += pr.second; // leave remaining volume
      reoptimize_flows_on_each_placement_ = prev_reopt;
      total_placed = placed_demand();
    }
  }
  return { total_placed, volume };
}

/* Rebalance existing placed demand such that each flow carries approximately
   `target_per_flow`. Internally removes and re-places the same total volume. */
std::pair<double,double> FlowPolicy::rebalance_demand(FlowGraph& fg,
                                                      NodeId src, NodeId dst,
                                                      FlowClass flowClass,
                                                      double target_per_flow) {
  double vol = placed_demand();
  remove_demand(fg);
  return place_demand(fg, src, dst, flowClass, vol, target_per_flow, std::nullopt);
}

/* Remove all placed flows for this policy from the FlowGraph and reset
   per-flow placed volumes. */
void FlowPolicy::remove_demand(FlowGraph& fg) {
  for (auto const& kv : flows_) {
    fg.remove(kv.first);
  }
  flows_.clear();
  best_path_cost_ = std::numeric_limits<Cost>::max();
}

/* Configure static paths to be used instead of dynamic SPF selection. If
   `max_flow_count` is not set, it is set to the number of provided paths. */
void FlowPolicy::set_static_paths(std::vector<std::tuple<NodeId, NodeId, PredDAG, Cost>> paths) {
  static_paths_ = std::move(paths);
  if (max_flow_count_.has_value() && static_cast<int>(static_paths_.size()) != *max_flow_count_) {
    throw std::invalid_argument("If set, max_flow_count must be equal to the number of static paths.");
  }
  if (!max_flow_count_.has_value()) {
    max_flow_count_ = static_cast<int>(static_paths_.size());
  }
}

} // namespace netgraph::core
