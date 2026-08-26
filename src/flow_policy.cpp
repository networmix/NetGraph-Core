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
#include <cstring>
#include <deque>
#include <limits>
#include <optional>

namespace netgraph::core {

/* Reject uses that would silently produce a wrong answer:
   - a FlowGraph wrapping a different graph than the policy routes on (SPF would run
     on one topology while flow is placed on another);
   - a (src, dst) pair different from the demand this policy already manages (the
     round-robin loop routes using src/dst from the existing flow records). */
void FlowPolicy::check_demand_target(const FlowGraph& fg, NodeId src, NodeId dst,
                                     const char* what) const {
  const auto* policy_graph = ctx_.graph.graph.get();
  if (policy_graph != nullptr && policy_graph != &fg.graph()) {
    throw std::invalid_argument(
        std::string("FlowPolicy::") + what +
        ": the FlowGraph wraps a different StrictMultiDiGraph than this policy; "
        "paths would be selected on one topology and placed on another");
  }
  if (!flows_.empty()) {
    const auto& existing = flows_.begin()->second;
    if (existing.src != src || existing.dst != dst) {
      throw std::invalid_argument(
          std::string("FlowPolicy::") + what +
          ": this policy already manages a demand for a different (src, dst) pair; "
          "use a separate FlowPolicy per demand or call remove_demand() first");
    }
  }
}

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
  // With static paths, flows are created directly from the pinned bundles in
  // place_demand_body and never reoptimized, so no dynamic bundle exists to hand
  // out. Callers reaching here with static paths configured get nothing.
  if (has_static_paths()) return std::nullopt;
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
  // dist.size() == num_nodes() for every SPF result, so the range check does
  // not depend on the SPF output; hoisting it lets the memo skip the call.
  if (dst < 0 || dst >= ctx_.graph.graph->num_nodes()) return std::nullopt;

  PredDAG dag;
  Cost dst_cost;
  const bool use_memo = (flow_placement_ == FlowPlacement::EqualBalanced);
  const auto stamp = fg.state_stamp();
  const bool with_residual = !opts.residual.empty();
  std::size_t hit_idx = spf_memo_.size();
  if (use_memo) {
    for (std::size_t i = 0; i < spf_memo_.size(); ++i) {
      const auto& e = spf_memo_[i];
      // min_flow's VALUE never reaches SPF in EB mode (the value-derived edge
      // mask is Proportional-only); only has_value() matters, via
      // require_residual. Keying on the value caused misses whenever rebalance
      // rounds adjusted the per-flow target against unchanged residuals.
      if (e.src != src || e.dst != dst || e.with_residual != with_residual ||
          e.has_min_flow != min_flow.has_value()) {
        continue;
      }
      // Fast path: same FlowGraph, no mutation since the entry was stored.
      // Content path: byte-identical residuals (rebalance remove+place
      // round-trips restore content while the version keeps advancing).
      const bool same = !with_residual || e.stamp == stamp ||
          (e.residual.size() == opts.residual.size() &&
           std::memcmp(e.residual.data(), opts.residual.data(),
                       opts.residual.size() * sizeof(Cap)) == 0);
      if (same) { hit_idx = i; break; }
    }
  }
  if (hit_idx < spf_memo_.size()) {
    dag = spf_memo_[hit_idx].dag;
    dst_cost = spf_memo_[hit_idx].dst_cost;
    // MRU: move the hit to the front so hot entries stay cheap to find.
    if (hit_idx != 0) {
      std::rotate(spf_memo_.begin(), spf_memo_.begin() + hit_idx,
                  spf_memo_.begin() + hit_idx + 1);
    }
  } else {
    auto res = ctx_.algorithms->spf(ctx_.graph, src, opts);
    dag = std::move(res.second);
    dst_cost = res.first[static_cast<std::size_t>(dst)];
    if (use_memo) {
      SpfMemoEntry entry;
      entry.src = src;
      entry.dst = dst;
      entry.with_residual = with_residual;
      entry.stamp = stamp;
      entry.residual.assign(opts.residual.begin(), opts.residual.end());
      entry.has_min_flow = min_flow.has_value();
      entry.dag = dag;
      entry.dst_cost = dst_cost;
      const std::size_t entry_bytes =
          entry.residual.size() * sizeof(Cap) +
          entry.dag.parent_offsets.size() * sizeof(std::int32_t) +
          entry.dag.parents.size() * (sizeof(NodeId) + sizeof(EdgeId)) +
          sizeof(SpfMemoEntry);
      const std::size_t cap = std::clamp<std::size_t>(
          kSpfMemoMaxBytes / std::max<std::size_t>(entry_bytes, 1), 1,
          kSpfMemoMaxEntries);
      spf_memo_.insert(spf_memo_.begin(), std::move(entry));
      if (spf_memo_.size() > cap) spf_memo_.resize(cap);
    }
  }
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
    const Cost absmax = max_path_cost_.value_or(std::numeric_limits<Cost>::max());
    // Apply the relative bound only when a best cost exists: multiplying the
    // INT64_MAX "no best yet" sentinel by the factor and casting back is
    // undefined behavior (the product exceeds the int64 range).
    Cost factor_bound = std::numeric_limits<Cost>::max();
    if (max_path_cost_factor_.has_value() &&
        best_path_cost_ != std::numeric_limits<Cost>::max()) {
      const double prod = static_cast<double>(best_path_cost_) * *max_path_cost_factor_;
      if (prod < static_cast<double>(std::numeric_limits<Cost>::max())) {
        factor_bound = static_cast<Cost>(prod);
      }
    }
    if (dst_cost > std::min<Cost>(absmax, factor_bound)) return std::nullopt;
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
  // Pinned means pinned: a static-path flow is never rerouted. (Deliberate
  // divergence from the original Python port, where reoptimization could
  // silently move a pinned flow onto an SPF path.) Must precede any
  // remove/re-place churn below.
  if (has_static_paths()) return nullptr;
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

  check_demand_target(fg, src, dst, "place_demand");

  auto [total_placed, remaining] = place_demand_body(fg, src, dst, flowClass, volume,
                                                     target_per_flow, min_flow);

  // For EQUAL_BALANCED placement, rebalance flows to maintain equal volumes.
  //
  // Iterative on purpose: the previous implementation recursed
  // place_demand -> rebalance_demand -> place_demand until balanced, and with many
  // pinned bundles of heterogeneous capacity the depth grows like
  // U * ln(imbalance / kMinFlow) -- a stack-overflow risk on worker threads. Each
  // round below performs exactly what one recursion level performed (re-place the
  // currently placed volume at the equal-share target), and the recursion's return
  // value telescoped to (placed_demand(), pre-rebalance leftover + volume lost in
  // rebalancing), which is reproduced after the loop.
  if (flow_placement_ == FlowPlacement::EqualBalanced && !flows_.empty()) {
    // Restore the reoptimize flag even if a round throws (bad_alloc is the only
    // realistic thrower here); otherwise the policy would stay permanently
    // non-reoptimizing.
    struct ReoptRestore {
      bool& flag;
      bool prev;
      ~ReoptRestore() { flag = prev; }
    } reopt_restore{reoptimize_flows_on_each_placement_, reoptimize_flows_on_each_placement_};
    reoptimize_flows_on_each_placement_ = false;
    const double pre_rounds_placed = placed_demand();
    bool rebalanced = false;
    // Backstop only; each round strictly reduces the imbalance in practice.
    constexpr int kMaxRebalanceRounds = 65536;
    for (int round = 0; round < kMaxRebalanceRounds && !flows_.empty(); ++round) {
      const double target_eq = placed_demand() / static_cast<double>(flows_.size());
      bool unbalanced = false;
      for (auto const& kv : flows_) {
        if (std::abs(target_eq - kv.second.placed_flow) >= kMinFlow) { unbalanced = true; break; }
      }
      if (!unbalanced) break;
      rebalanced = true;
      const double vol = placed_demand();
      remove_demand(fg);
      (void)place_demand_body(fg, src, dst, flowClass, vol, target_eq, std::nullopt);
    }
    if (rebalanced) {
      const double final_placed = placed_demand();
      remaining += pre_rounds_placed - final_placed;  // volume shed while rebalancing
      total_placed = final_placed;
    }
  }
  return { total_placed, remaining };
}

/* Placement core; see place_demand for the public contract. */
std::pair<double,double> FlowPolicy::place_demand_body(FlowGraph& fg,
                                                       NodeId src, NodeId dst,
                                                       FlowClass flowClass,
                                                       double volume,
                                                       std::optional<double> target_per_flow,
                                                       std::optional<double> min_flow) {
  const bool is_static = has_static_paths();

  // Compute target flow per flow-record.
  // target: the volume to place per flow (or globally if target_per_flow is unset).
  // per_target: refined target for EqualBalanced mode (considers source capacity).
  double target = target_per_flow.value_or(volume);
  double per_target = target;

  // EqualBalanced divisor: the number of flows the volume is split across. For a
  // pinned policy that is the number of USABLE bundles (a head-end hashes over up
  // LSPs only); dynamically it is the configured max_flow_count.
  int eb_divisor = 0;
  if (flow_placement_ == FlowPlacement::EqualBalanced) {
    if (is_static) eb_divisor = static_cast<int>(static_bundles_.size());
    else if (max_flow_count_.has_value()) eb_divisor = *max_flow_count_;
  }

  // For EqualBalanced, compute a per-flow target based on available source
  // capacity and the requested volume, divided by the number of flows.
  if (eb_divisor > 0) {
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
    // Per-flow target: min of requested volume and source capacity, per flow.
    double per_req = target / static_cast<double>(eb_divisor);
    double per_src = src_cap / static_cast<double>(eb_divisor);
    per_target = std::max(kMinFlow, std::min(per_req, per_src));
  }

  // Initialize flows if none exist yet.
  if (flows_.empty()) {
    if (is_static) {
      // Pinned bundles: one flow per usable bundle, in supply order, each bound to
      // its own (pruned) DAG and post-prune cost.
      if (src != static_src_ || dst != static_dst_) {
        throw std::invalid_argument(
            "Source and destination nodes of static paths do not match demand.");
      }
      for (auto const& b : static_bundles_) {
        FlowIndex idx{src, dst, flowClass, next_flow_id_++};
        flows_.emplace(idx, FlowRecord(idx, src, dst, b.dag, b.cost));
      }
    } else {
      // Dynamic paths: seed initial flows.
      int initial = min_flow_count_;
      if (max_flow_count_.has_value()) {
        initial = std::min(initial, *max_flow_count_);
      }
      auto min_req = (flow_placement_ == FlowPlacement::EqualBalanced && max_flow_count_.has_value())
                       ? std::optional<double>(per_target)
                       : min_flow;
      // Seeding places no flow, so residuals do not change between iterations and
      // every create_flow() here would recompute the identical SPF. Compute the
      // bundle once and copy it into each seeded flow.
      if (initial > 0) {
        if (auto pb = get_path_bundle(fg, src, dst, min_req)) {
          for (int i = 0; i < initial; ++i) {
            FlowIndex idx{src, dst, flowClass, next_flow_id_++};
            FlowRecord f(idx, src, dst, pb->first, pb->second);
            flows_.emplace(idx, std::move(f));
          }
        }
      }
    }
  }
  // Round-robin placement: iterate over flows and try to place volume on each.
  std::deque<FlowIndex> q;
  for (auto const& kv : flows_) q.push_back(kv.first);
  if (is_static) {
    // Pinned bundles have a documented precedence: bundle-supply order. Flow ids
    // are assigned in that order, so sorting makes the visit order deterministic
    // across platforms (unordered_map iteration order is not).
    std::sort(q.begin(), q.end(),
              [](const FlowIndex& a, const FlowIndex& b) { return a.flowId < b.flowId; });
  }
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
    if (flow_placement_ == FlowPlacement::EqualBalanced && !is_static) {
      if (auto pb = get_path_bundle(fg, f->src, f->dst, std::optional<double>(per_target))) {
        f->dag = std::move(pb->first);
        f->cost = pb->second;
      }
    }
    double need;
    if (target_per_flow.has_value()) {
      // When a per-flow target is specified (e.g., during rebalancing), cap by remaining per-flow target.
      need = std::max(0.0, target - f->placed_flow);
    } else if (eb_divisor > 0) {
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
    // Track recent placements. For a pinned policy, arm the window only once every
    // flow has been visited: with more bundles than the window, the early rounds of
    // small per-flow placements must not starve the bundles not yet visited.
    if (diminishing_returns_enabled_ &&
        (!is_static || iters >= static_cast<int>(flows_.size()))) {
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
    // A pinned policy neither grows its flow set nor reoptimizes: the pinned-ness
    // guard is explicit, never inferred from flow-count arithmetic.
    if (!is_static) {
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
    }
    if (iters >= max_total_iterations_) break;
  }

  // Reoptimize all flows after placement if enabled (no-op for pinned policies:
  // reoptimize_flow returns immediately when static paths are configured).
  if (reoptimize_flows_on_each_placement_) {
    for (auto& kv : flows_) {
      (void)reoptimize_flow(fg, kv.first, kMinFlow);
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
  // Must run before remove_demand() empties flows_, or the check inside
  // place_demand() would have nothing left to compare against and would
  // silently retarget this policy's volume onto a different node pair.
  check_demand_target(fg, src, dst, "rebalance_demand");
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

namespace {

/* Validate one pinned bundle against the graph, prune it against the policy's
   masks, and compute its post-prune min cost. Throws std::invalid_argument on a
   malformed bundle; returns nullopt for a structurally valid bundle that has no
   surviving src->dst walk under the masks (a DOWN LSP). */
std::optional<std::pair<PredDAG, Cost>>
prepare_static_bundle(const StrictMultiDiGraph& g, const PredDAG& dag,
                      NodeId src, NodeId dst,
                      std::span<const bool> node_mask,
                      std::span<const bool> edge_mask) {
  const auto N = static_cast<std::size_t>(g.num_nodes());
  const auto E = static_cast<std::size_t>(g.num_edges());
  const auto& off = dag.parent_offsets;
  const auto& par = dag.parents;
  const auto& via = dag.via_edges;

  // Shape: CSR offsets over N nodes, entry arrays sized by the final offset.
  if (off.size() != N + 1 || off.front() != 0) {
    throw std::invalid_argument("static path bundle: parent_offsets must have length num_nodes + 1 and start at 0");
  }
  for (std::size_t v = 0; v + 1 < off.size(); ++v) {
    if (off[v] > off[v + 1]) {
      throw std::invalid_argument("static path bundle: parent_offsets must be non-decreasing");
    }
  }
  const auto entries = static_cast<std::size_t>(off.back());
  if (par.size() != entries || via.size() != entries) {
    throw std::invalid_argument("static path bundle: parents/via_edges size must equal parent_offsets.back()");
  }

  // Entries: ids in range, and each via edge must actually connect parent -> child
  // in THIS graph. Without the endpoint check, a well-shaped DAG built for a
  // different graph (or with shuffled edge ids) would silently place flow on
  // arbitrary edges.
  const auto esrc = g.edge_src_view();
  const auto edst = g.edge_dst_view();
  for (std::size_t v = 0; v < N; ++v) {
    for (auto i = static_cast<std::size_t>(off[v]); i < static_cast<std::size_t>(off[v + 1]); ++i) {
      const auto pnode = par[i];
      const auto e = via[i];
      if (pnode < 0 || static_cast<std::size_t>(pnode) >= N) {
        throw std::invalid_argument("static path bundle: parent node id out of range");
      }
      if (e < 0 || static_cast<std::size_t>(e) >= E) {
        throw std::invalid_argument("static path bundle: via edge id out of range");
      }
      if (esrc[static_cast<std::size_t>(e)] != pnode ||
          edst[static_cast<std::size_t>(e)] != static_cast<NodeId>(v)) {
        throw std::invalid_argument(
            "static path bundle: via edge does not connect its parent to its node in "
            "this graph (bundle built for a different graph?)");
      }
    }
  }

  // Acyclicity (Kahn over parent -> child entries).
  {
    std::vector<int> indeg(N, 0);
    for (std::size_t v = 0; v < N; ++v) {
      indeg[v] = static_cast<int>(off[v + 1] - off[v]);
    }
    std::vector<NodeId> stack;
    stack.reserve(N);
    for (std::size_t v = 0; v < N; ++v) {
      if (indeg[v] == 0) stack.push_back(static_cast<NodeId>(v));
    }
    std::size_t seen = 0;
    // Child adjacency: parent -> list of children, derived on the fly.
    std::vector<std::vector<NodeId>> children(N);
    for (std::size_t v = 0; v < N; ++v) {
      for (auto i = static_cast<std::size_t>(off[v]); i < static_cast<std::size_t>(off[v + 1]); ++i) {
        children[static_cast<std::size_t>(par[i])].push_back(static_cast<NodeId>(v));
      }
    }
    while (!stack.empty()) {
      auto u = stack.back(); stack.pop_back();
      ++seen;
      for (auto v : children[static_cast<std::size_t>(u)]) {
        if (--indeg[static_cast<std::size_t>(v)] == 0) stack.push_back(v);
      }
    }
    if (seen != N) {
      throw std::invalid_argument("static path bundle: predecessor structure contains a cycle");
    }
  }

  // Structural src->dst connectivity ignoring masks: a bundle that never had an
  // src->dst walk is malformed for this demand, distinct from being down.
  auto backward_reaches_src = [&](auto&& admit_entry) {
    std::vector<char> reach(N, 0);
    std::vector<NodeId> bfs;
    bfs.push_back(dst);
    reach[static_cast<std::size_t>(dst)] = 1;
    for (std::size_t head = 0; head < bfs.size(); ++head) {
      const auto v = static_cast<std::size_t>(bfs[head]);
      for (auto i = static_cast<std::size_t>(off[v]); i < static_cast<std::size_t>(off[v + 1]); ++i) {
        if (!admit_entry(i)) continue;
        const auto pnode = par[i];
        if (!reach[static_cast<std::size_t>(pnode)]) {
          reach[static_cast<std::size_t>(pnode)] = 1;
          bfs.push_back(pnode);
        }
      }
    }
    return reach[static_cast<std::size_t>(src)] != 0;
  };
  if (!backward_reaches_src([](std::size_t) { return true; })) {
    throw std::invalid_argument("static path bundle: no src->dst walk exists in the bundle");
  }

  // Prune against the policy's masks (failure exclusions), then re-check
  // connectivity: no surviving walk means the LSP is DOWN.
  const bool use_nm = !node_mask.empty();
  const bool use_em = !edge_mask.empty();
  auto entry_up = [&](std::size_t i) {
    const auto pnode = static_cast<std::size_t>(par[i]);
    const auto e = static_cast<std::size_t>(via[i]);
    if (use_em && !edge_mask[e]) return false;
    if (use_nm && !node_mask[pnode]) return false;
    return true;
  };
  auto node_up = [&](NodeId v) { return !use_nm || node_mask[static_cast<std::size_t>(v)]; };
  if (!node_up(src) || !node_up(dst)) return std::nullopt;

  PredDAG pruned;
  pruned.parent_offsets.assign(N + 1, 0);
  for (std::size_t v = 0; v < N; ++v) {
    std::int32_t c = 0;
    if (node_up(static_cast<NodeId>(v))) {
      for (auto i = static_cast<std::size_t>(off[v]); i < static_cast<std::size_t>(off[v + 1]); ++i) {
        if (entry_up(i)) ++c;
      }
    }
    pruned.parent_offsets[v + 1] = pruned.parent_offsets[v] + c;
  }
  pruned.parents.resize(static_cast<std::size_t>(pruned.parent_offsets.back()));
  pruned.via_edges.resize(static_cast<std::size_t>(pruned.parent_offsets.back()));
  for (std::size_t v = 0; v < N; ++v) {
    auto w = static_cast<std::size_t>(pruned.parent_offsets[v]);
    if (!node_up(static_cast<NodeId>(v))) continue;
    for (auto i = static_cast<std::size_t>(off[v]); i < static_cast<std::size_t>(off[v + 1]); ++i) {
      if (!entry_up(i)) continue;
      pruned.parents[w] = par[i];
      pruned.via_edges[w] = via[i];
      ++w;
    }
  }

  const auto& poff = pruned.parent_offsets;
  {
    std::vector<char> reach(N, 0);
    std::vector<NodeId> bfs;
    bfs.push_back(dst);
    reach[static_cast<std::size_t>(dst)] = 1;
    for (std::size_t head = 0; head < bfs.size(); ++head) {
      const auto v = static_cast<std::size_t>(bfs[head]);
      for (auto i = static_cast<std::size_t>(poff[v]); i < static_cast<std::size_t>(poff[v + 1]); ++i) {
        const auto pnode = pruned.parents[i];
        if (!reach[static_cast<std::size_t>(pnode)]) {
          reach[static_cast<std::size_t>(pnode)] = 1;
          bfs.push_back(pnode);
        }
      }
    }
    if (!reach[static_cast<std::size_t>(src)]) return std::nullopt;  // DOWN
  }

  // Cost: min-cost src->dst walk within the PRUNED bundle, so the flow reports the
  // metric of a walk it can actually use. DP in a topological order of the pruned
  // DAG (acyclicity established above; pruning cannot introduce cycles).
  Cost best = std::numeric_limits<Cost>::max();
  {
    const auto costs = g.cost_view();
    std::vector<Cost> dist(N, std::numeric_limits<Cost>::max());
    dist[static_cast<std::size_t>(src)] = 0;
    std::vector<int> indeg(N, 0);
    std::vector<std::vector<std::pair<NodeId, EdgeId>>> children(N);
    for (std::size_t v = 0; v < N; ++v) {
      indeg[v] = static_cast<int>(poff[v + 1] - poff[v]);
      for (auto i = static_cast<std::size_t>(poff[v]); i < static_cast<std::size_t>(poff[v + 1]); ++i) {
        children[static_cast<std::size_t>(pruned.parents[i])].emplace_back(
            static_cast<NodeId>(v), pruned.via_edges[i]);
      }
    }
    std::vector<NodeId> stack;
    for (std::size_t v = 0; v < N; ++v) {
      if (indeg[v] == 0) stack.push_back(static_cast<NodeId>(v));
    }
    while (!stack.empty()) {
      auto u = stack.back(); stack.pop_back();
      const auto du = dist[static_cast<std::size_t>(u)];
      for (auto [v, e] : children[static_cast<std::size_t>(u)]) {
        if (du != std::numeric_limits<Cost>::max()) {
          const Cost cand = du + costs[static_cast<std::size_t>(e)];
          if (cand < dist[static_cast<std::size_t>(v)]) dist[static_cast<std::size_t>(v)] = cand;
        }
        if (--indeg[static_cast<std::size_t>(v)] == 0) stack.push_back(v);
      }
    }
    best = dist[static_cast<std::size_t>(dst)];
  }
  return std::make_optional(std::make_pair(std::move(pruned), best));
}

} // namespace

/* See flow_policy.hpp for the contract. */
void FlowPolicy::set_static_paths(NodeId src, NodeId dst, std::vector<PredDAG> bundles) {
  const auto* g = ctx_.graph.graph.get();
  if (g == nullptr) {
    throw std::invalid_argument("FlowPolicy::set_static_paths: policy has no graph");
  }
  if (!flows_.empty()) {
    throw std::invalid_argument(
        "FlowPolicy::set_static_paths: policy already holds flows; call remove_demand() first");
  }
  if (bundles.empty()) {
    throw std::invalid_argument("FlowPolicy::set_static_paths: bundles must be non-empty");
  }
  if (shortest_path_) {
    throw std::invalid_argument(
        "FlowPolicy::set_static_paths: incompatible with shortest_path=true (single-"
        "augmentation IP semantics contradict pinned multi-LSP placement)");
  }
  const auto N = g->num_nodes();
  if (src < 0 || src >= N || dst < 0 || dst >= N || src == dst) {
    throw std::invalid_argument("FlowPolicy::set_static_paths: src/dst out of range or equal");
  }
  // Validate the USER's max_flow_count against the supplied bundle count. The
  // configured value is never overwritten: the usable (up) count U drives the
  // per-flow math at placement time, and N - flow_count() is the down-LSP count.
  if (max_flow_count_.has_value() && *max_flow_count_ != static_cast<int>(bundles.size())) {
    throw std::invalid_argument("If set, max_flow_count must be equal to the number of static paths.");
  }

  std::vector<StaticBundle> usable;
  usable.reserve(bundles.size());
  for (auto const& dag : bundles) {
    if (auto prepared = prepare_static_bundle(*g, dag, src, dst, node_mask_, edge_mask_)) {
      usable.push_back(StaticBundle{std::move(prepared->first), prepared->second});
    }
  }
  static_src_ = src;
  static_dst_ = dst;
  static_supplied_count_ = static_cast<int>(bundles.size());
  static_bundles_ = std::move(usable);
}

} // namespace netgraph::core
