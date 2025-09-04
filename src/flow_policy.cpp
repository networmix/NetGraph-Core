/*
  FlowPolicy — orchestrates path selection and flow placement for a single demand.

  Responsibilities:
    - Select shortest-path predecessors (SPF) per policy via `shortest_paths`.
    - Create and track flows with their current DAG and path cost.
    - Place volumes using `FlowPlacement` (Proportional or EqualBalanced).
    - Optionally re-optimize flows and enforce path-cost constraints.

  Notes:
    - "Residual-aware" SPF considers current remaining capacity (or a per-flow
      threshold) when exploring edges.
    - When `selection.multipath == false`, we lock exactly one edge per (u,v)
      neighbor pair across all flows to keep single-path behavior consistent.
*/
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/constants.hpp"

#include <algorithm>
#include <deque>
#include <limits>
#include <optional>
#include <unordered_set>

namespace netgraph::core {

double FlowPolicy::placed_demand() const {
  double s = 0.0;
  for (auto const& kv : flows_) s += kv.second.placed_flow;
  return s;
}

/* Compute an SPF predecessor DAG and its destination cost under the current
   policy configuration. Honors optional static paths, residual/edge masks,
   and single-path edge locking for (u,v) neighbor groups. */
std::optional<std::pair<PredDAG, Cost>> FlowPolicy::get_path_bundle(const FlowGraph& fg,
                                                                    NodeId src, NodeId dst,
                                                                    std::optional<double> min_flow) {
  // Static path handling
  if (!static_paths_.empty()) {
    for (auto const& t : static_paths_) {
      if (std::get<0>(t) == src && std::get<1>(t) == dst) {
        return std::make_optional(std::make_pair(std::get<2>(t), std::get<3>(t)));
      }
    }
    return std::nullopt;
  }
  if (path_alg_ != PathAlg::SPF) return std::nullopt;
  // Use configured selection; EqualBalanced behavior is achieved by placement and rebalancing.
  EdgeSelection sel = selection_;
  const bool single_uv = !selection_.multipath; // interpret as: single edge per (u,v) group
  if (single_uv) {
    // Allow multiple parents in DAG, we'll compress to one edge per (u,v) parent group below
    sel.multipath = true;
  }
  // Require residual-aware SPF when policy needs residual OR when EqualBalanced
  // provides a per-flow threshold (min_flow). Equal-balanced seeding enforces
  // a minimum deliverable volume per flow when configured.
  const bool require_residual = (sel.require_capacity || (flow_placement_ == FlowPlacement::EqualBalanced && min_flow.has_value()));
  std::pair<std::vector<Cost>, PredDAG> res;
  if (require_residual) {
    // Build a combined edge mask if needed:
    // - threshold: enforce per-edge residual >= min_flow when provided
    // - locking: for Proportional, allow-only locked edges; for EqualBalanced, exclude previously used edges
    const auto residual = fg.residual_view();
    std::vector<unsigned char> em; const bool* edge_mask_ptr = nullptr;
    bool need_mask = min_flow.has_value() || (single_uv && !locked_uv_edge_.empty());
    if (need_mask) {
      em.assign(residual.size(), 1u);
      if (min_flow.has_value()) {
        double thr = *min_flow;
        for (std::size_t i=0;i<residual.size();++i) em[i] = static_cast<unsigned char>(static_cast<double>(residual[i]) >= thr);
      }
      if (single_uv && !locked_uv_edge_.empty()) {
        const auto& g = fg.graph();
        auto row = g.row_offsets_view();
        auto col = g.col_indices_view();
        auto aei = g.adj_edge_index_view();
        for (std::int32_t u = 0; u < g.num_nodes(); ++u) {
          auto s = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
          auto e = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
          std::size_t i = s;
          while (i < e) {
            auto v = static_cast<std::int32_t>(col[i]);
            std::uint64_t key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(u))<<32) | static_cast<std::uint32_t>(v);
            auto it = locked_uv_edge_.find(key);
            if (it != locked_uv_edge_.end()) {
              EdgeId used = it->second;
              std::size_t j = i;
              while (j < e && col[j] == v) {
                auto eid = static_cast<EdgeId>(aei[j]);
                if (flow_placement_ == FlowPlacement::Proportional) {
                  if (eid != used) em[static_cast<std::size_t>(eid)] = 0u; // allow-only locked
                } else { // EqualBalanced: exclude previously used edge to diversify flows
                  if (eid == used) em[static_cast<std::size_t>(eid)] = 0u;
                }
                ++j;
              }
              i = j;
            } else { std::size_t j = i; while (j < e && col[j] == v) ++j; i = j; }
          }
        }
      }
      edge_mask_ptr = reinterpret_cast<const bool*>(em.data());
    }
    res = shortest_paths(fg.graph(), src, dst, sel, residual, /*node_mask*/ nullptr, edge_mask_ptr);
  } else {
    // Non-residual SPF. Apply lock if needed for single-path
    if (single_uv && !locked_uv_edge_.empty()) {
      const auto& g = fg.graph();
      auto row = g.row_offsets_view();
      auto col = g.col_indices_view();
      auto aei = g.adj_edge_index_view();
      std::vector<unsigned char> em(static_cast<std::size_t>(g.num_edges()), 1u);
      for (std::int32_t u = 0; u < g.num_nodes(); ++u) {
        auto s = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
        auto e = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
        std::size_t i = s;
        while (i < e) {
          auto v = static_cast<std::int32_t>(col[i]);
          std::uint64_t key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(u))<<32) | static_cast<std::uint32_t>(v);
          auto it = locked_uv_edge_.find(key);
          if (it != locked_uv_edge_.end()) {
            EdgeId used = it->second;
            std::size_t j = i;
            while (j < e && col[j] == v) {
              auto eid = static_cast<EdgeId>(aei[j]);
              if (flow_placement_ == FlowPlacement::Proportional) {
                if (eid != used) em[static_cast<std::size_t>(eid)] = 0u; // allow-only locked
              } else { // EqualBalanced
                if (eid == used) em[static_cast<std::size_t>(eid)] = 0u; // exclude used
              }
              ++j;
            }
            i = j;
          } else { std::size_t j = i; while (j < e && col[j] == v) ++j; i = j; }
        }
      }
      res = shortest_paths(fg.graph(), src, dst, sel, std::span<const Cap>(), /*node_mask*/ nullptr, reinterpret_cast<const bool*>(em.data()));
    } else {
      res = shortest_paths(fg.graph(), src, dst, sel);
    }
  }
  const auto& dist = res.first;
  PredDAG dag = std::move(res.second);
  // If we want single edge per (u,v) group, compress the DAG accordingly while
  // keeping all distinct parents of each node.
  if (single_uv) {
    PredDAG ndag; ndag.parent_offsets.assign(dag.parent_offsets.size(), 0);
    std::vector<NodeId> nparents;
    std::vector<EdgeId> nvia;
    const std::vector<Cap>* residual_ptr = nullptr;
    std::vector<Cap> res_copy;
    if (require_residual) {
      // Snapshot residual
      auto rv = fg.residual_view();
      res_copy.assign(rv.begin(), rv.end());
      residual_ptr = &res_copy;
    }
    for (std::size_t v = 0; v + 1 < dag.parent_offsets.size(); ++v) {
      auto s = static_cast<std::size_t>(dag.parent_offsets[v]);
      auto e = static_cast<std::size_t>(dag.parent_offsets[v+1]);
      if (e <= s) { ndag.parent_offsets[v+1] = static_cast<std::int32_t>(nparents.size()); continue; }
      // Group by parent u
      std::unordered_map<NodeId, EdgeId> chosen;
      std::unordered_map<NodeId, double> best_resid;
      for (std::size_t i = s; i < e; ++i) {
        NodeId u = dag.parents[i];
        EdgeId eid = dag.via_edges[i];
        if (selection_.tie_break == EdgeTieBreak::PreferHigherResidual && residual_ptr) {
          double r = static_cast<double>((*residual_ptr)[static_cast<std::size_t>(eid)]);
          auto it = best_resid.find(u);
          if (it == best_resid.end() || r > it->second + 1e-18 || (std::abs(r - it->second) <= 1e-18 && eid < chosen[u])) {
            best_resid[u] = r; chosen[u] = eid;
          }
        } else {
          // Deterministic: keep smallest edge id per parent
          auto it = chosen.find(u);
          if (it == chosen.end() || eid < it->second) { chosen[u] = eid; }
        }
      }
      for (auto const& kv : chosen) { nparents.push_back(kv.first); nvia.push_back(kv.second); }
      ndag.parent_offsets[v+1] = static_cast<std::int32_t>(nparents.size());
    }
    ndag.parents = std::move(nparents);
    ndag.via_edges = std::move(nvia);
    dag = std::move(ndag);
  }
  if (dst < 0 || static_cast<std::size_t>(dst) >= dist.size()) return std::nullopt;
  Cost dst_cost = dist[static_cast<std::size_t>(dst)];
  if (best_path_cost_ == 0 || dst_cost < best_path_cost_) best_path_cost_ = dst_cost;
  // If policy is in shortest_path mode, disallow advancing to higher-cost tiers
  // within the same placement session. Only paths with cost equal to the best
  // discovered cost are allowed.
  if (shortest_path_ && dst_cost > best_path_cost_) {
    return std::nullopt;
  }
  if (max_path_cost_.has_value() || max_path_cost_factor_.has_value()) {
    double maxf = max_path_cost_factor_.value_or(1.0);
    Cost absmax = max_path_cost_.value_or(std::numeric_limits<Cost>::max());
    if (dst_cost > std::min<Cost>(absmax, static_cast<Cost>(static_cast<double>(best_path_cost_) * maxf))) return std::nullopt;
  }
  // Ensure there is at least one predecessor for dst
  if (static_cast<std::size_t>(dst) >= dag.parent_offsets.size()-1) return std::nullopt;
  if (dag.parent_offsets[static_cast<std::size_t>(dst)] == dag.parent_offsets[static_cast<std::size_t>(dst)+1]) return std::nullopt;
  // For single-path selection with proportional placement, require the
  // destination to have a unique min-cost predecessor; otherwise, bail in
  // ambiguous multi-parent cases.
  if (single_uv && flow_placement_ == FlowPlacement::Proportional) {
    // Re-run with multipath=true to detect multiple equal-cost predecessors of dst
    EdgeSelection probe = sel; probe.multipath = true;
    std::pair<std::vector<Cost>, PredDAG> probe_res;
    if (require_residual) {
      const auto residual = fg.residual_view();
      probe_res = shortest_paths(fg.graph(), src, dst, probe, residual);
    } else {
      probe_res = shortest_paths(fg.graph(), src, dst, probe);
    }
    auto s = static_cast<std::size_t>(probe_res.second.parent_offsets[static_cast<std::size_t>(dst)]);
    auto e = static_cast<std::size_t>(probe_res.second.parent_offsets[static_cast<std::size_t>(dst)+1]);
    // Count unique parent nodes among entries s..e-1
    std::unordered_set<NodeId> uniq_parents;
    for (std::size_t i = s; i < e; ++i) uniq_parents.insert(probe_res.second.parents[i]);
    if (uniq_parents.size() > 1) return std::nullopt;
  }
  // If single-path selection, record the chosen edge per (u,v) along destination parents
  if (!sel.multipath) {
    const auto& off = dag.parent_offsets;
    const auto& parents = dag.parents;
    const auto& vias = dag.via_edges;
    // For each node that has exactly one parent entry, we can lock (u->v)
    const auto& g = fg.graph(); (void)g;
    for (std::size_t v = 0; v + 1 < off.size(); ++v) {
      auto s = static_cast<std::size_t>(off[v]);
      auto e = static_cast<std::size_t>(off[v+1]);
      if (e == s + 1) {
        NodeId u = parents[s];
        EdgeId eid = vias[s];
        std::uint64_t key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(u))<<32) | static_cast<std::uint32_t>(static_cast<NodeId>(v));
        locked_uv_edge_[key] = eid;
      }
    }
  }
  return std::make_optional(std::make_pair(dag, dst_cost));
}

// Note: a previous `needs_residual` helper became redundant after inlining checks.

/* Create a new flow using the current path bundle. Returns nullptr if no
   admissible path is available given constraints. */
FlowRecord* FlowPolicy::create_flow(FlowGraph& fg, NodeId src, NodeId dst, std::int32_t flowClass,
                              std::optional<double> min_flow) {
  FlowIndex idx{src, dst, flowClass, next_flow_id_++};
  auto pb = get_path_bundle(fg, src, dst, min_flow);
  if (!pb.has_value()) return nullptr;
  auto [dag, cost] = std::move(pb.value());
  FlowRecord f(idx, src, dst, std::move(dag), cost);
  auto [it, ok] = flows_.emplace(idx, std::move(f));
  (void)ok;
  return &it->second;
}

/* Re-select a path for an existing flow, requesting at least (current+headroom)
   residual. On failure, restores the flow on its previous DAG. */
FlowRecord* FlowPolicy::reoptimize_flow(FlowGraph& fg, const FlowIndex& idx, double headroom) {
  auto it = flows_.find(idx);
  if (it == flows_.end()) return nullptr;
  FlowRecord& cur = it->second;
  const double current = cur.placed_flow;
  const double new_min = current + headroom;
  // Temporarily remove this flow
  fg.remove(idx);
  auto pb = get_path_bundle(fg, cur.src, cur.dst, new_min);
  if (!pb.has_value()) {
    // restore on original DAG
    Flow placed = fg.place(idx, cur.src, cur.dst, cur.dag, current, flow_placement_, false);
    cur.placed_flow = placed; // may be slightly less if capacity changed; acceptable
    return nullptr;
  }
  auto [dag, cost] = std::move(pb.value());
  cur.dag = std::move(dag);
  cur.cost = cost;
  Flow placed = fg.place(idx, cur.src, cur.dst, cur.dag, current, flow_placement_, false);
  cur.placed_flow = placed;
  return &cur;
}

/* Place `volume` of demand according to the policy. When `target_per_flow`
   is provided (e.g., during rebalancing), each flow aims for that target.
   Returns (total_placed, leftover). */
std::pair<double,double> FlowPolicy::place_demand(FlowGraph& fg,
                                                  NodeId src, NodeId dst,
                                                  std::int32_t flowClass,
                                                  double volume,
                                                  std::optional<double> target_per_flow,
                                                  std::optional<double> min_flow) {
  (void)min_flow; // reserved for future use; avoids unused parameter warning
  // prune missing flows: nothing to do at the ledger level; policies own flows
  double target = target_per_flow.value_or(volume);
  double per_target = target;
  // Behavior: do not seed equal-balanced flows when selection is
  // non-capacity-aware min-cost (ALL_MIN_COST) and no static paths are given.
  if (flow_placement_ == FlowPlacement::EqualBalanced && !selection_.require_capacity && static_paths_.empty()) {
    return { 0.0, volume };
  }
  // Allow K flows to reuse paths and to traverse multi-hop shortest tiers.
  if (flow_placement_ == FlowPlacement::EqualBalanced && max_flow_count_.has_value()) {
    // ECMP gating: in shortest_path mode, require at least two distinct
    // shortest-path next hops from src that can reach dst; otherwise place 0.
    if (shortest_path_) {
      EdgeSelection sel_ecmp = selection_;
      sel_ecmp.multipath = true;
      sel_ecmp.require_capacity = true;
      auto [dist0, dag0] = shortest_paths(fg.graph(), src, dst, sel_ecmp, fg.residual_view());
      const int N0 = fg.graph().num_nodes();
      // Build children adjacency and mark which nodes can reach dst via DAG
      std::vector<std::vector<int>> children0(static_cast<std::size_t>(N0));
      for (int v = 0; v < N0; ++v) {
        std::size_t s = static_cast<std::size_t>(dag0.parent_offsets[static_cast<std::size_t>(v)]);
        std::size_t e = static_cast<std::size_t>(dag0.parent_offsets[static_cast<std::size_t>(v+1)]);
        for (std::size_t i = s; i < e; ++i) {
          int u = static_cast<int>(dag0.parents[i]);
          children0[static_cast<std::size_t>(u)].push_back(v);
        }
      }
      std::vector<char> can_reach(static_cast<std::size_t>(N0), 0);
      if (dst >= 0 && dst < N0) {
        // reverse BFS on parent links to mark nodes that can reach dst
        std::deque<int> dq; dq.push_back(dst); can_reach[static_cast<std::size_t>(dst)] = 1;
        while (!dq.empty()) {
          int v = dq.front(); dq.pop_front();
          std::size_t s = static_cast<std::size_t>(dag0.parent_offsets[static_cast<std::size_t>(v)]);
          std::size_t e = static_cast<std::size_t>(dag0.parent_offsets[static_cast<std::size_t>(v+1)]);
          for (std::size_t i = s; i < e; ++i) {
            int u = static_cast<int>(dag0.parents[i]);
            if (!can_reach[static_cast<std::size_t>(u)]) { can_reach[static_cast<std::size_t>(u)] = 1; dq.push_back(u); }
          }
        }
      }
      int next_hops = 0;
      if (src >= 0 && src < N0) {
        for (int v : children0[static_cast<std::size_t>(src)]) {
          if (can_reach[static_cast<std::size_t>(v)]) ++next_hops;
        }
      }
      if (next_hops < 2) {
        return { 0.0, volume };
      }
    }
    // Allow K flows to reuse paths
    // Derive a reasonable per-flow target using request and src egress residual
    const auto& g = fg.graph();
    auto row = g.row_offsets_view();
    auto aei = g.adj_edge_index_view();
    auto residual = fg.residual_view();
    double src_cap = 0.0;
    if (src >= 0 && src < g.num_nodes()) {
      auto s = static_cast<std::size_t>(row[static_cast<std::size_t>(src)]);
      auto e = static_cast<std::size_t>(row[static_cast<std::size_t>(src)+1]);
      for (std::size_t j = s; j < e; ++j) {
        auto eid = static_cast<std::size_t>(aei[j]);
        src_cap += static_cast<double>(residual[eid]);
      }
    }
    double per_req = target / static_cast<double>(*max_flow_count_);
    double per_src = src_cap / static_cast<double>(*max_flow_count_);
    per_target = std::max(kMinFlow, std::min(per_req, per_src));
  }
  if (flows_.empty()) {
    if (!static_paths_.empty()) {
      if (max_flow_count_.has_value() && static_cast<int>(static_paths_.size()) != *max_flow_count_) {
        throw std::invalid_argument("If set, max_flow_count must be equal to the number of static paths.");
      }
      // Limit to matching src/dst
      for (auto const& t : static_paths_) {
        if (std::get<0>(t) == src && std::get<1>(t) == dst) {
          create_flow(fg, src, dst, flowClass, std::nullopt);
        } else {
          throw std::invalid_argument("Source and destination nodes of static paths do not match demand.");
        }
      }
    } else {
      int initial = min_flow_count_;
      if (flow_placement_ == FlowPlacement::EqualBalanced && max_flow_count_.has_value()) {
        // For equal-balanced placement, eagerly seed the requested number of flows
        // up to the configured maximum to ensure balanced distribution from the start.
        initial = std::min(min_flow_count_, *max_flow_count_);
      }
      for (int i=0;i<initial;++i) {
        // Seed flows using per-target as minimum per-edge residual requirement
        create_flow(fg, src, dst, flowClass,
                    (flow_placement_ == FlowPlacement::EqualBalanced && max_flow_count_.has_value()) ? std::optional<double>(per_target) : std::nullopt);
      }
    }
  }
  std::deque<FlowIndex> q;
  for (auto const& kv : flows_) q.push_back(kv.first);
  double total_placed = 0.0;
  int no_progress = 0;
  int iters = 0;
  // Diminishing-returns tracking
  std::deque<double> recent;
  const double initial_request = volume;
  while (volume >= kMinFlow && !q.empty()) {
    FlowIndex cur_idx = q.front(); q.pop_front();
    auto it_cur = flows_.find(cur_idx);
    if (it_cur == flows_.end()) continue;
    FlowRecord* f = &it_cur->second;
    // must have a DAG to place; skip otherwise
    if (f->dag.parent_offsets.empty()) { ++no_progress; if (no_progress>=max_no_progress_iterations_) break; continue; }
    // EqualBalanced with non-capacity-aware selection (ALL_MIN_COST) does not
    // place flow (guard to avoid unrealistic equal balancing without capacity awareness).
    if (flow_placement_ == FlowPlacement::EqualBalanced && !selection_.require_capacity && static_paths_.empty()) {
      ++no_progress; if (no_progress>=max_no_progress_iterations_) break; continue;
    }
    // For EqualBalanced (non-static), refresh DAG each round with per-flow target min_flow
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
      // During initial equal-balanced placement without per-target, allow each flow to place up to remaining volume.
      need = volume;
    } else {
      // Default behavior uses the global target amount.
      need = target;
    }
    const double request = std::min(need, volume);
    Flow placed = fg.place(f->index, f->src, f->dst, f->dag, request, flow_placement_, /*shortest_path*/ shortest_path_);
    f->placed_flow += placed;
    volume -= placed;
    total_placed += placed;
    ++iters;
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
    if (flow_placement_ == FlowPlacement::EqualBalanced && max_flow_count_.has_value()) {
      // Keep adding flows until we reach the configured maximum for equal-balanced.
      if (static_cast<int>(flows_.size()) < *max_flow_count_) {
        if (auto* nf = create_flow(fg, src, dst, flowClass, std::optional<double>(per_target))) q.push_back(nf->index);
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
    if (iters > max_total_iterations_) break;
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
                                                      std::int32_t flowClass,
                                                      double target_per_flow) {
  double vol = placed_demand();
  remove_demand(fg);
  return place_demand(fg, src, dst, flowClass, vol, target_per_flow, std::nullopt);
}

/* Remove all placed flows for this policy from the FlowGraph and reset
   per-flow placed volumes. */
void FlowPolicy::remove_demand(FlowGraph& fg) {
  for (auto const& kv : flows_) fg.remove(kv.first);
  for (auto& kv : flows_) kv.second.placed_flow = 0.0;
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
