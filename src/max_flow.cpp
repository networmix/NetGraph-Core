/*
  calc_max_flow — incremental max-flow via repeated SPF tiers over residuals.

  Uses FlowState to track residuals and place flow along SPF predecessor DAGs.
  Accumulates total flow, optional per-edge flows, cost distribution, and
  derives a min-cut from the final residual graph.

  Algorithm: Iterative augmentation using Successive Shortest Path on residual graphs,
  pushing flow across full ECMP/WCMP DAGs at each step.
*/
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/flow_state.hpp"
#include "netgraph/core/constants.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <future>
#include <limits>
#include <memory>
#include <optional>
#include <queue>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace netgraph::core {

namespace {

std::size_t thread_budget_from_env(const char* env_name, std::size_t candidate_count) {
  if (candidate_count <= 1) {
    return 1;
  }

  const char* env = std::getenv(env_name);
  if (env != nullptr && env[0] != '\0') {
    char* end = nullptr;
    unsigned long parsed = std::strtoul(env, &end, 10);
    if (end != env) {
      if (parsed == 0ul) {
        return 1;
      }
      return std::min<std::size_t>(static_cast<std::size_t>(parsed), candidate_count);
    }
  }

  const auto hw = std::thread::hardware_concurrency();
  const auto budget = hw > 0 ? static_cast<std::size_t>(hw) : 1u;
  return std::min<std::size_t>(budget, candidate_count);
}

std::optional<std::pair<EdgeId, Flow>>
evaluate_sensitivity_candidate(const StrictMultiDiGraph& g,
                               NodeId src,
                               NodeId dst,
                               FlowPlacement placement,
                               bool shortest_path,
                               bool require_capacity,
                               Flow baseline_flow,
                               std::span<const bool> node_mask,
                               EdgeId eid,
                               bool* local_mask,
                               std::size_t mask_size) {
  local_mask[static_cast<std::size_t>(eid)] = false;

  auto [new_flow, _] = calc_max_flow(
      g, src, dst, placement,
      shortest_path,
      require_capacity,
      /*with_edge_flows=*/false,
      /*with_reachable=*/false,
      /*with_residuals=*/false,
      node_mask, std::span<const bool>(local_mask, mask_size));

  local_mask[static_cast<std::size_t>(eid)] = true;

  double delta = baseline_flow - new_flow;
  if (delta > kMinFlow) {
    return std::pair<EdgeId, Flow>{eid, delta};
  }
  return std::nullopt;
}

} // namespace

std::pair<Flow, FlowSummary>
calc_max_flow(const StrictMultiDiGraph& g, NodeId src, NodeId dst,
              FlowPlacement placement, bool shortest_path,
              bool require_capacity,
              bool with_edge_flows,
              bool with_reachable,
              bool with_residuals,
              std::span<const bool> node_mask, std::span<const bool> edge_mask) {
  FlowSummary summary;
  const auto N = g.num_nodes();
  if (src < 0 || src >= N || dst < 0 || dst >= N || src == dst) {
    return {0.0, std::move(summary)};
  }
  if (!node_mask.empty() && node_mask.size() != static_cast<std::size_t>(g.num_nodes())) {
    throw std::invalid_argument("calc_max_flow: node_mask length mismatch");
  }
  if (!edge_mask.empty() && edge_mask.size() != static_cast<std::size_t>(g.num_edges())) {
    throw std::invalid_argument("calc_max_flow: edge_mask length mismatch");
  }
  const bool use_node_mask = !node_mask.empty();
  const bool use_edge_mask = !edge_mask.empty();
  if ((use_node_mask && !node_mask[static_cast<std::size_t>(src)]) ||
      (use_node_mask && !node_mask[static_cast<std::size_t>(dst)])) {
    return {0.0, std::move(summary)};
  }
  // Use FlowState for residual tracking and per-edge placement
  FlowState fs(g);
  Flow total = static_cast<Flow>(0.0);
  std::vector<std::pair<Cost,Flow>> cost_dist; // (cost, flow)

  // Iterate tiers: SPF over current residual or costs only
  // require_capacity=true: Require edges to have capacity, exclude saturated links (SDN/TE behavior)
  // require_capacity=false: Routes based on costs only, ignore capacity (IP/IGP behavior)
  while (true) {
    EdgeSelection sel;
    sel.multi_edge = true;
    sel.require_capacity = require_capacity;
    sel.tie_break = EdgeTieBreak::Deterministic;
    auto [dist, dag] = shortest_paths(
        g, src, dst,
        /*multipath=*/true,
        sel,
        require_capacity ? fs.residual_view() : std::span<const Cap>{},
        use_node_mask ? node_mask : std::span<const bool>{},
        use_edge_mask ? edge_mask : std::span<const bool>{});

    // No path if t has no parents in DAG
    if (static_cast<std::size_t>(dst) >= dag.parent_offsets.size() - 1 ||
        dag.parent_offsets[static_cast<std::size_t>(dst)] == dag.parent_offsets[static_cast<std::size_t>(dst) + 1]) {
      break;
    }

    Cost path_cost = dist[static_cast<std::size_t>(dst)];
    Flow placed = fs.place_on_dag(src, dst, dag, std::numeric_limits<double>::infinity(), placement);

    if (placed < kMinFlow) break;
    total += placed;
    // Merge by exact cost (integer)
    bool merged = false;
    for (auto& pr : cost_dist) { if (pr.first == path_cost) { pr.second += placed; merged = true; break; } }
    if (!merged) cost_dist.emplace_back(path_cost, placed);
    if (shortest_path) break;
  }

  // Completion phase. The tier loop above augments only along forward SPF DAGs
  // and can never cancel an earlier placement, so it may terminate below the
  // true maximum (its own min-cut then contradicts total_flow). Finish with
  // BFS augmentation on the full residual graph, where traversing an edge
  // backwards returns previously placed flow. Only meaningful for Proportional
  // max-flow semantics: EqualBalanced models ECMP admission and
  // require_capacity=false models fixed-cost IP routing, both of which are
  // placement models rather than max-flow computations.
  if (placement == FlowPlacement::Proportional && !shortest_path && require_capacity) {
    const auto row = g.row_offsets_view();
    const auto col = g.col_indices_view();
    const auto aei = g.adj_edge_index_view();
    const auto irow = g.in_row_offsets_view();
    const auto icol = g.in_col_indices_view();
    const auto iaei = g.in_adj_edge_index_view();
    const auto costv = g.cost_view();
    const auto residual = fs.residual_view();
    const auto eflow = fs.edge_flow_view();
    std::vector<std::int32_t> prev_node(static_cast<std::size_t>(N));
    std::vector<EdgeId> prev_edge(static_cast<std::size_t>(N));
    std::vector<char> prev_fwd(static_cast<std::size_t>(N));
    std::vector<char> seen(static_cast<std::size_t>(N));
    std::vector<std::pair<EdgeId, Flow>> fwd_deltas, rev_deltas;
    while (true) {
      std::fill(seen.begin(), seen.end(), 0);
      std::queue<std::int32_t> bfs;
      seen[static_cast<std::size_t>(src)] = 1;
      bfs.push(src);
      bool found = false;
      while (!bfs.empty() && !found) {
        auto u = bfs.front(); bfs.pop();
        auto us = static_cast<std::size_t>(u);
        // Forward residual arcs u -> v.
        for (auto j = static_cast<std::size_t>(row[us]); j < static_cast<std::size_t>(row[us + 1]); ++j) {
          auto v = static_cast<std::size_t>(col[j]);
          auto eid = static_cast<std::size_t>(aei[j]);
          if (use_edge_mask && !edge_mask[eid]) continue;
          if (use_node_mask && !node_mask[v]) continue;
          if (seen[v] || residual[eid] <= kMinCap) continue;
          seen[v] = 1; prev_node[v] = u; prev_edge[v] = static_cast<EdgeId>(eid); prev_fwd[v] = 1;
          if (static_cast<NodeId>(v) == dst) { found = true; break; }
          bfs.push(static_cast<std::int32_t>(v));
        }
        if (found) break;
        // Reverse residual arcs u -> w along edges w -> u carrying flow.
        for (auto j = static_cast<std::size_t>(irow[us]); j < static_cast<std::size_t>(irow[us + 1]); ++j) {
          auto w = static_cast<std::size_t>(icol[j]);
          auto eid = static_cast<std::size_t>(iaei[j]);
          if (use_edge_mask && !edge_mask[eid]) continue;
          if (use_node_mask && !node_mask[w]) continue;
          if (seen[w] || eflow[eid] <= kMinFlow) continue;
          seen[w] = 1; prev_node[w] = u; prev_edge[w] = static_cast<EdgeId>(eid); prev_fwd[w] = 0;
          if (static_cast<NodeId>(w) == dst) { found = true; break; }
          bfs.push(static_cast<std::int32_t>(w));
        }
      }
      if (!found) break;
      // Bottleneck and marginal path cost (forward costs minus cancelled costs).
      double bottleneck = std::numeric_limits<double>::infinity();
      Cost path_cost = 0;
      for (auto v = dst; v != src; v = prev_node[static_cast<std::size_t>(v)]) {
        auto vs = static_cast<std::size_t>(v);
        auto eid = static_cast<std::size_t>(prev_edge[vs]);
        if (prev_fwd[vs]) {
          bottleneck = std::min(bottleneck, static_cast<double>(residual[eid]));
          path_cost += costv[eid];
        } else {
          bottleneck = std::min(bottleneck, static_cast<double>(eflow[eid]));
          path_cost -= costv[eid];
        }
      }
      if (!(bottleneck >= kMinFlow)) break;
      fwd_deltas.clear(); rev_deltas.clear();
      for (auto v = dst; v != src; v = prev_node[static_cast<std::size_t>(v)]) {
        auto vs = static_cast<std::size_t>(v);
        if (prev_fwd[vs]) fwd_deltas.emplace_back(prev_edge[vs], static_cast<Flow>(bottleneck));
        else rev_deltas.emplace_back(prev_edge[vs], static_cast<Flow>(bottleneck));
      }
      fs.apply_deltas(fwd_deltas, /*add=*/true);
      fs.apply_deltas(rev_deltas, /*add=*/false);
      total += static_cast<Flow>(bottleneck);
      bool merged2 = false;
      for (auto& pr : cost_dist) { if (pr.first == path_cost) { pr.second += bottleneck; merged2 = true; break; } }
      if (!merged2) cost_dist.emplace_back(path_cost, static_cast<Flow>(bottleneck));
    }
  }

  summary.total_flow = total;
  if (with_edge_flows) {
    auto ef = fs.edge_flow_view();
    summary.edge_flows.assign(ef.begin(), ef.end());
  }
  if (with_residuals) {
    auto res = fs.residual_view();
    summary.residual_capacity.resize(static_cast<std::size_t>(g.num_edges()));
    for (std::size_t i = 0; i < summary.residual_capacity.size(); ++i) {
      summary.residual_capacity[i] = static_cast<Cap>(res[i]);
    }
  }
  if (!cost_dist.empty()) {
    std::sort(cost_dist.begin(), cost_dist.end(), [](auto const& a, auto const& b){ return a.first < b.first; });
    summary.costs.reserve(cost_dist.size());
    summary.flows.reserve(cost_dist.size());
    for (auto const& pr : cost_dist) { summary.costs.push_back(pr.first); summary.flows.push_back(pr.second); }
  }
  // Min-Cut: Set of edges crossing from the reachable set (S) to the unreachable set (T).
  // Computed by finding all nodes reachable from src in the residual graph.
  if (total >= kMinFlow) {
    auto mc = fs.compute_min_cut(src, node_mask, edge_mask);
    summary.min_cut = mc;
    if (with_reachable) {
      summary.reachable_nodes.assign(static_cast<std::size_t>(g.num_nodes()), 0u);
      auto residual = fs.residual_view();
      auto eflows = fs.edge_flow_view();
      const bool reach_use_node_mask = use_node_mask;
      const bool reach_use_edge_mask = use_edge_mask;
      const auto N = static_cast<std::size_t>(g.num_nodes());
      std::vector<std::int32_t> stack;
      stack.reserve(N);
      stack.push_back(src);
      while (!stack.empty()) {
        auto n = static_cast<std::size_t>(stack.back());
        stack.pop_back();
        if (summary.reachable_nodes[n]) continue;
        if (reach_use_node_mask && !node_mask[n]) continue;
        summary.reachable_nodes[n] = 1u;
        auto ro = g.row_offsets_view();
        auto ci = g.col_indices_view();
        auto ae = g.adj_edge_index_view();
        auto s = static_cast<std::size_t>(ro[n]);
        auto e = static_cast<std::size_t>(ro[n+1]);
        for (std::size_t p = s; p < e; ++p) {
          auto v = static_cast<std::size_t>(ci[p]);
          auto eid = static_cast<std::size_t>(ae[p]);
          if (reach_use_edge_mask && !edge_mask[eid]) continue;
          if (reach_use_node_mask && !node_mask[v]) continue;
          if (residual[eid] > kMinCap && !summary.reachable_nodes[v]) {
            stack.push_back(static_cast<std::int32_t>(v));
          }
        }
        auto iro = g.in_row_offsets_view();
        auto ici = g.in_col_indices_view();
        auto iae = g.in_adj_edge_index_view();
        auto rs = static_cast<std::size_t>(iro[n]);
        auto re = static_cast<std::size_t>(iro[n+1]);
        for (std::size_t p = rs; p < re; ++p) {
          auto u = static_cast<std::size_t>(ici[p]);
          auto eid = static_cast<std::size_t>(iae[p]);
          if (reach_use_edge_mask && !edge_mask[eid]) continue;
          if (reach_use_node_mask && !node_mask[u]) continue;
          if (eflows[eid] > kMinFlow && !summary.reachable_nodes[u]) {
            stack.push_back(static_cast<std::int32_t>(u));
          }
        }
      }
    }
  }
  return {summary.total_flow, std::move(summary)};
}

std::vector<FlowSummary>
batch_max_flow(const StrictMultiDiGraph& g,
               const std::vector<std::pair<NodeId,NodeId>>& pairs,
               FlowPlacement placement, bool shortest_path,
               bool require_capacity,
               bool with_edge_flows,
               bool with_reachable,
               bool with_residuals,
               const std::vector<std::span<const bool>>& node_masks,
               const std::vector<std::span<const bool>>& edge_masks) {
  // Validate all masks upfront
  for (const auto& nm : node_masks) {
    if (!nm.empty() && nm.size() != static_cast<std::size_t>(g.num_nodes())) {
      throw std::invalid_argument("batch_max_flow: node_mask length mismatch");
    }
  }
  for (const auto& em : edge_masks) {
    if (!em.empty() && em.size() != static_cast<std::size_t>(g.num_edges())) {
      throw std::invalid_argument("batch_max_flow: edge_mask length mismatch");
    }
  }

  std::vector<FlowSummary> out(pairs.size());
  auto run_range = [&](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      auto pr = pairs[i];
      std::span<const bool> nm = (i < node_masks.size() ? node_masks[i] : std::span<const bool>{});
      std::span<const bool> em = (i < edge_masks.size() ? edge_masks[i] : std::span<const bool>{});
      auto [val, summary] = calc_max_flow(g, pr.first, pr.second,
                                          placement, shortest_path,
                                          require_capacity,
                                          with_edge_flows,
                                          with_reachable,
                                          with_residuals,
                                          nm, em);
      out[i] = std::move(summary);
    }
  };
  // Pairs are independent over an immutable graph; evaluate chunks in parallel.
  // Thread count from NGRAPH_CORE_BATCH_THREADS env or hardware concurrency.
  const auto thread_budget = thread_budget_from_env("NGRAPH_CORE_BATCH_THREADS", pairs.size());
  if (thread_budget <= 1) {
    run_range(0, pairs.size());
    return out;
  }
  // Claim pairs one at a time from a shared counter rather than handing each worker a
  // fixed contiguous chunk. Per-pair cost varies by orders of magnitude (a dense
  // cross-fabric pair versus a disconnected one), so static chunking leaves workers
  // idle whenever the expensive pairs happen to land in the same chunk.
  std::atomic<std::size_t> next_pair{0};
  auto run_dynamic = [&]() {
    for (;;) {
      const std::size_t i = next_pair.fetch_add(1, std::memory_order_relaxed);
      if (i >= pairs.size()) break;
      run_range(i, i + 1);
    }
  };
  std::vector<std::future<void>> futures;
  futures.reserve(thread_budget);
  for (std::size_t w = 0; w < thread_budget; ++w) {
    futures.emplace_back(std::async(std::launch::async, run_dynamic));
  }
  for (auto& f : futures) f.get();
  return out;
}

std::vector<std::pair<EdgeId, Flow>>
sensitivity_analysis(const StrictMultiDiGraph& g, NodeId src, NodeId dst,
                     FlowPlacement placement, bool shortest_path,
                     bool require_capacity,
                     std::span<const bool> node_mask,
                     std::span<const bool> edge_mask) {
  // Step 1: Baseline analysis to identify saturated edges
  // Uses shortest_path mode to match the routing semantics being analyzed.
  auto [baseline_flow, summary] = calc_max_flow(
      g, src, dst, placement,
      shortest_path,
      require_capacity,
      /*with_edge_flows=*/false,
      /*with_reachable=*/false,
      /*with_residuals=*/true,
      node_mask, edge_mask);

  if (baseline_flow < kMinFlow || summary.residual_capacity.empty()) {
    return {};
  }

  std::vector<EdgeId> candidates;
  // Identify saturated edges (residual <= kMinCap)
  // Only consider edges that are not already masked out
  for (EdgeId eid = 0; eid < static_cast<EdgeId>(summary.residual_capacity.size()); ++eid) {
    if (!edge_mask.empty() && !edge_mask[eid]) continue; // Already masked
    if (summary.residual_capacity[eid] <= kMinCap) {
      candidates.push_back(eid);
    }
  }

  if (candidates.empty()) {
    return {};
  }

  // Prepare mutable mask buffer
  const auto N = static_cast<std::size_t>(g.num_edges());
  std::unique_ptr<bool[]> test_mask_buf(new bool[N]);
  if (edge_mask.empty()) {
    std::fill(test_mask_buf.get(), test_mask_buf.get() + N, true);
  } else {
    std::copy(edge_mask.begin(), edge_mask.end(), test_mask_buf.get());
  }
  std::vector<std::pair<EdgeId, Flow>> results;
  results.reserve(candidates.size());

  // Step 2: Iterate candidates, testing flow reduction when each is removed
  const auto thread_budget = thread_budget_from_env("NGRAPH_CORE_SENSITIVITY_THREADS", candidates.size());
  if (thread_budget <= 1) {
    for (EdgeId eid : candidates) {
      auto maybe_result = evaluate_sensitivity_candidate(
          g, src, dst, placement,
          shortest_path, require_capacity,
          baseline_flow, node_mask,
          eid, test_mask_buf.get(), N);
      if (maybe_result.has_value()) {
        results.push_back(*maybe_result);
      }
    }
    return results;
  }

  const auto chunk_size = (candidates.size() + thread_budget - 1) / thread_budget;
  std::vector<std::future<std::vector<std::pair<EdgeId, Flow>>>> futures;
  futures.reserve(thread_budget);

  for (std::size_t begin = 0; begin < candidates.size(); begin += chunk_size) {
    const std::size_t end = std::min<std::size_t>(begin + chunk_size, candidates.size());
    futures.emplace_back(std::async(std::launch::async, [&, begin, end]() {
      std::unique_ptr<bool[]> local_mask_buf(new bool[N]);
      if (!edge_mask.empty()) {
        std::copy(edge_mask.begin(), edge_mask.end(), local_mask_buf.get());
      } else {
        std::fill(local_mask_buf.get(), local_mask_buf.get() + N, true);
      }

      std::vector<std::pair<EdgeId, Flow>> local_results;
      local_results.reserve(end - begin);
      for (std::size_t idx = begin; idx < end; ++idx) {
        auto maybe_result = evaluate_sensitivity_candidate(
            g, src, dst, placement,
            shortest_path, require_capacity,
            baseline_flow, node_mask,
            candidates[idx], local_mask_buf.get(), N);
        if (maybe_result.has_value()) {
          local_results.push_back(*maybe_result);
        }
      }
      return local_results;
    }));
  }

  for (auto& future : futures) {
    auto local_results = future.get();
    results.insert(results.end(), local_results.begin(), local_results.end());
  }

  return results;
}

} // namespace netgraph::core
