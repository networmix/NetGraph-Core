/* Path enumeration from PredDAG (resolve_to_paths). */
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/profiling.hpp"

#include <algorithm>
#include <optional>
#include <stack>
#include <utility>
#include <vector>

namespace netgraph::core {

static inline void group_parents(const PredDAG& dag, NodeId v,
                                 std::vector<std::pair<NodeId, std::vector<EdgeId>>>& out) {
  out.clear();
  const auto start = static_cast<std::size_t>(dag.parent_offsets[static_cast<std::size_t>(v)]);
  const auto end   = static_cast<std::size_t>(dag.parent_offsets[static_cast<std::size_t>(v) + 1]);
  if (start >= end) return;
  // Collect entries, grouping by parent id
  // PredDAG does not guarantee sorted by parent; aggregate via map-like linear scan
  for (std::size_t i = start; i < end; ++i) {
    auto p = dag.parents[i];
    auto e = dag.via_edges[i];
    // find or append group
    bool found = false;
    for (auto& pr : out) {
      if (pr.first == p) { pr.second.push_back(e); found = true; break; }
    }
    if (!found) out.emplace_back(p, std::vector<EdgeId>{e});
  }
}

std::vector<std::vector<std::pair<NodeId, std::vector<EdgeId>>>>
resolve_to_paths(const PredDAG& dag, NodeId src, NodeId dst,
                 bool split_parallel_edges,
                 std::optional<std::int64_t> max_paths) {
  std::vector<std::vector<std::pair<NodeId, std::vector<EdgeId>>>> paths;
  if (src == dst) {
    // Trivial path: ((src, ()))
    std::vector<std::pair<NodeId, std::vector<EdgeId>>> p;
    p.emplace_back(src, std::vector<EdgeId>{});
    paths.push_back(std::move(p));
    return paths;
  }
  if (static_cast<std::size_t>(dst) >= dag.parent_offsets.size() - 1) return paths;
  if (dag.parent_offsets[static_cast<std::size_t>(dst)] == dag.parent_offsets[static_cast<std::size_t>(dst) + 1]) return paths;

  // Iterative DFS stack: each frame holds current node and index into its parent-groups.
  struct Frame { NodeId node; std::size_t idx; std::vector<std::pair<NodeId, std::vector<EdgeId>>> groups; };
  std::vector<Frame> stack;
  stack.reserve(16);
  // start from dst
  Frame start; start.node = dst; start.idx = 0; group_parents(dag, dst, start.groups);
  stack.push_back(std::move(start));

  std::vector<std::pair<NodeId, std::vector<EdgeId>>> current; // reversed path accum

  while (!stack.empty()) {
    auto& top = stack.back();
    if (top.idx >= top.groups.size()) {
      // backtrack
      stack.pop_back();
      if (!current.empty()) current.pop_back();
      continue;
    }
    auto [parent, edges] = top.groups[top.idx++];
    current.emplace_back(top.node, std::move(edges));
    if (parent == src) {
      // reached src; build forward segments: for each hop prev->next store (next, edges)
      std::vector<std::pair<NodeId, std::vector<EdgeId>>> segments;
      segments.reserve(current.size());
      for (auto it = current.rbegin(); it != current.rend(); ++it) {
        segments.emplace_back(it->first, it->second);
      }
      // Build path tuples: (src, edges for src->n1), (n1, edges for n1->n2), ..., (dst, ())
      std::vector<std::pair<NodeId, std::vector<EdgeId>>> path;
      path.reserve(segments.size() + 1);
      if (!segments.empty()) {
        // src element
        path.emplace_back(src, segments[0].second);
        // intermediate elements (attach next hop's edges to current node)
        for (std::size_t j = 1; j < segments.size(); ++j) {
          path.emplace_back(segments[j - 1].first, segments[j].second);
        }
        // dst element
        path.emplace_back(segments.back().first, std::vector<EdgeId>{});
      } else {
        // Degenerate: src==dst should be handled earlier, but keep form
        path.emplace_back(src, std::vector<EdgeId>{});
      }
      if (!split_parallel_edges) {
        paths.push_back(std::move(path));
      } else {
        // expand cartesian product over edge sets excluding last (dst has empty edges)
        // collect ranges
        std::vector<std::size_t> idxs(path.size(), 0);
        // Enumerate over all elements except the final dst (which has empty edges)
        const std::size_t start_i = 0;
        const std::size_t end_i = path.size() - 2; // last index before dst
        // initialize counters
        bool done = false;
        while (!done) {
          // build one concrete path
          std::vector<std::pair<NodeId, std::vector<EdgeId>>> concrete;
          concrete.reserve(path.size());
          for (std::size_t i = start_i; i <= end_i; ++i) {
            const auto& node = path[i].first;
            const auto& eds = path[i].second;
            if (!eds.empty()) {
              std::size_t sel = std::min(idxs[i], eds.size() - 1);
              concrete.emplace_back(node, std::vector<EdgeId>{ eds[sel] });
            } else {
              concrete.emplace_back(node, std::vector<EdgeId>{});
            }
          }
          // append dst with empty edges
          concrete.emplace_back(path.back().first, std::vector<EdgeId>{});
          paths.push_back(std::move(concrete));
          if (max_paths && static_cast<std::int64_t>(paths.size()) >= *max_paths) return paths;
          // increment counters (mixed radix)
          std::ptrdiff_t k = static_cast<std::ptrdiff_t>(end_i);
          while (k >= static_cast<std::ptrdiff_t>(start_i)) {
            if (path[static_cast<std::size_t>(k)].second.empty()) { --k; continue; }
            idxs[static_cast<std::size_t>(k)]++;
            if (idxs[static_cast<std::size_t>(k)] < path[static_cast<std::size_t>(k)].second.size()) break;
            idxs[static_cast<std::size_t>(k)] = 0;
            if (k == static_cast<std::ptrdiff_t>(start_i)) { done = true; break; }
            --k;
          }
          if (k < static_cast<std::ptrdiff_t>(start_i)) done = true;
        }
      }
      if (max_paths && static_cast<std::int64_t>(paths.size()) >= *max_paths) return paths;
      current.pop_back();
      continue;
    }
    // descend
    Frame next; next.node = parent; next.idx = 0; group_parents(dag, parent, next.groups);
    if (next.groups.empty()) {
      // dead end, backtrack
      current.pop_back();
      continue;
    }
    stack.push_back(std::move(next));
  }

  return paths;
}

} // namespace netgraph::core

/*
  Dijkstra shortest path algorithm with capacity-aware tie-breaking.

  Features:
    - Residual-aware traversal: uses dynamic residuals if provided, otherwise static capacity
    - Multipath mode: collects all equal-cost predecessor edges per node
    - Single-path mode: selects one edge per adjacency using:
      * Edge-level tie-breaking for parallel edges (PreferHigherResidual or Deterministic)
      * Node-level tie-breaking for equal-cost nodes (prefers higher bottleneck capacity)
    - Early exit when specific destination is reached
*/
#include "netgraph/core/shortest_paths.hpp"
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <queue>
#include <tuple>
#include <utility>
#include <vector>
#include "netgraph/core/constants.hpp"

namespace netgraph::core {

namespace {
static std::pair<std::vector<Cost>, PredDAG>
shortest_paths_core(const StrictMultiDiGraph& g, NodeId src,
                    std::optional<NodeId> dst,
                    bool multipath_arg,
                    const EdgeSelection& selection,
                    std::span<const Cap> residual,
                    std::span<const bool> node_mask,
                    std::span<const bool> edge_mask) {
  NGRAPH_PROFILE_SCOPE("shortest_paths_core");
  const auto N = g.num_nodes();
  const auto row = g.row_offsets_view();
  const auto col = g.col_indices_view();
  const auto aei = g.adj_edge_index_view();
  const auto cost = g.cost_view();
  const auto cap  = g.capacity_view();

  // Initialize distance array to infinity (max value).
  std::vector<Cost> dist(static_cast<std::size_t>(N), std::numeric_limits<Cost>::max());

  // Track minimum residual capacity along the path to each node (for node-level tie-breaking).
  // When distances are equal, prefer paths with higher bottleneck capacity.
  std::vector<Cap> min_residual_to_node(static_cast<std::size_t>(N), static_cast<Cap>(0));

  const bool use_node_mask = (node_mask.size() == static_cast<std::size_t>(g.num_nodes()));
  const bool use_edge_mask = (edge_mask.size() == static_cast<std::size_t>(g.num_edges()));
  const bool src_allowed = (src >= 0 && src < N && (!use_node_mask || node_mask[static_cast<std::size_t>(src)]));
  if (src_allowed) {
    dist[static_cast<std::size_t>(src)] = static_cast<Cost>(0);
    min_residual_to_node[static_cast<std::size_t>(src)] = std::numeric_limits<Cap>::max();
  }

  // pred_lists[v] stores predecessors for node v as (parent_node, [edges_from_parent]).
  // In multipath mode, multiple parents with equal-cost paths are retained.
  std::vector<std::vector<std::pair<NodeId, std::vector<EdgeId>>>> pred_lists(static_cast<std::size_t>(N));
  if (src_allowed) {
    pred_lists[static_cast<std::size_t>(src)] = {};
  } else {
    // Source is out of range or masked out: no traversal, return empty DAG.
    PredDAG dag;
    dag.parent_offsets.assign(static_cast<std::size_t>(N + 1), 0);
    return {std::move(dist), std::move(dag)};
  }

  // Priority queue for Dijkstra with capacity-aware node-level tie-breaking.
  // QItem is (cost, -residual, node). Negated residual ensures higher capacity = higher priority.
  // Lexicographic ordering: cost (minimize) -> residual (maximize) -> node (deterministic).
  // This naturally distributes flows across equal-cost paths based on available capacity.
  using QItem = std::tuple<Cost, Cap, NodeId>;
  auto cmp = [](const QItem& a, const QItem& b) { return a > b; };
  std::priority_queue<QItem, std::vector<QItem>, decltype(cmp)> pq(cmp);
  pq.emplace(static_cast<Cost>(0), -std::numeric_limits<Cap>::max(), src);
  Cost best_dst_cost = std::numeric_limits<Cost>::max();
  bool have_best_dst = false;
  const bool early_exit = dst.has_value();
  const NodeId dst_node = dst.value_or(-1);

  const bool has_residual = (residual.size() == static_cast<std::size_t>(g.num_edges()));
  const bool require_cap = selection.require_capacity || has_residual;
  const bool multipath = multipath_arg;

  while (!pq.empty()) {
    // Extract min-cost node from priority queue.
    // Structured binding: auto [d_u, neg_res_u, u] = ... destructures the tuple.
    auto [d_u, neg_res_u, u] = pq.top(); pq.pop();
    if (u < 0 || u >= N) continue;
    // Skip stale entries (node already processed at a lower cost).
    if (d_u > dist[static_cast<std::size_t>(u)]) continue;
    // Skip residual-stale entries in single-path mode (same cost but outdated residual).
    if (!multipath && d_u == dist[static_cast<std::size_t>(u)] &&
        -neg_res_u < min_residual_to_node[static_cast<std::size_t>(u)] - kEpsilon) continue;

    // Early exit optimization: record when we first reach destination.
    if (early_exit && u == dst_node && !have_best_dst) { best_dst_cost = d_u; have_best_dst = true; }
    if (early_exit && u == dst_node) {
      if (pq.empty() || std::get<0>(pq.top()) > best_dst_cost) break; else continue;
    }

    // Iterate over u's outgoing edges using CSR row offsets.
    auto start = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto end   = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
    std::size_t i = start;
    // Process edges grouped by destination node.
    // Multigraph may have multiple edges (u, v); CSR clusters them together.
    while (i < end) {
      NodeId v = col[i];
      // Skip masked nodes (skip entire neighbor group).
      if (use_node_mask && !node_mask[static_cast<std::size_t>(v)]) {
        std::size_t j_skip = i; while (j_skip < end && col[j_skip] == v) ++j_skip; i = j_skip; continue;
      }

      // Select best edge(s) from u to v according to policy.
      Cost min_edge_cost = std::numeric_limits<Cost>::max();
      std::vector<EdgeId> selected_edges;
      double best_rem_for_min_cost = -1.0;
      std::size_t j = i;
      int best_edge_id = -1;

      // Scan all parallel edges from u to v (they are consecutive in CSR).
      for (; j < end && col[j] == v; ++j) {
        auto e = static_cast<std::size_t>(aei[j]);
        if (use_edge_mask && !edge_mask[e]) continue;
        const Cap rem = has_residual ? residual[e] : cap[e];
        if (require_cap && rem < kMinCap) continue;
        const Cost ecost = static_cast<Cost>(cost[e]);
        if (ecost < min_edge_cost) {
          min_edge_cost = ecost;
          selected_edges.clear();
          if (selection.multi_edge) {
            selected_edges.push_back(static_cast<EdgeId>(aei[j]));
          } else {
            best_edge_id = static_cast<int>(e);
            best_rem_for_min_cost = static_cast<double>(rem);
          }
        } else if (ecost == min_edge_cost) {
          if (selection.multi_edge) {
            selected_edges.push_back(static_cast<EdgeId>(aei[j]));
          } else {
            // tie-break among equal-cost edges
            if (selection.tie_break == EdgeTieBreak::PreferHigherResidual) {
              if (static_cast<double>(rem) > best_rem_for_min_cost + kEpsilon) {
                best_edge_id = static_cast<int>(e);
                best_rem_for_min_cost = static_cast<double>(rem);
              } else if (std::abs(static_cast<double>(rem) - best_rem_for_min_cost) <= kEpsilon) {
                // further tie-break deterministically by smaller edge id
                if (best_edge_id < 0 || static_cast<int>(e) < best_edge_id) {
                  best_edge_id = static_cast<int>(e);
                }
              }
            } else {
              // Deterministic: smallest edge id
              if (best_edge_id < 0 || static_cast<int>(e) < best_edge_id) {
                best_edge_id = static_cast<int>(e);
              }
            }
          }
        }
      }
      if (!selection.multi_edge && best_edge_id >= 0) {
        selected_edges.clear();
        selected_edges.push_back(static_cast<EdgeId>(best_edge_id));
      }
      // Update distance and predecessors if we found a better path (or equal-cost with better capacity).
      if (!selected_edges.empty()) {
        Cost new_cost = static_cast<Cost>(d_u + min_edge_cost);
        auto v_idx = static_cast<std::size_t>(v);

        // Compute bottleneck capacity along path to v through u for node-level tie-breaking.
        // Uses dynamic residuals if provided, otherwise static capacity.
        // This allows tie-breaking by capacity even in cost-only routing mode.
        Cap max_edge_residual = static_cast<Cap>(0);
        for (auto edge_id : selected_edges) {
          const Cap rem = has_residual ? residual[static_cast<std::size_t>(edge_id)]
                                        : cap[static_cast<std::size_t>(edge_id)];
          if (rem > max_edge_residual) max_edge_residual = rem;
        }
        Cap path_residual = std::min(min_residual_to_node[static_cast<std::size_t>(u)], max_edge_residual);

        // Relaxation: found shorter path to v, or equal-cost path with better capacity (single-path mode).
        if (new_cost < dist[v_idx] ||
            (!multipath && new_cost == dist[v_idx] &&
             path_residual > min_residual_to_node[v_idx] + kEpsilon)) {
          dist[v_idx] = new_cost;
          min_residual_to_node[v_idx] = path_residual;
          pred_lists[v_idx].clear();
          pred_lists[v_idx].push_back({u, std::move(selected_edges)});
          pq.emplace(new_cost, -path_residual, v);  // Negate residual for max-heap behavior
        }
        // Multipath: found equal-cost alternative path to v.
        else if (multipath && new_cost == dist[v_idx]) {
          pred_lists[v_idx].push_back({u, std::move(selected_edges)});
          // Note: In multipath mode, we don't update min_residual_to_node because
          // we're collecting all equal-cost paths, not choosing based on residual.
        }
      }
      i = j;  // Advance to next neighbor group
    }
    if (have_best_dst) { if (pq.empty() || std::get<0>(pq.top()) > best_dst_cost) break; }
  }

  // Convert pred_lists to PredDAG using CSR-like layout.
  // parent_offsets[v]:parent_offsets[v+1] gives the range in parents/via_edges for node v.
  PredDAG dag;
  dag.parent_offsets.assign(static_cast<std::size_t>(N+1), 0);

  // Step 1: Count total predecessor entries per node.
  for (std::int32_t v=0; v<N; ++v) {
    std::size_t c=0;
    for (auto const& pe : pred_lists[static_cast<std::size_t>(v)]) c += pe.second.size();
    dag.parent_offsets[static_cast<std::size_t>(v+1)] = static_cast<std::int32_t>(c);
  }

  // Step 2: Convert counts to cumulative offsets (prefix sum).
  for (std::size_t k=1; k<dag.parent_offsets.size(); ++k)
    dag.parent_offsets[k] += dag.parent_offsets[k-1];

  // Step 3: Fill parents and via_edges arrays.
  dag.parents.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  dag.via_edges.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  for (std::int32_t v=0; v<N; ++v) {
    auto base = static_cast<std::size_t>(dag.parent_offsets[static_cast<std::size_t>(v)]);
    std::size_t k = 0;
    for (auto const& pe : pred_lists[static_cast<std::size_t>(v)]) {
      NodeId p = pe.first;
      for (auto e : pe.second) {
        dag.parents[base+k] = p;
        dag.via_edges[base+k] = e;
        ++k;
      }
    }
  }
  return {std::move(dist), std::move(dag)};
}
} // namespace

std::pair<std::vector<Cost>, PredDAG>
shortest_paths(const StrictMultiDiGraph& g, NodeId src,
               std::optional<NodeId> dst,
               bool multipath,
               const EdgeSelection& selection,
               std::span<const Cap> residual,
               std::span<const bool> node_mask,
               std::span<const bool> edge_mask) {
  if (!node_mask.empty() && node_mask.size() != static_cast<std::size_t>(g.num_nodes())) {
    throw std::invalid_argument("shortest_paths: node_mask length mismatch");
  }
  if (!edge_mask.empty() && edge_mask.size() != static_cast<std::size_t>(g.num_edges())) {
    throw std::invalid_argument("shortest_paths: edge_mask length mismatch");
  }
  if (!residual.empty() && residual.size() != static_cast<std::size_t>(g.num_edges())) {
    throw std::invalid_argument("shortest_paths: residual length mismatch");
  }
  return shortest_paths_core(g, src, dst, multipath, selection, residual, node_mask, edge_mask);
}

} // namespace netgraph::core
