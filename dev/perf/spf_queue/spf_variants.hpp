// Experimental copies of the production forward SPF loop (src/shortest_paths.cpp,
// shortest_paths_core) with the priority queue and the scratch workspace made
// pluggable. The relaxation, selection, tie-break and predecessor rules are kept
// verbatim so that any observable difference is attributable to the queue or
// the workspace policy. Research code only; not part of the library.
#pragma once

#include "netgraph/core/constants.hpp"
#include "netgraph/core/shortest_paths.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace spfx {

using namespace netgraph::core;
constexpr Cost INF = std::numeric_limits<Cost>::max();
using QItem = std::tuple<Cost, Cap, NodeId>;

// Reference queue: a binary heap over (cost, -bottleneck, node), the exact key
// std::priority_queue uses in production. Keys are unique per push (a node is
// re-pushed only with a strictly smaller cost or a strictly larger bottleneck), so
// any correct binary heap pops in the same order as std::priority_queue.
struct HeapQueue {
  static constexpr const char* name = "heap";
  std::vector<QItem> h;
  static bool gt(const QItem& a, const QItem& b) { return a > b; }
  void begin(Cost, Cost) { h.clear(); }
  void end() { h.clear(); }
  bool empty() const { return h.empty(); }
  Cost top_cost() { return std::get<0>(h.front()); }
  QItem pop() {
    std::pop_heap(h.begin(), h.end(), gt);
    QItem t = h.back(); h.pop_back(); return t;
  }
  void push(Cost c, Cap nr, NodeId v) {
    h.emplace_back(c, nr, v);
    std::push_heap(h.begin(), h.end(), gt);
  }
};

// Monotone bucket queue with ordered buckets. Dijkstra pops keys in
// non-decreasing cost order and every push is >= the cost being popped, so with
// W = max_cost/stride + 1 cyclic buckets each bucket holds exactly one cost value
// at any time. The active bucket is a binary heap on (-bottleneck, node); pushes
// into it (zero-cost edges) keep the heap property. The pop sequence is therefore
// identical to the global heap's (cost, -bottleneck, node) order, including
// stale-entry handling, zero-cost edges, single-path re-pushes and early exit.
template <bool kPow2 = false>
struct BucketHeapQueueT {
  static constexpr const char* name = kPow2 ? "bucket2" : "bucket";
  struct Item { Cap neg_res; NodeId node; };
  static bool gt(const Item& a, const Item& b) {
#ifdef SPFX_BAD_ORDER
    return a.node > b.node;   // teeth check: ignore bottleneck ordering
#else
    return a.neg_res != b.neg_res ? a.neg_res > b.neg_res : a.node > b.node;
#endif
  }
  std::vector<std::vector<Item>> buckets;
  std::vector<std::uint64_t> occ;
  std::size_t W = 0, cur_b = 0, pending = 0;
  Cost cur = 0, stride = 1;

  void begin(Cost max_cost, Cost stride_) {
    stride = stride_ > 0 ? stride_ : 1;
    W = static_cast<std::size_t>(max_cost / stride) + 1;
    if constexpr (kPow2) W = std::bit_ceil(W);   // modulo becomes a mask
    if (buckets.size() < W) buckets.resize(W);
    occ.assign((W + 63) / 64, 0);
    cur = 0; cur_b = 0; pending = 0;
  }
  void end() {
    // Early exit can leave items behind; clear only occupied buckets.
    for (std::size_t w = 0; w < occ.size(); ++w) {
      std::uint64_t word = occ[w];
      while (word) {
        const std::size_t b = (w << 6) + static_cast<std::size_t>(std::countr_zero(word));
        buckets[b].clear();
        word &= word - 1;
      }
      occ[w] = 0;
    }
    pending = 0;
  }
  bool empty() const { return pending == 0; }
  std::size_t bucket_of(Cost c) const {
    const Cost q = stride == 1 ? c : c / stride;
    if constexpr (kPow2) return static_cast<std::size_t>(q) & (W - 1);
    else return static_cast<std::size_t>(q % static_cast<Cost>(W));
  }
  void push(Cost c, Cap nr, NodeId v) {
    const std::size_t b = bucket_of(c);
    auto& vec = buckets[b];
    if (vec.empty()) occ[b >> 6] |= (std::uint64_t{1} << (b & 63));
    vec.push_back({nr, v});
    if (b == cur_b) std::push_heap(vec.begin(), vec.end(), gt);
    ++pending;
  }
  void settle_active() {
    if (!buckets[cur_b].empty()) return;
    std::size_t i = cur_b + 1; if (i == W) i = 0;
    std::size_t w = i >> 6;
    std::uint64_t word = occ[w] & (~std::uint64_t{0} << (i & 63));
    const std::size_t nwords = occ.size();
    while (!word) { w = (w + 1 == nwords) ? 0 : w + 1; word = occ[w]; }
    const std::size_t f = (w << 6) + static_cast<std::size_t>(std::countr_zero(word));
    const std::size_t steps = f > cur_b ? f - cur_b : f + W - cur_b;
    cur += static_cast<Cost>(steps) * stride; cur_b = f;
    std::make_heap(buckets[f].begin(), buckets[f].end(), gt);
  }
  Cost top_cost() { settle_active(); return cur; }
  QItem pop() {
    settle_active();
    auto& vec = buckets[cur_b];
    std::pop_heap(vec.begin(), vec.end(), gt);
    const Item it = vec.back(); vec.pop_back(); --pending;
    if (vec.empty()) occ[cur_b >> 6] &= ~(std::uint64_t{1} << (cur_b & 63));
    return {cur, it.neg_res, it.node};
  }
};
using BucketHeapQueue = BucketHeapQueueT<false>;
using BucketPow2Queue = BucketHeapQueueT<true>;

// Variant of the bucket queue that sorts a bucket once when it becomes active and
// consumes it linearly; same-cost pushes that arrive while the bucket is active
// (zero-cost edges) go to a small side heap, and pop takes the smaller head.
struct BucketSortQueue {
  static constexpr const char* name = "bsort";
  using Item = BucketHeapQueue::Item;
  static bool lt(const Item& a, const Item& b) {
    return a.neg_res != b.neg_res ? a.neg_res < b.neg_res : a.node < b.node;
  }
  static bool gt(const Item& a, const Item& b) { return lt(b, a); }
  std::vector<std::vector<Item>> buckets;
  std::vector<std::uint64_t> occ;
  std::vector<Item> side;      // heap of same-cost pushes into the active bucket
  std::size_t W = 0, cur_b = 0, pending = 0, head = 0;
  Cost cur = 0, stride = 1;

  void begin(Cost max_cost, Cost stride_) {
    stride = stride_ > 0 ? stride_ : 1;
    W = static_cast<std::size_t>(max_cost / stride) + 1;
    if (buckets.size() < W) buckets.resize(W);
    occ.assign((W + 63) / 64, 0);
    side.clear(); cur = 0; cur_b = 0; pending = 0; head = 0;
  }
  void end() {
    for (std::size_t w = 0; w < occ.size(); ++w) {
      std::uint64_t word = occ[w];
      while (word) {
        const std::size_t b = (w << 6) + static_cast<std::size_t>(std::countr_zero(word));
        buckets[b].clear();
        word &= word - 1;
      }
      occ[w] = 0;
    }
    buckets[cur_b].clear(); side.clear(); pending = 0; head = 0;
  }
  bool empty() const { return pending == 0; }
  std::size_t bucket_of(Cost c) const {
    return static_cast<std::size_t>((c / stride) % static_cast<Cost>(W));
  }
  bool active_empty() const { return head >= buckets[cur_b].size() && side.empty(); }
  void push(Cost c, Cap nr, NodeId v) {
    const std::size_t b = bucket_of(c);
    if (b == cur_b) {
      side.push_back({nr, v});
      std::push_heap(side.begin(), side.end(), gt);
    } else {
      auto& vec = buckets[b];
      if (vec.empty()) occ[b >> 6] |= (std::uint64_t{1} << (b & 63));
      vec.push_back({nr, v});
    }
    ++pending;
  }
  void settle_active() {
    if (!active_empty()) return;
    buckets[cur_b].clear(); head = 0;
    std::size_t i = cur_b + 1; if (i == W) i = 0;
    std::size_t w = i >> 6;
    std::uint64_t word = occ[w] & (~std::uint64_t{0} << (i & 63));
    const std::size_t nwords = occ.size();
    while (!word) { w = (w + 1 == nwords) ? 0 : w + 1; word = occ[w]; }
    const std::size_t f = (w << 6) + static_cast<std::size_t>(std::countr_zero(word));
    const std::size_t steps = f > cur_b ? f - cur_b : f + W - cur_b;
    cur += static_cast<Cost>(steps) * stride; cur_b = f;
    occ[f >> 6] &= ~(std::uint64_t{1} << (f & 63));   // active bucket is tracked by head/side
    auto& vec = buckets[f];
    std::sort(vec.begin(), vec.end(), lt);
  }
  Cost top_cost() { settle_active(); return cur; }
  QItem pop() {
    settle_active();
    auto& vec = buckets[cur_b];
    Item it;
    if (head < vec.size() && (side.empty() || lt(vec[head], side.front()))) {
      it = vec[head++];
    } else {
      std::pop_heap(side.begin(), side.end(), gt);
      it = side.back(); side.pop_back();
    }
    --pending;
    return {cur, it.neg_res, it.node};
  }
};

enum class WsMode { Fresh, ReuseFull, ReuseSparse, ReuseHybrid };
// ReuseHybrid tracks touched nodes like ReuseSparse but chooses the reset, output
// and DAG-conversion strategy at the end of the query from the touched count.

// Scratch arrays of the production loop. Fresh: allocate per query as production
// does (including reserve(E) on the entry arrays). ReuseFull: keep the vectors and
// reset all N entries. ReuseSparse: keep the vectors and reset only touched nodes.
struct Workspace {
  std::vector<Cost> dist;
  std::vector<Cap> minres;
  std::vector<std::int32_t> pred_head, pred_tail;
  std::vector<char> settled;
  std::vector<NodeId> ent_parent;
  std::vector<EdgeId> ent_edge;
  std::vector<std::int32_t> ent_next;
  std::vector<NodeId> touched;
  WsMode mode = WsMode::Fresh;

  bool tracking() const { return mode == WsMode::ReuseSparse || mode == WsMode::ReuseHybrid; }
  // Decisions taken after the search from the touched fraction (hybrid only).
  bool sparse_reset() const {
    return mode == WsMode::ReuseSparse || (mode == WsMode::ReuseHybrid && touched.size() * 4 < dist.size());
  }
  void begin(std::int32_t n, std::int32_t e) {
    const auto N = static_cast<std::size_t>(n);
    if (tracking()) {
      if (dist.size() != N) {
        dist.assign(N, INF); minres.assign(N, 0); pred_head.assign(N, -1);
        pred_tail.assign(N, -1); settled.assign(N, 0);
      }
      touched.clear();
    } else {
      dist.assign(N, INF); minres.assign(N, 0); pred_head.assign(N, -1);
      pred_tail.assign(N, -1); settled.assign(N, 0);
    }
    ent_parent.clear(); ent_edge.clear(); ent_next.clear();
    if (mode == WsMode::Fresh) {
      const auto E = static_cast<std::size_t>(e);
      ent_parent.reserve(E); ent_edge.reserve(E); ent_next.reserve(E);
    }
  }
  inline void touch(NodeId v) { if (tracking()) touched.push_back(v); }
  void end() {
    if (!tracking()) return;
    if (sparse_reset()) {
      for (auto v : touched) {
        const auto i = static_cast<std::size_t>(v);
        dist[i] = INF; minres[i] = 0; pred_head[i] = -1; pred_tail[i] = -1; settled[i] = 0;
      }
    } else {
      std::fill(dist.begin(), dist.end(), INF); std::fill(minres.begin(), minres.end(), 0);
      std::fill(pred_head.begin(), pred_head.end(), -1); std::fill(pred_tail.begin(), pred_tail.end(), -1);
      std::fill(settled.begin(), settled.end(), 0);
    }
    touched.clear();
  }
};

// Copy of shortest_paths_core with the queue and workspace injected.
template <class Queue>
std::pair<std::vector<Cost>, PredDAG>
spf_variant(const StrictMultiDiGraph& g, NodeId src,
            std::optional<NodeId> dst,
            bool multipath_arg,
            const EdgeSelection& selection,
            std::span<const Cap> residual,
            std::span<const bool> node_mask,
            std::span<const bool> edge_mask,
            Queue& pq, Workspace& ws, Cost max_cost, Cost stride) {
  const auto N = g.num_nodes();
  const auto row = g.row_offsets_view();
  const auto col = g.col_indices_view();
  const auto aei = g.adj_edge_index_view();
  const auto cost = g.cost_view();
  const auto cap  = g.capacity_view();

  ws.begin(N, g.num_edges());
  auto& dist = ws.dist;
  auto& min_residual_to_node = ws.minres;
  auto& pred_head = ws.pred_head;
  auto& pred_tail = ws.pred_tail;
  auto& ent_parent = ws.ent_parent;
  auto& ent_edge = ws.ent_edge;
  auto& ent_next = ws.ent_next;
  auto& settled = ws.settled;

  auto finish_dist = [&]() {
    std::vector<Cost> out;
    if (ws.mode == WsMode::Fresh) { out = std::move(dist); dist.clear(); }
    else if (ws.mode == WsMode::ReuseFull || !ws.sparse_reset()) { out = dist; }
    else {
      out.assign(static_cast<std::size_t>(N), INF);
      for (auto v : ws.touched) out[static_cast<std::size_t>(v)] = dist[static_cast<std::size_t>(v)];
    }
    return out;
  };

  const bool use_node_mask = (node_mask.size() == static_cast<std::size_t>(g.num_nodes()));
  const bool use_edge_mask = (edge_mask.size() == static_cast<std::size_t>(g.num_edges()));
  const bool src_allowed = (src >= 0 && src < N && (!use_node_mask || node_mask[static_cast<std::size_t>(src)]));
  if (src_allowed) {
    dist[static_cast<std::size_t>(src)] = static_cast<Cost>(0);
    min_residual_to_node[static_cast<std::size_t>(src)] = std::numeric_limits<Cap>::max();
    ws.touch(src);
  }

  auto pred_clear = [&](std::size_t v){ pred_head[v] = -1; pred_tail[v] = -1; };
  auto pred_append = [&](std::size_t v, NodeId p, EdgeId e){
    const auto idx = static_cast<std::int32_t>(ent_parent.size());
    ent_parent.push_back(p); ent_edge.push_back(e); ent_next.push_back(-1);
    if (pred_head[v] < 0) { pred_head[v] = idx; }
    else { ent_next[static_cast<std::size_t>(pred_tail[v])] = idx; }
    pred_tail[v] = idx;
  };
  if (!src_allowed) {
    PredDAG dag;
    dag.parent_offsets.assign(static_cast<std::size_t>(N + 1), 0);
    auto out = finish_dist();
    ws.end();
    return {std::move(out), std::move(dag)};
  }

  pq.begin(max_cost, stride);
  pq.push(static_cast<Cost>(0), -std::numeric_limits<Cap>::max(), src);
  Cost best_dst_cost = std::numeric_limits<Cost>::max();
  bool have_best_dst = false;
  const bool early_exit = dst.has_value();
  const NodeId dst_node = dst.value_or(-1);

  const bool has_residual = (residual.size() == static_cast<std::size_t>(g.num_edges()));
  const bool require_cap = selection.require_capacity || has_residual;
  const bool multipath = multipath_arg;

  std::vector<EdgeId> sel_buf; sel_buf.reserve(16);
  while (!pq.empty()) {
    auto [d_u, neg_res_u, u] = pq.pop();
    if (u < 0 || u >= N) continue;
    if (d_u > dist[static_cast<std::size_t>(u)]) continue;
    if (!multipath && d_u == dist[static_cast<std::size_t>(u)] &&
        -neg_res_u < min_residual_to_node[static_cast<std::size_t>(u)] - kEpsilon) continue;
    settled[static_cast<std::size_t>(u)] = 1;

    if (early_exit && u == dst_node && !have_best_dst) { best_dst_cost = d_u; have_best_dst = true; }
    if (early_exit && u == dst_node) {
      if (pq.empty() || pq.top_cost() > best_dst_cost) break; else continue;
    }

    auto start = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto end   = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
    std::size_t i = start;
    while (i < end) {
      NodeId v = col[i];
      if (use_node_mask && !node_mask[static_cast<std::size_t>(v)]) {
        std::size_t j_skip = i; while (j_skip < end && col[j_skip] == v) ++j_skip; i = j_skip; continue;
      }
      Cost min_edge_cost = std::numeric_limits<Cost>::max();
      std::vector<EdgeId>& selected_edges = sel_buf; selected_edges.clear();
      double best_rem_for_min_cost = -1.0;
      std::size_t j = i;
      int best_edge_id = -1;
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
            if (selection.tie_break == EdgeTieBreak::PreferHigherResidual) {
              if (static_cast<double>(rem) > best_rem_for_min_cost + kEpsilon) {
                best_edge_id = static_cast<int>(e);
                best_rem_for_min_cost = static_cast<double>(rem);
              } else if (std::abs(static_cast<double>(rem) - best_rem_for_min_cost) <= kEpsilon) {
                if (best_edge_id < 0 || static_cast<int>(e) < best_edge_id) {
                  best_edge_id = static_cast<int>(e);
                }
              }
            } else {
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
      if (!selected_edges.empty()) {
        Cost new_cost = static_cast<Cost>(d_u + min_edge_cost);
        auto v_idx = static_cast<std::size_t>(v);
        Cap max_edge_residual = static_cast<Cap>(0);
        for (auto edge_id : selected_edges) {
          const Cap rem = has_residual ? residual[static_cast<std::size_t>(edge_id)]
                                        : cap[static_cast<std::size_t>(edge_id)];
          if (rem > max_edge_residual) max_edge_residual = rem;
        }
        Cap path_residual = std::min(min_residual_to_node[static_cast<std::size_t>(u)], max_edge_residual);
        if (new_cost < dist[v_idx] ||
            (!multipath && new_cost == dist[v_idx] && !settled[v_idx] &&
             path_residual > min_residual_to_node[v_idx] + kEpsilon)) {
          if (dist[v_idx] == INF) ws.touch(v);
          dist[v_idx] = new_cost;
          min_residual_to_node[v_idx] = path_residual;
          pred_clear(v_idx);
          for (auto sel_e : selected_edges) pred_append(v_idx, u, sel_e);
          pq.push(new_cost, -path_residual, v);
        }
        else if (multipath && new_cost == dist[v_idx] && !settled[v_idx]) {
          for (auto sel_e : selected_edges) pred_append(v_idx, u, sel_e);
        }
      }
      i = j;
    }
    if (have_best_dst) { if (pq.empty() || pq.top_cost() > best_dst_cost) break; }
  }
  pq.end();

  PredDAG dag;
  dag.parent_offsets.assign(static_cast<std::size_t>(N+1), 0);
  if (ws.tracking() && ws.sparse_reset()) {
    for (auto v : ws.touched) {
      std::size_t c = 0;
      for (std::int32_t i = pred_head[static_cast<std::size_t>(v)]; i >= 0; i = ent_next[static_cast<std::size_t>(i)]) ++c;
      dag.parent_offsets[static_cast<std::size_t>(v+1)] = static_cast<std::int32_t>(c);
    }
  } else {
    for (std::int32_t v=0; v<N; ++v) {
      std::size_t c=0;
      for (std::int32_t i = pred_head[static_cast<std::size_t>(v)]; i >= 0; i = ent_next[static_cast<std::size_t>(i)]) ++c;
      dag.parent_offsets[static_cast<std::size_t>(v+1)] = static_cast<std::int32_t>(c);
    }
  }
  for (std::size_t k=1; k<dag.parent_offsets.size(); ++k)
    dag.parent_offsets[k] += dag.parent_offsets[k-1];
  dag.parents.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  dag.via_edges.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  auto fill = [&](std::int32_t v) {
    auto base = static_cast<std::size_t>(dag.parent_offsets[static_cast<std::size_t>(v)]);
    std::size_t k = 0;
    for (std::int32_t i = pred_head[static_cast<std::size_t>(v)]; i >= 0; i = ent_next[static_cast<std::size_t>(i)]) {
      dag.parents[base+k] = ent_parent[static_cast<std::size_t>(i)];
      dag.via_edges[base+k] = ent_edge[static_cast<std::size_t>(i)];
      ++k;
    }
  };
  if (ws.tracking() && ws.sparse_reset()) { for (auto v : ws.touched) fill(v); }
  else { for (std::int32_t v=0; v<N; ++v) fill(v); }
  auto out = finish_dist();
  ws.end();
  return {std::move(out), std::move(dag)};
}

// Graph-level dispatch metadata: maximum edge cost and the gcd of all costs.
struct CostProfile { Cost max_cost = 0; Cost stride = 1; };
inline CostProfile profile_costs(const StrictMultiDiGraph& g) {
  CostProfile p; Cost gcd = 0;
  for (auto c : g.cost_view()) {
    if (c > p.max_cost) p.max_cost = c;
    Cost a = gcd, b = c; while (b) { Cost t = a % b; a = b; b = t; } gcd = a;
  }
  p.stride = gcd > 0 ? gcd : 1;
  return p;
}

} // namespace spfx
