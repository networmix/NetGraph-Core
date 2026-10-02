/* Path enumeration from PredDAG (resolve_to_paths). */
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/constants.hpp"
#include "netgraph/core/profiling.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
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

std::pair<std::vector<Cost>, PredDAG>
path_to_pred_dag(const StrictMultiDiGraph& g,
                 std::span<const NodeId> nodes,
                 std::span<const EdgeId> edges) {
  const auto cost_view = g.cost_view();
  std::vector<Cost> dist(static_cast<std::size_t>(g.num_nodes()),
                         std::numeric_limits<Cost>::max());
  PredDAG dag;
  dag.parent_offsets.assign(static_cast<std::size_t>(g.num_nodes() + 1), 0);
  if (!nodes.empty()) {
    dist[static_cast<std::size_t>(nodes.front())] = 0;
    for (std::size_t i = 1; i < nodes.size(); ++i) {
      auto u = nodes[i - 1]; auto v = nodes[i]; auto e = edges[i - 1];
      dist[static_cast<std::size_t>(v)] =
          dist[static_cast<std::size_t>(u)] + cost_view[static_cast<std::size_t>(e)];
      dag.parent_offsets[static_cast<std::size_t>(v + 1)] = 1;
    }
    for (std::size_t v = 1; v < dag.parent_offsets.size(); ++v)
      dag.parent_offsets[v] += dag.parent_offsets[v - 1];
    dag.parents.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
    dag.via_edges.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
    for (std::size_t i = 1; i < nodes.size(); ++i) {
      auto v = nodes[i];
      auto base = static_cast<std::size_t>(dag.parent_offsets[static_cast<std::size_t>(v)]);
      dag.parents[base] = nodes[i - 1];
      dag.via_edges[base] = edges[i - 1];
    }
  }
  return {std::move(dist), std::move(dag)};
}

PredDAG make_path_dag(const StrictMultiDiGraph& g, std::span<const EdgeId> edges) {
  if (edges.empty()) {
    throw std::invalid_argument("make_path_dag: edges must be non-empty");
  }
  const auto E = static_cast<std::size_t>(g.num_edges());
  const auto esrc = g.edge_src_view();
  const auto edst = g.edge_dst_view();
  std::vector<NodeId> nodes;
  nodes.reserve(edges.size() + 1);
  std::vector<char> seen(static_cast<std::size_t>(g.num_nodes()), 0);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto e = edges[i];
    if (e < 0 || static_cast<std::size_t>(e) >= E) {
      throw std::invalid_argument("make_path_dag: edge id out of range");
    }
    const auto u = esrc[static_cast<std::size_t>(e)];
    const auto v = edst[static_cast<std::size_t>(e)];
    if (i == 0) {
      nodes.push_back(u);
      seen[static_cast<std::size_t>(u)] = 1;
    } else if (u != nodes.back()) {
      throw std::invalid_argument(
          "make_path_dag: edges are not contiguous (edge source does not match the "
          "previous edge's destination)");
    }
    if (seen[static_cast<std::size_t>(v)]) {
      throw std::invalid_argument("make_path_dag: path revisits a node (must be a simple path)");
    }
    seen[static_cast<std::size_t>(v)] = 1;
    nodes.push_back(v);
  }
  return path_to_pred_dag(g, nodes, edges).second;
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
  // on_path[v] marks nodes on the current DFS stack. A well-formed PredDAG is
  // acyclic, but zero-cost edges (or a caller-supplied DAG) can contain cycles;
  // without this guard the enumeration below would walk them forever.
  std::vector<char> on_path(dag.parent_offsets.size() > 0 ? dag.parent_offsets.size() - 1 : 0, 0);
  // start from dst
  Frame start; start.node = dst; start.idx = 0; group_parents(dag, dst, start.groups);
  stack.push_back(std::move(start));
  on_path[static_cast<std::size_t>(dst)] = 1;

  std::vector<std::pair<NodeId, std::vector<EdgeId>>> current; // reversed path accum

  while (!stack.empty()) {
    auto& top = stack.back();
    if (top.idx >= top.groups.size()) {
      // backtrack
      on_path[static_cast<std::size_t>(top.node)] = 0;
      stack.pop_back();
      if (!current.empty()) current.pop_back();
      continue;
    }
    auto [parent, edges] = top.groups[top.idx++];
    // Skip parents already on the current path (cycle in the input DAG).
    if (parent != src && on_path[static_cast<std::size_t>(parent)]) continue;
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
    on_path[static_cast<std::size_t>(parent)] = 1;
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

  Frontier queue. The reference queue is a binary heap on (cost, -bottleneck,
  node). When the graph's costs fit a bounded number of buckets
  (W = max_cost / gcd(costs) + 1 <= kMaxBucketWidth) the search uses a monotone
  bucket queue instead: Dijkstra pops costs in non-decreasing order and every
  push is at least the cost being popped, so with W cyclic buckets each bucket
  holds exactly one cost value at a time. The active bucket is itself a small
  heap on (-bottleneck, node), so the pop sequence is identical to the heap's,
  including zero-cost pushes into the active bucket, single-path re-pushes,
  stale entries and destination early exit. The two queues therefore produce
  bit-identical results; ShortestPaths.QueueDifferential_* enforce that.

  Scratch state lives in a thread-local workspace that is reset after each
  search: only the touched nodes when few were touched, the whole arrays
  otherwise. This removes the O(N) per-query initialisation that dominated small
  destination-limited searches on large graphs. A workspace holds no pointer
  into any graph, residual or mask; see ForwardWorkspace for what it retains
  between searches and when it shrinks.
*/
namespace netgraph::core {

namespace {

constexpr Cost kInf = std::numeric_limits<Cost>::max();
constexpr std::size_t kMaxBucketWidth = std::size_t{1} << 16;

using QItem = std::tuple<Cost, Cap, NodeId>;

// Reference queue: binary heap on (cost, -bottleneck, node). Keys are unique per
// push (a node is re-pushed only with a strictly smaller cost or a strictly
// larger bottleneck), so the pop order is fully determined by the key.
struct HeapQueue {
  // Functor, not a function pointer: std::push_heap / pop_heap inline it.
  struct Greater { bool operator()(const QItem& a, const QItem& b) const { return a > b; } };
  std::vector<QItem> h;
  void begin(std::size_t, Cost) { h.clear(); }
  void end() { h.clear(); }
  void drop_storage() { std::vector<QItem> fresh; h.swap(fresh); }
  bool empty() const { return h.empty(); }
  Cost top_cost() const { return std::get<0>(h.front()); }
  QItem pop() {
    std::pop_heap(h.begin(), h.end(), Greater{});
    QItem t = h.back(); h.pop_back(); return t;
  }
  void push(Cost c, Cap nr, NodeId v) {
    h.emplace_back(c, nr, v);
    std::push_heap(h.begin(), h.end(), Greater{});
  }
};

// Monotone bucket queue with ordered buckets (see the file comment).
struct BucketQueue {
  struct Item { Cap neg_res; NodeId node; };
  struct Greater {
    bool operator()(const Item& a, const Item& b) const {
      return a.neg_res != b.neg_res ? a.neg_res > b.neg_res : a.node > b.node;
    }
  };
  std::vector<std::vector<Item>> buckets;
  std::vector<std::uint64_t> occ;     // one bit per non-empty bucket; all zero between searches
  std::vector<std::uint32_t> dirty;   // occupancy words set during the current search
  std::size_t W = 0, cur_b = 0, pending = 0;
  std::size_t item_capacity = 0;      // sum of the buckets' retained capacities
  std::size_t peak_pending = 0;       // largest frontier of the current search
  Cost cur = 0, stride = 1;
  // A frontier of exactly one item is held inline rather than in a bucket, so
  // chain-like stretches pay no bucket bookkeeping; a second push promotes it.
  bool inline_item = false;
  QItem single;

  void begin(std::size_t width, Cost stride_) {
    stride = stride_;
    W = width;
    if (buckets.size() < W) buckets.resize(W);
    const std::size_t words = (W + 63) / 64;
    if (occ.size() < words) occ.resize(words, 0);   // existing words are already zero
    cur = 0; cur_b = 0; pending = 0; peak_pending = 0; inline_item = false;
  }
  // Early exit can leave items behind: clear the buckets of the words this
  // search set (O(buckets used), not O(W)). Retained item storage is dropped
  // when it exceeds a small multiple of what this search needed, so the
  // workspace does not accumulate the high-water marks of every bucket ever
  // used across differently shaped graphs.
  void end() {
    for (auto w : dirty) {
      std::uint64_t word = occ[w];
      while (word) {
        const std::size_t b = (static_cast<std::size_t>(w) << 6) + static_cast<std::size_t>(std::countr_zero(word));
        buckets[b].clear();
        word &= word - 1;
      }
      occ[w] = 0;
    }
    dirty.clear();
    pending = 0; inline_item = false;
    if (item_capacity > 4 * peak_pending + 4096) drop_storage();
  }
  void drop_storage() {
    for (auto& b : buckets) { std::vector<Item> fresh; b.swap(fresh); }
    item_capacity = 0;
  }
  bool empty() const { return pending == 0; }
  std::size_t bucket_of(Cost c) const {
    const Cost q = stride == 1 ? c : c / stride;
    return static_cast<std::size_t>(q % static_cast<Cost>(W));
  }
  void push_bucket(Cost c, Cap nr, NodeId v) {
    const std::size_t b = bucket_of(c);
    auto& vec = buckets[b];
    if (vec.empty()) {
      const std::size_t w = b >> 6;
      if (occ[w] == 0) dirty.push_back(static_cast<std::uint32_t>(w));
      occ[w] |= (std::uint64_t{1} << (b & 63));
    }
    const std::size_t cap_before = vec.capacity();
    vec.push_back({nr, v});
    if (vec.capacity() != cap_before) item_capacity += vec.capacity() - cap_before;
    if (b == cur_b) std::push_heap(vec.begin(), vec.end(), Greater{});
    ++pending;
    if (pending > peak_pending) peak_pending = pending;
  }
  void push(Cost c, Cap nr, NodeId v) {
    if (pending == 0) {
      single = {c, nr, v}; inline_item = true; pending = 1;
      if (peak_pending == 0) peak_pending = 1;
      return;
    }
    if (inline_item) {
      // Promote the inline item into the buckets. cur is the cost of the last
      // popped item, so every pending cost lies in [cur, cur + max_cost] and
      // the cyclic bucket invariant holds from cur_b = bucket_of(cur).
      inline_item = false; pending = 0;
      cur_b = bucket_of(cur);
      push_bucket(std::get<0>(single), std::get<1>(single), std::get<2>(single));
    }
    push_bucket(c, nr, v);
  }
  // Advance to the next non-empty bucket when the active one is exhausted.
  // Precondition: pending > 0 and no inline item.
  void settle_active() {
    if (!buckets[cur_b].empty()) return;
    std::size_t i = cur_b + 1; if (i == W) i = 0;
    std::size_t w = i >> 6;
    std::uint64_t word = occ[w] & (~std::uint64_t{0} << (i & 63));
    const std::size_t nwords = (W + 63) / 64;
    while (!word) { w = (w + 1 == nwords) ? 0 : w + 1; word = occ[w]; }
    const std::size_t f = (w << 6) + static_cast<std::size_t>(std::countr_zero(word));
    const std::size_t steps = f > cur_b ? f - cur_b : f + W - cur_b;
    cur += static_cast<Cost>(steps) * stride; cur_b = f;
    std::make_heap(buckets[f].begin(), buckets[f].end(), Greater{});
  }
  Cost top_cost() {
    if (inline_item) return std::get<0>(single);
    settle_active(); return cur;
  }
  QItem pop() {
    if (inline_item) {
      inline_item = false; pending = 0;
      cur = std::get<0>(single);   // cur_b keeps the last activation; promotion recomputes it
      return single;
    }
    settle_active();
    auto& vec = buckets[cur_b];
    std::pop_heap(vec.begin(), vec.end(), Greater{});
    const Item it = vec.back(); vec.pop_back(); --pending;
    if (vec.empty()) occ[cur_b >> 6] &= ~(std::uint64_t{1} << (cur_b & 63));
    return {cur, it.neg_res, it.node};
  }
};

struct BucketPlan { bool eligible = false; std::size_t width = 0; Cost stride = 1; };

BucketPlan plan_buckets(const StrictMultiDiGraph& g) {
  BucketPlan p;
  p.stride = g.cost_gcd() > 0 ? g.cost_gcd() : 1;
  const Cost w = g.max_cost() / p.stride + 1;   // max_cost < 2^62, no overflow
  p.eligible = w <= static_cast<Cost>(kMaxBucketWidth);
  p.width = p.eligible ? static_cast<std::size_t>(w) : 0;
  return p;
}

SpfQueue env_queue_override() {
  static const SpfQueue q = [] {
    const char* v = std::getenv("NGRAPH_CORE_SPF_QUEUE");
    if (!v) return SpfQueue::Auto;
    const std::string s(v);
    if (s == "heap") return SpfQueue::Heap;
    if (s == "bucket") return SpfQueue::Bucket;
    return SpfQueue::Auto;
  }();
  return q;
}

bool use_bucket_queue(const StrictMultiDiGraph& g, SpfQueue requested, BucketPlan& plan) {
  const SpfQueue q = requested == SpfQueue::Auto ? env_queue_override() : requested;
  if (q == SpfQueue::Heap) return false;
  plan = plan_buckets(g);
  return plan.eligible;
}

// Shared reset policy: sparse when few nodes were touched, wholesale otherwise.
// The same predicate decides how outputs are materialised, so the two agree.
template <class WS>
bool sparse_reset(const WS& ws) { return ws.touched.size() * 4 < ws.dist.size(); }

// Resizes a fixed-size scratch array; returns true when a much larger buffer
// was released because the graph shrank.
template <class Vec>
bool ensure_sized(Vec& v, std::size_t n, typename Vec::value_type init) {
  if (v.size() == n) return false;
  bool released = false;
  if (n * 4 < v.capacity()) { Vec fresh; v.swap(fresh); released = true; }
  v.assign(n, init);
  return released;
}

template <class Vec>
void drop(Vec& v) { Vec fresh; v.swap(fresh); }

// Forward search scratch. Invariant between searches: dist == kInf,
// min_res == 0, pred_head/pred_tail == -1, pred_count == 0, settled == 0 for
// every node; the entry arrays, touched list and queues are empty.
//
// Retention: the node arrays are sized to the current graph and released when
// a graph less than a quarter the size is searched; at that point the entry
// arrays, touched list, heap and bucket storage are released too. Otherwise
// those growable buffers keep their high-water capacity, except that the
// bucket queue drops its item storage when it exceeds a small multiple of the
// last search's peak frontier. The searches are synchronous and do not call
// back into user code, so a workspace is never re-entered on its thread.
struct ForwardWorkspace {
  std::vector<Cost> dist;
  std::vector<Cap> min_res;
  std::vector<std::int32_t> pred_head, pred_tail, pred_count;
  std::vector<char> settled;
  std::vector<NodeId> ent_parent;
  std::vector<EdgeId> ent_edge;
  std::vector<std::int32_t> ent_next;
  std::vector<NodeId> touched;
  std::vector<EdgeId> sel_buf;
  HeapQueue heap;
  BucketQueue bucket;

  void acquire(std::size_t n) {
    const bool shrank = ensure_sized(dist, n, kInf);
    ensure_sized(min_res, n, static_cast<Cap>(0));
    ensure_sized(pred_head, n, -1);
    ensure_sized(pred_tail, n, -1);
    ensure_sized(pred_count, n, 0);
    ensure_sized(settled, n, 0);
    if (shrank) {
      drop(ent_parent); drop(ent_edge); drop(ent_next); drop(touched); drop(sel_buf);
      heap.drop_storage(); bucket.drop_storage();
    }
    ent_parent.clear(); ent_edge.clear(); ent_next.clear(); touched.clear();
  }
  void release() {
    if (sparse_reset(*this)) {
      for (auto v : touched) {
        const auto i = static_cast<std::size_t>(v);
        dist[i] = kInf; min_res[i] = 0; pred_head[i] = -1; pred_tail[i] = -1; pred_count[i] = 0; settled[i] = 0;
      }
    } else {
      std::fill(dist.begin(), dist.end(), kInf);
      std::fill(min_res.begin(), min_res.end(), static_cast<Cap>(0));
      std::fill(pred_head.begin(), pred_head.end(), -1);
      std::fill(pred_tail.begin(), pred_tail.end(), -1);
      std::fill(pred_count.begin(), pred_count.end(), 0);
      std::fill(settled.begin(), settled.end(), 0);
    }
    touched.clear();
    heap.end(); bucket.end();
  }
};

// Reverse search scratch; same invariant. Successor lists per tail node u are
// a flat intrusive list of (v, e) entries (succ_head/succ_tail into ent_*).
// has_in_entry[v]: some DAG entry u -> v exists; a fanout edge may not leave
// such a node.
struct ReverseWorkspace {
  std::vector<Cost> dist;
  std::vector<Cap> min_res;
  std::vector<std::int32_t> succ_head, succ_tail;
  std::vector<char> settled;
  std::vector<char> has_in_entry;
  std::vector<NodeId> ent_node;
  std::vector<EdgeId> ent_edge;
  std::vector<std::int32_t> ent_next;
  std::vector<NodeId> touched;
  std::vector<EdgeId> sel_buf;
  HeapQueue heap;
  BucketQueue bucket;

  void acquire(std::size_t n) {
    const bool shrank = ensure_sized(dist, n, kInf);
    ensure_sized(min_res, n, static_cast<Cap>(0));
    ensure_sized(succ_head, n, -1);
    ensure_sized(succ_tail, n, -1);
    ensure_sized(settled, n, 0);
    ensure_sized(has_in_entry, n, 0);
    if (shrank) {
      drop(ent_node); drop(ent_edge); drop(ent_next); drop(touched); drop(sel_buf);
      heap.drop_storage(); bucket.drop_storage();
    }
    ent_node.clear(); ent_edge.clear(); ent_next.clear(); touched.clear();
  }
  void release() {
    if (sparse_reset(*this)) {
      for (auto v : touched) {
        const auto i = static_cast<std::size_t>(v);
        dist[i] = kInf; min_res[i] = 0; succ_head[i] = -1; succ_tail[i] = -1;
        settled[i] = 0; has_in_entry[i] = 0;
      }
    } else {
      std::fill(dist.begin(), dist.end(), kInf);
      std::fill(min_res.begin(), min_res.end(), static_cast<Cap>(0));
      std::fill(succ_head.begin(), succ_head.end(), -1);
      std::fill(succ_tail.begin(), succ_tail.end(), -1);
      std::fill(settled.begin(), settled.end(), 0);
      std::fill(has_in_entry.begin(), has_in_entry.end(), 0);
    }
    touched.clear();
    heap.end(); bucket.end();
  }
};

thread_local ForwardWorkspace tls_forward;
thread_local ReverseWorkspace tls_reverse;

// Restores the workspace invariant on every exit path, including exceptions.
// Callers must record a node in `touched` before writing its scratch entries,
// so that a throwing push_back cannot leave an unrecorded modification behind.
template <class WS>
struct WorkspaceLease {
  WS& ws;
  ~WorkspaceLease() { ws.release(); }
};

// The search loop, shared by both queues. Precondition: src is allowed and
// initialised in ws (dist 0, min_res max, touched).
template <class Queue>
void run_forward(Queue& pq, ForwardWorkspace& ws, const StrictMultiDiGraph& g, NodeId src,
                 std::optional<NodeId> dst, bool multipath, const EdgeSelection& selection,
                 std::span<const Cap> residual, std::span<const bool> node_mask,
                 std::span<const bool> edge_mask, const BucketPlan& plan) {
  const auto N = g.num_nodes();
  const auto row = g.row_offsets_view();
  const auto col = g.col_indices_view();
  const auto aei = g.adj_edge_index_view();
  const auto cost = g.cost_view();
  const auto cap  = g.capacity_view();
  // Raw pointers into the fixed-size arrays: they are not resized during the
  // search, and hoisting them keeps the compiler from reloading the vector
  // bases around every store (the workspace is a thread-local object).
  Cost* const dist = ws.dist.data();
  Cap* const min_residual_to_node = ws.min_res.data();
  std::int32_t* const pred_head = ws.pred_head.data();
  std::int32_t* const pred_tail = ws.pred_tail.data();
  std::int32_t* const pred_count = ws.pred_count.data();
  char* const settled = ws.settled.data();
  auto& ent_parent = ws.ent_parent;
  auto& ent_edge = ws.ent_edge;
  auto& ent_next = ws.ent_next;

  const bool use_node_mask = (node_mask.size() == static_cast<std::size_t>(N));
  const bool use_edge_mask = (edge_mask.size() == static_cast<std::size_t>(g.num_edges()));
  const bool has_residual = (residual.size() == static_cast<std::size_t>(g.num_edges()));
  const bool require_cap = selection.require_capacity || has_residual;

  // Predecessor storage as a flat intrusive list: pred_head/pred_tail index into
  // the ent_* arrays, whose entries are (parent, via_edge) pairs appended in
  // discovery order. pred_count tracks the live length of each list so the
  // result can be sized without walking the lists twice.
  auto pred_clear = [&](std::size_t v){ pred_head[v] = -1; pred_tail[v] = -1; pred_count[v] = 0; };
  auto pred_append = [&](std::size_t v, NodeId p, EdgeId e){
    const auto idx = static_cast<std::int32_t>(ent_parent.size());
    ent_parent.push_back(p); ent_edge.push_back(e); ent_next.push_back(-1);
    if (pred_head[v] < 0) { pred_head[v] = idx; }
    else { ent_next[static_cast<std::size_t>(pred_tail[v])] = idx; }
    pred_tail[v] = idx;
    ++pred_count[v];
  };

  // Queue items are (cost, -residual, node): cost (minimize) -> residual
  // (maximize) -> node (deterministic). This naturally distributes flows across
  // equal-cost paths based on available capacity.
  pq.begin(plan.width, plan.stride);
  pq.push(static_cast<Cost>(0), -std::numeric_limits<Cap>::max(), src);
  Cost best_dst_cost = kInf;
  bool have_best_dst = false;
  const bool early_exit = dst.has_value();
  const NodeId dst_node = dst.value_or(-1);

  // settled[v] is set once v is popped at its final distance. Equal-cost
  // predecessor updates are only accepted while v is unsettled: with positive
  // edge costs every equal-cost parent is discovered before v settles, and with
  // zero-cost edges this guard is what keeps the PredDAG acyclic (previously a
  // zero-cost pair u<->v recorded each node as the other's parent).
  std::vector<EdgeId>& selected_edges = ws.sel_buf;
  while (!pq.empty()) {
    auto [d_u, neg_res_u, u] = pq.pop();
    if (u < 0 || u >= N) continue;
    // Skip stale entries (node already processed at a lower cost).
    if (d_u > dist[static_cast<std::size_t>(u)]) continue;
    // Skip residual-stale entries in single-path mode (same cost but outdated residual).
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
    // Process edges grouped by destination node (parallel edges are consecutive in CSR).
    // A settled neighbour cannot be improved (costs are non-negative) and no
    // longer accepts equal-cost or capacity-tie updates, so its whole parallel
    // group is skipped before any edge is examined.
    while (i < end) {
      NodeId v = col[i];
      if (settled[static_cast<std::size_t>(v)] || (use_node_mask && !node_mask[static_cast<std::size_t>(v)])) {
        std::size_t j_skip = i; while (j_skip < end && col[j_skip] == v) ++j_skip; i = j_skip; continue;
      }

      // Select best edge(s) from u to v according to policy.
      Cost min_edge_cost = kInf;
      selected_edges.clear();
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
      if (!selected_edges.empty()) {
        Cost new_cost = static_cast<Cost>(d_u + min_edge_cost);
        auto v_idx = static_cast<std::size_t>(v);
        // Bottleneck capacity along the path to v through u (node-level tie-breaking).
        Cap max_edge_residual = static_cast<Cap>(0);
        for (auto edge_id : selected_edges) {
          const Cap rem = has_residual ? residual[static_cast<std::size_t>(edge_id)]
                                        : cap[static_cast<std::size_t>(edge_id)];
          if (rem > max_edge_residual) max_edge_residual = rem;
        }
        Cap path_residual = std::min(min_residual_to_node[static_cast<std::size_t>(u)], max_edge_residual);
        // Relaxation: shorter path, or equal-cost path with better capacity (single-path mode).
        if (new_cost < dist[v_idx] ||
            (!multipath && new_cost == dist[v_idx] && !settled[v_idx] &&
             path_residual > min_residual_to_node[v_idx] + kEpsilon)) {
          if (dist[v_idx] == kInf) ws.touched.push_back(v);
          dist[v_idx] = new_cost;
          min_residual_to_node[v_idx] = path_residual;
          pred_clear(v_idx);
          for (auto sel_e : selected_edges) pred_append(v_idx, u, sel_e);
          pq.push(new_cost, -path_residual, v);
        }
        // Multipath: equal-cost alternative (only while v is unsettled; see above).
        // min_residual_to_node is not updated here: all equal-cost paths are kept,
        // none is chosen by residual.
        else if (multipath && new_cost == dist[v_idx] && !settled[v_idx]) {
          for (auto sel_e : selected_edges) pred_append(v_idx, u, sel_e);
        }
      }
      i = j;
    }
    if (have_best_dst) { if (pq.empty() || pq.top_cost() > best_dst_cost) break; }
  }
}

// Materialise (distances, PredDAG) from the forward workspace. Per-node
// predecessor lists are independent, so visiting only touched nodes yields the
// same arrays as visiting every node.
std::pair<std::vector<Cost>, PredDAG> finish_forward(ForwardWorkspace& ws, std::int32_t N) {
  const bool sparse = sparse_reset(ws);
  std::vector<Cost> dist;
  if (sparse) {
    dist.assign(static_cast<std::size_t>(N), kInf);
    for (auto v : ws.touched) dist[static_cast<std::size_t>(v)] = ws.dist[static_cast<std::size_t>(v)];
  } else {
    dist = ws.dist;
  }
  PredDAG dag;
  dag.parent_offsets.assign(static_cast<std::size_t>(N+1), 0);
  if (sparse) {
    for (auto v : ws.touched) dag.parent_offsets[static_cast<std::size_t>(v) + 1] = ws.pred_count[static_cast<std::size_t>(v)];
  } else {
    for (std::int32_t v = 0; v < N; ++v) dag.parent_offsets[static_cast<std::size_t>(v) + 1] = ws.pred_count[static_cast<std::size_t>(v)];
  }
  for (std::size_t k = 1; k < dag.parent_offsets.size(); ++k)
    dag.parent_offsets[k] += dag.parent_offsets[k-1];
  dag.parents.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  dag.via_edges.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  auto fill = [&](std::int32_t v) {
    auto base = static_cast<std::size_t>(dag.parent_offsets[static_cast<std::size_t>(v)]);
    std::size_t k = 0;
    for (std::int32_t i = ws.pred_head[static_cast<std::size_t>(v)]; i >= 0; i = ws.ent_next[static_cast<std::size_t>(i)]) {
      dag.parents[base+k] = ws.ent_parent[static_cast<std::size_t>(i)];
      dag.via_edges[base+k] = ws.ent_edge[static_cast<std::size_t>(i)];
      ++k;
    }
  };
  if (sparse) { for (auto v : ws.touched) fill(v); }
  else { for (std::int32_t v = 0; v < N; ++v) fill(v); }
  return {std::move(dist), std::move(dag)};
}

std::pair<std::vector<Cost>, PredDAG>
shortest_paths_core(const StrictMultiDiGraph& g, NodeId src,
                    std::optional<NodeId> dst,
                    bool multipath,
                    const EdgeSelection& selection,
                    std::span<const Cap> residual,
                    std::span<const bool> node_mask,
                    std::span<const bool> edge_mask,
                    SpfQueue queue) {
  NGRAPH_PROFILE_SCOPE("shortest_paths_core");
  const auto N = g.num_nodes();
  const bool use_node_mask = (node_mask.size() == static_cast<std::size_t>(N));
  const bool src_allowed = (src >= 0 && src < N && (!use_node_mask || node_mask[static_cast<std::size_t>(src)]));
  if (!src_allowed) {
    // Source is out of range or masked out: no traversal, return empty DAG.
    PredDAG dag;
    dag.parent_offsets.assign(static_cast<std::size_t>(N + 1), 0);
    return {std::vector<Cost>(static_cast<std::size_t>(N), kInf), std::move(dag)};
  }
  ForwardWorkspace& ws = tls_forward;
  ws.acquire(static_cast<std::size_t>(N));
  WorkspaceLease<ForwardWorkspace> lease{ws};
  ws.touched.push_back(src);   // before the writes: see WorkspaceLease
  ws.dist[static_cast<std::size_t>(src)] = static_cast<Cost>(0);
  ws.min_res[static_cast<std::size_t>(src)] = std::numeric_limits<Cap>::max();
  BucketPlan plan;
  if (use_bucket_queue(g, queue, plan)) {
    run_forward(ws.bucket, ws, g, src, dst, multipath, selection, residual, node_mask, edge_mask, plan);
  } else {
    run_forward(ws.heap, ws, g, src, dst, multipath, selection, residual, node_mask, edge_mask, plan);
  }
  return finish_forward(ws, N);
}

} // namespace

std::pair<std::vector<Cost>, PredDAG>
shortest_paths(const StrictMultiDiGraph& g, NodeId src,
               std::optional<NodeId> dst,
               bool multipath,
               const EdgeSelection& selection,
               std::span<const Cap> residual,
               std::span<const bool> node_mask,
               std::span<const bool> edge_mask,
               SpfQueue queue) {
  if (!node_mask.empty() && node_mask.size() != static_cast<std::size_t>(g.num_nodes())) {
    throw std::invalid_argument("shortest_paths: node_mask length mismatch");
  }
  if (!edge_mask.empty() && edge_mask.size() != static_cast<std::size_t>(g.num_edges())) {
    throw std::invalid_argument("shortest_paths: edge_mask length mismatch");
  }
  if (!residual.empty() && residual.size() != static_cast<std::size_t>(g.num_edges())) {
    throw std::invalid_argument("shortest_paths: residual length mismatch");
  }
  return shortest_paths_core(g, src, dst, multipath, selection, residual, node_mask, edge_mask, queue);
}

} // namespace netgraph::core

/*
  Reverse Dijkstra: shortest paths from every node *to* one destination.

  Mirrors the forward search over the in-adjacency (in_row_offsets /
  in_col_indices / in_adj_edge_index): the queue settles nodes by their distance
  to dst, edge selection runs per (u -> v) parallel group exactly as in the
  forward variant, and node-level tie-breaking in single-path mode prefers the
  higher bottleneck capacity toward dst. Successor entries (v, e) are kept per
  tail node u while u is unsettled, then bucketed by v into the forward PredDAG
  layout that placement consumes. The same two queues apply.
*/
namespace netgraph::core {

namespace {

template <class Queue>
void run_reverse(Queue& pq, ReverseWorkspace& ws, const StrictMultiDiGraph& g, NodeId target,
                 bool multipath, const EdgeSelection& selection,
                 std::span<const Cap> residual, std::span<const bool> node_mask,
                 std::span<const bool> edge_mask, const BucketPlan& plan) {
  const auto N = g.num_nodes();
  const auto E = g.num_edges();
  const auto irow = g.in_row_offsets_view();
  const auto icol = g.in_col_indices_view();
  const auto iaei = g.in_adj_edge_index_view();
  const auto cost = g.cost_view();
  const auto cap  = g.capacity_view();
  Cost* const dist = ws.dist.data();
  Cap* const min_residual_to_dst = ws.min_res.data();
  std::int32_t* const succ_head = ws.succ_head.data();
  std::int32_t* const succ_tail = ws.succ_tail.data();
  char* const settled = ws.settled.data();
  char* const has_in_entry = ws.has_in_entry.data();
  auto& ent_node = ws.ent_node;
  auto& ent_edge = ws.ent_edge;
  auto& ent_next = ws.ent_next;

  const bool use_node_mask = (node_mask.size() == static_cast<std::size_t>(N));
  const bool use_edge_mask = (edge_mask.size() == static_cast<std::size_t>(E));
  const bool has_residual = (residual.size() == static_cast<std::size_t>(E));
  const bool require_cap = selection.require_capacity || has_residual;

  auto succ_clear = [&](std::size_t u){ succ_head[u] = -1; succ_tail[u] = -1; };
  auto succ_append = [&](std::size_t u, NodeId v, EdgeId e){
    const auto idx = static_cast<std::int32_t>(ent_node.size());
    ent_node.push_back(v); ent_edge.push_back(e); ent_next.push_back(-1);
    if (succ_head[u] < 0) { succ_head[u] = idx; }
    else { ent_next[static_cast<std::size_t>(succ_tail[u])] = idx; }
    succ_tail[u] = idx;
  };

  pq.begin(plan.width, plan.stride);
  pq.push(static_cast<Cost>(0), -std::numeric_limits<Cap>::max(), target);
  std::vector<EdgeId>& selected_edges = ws.sel_buf;

  while (!pq.empty()) {
    auto [d_v, neg_res_v, v] = pq.pop();
    if (v < 0 || v >= N) continue;
    if (d_v > dist[static_cast<std::size_t>(v)]) continue;
    if (!multipath && d_v == dist[static_cast<std::size_t>(v)] &&
        -neg_res_v < min_residual_to_dst[static_cast<std::size_t>(v)] - kEpsilon) continue;
    settled[static_cast<std::size_t>(v)] = 1;

    // In-edges of v, clustered by tail node u (edges are sorted by (src, dst)).
    auto start = static_cast<std::size_t>(irow[static_cast<std::size_t>(v)]);
    auto end   = static_cast<std::size_t>(irow[static_cast<std::size_t>(v)+1]);
    std::size_t i = start;
    while (i < end) {
      NodeId u = icol[i];
      if (settled[static_cast<std::size_t>(u)] || (use_node_mask && !node_mask[static_cast<std::size_t>(u)])) {
        std::size_t j_skip = i; while (j_skip < end && icol[j_skip] == u) ++j_skip; i = j_skip; continue;
      }
      Cost min_edge_cost = kInf;
      selected_edges.clear();
      double best_rem_for_min_cost = -1.0;
      std::size_t j = i;
      int best_edge_id = -1;
      for (; j < end && icol[j] == u; ++j) {
        auto e = static_cast<std::size_t>(iaei[j]);
        if (use_edge_mask && !edge_mask[e]) continue;
        const Cap rem = has_residual ? residual[e] : cap[e];
        if (require_cap && rem < kMinCap) continue;
        const Cost ecost = static_cast<Cost>(cost[e]);
        if (ecost < min_edge_cost) {
          min_edge_cost = ecost;
          selected_edges.clear();
          if (selection.multi_edge) {
            selected_edges.push_back(static_cast<EdgeId>(iaei[j]));
          } else {
            best_edge_id = static_cast<int>(e);
            best_rem_for_min_cost = static_cast<double>(rem);
          }
        } else if (ecost == min_edge_cost) {
          if (selection.multi_edge) {
            selected_edges.push_back(static_cast<EdgeId>(iaei[j]));
          } else if (selection.tie_break == EdgeTieBreak::PreferHigherResidual) {
            if (static_cast<double>(rem) > best_rem_for_min_cost + kEpsilon) {
              best_edge_id = static_cast<int>(e);
              best_rem_for_min_cost = static_cast<double>(rem);
            } else if (std::abs(static_cast<double>(rem) - best_rem_for_min_cost) <= kEpsilon) {
              if (best_edge_id < 0 || static_cast<int>(e) < best_edge_id) best_edge_id = static_cast<int>(e);
            }
          } else {
            if (best_edge_id < 0 || static_cast<int>(e) < best_edge_id) best_edge_id = static_cast<int>(e);
          }
        }
      }
      if (!selection.multi_edge && best_edge_id >= 0) {
        selected_edges.clear();
        selected_edges.push_back(static_cast<EdgeId>(best_edge_id));
      }
      if (!selected_edges.empty()) {
        const Cost new_cost = static_cast<Cost>(d_v + min_edge_cost);
        const auto u_idx = static_cast<std::size_t>(u);
        Cap max_edge_residual = static_cast<Cap>(0);
        for (auto edge_id : selected_edges) {
          const Cap rem = has_residual ? residual[static_cast<std::size_t>(edge_id)]
                                       : cap[static_cast<std::size_t>(edge_id)];
          if (rem > max_edge_residual) max_edge_residual = rem;
        }
        const Cap path_residual = std::min(min_residual_to_dst[static_cast<std::size_t>(v)], max_edge_residual);
        if (new_cost < dist[u_idx] ||
            (!multipath && new_cost == dist[u_idx] && !settled[u_idx] &&
             path_residual > min_residual_to_dst[u_idx] + kEpsilon)) {
          if (dist[u_idx] == kInf) ws.touched.push_back(u);
          dist[u_idx] = new_cost;
          min_residual_to_dst[u_idx] = path_residual;
          succ_clear(u_idx);
          for (auto sel_e : selected_edges) succ_append(u_idx, v, sel_e);
          has_in_entry[static_cast<std::size_t>(v)] = 1;
          pq.push(new_cost, -path_residual, u);
        } else if (multipath && new_cost == dist[u_idx] && !settled[u_idx]) {
          for (auto sel_e : selected_edges) succ_append(u_idx, v, sel_e);
          has_in_entry[static_cast<std::size_t>(v)] = 1;
        }
      }
      i = j;
    }
  }
}

std::pair<std::vector<Cost>, PredDAG>
shortest_paths_to_core(const StrictMultiDiGraph& g, NodeId target,
                       bool multipath,
                       const EdgeSelection& selection,
                       std::span<const Cap> residual,
                       std::span<const bool> node_mask,
                       std::span<const bool> edge_mask,
                       std::span<const EdgeId> fanout_edges,
                       SpfQueue queue) {
  NGRAPH_PROFILE_SCOPE("shortest_paths_to_core");
  const auto N = g.num_nodes();
  const auto E = g.num_edges();
  const auto esrc = g.edge_src_view();
  const auto edst = g.edge_dst_view();
  const auto cost = g.cost_view();
  const auto cap  = g.capacity_view();
  const bool use_node_mask = (node_mask.size() == static_cast<std::size_t>(N));
  const bool use_edge_mask = (edge_mask.size() == static_cast<std::size_t>(E));
  const bool has_residual = (residual.size() == static_cast<std::size_t>(E));
  const bool require_cap = selection.require_capacity || has_residual;
  const bool target_allowed = (target >= 0 && target < N &&
                               (!use_node_mask || node_mask[static_cast<std::size_t>(target)]));

  ReverseWorkspace& ws = tls_reverse;
  ws.acquire(static_cast<std::size_t>(N));
  WorkspaceLease<ReverseWorkspace> lease{ws};

  if (target_allowed) {
    ws.touched.push_back(target);   // before the writes: see WorkspaceLease
    ws.dist[static_cast<std::size_t>(target)] = static_cast<Cost>(0);
    ws.min_res[static_cast<std::size_t>(target)] = std::numeric_limits<Cap>::max();
    BucketPlan plan;
    if (use_bucket_queue(g, queue, plan)) {
      run_reverse(ws.bucket, ws, g, target, multipath, selection, residual, node_mask, edge_mask, plan);
    } else {
      run_reverse(ws.heap, ws, g, target, multipath, selection, residual, node_mask, edge_mask, plan);
    }
  }

  auto succ_append = [&](std::size_t u, NodeId v, EdgeId e){
    const auto idx = static_cast<std::int32_t>(ws.ent_node.size());
    ws.ent_node.push_back(v); ws.ent_edge.push_back(e); ws.ent_next.push_back(-1);
    if (ws.succ_head[u] < 0) { ws.succ_head[u] = idx; }
    else { ws.ent_next[static_cast<std::size_t>(ws.succ_tail[u])] = idx; }
    ws.succ_tail[u] = idx;
  };

  // Forced fan-out entries (see the header).
  for (auto e_raw : fanout_edges) {
    if (e_raw < 0 || e_raw >= E) {
      throw std::invalid_argument("shortest_paths_to: fanout edge id out of range");
    }
    const auto e = static_cast<std::size_t>(e_raw);
    const NodeId u = esrc[e];
    const NodeId v = edst[e];
    if (ws.has_in_entry[static_cast<std::size_t>(u)]) {
      throw std::invalid_argument(
          "shortest_paths_to: a fanout edge must leave a node with no incoming DAG entry "
          "(otherwise the DAG could contain a cycle)");
    }
    if (use_edge_mask && !edge_mask[e]) continue;
    if (use_node_mask && (!node_mask[static_cast<std::size_t>(u)] || !node_mask[static_cast<std::size_t>(v)])) continue;
    const Cap rem = has_residual ? residual[e] : cap[e];
    if (require_cap && rem < kMinCap) continue;
    if (ws.dist[static_cast<std::size_t>(v)] == kInf) continue;
    bool present = false;
    for (std::int32_t i = ws.succ_head[static_cast<std::size_t>(u)]; i >= 0; i = ws.ent_next[static_cast<std::size_t>(i)]) {
      if (ws.ent_edge[static_cast<std::size_t>(i)] == e_raw) { present = true; break; }
    }
    if (present) continue;
    if (ws.dist[static_cast<std::size_t>(u)] == kInf) ws.touched.push_back(u);
    succ_append(static_cast<std::size_t>(u), v, e_raw);
    if (ws.dist[static_cast<std::size_t>(u)] == kInf) {
      ws.dist[static_cast<std::size_t>(u)] = static_cast<Cost>(cost[e]) + ws.dist[static_cast<std::size_t>(v)];
    }
  }

  // Bucket successor entries (u -> v via e) by v into the forward PredDAG
  // layout. Entries within a v range are ordered by ascending tail node, so a
  // sparse pass must visit the touched tails in ascending order.
  const bool sparse = sparse_reset(ws);
  std::vector<Cost> dist;
  if (sparse) {
    std::sort(ws.touched.begin(), ws.touched.end());
    dist.assign(static_cast<std::size_t>(N), kInf);
    for (auto v : ws.touched) dist[static_cast<std::size_t>(v)] = ws.dist[static_cast<std::size_t>(v)];
  } else {
    dist = ws.dist;
  }
  PredDAG dag;
  dag.parent_offsets.assign(static_cast<std::size_t>(N+1), 0);
  auto count = [&](std::int32_t u) {
    for (std::int32_t i = ws.succ_head[static_cast<std::size_t>(u)]; i >= 0; i = ws.ent_next[static_cast<std::size_t>(i)]) {
      dag.parent_offsets[static_cast<std::size_t>(ws.ent_node[static_cast<std::size_t>(i)]) + 1] += 1;
    }
  };
  if (sparse) { for (auto u : ws.touched) count(u); }
  else { for (std::int32_t u = 0; u < N; ++u) count(u); }
  for (std::size_t k = 1; k < dag.parent_offsets.size(); ++k) dag.parent_offsets[k] += dag.parent_offsets[k-1];
  dag.parents.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  dag.via_edges.resize(static_cast<std::size_t>(dag.parent_offsets.back()));
  std::vector<std::int32_t> cursor(dag.parent_offsets.begin(), dag.parent_offsets.end() - 1);
  auto fill = [&](std::int32_t u) {
    for (std::int32_t i = ws.succ_head[static_cast<std::size_t>(u)]; i >= 0; i = ws.ent_next[static_cast<std::size_t>(i)]) {
      const auto v = static_cast<std::size_t>(ws.ent_node[static_cast<std::size_t>(i)]);
      const auto pos = static_cast<std::size_t>(cursor[v]++);
      dag.parents[pos] = u;
      dag.via_edges[pos] = ws.ent_edge[static_cast<std::size_t>(i)];
    }
  };
  if (sparse) { for (auto u : ws.touched) fill(u); }
  else { for (std::int32_t u = 0; u < N; ++u) fill(u); }
  return {std::move(dist), std::move(dag)};
}

} // namespace

std::pair<std::vector<Cost>, PredDAG>
shortest_paths_to(const StrictMultiDiGraph& g, NodeId dst,
                  bool multipath,
                  const EdgeSelection& selection,
                  std::span<const Cap> residual,
                  std::span<const bool> node_mask,
                  std::span<const bool> edge_mask,
                  std::span<const EdgeId> fanout_edges,
                  SpfQueue queue) {
  if (!node_mask.empty() && node_mask.size() != static_cast<std::size_t>(g.num_nodes())) {
    throw std::invalid_argument("shortest_paths_to: node_mask length mismatch");
  }
  if (!edge_mask.empty() && edge_mask.size() != static_cast<std::size_t>(g.num_edges())) {
    throw std::invalid_argument("shortest_paths_to: edge_mask length mismatch");
  }
  if (!residual.empty() && residual.size() != static_cast<std::size_t>(g.num_edges())) {
    throw std::invalid_argument("shortest_paths_to: residual length mismatch");
  }
  return shortest_paths_to_core(g, dst, multipath, selection, residual, node_mask, edge_mask, fanout_edges, queue);
}

} // namespace netgraph::core
