/*
  The one traversal core of every search (detail/traversal.hpp), the public
  spf, and the constructors of the search options.

  Distances are int64 path costs: every cost is positive and the sum of all costs is below 2^62,
  so no simple path cost overflows and kUnreachable (INT64_MAX) never collides.

  Both queues pop in (distance, node) order: the heap by its key, the cyclic buckets because
  each holds one distance at a time and the active bucket is sorted by node. After popping
  distance d, pending distances lie in [d, d + max_cost], all multiples of cost_gcd. A node is
  pushed only on a strict improvement, so it appears at most once per bucket. Equal-distance pop
  order cannot change the result: positive costs make every tight predecessor of a node
  strictly nearer, so it is finalised first.
*/
#include "netgraph/core/search.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <functional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

#include "netgraph/core/detail/edge_amounts.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/graph_access.hpp"
#include "netgraph/core/detail/hooks_flag.hpp"
#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/traversal.hpp"

namespace netgraph::core::detail {

int compare(const Priority& a, const Priority& b) noexcept {
  if (a.infinite || b.infinite) return a.infinite == b.infinite ? 0 : (a.infinite ? 1 : -1);
  if (a.den == 1 && b.den == 1) return a.num == b.num ? 0 : (a.num < b.num ? -1 : 1);
  const U128 left = mul_u64(static_cast<std::uint64_t>(a.num), static_cast<std::uint64_t>(b.den));
  const U128 right = mul_u64(static_cast<std::uint64_t>(b.num), static_cast<std::uint64_t>(a.den));
  if (left == right) return 0;
  return left < right ? -1 : 1;
}

Priority priority_value(EdgePriority priority, Amount residual, Amount capacity) noexcept {
  if (priority == EdgePriority::Capacity) return {capacity, 1, false};
  if (priority == EdgePriority::Residual) return {residual, 1, false};
  return {residual, capacity, false};
}

Priority PriorityValues::of(EdgeId e) const noexcept {
  const auto ei = static_cast<std::size_t>(e);
  const Amount cap = g->capacity_units()[ei];
  if (residual.empty()) return {cap, 1, false};
  return priority_value(priority, residual[ei], cap);
}

PriorityFn PriorityFn::of(const PriorityValues& values) noexcept {
  return {&values, [](const void* ctx, EdgeId e) noexcept { return static_cast<const PriorityValues*>(ctx)->of(e); }};
}

std::int64_t bucket_width(const Graph& g) noexcept {
  const Cost stride = std::max<Cost>(1, g.cost_gcd());
  return g.max_cost() / stride + 1;  // max_cost < 2^62: no overflow
}

bool uses_buckets(const Graph& g, SpfQueue queue) noexcept {
  return queue == SpfQueue::Auto && bucket_width(g) <= kMaxBucketWidth;
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
QueueCounts& thread_queue_counts() noexcept {
  thread_local QueueCounts counts;
  return counts;
}
#endif

namespace {

constexpr Cost kInf = kUnreachable;

// Node marks of the workspace (one byte per node, zero between leases).
enum : std::uint8_t {
  kTarget = 1,      // a live target of the current search, or a destination of a support route
  kExcluded = 2,    // a KSP spur's excluded prefix node
  kKeep = 4,        // a node of the pruned DAG (the walk back from the targets) or of a support route
  kBannedHead = 8,  // a KSP Collapse ban: the spur node's transition to this node
};

// Releases a buffer whose retained capacity exceeds a small multiple of what the last use
// needed, so a thread's workspace does not keep the high-water mark of every graph it saw.
template <class Vec>
void trim(Vec& v, std::size_t used) noexcept {
  if (v.capacity() > 4 * used + 65536) {
    Vec fresh;
    v.swap(fresh);
  }
}

// The reference queue: a binary heap on (distance, node).
class HeapQueue {
public:
  void begin(std::span<const NodeId> roots) {
    h_.clear();
    peak_ = 0;
    for (NodeId r : roots) push(0, r);
  }
  void push(Cost c, NodeId v) {
    h_.emplace_back(c, v);  // strong guarantee: a failed growth leaves the heap as it was
    std::push_heap(h_.begin(), h_.end(), std::greater<> {});
    if (h_.size() > peak_) peak_ = h_.size();
  }
  bool pop(Cost& c, NodeId& v) noexcept {
    if (h_.empty()) return false;
    std::pop_heap(h_.begin(), h_.end(), std::greater<> {});
    c = h_.back().first;
    v = h_.back().second;
    h_.pop_back();
    return true;
  }
  void end() noexcept {
    h_.clear();
    trim(h_, peak_);
  }

private:
  std::vector<std::pair<Cost, NodeId>> h_ {};
  std::size_t peak_ {0};
};

// Monotone cyclic buckets (see the file comment): W' = bit_ceil(W) >= 64 buckets, each an
// intrusive list into one item pool, an occupancy bit per bucket and the words set during the
// search; the active bucket is copied out and sorted by node. A frontier of exactly one item is
// held inline (chain-like stretches pay no bucket bookkeeping); a second push promotes it.
// Between searches every head is -1 and every occupancy bit 0.
class BucketQueue {
public:
  void begin(std::int64_t width, Cost stride, std::size_t num_nodes, std::span<const NodeId> roots) {
    const std::size_t w2 = std::max<std::size_t>(64, std::bit_ceil(static_cast<std::size_t>(width)));
    const std::size_t words = w2 / 64;
    // Each growth keeps the invariant of the existing part (heads -1, bits 0), and none can leave
    // the queue unusable: every later use reads only up to the sizes reached.
    if (head_.size() < w2) head_.resize(w2, -1);
    if (occ_.size() < words) occ_.resize(words, 0);
    if (dirty_flag_.size() < words) dirty_flag_.resize(words, 0);
    if (dirty_.capacity() < words) dirty_.reserve(words);       // one entry per word: never grows
    if (active_.capacity() < num_nodes) active_.reserve(num_nodes);  // a bucket holds a node once
    mask_ = w2 - 1;
    words_ = words;
    stride_ = stride;
    shift_ = std::has_single_bit(static_cast<std::uint64_t>(stride)) ? std::countr_zero(static_cast<std::uint64_t>(stride)) : -1;
    cur_ = 0;
    cur_b_ = 0;
    pending_ = 0;
    held_ = false;
    items_.clear();
    peak_ = 0;
    active_.assign(roots.begin(), roots.end());  // distance 0, sorted; within the reserved capacity
    pos_ = 0;
  }
  void push(Cost c, NodeId v) {
    if (!held_ && pending_ == 0 && pos_ == active_.size()) {
      held_ = true;
      held_cost_ = c;
      held_node_ = v;
      return;
    }
    if (held_) {
      put(held_cost_, held_node_);
      held_ = false;
    }
    put(c, v);
  }
  bool pop(Cost& c, NodeId& v) noexcept {
    if (pos_ < active_.size()) {
      c = cur_;
      v = active_[pos_++];
      return true;
    }
    if (held_) {
      held_ = false;
      cur_ = held_cost_;
      cur_b_ = index(cur_);
      c = cur_;
      v = held_node_;
      return true;
    }
    if (pending_ == 0) return false;
    activate();
    c = cur_;
    v = active_[pos_++];
    return true;
  }
  void end() noexcept {
    for (const std::uint32_t w : dirty_) {
      std::uint64_t word = occ_[w];
      while (word != 0) {
        head_[(static_cast<std::size_t>(w) << 6) + static_cast<std::size_t>(std::countr_zero(word))] = -1;
        word &= word - 1;
      }
      occ_[w] = 0;
      dirty_flag_[w] = 0;
    }
    dirty_.clear();
    items_.clear();
    active_.clear();
    pos_ = 0;
    pending_ = 0;
    held_ = false;
    trim(items_, peak_);
  }

private:
  struct Item {
    NodeId node;
    std::int32_t next;
  };
  [[nodiscard]] std::size_t index(Cost c) const noexcept {
    const Cost q = shift_ >= 0 ? (c >> shift_) : c / stride_;
    return static_cast<std::size_t>(q) & mask_;
  }
  void put(Cost c, NodeId v) {
    const std::size_t b = index(c);
    items_.push_back({v, head_[b]});  // the only allocation; a failure changes nothing
    head_[b] = static_cast<std::int32_t>(items_.size() - 1);
    const std::size_t w = b >> 6;
    if (dirty_flag_[w] == 0) {
      dirty_flag_[w] = 1;
      dirty_.push_back(static_cast<std::uint32_t>(w));  // reserved for every word
    }
    occ_[w] |= std::uint64_t {1} << (b & 63);
    ++pending_;
    if (items_.size() > peak_) peak_ = items_.size();
  }
  // The next non-empty bucket after the active one, cyclically (pending_ > 0): its items become
  // the active bucket, sorted by node.
  void activate() noexcept {
    std::size_t i = (cur_b_ + 1) & mask_;
    std::size_t w = i >> 6;
    std::uint64_t word = occ_[w] & (~std::uint64_t {0} << (i & 63));
    while (word == 0) {
      w = w + 1 == words_ ? 0 : w + 1;
      word = occ_[w];
    }
    const std::size_t f = (w << 6) + static_cast<std::size_t>(std::countr_zero(word));
    cur_ += static_cast<Cost>((f - cur_b_) & mask_) * stride_;
    cur_b_ = f;
    active_.clear();
    pos_ = 0;
    for (std::int32_t k = head_[f]; k >= 0; k = items_[static_cast<std::size_t>(k)].next) {
      active_.push_back(items_[static_cast<std::size_t>(k)].node);  // within the reserved capacity
    }
    head_[f] = -1;
    occ_[w] &= ~(std::uint64_t {1} << (f & 63));
    pending_ -= active_.size();
    if (active_.size() > 1) std::sort(active_.begin(), active_.end());
  }

  std::vector<std::int32_t> head_ {};
  std::vector<std::uint64_t> occ_ {};
  std::vector<std::uint8_t> dirty_flag_ {};
  std::vector<std::uint32_t> dirty_ {};
  std::vector<Item> items_ {};
  std::vector<NodeId> active_ {};
  std::size_t pos_ {0};
  std::size_t mask_ {0};
  std::size_t words_ {0};
  Cost stride_ {1};
  int shift_ {0};
  Cost cur_ {0};
  std::size_t cur_b_ {0};
  std::size_t pending_ {0};
  std::size_t peak_ {0};
  bool held_ {false};
  Cost held_cost_ {0};
  NodeId held_node_ {0};
};

struct PoolEntry {
  EdgeId edge;
  std::int32_t next;
};

// The per-thread search workspace. Between leases: dist == kInf, head == -1, loc == -1 and
// flags == 0 for every node; touched, marked, the pool and both queues empty. touched and marked
// have the capacity of the node arrays, so recording a node never allocates; scratch buffers are
// rewritten before every use.
struct Workspace {
  std::vector<Cost> dist {};
  std::vector<std::int32_t> head {};
  std::vector<std::int32_t> loc {};
  std::vector<std::uint8_t> flags {};
  std::vector<NodeId> touched {};  // the nodes whose dist (and head) a lease wrote
  std::vector<NodeId> marked {};   // the nodes whose flags (and loc) a lease wrote
  std::vector<PoolEntry> pool {};  // the tight-edge lists
  std::size_t pool_peak {0};
  HeapQueue heap {};
  BucketQueue buckets {};
  // scratch
  std::vector<EdgeId> list {};
  std::vector<EdgeId> pruned {};
  std::vector<NodeId> kept {};
  std::vector<NodeId> stack {};
  std::vector<std::int32_t> off_out {};
  std::vector<std::int32_t> off_in {};
  std::vector<std::int32_t> idx_out {};
  std::vector<std::int32_t> idx_in {};
  std::vector<std::int32_t> hops {};
  std::vector<std::uint8_t> keep_edge {};
  std::vector<std::uint8_t> has_best {};
  std::vector<Priority> best {};
  bool leased {false};
};

thread_local Workspace tls_workspace;

// The lease of the thread's workspace for one search (and its extraction): sizes the node arrays
// for the graph, refuses an overlapping lease, and restores the workspace's invariant on every
// exit, exceptions included. A failed acquisition changes nothing (the arrays are replaced
// together or not at all). The workspace holds no pointer into a graph, a residual or an
// Eligibility between leases.
class Lease {
public:
  explicit Lease(std::size_t n) : ws_(tls_workspace) {
    if (ws_.leased) throw std::logic_error("search: the thread's workspace is already leased (overlapping searches)");
    if (ws_.dist.size() < n || ws_.dist.size() > 4 * n + 4096) {
      std::vector<Cost> dist(n, kInf);
      std::vector<std::int32_t> head(n, -1);
      std::vector<std::int32_t> loc(n, -1);
      std::vector<std::uint8_t> flags(n, 0);
      std::vector<NodeId> touched;
      std::vector<NodeId> marked;
      touched.reserve(n);
      marked.reserve(n);
      ws_.dist.swap(dist);
      ws_.head.swap(head);
      ws_.loc.swap(loc);
      ws_.flags.swap(flags);
      ws_.touched.swap(touched);
      ws_.marked.swap(marked);
    }
    ws_.leased = true;
  }
  Lease(const Lease&) = delete;
  Lease& operator=(const Lease&) = delete;
  ~Lease() { release(); }
  [[nodiscard]] Workspace& ws() noexcept { return ws_; }

private:
  void release() noexcept {
    auto& w = ws_;
    const std::size_t n = w.dist.size();
    if (w.touched.size() * 4 >= n) {
      std::fill(w.dist.begin(), w.dist.end(), kInf);
      std::fill(w.head.begin(), w.head.end(), -1);
    } else {
      for (NodeId v : w.touched) {
        w.dist[static_cast<std::size_t>(v)] = kInf;
        w.head[static_cast<std::size_t>(v)] = -1;
      }
    }
    if (w.marked.size() * 4 >= n) {
      std::fill(w.flags.begin(), w.flags.end(), std::uint8_t {0});
      std::fill(w.loc.begin(), w.loc.end(), -1);
    } else {
      for (NodeId v : w.marked) {
        w.flags[static_cast<std::size_t>(v)] = 0;
        w.loc[static_cast<std::size_t>(v)] = -1;
      }
    }
    w.touched.clear();
    w.marked.clear();
    w.pool_peak = std::max(w.pool_peak, w.pool.size());
    w.pool.clear();
    trim(w.pool, w.pool_peak);
    w.pool_peak = 0;
    w.heap.end();
    w.buckets.end();
    w.leased = false;
  }
  Workspace& ws_;
};

void mark(Workspace& ws, NodeId v, std::uint8_t bit) noexcept {
  auto& f = ws.flags[static_cast<std::size_t>(v)];
  if (f == 0) ws.marked.push_back(v);  // within the reserved capacity
  f = static_cast<std::uint8_t>(f | bit);
}

// A KSP spur is the node where an alternative path leaves an earlier path's prefix. Spur rules
// apply to forward searches only: excluded prefix nodes (kExcluded) prevent cycles. At the spur,
// multi_edge = true (Expand) bans the listed edges; false (Collapse) bans every edge to a head
// marked kBannedHead.
struct Spur {
  NodeId node {-1};
  std::span<const EdgeId> banned_edges {};
  bool ban_heads {false};
};

// One search over validated inputs. `roots`: allowed, sorted, unique, non-empty; with `targets`
// the live targets carry kTarget and the search stops at the first one popped.
struct SearchIn {
  const Graph* g {nullptr};
  Gate gate {};
  Direction direction {Direction::Forward};
  SpfQueue queue {SpfQueue::Auto};
  std::span<const NodeId> roots {};
  bool targets {false};
  bool lists {true};
  const Spur* spur {nullptr};
};

template <Direction D, class Q>
Cost loop(Workspace& ws, Q& q, const SearchIn& in) {
  const Graph& g = *in.g;
  const Cost* const cost = g.cost().data();
  const Amount* const cap = g.capacity_units().data();
  const NodeId* const edst = g.edge_dst().data();
  const std::int32_t* const row = g.row_offsets().data();
  const std::int32_t* const in_row = g.in_row_offsets().data();
  const EdgeId* const in_adj = g.in_adj_edge_index().data();
  const NodeId* const in_col = detail::GraphAccess::in_col_indices(g).data();  // in_col[s] = esrc[in_adj[s]]
  const std::uint8_t* const allowed = in.gate.elig != nullptr ? in.gate.elig->edge_allowed().data() : nullptr;
  const Amount* const res = in.gate.residual.empty() ? nullptr : in.gate.residual.data();
  const Amount minr = in.gate.min_residual;
  // The node arrays are not resized during a lease: raw pointers keep the compiler from
  // reloading the vector bases around every store into the thread-local workspace.
  Cost* const dist = ws.dist.data();
  std::int32_t* const head = ws.head.data();
  const std::uint8_t* const flags = ws.flags.data();
  const bool targets = in.targets;
  const bool lists = in.lists;
  const Spur* const spur = in.spur;
  auto candidate = [&](std::size_t e) noexcept {
    return (allowed == nullptr || allowed[e] != 0) && cap[e] >= 1 && (res == nullptr || res[e] >= minr);
  };
  // Records a node in `touched` before its first write (never allocates) and a tight edge in the
  // pool before it is linked (a failed growth changes nothing).
  auto relax = [&](NodeId x, Cost nd, EdgeId e) {
    const auto xi = static_cast<std::size_t>(x);
    const Cost dx = dist[xi];
    if (nd < dx) {
      if (lists) ws.pool.push_back({e, -1});
      if (dx == kInf) ws.touched.push_back(x);
      if (lists) head[xi] = static_cast<std::int32_t>(ws.pool.size() - 1);
      dist[xi] = nd;
      q.push(nd, x);
    } else if (lists && nd == dx) {
      ws.pool.push_back({e, head[xi]});
      head[xi] = static_cast<std::int32_t>(ws.pool.size() - 1);
    }
  };
  Cost d = 0;
  NodeId u = 0;
  while (q.pop(d, u)) {
    const auto ui = static_cast<std::size_t>(u);
    if (d != dist[ui]) continue;  // a stale entry: the node improved after this push
    if (targets && (flags[ui] & kTarget) != 0) return d;  // the nearest targets' distance
    if constexpr (D == Direction::Forward) {
      const bool at_spur = spur != nullptr && u == spur->node;
      for (auto e = static_cast<std::size_t>(row[ui]); e < static_cast<std::size_t>(row[ui + 1]); ++e) {
        if (!candidate(e)) continue;
        const NodeId v = edst[e];
        if (spur != nullptr) {
          const std::uint8_t fv = flags[static_cast<std::size_t>(v)];
          if ((fv & kExcluded) != 0) continue;
          if (at_spur) {
            if (spur->ban_heads ? (fv & kBannedHead) != 0
                                : std::find(spur->banned_edges.begin(), spur->banned_edges.end(),
                                            static_cast<EdgeId>(e)) != spur->banned_edges.end()) {
              continue;
            }
          }
        }
        relax(v, d + cost[e], static_cast<EdgeId>(e));
      }
    } else {
      for (auto s = static_cast<std::size_t>(in_row[ui]); s < static_cast<std::size_t>(in_row[ui + 1]); ++s) {
        const auto e = static_cast<std::size_t>(in_adj[s]);
        if (!candidate(e)) continue;
        relax(in_col[s], d + cost[e], static_cast<EdgeId>(e));
      }
    }
  }
  return kInf;
}

// Runs one search into the leased workspace; returns the nearest targets' distance d* (kInf
// without targets or when none is reached).
Cost run(Workspace& ws, const SearchIn& in) {
  for (NodeId r : in.roots) {
    ws.touched.push_back(r);  // before the write; within the reserved capacity
    ws.dist[static_cast<std::size_t>(r)] = 0;
  }
  const Graph& g = *in.g;
  const bool forward = in.direction == Direction::Forward;
  if (uses_buckets(g, in.queue)) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
    ++thread_queue_counts().buckets;
#endif
    ws.buckets.begin(bucket_width(g), std::max<Cost>(1, g.cost_gcd()), ws.dist.size(), in.roots);
    return forward ? loop<Direction::Forward>(ws, ws.buckets, in) : loop<Direction::Reverse>(ws, ws.buckets, in);
  }
  ws.heap.begin(in.roots);
  return forward ? loop<Direction::Forward>(ws, ws.heap, in) : loop<Direction::Reverse>(ws, ws.heap, in);
}

// How a list's entries are read after a search: the list of a node holds its tight in-edges
// (Forward: the other endpoint is the edge's source) or its tight out-edges toward the roots
// (Reverse: the other endpoint is the edge's destination); Collapse keeps one member per pair.
struct Emit {
  const Graph* g {nullptr};
  bool forward {true};
  bool collapse {false};
  bool prefer {false};
  PriorityFn prio {};
};

// Appends node x's list to `out` with Collapse applied. The entries of one pair are consecutive
// in the list (each node's scan inserts them one after another; they are the pair's cheapest
// candidates, consecutive in both adjacency orders), so one pass groups them by the other endpoint.
void read_list(const Workspace& ws, const Emit& em, NodeId x, std::vector<EdgeId>& out) {
  const NodeId* const other = em.forward ? em.g->edge_src().data() : em.g->edge_dst().data();
  std::int32_t k = ws.head[static_cast<std::size_t>(x)];
  if (!em.collapse) {
    for (; k >= 0; k = ws.pool[static_cast<std::size_t>(k)].next) out.push_back(ws.pool[static_cast<std::size_t>(k)].edge);
    return;
  }
  while (k >= 0) {
    EdgeId best = ws.pool[static_cast<std::size_t>(k)].edge;
    const NodeId y = other[static_cast<std::size_t>(best)];
    k = ws.pool[static_cast<std::size_t>(k)].next;
    for (; k >= 0; k = ws.pool[static_cast<std::size_t>(k)].next) {
      const EdgeId e = ws.pool[static_cast<std::size_t>(k)].edge;
      if (other[static_cast<std::size_t>(e)] != y) break;
      if (em.prefer) {
        const int c = compare(em.prio(e), em.prio(best));
        if (c > 0 || (c == 0 && e < best)) best = e;
      } else if (e < best) {
        best = e;
      }
    }
    out.push_back(best);
  }
}

// Sorts a short run ascending (a node's list).
void sort_small(std::vector<EdgeId>& v, std::size_t from) {
  if (v.size() - from > 1) std::sort(v.begin() + static_cast<std::ptrdiff_t>(from), v.end());
}

// The full search's DAG from the lists (one dense pass over the nodes; a full search costs O(N)
// for its costs anyway): Forward lists are per child, so the entries come grouped by child and
// the ascending order by a counting pass over the parents; Reverse lists are per parent, so the
// ascending order comes directly and the grouping by a counting pass over the children.
PredDAG emit_full(Workspace& ws, const Emit& em) {
  const Graph& g = *em.g;
  const auto n = static_cast<std::size_t>(g.num_nodes());
  const auto esrc = g.edge_src();
  const auto edst = g.edge_dst();
  std::vector<NodeId> nodes;
  std::vector<std::int32_t> offsets {0};
  std::vector<NodeId> parents;
  std::vector<EdgeId> via;
  std::vector<EdgeId> asc;
  auto& list = ws.list;
  list.clear();
  if (em.forward) {
    for (std::size_t v = 0; v < n; ++v) {
      if (ws.head[v] < 0) continue;
      const std::size_t from = list.size();
      read_list(ws, em, static_cast<NodeId>(v), list);
      sort_small(list, from);
      nodes.push_back(static_cast<NodeId>(v));
      offsets.push_back(static_cast<std::int32_t>(list.size()));
    }
    via = list;
    parents.resize(via.size());
    for (std::size_t k = 0; k < via.size(); ++k) parents[k] = esrc[static_cast<std::size_t>(via[k])];
    // Ascending EdgeId = grouped by source (EdgeIds rank by source first), ascending within one.
    std::vector<std::int32_t> slot(n + 1, 0);
    for (NodeId p : parents) ++slot[static_cast<std::size_t>(p) + 1];
    for (std::size_t v = 1; v <= n; ++v) slot[v] += slot[v - 1];
    asc.resize(via.size());
    for (std::size_t k = 0; k < via.size(); ++k) {
      asc[static_cast<std::size_t>(slot[static_cast<std::size_t>(parents[k])]++)] = via[k];
    }
    for (std::size_t v = 0, from = 0; v < n; ++v) {
      const auto to = static_cast<std::size_t>(slot[v]);
      if (to - from > 1) std::sort(asc.begin() + static_cast<std::ptrdiff_t>(from), asc.begin() + static_cast<std::ptrdiff_t>(to));
      from = to;
    }
  } else {
    const auto row = g.row_offsets();
    for (std::size_t u = 0; u < n; ++u) {
      if (ws.head[u] < 0) continue;
      if (!em.collapse && row[u + 1] - row[u] <= 64) {
        // The list holds distinct out-edges of u, ids in [row[u], row[u + 1]): one bit per id gives
        // them in ascending order without a sort.
        std::uint64_t bits = 0;
        for (std::int32_t k = ws.head[u]; k >= 0; k = ws.pool[static_cast<std::size_t>(k)].next) {
          bits |= std::uint64_t {1} << (ws.pool[static_cast<std::size_t>(k)].edge - row[u]);
        }
        for (; bits != 0; bits &= bits - 1) list.push_back(row[u] + std::countr_zero(bits));
        continue;
      }
      const std::size_t from = list.size();
      read_list(ws, em, static_cast<NodeId>(u), list);
      sort_small(list, from);
    }
    asc = list;  // grouped by source in ascending NodeId, ascending within each: ascending EdgeId
    std::vector<std::int32_t> slot(n + 1, 0);
    for (EdgeId e : asc) ++slot[static_cast<std::size_t>(edst[static_cast<std::size_t>(e)]) + 1];
    for (std::size_t v = 1; v <= n; ++v) slot[v] += slot[v - 1];
    via.resize(asc.size());
    for (EdgeId e : asc) via[static_cast<std::size_t>(slot[static_cast<std::size_t>(edst[static_cast<std::size_t>(e)])]++)] = e;
    parents.resize(via.size());
    for (std::size_t k = 0; k < via.size(); ++k) parents[k] = esrc[static_cast<std::size_t>(via[k])];
    for (std::size_t v = 0, from = 0; v < n; ++v) {
      const auto to = static_cast<std::size_t>(slot[v]);  // the end of child v's entries
      if (to > from) {
        nodes.push_back(static_cast<NodeId>(v));
        offsets.push_back(static_cast<std::int32_t>(to));
      }
      from = to;
    }
  }
  return PredDAGBuilder::of_parts(g, std::move(nodes), std::move(offsets), std::move(parents), std::move(via),
                                  std::move(asc));
}

// The walk back from the targets at d* along the lists (anycast pruning): every tight edge on a
// walk from a root to a nearest target, Collapse applied, into ws.pruned (any order); the nodes on
// those walks carry kKeep and are listed in ws.kept.
void prune(Workspace& ws, const Emit& em, std::span<const NodeId> targets, Cost dstar) {
  const NodeId* const other = em.forward ? em.g->edge_src().data() : em.g->edge_dst().data();
  ws.pruned.clear();
  ws.kept.clear();
  ws.stack.clear();
  for (NodeId t : targets) {
    if (ws.dist[static_cast<std::size_t>(t)] != dstar) continue;
    mark(ws, t, kKeep);
    ws.kept.push_back(t);
    ws.stack.push_back(t);
  }
  while (!ws.stack.empty()) {
    const NodeId x = ws.stack.back();
    ws.stack.pop_back();
    const std::size_t from = ws.pruned.size();
    read_list(ws, em, x, ws.pruned);
    for (std::size_t k = from; k < ws.pruned.size(); ++k) {
      const NodeId y = other[static_cast<std::size_t>(ws.pruned[k])];
      if ((ws.flags[static_cast<std::size_t>(y)] & kKeep) != 0) continue;
      mark(ws, y, kKeep);
      ws.kept.push_back(y);
      ws.stack.push_back(y);
    }
  }
}

// Choose one path from tight DAG `edges` (ws.pruned, any order) and `nodes` (ws.kept, kKeep).
// pot(v) is distance from the starts: 0 there, `total` at the ends. Increasing pot orders
// every predecessor before its child (topological order).
// Under PreferHigherPriority (`prio`), first maximise the smallest edge priority on a walk
// (its bottleneck), capped by `cap` for a KSP prefix; call this B*. Keep edges with priority >= B*.
// Then compute minimum remaining hops h from the ends. Start where h is smallest and repeatedly
// take the lowest EdgeId with h(v) = h(u) - 1. Decreasing h prevents cycles; choosing the smallest
// eligible id at each position gives the smallest EdgeId sequence among minimum-hop paths.
struct Walk {
  std::vector<EdgeId> edges {};  // in walk order
  Cost cost {0};
  Priority bottleneck {0, 1, true};  // the minimum priority along the walk (prio given)
};

template <class Pot>
Walk extract(Workspace& ws, const Graph& g, Pot pot, Cost total, const PriorityFn* prio, const Priority* cap) {
  auto& nodes = ws.kept;
  auto& edges = ws.pruned;
  const auto esrc = g.edge_src();
  const auto edst = g.edge_dst();
  const auto cost = g.cost();
  std::sort(nodes.begin(), nodes.end(), [&](NodeId a, NodeId b) {
    const Cost pa = pot(a);
    const Cost pb = pot(b);
    return pa != pb ? pa < pb : a < b;
  });
  std::sort(edges.begin(), edges.end());
  const std::size_t k = nodes.size();
  for (std::size_t i = 0; i < k; ++i) ws.loc[static_cast<std::size_t>(nodes[i])] = static_cast<std::int32_t>(i);
  auto local = [&](NodeId v) { return static_cast<std::size_t>(ws.loc[static_cast<std::size_t>(v)]); };
  // Local CSRs over the edges in ascending EdgeId: out-lists by tail, in-lists by head.
  ws.off_out.assign(k + 1, 0);
  ws.off_in.assign(k + 1, 0);
  for (EdgeId e : edges) {
    ++ws.off_out[local(esrc[static_cast<std::size_t>(e)]) + 1];
    ++ws.off_in[local(edst[static_cast<std::size_t>(e)]) + 1];
  }
  for (std::size_t i = 1; i <= k; ++i) {
    ws.off_out[i] += ws.off_out[i - 1];
    ws.off_in[i] += ws.off_in[i - 1];
  }
  ws.idx_out.resize(edges.size());
  ws.idx_in.resize(edges.size());
  ws.hops.assign(k, 0);  // reused as the fill cursors, then as h
  for (std::size_t p = 0; p < edges.size(); ++p) {
    const std::size_t t = local(esrc[static_cast<std::size_t>(edges[p])]);
    ws.idx_out[static_cast<std::size_t>(ws.off_out[t] + ws.hops[t]++)] = static_cast<std::int32_t>(p);
  }
  std::fill(ws.hops.begin(), ws.hops.end(), 0);
  for (std::size_t p = 0; p < edges.size(); ++p) {
    const std::size_t h = local(edst[static_cast<std::size_t>(edges[p])]);
    ws.idx_in[static_cast<std::size_t>(ws.off_in[h] + ws.hops[h]++)] = static_cast<std::int32_t>(p);
  }
  ws.keep_edge.assign(edges.size(), 1);
  if (prio != nullptr) {
    ws.best.assign(k, Priority {});
    ws.has_best.assign(k, 0);
    for (std::size_t i = 0; i < k && pot(nodes[i]) == 0; ++i) {
      ws.best[i] = Priority {0, 1, true};
      ws.has_best[i] = 1;
    }
    for (std::size_t i = 0; i < k; ++i) {
      if (ws.has_best[i] == 0) continue;
      for (auto s = ws.off_out[i]; s < ws.off_out[i + 1]; ++s) {
        const EdgeId e = edges[static_cast<std::size_t>(ws.idx_out[static_cast<std::size_t>(s)])];
        Priority c = (*prio)(e);
        if (compare(ws.best[i], c) < 0) c = ws.best[i];
        const std::size_t j = local(edst[static_cast<std::size_t>(e)]);
        if (ws.has_best[j] == 0 || compare(c, ws.best[j]) > 0) {
          ws.best[j] = c;
          ws.has_best[j] = 1;
        }
      }
    }
    std::optional<Priority> b_star;
    for (std::size_t i = k; i-- > 0 && pot(nodes[i]) == total;) {
      if (ws.has_best[i] != 0 && (!b_star.has_value() || compare(ws.best[i], *b_star) > 0)) b_star = ws.best[i];
    }
    if (!b_star.has_value()) throw std::logic_error("canonical path: no walk reaches an end of the tight DAG");
    if (cap != nullptr && compare(*cap, *b_star) < 0) b_star = *cap;
    for (std::size_t p = 0; p < edges.size(); ++p) ws.keep_edge[p] = compare((*prio)(edges[p]), *b_star) >= 0 ? 1 : 0;
  }
  // h: the hop count to the nearest end inside the (restricted) DAG; -1 where none is reached.
  auto& h = ws.hops;
  std::fill(h.begin(), h.end(), -1);
  auto& queue = ws.stack;
  queue.clear();
  for (std::size_t i = k; i-- > 0 && pot(nodes[i]) == total;) {
    h[i] = 0;
    queue.push_back(nodes[i]);
  }
  for (std::size_t qi = 0; qi < queue.size(); ++qi) {
    const std::size_t x = local(queue[qi]);
    for (auto s = ws.off_in[x]; s < ws.off_in[x + 1]; ++s) {
      const auto p = static_cast<std::size_t>(ws.idx_in[static_cast<std::size_t>(s)]);
      if (ws.keep_edge[p] == 0) continue;
      const std::size_t y = local(esrc[static_cast<std::size_t>(edges[p])]);
      if (h[y] < 0) {
        h[y] = h[x] + 1;
        queue.push_back(nodes[y]);
      }
    }
  }
  // The first qualifying out-edge of a node (out-lists ascending: its smallest), or -1.
  auto next_edge = [&](std::size_t x) -> EdgeId {
    for (auto s = ws.off_out[x]; s < ws.off_out[x + 1]; ++s) {
      const auto p = static_cast<std::size_t>(ws.idx_out[static_cast<std::size_t>(s)]);
      if (ws.keep_edge[p] == 0) continue;
      const EdgeId e = edges[p];
      if (h[local(edst[static_cast<std::size_t>(e)])] == h[x] - 1) return e;
    }
    return -1;
  };
  std::int32_t hmin = -1;
  for (std::size_t i = 0; i < k && pot(nodes[i]) == 0; ++i) {
    if (h[i] >= 0 && (hmin < 0 || h[i] < hmin)) hmin = h[i];
  }
  if (hmin <= 0) throw std::logic_error("canonical path: no start reaches an end of the restricted DAG");
  EdgeId pick = -1;
  for (std::size_t i = 0; i < k && pot(nodes[i]) == 0; ++i) {
    if (h[i] != hmin) continue;
    const EdgeId e = next_edge(i);
    if (e >= 0 && (pick < 0 || e < pick)) pick = e;
  }
  Walk out;
  out.edges.reserve(static_cast<std::size_t>(hmin));
  while (pick >= 0) {
    out.edges.push_back(pick);
    out.cost += cost[static_cast<std::size_t>(pick)];
    if (prio != nullptr) {
      const Priority p = (*prio)(pick);
      if (compare(p, out.bottleneck) < 0) out.bottleneck = p;
    }
    const std::size_t v = local(edst[static_cast<std::size_t>(pick)]);
    if (h[v] == 0) break;
    pick = next_edge(v);
  }
  if (out.edges.size() != static_cast<std::size_t>(hmin) || out.cost != total) {
    throw std::logic_error("canonical path: the greedy walk is not a shortest walk of the DAG");
  }
  return out;
}

// The live (allowed) members of a sorted node set.
std::vector<NodeId> live(std::span<const NodeId> set, const Eligibility* elig) {
  std::vector<NodeId> out;
  out.reserve(set.size());
  for (NodeId v : set) {
    if (elig == nullptr || elig->node_ok(v)) out.push_back(v);
  }
  return out;
}

Emit emit_of(const Graph& g, Direction direction, const PathSelection& selection, const PriorityFn& prio) {
  Emit em;
  em.g = &g;
  em.forward = direction == Direction::Forward;
  em.collapse = !selection.multi_edge;
  em.prefer = selection.tie_break == EdgeTieBreak::PreferHigherPriority;
  em.prio = prio;
  return em;
}

// One KSP spur search: the canonical walk from `spur` to `dst` over the
// candidates of `gate`, without the excluded prefix nodes, with the spur node's banned
// transitions (Expand: the banned edges; Collapse: every edge to a banned head, no re-selection),
// B* capped by `cap` (the root prefix's bottleneck). nullopt when dst is unreachable.
std::optional<Walk> spur_walk(const Graph& g, const Gate& gate, const PathSelection& selection, const PriorityFn& prio,
                              NodeId spur, NodeId dst, std::span<const NodeId> excluded, std::span<const EdgeId> banned,
                              const Priority* cap) {
  Lease lease(static_cast<std::size_t>(g.num_nodes()));
  Workspace& ws = lease.ws();
  mark(ws, dst, kTarget);
  for (NodeId v : excluded) mark(ws, v, kExcluded);
  Spur rule;
  rule.node = spur;
  if (selection.multi_edge) {
    rule.banned_edges = banned;
  } else {
    rule.ban_heads = true;
    for (EdgeId e : banned) mark(ws, g.edge_dst()[static_cast<std::size_t>(e)], kBannedHead);
  }
  const NodeId roots[1] = {spur};
  SearchIn in;
  in.g = &g;
  in.gate = gate;
  in.roots = roots;
  in.targets = true;
  in.spur = &rule;
  const Cost dstar = run(ws, in);
  if (dstar == kInf) return std::nullopt;
  const Emit em = emit_of(g, Direction::Forward, selection, prio);
  const NodeId targets[1] = {dst};
  prune(ws, em, targets, dstar);
  const bool prefer = selection.tie_break == EdgeTieBreak::PreferHigherPriority;
  return extract(ws, g, [&](NodeId v) { return ws.dist[static_cast<std::size_t>(v)]; }, dstar,
                 prefer ? &prio : nullptr, cap);
}


// The walks of resolve_to_paths (validated inputs): see the public function.
std::vector<Path> enumerate_paths(const Graph& g, const PredDAG& dag, NodeId src, NodeId dst, bool split,
                                  std::optional<std::int64_t> max_paths) {
  const auto children = PredDAGBuilder::member_nodes(dag);
  const auto offsets = PredDAGBuilder::member_offsets(dag);
  const auto parents = dag.parents();
  const auto via = dag.via_edges();
  const auto cost = g.cost();
  std::vector<Path> out;
  // The entries of a child: [offsets[i], offsets[i + 1]) for children[i], none for another node.
  auto entries_of = [&](NodeId v) -> std::pair<std::int32_t, std::int32_t> {
    const auto it = std::lower_bound(children.begin(), children.end(), v);
    if (it == children.end() || *it != v) return {0, 0};
    const auto i = static_cast<std::size_t>(it - children.begin());
    return {offsets[i], offsets[i + 1]};
  };
  Lease lease(static_cast<std::size_t>(g.num_nodes()));
  Workspace& ws = lease.ws();
  // The nodes that reach dst inside the DAG (the walk back over the entries), local index by
  // discovery.
  auto& kept = ws.kept;
  kept.clear();
  ws.stack.clear();
  mark(ws, dst, kKeep);
  kept.push_back(dst);
  ws.stack.push_back(dst);
  while (!ws.stack.empty()) {
    const NodeId v = ws.stack.back();
    ws.stack.pop_back();
    const auto [lo, hi] = entries_of(v);
    for (auto k = lo; k < hi; ++k) {
      const NodeId u = parents[static_cast<std::size_t>(k)];
      if ((ws.flags[static_cast<std::size_t>(u)] & kKeep) != 0) continue;
      mark(ws, u, kKeep);
      kept.push_back(u);
      ws.stack.push_back(u);
    }
  }
  if ((ws.flags[static_cast<std::size_t>(src)] & kKeep) == 0) return out;  // src is not on a walk to dst
  for (std::size_t i = 0; i < kept.size(); ++i) ws.loc[static_cast<std::size_t>(kept[i])] = static_cast<std::int32_t>(i);
  auto local = [&](NodeId v) { return static_cast<std::size_t>(ws.loc[static_cast<std::size_t>(v)]); };
  // The hops: per child, its entries come in ascending via, so grouped by parent (EdgeIds rank by
  // source first); a hop is one parent's run (its members ascending, its cost the cheapest, its
  // key the lowest id), or one member when split.
  struct Hop {
    std::int32_t tail;
    std::int32_t head;
    std::int32_t first;  // the members are via[first .. first + count)
    std::int32_t count;
    Cost cost;
    EdgeId key;
  };
  std::vector<Hop> hops;
  for (NodeId v : kept) {
    const auto [lo, hi] = entries_of(v);
    for (auto k = lo; k < hi;) {
      const NodeId u = parents[static_cast<std::size_t>(k)];
      auto end = k + 1;
      while (end < hi && parents[static_cast<std::size_t>(end)] == u) ++end;
      if (split) {
        for (auto m = k; m < end; ++m) {
          const EdgeId e = via[static_cast<std::size_t>(m)];
          hops.push_back({static_cast<std::int32_t>(local(u)), static_cast<std::int32_t>(local(v)), m, 1,
                          cost[static_cast<std::size_t>(e)], e});
        }
      } else {
        Cost c = cost[static_cast<std::size_t>(via[static_cast<std::size_t>(k)])];
        for (auto m = k + 1; m < end; ++m) c = std::min(c, cost[static_cast<std::size_t>(via[static_cast<std::size_t>(m)])]);
        hops.push_back({static_cast<std::int32_t>(local(u)), static_cast<std::int32_t>(local(v)), k, end - k, c,
                        via[static_cast<std::size_t>(k)]});
      }
      k = end;
    }
  }
  const std::size_t n = kept.size();
  std::vector<std::int32_t> out_off(n + 1, 0);
  for (const auto& h : hops) ++out_off[static_cast<std::size_t>(h.tail) + 1];
  for (std::size_t i = 1; i <= n; ++i) out_off[i] += out_off[i - 1];
  std::vector<std::int32_t> out_hop(hops.size());
  {
    std::vector<std::int32_t> fill(out_off.begin(), out_off.end() - 1);
    for (std::size_t h = 0; h < hops.size(); ++h) {
      out_hop[static_cast<std::size_t>(fill[static_cast<std::size_t>(hops[h].tail)]++)] = static_cast<std::int32_t>(h);
    }
  }
  // A completion is a walk suffix from a node to dst. Per node: the completions found so far
  // (sorted) and a heap holding the next completion through each hop.
  struct Elem {
    Cost cost;
    std::int64_t hops;
    std::int32_t hop;   // -1: the empty completion of dst
    std::int64_t rank;  // the completion of the hop's head this one extends
  };
  struct Next {
    Cost cost;
    std::int64_t hops;
    EdgeId key;
    std::int32_t hop;
    std::int64_t rank;
  };
  auto after = [](const Next& a, const Next& b) {  // a pops after b: (cost, hops, key) ascending
    if (a.cost != b.cost) return a.cost > b.cost;
    if (a.hops != b.hops) return a.hops > b.hops;
    return a.key > b.key;
  };
  std::vector<std::vector<Elem>> list(n);
  std::vector<std::vector<Next>> heap(n);
  std::vector<std::uint8_t> ready(n, 0);
  const std::size_t t = local(dst);
  list[t].push_back({0, 0, -1, -1});
  ready[t] = 1;  // a walk ends at dst: its out-hops never extend a completion
  auto offer = [&](std::size_t x, std::int32_t h, std::int64_t rank) {
    const auto& s = list[static_cast<std::size_t>(hops[static_cast<std::size_t>(h)].head)][static_cast<std::size_t>(rank)];
    heap[x].push_back({hops[static_cast<std::size_t>(h)].cost + s.cost, s.hops + 1, hops[static_cast<std::size_t>(h)].key, h, rank});
    std::push_heap(heap[x].begin(), heap[x].end(), after);
  };
  // get(x, i): whether node x has an i-th completion (computing it if needed), with an explicit
  // frame stack. Stages: 0 entry, 1 the next hop to initialise, 2 back from initialising one,
  // 3 pop the next completion, 4 back from extending the popped hop.
  struct Frame {
    std::size_t x;
    std::int64_t i;
    int stage;
    std::int32_t k;
    std::int64_t r;
  };
  std::vector<Frame> frames;
  auto get = [&](std::size_t x0, std::int64_t i0) {
    frames.assign(1, Frame {x0, i0, 0, 0, 0});
    bool ret = false;
    while (!frames.empty()) {
      const std::size_t top = frames.size() - 1;
      Frame f = frames[top];
      switch (f.stage) {
        case 0:
          if (static_cast<std::int64_t>(list[f.x].size()) > f.i) {
            ret = true;
            frames.pop_back();
          } else if (ready[f.x] == 0) {
            frames[top].stage = 1;
            frames[top].k = out_off[f.x];
          } else {
            frames[top].stage = 3;
          }
          break;
        case 1:
          if (f.k == out_off[f.x + 1]) {
            ready[f.x] = 1;
            frames[top].stage = 3;
          } else {
            frames[top].stage = 2;
            frames.push_back({static_cast<std::size_t>(hops[static_cast<std::size_t>(out_hop[static_cast<std::size_t>(f.k)])].head), 0,
                              0, 0, 0});
          }
          break;
        case 2:
          if (ret) offer(f.x, out_hop[static_cast<std::size_t>(f.k)], 0);
          frames[top].k = f.k + 1;
          frames[top].stage = 1;
          break;
        case 3: {
          auto& hx = heap[f.x];
          if (hx.empty()) {
            ret = false;
            frames.pop_back();
            break;
          }
          std::pop_heap(hx.begin(), hx.end(), after);
          const Next c = hx.back();
          hx.pop_back();
          list[f.x].push_back({c.cost, c.hops, c.hop, c.rank});
          frames[top].stage = 4;
          frames[top].k = c.hop;
          frames[top].r = c.rank + 1;
          frames.push_back({static_cast<std::size_t>(hops[static_cast<std::size_t>(c.hop)].head), c.rank + 1, 0, 0, 0});
          break;
        }
        default:  // 4
          if (ret) offer(f.x, f.k, f.r);
          frames[top].stage = 0;
          break;
      }
    }
    return ret;
  };
  const std::size_t s = local(src);
  for (std::int64_t j = 0; !max_paths.has_value() || j < *max_paths; ++j) {
    if (!get(s, j)) break;
    std::vector<std::vector<EdgeId>> walk;
    std::size_t x = s;
    std::int64_t idx = j;
    while (x != t) {
      const Elem& e = list[x][static_cast<std::size_t>(idx)];
      const Hop& h = hops[static_cast<std::size_t>(e.hop)];
      walk.emplace_back(via.begin() + h.first, via.begin() + h.first + h.count);
      x = static_cast<std::size_t>(h.head);
      idx = e.rank;
    }
    out.push_back(PathBuilder::of_hops(g, std::move(walk)));
  }
  return out;
}

} // namespace

std::optional<SearchRoute> route_search(const Graph& g, const Gate& gate, const PathSelection& selection, NodeId src,
                                        std::span<const NodeId> dests) {
  if (dests.empty()) return std::nullopt;
  const PriorityValues values {&g, selection.priority, gate.residual};
  const PriorityFn prio = PriorityFn::of(values);
  Lease lease(static_cast<std::size_t>(g.num_nodes()));
  Workspace& ws = lease.ws();
  for (NodeId t : dests) mark(ws, t, kTarget);
  const NodeId roots[1] = {src};
  SearchIn in;
  in.g = &g;
  in.gate = gate;
  in.roots = roots;
  in.targets = true;
  const Cost dstar = run(ws, in);
  if (dstar == kInf) return std::nullopt;
  const Emit em = emit_of(g, Direction::Forward, selection, prio);
  prune(ws, em, dests, dstar);
  SearchRoute out;
  out.cost = dstar;
  if (selection.multipath) {
    out.edges = ws.pruned;
  } else {
    const bool prefer = selection.tie_break == EdgeTieBreak::PreferHigherPriority;
    auto walk = extract(ws, g, [&](NodeId v) { return ws.dist[static_cast<std::size_t>(v)]; }, dstar,
                        prefer ? &prio : nullptr, nullptr);
    out.edges = std::move(walk.edges);
  }
  std::sort(out.edges.begin(), out.edges.end());
  return out;
}

SearchRoute support_route(const Graph& g, std::span<const EdgeId> support, NodeId src, std::span<const NodeId> dests,
                          const PathSelection& selection, PriorityFn view) {
  const auto row = g.row_offsets();
  const auto edst = g.edge_dst();
  const auto esrc = g.edge_src();
  const auto cost = g.cost();
  Lease lease(static_cast<std::size_t>(g.num_nodes()));
  Workspace& ws = lease.ws();
  for (NodeId t : dests) mark(ws, t, kTarget);
  // The part of the support reachable from src: every walk of a support from one node has one
  // cost to each node it reaches (the support is a reverse multipath DAG), so the order of the
  // walk assigns each node its distance from src exactly.
  auto& route = ws.pruned;
  route.clear();
  ws.kept.clear();
  ws.stack.clear();
  mark(ws, src, kKeep);
  ws.kept.push_back(src);
  ws.touched.push_back(src);
  ws.dist[static_cast<std::size_t>(src)] = 0;
  ws.stack.push_back(src);
  Cost total = kInf;
  while (!ws.stack.empty()) {
    const NodeId u = ws.stack.back();
    ws.stack.pop_back();
    const auto ui = static_cast<std::size_t>(u);
    if ((ws.flags[ui] & kTarget) != 0) total = ws.dist[ui];
    auto it = std::lower_bound(support.begin(), support.end(), row[ui]);
    for (; it != support.end() && *it < row[ui + 1]; ++it) {
      route.push_back(*it);
      const NodeId v = edst[static_cast<std::size_t>(*it)];
      const auto vi = static_cast<std::size_t>(v);
      if ((ws.flags[vi] & kKeep) != 0) continue;
      mark(ws, v, kKeep);
      ws.kept.push_back(v);
      ws.touched.push_back(v);
      ws.dist[vi] = ws.dist[ui] + cost[static_cast<std::size_t>(*it)];
      ws.stack.push_back(v);
    }
  }
  SearchRoute out;
  if (total == kInf) return out;  // no destination reached (the roster excludes such a source)
  std::sort(route.begin(), route.end());
  const bool prefer = selection.tie_break == EdgeTieBreak::PreferHigherPriority;
  if (!selection.multi_edge) {
    // One member per (u, v) by select_adjacency on the view (no gating: the support's members are
    // all statically eligible). The members of a pair in the support have one cost, so they are
    // consecutive in id order (EdgeIds rank by (src, cost, dst, input)).
    std::size_t w = 0;
    for (std::size_t i = 0; i < route.size();) {
      const auto ei = static_cast<std::size_t>(route[i]);
      EdgeId best = route[i];
      std::size_t j = i + 1;
      for (; j < route.size() && esrc[static_cast<std::size_t>(route[j])] == esrc[ei] &&
             edst[static_cast<std::size_t>(route[j])] == edst[ei];
           ++j) {
        if (replaces(g, route[j], best, prefer, view)) best = route[j];
      }
      route[w++] = best;
      i = j;
    }
    route.resize(w);
  }
  out.cost = total;
  if (selection.multipath) {
    out.edges = route;
    return out;
  }
  auto walk = extract(ws, g, [&](NodeId v) { return ws.dist[static_cast<std::size_t>(v)]; }, total,
                      prefer ? &view : nullptr, nullptr);
  out.edges = std::move(walk.edges);
  std::sort(out.edges.begin(), out.edges.end());
  return out;
}

PredDAG static_support(const Graph& g, std::span<const NodeId> roots, Direction direction, const Eligibility* elig) {
  const auto live_roots = live(roots, elig);
  if (live_roots.empty()) return PredDAGBuilder::empty(g);
  Lease lease(static_cast<std::size_t>(g.num_nodes()));
  Workspace& ws = lease.ws();
  SearchIn in;
  in.g = &g;
  in.gate.elig = elig;
  in.direction = direction;
  in.roots = live_roots;
  (void)run(ws, in);
  Emit em;
  em.g = &g;
  em.forward = direction == Direction::Forward;
  return emit_full(ws, em);
}

std::vector<Cost> static_costs(const Graph& g, const Eligibility* elig, std::span<const NodeId> dests,
                               std::span<const NodeId> sources) {
  std::vector<Cost> out(sources.size(), kInf);
  const auto roots = live(dests, elig);
  if (roots.empty()) return out;
  Lease lease(static_cast<std::size_t>(g.num_nodes()));
  Workspace& ws = lease.ws();
  SearchIn in;
  in.g = &g;
  in.gate.elig = elig;
  in.direction = Direction::Reverse;
  in.roots = roots;
  in.lists = false;
  (void)run(ws, in);
  for (std::size_t i = 0; i < sources.size(); ++i) out[i] = ws.dist[static_cast<std::size_t>(sources[i])];
  return out;
}

void check_ceiling(const std::optional<CostCeiling>& ceiling) {
  if (!ceiling.has_value()) return;
  if (const auto* a = std::get_if<Absolute>(&*ceiling)) {
    if (a->cost < 0) throw std::invalid_argument("max_path_cost Absolute cost must be >= 0 (got " + std::to_string(a->cost) + ")");
    return;
  }
  const double f = std::get<Stretch>(*ceiling).factor;
  if (!(f >= 1.0) || !std::isfinite(f)) {
    throw std::invalid_argument("max_path_cost Stretch factor must be finite and in [1, inf) (got " + std::to_string(f) + ")");
  }
}

} // namespace netgraph::core::detail

namespace netgraph::core {

// The selection of a search (options.hpp): its two enum-typed fields inside their enumerations.
PathSelection::PathSelection(EdgeTieBreak tie_break_, EdgePriority priority_, bool multipath_, bool multi_edge_)
    : tie_break(tie_break_), priority(priority_), multipath(multipath_), multi_edge(multi_edge_) {
  detail::check_member(tie_break, "tie_break");
  detail::check_member(priority, "priority");
}

SpfOptions::SpfOptions(Direction direction_, PathSelection selection_) : direction(direction_), selection(selection_) {
  detail::check_member(direction, "direction");
}

SpfResult spf(const Graph& g, std::span<const NodeId> roots, const SpfOptions& opts, SpfQueue queue) {
  using namespace detail;
  check_member(opts.direction, "direction");
  check_selection(opts.selection);
  check_member(queue, "queue");
  const auto root_set = node_set(roots, g.num_nodes(), "roots");
  std::vector<NodeId> target_set;
  if (opts.targets.has_value()) {
    target_set = node_set(*opts.targets, g.num_nodes(), "targets");
    check_disjoint(root_set, target_set, "roots", "targets");
  } else if (!opts.selection.multipath) {
    throw std::invalid_argument("multipath = false requires targets (the single path joins a root to a target)");
  }
  const std::span<const Amount> residual = opts.residual.value_or(std::span<const Amount> {});
  if (opts.residual.has_value()) check_edge_amounts(g, residual, "residual");
  if (opts.min_residual < 1) {
    throw std::invalid_argument("min_residual must be at least 1 unit (got " + std::to_string(opts.min_residual) + ")");
  }
  if (opts.elig != nullptr && opts.elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  const PriorityValues values {&g, opts.selection.priority, residual};
  const PriorityFn prio = PriorityFn::of(values);
  const Emit em = emit_of(g, opts.direction, opts.selection, prio);
  const auto live_roots = live(root_set, opts.elig);
  SearchIn in;
  in.g = &g;
  in.gate = Gate {opts.elig, residual, opts.min_residual};
  in.direction = opts.direction;
  in.queue = queue;
  in.roots = live_roots;
  if (!opts.targets.has_value()) {
    const auto n = static_cast<std::size_t>(g.num_nodes());
    if (live_roots.empty()) {
      return SpfResult {std::vector<Cost>(n, kUnreachable), PredDAGBuilder::empty(g), std::nullopt};
    }
    Lease lease(n);
    Workspace& ws = lease.ws();
    (void)run(ws, in);
    std::vector<Cost> costs(ws.dist.begin(), ws.dist.begin() + static_cast<std::ptrdiff_t>(n));
    PredDAG dag = emit_full(ws, em);
    return SpfResult {std::move(costs), std::move(dag), std::nullopt};
  }
  const auto live_targets = live(target_set, opts.elig);
  if (live_roots.empty() || live_targets.empty()) return SpfResult {std::nullopt, PredDAGBuilder::empty(g), kUnreachable};
  Lease lease(static_cast<std::size_t>(g.num_nodes()));
  Workspace& ws = lease.ws();
  for (NodeId t : live_targets) mark(ws, t, kTarget);
  in.targets = true;
  const Cost dstar = run(ws, in);
  if (dstar == kInf) return SpfResult {std::nullopt, PredDAGBuilder::empty(g), kUnreachable};
  prune(ws, em, live_targets, dstar);
  if (opts.selection.multipath) return SpfResult {std::nullopt, PredDAGBuilder::of_edges(g, ws.pruned), dstar};
  // The walks run from the walk starts (Forward: the roots; Reverse: the targets) to the ends;
  // the potential is the forward distance from the starts (Reverse: d* - the distance to the roots).
  const bool prefer = opts.selection.tie_break == EdgeTieBreak::PreferHigherPriority;
  Walk walk;
  if (opts.direction == Direction::Forward) {
    walk = extract(ws, g, [&](NodeId v) { return ws.dist[static_cast<std::size_t>(v)]; }, dstar,
                   prefer ? &prio : nullptr, nullptr);
  } else {
    walk = extract(ws, g, [&](NodeId v) { return dstar - ws.dist[static_cast<std::size_t>(v)]; }, dstar,
                   prefer ? &prio : nullptr, nullptr);
  }
  return SpfResult {std::nullopt, PredDAGBuilder::of_edges(g, walk.edges), dstar};
}


KspOptions::KspOptions(std::int64_t k_, PathSelection selection_) : k(k_), selection(selection_) {}

std::vector<std::pair<Cost, Path>> ksp(const Graph& g, NodeId src, NodeId dst, const KspOptions& opts) {
  using namespace detail;
  check_selection(opts.selection);
  const NodeId s1[1] = {src};
  const NodeId t1[1] = {dst};
  check_range(s1, g.num_nodes(), "src");
  check_range(t1, g.num_nodes(), "dst");
  check_disjoint(s1, t1, "src", "dst");
  if (opts.k < 1) throw std::invalid_argument("k must be >= 1 (got " + std::to_string(opts.k) + ")");
  check_ceiling(opts.max_path_cost);
  const std::span<const Amount> residual = opts.residual.value_or(std::span<const Amount> {});
  if (opts.residual.has_value()) check_edge_amounts(g, residual, "residual");
  if (opts.elig != nullptr && opts.elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  std::vector<std::pair<Cost, Path>> out;
  if (opts.elig != nullptr && (!opts.elig->node_ok(src) || !opts.elig->node_ok(dst))) return out;
  // The ceiling: Absolute(c), or Stretch(f) times the static shortest cost of the member graph
  // under the Eligibility (residual ignored, as for the CapacityAware ceiling); no static route
  // gives [].
  const auto* absolute = opts.max_path_cost.has_value() ? std::get_if<Absolute>(&*opts.max_path_cost) : nullptr;
  const auto* stretch = opts.max_path_cost.has_value() ? std::get_if<Stretch>(&*opts.max_path_cost) : nullptr;
  Cost base = 0;
  if (stretch != nullptr) {
    base = static_costs(g, opts.elig, t1, s1)[0];
    if (base == kUnreachable) return out;
  }
  auto within = [&](Cost c) {
    if (absolute != nullptr) return c <= absolute->cost;
    if (stretch != nullptr) return stretch_le(c, base, stretch->factor);
    return true;
  };
  const Gate gate {opts.elig, residual, 1};  // KSP gates a supplied residual at 1 unit
  const PriorityValues values {&g, opts.selection.priority, residual};
  const PriorityFn prio = PriorityFn::of(values);
  const bool prefer = opts.selection.tie_break == EdgeTieBreak::PreferHigherPriority;
  // Yen's method finds alternatives by keeping a prefix and replacing its suffix at a spur.
  // The first path and each spur suffix use the canonical search with that spur's exclusions.
  // Match prefixes by EdgeId sequence. Ban the next edge of every accepted matching path
  // (Expand), or its entire (u, v) transition (Collapse); exclude earlier prefix vertices to
  // prevent cycles. Store candidates once in a set ordered by (cost, decreasing bottleneck
  // under PreferHigherPriority, hops, EdgeId sequence).
  struct Cand {
    Cost cost {0};
    Priority bottleneck {0, 1, true};
    std::vector<EdgeId> edges {};
  };
  struct CandLess {
    bool prefer;
    bool operator()(const Cand& a, const Cand& b) const {
      if (a.cost != b.cost) return a.cost < b.cost;
      if (prefer) {
        const int c = compare(a.bottleneck, b.bottleneck);
        if (c != 0) return c > 0;
      }
      if (a.edges.size() != b.edges.size()) return a.edges.size() < b.edges.size();
      return a.edges < b.edges;
    }
  };
  auto first = spur_walk(g, gate, opts.selection, prio, src, dst, {}, {}, nullptr);
  if (!first.has_value() || !within(first->cost)) return out;
  std::vector<Cand> accepted;
  accepted.push_back(Cand {first->cost, first->bottleneck, std::move(first->edges)});
  std::set<Cand, CandLess> candidates(CandLess {prefer});
  const auto edst = g.edge_dst();
  const auto cost = g.cost();
  std::vector<NodeId> nodes;
  std::vector<EdgeId> banned;
  while (static_cast<std::int64_t>(accepted.size()) < opts.k) {
    const std::vector<EdgeId> last = accepted.back().edges;
    nodes.assign(1, src);
    for (EdgeId e : last) nodes.push_back(edst[static_cast<std::size_t>(e)]);
    Cost prefix_cost = 0;
    Priority prefix_bottleneck {0, 1, true};
    for (std::size_t j = 0; j < last.size(); ++j) {
      banned.clear();
      for (const auto& p : accepted) {
        if (p.edges.size() > j && std::equal(last.begin(), last.begin() + static_cast<std::ptrdiff_t>(j), p.edges.begin())) {
          banned.push_back(p.edges[j]);
        }
      }
      const std::span<const NodeId> excluded(nodes.data(), j);
      auto sp = spur_walk(g, gate, opts.selection, prio, nodes[j], dst, excluded, banned,
                          prefer && j > 0 ? &prefix_bottleneck : nullptr);
      if (sp.has_value()) {
        Cand c;
        c.cost = prefix_cost + sp->cost;
        c.bottleneck = compare(sp->bottleneck, prefix_bottleneck) < 0 ? sp->bottleneck : prefix_bottleneck;
        c.edges.assign(last.begin(), last.begin() + static_cast<std::ptrdiff_t>(j));
        c.edges.insert(c.edges.end(), sp->edges.begin(), sp->edges.end());
        // Never an accepted path: every accepted path with this prefix has its next edge (Expand)
        // or its next transition (Collapse) banned at the spur, and the candidate's next edge is
        // the spur walk's first. Only a pending candidate can be equal; the set keeps one.
        candidates.insert(std::move(c));
      }
      prefix_cost += cost[static_cast<std::size_t>(last[j])];
      if (prefer) {
        const Priority pj = prio(last[j]);
        if (compare(pj, prefix_bottleneck) < 0) prefix_bottleneck = pj;
      }
    }
    if (candidates.empty()) break;
    auto best = candidates.begin();
    if (!within(best->cost)) break;  // every later path costs at least as much
    accepted.push_back(*best);
    candidates.erase(best);
  }
  out.reserve(accepted.size());
  for (auto& p : accepted) out.emplace_back(p.cost, PathBuilder::of_walk(g, p.edges));
  return out;
}

std::vector<Path> resolve_to_paths(const Graph& g, const PredDAG& dag, NodeId src, NodeId dst, bool split_parallel_edges,
                                   std::optional<std::int64_t> max_paths) {
  using namespace detail;
  if (dag.graph_uid() != g.uid()) throw std::invalid_argument("foreign PredDAG: built for another graph");
  const NodeId s1[1] = {src};
  const NodeId t1[1] = {dst};
  check_range(s1, g.num_nodes(), "src");
  check_range(t1, g.num_nodes(), "dst");
  check_disjoint(s1, t1, "src", "dst");
  if (max_paths.has_value() && *max_paths < 1) {
    throw std::invalid_argument("max_paths must be >= 1 (got " + std::to_string(*max_paths) + ")");
  }
  return enumerate_paths(g, dag, src, dst, split_parallel_edges, max_paths);
}

} // namespace netgraph::core
