/*
  Graph: the one graph type on integer units (APPROACH 2.1, 3).

  Construction validates the inputs in the documented order, then builds layout B' by stable
  counting passes (each distributes items into buckets numbered by a small integer key and keeps
  the earlier order inside a bucket). The cost rank of an edge is the position of its cost among
  the K distinct costs, sorted. The id order is (src, cost, dst, input), whose rank is the public
  EdgeId, so the out-edges of u are the id range [row_offsets[u], row_offsets[u + 1]); the
  in-edges of v are listed in (cost, src, input) order. An edgeless graph is its validation and
  two zero-filled offset arrays.
  * N K >= E: three passes (dst, cost rank, src; least significant first) give the id order, and
    two over the ids (cost rank, then dst) the reverse lists: within one (dst, cost) group the
    id order is (src, input).
  * N K < E (combined keys): a pass by dst, then one by (src, rank), give the id order, and one
    pass over the ids by (dst, rank) the reverse lists; the two (N K + 1)-entry histograms live
    in the not-yet-filled capacity and cost arrays.
  The offset arrays serve as their own cursors; a pass that advances them is followed by one
  that walks back, so no pass over the nodes restores them.

  The ranks: each edge first gets its cost's first-appearance ordinal from an open-addressing
  table (fixed mixing, linear probing, load at most 1/4, doubled from 64 up to 2^18 slots, each
  growth rebuilt from the distinct costs in their order); an
  edge whose cost equals the previous edge's reuses that ordinal without a lookup. The K
  distinct costs are then radix-sorted. Hashing stops on the 65,537th distinct cost, or before
  a probe (one table slot examined, the terminating one included, counted across the lookups
  and the growths) would exceed 8E + 8 x 65,536; the ranks then come from a radix sort of all
  costs (16-bit digits). Both paths give the same ranks.

  Memory: every retained array is allocated first, at its final size; until it is filled, each
  serves as scratch. The ranks go into in_col_indices_; the hash table borrows capacity_ (keys) and
  src_ (ordinals) while it has at most E slots, the distinct costs cost_, their sort dst_,
  input_index_ and in_adj_edge_index_; the passes use in_adj_edge_index_ (the middle forward
  pass), capacity_ and cost_ (the K rank cursors or the combined histograms), and the fallback
  input_index_ and in_adj_edge_index_. The only other buffers are a table larger than E slots (its
  own vectors, the old and new ones briefly together at a growth) and the fallback's 2^16 digit
  counts. A table of at most E slots frees no separate buffer that the allocator could keep cached.
*/
#include "netgraph/core/graph.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#include "netgraph/core/detail/graph_access.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/detail/uid.hpp"

namespace netgraph::core {

namespace detail {

std::uint64_t next_uid() noexcept {
  static std::atomic<std::uint64_t> counter {0};
  return counter.fetch_add(1, std::memory_order_relaxed) + 1;
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
ConstructionCounts& thread_construction_counts() noexcept {
  thread_local ConstructionCounts counts;
  return counts;
}
#endif

} // namespace detail

namespace {

using detail::next_uid;

// The exponent k of a resolution q = 2^k with k in [-32, 32].
int resolution_exponent(double q) {
  if (!std::isfinite(q) || q <= 0.0) {
    throw std::invalid_argument("resolution must be a power of two in [2^-32, 2^32]");
  }
  int exp = 0;
  const double mantissa = std::frexp(q, &exp);  // q = mantissa * 2^exp, mantissa in [0.5, 1)
  const int k = exp - 1;
  if (mantissa != 0.5 || k < -32 || k > 32) {
    throw std::invalid_argument("resolution must be a power of two in [2^-32, 2^32]");
  }
  return k;
}

// The cost ranks ---------------------------------------------------------------------------


constexpr std::int64_t kMaxHashedCosts = 65536;  // L

// The cumulative probe budget: take() runs before every table slot is examined and refuses the
// probe that would exceed the limit.
struct ProbeBudget {
  std::int64_t used {0};
  std::int64_t limit {0};
  bool take() noexcept {
    if (used == limit) return false;
    ++used;
    return true;
  }
};

constexpr std::size_t kNoSlot = ~std::size_t {0};

// The distinct-cost table: key 0 marks an empty slot (costs are >= 1 after validation).
struct CostTable {
  std::span<Cost> key;
  std::span<std::int32_t> ord;  // the key's first-appearance ordinal
  std::size_t mask {0};

  // The slot holding c, or the empty slot where c belongs; kNoSlot when the budget refuses a
  // probe. Every slot examined is one probe, the terminating one included.
  std::size_t find(Cost c, ProbeBudget& b) const noexcept {
    std::size_t s = static_cast<std::size_t>(detail::cost_hash(c)) & mask;
    for (;;) {
      if (!b.take()) return kNoSlot;
      if (key[s] == 0 || key[s] == c) return s;
      s = (s + 1) & mask;
    }
  }
};

// Where the table lives: the first `slots` entries of two retained arrays not yet filled (E slots
// each), so the table allocates nothing and leaves no freed block behind (a freed block the
// allocator keeps cached would count in the process's resident peak); only a table larger than E
// slots gets vectors of its own. A new table is empty: growth rebuilds it from the distinct costs.
struct TableStore {
  std::span<Cost> key_room;
  std::span<std::int32_t> ord_room;
  std::vector<Cost> key_own {};
  std::vector<std::int32_t> ord_own {};

  CostTable make(std::size_t slots) {
    if (slots <= key_room.size()) {
      const std::span<Cost> key = key_room.first(slots);
      std::fill(key.begin(), key.end(), Cost {0});
      return {key, ord_room.first(slots), slots - 1};
    }
    key_own.assign(slots, 0);
    ord_own.resize(slots);
    return {key_own, ord_own, slots - 1};
  }
};

// The retained arrays the hash path borrows before they are filled (each E entries, K <= E).
struct HashRooms {
  std::span<Cost> keys;             // the table's keys while it has at most E slots
  std::span<std::int32_t> ords;     // the table's ordinals, the same
  std::span<Cost> distinct;         // the distinct costs in first-appearance order
  std::span<std::int32_t> a, b, c;  // the rank step: ordinal -> rank, and its two position arrays
};

enum class RankPath { Hashed, DistinctLimit, ExhaustedCollection, ExhaustedGrowth };

struct RankOutcome {
  RankPath path {RankPath::Hashed};
  std::int64_t distinct {0};
  std::int64_t table_slots {0};
  std::int64_t probes {0};
  std::int64_t budget {0};
};

// Ranks from a sort of all E costs: a stable LSD radix sort of the input positions by cost,
// 16 bits per pass, then one scan of the sorted positions. Digits above the largest cost's bit
// length are not visited, and a pass whose digit is the same for every cost is skipped.
// `a` and `b` (E slots each) hold the positions; `count` has 2^16 entries.
template <int kBits>
std::int64_t ranks_by_sort(std::span<const Cost> cost, Cost max_cost, std::span<std::int32_t> rank,
                           std::span<std::int32_t> a, std::span<std::int32_t> b,
                           std::span<std::int32_t> count) {
  const std::size_t m = cost.size();
  if (m == 0) return 0;
  constexpr std::uint64_t kMask = (std::uint64_t {1} << kBits) - 1;
  std::iota(a.begin(), a.end(), 0);
  const auto top = static_cast<std::uint64_t>(max_cost);
  for (int shift = 0; shift < 64 && (top >> shift) != 0; shift += kBits) {
    auto digit = [&](std::int32_t i) {
      const auto c = static_cast<std::uint64_t>(cost[static_cast<std::size_t>(i)]);
      return static_cast<std::size_t>((c >> shift) & kMask);
    };
    std::fill(count.begin(), count.end(), 0);
    for (std::size_t j = 0; j < m; ++j) ++count[digit(static_cast<std::int32_t>(j))];
    if (count[digit(0)] == static_cast<std::int32_t>(m)) continue;  // one digit value
    for (std::size_t d = 0, sum = 0; d < count.size(); ++d) {
      const auto c = static_cast<std::size_t>(count[d]);
      count[d] = static_cast<std::int32_t>(sum);
      sum += c;
    }
    for (std::size_t j = 0; j < m; ++j) b[static_cast<std::size_t>(count[digit(a[j])]++)] = a[j];
    std::swap(a, b);
  }
  std::int32_t r = 0;
  rank[static_cast<std::size_t>(a[0])] = 0;
  for (std::size_t j = 1; j < m; ++j) {
    const auto p = static_cast<std::size_t>(a[j]);
    if (cost[p] != cost[static_cast<std::size_t>(a[j - 1])]) ++r;
    rank[p] = r;
  }
  return std::int64_t {r} + 1;
}

// Ordinals into `rank` (input coordinates), then ranks. The table, the distinct costs and the rank
// step live in `room` (retained arrays not yet filled); a table larger than E slots gets vectors of
// its own, freed on return. Growth rebuilds the bigger table from the distinct costs in their
// order, counting every probe. Any path other than Hashed leaves `rank` partly written for the
// caller to overwrite.
RankOutcome ranks_by_hash(std::span<const Cost> cost, std::span<std::int32_t> rank, const HashRooms& room) {
  const std::size_t m = cost.size();
  RankOutcome out;
  ProbeBudget b {0, 8 * static_cast<std::int64_t>(m) + 8 * kMaxHashedCosts};
  out.budget = b.limit;
  if (m == 0) return out;  // no cost, no table
  std::size_t k = 0;       // the distinct costs so far: room.distinct[0, k)
  {
    TableStore store {room.keys, room.ords};
    CostTable t = store.make(64);
    auto stop = [&](RankPath path) {
      out.path = path;
      out.probes = b.used;
      out.table_slots = static_cast<std::int64_t>(t.key.size());
      return out;
    };
    Cost prev = 0;  // a run of equal costs reuses the previous ordinal without a lookup
    std::int32_t prev_ord = 0;
    for (std::size_t i = 0; i < m; ++i) {
      const Cost c = cost[i];
      if (c == prev) {
        rank[i] = prev_ord;
        continue;
      }
      std::size_t s = t.find(c, b);
      if (s == kNoSlot) return stop(RankPath::ExhaustedCollection);
      if (t.key[s] == 0) {
        if (static_cast<std::int64_t>(k) == kMaxHashedCosts) return stop(RankPath::DistinctLimit);
        if (4 * (k + 1) > t.key.size()) {  // keep the load <= 1/4: grow first
          t = store.make(t.key.size() * 2);
          for (std::size_t j = 0; j < k; ++j) {
            const std::size_t d = t.find(room.distinct[j], b);
            if (d == kNoSlot) return stop(RankPath::ExhaustedGrowth);
            t.key[d] = room.distinct[j];
            t.ord[d] = static_cast<std::int32_t>(j);
          }
          s = t.find(c, b);
          if (s == kNoSlot) return stop(RankPath::ExhaustedGrowth);
        }
        t.key[s] = c;
        t.ord[s] = static_cast<std::int32_t>(k);
        room.distinct[k++] = c;
      }
      rank[i] = t.ord[s];
      prev = c;
      prev_ord = t.ord[s];
    }
    out.table_slots = static_cast<std::int64_t>(t.key.size());
  }
  out.probes = b.used;
  // Ordinal -> rank without the table: a radix sort of the K distinct costs.
  const std::span<const Cost> by_ord = room.distinct.first(k);
  const std::span<std::int32_t> rank_of_ord = room.a.first(k);
  std::array<std::int32_t, std::size_t {1} << 8> count {};
  ranks_by_sort<8>(by_ord, *std::max_element(by_ord.begin(), by_ord.end()), rank_of_ord, room.b.first(k),
                   room.c.first(k), count);
  for (std::size_t i = 0; i < m; ++i) rank[i] = rank_of_ord[static_cast<std::size_t>(rank[i])];
  out.distinct = static_cast<std::int64_t>(k);
  return out;
}

// Offsets from cursors that a scatter advanced to the ends: off[k] = start(k + 1) becomes
// off[k] = start(k), off[0] = 0 (the last entry of an (N + 1)-array keeps its value).
// Counts in off[0, n) become inclusive sums (each node's end), off[n] the total: a scatter that
// decrements these cursors, walking its input backwards, leaves the starts (stable, no shift).
// The sum stays in a register (a sum through memory chains every node on the previous store).
void counts_to_ends(std::span<std::int32_t> off) noexcept {
  std::int32_t sum = 0;
  for (std::size_t k = 0; k + 1 < off.size(); ++k) {
    sum += off[k];
    off[k] = sum;
  }
  off.back() = sum;
}

// Counts in off[1, n] (node v's at v + 1) become the starts, off[0] = 0 (a running sum as above).
void counts_to_starts(std::span<std::int32_t> off) noexcept {
  std::int32_t sum = 0;
  for (std::size_t k = 1; k < off.size(); ++k) {
    sum += off[k];
    off[k] = sum;
  }
}

template <typename T>
void ends_to_starts(std::span<T> off) noexcept {
  for (std::size_t k = off.size(); k-- > 1;) off[k] = off[k - 1];
  if (!off.empty()) off[0] = 0;
}

} // namespace

Amount to_units(double x, double resolution) {
  const int k = resolution_exponent(resolution);
  if (std::isnan(x)) {
    throw std::invalid_argument("to_units: NaN has no unit value");
  }
  if (std::isinf(x)) {
    throw std::invalid_argument("to_units: an infinite amount is not finite and has no unit value");
  }
  if (x < 0.0) {
    throw std::invalid_argument("to_units: a negative amount has no unit value");
  }
  // x / q is exact (q is a power of two) unless it overflows to infinity.
  const double y = std::ldexp(x, -k);
  constexpr double kTwo52 = 4503599627370496.0;  // 2^52: doubles at or above it are integers
  constexpr double kLimit = 4611686018427387904.0;  // 2^62
  double r;
  if (y >= kTwo52) {
    r = y;
  } else {
    const double f = std::floor(y);
    r = (y - f >= 0.5) ? f + 1.0 : f;  // y - f is exact below 2^52: round half up
  }
  if (!(r < kLimit)) {
    throw std::invalid_argument(
        "to_units: the amount is at or above the bound of the units domain (2^62 units)");
  }
  return static_cast<Amount>(r);
}

double as_float(Amount x, double resolution) {
  const int k = resolution_exponent(resolution);
  return std::ldexp(static_cast<double>(x), k);
}

Graph::Graph(std::int32_t num_nodes,
             std::span<const NodeId> src,
             std::span<const NodeId> dst,
             std::span<const Amount> capacity_units,
             std::span<const Cost> cost,
             double resolution) {
  (void)resolution_exponent(resolution);
  // The unit values.
  for (std::size_t i = 0; i < capacity_units.size(); ++i) {
    const Amount c = capacity_units[i];
    if (c < 0 || c >= kDomainLimit) {
      throw std::invalid_argument("capacity_units[" + std::to_string(i) + "] = " + std::to_string(c) +
                                  " is outside the units domain [0, 2^62)");
    }
  }
  // The array shapes.
  if (num_nodes < 0) {
    throw std::invalid_argument("num_nodes must be >= 0 (got " + std::to_string(num_nodes) + ")");
  }
  const std::size_t m = src.size();
  if (dst.size() != m || capacity_units.size() != m || cost.size() != m) {
    throw std::invalid_argument("src, dst, capacity_units and cost must have the same length (got " +
                                std::to_string(src.size()) + ", " + std::to_string(dst.size()) + ", " +
                                std::to_string(capacity_units.size()) + ", " +
                                std::to_string(cost.size()) + ")");
  }
  if (m > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    throw std::invalid_argument("the number of edges exceeds INT32_MAX (array length)");
  }
  // Zero or negative costs (rule 28).
  Cost max_cost = 0;
  Cost cost_gcd = 0;
  for (std::size_t i = 0; i < m; ++i) {
    if (cost[i] < 1) {
      throw std::invalid_argument("cost[" + std::to_string(i) + "] = " + std::to_string(cost[i]) +
                                  ": costs must be positive (cost >= 1)");
    }
    max_cost = std::max(max_cost, cost[i]);
    // A running GCD of 1 stays 1 (gcd(1, c) = 1): the remaining edges are only checked.
    if (cost_gcd != 1) cost_gcd = std::gcd(cost_gcd, cost[i]);
  }
  // The remaining bounds.
  for (std::size_t i = 0; i < m; ++i) {
    if (src[i] < 0 || src[i] >= num_nodes || dst[i] < 0 || dst[i] >= num_nodes) {
      throw std::invalid_argument("edge " + std::to_string(i) + ": endpoint (" + std::to_string(src[i]) +
                                  ", " + std::to_string(dst[i]) + ") out of node range [0, " +
                                  std::to_string(num_nodes) + ")");
    }
  }
  constexpr auto kLimit = static_cast<std::uint64_t>(kDomainLimit);
  std::uint64_t cap_sum = 0;
  std::uint64_t cost_sum = 0;
  for (std::size_t i = 0; i < m; ++i) {
    cap_sum += static_cast<std::uint64_t>(capacity_units[i]);  // each < 2^62: no wrap before the check
    if (cap_sum >= kLimit) {
      throw std::invalid_argument("the sum of all capacities reaches the bound of the units domain (2^62 units)");
    }
    cost_sum += static_cast<std::uint64_t>(cost[i]);
    if (cost_sum >= kLimit) {
      throw std::invalid_argument(
          "the sum of all edge costs reaches the bound 2^62: path costs would overflow int64");
    }
  }

  uid_ = next_uid();
  resolution_ = resolution;
  num_nodes_ = num_nodes;
  max_cost_ = max_cost;
  cost_gcd_ = cost_gcd;
  const auto n = static_cast<std::size_t>(num_nodes);

  // Every retained array first (36 E + 8 (N + 1) bytes); until it is filled, each serves as
  // scratch. The cost ranks go into the copy's array (by input position); the hash table, the
  // distinct costs and the rank step borrow the others, so the hash path allocates nothing of its
  // own while the table has at most E slots (no freed block stays behind in the allocator's cache).
  capacity_.resize(m);
  cost_.resize(m);
  src_.resize(m);
  dst_.resize(m);
  input_index_.resize(m);
  in_adj_edge_index_.resize(m);
  in_col_indices_.resize(m);
  row_offsets_.assign(n + 1, 0);
  in_row_offsets_.assign(n + 1, 0);
  const std::span<std::int32_t> rank(in_col_indices_);
  const RankOutcome ro =
      ranks_by_hash(cost, rank, HashRooms {capacity_, src_, cost_, dst_, input_index_, in_adj_edge_index_});
  const std::span<std::int32_t> order(input_index_);      // forward passes 1 and 3 (final)
  const std::span<std::int32_t> mid(in_adj_edge_index_);  // forward pass 2
  const std::span<Amount> bucket(capacity_);              // the K rank cursors (K <= E)

  // Any other path: the ranks from a sort of all costs (the same ranks).
  std::int64_t k = ro.distinct;
  if (ro.path != RankPath::Hashed) {
    std::vector<std::int32_t> count(std::size_t {1} << 16);
    k = ranks_by_sort<16>(cost, max_cost, rank, input_index_, in_adj_edge_index_, count);
  }
  const auto kk = static_cast<std::size_t>(k);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  {
    auto& c = detail::thread_construction_counts();
    ++(ro.path == RankPath::Hashed ? c.hashed : c.fallback);
    if (ro.path == RankPath::DistinctLimit) ++c.distinct_limit;
    if (ro.path == RankPath::ExhaustedCollection) ++c.exhausted_collection;
    if (ro.path == RankPath::ExhaustedGrowth) ++c.exhausted_growth;
    c.last_distinct = k;
    c.last_table_slots = ro.table_slots;
    c.last_probes = ro.probes;
    c.last_budget = ro.budget;
  }
#endif
  // An edgeless graph is its validation and its two zero-filled offset arrays (APPROACH 3).
  if (m == 0) return;

  // The passes. With few nodes times few costs (N K < E), the cost rank rides with an endpoint as
  // one key: forward by dst, then by (src, rank); reverse by (dst, rank) over the ids (proposal P4).
  // Otherwise the five single-key passes: forward by dst, rank, src; reverse by rank, then dst.
  // Both give the same orders. P4's two (N K + 1)-entry histograms live in the not-yet-filled
  // capacity and cost arrays (N K + 1 <= E), so its work is O(E).
  const bool fused = n * kk < m;
  const std::span<Amount> by_dst_rank(capacity_);  // P4: (dst, rank) histogram, then cursors
  const std::span<Cost> by_src_rank(cost_);        // P4: (src, rank) histogram, then cursors
  if (fused) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
    ++detail::thread_construction_counts().fused;
#endif
    const std::size_t nk = n * kk;
    std::fill(by_dst_rank.begin(), by_dst_rank.begin() + static_cast<std::ptrdiff_t>(nk + 1), 0);
    std::fill(by_src_rank.begin(), by_src_rank.begin() + static_cast<std::ptrdiff_t>(nk + 1), 0);
    // One counter update per loop: two in one loop ran several times slower (NOTES section 1).
    for (std::size_t i = 0; i < m; ++i) ++by_dst_rank[static_cast<std::size_t>(dst[i]) * kk + static_cast<std::size_t>(rank[i]) + 1];
    for (std::size_t i = 0; i < m; ++i) ++by_src_rank[static_cast<std::size_t>(src[i]) * kk + static_cast<std::size_t>(rank[i]) + 1];
    for (std::size_t x = 0, sum = 0; x <= nk; ++x) sum = static_cast<std::size_t>(by_dst_rank[x] += static_cast<Amount>(sum));
    for (std::size_t x = 0, sum = 0; x <= nk; ++x) sum = static_cast<std::size_t>(by_src_rank[x] += static_cast<Cost>(sum));
    for (std::size_t v = 0; v < n; ++v) {
      in_row_offsets_[v] = static_cast<std::int32_t>(by_dst_rank[(v + 1) * kk]);  // ends: pass 1's cursors
      row_offsets_[v] = static_cast<std::int32_t>(by_src_rank[v * kk]);
    }
    in_row_offsets_[n] = row_offsets_[n] = static_cast<std::int32_t>(m);
  } else {
    for (std::size_t i = 0; i < m; ++i) ++in_row_offsets_[static_cast<std::size_t>(dst[i]) + 1];
    counts_to_starts(std::span<std::int32_t>(in_row_offsets_));
  }

  // Forward pass 1, by dst, from input order. The reverse offsets serve as cursors, so no pass
  // over the nodes restores them: P4 walks the input backwards from the ends to the starts; the
  // single-key passes walk it forwards from the starts to the ends, and reverse pass 2 walks back.
  if (fused) {
    for (std::size_t i = m; i-- > 0;) {
      order[static_cast<std::size_t>(--in_row_offsets_[static_cast<std::size_t>(dst[i])])] = static_cast<std::int32_t>(i);
    }
  } else {
    for (std::size_t i = 0; i < m; ++i) {
      order[static_cast<std::size_t>(in_row_offsets_[static_cast<std::size_t>(dst[i])]++)] = static_cast<std::int32_t>(i);
    }
  }

  if (fused) {
    // Forward pass 2, by (src, rank), into `mid`: the final id order; the two arrays trade places.
    for (std::size_t j = 0; j < m; ++j) {
      const auto i = static_cast<std::size_t>(order[j]);
      const std::size_t key = static_cast<std::size_t>(src[i]) * kk + static_cast<std::size_t>(rank[i]);
      mid[static_cast<std::size_t>(by_src_rank[key]++)] = static_cast<std::int32_t>(i);
    }
    input_index_.swap(in_adj_edge_index_);
    // Reverse, the ids by (dst, rank): within a key the id order gives (src, input).
    for (std::size_t e = 0; e < m; ++e) {
      const auto i = static_cast<std::size_t>(input_index_[e]);
      const std::size_t key = static_cast<std::size_t>(dst[i]) * kk + static_cast<std::size_t>(rank[i]);
      in_adj_edge_index_[static_cast<std::size_t>(by_dst_rank[key]++)] = static_cast<std::int32_t>(e);
    }
  } else {
    // Forward pass 2, by cost rank; the cursors stay for reverse pass 1.
    std::fill(bucket.begin(), bucket.begin() + static_cast<std::ptrdiff_t>(kk), 0);
    for (std::size_t i = 0; i < m; ++i) ++bucket[static_cast<std::size_t>(rank[i])];
    for (std::size_t r = 0, sum = 0; r < kk; ++r) {
      const auto c = static_cast<std::size_t>(bucket[r]);
      bucket[r] = static_cast<Amount>(sum);
      sum += c;
    }
    for (std::size_t j = 0; j < m; ++j) {
      const auto i = static_cast<std::size_t>(order[j]);
      mid[static_cast<std::size_t>(bucket[static_cast<std::size_t>(rank[i])]++)] = static_cast<std::int32_t>(i);
    }
    ends_to_starts(bucket.first(kk));

    // Forward pass 3, by src; the src histogram is the forward offsets.
    for (std::size_t i = 0; i < m; ++i) ++row_offsets_[static_cast<std::size_t>(src[i])];
    counts_to_ends(std::span<std::int32_t>(row_offsets_));
    for (std::size_t j = m; j-- > 0;) {
      const auto i = static_cast<std::size_t>(mid[j]);
      order[static_cast<std::size_t>(--row_offsets_[static_cast<std::size_t>(src[i])])] = static_cast<std::int32_t>(i);
    }

    // Reverse pass 1, by cost rank over the ids (the id order is (src, cost, dst, input)).
    const std::span<std::int32_t> by_cost(in_adj_edge_index_);  // free again after pass 3
    for (std::size_t e = 0; e < m; ++e) {
      const auto r = static_cast<std::size_t>(rank[static_cast<std::size_t>(order[e])]);
      by_cost[static_cast<std::size_t>(bucket[r]++)] = static_cast<std::int32_t>(e);
    }
  }

  // The per-edge values, permuted once (the cursors and the ranks are no longer read).
  for (std::size_t e = 0; e < m; ++e) {
    const auto i = static_cast<std::size_t>(input_index_[e]);
    capacity_[e] = capacity_units[i];
    cost_[e] = cost[i];
    src_[e] = src[i];
    dst_[e] = dst[i];
  }

  if (!fused) {
    // Reverse pass 2, by dst, into the ranks' array, walked backwards: the reverse offsets hold
    // the ends after pass 1 and return to the starts. Then the two slot arrays trade places.
    const std::span<std::int32_t> slots(in_col_indices_);
    const std::span<std::int32_t> by_cost(in_adj_edge_index_);
    for (std::size_t j = m; j-- > 0;) {
      const auto e = static_cast<std::size_t>(by_cost[j]);
      slots[static_cast<std::size_t>(--in_row_offsets_[static_cast<std::size_t>(dst_[e])])] = static_cast<std::int32_t>(e);
    }
    in_adj_edge_index_.swap(in_col_indices_);
  }
  // The source copy, over the ranks' (or the free) array.
  for (std::size_t s = 0; s < m; ++s) in_col_indices_[s] = src_[static_cast<std::size_t>(in_adj_edge_index_[s])];
}

Graph::Graph(const Graph& other)
    : uid_(next_uid()),
      resolution_(other.resolution_),
      num_nodes_(other.num_nodes_),
      max_cost_(other.max_cost_),
      cost_gcd_(other.cost_gcd_),
      capacity_(other.capacity_),
      cost_(other.cost_),
      src_(other.src_),
      dst_(other.dst_),
      input_index_(other.input_index_),
      row_offsets_(other.row_offsets_),
      in_row_offsets_(other.in_row_offsets_),
      in_adj_edge_index_(other.in_adj_edge_index_),
      in_col_indices_(other.in_col_indices_) {}

// No allocation: the source keeps empty vectors, and its offset views read the static zero.
Graph::Graph(Graph&& other) noexcept
    : uid_(std::exchange(other.uid_, 0)),
      resolution_(other.resolution_),
      num_nodes_(std::exchange(other.num_nodes_, 0)),
      max_cost_(std::exchange(other.max_cost_, 0)),
      cost_gcd_(std::exchange(other.cost_gcd_, 0)),
      capacity_(std::move(other.capacity_)),
      cost_(std::move(other.cost_)),
      src_(std::move(other.src_)),
      dst_(std::move(other.dst_)),
      input_index_(std::move(other.input_index_)),
      row_offsets_(std::move(other.row_offsets_)),
      in_row_offsets_(std::move(other.in_row_offsets_)),
      in_adj_edge_index_(std::move(other.in_adj_edge_index_)),
      in_col_indices_(std::move(other.in_col_indices_)) {}

} // namespace netgraph::core
