/*
  Graph: the one graph type on integer units.

  Validate inputs in the documented order, then build edge ids and CSR arrays. N is the node
  count, E the edge count, and K the number of distinct costs. Stable counting passes group by
  a small integer key while preserving order within each group. A cost's rank is its position
  among sorted distinct costs. EdgeIds sort by (source, cost, destination, input); u's out-edges occupy
  [outgoing_offsets[u], outgoing_offsets[u + 1]); v's in-edges sort by (cost, source, input).
  An edgeless graph is its validation and two zero-filled offset arrays.
  * N K >= E: three passes (destination, cost rank, source; least significant first) give the id order, and
    two over the ids (cost rank, then destination) the reverse lists: within one (destination, cost) group the
    id order is (source, input).
  * N K < E (combined keys): a pass by destination, then one by (source, rank), give the id order, and one
    pass over the ids by (destination, rank) the reverse lists; the two (N K + 1)-entry histograms live
    in the not-yet-filled capacity and cost arrays.
  The offset arrays serve as their own cursors; a pass that advances them is followed by one
  that walks back, so no pass over the nodes restores them.

  Assign each cost its first-appearance index (ordinal) with a hash table. Resolve collisions
  by trying successive slots (linear probing). Use fixed mixing, at most 1/4 occupancy, and
  double from 64 to 2^18 slots, rebuilding in distinct-cost order. Consecutive equal costs reuse
  the previous ordinal without lookup. Radix-sort the K distinct costs, processing digits from
  least to most significant. Stop hashing at the 65,537th distinct cost, or before examining
  more than 8E + 8 x 65,536 slots. Count every examined slot, including the final slot, across
  lookups and growths. On either limit, radix-sort all costs in 16-bit digits.
  Both methods produce identical ranks.

  Allocate all retained arrays at final size, then use them as scratch before filling them.
  Store ranks in in_col_indices_. While the hash table has <= E slots, capacity_ holds its keys
  and src_ its ordinals. cost_ holds distinct costs; dst_, input_index_ and in_adj_edge_index_
  serve their sort. Counting passes use in_adj_edge_index_ for the middle forward pass, and
  capacity_/cost_ for K rank cursors or combined histograms. The fallback uses input_index_ and
  in_adj_edge_index_. Extra buffers are limited to tables > E slots (old and new coexist during
  growth) and the fallback's 2^16 digit counts. Smaller tables free no separate buffer that the
  allocator could retain.
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
#include "netgraph/core/detail/hooks_flag.hpp"
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
  std::span<std::int32_t> ordinals;  // the key's first-appearance ordinal
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
struct CostTableStorage {
  std::span<Cost> borrowed_keys;
  std::span<std::int32_t> borrowed_ordinals;
  std::vector<Cost> owned_keys {};
  std::vector<std::int32_t> owned_ordinals {};

  CostTable make(std::size_t slots) {
    if (slots <= borrowed_keys.size()) {
      const std::span<Cost> key = borrowed_keys.first(slots);
      std::fill(key.begin(), key.end(), Cost {0});
      return {key, borrowed_ordinals.first(slots), slots - 1};
    }
    owned_keys.assign(slots, 0);
    owned_ordinals.resize(slots);
    return {owned_keys, owned_ordinals, slots - 1};
  }
};

// The retained arrays the hash path borrows before they are filled (each E entries, K <= E).
struct CostRankBuffers {
  std::span<Cost> keys;             // the table's keys while it has at most E slots
  std::span<std::int32_t> ordinals;     // the table's ordinals, the same
  std::span<Cost> distinct;         // the distinct costs in first-appearance order
  std::span<std::int32_t> ordinal_to_rank, sort_positions, sort_scratch;  // the rank step: ordinal -> rank, and its two position arrays
};

enum class RankPath { Hashed, DistinctLimit, ExhaustedCollection, ExhaustedGrowth };

struct RankOutcome {
  RankPath path {RankPath::Hashed};
  std::int64_t distinct {0};
  std::int64_t table_slots {0};
  std::int64_t probes {0};
  std::int64_t budget {0};
};

// Ranks from a stable LSD radix sort of the positions of `cost` by cost, kBits bits per pass (all
// E costs and 16 bits on the fallback path, the K distinct costs and 8 bits on the hashed path),
// then one scan of the sorted positions. Digits above the largest cost's bit length are not
// visited, and a pass whose digit is the same for every cost is skipped. `a` and `b`
// (cost.size() slots each) hold the positions; `count` has 2^kBits entries. `cost` is never
// empty: ranks_by_hash returns Hashed for an edgeless graph before calling this, a successful
// collection of positive costs keeps at least one, and the fallback sort of all costs runs only
// after hashing failed, which needs an edge.
template <int kBits>
std::int64_t ranks_by_sort(std::span<const Cost> cost, Cost max_edge_cost, std::span<std::int32_t> rank,
                           std::span<std::int32_t> a, std::span<std::int32_t> b,
                           std::span<std::int32_t> count) {
  const std::size_t m = cost.size();
  constexpr std::uint64_t kMask = (std::uint64_t {1} << kBits) - 1;
  std::iota(a.begin(), a.end(), 0);
  const auto top = static_cast<std::uint64_t>(max_edge_cost);
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
RankOutcome ranks_by_hash(std::span<const Cost> cost, std::span<std::int32_t> rank, const CostRankBuffers& room) {
  const std::size_t m = cost.size();
  RankOutcome out;
  ProbeBudget b {0, 8 * static_cast<std::int64_t>(m) + 8 * kMaxHashedCosts};
  out.budget = b.limit;
  if (m == 0) return out;  // no cost, no table
  std::size_t k = 0;       // the distinct costs so far: room.distinct[0, k)
  {
    CostTableStorage store {room.keys, room.ordinals};
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
            t.ordinals[d] = static_cast<std::int32_t>(j);
          }
          s = t.find(c, b);
          if (s == kNoSlot) return stop(RankPath::ExhaustedGrowth);
        }
        t.key[s] = c;
        t.ordinals[s] = static_cast<std::int32_t>(k);
        room.distinct[k++] = c;
      }
      rank[i] = t.ordinals[s];
      prev = c;
      prev_ord = t.ordinals[s];
    }
    out.table_slots = static_cast<std::int64_t>(t.key.size());
  }
  out.probes = b.used;
  // Ordinal -> rank without the table: a radix sort of the K distinct costs.
  const std::span<const Cost> by_ord = room.distinct.first(k);
  const std::span<std::int32_t> rank_of_ord = room.ordinal_to_rank.first(k);
  std::array<std::int32_t, std::size_t {1} << 8> count {};
  ranks_by_sort<8>(by_ord, *std::max_element(by_ord.begin(), by_ord.end()), rank_of_ord, room.sort_positions.first(k),
                   room.sort_scratch.first(k), count);
  for (std::size_t i = 0; i < m; ++i) rank[i] = rank_of_ord[static_cast<std::size_t>(rank[i])];
  out.distinct = static_cast<std::int64_t>(k);
  return out;
}

// Convert counts in off[0, n) to cumulative ends; off[n] becomes the total. Scattering input
// backwards while decrementing these cursors leaves group starts and preserves order within
// each group. Forward scattering leaves ends; ends_to_starts below shifts those back to starts
// (off[k] = off[k - 1] for every k >= 1, then off[0] = 0). Keep the running sum in a register so
// each iteration does not wait for the previous store.
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
  // x / q is exact (q is a power of two) unless it overflows to infinity or falls below the
  // normal range (below 2^-1022 it may lose low bits; it then rounds to 0 units either way).
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
             std::span<const NodeId> edge_sources,
             std::span<const NodeId> edge_destinations,
             std::span<const Amount> edge_capacity_units,
             std::span<const Cost> edge_costs,
             double resolution) {
  (void)resolution_exponent(resolution);
  // The unit values.
  for (std::size_t i = 0; i < edge_capacity_units.size(); ++i) {
    const Amount c = edge_capacity_units[i];
    if (c < 0 || c >= kDomainLimit) {
      throw std::invalid_argument("edge_capacity_units[" + std::to_string(i) + "] = " + std::to_string(c) +
                                  " is outside the units domain [0, 2^62)");
    }
  }
  // The array shapes.
  if (num_nodes < 0) {
    throw std::invalid_argument("num_nodes must be >= 0 (got " + std::to_string(num_nodes) + ")");
  }
  const std::size_t m = edge_sources.size();
  if (edge_destinations.size() != m || edge_capacity_units.size() != m || edge_costs.size() != m) {
    throw std::invalid_argument("edge_sources, edge_destinations, edge_capacity_units and edge_costs must have the same length (got " +
                                std::to_string(edge_sources.size()) + ", " + std::to_string(edge_destinations.size()) + ", " +
                                std::to_string(edge_capacity_units.size()) + ", " +
                                std::to_string(edge_costs.size()) + ")");
  }
  if (m > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    throw std::invalid_argument("the number of edges exceeds INT32_MAX (array length)");
  }
  // Zero or negative costs (costs are positive).
  Cost max_edge_cost = 0;
  Cost edge_cost_gcd = 0;
  for (std::size_t i = 0; i < m; ++i) {
    if (edge_costs[i] < 1) {
      throw std::invalid_argument("edge_costs[" + std::to_string(i) + "] = " + std::to_string(edge_costs[i]) +
                                  ": costs must be positive (edge_costs >= 1)");
    }
    max_edge_cost = std::max(max_edge_cost, edge_costs[i]);
    // A running GCD of 1 stays 1 (gcd(1, c) = 1): the remaining edges are only checked.
    if (edge_cost_gcd != 1) edge_cost_gcd = std::gcd(edge_cost_gcd, edge_costs[i]);
  }
  // The remaining bounds.
  for (std::size_t i = 0; i < m; ++i) {
    if (edge_sources[i] < 0 || edge_sources[i] >= num_nodes || edge_destinations[i] < 0 || edge_destinations[i] >= num_nodes) {
      throw std::invalid_argument("edge " + std::to_string(i) + ": endpoint (" + std::to_string(edge_sources[i]) +
                                  ", " + std::to_string(edge_destinations[i]) + ") out of node range [0, " +
                                  std::to_string(num_nodes) + ")");
    }
  }
  constexpr auto kLimit = static_cast<std::uint64_t>(kDomainLimit);
  std::uint64_t cap_sum = 0;
  std::uint64_t cost_sum = 0;
  for (std::size_t i = 0; i < m; ++i) {
    cap_sum += static_cast<std::uint64_t>(edge_capacity_units[i]);  // each < 2^62: no wrap before the check
    if (cap_sum >= kLimit) {
      throw std::invalid_argument("the sum of all capacities reaches the bound of the units domain (2^62 units)");
    }
    cost_sum += static_cast<std::uint64_t>(edge_costs[i]);
    if (cost_sum >= kLimit) {
      throw std::invalid_argument(
          "the sum of all edge costs reaches the bound 2^62: path costs would overflow int64");
    }
  }

  uid_ = next_uid();
  resolution_ = resolution;
  num_nodes_ = num_nodes;
  max_edge_cost_ = max_edge_cost;
  edge_cost_gcd_ = edge_cost_gcd;
  const auto n = static_cast<std::size_t>(num_nodes);

  // Every retained array first (36 E + 8 (N + 1) bytes); until it is filled, each serves as
  // scratch. The cost ranks go into the copy's array (by input position); the hash table, the
  // distinct costs and the rank step borrow the others, so the hash path allocates nothing of its
  // own while the table has at most E slots (no freed block stays behind in the allocator's cache).
  edge_capacity_units_.resize(m);
  edge_costs_.resize(m);
  edge_sources_.resize(m);
  edge_destinations_.resize(m);
  input_positions_.resize(m);
  incoming_edges_.resize(m);
  incoming_edge_sources_.resize(m);
  outgoing_offsets_.assign(n + 1, 0);
  incoming_offsets_.assign(n + 1, 0);
  const std::span<std::int32_t> rank(incoming_edge_sources_);
  const RankOutcome ro =
      ranks_by_hash(edge_costs, rank, CostRankBuffers {edge_capacity_units_, edge_sources_, edge_costs_, edge_destinations_, input_positions_, incoming_edges_});
  const std::span<std::int32_t> order(input_positions_);      // forward passes 1 and 3 (final)
  const std::span<std::int32_t> mid(incoming_edges_);  // forward pass 2
  const std::span<Amount> bucket(edge_capacity_units_);              // the K rank cursors (K <= E)

  // Any other path: the ranks from a sort of all costs (the same ranks).
  std::int64_t k = ro.distinct;
  if (ro.path != RankPath::Hashed) {
    std::vector<std::int32_t> count(std::size_t {1} << 16);
    k = ranks_by_sort<16>(edge_costs, max_edge_cost, rank, input_positions_, incoming_edges_, count);
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
  // An edgeless graph is its validation and its two zero-filled offset arrays.
  if (m == 0) return;

  // The passes. With few nodes times few costs (N K < E), the cost rank rides with an endpoint as
  // one key: forward by destination, then by (source, rank); reverse by (destination, rank) over the ids.
  // Otherwise the five single-key passes: forward by destination, rank, source; reverse by rank, then destination.
  // Both give the same orders. The combined-key passes keep their two (N K + 1)-entry histograms
  // in the not-yet-filled capacity and cost arrays (N K + 1 <= E), so their work is O(E).
  const bool fused = n * kk < m;
  const std::span<Amount> by_dst_rank(edge_capacity_units_);  // fused: (destination, rank) histogram, then cursors
  const std::span<Cost> by_src_rank(edge_costs_);        // fused: (source, rank) histogram, then cursors
  if (fused) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
    ++detail::thread_construction_counts().fused;
#endif
    const std::size_t nk = n * kk;
    std::fill(by_dst_rank.begin(), by_dst_rank.begin() + static_cast<std::ptrdiff_t>(nk + 1), 0);
    std::fill(by_src_rank.begin(), by_src_rank.begin() + static_cast<std::ptrdiff_t>(nk + 1), 0);
    // Keep the counter updates in separate loops for speed.
    for (std::size_t i = 0; i < m; ++i) ++by_dst_rank[static_cast<std::size_t>(edge_destinations[i]) * kk + static_cast<std::size_t>(rank[i]) + 1];
    for (std::size_t i = 0; i < m; ++i) ++by_src_rank[static_cast<std::size_t>(edge_sources[i]) * kk + static_cast<std::size_t>(rank[i]) + 1];
    for (std::size_t x = 0, sum = 0; x <= nk; ++x) sum = static_cast<std::size_t>(by_dst_rank[x] += static_cast<Amount>(sum));
    for (std::size_t x = 0, sum = 0; x <= nk; ++x) sum = static_cast<std::size_t>(by_src_rank[x] += static_cast<Cost>(sum));
    for (std::size_t v = 0; v < n; ++v) {
      incoming_offsets_[v] = static_cast<std::int32_t>(by_dst_rank[(v + 1) * kk]);  // ends: pass 1's cursors
      outgoing_offsets_[v] = static_cast<std::int32_t>(by_src_rank[v * kk]);
    }
    incoming_offsets_[n] = outgoing_offsets_[n] = static_cast<std::int32_t>(m);
  } else {
    for (std::size_t i = 0; i < m; ++i) ++incoming_offsets_[static_cast<std::size_t>(edge_destinations[i]) + 1];
    counts_to_starts(std::span<std::int32_t>(incoming_offsets_));
  }

  // Forward pass 1, by destination, from input order. The reverse offsets serve as cursors, so no pass
  // over the nodes restores them: the combined-key passes walk the input backwards from the ends
  // to the starts; the single-key passes walk it forwards from the starts to the ends, and
  // reverse pass 2 walks back.
  if (fused) {
    for (std::size_t i = m; i-- > 0;) {
      order[static_cast<std::size_t>(--incoming_offsets_[static_cast<std::size_t>(edge_destinations[i])])] = static_cast<std::int32_t>(i);
    }
  } else {
    for (std::size_t i = 0; i < m; ++i) {
      order[static_cast<std::size_t>(incoming_offsets_[static_cast<std::size_t>(edge_destinations[i])]++)] = static_cast<std::int32_t>(i);
    }
  }

  if (fused) {
    // Forward pass 2, by (source, rank), into `mid`: the final id order; the two arrays trade places.
    for (std::size_t j = 0; j < m; ++j) {
      const auto i = static_cast<std::size_t>(order[j]);
      const std::size_t key = static_cast<std::size_t>(edge_sources[i]) * kk + static_cast<std::size_t>(rank[i]);
      mid[static_cast<std::size_t>(by_src_rank[key]++)] = static_cast<std::int32_t>(i);
    }
    input_positions_.swap(incoming_edges_);
    // Reverse, the ids by (destination, rank): within a key the id order gives (source, input).
    for (std::size_t e = 0; e < m; ++e) {
      const auto i = static_cast<std::size_t>(input_positions_[e]);
      const std::size_t key = static_cast<std::size_t>(edge_destinations[i]) * kk + static_cast<std::size_t>(rank[i]);
      incoming_edges_[static_cast<std::size_t>(by_dst_rank[key]++)] = static_cast<std::int32_t>(e);
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

    // Forward pass 3, by source; the source histogram is the forward offsets.
    for (std::size_t i = 0; i < m; ++i) ++outgoing_offsets_[static_cast<std::size_t>(edge_sources[i])];
    counts_to_ends(std::span<std::int32_t>(outgoing_offsets_));
    for (std::size_t j = m; j-- > 0;) {
      const auto i = static_cast<std::size_t>(mid[j]);
      order[static_cast<std::size_t>(--outgoing_offsets_[static_cast<std::size_t>(edge_sources[i])])] = static_cast<std::int32_t>(i);
    }

    // Reverse pass 1, by cost rank over the ids (the id order is (source, cost, destination, input)).
    const std::span<std::int32_t> by_cost(incoming_edges_);  // free again after pass 3
    for (std::size_t e = 0; e < m; ++e) {
      const auto r = static_cast<std::size_t>(rank[static_cast<std::size_t>(order[e])]);
      by_cost[static_cast<std::size_t>(bucket[r]++)] = static_cast<std::int32_t>(e);
    }
  }

  // The per-edge values, permuted once (the cursors and the ranks are no longer read).
  for (std::size_t e = 0; e < m; ++e) {
    const auto i = static_cast<std::size_t>(input_positions_[e]);
    edge_capacity_units_[e] = edge_capacity_units[i];
    edge_costs_[e] = edge_costs[i];
    edge_sources_[e] = edge_sources[i];
    edge_destinations_[e] = edge_destinations[i];
  }

  if (!fused) {
    // Reverse pass 2, by destination, into the ranks' array, walked backwards: the reverse offsets hold
    // the ends after pass 1 and return to the starts. Then the two slot arrays trade places.
    const std::span<std::int32_t> slots(incoming_edge_sources_);
    const std::span<std::int32_t> by_cost(incoming_edges_);
    for (std::size_t j = m; j-- > 0;) {
      const auto e = static_cast<std::size_t>(by_cost[j]);
      slots[static_cast<std::size_t>(--incoming_offsets_[static_cast<std::size_t>(edge_destinations_[e])])] = static_cast<std::int32_t>(e);
    }
    incoming_edges_.swap(incoming_edge_sources_);
  }
  // The source copy, over the ranks' (or the free) array.
  for (std::size_t s = 0; s < m; ++s) incoming_edge_sources_[s] = edge_sources_[static_cast<std::size_t>(incoming_edges_[s])];
}

Graph::Graph(const Graph& other)
    : uid_(next_uid()),
      resolution_(other.resolution_),
      num_nodes_(other.num_nodes_),
      max_edge_cost_(other.max_edge_cost_),
      edge_cost_gcd_(other.edge_cost_gcd_),
      edge_capacity_units_(other.edge_capacity_units_),
      edge_costs_(other.edge_costs_),
      edge_sources_(other.edge_sources_),
      edge_destinations_(other.edge_destinations_),
      input_positions_(other.input_positions_),
      outgoing_offsets_(other.outgoing_offsets_),
      incoming_offsets_(other.incoming_offsets_),
      incoming_edges_(other.incoming_edges_),
      incoming_edge_sources_(other.incoming_edge_sources_) {}

// No allocation: the source keeps empty vectors, and its offset views read the static zero.
Graph::Graph(Graph&& other) noexcept
    : uid_(std::exchange(other.uid_, 0)),
      resolution_(other.resolution_),
      num_nodes_(std::exchange(other.num_nodes_, 0)),
      max_edge_cost_(std::exchange(other.max_edge_cost_, 0)),
      edge_cost_gcd_(std::exchange(other.edge_cost_gcd_, 0)),
      edge_capacity_units_(std::move(other.edge_capacity_units_)),
      edge_costs_(std::move(other.edge_costs_)),
      edge_sources_(std::move(other.edge_sources_)),
      edge_destinations_(std::move(other.edge_destinations_)),
      input_positions_(std::move(other.input_positions_)),
      outgoing_offsets_(std::move(other.outgoing_offsets_)),
      incoming_offsets_(std::move(other.incoming_offsets_)),
      incoming_edges_(std::move(other.incoming_edges_)),
      incoming_edge_sources_(std::move(other.incoming_edge_sources_)) {}

} // namespace netgraph::core
