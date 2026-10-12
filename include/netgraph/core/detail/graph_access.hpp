/* Internal graph access for reverse searches and tests. The stored reverse-source array
   avoids fetching each source through a second array; tests also read retained array bytes.
   Everything is inline, with no Python name. Construction counters exist only in private-hooks
   builds (their accessor, the one library symbol, is defined in graph.cpp) and let C++ tests
   observe which construction branches ran. */
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>

#include "netgraph/core/detail/hooks_flag.hpp"
#include "netgraph/core/graph.hpp"

namespace netgraph::core::detail {

// The fixed hash of the distinct-cost table (splitmix64's finalizer; the table index is its low
// bits). The construction uses it; the tests build colliding costs from it.
[[nodiscard]] inline std::uint64_t cost_hash(Cost c) noexcept {
  auto x = static_cast<std::uint64_t>(c);
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}

struct GraphAccess {
  // in_col_indices[slot] = edge_src()[in_adj_edge_index()[slot]], filled before the graph is
  // published and never written after.
  [[nodiscard]] static std::span<const NodeId> in_col_indices(const Graph& g) noexcept { return g.in_col_indices_; }

  // The retained arrays' heap payload: capacity() and size() times the element size, summed.
  struct Retained {
    std::size_t capacity_bytes {0};
    std::size_t size_bytes {0};
    std::size_t arrays {0};
  };
  [[nodiscard]] static Retained retained(const Graph& g) noexcept {
    Retained r;
    auto add = [&r](const auto& v) {
      using T = typename std::remove_cvref_t<decltype(v)>::value_type;
      r.capacity_bytes += v.capacity() * sizeof(T);
      r.size_bytes += v.size() * sizeof(T);
      ++r.arrays;
    };
    add(g.capacity_);
    add(g.cost_);
    add(g.src_);
    add(g.dst_);
    add(g.input_index_);
    add(g.row_offsets_);
    add(g.in_row_offsets_);
    add(g.in_adj_edge_index_);
    add(g.in_col_indices_);
    return r;
  }
};

#if NETGRAPH_CORE_PRIVATE_HOOKS
// Test evidence (C++ tests, hooks builds only): this thread's construction counts so far, and
// what its last construction did (the last_* fields).
struct ConstructionCounts {
  std::int64_t hashed {0};                // constructions whose ranks came from the hash table
  std::int64_t fallback {0};              // constructions whose ranks came from the sort of all costs
  std::int64_t distinct_limit {0};        // fallbacks at the (L + 1)th distinct cost
  std::int64_t exhausted_collection {0};  // fallbacks: the probe budget ran out at an input lookup
  std::int64_t exhausted_growth {0};      // fallbacks: the probe budget ran out during a table growth
  std::int64_t last_distinct {0};         // K of the last construction
  std::int64_t last_table_slots {0};
  std::int64_t last_probes {0};           // probes examined by the last construction (<= last_budget)
  std::int64_t last_budget {0};           // its B = 8E + 8 x 65,536
  std::int64_t fused {0};                 // constructions taking the combined-key passes (N K < E)
};
ConstructionCounts& thread_construction_counts() noexcept;
#endif

} // namespace netgraph::core::detail
