/* Internal access to a ShortestPathCache: its mode (the evaluators' rule for a supplied cache) and the
   frozen statistics of one public call. */
#pragma once

#include <atomic>
#include <cstdint>
#include <memory>
#include <span>
#include <type_traits>
#include <utility>

#include "netgraph/core/shortest_path_cache.hpp"

namespace netgraph::core::detail {

class FrozenCacheHitTally;

struct DagCacheAccess {
  // Whether the cache is frozen (a frozen cache may be shared by every worker).
  [[nodiscard]] static bool frozen(const ShortestPathCache& cache) noexcept { return cache.frozen_; }
  // A lookup inside a public call (FrozenCacheHitTally): a frozen hit is recorded in `tally` instead of
  // the cache, a frozen miss is counted by the call's boundary (frozen_call); a building cache
  // counts at once, as ShortestPathCache::lookup does.
  [[nodiscard]] static std::shared_ptr<const PredecessorDAG> lookup(ShortestPathCache& cache, const Graph& graph,
                                                             std::span<const NodeId> roots, Direction direction,
                                                             const Eligibility* eligibility, FrozenCacheHitTally& tally) {
    return cache.lookup_in_call(graph, roots, direction, eligibility, &tally);
  }
  static void add_hits(ShortestPathCache& cache, std::int64_t n) noexcept {
    cache.hits_.fetch_add(n, std::memory_order_relaxed);
  }
  static void add_miss(ShortestPathCache& cache) noexcept { cache.misses_.fetch_add(1, std::memory_order_relaxed); }
};

// Collect one public call's frozen-cache hits here; publish them only after all fallible work,
// including the Python result, succeeds. Failure discards the hits, except that ShortestPathCacheMiss
// records one miss regardless of its workers' lookups. A building cache keeps its statistics
// even on failure. Workers share this call's atomic tally; separate calls have separate
// tallies, so concurrent successful calls keep all their increments.
class FrozenCacheHitTally {
public:
  explicit FrozenCacheHitTally(ShortestPathCache* cache) noexcept : cache_(cache) {}
  FrozenCacheHitTally(const FrozenCacheHitTally&) = delete;
  FrozenCacheHitTally& operator=(const FrozenCacheHitTally&) = delete;
  ~FrozenCacheHitTally() = default;

  [[nodiscard]] const ShortestPathCache* cache() const noexcept { return cache_; }
  void hit() noexcept { hits_.fetch_add(1, std::memory_order_relaxed); }
  // The call succeeded: its hits reach the cache (one atomic addition, cannot fail).
  void publish() noexcept {
    const auto n = hits_.exchange(0, std::memory_order_relaxed);
    if (cache_ != nullptr && n != 0) DagCacheAccess::add_hits(*cache_, n);
  }
  // The call raises ShortestPathCacheMiss: one miss; its hits are dropped.
  void miss() noexcept {
    if (cache_ != nullptr) DagCacheAccess::add_miss(*cache_);
  }

private:
  ShortestPathCache* cache_;
  std::atomic<std::int64_t> hits_ {0};
};

// Runs body(tally) as one public call on `cache` (the supplied cache, or null): publishes the
// recorded hits once body has returned (its last failing step: body builds the call's result,
// the Python object included) and returns body's value, which then moves without failing; when
// body throws, nothing is published, except one miss for a ShortestPathCacheMiss.
template <class Body>
auto frozen_call(ShortestPathCache* cache, Body&& body) -> decltype(std::forward<Body>(body)(std::declval<FrozenCacheHitTally&>())) {
  using Result = decltype(std::forward<Body>(body)(std::declval<FrozenCacheHitTally&>()));
  static_assert(std::is_nothrow_move_constructible_v<Result>, "the result moves out after the publication");
  FrozenCacheHitTally tally(cache);
  try {
    Result result = std::forward<Body>(body)(tally);
    tally.publish();
    return result;
  } catch (const ShortestPathCacheMiss&) {
    tally.miss();
    throw;
  }
}

} // namespace netgraph::core::detail
