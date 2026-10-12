/* Internal access to a DagCache: its mode (the evaluators' rule for a supplied cache) and the
   frozen statistics of one public call. */
#pragma once

#include <atomic>
#include <cstdint>
#include <memory>
#include <span>
#include <type_traits>
#include <utility>

#include "netgraph/core/dag_cache.hpp"

namespace netgraph::core::detail {

class FrozenTally;

struct DagCacheAccess {
  // Whether the cache is frozen (a frozen cache may be shared by every worker).
  [[nodiscard]] static bool frozen(const DagCache& cache) noexcept { return cache.frozen_; }
  // A lookup inside a public call (FrozenTally): a frozen hit is recorded in `tally` instead of
  // the cache, a frozen miss is counted by the call's boundary (frozen_call); a building cache
  // counts at once, as DagCache::lookup does.
  [[nodiscard]] static std::shared_ptr<const PredDAG> lookup(DagCache& cache, const Graph& g,
                                                             std::span<const NodeId> roots, Direction direction,
                                                             const Eligibility* elig, FrozenTally& tally) {
    return cache.lookup_in_call(g, roots, direction, elig, &tally);
  }
  static void add_hits(DagCache& cache, std::int64_t n) noexcept {
    cache.hits_.fetch_add(n, std::memory_order_relaxed);
  }
  static void add_miss(DagCache& cache) noexcept { cache.misses_.fetch_add(1, std::memory_order_relaxed); }
};

// Collect one public call's frozen-cache hits here; publish them only after all fallible work,
// including the Python result, succeeds. Failure discards the hits, except that DagCacheMiss
// records one miss regardless of its workers' lookups. A building cache keeps its statistics
// even on failure. Workers share this call's atomic tally; separate calls have separate
// tallies, so concurrent successful calls keep all their increments.
class FrozenTally {
public:
  explicit FrozenTally(DagCache* cache) noexcept : cache_(cache) {}
  FrozenTally(const FrozenTally&) = delete;
  FrozenTally& operator=(const FrozenTally&) = delete;
  ~FrozenTally() = default;

  [[nodiscard]] const DagCache* cache() const noexcept { return cache_; }
  void hit() noexcept { hits_.fetch_add(1, std::memory_order_relaxed); }
  // The call succeeded: its hits reach the cache (one atomic addition, cannot fail).
  void publish() noexcept {
    const auto n = hits_.exchange(0, std::memory_order_relaxed);
    if (cache_ != nullptr && n != 0) DagCacheAccess::add_hits(*cache_, n);
  }
  // The call raises DagCacheMiss: one miss; its hits are dropped.
  void miss() noexcept {
    if (cache_ != nullptr) DagCacheAccess::add_miss(*cache_);
  }

private:
  DagCache* cache_;
  std::atomic<std::int64_t> hits_ {0};
};

// Runs body(tally) as one public call on `cache` (the supplied cache, or null): publishes the
// recorded hits once body has returned (its last failing step: body builds the call's result,
// the Python object included) and returns body's value, which then moves without failing; when
// body throws, nothing is published, except one miss for a DagCacheMiss.
template <class Body>
auto frozen_call(DagCache* cache, Body&& body) -> decltype(std::forward<Body>(body)(std::declval<FrozenTally&>())) {
  using Result = decltype(std::forward<Body>(body)(std::declval<FrozenTally&>()));
  static_assert(std::is_nothrow_move_constructible_v<Result>, "the result moves out after the publication");
  FrozenTally tally(cache);
  try {
    Result result = std::forward<Body>(body)(tally);
    tally.publish();
    return result;
  } catch (const DagCacheMiss&) {
    tally.miss();
    throw;
  }
}

} // namespace netgraph::core::detail
