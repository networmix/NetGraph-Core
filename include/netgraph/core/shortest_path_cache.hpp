/* The ShortestPathCache: building and frozen modes, the dedicated frozen-miss exception. */
#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <map>
#include <memory>
#include <span>
#include <string>
#include <tuple>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

namespace detail {
struct DagCacheAccess;
class FrozenCacheHitTally;
class CacheGate;
} // namespace detail

// A frozen-cache lookup with no matching key throws ShortestPathCacheMiss. The shared message can be
// copied without allocation. Keep what() out of line: this gives the exception one vtable and
// type_info, so the binding can reliably catch it.
class ShortestPathCacheMiss : public std::exception {
public:
  explicit ShortestPathCacheMiss(const std::string& what) : what_(std::make_shared<const std::string>(what)) {}
  [[nodiscard]] const char* what() const noexcept override;

private:
  std::shared_ptr<const std::string> what_;
};

// A cache of static shortest-path supports bound to one graph uid. A support contains every
// statically eligible edge, whatever its capacity, on a shortest path from the roots (Forward)
// or to them (Reverse), ignoring residual capacity. It keeps every equal-cost path and every
// parallel edge. The key is (sorted root set, direction, Eligibility uid; 0 for no Eligibility).
// Values are immutable shared PredecessorDAGs; reverse searches keep the original edge orientation.
// Values remain valid after clear() or cache destruction. The route choice is applied when a policy
// creates a flow and is never cached.
//
// Building mode (a new cache): one exclusive lease per call. A hit returns the value; a miss
// runs the search, inserts and returns (a failed search counts one miss and inserts nothing); an
// invalid query (a graph or Eligibility of another graph, an empty root set, a root out of range,
// a direction outside its enumeration) throws std::invalid_argument and counts nothing. clear()
// empties the cache; the hit and miss counters persist.
//
// Frozen mode (after freeze()): read-only and shareable; lookups from several threads at once are
// safe. A hit returns the value; a miss counts one miss and throws ShortestPathCacheMiss (nothing is
// inserted); clear() throws std::invalid_argument; freeze() again has no effect; thaw() returns to
// building mode after outstanding calls release their leases. Lifecycle calls take exclusive
// admission, in FIFO order with lookups, so later building calls and frozen readers cannot
// starve them. Callers keep the cache alive throughout every call; returned values own their data.
//
// Each statistic is atomic, but reading several fields need not give one consistent snapshot.
// `bytes` counts the entries' stored DAG arrays; `live_values` counts values alive anywhere,
// including in the cache. A building cache keeps completed entries and lookup statistics even
// if the enclosing public call fails. Entries are always complete; keeping them changes no
// later result. This is the exception to public calls' all-or-nothing state guarantee.
// For a frozen cache, hits count only after the whole public call succeeds, including its
// Python result. A failed call counts nothing, except one miss if it raises ShortestPathCacheMiss.
// Direct lookup() follows the same rule: one hit on success, one miss on ShortestPathCacheMiss.
class ShortestPathCache {
public:
  struct Stats {
    std::int64_t hits {0};
    std::int64_t misses {0};
    std::int64_t entries {0};
    std::int64_t bytes {0};
    std::int64_t live_values {0};
  };

  explicit ShortestPathCache(const Graph& graph);
  ShortestPathCache(const ShortestPathCache&) = delete;
  ShortestPathCache& operator=(const ShortestPathCache&) = delete;
  ~ShortestPathCache();

  [[nodiscard]] std::uint64_t graph_uid() const noexcept { return graph_uid_; }
  [[nodiscard]] std::shared_ptr<const PredecessorDAG> lookup(const Graph& graph, std::span<const NodeId> roots,
                                                      Direction direction, const Eligibility* eligibility);
  void freeze();
  void thaw();
  void clear();
  [[nodiscard]] Stats stats() const noexcept;

private:
  friend struct detail::DagCacheAccess;
  // Caller owns the lifecycle lease, or this cache is exclusive to one worker.
  void clear_contents() noexcept;
  // lookup() inside a public call: with a tally of this cache, a frozen hit is recorded in it and a
  // frozen miss counts nothing here (the call's boundary counts it).
  std::shared_ptr<const PredecessorDAG> lookup_in_call(const Graph& graph, std::span<const NodeId> roots, Direction direction,
                                                const Eligibility* eligibility, detail::FrozenCacheHitTally* tally);
  std::unique_ptr<detail::CacheGate> gate_;
  using Key = std::tuple<std::vector<NodeId>, Direction, std::uint64_t>;
  // A lookup's key over borrowed roots (a validated sorted, duplicate-free span): a hit needs no copy.
  using KeyView = std::tuple<std::span<const NodeId>, Direction, std::uint64_t>;
  // The order of Key (lexicographic roots, then direction, then Eligibility uid), also between a Key
  // and a KeyView.
  struct KeyLess {
    using is_transparent = void;
    static bool less(std::span<const NodeId> ra, Direction da, std::uint64_t ua, std::span<const NodeId> rb,
                     Direction db, std::uint64_t ub) noexcept {
      const std::size_t n = std::min(ra.size(), rb.size());
      for (std::size_t i = 0; i < n; ++i) {
        if (ra[i] != rb[i]) return ra[i] < rb[i];
      }
      if (ra.size() != rb.size()) return ra.size() < rb.size();
      return da != db ? da < db : ua < ub;
    }
    template <class A, class B>
    bool operator()(const A& a, const B& b) const noexcept {
      return less(std::get<0>(a), std::get<1>(a), std::get<2>(a), std::get<0>(b), std::get<1>(b), std::get<2>(b));
    }
  };
  std::uint64_t graph_uid_ {0};
  std::map<Key, std::shared_ptr<const PredecessorDAG>, KeyLess> entries_ {};
  std::atomic<std::int64_t> hits_ {0};
  std::atomic<std::int64_t> misses_ {0};
  std::atomic<std::int64_t> num_entries_ {0};
  std::atomic<std::int64_t> bytes_ {0};
  // Shared with every value's deleter, so it outlives the cache while values live.
  std::shared_ptr<std::atomic<std::int64_t>> live_ {std::make_shared<std::atomic<std::int64_t>>(0)};
};

} // namespace netgraph::core
