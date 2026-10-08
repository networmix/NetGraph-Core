/* The DagCache of APPROACH 10: building and frozen modes, the dedicated frozen-miss exception. */
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
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

namespace detail {
struct DagCacheAccess;
class FrozenTally;
} // namespace detail

// A lookup in a frozen DagCache whose key is absent (APPROACH 10): a dedicated exception, not a
// std::runtime_error, so no evaluator converts it into a status (a sensitivity entry's
// `unevaluable` comes from a std::runtime_error); every evaluator propagates it unchanged. Its
// message is shared, so copying it never allocates. what() is defined in dag_cache.cpp: as the
// key function it gives the type one vtable and one type_info. An all-inline class gets a weak
// copy per translation unit, and the LTO build of the extension kept two, so the binding's catch
// of DagCacheMiss missed the thrown one (NOTES R26).
class DagCacheMiss : public std::exception {
public:
  explicit DagCacheMiss(const std::string& what) : what_(std::make_shared<const std::string>(what)) {}
  [[nodiscard]] const char* what() const noexcept override;

private:
  std::shared_ptr<const std::string> what_;
};

// A cache of shortest-path supports bound to one graph (by uid). The key is (the roots as a
// sorted set, the direction, the Eligibility uid; no Eligibility is uid 0); the value is the full
// static support of that search (multipath = true, multi_edge = true, no residual: every tight
// edge of the statically eligible member graph over the nodes it reaches, forward-oriented for a
// reverse search) as an immutable shared PredDAG, valid after clear() and after the cache is
// destroyed. Selection is never cached: a policy projects it at creation (APPROACH 9.1, 10).
//
// Building mode (a new cache): one writer, no concurrent readers. A hit returns the value; a miss
// runs the search, inserts and returns (a failed search counts one miss and inserts nothing); an
// invalid query (a graph or Eligibility of another graph, an empty root set, a root out of range,
// a direction outside its enumeration) throws std::invalid_argument and counts nothing. clear()
// empties the cache; the hit and miss counters persist.
//
// Frozen mode (after freeze()): read-only and shareable; lookups from several threads at once are
// safe. A hit returns the value; a miss counts one miss and throws DagCacheMiss (nothing is
// inserted); clear() throws std::invalid_argument; freeze() again has no effect; thaw() returns to
// building mode once the caller has joined every reader (no lookup may run concurrently with
// freeze(), thaw() or destruction).
//
// Statistics are atomic per field, not a consistent snapshot: `bytes` is the retained payload of
// the entries (the DAGs' stored arrays), `live_values` the values alive anywhere, entries
// included. A supplied building cache whose public call fails later keeps the complete entries and
// statistics of its successful lookups (APPROACH 13, limit 2). A frozen cache's statistics follow
// the public call that reads it (a policy constructor, an evaluator; NOTES R27): its hits are
// counted when the call succeeds, its Python result included, and a failed call changes no
// statistic, except that a call raising DagCacheMiss counts one miss. A direct lookup() is such a
// call: a hit counts one hit, a miss one miss.
class DagCache {
public:
  struct Stats {
    std::int64_t hits {0};
    std::int64_t misses {0};
    std::int64_t entries {0};
    std::int64_t bytes {0};
    std::int64_t live_values {0};
  };

  explicit DagCache(const Graph& graph);
  DagCache(const DagCache&) = delete;
  DagCache& operator=(const DagCache&) = delete;
  ~DagCache() = default;

  [[nodiscard]] std::uint64_t graph_uid() const noexcept { return graph_uid_; }
  [[nodiscard]] std::shared_ptr<const PredDAG> lookup(const Graph& g, std::span<const NodeId> roots,
                                                      Direction direction, const Eligibility* elig);
  void freeze() noexcept;
  void thaw() noexcept;
  void clear();
  [[nodiscard]] Stats stats() const noexcept;

private:
  friend struct detail::DagCacheAccess;
  // lookup() inside a public call: with a tally of this cache, a frozen hit is recorded in it and a
  // frozen miss counts nothing here (the call's boundary counts it).
  std::shared_ptr<const PredDAG> lookup_in_call(const Graph& g, std::span<const NodeId> roots, Direction direction,
                                                const Eligibility* elig, detail::FrozenTally* tally);
  bool frozen_ {false};
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
  std::map<Key, std::shared_ptr<const PredDAG>, KeyLess> entries_ {};
  std::atomic<std::int64_t> hits_ {0};
  std::atomic<std::int64_t> misses_ {0};
  std::atomic<std::int64_t> num_entries_ {0};
  std::atomic<std::int64_t> bytes_ {0};
  // Shared with every value's deleter, so it outlives the cache while values live.
  std::shared_ptr<std::atomic<std::int64_t>> live_ {std::make_shared<std::atomic<std::int64_t>>(0)};
};

} // namespace netgraph::core
