/* The DagCache of APPROACH 10 in building mode. */
#pragma once

#include <atomic>
#include <cstdint>
#include <map>
#include <memory>
#include <span>
#include <tuple>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// A cache of shortest-path supports bound to one graph (by uid). The key is (the roots as a
// sorted set, the direction, the Eligibility uid; no Eligibility is uid 0); the value is the full
// static support of that search (multipath = true, multi_edge = true, no residual: every tight
// edge of the statically eligible member graph over the nodes it reaches, forward-oriented for a
// reverse search) as an immutable shared PredDAG, valid after clear() and after the cache is
// destroyed. Selection is never cached: a policy projects it at creation (APPROACH 9.1, 10).
//
// Building mode: one writer, no concurrent readers. A hit returns the value; a miss runs the
// search, inserts and returns (a failed search counts one miss and inserts nothing); an invalid
// query (a graph or Eligibility of another graph, an empty root set, a root out of range) throws
// std::invalid_argument and counts nothing. clear() empties the cache; the hit and miss counters
// persist. Statistics are atomic per field, not a consistent snapshot: `bytes` is the retained
// payload of the entries (the DAGs' stored arrays), `live_values` the values alive anywhere,
// entries included. A supplied cache whose public call fails later keeps the complete entries
// and statistics of its successful lookups (APPROACH 13, limit 2).
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
  void clear() noexcept;
  [[nodiscard]] Stats stats() const noexcept;

private:
  using Key = std::tuple<std::vector<NodeId>, Direction, std::uint64_t>;
  std::uint64_t graph_uid_ {0};
  std::map<Key, std::shared_ptr<const PredDAG>> entries_ {};
  std::atomic<std::int64_t> hits_ {0};
  std::atomic<std::int64_t> misses_ {0};
  std::atomic<std::int64_t> num_entries_ {0};
  std::atomic<std::int64_t> bytes_ {0};
  // Shared with every value's deleter, so it outlives the cache while values live.
  std::shared_ptr<std::atomic<std::int64_t>> live_ {std::make_shared<std::atomic<std::int64_t>>(0)};
};

} // namespace netgraph::core
