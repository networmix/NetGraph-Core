/*
  The DagCache: a map from (roots, direction, Eligibility uid) to the full static
  support of that search, with atomic statistics; building mode (a miss searches and inserts) and
  frozen mode (read-only: concurrent hits read the map and bump atomic counters, a miss throws).
*/
#include "netgraph/core/dag_cache.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

#include "netgraph/core/detail/cache_access.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/traversal.hpp"

namespace netgraph::core {

namespace {

// Decrements the live-value counter when the last holder of a value lets it go.
struct LiveDeleter {
  std::shared_ptr<std::atomic<std::int64_t>> live;
  void operator()(const PredDAG* dag) const noexcept {
    live->fetch_sub(1, std::memory_order_relaxed);
    delete dag;
  }
};

} // namespace

// The key function of DagCacheMiss (dag_cache.hpp): its vtable and type_info live here only.
const char* DagCacheMiss::what() const noexcept { return what_->c_str(); }

DagCache::DagCache(const Graph& graph) : graph_uid_(graph.uid()) {}

std::shared_ptr<const PredDAG> DagCache::lookup(const Graph& g, std::span<const NodeId> roots,
                                                Direction direction, const Eligibility* elig) {
  return lookup_in_call(g, roots, direction, elig, nullptr);
}

std::shared_ptr<const PredDAG> DagCache::lookup_in_call(const Graph& g, std::span<const NodeId> roots,
                                                        Direction direction, const Eligibility* elig,
                                                        detail::FrozenTally* tally) {
  // Inside a public call a frozen cache's statistics wait for the call's boundary.
  const bool staged = frozen_ && tally != nullptr && tally->cache() == this;
  // Invalid queries throw before anything is counted.
  if (g.uid() != graph_uid_) throw std::invalid_argument("DagCache: the graph is not the cache's graph (foreign graph)");
  if (elig != nullptr && elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  // The roots as a sorted set: borrowed when the caller's span already is one (a hit then copies
  // nothing), else normalised (node_set's checks either way, in its order).
  std::vector<NodeId> normalised;
  std::span<const NodeId> key_roots = roots;
  if (roots.empty() || !std::is_sorted(roots.begin(), roots.end()) ||
      std::adjacent_find(roots.begin(), roots.end()) != roots.end()) {
    normalised = detail::node_set(roots, g.num_nodes(), "roots");
    key_roots = normalised;
  } else {
    detail::check_range(roots, g.num_nodes(), "roots");
  }
  detail::check_member(direction, "direction");
  const std::uint64_t elig_uid = elig != nullptr ? elig->uid() : 0;
  const auto it = entries_.find(KeyView {key_roots, direction, elig_uid});
  if (it != entries_.end()) {
    if (staged) {
      tally->hit();
    } else {
      hits_.fetch_add(1, std::memory_order_relaxed);
    }
    return it->second;
  }
  if (frozen_) {
    // A frozen miss counts one miss (the lookup found no entry; inside a call, the boundary
    // counts it) and inserts nothing; the exception is built first, so a failed allocation there
    // counts nothing.
    std::string roots_text;
    for (NodeId r : key_roots) roots_text += (roots_text.empty() ? "" : ", ") + std::to_string(r);
    const DagCacheMiss miss("DagCache: the frozen cache has no entry for roots {" + roots_text + "}, direction " +
                            (direction == Direction::Forward ? "Forward" : "Reverse") + ", Eligibility uid " +
                            std::to_string(elig_uid) + " (a frozen cache never searches)");
    if (!staged) misses_.fetch_add(1, std::memory_order_relaxed);
    throw miss;
  }
  // The entry's owned key first (the borrowed roots copied): a failed copy counts nothing, like a
  // failed normalisation. Then a failed search (or a failed insertion) inserts nothing and the miss
  // stays counted.
  Key key {normalised.empty() ? std::vector<NodeId>(key_roots.begin(), key_roots.end()) : std::move(normalised),
           direction, elig_uid};
  misses_.fetch_add(1, std::memory_order_relaxed);
  auto value = std::make_unique<PredDAG>(detail::static_support(g, std::get<0>(key), direction, elig));
  const auto bytes = static_cast<std::int64_t>(detail::PredDAGBuilder::payload_bytes(*value));
  live_->fetch_add(1, std::memory_order_relaxed);  // the deleter runs even if the control block fails
  std::shared_ptr<const PredDAG> shared(value.release(), LiveDeleter {live_});
  entries_.emplace(std::move(key), shared);  // strong guarantee: a throw leaves the map unchanged
  num_entries_.fetch_add(1, std::memory_order_relaxed);
  bytes_.fetch_add(bytes, std::memory_order_relaxed);
  return shared;
}

void DagCache::freeze() noexcept { frozen_ = true; }

void DagCache::thaw() noexcept { frozen_ = false; }

void DagCache::clear() {
  if (frozen_) {
    throw std::invalid_argument(
        "DagCache.clear is rejected while the cache is frozen (join every reader, then thaw() and clear())");
  }
  entries_.clear();
  num_entries_.store(0, std::memory_order_relaxed);
  bytes_.store(0, std::memory_order_relaxed);
}

DagCache::Stats DagCache::stats() const noexcept {
  Stats s;
  s.hits = hits_.load(std::memory_order_relaxed);
  s.misses = misses_.load(std::memory_order_relaxed);
  s.entries = num_entries_.load(std::memory_order_relaxed);
  s.bytes = bytes_.load(std::memory_order_relaxed);
  s.live_values = live_->load(std::memory_order_relaxed);
  return s;
}

} // namespace netgraph::core
