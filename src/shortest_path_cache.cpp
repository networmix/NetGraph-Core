/*
  The ShortestPathCache: a map from (roots, direction, Eligibility uid) to the full static
  support of that search, with atomic statistics; building mode (a miss searches and inserts) and
  frozen mode (read-only: concurrent hits read the map and bump atomic counters, a miss throws).
*/
#include "netgraph/core/shortest_path_cache.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <utility>

#include "netgraph/core/detail/cache_access.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/route_builders.hpp"
#include "netgraph/core/detail/search_detail.hpp"

namespace netgraph::core {

namespace {

// Decrements the live-value counter when the last holder of a value lets it go.
struct LiveDeleter {
  std::shared_ptr<std::atomic<std::int64_t>> live;
  void operator()(const PredecessorDAG* dag) const noexcept {
    live->fetch_sub(1, std::memory_order_relaxed);
    delete dag;
  }
};

} // namespace

// The key function of ShortestPathCacheMiss (shortest_path_cache.hpp): its vtable and type_info live here only.
const char* ShortestPathCacheMiss::what() const noexcept { return what_->c_str(); }

ShortestPathCache::ShortestPathCache(const Graph& graph)
    : gate_(std::make_unique<detail::CacheGate>()), graph_uid_(graph.uid()) {}

ShortestPathCache::~ShortestPathCache() = default;

std::shared_ptr<const PredecessorDAG> ShortestPathCache::lookup(const Graph& graph, std::span<const NodeId> roots,
                                                Direction direction, const Eligibility* eligibility) {
  if (graph.uid() != graph_uid_) throw std::invalid_argument("ShortestPathCache: the graph is not the cache's graph (foreign graph)");
  if (eligibility != nullptr && eligibility->graph_uid() != graph.uid()) throw std::invalid_argument("foreign Eligibility: built for another graph");
  // Copy borrowed roots before waiting. Small direct hits need no heap allocation.
  std::array<NodeId, 8> small;
  std::vector<NodeId> large;
  std::span<NodeId> owned;
  if (roots.size() <= small.size()) {
    owned = std::span(small).first(roots.size());
    std::copy(roots.begin(), roots.end(), owned.begin());
  } else {
    large.assign(roots.begin(), roots.end());
    owned = large;
  }
  if (owned.empty()) (void)detail::node_set(owned, graph.num_nodes(), "roots");
  detail::check_range(owned, graph.num_nodes(), "roots");
  std::sort(owned.begin(), owned.end());
  const auto count = static_cast<std::size_t>(std::unique(owned.begin(), owned.end()) - owned.begin());
  detail::check_member(direction, "direction");
  const detail::CacheGate::Lease lease(*gate_);
  return lookup_in_call(graph, owned.first(count), direction, eligibility, nullptr);
}

std::shared_ptr<const PredecessorDAG> ShortestPathCache::lookup_in_call(const Graph& graph, std::span<const NodeId> roots,
                                                        Direction direction, const Eligibility* eligibility,
                                                        detail::FrozenCacheHitTally* tally) {
  // Inside a public call a frozen cache's statistics wait for the call's boundary.
  const bool staged = gate_->frozen() && tally != nullptr && tally->cache() == this;
  // Invalid queries throw before anything is counted.
  if (graph.uid() != graph_uid_) throw std::invalid_argument("ShortestPathCache: the graph is not the cache's graph (foreign graph)");
  if (eligibility != nullptr && eligibility->graph_uid() != graph.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  // The roots as a sorted set: borrowed when the caller's span already is one (a hit then copies
  // nothing), else normalised (node_set's checks either way, in its order).
  std::vector<NodeId> normalised;
  std::span<const NodeId> key_roots = roots;
  if (roots.empty() || !std::is_sorted(roots.begin(), roots.end()) ||
      std::adjacent_find(roots.begin(), roots.end()) != roots.end()) {
    normalised = detail::node_set(roots, graph.num_nodes(), "roots");
    key_roots = normalised;
  } else {
    detail::check_range(roots, graph.num_nodes(), "roots");
  }
  detail::check_member(direction, "direction");
  const std::uint64_t elig_uid = eligibility != nullptr ? eligibility->uid() : 0;
  const auto it = entries_.find(KeyView {key_roots, direction, elig_uid});
  if (it != entries_.end()) {
    if (staged) {
      tally->hit();
    } else {
      hits_.fetch_add(1, std::memory_order_relaxed);
    }
    return it->second;
  }
  if (gate_->frozen()) {
    // A frozen miss counts one miss (the lookup found no entry; inside a call, the boundary
    // counts it) and inserts nothing; the exception is built first, so a failed allocation there
    // counts nothing.
    std::string roots_text;
    for (NodeId r : key_roots) roots_text += (roots_text.empty() ? "" : ", ") + std::to_string(r);
    const ShortestPathCacheMiss miss("ShortestPathCache: the frozen cache has no entry for roots {" + roots_text + "}, direction " +
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
  auto value = std::make_unique<PredecessorDAG>(detail::static_support(graph, std::get<0>(key), direction, eligibility));
  const auto bytes = static_cast<std::int64_t>(detail::PredecessorDAGBuilder::payload_bytes(*value));
  live_->fetch_add(1, std::memory_order_relaxed);  // the deleter runs even if the control block fails
  std::shared_ptr<const PredecessorDAG> shared(value.release(), LiveDeleter {live_});
  entries_.emplace(std::move(key), shared);  // strong guarantee: a throw leaves the map unchanged
  num_entries_.fetch_add(1, std::memory_order_relaxed);
  bytes_.fetch_add(bytes, std::memory_order_relaxed);
  return shared;
}

void ShortestPathCache::freeze() {
  const detail::CacheGate::Lease lease(*gate_, true);
  gate_->freeze(true);
}

void ShortestPathCache::thaw() {
  const detail::CacheGate::Lease lease(*gate_, true);
  gate_->freeze(false);
}

void ShortestPathCache::clear() {
  const detail::CacheGate::Lease lease(*gate_, true);
  if (gate_->frozen()) {
    throw std::invalid_argument(
        "ShortestPathCache.clear is rejected while the cache is frozen (thaw() it, then clear())");
  }
  clear_contents();
}

void ShortestPathCache::clear_contents() noexcept {
  entries_.clear();
  num_entries_.store(0, std::memory_order_relaxed);
  bytes_.store(0, std::memory_order_relaxed);
}

ShortestPathCache::Stats ShortestPathCache::stats() const noexcept {
  Stats s;
  s.hits = hits_.load(std::memory_order_relaxed);
  s.misses = misses_.load(std::memory_order_relaxed);
  s.entries = num_entries_.load(std::memory_order_relaxed);
  s.bytes = bytes_.load(std::memory_order_relaxed);
  s.live_values = live_->load(std::memory_order_relaxed);
  return s;
}

} // namespace netgraph::core
