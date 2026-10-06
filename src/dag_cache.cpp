/*
  The DagCache of APPROACH 10 in building mode: a map from (roots, direction, Eligibility uid)
  to the full static support of that search, with atomic statistics.
*/
#include "netgraph/core/dag_cache.hpp"

#include <stdexcept>
#include <string>
#include <utility>

#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/search.hpp"

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

DagCache::DagCache(const Graph& graph) : graph_uid_(graph.uid()) {}

std::shared_ptr<const PredDAG> DagCache::lookup(const Graph& g, std::span<const NodeId> roots,
                                                Direction direction, const Eligibility* elig) {
  // Invalid queries throw before anything is counted.
  if (g.uid() != graph_uid_) throw std::invalid_argument("DagCache: the graph is not the cache's graph (foreign graph)");
  if (elig != nullptr && elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  auto key_roots = detail::node_set(roots, g.num_nodes(), "roots");
  detail::check_member(direction, "direction");
  Key key {std::move(key_roots), direction, elig != nullptr ? elig->uid() : 0};
  const auto it = entries_.find(key);
  if (it != entries_.end()) {
    hits_.fetch_add(1, std::memory_order_relaxed);
    return it->second;
  }
  misses_.fetch_add(1, std::memory_order_relaxed);
  // A failed search (or a failed insertion) inserts nothing; the miss stays counted.
  auto value = std::make_unique<PredDAG>(detail::static_support(g, std::get<0>(key), direction, elig));
  const auto bytes = static_cast<std::int64_t>(detail::PredDAGBuilder::payload_bytes(*value));
  live_->fetch_add(1, std::memory_order_relaxed);  // the deleter runs even if the control block fails
  std::shared_ptr<const PredDAG> shared(value.release(), LiveDeleter {live_});
  entries_.emplace(std::move(key), shared);  // strong guarantee: a throw leaves the map unchanged
  num_entries_.fetch_add(1, std::memory_order_relaxed);
  bytes_.fetch_add(bytes, std::memory_order_relaxed);
  return shared;
}

void DagCache::clear() noexcept {
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
