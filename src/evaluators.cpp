/*
  capacity() is place(+inf) on a fresh ledger. batch_capacity() runs it per pair and
  sensitivity() per candidate deletion on `workers` threads; results follow input order and
  the lowest failing index supplies the exception. Global repairs the baseline flow after
  each deletion. CostOnly skips an edge that no roster source's reachable part of the static
  shortest-path support contains: removing it changes no route or amount. Other deletions use
  fresh ledgers.
*/
#include "netgraph/core/evaluators.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>

#include "netgraph/core/detail/cache_access.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/evaluator_checks.hpp"
#include "netgraph/core/detail/flow_core.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"

namespace netgraph::core {

namespace {

// A ledger of g for one evaluation: the graph is borrowed (the ledger dies within the call).
std::shared_ptr<const Graph> borrowed(const Graph& g) {
  return std::shared_ptr<const Graph>(std::shared_ptr<const void> {}, &g);
}

// One evaluation inside a public call: a supplied frozen cache's hit goes to the call's tally.
PlacementResult evaluate(const Graph& g, std::span<const NodeId> sources,
                         std::span<const NodeId> destinations, const FlowPolicyConfig& cfg, const Eligibility* elig,
                         DagCache* cache, const OutputOptions& outputs, detail::FrozenTally& tally) {
  FlowLedger fg(borrowed(g));
  FlowPolicy policy = detail::PolicyAccess::make(fg, sources, destinations, cfg, elig, cache, tally);
  return policy.place(std::nullopt, outputs);
}

void check_common(const Graph& g, const Eligibility* elig, const DagCache* cache, std::int32_t workers) {
  if (workers < 1) throw std::invalid_argument("workers must be >= 1 (got " + std::to_string(workers) + ")");
  if (elig != nullptr && elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  if (cache != nullptr && cache->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign DagCache: bound to another graph");
  }
  if (cache != nullptr && workers != 1 && !detail::DagCacheAccess::frozen(*cache)) {
    throw std::invalid_argument(
        "a supplied building DagCache has one writer: it is rejected unless workers == 1 (got workers = " +
        std::to_string(workers) + "; a frozen cache may be shared by every worker)");
  }
}

// The number of threads run_indexed uses for n tasks.
std::size_t worker_threads(std::size_t n, std::int32_t workers) {
  return std::max<std::size_t>(1, static_cast<std::size_t>(std::min<std::int64_t>(workers, static_cast<std::int64_t>(n))));
}

// Runs task(i, cache, t) for i in [0, n) on worker_threads(n, workers) threads (t the worker's
// index), each with its own private building cache unless `shared` is given (a supplied cache:
// building with workers == 1, or frozen and shared by every worker). Indices are taken in increasing
// order and a taken index always runs, so after an exception every lower index has run: the
// lowest failing index is the one a sequential run would report, and its exception is rethrown.
// Everything a worker needs (the result slots, the private caches) is allocated here, before
// any thread starts, so a failure throws in the caller; the worker body cannot throw (an
// exception escaping a std::thread terminates the process): a task's exception is caught into
// its slot, and every launched thread is joined before anything is rethrown.
template <class Task>
void run_indexed(std::size_t n, std::int32_t workers, const Graph& g, DagCache* shared, Task task) {
  if (n == 0) return;
  const auto threads = worker_threads(n, workers);
  std::vector<std::exception_ptr> errors(n);
  std::vector<std::unique_ptr<DagCache>> own;
  if (shared == nullptr) {
    own.reserve(threads);
    for (std::size_t t = 0; t < threads; ++t) own.push_back(std::make_unique<DagCache>(g));
  }
  std::atomic<std::size_t> next {0};
  std::atomic<bool> abort {false};
  auto worker = [&](std::size_t t) noexcept {
    DagCache& cache = shared != nullptr ? *shared : *own[t];
    while (!abort.load(std::memory_order_relaxed)) {
      const std::size_t i = next.fetch_add(1, std::memory_order_relaxed);
      if (i >= n) break;
      try {
        task(i, cache, t);
      } catch (...) {
        errors[i] = std::current_exception();
        abort.store(true, std::memory_order_relaxed);
      }
    }
  };
  if (threads == 1) {
    worker(0);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(threads);
    try {
      for (std::size_t t = 0; t < threads; ++t) pool.emplace_back(worker, t);
    } catch (...) {
      abort.store(true);
      for (auto& th : pool) th.join();
      throw;
    }
    for (auto& th : pool) th.join();
  }
  for (const auto& e : errors) {
    if (e) std::rethrow_exception(e);
  }
}

// CostOnly sensitivity: the edges of the full static support (every edge on a static shortest
// path to the destination set) reachable from the roster sources through support edges. This is
// not only the baseline's selected routes. A source's route is read from the part of this set
// reachable from it, then selected on the priority values of the working residual when that
// source's flow is created; earlier creations may already have used capacity, even on a fresh
// ledger. Deleting an edge outside the set leaves every reachable node's distance and support
// edges as they were (distances only grow, and no shortest path of a reachable node used the
// edge). So, creation by creation in their order, each sees the same reachable support and the
// same residual, and every route and placement is unchanged: the loss is exactly 0 and the
// deletion is not re-evaluated.
std::vector<std::uint8_t> reachable_support(const Graph& g, const detail::PolicyEngine& pe) {
  const auto num_edges = static_cast<std::size_t>(g.num_edges());
  std::vector<std::uint8_t> in_support(num_edges, 0);
  for (EdgeId e : detail::PredDAGBuilder::ascending(*pe.support)) in_support[static_cast<std::size_t>(e)] = 1;
  std::vector<std::uint8_t> reached(num_edges, 0);
  std::vector<std::uint8_t> seen(static_cast<std::size_t>(g.num_nodes()), 0);
  const auto row = g.row_offsets();
  const auto edge_dst = g.edge_dst();
  std::vector<NodeId> stack(pe.roster.begin(), pe.roster.end());
  for (NodeId s : stack) seen[static_cast<std::size_t>(s)] = 1;
  while (!stack.empty()) {
    const auto u = static_cast<std::size_t>(stack.back());
    stack.pop_back();
    for (auto e = static_cast<std::size_t>(row[u]); e < static_cast<std::size_t>(row[u + 1]); ++e) {
      if (in_support[e] == 0) continue;
      reached[e] = 1;
      const auto v = static_cast<std::size_t>(edge_dst[e]);
      if (seen[v] == 0) {
        seen[v] = 1;
        stack.push_back(static_cast<NodeId>(v));
      }
    }
  }
  return reached;
}

} // namespace

namespace detail {

void check_batch_arguments(const Graph& g, const CapacityOptions& opts, std::int32_t workers) {
  check_common(g, opts.elig, opts.cache, workers);
  if (std::holds_alternative<Pinned>(opts.cfg.routing)) {
    throw std::invalid_argument("Pinned routing is rejected in batch_capacity (bundles bind their endpoints)");
  }
  validate_config(opts.cfg, 2, 2);  // Pinned (the endpoint-count rule) is rejected above
}

void check_batch_pair(const Graph& g, std::span<const NodeId> sources,
                      std::span<const NodeId> destinations, const Eligibility* elig, std::size_t index) {
  const auto s = node_set(sources, g.num_nodes(), "sources");
  const auto t = node_set(destinations, g.num_nodes(), "destinations");
  check_disjoint(s, t, "sources", "destinations");
  if (elig != nullptr && elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: per_pair_elig[" + std::to_string(index) +
                                "] built for another graph");
  }
}

PlacementResult capacity_in_call(const Graph& g, std::span<const NodeId> sources,
                                 std::span<const NodeId> destinations, const CapacityOptions& opts, FrozenTally& tally) {
  return evaluate(g, sources, destinations, opts.cfg, opts.elig, opts.cache, opts.outputs, tally);
}

std::vector<PlacementResult> batch_capacity_in_call(
    const Graph& g, const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>>& pairs,
    const CapacityOptions& opts, std::span<const Eligibility* const> per_pair_elig, std::int32_t workers,
    FrozenTally& tally) {
  check_batch_arguments(g, opts, workers);
  if (!per_pair_elig.empty() && per_pair_elig.size() != pairs.size()) {
    throw std::invalid_argument("per_pair_elig has " + std::to_string(per_pair_elig.size()) + " entries for " +
                                std::to_string(pairs.size()) + " pairs");
  }
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    check_batch_pair(g, pairs[i].first, pairs[i].second, per_pair_elig.empty() ? nullptr : per_pair_elig[i], i);
  }
  std::vector<std::optional<PlacementResult>> results(pairs.size());
  run_indexed(pairs.size(), workers, g, opts.cache, [&](std::size_t i, DagCache& cache, std::size_t /*t*/) {
    const Eligibility* extra = per_pair_elig.empty() ? nullptr : per_pair_elig[i];
    std::optional<Eligibility> both;
    const Eligibility* elig = opts.elig;
    if (extra != nullptr) {
      if (elig == nullptr) {
        elig = extra;
      } else {
        both = Eligibility::intersect(*opts.elig, *extra);
        elig = &*both;
      }
    }
    results[i] = evaluate(g, pairs[i].first, pairs[i].second, opts.cfg, elig, &cache, opts.outputs, tally);
  });
  std::vector<PlacementResult> out;
  out.reserve(results.size());
  for (auto& r : results) out.push_back(std::move(*r));
  return out;
}

SensitivityResult sensitivity_in_call(const Graph& g, std::span<const NodeId> sources,
                                      std::span<const NodeId> destinations, const SensitivityOptions& opts,
                                      std::int32_t workers, FrozenTally& tally) {
  check_common(g, opts.elig, opts.cache, workers);
  const bool global = std::holds_alternative<Global>(opts.cfg.routing);
  const bool cost_only = std::holds_alternative<CostOnly>(opts.cfg.routing);
  SensitivityResult out;
  std::vector<std::uint8_t> on_route;  // CostOnly: the baseline's reachable support
  detail::DeletionCapture capture;     // Global: the baseline solve's final state
  {
    std::unique_ptr<DagCache> own;
    DagCache* cache = opts.cache;
    if (cache == nullptr) {
      own = std::make_unique<DagCache>(g);
      cache = own.get();
    }
    FlowLedger fg(borrowed(g));
    FlowPolicy policy = PolicyAccess::make(fg, sources, destinations, opts.cfg, opts.elig, cache, tally);
    auto& pe = PolicyAccess::engine(policy);
    if (global) pe.capture = &capture;
    out.baseline = policy.place(std::nullopt, opts.outputs);
    pe.capture = nullptr;
    if (cost_only) on_route = reachable_support(g, pe);
  }
  std::vector<EdgeId> candidates;
  if (global) {
    // The baseline's positive-flow edges; no state, and no candidate, when no roster source
    // reaches an allowed destination (nothing is placed).
    if (capture.state.has_value()) {
      if (capture.state->value() != out.baseline.placed) {
        throw std::logic_error("sensitivity: the captured state differs from the baseline");
      }
      const auto flows = capture.state->edge_flows();
      for (EdgeId e = 0; e < g.num_edges(); ++e) {
        if (flows[static_cast<std::size_t>(e)] > 0) candidates.push_back(e);
      }
    }
  } else {
    const auto capacity_units = g.capacity_units();
    for (EdgeId e = 0; e < g.num_edges(); ++e) {
      if ((opts.elig == nullptr || opts.elig->edge_ok(e)) && capacity_units[static_cast<std::size_t>(e)] >= 1) {
        candidates.push_back(e);
      }
    }
  }
  const Amount base = out.baseline.placed;
  out.entries.resize(candidates.size());
  // Global: one resolver per worker (its copy of the state and its working storage), made here.
  std::vector<detail::DeletionResolver> resolvers;
  if (global && !candidates.empty()) {
    const auto threads = worker_threads(candidates.size(), workers);
    resolvers.reserve(threads);
    for (std::size_t t = 0; t < threads; ++t) resolvers.emplace_back(*capture.state);
  }
  run_indexed(candidates.size(), workers, g, opts.cache, [&](std::size_t i, DagCache& cache, std::size_t t) {
    const EdgeId e = candidates[i];
    auto& entry = out.entries[i];
    entry.edge = e;
    // A failed deletion fails the call: run_indexed raises the lowest failing candidate's exception.
    if (global) {
      entry.loss = base - resolvers[t].resolve(e).value;
      return;
    }
    if (cost_only && on_route[static_cast<std::size_t>(e)] == 0) {
      entry.loss = 0;
    } else {
      const EdgeId one[1] = {e};
      Eligibility deleted = g.eligibility({}, one);
      if (opts.elig != nullptr) deleted = Eligibility::intersect(*opts.elig, deleted);
      const auto r = evaluate(g, sources, destinations, opts.cfg, &deleted, &cache, OutputOptions {}, tally);
      entry.loss = base - r.placed;
    }
    if (&cache != opts.cache) cache.clear();
  });
  return out;
}

std::int32_t default_workers(bool sensitivity) noexcept {
  const char* v = sensitivity ? std::getenv("NGRAPH_CORE_SENSITIVITY_THREADS") : std::getenv("NGRAPH_CORE_BATCH_THREADS");
  if (v == nullptr) return 1;
  char* end = nullptr;
  const long n = std::strtol(v, &end, 10);
  if (end == v || *end != '\0' || n < 1 || n > 4096) return 1;
  return static_cast<std::int32_t>(n);
}

} // namespace detail

// Each evaluator is one public call on a supplied frozen cache's statistics.
PlacementResult capacity(const Graph& g, std::span<const NodeId> sources,
                         std::span<const NodeId> destinations, const CapacityOptions& opts) {
  return detail::frozen_call(opts.cache, [&](detail::FrozenTally& tally) {
    return detail::capacity_in_call(g, sources, destinations, opts, tally);
  });
}

std::vector<PlacementResult> batch_capacity(
    const Graph& g, const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>>& pairs,
    const CapacityOptions& opts, std::span<const Eligibility* const> per_pair_elig, std::int32_t workers) {
  return detail::frozen_call(opts.cache, [&](detail::FrozenTally& tally) {
    return detail::batch_capacity_in_call(g, pairs, opts, per_pair_elig, workers, tally);
  });
}

SensitivityResult sensitivity(const Graph& g, std::span<const NodeId> sources,
                              std::span<const NodeId> destinations, const SensitivityOptions& opts,
                              std::int32_t workers) {
  return detail::frozen_call(opts.cache, [&](detail::FrozenTally& tally) {
    return detail::sensitivity_in_call(g, sources, destinations, opts, workers, tally);
  });
}

} // namespace netgraph::core
