/*
  evaluate_demand() is a fresh policy's place(offer) on a fresh ledger (absent: +inf; T1 and T2:
  evaluate(offer, design)). evaluate_pairs() runs it per pair and evaluate_edge_failures() per candidate
  deletion on `workers` threads; results follow
  input order and the lowest failing index supplies the exception. Unrestricted CentralizedTE with the
  optimized source allocation deletes only the edges its baseline loads. At a maximum baseline it
  derives scalar losses from residual component certificates or the first deletion repair,
  resetting each worker's private state independently, and solves afresh at a binding
  offer. ShortestPaths skips an edge that no
  roster source's reachable part of the static shortest-path support contains: removing it
  changes no route or amount. Other deletions (the equal source allocation of CentralizedTE among them: a deletion
  can raise its maximum; CentralizedTE with a support; T1 and T2, which set up or lay out again; Pinned,
  which replays its program under each deletion) use fresh ledgers. Only ShortestPaths reads a
  ShortestPathCache.
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
#include "netgraph/core/detail/min_cost_flow.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/route_builders.hpp"

namespace netgraph::core {

namespace {

// A ledger of graph for one evaluation: the graph is borrowed (the ledger dies within the call).
std::shared_ptr<const Graph> borrowed(const Graph& graph) {
  return std::shared_ptr<const Graph>(std::shared_ptr<const void> {}, &graph);
}

// The call an evaluator makes on its fresh policy: place(offer), or for T1 and T2
// evaluate(offer, design) (`full`).
struct EvaluatorCall {
  std::optional<Amount> offer {};   // nullopt: +inf (no offer, or the first-hop experiment)
  std::optional<Amount> design {};  // the tunnel models' design D (nullopt: unlimited)
  bool full {false};
};

bool in_domain(Amount x) noexcept { return x >= 0 && x < kDomainLimit; }

// The evaluators' amounts, after the policy's construction checked the configuration: an
// explicitly given offer (finite, in its domain); the design demand (its model, its domain, T1's
// finite); T1's missing design; then the absent values resolved and the missing Lossy offer.
EvaluatorCall evaluator_call(const FlowPolicyConfig& config, const std::optional<std::optional<Amount>>& offer,
                             const std::optional<std::optional<Amount>>& design) {
  const auto& r = config.routing;
  if (offer.has_value() && (!offer->has_value() || !in_domain(**offer))) {
    throw std::invalid_argument("offer must be an integer in [0, 2^62)");
  }
  const std::string model = detail::model_class(r);
  const bool t1 = std::holds_alternative<DistributedTE>(r);
  const auto* gl = std::get_if<CentralizedTE>(&r);
  const bool t2 = gl != nullptr && gl->tunnel_count.has_value();
  if (design.has_value()) {
    if (!t1 && !t2) throw std::invalid_argument(model + " does not accept design_demand");
    if (design->has_value() && !in_domain(**design)) {
      throw std::invalid_argument("design_demand must be an integer in [0, 2^62) or permitted +inf");
    }
    if (t1 && !design->has_value()) throw std::invalid_argument("DistributedTE requires a finite design demand");
  }
  if (t1 && !design.has_value()) throw std::invalid_argument("DistributedTE evaluation requires design_demand");
  EvaluatorCall call;
  if (offer.has_value()) call.offer = **offer;
  if (t1 || t2) {
    call.full = true;
    // T2's absent design is the offer when it is finite, else unlimited.
    call.design = design.has_value() ? *design : call.offer;
  }
  if (config.admission == Admission::Lossy && !offer.has_value()) {
    if (t2 && call.design.has_value()) {
      call.offer = call.design;
    } else if (t1 || t2 || std::holds_alternative<Pinned>(r)) {
      throw std::invalid_argument(model + " LOSSY evaluation requires a finite offer");
    }
    // ShortestPaths Equal/Weighted: the first-hop experiment (place(+inf)).
  }
  return call;
}

// The evaluation of a fresh policy.
PlacementResult run(FlowPolicy& policy, const EvaluatorCall& call, const OutputOptions& outputs) {
  if (call.full) return policy.evaluate(call.offer, std::optional<std::optional<Amount>>(std::in_place, call.design), outputs);
  return policy.place(call.offer, outputs);
}

// One evaluation inside a public call: a supplied frozen cache's hit goes to the call's tally.
PlacementResult evaluate(const Graph& graph, std::span<const NodeId> sources, std::span<const NodeId> destinations,
                         const EvaluationOptions& opts, const Eligibility* eligibility, ShortestPathCache* cache,
                         const OutputOptions& outputs, detail::FrozenCacheHitTally& tally, const detail::TunnelRoster* roster = nullptr) {
  FlowLedger fg(borrowed(graph));
  FlowPolicy policy = detail::PolicyAccess::make(fg, sources, destinations, opts.config, eligibility, cache, tally, roster);
  return run(policy, evaluator_call(policy.config(), opts.offer, opts.design_demand), outputs);
}

void check_common(const Graph& graph, const Eligibility* eligibility, const ShortestPathCache* cache, std::int32_t workers) {
  if (workers < 1) throw std::invalid_argument("workers must be >= 1 (got " + std::to_string(workers) + ")");
  if (eligibility != nullptr && eligibility->graph_uid() != graph.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  if (cache != nullptr && cache->graph_uid() != graph.uid()) {
    throw std::invalid_argument("foreign ShortestPathCache: bound to another graph");
  }
}

// Called only while the evaluator's single supplied-cache lease is held.
void check_cache_mode(const ShortestPathCache* cache, std::int32_t workers) {
  if (cache != nullptr && workers != 1 && !detail::DagCacheAccess::frozen(*cache)) {
    throw std::invalid_argument(
        "a supplied building ShortestPathCache has one writer: it is rejected unless workers == 1 (got workers = " +
        std::to_string(workers) + "; a frozen cache may be shared by every worker)");
  }
}

// The number of threads run_indexed uses for n tasks.
std::size_t worker_threads(std::size_t n, std::int32_t workers) {
  return std::max<std::size_t>(1, static_cast<std::size_t>(std::min<std::int64_t>(workers, static_cast<std::int64_t>(n))));
}

// A worker's cache for one evaluation: only ShortestPaths reads one (the other models reject a cache,
// and the caller's was validated before the workers started).
ShortestPathCache* worker_cache(const FlowPolicyConfig& config, ShortestPathCache& cache) noexcept {
  return std::holds_alternative<ShortestPaths>(config.routing) ? &cache : nullptr;
}

// Increasing indices ensure every task below a failing index has already been taken.
// Taken tasks finish even after an abort. Each exception has its own slot, and every
// launched thread is joined before the lowest-index exception is rethrown.
// The caller allocates result slots before entering; error and thread storage is private here.
template <class Task>
void run_tasks(std::size_t n, std::int32_t workers, Task task) {
  if (n == 0) return;
  const auto threads = worker_threads(n, workers);
  std::vector<std::exception_ptr> errors(n);
  std::atomic<std::size_t> next {0};
  std::atomic<bool> abort {false};
#if NETGRAPH_CORE_PRIVATE_HOOKS
  const auto* caller_probe = detail::current_concurrent_probe();
#endif
  auto worker = [&](std::size_t t) noexcept {
#if NETGRAPH_CORE_PRIVATE_HOOKS
    const detail::ScopedConcurrentProbe inherited(caller_probe);
#endif
    while (!abort.load(std::memory_order_relaxed)) {
      const std::size_t i = next.fetch_add(1, std::memory_order_relaxed);
      if (i >= n) break;
      try {
        detail::concurrent_event(detail::ConcurrentStage::WorkerBegin, &i);
        task(i, t);
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
      for (std::size_t t = 0; t < threads; ++t) {
        detail::concurrent_event(detail::ConcurrentStage::ThreadLaunch, &t);
        pool.emplace_back(worker, t);
      }
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

// Evaluation workers borrow the supplied cache or own one private building cache each.
// First-pass roster checks use run_tasks directly and need no cache or lease.
template <class Task>
void run_indexed(std::size_t n, std::int32_t workers, const Graph& graph, ShortestPathCache* shared, Task task) {
  if (n == 0) return;
  std::vector<std::unique_ptr<ShortestPathCache>> own;
  if (shared == nullptr) {
    const auto threads = worker_threads(n, workers);
    own.reserve(threads);
    for (std::size_t t = 0; t < threads; ++t) own.push_back(std::make_unique<ShortestPathCache>(graph));
  }
  run_tasks(n, workers, [&](std::size_t i, std::size_t t) {
    task(i, shared != nullptr ? *shared : *own[t], t);
  });
}

// ShortestPaths sensitivity: the edges of the full static support (every edge on a static shortest
// path to the destination set) reachable from the roster sources through support edges. This is
// not only the baseline's selected routes. Each source's route is read from the part of this set
// reachable from it, then selected on the priority values of the residual at the preparing call's
// entry (a fresh ledger's capacities). Deleting an edge outside the set leaves every reachable
// node's distance and support edges as they were (distances only grow, and no shortest path of a
// reachable node used the edge), and the roster as it was (a roster source keeps its support
// edges; a deletion never connects a source). The support ignores capacity, so this holds for
// edges of capacity 0 too. So every source's route, the pooled program (or the controller's
// support) built from them, and the view they are read on are unchanged, and so is every amount:
// the loss is exactly 0 and the deletion is not re-evaluated.
std::vector<std::uint8_t> reachable_support(const Graph& graph, const detail::PolicyEngine& pe) {
  const auto num_edges = static_cast<std::size_t>(graph.num_edges());
  std::vector<std::uint8_t> in_support(num_edges, 0);
  for (EdgeId e : detail::PredecessorDAGBuilder::ascending(*pe.support)) in_support[static_cast<std::size_t>(e)] = 1;
  std::vector<std::uint8_t> reached(num_edges, 0);
  std::vector<std::uint8_t> seen(static_cast<std::size_t>(graph.num_nodes()), 0);
  const auto row = graph.outgoing_offsets();
  const auto edge_destinations = graph.edge_destinations();
  std::vector<NodeId> stack(pe.reachable_sources.begin(), pe.reachable_sources.end());
  for (NodeId s : stack) seen[static_cast<std::size_t>(s)] = 1;
  while (!stack.empty()) {
    const auto u = static_cast<std::size_t>(stack.back());
    stack.pop_back();
    for (auto e = static_cast<std::size_t>(row[u]); e < static_cast<std::size_t>(row[u + 1]); ++e) {
      if (in_support[e] == 0) continue;
      reached[e] = 1;
      const auto v = static_cast<std::size_t>(edge_destinations[e]);
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

#if NETGRAPH_CORE_PRIVATE_HOOKS
namespace {
thread_local const SensitivityProbe* edge_failure_probe = nullptr;
}
const SensitivityProbe* current_edge_failure_probe() noexcept { return edge_failure_probe; }
void set_edge_failure_probe(const SensitivityProbe* probe) noexcept { edge_failure_probe = probe; }
#endif

std::optional<TunnelRoster> check_demand_arguments(const Graph& graph, std::span<const NodeId> sources,
    std::span<const NodeId> destinations, const EvaluationOptions& options, ArgumentWork search) {
  auto roster = check_policy_construction(graph, sources, destinations, options.config, options.eligibility, options.cache, search);
  (void)evaluator_call(options.config, options.offer, options.design_demand);
  return roster;
}

std::optional<TunnelRoster> check_edge_failure_arguments(const Graph& graph, std::span<const NodeId> sources,
    std::span<const NodeId> destinations, const EvaluationOptions& options, std::int32_t workers, ArgumentWork search) {
  check_common(graph, options.eligibility, options.cache, workers);
  return check_demand_arguments(graph, sources, destinations, options, search);
}

void check_batch_arguments(const Graph& graph, const EvaluationOptions& options, std::int32_t workers) {
  check_common(graph, options.eligibility, options.cache, workers);
  if (std::holds_alternative<Pinned>(options.config.routing)) throw std::invalid_argument("evaluate_pairs does not accept Pinned");
  if (const auto* global = std::get_if<CentralizedTE>(&options.config.routing); global != nullptr && global->support &&
                                                                       global->support->graph_uid() != graph.uid()) {
    throw std::invalid_argument("support belongs to another graph");
  }
  check_config(options.config, options.cache != nullptr, graph.num_nodes(), graph.num_edges(), std::nullopt);
  (void)evaluator_call(options.config, options.offer, options.design_demand);
}

void check_batch_pair(const Graph& graph, const FlowPolicyConfig& config, std::span<const NodeId> sources,
                      std::span<const NodeId> destinations, const Eligibility* eligibility, std::size_t index) {
  const auto s = node_set(sources, graph.num_nodes(), "sources");
  const auto t = node_set(destinations, graph.num_nodes(), "destinations");
  check_disjoint(s, t, "sources", "destinations");
  if (eligibility != nullptr && eligibility->graph_uid() != graph.uid()) {
    throw std::invalid_argument("foreign Eligibility: per_pair_eligibility[" + std::to_string(index) +
                                "] built for another graph");
  }
  check_source_count(config, s.size());
}

std::vector<TunnelRoster> prepare_batch_rosters(const Graph& graph,
    std::span<const std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs,
    const EvaluationOptions& options, std::span<const Eligibility* const> per_pair, std::int32_t workers,
    ArgumentWork search) {
  const auto* centralized = std::get_if<CentralizedTE>(&options.config.routing);
  if (pairs.empty() || !(std::holds_alternative<DistributedTE>(options.config.routing) ||
                        (centralized && centralized->tunnel_count.has_value()))) return {};
  std::vector<TunnelRoster> rosters;
  auto work = [&] {
    rosters.resize(pairs.size());
    run_tasks(pairs.size(), workers, [&](std::size_t i, std::size_t /*t*/) {
      const auto* common = options.eligibility;
      const auto* eligibility = per_pair.empty() ? nullptr : per_pair[i];
      std::optional<Eligibility> both;
      if (common && eligibility) both = Eligibility::intersect(*common, *eligibility);
      const auto* mask = both ? &*both : (common ? common : eligibility);
      rosters[i] = *check_policy_construction(graph, pairs[i].first, pairs[i].second, options.config, mask, nullptr);
    });
  };
  if (search != nullptr) search([](void* context) { (*static_cast<decltype(work)*>(context))(); }, &work);
  else work();
  return rosters;
}

PlacementResult evaluate_demand_in_call(const Graph& graph, std::span<const NodeId> sources,
                                 std::span<const NodeId> destinations, const EvaluationOptions& options, FrozenCacheHitTally& tally, const TunnelRoster* roster) {
  return evaluate(graph, sources, destinations, options, options.eligibility, options.cache, options.outputs, tally, roster);
}

std::vector<PlacementResult> evaluate_pairs_in_call(
    const Graph& graph, const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>>& pairs,
    const EvaluationOptions& options, std::span<const Eligibility* const> per_pair_eligibility, std::int32_t workers,
    FrozenCacheHitTally& tally, std::span<const TunnelRoster> rosters) {
  check_cache_mode(options.cache, workers);
  std::vector<std::optional<PlacementResult>> results(pairs.size());
  run_indexed(pairs.size(), workers, graph, options.cache, [&](std::size_t i, ShortestPathCache& cache, std::size_t /*t*/) {
    const Eligibility* extra = per_pair_eligibility.empty() ? nullptr : per_pair_eligibility[i];
    std::optional<Eligibility> both;
    const Eligibility* elig = options.eligibility;
    if (extra != nullptr) {
      if (elig == nullptr) {
        elig = extra;
      } else {
        both = Eligibility::intersect(*options.eligibility, *extra);
        elig = &*both;
      }
    }
    results[i] = evaluate(graph, pairs[i].first, pairs[i].second, options, elig, worker_cache(options.config, cache),
                          options.outputs, tally, rosters.empty() ? nullptr : &rosters[i]);
  });
  std::vector<PlacementResult> out;
  out.reserve(results.size());
  for (auto& r : results) out.push_back(std::move(*r));
  return out;
}

SensitivityResult evaluate_edge_failures_in_call(const Graph& graph, std::span<const NodeId> sources,
                                      std::span<const NodeId> destinations, const EvaluationOptions& options,
                                      std::int32_t workers, FrozenCacheHitTally& tally, const TunnelRoster* roster) {
  check_cache_mode(options.cache, workers);
  // Unrestricted CentralizedTE with the optimized source allocation skips the links its baseline leaves without
  // flow (deleting one keeps the baseline feasible, so the delivery stays). When the baseline is a maximum flow (no
  // offer, or a finite offer above what it delivers) each loss uses a certificate or the first
  // repair from an independent baseline reset; no partial allocation is published. When the offer
  // binds (the baseline delivers it, perhaps below the maximum) each deletion is solved afresh at the offer. The
  // equal source allocation, a support and T2 solve every allowed deletion afresh.
  const auto* global_routing = std::get_if<CentralizedTE>(&options.config.routing);
  const bool global = global_routing != nullptr && !global_routing->support && !global_routing->tunnel_count.has_value() &&
                      global_routing->source_allocation == SourceAllocation::Optimized;
  const bool shortest_paths = std::holds_alternative<ShortestPaths>(options.config.routing);
  SensitivityResult out;
  std::vector<std::uint8_t> on_route;  // ShortestPaths: the baseline's reachable support
  detail::DeletionCapture capture;     // optimized CentralizedTE: the baseline solve's final state
  {
    std::unique_ptr<ShortestPathCache> own;
    ShortestPathCache* cache = options.cache;
    if (cache == nullptr && shortest_paths) {
      own = std::make_unique<ShortestPathCache>(graph);
      cache = own.get();
    }
    FlowLedger fg(borrowed(graph));
    FlowPolicy policy = PolicyAccess::make(fg, sources, destinations, options.config, options.eligibility, cache, tally, roster);
    auto& pe = PolicyAccess::engine(policy);
    if (global) pe.capture = &capture;
    out.baseline = run(policy, evaluator_call(policy.config(), options.offer, options.design_demand), options.outputs);
    pe.capture = nullptr;
    if (shortest_paths) on_route = reachable_support(graph, pe);
  }
  std::vector<EdgeId> candidates;
  std::vector<Amount> baseline_flows;
  if (global) {
    // The baseline's positive-flow edges; no state, and no candidate, when no roster source
    // reaches an allowed destination (nothing is placed).
    if (capture.state.has_value()) {
      if (capture.state->value() != out.baseline.call_delivered) {
        throw std::logic_error("evaluate_edge_failures: the captured state differs from the baseline");
      }
      baseline_flows = capture.state->edge_flows();
      for (EdgeId e = 0; e < graph.num_edges(); ++e) {
        if (baseline_flows[static_cast<std::size_t>(e)] > 0) candidates.push_back(e);
      }
    }
  } else {
    // Every allowed edge: an edge of capacity 0 can carry a route (and a roster source), and a
    // deletion can raise an equal source allocation's total (its source then leaves the roster).
    for (EdgeId e = 0; e < graph.num_edges(); ++e) {
      if (options.eligibility == nullptr || options.eligibility->allows_edge(e)) candidates.push_back(e);
    }
  }
  const Amount base = out.baseline.call_delivered;
  const bool repair = global && (!options.offer.has_value() || base < **options.offer);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  const auto* probe = current_edge_failure_probe();
  const bool full_repair = probe != nullptr && probe->full_repair;
  const bool classify = !full_repair && (probe == nullptr || probe->cut_edges);
#else
  constexpr bool full_repair = false;
  constexpr bool classify = true;
#endif
  out.entries.resize(candidates.size());
  std::vector<std::uint8_t> cut_edges;
  if (repair && !candidates.empty() && classify) cut_edges = detail::deletion_cut_edges(*capture.state);
  const bool needs_resolver = repair && std::any_of(candidates.begin(), candidates.end(), [&](EdgeId e) {
    return cut_edges.empty() || cut_edges[static_cast<std::size_t>(e)] == 0;
  });
  // CentralizedTE: one resolver per worker (its copy of the state and its working storage), made here.
  std::vector<detail::DeletionResolver> resolvers;
  if (needs_resolver) {
    const auto threads = worker_threads(candidates.size(), workers);
    resolvers.reserve(threads);
    for (std::size_t t = 0; t < threads; ++t) resolvers.emplace_back(*capture.state);
  }
  run_indexed(candidates.size(), workers, graph, options.cache, [&](std::size_t i, ShortestPathCache& cache, std::size_t t) {
    const EdgeId e = candidates[i];
    auto& entry = out.entries[i];
    entry.edge = e;
    // A failed deletion fails the call: run_indexed raises the lowest failing candidate's exception.
    if (repair) {
      const auto k = static_cast<std::size_t>(e);
      const bool cut = !cut_edges.empty() && cut_edges[k] != 0;
#if NETGRAPH_CORE_PRIVATE_HOOKS
      if (probe != nullptr && probe->observe != nullptr) {
        probe->observe(probe->context, e, full_repair ? SensitivityWork::FullRepair : cut ? SensitivityWork::Cut : SensitivityWork::FirstRepair);
      }
#endif
      if (cut) entry.delivery_reduction = baseline_flows[k];
      else if (full_repair) entry.delivery_reduction = base - resolvers[t].resolve(e).value;
      else entry.delivery_reduction = resolvers[t].loss(e);
      return;
    }
    if (shortest_paths && on_route[static_cast<std::size_t>(e)] == 0) {
      entry.delivery_reduction = 0;
    } else {
      const EdgeId one[1] = {e};
      Eligibility deleted = graph.eligibility({}, one);
      if (options.eligibility != nullptr) deleted = Eligibility::intersect(*options.eligibility, deleted);
      const auto r = evaluate(graph, sources, destinations, options, &deleted, worker_cache(options.config, cache),
                              OutputOptions {}, tally);
      entry.delivery_reduction = base - r.call_delivered;
    }
    if (&cache != options.cache) detail::DagCacheAccess::clear_private(cache);
  });
  return out;
}

std::int32_t default_workers(bool edge_failures) noexcept {
  const char* v = edge_failures ? std::getenv("NGRAPH_CORE_SENSITIVITY_THREADS") : std::getenv("NGRAPH_CORE_BATCH_THREADS");
  if (v == nullptr) return 1;
  char* end = nullptr;
  const long n = std::strtol(v, &end, 10);
  if (end == v || *end != '\0' || n < 1 || n > 4096) return 1;
  return static_cast<std::int32_t>(n);
}

} // namespace detail

// Each evaluator is one public call on a supplied frozen cache's statistics.
PlacementResult evaluate_demand(const Graph& graph, std::span<const NodeId> input_sources,
                         std::span<const NodeId> input_destinations, const EvaluationOptions& input_options) {
  const EvaluationOptions options = input_options;
  const std::vector<NodeId> sources(input_sources.begin(), input_sources.end());
  const std::vector<NodeId> destinations(input_destinations.begin(), input_destinations.end());
  const auto roster = detail::check_demand_arguments(graph, sources, destinations, options);
  return detail::frozen_call(options.cache, [&](detail::FrozenCacheHitTally& tally) {
    return detail::evaluate_demand_in_call(graph, sources, destinations, options, tally, roster ? &*roster : nullptr);
  });
}

std::vector<PlacementResult> evaluate_pairs(
    const Graph& graph, const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>>& input_pairs,
    const EvaluationOptions& input_options, std::span<const Eligibility* const> input_per_pair, std::int32_t workers) {
  const EvaluationOptions options = input_options;
  const auto pairs = input_pairs;
  const std::vector<const Eligibility*> per_pair_eligibility(input_per_pair.begin(), input_per_pair.end());
  detail::check_batch_arguments(graph, options, workers);
  if (!per_pair_eligibility.empty() && per_pair_eligibility.size() != pairs.size()) {
    throw std::invalid_argument("per_pair_eligibility has " + std::to_string(per_pair_eligibility.size()) + " entries for " +
                                std::to_string(pairs.size()) + " pairs");
  }
  std::exception_ptr pending;
  std::size_t valid = 0;
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    try {
      detail::check_batch_pair(graph, options.config, pairs[i].first, pairs[i].second,
          per_pair_eligibility.empty() ? nullptr : per_pair_eligibility[i], i);
      ++valid;
    } catch (...) {
      pending = std::current_exception();
      break;
    }
  }
  // Searches of the valid prefix precede the later pair's attached-check error.
  // The search workers themselves report the lowest failing pair and are joined here.
  const auto rosters = detail::prepare_batch_rosters(graph, std::span(pairs).first(valid), options, per_pair_eligibility, workers);
  if (pending) std::rethrow_exception(pending);
  return detail::frozen_call(options.cache, [&](detail::FrozenCacheHitTally& tally) {
    return detail::evaluate_pairs_in_call(graph, pairs, options, per_pair_eligibility, workers, tally, rosters);
  });
}

SensitivityResult evaluate_edge_failures(const Graph& graph, std::span<const NodeId> input_sources,
                              std::span<const NodeId> input_destinations, const EvaluationOptions& input_options,
                              std::int32_t workers) {
  const EvaluationOptions options = input_options;
  const std::vector<NodeId> sources(input_sources.begin(), input_sources.end());
  const std::vector<NodeId> destinations(input_destinations.begin(), input_destinations.end());
  const auto roster = detail::check_edge_failure_arguments(graph, sources, destinations, options, workers);
  return detail::frozen_call(options.cache, [&](detail::FrozenCacheHitTally& tally) {
    return detail::evaluate_edge_failures_in_call(graph, sources, destinations, options, workers, tally, roster ? &*roster : nullptr);
  });
}

} // namespace netgraph::core
