/* The evaluators of APPROACH 11: capacity (a FlowPolicy placement of +inf on a fresh ledger; the
   exact maximum flow is capacity() with routing Global and placement MaxFlow), batch_capacity
   (capacity per pair) and sensitivity (the loss of deleting each candidate edge). */
#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// The configuration, the static Eligibility, the DagCache and the optional outputs of one
// evaluation (the Eligibility and the cache outlive the call).
struct CapacityOptions {
  FlowPolicyConfig cfg;
  const Eligibility* elig {nullptr};
  DagCache* cache {nullptr};
  OutputOptions outputs {};
};
// The same fields; `outputs` apply to the baseline.
struct SensitivityOptions {
  FlowPolicyConfig cfg;
  const Eligibility* elig {nullptr};
  DagCache* cache {nullptr};
  OutputOptions outputs {};
};

enum class SensitivityStatus { Evaluated, Unevaluable };

// One candidate edge: loss = baseline placed - placed without the edge (signed, in units, never
// thresholded); an unevaluable deletion has no loss and the model exception's message in reason.
struct SensitivityEntry {
  EdgeId edge {0};
  std::optional<Amount> loss {};
  SensitivityStatus status {SensitivityStatus::Evaluated};
  std::string reason {};
};

struct SensitivityResult {
  PlacementResult baseline {};
  std::vector<SensitivityEntry> entries {};  // ascending EdgeId; an absent edge has loss 0
};

// Each evaluator is one public call on a supplied frozen cache's statistics (dag_cache.hpp): its
// hits count when it returns; a failed call changes no statistic, except one miss when it raises
// DagCacheMiss.

// place(+inf) of a FlowPolicy with opts.cfg on a fresh ledger of g (its capacities).
[[nodiscard]] PlacementResult capacity(const Graph& g, std::span<const NodeId> sources,
                                       std::span<const NodeId> destinations, const CapacityOptions& opts);

// capacity() per pair, in pair order, with opts.elig intersected with per_pair_elig[i] (null
// entries allowed; an empty span means none: the Python binding checks the length of every
// sequence it is given). Validated before any worker starts, the arguments first (workers >= 1,
// opts.elig's and the cache's graph, a supplied building cache only with workers == 1: it has
// one writer, while a frozen cache is shared by every worker; Pinned rejected: bundles bind
// endpoints; the per_pair_elig length; the
// configuration), then each pair completely in pair order (its endpoint sets, their overlap,
// its Eligibility's graph), so the lowest failing pair index is reported. Without a supplied
// cache each worker owns a private building cache, allocated before any thread starts. A worker
// exception aborts the batch and the exception of the lowest pair index is raised; an empty
// batch returns an empty list. A frozen cache's DagCacheMiss propagates unchanged; a per-pair
// intersection that differs from both operands has a fresh uid, absent from a frozen cache
// (APPROACH 10: already effective per-pair values with a null common value, or a common value
// containing them, reuse the prepopulated operand).
[[nodiscard]] std::vector<PlacementResult> batch_capacity(
    const Graph& g, const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>>& pairs,
    const CapacityOptions& opts, std::span<const Eligibility* const> per_pair_elig, std::int32_t workers);

// The baseline capacity() and one entry per candidate edge: under Global the edges with
// positive flow in the baseline; under every other model every edge opts.elig allows with at
// least 1 unit of static capacity. A loss is the baseline's placed minus the placed amount of
// capacity() with the edge deleted (opts.elig without it; the same configuration and cache).
// Global evaluates a deletion by the flow core's deletion re-solve from the baseline solve's
// final state (APPROACH 2.2: the minimum-cost flow of the maximum value without the edge, that
// value exactly), on one copy of the state per worker; an exact deletion on valid inputs always
// has an optimum, so any failure there propagates. CostOnly gives an edge outside every roster
// source's reachable part of the static support loss 0 without a placement, a result-equivalent
// predicate (no route uses the edge, so the deletion changes no route and no amount); every
// other deletion is placed on a fresh ledger, where a model exception (std::runtime_error) makes
// the entry unevaluable with its message as the reason. Every other exception propagates (the
// lowest candidate's). Workers and the cache
// rule as batch_capacity; a worker's private cache is cleared per candidate (fresh Eligibility
// uids never hit). CostOnly is the one model reading the cache: with a frozen cache a
// placed deletion raises DagCacheMiss (a fresh uid), which propagates unchanged (it is not a
// std::runtime_error, so never `unevaluable`); the other models read no cache.
[[nodiscard]] SensitivityResult sensitivity(const Graph& g, std::span<const NodeId> sources,
                                            std::span<const NodeId> destinations, const SensitivityOptions& opts,
                                            std::int32_t workers);

} // namespace netgraph::core
