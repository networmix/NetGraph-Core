/* The evaluators: evaluate_demand (a fresh FlowPolicy's evaluation on a fresh ledger, at the offer;
   the exact maximum flow is evaluate_demand of CentralizedTE with the optimized source allocation),
   evaluate_pairs (evaluate_demand per pair) and evaluate_edge_failures (the delivery reduction of
   deleting each candidate edge). */
#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/shortest_path_cache.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// The options of every evaluator: the configuration, the static Eligibility, the ShortestPathCache, the
// optional outputs of each result (evaluate_edge_failures: of its baseline), `design_demand` (T1's design,
// required and finite; T2's, absent meaning the offer when it is finite, else unlimited; every
// other model rejects one) and `offer`, a finite offer for every model (absent: the uncapped
// Lossless question, or the first-hop experiment of ShortestPaths Equal/Weighted Lossy; an inner
// nullopt is an explicit +inf, which the evaluators reject). The Eligibility and the cache outlive
// the call.
// | model                                | Lossless                                     | Lossy |
// | ShortestPaths Equal/Weighted         | no offer: the uncapped bound; finite: up to it | finite: that offer; none: the first-hop experiment |
// | ShortestPaths Optimized, CentralizedTE | no offer: the uncapped optimum; finite: up to it | rejected |
// | T1                                   | design required; no offer: the tunnels' lossless bound; finite: up to it | design and offer required |
// | T2                                   | design absent: the offer if finite, else unlimited; no offer: the bound at that design | finite offer; none: the design if finite |
// | Pinned                               | no offer: the fixed-program bound; finite: up to it | offer required |
struct EvaluationOptions {
  FlowPolicyConfig config;
  const Eligibility* eligibility {nullptr};
  ShortestPathCache* cache {nullptr};
  OutputOptions outputs {};
  std::optional<std::optional<Amount>> design_demand {};
  std::optional<std::optional<Amount>> offer {};
};

// One candidate edge: delivery_reduction = the baseline's call_delivered - call_delivered without
// the edge (signed, in units, never thresholded).
struct SensitivityEntry {
  EdgeId edge {0};
  Amount delivery_reduction {0};
};

struct SensitivityResult {
  PlacementResult baseline {};
  std::vector<SensitivityEntry> entries {};  // ascending EdgeId; an absent edge has delivery_reduction 0
};

// Each evaluator is one public call on a supplied frozen cache's statistics (shortest_path_cache.hpp): its
// hits count when it returns; a failed call changes no statistic, except one miss when it raises
// ShortestPathCacheMiss.

// The evaluation of a fresh FlowPolicy with options.config on a fresh ledger of graph (its capacities):
// place(offer) (absent: +inf), or for T1 and T2 evaluate(offer, design). Checked in this order:
// the policy's construction (the configuration, in the model matrix's order); then an explicitly given offer (`offer must
// be an integer in [0, 2^62)`, +inf included); then the design demand: only T1 and T2 take one
// (`<model> does not accept design_demand`), in [0, 2^62) or +inf (`design_demand must be an
// integer in [0, 2^62) or permitted +inf`), T1's finite (`DistributedTE requires a finite design
// demand`); then T1's missing design (`DistributedTE evaluation requires design_demand`); then the
// missing Lossy offer: T1 and Pinned require one, T2 takes its finite design instead
// (`<DistributedTE or CentralizedTE or Pinned> LOSSY evaluation requires a finite offer`).
[[nodiscard]] PlacementResult evaluate_demand(const Graph& graph, std::span<const NodeId> sources,
                                       std::span<const NodeId> destinations, const EvaluationOptions& options);

// Run evaluate_demand() per pair, returning results in pair order. Intersect options.eligibility with
// per_pair_eligibility[i]; null entries are allowed, and an empty span supplies no per-pair values.
// The Python binding checks the length of every supplied sequence.
// Before evaluation workers start, validate in order: workers >= 1;
// options.eligibility's and the cache's graph; Pinned (rejected: a program has fixed endpoints); the configuration;
// offer and design_demand (the rules of evaluate_demand); and per_pair_eligibility's length. Then validate each pair
// completely in index order: endpoint sets, overlap, Eligibility's graph, the pooled single-path rule for its sources
// and, for the tunnel models, the tunnel limit. This reports the lowest failing pair index. Last, under the cache's
// lease: a supplied building cache requires workers == 1. First-pass tunnel searches run in parallel
// over the checked prefix, before any lease; their immutable rosters are retained for evaluation.
// A supplied building cache requires workers == 1 because it allows one writer. Workers share
// a supplied frozen cache, and the one worker uses a supplied building cache; without a supplied
// cache each worker gets a private building cache before threads start.
// A worker exception aborts the batch; rethrow the exception from the lowest pair index.
// An empty batch returns an empty list. ShortestPathCacheMiss always propagates unchanged.
// An intersection equal to an operand reuses its uid and any prepopulated cache entry.
// Otherwise it gets a fresh uid, absent from a frozen cache. To reuse entries, supply already
// effective per-pair values with no common Eligibility, or one that contains them.
[[nodiscard]] std::vector<PlacementResult> evaluate_pairs(
    const Graph& graph, const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>>& pairs,
    const EvaluationOptions& options, std::span<const Eligibility* const> per_pair_eligibility, std::int32_t workers);

// Return the baseline evaluate_demand() and one entry per candidate edge. Unrestricted CentralizedTE
// with the optimized source allocation: the candidates carry positive baseline flow (deleting any
// other edge keeps the baseline feasible, at any offer); otherwise every edge options.eligibility
// allows, whatever its capacity (deleting an edge of capacity 0 can change routes and rosters, so
// its delivery_reduction can be nonzero; with the equal source allocation, a deletion that removes a
// source from the roster can raise the total, a negative delivery_reduction). delivery_reduction is
// the baseline's call_delivered minus the call_delivered with that edge deleted, using the same
// config, offer, design and cache. When the unrestricted optimized CentralizedTE baseline is a
// maximum flow (no offer, or an offer above what it delivers), each deletion starts from the
// baseline solve's final state and repairs it, using one state copy per worker: the exact maximum
// flow at minimum cost. At an offer the baseline delivers (it binds, perhaps below the maximum),
// each deletion is solved afresh at the offer. CentralizedTE with a support solves every deletion afresh;
// T1 sets up again for the same design under each deletion, T2 redoes its sizing, count selection
// and objectives; Pinned replays the same program under
// each deletion (a fixed-route failure question).
// ShortestPaths skips a deletion when no roster source's reachable part of the static shortest-path
// support contains the edge: no route uses it, so its delivery_reduction is 0. All other deletions use a fresh
// ledger. A failed deletion fails the call (no partial result), raising the exception of the
// lowest failing candidate index.
// Worker and cache rules are as in evaluate_pairs. Clear each worker's private cache per
// candidate because deletion Eligibility uids are fresh and cannot reuse entries.
// Only ShortestPaths reads the cache. A deletion that needs placement with a frozen cache raises
// ShortestPathCacheMiss for its fresh uid.
[[nodiscard]] SensitivityResult evaluate_edge_failures(const Graph& graph, std::span<const NodeId> sources,
                                            std::span<const NodeId> destinations, const EvaluationOptions& options,
                                            std::int32_t workers);

} // namespace netgraph::core
