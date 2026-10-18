/* The evaluators: capacity (a FlowPolicy placement of +inf on
   a fresh ledger; the exact maximum flow is capacity() with routing Global and the optimized
   origin), batch_capacity (capacity per pair) and sensitivity (the loss of deleting each
   candidate edge). */
#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// The options of every evaluator: the configuration, the static Eligibility, the DagCache, the
// optional outputs of each result (sensitivity: of its baseline) and `offered`, the offer of a
// Pinned Lossy evaluation (absent: none; the inner nullopt is +inf, which its rule rejects). The
// Eligibility and the cache outlive the call.
// | model                         | Lossless                         | Lossy |
// | CostOnly EqualBalanced/Weighted | the uncapped bound (no offered) | the first-hop experiment (no offered) |
// | CostOnly MaxFlow, Global       | the uncapped optimum (no offered) | rejected |
// | CapacityAware                  | the uncapped placement (no offered) | rejected |
// | Pinned                         | the fixed-program bound (no offered) | offered required and finite |
struct QueryOptions {
  FlowPolicyConfig config;
  const Eligibility* elig {nullptr};
  DagCache* cache {nullptr};
  OutputOptions outputs {};
  std::optional<std::optional<Amount>> offered {};
};

// One candidate edge: loss = baseline placed - placed without the edge (signed, in units, never
// thresholded).
struct SensitivityEntry {
  EdgeId edge {0};
  Amount loss {0};
};

struct SensitivityResult {
  PlacementResult baseline {};
  std::vector<SensitivityEntry> entries {};  // ascending EdgeId; an absent edge has loss 0
};

// Each evaluator is one public call on a supplied frozen cache's statistics (dag_cache.hpp): its
// hits count when it returns; a failed call changes no statistic, except one miss when it raises
// DagCacheMiss.

// place(+inf) of a FlowPolicy with opts.config on a fresh ledger of g (its capacities); Pinned
// Lossy: place(offered). `offered` is checked after the policy's construction: an amount in [0,
// 2^62) or +inf; a Lossless model rejects it (`<model> LOSSLESS capacity does not accept
// offered`), CostOnly Lossy too (`CostOnly LOSSY capacity uses first-hop injection`), and Pinned
// Lossy requires it finite (`Pinned LOSSY capacity requires finite offered`).
[[nodiscard]] PlacementResult capacity(const Graph& g, std::span<const NodeId> sources,
                                       std::span<const NodeId> destinations, const QueryOptions& opts);

// Run capacity() per pair, returning results in pair order. Intersect opts.elig with
// per_pair_elig[i]; null entries are allowed, and an empty span supplies no per-pair values.
// The Python binding checks the length of every supplied sequence.
// Before starting workers, validate in order: workers >= 1; opts.elig's and the cache's graph;
// a supplied building cache with workers > 1; Pinned (rejected: a program has fixed endpoints);
// the configuration; offered (the rule of capacity); and per_pair_elig's length. Then validate each pair completely in index
// order: endpoint sets, overlap, Eligibility's graph, the pooled single-path rule for its sources.
// This reports the lowest failing pair index.
// A supplied building cache requires workers == 1 because it allows one writer. Workers share
// a supplied frozen cache, and the one worker uses a supplied building cache; without a supplied
// cache each worker gets a private building cache before threads start.
// A worker exception aborts the batch; rethrow the exception from the lowest pair index.
// An empty batch returns an empty list. DagCacheMiss always propagates unchanged.
// An intersection equal to an operand reuses its uid and any prepopulated cache entry.
// Otherwise it gets a fresh uid, absent from a frozen cache. To reuse entries, supply already
// effective per-pair values with no common Eligibility, or one that contains them.
[[nodiscard]] std::vector<PlacementResult> batch_capacity(
    const Graph& g, const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>>& pairs,
    const QueryOptions& opts, std::span<const Eligibility* const> per_pair_elig, std::int32_t workers);

// Return baseline capacity() and one entry per candidate edge. Unrestricted Global candidates with
// the optimized origin carry positive baseline flow; otherwise every edge opts.elig allows, whatever
// its capacity (deleting an edge of capacity 0 can change routes and rosters, so its loss can be
// nonzero; with the equal origin, a deletion that removes a source from the roster can raise the
// total, a negative loss). Loss is baseline placed minus placed with that edge deleted, using the
// same config, offered and cache. Unrestricted optimized Global starts from the baseline solve's
// final state and repairs it after deleting the edge, using one state copy per worker. The result
// is the exact maximum flow at minimum cost. Global with a support solves every deletion afresh;
// Pinned replays the same program under each deletion (a fixed-route failure question).
// CostOnly skips a deletion when no roster source's reachable part of the static shortest-path
// support contains the edge: no route uses it, so loss is 0. All other deletions use a fresh
// ledger. A failed deletion fails the call (no partial result), raising the exception of the
// lowest failing candidate index.
// Worker and cache rules are as in batch_capacity. Clear each worker's private cache per
// candidate because deletion Eligibility uids are fresh and cannot reuse entries.
// Only CostOnly reads the cache. A deletion that needs placement with a frozen cache raises
// DagCacheMiss for its fresh uid.
[[nodiscard]] SensitivityResult sensitivity(const Graph& g, std::span<const NodeId> sources,
                                            std::span<const NodeId> destinations, const QueryOptions& opts,
                                            std::int32_t workers);

} // namespace netgraph::core
