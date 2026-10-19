/* The evaluators' argument checks, worker default and in-call forms (internal), shared with the
   Python binding so that both report the same first error and publish a frozen cache's
   statistics at the same boundary. */
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/detail/cache_access.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

// evaluate_pairs's validation in its two parts: the arguments (workers, foreign objects, the
// cache rule, Pinned, the configuration), then one pair (its endpoint sets, their overlap, its
// Eligibility's graph, the pooled single-path rule for its source count), the pairs in index
// order.
void check_batch_arguments(const Graph& graph, const EvaluationOptions& options, std::int32_t workers);
void check_batch_pair(const Graph& graph, const FlowPolicyConfig& config, std::span<const NodeId> sources,
                      std::span<const NodeId> destinations, const Eligibility* eligibility, std::size_t index);

// The evaluators' default worker count from the environment (NGRAPH_CORE_BATCH_THREADS for
// evaluate_pairs, NGRAPH_CORE_SENSITIVITY_THREADS for evaluate_edge_failures): an integer in [1, 4096],
// else 1. An explicit workers >= 1 is always honoured; the environment supplies defaults only.
std::int32_t default_workers(bool edge_failures) noexcept;

// The evaluators inside a public call: as evaluate_demand, evaluate_pairs and evaluate_edge_failures, with a
// supplied frozen cache's hits recorded in `tally` and nothing published. The public evaluators
// are frozen_call over these; the Python binding runs frozen_call over these and the construction
// of its result object, so the hits count only once that object exists.
PlacementResult evaluate_demand_in_call(const Graph& graph, std::span<const NodeId> sources,
                                 std::span<const NodeId> destinations, const EvaluationOptions& options, FrozenCacheHitTally& tally);
std::vector<PlacementResult> evaluate_pairs_in_call(
    const Graph& graph, const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>>& pairs,
    const EvaluationOptions& options, std::span<const Eligibility* const> per_pair_eligibility, std::int32_t workers,
    FrozenCacheHitTally& tally);
SensitivityResult evaluate_edge_failures_in_call(const Graph& graph, std::span<const NodeId> sources,
                                      std::span<const NodeId> destinations, const EvaluationOptions& options,
                                      std::int32_t workers, FrozenCacheHitTally& tally);

} // namespace netgraph::core::detail
