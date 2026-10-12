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

// batch_capacity's validation in its two parts: the arguments (workers, foreign objects, the
// cache rule, Pinned, the configuration), then one pair (its endpoint sets, their overlap, its
// Eligibility's graph), the pairs in index order.
void check_batch_arguments(const Graph& g, const CapacityOptions& opts, std::int32_t workers);
void check_batch_pair(const Graph& g, std::span<const NodeId> sources,
                      std::span<const NodeId> destinations, const Eligibility* elig, std::size_t index);

// The evaluators' default worker count from the environment (NGRAPH_CORE_BATCH_THREADS for
// batch_capacity, NGRAPH_CORE_SENSITIVITY_THREADS for sensitivity): an integer in [1, 4096],
// else 1. An explicit workers >= 1 is always honoured; the environment supplies defaults only.
std::int32_t default_workers(bool sensitivity) noexcept;

// The evaluators inside a public call: as capacity, batch_capacity and sensitivity, with a
// supplied frozen cache's hits recorded in `tally` and nothing published. The public evaluators
// are frozen_call over these; the Python binding runs frozen_call over these and the construction
// of its result object, so the hits count only once that object exists.
PlacementResult capacity_in_call(const Graph& g, std::span<const NodeId> sources,
                                 std::span<const NodeId> destinations, const CapacityOptions& opts, FrozenTally& tally);
std::vector<PlacementResult> batch_capacity_in_call(
    const Graph& g, const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>>& pairs,
    const CapacityOptions& opts, std::span<const Eligibility* const> per_pair_elig, std::int32_t workers,
    FrozenTally& tally);
SensitivityResult sensitivity_in_call(const Graph& g, std::span<const NodeId> sources,
                                      std::span<const NodeId> destinations, const SensitivityOptions& opts,
                                      std::int32_t workers, FrozenTally& tally);

} // namespace netgraph::core::detail
