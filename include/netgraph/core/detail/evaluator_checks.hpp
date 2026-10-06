/* The evaluators' argument checks and worker default (internal), shared with the Python
   binding so that both report the same first error (APPROACH 4 and 11). */
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>

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
// else 1. An explicit workers >= 1 is always honoured; the environment supplies defaults only
// (APPROACH 11).
std::int32_t default_workers(bool sensitivity) noexcept;

} // namespace netgraph::core::detail
