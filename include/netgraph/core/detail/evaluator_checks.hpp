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
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

#if NETGRAPH_CORE_PRIVATE_HOOKS
// Native-only controls for comparing edge-failure shortcuts with complete repairs.
// The caller owns the probe through all worker joins; observers may not use Python.
enum class SensitivityWork { FullRepair, FirstRepair, Cut };
struct SensitivityProbe {
  bool full_repair {false};
  bool cut_edges {true};
  void (*observe)(void*, EdgeId, SensitivityWork) {nullptr};
  void* context {nullptr};
};
[[nodiscard]] const SensitivityProbe* current_edge_failure_probe() noexcept;
void set_edge_failure_probe(const SensitivityProbe* probe) noexcept;
struct ScopedSensitivityProbe {
  const SensitivityProbe* previous;
  explicit ScopedSensitivityProbe(const SensitivityProbe* probe) noexcept : previous(current_edge_failure_probe()) { set_edge_failure_probe(probe); }
  ~ScopedSensitivityProbe() { set_edge_failure_probe(previous); }
  ScopedSensitivityProbe(const ScopedSensitivityProbe&) = delete;
  ScopedSensitivityProbe& operator=(const ScopedSensitivityProbe&) = delete;
};
#endif

// evaluate_pairs's validation in its two parts: the arguments (workers, foreign objects, the
// Pinned, the configuration), then one pair (its endpoint sets, their overlap, its
// Eligibility's graph, the pooled single-path rule for its source count), the pairs in index
// order. Graph searches for the valid prefix run in one parallel native phase,
// before any later pair's check error is rethrown. Their rosters are reused by workers.
void check_batch_arguments(const Graph& graph, const EvaluationOptions& options, std::int32_t workers);
void check_batch_pair(const Graph& graph, const FlowPolicyConfig& config, std::span<const NodeId> sources,
                      std::span<const NodeId> destinations, const Eligibility* eligibility, std::size_t index);
std::vector<TunnelRoster> prepare_batch_rosters(const Graph& graph,
    std::span<const std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs,
    const EvaluationOptions& options, std::span<const Eligibility* const> per_pair, std::int32_t workers,
    ArgumentWork search = nullptr);
std::optional<TunnelRoster> check_demand_arguments(const Graph& graph, std::span<const NodeId> sources,
    std::span<const NodeId> destinations, const EvaluationOptions& options, ArgumentWork search = nullptr);
std::optional<TunnelRoster> check_edge_failure_arguments(const Graph& graph, std::span<const NodeId> sources,
    std::span<const NodeId> destinations, const EvaluationOptions& options, std::int32_t workers, ArgumentWork search = nullptr);

// The evaluators' default worker count from the environment (NGRAPH_CORE_BATCH_THREADS for
// evaluate_pairs, NGRAPH_CORE_SENSITIVITY_THREADS for evaluate_edge_failures): an integer in [1, 4096],
// else 1. An explicit workers >= 1 is always honoured; the environment supplies defaults only.
std::int32_t default_workers(bool edge_failures) noexcept;

// The evaluators inside a public call, after all argument/immutable checks above and with the
// supplied cache's lease held: as evaluate_demand, evaluate_pairs and evaluate_edge_failures, with a
// supplied frozen cache's hits recorded in `tally` and nothing published. The public evaluators
// are frozen_call over these; the Python binding allocates its result shells first, then runs
// frozen_call over these and fills the shells' native values before the hits are published.
// A nonempty roster span has one checked tunnel roster per pair, for the same
// config/endpoints/effective mask. It lives through every worker join; a single
// roster pointer applies only to the matching demand or evaluate_edge_failures baseline.
PlacementResult evaluate_demand_in_call(const Graph& graph, std::span<const NodeId> sources,
                                 std::span<const NodeId> destinations, const EvaluationOptions& options, FrozenCacheHitTally& tally, const TunnelRoster* roster = nullptr);
std::vector<PlacementResult> evaluate_pairs_in_call(
    const Graph& graph, const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>>& pairs,
    const EvaluationOptions& options, std::span<const Eligibility* const> per_pair_eligibility, std::int32_t workers,
    FrozenCacheHitTally& tally, std::span<const TunnelRoster> rosters = {});
SensitivityResult evaluate_edge_failures_in_call(const Graph& graph, std::span<const NodeId> sources,
                                      std::span<const NodeId> destinations, const EvaluationOptions& options,
                                      std::int32_t workers, FrozenCacheHitTally& tally, const TunnelRoster* roster = nullptr);

} // namespace netgraph::core::detail
