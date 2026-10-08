/* The state of a FlowPolicy (internal; APPROACH 8 and 9). The transaction lives in
   src/flow_policy.cpp; the evaluators, the private bindings and the tests read this state
   through PolicyAccess. */
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/pred_dag.hpp"

namespace netgraph::core::detail {

struct DeletionCapture;

// A live flow of the policy between calls. Its amounts are its ledger record (delivered = the
// cumulative T of a Lossless flow; edges; Lossy offered and drops); the plan is frozen at the
// flow's creation (and replaced by an adopted reroute). Global flows have no plan.
struct LiveFlow {
  std::int64_t id {-1};
  NodeId src {0};
  std::shared_ptr<const KernelPlan> plan {};
  std::shared_ptr<const PredDAG> bundle {};
  Cost cost {kUnreachable};
  std::int32_t bundle_index {-1};  // Pinned: the supplied bundle
};

// One supplied Pinned bundle prepared at construction: its usable route (static eligibility,
// capacity >= 1 unit, pruned to the src -> dst walks) as the plan of every static-membership
// kernel, its cost, and the record bundle (the supplied DAG itself when nothing was pruned).
struct PreparedPinnedBundle {
  std::shared_ptr<const KernelPlan> route_plan {};
  Cost cost {kUnreachable};
  std::shared_ptr<const PredDAG> route_dag {};
};

struct PolicyEngine {
  PolicyEngine(FlowLedger& ledger, FlowPolicyConfig cfg) : fg(ledger), config(std::move(cfg)) {}

  FlowLedger& fg;
  std::shared_ptr<const Graph> graph {};
  FlowPolicyConfig config;
  const Eligibility* elig {nullptr};
  std::vector<NodeId> sources {};       // the normalised sets
  std::vector<NodeId> destinations {};
  std::vector<NodeId> dests_ok {};      // the allowed destinations (the records' destination set)
  std::vector<NodeId> roster {};        // fixed at construction, ascending
  std::vector<Cost> baseline {};        // CapacityAware: the static shortest cost per roster source
  std::shared_ptr<const PredDAG> support {};       // CostOnly: the full static support
  std::vector<PreparedPinnedBundle> bundles {};    // Pinned
  std::vector<std::int32_t> usable {};             // Pinned: indices of usable bundles
  std::int64_t owner {0};
  // Epoch state.
  std::vector<std::vector<LiveFlow>> cohorts {};  // per roster index, in FlowId order
  std::vector<std::pair<Cost, Amount>> curve {};  // Global: the cumulative marginal curve (9.1)
  Counter128 requested {};
  bool infinite {false};
  // Global: a capture request for the next solve's final state (the exact sensitivity model's
  // baseline, APPROACH 11); set and cleared by the evaluator.
  DeletionCapture* capture {nullptr};
};

// The internal access to a policy's engine (the pattern of LedgerAccess): FlowPolicy's public
// interface is the contract's; everything else reaches the state here.
struct PolicyAccess {
  [[nodiscard]] static const PolicyEngine& engine(const FlowPolicy& policy) noexcept { return *policy.engine_; }
  [[nodiscard]] static PolicyEngine& engine(FlowPolicy& policy) noexcept { return *policy.engine_; }
  // A policy constructed inside another public call, whose tally records a supplied frozen
  // cache's hits (NOTES R27); returned as a prvalue (FlowPolicy does not move).
  [[nodiscard]] static FlowPolicy make(FlowLedger& ledger, std::span<const NodeId> sources,
                                       std::span<const NodeId> destinations, FlowPolicyConfig config,
                                       const Eligibility* elig, DagCache* cache, FrozenTally& tally) {
    return FlowPolicy(ledger, sources, destinations, std::move(config), elig, cache, tally);
  }
};

// The configuration checks of FlowPolicy (APPROACH 9; enumeration values first, rule 13),
// shared with the evaluators' argument checks.
void validate_config(const FlowPolicyConfig& cfg, std::size_t num_sources, std::size_t num_destinations);

} // namespace netgraph::core::detail
