/* The state of a FlowPolicy (internal). The transaction lives in
   src/flow_policy.cpp; the evaluators, the private bindings and the tests read this state
   through PolicyAccess. */
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/detail/forwarding.hpp"
#include "netgraph/core/detail/route_kernels.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/forwarding_program.hpp"

namespace netgraph::core::detail {

struct DeletionCapture;
// A detached argument search may be run by the binding; the callback and its
// context are native stack values, invoked synchronously before any lock.
using ArgumentWork = void (*)(void (*)(void*), void*);
struct TunnelRoster {
  std::vector<NodeId> heads;
  Cost c0 {0};
};
enum class PolicyOperation { Place, Evaluate, Replace, Forward };

// A live flow of the policy between calls. Its amounts are its ledger record (delivered; edge flows;
// Lossy offered and edge drops; T1's holds). The flows of ShortestPaths and CentralizedTE are one per roster
// source (ShortestPaths's sources share the epoch's program); a T1, T2 or Pinned policy's are one per
// program entry (`entry_id`).
struct LiveFlow {
  std::int64_t record_counter {-1};
  NodeId source {0};
  std::shared_ptr<const PredecessorDAG> route {};
  Cost min_route_cost {kUnreachable};
  std::int32_t entry_id {-1};
};

// The state accumulated since construction or the last remove(), replace() or forward() (the
// policy epoch): a value, so a full-offer call builds its own and publishes it whole.
struct PolicyEpoch {
  // The live flows: per roster index (ShortestPaths and CentralizedTE: one per source), or per program entry
  // (T1, T2, Pinned).
  std::vector<std::vector<LiveFlow>> record_groups {};
  // The installed program (null before the first preparation) and, for T1 and T2, its setup summary
  // (shared, never changed: a policy without tunnels carries a null pointer, not a TunnelPlan).
  std::shared_ptr<const ForwardingProgram> program {};
  std::shared_ptr<const TunnelPlan> plan {};
  // ShortestPaths Optimized: the union of the roster's routes, the controller's support (ascending).
  std::vector<EdgeId> support {};
  // Unrestricted CentralizedTE with the optimized source allocation: the cumulative marginal curve of the last solve
  // (absent before one, and after a forward).
  std::optional<std::vector<std::pair<Cost, Amount>>> marginal_cost_curve {};
  // Pinned, Lossy: per view head, the shares it lost itself (black-holed, no surviving entry).
  std::vector<Counter128> source_drops_by_source {};
  // The by_tunnel rows of the state (a zero call reports them).
  std::vector<TunnelFlow> by_tunnel {};
  Counter128 requested {};
  Counter128 offered {};  // Lossy: the traffic sent (with an empty roster: lost at the sources)
  bool infinite {false};
  bool prepared {false};  // the epoch's first positive call (or a full-offer call) has run
};

struct PolicyEngine {
  PolicyEngine(FlowLedger& ledger, FlowPolicyConfig&& config) : ledger(ledger), config(std::move(config)) {}

  FlowLedger& ledger;
  std::shared_ptr<const Graph> graph {};
  FlowPolicyConfig config;       // validated, source allocation and weight basis resolved
  // The model's route choice: ShortestPaths' fields, DistributedTE's tie rule and parallel_edges as
  // a single path, the default for CentralizedTE and Pinned (unread).
  RouteChoice route_choice {};
  SourceAllocation source_allocation {SourceAllocation::Equal};  // source_allocation_of(routing), or Equal for Pinned (unread)
  const Eligibility* eligibility {nullptr};
  // T2 on a support: the policy's mask with every edge outside the support excluded too (the
  // one allowed set of every T2 step); `eligibility` points to it.
  std::unique_ptr<const Eligibility> support_mask {};
  std::vector<NodeId> sources {};       // the normalised sets
  std::vector<NodeId> destinations {};
  std::vector<NodeId> allowed_destinations {};      // the allowed destinations (the records' destination set)
  std::vector<NodeId> reachable_sources {};        // fixed at construction, ascending
  Cost c0 {0};                          // T1 with a CostStretch ceiling: the smallest static shortest cost
                                        // over the roster heads (the pool baseline)
  std::shared_ptr<const PredecessorDAG> support {};  // ShortestPaths: the full static support
  // CentralizedTE with a support: the trimmed support (ascending), shared with its programs' replay basis.
  std::shared_ptr<const std::vector<EdgeId>> centralized_te_support {};
  // Unrestricted CentralizedTE, T1 and T2: the mask's excluded nodes and edges (ascending), their programs'
  // replay basis.
  std::vector<NodeId> excluded_nodes {};
  std::vector<EdgeId> excluded_edges {};
  // Pinned: its program as this policy forwards it (the replay under the mask, the black holes),
  // and whether any head can place traffic (some head with a positive injection weight and a
  // surviving entry of positive weight with a usable route).
  ForwardingView view {};
  bool can_forward {false};
  std::int64_t policy_id {0};
  PolicyEpoch epoch {};
  // CentralizedTE: a capture request for the next solve's final state (the exact sensitivity model's
  // baseline); set and cleared by the evaluator.
  DeletionCapture* capture {nullptr};
};

// The internal access to a policy's engine (the pattern of LedgerAccess): FlowPolicy's public
// interface is the documented API; everything else reaches the state here.
struct PolicyAccess {
  [[nodiscard]] static bool constructed(const FlowPolicy& policy) noexcept { return policy.ready_.load(std::memory_order_acquire); }
  [[nodiscard]] static const PolicyEngine& engine(const FlowPolicy& policy) noexcept { return *policy.engine_; }
  [[nodiscard]] static PolicyEngine& engine(FlowPolicy& policy) noexcept { return *policy.engine_; }
  // Argument-only construction work, performed before a cache lease or ledger gate.
  [[nodiscard]] static std::unique_ptr<PolicyEngine> prepare(FlowLedger& ledger, std::span<const NodeId> sources,
      std::span<const NodeId> destinations, FlowPolicyConfig config, const Eligibility* eligibility, ShortestPathCache* cache,
      const TunnelRoster* roster = nullptr, ArgumentWork search = nullptr);
  // An unbound, destructible shell; installation is the only publishing step.
  [[nodiscard]] static std::unique_ptr<FlowPolicy> shell();
  static void construct(FlowPolicy& policy, std::unique_ptr<PolicyEngine> prepared, ShortestPathCache* cache);
  // A policy constructed inside another public call, whose tally records a supplied frozen
  // cache's hits; returned as a prvalue (FlowPolicy does not move).
  [[nodiscard]] static FlowPolicy make(FlowLedger& ledger, std::span<const NodeId> sources,
                                       std::span<const NodeId> destinations, FlowPolicyConfig config,
                                       const Eligibility* eligibility, ShortestPathCache* cache, FrozenCacheHitTally& tally, const TunnelRoster* roster = nullptr) {
    return FlowPolicy(ledger, sources, destinations, std::move(config), eligibility, cache, tally, roster);
  }
};

// The configuration checks of FlowPolicy, shared with the evaluators' argument checks, in
// FlowPolicyConfig's order (flow_policy.hpp): the enumeration values, the model matrix with a
// supplied ShortestPathCache last (`cache` says whether one is supplied), the tunnel models'
// parameters and edge_weights' length and entries against `num_edges`. With `num_sources` (the
// normalized input sources) the two source-count rules run at their places in the matrix (T2's
// several-sources rule first, the pooled single-path rule after the ineffective tie rule);
// evaluate_pairs checks them per pair with check_source_count. Then, with the roster size, the
// tunnel limit (check_tunnel_limit).
void check_config(const FlowPolicyConfig& config, bool cache, std::int32_t num_edges, std::optional<std::size_t> num_sources);
void check_source_count(const FlowPolicyConfig& config, std::size_t num_sources);
// The policy's snapshot: ShortestPaths' absent source_allocation and an absent weight_basis of
// Weighted resolved.
[[nodiscard]] FlowPolicyConfig resolved_config(FlowPolicyConfig config);
// The resolved source allocation of a computing model; none for Pinned (its program has one).
[[nodiscard]] std::optional<SourceAllocation> source_allocation_of(const Routing& r) noexcept;
void check_tunnel_limit(const FlowPolicyConfig& cfg, std::size_t num_heads);
// Construction's immutable checks, shared with evaluator preflight. Endpoints are
// normalized here before graph ownership and model checks.
std::optional<TunnelRoster> check_policy_construction(const Graph& graph, std::span<const NodeId> sources,
    std::span<const NodeId> destinations, const FlowPolicyConfig& config,
    const Eligibility* eligibility, const ShortestPathCache* cache, ArgumentWork search = nullptr);
// place()'s first check, before its offer is read (the binding runs it before converting it): T1
// and T2 raise `<model> does not support incremental placement`.
void check_place_operation(const PolicyEngine& pe);
// forward()'s first check, before its offer is read (the binding runs it before converting
// it): `Pinned forwards with replace`. The installed program is checked under the
// gate, after every argument check.
void check_forward_operation(const PolicyEngine& pe);
void check_policy_offer(const PolicyEngine& pe, PolicyOperation operation, std::optional<Amount> offer);
void check_policy_design_model(const PolicyEngine& pe);
// All argument/immutable checks of a scalar call. No epoch or ledger state is read.
void check_policy_arguments(const PolicyEngine& pe, PolicyOperation operation, std::optional<Amount> offer,
    const std::optional<std::optional<Amount>>& design = {}, const Eligibility* black_holes = nullptr);
void check_list_arguments(FlowLedger& ledger, std::span<const Action> actions, BatchMode mode, const OutputOptions& outputs);
// The model's name in messages: ShortestPaths, CentralizedTE, T1 (DistributedTE), T2 (CentralizedTE with tunnel_count) or
// Pinned.
[[nodiscard]] const char* model_name(const Routing& r) noexcept;
// The model's class name, T2 included (CentralizedTE): the cache rule and the evaluators' texts.
[[nodiscard]] const char* model_class(const Routing& r) noexcept;

} // namespace netgraph::core::detail
