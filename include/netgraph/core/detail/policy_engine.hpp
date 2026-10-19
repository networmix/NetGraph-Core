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
#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/program.hpp"

namespace netgraph::core::detail {

struct DeletionCapture;

// A live flow of the policy between calls. Its amounts are its ledger record (delivered; edges;
// Lossy offered and drops; T1's holds). The flows of CostOnly and Global are one per roster
// source (CostOnly's sources share the epoch's program); a T1, T2 or Pinned policy's are one per
// program entry (`entry`).
struct LiveFlow {
  std::int64_t id {-1};
  NodeId src {0};
  std::shared_ptr<const PredDAG> bundle {};
  Cost cost {kUnreachable};
  std::int32_t entry {-1};
};

// The state accumulated since construction or the last remove(), replace() or forward() (the
// policy epoch): a value, so a full-offer call builds its own and publishes it whole.
struct Epoch {
  // The live flows: per roster index (CostOnly and Global: one per source), or per program entry
  // (T1, T2, Pinned).
  std::vector<std::vector<LiveFlow>> cohorts {};
  // The installed program (null before the first preparation) and, for T1 and T2, its setup summary
  // (shared, never changed: a policy without tunnels carries a null pointer, not a TunnelPlan).
  std::shared_ptr<const ForwardingProgram> program {};
  std::shared_ptr<const TunnelPlan> plan {};
  // CostOnly MaxFlow: the union of the roster's routes, the controller's support (ascending).
  std::vector<EdgeId> support {};
  // Unrestricted Global with the optimized origin: the cumulative marginal curve of the last solve
  // (absent before one, and after a forward).
  std::optional<std::vector<std::pair<Cost, Amount>>> curve {};
  // Pinned, Lossy: per view head, the shares it lost itself (black-holed, no surviving entry).
  std::vector<Counter128> head_lost {};
  // The by_tunnel rows of the state (a zero call reports them).
  std::vector<TunnelFlow> by_tunnel {};
  Counter128 requested {};
  Counter128 offered {};  // Lossy: the traffic sent (with an empty roster: lost at the sources)
  bool infinite {false};
  bool prepared {false};  // the epoch's first positive call (or a full-offer call) has run
};

struct PolicyEngine {
  PolicyEngine(FlowLedger& ledger, FlowPolicyConfig cfg) : fg(ledger), config(std::move(cfg)) {}

  FlowLedger& fg;
  std::shared_ptr<const Graph> graph {};
  FlowPolicyConfig config;       // validated, origin, selection and weights resolved
  PathSelection selection {};    // *config.selection, or the default for Global and Pinned (unread)
  Origin origin {Origin::Equal};  // *config.origin, or Equal for Pinned (unread)
  const Eligibility* elig {nullptr};
  std::vector<NodeId> sources {};       // the normalised sets
  std::vector<NodeId> destinations {};
  std::vector<NodeId> dests_ok {};      // the allowed destinations (the records' destination set)
  std::vector<NodeId> roster {};        // fixed at construction, ascending
  Cost c0 {0};                          // T1 with a Stretch ceiling: the smallest static shortest cost
                                        // over the roster heads (the pool baseline)
  std::shared_ptr<const PredDAG> support {};  // CostOnly: the full static support
  // Global with a support: the trimmed support (ascending), shared with its programs' replay basis.
  std::shared_ptr<const std::vector<EdgeId>> global_support {};
  // Unrestricted Global, T1 and T2: the mask's excluded nodes and edges (ascending), their programs'
  // replay basis.
  std::vector<NodeId> excluded_nodes {};
  std::vector<EdgeId> excluded_edges {};
  // Pinned: its program as this policy forwards it (the replay under the mask, the black holes),
  // and whether any head can place traffic (some head with a positive injection weight and a
  // surviving entry of positive weight with a usable route).
  ForwardingView view {};
  bool active {false};
  std::int64_t owner {0};
  Epoch epoch {};
  // Global: a capture request for the next solve's final state (the exact sensitivity model's
  // baseline); set and cleared by the evaluator.
  DeletionCapture* capture {nullptr};
};

// The internal access to a policy's engine (the pattern of LedgerAccess): FlowPolicy's public
// interface is the documented API; everything else reaches the state here.
struct PolicyAccess {
  [[nodiscard]] static const PolicyEngine& engine(const FlowPolicy& policy) noexcept { return *policy.engine_; }
  [[nodiscard]] static PolicyEngine& engine(FlowPolicy& policy) noexcept { return *policy.engine_; }
  // A policy constructed inside another public call, whose tally records a supplied frozen
  // cache's hits; returned as a prvalue (FlowPolicy does not move).
  [[nodiscard]] static FlowPolicy make(FlowLedger& ledger, std::span<const NodeId> sources,
                                       std::span<const NodeId> destinations, FlowPolicyConfig config,
                                       const Eligibility* elig, DagCache* cache, FrozenTally& tally) {
    return FlowPolicy(ledger, sources, destinations, std::move(config), elig, cache, tally);
  }
};

// The configuration checks of FlowPolicy, shared with the evaluators' argument checks, in
// FlowPolicyConfig's order (the enumeration values, the model's row, the DagCache: `cache` says
// whether one is supplied, the tunnel models' parameters, edge_weights' length and entries against
// `num_edges`);
// returns its selection, origin and weights with an absent one resolved (nothing else is copied).
// Then, with the normalized endpoint counts, the pooled single-path rule (check_source_count: an
// absent selection is PathSelection(), multipath absent, so it reads the configuration resolved or
// not alike) and T2's several-sources rule, and with the roster size the tunnel limit
// (check_tunnel_limit).
struct ResolvedAxes {
  std::optional<PathSelection> selection;
  std::optional<Origin> origin;
  std::optional<WeightBasis> weights;
};
[[nodiscard]] ResolvedAxes resolve_config(const FlowPolicyConfig& cfg, bool cache, std::int32_t num_edges);
void check_source_count(const FlowPolicyConfig& cfg, std::size_t num_sources, std::size_t num_destinations);
void check_tunnel_limit(const FlowPolicyConfig& cfg, std::size_t num_heads);
// forward()'s first checks, before its offer is read (the binding runs them before converting
// it): `Pinned forwards with replace`, then `forward requires an installed program`.
void check_forward_operation(const PolicyEngine& pe);
// The model's name in messages: CostOnly, Global, T1 (CapacityAware), T2 (Global with tunnels) or
// Pinned.
[[nodiscard]] const char* model_name(const Routing& r) noexcept;

} // namespace netgraph::core::detail
