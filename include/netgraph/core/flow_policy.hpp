/* The policy engine FlowPolicy (APPROACH 8 and 9). */
#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

namespace detail {
struct PolicyEngine;
struct PolicyAccess;
class FrozenTally;
} // namespace detail

// CapacityAware creations per call (stage 1 and accepted growth, across sources) before the
// call stops with limit_flows (APPROACH 9.2): an engine constant.
inline constexpr std::int64_t kMaxFlowsPerCall = 10000;

// The diminishing-returns guard of CapacityAware growth: stop when at least `window` growth
// steps were accepted in the call and the last `window` gains sum below `fraction` times the
// source's offer (finite offers only). window >= 1, fraction in (0, 1).
struct DiminishingReturns {
  std::int32_t window {1};
  double fraction {0.5};
};

// The routing modes (APPROACH 9.1): a field exists only in the variant it applies to.
// CostOnly: IGP forwarding, one flow per source on the static shortest-path DAG (or its
//   canonical single path), never refreshed, grown or reoptimised.
// CapacityAware: distributed CSPF tunnels: flow_count flows created by the canonical search on
//   the selection view, grown up to max_flow_count (absent: no ceiling), gated by
//   min_flow_bandwidth (units) and the max_path_cost ceiling, optionally reoptimised and guarded
//   by `diminishing` (absent: off).
// Global: the omniscient controller, the minimum-cost flow of the cumulative value over the
//   whole eligible graph (a replacement kernel; placement MaxFlow only).
// Pinned: explicit route structures, one flow per usable bundle in supply order (single
//   endpoints; a Path converts with to_dag()).
struct CostOnly {};
struct CapacityAware {
  std::int32_t flow_count {1};
  std::optional<std::int32_t> max_flow_count {};
  std::optional<CostCeiling> max_path_cost {};
  Amount min_flow_bandwidth {0};
  bool reoptimize {false};
  std::optional<DiminishingReturns> diminishing {};
};
struct Global {};
struct Pinned {
  std::vector<std::shared_ptr<const PredDAG>> bundles {};
};
using Routing = std::variant<CostOnly, CapacityAware, Global, Pinned>;

// The four outer axes are explicit (no defaults). Global and Pinned read selection.priority only
// for the Weighted kernel (membership and weights) and ignore the rest; CostOnly and
// CapacityAware read all four fields (APPROACH 9). The constructor rejects a placement,
// admission or selection field outside its enumeration (std::invalid_argument naming the field;
// rule 13); the fields stay assignable, so the policy and the evaluators check them again with
// the combination rules.
struct FlowPolicyConfig {
  FlowPolicyConfig(Routing routing_, FlowPlacement placement_, Admission admission_, PathSelection selection_);
  Routing routing;
  FlowPlacement placement;
  Admission admission;
  PathSelection selection;
};

// Optional outputs of a placement (absent unless requested): edge_flows (this policy's own
// allocation per edge after the call), residual and bottleneck (the ledger's after the call),
// and on Global the min_cut of the solve when bottleneck is requested.
struct OutputOptions {
  bool edge_flows {false};
  bool residual {false};
  bool bottleneck {false};
};

// The local outcome of a source and the call's termination (APPROACH 9.4).
enum class Termination { Complete, LimitFlows, LimitDiminishing, NoRoute, Exhausted };

// The result of one placement call: the deltas of this call (placed, remaining = request -
// placed, nullopt for +inf; drops = this call's transit drops), the termination and one local
// outcome per roster source (ascending NodeId), the two flags of the 9.4 table; Global adds the
// marginal cost curve of its cumulative solve. The optional outputs are present when requested.
struct PlacementResult {
  Amount placed {0};
  std::optional<Amount> remaining {};
  Termination termination {Termination::Complete};
  std::vector<std::pair<NodeId, Termination>> local_outcomes {};
  bool routes_changed {false};
  bool weights_changed {false};
  std::vector<std::pair<EdgeId, Amount>> drops {};
  std::optional<std::vector<std::pair<Cost, Amount>>> cost_distribution {};
  std::optional<std::vector<Amount>> edge_flows {};
  std::optional<std::vector<Amount>> residual {};
  std::optional<std::vector<EdgeId>> bottleneck {};
  // Global with bottleneck requested: the full-residual cut when the solve was maximal on its
  // view (an empty roster's maximal flow 0 has the empty cut); nullopt otherwise: not requested,
  // not Global, a solve stopped by its cap, or a zero call (nothing is solved).
  std::optional<std::vector<EdgeId>> min_cut {};
};

// The one placement engine (APPROACH 9). A policy binds to its ledger at construction and is the
// only way flows enter it: place() runs the transaction of 9.2 and commits one atomic exchange;
// remove() reverts every flow of the policy, retires the ids and opens a new epoch; destroying a
// policy without remove() leaves its records in the ledger as orphans.
//
// Validation (std::invalid_argument, before anything binds): the endpoint sets (non-empty, ids
// in range, disjoint), a foreign Eligibility or DagCache, the configuration (enumeration values;
// Lossy with MaxFlow, Global or CapacityAware; Global without MaxFlow; Weighted with Headroom;
// flow_count >= 1; max_flow_count >= flow_count; min_flow_bandwidth >= 0; the ceiling; the
// diminishing guard; Pinned: single endpoints, a non-empty bundle list of this graph without
// duplicates, every bundle with a src -> dst walk before pruning). place() validates the
// request ([0, 2^62) units or +inf; a finite request after +inf in the epoch is rejected, also
// 0) before any change.
//
// Lifetime: the ledger, the Eligibility and the DagCache outlive the policy (the policy holds
// pointers to them). The roster (the eligible sources that statically reach an allowed
// destination) and, for Pinned, the usable bundles are fixed at construction.
class FlowPolicy {
public:
  FlowPolicy(FlowLedger& ledger, std::span<const NodeId> sources, std::span<const NodeId> destinations,
             FlowPolicyConfig config, const Eligibility* elig = nullptr, DagCache* cache = nullptr);
  FlowPolicy(const FlowPolicy&) = delete;
  FlowPolicy& operator=(const FlowPolicy&) = delete;
  FlowPolicy(FlowPolicy&&) = delete;
  FlowPolicy& operator=(FlowPolicy&&) = delete;
  ~FlowPolicy();

  // volume: units, nullopt is +inf.
  PlacementResult place(std::optional<Amount> volume, const OutputOptions& outputs = {});
  void remove();
  // The sum of delivered over this policy's live records.
  [[nodiscard]] Amount placed_demand() const;
  // This policy's live records as snapshots in FlowId order.
  [[nodiscard]] std::vector<FlowRecord> records() const;
  // One of this policy's live records; a FlowId of another ledger or policy, a retired or an
  // unknown id throws std::invalid_argument.
  [[nodiscard]] FlowRecord record(const FlowId& flow_id) const;
  // Pinned: the number of bundles with a src -> dst walk after static pruning; 0 otherwise.
  [[nodiscard]] std::int32_t usable_bundles() const noexcept;
  // The configuration (a copy taken at construction).
  [[nodiscard]] const FlowPolicyConfig& config() const noexcept;

private:
  friend struct detail::PolicyAccess;
  // Construction inside another public call (an evaluator): a supplied frozen cache's hits are
  // recorded in that call's tally (NOTES R27).
  FlowPolicy(FlowLedger& ledger, std::span<const NodeId> sources, std::span<const NodeId> destinations,
             FlowPolicyConfig config, const Eligibility* elig, DagCache* cache, detail::FrozenTally& tally);
  std::unique_ptr<detail::PolicyEngine> engine_;
};

} // namespace netgraph::core
