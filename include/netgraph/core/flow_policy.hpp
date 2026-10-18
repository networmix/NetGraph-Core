/* The policy engine FlowPolicy. */
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

// CapacityAware creations per call (stage 1 and accepted growth, across sources): reaching it
// stops creation, and a source whose offer is not met then ends with LimitFlows. An engine
// constant.
inline constexpr std::int64_t kMaxFlowsPerCall = 10000;

// The diminishing-returns guard of CapacityAware growth: stop when at least `window` growth
// steps were accepted in the call and the last `window` gains sum below `fraction` times the
// source's offer (finite offers only). window >= 1, fraction in (0, 1).
struct DiminishingReturns {
  std::int32_t window {1};
  double fraction {0.5};
};

// Each routing mode has only the fields it uses.
// CostOnly: hop-by-hop routing on each roster source's static shortest routes (all shortest
//   paths, with no directed cycles, edges of capacity 0 included, or the canonical single
//   path), prepared once per epoch and never refreshed, grown or reoptimised. The sources of
//   the demand are pooled: a router shared by several sources adds their traffic and splits it
//   once. EqualBalanced and Weighted forward the equal origin; MaxFlow is a controller
//   restricted to the union of those routes.
// CapacityAware: flow_count flows per source found by shortest-path searches that require
//   enough unused capacity on each edge. Grow up to max_flow_count (absent: no ceiling), subject
//   to min_flow_bandwidth (units) and max_path_cost. Rerouting is optional; `diminishing`
//   stops growth when recent gains are too small (absent: off). Its sources are not pooled.
// Global: the ideal controller, a minimum-cost maximum flow over the whole eligible graph of
//   what the policy already delivers plus this call's request (placed_demand() + volume: what an
//   earlier call could not place is not requested again).
// Pinned: one flow per usable supplied route bundle (a walk survives the Eligibility, whatever
//   the capacities), in supply order; single endpoints only. Convert a Path with to_dag().
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

// One demand's settings. A field that does not apply to the routing mode stays absent; a policy
// or an evaluator checks the combination when it uses the configuration (the fields remain
// assignable), in this order: the enumeration values, then the model's row of the matrix below,
// its columns left to right (placement, origin, admission, selection), then the DagCache, then
// CapacityAware's own parameters, then the pooled single-path rule.
//
// | routing       | placement               | origin (default / allowed) | admission         | selection |
// | CostOnly      | EqualBalanced, Weighted | Equal / Equal              | Lossless, Lossy   | all four fields |
// | CostOnly      | MaxFlow                 | Optimized / both           | Lossless          | all four fields |
// | Global        | none                    | Optimized / both           | Lossless          | none |
// | CapacityAware | required (all three)    | Equal / Equal              | Lossless          | all four fields |
// | Pinned        | required (all three)    | none                       | Lossless; Lossy   | priority (Weighted) |
//
// An absent origin resolves to the default and an absent selection to PathSelection() in the
// policy's snapshot (config()). Only CostOnly reads a DagCache; another model rejects one. Weighted takes its weights from
// selection.priority (Capacity or Residual; Headroom is rejected). Pooled single-path rule:
// CostOnly with more than one input source (normalized, before any Eligibility) rejects
// multipath = false together with PreferHigherPriority, whose bottleneck choice would differ per
// source at a shared router.
struct FlowPolicyConfig {
  Routing routing {CostOnly {}};
  std::optional<FlowPlacement> placement {};
  Admission admission {Admission::Lossless};
  std::optional<PathSelection> selection {};
  std::optional<Origin> origin {};
};

// Optional outputs of a placement (absent unless requested): edge_flows (this policy's own
// allocation per edge after the call), residual and bottleneck (the ledger's after the call),
// and min_cut (Global with the optimized origin; see PlacementResult::min_cut).
struct OutputOptions {
  bool edge_flows {false};
  bool residual {false};
  bool bottleneck {false};
  bool min_cut {false};
};

// The call's outcome, the first that applies: Complete for a zero request; NoRoute for a
// positive request with an empty roster (no eligible source statically reaches an allowed
// destination); NoActiveWeight when routes exist but none has a positive weight (Weighted over
// edges of capacity 0 only); a Lossy finite request: Complete (delivery and losses account for
// the whole offer, which does not mean nothing was lost); otherwise Complete when the whole
// request was delivered, else Exhausted (an infinite request is Exhausted). A controller whose
// maximum is 0 over existing routes is Exhausted, not NoRoute.
// CapacityAware alone reduces per-source outcomes instead: a source is served when its offer is
// 0 or delivered; an unserved source's outcome is LimitFlows (creation or growth hit
// kMaxFlowsPerCall, or growth hit max_flow_count > flow_count), else LimitDiminishing (the
// diminishing guard), else NoRoute (no flows), else Exhausted; the call takes the highest
// unserved outcome in that order, Complete when every source is served.
enum class Termination { Complete, LimitFlows, LimitDiminishing, NoRoute, Exhausted, NoActiveWeight };

// One normalized input source's forwarding in the epoch (zero rows included: a masked or
// unreachable source): offered (Lossy only: its injections), admitted (Lossless: what it
// delivers; Lossy: offered), delivered, source_loss (offered shares lost at the source for lack of
// a usable route) and dropped (its losses on links). Lossy: offered = delivered + source_loss +
// dropped.
struct SourceFlow {
  NodeId source {0};
  std::optional<Counter128> offered {};
  Counter128 admitted {};
  Amount delivered {0};
  Counter128 source_loss {};
  Counter128 dropped {};
  friend bool operator==(const SourceFlow&, const SourceFlow&) = default;
};

// The result of one call: this call's values, then the demand's state after it (totals of the
// epoch), then the optional outputs.
// This call: placed (units delivered by the call), remaining = request - placed (nullopt for
// +inf), delivery_delta (new minus old delivery: placed for a placement), call_drops (this call's
// drops per edge, ascending EdgeId), call_source_loss, termination and the change flags.
// routes_changed means a group of flows was created (even at zero delivery: CostOnly prepares
// every roster source's route at the epoch's first positive call), a retained flow was rerouted,
// growth added a flow, or Global changed a source's used-edge set; weights_changed means a
// rerouted Weighted flow captured its weights again. A candidate equal to the current route is not
// a change.
// State: requested (the sum of the requests; nullopt after +inf), offered (Lossy only: the
// traffic actually sent; the first-hop injection of a Lossy +inf), admitted (Lossless: the total
// carried; Lossy: offered), delivered, unadmitted (Lossless: requested - delivered, nullopt after
// +inf; Lossy: 0), source_loss and drops (per edge). Lossy: offered = source_loss + sum(drops) +
// delivered. by_source has one row per normalized input source in ascending NodeId; the rows sum
// to the totals except for Lossy with an empty roster, where the whole offer is source loss of
// the demand and every row is zero. total_cost (controllers: Global and CostOnly MaxFlow) is the
// exact sum of cost * flow over the policy's allocation (below 2^124).
// cost_distribution: Global with the optimized origin only, the cumulative marginal cost curve:
// (cost, units) pairs in ascending cost, where delivering `units` more units raised the minimum
// total cost of the flow by `cost` per unit (more delivery can reroute flow already placed, so the
// pairs are not the costs of the paths in the final allocation); empty for an empty roster, and
// absent before the epoch's first positive call.
struct PlacementResult {
  Amount placed {0};
  std::optional<Amount> remaining {};
  Amount delivery_delta {0};
  std::vector<std::pair<EdgeId, Amount>> call_drops {};
  Amount call_source_loss {0};
  Termination termination {Termination::Complete};
  bool routes_changed {false};
  bool weights_changed {false};
  std::optional<Counter128> requested {};
  std::optional<Counter128> offered {};
  Counter128 admitted {};
  Amount delivered {0};
  std::optional<Counter128> unadmitted {};
  Counter128 source_loss {};
  std::vector<std::pair<EdgeId, Counter128>> drops {};
  std::vector<SourceFlow> by_source {};
  std::optional<Counter128> total_cost {};
  std::optional<std::vector<std::pair<Cost, Amount>>> cost_distribution {};
  std::optional<std::vector<Amount>> edge_flows {};
  std::optional<std::vector<Amount>> residual {};
  std::optional<std::vector<EdgeId>> bottleneck {};
  // Requested by OutputOptions::min_cut; Global with the optimized origin only, when the solve
  // reached the maximum flow: the edges leaving the set of nodes the roster sources reach through
  // residual arcs (free capacity plus this policy's own load, reverse arcs that cancel flow
  // included), whatever their free capacity (an allowed edge of capacity 0 is listed; an excluded
  // edge is not an edge here), ascending EdgeId. That is the minimum cut closest to the sources; an
  // empty roster returns an empty cut. nullopt otherwise (another model, the equal origin, a call
  // below the maximum, a zero call).
  std::optional<std::vector<EdgeId>> min_cut {};
};

// A policy binds to its ledger at construction. place() validates, then computes on working
// copies and commits all changes together, or none. CostOnly prepares its routes at the epoch's
// first positive call (every roster source's shortest routes on the ledger's residual at that
// call's entry) and keeps them for the epoch; EqualBalanced and Weighted forward the equal origin
// over the pooled routes: Lossless admits the largest total a in [admitted, admitted + volume]
// whose no-drop forwarding fits (a refused unit is not offered again), Lossy forwards the new
// round-robin injections on the residual at call entry (earlier drops are not offered again).
// The controllers (Global, CostOnly MaxFlow) re-solve admitted + volume from zero on the
// residual plus their own allocation; the equal origin admits the largest total whose round-robin
// supplies are all carried. CapacityAware and Pinned run per source (stages of creation,
// reoptimisation, amounts and growth). remove() releases every policy flow, retires their ids and
// starts a new epoch: the period since construction or the last remove() over which requests and
// flows accumulate. Destroying a policy without remove() leaves its records in the ledger as
// orphans.
//
// Validation (std::invalid_argument, before anything binds), in this order: the endpoint sets
// (non-empty, ids in range, disjoint), a foreign Eligibility or DagCache, the configuration
// (FlowPolicyConfig's order), then Pinned's bundles (a non-empty list of this graph without null
// or duplicate bundles, every bundle with a src -> dst walk before pruning). place() validates
// the request ([0, 2^62) units or +inf; a finite request after +inf in the epoch is rejected,
// also 0) before any change. A solve of the flow core (Global, CostOnly MaxFlow, a MaxFlow flow)
// whose network exceeds the core's int32 numbering throws std::invalid_argument, and the call
// changes nothing: more than INT32_MAX - 2 nodes, or more than 2^30 - 1 edges, sources and
// destinations together.
//
// Lifetime: the ledger, the Eligibility and the DagCache outlive the policy (the policy holds
// pointers to them). The roster (the eligible sources that statically reach an allowed
// destination, over allowed edges whatever their capacity) and, for Pinned, the usable bundles
// are fixed at construction.
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
  // The configuration as validated at construction: a copy whose absent origin and selection are
  // resolved to the model's defaults (FlowPolicyConfig).
  [[nodiscard]] const FlowPolicyConfig& config() const noexcept;

private:
  friend struct detail::PolicyAccess;
  // Construction inside another public call (an evaluator): a supplied frozen cache's hits are
  // recorded in that call's tally.
  FlowPolicy(FlowLedger& ledger, std::span<const NodeId> sources, std::span<const NodeId> destinations,
             FlowPolicyConfig config, const Eligibility* elig, DagCache* cache, detail::FrozenTally& tally);
  std::unique_ptr<detail::PolicyEngine> engine_;
};

} // namespace netgraph::core
