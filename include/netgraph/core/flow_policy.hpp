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
// CostOnly: one flow per source on the static shortest-path DAG (all shortest paths, with no
//   directed cycles), or its canonical single path; never refreshed, grown or reoptimised.
// CapacityAware: flow_count flows found by shortest-path searches that require enough unused
//   capacity on each edge. Grow up to max_flow_count (absent: no ceiling), subject to
//   min_flow_bandwidth (units) and max_path_cost. Rerouting is optional; `diminishing`
//   stops growth when recent gains are too small (absent: off).
// Global: replace the policy's allocation with a minimum-cost flow of what it already delivers
//   plus this call's request (placed_demand() + volume: what an earlier call could not place is
//   not requested again) over the whole eligible graph (placement MaxFlow only).
// Pinned: one flow per usable supplied route bundle, in supply order; single endpoints only.
//   Convert a Path with to_dag().
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

// routing, placement, admission and selection are required; none has a default.
// Global and Pinned use selection.priority only for Weighted membership and weights, and
// ignore the other selection fields. CostOnly and CapacityAware use all four selection fields.
// Setting an ignored field is not an error. Invalid placement, admission or selection enum
// values throw std::invalid_argument naming the field. Fields remain assignable, so policies
// and evaluators check them again, then check the supported combinations.
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

// A source is served if its offer is 0, Lossless delivers the offer, or Lossy forwards the
// entire finite offer (including drops). Choose its first matching outcome:
// Complete: served.
// LimitFlows: creation/growth hit kMaxFlowsPerCall, or growth hit max_flow_count > flow_count.
// LimitDiminishing: the diminishing-returns guard stopped growth.
// NoRoute: no flows remain; for Global, maximum flow is 0 and nothing was ever placed.
// Exhausted: any other unserved case, including full flows with no growth, a full fixed-size
// group of flows, Global at maximum flow, or an infinite Lossy offer.
// Global and Pinned act as one source: every local outcome equals the call's termination.
// The call is Complete for request 0 or when all sources are served, and NoRoute for a positive
// request with an empty roster: no eligible source statically reaches an allowed destination.
// Otherwise, take the highest unserved outcome: LimitFlows > LimitDiminishing > Exhausted >
// NoRoute.
enum class Termination { Complete, LimitFlows, LimitDiminishing, NoRoute, Exhausted };

// Per-call changes: placed, remaining = request - placed (nullopt for +inf), and transit drops.
// Include the termination and one local outcome per roster source in ascending NodeId.
// Global also returns its cumulative marginal cost curve: (cost, units) pairs in ascending cost,
// where delivering `units` more units raises the minimum total cost of the flow by `cost` per
// unit. Delivering more can reroute flow already placed, so the pairs are not the costs of the
// paths in the final allocation. Optional outputs appear only if asked.
// Both change flags start false and report committed changes only. routes_changed means a
// source's group of flows was created (even at zero delivery), a retained flow was rerouted,
// growth added a flow, or Global changed a source's used-edge set. weights_changed means a
// rerouted Weighted flow captured its weights again; Global never sets it. A candidate equal
// to the current route is not a change.
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
  // Global, when bottleneck is requested: positive-view edges leaving the set reachable from
  // the real sources through residual arcs, including reverse arcs that cancel flow. A maximal
  // solve returns this minimum cut; an empty roster returns an empty cut for its maximal flow
  // of 0. nullopt if not requested, not Global, stopped by the cap, or a zero call (no solve).
  std::optional<std::vector<EdgeId>> min_cut {};
};

// A policy binds to its ledger at construction. place() validates, then Global re-solves what
// it already delivers plus the request; other modes create, reroute, fill and grow flows per roster source using
// one working copy of unused capacity. All changes commit together, or none do.
// remove() releases every policy flow, retires their ids and starts a new epoch: the period
// since construction or the last remove() over which requests and flows accumulate. Destroying
// a policy without remove() leaves its records in the ledger as orphans.
//
// Validation (std::invalid_argument, before anything binds): the endpoint sets (non-empty, ids
// in range, disjoint), a foreign Eligibility or DagCache, the configuration (enumeration values;
// Lossy with MaxFlow, Global or CapacityAware; Global without MaxFlow; Weighted with Headroom;
// flow_count >= 1; max_flow_count >= flow_count; min_flow_bandwidth in [0, 2^62) units; the
// ceiling; the diminishing guard; Pinned: single endpoints, a non-empty bundle list of this graph
// without null or duplicate bundles, every bundle with a src -> dst walk before pruning). The order between
// these three groups of checks is not specified: when defects of two groups meet, the one
// reported may differ between this constructor and a language binding. place() validates the request ([0, 2^62) units or +inf; a
// finite request after +inf in the epoch is rejected, also 0) before any change. A solve of the
// flow core (Global, or a MaxFlow flow) whose network exceeds the core's int32 numbering throws
// std::invalid_argument, and the call changes nothing: more than INT32_MAX - 2 nodes, or more than
// 2^30 - 1 edges, sources and destinations together (Global counts the graph's nodes, the
// eligible edges, the roster and the eligible destinations; a MaxFlow flow the nodes of its
// member edges and endpoints, its member edges, its source and its destinations).
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
  // recorded in that call's tally.
  FlowPolicy(FlowLedger& ledger, std::span<const NodeId> sources, std::span<const NodeId> destinations,
             FlowPolicyConfig config, const Eligibility* elig, DagCache* cache, detail::FrozenTally& tally);
  std::unique_ptr<detail::PolicyEngine> engine_;
};

} // namespace netgraph::core
