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
#include "netgraph/core/program.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

namespace detail {
struct PolicyEngine;
struct PolicyAccess;
class FrozenTally;
} // namespace detail

// The tunnel count of the tunnel models, per head: `count` in [1, 10000] (checked where a policy
// or an evaluator uses it); `up_to` (T2 only; T1 rejects it): the smallest sufficient count up to
// `count`.
struct Tunnels {
  std::int64_t count {1};
  bool up_to {false};
  friend bool operator==(const Tunnels&, const Tunnels&) = default;
};

// Each routing mode has only the fields it uses.
// CostOnly: hop-by-hop routing on each roster source's static shortest routes (all shortest
//   paths, with no directed cycles, edges of capacity 0 included, or the canonical single
//   path), prepared once per epoch and fixed for it. The sources of the demand are pooled: a
//   router shared by several sources adds their traffic and splits it once. EqualBalanced and Weighted forward the equal origin; MaxFlow is a controller
//   restricted to the union of those routes. Its program: one split route per roster source,
//   pooled into one split context (EqualBalanced, Weighted), or the controller's paths (MaxFlow).
// CapacityAware (T1): every roster head sets up its own tunnels.count strict tunnels by CSPF, no
//   controller. A design volume D is split RR(D, H, h) over the heads and each head's share
//   RR(d, N, j) over its tunnels; in round order (round j, heads ascending: index j * H + i) each
//   tunnel keeps the allowed links whose unreserved capacity U (the ledger's, this policy's own
//   holds released, minus this setup's earlier up tunnels) is at least max(request,
//   min_unreserved_bandwidth) (a request of 0 with a minimum of 0 has no filter), collapses
//   parallel links (multi_edge = false), and takes the single path in the selection's order on
//   static costs (priorities read U: Residual U, Headroom U / C, Capacity C); it is Up, holding its
//   request on every link, unless no path passes the filter (NoBandwidth) or the cheapest exceeds
//   max_path_cost (CostCeiling): Absolute inclusive, Stretch(f) times c0, the smallest static
//   shortest cost over the roster heads under the mask, zero-capacity links included. Forwarding
//   is the program's: round robin over the heads, then over a head's up tunnels; a head without
//   one loses its share. Its program: those tunnels, weight 1 up and 0 down.
// Global: the ideal controller, a minimum-cost maximum flow over the whole eligible graph, or
//   exactly over `support` (a DAG of the policy's graph, trimmed to the walks from the sources to
//   the allowed destinations under the Eligibility), of what the policy already delivers plus this
//   call's request (placed_demand() + volume: what an earlier call could not place is not
//   requested again). Its program: the solved flow cut into paths per source, each path a tunnel
//   weighted by its amount.
// Pinned: forwarding over a given program, computing nothing; the program's endpoints are the
//   policy's. Its Eligibility repairs failures (a tunnel touching an excluded element is Failed, a
//   split route loses its excluded members, the roster is rebuilt from the program's basis); the
//   nodes and edges `black_holes` excludes are black holes (stale forwarding: the program is
//   forwarded unchanged and traffic reaching one is dropped there). One record per program entry.
struct CostOnly {};
struct CapacityAware {
  Tunnels tunnels {};
  std::optional<CostCeiling> max_path_cost {};
  Amount min_unreserved_bandwidth {0};
};
struct Global {
  std::shared_ptr<const PredDAG> support {};
};
struct Pinned {
  std::shared_ptr<const ForwardingProgram> program {};
  std::shared_ptr<const Eligibility> black_holes {};
};
using Routing = std::variant<CostOnly, CapacityAware, Global, Pinned>;

// One demand's settings. A field that does not apply to the routing mode stays absent; a policy
// or an evaluator checks the combination when it uses the configuration (the fields remain
// assignable), in this order: the enumeration values, then the model's row of the matrix below
// (T1: an "up to" count first, it is not T1's row), its columns left to right (placement, origin,
// admission, selection, weights with edge_weights), then the DagCache, then T1's own parameters
// (the tunnel count, max_path_cost, min_unreserved_bandwidth), then edge_weights' length and
// entries, then the roster and count rules (the pooled single-path rule, T1's tunnel limit).
//
// | routing           | placement            | origin (default / allowed) | admission       | selection                          | weights |
// | CostOnly          | EqualBalanced        | Equal / Equal              | Lossless, Lossy | all four fields                    | none |
// | CostOnly          | Weighted             | Equal / Equal              | Lossless, Lossy | all four fields                    | Capacity (default), Explicit, ResidualSnapshot |
// | CostOnly          | MaxFlow              | Optimized / both           | Lossless        | all four fields                    | none |
// | Global            | none                 | Optimized / both           | Lossless        | none                               | none |
// | T1: CapacityAware | none                 | Equal / Equal              | Lossless, Lossy | tie_break, priority, multi_edge    | none |
// | Pinned            | none                 | none (the program has it)  | Lossless, Lossy | none                               | none |
//
// An absent origin resolves to the default and an absent selection to PathSelection() in the
// policy's snapshot (config()); an absent weights of CostOnly Weighted resolves to Capacity. Only
// CostOnly reads a DagCache; another model rejects one. edge_weights (one integer per edge in
// [0, 2^62), shared by every configuration copied from this one) is given exactly when weights is
// Explicit. T1 rejects multipath = true (a tunnel is one path; absent and false are accepted).
// Pooled single-path rule: CostOnly with more than one input source (normalized, before any
// Eligibility) rejects multipath = false together with PreferHigherPriority, whose bottleneck
// choice would differ per source at a shared router. T1's tunnel limit: roster heads times the
// tunnel count at most 2^30 - 1.
struct FlowPolicyConfig {
  Routing routing {CostOnly {}};
  std::optional<FlowPlacement> placement {};
  Admission admission {Admission::Lossless};
  std::optional<PathSelection> selection {};
  std::optional<Origin> origin {};
  std::optional<WeightBasis> weights {};
  std::shared_ptr<const std::vector<Amount>> edge_weights {};
};

// Optional outputs of a placement (absent unless requested): edge_flows (this policy's own
// allocation per edge after the call), residual and bottleneck (the ledger's after the call),
// min_cut (Global with the optimized origin; see PlacementResult::min_cut) and
// reservation_edge_flows (this policy's strict holds per edge after the call, zeros for a model
// without holds).
struct OutputOptions {
  bool edge_flows {false};
  bool residual {false};
  bool bottleneck {false};
  bool min_cut {false};
  bool reservation_edge_flows {false};
};

// The call's outcome, the first that applies: Complete for a zero request; NoRoute for a
// positive request with an empty roster (no eligible source statically reaches an allowed
// destination); NoActiveWeight when the program has an entry with a path (an up tunnel, a
// controller path or a split route) but no roster source has both a positive injection weight and
// an entry of positive weight with a usable route (Weighted over edges of capacity 0 only, all
// entry weights 0, every remaining captured weight 0, every up tunnel failed); a program whose
// tunnels are all down is not this case. Then a Lossy finite request: Complete (delivery and
// losses account for the whole offer, which does not mean nothing was lost); otherwise Complete
// when the whole request was delivered, else Exhausted (an infinite request is Exhausted). A
// controller whose maximum is 0 over existing routes is Exhausted, not NoRoute.
enum class Termination { Complete, NoRoute, Exhausted, NoActiveWeight };

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

// One program tunnel's forwarding in the epoch: its entry identifier and head, its status in this
// forwarding (Up, Failed, Blackholed), assigned (its share of the injections: Lossy the shares sent
// to it, Lossless what it carries), delivered and dropped (its losses on links).
struct TunnelFlow {
  std::int32_t index {0};
  NodeId head {0};
  TunnelStatus status {TunnelStatus::Up};
  Counter128 assigned {};
  Amount delivered {0};
  Counter128 dropped {};
  friend bool operator==(const TunnelFlow&, const TunnelFlow&) = default;
};

// The tunnel models' setup summary (PlacementResult::plan): count (the tunnels per head),
// design_demand (nullopt: unlimited), head_design ((head, RR(D, H, i)) for every roster head in
// ascending NodeId; empty for an unlimited design) and reserved_total (T1's strict holds, counted
// once per up tunnel: the sum of their requests).
struct TunnelPlan {
  std::int64_t count {0};
  std::optional<Amount> design_demand {};
  std::vector<std::pair<NodeId, Amount>> head_design {};
  Amount reserved_total {0};
  friend bool operator==(const TunnelPlan&, const TunnelPlan&) = default;
};

// The result of one call: this call's values, then the demand's state after it (totals of the
// epoch), then the optional outputs.
// This call: placed (units delivered by the call; replace, forward and evaluate: the full new
// delivery), remaining = request - placed (nullopt for +inf), delivery_delta (new minus old
// delivery; placed for a placement, may be negative for replace), call_drops (this call's drops
// per edge, ascending EdgeId), call_source_loss, termination and the change flags, which compare
// the program before and after the call ("no program" is the empty program; two programs without
// an entry are equal): routes_changed the ordered heads, paths, split tables and duplicates;
// weights_changed the injection (the source weights) and every entry and split-table weight (an
// absent entry weighs 0); reservations_changed every tunnel's requested and reserved amounts by
// index (an absent tunnel holds (0, 0); only T1 sets it). An offer-only change to a program that
// is not recomputed (forward, Pinned, CostOnly increments) sets none of them; a controller's
// re-solve that changes a path's amount is a weight change; an identical replace sets none.
// State: program (the installed program after the call; evaluate: the computed one; null before
// the epoch's first preparation), plan (T1: the setup of that program; nullopt for every other
// model, Pinned over a T1 program too), requested (the sum of the requests;
// nullopt after +inf), offered (Lossy only: the traffic actually sent; the first-hop injection of a
// Lossy +inf), admitted (Lossless: the total carried; Lossy: offered), delivered, unadmitted
// (Lossless: requested - delivered, nullopt after +inf; Lossy: 0), source_loss and drops (per
// edge). Lossy: offered = source_loss + sum(drops) + delivered. by_source has one row per
// normalized input source in ascending NodeId; the rows sum to the totals except for Lossy with an
// empty roster, where the whole offer is source loss of the demand and every row is zero.
// by_tunnel has one row per program tunnel in ascending index (empty without tunnels); it is an
// alternative view of the same traffic, never added to by_source. total_cost (controllers: Global
// and CostOnly MaxFlow) is the exact sum of cost * flow over the policy's allocation (below 2^124;
// nullopt for the tunnel models and Pinned).
// cost_distribution: unrestricted Global with the optimized origin only, the cumulative marginal
// cost curve: (cost, units) pairs in ascending cost, where delivering `units` more units raised the
// minimum total cost of the flow by `cost` per unit (more delivery can reroute flow already placed,
// so the pairs are not the costs of the paths in the final allocation); empty for an empty roster,
// and absent before the epoch's first positive call. Optional outputs describe the same state as
// the rest of the result (evaluate: the evaluated state).
struct PlacementResult {
  Amount placed {0};
  std::optional<Amount> remaining {};
  Amount delivery_delta {0};
  std::vector<std::pair<EdgeId, Amount>> call_drops {};
  Amount call_source_loss {0};
  Termination termination {Termination::Complete};
  bool routes_changed {false};
  bool weights_changed {false};
  bool reservations_changed {false};
  std::shared_ptr<const ForwardingProgram> program {};
  std::optional<TunnelPlan> plan {};
  std::optional<Counter128> requested {};
  std::optional<Counter128> offered {};
  Counter128 admitted {};
  Amount delivered {0};
  std::optional<Counter128> unadmitted {};
  Counter128 source_loss {};
  std::vector<std::pair<EdgeId, Counter128>> drops {};
  std::vector<SourceFlow> by_source {};
  std::vector<TunnelFlow> by_tunnel {};
  std::optional<Counter128> total_cost {};
  std::optional<std::vector<std::pair<Cost, Amount>>> cost_distribution {};
  std::optional<std::vector<Amount>> edge_flows {};
  std::optional<std::vector<Amount>> residual {};
  std::optional<std::vector<EdgeId>> bottleneck {};
  std::optional<std::vector<Amount>> reservation_edge_flows {};
  // Requested by OutputOptions::min_cut; unrestricted Global with the optimized origin only, when
  // the solve reached the maximum flow: the edges leaving the set of nodes the roster sources reach
  // through residual arcs (free capacity plus this policy's own load, reverse arcs that cancel flow
  // included), whatever their free capacity (an allowed edge of capacity 0 is listed; an excluded
  // edge is not an edge here), ascending EdgeId. That is the minimum cut closest to the sources; an
  // empty roster returns an empty cut. nullopt otherwise (another model, the equal origin, a
  // support, a call below the maximum, a zero call, forward).
  std::optional<std::vector<EdgeId>> min_cut {};
};

// A policy binds to its ledger at construction. Every call validates, then computes on working
// copies and publishes all changes together, or none.
// place(v) adds v to the epoch: CostOnly prepares its program at the epoch's first positive call
// (every roster source's shortest routes on the ledger's residual at that call's entry) and keeps
// it for the epoch; EqualBalanced and Weighted forward the equal origin over the pooled routes:
// Lossless admits the largest total a in [admitted, admitted + volume] whose no-drop forwarding fits
// the residual plus the policy's own load (a refused unit is not offered again), Lossy forwards
// the new round-robin injections on the residual at call entry (earlier drops are not offered
// again). The controllers (Global, CostOnly MaxFlow) re-solve admitted + volume from zero on the
// residual plus their own allocation and export the solve as their program; the equal origin admits
// the largest total whose round-robin supplies are all carried. Pinned forwards its program the
// same way (Lossless: the largest total in [admitted, admitted + volume]; Lossy: the round-robin
// or captured-weight increments of the injection, each head splitting its increment afresh).
// T1 has no place (rejected, also for 0).
// evaluate(offer, design_demand) is a full, independent computation and forwarding of the offer
// from the policy's own released view (the residual plus its own carried load; T1 sets up on the
// unreserved capacity with its own holds released): it publishes nothing, creates no record and
// leaves the epoch as it was. replace(offer, design_demand) is the same, published atomically: the
// old records and holds are retired, new ones created, and a new epoch starts at this offer
// (delivery may go down). design_demand (T1 only; every other model rejects it): absent means D =
// offer, an inner nullopt +inf; T1 needs a finite D. forward(offer, black_holes) forwards a full
// new offer over the installed program with no route computation and no setup (Lossless: the
// fixed program's lossless bound; Lossy: the offer's injections) on the released view, starts a
// new epoch like replace, and keeps the program and T1's holds; with black_holes, the nodes and
// edges it excludes are black holes for that call. A Lossy replace, forward or evaluate needs a
// finite offer; a Lossy +inf is CostOnly place's first-hop experiment only. remove() releases
// every policy flow and hold, retires their ids, clears the program and starts a new epoch: the
// period since construction, the last remove(), replace() or forward() over which requests and
// flows accumulate. Destroying a policy without remove() leaves its records and holds in the
// ledger as orphans.
//
// Validation (std::invalid_argument, before anything binds), in this order: the endpoint sets
// (non-empty, ids in range, disjoint), the ownership of the Eligibility, the program, the support,
// the DagCache and the black holes by the ledger's graph, the configuration (FlowPolicyConfig's
// order), then the roster and count rules (the pooled single-path rule; T1's tunnel limit; a
// Pinned program's endpoints equal to the policy's). A call validates the operation for the model
// (T1's place), an installed program (forward), the request ([0, 2^62) units or +inf), a Lossy
// full offer's finiteness, the design demand (the model accepts one, [0, 2^62) or +inf, T1
// finite), a finite place after +inf in the epoch (rejected, also 0) and the black holes' graph
// (forward) before any change. A solve of the flow core (Global, CostOnly MaxFlow, a MaxFlow flow) whose
// network exceeds the core's int32 numbering throws std::invalid_argument, and the call changes
// nothing: more than INT32_MAX - 2 nodes, or more than 2^30 - 1 edges, sources and destinations
// together. A Lossy forwarding whose losses wait on each other in a cycle that alternates between
// two states or does not settle within its round limit throws std::invalid_argument; a failed
// check of an iterated forwarding throws std::runtime_error; neither publishes anything.
//
// Lifetime: the ledger, the Eligibility and the DagCache outlive the policy (the policy holds
// pointers to them); a Pinned policy shares its program and black holes. The roster (the eligible
// sources that statically reach an allowed destination, over allowed edges whatever their
// capacity, T1's heads; Pinned: the program's roster heads that reach an allowed destination over
// its basis under the Eligibility) is fixed at construction.
class FlowPolicy {
public:
  FlowPolicy(FlowLedger& ledger, std::span<const NodeId> sources, std::span<const NodeId> destinations,
             FlowPolicyConfig config, const Eligibility* elig = nullptr, DagCache* cache = nullptr);
  FlowPolicy(const FlowPolicy&) = delete;
  FlowPolicy& operator=(const FlowPolicy&) = delete;
  FlowPolicy(FlowPolicy&&) = delete;
  FlowPolicy& operator=(FlowPolicy&&) = delete;
  ~FlowPolicy();

  // volume, offer and an inner design_demand: units, nullopt is +inf; an absent design_demand is
  // the offer.
  PlacementResult place(std::optional<Amount> volume, const OutputOptions& outputs = {});
  [[nodiscard]] PlacementResult evaluate(std::optional<Amount> offer,
                                         const std::optional<std::optional<Amount>>& design_demand = std::nullopt,
                                         const OutputOptions& outputs = {}) const;
  PlacementResult replace(std::optional<Amount> offer,
                          const std::optional<std::optional<Amount>>& design_demand = std::nullopt,
                          const OutputOptions& outputs = {});
  // black_holes: an Eligibility of the policy's graph, or null for none (Pinned: rejected, its
  // black holes are part of its configuration).
  PlacementResult forward(std::optional<Amount> offer, const Eligibility* black_holes = nullptr,
                          const OutputOptions& outputs = {});
  void remove();
  // The sum of delivered over this policy's live records.
  [[nodiscard]] Amount placed_demand() const;
  // This policy's live records as snapshots in FlowId order.
  [[nodiscard]] std::vector<FlowRecord> records() const;
  // One of this policy's live records; a FlowId of another ledger or policy, a retired or an
  // unknown id throws std::invalid_argument.
  [[nodiscard]] FlowRecord record(const FlowId& flow_id) const;
  // The installed program: null before the epoch's first preparation and after remove().
  [[nodiscard]] std::shared_ptr<const ForwardingProgram> program() const noexcept;
  // The configuration as validated at construction: a copy whose absent origin, selection and
  // weights are resolved to the model's defaults (FlowPolicyConfig).
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
