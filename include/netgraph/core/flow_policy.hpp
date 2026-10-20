/* The policy engine FlowPolicy. */
#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include "netgraph/core/shortest_path_cache.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/forwarding_program.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

namespace detail {
struct PolicyEngine;
struct PolicyAccess;
class FrozenCacheHitTally;
} // namespace detail

// The tunnel count of the tunnel models, per head: `per_head` in [1, 10000] (checked where a policy
// or an evaluator uses it); `up_to` (T2 only; T1 rejects it): one common count k in 1..per_head for
// every head, the smallest whose program carries the live heads' design shares (else the largest
// live-share score, then the smallest count; no promise of whole-demand lossless delivery).
struct TunnelCount {
  std::int64_t per_head {1};
  bool up_to {false};
  friend bool operator==(const TunnelCount&, const TunnelCount&) = default;
};

// Each routing mode has only the fields it uses.
// ShortestPaths: hop-by-hop routing on each roster source's static shortest routes (all shortest
//   paths, with no directed cycles, edges of capacity 0 included, or the canonical single
//   path), prepared once per epoch and fixed for it. The sources of the demand are pooled: a
//   router shared by several sources adds their traffic and splits it once. Equal and Weighted forward the equal source allocation; Optimized is a controller
//   restricted to the union of those routes. Its program: one split route per roster source,
//   pooled into one split context (Equal, Weighted), or the controller's paths (Optimized). Its
//   route choice: tie_break (RouteTieBreak), multipath (false: the canonical single path) and
//   parallel_edges (false: one member of each (u, v) pair).
// DistributedTE (T1): every roster head sets up its own tunnel_count.per_head strict tunnels by CSPF, no
//   controller. A design volume D is split RR(D, H, h) over the heads and each head's share
//   RR(d, N, j) over its tunnels; in round order (round j, heads ascending: index j * H + i) each
//   tunnel keeps the allowed links whose unreserved capacity U (the ledger's, this policy's own
//   holds released, minus this setup's earlier up tunnels) is at least max(request,
//   min_unreserved_bandwidth) (a request of 0 with a minimum of 0 has no filter), collapses
//   parallel links (parallel_edges = false), and takes the single path in the order of tie_break
//   on static costs (scores read U: HigherAvailableCapacity U, HigherHeadroom U / C, HigherCapacity
//   C); it is Up, holding its
//   request on every link, unless no path passes the filter (InsufficientBandwidth) or the cheapest
//   exceeds cost_ceiling (CostCeilingExceeded): AbsoluteCost inclusive, CostStretch(f) times c0, the
//   smallest static
//   shortest cost over the roster heads under the mask, zero-capacity links included. Forwarding
//   is the program's: round robin over the heads, then over a head's up tunnels; a head without
//   one loses its share. Its program: those tunnels, weight 1 up and 0 down.
// CentralizedTE: the ideal controller, a minimum-cost maximum flow over the whole eligible graph, or
//   exactly over `support` (a DAG of the policy's graph, trimmed to the walks from the sources to
//   the allowed destinations under the Eligibility), of what the policy already delivers plus this
//   call's request (what the policy delivers plus the offer: what an earlier call could not place is not
//   requested again). Its program: the solved flow cut into paths per source, each path a tunnel
//   weighted by its amount.
// CentralizedTE with `tunnel_count` (T2; with `support`, the allowed links are the support's edges
//   intersected with the Eligibility, in every step below): a controller lays out, for every live
//   roster head (one reaching an allowed destination over links with free capacity V >= 1, V the
//   residual with this policy's own traffic released; the other heads are dead: they keep their
//   roster position and lose their share), k single-path tunnels with one common whole-unit planned
//   bandwidth: first the largest common size U (quota flows of K = a k path units, a live heads),
//   then for a finite design D the working size u = min(U, max(1, max live ceil(RR(D, H, i) / k)))
//   (unlimited: U), then the cheapest routes (the route-cost sum J), then with prefer_headroom
//   the smallest headroom penalty (each tunnel on a link pays floor(2^32 (B - V + j u) / B) for its
//   slot j, B the static capacity), then one fixed unit-by-unit augmentation (fewest arcs, then the
//   smallest arc keys). k = tunnel_count.per_head, or with up_to the smallest count in 1..per_head whose
//   program carries the live part of D (else the best, then the smallest count). U = 0 with live
//   heads, or no live head, gives a program without tunnels that keeps the roster. Forwarding is
//   the program's: round robin over the heads, then over a head's tunnels. Its program: the
//   tunnels (index j * H + i, weight 1, planned u); the plan names U, u, s* and the objectives.
// Pinned: forwarding over a given program, computing nothing; the program's endpoints are the
//   policy's. Its Eligibility repairs failures (a tunnel touching an excluded element is Failed, a
//   split route loses its excluded members, the roster is rebuilt from the program's basis); the
//   nodes and edges `black_holes` excludes are black holes (stale forwarding: the program is
//   forwarded unchanged and traffic reaching one is dropped there). One record per program entry.
struct ShortestPaths {
  RouteAllocation allocation {RouteAllocation::Equal};
  std::optional<SourceAllocation> source_allocation {};
  std::optional<WeightBasis> weight_basis {};
  std::shared_ptr<const std::vector<Amount>> edge_weights {};
  RouteTieBreak tie_break {RouteTieBreak::FewerHops};
  bool multipath {true};
  bool parallel_edges {true};
};
struct DistributedTE {
  TunnelCount tunnel_count {};
  std::optional<CostCeiling> cost_ceiling {};
  Amount min_unreserved_bandwidth {0};
  RouteTieBreak tie_break {RouteTieBreak::FewerHops};
  bool parallel_edges {true};
};
// CentralizedTE requires its source allocation: there is no default constructor.
struct CentralizedTE {
  explicit CentralizedTE(SourceAllocation source_allocation_) : source_allocation(source_allocation_) {}
  SourceAllocation source_allocation;
  std::shared_ptr<const PredecessorDAG> support {};
  std::optional<TunnelCount> tunnel_count {};
  bool prefer_headroom {false};
  friend bool operator==(const CentralizedTE&, const CentralizedTE&) = default;
};
struct Pinned {
  std::shared_ptr<const ForwardingProgram> program {};
  std::shared_ptr<const Eligibility> black_holes {};
};
using Routing = std::variant<ShortestPaths, DistributedTE, CentralizedTE, Pinned>;

// One demand's settings: the routing model, which holds every option it reads, and the admission.
// A policy or an evaluator checks the configuration when it uses it (the fields remain
// assignable), in this order: the enumeration values; then the model matrix below, in the order
// route and source allocation (ShortestPaths EQUAL or WEIGHTED takes only the equal source
// allocation; T2 with several input sources only the equal one), admission, the route-choice
// controls (the ineffective tie rule, then the pooled single-path rule), prefer_headroom, the
// weight combinations and, last, a ShortestPathCache on a model that does not read it; then the
// tunnel count, an "up to" count for DistributedTE, the ceilings, min_unreserved_bandwidth and
// edge_weights' length and entries; then the tunnel limit; then a Pinned program's endpoints.
//
// | model                             | allocation    | source allocation (default / allowed)  | admission       | route choice                            | weights |
// | ShortestPaths                     | Equal         | Equal / Equal                          | Lossless, Lossy | tie_break, multipath, parallel_edges    | none |
// | ShortestPaths                     | Weighted      | Equal / Equal                          | Lossless, Lossy | tie_break, multipath, parallel_edges    | Capacity (default), Explicit, ResidualSnapshot |
// | ShortestPaths                     | Optimized     | Optimized / both                       | Lossless        | tie_break, multipath, parallel_edges    | none |
// | CentralizedTE                     | optimal       | required / both                        | Lossless        | none                                    | none |
// | DistributedTE (T1)                | tunnels       | Equal, fixed                           | Lossless, Lossy | tie_break, parallel_edges (single path) | none |
// | CentralizedTE with tunnels (T2)   | tunnels       | required; Optimized only with one source | Lossless, Lossy | none                                  | none |
// | Pinned                            | the program's | the program's                          | Lossless, Lossy | none                                    | none |
//
// A Higher* tie_break of ShortestPaths with multipath and parallel_edges both true has nothing to
// choose and is rejected. The policy's snapshot (config()) resolves ShortestPaths' absent
// source_allocation (Optimized for Optimized, else Equal) and an absent weight_basis of Weighted
// (Capacity). edge_weights (one integer per edge in [0, 2^62), shared by every configuration
// copied from this one) is given exactly when weight_basis is Explicit. T2's headroom objective is
// CentralizedTE's prefer_headroom (true needs tunnels); a tie rule never enables a controller
// objective. Pooled single-path rule: ShortestPaths with more than one input source (normalized,
// before any Eligibility) and multipath = false requires tie_break = FewerHops. The tunnel limit
// (T1, T2): roster heads times the tunnel count (the largest count for up_to) at most 2^30 - 1.
struct FlowPolicyConfig {
  Routing routing {ShortestPaths {}};
  Admission admission {Admission::Lossless};
};

// Optional outputs of a placement (absent unless requested): edge_flows (this policy's own
// allocation per edge after the call), residual and saturated_boundary_edges (the ledger's after
// the call), min_cut (CentralizedTE with the optimized source allocation; see
// PlacementResult::min_cut) and edge_reservations (this policy's strict holds per edge after the
// call, zeros for a model without holds).
struct OutputOptions {
  bool edge_flows {false};
  bool residual {false};
  bool saturated_boundary_edges {false};
  bool min_cut {false};
  bool edge_reservations {false};
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
// delivers; Lossy: offered), delivered, source_drops (offered shares lost at the source for lack
// of a usable route) and total_edge_drops (its losses on links). Lossy: offered = delivered +
// source_drops + total_edge_drops.
struct SourceFlow {
  NodeId source {0};
  std::optional<Counter128> offered {};
  Counter128 admitted {};
  Amount delivered {0};
  Counter128 source_drops {};
  Counter128 total_edge_drops {};
  friend bool operator==(const SourceFlow&, const SourceFlow&) = default;
};

// One program tunnel's forwarding in the epoch: its entry identifier and head, its status in this
// forwarding (Up, Failed, Blackholed), assigned (its share of the injections: Lossy the shares sent
// to it, Lossless what it carries), delivered and total_edge_drops (its losses on links).
struct TunnelFlow {
  std::int32_t entry_id {0};
  NodeId head {0};
  TunnelStatus status {TunnelStatus::Up};
  Counter128 assigned {};
  Amount delivered {0};
  Counter128 total_edge_drops {};
  friend bool operator==(const TunnelFlow&, const TunnelFlow&) = default;
};

// The tunnel models' setup summary (PlacementResult::plan): tunnels_per_head (T2 up to: the
// chosen count), design_demand (nullopt: unlimited), design_by_head ((head, RR(D, H, i)) for every
// roster head in ascending NodeId; empty for an unlimited design) and total_reserved_bandwidth
// (T1's strict holds, counted once per up tunnel: the sum of their requests; 0 for T2). T2 only
// (nullopt for T1): heads_with_capacity (ascending), insufficient_resolution (heads with capacity
// exist but U = 0: the program has no tunnels), max_fractional_bandwidth_per_tunnel (the exact
// fluid optimum s* as a reduced fraction (numerator, denominator); 0 / 1 without such a head),
// max_bandwidth_per_tunnel (U), planned_bandwidth_per_tunnel (u), total_planned_bandwidth (K u),
// total_tunnel_cost (J, below 2^92), headroom_penalty (Phi, below 2^94; with prefer_headroom only)
// and max_lossless_offer (F: the program's lossless maximum on the call's view, a dead head's
// share blocking).
struct TunnelPlan {
  std::int64_t tunnels_per_head {0};
  std::optional<Amount> design_demand {};
  std::vector<std::pair<NodeId, Amount>> design_by_head {};
  Amount total_reserved_bandwidth {0};
  std::optional<std::vector<NodeId>> heads_with_capacity {};
  std::optional<bool> insufficient_resolution {};
  std::optional<std::pair<Amount, Amount>> max_fractional_bandwidth_per_tunnel {};
  std::optional<Amount> max_bandwidth_per_tunnel {};
  std::optional<Amount> planned_bandwidth_per_tunnel {};
  std::optional<Amount> total_planned_bandwidth {};
  std::optional<Counter128> total_tunnel_cost {};
  std::optional<Counter128> headroom_penalty {};
  std::optional<Amount> max_lossless_offer {};
  friend bool operator==(const TunnelPlan&, const TunnelPlan&) = default;
};

// The result of one call: this call's values, then the demand's state after it (totals of the
// epoch), then the optional outputs.
// This call: call_delivered (units delivered by the call; replace, forward and evaluate: the full
// new delivery), call_undelivered = request - call_delivered (nullopt for +inf), delivery_delta
// (new minus old delivery; call_delivered for a placement, may be negative for replace),
// call_edge_drops (this call's drops per edge, ascending EdgeId), call_source_drops, termination
// and the change flags, which compare
// the program before and after the call ("no program" is the empty program; two programs without
// an entry are equal): routes_changed the ordered heads, paths, split tables and duplicates;
// weights_changed the injection (the source weights) and every entry and split-table weight (an
// absent entry weighs 0); reservations_changed every tunnel's requested and reserved amounts by
// index (an absent tunnel holds (0, 0); only T1 sets it: T2's planned bandwidth is no
// reservation). An offer-only change to a program that
// is not recomputed (forward, Pinned, ShortestPaths increments) sets none of them; a controller's
// re-solve that changes a path's amount is a weight change; an identical replace sets none.
// State: program (the installed program after the call; evaluate: the computed one; null before
// the epoch's first preparation), plan (T1 and T2: the setup of that program; nullopt for every
// other model, Pinned over a tunnel program too), requested (the sum of the requests;
// nullopt after +inf), offered (Lossy only: the traffic actually sent; the first-hop injection of a
// Lossy +inf), admitted (Lossless: the total carried; Lossy: offered), delivered, unadmitted
// (Lossless: requested - delivered, nullopt after +inf; Lossy: 0), source_drops and edge_drops
// (per edge). Lossy: offered = source_drops + sum(edge_drops) + delivered. by_source has one row per
// normalized input source in ascending NodeId; the rows sum to the totals except for Lossy with an
// empty roster, where the whole offer is source loss of the demand and every row is zero.
// by_tunnel has one row per program tunnel in ascending index (empty without tunnels); it is an
// alternative view of the same traffic, never added to by_source. total_cost (controllers: CentralizedTE
// and ShortestPaths Optimized) is the exact sum of cost * flow over the policy's allocation (below 2^124;
// nullopt for the tunnel models and Pinned).
// marginal_cost_curve: unrestricted CentralizedTE with the optimized source allocation only, the cumulative marginal
// cost curve: (cost, units) pairs in ascending cost, where delivering `units` more units raised the
// minimum total cost of the flow by `cost` per unit (more delivery can reroute flow already placed,
// so the pairs are not the costs of the paths in the final allocation); empty for an empty roster,
// and absent before the epoch's first positive call. Optional outputs describe the same state as
// the rest of the result (evaluate: the evaluated state).
struct PlacementResult {
  Amount call_delivered {0};
  std::optional<Amount> call_undelivered {};
  Amount delivery_delta {0};
  std::vector<std::pair<EdgeId, Amount>> call_edge_drops {};
  Amount call_source_drops {0};
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
  Counter128 source_drops {};
  std::vector<std::pair<EdgeId, Counter128>> edge_drops {};
  std::vector<SourceFlow> by_source {};
  std::vector<TunnelFlow> by_tunnel {};
  std::optional<Counter128> total_cost {};
  std::optional<std::vector<std::pair<Cost, Amount>>> marginal_cost_curve {};
  std::optional<std::vector<Amount>> edge_flows {};
  std::optional<std::vector<Amount>> residual {};
  std::optional<std::vector<EdgeId>> saturated_boundary_edges {};
  std::optional<std::vector<Amount>> edge_reservations {};
  // Requested by OutputOptions::min_cut; unrestricted CentralizedTE with the optimized source allocation only, when
  // the solve reached the maximum flow: the edges leaving the set of nodes the roster sources reach
  // through residual arcs (free capacity plus this policy's own load, reverse arcs that cancel flow
  // included), whatever their free capacity (an allowed edge of capacity 0 is listed; an excluded
  // edge is not an edge here), ascending EdgeId. That is the minimum cut closest to the sources; an
  // empty roster returns an empty cut. nullopt otherwise (another model, the equal source allocation, a
  // support, a call below the maximum, a zero call, forward).
  std::optional<std::vector<EdgeId>> min_cut {};
};

// A policy binds to its ledger at construction. Every call validates, then computes on working
// copies and publishes all changes together, or none.
// place(offer) adds offer to the epoch: ShortestPaths prepares its program at the epoch's first positive call
// (every roster source's shortest routes on the ledger's residual at that call's entry) and keeps
// it for the epoch; Equal and Weighted forward the equal source allocation over the pooled routes:
// Lossless admits the largest total a in [admitted, admitted + offer] whose no-drop forwarding fits
// the residual plus the policy's own load (a refused unit is not offered again), Lossy forwards
// the new round-robin injections on the residual at call entry (earlier drops are not offered
// again). The controllers (CentralizedTE, ShortestPaths Optimized) re-solve admitted + offer from zero on the
// residual plus their own allocation and export the solve as their program; the equal source allocation admits
// the largest total whose round-robin supplies are all carried. Pinned forwards its program the
// same way (Lossless: the largest total in [admitted, admitted + offer]; Lossy: the round-robin
// or captured-weight increments of the injection, each head splitting its increment afresh).
// T1 and T2 have no place (rejected, also for 0).
// evaluate(offer, design_demand) is a full, independent computation and forwarding of the offer
// from the policy's own released view (the residual plus its own carried load; T1 sets up on the
// unreserved capacity with its own holds released): it publishes nothing, creates no record and
// leaves the epoch as it was. replace(offer, design_demand) is the same, published atomically: the
// old records and holds are retired, new ones created, and a new epoch starts at this offer
// (delivery may go down). design_demand (T1 and T2; every other model rejects it): absent means D
// = offer, an inner nullopt +inf; T1 needs a finite D, T2 takes an unlimited one (u = U). T2 lays
// out its tunnels on the residual with its own traffic released. forward(offer, black_holes) forwards a full
// new offer over the installed program with no route computation and no setup (Lossless: the
// fixed program's lossless bound; Lossy: the offer's injections) on the released view, starts a
// new epoch like replace, and keeps the program, its plan and T1's holds; with black_holes, the nodes and
// edges it excludes are black holes for that call. A Lossy replace, forward or evaluate needs a
// finite offer; a Lossy +inf is ShortestPaths place's first-hop experiment only. remove() releases
// every policy flow and hold, retires their ids, clears the program and starts a new epoch: the
// period since construction, the last remove(), replace() or forward() over which requests and
// flows accumulate. Destroying a policy without remove() leaves its records and holds in the
// ledger as orphans.
//
// Validation (std::invalid_argument, before anything binds), in this order: the endpoint sets
// (non-empty, ids in range, disjoint), the ownership of the Eligibility, the program, the support,
// the ShortestPathCache and the black holes by the ledger's graph, the configuration (FlowPolicyConfig's
// order), then the roster and count rules (the pooled single-path rule; T2's several-sources
// rule; the tunnel limit; a Pinned program's endpoints equal to the policy's). A call validates the
// operation for the model (T1's and T2's place), an installed program (forward), the request ([0,
// 2^62) units or +inf), a Lossy full offer's finiteness, the design demand (the model accepts one,
// [0, 2^62) or +inf, T1 finite), a finite place after +inf in the epoch (rejected, also 0) and the black holes' graph
// (forward) before any change. A solve of the flow core (CentralizedTE, ShortestPaths Optimized, a Optimized flow) whose
// network exceeds the core's int32 numbering throws std::invalid_argument, and the call changes
// nothing: more than INT32_MAX - 2 nodes, or more than 2^30 - 1 edges, sources and destinations
// together. A Lossy forwarding whose losses wait on each other in a cycle that alternates between
// two states or does not settle within its round limit throws std::invalid_argument; a failed
// check of an iterated forwarding throws std::runtime_error; neither publishes anything.
//
// Lifetime: the ledger, the Eligibility and the ShortestPathCache outlive the policy (the policy holds
// pointers to them); a Pinned policy shares its program and black holes. The roster (the eligible
// sources that statically reach an allowed destination, over allowed edges whatever their
// capacity, the tunnel models' heads; Pinned: the program's roster heads that reach an allowed destination over
// its basis under the Eligibility) is fixed at construction.
class FlowPolicy {
public:
  FlowPolicy(FlowLedger& ledger, std::span<const NodeId> sources, std::span<const NodeId> destinations,
             FlowPolicyConfig config, const Eligibility* eligibility = nullptr, ShortestPathCache* cache = nullptr);
  FlowPolicy(const FlowPolicy&) = delete;
  FlowPolicy& operator=(const FlowPolicy&) = delete;
  FlowPolicy(FlowPolicy&&) = delete;
  FlowPolicy& operator=(FlowPolicy&&) = delete;
  ~FlowPolicy();

  // offer and an inner design_demand: units, nullopt is +inf; an absent design_demand is
  // the offer.
  PlacementResult place(std::optional<Amount> offer, const OutputOptions& outputs = {});
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
  [[nodiscard]] Amount delivered() const;
  // This policy's live records as snapshots in FlowId order.
  [[nodiscard]] std::vector<FlowRecord> records() const;
  // One of this policy's live records; a FlowId of another ledger or policy, a retired or an
  // unknown id throws std::invalid_argument.
  [[nodiscard]] FlowRecord record(const FlowId& flow_id) const;
  // The installed program: null before the epoch's first preparation and after remove().
  [[nodiscard]] std::shared_ptr<const ForwardingProgram> program() const noexcept;
  // The configuration as validated at construction: a copy whose absent source_allocation and
  // weight_basis are resolved to the model's defaults (FlowPolicyConfig).
  [[nodiscard]] const FlowPolicyConfig& config() const noexcept;

private:
  friend struct detail::PolicyAccess;
  // Construction inside another public call (an evaluator): a supplied frozen cache's hits are
  // recorded in that call's tally.
  FlowPolicy(FlowLedger& ledger, std::span<const NodeId> sources, std::span<const NodeId> destinations,
             FlowPolicyConfig config, const Eligibility* eligibility, ShortestPathCache* cache, detail::FrozenCacheHitTally& tally);
  std::unique_ptr<detail::PolicyEngine> engine_;
};

// ---------------------------------------------------------------- several demands at once

// How apply_actions runs its list: Ordered, each action as its scalar call on the state the earlier
// ones left; Simultaneous, every action forwarded at once with the loss rule across demands.
enum class BatchMode { Ordered, Simultaneous };

// The actions of a list. Offers are as in FlowPolicy::place (nullopt: +inf); Replace's design demand
// as in FlowPolicy::replace (outer nullopt absent, inner nullopt +inf); Forward's black holes an
// Eligibility of the policy's graph or null. They borrow the policy and the black holes, which must
// outlive the call.
struct Place {
  std::reference_wrapper<FlowPolicy> policy;
  std::optional<Amount> offer;
};
struct Replace {
  std::reference_wrapper<FlowPolicy> policy;
  std::optional<Amount> offer;
  std::optional<std::optional<Amount>> design_demand {};
};
struct Forward {
  std::reference_wrapper<FlowPolicy> policy;
  std::optional<Amount> offer;
  const Eligibility* black_holes {nullptr};
};
using Action = std::variant<Place, Replace, Forward>;

// Several policies of one ledger in one operation, validated whole, then computed, then published
// once (one atomic exchange of every policy's records, then their epochs): an exception leaves the
// ledger and every policy unchanged. Results in list order, each the scalar call's fields for its
// policy.
//
// Ordered: each action runs in list order on a private working state (the ledger's carried loads and
// holds as the earlier actions left them) exactly as its scalar call would; its outputs (residual,
// saturated_boundary_edges) read the state after it; new records take FlowIds in list order, then in
// each call's draft order. Simultaneous: Replace of ShortestPaths Equal or Weighted (its program
// prepared on the view) and of Pinned, Forward of any policy with an installed program, every one
// Lossy with a finite offer, forwarded at once on one view (the ledger's residual plus every
// participant's old carried load; outside owners and orphans fixed, holds never added) by the loss
// rule across demands: at a link that cannot carry its arrivals, highest averages over the demands
// (ties by policy creation order, not list position), then within a demand over its sources, then
// its entries; a pooled ShortestPaths demand splits once at each router and attributes each member's
// amount to its sources; each participant's black holes are its own (what reaches one is dropped
// there); losses that wait on each other in a cycle across demands are iterated as for one program,
// with its two errors. Every entry reports the one residual after the list and its own
// saturated_boundary_edges on it. An empty list returns nothing after checking the ledger and mode.
//
// Validation (std::invalid_argument, before anything is computed): the ledger (a moved-from one:
// std::logic_error), `mode` a BatchMode value (`mode: <v> is not a value of BatchMode`); then each
// action in list order, all of its checks before the next one's: `list policy belongs to another
// ledger`; `duplicate policy in list` (at the second occurrence); Simultaneous only, in this order,
// `SIMULTANEOUS does not accept Place`, `SIMULTANEOUS Replace requires ShortestPaths EQUAL or WEIGHTED,
// or Pinned`, `SIMULTANEOUS requires admission=LOSSY`, `SIMULTANEOUS requires finite offers`; then
// the scalar call's own checks in their order. A computation error (the loss iteration's errors, a
// solver or forwarding size, the record identity limit counted for the whole list) also changes
// nothing.
std::vector<PlacementResult> apply_actions(FlowLedger& ledger, std::span<const Action> actions,
                                           BatchMode mode = BatchMode::Ordered, const OutputOptions& outputs = {});

} // namespace netgraph::core
