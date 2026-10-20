/*
  The policy engine FlowPolicy: one transaction for every routing mode and operation, published
  through the ledger as one atomic exchange per call.

  A call computes from the epoch it continues (place) or a fresh one (evaluate, replace and
  forward, on the residual with the policy's own load released) into a staged value: the result,
  the drafts of the commit and the epoch after it. place, replace and forward publish it (the
  ledger commit, then the epoch); evaluate discards it. The change flags compare the program the
  policy had installed with the call's.

  * ShortestPaths Equal and Weighted (place_pooled): the epoch's first positive call (or a
    replace) prepares every roster source's route on the working residual at its entry and joins
    them into the program's one pooled split context, a plan rooted at the whole roster (a router
    reached by several sources splits their added traffic once). Weighted weights come from the
    configured basis (capacity, explicit edge weights, or the free capacity read then and frozen).
    Lossless admits the largest total of the equal source allocation whose no-drop forwarding fits V = W +
    own; Lossy forwards the round-robin increments of the offer on W. The amounts are attributed
    back to the sources (kernels' forward()).
  * The controllers (place_optimized_flow): unrestricted CentralizedTE over every eligible edge, CentralizedTE with a
    support over that support trimmed, ShortestPaths Optimized over the union of the roster's routes
    (prepared as above). A minimum-cost flow of the cumulative value on V = W + own, from one
    virtual source (the optimized source allocation) or with source i's arc capped at its round-robin share of
    the largest total whose shares all fit (the equal source allocation, integer bisection), cut into paths
    per source (decompose()): the program, each path a tunnel weighted by its amount; each source's
    record is the sum of its paths.
  * Pinned (place_pinned): forwarding its program as the policy's view replays it under its mask
    with its black holes, through the forwarding engine (forwarding.hpp), one record per program
    entry (forward_entries); forward() (forward_installed) forwards a computing model's installed
    program the same way, with records per roster source.
  * T1 (apply_distributed_te_call): replace and evaluate set up the tunnels (distributed_te.cpp) on the unreserved view
    with the policy's own holds released, build the tunnel program and forward it like a Pinned
    program, one record per entry, each up tunnel's record holding its request on its path.
  * T2 (place_t2): replace and evaluate lay out the tunnels (centralized_te_tunnels.cpp) on the residual with
    the policy's own traffic released, build the tunnel program (dead heads keep their roster
    place), take its lossless bound as the plan's forwarding bound and forward it the same way; a
    plan books nothing.
  * The tunnel models' forward() (forward_tunnels) forwards the installed program the same way,
    keeping its plan and T1's holds.

  apply_actions (several policies of one ledger): the list is validated whole, then computed (ORDERED:
  each action's call on a working copy of the ledger's totals; SIMULTANEOUS: every participant
  prepared on one view, forwarded by one joint evaluation, finished as its call), then published
  as one exchange of every owner's records followed by their epochs.

  Work on copies and publish all changes together; an exception leaves ledger and policy intact.
  Store only changes to ledger residual in per-thread slots. Unrestricted CentralizedTE reads a dense
  residual and solves over every edge (its view and flows per EdgeId); the other controllers
  gather their view over their support and keep their flows aligned with it. T1's setup reads a
  dense unreserved view; ShortestPaths preparation reads the working residual edge by edge.

  The transaction trace (the call's intermediate values, keyed as the tests' model of the policy
  keys them) exists in a build with the private hooks only: every
  trace block is guarded by kTraced and is compile-time dead without them.
*/
#include "netgraph/core/flow_policy.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

#include "netgraph/core/detail/cache_access.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/min_cost_flow.hpp"
#include "netgraph/core/detail/forwarding.hpp"
#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/route_kernels.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/route_builders.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/detail/trace.hpp"
#include "netgraph/core/detail/search_detail.hpp"
#include "netgraph/core/detail/tunnel_setup.hpp"

namespace netgraph::core {

namespace detail {

namespace {

// ================================================================ working residual

// The ledger's totals a call reads: the ledger's own (a scalar call), or a list's working copy
// (apply_actions: ORDERED's totals after the earlier actions, SIMULTANEOUS's with every
// participant's load released), which outlives every reader. `reservation_baseline` is the
// ledger's (the capacities when it has none); an empty `held` reads as all zero, as the ledger's
// unstaged holds do.
struct LedgerTotals {
  std::span<const Amount> service_baseline;
  std::span<const Amount> carried;
  std::span<const Amount> reservation_baseline;
  std::span<const Amount> held;
};

LedgerTotals totals_of(const FlowLedger& ledger) {
  return {LedgerAccess::baseline(ledger), LedgerAccess::carried(ledger), LedgerAccess::reservation_baseline(ledger),
          LedgerAccess::reserved(ledger)};
}

// W = service baseline - carried() + delta; delta holds this call's releases and new allocations.
// Each call advances a scratch epoch (a numeric stamp, separate from the policy epoch).
// Read a slot only when its stamp matches; later calls ignore stale slots even after a throw.
// Clear all stamps only on wraparound (after 2^32 calls on one thread, which a long-running
// process reaches). Slots share one per-thread array, so failed growth leaves the old array
// intact.
class ServiceCapacityWorkspace {
public:
  explicit ServiceCapacityWorkspace(const LedgerTotals& totals)
      : storage_(storage()), service_baseline_(totals.service_baseline), carried_(totals.carried) {
    if (storage_.slot.size() < service_baseline_.size()) storage_.slot.resize(service_baseline_.size());
    if (++storage_.generation == 0) {
      for (auto& x : storage_.slot) x.stamp = 0;
      storage_.generation = 1;
    }
  }
  ServiceCapacityWorkspace(const ServiceCapacityWorkspace&) = delete;
  ServiceCapacityWorkspace& operator=(const ServiceCapacityWorkspace&) = delete;

  [[nodiscard]] Amount at(EdgeId e) const noexcept {
    const auto ei = static_cast<std::size_t>(e);
    const auto& x = storage_.slot[ei];
    return service_baseline_[ei] - carried_[ei] + (x.stamp == storage_.generation ? x.delta : 0);
  }
  void add(EdgeId e, Amount v) noexcept {
    auto& x = storage_.slot[static_cast<std::size_t>(e)];
    if (x.stamp != storage_.generation) {
      x.stamp = storage_.generation;
      x.delta = 0;
    }
    x.delta += v;
  }
  // Unused capacity for each edge of `members`, in the same order (member-aligned).
  void gather(std::span<const EdgeId> members, std::vector<Amount>& out) const {
    out.resize(members.size());
    for (std::size_t m = 0; m < members.size(); ++m) out[m] = at(members[m]);
  }
  // The dense working residual (searches, CentralizedTE, traces).
  void dense(std::vector<Amount>& out) const {
    out.resize(service_baseline_.size());
    for (std::size_t e = 0; e < out.size(); ++e) {
      const auto& x = storage_.slot[e];
      out[e] = service_baseline_[e] - carried_[e] + (x.stamp == storage_.generation ? x.delta : 0);
    }
  }
  // The policy's own carried load released (the view of a full-offer call): its live records'
  // edges added back.
  void release_records(std::span<const FlowRecord* const> records) noexcept {
    for (const FlowRecord* r : records) {
      for (const auto& [e, x] : r->edge_flows) add(e, x);
    }
  }

private:
  struct Slot {
    Amount delta {0};
    std::uint32_t stamp {0};
  };
  struct Storage {
    std::vector<Slot> slot;
    std::uint32_t generation {0};
  };
  static Storage& storage() {
    thread_local Storage st;
    return st;
  }
  Storage& storage_;
  std::span<const Amount> service_baseline_;
  std::span<const Amount> carried_;
};

// ================================================================ small helpers

[[noreturn]] void invalid(const std::string& what) { throw std::invalid_argument(what); }

// The tunnel count per head and the tunnels of one program (H x N, int32 entry indices).
constexpr std::int64_t kMaxTunnelCount = 10000;
constexpr std::int64_t kMaxTunnels = (std::int64_t {1} << 30) - 1;

// The policy's routing mode (the tag of the Routing variant).
template <class R>
bool is(const PolicyEngine& pe) noexcept {
  return std::holds_alternative<R>(pe.config.routing);
}

std::string str(std::int64_t v) { return std::to_string(v); }

// The member-aligned amounts of a record's (EdgeId, amount) pairs (both ascending).
template <class T>
void align(std::span<const EdgeId> members, const std::vector<std::pair<EdgeId, T>>& pairs, std::vector<T>& out) {
  out.assign(members.size(), T {});
  std::size_t m = 0;
  for (const auto& [e, x] : pairs) {
    while (m < members.size() && members[m] < e) ++m;
    if (m == members.size() || members[m] != e) {
      throw std::logic_error("policy: a record's allocation lies outside its flow's split graph");
    }
    out[m] = x;
  }
}

// A ShortestPaths preparation's inputs: the policy, its graph and the call's working residual.
struct PolicyCallContext {
  const PolicyEngine& engine;
  const Graph& graph;
  ServiceCapacityWorkspace& service_capacity;
};

// ================================================================ routes and plans

struct RouteSupport {
  std::vector<EdgeId> edges {};  // ascending
  Cost cost {kUnreachable};
};

// The edge scores of the creation view W, as the route choice reads them (ShortestPaths). The
// PriorityReader of `as_priority_reader` points at this view: it must outlive every call through it.
struct WorkingPriority {
  const ServiceCapacityWorkspace* service_capacity;
  std::span<const Amount> capacity;
  RouteTieBreak tie_break;
  [[nodiscard]] PriorityReader as_priority_reader() const noexcept {
    return {this, [](const void* ctx, EdgeId e) noexcept {
              const auto& v = *static_cast<const WorkingPriority*>(ctx);
              return priority_value(v.tie_break, v.service_capacity->at(e), v.capacity[static_cast<std::size_t>(e)]);
            }};
  }
};

// The record route of a route: `share` itself when it is the whole route, else a new DAG.
std::shared_ptr<const PredecessorDAG> record_route_of(const Graph& graph, std::span<const EdgeId> route,
                                         const std::shared_ptr<const PredecessorDAG>& share) {
  if (share && PredecessorDAGBuilder::ascending(*share).size() == route.size()) return share;
  return std::make_shared<const PredecessorDAG>(PredecessorDAGBuilder::of_edges(graph, route));
}

// ShortestPaths: keep static shortest-path support reachable from source, then apply the route choice on W.
// For parallel_edges = false (Collapse), choose one edge per (u, v) by (cost, higher tie value under
// a HIGHER_* tie_break, lowest EdgeId). For multipath = false, choose the canonical single
// path in stages. Both use W's tie values without excluding edges for low residual.
RouteSupport shortest_path_support(const PolicyCallContext& c, NodeId source) {
  const auto& g = c.graph;
  const auto& pe = c.engine;
  const auto& choice = pe.route_choice;
  const WorkingPriority view {&c.service_capacity, g.edge_capacity_units(), choice.tie_break};
  auto r = support_route(g, PredecessorDAGBuilder::ascending(*pe.support), source, pe.allowed_destinations, choice, view.as_priority_reader());
  RouteSupport out;
  out.edges = std::move(r.edges);
  if (choice.single_path() && !out.edges.empty()) out.cost = r.cost;
  return out;
}

// ShortestPaths Weighted's static weight source (Capacity or Explicit; ResidualSnapshot reads the view
// at preparation and is made there).
WeightSource shortest_path_weight_source(const PolicyEngine& pe) noexcept {
  const auto& sp = std::get<ShortestPaths>(pe.config.routing);
  if (sp.weight_basis == WeightBasis::Explicit) return {WeightBasis::Explicit, *sp.edge_weights};
  return {sp.weight_basis.value_or(WeightBasis::Capacity), {}};
}

// ShortestPaths's creation products in one pass (prepare_shortest_paths, identical to
// shortest_path_support, record_route_of and the plan where it applies); nullopt elsewhere.
std::optional<PreparedRouteResult> creation_products(const PolicyCallContext& c, NodeId source) {
  const auto& pe = c.engine;
  const auto& choice = pe.route_choice;
  const WorkingPriority view {&c.service_capacity, c.graph.edge_capacity_units(), choice.tie_break};
  return prepare_shortest_paths(c.graph, pe.support, source, pe.allowed_destinations, choice, view.as_priority_reader(), std::get<ShortestPaths>(pe.config.routing).allocation,
                            pe.config.admission, shortest_path_weight_source(pe));
}

// ---------------------------------------------------------------- outputs

// The optional outputs of a call whose policy records become `after`, computed before the commit
// (so nothing allocates once the ledger changed): edge_flows = the policy's own allocation after
// the call, residual = the ledger's after it (the totals' service baseline - carried + own before -
// own after), saturated_boundary_edges on that residual, edge_reservations = the policy's own holds
// after the call.
// `totals`: the totals the call read (nullptr: the ledger's, read only for these two outputs).
void prepare_outputs(const PolicyEngine& pe, const LedgerTotals* totals, const OutputOptions& outputs,
                     std::span<const FlowRecord* const> after, PlacementResult& out) {
  const auto m = static_cast<std::size_t>(pe.graph->num_edges());
  if (outputs.edge_reservations) {
    std::vector<Amount> held(m, 0);
    for (const FlowRecord* r : after) {
      for (const auto& [e, x] : r->edge_reservations) held[static_cast<std::size_t>(e)] += x;
    }
    out.edge_reservations = std::move(held);
  }
  if (!outputs.edge_flows && !outputs.residual && !outputs.saturated_boundary_edges) return;
  std::vector<Amount> own(m, 0);
  for (const FlowRecord* r : after) {
    for (const auto& [e, x] : r->edge_flows) own[static_cast<std::size_t>(e)] += x;
  }
  if (outputs.residual || outputs.saturated_boundary_edges) {
    const LedgerTotals t = totals != nullptr ? *totals : totals_of(pe.ledger);
    std::vector<Amount> res(m);
    for (std::size_t e = 0; e < m; ++e) res[e] = t.service_baseline[e] - t.carried[e];
    for (const auto& cohort : pe.epoch.record_groups) {
      for (const auto& s : cohort) {
        for (const auto& [e, x] : LedgerAccess::live(pe.ledger, s.record_counter).edge_flows) res[static_cast<std::size_t>(e)] += x;
      }
    }
    for (std::size_t e = 0; e < m; ++e) res[e] -= own[e];
    if (outputs.saturated_boundary_edges) out.saturated_boundary_edges = saturated_boundary_edges_of(*pe.graph, res, pe.sources, pe.eligibility);
    if (outputs.residual) out.residual = std::move(res);
  }
  if (outputs.edge_flows) out.edge_flows = std::move(own);
}

bool any_output(const OutputOptions& outputs) noexcept {
  return outputs.edge_flows || outputs.residual || outputs.saturated_boundary_edges || outputs.edge_reservations;
}

// ---------------------------------------------------------------- the demand's state

// The result's 128-bit amounts (Counter128 values stay below 2^127).
Counter128 c128(Amount x) noexcept { return Counter128 {0, static_cast<std::uint64_t>(x)}; }

// a += b, or std::invalid_argument("<what> limit exceeded") when the sum reaches 2^127.
void accumulate(Counter128& a, const Counter128& b, const char* what) {
  const std::uint64_t lo = a.lo + b.lo;
  const std::uint64_t hi = a.hi + b.hi + (lo < a.lo ? 1U : 0U);
  if (hi >= (std::uint64_t {1} << 63)) throw std::invalid_argument(std::string(what) + " limit exceeded");
  a = Counter128 {hi, lo};
}

// a - b, for a >= b.
Counter128 minus(const Counter128& a, const Counter128& b) noexcept {
  return Counter128 {a.hi - b.hi - (a.lo < b.lo ? 1U : 0U), a.lo - b.lo};
}

// The demand's state after a call: its records after the call (`after`: one per flow, per roster
// source or per program entry; the record's source names its row), the shares a Pinned head lost
// itself (`lost`: per source, outside every record) and the epoch's counters after it. Lossy: a
// row's source loss is what its source was offered and neither delivered nor dropped; the
// demand's is its offered total minus delivery and drops, which with an empty roster is the
// whole offer while every row stays zero.
void fill_state(const PolicyEngine& pe, std::span<const FlowRecord* const> after, const Counter128& requested,
                bool infinite, const Counter128& offered, PlacementResult& out,
                std::span<const std::pair<NodeId, Counter128>> lost = {}) {
  const bool lossy = pe.config.admission == Admission::Lossy;
  if (!infinite) out.requested = requested;
  out.by_source.resize(pe.sources.size());
  for (std::size_t i = 0; i < pe.sources.size(); ++i) {
    out.by_source[i].source = pe.sources[i];
    if (lossy) out.by_source[i].offered = Counter128 {};
  }
  // The records' drops, gathered for one merge (one record's are already ascending and unique).
  std::vector<std::pair<EdgeId, Counter128>> drops;
  std::size_t dropping = 0;
  if (lossy) {
    std::size_t n = 0;
    for (const FlowRecord* r : after) {
      n += r->edge_drops.size();
      dropping += r->edge_drops.empty() ? 0 : 1;
    }
    drops.reserve(n);
  }
  for (const FlowRecord* r : after) {
    const auto it = std::lower_bound(pe.sources.begin(), pe.sources.end(), r->sources.front());
    auto& row = out.by_source[static_cast<std::size_t>(it - pe.sources.begin())];
    row.delivered += r->delivered;
    out.delivered += r->delivered;
    if (!lossy) continue;
    accumulate(*row.offered, r->offered.value_or(Counter128 {}), "offered");
    for (const auto& [e, x] : r->edge_drops) {
      accumulate(row.total_edge_drops, x, "drops");
      drops.emplace_back(e, x);
    }
  }
  for (const auto& [s, x] : lost) {
    const auto it = std::lower_bound(pe.sources.begin(), pe.sources.end(), s);
    accumulate(*out.by_source[static_cast<std::size_t>(it - pe.sources.begin())].offered, x, "offered");
  }
  if (!lossy) {
    for (auto& row : out.by_source) row.admitted = c128(row.delivered);
    out.admitted = c128(out.delivered);
    if (!infinite) out.unadmitted = minus(requested, out.admitted);
    return;
  }
  for (auto& row : out.by_source) {
    row.admitted = *row.offered;
    row.source_drops = minus(minus(*row.offered, c128(row.delivered)), row.total_edge_drops);
  }
  if (dropping > 1) std::sort(drops.begin(), drops.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  out.edge_drops.reserve(drops.size());
  Counter128 dropped {};
  for (const auto& [e, x] : drops) {
    if (out.edge_drops.empty() || out.edge_drops.back().first != e) out.edge_drops.emplace_back(e, Counter128 {});
    accumulate(out.edge_drops.back().second, x, "drops");
    accumulate(dropped, x, "drops");
  }
  out.offered = offered;
  out.admitted = offered;
  out.unadmitted = Counter128 {};
  out.source_drops = minus(minus(offered, c128(out.delivered)), dropped);
}

// The exact sum of cost * flow over the records' edges (a controller's total cost; below 2^124:
// the costs and the capacities each sum below 2^62).
Counter128 total_cost_of(const Graph& graph, std::span<const FlowRecord* const> after) {
  const auto cost = graph.edge_costs();
  U128 sum {};
  for (const FlowRecord* r : after) {
    for (const auto& [e, x] : r->edge_flows) {
      const U128 p = mul_u64(static_cast<std::uint64_t>(cost[static_cast<std::size_t>(e)]), static_cast<std::uint64_t>(x));
      const std::uint64_t lo = sum.lo + p.lo;
      sum = U128 {sum.hi + p.hi + (lo < sum.lo ? 1U : 0U), lo};
    }
  }
  return Counter128 {sum.hi, sum.lo};
}

// The records of the policy's live flows (a call that commits nothing).
std::vector<const FlowRecord*> current_records(const PolicyEngine& pe) {
  std::vector<const FlowRecord*> out;
  for (const auto& cohort : pe.epoch.record_groups) {
    for (const auto& s : cohort) out.push_back(&LedgerAccess::live(pe.ledger, s.record_counter));
  }
  return out;
}

// The drafts' records as a list, in per-thread storage reused by every call (valid until this
// thread's next draft_records; it grows before any commit).
std::span<const FlowRecord* const> draft_records(const std::vector<RecordDraft>& drafts) {
  thread_local std::vector<const FlowRecord*> out;
  out.clear();
  for (const auto& d : drafts) out.push_back(&d.record);
  return out;
}

// T2: CentralizedTE with a tunnel count.
bool is_t2(const Routing& r) noexcept {
  const auto* gl = std::get_if<CentralizedTE>(&r);
  return gl != nullptr && gl->tunnel_count.has_value();
}
bool t2(const PolicyEngine& pe) noexcept { return is_t2(pe.config.routing); }

// The tunnel models: T1 and T2 (one record per program entry; a full-offer call is their only
// computation).
bool tunnel_model(const PolicyEngine& pe) noexcept { return is<DistributedTE>(pe) || t2(pe); }

// CentralizedTE without tunnels (with or without a support), or ShortestPaths Optimized: a controller.
bool uses_optimized_flow(const PolicyEngine& pe) noexcept {
  return (is<CentralizedTE>(pe) && !t2(pe)) || (is<ShortestPaths>(pe) && std::get<ShortestPaths>(pe.config.routing).allocation == RouteAllocation::Optimized);
}

// Unrestricted CentralizedTE: the controller over every eligible edge (dense arrays, the curve, the cut).
bool is_unrestricted_centralized_te(const PolicyEngine& pe) noexcept {
  return uses_optimized_flow(pe) && is<CentralizedTE>(pe) && !std::get<CentralizedTE>(pe.config.routing).support;
}

// ================================================================ the transaction

Trace route_trace(std::span<const EdgeId> edges) { return Trace::ints(edges); }

std::string key(NodeId v) { return std::to_string(v); }
// ================================================================ calls

// What a call computes. place(v) continues the epoch `cur` (the policy's); evaluate, replace and
// forward start from a fresh epoch (forward keeps the installed program) on the released view.
enum class PolicyOperation { Place, Evaluate, Replace, Forward };

struct PolicyCall {
  PolicyOperation operation {PolicyOperation::Place};
  std::optional<Amount> offer {};          // nullopt: +inf
  std::optional<std::optional<Amount>> design {};  // evaluate and replace: absent, or a design (nullopt: +inf)
  const Eligibility* black_holes {nullptr};   // forward's black holes
  const OutputOptions* outputs {nullptr};
  Trace* trace {nullptr};
  const LedgerTotals* totals {nullptr};  // the totals the call reads (nullptr: the ledger's)
  [[nodiscard]] bool is_full_offer() const noexcept { return operation != PolicyOperation::Place; }
};

LedgerTotals totals_of(const PolicyEngine& pe, const PolicyCall& rq) { return rq.totals != nullptr ? *rq.totals : totals_of(pe.ledger); }

// A computed call, published or discarded whole: its result, the drafts of its commit (a place:
// every live record plus the new ones; a full-offer call: only new ones, the live ones retired),
// the epoch after it, and per new draft (in draft order) its flow in that epoch.
struct PendingPolicyUpdate {
  PlacementResult result {};
  std::vector<RecordDraft> drafts {};
  PolicyEpoch next_epoch {};
  std::vector<std::pair<std::size_t, std::size_t>> fresh {};
};

// The requested total after the call (a full-offer call: its offer alone) and whether it is
// infinite.
Counter128 requested_after(const PolicyEpoch& cur, std::optional<Amount> v) noexcept {
  Counter128 r = cur.requested;
  if (v.has_value()) (void)r.add(*v);  // checked at the call's entry
  return r;
}

// The epoch a full-offer call starts: one empty cohort per roster source (or program entry);
// forward keeps the installed program and ShortestPaths Optimized's support.
PolicyEpoch fresh_epoch(const PolicyEngine& pe, bool keep_program) {
  PolicyEpoch e;
  e.record_groups.resize(pe.epoch.record_groups.size());
  if (keep_program) {
    e.program = pe.epoch.program;
    e.plan = pe.epoch.plan;
    e.support = pe.epoch.support;
  }
  return e;
}

// The termination of a forwarding call (Termination's rule after the empty roster): rule (3) when
// the program has an entry with a path (`has_path`) and no head can use one (`active` false).
Termination forwarding_termination(bool active, bool has_path, bool lossy, std::optional<Amount> v,
                                   Amount call_delivered) noexcept {
  if (has_path && !active) return Termination::NoActiveWeight;
  if (lossy) return v.has_value() ? Termination::Complete : Termination::Exhausted;
  return v.has_value() && call_delivered == *v ? Termination::Complete : Termination::Exhausted;
}

// ---------------------------------------------------------------- the controllers' paths

// decompose() (program_data.hpp): each node's out-edge cursor only advances, remaining flow never
// increasing.
// decompose()'s per-thread scratch: per node its cursor and its range's end (the positions of its
// positive out-edges; cursor -1: none), the destination marks, its roster slot (-1: none), the
// nodes whose cursor is set, the walk and a source's path pieces. A call resets what it set (also when it throws), so it costs O(the nodes
// it walks and the flow's positions), not O(N).
struct DecomposeScratch {
  std::vector<std::int64_t> cur;
  std::vector<std::size_t> end;
  std::vector<std::uint8_t> is_dst;
  std::vector<std::int32_t> slot;
  std::vector<NodeId> touched;
  std::vector<std::size_t> walk;
  std::vector<std::pair<EdgeId, Amount>> acc;
};

} // namespace

Decomposition decompose(const Graph& graph, std::span<const EdgeId> edges, std::vector<Amount> rem,
                        std::span<const NodeId> reachable_sources, std::span<const NodeId> destinations) {
  const bool dense = edges.empty() && rem.size() == static_cast<std::size_t>(graph.num_edges());
  const auto edge_sources = graph.edge_sources();
  const auto edge_destinations = graph.edge_destinations();
  // Position p of the flows: EdgeId p (dense) or edges[p]. Out-edges have consecutive EdgeIds, so
  // a node's positions are one range.
  auto edge_at = [&](std::size_t p) { return dense ? static_cast<EdgeId>(p) : edges[p]; };
  thread_local DecomposeScratch sc;
  const auto n = static_cast<std::size_t>(graph.num_nodes());
  // Each grown on its own (a failed growth leaves the others' sizes as they were).
  if (sc.cur.size() < n) sc.cur.resize(n, -1);
  if (sc.end.size() < n) sc.end.resize(n, 0);
  if (sc.is_dst.size() < n) sc.is_dst.resize(n, 0);
  if (sc.slot.size() < n) sc.slot.resize(n, -1);
  struct Reset {
    DecomposeScratch& sc;
    std::span<const NodeId> dests;
    std::span<const NodeId> roster;
    ~Reset() {
      for (NodeId u : sc.touched) sc.cur[static_cast<std::size_t>(u)] = -1;
      sc.touched.clear();
      for (NodeId d : dests) sc.is_dst[static_cast<std::size_t>(d)] = 0;
      for (NodeId s : roster) sc.slot[static_cast<std::size_t>(s)] = -1;
    }
  } reset {sc, destinations, reachable_sources};
  for (NodeId d : destinations) sc.is_dst[static_cast<std::size_t>(d)] = 1;
  for (std::size_t i = 0; i < reachable_sources.size(); ++i) sc.slot[static_cast<std::size_t>(reachable_sources[i])] = static_cast<std::int32_t>(i);
  // A node without positive out-edges has an empty range.
  auto cursor = [&](NodeId u) -> std::int64_t& {
    auto& c = sc.cur[static_cast<std::size_t>(u)];
    if (c < 0) {
      sc.touched.push_back(u);
      c = 0;
      sc.end[static_cast<std::size_t>(u)] = 0;
    }
    return c;
  };
  Decomposition out;
  out.shares.resize(reachable_sources.size());
  out.supplies.resize(reachable_sources.size());
  // One pass over the flow's positions before any walk: each node's range of positive out-edges
  // (the zero positions around them are skipped by the walks anyway); each source's supply, its
  // net outflow (transit through it from another source is not supply; a walk through another
  // source lowers its in- and outflow alike, as sources are never destinations); the number of
  // positive positions.
  std::size_t positive = 0;
  for (std::size_t p = 0; p < rem.size(); ++p) {
    if (rem[p] == 0) continue;
    ++positive;
    const auto e = static_cast<std::size_t>(edge_at(p));
    const auto u = static_cast<std::size_t>(edge_sources[e]);
    if (sc.cur[u] < 0) {
      sc.touched.push_back(edge_sources[e]);
      sc.cur[u] = static_cast<std::int64_t>(p);
    }
    sc.end[u] = p + 1;
    const auto from = sc.slot[u];
    const auto to = sc.slot[static_cast<std::size_t>(edge_destinations[e])];
    if (from >= 0) out.supplies[static_cast<std::size_t>(from)] += rem[p];
    if (to >= 0) out.supplies[static_cast<std::size_t>(to)] -= rem[p];
  }
  // A one-source flow is all that source's share: the walks cover the flow (checked below), so the
  // sum of its paths is the flow's positive positions, ascending (the per-path sum's degenerate
  // form).
  const bool one = reachable_sources.size() == 1;
  // Every path zeroes a position or ends a source's supply.
  out.paths.reserve(positive + reachable_sources.size());
  if (one) {
    out.shares[0].reserve(positive);
    for (std::size_t p = 0; p < rem.size(); ++p) {
      if (rem[p] != 0) out.shares[0].emplace_back(edge_at(p), rem[p]);
    }
  }
  auto& walk = sc.walk;
  auto& acc = sc.acc;
  for (std::size_t i = 0; i < reachable_sources.size(); ++i) {
    const NodeId s = reachable_sources[i];
    Amount left = out.supplies[i];
    acc.clear();
    while (left > 0) {
      walk.clear();
      NodeId u = s;
      while (walk.empty() || sc.is_dst[static_cast<std::size_t>(u)] == 0) {
        auto& p = cursor(u);
        const auto stop = sc.end[static_cast<std::size_t>(u)];
        while (static_cast<std::size_t>(p) < stop && rem[static_cast<std::size_t>(p)] == 0) ++p;
        if (static_cast<std::size_t>(p) == stop) {
          throw std::logic_error("CentralizedTE decomposition: flow conservation broken (no continuation)");
        }
        walk.push_back(static_cast<std::size_t>(p));
        u = edge_destinations[static_cast<std::size_t>(edge_at(static_cast<std::size_t>(p)))];
        if (walk.size() > n) throw std::logic_error("CentralizedTE decomposition: a cycle with positive flow");
      }
      Amount amount = left;
      for (std::size_t p : walk) amount = std::min(amount, rem[p]);
      OptimizedFlowPath path;
      path.roster_index = static_cast<std::int32_t>(i);
      path.weight = amount;
      path.edges.reserve(walk.size());
      for (std::size_t p : walk) {
        rem[p] -= amount;
        path.edges.push_back(edge_at(p));
        if (!one) acc.emplace_back(edge_at(p), amount);
      }
      out.paths.push_back(std::move(path));
      left -= amount;
    }
    if (one) continue;
    std::sort(acc.begin(), acc.end());
    auto& res = out.shares[i];
    for (const auto& [e, x] : acc) {
      if (!res.empty() && res.back().first == e) {
        res.back().second += x;
      } else {
        res.emplace_back(e, x);
      }
    }
  }
  for (Amount x : rem) {
    if (x != 0) throw std::logic_error("CentralizedTE decomposition: the walks do not cover the flow");
  }
  return out;
}

namespace {

// ---------------------------------------------------------------- ShortestPaths preparation

// A ShortestPaths demand's routes for the epoch, prepared at its first positive call (or a replace)
// on the working residual at that call's entry: per roster source its route (the record's route
// and route cost), and its program (Equal, Weighted: one split route per source in one
// plan rooted at the roster) or the routes' union (Optimized, the controller's support).
struct PreparedShortestPaths {
  std::vector<LiveFlow> flows {};  // per roster source, new (id -1)
  std::shared_ptr<const ForwardingProgram> program {};
  std::vector<EdgeId> edges {};    // ascending: the union of the routes
};

PreparedShortestPaths prepare_shortest_path_program(PolicyCallContext& c) {
  const auto& pe = c.engine;
  const auto& g = c.graph;
  const std::size_t k = pe.reachable_sources.size();
  const RouteAllocation allocation = std::get<ShortestPaths>(pe.config.routing).allocation;
  PreparedShortestPaths out;
  out.flows.resize(k);
  std::vector<std::vector<EdgeId>> routes(k);
  std::shared_ptr<const PreparedRoute> single;  // one source's plan from its one-pass creation
  for (std::size_t i = 0; i < k; ++i) {
    const NodeId s = pe.reachable_sources[i];
    auto& f = out.flows[i];
    f.source = s;
    if (auto made = creation_products(c, s)) {
      f.route = std::move(made->route);
      f.min_route_cost = made->min_route_cost;
      if (k == 1 && allocation != RouteAllocation::Optimized) {
        single = std::make_shared<const PreparedRoute>(std::move(made->prepared_route));
      } else {
        routes[i] = std::move(made->prepared_route.support_edges);
      }
    } else {
      RouteSupport r = shortest_path_support(c, s);
      f.route = record_route_of(g, r.edges, pe.support);
      f.min_route_cost = r.cost != kUnreachable ? r.cost : min_route_cost(g, *f.route, s, pe.allowed_destinations);
      routes[i] = std::move(r.edges);
    }
  }
  // One split-routing source: its one-pass creation plan (the pooled builder's plan of its route,
  // rooted at it) is the program; its route is the union.
  if (single) {
    std::vector<ProgramEntry> entries(1);
    entries[0].route = out.flows[0].route;
    entries[0].min_route_cost = out.flows[0].min_route_cost;
    out.program = shortest_path_program(g, pe.sources, pe.destinations, pe.reachable_sources, std::move(single), std::move(entries));
    return out;
  }
  for (const auto& r : routes) out.edges.insert(out.edges.end(), r.begin(), r.end());
  std::sort(out.edges.begin(), out.edges.end());
  out.edges.erase(std::unique(out.edges.begin(), out.edges.end()), out.edges.end());
  // One continuation per router: every source's route leaves each of its routers by all the
  // union's edges there. The routes are read on one view under one mask, so a router's
  // continuation depends on the router alone (a single path's suffix is the canonical path of
  // its first node; the rejected multipath = false with a HIGHER_* tie_break is the exception).
  // One route is its own union: nothing to check.
  const auto row = g.outgoing_offsets();
  const auto edge_sources = g.edge_sources();
  for (std::size_t i = 0; k > 1 && i < k; ++i) {
    const auto& r = routes[i];
    for (std::size_t a = 0; a < r.size();) {
      const auto u = static_cast<std::size_t>(edge_sources[static_cast<std::size_t>(r[a])]);
      std::size_t b = a;
      while (b < r.size() && static_cast<std::size_t>(edge_sources[static_cast<std::size_t>(r[b])]) == u) ++b;
      const auto lo = std::lower_bound(out.edges.begin(), out.edges.end(), row[u]);
      const auto hi = std::lower_bound(lo, out.edges.end(), row[u + 1]);
      if (static_cast<std::size_t>(hi - lo) != b - a) {
        throw std::logic_error("ShortestPaths: two roster sources continue differently at node " + std::to_string(u));
      }
      a = b;
    }
  }
  if (allocation == RouteAllocation::Optimized) return out;
  const auto support = PredecessorDAGBuilder::ascending(*pe.support);
  const auto dag = std::equal(out.edges.begin(), out.edges.end(), support.begin(), support.end())
                       ? pe.support
                       : std::make_shared<const PredecessorDAG>(PredecessorDAGBuilder::of_edges(g, out.edges));
  // Weighted weights: the static source, or the free capacity read now (ResidualSnapshot, frozen
  // for the epoch).
  std::vector<Amount> creation;
  WeightSource weights = shortest_path_weight_source(pe);
  if (allocation == RouteAllocation::Weighted && weights.basis == WeightBasis::ResidualSnapshot) {
    creation.resize(out.edges.size());
    for (std::size_t i = 0; i < out.edges.size(); ++i) creation[i] = c.service_capacity.at(out.edges[i]);
    weights.values = creation;
  }
  std::vector<ProgramEntry> entries(k);
  for (std::size_t i = 0; i < k; ++i) {
    entries[i].route = out.flows[i].route;
    entries[i].min_route_cost = out.flows[i].min_route_cost;
  }
  out.program = pooled_dag_program(g, pe.sources, pe.destinations, pe.reachable_sources, *dag, allocation, weights,
                                   std::move(entries));
  return out;
}

Trace routes_trace(const PreparedShortestPaths& p) {
  Trace t = Trace::map();
  for (const auto& f : p.flows) t.set(key(f.source), route_trace(PredecessorDAGBuilder::ascending(*f.route)));
  return t;
}

// A record of one roster source of ShortestPaths or CentralizedTE (its route, its route cost).
FlowRecord source_record(const PolicyEngine& pe, const LiveFlow& f) {
  FlowRecord r;
  r.sources = {f.source};
  r.destinations = pe.allowed_destinations;
  r.admission = pe.config.admission;
  r.route = f.route;
  r.min_route_cost = f.min_route_cost;
  return r;
}

// `base` plus `add` (both ascending by EdgeId, positive entries): the record's cumulative edge
// amounts (Amount) or drops (Counter128).
std::vector<std::pair<EdgeId, Amount>> plus(const std::vector<std::pair<EdgeId, Amount>>& base,
                                            const std::vector<std::pair<EdgeId, Amount>>& add) {
  std::vector<std::pair<EdgeId, Amount>> out;
  out.reserve(base.size() + add.size());
  std::size_t i = 0;
  for (const auto& [e, x] : add) {
    while (i < base.size() && base[i].first < e) out.push_back(base[i++]);
    if (i < base.size() && base[i].first == e) {
      out.emplace_back(e, base[i++].second + x);
    } else {
      out.emplace_back(e, x);
    }
  }
  out.insert(out.end(), base.begin() + static_cast<std::ptrdiff_t>(i), base.end());
  return out;
}

std::vector<std::pair<EdgeId, Counter128>> plus(const std::vector<std::pair<EdgeId, Counter128>>& base,
                                                const std::vector<std::pair<EdgeId, Amount>>& add) {
  std::vector<std::pair<EdgeId, Counter128>> out;
  out.reserve(base.size() + add.size());
  std::size_t i = 0;
  for (const auto& [e, x] : add) {
    while (i < base.size() && base[i].first < e) out.push_back(base[i++]);
    Counter128 sum = i < base.size() && base[i].first == e ? base[i++].second : Counter128 {};
    accumulate(sum, c128(x), "drops");
    out.emplace_back(e, sum);
  }
  out.insert(out.end(), base.begin() + static_cast<std::ptrdiff_t>(i), base.end());
  return out;
}

// The positive entries of a member-aligned array as (EdgeId, amount) pairs.
std::vector<std::pair<EdgeId, Amount>> positive(std::span<const EdgeId> members, std::span<const Amount> a) {
  std::vector<std::pair<EdgeId, Amount>> out;
  out.reserve(static_cast<std::size_t>(std::count_if(a.begin(), a.end(), [](Amount x) { return x != 0; })));
  for (std::size_t m = 0; m < members.size(); ++m) {
    if (a[m] != 0) out.emplace_back(members[m], a[m]);
  }
  return out;
}

// The sum of several ascending (EdgeId, amount) lists.
std::vector<std::pair<EdgeId, Amount>> merged(std::vector<std::pair<EdgeId, Amount>> all) {
  std::sort(all.begin(), all.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  std::vector<std::pair<EdgeId, Amount>> out;
  out.reserve(all.size());
  for (const auto& [e, x] : all) {
    if (!out.empty() && out.back().first == e) {
      out.back().second += x;
    } else {
      out.emplace_back(e, x);
    }
  }
  return out;
}

// One record per roster source (ShortestPaths, CentralizedTE): the drafts in roster order, kept ids when
// the epoch already has them, else new (with their place in the next epoch).
void source_drafts(const PolicyEngine& pe, const std::vector<LiveFlow>& flows, PendingPolicyUpdate& s) {
  const std::size_t k = flows.size();
  s.drafts.resize(k);
  s.next_epoch.record_groups.assign(k, {});
  for (std::size_t i = 0; i < k; ++i) {
    if (flows[i].record_counter >= 0) {
      s.drafts[i].record_counter = flows[i].record_counter;
    } else {
      s.fresh.emplace_back(i, 0);
    }
    s.drafts[i].record = source_record(pe, flows[i]);
    s.next_epoch.record_groups[i].push_back(flows[i]);
  }
}

// The programs of a call and the flags they give. Only T1's own setup books holds, so only it
// reports a reservation change; Pinned over a T1 program keeps the imported amounts and books none.
void set_flags(PlacementResult& out, const ForwardingProgram* before, const ForwardingProgram* after,
               bool reserves) noexcept {
  const auto change = compare_programs(before, after, reserves);
  out.routes_changed = change.routes;
  out.weights_changed = change.weights;
  out.reservations_changed = change.reservations;
}

// The by_tunnel rows of a controller's program: every path up, assigned = delivered = its amount.
std::vector<TunnelFlow> optimized_flow_tunnel_rows(const ForwardingProgram& program) {
  std::vector<TunnelFlow> rows;
  rows.reserve(program.tunnels().size());
  for (const auto& t : program.tunnels()) rows.push_back(TunnelFlow {t.entry_id, t.head, TunnelStatus::Up, c128(t.weight), t.weight, {}});
  return rows;
}

// ---------------------------------------------------------------- ShortestPaths Equal and Weighted

PendingPolicyUpdate place_pooled(const PolicyEngine& pe, const PolicyEpoch& cur, ServiceCapacityWorkspace& W, const PolicyCall& rq,
                    const ForwardingProgram* before_program) {
  const auto& g = *pe.graph;
  const auto v = rq.offer;
  Trace* trace = rq.trace;
  PolicyCallContext c {pe, g, W};
  const std::size_t k = pe.reachable_sources.size();
  const bool lossy = pe.config.admission == Admission::Lossy;
  const bool created = cur.record_groups.front().empty();
  PreparedShortestPaths prep;
  std::vector<LiveFlow> flows(k);
  if (created) {
    prep = prepare_shortest_path_program(c);
    flows = prep.flows;
  } else {
    for (std::size_t i = 0; i < k; ++i) flows[i] = cur.record_groups[i].front();
  }
  const auto program = created ? prep.program : cur.program;
  const PreparedRoute& plan = *ProgramAccess::data(*program).split_contexts.front().prepared_route;
  const auto members = std::span<const EdgeId>(plan.forwarding_edges);
  std::vector<const FlowRecord*> before(k, nullptr);
  if (!created) {
    for (std::size_t i = 0; i < k; ++i) before[i] = &LedgerAccess::live(pe.ledger, flows[i].record_counter);
  }
  PendingPolicyUpdate s;
  source_drafts(pe, flows, s);
  auto& drafts = s.drafts;
  // A one-source program is its own attribution.
  Attribution attr;
  Attribution* per_source = k > 1 ? &attr : nullptr;
  std::vector<Amount> inject(k, 0);
  Forwarding fw;
  Amount call_delivered = 0;
  Counter128 offered = cur.offered;
  Trace step_tr = Trace::map();
  if (!lossy) {
    // The largest total in [admitted, admitted + v] whose no-drop forwarding fits V = W + own.
    std::vector<Amount> own(members.size(), 0);
    std::vector<Amount> part;
    Amount admitted = 0;
    for (const FlowRecord* r : before) {
      if (r == nullptr) continue;
      admitted += r->delivered;
      align(members, r->edge_flows, part);
      for (std::size_t m = 0; m < members.size(); ++m) own[m] += part[m];
    }
    std::vector<Amount> view;
    W.gather(members, view);
    for (std::size_t m = 0; m < members.size(); ++m) view[m] += own[m];
    auto step = lossless_step(plan, admitted, view, own, v);
    call_delivered = step.admitted - admitted;
    if (per_source == nullptr) {
      // One source holds the whole admitted allocation (the attribution's degenerate form): it
      // delivers all of it, with no second forwarding.
      auto& r = drafts.front().record;
      r.edge_flows = positive(members, step.edge_flows);
      r.delivered = step.admitted;
    } else {
      // Several sources: forward the admitted total again to attribute it.
      round_robin(step.admitted, inject);
      fw = forward(plan, Injection {inject, false}, {}, per_source);
      if (fw.edge_flows != step.edge_flows) {
        throw std::logic_error("ShortestPaths: the attributed forwarding differs from the admitted allocation");
      }
      for (std::size_t i = 0; i < k; ++i) {
        auto& r = drafts[i].record;
        for (const auto& x : attr.entries[i]) {
          if (x.carried != 0) r.edge_flows.emplace_back(members[static_cast<std::size_t>(x.edge_position)], x.carried);
        }
        r.delivered = attr.delivered[i];
      }
    }
    if (kTraced && trace) step_tr = Trace::ints(std::vector<Amount> {admitted, step.admitted});
  } else {
    // The new round-robin injections of the offer (or the first-hop experiment) on W.
    Injection in;
    if (v.has_value()) {
      round_robin_increments(cur.offered, *v, inject);
      in.per_root = inject;
    } else {
      in.first_hop = true;
    }
    std::vector<Amount> w;
    W.gather(members, w);
    fw = forward(plan, in, w, per_source);
    call_delivered = fw.delivered;
    for (std::size_t i = 0; i < k; ++i) {
      auto& r = drafts[i].record;
      std::vector<std::pair<EdgeId, Amount>> carried;
      std::vector<std::pair<EdgeId, Amount>> dropped;
      Amount offer = fw.offered;
      Amount delivered = fw.delivered;
      if (per_source == nullptr) {
        carried = positive(members, fw.edge_flows);
        dropped = positive(members, fw.edge_drops);
      } else {
        for (const auto& x : attr.entries[i]) {
          const EdgeId e = members[static_cast<std::size_t>(x.edge_position)];
          if (x.carried != 0) carried.emplace_back(e, x.carried);
          if (x.drops != 0) dropped.emplace_back(e, x.drops);
        }
        offer = attr.offered[i];
        delivered = attr.delivered[i];
      }
      // The counters' 2^127 checks: a drop counter first, then the offered ones.
      const FlowRecord* b = before[i];
      r.edge_drops = plus(b != nullptr ? b->edge_drops : std::vector<std::pair<EdgeId, Counter128>> {}, dropped);
      r.offered = b != nullptr ? b->offered.value_or(Counter128 {}) : Counter128 {};
      accumulate(*r.offered, c128(offer), "offered");
      r.delivered = (b != nullptr ? b->delivered : 0) + delivered;
      r.edge_flows = b != nullptr ? plus(b->edge_flows, carried) : std::move(carried);
    }
    accumulate(offered, c128(fw.offered), "offered");
    if (kTraced && trace) {
      step_tr = Trace::map();
      for (std::size_t i = 0; i < k; ++i) {
        step_tr.set(key(pe.reachable_sources[i]), v.has_value() ? Trace::integer(inject[i]) : Trace::str("inf"));
      }
    }
  }
  const Counter128 requested = requested_after(cur, v);
  PlacementResult& out = s.result;
  out.call_delivered = call_delivered;
  if (v.has_value()) out.call_undelivered = *v - call_delivered;
  out.delivery_delta = call_delivered;
  if (lossy) out.call_edge_drops = positive(members, fw.edge_drops);
  out.call_source_drops = fw.source_drops;
  if (members.empty()) {
    out.termination = Termination::NoActiveWeight;
  } else if (lossy) {
    out.termination = v.has_value() ? Termination::Complete : Termination::Exhausted;
  } else {
    out.termination = v.has_value() && call_delivered == *v ? Termination::Complete : Termination::Exhausted;
  }
  if (created) set_flags(out, before_program, program.get(), false);
  out.program = program;
  fill_state(pe, draft_records(drafts), requested, cur.infinite || !v.has_value(), offered, out);
  if (kTraced && trace) {
    trace->set("prepare", created ? routes_trace(prep) : Trace::null());
    trace->set(lossy ? "offers" : "T", std::move(step_tr));
  }
  s.next_epoch.program = program;
  s.next_epoch.requested = requested;
  s.next_epoch.offered = offered;
  s.next_epoch.prepared = true;
  s.next_epoch.infinite = cur.infinite || !v.has_value();
  return s;
}

// ---------------------------------------------------------------- the controllers

// One solve of a controller from zero on the view V with value cap `cap` (nullopt: uncapped);
// `caps` (the equal source allocation) caps source i's arc at caps[i]. Unrestricted CentralizedTE solves over every
// eligible edge, its V and flows per EdgeId; CentralizedTE with a support and ShortestPaths Optimized over
// `support`, their V and flows aligned with it.
struct OptimizedFlowResult {
  Amount flow_value {0};
  std::vector<Amount> edge_flows {};
  MinCostFlowResult payload {};  // unrestricted CentralizedTE: the flow core's payload (the curve, the cut, the trace's values)
};

OptimizedFlowResult solve(const PolicyEngine& pe, std::span<const EdgeId> support, std::span<const Amount> V,
            std::optional<Amount> cap, std::span<const Amount> caps) {
  const auto& g = *pe.graph;
  OptimizedFlowResult out;
  if (is_unrestricted_centralized_te(pe)) {
    out.payload = solve_min_cost_flow(g, V, pe.reachable_sources, pe.allowed_destinations, pe.eligibility, cap, caps.empty() ? pe.capture : nullptr, caps);
    out.flow_value = out.payload.flow_value;
    out.edge_flows = out.payload.edge_flows;
    return out;
  }
  auto r = solve_subset(g, support, V, pe.reachable_sources, pe.allowed_destinations, cap, caps);
  out.flow_value = r.flow_value;
  out.edge_flows = std::move(r.edge_flows);
  return out;
}

// The equal source allocation: the largest total T in [lo, hi] whose round-robin supplies RR(T, H, i) are all
// carried (a solve capped at T with source i's arc capped at its share reaches T), given that lo's
// are. Feasible totals form a prefix (each share never decreases in T, and a smaller supply
// vector stays feasible), so integer bisection is exact; T's solve from zero is returned, so the
// result does not depend on the probes.
OptimizedFlowResult solve_equal(const PolicyEngine& pe, std::span<const EdgeId> support, std::span<const Amount> V, Amount lo,
                  Amount hi, Trace* probes) {
  std::vector<Amount> caps(pe.reachable_sources.size());
  auto attempt = [&](Amount t) -> std::optional<OptimizedFlowResult> {
    std::fill(caps.begin(), caps.end(), 0);
    round_robin(t, caps);
    auto s = solve(pe, support, V, t, caps);
    if (kTraced && probes != nullptr) probes->push(Trace::integer(t));
    if (s.flow_value != t) return std::nullopt;
    return s;
  };
  if (auto s = attempt(hi)) return std::move(*s);
  std::optional<OptimizedFlowResult> best;
  while (hi - lo > 1) {
    const Amount mid = lo + (hi - lo) / 2;
    if (auto s = attempt(mid)) {
      lo = mid;
      best = std::move(s);
    } else {
      hi = mid;
    }
  }
  if (!best) best = attempt(lo);
  if (!best) throw std::logic_error("equal source allocation: the admitted total is no longer carried");
  return std::move(*best);
}

// The upper end of the equal search: the view summed over the edges leaving the roster sources
// (unrestricted CentralizedTE: the eligible ones, V per EdgeId; otherwise the support's, V aligned with
// it); the supplies leave through them.
Amount source_outgoing_capacity_bound(const PolicyEngine& pe, std::span<const EdgeId> support, std::span<const Amount> V) {
  const auto row = pe.graph->outgoing_offsets();
  const bool dense = is_unrestricted_centralized_te(pe);
  Amount sum = 0;
  for (NodeId s : pe.reachable_sources) {
    const auto si = static_cast<std::size_t>(s);
    if (dense) {
      for (auto e = row[si]; e < row[si + 1]; ++e) {
        if (pe.eligibility == nullptr || pe.eligibility->allows_edge(e)) sum += V[static_cast<std::size_t>(e)];
      }
    } else {
      for (auto it = std::lower_bound(support.begin(), support.end(), row[si]);
           it != support.end() && *it < row[si + 1]; ++it) {
        sum += V[static_cast<std::size_t>(it - support.begin())];
      }
    }
  }
  return sum;
}

// The program basis of a controller: unrestricted CentralizedTE the topology under its mask; CentralizedTE with
// a support that support; ShortestPaths Optimized each source's route.
ReplayBasis optimized_flow_replay_basis(const PolicyEngine& pe, const std::vector<LiveFlow>& flows) {
  ReplayBasis b;
  if (is_unrestricted_centralized_te(pe)) {
    b.uses_topology_basis = true;
    b.excluded_nodes = pe.excluded_nodes;
    b.excluded_edges = pe.excluded_edges;
  } else if (is<CentralizedTE>(pe)) {
    b.common_support = pe.centralized_te_support;
  } else {
    b.source_supports.reserve(flows.size());
    for (const auto& f : flows) b.source_supports.push_back(f.route);
  }
  return b;
}

PendingPolicyUpdate place_optimized_flow(const PolicyEngine& pe, const PolicyEpoch& cur, ServiceCapacityWorkspace& W, const PolicyCall& rq,
                        const ForwardingProgram* before_program) {
  const auto& g = *pe.graph;
  const auto v = rq.offer;
  Trace* trace = rq.trace;
  const bool dense = is_unrestricted_centralized_te(pe);
  const bool global = is<CentralizedTE>(pe);
  const bool equal = pe.source_allocation == SourceAllocation::Equal;
  const std::size_t k = pe.reachable_sources.size();
  const bool created = cur.record_groups.front().empty();
  // ShortestPaths Optimized: the support of the epoch (its routes, prepared at its first call).
  PreparedShortestPaths prep;
  if (!global && created) {
    PolicyCallContext c {pe, g, W};
    prep = prepare_shortest_path_program(c);
  }
  const std::span<const EdgeId> support = dense      ? std::span<const EdgeId> {}
                                          : global   ? std::span<const EdgeId>(*pe.centralized_te_support)
                                          : created  ? std::span<const EdgeId>(prep.edges)
                                                     : std::span<const EdgeId>(cur.support);
  std::vector<LiveFlow> flows(k);
  for (std::size_t i = 0; i < k; ++i) {
    if (!created) {
      flows[i] = cur.record_groups[i].front();
    } else if (!global) {
      flows[i] = prep.flows[i];
    } else {
      flows[i].source = pe.reachable_sources[i];
    }
  }
  // V = W + own: the working residual plus the policy's own allocation (its live records, or none
  // in a full-offer call, whose view already holds them released), per EdgeId (unrestricted
  // CentralizedTE) or aligned with the support (a record's edges lie in it, both ascending).
  std::vector<Amount> V;
  if (dense) {
    W.dense(V);
  } else {
    W.gather(support, V);
  }
  Amount P = 0;
  for (const auto& f : flows) {
    if (f.record_counter < 0) continue;
    const FlowRecord& rec = LedgerAccess::live(pe.ledger, f.record_counter);
    P += rec.delivered;
    std::size_t p = 0;
    for (const auto& [e, x] : rec.edge_flows) {
      if (dense) {
        V[static_cast<std::size_t>(e)] += x;
        continue;
      }
      while (p < support.size() && support[p] < e) ++p;
      if (p == support.size() || support[p] != e) throw std::logic_error("controller: a record edge outside the support");
      V[p] += x;
    }
  }
  const std::optional<Amount> cap = v.has_value() ? std::optional<Amount>(P + *v) : std::nullopt;
  if (pe.capture != nullptr) {
    if (equal || !dense) throw std::logic_error("a deletion capture of a controller other than optimized CentralizedTE");
    pe.capture->graph = pe.graph;
  }
  Trace probes = Trace::list();
  OptimizedFlowResult sol;
  if (equal) {
    const Amount top = source_outgoing_capacity_bound(pe, support, V);
    sol = solve_equal(pe, support, V, P, cap.has_value() ? std::min(*cap, top) : top, &probes);
  } else {
    sol = solve(pe, support, V, cap, {});
  }
  if (sol.flow_value < P) throw std::logic_error("controller: a call reduced delivered()");
  // The program: each source's supply cut into paths; each source's share is the sum of its
  // paths.
  // The solved flows move into the decomposition (a trace reads them afterwards: a copy).
  const bool traced = kTraced && trace != nullptr;
  auto dec = decompose(g, support, traced ? sol.edge_flows : std::move(sol.edge_flows), pe.reachable_sources, pe.allowed_destinations);
  PendingPolicyUpdate s;
  // A supplied support, which a share that uses every one of its edges reuses as its route (a
  // PredecessorDAG is fixed by its graph and edge set: the same value).
  const auto* gl = std::get_if<CentralizedTE>(&pe.config.routing);
  const std::shared_ptr<const PredecessorDAG> supplied = gl != nullptr ? gl->support : nullptr;
  for (std::size_t i = 0; i < k; ++i) {
    const NodeId src = pe.reachable_sources[i];
    auto& f = flows[i];
    if (global) {
      // CentralizedTE's record route is the edges its share uses.
      const auto& share = dec.shares[i];
      auto same = [&](std::span<const EdgeId> edges) {
        return std::equal(share.begin(), share.end(), edges.begin(), edges.end(),
                          [](const auto& ex, EdgeId e) { return ex.first == e; });
      };
      if (f.route && same(PredecessorDAGBuilder::ascending(*f.route))) continue;
      if (supplied && same(PredecessorDAGBuilder::ascending(*supplied))) {
        f.route = supplied;
      } else {
        std::vector<EdgeId> route;
        route.reserve(share.size());
        for (const auto& [e, x] : share) route.push_back(e);
        f.route = std::make_shared<const PredecessorDAG>(PredecessorDAGBuilder::of_edges(g, route));
      }
      f.min_route_cost = share.empty() ? kUnreachable : min_route_cost(g, *f.route, src, pe.allowed_destinations);
    }
  }
  source_drafts(pe, flows, s);
  for (std::size_t i = 0; i < k; ++i) {
    s.drafts[i].record.delivered = dec.supplies[i];
    s.drafts[i].record.edge_flows = std::move(dec.shares[i]);
  }
  const auto& supplies = dec.supplies;
  auto program = optimized_flow_program(g, pe.sources, pe.destinations, pe.reachable_sources, std::move(dec.paths),
                                    equal ? std::span<const Amount> {} : std::span<const Amount>(supplies),
                                    optimized_flow_replay_basis(pe, flows));
  const Counter128 requested = requested_after(cur, v);
  PlacementResult& out = s.result;
  out.call_delivered = sol.flow_value - P;
  if (v.has_value()) out.call_undelivered = *v - out.call_delivered;
  out.delivery_delta = out.call_delivered;
  out.termination = v.has_value() && out.call_delivered == *v ? Termination::Complete : Termination::Exhausted;
  set_flags(out, before_program, program.get(), false);
  out.program = program;
  out.by_tunnel = optimized_flow_tunnel_rows(*program);
  const auto after = draft_records(s.drafts);
  fill_state(pe, after, requested, cur.infinite || !v.has_value(), Counter128 {}, out);
  out.total_cost = total_cost_of(g, after);
  const bool curve = dense && !equal;
  if (curve) out.marginal_cost_curve = sol.payload.marginal_cost_curve;
  if (curve && rq.outputs != nullptr && rq.outputs->min_cut) out.min_cut = sol.payload.min_cut;
  if (kTraced && trace) {
    if (!global) trace->set("prepare", created ? routes_trace(prep) : Trace::null());
    if (equal) trace->set("probes", std::move(probes));
    if (dense) trace->set("V", Trace::ints(V));
    trace->set("cap", cap.has_value() ? Trace::integer(*cap) : Trace::null());
    Trace fc = Trace::map();
    fc.set("placed", Trace::integer(sol.flow_value));
    if (dense) {
      const auto& r = sol.payload;
      fc.set("cost_distribution", Trace::pairs(r.marginal_cost_curve));
      fc.set("edge_flows", Trace::ints(r.edge_flows));
      fc.set("min_cut", r.min_cut.has_value() ? Trace::ints(*r.min_cut) : Trace::null());
      fc.set("potentials", Trace::ints(r.potentials));
      Trace counters = Trace::map();
      counters.set("dijkstra_runs", Trace::integer(r.counters.dijkstra_runs));
      counters.set("dijkstra_pops", Trace::integer(r.counters.dijkstra_pops));
      counters.set("settled", Trace::integer(r.counters.settled));
      counters.set("arc_scans", Trace::integer(r.counters.arc_scans));
      counters.set("dinic_phases", Trace::integer(r.counters.dinic_phases));
      counters.set("augmentations", Trace::integer(r.counters.augmentations));
      fc.set("counters", std::move(counters));
    } else {
      std::vector<Amount> flows_dense(static_cast<std::size_t>(g.num_edges()), 0);
      for (std::size_t p = 0; p < support.size(); ++p) flows_dense[static_cast<std::size_t>(support[p])] = sol.edge_flows[p];
      fc.set("edge_flows", Trace::ints(flows_dense));
    }
    Trace supply = Trace::map();
    Trace shares_tr = Trace::map();
    for (std::size_t i = 0; i < k; ++i) {
      supply.set(key(pe.reachable_sources[i]), Trace::integer(supplies[i]));
      shares_tr.set(key(pe.reachable_sources[i]), Trace::pairs(s.drafts[i].record.edge_flows));
    }
    fc.set("supply", std::move(supply));
    trace->set("flow_core", std::move(fc));
    trace->set("shares", std::move(shares_tr));
  }
  if (curve) s.next_epoch.marginal_cost_curve = *out.marginal_cost_curve;  // the cumulative curve a later zero call reports
  if (!global) s.next_epoch.support = created ? prep.edges : cur.support;
  s.next_epoch.program = program;
  s.next_epoch.by_tunnel = out.by_tunnel;
  s.next_epoch.requested = requested;
  s.next_epoch.prepared = true;
  s.next_epoch.infinite = cur.infinite || !v.has_value();
  return s;
}

// ---------------------------------------------------------------- forwarding a program

// The by_tunnel rows of a forwarding: per tunnel its status in the view and its amounts
// (assigned, delivered and dropped amounts per program entry).
std::vector<TunnelFlow> tunnel_rows(const ForwardingView& view, std::span<const Counter128> assigned,
                                    std::span<const Amount> delivered, std::span<const Counter128> dropped) {
  const ProgramData& p = *view.data;
  std::vector<TunnelFlow> rows;
  rows.reserve(p.tunnels.size());
  for (std::size_t t = 0; t < p.tunnels.size(); ++t) {
    const auto e = static_cast<std::size_t>(p.tunnel_entry_positions[t]);
    rows.push_back(TunnelFlow {p.tunnels[t].entry_id, p.tunnels[t].head, view.status(e), assigned[e], delivered[e], dropped[e]});
  }
  return rows;
}

// The free capacity of the forwarding's edges in the call's working residual.
AvailableCapacityFn available_capacity_reader(const ServiceCapacityWorkspace& W) {
  return [&W](std::span<const EdgeId> edges, std::span<Amount> out) {
    for (std::size_t i = 0; i < edges.size(); ++i) out[i] = W.at(edges[i]);
  };
}

// Some head of the view can place traffic: a positive injection weight and a surviving entry of
// positive weight whose route has a member.
bool view_active(const ForwardingView& v) {
  const auto& p = *v.data;
  for (std::size_t h = 0; h < v.num_sources(); ++h) {
    if (!v.source_weights().empty() && v.source_weights()[h] == 0) continue;
    for (const std::int32_t ei : v.entries_for_source(h)) {
      const auto e = static_cast<std::size_t>(ei);
      const auto& entry = p.entries[e];
      if (entry.weight == 0) continue;
      if (entry.tunnel_position >= 0) return true;
      const auto r = v.context_source_index(e);
      if (r >= 0 && v.prepared_route(static_cast<std::size_t>(entry.split_context_position))->source_positions[static_cast<std::size_t>(r)] >= 0) return true;
    }
  }
  return false;
}

// A program forwarded as `view` sees it, one record per program entry (Pinned: its replay, every
// call; T1: its own program, after setup or installed): place, evaluate, replace and T1's
// forward. `active`: some head can place traffic (view_active). With `holds` (a T1 policy's own
// program), each up tunnel's record holds its request on every edge of its path; `plan` is the
// result's and the epoch's. `given` (Lossy, a full offer): the forwarding already computed with the
// heads' shares of the offer (a SIMULTANEOUS list's joint evaluation), moved from.
PendingPolicyUpdate forward_entries(const PolicyEngine& pe, const PolicyEpoch& cur, ServiceCapacityWorkspace& W, const PolicyCall& rq,
                       const ForwardingProgram* before_program, const ForwardingView& view, bool active, bool holds,
                       std::shared_ptr<const TunnelPlan> plan, ProgramFlow* given = nullptr) {
  const auto v = rq.offer;
  const ProgramData& p = *view.data;
  const std::size_t n = p.entries.size();
  const std::size_t k = view.num_sources();
  const bool lossy = pe.config.admission == Admission::Lossy;
  const bool created = !cur.prepared;
  // The entries' live records (none before the first call).
  std::vector<const FlowRecord*> before;
  if (!created) {
    before.resize(n);
    for (std::size_t e = 0; e < n; ++e) before[e] = &LedgerAccess::live(pe.ledger, cur.record_groups[e].front().record_counter);
  }
  ProgramFlow flow;
  Amount call_delivered = 0;
  Counter128 offered = cur.offered;
  if (!lossy) {
    // The largest total in [admitted, admitted + v] on the residual plus the policy's own load.
    Amount admitted = 0;
    for (const FlowRecord* r : before) {
      if (r == nullptr) continue;
      admitted += r->delivered;
    }
    if (!created) W.release_records(before);
    const std::optional<Amount> hi = v.has_value() ? std::optional<Amount>(admitted + *v) : std::nullopt;
    flow = detail::lossless(view, available_capacity_reader(W), admitted, hi);
    call_delivered = flow.admitted - admitted;
  } else if (given != nullptr) {
    flow = std::move(*given);
    call_delivered = flow.delivered;
    flow.source_drops += *v - flow.offered;
    accumulate(offered, c128(*v), "offered");
  } else {
    std::vector<Amount> inject(k, 0);
    if (rq.is_full_offer()) {
      source_shares(view, *v, inject);
    } else {
      source_offer_increments(view, cur.offered, *v, inject);
    }
    flow = detail::lossy(view, available_capacity_reader(W), inject);
    call_delivered = flow.delivered;
    // The whole offer is sent: what no head receives (every remaining captured weight 0) is lost at
    // the sources, outside every row.
    flow.source_drops += *v - flow.offered;
    accumulate(offered, c128(*v), "offered");
  }
  PendingPolicyUpdate s;
  s.drafts.resize(n);
  s.next_epoch.record_groups.assign(n, {});
  // The by_tunnel amounts per entry, kept only when the program has tunnels (only they have rows).
  const bool rows = !p.tunnels.empty();
  std::vector<Counter128> assigned(rows ? n : 0);
  std::vector<Amount> delivered(rows ? n : 0, 0);
  std::vector<Counter128> dropped(rows ? n : 0);
  std::vector<std::pair<EdgeId, Amount>> call_edge_drops;
  for (std::size_t e = 0; e < n; ++e) {
    const auto& entry = p.entries[e];
    auto& f = flow.entries[e];
    LiveFlow live;
    live.source = entry.head;
    live.route = view.route(e);
    live.min_route_cost = view.min_route_cost(e);
    live.entry_id = entry.entry_id;
    if (!created) {
      live.record_counter = cur.record_groups[e].front().record_counter;
      s.drafts[e].record_counter = live.record_counter;
    } else {
      s.fresh.emplace_back(e, 0);
    }
    s.next_epoch.record_groups[e].push_back(live);
    FlowRecord& r = s.drafts[e].record;
    r.sources = {entry.head};
    r.destinations = pe.allowed_destinations;
    r.admission = pe.config.admission;
    r.route = live.route;
    r.min_route_cost = live.min_route_cost;
    r.entry_id = entry.entry_id;
    if (holds && entry.tunnel_position >= 0) {
      const auto& t = p.tunnels[static_cast<std::size_t>(entry.tunnel_position)];
      if (t.reserved_bandwidth != 0) {
        const auto path = path_of(t);
        r.edge_reservations.reserve(path.size());
        for (EdgeId x : path) r.edge_reservations.emplace_back(x, t.reserved_bandwidth);
        std::sort(r.edge_reservations.begin(), r.edge_reservations.end());
      }
    }
    const FlowRecord* b = created ? nullptr : before[e];
    if (!lossy) {
      r.edge_flows = std::move(f.edge_flows);
      r.delivered = f.delivered;
      if (rows) {
        assigned[e] = c128(f.delivered);
        delivered[e] = f.delivered;
      }
      continue;
    }
    r.edge_drops = plus(b != nullptr ? b->edge_drops : std::vector<std::pair<EdgeId, Counter128>> {}, f.edge_drops);
    r.offered = b != nullptr ? b->offered.value_or(Counter128 {}) : Counter128 {};
    accumulate(*r.offered, c128(f.assigned), "offered");
    r.delivered = (b != nullptr ? b->delivered : 0) + f.delivered;
    r.edge_flows = b != nullptr ? plus(b->edge_flows, f.edge_flows) : std::move(f.edge_flows);
    call_edge_drops.insert(call_edge_drops.end(), f.edge_drops.begin(), f.edge_drops.end());
    Counter128 drops;
    for (const auto& [x, d] : r.edge_drops) accumulate(drops, d, "drops");
    if (rows) {
      assigned[e] = *r.offered;
      delivered[e] = r.delivered;
      dropped[e] = drops;
    }
  }
  // The shares the heads lost themselves (Lossy): per view head, cumulative.
  std::vector<std::pair<NodeId, Counter128>> lost;
  if (lossy) {
    s.next_epoch.source_drops_by_source = cur.source_drops_by_source.empty() ? std::vector<Counter128>(k) : cur.source_drops_by_source;
    for (std::size_t h = 0; h < k; ++h) {
      accumulate(s.next_epoch.source_drops_by_source[h], c128(flow.source_drops_by_source[h]), "offered");
      if (!s.next_epoch.source_drops_by_source[h].is_zero()) lost.emplace_back(view.source(h), s.next_epoch.source_drops_by_source[h]);
    }
  }
  const Counter128 requested = requested_after(cur, v);
  PlacementResult& out = s.result;
  out.call_delivered = call_delivered;
  if (v.has_value()) out.call_undelivered = *v - call_delivered;
  out.delivery_delta = call_delivered;
  if (lossy) out.call_edge_drops = merged(std::move(call_edge_drops));
  out.call_source_drops = flow.source_drops;
  // No head reaches a destination (every entry failed): NO_ROUTE, the records kept per entry.
  out.termination = k == 0 ? Termination::NoRoute : forwarding_termination(active, p.has_path(), lossy, v, call_delivered);
  set_flags(out, before_program, view.program.get(), is<DistributedTE>(pe));
  out.program = view.program;
  if (plan) out.plan = *plan;
  out.by_tunnel = tunnel_rows(view, assigned, delivered, dropped);
  fill_state(pe, draft_records(s.drafts), requested, cur.infinite || !v.has_value(), offered, out, lost);
  s.next_epoch.program = view.program;
  s.next_epoch.plan = std::move(plan);
  s.next_epoch.by_tunnel = out.by_tunnel;
  s.next_epoch.requested = requested;
  s.next_epoch.offered = offered;
  s.next_epoch.prepared = true;
  s.next_epoch.infinite = cur.infinite || !v.has_value();
  return s;
}

// Pinned: forwarding its program as the policy replays it (place, evaluate, replace); no holds.
PendingPolicyUpdate place_pinned(const PolicyEngine& pe, const PolicyEpoch& cur, ServiceCapacityWorkspace& W, const PolicyCall& rq,
                    const ForwardingProgram* before_program) {
  return forward_entries(pe, cur, W, rq, before_program, pe.view, pe.can_forward, false, nullptr);
}

// T1's setup view U: the reservation baseline minus every hold, this policy's own holds released
// (never above static capacity: R0 <= capacity and the own holds are part of S).
std::vector<Amount> reservation_view(const PolicyEngine& pe, const LedgerTotals& totals) {
  const auto r0 = totals.reservation_baseline;
  const auto held = totals.held;
  std::vector<Amount> U(r0.begin(), r0.end());
  if (!held.empty()) {
    for (std::size_t e = 0; e < U.size(); ++e) U[e] -= held[e];
  }
  for (const FlowRecord* r : current_records(pe)) {
    for (const auto& [e, x] : r->edge_reservations) U[static_cast<std::size_t>(e)] += x;
  }
  return U;
}

// The design volume of a tunnel model's call: the explicit design, or the offer (nullopt:
// unlimited; T1's is finite, checked by check_call).
std::optional<Amount> design_of(const PolicyCall& rq) noexcept { return rq.design.has_value() ? *rq.design : rq.offer; }

// The plan of a T1 setup for the design D: the tunnel count, D, the heads' shares and the holds
// of the up tunnels counted once each.
TunnelPlan distributed_te_plan(const PolicyEngine& pe, Amount D, const ForwardingProgram* program) {
  TunnelPlan plan;
  plan.tunnels_per_head = std::get<DistributedTE>(pe.config.routing).tunnel_count.per_head;
  plan.design_demand = D;
  plan.design_by_head = design_by_head(pe.reachable_sources, D);
  if (program != nullptr) {
    for (const auto& t : program->tunnels()) plan.total_reserved_bandwidth += t.reserved_bandwidth;
  }
  return plan;
}

// The tunnel models' forward(): the installed program forwarded with its plan (and the call's
// black holes), records per program entry, T1's holds kept.
PendingPolicyUpdate forward_tunnels(const PolicyEngine& pe, const PolicyEpoch& cur, ServiceCapacityWorkspace& W, const PolicyCall& rq,
                       const ForwardingProgram* before_program) {
  const ForwardingView view = make_forwarding_view(*pe.graph, pe.epoch.program, nullptr, pe.allowed_destinations, rq.black_holes, true);
  return forward_entries(pe, cur, W, rq, before_program, view, view_active(view), is<DistributedTE>(pe), pe.epoch.plan);
}

// T1: replace and evaluate set up the tunnels for the design on the unreserved view and forward
// the new program, one record per entry.
PendingPolicyUpdate apply_distributed_te_call(const PolicyEngine& pe, const PolicyEpoch& cur, ServiceCapacityWorkspace& W, const PolicyCall& rq,
                const ForwardingProgram* before_program) {
  const auto& g = *pe.graph;
  const Amount D = *design_of(rq);
  const auto& ca = std::get<DistributedTE>(pe.config.routing);
  auto tunnels = setup_distributed_te_tunnels(g, pe.eligibility, pe.reachable_sources, pe.allowed_destinations, ca, pe.c0, D, reservation_view(pe, totals_of(pe, rq)));
  const auto program = distributed_te_program(g, pe.sources, pe.destinations, pe.reachable_sources, std::move(tunnels),
                                      pe.excluded_nodes, pe.excluded_edges);
  const ForwardingView view = make_forwarding_view(g, program, nullptr, pe.allowed_destinations, nullptr, true);
  return forward_entries(pe, cur, W, rq, before_program, view, view_active(view), true,
                         std::make_shared<const TunnelPlan>(distributed_te_plan(pe, D, program.get())));
}

// The plan of a T2 layout (zeros without a tunnel) for the design D (nullopt: unlimited), with its
// K tunnels and its program's forwarding bound F.
TunnelPlan t2_plan(const PolicyEngine& pe, std::optional<Amount> D, const T2Layout& lay, std::size_t K, Amount F) {
  TunnelPlan plan;
  plan.tunnels_per_head = lay.count;
  plan.design_demand = D;
  if (D.has_value()) plan.design_by_head = design_by_head(pe.reachable_sources, *D);
  plan.heads_with_capacity = lay.live;
  plan.insufficient_resolution = !lay.live.empty() && lay.common == 0;
  plan.max_fractional_bandwidth_per_tunnel = std::pair<Amount, Amount> {lay.fluid_num, lay.fluid_den};
  plan.max_bandwidth_per_tunnel = lay.common;
  plan.planned_bandwidth_per_tunnel = lay.planned;
  plan.total_planned_bandwidth = static_cast<Amount>(K) * lay.planned;
  plan.total_tunnel_cost = lay.route_cost;
  if (std::get<CentralizedTE>(pe.config.routing).prefer_headroom) plan.headroom_penalty = lay.headroom_penalty;
  plan.max_lossless_offer = F;
  return plan;
}

// T2: replace and evaluate lay out the tunnels for the design on the released residual and
// forward the new program, one record per entry; the plan's forwarding bound is the program's
// lossless maximum on that view.
PendingPolicyUpdate place_t2(const PolicyEngine& pe, const PolicyEpoch& cur, ServiceCapacityWorkspace& W, const PolicyCall& rq,
                const ForwardingProgram* before_program) {
  const auto& g = *pe.graph;
  const std::optional<Amount> D = design_of(rq);
  std::vector<Amount> V;
  W.dense(V);
  const auto& gl = std::get<CentralizedTE>(pe.config.routing);
  auto lay = t2_layout(g, pe.eligibility, pe.reachable_sources, pe.allowed_destinations, *gl.tunnel_count, gl.prefer_headroom, D, V);
  const auto program = distributed_te_program(g, pe.sources, pe.destinations, pe.reachable_sources, std::move(lay.tunnels),
                                      pe.excluded_nodes, pe.excluded_edges);
  const ForwardingView view = make_forwarding_view(g, program, nullptr, pe.allowed_destinations, nullptr, true);
  const Amount F = detail::lossless(view, available_capacity_reader(W), 0, std::nullopt).admitted;
  return forward_entries(pe, cur, W, rq, before_program, view, view_active(view), false,
                         std::make_shared<const TunnelPlan>(t2_plan(pe, D, lay, program->tunnels().size(), F)));
}

// The view forward() forwards a computing model's installed program on (its roster, the call's
// black holes).
ForwardingView installed_view(const PolicyEngine& pe, const PolicyCall& rq) {
  return make_forwarding_view(*pe.graph, pe.epoch.program, nullptr, pe.allowed_destinations, rq.black_holes);
}

// A SIMULTANEOUS ShortestPaths replace (EQUAL or WEIGHTED): its program prepared on the list's view,
// its per-source flows (routes and costs, new ids) and the program installed before, for the flags.
struct PreparedReplace {
  std::shared_ptr<const ForwardingProgram> program {};
  std::vector<LiveFlow> flows {};
  const ForwardingProgram* before {nullptr};
};

// A forwarding of a computing model's program as `view` (its roster) sees it, finished into records
// per roster source: forward() of the installed program, or (`replaced`) a SIMULTANEOUS ShortestPaths
// replace of its newly prepared one (its routes and costs, its flags, its termination).
PendingPolicyUpdate finish_by_source(const PolicyEngine& pe, const PolicyCall& rq, const ForwardingView& view, ProgramFlow& flow,
                                     const PreparedReplace* replaced) {
  const auto& g = *pe.graph;
  const auto v = rq.offer;
  const auto& program = replaced != nullptr ? replaced->program : pe.epoch.program;
  const ProgramData& p = *view.data;
  const std::size_t k = pe.reachable_sources.size();
  const std::size_t n = p.entries.size();
  const bool lossy = pe.config.admission == Admission::Lossy;
  PendingPolicyUpdate s;
  s.next_epoch = fresh_epoch(pe, replaced == nullptr);
  s.next_epoch.program = program;
  std::vector<LiveFlow> flows(k);
  // Per roster source: its entries' carried, drops and deliveries summed (the view's heads are the
  // roster: nothing is replayed).
  std::vector<std::vector<std::pair<EdgeId, Amount>>> carried(k);
  std::vector<std::vector<std::pair<EdgeId, Amount>>> drops(k);
  std::vector<Amount> got(k, 0);
  for (std::size_t e = 0; e < n; ++e) {
    const auto slot = static_cast<std::size_t>(p.entries[e].roster_index);
    const auto& f = flow.entries[e];
    carried[slot].insert(carried[slot].end(), f.edge_flows.begin(), f.edge_flows.end());
    drops[slot].insert(drops[slot].end(), f.edge_drops.begin(), f.edge_drops.end());
    got[slot] += f.delivered;
  }
  for (std::size_t i = 0; i < k; ++i) {
    flows[i] = replaced != nullptr ? replaced->flows[i] : pe.epoch.record_groups[i].front();
    flows[i].record_counter = -1;
    carried[i] = merged(std::move(carried[i]));
    drops[i] = merged(std::move(drops[i]));
    if (is<CentralizedTE>(pe)) {
      // CentralizedTE's record route is the edges its share uses.
      std::vector<EdgeId> route;
      route.reserve(carried[i].size());
      for (const auto& [e, x] : carried[i]) route.push_back(e);
      flows[i].route = std::make_shared<const PredecessorDAG>(PredecessorDAGBuilder::of_edges(g, route));
      flows[i].min_route_cost = route.empty() ? kUnreachable : min_route_cost(g, *flows[i].route, pe.reachable_sources[i], pe.allowed_destinations);
    }
  }
  source_drafts(pe, flows, s);
  std::vector<std::pair<EdgeId, Amount>> call_edge_drops;
  for (std::size_t i = 0; i < k; ++i) {
    FlowRecord& r = s.drafts[i].record;
    r.edge_flows = carried[i];
    r.delivered = got[i];
    if (!lossy) continue;
    r.edge_drops = plus(std::vector<std::pair<EdgeId, Counter128>> {}, drops[i]);
    r.offered = c128(flow.offered_by_source[i]);
    call_edge_drops.insert(call_edge_drops.end(), drops[i].begin(), drops[i].end());
  }
  std::vector<Counter128> assigned(n);
  std::vector<Amount> delivered(n, 0);
  std::vector<Counter128> dropped(n);
  for (std::size_t e = 0; e < n; ++e) {
    const auto& f = flow.entries[e];
    assigned[e] = c128(lossy ? f.assigned : f.delivered);
    delivered[e] = f.delivered;
    for (const auto& [x, d] : f.edge_drops) accumulate(dropped[e], c128(d), "drops");
  }
  const Counter128 requested = v.has_value() ? c128(*v) : Counter128 {};
  const Counter128 offered = lossy ? c128(flow.offered) : Counter128 {};
  PlacementResult& out = s.result;
  out.call_delivered = lossy ? flow.delivered : flow.admitted;
  if (v.has_value()) out.call_undelivered = *v - out.call_delivered;
  out.call_edge_drops = merged(std::move(call_edge_drops));
  out.call_source_drops = flow.source_drops;
  bool active = false;
  for (std::size_t h = 0; h < view.num_sources() && !active; ++h) {
    for (const std::int32_t e : view.entries_for_source(h)) active = active || p.entries[static_cast<std::size_t>(e)].weight > 0;
  }
  if (p.pools_source_inflows) active = !view.prepared_route(0)->forwarding_edges.empty();
  out.termination = forwarding_termination(active, p.has_path(), lossy, v, out.call_delivered);
  if (replaced != nullptr) {
    // A replace's own rule (place_pooled's): no member, no active weight; else the Lossy offer
    // is sent whole.
    out.termination = view.prepared_route(0)->forwarding_edges.empty() ? Termination::NoActiveWeight : Termination::Complete;
    set_flags(out, replaced->before, program.get(), false);
  }
  out.program = program;
  out.by_tunnel = tunnel_rows(view, assigned, delivered, dropped);
  const auto after = draft_records(s.drafts);
  fill_state(pe, after, requested, !v.has_value(), offered, out);
  if (uses_optimized_flow(pe)) out.total_cost = total_cost_of(g, after);
  s.next_epoch.by_tunnel = out.by_tunnel;
  s.next_epoch.requested = requested;
  s.next_epoch.offered = offered;
  s.next_epoch.prepared = true;
  s.next_epoch.infinite = !v.has_value();
  return s;
}

// forward(): a computing model's installed program forwarded with no route computation, on the
// released view, records per roster source (ShortestPaths and CentralizedTE).
PendingPolicyUpdate forward_installed(const PolicyEngine& pe, ServiceCapacityWorkspace& W, const PolicyCall& rq) {
  const auto v = rq.offer;
  const ForwardingView view = installed_view(pe, rq);
  ProgramFlow flow;
  if (pe.config.admission != Admission::Lossy) {
    flow = detail::lossless(view, available_capacity_reader(W), 0, v);
  } else {
    std::vector<Amount> inject(pe.reachable_sources.size(), 0);
    source_shares(view, *v, inject);
    flow = detail::lossy(view, available_capacity_reader(W), inject);
  }
  return finish_by_source(pe, rq, view, flow, nullptr);
}

// ---------------------------------------------------------------- construction

} // namespace

const char* model_name(const Routing& r) noexcept {
  if (std::holds_alternative<ShortestPaths>(r)) return "ShortestPaths";
  if (is_t2(r)) return "CentralizedTE with tunnel_count";
  if (std::holds_alternative<CentralizedTE>(r)) return "CentralizedTE";
  if (std::holds_alternative<DistributedTE>(r)) return "DistributedTE";
  return "Pinned";
}

const char* model_class(const Routing& r) noexcept {
  if (std::holds_alternative<ShortestPaths>(r)) return "ShortestPaths";
  if (std::holds_alternative<CentralizedTE>(r)) return "CentralizedTE";
  if (std::holds_alternative<DistributedTE>(r)) return "DistributedTE";
  return "Pinned";
}

namespace {

// T2's several-sources rule: OPTIMIZED with one input source is the equal calculation; with several
// it would be another (strongly NP-hard) problem.
void check_t2_sources(const Routing& r, std::size_t num_sources) {
  const auto* gl = std::get_if<CentralizedTE>(&r);
  if (gl != nullptr && gl->tunnel_count.has_value() && num_sources > 1 && gl->source_allocation == SourceAllocation::Optimized) {
    invalid("CentralizedTE with tunnel_count and several sources requires source_allocation=EQUAL");
  }
}

// The pooled single-path rule: a bottleneck choice would differ per source at a shared router.
void check_pooled_single_path(const Routing& r, std::size_t num_sources) {
  const auto* sp = std::get_if<ShortestPaths>(&r);
  if (sp != nullptr && num_sources > 1 && !sp->multipath && detail::RouteChoice {sp->tie_break}.prefers_score()) {
    invalid("ShortestPaths with multiple input sources requires tie_break=FEWER_HOPS when multipath=False");
  }
}

} // namespace

void check_config(const FlowPolicyConfig& config, bool cache, std::int32_t num_edges, std::optional<std::size_t> num_sources) {
  const auto& r = config.routing;
  if (r.valueless_by_exception()) invalid("routing: the variant holds no routing mode");
  // The enumeration values.
  check_member(config.admission, "admission");
  const auto* sp = std::get_if<ShortestPaths>(&r);
  const auto* gl = std::get_if<CentralizedTE>(&r);
  const auto* t1 = std::get_if<DistributedTE>(&r);
  if (sp != nullptr) {
    check_member(sp->allocation, "allocation");
    if (sp->source_allocation.has_value()) check_member(*sp->source_allocation, "source_allocation");
    if (sp->weight_basis.has_value()) check_member(*sp->weight_basis, "weight_basis");
    check_member(sp->tie_break, "tie_break");
  }
  if (gl != nullptr) check_member(gl->source_allocation, "source_allocation");
  if (t1 != nullptr) check_member(t1->tie_break, "tie_break");
  const bool t2 = is_t2(r);
  const bool lossy = config.admission == Admission::Lossy;
  // The model matrix in its order: route and source allocation, admission, the route-choice
  // controls, prefer_headroom, the weight combinations, then a cache on a model that does not read it.
  if (sp != nullptr) {
    const bool optimized = sp->allocation == RouteAllocation::Optimized;
    if (!optimized && sp->source_allocation == SourceAllocation::Optimized) {
      invalid("ShortestPaths with EQUAL or WEIGHTED requires source_allocation=EQUAL");
    }
    if (optimized && lossy) invalid("ShortestPaths OPTIMIZED requires admission=LOSSLESS");
    if (detail::RouteChoice {sp->tie_break, sp->multipath, sp->parallel_edges}.tie_is_ineffective()) {
      invalid("ShortestPaths tie_break requires multipath=False or parallel_edges=False");
    }
    if (num_sources.has_value()) check_pooled_single_path(r, *num_sources);
    const bool weighted = sp->allocation == RouteAllocation::Weighted;
    if (sp->weight_basis.has_value() && !weighted) invalid("weight_basis requires allocation=WEIGHTED");
    if (sp->weight_basis == WeightBasis::Explicit && !sp->edge_weights) invalid("weight_basis=EXPLICIT requires edge_weights");
    if (sp->edge_weights && sp->weight_basis != WeightBasis::Explicit) invalid("edge_weights require weight_basis=EXPLICIT");
  }
  if (gl != nullptr) {
    if (num_sources.has_value()) check_t2_sources(r, *num_sources);
    if (!t2 && lossy) invalid("CentralizedTE without tunnel_count requires admission=LOSSLESS");
    if (gl->prefer_headroom && !t2) invalid("prefer_headroom requires tunnel_count");
  }
  if (cache && sp == nullptr) invalid(std::string(model_class(r)) + " does not use ShortestPathCache");
  // Then the tunnel count, an "up to" count for DistributedTE, the ceilings, the search minimum and
  // edge_weights' length and entries (they need the graph).
  if (t2 && (gl->tunnel_count->per_head < 1 || gl->tunnel_count->per_head > kMaxTunnelCount)) invalid("tunnel count must be in [1, 10000]");
  if (t1 != nullptr) {
    if (t1->tunnel_count.per_head < 1 || t1->tunnel_count.per_head > kMaxTunnelCount) invalid("tunnel count must be in [1, 10000]");
    if (t1->tunnel_count.up_to) invalid("DistributedTE requires an exact tunnel count");
    check_ceiling(t1->cost_ceiling);
    if (t1->min_unreserved_bandwidth < 0 || t1->min_unreserved_bandwidth >= kDomainLimit) {
      invalid("min_unreserved_bandwidth must be an integer in [0, 2^62)");
    }
  }
  if (sp != nullptr && sp->allocation == RouteAllocation::Weighted && sp->edge_weights) {
    const auto& w = *sp->edge_weights;
    if (w.size() != static_cast<std::size_t>(num_edges)) invalid("edge_weights: expected " + str(num_edges) + " entries");
    for (std::size_t e = 0; e < w.size(); ++e) {
      if (w[e] < 0 || w[e] >= kDomainLimit) invalid("edge_weights[" + std::to_string(e) + "] must be an integer in [0, 2^62)");
    }
  }
}

void check_source_count(const FlowPolicyConfig& config, std::size_t num_sources) {
  check_pooled_single_path(config.routing, num_sources);
  check_t2_sources(config.routing, num_sources);
}

FlowPolicyConfig resolved_config(FlowPolicyConfig config) {
  if (auto* sp = std::get_if<ShortestPaths>(&config.routing)) {
    if (!sp->source_allocation.has_value()) {
      sp->source_allocation = sp->allocation == RouteAllocation::Optimized ? SourceAllocation::Optimized : SourceAllocation::Equal;
    }
    if (sp->allocation == RouteAllocation::Weighted && !sp->weight_basis.has_value()) sp->weight_basis = WeightBasis::Capacity;
  }
  return config;
}

std::optional<SourceAllocation> source_allocation_of(const Routing& r) noexcept {
  if (const auto* sp = std::get_if<ShortestPaths>(&r)) {
    return sp->source_allocation.value_or(sp->allocation == RouteAllocation::Optimized ? SourceAllocation::Optimized
                                                                                          : SourceAllocation::Equal);
  }
  if (const auto* gl = std::get_if<CentralizedTE>(&r)) return gl->source_allocation;
  if (std::holds_alternative<DistributedTE>(r)) return SourceAllocation::Equal;
  return std::nullopt;
}

void check_tunnel_limit(const FlowPolicyConfig& cfg, std::size_t num_heads) {
  const auto* t1 = std::get_if<DistributedTE>(&cfg.routing);
  const auto* gl = std::get_if<CentralizedTE>(&cfg.routing);
  if (t1 == nullptr && (gl == nullptr || !gl->tunnel_count.has_value())) return;
  // The count (the largest for "up to") is in [1, 10000] (check_config) and H < 2^31: the
  // product fits 64 bits.
  const std::int64_t count = t1 != nullptr ? t1->tunnel_count.per_head : gl->tunnel_count->per_head;
  const auto heads = static_cast<std::int64_t>(num_heads);
  if (heads * count > kMaxTunnels) {
    invalid("tunnel limit exceeded: " + str(heads) + " heads x " + str(count) + " tunnels > " + str(kMaxTunnels));
  }
}

} // namespace detail

// ================================================================ FlowPolicy

using detail::LedgerAccess;
using detail::PolicyEngine;

namespace {

// The engine of a new policy, bound to its ledger as the last step. A supplied frozen cache's hit
// is recorded in `tally`: the public constructor publishes it once construction has
// succeeded, an evaluator once its whole call has.
std::unique_ptr<PolicyEngine> make_engine(FlowLedger& ledger, std::span<const NodeId> sources,
                                          std::span<const NodeId> destinations, FlowPolicyConfig config,
                                          const Eligibility* eligibility, ShortestPathCache* cache, detail::FrozenCacheHitTally& tally) {
  auto pe = std::make_unique<PolicyEngine>(ledger, std::move(config));
  pe->graph = ledger.graph();
  const auto& g = *pe->graph;
  // The endpoints, then the ownership of the mask, the program, the support, the cache and the
  // black holes, then the configuration.
  pe->sources = detail::node_set(sources, g.num_nodes(), "sources");
  pe->destinations = detail::node_set(destinations, g.num_nodes(), "destinations");
  detail::check_disjoint(pe->sources, pe->destinations, "sources", "destinations");
  if (eligibility != nullptr && eligibility->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  const auto* pinned = std::get_if<Pinned>(&pe->config.routing);
  const auto* global = std::get_if<CentralizedTE>(&pe->config.routing);
  if (pinned != nullptr) {
    if (!pinned->program) throw std::invalid_argument("Pinned requires a program");
    if (pinned->program->graph_uid() != g.uid()) throw std::invalid_argument("program belongs to another graph");
  }
  if (global != nullptr && global->support && global->support->graph_uid() != g.uid()) {
    throw std::invalid_argument("support belongs to another graph");
  }
  if (cache != nullptr && cache->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign ShortestPathCache: bound to another graph");
  }
  if (pinned != nullptr && pinned->black_holes && pinned->black_holes->graph_uid() != g.uid()) {
    throw std::invalid_argument("black_holes belongs to another graph");
  }
  pe->eligibility = eligibility;
  detail::check_config(pe->config, cache != nullptr, g.num_edges(), pe->sources.size());
  pe->config = detail::resolved_config(std::move(pe->config));
  if (global != nullptr && global->support && global->tunnel_count.has_value()) {
    // T2 on a support: the support's edges intersected with the mask are the allowed links of every
    // step (roster, sizing, objectives, decomposition, the programs' replay basis).
    std::vector<NodeId> xn;
    std::vector<EdgeId> xe;
    const auto kept = detail::PredecessorDAGBuilder::ascending(*global->support);
    for (NodeId v = 0; v < g.num_nodes(); ++v) {
      if (eligibility != nullptr && !eligibility->allows_node(v)) xn.push_back(v);
    }
    for (EdgeId e = 0; e < g.num_edges(); ++e) {
      const bool in_support = std::binary_search(kept.begin(), kept.end(), e);
      if (!in_support || (eligibility != nullptr && !eligibility->allows_edge(e))) xe.push_back(e);
    }
    pe->support_mask = std::make_unique<const Eligibility>(g.eligibility(xn, xe));
    eligibility = pe->support_mask.get();
    pe->eligibility = eligibility;
  }
  if (pinned != nullptr) {
    const auto& prog = *pinned->program;
    if (!std::equal(prog.sources().begin(), prog.sources().end(), pe->sources.begin(), pe->sources.end()) ||
        !std::equal(prog.destinations().begin(), prog.destinations().end(), pe->destinations.begin(),
                    pe->destinations.end())) {
      throw std::invalid_argument("program endpoints do not match policy endpoints");
    }
  }
  if (const auto* sp = std::get_if<ShortestPaths>(&pe->config.routing)) {
    pe->route_choice = detail::RouteChoice {sp->tie_break, sp->multipath, sp->parallel_edges};
  } else if (const auto* ca = std::get_if<DistributedTE>(&pe->config.routing)) {
    pe->route_choice = detail::RouteChoice {ca->tie_break, false, ca->parallel_edges};
  }
  if (const auto s = detail::source_allocation_of(pe->config.routing)) pe->source_allocation = *s;
  auto node_allowed = [&](NodeId v) { return eligibility == nullptr || eligibility->allows_node(v); };
  for (NodeId d : pe->destinations) {
    if (node_allowed(d)) pe->allowed_destinations.push_back(d);
  }
  if (pinned != nullptr) {
    // The program as this policy forwards it: replayed under the mask, the black holes laid over.
    pe->view = detail::make_forwarding_view(g, pinned->program, eligibility, pe->allowed_destinations, pinned->black_holes.get(), true);
    for (std::size_t h = 0; h < pe->view.num_sources(); ++h) pe->reachable_sources.push_back(pe->view.source(h));
    pe->can_forward = detail::view_active(pe->view);
    pe->epoch.record_groups.resize(detail::ProgramAccess::data(*pinned->program).entries.size());
    pe->policy_id = LedgerAccess::bind(ledger);
    return pe;
  }
  if (detail::is<ShortestPaths>(*pe)) {
    pe->support = cache != nullptr
                      ? detail::DagCacheAccess::lookup(*cache, g, pe->destinations, Direction::Reverse, eligibility, tally)
                      : std::make_shared<const PredecessorDAG>(detail::static_support(g, pe->destinations, Direction::Reverse, eligibility));
    const auto edges = detail::PredecessorDAGBuilder::ascending(*pe->support);
    const auto row = g.outgoing_offsets();
    for (NodeId s : pe->sources) {
      if (!node_allowed(s)) continue;
      const auto it = std::lower_bound(edges.begin(), edges.end(), row[static_cast<std::size_t>(s)]);
      if (it != edges.end() && *it < row[static_cast<std::size_t>(s) + 1]) pe->reachable_sources.push_back(s);
    }
  } else if (global != nullptr && global->support && !global->tunnel_count.has_value()) {
    // The support trimmed to the walks from the sources to the allowed destinations under the mask;
    // the roster is the sources that reach one inside it.
    std::vector<EdgeId> trimmed;
    for (NodeId s : pe->sources) {
      if (!node_allowed(s)) continue;
      auto part = detail::prune(g, *global->support, s, pe->allowed_destinations, eligibility);
      if (part.empty()) continue;
      pe->reachable_sources.push_back(s);
      trimmed.insert(trimmed.end(), part.begin(), part.end());
    }
    std::sort(trimmed.begin(), trimmed.end());
    trimmed.erase(std::unique(trimmed.begin(), trimmed.end()), trimmed.end());
    pe->centralized_te_support = std::make_shared<const std::vector<EdgeId>>(std::move(trimmed));
  } else {
    // The static costs fix the roster (unrestricted CentralizedTE, T1 and T2); T1's CostStretch baseline c0
    // is their minimum over the roster heads (zero-capacity edges included, before any filter).
    const auto dist = detail::static_costs(g, eligibility, pe->allowed_destinations, pe->sources);
    pe->c0 = kUnreachable;
    for (std::size_t i = 0; i < pe->sources.size(); ++i) {
      const NodeId s = pe->sources[i];
      if (!node_allowed(s) || dist[i] == kUnreachable) continue;
      pe->reachable_sources.push_back(s);
      pe->c0 = std::min(pe->c0, dist[i]);
    }
    detail::check_tunnel_limit(pe->config, pe->reachable_sources.size());
    if (eligibility != nullptr && eligibility->uid() != 0) {
      // Its programs' replay basis: the topology under this mask.
      for (NodeId v = 0; v < g.num_nodes(); ++v) {
        if (!eligibility->allows_node(v)) pe->excluded_nodes.push_back(v);
      }
      for (EdgeId e = 0; e < g.num_edges(); ++e) {
        if (!eligibility->allows_edge(e)) pe->excluded_edges.push_back(e);
      }
    }
  }
  // The tunnel models' cohorts are per program entry, sized by their first call.
  if (!detail::tunnel_model(*pe)) pe->epoch.record_groups.resize(pe->reachable_sources.size());
  pe->policy_id = LedgerAccess::bind(ledger);  // last: nothing is bound when construction fails
  return pe;
}

} // namespace

// A public call of its own: a supplied frozen cache counts the hit once the policy exists.
FlowPolicy::FlowPolicy(FlowLedger& ledger, std::span<const NodeId> sources, std::span<const NodeId> destinations,
                       FlowPolicyConfig config, const Eligibility* eligibility, ShortestPathCache* cache)
    : engine_(detail::frozen_call(cache, [&](detail::FrozenCacheHitTally& tally) {
        return make_engine(ledger, sources, destinations, std::move(config), eligibility, cache, tally);
      })) {}

FlowPolicy::FlowPolicy(FlowLedger& ledger, std::span<const NodeId> sources, std::span<const NodeId> destinations,
                       FlowPolicyConfig config, const Eligibility* eligibility, ShortestPathCache* cache, detail::FrozenCacheHitTally& tally)
    : engine_(make_engine(ledger, sources, destinations, std::move(config), eligibility, cache, tally)) {}

// The destruction law: the records stay in the ledger as orphans.
FlowPolicy::~FlowPolicy() { LedgerAccess::release(engine_->ledger, engine_->policy_id); }

namespace detail {
namespace {

// The program a computing model installs for an empty roster: an empty one. Pinned never
// comes here: its forwarding with no head is place_pinned's.
std::shared_ptr<const ForwardingProgram> empty_roster_program(const PolicyEngine& pe) {
  return empty_program(pe.graph->uid(), pe.sources, pe.destinations, uses_optimized_flow(pe) && pe.source_allocation == SourceAllocation::Optimized);
}

// The current state of a zero call: its records, its program and rows, every flag false; before
// the epoch's first preparation all zeros (no records, no program, no curve). A zero call solves
// nothing, so it claims no cut.
PlacementResult zero_state(const PolicyEngine& pe, const LedgerTotals* totals, const OutputOptions& outputs) {
  PlacementResult out;
  out.call_undelivered = 0;
  const auto& ep = pe.epoch;
  const auto records = current_records(pe);
  std::vector<std::pair<NodeId, Counter128>> lost;
  for (std::size_t h = 0; h < ep.source_drops_by_source.size(); ++h) {
    if (!ep.source_drops_by_source[h].is_zero()) lost.emplace_back(pe.view.source(h), ep.source_drops_by_source[h]);
  }
  fill_state(pe, records, ep.requested, ep.infinite, ep.offered, out, lost);
  if (uses_optimized_flow(pe)) out.total_cost = total_cost_of(*pe.graph, records);
  if (is_unrestricted_centralized_te(pe) && pe.source_allocation == SourceAllocation::Optimized && ep.prepared && ep.marginal_cost_curve.has_value()) out.marginal_cost_curve = *ep.marginal_cost_curve;
  out.program = ep.program;
  out.by_tunnel = ep.by_tunnel;
  prepare_outputs(pe, totals, outputs, records, out);
  return out;
}

// A computing model's call when no source reaches an allowed destination (Pinned with no head
// forwards through place_pinned, keeping its per-entry records): nothing is forwarded. A Lossy
// offer is lost before it enters (the demand's source loss; no source's row has it). CentralizedTE's flow
// is 0 and maximal: its curve is empty and its cut is empty (no roster source, so the closure is
// empty).
PendingPolicyUpdate empty_roster(const PolicyEngine& pe, const PolicyEpoch& cur, const PolicyCall& rq, const ForwardingProgram* before) {
  const bool lossy = pe.config.admission == Admission::Lossy;
  const bool optimized_global = is_unrestricted_centralized_te(pe) && pe.source_allocation == SourceAllocation::Optimized;
  const auto v = rq.offer;
  PendingPolicyUpdate s;
  s.next_epoch.record_groups.resize(cur.record_groups.size());
  s.next_epoch.requested = requested_after(cur, v);
  s.next_epoch.offered = cur.offered;
  if (lossy && v.has_value()) (void)s.next_epoch.offered.add(*v);  // offered <= requested
  s.next_epoch.prepared = true;
  s.next_epoch.infinite = cur.infinite || !v.has_value();
  s.next_epoch.program = cur.program ? cur.program : empty_roster_program(pe);
  // The tunnel models: forward keeps the installed plan; a setup plans no head.
  if (rq.operation == PolicyOperation::Forward) {
    s.next_epoch.plan = cur.plan;
  } else if (is<DistributedTE>(pe)) {
    s.next_epoch.plan = std::make_shared<const TunnelPlan>(distributed_te_plan(pe, *design_of(rq), nullptr));
  } else if (t2(pe)) {
    T2Layout none;
    none.count = std::get<CentralizedTE>(pe.config.routing).tunnel_count->up_to ? 1 : std::get<CentralizedTE>(pe.config.routing).tunnel_count->per_head;
    s.next_epoch.plan = std::make_shared<const TunnelPlan>(t2_plan(pe, design_of(rq), none, 0, 0));
  }
  PlacementResult& out = s.result;
  out.call_undelivered = v;
  out.delivery_delta = 0;
  out.termination = Termination::NoRoute;
  if (lossy && v.has_value()) out.call_source_drops = *v;
  set_flags(out, before, s.next_epoch.program.get(), is<DistributedTE>(pe));
  out.program = s.next_epoch.program;
  if (s.next_epoch.plan) out.plan = *s.next_epoch.plan;
  fill_state(pe, {}, s.next_epoch.requested, s.next_epoch.infinite, s.next_epoch.offered, out);
  if (uses_optimized_flow(pe)) out.total_cost = Counter128 {};
  if (optimized_global) {
    out.marginal_cost_curve.emplace();
    if (rq.outputs != nullptr && rq.outputs->min_cut) out.min_cut.emplace();
    s.next_epoch.marginal_cost_curve.emplace();
  }
  if (kTraced && rq.trace) rq.trace->set("roster", Trace::list());
  return s;
}

} // namespace

void check_place_operation(const PolicyEngine& pe) {
  if (tunnel_model(pe)) invalid(std::string(model_name(pe.config.routing)) + " does not support incremental placement");
}

void check_forward_operation(const PolicyEngine& pe) {
  if (is<Pinned>(pe)) invalid("Pinned forwards with replace");
  if (!pe.epoch.program) invalid("forward requires an installed program");
}

namespace {

// The call's validation (before any change, also before the zero shortcut), in this order: the
// operation for the model (T1's and T2's place) and an installed program (forward), the offer, a
// Lossy full offer's finiteness, the design demand (the model, its domain, T1's finiteness), a
// finite place after +inf, the requested total's limit, the black holes' graph.
void check_call(const PolicyEngine& pe, const PolicyCall& rq) {
  const bool lossy = pe.config.admission == Admission::Lossy;
  const bool t1 = is<DistributedTE>(pe);
  if (rq.operation == PolicyOperation::Place) check_place_operation(pe);
  if (rq.operation == PolicyOperation::Forward) check_forward_operation(pe);
  if (rq.offer.has_value() && (*rq.offer < 0 || *rq.offer >= kDomainLimit)) {
    invalid("offer must be an integer in [0, 2^62) or permitted +inf");
  }
  // A Lossy full offer is finite; ShortestPaths's place(+inf) is its first-hop experiment.
  if (lossy && !rq.offer.has_value() && (rq.is_full_offer() || is<Pinned>(pe))) invalid("LOSSY forwarding requires a finite offer");
  // The design demand: the tunnel models' alone, in its domain, T1's finite (absent: the offer).
  if (rq.design.has_value()) {
    if (!tunnel_model(pe)) invalid(std::string(model_name(pe.config.routing)) + " does not accept design_demand");
    if (rq.design->has_value() && (**rq.design < 0 || **rq.design >= kDomainLimit)) {
      invalid("design_demand must be an integer in [0, 2^62) or permitted +inf");
    }
  }
  if (t1 && rq.operation != PolicyOperation::Forward && (rq.design.has_value() ? !rq.design->has_value() : !rq.offer.has_value())) {
    invalid("DistributedTE requires a finite design demand");
  }
  if (rq.operation == PolicyOperation::Place) {
    if (pe.epoch.infinite && rq.offer.has_value()) invalid("finite placement after infinite placement requires a new epoch");
    // The only check of the cumulative request: nothing changes the requested total between here
    // and the call's publication, so the later additions repeat this sum and cannot fail.
    if (rq.offer.has_value()) {
      Counter128 probe = pe.epoch.requested;
      if (!probe.add(*rq.offer)) invalid("requested limit exceeded");
    }
  }
  if (rq.operation == PolicyOperation::Forward && rq.black_holes != nullptr && rq.black_holes->graph_uid() != pe.graph->uid()) {
    invalid("black_holes belongs to another graph");
  }
}

// The computation of a call: nothing published. A full-offer call works on the residual with
// the policy's own load released, from a fresh epoch.
PendingPolicyUpdate compute(const PolicyEngine& pe, const PolicyCall& rq) {
  const LedgerTotals totals = totals_of(pe, rq);
  ServiceCapacityWorkspace W(totals);
  if (rq.is_full_offer()) W.release_records(current_records(pe));
  const PolicyEpoch fresh = rq.is_full_offer() ? fresh_epoch(pe, rq.operation == PolicyOperation::Forward) : PolicyEpoch {};
  const PolicyEpoch& cur = rq.is_full_offer() ? fresh : pe.epoch;
  // The flags compare the program the policy has installed now.
  const ForwardingProgram* before = pe.epoch.program.get();
  PendingPolicyUpdate s;
  if (pe.reachable_sources.empty() && !is<Pinned>(pe)) {
    // No roster source (a computing model): nothing is computed. A Pinned policy whose heads all
    // fail keeps its per-entry records and forwards nothing (place_pinned).
    s = empty_roster(pe, cur, rq, before);
  } else if (tunnel_model(pe) && rq.operation == PolicyOperation::Forward) {
    s = forward_tunnels(pe, cur, W, rq, before);
  } else if (is<DistributedTE>(pe)) {
    s = apply_distributed_te_call(pe, cur, W, rq, before);
  } else if (t2(pe)) {
    s = place_t2(pe, cur, W, rq, before);
  } else if (rq.operation == PolicyOperation::Forward) {
    s = forward_installed(pe, W, rq);
  } else if (is<Pinned>(pe)) {
    s = place_pinned(pe, cur, W, rq, before);
  } else if (uses_optimized_flow(pe)) {
    s = place_optimized_flow(pe, cur, W, rq, before);
  } else {
    s = place_pooled(pe, cur, W, rq, before);
  }
  // A full-offer call: its call_delivered is the whole new delivery, compared with the old one.
  if (rq.is_full_offer()) {
    Amount old = 0;
    for (const FlowRecord* r : current_records(pe)) old += r->delivered;
    s.result.call_delivered = s.result.delivered;
    s.result.delivery_delta = s.result.delivered - old;
    if (rq.offer.has_value()) s.result.call_undelivered = *rq.offer - s.result.call_delivered;
  }
  // The first termination row: a request of 0 is COMPLETE, also for a full-offer call, which
  // still prepared and publishes what it computed.
  if (rq.offer == Amount {0}) s.result.termination = Termination::Complete;
  if (rq.outputs != nullptr && any_output(*rq.outputs)) prepare_outputs(pe, rq.totals, *rq.outputs, draft_records(s.drafts), s.result);
  return s;
}

// The publication: the ledger's one atomic exchange (a place commits every live record and the new
// ones; a full-offer call retires the live ones and creates the new ones), then the epoch.
void publish(PolicyEngine& pe, PendingPolicyUpdate& s, bool full) {
  const auto ids = full ? LedgerAccess::exchange(pe.ledger, pe.policy_id, s.drafts) : LedgerAccess::commit(pe.ledger, pe.policy_id, s.drafts);
  // Committed: the epoch state follows (nothing below throws).
  for (std::size_t n = 0; n < s.fresh.size(); ++n) s.next_epoch.record_groups[s.fresh[n].first][s.fresh[n].second].record_counter = ids[n].second;
  pe.epoch = std::move(s.next_epoch);
}

// A publishing call (place, replace, forward; evaluate computes without it). `before_publish`
// (optional) runs once the call's result and trace are final and before anything is published (the
// ledger commit, the epoch state): an exception from it aborts the call with nothing changed (the
// Python binding builds its result objects there).
PlacementResult run_call(PolicyEngine& pe, const PolicyCall& rq, const std::function<void()>* before_publish) {
  check_call(pe, rq);
  if (rq.operation == PolicyOperation::Place && rq.offer.has_value() && *rq.offer == 0) {
    PlacementResult out = zero_state(pe, rq.totals, *rq.outputs);
    if (before_publish != nullptr && *before_publish) (*before_publish)();
    return out;
  }
  PendingPolicyUpdate s = compute(pe, rq);
  if (before_publish != nullptr && *before_publish) (*before_publish)();
  PlacementResult out = std::move(s.result);
  publish(pe, s, rq.is_full_offer());
  return out;
}

// ================================================================ lists (apply_actions)

// One action of a list, validated: its policy's engine and the scalar call it stands for.
struct ListAction {
  PolicyEngine* engine {nullptr};
  PolicyCall call {};
};

// The list's validation, before anything is computed: the ledger, the mode, then each
// action in list order, all of its checks before the next one's (ownership, the duplicate rule,
// the SIMULTANEOUS rules, the scalar call's own checks).
std::vector<ListAction> validate_list(FlowLedger& ledger, std::span<const Action> actions, BatchMode mode, const OutputOptions& outputs) {
  (void)LedgerAccess::baseline(ledger);  // a moved-from ledger: std::logic_error
  switch (mode) {
    case BatchMode::Ordered:
    case BatchMode::Simultaneous: break;
    default: not_a_member("mode", static_cast<int>(mode), "BatchMode");
  }
  std::vector<ListAction> out;
  out.reserve(actions.size());
  std::set<std::int64_t> seen;
  for (const Action& a : actions) {
    ListAction x;
    FlowPolicy* policy = nullptr;
    std::visit(
        [&](const auto& act) {
          using T = std::decay_t<decltype(act)>;
          policy = &act.policy.get();
          x.call.offer = act.offer;
          x.call.outputs = &outputs;
          if constexpr (std::is_same_v<T, Place>) {
            x.call.operation = PolicyOperation::Place;
          } else if constexpr (std::is_same_v<T, Replace>) {
            x.call.operation = PolicyOperation::Replace;
            x.call.design = act.design_demand;
          } else {
            x.call.operation = PolicyOperation::Forward;
            x.call.black_holes = act.black_holes;
          }
        },
        a);
    PolicyEngine& pe = PolicyAccess::engine(*policy);
    x.engine = &pe;
    if (&pe.ledger != &ledger) invalid("list policy belongs to another ledger");
    if (!seen.insert(pe.policy_id).second) invalid("duplicate policy in list");
    if (mode == BatchMode::Simultaneous) {
      if (x.call.operation == PolicyOperation::Place) invalid("SIMULTANEOUS does not accept Place");
      if (x.call.operation == PolicyOperation::Replace && !is<Pinned>(pe) &&
          !(is<ShortestPaths>(pe) && std::get<ShortestPaths>(pe.config.routing).allocation != RouteAllocation::Optimized)) {
        invalid("SIMULTANEOUS Replace requires ShortestPaths EQUAL or WEIGHTED, or Pinned");
      }
      if (pe.config.admission != Admission::Lossy) invalid("SIMULTANEOUS requires admission=LOSSY");
      if (!x.call.offer.has_value()) invalid("SIMULTANEOUS requires finite offers");
    }
    check_call(pe, x.call);
    out.push_back(x);
  }
  return out;
}

// The list's one publication: one atomic exchange of every publishing policy's records (new ids in
// list order, then draft order), then their epochs (nothing after the exchange throws).
void publish_list(FlowLedger& ledger, std::span<ListAction> list, std::span<PendingPolicyUpdate> pending, const std::vector<std::uint8_t>& publishes) {
  std::vector<OwnerDrafts> owners;
  std::vector<std::size_t> which;
  for (std::size_t i = 0; i < list.size(); ++i) {
    if (publishes[i] == 0) continue;
    owners.push_back({list[i].engine->policy_id, pending[i].drafts, list[i].call.is_full_offer()});
    which.push_back(i);
  }
  std::vector<std::vector<FlowId>> created(owners.size());
  LedgerAccess::commit_batch(ledger, owners, created);
  for (std::size_t k = 0; k < which.size(); ++k) {
    PendingPolicyUpdate& s = pending[which[k]];
    for (std::size_t n = 0; n < s.fresh.size(); ++n) s.next_epoch.record_groups[s.fresh[n].first][s.fresh[n].second].record_counter = created[k][n].second;
    list[which[k]].engine->epoch = std::move(s.next_epoch);
  }
}

// ORDERED: each action as its scalar call on the working totals the earlier ones left (their
// carried loads and holds exchanged for the drafts), its outputs after it.
std::vector<PlacementResult> run_ordered(FlowLedger& ledger, std::span<ListAction> list, const OutputOptions& outputs) {
  const auto c = LedgerAccess::carried(ledger);
  const auto h = LedgerAccess::reserved(ledger);
  std::vector<Amount> carried(c.begin(), c.end());
  std::vector<Amount> held(h.begin(), h.end());
  LedgerTotals totals {LedgerAccess::baseline(ledger), carried, LedgerAccess::reservation_baseline(ledger), held};
  std::vector<PendingPolicyUpdate> pending(list.size());
  std::vector<std::uint8_t> publishes(list.size(), 0);
  std::vector<PlacementResult> results(list.size());
  for (std::size_t i = 0; i < list.size(); ++i) {
    PolicyEngine& pe = *list[i].engine;
    PolicyCall& rq = list[i].call;
    rq.totals = &totals;
    if (rq.operation == PolicyOperation::Place && rq.offer == Amount {0}) {
      results[i] = zero_state(pe, &totals, outputs);
      continue;
    }
    pending[i] = compute(pe, rq);
    publishes[i] = 1;
    // The working totals after it: the policy's live records out, its drafts in.
    auto exchange = [&](const FlowRecord& r, Amount sign) {
      for (const auto& [e, x] : r.edge_flows) carried[static_cast<std::size_t>(e)] += sign * x;
      if (r.edge_reservations.empty()) return;
      if (held.empty()) {
        held.assign(carried.size(), 0);
        totals.held = held;
      }
      for (const auto& [e, x] : r.edge_reservations) held[static_cast<std::size_t>(e)] += sign * x;
    };
    for (const FlowRecord* r : current_records(pe)) exchange(*r, -1);
    for (const auto& d : pending[i].drafts) exchange(d.record, 1);
    results[i] = std::move(pending[i].result);
  }
  publish_list(ledger, list, pending, publishes);
  return results;
}

// SIMULTANEOUS: every participant prepared and forwarded on one view (the ledger's residual plus
// every participant's old carried load) by one joint evaluation, demands in creation order; each
// finished as its scalar call finishes; the outputs after the whole list.
std::vector<PlacementResult> run_simultaneous(FlowLedger& ledger, std::span<ListAction> list, const OutputOptions& outputs) {
  const std::size_t n = list.size();
  const auto c = LedgerAccess::carried(ledger);
  std::vector<Amount> carried(c.begin(), c.end());
  for (const auto& x : list) {
    for (const FlowRecord* r : current_records(*x.engine)) {
      for (const auto& [e, a] : r->edge_flows) carried[static_cast<std::size_t>(e)] -= a;
    }
  }
  const LedgerTotals totals {LedgerAccess::baseline(ledger), carried, LedgerAccess::reservation_baseline(ledger), LedgerAccess::reserved(ledger)};
  // One workspace for the whole computation: a second one on this thread would void its deltas.
  ServiceCapacityWorkspace W(totals);
  enum class Kind { EmptyRoster, Entries, BySource, ShortestPathsReplace };
  struct Participant {
    Kind kind {Kind::EmptyRoster};
    PolicyEpoch fresh {};
    std::optional<ForwardingView> view {};  // its own (a Pinned replace forwards the policy's)
    PreparedShortestPaths prep {};
    std::vector<Amount> inject {};
    ProgramFlow flow {};
    [[nodiscard]] const ForwardingView* forwarded(const PolicyEngine& pe) const noexcept { return view ? &*view : &pe.view; }
  };
  std::vector<Participant> parts(n);
  std::vector<std::size_t> forwarding;
  for (std::size_t i = 0; i < n; ++i) {
    PolicyEngine& pe = *list[i].engine;
    const PolicyCall& rq = list[i].call;
    auto& part = parts[i];
    part.fresh = fresh_epoch(pe, rq.operation == PolicyOperation::Forward);
    if (pe.reachable_sources.empty() && !is<Pinned>(pe)) continue;
    if (rq.operation == PolicyOperation::Replace && is<Pinned>(pe)) {
      part.kind = Kind::Entries;
    } else if (rq.operation == PolicyOperation::Replace) {
      PolicyCallContext ctx {pe, *pe.graph, W};
      part.prep = prepare_shortest_path_program(ctx);
      part.view.emplace(make_forwarding_view(*pe.graph, part.prep.program, nullptr, pe.allowed_destinations, nullptr));
      part.kind = Kind::ShortestPathsReplace;
    } else if (tunnel_model(pe)) {
      part.view.emplace(make_forwarding_view(*pe.graph, pe.epoch.program, nullptr, pe.allowed_destinations, rq.black_holes, true));
      part.kind = Kind::Entries;
    } else {
      part.view.emplace(installed_view(pe, rq));
      part.kind = Kind::BySource;
    }
    const ForwardingView& view = *part.forwarded(pe);
    part.inject.assign(view.num_sources(), 0);
    source_shares(view, *rq.offer, part.inject);
    forwarding.push_back(i);
  }
  // The joint forwarding: the demands in their policies' creation order.
  std::sort(forwarding.begin(), forwarding.end(),
            [&](std::size_t a, std::size_t b) { return list[a].engine->policy_id < list[b].engine->policy_id; });
  std::vector<ForwardingDemand> demands;
  demands.reserve(forwarding.size());
  for (auto i : forwarding) demands.push_back({parts[i].forwarded(*list[i].engine), parts[i].inject});
  auto flows = detail::lossy(demands, available_capacity_reader(W));
  for (std::size_t k = 0; k < forwarding.size(); ++k) parts[forwarding[k]].flow = std::move(flows[k]);
  // Each participant finished as its scalar call finishes.
  std::vector<PendingPolicyUpdate> pending(n);
  for (std::size_t i = 0; i < n; ++i) {
    PolicyEngine& pe = *list[i].engine;
    const PolicyCall& rq = list[i].call;
    auto& part = parts[i];
    const ForwardingProgram* before = pe.epoch.program.get();
    PendingPolicyUpdate& s = pending[i];
    switch (part.kind) {
      case Kind::EmptyRoster: s = empty_roster(pe, part.fresh, rq, before); break;
      case Kind::Entries:
        if (rq.operation == PolicyOperation::Replace) {
          s = forward_entries(pe, part.fresh, W, rq, before, pe.view, pe.can_forward, false, nullptr, &part.flow);
        } else {
          s = forward_entries(pe, part.fresh, W, rq, before, *part.view, view_active(*part.view), is<DistributedTE>(pe), pe.epoch.plan,
                              &part.flow);
        }
        break;
      case Kind::BySource: s = finish_by_source(pe, rq, *part.view, part.flow, nullptr); break;
      case Kind::ShortestPathsReplace: {
        const PreparedReplace replaced {part.prep.program, part.prep.flows, before};
        s = finish_by_source(pe, rq, *part.view, part.flow, &replaced);
        break;
      }
    }
    // A full-offer call: its call_delivered is the whole new delivery, compared with the old one.
    Amount old = 0;
    for (const FlowRecord* r : current_records(pe)) old += r->delivered;
    s.result.call_delivered = s.result.delivered;
    s.result.delivery_delta = s.result.delivered - old;
    s.result.call_undelivered = *rq.offer - s.result.call_delivered;
    if (*rq.offer == 0) s.result.termination = Termination::Complete;
  }
  // The outputs after the whole list: one residual for every entry, each entry's boundary from its
  // own sources and mask on it.
  if (any_output(outputs)) {
    std::vector<Amount> after(carried);
    for (const auto& s : pending) {
      for (const auto& d : s.drafts) {
        for (const auto& [e, x] : d.record.edge_flows) after[static_cast<std::size_t>(e)] += x;
      }
    }
    const auto baseline = LedgerAccess::baseline(ledger);
    std::vector<Amount> residual(after.size());
    for (std::size_t e = 0; e < residual.size(); ++e) residual[e] = baseline[e] - after[e];
    for (std::size_t i = 0; i < n; ++i) {
      const PolicyEngine& pe = *list[i].engine;
      PlacementResult& out = pending[i].result;
      const auto m = residual.size();
      if (outputs.edge_reservations) {
        std::vector<Amount> held(m, 0);
        for (const auto& d : pending[i].drafts) {
          for (const auto& [e, x] : d.record.edge_reservations) held[static_cast<std::size_t>(e)] += x;
        }
        out.edge_reservations = std::move(held);
      }
      if (outputs.edge_flows) {
        std::vector<Amount> own(m, 0);
        for (const auto& d : pending[i].drafts) {
          for (const auto& [e, x] : d.record.edge_flows) own[static_cast<std::size_t>(e)] += x;
        }
        out.edge_flows = std::move(own);
      }
      if (outputs.saturated_boundary_edges) out.saturated_boundary_edges = saturated_boundary_edges_of(*pe.graph, residual, pe.sources, pe.eligibility);
      if (outputs.residual) out.residual = residual;
    }
  }
  std::vector<PlacementResult> results(n);
  for (std::size_t i = 0; i < n; ++i) results[i] = std::move(pending[i].result);
  publish_list(ledger, list, pending, std::vector<std::uint8_t>(n, 1));
  return results;
}

} // namespace
} // namespace detail

PlacementResult FlowPolicy::place(std::optional<Amount> offer, const OutputOptions& outputs) {
  return detail::run_call(*engine_, detail::PolicyCall {detail::PolicyOperation::Place, offer, {}, nullptr, &outputs, nullptr}, nullptr);
}

PlacementResult FlowPolicy::evaluate(std::optional<Amount> offer, const std::optional<std::optional<Amount>>& design_demand,
                                     const OutputOptions& outputs) const {
  const detail::PolicyCall rq {detail::PolicyOperation::Evaluate, offer, design_demand, nullptr, &outputs, nullptr};
  detail::check_call(*engine_, rq);
  return std::move(detail::compute(*engine_, rq).result);
}

PlacementResult FlowPolicy::replace(std::optional<Amount> offer, const std::optional<std::optional<Amount>>& design_demand,
                                    const OutputOptions& outputs) {
  return detail::run_call(*engine_, detail::PolicyCall {detail::PolicyOperation::Replace, offer, design_demand, nullptr, &outputs, nullptr},
                          nullptr);
}

PlacementResult FlowPolicy::forward(std::optional<Amount> offer, const Eligibility* black_holes,
                                    const OutputOptions& outputs) {
  return detail::run_call(*engine_, detail::PolicyCall {detail::PolicyOperation::Forward, offer, {}, black_holes, &outputs, nullptr},
                          nullptr);
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
PlacementResult detail::place_traced(FlowPolicy& policy, std::optional<Amount> offer, const OutputOptions& outputs,
                                     Trace& trace, const std::function<void()>* before_publish) {
  trace = Trace::map();
  return run_call(PolicyAccess::engine(policy), PolicyCall {PolicyOperation::Place, offer, {}, nullptr, &outputs, &trace}, before_publish);
}
#endif

std::vector<PlacementResult> apply_actions(FlowLedger& ledger, std::span<const Action> actions, BatchMode mode, const OutputOptions& outputs) {
  auto list = detail::validate_list(ledger, actions, mode, outputs);
  if (list.empty()) return {};
  return mode == BatchMode::Ordered ? detail::run_ordered(ledger, list, outputs) : detail::run_simultaneous(ledger, list, outputs);
}

void FlowPolicy::remove() {
  auto& pe = *engine_;
  (void)LedgerAccess::detach(pe.ledger, pe.policy_id);
  // A new epoch (nothing below throws): the same number of cohorts, every one empty.
  for (auto& c : pe.epoch.record_groups) c.clear();
  pe.epoch.program.reset();
  pe.epoch.plan.reset();
  pe.epoch.support.clear();
  pe.epoch.marginal_cost_curve.reset();
  pe.epoch.source_drops_by_source.clear();
  pe.epoch.by_tunnel.clear();
  pe.epoch.requested = Counter128 {};
  pe.epoch.offered = Counter128 {};
  pe.epoch.infinite = false;
  pe.epoch.prepared = false;
}

Amount FlowPolicy::delivered() const {
  const auto& pe = *engine_;
  Amount total = 0;
  for (const auto& c : pe.epoch.record_groups) {
    for (const auto& s : c) total += LedgerAccess::live(pe.ledger, s.record_counter).delivered;
  }
  return total;
}

std::vector<FlowRecord> FlowPolicy::records() const {
  const auto& pe = *engine_;
  std::vector<std::int64_t> ids;
  for (const auto& c : pe.epoch.record_groups) {
    for (const auto& s : c) ids.push_back(s.record_counter);
  }
  std::sort(ids.begin(), ids.end());
  std::vector<FlowRecord> out;
  out.reserve(ids.size());
  for (std::int64_t id : ids) out.push_back(LedgerAccess::live(pe.ledger, id));
  return out;
}

FlowRecord FlowPolicy::record(const FlowId& flow_id) const {
  const auto& pe = *engine_;
  const FlowLedger& fg = pe.ledger;
  // The ledger's checks first (another ledger, a retired or unknown id), then the owner.
  FlowRecord rec = fg.record(flow_id);
  if (rec.policy_id != pe.policy_id) {
    throw std::invalid_argument("FlowId (" + std::to_string(flow_id.first) + ", " + std::to_string(flow_id.second) +
                                ") is not a record of this policy");
  }
  return rec;
}

std::shared_ptr<const ForwardingProgram> FlowPolicy::program() const noexcept { return engine_->epoch.program; }

const FlowPolicyConfig& FlowPolicy::config() const noexcept { return engine_->config; }

} // namespace netgraph::core
