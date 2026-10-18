/*
  The policy engine FlowPolicy: one transaction for every routing mode and operation, published
  through the ledger as one atomic exchange per call.

  A call computes from the epoch it continues (place) or a fresh one (evaluate, replace and
  forward, on the residual with the policy's own load released) into a staged value: the result,
  the drafts of the commit and the epoch after it. place, replace and forward publish it (the
  ledger commit, then the epoch); evaluate discards it. The change flags compare the program the
  policy had installed with the call's.

  * CostOnly EqualBalanced and Weighted (place_pooled): the epoch's first positive call (or a
    replace) prepares every roster source's route on the working residual at its entry and joins
    them into the program's one pooled split context, a plan rooted at the whole roster (a router
    reached by several sources splits their added traffic once). Weighted weights come from the
    configured basis (capacity, explicit edge weights, or the free capacity read then and frozen).
    Lossless admits the largest total of the equal origin whose no-drop forwarding fits V = W +
    own; Lossy forwards the round-robin increments of the offer on W. The amounts are attributed
    back to the sources (kernels' forward()).
  * The controllers (place_controller): unrestricted Global over every eligible edge, Global with a
    support over that support trimmed, CostOnly MaxFlow over the union of the roster's routes
    (prepared as above). A minimum-cost flow of the cumulative value on V = W + own, from one
    virtual source (the optimized origin) or with source i's arc capped at its round-robin share of
    the largest total whose shares all fit (the equal origin, integer bisection), cut into paths
    per source (decompose()): the program, each path a tunnel weighted by its amount; each source's
    record is the sum of its paths.
  * Pinned (place_pinned): forwarding its program as the policy's view replays it under its mask
    with its black holes, through the forwarding engine (forwarding.hpp), one record per program
    entry; forward() (forward_installed) forwards a computing model's installed program the same
    way, with records per roster source.
  * CapacityAware (place_flows) visits the roster in ascending NodeId. Round robin gives each
    source an equal integer share of the cumulative request, with leftover units going to earlier
    sources. Offer each source the increase in its share, then run stages 1 creation, 2
    reoptimisation, 3 amounts and 4 growth on shared unused capacity W (the working residual). One
    `created` counter covers the call. It has no program.

  Work on copies and publish all changes together; an exception leaves ledger and policy intact.
  Store only changes to ledger residual in per-thread slots. Unrestricted Global reads a dense
  residual and solves over every edge (its view and flows per EdgeId); the other controllers
  gather their view over their support and keep their flows aligned with it. CapacityAware reads a
  dense residual copy; CostOnly preparation reads the working residual edge by edge.

  The transaction trace (the call's intermediate values, keyed as the tests' model of the policy
  keys them) exists in a build with the private hooks only: every
  trace block is guarded by kTraced and is compile-time dead without them.
*/
#include "netgraph/core/flow_policy.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "netgraph/core/detail/cache_access.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/flow_core.hpp"
#include "netgraph/core/detail/forwarding.hpp"
#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/detail/trace.hpp"
#include "netgraph/core/detail/traversal.hpp"

namespace netgraph::core {

namespace detail {

#if NETGRAPH_CORE_PRIVATE_HOOKS
std::atomic<std::int64_t>& max_flows_per_call() noexcept {
  static std::atomic<std::int64_t> cap {kMaxFlowsPerCall};
  return cap;
}
#endif

namespace {

// ================================================================ working residual

// W = baseline - allocated + delta; delta holds this call's releases and new allocations.
// Each call advances a scratch epoch (a numeric stamp, separate from the policy epoch).
// Read a slot only when its stamp matches; later calls ignore stale slots even after a throw.
// Clear all stamps only on wraparound (after 2^32 calls on one thread, which a long-running
// process reaches). Slots share one per-thread array, so failed growth leaves the old array
// intact.
class Working {
public:
  explicit Working(const FlowLedger& fg)
      : s_(storage()), base_(LedgerAccess::baseline(fg)), held_(LedgerAccess::allocated(fg)) {
    if (s_.slot.size() < base_.size()) s_.slot.resize(base_.size());
    if (++s_.epoch == 0) {
      for (auto& x : s_.slot) x.stamp = 0;
      s_.epoch = 1;
    }
  }
  Working(const Working&) = delete;
  Working& operator=(const Working&) = delete;

  [[nodiscard]] Amount at(EdgeId e) const noexcept {
    const auto ei = static_cast<std::size_t>(e);
    const auto& x = s_.slot[ei];
    return base_[ei] - held_[ei] + (x.stamp == s_.epoch ? x.delta : 0);
  }
  void add(EdgeId e, Amount v) noexcept {
    auto& x = s_.slot[static_cast<std::size_t>(e)];
    if (x.stamp != s_.epoch) {
      x.stamp = s_.epoch;
      x.delta = 0;
    }
    x.delta += v;
  }
  // Unused capacity for each edge of `members`, in the same order (member-aligned).
  void gather(std::span<const EdgeId> members, std::vector<Amount>& out) const {
    out.resize(members.size());
    for (std::size_t m = 0; m < members.size(); ++m) out[m] = at(members[m]);
  }
  // The dense working residual (searches, Global, traces).
  void dense(std::vector<Amount>& out) const {
    out.resize(base_.size());
    for (std::size_t e = 0; e < out.size(); ++e) {
      const auto& x = s_.slot[e];
      out[e] = base_[e] - held_[e] + (x.stamp == s_.epoch ? x.delta : 0);
    }
  }
  // The policy's own carried load released (the view of a full-offer call): its live records'
  // edges added back.
  void release_records(std::span<const FlowRecord* const> records) noexcept {
    for (const FlowRecord* r : records) {
      for (const auto& [e, x] : r->edges) add(e, x);
    }
  }
  // W += a (release) or W -= a (consume) over member-aligned amounts.
  void release(std::span<const EdgeId> members, std::span<const Amount> a) noexcept {
    for (std::size_t m = 0; m < members.size(); ++m) {
      if (a[m] != 0) add(members[m], a[m]);
    }
  }
  void consume(std::span<const EdgeId> members, std::span<const Amount> a) noexcept {
    for (std::size_t m = 0; m < members.size(); ++m) {
      if (a[m] != 0) add(members[m], -a[m]);
    }
  }

private:
  struct Slot {
    Amount delta {0};
    std::uint32_t stamp {0};
  };
  struct Storage {
    std::vector<Slot> slot;
    std::uint32_t epoch {0};
  };
  static Storage& storage() {
    thread_local Storage st;
    return st;
  }
  Storage& s_;
  std::span<const Amount> base_;
  std::span<const Amount> held_;
};

// ================================================================ small helpers

[[noreturn]] void invalid(const std::string& what) { throw std::invalid_argument(what); }

// The policy's routing mode (the tag of the Routing variant).
template <class R>
bool is(const PolicyEngine& pe) noexcept {
  return std::holds_alternative<R>(pe.config.routing);
}

std::string str(std::int64_t v) { return std::to_string(v); }

const char* outcome_name(Termination t) noexcept {
  switch (t) {
    case Termination::Complete: return "complete";
    case Termination::LimitFlows: return "limit_flows";
    case Termination::LimitDiminishing: return "limit_diminishing";
    case Termination::NoRoute: return "no_route";
    case Termination::Exhausted: return "exhausted";
    case Termination::NoActiveWeight: return "no_active_weight";
  }
  return "?";
}

// The reduction rank of the call's termination (higher wins among unserved sources).
int rank(Termination t) noexcept {
  switch (t) {
    case Termination::LimitFlows: return 4;
    case Termination::LimitDiminishing: return 3;
    case Termination::Exhausted: return 2;
    case Termination::NoRoute: return 1;
    case Termination::Complete:
    case Termination::NoActiveWeight: return 0;  // not an outcome of the reduction
  }
  return 0;
}

Trace amount_or_inf(std::optional<Amount> v) { return v.has_value() ? Trace::integer(*v) : Trace::str("inf"); }

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

// Round robin of a cumulative 128-bit amount (the request, or a Lossy demand's offer): the
// per-source increment of the share when the cumulative amount goes from R to R + v over k
// sources, offer(i, increment) for i in [0, k). It depends only on r = R mod k:
// floor((r + v) / k) + [i < (r + v) mod k] - [i < r].
template <class Offer>
void cumulative_offers(const Counter128& requested, Amount v, std::size_t k, Offer&& offer) {
  const auto kk = static_cast<std::uint64_t>(k);
  const std::uint64_t two64 = ((std::numeric_limits<std::uint64_t>::max() % kk) + 1) % kk;
  const std::uint64_t r = ((requested.hi % kk) * two64 % kk + requested.lo % kk) % kk;
  const auto total = static_cast<std::uint64_t>(v) + r;
  const auto base = static_cast<Amount>(total / kk);
  const std::uint64_t extra = total % kk;
  for (std::size_t i = 0; i < k; ++i) offer(i, base + (i < extra ? 1 : 0) - (i < r ? 1 : 0));
}

// ================================================================ working flows

struct WorkFlow {
  std::int64_t id {-1};  // -1: created in this call
  NodeId src {0};
  std::shared_ptr<const KernelPlan> plan {};
  std::shared_ptr<const PredDAG> bundle {};
  Cost cost {kUnreachable};
  Amount T {0};                 // Lossless: the cumulative amount (= delivered)
  std::vector<Amount> alloc {}; // member-aligned with plan->members (cumulative carried for Lossy)
  std::optional<Counter128> offered {};  // Lossy
  Amount delivered {0};                  // Lossy (Lossless: T)
  std::vector<Counter128> drops {};      // Lossy, member-aligned, cumulative
  std::int64_t seq {0};                  // creation order within the call (new flows)
  [[nodiscard]] bool is_new() const noexcept { return id < 0; }
  [[nodiscard]] std::span<const EdgeId> members() const noexcept { return plan->members; }
};

struct SourceRun {
  std::size_t ri {0};
  NodeId src {0};
  std::optional<Amount> offer {};  // nullopt: +inf
  Amount P {0};
  std::optional<Amount> B {};
  Amount T_entry {0};
  std::vector<WorkFlow> flows {};
  bool capped {false};
  std::optional<Termination> limit {};
  std::vector<Amount> gains {};
  Amount injected {0};
  Trace* tr {nullptr};  // this source's trace map
};

struct Call {
  const PolicyEngine& pe;
  const Graph& g;
  Working& W;
  std::int64_t cap {kMaxFlowsPerCall};
  std::int64_t created {0};
  bool routes_changed {false};
  bool weights_changed {false};
  std::int64_t seq {0};
  std::vector<std::pair<EdgeId, Amount>> drops {};
  Amount offered {0};      // Lossy: this call's offers (the first-hop injection of +inf)
  Amount source_loss {0};  // Lossy: this call's shares of flows without a member
  Trace* trace {nullptr};
};

// ================================================================ routes and plans

struct Route {
  std::vector<EdgeId> edges {};  // ascending
  Cost cost {kUnreachable};
};

// A CapacityAware flow's plan on its route (the record's bundle) with its membership captured on
// the creation view (Weighted with priority Residual: `dense_view` when the search read one, else
// the working residual W); with `cost`, the route's minimum walk cost.
std::shared_ptr<const KernelPlan> plan_on(const Call& c, NodeId src, const PredDAG& route,
                                          const std::vector<Amount>* dense_view, Cost* cost) {
  const auto& pe = c.pe;
  std::vector<Amount> creation;
  // CapacityAware's Weighted weights come from its selection's priority (Capacity or Residual).
  WeightSource weights {WeightBasis::Capacity, {}};
  if (*pe.config.placement == FlowPlacement::Weighted && pe.selection.priority == EdgePriority::Residual) {
    const auto edges = PredDAGBuilder::ascending(route);
    creation.resize(edges.size());
    for (std::size_t i = 0; i < edges.size(); ++i) {
      creation[i] = dense_view != nullptr ? (*dense_view)[static_cast<std::size_t>(edges[i])] : c.W.at(edges[i]);
    }
    weights = WeightSource {WeightBasis::ResidualSnapshot, creation};
  }
  return std::make_shared<const KernelPlan>(make_plan(c.g, *pe.config.placement, pe.config.admission, weights, src,
                                                      pe.dests_ok, route, cost));
}

// The priority values of the creation view W, as the selection reads them (CostOnly). The
// PriorityFn of `on_w` points at this view: it must outlive every call through that PriorityFn.
struct WorkingPriority {
  const Working* W;
  std::span<const Amount> capacity;
  EdgePriority priority;
  [[nodiscard]] PriorityFn on_w() const noexcept {
    return {this, [](const void* ctx, EdgeId e) noexcept {
              const auto& v = *static_cast<const WorkingPriority*>(ctx);
              return priority_value(v.priority, v.W->at(e), v.capacity[static_cast<std::size_t>(e)]);
            }};
  }
};

// The record bundle of a route: `share` itself when it is the whole route, else a new DAG.
std::shared_ptr<const PredDAG> bundle_of(const Graph& g, std::span<const EdgeId> route,
                                         const std::shared_ptr<const PredDAG>& share) {
  if (share && PredDAGBuilder::ascending(*share).size() == route.size()) return share;
  return std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, route));
}

// CostOnly: keep static shortest-path support reachable from src, then apply selection on W.
// For multi_edge = false (Collapse), choose one edge per (u, v) by (cost, higher priority under
// PreferHigherPriority, lowest EdgeId). For multipath = false, choose the canonical single
// path in stages. Both use W's priority values without excluding edges for low residual.
Route cost_only_route(const Call& c, NodeId src) {
  const auto& g = c.g;
  const auto& pe = c.pe;
  const auto& sel = pe.selection;
  const WorkingPriority view {&c.W, g.capacity_units(), sel.priority};
  auto r = support_route(g, PredDAGBuilder::ascending(*pe.support), src, pe.dests_ok, sel, view.on_w());
  Route out;
  out.edges = std::move(r.edges);
  if (!sel.multipath && !out.edges.empty()) out.cost = r.cost;
  return out;
}

// CostOnly Weighted's static weight source (Capacity or Explicit; ResidualSnapshot reads the view
// at preparation and is made there).
WeightSource cost_only_weights(const PolicyEngine& pe) noexcept {
  if (pe.config.weights == WeightBasis::Explicit) return {WeightBasis::Explicit, *pe.config.edge_weights};
  return {pe.config.weights.value_or(WeightBasis::Capacity), {}};
}

// CostOnly's creation products in one pass (kernels' cost_only_creation, identical to
// cost_only_route, bundle_of and the plan where it applies); nullopt elsewhere.
std::optional<CreationProducts> creation_products(const Call& c, NodeId src) {
  const auto& pe = c.pe;
  const auto& sel = pe.selection;
  const WorkingPriority view {&c.W, c.g.capacity_units(), sel.priority};
  return cost_only_creation(c.g, pe.support, src, pe.dests_ok, sel, view.on_w(), *pe.config.placement,
                            pe.config.admission, cost_only_weights(pe));
}

// CapacityAware: the canonical route on the selection view `view` (dense): static
// eligibility, residual >= max(1, min_flow_bandwidth), Collapse on the view, the tight DAG or the
// staged single path toward the allowed destinations, then the max_path_cost ceiling.
std::optional<Route> search(const Call& c, std::size_t ri, std::span<const Amount> view) {
  const auto& pe = c.pe;
  const auto& g = c.g;
  const auto& ca = std::get<CapacityAware>(pe.config.routing);
  const auto& sel = pe.selection;
  const Gate gate {pe.elig, view, std::max<Amount>(1, ca.min_flow_bandwidth)};
  auto found = route_search(g, gate, sel, pe.roster[ri], pe.dests_ok);
  if (!found.has_value()) return std::nullopt;
  Route out;
  out.edges = std::move(found->edges);
  out.cost = found->cost;
  if (ca.max_path_cost.has_value()) {
    if (const auto* a = std::get_if<Absolute>(&*ca.max_path_cost)) {
      if (out.cost > a->cost) return std::nullopt;
    } else {
      // The roster holds the sources with a static route, so the baseline is a positive cost.
      if (!stretch_le(out.cost, pe.baseline[ri], std::get<Stretch>(*ca.max_path_cost).factor)) return std::nullopt;
    }
  }
  return out;
}

// ================================================================ amounts of one flow

// T_max of a flow on a member-aligned view (the trace's bounds): MaxFlow's uncapped step from 0,
// else the exact search.
Amount flow_t_max(const Graph& g, const KernelPlan& plan, std::span<const Amount> view_m) {
  if (plan.members.empty()) return 0;
  if (plan.placement == FlowPlacement::MaxFlow) return lossless_step(g, plan, 0, view_m, {}, std::nullopt).total;
  return t_max(plan, view_m).value;
}

// ================================================================ the cohort (EqualBalanced)

// Split T by round robin over the source's EqualBalanced flows in their order. Allocate each
// share using that flow's plan, then sum edge amounts across the union of their members.
class CohortEval {
public:
  CohortEval(const std::vector<const KernelPlan*>& plans, const Working& W) : plans_(plans) {
    for (const auto* p : plans_) uni_.insert(uni_.end(), p->members.begin(), p->members.end());
    std::sort(uni_.begin(), uni_.end());
    uni_.erase(std::unique(uni_.begin(), uni_.end()), uni_.end());
    idx_.resize(plans_.size());
    for (std::size_t i = 0; i < plans_.size(); ++i) {
      idx_[i].resize(plans_[i]->members.size());
      for (std::size_t m = 0; m < plans_[i]->members.size(); ++m) {
        idx_[i][m] = static_cast<std::int32_t>(std::lower_bound(uni_.begin(), uni_.end(), plans_[i]->members[m]) -
                                               uni_.begin());
      }
    }
    W.gather(uni_, view_);
    acc_.assign(uni_.size(), 0);
    // The search's upper bound: the view summed over the source members of every flow (each
    // plan has the one root, the source).
    std::vector<std::uint8_t> src_member(uni_.size(), 0);
    for (std::size_t i = 0; i < plans_.size(); ++i) {
      const auto p = plans_[i]->root_pos.front();
      if (p < 0) continue;
      for (auto k = plans_[i]->out_begin[static_cast<std::size_t>(p)]; k < plans_[i]->out_end[static_cast<std::size_t>(p)]; ++k) {
        src_member[static_cast<std::size_t>(idx_[i][static_cast<std::size_t>(k)])] = 1;
      }
    }
    for (std::size_t u = 0; u < uni_.size(); ++u) {
      if (src_member[u] != 0) hi_ += view_[u];
    }
  }
  [[nodiscard]] Amount upper() const noexcept { return hi_; }

  // Whether the cohort allocation at total T fits the view.
  bool fits(Amount T) {
    std::fill(acc_.begin(), acc_.end(), 0);
    const auto k = plans_.size();
    std::vector<Amount> parts(k);
    round_robin(T, parts);
    for (std::size_t i = 0; i < k; ++i) {
      if (parts[i] == 0) continue;
      tmp_.resize(plans_[i]->members.size());
      allocate(*plans_[i], parts[i], tmp_);
      for (std::size_t m = 0; m < tmp_.size(); ++m) {
        const auto u = static_cast<std::size_t>(idx_[i][m]);
        acc_[u] += tmp_[m];
        if (acc_[u] > view_[u]) return false;
      }
    }
    return true;
  }
  // The largest T in [lo, hi] that fits, given fits(lo) (bisection, exact; at most 62 probes).
  Amount largest(Amount lo, Amount hi) {
    while (lo < hi) {
      const Amount mid = lo + (hi - lo + 1) / 2;
      if (fits(mid)) {
        lo = mid;
      } else {
        hi = mid - 1;
      }
    }
    return lo;
  }
  // The full T_max (the trace's `bound`).
  Amount t_max() { return largest(0, hi_); }

private:
  const std::vector<const KernelPlan*>& plans_;
  std::vector<EdgeId> uni_ {};
  std::vector<std::vector<std::int32_t>> idx_ {};
  std::vector<Amount> view_ {};
  std::vector<Amount> acc_ {};
  std::vector<Amount> tmp_ {};
  Amount hi_ {0};
};

// Splits T over the cohort by round robin; a flow whose share changed (or is new) gets the
// allocation of its new share (the amounts stage of EqualBalanced).
void set_cohort(std::vector<WorkFlow>& flows, Amount T) {
  std::vector<Amount> parts(flows.size());
  round_robin(T, parts);
  for (std::size_t i = 0; i < flows.size(); ++i) {
    auto& f = flows[i];
    if (parts[i] != f.T || f.is_new()) {
      std::vector<Amount> a(f.plan->members.size());
      allocate(*f.plan, parts[i], a);
      f.alloc = std::move(a);
    }
    f.T = parts[i];
  }
}

std::vector<const KernelPlan*> plans_of(const std::vector<WorkFlow>& flows) {
  std::vector<const KernelPlan*> out;
  out.reserve(flows.size() + 1);
  for (const auto& f : flows) out.push_back(f.plan.get());
  return out;
}

// ---------------------------------------------------------------- outputs

// The optional outputs of a call whose policy allocation becomes `after` (the records' edge
// lists), computed before the commit (so nothing allocates once the ledger changed):
// edge_flows = the policy's own allocation after the call, residual = the ledger's after it
// (now + own before - own after), bottleneck on that residual.
void prepare_outputs(const PolicyEngine& pe, const OutputOptions& outputs,
                     const std::vector<const std::vector<std::pair<EdgeId, Amount>>*>& after, PlacementResult& out) {
  if (!outputs.edge_flows && !outputs.residual && !outputs.bottleneck) return;
  const auto m = static_cast<std::size_t>(pe.graph->num_edges());
  std::vector<Amount> own(m, 0);
  for (const auto* edges : after) {
    for (const auto& [e, x] : *edges) own[static_cast<std::size_t>(e)] += x;
  }
  if (outputs.residual || outputs.bottleneck) {
    std::vector<Amount> res;
    LedgerAccess::residual_into(pe.fg, res);
    for (const auto& cohort : pe.epoch.cohorts) {
      for (const auto& s : cohort) {
        for (const auto& [e, x] : LedgerAccess::live(pe.fg, s.id).edges) res[static_cast<std::size_t>(e)] += x;
      }
    }
    for (std::size_t e = 0; e < m; ++e) res[e] -= own[e];
    if (outputs.bottleneck) out.bottleneck = bottleneck_of(*pe.graph, res, pe.sources, pe.elig);
    if (outputs.residual) out.residual = std::move(res);
  }
  if (outputs.edge_flows) out.edge_flows = std::move(own);
}

// The current records' edge lists (a call that commits nothing).
std::vector<const std::vector<std::pair<EdgeId, Amount>>*> current_edges(const PolicyEngine& pe) {
  std::vector<const std::vector<std::pair<EdgeId, Amount>>*> out;
  for (const auto& cohort : pe.epoch.cohorts) {
    for (const auto& s : cohort) out.push_back(&LedgerAccess::live(pe.fg, s.id).edges);
  }
  return out;
}

bool any_output(const OutputOptions& outputs) noexcept {
  return outputs.edge_flows || outputs.residual || outputs.bottleneck;
}

std::vector<const std::vector<std::pair<EdgeId, Amount>>*> draft_edges(const std::vector<Draft>& drafts) {
  std::vector<const std::vector<std::pair<EdgeId, Amount>>*> out;
  out.reserve(drafts.size());
  for (const auto& d : drafts) out.push_back(&d.record.edges);
  return out;
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
      n += r->drops.size();
      dropping += r->drops.empty() ? 0 : 1;
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
    for (const auto& [e, x] : r->drops) {
      accumulate(row.dropped, x, "drops");
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
    row.source_loss = minus(minus(*row.offered, c128(row.delivered)), row.dropped);
  }
  if (dropping > 1) std::sort(drops.begin(), drops.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  out.drops.reserve(drops.size());
  Counter128 dropped {};
  for (const auto& [e, x] : drops) {
    if (out.drops.empty() || out.drops.back().first != e) out.drops.emplace_back(e, Counter128 {});
    accumulate(out.drops.back().second, x, "drops");
    accumulate(dropped, x, "drops");
  }
  out.offered = offered;
  out.admitted = offered;
  out.unadmitted = Counter128 {};
  out.source_loss = minus(minus(offered, c128(out.delivered)), dropped);
}

// The exact sum of cost * flow over the records' edges (a controller's total cost; below 2^124:
// the costs and the capacities each sum below 2^62).
Counter128 total_cost_of(const Graph& g, std::span<const FlowRecord* const> after) {
  const auto cost = g.cost();
  U128 sum {};
  for (const FlowRecord* r : after) {
    for (const auto& [e, x] : r->edges) {
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
  for (const auto& cohort : pe.epoch.cohorts) {
    for (const auto& s : cohort) out.push_back(&LedgerAccess::live(pe.fg, s.id));
  }
  return out;
}

// The drafts' records as a list, in per-thread storage reused by every call (valid until this
// thread's next draft_records; it grows before any commit).
std::span<const FlowRecord* const> draft_records(const std::vector<Draft>& drafts) {
  thread_local std::vector<const FlowRecord*> out;
  out.clear();
  for (const auto& d : drafts) out.push_back(&d.record);
  return out;
}

// Global (with or without a support), or CostOnly MaxFlow: a controller.
bool controller(const PolicyEngine& pe) noexcept {
  return is<Global>(pe) || (is<CostOnly>(pe) && *pe.config.placement == FlowPlacement::MaxFlow);
}

// Unrestricted Global: the controller over every eligible edge (dense arrays, the curve, the cut).
bool dense_global(const PolicyEngine& pe) noexcept {
  return is<Global>(pe) && !std::get<Global>(pe.config.routing).support;
}

// ================================================================ the transaction

Trace route_trace(std::span<const EdgeId> edges) { return Trace::ints(edges); }

Trace alloc_trace(const WorkFlow& f) {
  Trace t = Trace::list();
  const auto members = f.members();
  for (std::size_t m = 0; m < members.size(); ++m) {
    if (f.alloc[m] == 0) continue;
    Trace p = Trace::list();
    p.push(Trace::integer(members[m]));
    p.push(Trace::integer(f.alloc[m]));
    t.push(std::move(p));
  }
  return t;
}

// A new working flow on `bundle` (its record's bundle) with the route's cost (kUnreachable:
// computed with the plan) and its plan captured on the creation view, unless the route's static
// plan is given.
WorkFlow new_flow(Call& c, NodeId src, std::shared_ptr<const PredDAG> bundle, Cost cost,
                  const std::vector<Amount>* dense_view, std::shared_ptr<const KernelPlan> plan = nullptr) {
  WorkFlow f;
  f.src = src;
  f.bundle = std::move(bundle);
  f.cost = cost;
  f.plan = plan ? std::move(plan) : plan_on(c, src, *f.bundle, dense_view, cost == kUnreachable ? &f.cost : nullptr);
  f.alloc.assign(f.plan->members.size(), 0);
  if (c.pe.config.admission == Admission::Lossy) {
    f.offered = Counter128 {};
    f.drops.assign(f.plan->members.size(), Counter128 {});
  }
  return f;
}

// ---------------------------------------------------------------- stage 1: creation

void stage1(Call& c, SourceRun& st) {
  const auto& pe = c.pe;
  Trace tr = Trace::map();
  {
    const auto& ca = std::get<CapacityAware>(pe.config.routing);
    std::vector<Amount> S;
    c.W.dense(S);
    std::optional<Amount> rem = st.offer;
    std::optional<Amount> target;
    if (*pe.config.placement == FlowPlacement::EqualBalanced && st.offer.has_value()) {
      target = (st.P + *st.offer) / ca.flow_count;
    }
    Trace reservations = Trace::list();
    const std::int64_t budget = std::min<std::int64_t>(ca.flow_count, std::max<std::int64_t>(0, c.cap - c.created));
    bool stopped = false;
    std::vector<Amount> view_m;
    std::vector<Amount> zeros;
    for (std::int64_t i = 0; i < budget; ++i) {
      auto route = search(c, st.ri, S);
      if (!route.has_value()) {
        stopped = true;
        break;
      }
      WorkFlow f = new_flow(c, st.src, bundle_of(c.g, route->edges, nullptr), route->cost, &S);
      const auto members = f.members();
      view_m.resize(members.size());
      for (std::size_t m = 0; m < members.size(); ++m) view_m[m] = S[static_cast<std::size_t>(members[m])];
      zeros.assign(members.size(), 0);
      std::optional<Amount> request = *pe.config.placement == FlowPlacement::EqualBalanced ? target : rem;
      auto step = lossless_step(c.g, *f.plan, 0, view_m, zeros, request);
      const Amount t = step.total;
      if (*pe.config.placement != FlowPlacement::EqualBalanced && rem.has_value()) *rem -= t;
      if (kTraced && c.trace) {
        Trace row = Trace::map();
        row.set("route", route_trace(route->edges));
        row.set("S", Trace::ints(S));
        row.set("T_bound", Trace::integer(flow_t_max(c.g, *f.plan, view_m)));
        row.set("reserved", Trace::integer(t));
        reservations.push(std::move(row));
      }
      for (std::size_t m = 0; m < members.size(); ++m) S[static_cast<std::size_t>(members[m])] -= step.alloc[m];
      f.seq = ++c.seq;
      st.flows.push_back(std::move(f));
      ++c.created;
    }
    st.capped = budget < ca.flow_count && !stopped;
    if (kTraced && c.trace) {
      tr.set("eb_target", target.has_value() ? Trace::integer(*target) : Trace::null());
      tr.set("reservations", std::move(reservations));
      tr.set("created_cap_blocked", Trace::boolean(st.capped));
    }
  }
  if (!st.flows.empty()) c.routes_changed = true;
  if (kTraced && st.tr) st.tr->set("stage1", std::move(tr));
}

// ---------------------------------------------------------------- stage 2: reoptimisation

void stage2(Call& c, SourceRun& st) {
  Trace rows = Trace::list();
  std::vector<Amount> Wi;
  std::vector<Amount> view_m;
  std::vector<Amount> zeros;
  for (auto& f : st.flows) {
    if (f.is_new()) continue;
    c.W.dense(Wi);
    const auto members = f.members();
    for (std::size_t m = 0; m < members.size(); ++m) Wi[static_cast<std::size_t>(members[m])] += f.alloc[m];
    auto cand = search(c, st.ri, Wi);
    Trace row = Trace::map();
    if (kTraced && c.trace) {
      row.set("flow", Trace::ints(std::vector<std::int64_t> {static_cast<std::int64_t>(c.pe.fg.uid()), f.id}));
      row.set("T", Trace::integer(f.T));
      row.set("candidate", cand.has_value() ? route_trace(cand->edges) : Trace::null());
    }
    if (!cand.has_value() || cand->edges == f.plan->route) {
      if (kTraced && c.trace) {
        row.set("decision", Trace::str("kept"));
        rows.push(std::move(row));
      }
      continue;
    }
    WorkFlow trial = new_flow(c, st.src, bundle_of(c.g, cand->edges, nullptr), cand->cost, &Wi);
    const auto tm = trial.members();
    view_m.resize(tm.size());
    for (std::size_t m = 0; m < tm.size(); ++m) view_m[m] = Wi[static_cast<std::size_t>(tm[m])];
    zeros.assign(tm.size(), 0);
    // The candidate's allocation at T_i (adopted iff it delivers T_i exactly).
    auto step = lossless_step(c.g, *trial.plan, 0, view_m, zeros, f.T);
    if (step.total < f.T) {
      if (kTraced && c.trace) {
        row.set("decision", Trace::str("rejected: delivers less than T_i"));
        row.set("bound", Trace::integer(flow_t_max(c.g, *trial.plan, view_m)));
        rows.push(std::move(row));
      }
      continue;
    }
    // Adopted: a replacement at T_i; W := W_i - alloc_c.
    c.W.release(members, f.alloc);
    trial.alloc = std::move(step.alloc);
    c.W.consume(tm, trial.alloc);
    f.plan = std::move(trial.plan);
    f.bundle = std::move(trial.bundle);
    f.cost = trial.cost;
    f.alloc = std::move(trial.alloc);
    c.routes_changed = true;
    if (*c.pe.config.placement == FlowPlacement::Weighted) c.weights_changed = true;
    if (kTraced && c.trace) {
      row.set("decision", Trace::str("adopted"));
      row.set("alloc", alloc_trace(f));
      rows.push(std::move(row));
    }
  }
  if (kTraced && st.tr) st.tr->set("stage2", std::move(rows));
}

// ---------------------------------------------------------------- stage 3: amounts

// The member-aligned view of stage 3, per thread: every use gathers it whole first, and no use
// nests another, so its contents never outlive one flow's step.
std::vector<Amount>& view_scratch() {
  thread_local std::vector<Amount> v;
  return v;
}

void stage3_seq(Call& c, SourceRun& st) {
  Trace rows = Trace::list();
  std::vector<Amount>& view_m = view_scratch();
  for (auto& f : st.flows) {
    const auto members = f.members();
    c.W.release(members, f.alloc);  // V = W + own
    c.W.gather(members, view_m);
    const Amount T0 = f.T;
    auto step = lossless_step(c.g, *f.plan, f.T, view_m, f.alloc, st.B);
    f.T = step.total;
    f.alloc = std::move(step.alloc);
    c.W.consume(members, f.alloc);  // W = V - alloc
    if (st.B.has_value()) *st.B -= f.T - T0;
    if (kTraced && c.trace) {
      Trace row = Trace::map();
      row.set("T", Trace::ints(std::vector<Amount> {T0, f.T}));
      row.set("curve", f.plan->placement == FlowPlacement::MaxFlow ? Trace::pairs(step.cost_distribution) : Trace::null());
      row.set("alloc", alloc_trace(f));
      rows.push(std::move(row));
    }
  }
  if (kTraced && st.tr) st.tr->set("stage3", std::move(rows));
}

void stage3_eb(Call& c, SourceRun& st) {
  auto& flows = st.flows;
  Amount T = 0;
  for (const auto& f : flows) T += f.T;
  for (const auto& f : flows) c.W.release(f.members(), f.alloc);  // V = W + cohort allocation
  Amount t_new = T;
  std::optional<Amount> bound;
  if (flows.size() == 1) {
    auto& f = flows.front();
    std::vector<Amount>& view_m = view_scratch();
    c.W.gather(f.members(), view_m);
    auto step = lossless_step(c.g, *f.plan, f.T, view_m, f.alloc, st.B);
    t_new = step.total;
    f.alloc = std::move(step.alloc);
    f.T = t_new;
    if (kTraced && c.trace) bound = flow_t_max(c.g, *f.plan, view_m);
  } else {
    const auto plans = plans_of(flows);
    CohortEval cohort(plans, c.W);
    const Amount target = st.B.has_value() ? T + *st.B : cohort.upper();
    if (target > T) {
      t_new = cohort.fits(target) ? target : cohort.largest(T, target - 1);
    }
    if (kTraced && c.trace) bound = cohort.t_max();
    set_cohort(flows, t_new);
  }
  for (const auto& f : flows) c.W.consume(f.members(), f.alloc);
  if (st.B.has_value()) *st.B -= t_new - T;
  if (kTraced && st.tr) {
    Trace row = Trace::map();
    row.set("T", Trace::ints(std::vector<Amount> {T, t_new}));
    row.set("bound", bound.has_value() ? Trace::integer(*bound) : Trace::null());
    Trace shares = Trace::list();
    for (const auto& f : flows) shares.push(Trace::integer(f.T));
    row.set("shares", std::move(shares));
    st.tr->set("stage3", std::move(row));
  }
}

// ---------------------------------------------------------------- stage 4: growth

// fraction * v as the reduced rational "num/den" (or "num"): the exact text the transaction trace
// records for the diminishing guard's threshold.
std::string fraction_text(double fraction, Amount v) {
  const auto [m, x] = dyadic(fraction);
  U128 num = mul_u64(m, static_cast<std::uint64_t>(v));
  int k = -x;  // den = 2^k
  while (k > 0 && (num.lo & 1U) == 0 && !(num == U128 {0, 0})) {
    num = U128 {num.hi >> 1, (num.lo >> 1) | (num.hi << 63)};
    --k;
  }
  if (num == U128 {0, 0}) return "0";
  auto dec = [](U128 n) {
    std::string out;
    while (!(n == U128 {0, 0})) {
      const std::uint64_t hi_q = n.hi / 10;
      const std::uint64_t hi_r = n.hi % 10;
      const std::uint64_t lo_q = div_u128(U128 {hi_r, n.lo}, 10);
      const std::uint64_t lo_r = n.lo - lo_q * 10;  // the remainder of (hi_r * 2^64 + lo) / 10
      out.push_back(static_cast<char>('0' + lo_r));
      n = U128 {hi_q, lo_q};
    }
    std::reverse(out.begin(), out.end());
    return out;
  };
  if (k == 0) return dec(num);
  if (k <= 126) return dec(num) + "/" + dec(k >= 64 ? U128 {std::uint64_t {1} << (k - 64), 0} : U128 {0, std::uint64_t {1} << k});
  return dec(num) + "/2**" + std::to_string(k);
}

// Compare gains.sum < fraction * v exactly, reading fraction as its exact stored double value;
// v < 2^62.
bool below_fraction(Amount recent, double fraction, Amount v) {
  const auto [m, x] = dyadic(fraction);  // x < 0 since fraction < 1
  const U128 rhs = mul_u64(m, static_cast<std::uint64_t>(v));  // below 2^115
  // recent < rhs * 2^x  <=>  recent * 2^-x < rhs; with -x large the left side exceeds 2^128
  // unless recent is 0.
  if (recent == 0) return !(rhs == U128 {0, 0});
  const U128 r {0, static_cast<std::uint64_t>(recent)};
  if (bits(r) - x > 127) return false;
  return shl(r, -x) < rhs;
}

void stage4(Call& c, SourceRun& st) {
  const auto& pe = c.pe;
  const auto& ca = std::get<CapacityAware>(pe.config.routing);
  const auto v = st.offer;
  Trace steps = Trace::list();
  std::int32_t accepted = 0;
  std::vector<Amount> Wd;
  std::vector<Amount> view_m;
  std::vector<Amount> zeros;
  auto stop = [&](const char* why) {
    if (kTraced && c.trace) {
      Trace s = Trace::map();
      s.set("stop", Trace::str(why));
      steps.push(std::move(s));
    }
  };
  while (!st.B.has_value() || *st.B > 0) {
    if (c.created >= c.cap) {
      st.limit = Termination::LimitFlows;
      stop("kMaxFlowsPerCall");
      break;
    }
    if (ca.max_flow_count.has_value() && static_cast<std::int64_t>(st.flows.size()) >= *ca.max_flow_count) {
      if (*ca.max_flow_count > ca.flow_count) st.limit = Termination::LimitFlows;
      if (kTraced && c.trace) {
        Trace s = Trace::map();
        s.set("stop", Trace::str("max_flow_count"));
        s.set("k", Trace::integer(static_cast<std::int64_t>(st.flows.size())));
        steps.push(std::move(s));
      }
      break;
    }
    c.W.dense(Wd);
    auto cand = search(c, st.ri, Wd);
    if (!cand.has_value()) {
      stop("no candidate");
      break;
    }
    if (pe.selection.multipath) {
      // A multipath DAG already spans every equal-cost path, so a second
      // flow on a structurally equal DAG gains only round-robin remainders: growth stops. A
      // repeated single path is a legitimate tunnel.
      bool equal = false;
      for (const auto& f : st.flows) equal = equal || f.plan->route == cand->edges;
      if (equal) {
        if (kTraced && c.trace) {
          Trace s = Trace::map();
          s.set("stop", Trace::str("candidate equals a cohort route"));
          s.set("route", route_trace(cand->edges));
          steps.push(std::move(s));
        }
        break;
      }
    }
    WorkFlow f = new_flow(c, st.src, bundle_of(c.g, cand->edges, nullptr), cand->cost, &Wd);
    Trace step = Trace::map();
    if (kTraced && c.trace) {
      step.set("route", route_trace(cand->edges));
      step.set("W", Trace::ints(Wd));
    }
    Amount gain = 0;
    if (*pe.config.placement == FlowPlacement::EqualBalanced) {
      Amount T = 0;
      for (const auto& x : st.flows) T += x.T;
      for (const auto& x : st.flows) c.W.release(x.members(), x.alloc);  // V = W + cohort allocation
      auto plans = plans_of(st.flows);
      plans.push_back(f.plan.get());
      CohortEval trial(plans, c.W);
      const Amount target = v.has_value() ? st.T_entry + *v : trial.upper();
      Amount t_cand = T;  // rejected unless the trial cohort carries more than T
      if (target > T && trial.fits(T + 1)) {
        t_cand = trial.fits(target) ? target : trial.largest(T + 1, target - 1);
      }
      if (kTraced && c.trace) {
        const Amount bound = trial.t_max();
        step.set("T", Trace::integer(T));
        step.set("T_candidate", Trace::integer(v.has_value() ? std::min(st.T_entry + *v, bound) : bound));
        step.set("bound", Trace::integer(bound));
      }
      if (t_cand <= T) {
        for (const auto& x : st.flows) c.W.consume(x.members(), x.alloc);  // W restored
        if (kTraced && c.trace) {
          step.set("decision", Trace::str("rejected: cohort total not larger"));
          steps.push(std::move(step));
        }
        break;
      }
      f.seq = ++c.seq;
      st.flows.push_back(std::move(f));
      set_cohort(st.flows, t_cand);
      for (const auto& x : st.flows) c.W.consume(x.members(), x.alloc);
      gain = t_cand - T;
      if (kTraced && c.trace) {
        Trace shares = Trace::list();
        for (const auto& x : st.flows) shares.push(Trace::integer(x.T));
        step.set("shares", std::move(shares));
      }
    } else {
      const auto members = f.members();
      c.W.gather(members, view_m);
      zeros.assign(members.size(), 0);
      auto fill = lossless_step(c.g, *f.plan, 0, view_m, zeros, st.B);
      if (kTraced && c.trace) {
        step.set("bound", Trace::integer(flow_t_max(c.g, *f.plan, view_m)));
        step.set("fill", Trace::integer(fill.total));
      }
      if (fill.total == 0) {
        if (kTraced && c.trace) {
          step.set("decision", Trace::str("rejected: zero fill"));
          steps.push(std::move(step));
        }
        break;
      }
      f.T = fill.total;
      f.alloc = std::move(fill.alloc);
      c.W.consume(members, f.alloc);
      f.seq = ++c.seq;
      gain = f.T;
      st.flows.push_back(std::move(f));
    }
    ++c.created;
    c.routes_changed = true;
    ++accepted;
    st.gains.push_back(gain);
    if (st.B.has_value()) *st.B -= gain;
    if (kTraced && c.trace) {
      step.set("decision", Trace::str("accepted"));
      step.set("gain", Trace::integer(gain));
      steps.push(std::move(step));
    }
    if (ca.diminishing.has_value() && v.has_value() && accepted >= ca.diminishing->window) {
      Amount recent = 0;
      const auto w = static_cast<std::size_t>(ca.diminishing->window);
      for (std::size_t i = st.gains.size() - w; i < st.gains.size(); ++i) recent += st.gains[i];
      if (below_fraction(recent, ca.diminishing->fraction, *v) && *st.B > 0) {
        st.limit = Termination::LimitDiminishing;
        if (kTraced && c.trace) {
          Trace s = Trace::map();
          s.set("stop", Trace::str("diminishing"));
          s.set("recent", Trace::integer(recent));
          s.set("threshold", Trace::str(fraction_text(ca.diminishing->fraction, *v)));
          steps.push(std::move(s));
        }
        break;
      }
    }
  }
  if (kTraced && st.tr) st.tr->set("stage4", std::move(steps));
}

// ---------------------------------------------------------------- per source

void run_source(Call& c, SourceRun& st) {
  const auto& pe = c.pe;
  if (st.flows.empty()) stage1(c, st);
  if (st.flows.empty()) return;
  if (std::get<CapacityAware>(pe.config.routing).reoptimize) stage2(c, st);
  if (*pe.config.placement == FlowPlacement::EqualBalanced) {
    stage3_eb(c, st);
  } else {
    stage3_seq(c, st);
  }
  stage4(c, st);
}

// A CapacityAware source's outcome (Lossless; Termination's reduction).
Termination local_outcome(const SourceRun& st) {
  Amount after = 0;
  for (const auto& f : st.flows) after += f.T;
  const Amount delivered = after - st.P;
  if (st.offer.has_value() && (*st.offer == 0 || delivered == *st.offer)) return Termination::Complete;
  if (st.limit == Termination::LimitFlows || st.capped) return Termination::LimitFlows;
  if (st.limit == Termination::LimitDiminishing) return Termination::LimitDiminishing;
  if (st.flows.empty()) return Termination::NoRoute;
  return Termination::Exhausted;
}

// The working copy of a live flow: its state and its record's amounts, member-aligned.
WorkFlow working_copy(const PolicyEngine& pe, const LiveFlow& s) {
  const FlowRecord& rec = LedgerAccess::live(pe.fg, s.id);
  WorkFlow f;
  f.id = s.id;
  f.src = s.src;
  f.plan = s.plan;
  f.bundle = s.bundle;
  f.cost = s.cost;
  align(f.plan->members, rec.edges, f.alloc);
  if (pe.config.admission == Admission::Lossy) {
    f.offered = rec.offered.value_or(Counter128 {});
    f.delivered = rec.delivered;
    align(f.plan->members, rec.drops, f.drops);
  } else {
    f.T = rec.delivered;
  }
  return f;
}

Draft draft_of(const PolicyEngine& pe, const WorkFlow& f) {
  Draft d;
  if (!f.is_new()) d.id = f.id;
  auto& r = d.record;
  r.sources = {f.src};
  r.destinations = pe.dests_ok;
  r.admission = pe.config.admission;
  r.bundle = f.bundle;
  r.route_cost = f.cost;
  const auto members = f.members();
  // Exactly sized, written by index (the nonzero members in member order).
  r.edges.resize(static_cast<std::size_t>(std::count_if(f.alloc.begin(), f.alloc.end(), [](Amount x) { return x != 0; })));
  for (std::size_t m = 0, k = 0; m < members.size(); ++m) {
    if (f.alloc[m] != 0) r.edges[k++] = {members[m], f.alloc[m]};
  }
  if (pe.config.admission == Admission::Lossy) {
    r.offered = f.offered;
    r.delivered = f.delivered;
    r.drops.resize(static_cast<std::size_t>(
        std::count_if(f.drops.begin(), f.drops.end(), [](const Counter128& x) { return !x.is_zero(); })));
    for (std::size_t m = 0, k = 0; m < members.size(); ++m) {
      if (!f.drops[m].is_zero()) r.drops[k++] = {members[m], f.drops[m]};
    }
  } else {
    r.delivered = f.T;
  }
  return d;
}

// The epoch state after the commit is built before it (every allocation), then swapped in.
std::vector<std::vector<LiveFlow>> next_cohorts(const std::vector<SourceRun>& runs, std::size_t num_roster) {
  std::vector<std::vector<LiveFlow>> out(num_roster);
  for (const auto& st : runs) {
    auto& c = out[st.ri];
    c.reserve(st.flows.size());
    for (const auto& f : st.flows) c.push_back(LiveFlow {f.id, f.src, f.plan, f.bundle, f.cost, -1});
  }
  return out;
}

std::string key(NodeId v) { return std::to_string(v); }
// ================================================================ calls

// What a call computes. place(v) continues the epoch `cur` (the policy's); evaluate, replace and
// forward start from a fresh epoch (forward keeps the installed program) on the released view.
enum class Op { Place, Evaluate, Replace, Forward };

struct Request {
  Op op {Op::Place};
  std::optional<Amount> v {};          // nullopt: +inf
  const Eligibility* dark {nullptr};   // forward's black holes
  const OutputOptions* outputs {nullptr};
  Trace* trace {nullptr};
  [[nodiscard]] bool full() const noexcept { return op != Op::Place; }
};

// A computed call, published or discarded whole: its result, the drafts of its commit (a place:
// every live record plus the new ones; a full-offer call: only new ones, the live ones retired),
// the epoch after it, and per new draft (in draft order) its flow in that epoch.
struct Staged {
  PlacementResult out {};
  std::vector<Draft> drafts {};
  Epoch next {};
  std::vector<std::pair<std::size_t, std::size_t>> fresh {};
};

// The requested total after the call (a full-offer call: its offer alone) and whether it is
// infinite.
Counter128 requested_after(const Epoch& cur, std::optional<Amount> v) noexcept {
  Counter128 r = cur.requested;
  if (v.has_value()) (void)r.add(*v);  // checked at the call's entry
  return r;
}

// The epoch a full-offer call starts: one empty cohort per roster source (or program entry);
// forward keeps the installed program and CostOnly MaxFlow's support.
Epoch fresh_epoch(const PolicyEngine& pe, bool keep_program) {
  Epoch e;
  e.cohorts.resize(pe.epoch.cohorts.size());
  if (keep_program) {
    e.program = pe.epoch.program;
    e.support = pe.epoch.support;
  }
  return e;
}

// The termination of a forwarding call (Termination's rule after the empty roster).
Termination forwarding_termination(bool active, bool lossy, std::optional<Amount> v, Amount placed) noexcept {
  if (!active) return Termination::NoActiveWeight;
  if (lossy) return v.has_value() ? Termination::Complete : Termination::Exhausted;
  return v.has_value() && placed == *v ? Termination::Complete : Termination::Exhausted;
}

// ---------------------------------------------------------------- CapacityAware

Staged place_flows(const PolicyEngine& pe, const Epoch& cur, Working& W, const Request& rq) {
  const auto& g = *pe.graph;
  const auto v = rq.v;
  Trace* trace = rq.trace;
  Call c {pe, g, W};
#if NETGRAPH_CORE_PRIVATE_HOOKS
  c.cap = max_flows_per_call().load(std::memory_order_relaxed);
#endif
  c.trace = trace;
  const std::size_t k = pe.roster.size();
  std::vector<SourceRun> runs(k);
  if (v.has_value()) cumulative_offers(cur.requested, *v, k, [&](std::size_t i, Amount x) { runs[i].offer = x; });
  Trace sources_tr = Trace::map();
  for (std::size_t i = 0; i < k; ++i) {
    auto& st = runs[i];
    st.ri = i;
    st.src = pe.roster[i];
    st.B = st.offer;
    st.flows.reserve(cur.cohorts[i].size());
    for (const auto& s : cur.cohorts[i]) st.flows.push_back(working_copy(pe, s));
    for (const auto& f : st.flows) {
      st.P += f.T;
      st.T_entry += f.T;
    }
  }
  for (auto& st : runs) {
    Trace* tr = nullptr;
    if (kTraced && trace) tr = &sources_tr.set(key(st.src), Trace::map());
    st.tr = tr;
    run_source(c, st);
    st.tr = nullptr;
  }
  // The commit: retained flows (by source, FlowId order) then new flows in creation order.
  Staged s;
  std::vector<std::pair<std::size_t, std::size_t>> fresh;
  for (std::size_t i = 0; i < k; ++i) {
    for (std::size_t j = 0; j < runs[i].flows.size(); ++j) {
      const auto& f = runs[i].flows[j];
      if (f.is_new()) {
        fresh.emplace_back(i, j);
      } else {
        s.drafts.push_back(draft_of(pe, f));
      }
    }
  }
  std::sort(fresh.begin(), fresh.end(), [&](const auto& a, const auto& b) {
    return runs[a.first].flows[a.second].seq < runs[b.first].flows[b.second].seq;
  });
  for (const auto& [i, j] : fresh) s.drafts.push_back(draft_of(pe, runs[i].flows[j]));
  s.fresh = std::move(fresh);
  const Counter128 requested = requested_after(cur, v);
  Amount placed = 0;
  std::vector<Termination> outcomes;
  outcomes.reserve(k);
  Termination term = Termination::Complete;
  for (const auto& st : runs) {
    Amount after = 0;
    for (const auto& f : st.flows) after += f.T;
    placed += after - st.P;
    outcomes.push_back(local_outcome(st));
    if (rank(outcomes.back()) > rank(term)) term = outcomes.back();
  }
  PlacementResult& out = s.out;
  out.placed = placed;
  if (v.has_value()) out.remaining = *v - placed;
  out.delivery_delta = placed;
  out.termination = term;
  out.routes_changed = c.routes_changed;
  out.weights_changed = c.weights_changed;
  fill_state(pe, draft_records(s.drafts), requested, cur.infinite || !v.has_value(), Counter128 {}, out);
  if (kTraced && trace) {
    Trace offers_tr = Trace::map();
    Trace outcomes_tr = Trace::map();
    for (std::size_t i = 0; i < k; ++i) {
      offers_tr.set(key(runs[i].src), amount_or_inf(runs[i].offer));
      outcomes_tr.set(key(runs[i].src), Trace::str(outcome_name(outcomes[i])));
    }
    trace->set("offers", std::move(offers_tr));
    trace->set("local_outcomes", std::move(outcomes_tr));
    trace->set("created", Trace::integer(c.created));
    trace->set("sources", std::move(sources_tr));
    std::vector<Amount> wf;
    W.dense(wf);
    trace->set("W_final", Trace::ints(wf));
  }
  s.next.cohorts = next_cohorts(runs, k);
  s.next.requested = requested;
  s.next.prepared = true;
  s.next.infinite = cur.infinite || !v.has_value();
  return s;
}

// ---------------------------------------------------------------- the controllers' paths

// The flow cut into paths (4.1), source by source in ascending NodeId: while a source has supply
// left, take the smallest EdgeId sequence with positive remaining flow, stopping at the first
// destination; assign min(supply left, smallest remaining edge flow on that walk), subtract it
// from the walk, record the walk as one path. Each node's out-edge cursor only advances:
// remaining flow never increases. `edges` (ascending) are the flows' edges (empty: every edge, the
// flows per EdgeId); `rem` is aligned with them. Per source also its share (the edges and amounts
// its paths use, ascending).
struct Decomposition {
  std::vector<ControllerPath> paths {};
  std::vector<std::vector<std::pair<EdgeId, Amount>>> shares {};
  std::vector<Amount> supplies {};
};

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

Decomposition decompose(const Graph& g, std::span<const EdgeId> edges, std::vector<Amount> rem,
                        std::span<const NodeId> roster, std::span<const NodeId> dests) {
  const bool dense = edges.empty() && rem.size() == static_cast<std::size_t>(g.num_edges());
  const auto edge_src = g.edge_src();
  const auto edge_dst = g.edge_dst();
  // Position p of the flows: EdgeId p (dense) or edges[p]. Out-edges have consecutive EdgeIds, so
  // a node's positions are one range.
  auto edge_at = [&](std::size_t p) { return dense ? static_cast<EdgeId>(p) : edges[p]; };
  thread_local DecomposeScratch sc;
  const auto n = static_cast<std::size_t>(g.num_nodes());
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
  } reset {sc, dests, roster};
  for (NodeId d : dests) sc.is_dst[static_cast<std::size_t>(d)] = 1;
  for (std::size_t i = 0; i < roster.size(); ++i) sc.slot[static_cast<std::size_t>(roster[i])] = static_cast<std::int32_t>(i);
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
  out.shares.resize(roster.size());
  out.supplies.resize(roster.size());
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
    const auto u = static_cast<std::size_t>(edge_src[e]);
    if (sc.cur[u] < 0) {
      sc.touched.push_back(edge_src[e]);
      sc.cur[u] = static_cast<std::int64_t>(p);
    }
    sc.end[u] = p + 1;
    const auto from = sc.slot[u];
    const auto to = sc.slot[static_cast<std::size_t>(edge_dst[e])];
    if (from >= 0) out.supplies[static_cast<std::size_t>(from)] += rem[p];
    if (to >= 0) out.supplies[static_cast<std::size_t>(to)] -= rem[p];
  }
  // A one-source flow is all that source's share: the walks cover the flow (checked below), so the
  // sum of its paths is the flow's positive positions, ascending (the per-path sum's degenerate
  // form).
  const bool one = roster.size() == 1;
  // Every path zeroes a position or ends a source's supply.
  out.paths.reserve(positive + roster.size());
  if (one) {
    out.shares[0].reserve(positive);
    for (std::size_t p = 0; p < rem.size(); ++p) {
      if (rem[p] != 0) out.shares[0].emplace_back(edge_at(p), rem[p]);
    }
  }
  auto& walk = sc.walk;
  auto& acc = sc.acc;
  for (std::size_t i = 0; i < roster.size(); ++i) {
    const NodeId s = roster[i];
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
          throw std::logic_error("Global decomposition: flow conservation broken (no continuation)");
        }
        walk.push_back(static_cast<std::size_t>(p));
        u = edge_dst[static_cast<std::size_t>(edge_at(static_cast<std::size_t>(p)))];
        if (walk.size() > n) throw std::logic_error("Global decomposition: a cycle with positive flow");
      }
      Amount amount = left;
      for (std::size_t p : walk) amount = std::min(amount, rem[p]);
      ControllerPath path;
      path.slot = static_cast<std::int32_t>(i);
      path.amount = amount;
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
    if (x != 0) throw std::logic_error("Global decomposition: the walks do not cover the flow");
  }
  return out;
}

// ---------------------------------------------------------------- CostOnly preparation

// A CostOnly demand's routes for the epoch, prepared at its first positive call (or a replace)
// on the working residual at that call's entry: per roster source its route (the record's bundle
// and route cost), and its program (EqualBalanced, Weighted: one split route per source in one
// plan rooted at the roster) or the routes' union (MaxFlow, the controller's support).
struct Prepared {
  std::vector<LiveFlow> flows {};  // per roster source, new (id -1)
  std::shared_ptr<const ForwardingProgram> program {};
  std::vector<EdgeId> edges {};    // ascending: the union of the routes
};

Prepared prepare_cost_only(Call& c) {
  const auto& pe = c.pe;
  const auto& g = c.g;
  const std::size_t k = pe.roster.size();
  const FlowPlacement placement = *pe.config.placement;
  Prepared out;
  out.flows.resize(k);
  std::vector<std::vector<EdgeId>> routes(k);
  std::shared_ptr<const KernelPlan> single;  // one source's plan from its one-pass creation
  for (std::size_t i = 0; i < k; ++i) {
    const NodeId s = pe.roster[i];
    auto& f = out.flows[i];
    f.src = s;
    if (auto made = creation_products(c, s)) {
      f.bundle = std::move(made->bundle);
      f.cost = made->cost;
      if (k == 1 && placement != FlowPlacement::MaxFlow) {
        single = std::make_shared<const KernelPlan>(std::move(made->plan));
      } else {
        routes[i] = std::move(made->plan.route);
      }
    } else {
      Route r = cost_only_route(c, s);
      f.bundle = bundle_of(g, r.edges, pe.support);
      f.cost = r.cost != kUnreachable ? r.cost : route_cost(g, *f.bundle, s, pe.dests_ok);
      routes[i] = std::move(r.edges);
    }
  }
  // The program: one entry per roster source with its record bundle (its route, the replay basis)
  // and route cost.
  auto program_of = [&](std::shared_ptr<const KernelPlan> plan) {
    std::vector<ProgramEntry> entries(k);
    for (std::size_t i = 0; i < k; ++i) {
      entries[i].bundle = out.flows[i].bundle;
      entries[i].cost = out.flows[i].cost;
    }
    return cost_only_program(g, pe.sources, pe.destinations, pe.roster, std::move(plan), std::move(entries));
  };
  // One split-routing source: its one-pass creation plan is the program (its route is the union).
  if (single) {
    out.program = program_of(std::move(single));
    return out;
  }
  for (const auto& r : routes) out.edges.insert(out.edges.end(), r.begin(), r.end());
  std::sort(out.edges.begin(), out.edges.end());
  out.edges.erase(std::unique(out.edges.begin(), out.edges.end()), out.edges.end());
  // One continuation per router: every source's route leaves each of its routers by all the
  // union's edges there. The routes are read on one view under one mask, so a router's
  // continuation depends on the router alone (a single path's suffix is the canonical path of
  // its first node; the rejected multipath = false with PreferHigherPriority is the exception).
  // One route is its own union: nothing to check.
  const auto row = g.row_offsets();
  const auto edge_src = g.edge_src();
  for (std::size_t i = 0; k > 1 && i < k; ++i) {
    const auto& r = routes[i];
    for (std::size_t a = 0; a < r.size();) {
      const auto u = static_cast<std::size_t>(edge_src[static_cast<std::size_t>(r[a])]);
      std::size_t b = a;
      while (b < r.size() && static_cast<std::size_t>(edge_src[static_cast<std::size_t>(r[b])]) == u) ++b;
      const auto lo = std::lower_bound(out.edges.begin(), out.edges.end(), row[u]);
      const auto hi = std::lower_bound(lo, out.edges.end(), row[u + 1]);
      if (static_cast<std::size_t>(hi - lo) != b - a) {
        throw std::logic_error("CostOnly: two roster sources continue differently at node " + std::to_string(u));
      }
      a = b;
    }
  }
  if (placement == FlowPlacement::MaxFlow) return out;
  const auto support = PredDAGBuilder::ascending(*pe.support);
  const auto dag = std::equal(out.edges.begin(), out.edges.end(), support.begin(), support.end())
                       ? pe.support
                       : std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, out.edges));
  // Weighted weights: the static source, or the free capacity read now (ResidualSnapshot, frozen
  // for the epoch).
  std::vector<Amount> creation;
  WeightSource weights = cost_only_weights(pe);
  if (placement == FlowPlacement::Weighted && weights.basis == WeightBasis::ResidualSnapshot) {
    creation.resize(out.edges.size());
    for (std::size_t i = 0; i < out.edges.size(); ++i) creation[i] = c.W.at(out.edges[i]);
    weights.values = creation;
  }
  out.program = program_of(std::make_shared<const KernelPlan>(
      make_plan(g, placement, pe.config.admission, weights, pe.roster, pe.dests_ok, *dag)));
  return out;
}

Trace routes_trace(const Prepared& p) {
  Trace t = Trace::map();
  for (const auto& f : p.flows) t.set(key(f.src), route_trace(PredDAGBuilder::ascending(*f.bundle)));
  return t;
}

// A record of one roster source of CostOnly or Global (its route, its route cost).
FlowRecord source_record(const PolicyEngine& pe, const LiveFlow& f) {
  FlowRecord r;
  r.sources = {f.src};
  r.destinations = pe.dests_ok;
  r.admission = pe.config.admission;
  r.bundle = f.bundle;
  r.route_cost = f.cost;
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

// One record per roster source (CostOnly, Global): the drafts in roster order, kept ids when
// the epoch already has them, else new (with their place in the next epoch).
void source_drafts(const PolicyEngine& pe, const std::vector<LiveFlow>& flows, Staged& s) {
  const std::size_t k = flows.size();
  s.drafts.resize(k);
  s.next.cohorts.assign(k, {});
  for (std::size_t i = 0; i < k; ++i) {
    if (flows[i].id >= 0) {
      s.drafts[i].id = flows[i].id;
    } else {
      s.fresh.emplace_back(i, 0);
    }
    s.drafts[i].record = source_record(pe, flows[i]);
    s.next.cohorts[i].push_back(flows[i]);
  }
}

// The programs of a call and the flags they give (10.4).
void set_flags(PlacementResult& out, const ForwardingProgram* before, const ForwardingProgram* after) noexcept {
  const auto change = compare_programs(before, after);
  out.routes_changed = change.routes;
  out.weights_changed = change.weights;
}

// The by_tunnel rows of a controller's program: every path up, assigned = delivered = its amount.
std::vector<TunnelFlow> controller_rows(const ForwardingProgram& program) {
  std::vector<TunnelFlow> rows;
  rows.reserve(program.tunnels().size());
  for (const auto& t : program.tunnels()) rows.push_back(TunnelFlow {t.index, t.head, TunnelStatus::Up, c128(t.weight), t.weight, {}});
  return rows;
}

// ---------------------------------------------------------------- CostOnly EqualBalanced and Weighted

Staged place_pooled(const PolicyEngine& pe, const Epoch& cur, Working& W, const Request& rq,
                    const ForwardingProgram* before_program) {
  const auto& g = *pe.graph;
  const auto v = rq.v;
  Trace* trace = rq.trace;
  Call c {pe, g, W};
  c.trace = trace;
  const std::size_t k = pe.roster.size();
  const bool lossy = pe.config.admission == Admission::Lossy;
  const bool created = cur.cohorts.front().empty();
  Prepared prep;
  std::vector<LiveFlow> flows(k);
  if (created) {
    prep = prepare_cost_only(c);
    flows = prep.flows;
  } else {
    for (std::size_t i = 0; i < k; ++i) flows[i] = cur.cohorts[i].front();
  }
  const auto program = created ? prep.program : cur.program;
  const KernelPlan& plan = *ProgramAccess::data(*program).contexts.front().plan;
  const auto members = std::span<const EdgeId>(plan.members);
  std::vector<const FlowRecord*> before(k, nullptr);
  if (!created) {
    for (std::size_t i = 0; i < k; ++i) before[i] = &LedgerAccess::live(pe.fg, flows[i].id);
  }
  Staged s;
  source_drafts(pe, flows, s);
  auto& drafts = s.drafts;
  // A one-source program is its own attribution.
  Attribution attr;
  Attribution* per_source = k > 1 ? &attr : nullptr;
  std::vector<Amount> inject(k, 0);
  Forwarding fw;
  Amount placed = 0;
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
      align(members, r->edges, part);
      for (std::size_t m = 0; m < members.size(); ++m) own[m] += part[m];
    }
    std::vector<Amount> view;
    W.gather(members, view);
    for (std::size_t m = 0; m < members.size(); ++m) view[m] += own[m];
    auto step = lossless_step(g, plan, admitted, view, own, v);
    placed = step.total - admitted;
    if (per_source == nullptr) {
      // One source holds the whole admitted allocation (the attribution's degenerate form): it
      // delivers all of it, with no second forwarding.
      auto& r = drafts.front().record;
      r.edges = positive(members, step.alloc);
      r.delivered = step.total;
    } else {
      // Several sources: forward the admitted total again to attribute it.
      round_robin(step.total, inject);
      fw = forward(plan, Injection {inject, false}, {}, per_source);
      if (fw.carried != step.alloc) {
        throw std::logic_error("CostOnly: the attributed forwarding differs from the admitted allocation");
      }
      for (std::size_t i = 0; i < k; ++i) {
        auto& r = drafts[i].record;
        for (const auto& x : attr.entries[i]) {
          if (x.carried != 0) r.edges.emplace_back(members[static_cast<std::size_t>(x.member)], x.carried);
        }
        r.delivered = attr.delivered[i];
      }
    }
    if (kTraced && trace) step_tr = Trace::ints(std::vector<Amount> {admitted, step.total});
  } else {
    // The new round-robin injections of the offer (or the first-hop experiment) on W.
    Injection in;
    if (v.has_value()) {
      cumulative_offers(cur.offered, *v, k, [&](std::size_t i, Amount x) { inject[i] = x; });
      in.per_root = inject;
    } else {
      in.first_hop = true;
    }
    std::vector<Amount> w;
    W.gather(members, w);
    fw = forward(plan, in, w, per_source);
    placed = fw.delivered;
    for (std::size_t i = 0; i < k; ++i) {
      auto& r = drafts[i].record;
      std::vector<std::pair<EdgeId, Amount>> carried;
      std::vector<std::pair<EdgeId, Amount>> dropped;
      Amount offer = fw.offered;
      Amount delivered = fw.delivered;
      if (per_source == nullptr) {
        carried = positive(members, fw.carried);
        dropped = positive(members, fw.drops);
      } else {
        for (const auto& x : attr.entries[i]) {
          const EdgeId e = members[static_cast<std::size_t>(x.member)];
          if (x.carried != 0) carried.emplace_back(e, x.carried);
          if (x.dropped != 0) dropped.emplace_back(e, x.dropped);
        }
        offer = attr.offered[i];
        delivered = attr.delivered[i];
      }
      // The counters' 2^127 checks: a drop counter first, then the offered ones.
      const FlowRecord* b = before[i];
      r.drops = plus(b != nullptr ? b->drops : std::vector<std::pair<EdgeId, Counter128>> {}, dropped);
      r.offered = b != nullptr ? b->offered.value_or(Counter128 {}) : Counter128 {};
      accumulate(*r.offered, c128(offer), "offered");
      r.delivered = (b != nullptr ? b->delivered : 0) + delivered;
      r.edges = b != nullptr ? plus(b->edges, carried) : std::move(carried);
    }
    accumulate(offered, c128(fw.offered), "offered");
    if (kTraced && trace) {
      step_tr = Trace::map();
      for (std::size_t i = 0; i < k; ++i) {
        step_tr.set(key(pe.roster[i]), v.has_value() ? Trace::integer(inject[i]) : Trace::str("inf"));
      }
    }
  }
  const Counter128 requested = requested_after(cur, v);
  PlacementResult& out = s.out;
  out.placed = placed;
  if (v.has_value()) out.remaining = *v - placed;
  out.delivery_delta = placed;
  if (lossy) out.call_drops = positive(members, fw.drops);
  out.call_source_loss = fw.source_loss;
  if (members.empty()) {
    out.termination = Termination::NoActiveWeight;
  } else if (lossy) {
    out.termination = v.has_value() ? Termination::Complete : Termination::Exhausted;
  } else {
    out.termination = v.has_value() && placed == *v ? Termination::Complete : Termination::Exhausted;
  }
  if (created) set_flags(out, before_program, program.get());
  out.program = program;
  fill_state(pe, draft_records(drafts), requested, cur.infinite || !v.has_value(), offered, out);
  if (kTraced && trace) {
    trace->set("prepare", created ? routes_trace(prep) : Trace::null());
    trace->set(lossy ? "offers" : "T", std::move(step_tr));
  }
  s.next.program = program;
  s.next.requested = requested;
  s.next.offered = offered;
  s.next.prepared = true;
  s.next.infinite = cur.infinite || !v.has_value();
  return s;
}

// ---------------------------------------------------------------- the controllers

// One solve of a controller from zero on the view V with value cap `cap` (nullopt: uncapped);
// `caps` (the equal origin) caps source i's arc at caps[i]. Unrestricted Global solves over every
// eligible edge, its V and flows per EdgeId; Global with a support and CostOnly MaxFlow over
// `support`, their V and flows aligned with it.
struct Solve {
  Amount placed {0};
  std::vector<Amount> flows {};
  FlowCoreResult payload {};  // unrestricted Global: the flow core's payload (the curve, the cut, the trace's values)
};

Solve solve(const PolicyEngine& pe, std::span<const EdgeId> support, std::span<const Amount> V,
            std::optional<Amount> cap, std::span<const Amount> caps) {
  const auto& g = *pe.graph;
  Solve out;
  if (dense_global(pe)) {
    out.payload = place_global_units(g, V, pe.roster, pe.dests_ok, pe.elig, cap, caps.empty() ? pe.capture : nullptr, caps);
    out.placed = out.payload.placed;
    out.flows = out.payload.edge_flows;
    return out;
  }
  auto r = solve_subset(g, support, V, pe.roster, pe.dests_ok, cap, caps);
  out.placed = r.placed;
  out.flows = std::move(r.flows);
  return out;
}

// The equal origin: the largest total T in [lo, hi] whose round-robin supplies RR(T, H, i) are all
// carried (a solve capped at T with source i's arc capped at its share reaches T), given that lo's
// are. Feasible totals form a prefix (each share never decreases in T, and a smaller supply
// vector stays feasible), so integer bisection is exact; T's solve from zero is returned, so the
// result does not depend on the probes.
Solve solve_equal(const PolicyEngine& pe, std::span<const EdgeId> support, std::span<const Amount> V, Amount lo,
                  Amount hi, Trace* probes) {
  std::vector<Amount> caps(pe.roster.size());
  auto attempt = [&](Amount t) -> std::optional<Solve> {
    std::fill(caps.begin(), caps.end(), 0);
    round_robin(t, caps);
    auto s = solve(pe, support, V, t, caps);
    if (kTraced && probes != nullptr) probes->push(Trace::integer(t));
    if (s.placed != t) return std::nullopt;
    return s;
  };
  if (auto s = attempt(hi)) return std::move(*s);
  std::optional<Solve> best;
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
  if (!best) throw std::logic_error("equal origin: the admitted total is no longer carried");
  return std::move(*best);
}

// The upper end of the equal search: the view summed over the edges leaving the roster sources
// (unrestricted Global: the eligible ones, V per EdgeId; otherwise the support's, V aligned with
// it); the supplies leave through them.
Amount roster_out(const PolicyEngine& pe, std::span<const EdgeId> support, std::span<const Amount> V) {
  const auto row = pe.graph->row_offsets();
  const bool dense = dense_global(pe);
  Amount sum = 0;
  for (NodeId s : pe.roster) {
    const auto si = static_cast<std::size_t>(s);
    if (dense) {
      for (auto e = row[si]; e < row[si + 1]; ++e) {
        if (pe.elig == nullptr || pe.elig->edge_ok(e)) sum += V[static_cast<std::size_t>(e)];
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

// The program basis of a controller: unrestricted Global the topology under its mask; Global with
// a support that support; CostOnly MaxFlow each source's route.
ControllerBasis controller_basis(const PolicyEngine& pe, const std::vector<LiveFlow>& flows) {
  ControllerBasis b;
  if (dense_global(pe)) {
    b.topology = true;
    b.excluded_nodes = pe.excluded_nodes;
    b.excluded_edges = pe.excluded_edges;
  } else if (is<Global>(pe)) {
    b.common = pe.global_support;
  } else {
    b.per_source.reserve(flows.size());
    for (const auto& f : flows) b.per_source.push_back(f.bundle);
  }
  return b;
}

Staged place_controller(const PolicyEngine& pe, const Epoch& cur, Working& W, const Request& rq,
                        const ForwardingProgram* before_program) {
  const auto& g = *pe.graph;
  const auto v = rq.v;
  Trace* trace = rq.trace;
  const bool dense = dense_global(pe);
  const bool global = is<Global>(pe);
  const bool equal = pe.origin == Origin::Equal;
  const std::size_t k = pe.roster.size();
  const bool created = cur.cohorts.front().empty();
  // CostOnly MaxFlow: the support of the epoch (its routes, prepared at its first call).
  Prepared prep;
  if (!global && created) {
    Call c {pe, g, W};
    prep = prepare_cost_only(c);
  }
  const std::span<const EdgeId> support = dense      ? std::span<const EdgeId> {}
                                          : global   ? std::span<const EdgeId>(*pe.global_support)
                                          : created  ? std::span<const EdgeId>(prep.edges)
                                                     : std::span<const EdgeId>(cur.support);
  std::vector<LiveFlow> flows(k);
  for (std::size_t i = 0; i < k; ++i) {
    if (!created) {
      flows[i] = cur.cohorts[i].front();
    } else if (!global) {
      flows[i] = prep.flows[i];
    } else {
      flows[i].src = pe.roster[i];
    }
  }
  // V = W + own: the working residual plus the policy's own allocation (its live records, or none
  // in a full-offer call, whose view already holds them released), per EdgeId (unrestricted
  // Global) or aligned with the support (a record's edges lie in it, both ascending).
  std::vector<Amount> V;
  if (dense) {
    W.dense(V);
  } else {
    W.gather(support, V);
  }
  Amount P = 0;
  for (const auto& f : flows) {
    if (f.id < 0) continue;
    const FlowRecord& rec = LedgerAccess::live(pe.fg, f.id);
    P += rec.delivered;
    std::size_t p = 0;
    for (const auto& [e, x] : rec.edges) {
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
    if (equal || !dense) throw std::logic_error("a deletion capture of a controller other than optimized Global");
    pe.capture->graph = pe.graph;
  }
  Trace probes = Trace::list();
  Solve sol;
  if (equal) {
    const Amount top = roster_out(pe, support, V);
    sol = solve_equal(pe, support, V, P, cap.has_value() ? std::min(*cap, top) : top, &probes);
  } else {
    sol = solve(pe, support, V, cap, {});
  }
  if (sol.placed < P) throw std::logic_error("controller: a call reduced placed_demand()");
  // The program: each source's supply cut into paths (4.1); each source's share is the sum of its
  // paths.
  // The solved flows move into the decomposition (a trace reads them afterwards: a copy).
  const bool traced = kTraced && trace != nullptr;
  auto dec = decompose(g, support, traced ? sol.flows : std::move(sol.flows), pe.roster, pe.dests_ok);
  Staged s;
  // A supplied support, which a share that uses every one of its edges reuses as its bundle (a
  // PredDAG is fixed by its graph and edge set: the same value).
  const auto* gl = std::get_if<Global>(&pe.config.routing);
  const std::shared_ptr<const PredDAG> supplied = gl != nullptr ? gl->support : nullptr;
  for (std::size_t i = 0; i < k; ++i) {
    const NodeId src = pe.roster[i];
    auto& f = flows[i];
    if (global) {
      // Global's record bundle is the edges its share uses.
      const auto& share = dec.shares[i];
      auto same = [&](std::span<const EdgeId> edges) {
        return std::equal(share.begin(), share.end(), edges.begin(), edges.end(),
                          [](const auto& ex, EdgeId e) { return ex.first == e; });
      };
      if (f.bundle && same(PredDAGBuilder::ascending(*f.bundle))) continue;
      if (supplied && same(PredDAGBuilder::ascending(*supplied))) {
        f.bundle = supplied;
      } else {
        std::vector<EdgeId> route;
        route.reserve(share.size());
        for (const auto& [e, x] : share) route.push_back(e);
        f.bundle = std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, route));
      }
      f.cost = share.empty() ? kUnreachable : route_cost(g, *f.bundle, src, pe.dests_ok);
    }
  }
  source_drafts(pe, flows, s);
  for (std::size_t i = 0; i < k; ++i) {
    s.drafts[i].record.delivered = dec.supplies[i];
    s.drafts[i].record.edges = std::move(dec.shares[i]);
  }
  const auto& supplies = dec.supplies;
  auto program = controller_program(g, pe.sources, pe.destinations, pe.roster, std::move(dec.paths),
                                    equal ? std::span<const Amount> {} : std::span<const Amount>(supplies),
                                    controller_basis(pe, flows));
  const Counter128 requested = requested_after(cur, v);
  PlacementResult& out = s.out;
  out.placed = sol.placed - P;
  if (v.has_value()) out.remaining = *v - out.placed;
  out.delivery_delta = out.placed;
  out.termination = v.has_value() && out.placed == *v ? Termination::Complete : Termination::Exhausted;
  set_flags(out, before_program, program.get());
  out.program = program;
  out.by_tunnel = controller_rows(*program);
  const auto after = draft_records(s.drafts);
  fill_state(pe, after, requested, cur.infinite || !v.has_value(), Counter128 {}, out);
  out.total_cost = total_cost_of(g, after);
  const bool curve = dense && !equal;
  if (curve) out.cost_distribution = sol.payload.cost_distribution;
  if (curve && rq.outputs != nullptr && rq.outputs->min_cut) out.min_cut = sol.payload.min_cut;
  if (kTraced && trace) {
    if (!global) trace->set("prepare", created ? routes_trace(prep) : Trace::null());
    if (equal) trace->set("probes", std::move(probes));
    if (dense) trace->set("V", Trace::ints(V));
    trace->set("cap", cap.has_value() ? Trace::integer(*cap) : Trace::null());
    Trace fc = Trace::map();
    fc.set("placed", Trace::integer(sol.placed));
    if (dense) {
      const auto& r = sol.payload;
      fc.set("cost_distribution", Trace::pairs(r.cost_distribution));
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
      for (std::size_t p = 0; p < support.size(); ++p) flows_dense[static_cast<std::size_t>(support[p])] = sol.flows[p];
      fc.set("edge_flows", Trace::ints(flows_dense));
    }
    Trace supply = Trace::map();
    Trace shares_tr = Trace::map();
    for (std::size_t i = 0; i < k; ++i) {
      supply.set(key(pe.roster[i]), Trace::integer(supplies[i]));
      shares_tr.set(key(pe.roster[i]), Trace::pairs(s.drafts[i].record.edges));
    }
    fc.set("supply", std::move(supply));
    trace->set("flow_core", std::move(fc));
    trace->set("shares", std::move(shares_tr));
  }
  if (curve) s.next.curve = *out.cost_distribution;  // the cumulative curve a later zero call reports
  if (!global) s.next.support = created ? prep.edges : cur.support;
  s.next.program = program;
  s.next.by_tunnel = out.by_tunnel;
  s.next.requested = requested;
  s.next.prepared = true;
  s.next.infinite = cur.infinite || !v.has_value();
  return s;
}

// ---------------------------------------------------------------- forwarding a program

// The by_tunnel rows of a forwarding: per tunnel its status in the view and its amounts
// (`assigned`, `delivered`, `dropped` per program entry).
std::vector<TunnelFlow> tunnel_rows(const ForwardingView& view, std::span<const Counter128> assigned,
                                    std::span<const Amount> delivered, std::span<const Counter128> dropped) {
  const ProgramData& p = *view.data;
  std::vector<TunnelFlow> rows;
  rows.reserve(p.tunnels.size());
  for (std::size_t t = 0; t < p.tunnels.size(); ++t) {
    const auto e = static_cast<std::size_t>(p.tunnel_entry[t]);
    rows.push_back(TunnelFlow {p.tunnels[t].index, p.tunnels[t].head, view.status(e), assigned[e], delivered[e], dropped[e]});
  }
  return rows;
}

// The free capacity of the forwarding's edges in the call's working residual.
FreeFn free_of(const Working& W) {
  return [&W](std::span<const EdgeId> edges, std::span<Amount> out) {
    for (std::size_t i = 0; i < edges.size(); ++i) out[i] = W.at(edges[i]);
  };
}

// Pinned: forwarding its program (place, evaluate, replace), one record per program entry.
Staged place_pinned(const PolicyEngine& pe, const Epoch& cur, Working& W, const Request& rq,
                    const ForwardingProgram* before_program) {
  const auto& g = *pe.graph;
  const auto v = rq.v;
  const ForwardingView& view = pe.view;
  const ProgramData& p = *view.data;
  const std::size_t n = p.entries.size();
  const std::size_t k = view.num_heads();
  const bool lossy = pe.config.admission == Admission::Lossy;
  const bool created = !cur.prepared;
  // The entries' live records (none before the first call).
  std::vector<const FlowRecord*> before;
  if (!created) {
    before.resize(n);
    for (std::size_t e = 0; e < n; ++e) before[e] = &LedgerAccess::live(pe.fg, cur.cohorts[e].front().id);
  }
  ProgramFlow flow;
  Amount placed = 0;
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
    flow = detail::lossless(g, view, free_of(W), admitted, hi);
    placed = flow.total - admitted;
  } else {
    std::vector<Amount> inject(k, 0);
    if (rq.full()) {
      head_shares(view, *v, inject);
    } else {
      head_increments(view, cur.offered, *v, inject);
    }
    flow = detail::lossy(g, view, free_of(W), inject);
    placed = flow.delivered;
    // The whole offer is sent: what no head receives (every remaining captured weight 0) is lost at
    // the sources, outside every row.
    flow.source_loss += *v - flow.offered;
    accumulate(offered, c128(*v), "offered");
  }
  Staged s;
  s.drafts.resize(n);
  s.next.cohorts.assign(n, {});
  // The by_tunnel amounts per entry, kept only when the program has tunnels (only they have rows).
  const bool rows = !p.tunnels.empty();
  std::vector<Counter128> assigned(rows ? n : 0);
  std::vector<Amount> delivered(rows ? n : 0, 0);
  std::vector<Counter128> dropped(rows ? n : 0);
  std::vector<std::pair<EdgeId, Amount>> call_drops;
  for (std::size_t e = 0; e < n; ++e) {
    const auto& entry = p.entries[e];
    auto& f = flow.entries[e];
    LiveFlow live;
    live.src = entry.head;
    live.bundle = view.bundle(e);
    live.cost = view.cost(e);
    live.entry = entry.index;
    // A split route's plan (its members and weights as forwarded); a tunnel has none.
    if (entry.context >= 0) live.plan = view.plan(static_cast<std::size_t>(entry.context));
    if (!created) {
      live.id = cur.cohorts[e].front().id;
      s.drafts[e].id = live.id;
    } else {
      s.fresh.emplace_back(e, 0);
    }
    s.next.cohorts[e].push_back(live);
    FlowRecord& r = s.drafts[e].record;
    r.sources = {entry.head};
    r.destinations = pe.dests_ok;
    r.admission = pe.config.admission;
    r.bundle = live.bundle;
    r.route_cost = live.cost;
    r.entry = entry.index;
    const FlowRecord* b = created ? nullptr : before[e];
    if (!lossy) {
      r.edges = std::move(f.carried);
      r.delivered = f.delivered;
      if (rows) {
        assigned[e] = c128(f.delivered);
        delivered[e] = f.delivered;
      }
      continue;
    }
    r.drops = plus(b != nullptr ? b->drops : std::vector<std::pair<EdgeId, Counter128>> {}, f.dropped);
    r.offered = b != nullptr ? b->offered.value_or(Counter128 {}) : Counter128 {};
    accumulate(*r.offered, c128(f.assigned), "offered");
    r.delivered = (b != nullptr ? b->delivered : 0) + f.delivered;
    r.edges = b != nullptr ? plus(b->edges, f.carried) : std::move(f.carried);
    call_drops.insert(call_drops.end(), f.dropped.begin(), f.dropped.end());
    Counter128 drops;
    for (const auto& [x, d] : r.drops) accumulate(drops, d, "drops");
    if (rows) {
      assigned[e] = *r.offered;
      delivered[e] = r.delivered;
      dropped[e] = drops;
    }
  }
  // The shares the heads lost themselves (Lossy): per view head, cumulative.
  std::vector<std::pair<NodeId, Counter128>> lost;
  if (lossy) {
    s.next.head_lost = cur.head_lost.empty() ? std::vector<Counter128>(k) : cur.head_lost;
    for (std::size_t h = 0; h < k; ++h) {
      accumulate(s.next.head_lost[h], c128(flow.head_lost[h]), "offered");
      if (!s.next.head_lost[h].is_zero()) lost.emplace_back(view.head(h), s.next.head_lost[h]);
    }
  }
  const Counter128 requested = requested_after(cur, v);
  PlacementResult& out = s.out;
  out.placed = placed;
  if (v.has_value()) out.remaining = *v - placed;
  out.delivery_delta = placed;
  if (lossy) out.call_drops = merged(std::move(call_drops));
  out.call_source_loss = flow.source_loss;
  // No head reaches a destination (every entry failed): NO_ROUTE, the records kept per entry.
  out.termination = k == 0 ? Termination::NoRoute : forwarding_termination(pe.active, lossy, v, placed);
  set_flags(out, before_program, view.program.get());
  out.program = view.program;
  out.by_tunnel = tunnel_rows(view, assigned, delivered, dropped);
  fill_state(pe, draft_records(s.drafts), requested, cur.infinite || !v.has_value(), offered, out, lost);
  s.next.program = view.program;
  s.next.by_tunnel = out.by_tunnel;
  s.next.requested = requested;
  s.next.offered = offered;
  s.next.prepared = true;
  s.next.infinite = cur.infinite || !v.has_value();
  return s;
}

// forward(): a computing model's installed program forwarded with no route computation, on the
// released view, records per roster source (CostOnly and Global).
Staged forward_installed(const PolicyEngine& pe, Working& W, const Request& rq) {
  const auto& g = *pe.graph;
  const auto v = rq.v;
  const auto& program = pe.epoch.program;
  const ForwardingView view = make_view(g, program, nullptr, pe.dests_ok, rq.dark);
  const ProgramData& p = *view.data;
  const std::size_t k = pe.roster.size();
  const std::size_t n = p.entries.size();
  const bool lossy = pe.config.admission == Admission::Lossy;
  ProgramFlow flow;
  if (!lossy) {
    flow = detail::lossless(g, view, free_of(W), 0, v);
  } else {
    std::vector<Amount> inject(k, 0);
    head_shares(view, *v, inject);
    flow = detail::lossy(g, view, free_of(W), inject);
  }
  Staged s;
  s.next = fresh_epoch(pe, true);
  std::vector<LiveFlow> flows(k);
  // Per roster source: its entries' carried, drops and deliveries summed (the view's heads are the
  // roster: nothing is replayed).
  std::vector<std::vector<std::pair<EdgeId, Amount>>> carried(k);
  std::vector<std::vector<std::pair<EdgeId, Amount>>> drops(k);
  std::vector<Amount> got(k, 0);
  for (std::size_t e = 0; e < n; ++e) {
    const auto slot = static_cast<std::size_t>(p.entries[e].slot);
    const auto& f = flow.entries[e];
    carried[slot].insert(carried[slot].end(), f.carried.begin(), f.carried.end());
    drops[slot].insert(drops[slot].end(), f.dropped.begin(), f.dropped.end());
    got[slot] += f.delivered;
  }
  for (std::size_t i = 0; i < k; ++i) {
    flows[i] = pe.epoch.cohorts[i].front();
    flows[i].id = -1;
    carried[i] = merged(std::move(carried[i]));
    drops[i] = merged(std::move(drops[i]));
    if (is<Global>(pe)) {
      // Global's record bundle is the edges its share uses.
      std::vector<EdgeId> route;
      route.reserve(carried[i].size());
      for (const auto& [e, x] : carried[i]) route.push_back(e);
      flows[i].bundle = std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, route));
      flows[i].cost = route.empty() ? kUnreachable : route_cost(g, *flows[i].bundle, pe.roster[i], pe.dests_ok);
    }
  }
  source_drafts(pe, flows, s);
  std::vector<std::pair<EdgeId, Amount>> call_drops;
  for (std::size_t i = 0; i < k; ++i) {
    FlowRecord& r = s.drafts[i].record;
    r.edges = carried[i];
    r.delivered = got[i];
    if (!lossy) continue;
    r.drops = plus(std::vector<std::pair<EdgeId, Counter128>> {}, drops[i]);
    r.offered = c128(flow.head_offered[i]);
    call_drops.insert(call_drops.end(), drops[i].begin(), drops[i].end());
  }
  std::vector<Counter128> assigned(n);
  std::vector<Amount> delivered(n, 0);
  std::vector<Counter128> dropped(n);
  for (std::size_t e = 0; e < n; ++e) {
    const auto& f = flow.entries[e];
    assigned[e] = c128(lossy ? f.assigned : f.delivered);
    delivered[e] = f.delivered;
    for (const auto& [x, d] : f.dropped) accumulate(dropped[e], c128(d), "drops");
  }
  const Counter128 requested = v.has_value() ? c128(*v) : Counter128 {};
  const Counter128 offered = lossy ? c128(flow.offered) : Counter128 {};
  PlacementResult& out = s.out;
  out.placed = lossy ? flow.delivered : flow.total;
  if (v.has_value()) out.remaining = *v - out.placed;
  out.call_drops = merged(std::move(call_drops));
  out.call_source_loss = flow.source_loss;
  bool active = false;
  for (std::size_t h = 0; h < view.num_heads() && !active; ++h) {
    for (const std::int32_t e : view.head_entries(h)) active = active || p.entries[static_cast<std::size_t>(e)].weight > 0;
  }
  if (p.pooled) active = !view.plan(0)->members.empty();
  out.termination = forwarding_termination(active, lossy, v, out.placed);
  out.program = program;
  out.by_tunnel = tunnel_rows(view, assigned, delivered, dropped);
  const auto after = draft_records(s.drafts);
  fill_state(pe, after, requested, !v.has_value(), offered, out);
  if (controller(pe)) out.total_cost = total_cost_of(g, after);
  s.next.by_tunnel = out.by_tunnel;
  s.next.requested = requested;
  s.next.offered = offered;
  s.next.prepared = true;
  s.next.infinite = !v.has_value();
  return s;
}

// ---------------------------------------------------------------- construction

} // namespace

const char* model_name(const Routing& r) noexcept {
  if (std::holds_alternative<CostOnly>(r)) return "CostOnly";
  if (std::holds_alternative<Global>(r)) return "Global";
  if (std::holds_alternative<CapacityAware>(r)) return "T1";
  return "Pinned";
}

ResolvedAxes resolve_config(const FlowPolicyConfig& cfg, bool cache, std::int32_t num_edges) {
  const auto& r = cfg.routing;
  if (r.valueless_by_exception()) invalid("routing: the variant holds no routing mode");
  // The enumeration values.
  if (cfg.placement.has_value()) check_member(*cfg.placement, "placement");
  check_member(cfg.admission, "admission");
  if (cfg.origin.has_value()) check_member(*cfg.origin, "origin");
  if (cfg.selection.has_value()) check_selection(*cfg.selection);
  if (cfg.weights.has_value()) check_member(*cfg.weights, "weights");
  // The model's row, its columns left to right: placement, origin, admission, selection, weights
  // (with edge_weights).
  ResolvedAxes out {cfg.selection, cfg.origin, cfg.weights};
  const bool cost_only = std::holds_alternative<CostOnly>(r);
  const bool global = std::holds_alternative<Global>(r);
  const bool capacity_aware = std::holds_alternative<CapacityAware>(r);
  const bool pinned = std::holds_alternative<Pinned>(r);
  const std::string model = model_name(r);
  if (global || pinned) {
    if (cfg.placement.has_value()) invalid(model + " does not accept placement");
  } else if (!cfg.placement.has_value()) {
    invalid(std::string(cost_only ? "CostOnly" : "CapacityAware") + " requires placement");
  }
  const bool max_flow = global || (cfg.placement.has_value() && *cfg.placement == FlowPlacement::MaxFlow);
  if (cost_only && !max_flow && cfg.origin == Origin::Optimized) {
    invalid("CostOnly with EQUAL_BALANCED or WEIGHTED requires origin=EQUAL");
  }
  if (capacity_aware && cfg.origin == Origin::Optimized) invalid("T1 requires origin=EQUAL");
  if (pinned && cfg.origin.has_value()) invalid("Pinned does not accept origin");
  if (!cfg.origin.has_value() && !pinned) {
    out.origin = cost_only && !max_flow ? Origin::Equal : capacity_aware ? Origin::Equal : Origin::Optimized;
  }
  if (cfg.admission == Admission::Lossy) {
    if (max_flow) invalid("MAX_FLOW requires admission=LOSSLESS");
    if (capacity_aware) invalid("CapacityAware requires admission=LOSSLESS");
  }
  if (global || pinned) {
    if (cfg.selection.has_value()) invalid(model + " does not accept selection");
  } else {
    if (!cfg.selection.has_value()) out.selection = PathSelection();
    check_kernel_combination(*cfg.placement, cfg.admission);
    // CapacityAware's Weighted flows read their weights from the priority (Capacity or Residual),
    // so Headroom (a ratio, not units) is rejected for them.
    if (capacity_aware && *cfg.placement == FlowPlacement::Weighted && out.selection->priority == EdgePriority::Headroom) {
      invalid("the Weighted kernel is rejected with priority Headroom (weights are units)");
    }
  }
  // The weights column: CostOnly Weighted reads them (Capacity by default); another CostOnly
  // placement names the placement it needs; every other model takes none.
  const bool weighted = cost_only && *cfg.placement == FlowPlacement::Weighted;
  if (weighted) {
    if (!cfg.weights.has_value()) out.weights = WeightBasis::Capacity;
    if (cfg.weights == WeightBasis::Explicit && !cfg.edge_weights) invalid("weights=EXPLICIT requires edge_weights");
    if (cfg.edge_weights && cfg.weights != WeightBasis::Explicit) invalid("edge_weights require weights=EXPLICIT");
  } else if (cost_only) {
    if (cfg.weights.has_value()) invalid("weights require CostOnly placement=WEIGHTED");
    if (cfg.edge_weights) invalid("edge_weights require weights=EXPLICIT");
  } else {
    if (cfg.weights.has_value()) invalid(model + " does not accept weights");
    if (cfg.edge_weights) invalid(model + " does not accept edge_weights");
  }
  if (cache && !cost_only) invalid(model + " does not use DagCache");
  if (const auto* ca = std::get_if<CapacityAware>(&r)) {
    if (ca->flow_count < 1) invalid("CapacityAware flow_count must be >= 1 (got " + str(ca->flow_count) + ")");
    if (ca->max_flow_count.has_value() && *ca->max_flow_count < ca->flow_count) {
      invalid("CapacityAware max_flow_count (" + str(*ca->max_flow_count) + ") must be >= flow_count (" +
              str(ca->flow_count) + ")");
    }
    if (ca->min_flow_bandwidth < 0 || ca->min_flow_bandwidth >= kDomainLimit) {
      invalid("CapacityAware min_flow_bandwidth must be in [0, 2^62) units (got " + str(ca->min_flow_bandwidth) + ")");
    }
    check_ceiling(ca->max_path_cost);
    if (ca->diminishing.has_value()) {
      const auto& d = *ca->diminishing;
      if (d.window < 1) invalid("diminishing window must be >= 1 (got " + str(d.window) + ")");
      if (!(d.fraction > 0.0 && d.fraction < 1.0)) {
        invalid("diminishing fraction must be in (0, 1) (got " + std::to_string(d.fraction) + ")");
      }
    }
  }
  // edge_weights' length and entries (they need the graph).
  if (weighted && cfg.edge_weights) {
    const auto& w = *cfg.edge_weights;
    if (w.size() != static_cast<std::size_t>(num_edges)) invalid("edge_weights: expected " + str(num_edges) + " entries");
    for (std::size_t e = 0; e < w.size(); ++e) {
      if (w[e] < 0 || w[e] >= kDomainLimit) invalid("edge_weights[" + std::to_string(e) + "] must be an integer in [0, 2^62)");
    }
  }
  return out;
}

void check_source_count(const FlowPolicyConfig& cfg, std::size_t num_sources, std::size_t num_destinations) {
  const auto& r = cfg.routing;
  (void)num_destinations;
  // The pooled single-path rule: a bottleneck choice would differ per source at a shared router.
  if (std::holds_alternative<CostOnly>(r) && num_sources > 1 && cfg.selection.has_value() &&
      !cfg.selection->multipath && cfg.selection->tie_break == EdgeTieBreak::PreferHigherPriority) {
    invalid("CostOnly with multiple input sources rejects multipath=false with PREFER_HIGHER_PRIORITY");
  }
}

} // namespace detail

// ================================================================ FlowPolicy

using detail::LedgerAccess;
using detail::PolicyEngine;

namespace {

// Some head can place traffic: a positive injection weight and a surviving entry of positive
// weight whose route has a member.
bool view_active(const detail::ForwardingView& v) {
  const auto& p = *v.data;
  for (std::size_t h = 0; h < v.num_heads(); ++h) {
    if (!v.weights().empty() && v.weights()[h] == 0) continue;
    for (const std::int32_t ei : v.head_entries(h)) {
      const auto e = static_cast<std::size_t>(ei);
      const auto& entry = p.entries[e];
      if (entry.weight == 0) continue;
      if (entry.tunnel >= 0) return true;
      const auto r = v.root_of(e);
      if (r >= 0 && v.plan(static_cast<std::size_t>(entry.context))->root_pos[static_cast<std::size_t>(r)] >= 0) return true;
    }
  }
  return false;
}

// The engine of a new policy, bound to its ledger as the last step. A supplied frozen cache's hit
// is recorded in `tally`: the public constructor publishes it once construction has
// succeeded, an evaluator once its whole call has.
std::unique_ptr<PolicyEngine> make_engine(FlowLedger& ledger, std::span<const NodeId> sources,
                                          std::span<const NodeId> destinations, FlowPolicyConfig config,
                                          const Eligibility* elig, DagCache* cache, detail::FrozenTally& tally) {
  auto pe = std::make_unique<PolicyEngine>(ledger, std::move(config));
  pe->graph = ledger.graph();
  const auto& g = *pe->graph;
  // The endpoints, then the ownership of the mask, the program, the support, the cache and the
  // black holes, then the configuration.
  pe->sources = detail::node_set(sources, g.num_nodes(), "sources");
  pe->destinations = detail::node_set(destinations, g.num_nodes(), "destinations");
  detail::check_disjoint(pe->sources, pe->destinations, "sources", "destinations");
  if (elig != nullptr && elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  const auto* pinned = std::get_if<Pinned>(&pe->config.routing);
  const auto* global = std::get_if<Global>(&pe->config.routing);
  if (pinned != nullptr) {
    if (!pinned->program) throw std::invalid_argument("Pinned requires a program");
    if (pinned->program->graph_uid() != g.uid()) throw std::invalid_argument("program belongs to another graph");
  }
  if (global != nullptr && global->support && global->support->graph_uid() != g.uid()) {
    throw std::invalid_argument("support belongs to another graph");
  }
  if (cache != nullptr && cache->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign DagCache: bound to another graph");
  }
  if (pinned != nullptr && pinned->black_holes && pinned->black_holes->graph_uid() != g.uid()) {
    throw std::invalid_argument("black_holes belongs to another graph");
  }
  pe->elig = elig;
  const auto axes = detail::resolve_config(pe->config, cache != nullptr, g.num_edges());
  pe->config.selection = axes.selection;
  pe->config.origin = axes.origin;
  pe->config.weights = axes.weights;
  detail::check_source_count(pe->config, pe->sources.size(), pe->destinations.size());
  if (pinned != nullptr) {
    const auto& prog = *pinned->program;
    if (!std::equal(prog.sources().begin(), prog.sources().end(), pe->sources.begin(), pe->sources.end()) ||
        !std::equal(prog.destinations().begin(), prog.destinations().end(), pe->destinations.begin(),
                    pe->destinations.end())) {
      throw std::invalid_argument("program endpoints do not match policy endpoints");
    }
  }
  if (pe->config.selection.has_value()) pe->selection = *pe->config.selection;
  if (pe->config.origin.has_value()) pe->origin = *pe->config.origin;
  auto node_ok = [&](NodeId v) { return elig == nullptr || elig->node_ok(v); };
  for (NodeId d : pe->destinations) {
    if (node_ok(d)) pe->dests_ok.push_back(d);
  }
  if (pinned != nullptr) {
    // The program as this policy forwards it: replayed under the mask, the black holes laid over.
    pe->view = detail::make_view(g, pinned->program, elig, pe->dests_ok, pinned->black_holes.get(), true);
    for (std::size_t h = 0; h < pe->view.num_heads(); ++h) pe->roster.push_back(pe->view.head(h));
    pe->active = view_active(pe->view);
    pe->epoch.cohorts.resize(detail::ProgramAccess::data(*pinned->program).entries.size());
    pe->owner = LedgerAccess::bind(ledger);
    return pe;
  }
  if (detail::is<CostOnly>(*pe)) {
    pe->support = cache != nullptr
                      ? detail::DagCacheAccess::lookup(*cache, g, pe->destinations, Direction::Reverse, elig, tally)
                      : std::make_shared<const PredDAG>(detail::static_support(g, pe->destinations, Direction::Reverse, elig));
    const auto edges = detail::PredDAGBuilder::ascending(*pe->support);
    const auto row = g.row_offsets();
    for (NodeId s : pe->sources) {
      if (!node_ok(s)) continue;
      const auto it = std::lower_bound(edges.begin(), edges.end(), row[static_cast<std::size_t>(s)]);
      if (it != edges.end() && *it < row[static_cast<std::size_t>(s) + 1]) pe->roster.push_back(s);
    }
  } else if (global != nullptr && global->support) {
    // The support trimmed to the walks from the sources to the allowed destinations under the mask;
    // the roster is the sources that reach one inside it.
    std::vector<EdgeId> trimmed;
    for (NodeId s : pe->sources) {
      if (!node_ok(s)) continue;
      auto part = detail::prune(g, *global->support, s, pe->dests_ok, elig);
      if (part.empty()) continue;
      pe->roster.push_back(s);
      trimmed.insert(trimmed.end(), part.begin(), part.end());
    }
    std::sort(trimmed.begin(), trimmed.end());
    trimmed.erase(std::unique(trimmed.begin(), trimmed.end()), trimmed.end());
    pe->global_support = std::make_shared<const std::vector<EdgeId>>(std::move(trimmed));
  } else {
    // The static costs fix the roster; only a Stretch ceiling reads them again, so only it keeps them.
    const auto* ca = std::get_if<CapacityAware>(&pe->config.routing);
    const bool keep_costs =
        ca != nullptr && ca->max_path_cost.has_value() && std::holds_alternative<Stretch>(*ca->max_path_cost);
    const auto dist = detail::static_costs(g, elig, pe->dests_ok, pe->sources);
    for (std::size_t i = 0; i < pe->sources.size(); ++i) {
      const NodeId s = pe->sources[i];
      if (!node_ok(s) || dist[i] == kUnreachable) continue;
      pe->roster.push_back(s);
      if (keep_costs) pe->baseline.push_back(dist[i]);
    }
    if (global != nullptr && elig != nullptr && elig->uid() != 0) {
      // Its programs' replay basis: the topology under this mask.
      for (NodeId v = 0; v < g.num_nodes(); ++v) {
        if (!elig->node_ok(v)) pe->excluded_nodes.push_back(v);
      }
      for (EdgeId e = 0; e < g.num_edges(); ++e) {
        if (!elig->edge_ok(e)) pe->excluded_edges.push_back(e);
      }
    }
  }
  pe->epoch.cohorts.resize(pe->roster.size());
  pe->owner = LedgerAccess::bind(ledger);  // last: nothing is bound when construction fails
  return pe;
}

} // namespace

// A public call of its own: a supplied frozen cache counts the hit once the policy exists.
FlowPolicy::FlowPolicy(FlowLedger& ledger, std::span<const NodeId> sources, std::span<const NodeId> destinations,
                       FlowPolicyConfig config, const Eligibility* elig, DagCache* cache)
    : engine_(detail::frozen_call(cache, [&](detail::FrozenTally& tally) {
        return make_engine(ledger, sources, destinations, std::move(config), elig, cache, tally);
      })) {}

FlowPolicy::FlowPolicy(FlowLedger& ledger, std::span<const NodeId> sources, std::span<const NodeId> destinations,
                       FlowPolicyConfig config, const Eligibility* elig, DagCache* cache, detail::FrozenTally& tally)
    : engine_(make_engine(ledger, sources, destinations, std::move(config), elig, cache, tally)) {}

// The destruction law: the records stay in the ledger as orphans.
FlowPolicy::~FlowPolicy() { LedgerAccess::release(engine_->fg, engine_->owner); }

namespace detail {
namespace {

// The program a computing model installs for an empty roster (R22): CostOnly and Global an empty
// one, CapacityAware none. Pinned never comes here: its forwarding with no head is place_pinned's.
std::shared_ptr<const ForwardingProgram> empty_roster_program(const PolicyEngine& pe) {
  if (is<CapacityAware>(pe)) return nullptr;
  return empty_program(pe.graph->uid(), pe.sources, pe.destinations, controller(pe) && pe.origin == Origin::Optimized);
}

// The current state of a zero call: its records, its program and rows, every flag false; before
// the epoch's first preparation all zeros (no records, no program, no curve). A zero call solves
// nothing, so it claims no cut.
PlacementResult zero_state(const PolicyEngine& pe, const OutputOptions& outputs) {
  PlacementResult out;
  out.remaining = 0;
  const auto& ep = pe.epoch;
  const auto records = current_records(pe);
  std::vector<std::pair<NodeId, Counter128>> lost;
  for (std::size_t h = 0; h < ep.head_lost.size(); ++h) {
    if (!ep.head_lost[h].is_zero()) lost.emplace_back(pe.view.head(h), ep.head_lost[h]);
  }
  fill_state(pe, records, ep.requested, ep.infinite, ep.offered, out, lost);
  if (controller(pe)) out.total_cost = total_cost_of(*pe.graph, records);
  if (dense_global(pe) && pe.origin == Origin::Optimized && ep.prepared && ep.curve.has_value()) out.cost_distribution = *ep.curve;
  out.program = ep.program;
  out.by_tunnel = ep.by_tunnel;
  prepare_outputs(pe, outputs, current_edges(pe), out);
  return out;
}

// A computing model's call when no source reaches an allowed destination (Pinned with no head
// forwards through place_pinned, keeping its per-entry records): nothing is forwarded. A Lossy
// offer is lost before it enters (the demand's source loss; no source's row has it). Global's flow
// is 0 and maximal: its curve is empty and its cut is empty (no roster source, so the closure is
// empty).
Staged empty_roster(const PolicyEngine& pe, const Epoch& cur, const Request& rq, const ForwardingProgram* before) {
  const bool lossy = pe.config.admission == Admission::Lossy;
  const bool optimized_global = dense_global(pe) && pe.origin == Origin::Optimized;
  const auto v = rq.v;
  Staged s;
  s.next.cohorts.resize(cur.cohorts.size());
  s.next.requested = requested_after(cur, v);
  s.next.offered = cur.offered;
  if (lossy && v.has_value()) (void)s.next.offered.add(*v);  // offered <= requested
  s.next.prepared = true;
  s.next.infinite = cur.infinite || !v.has_value();
  s.next.program = cur.program ? cur.program : empty_roster_program(pe);
  PlacementResult& out = s.out;
  out.remaining = v;
  out.delivery_delta = 0;
  out.termination = Termination::NoRoute;
  if (lossy && v.has_value()) out.call_source_loss = *v;
  set_flags(out, before, s.next.program.get());
  out.program = s.next.program;
  fill_state(pe, {}, s.next.requested, s.next.infinite, s.next.offered, out);
  if (controller(pe)) out.total_cost = Counter128 {};
  if (optimized_global) {
    out.cost_distribution.emplace();
    if (rq.outputs != nullptr && rq.outputs->min_cut) out.min_cut.emplace();
    s.next.curve.emplace();
  }
  if (kTraced && rq.trace) rq.trace->set("roster", Trace::list());
  return s;
}

} // namespace

void check_forward_operation(const PolicyEngine& pe) {
  if (is<Pinned>(pe)) invalid("Pinned forwards with replace");
  if (!pe.epoch.program) invalid("forward requires an installed program");
}

namespace {

// The call's validation (before any change, also before the zero shortcut), in this order: the
// operation for the model and an installed program (forward), the offer, a Lossy full offer's
// finiteness, a finite place after +inf, the requested total's limit, the black holes' graph.
void check_call(const PolicyEngine& pe, const Request& rq) {
  const bool lossy = pe.config.admission == Admission::Lossy;
  if (rq.op == Op::Forward) check_forward_operation(pe);
  if (rq.v.has_value() && (*rq.v < 0 || *rq.v >= kDomainLimit)) {
    invalid("offer must be an integer in [0, 2^62) or permitted +inf");
  }
  // A Lossy full offer is finite; CostOnly's place(+inf) is its first-hop experiment.
  if (lossy && !rq.v.has_value() && (rq.full() || is<Pinned>(pe))) invalid("LOSSY forwarding requires a finite offer");
  if (rq.op == Op::Place) {
    if (pe.epoch.infinite && rq.v.has_value()) invalid("finite placement after infinite placement requires a new epoch");
    // The only check of the cumulative request: nothing changes the requested total between here
    // and the call's publication, so the later additions repeat this sum and cannot fail.
    if (rq.v.has_value()) {
      Counter128 probe = pe.epoch.requested;
      if (!probe.add(*rq.v)) invalid("requested limit exceeded");
    }
  }
  if (rq.op == Op::Forward && rq.dark != nullptr && rq.dark->graph_uid() != pe.graph->uid()) {
    invalid("black_holes belongs to another graph");
  }
}

// The computation of a call: nothing published. A full-offer call works on the residual with
// the policy's own load released, from a fresh epoch.
Staged compute(const PolicyEngine& pe, const Request& rq) {
  Working W(pe.fg);
  if (rq.full()) W.release_records(current_records(pe));
  const Epoch fresh = rq.full() ? fresh_epoch(pe, rq.op == Op::Forward) : Epoch {};
  const Epoch& cur = rq.full() ? fresh : pe.epoch;
  // The flags compare the program the policy has installed now.
  const ForwardingProgram* before = pe.epoch.program.get();
  Staged s;
  if (pe.roster.empty() && !is<Pinned>(pe)) {
    // No roster source (a computing model): nothing is computed. A Pinned policy whose heads all
    // fail keeps its per-entry records and forwards nothing (place_pinned).
    s = empty_roster(pe, cur, rq, before);
  } else if (rq.op == Op::Forward) {
    s = forward_installed(pe, W, rq);
  } else if (is<Pinned>(pe)) {
    s = place_pinned(pe, cur, W, rq, before);
  } else if (controller(pe)) {
    s = place_controller(pe, cur, W, rq, before);
  } else if (is<CostOnly>(pe)) {
    s = place_pooled(pe, cur, W, rq, before);
  } else {
    s = place_flows(pe, cur, W, rq);
  }
  // A full-offer call: its placed is the whole new delivery, compared with the old one.
  if (rq.full()) {
    Amount old = 0;
    for (const FlowRecord* r : current_records(pe)) old += r->delivered;
    s.out.placed = s.out.delivered;
    s.out.delivery_delta = s.out.delivered - old;
    if (rq.v.has_value()) s.out.remaining = *rq.v - s.out.placed;
  }
  // The first termination row (14.2): a request of 0 is COMPLETE, also for a full-offer call, which
  // still prepared and publishes what it computed.
  if (rq.v == Amount {0}) s.out.termination = Termination::Complete;
  if (rq.outputs != nullptr && any_output(*rq.outputs)) prepare_outputs(pe, *rq.outputs, draft_edges(s.drafts), s.out);
  return s;
}

// The publication: the ledger's one atomic exchange (a place commits every live record and the new
// ones; a full-offer call retires the live ones and creates the new ones), then the epoch.
void publish(PolicyEngine& pe, Staged& s, bool full) {
  const auto ids = full ? LedgerAccess::exchange(pe.fg, pe.owner, s.drafts) : LedgerAccess::commit(pe.fg, pe.owner, s.drafts);
  // Committed: the epoch state follows (nothing below throws).
  for (std::size_t n = 0; n < s.fresh.size(); ++n) s.next.cohorts[s.fresh[n].first][s.fresh[n].second].id = ids[n].second;
  pe.epoch = std::move(s.next);
}

// `staged` (optional) runs once the call's result and trace are final and before anything is
// published (the ledger commit, the epoch state): an exception from it aborts the call with
// nothing changed (the Python binding builds its result objects there).
PlacementResult run_call(PolicyEngine& pe, const Request& rq, const std::function<void()>* staged) {
  check_call(pe, rq);
  if (rq.op == Op::Place && rq.v.has_value() && *rq.v == 0) {
    PlacementResult out = zero_state(pe, *rq.outputs);
    if (staged != nullptr && *staged) (*staged)();
    return out;
  }
  Staged s = compute(pe, rq);
  if (staged != nullptr && *staged) (*staged)();
  if (rq.op == Op::Evaluate) return std::move(s.out);
  PlacementResult out = std::move(s.out);
  publish(pe, s, rq.full());
  return out;
}

} // namespace
} // namespace detail

PlacementResult FlowPolicy::place(std::optional<Amount> volume, const OutputOptions& outputs) {
  return detail::run_call(*engine_, detail::Request {detail::Op::Place, volume, nullptr, &outputs, nullptr}, nullptr);
}

PlacementResult FlowPolicy::evaluate(std::optional<Amount> offer, const OutputOptions& outputs) const {
  const detail::Request rq {detail::Op::Evaluate, offer, nullptr, &outputs, nullptr};
  detail::check_call(*engine_, rq);
  return std::move(detail::compute(*engine_, rq).out);
}

PlacementResult FlowPolicy::replace(std::optional<Amount> offer, const OutputOptions& outputs) {
  return detail::run_call(*engine_, detail::Request {detail::Op::Replace, offer, nullptr, &outputs, nullptr}, nullptr);
}

PlacementResult FlowPolicy::forward(std::optional<Amount> offer, const Eligibility* black_holes,
                                    const OutputOptions& outputs) {
  return detail::run_call(*engine_, detail::Request {detail::Op::Forward, offer, black_holes, &outputs, nullptr}, nullptr);
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
PlacementResult detail::place_traced(FlowPolicy& policy, std::optional<Amount> volume, const OutputOptions& outputs,
                                     Trace& trace, const std::function<void()>* staged) {
  trace = Trace::map();
  return run_call(PolicyAccess::engine(policy), Request {Op::Place, volume, nullptr, &outputs, &trace}, staged);
}

PlacementResult detail::call_staged(FlowPolicy& policy, CallKind kind, std::optional<Amount> offer,
                                    const Eligibility* black_holes, const OutputOptions& outputs,
                                    const std::function<void()>* staged) {
  const Op op = kind == CallKind::Place      ? Op::Place
                : kind == CallKind::Evaluate ? Op::Evaluate
                : kind == CallKind::Replace  ? Op::Replace
                                             : Op::Forward;
  return run_call(PolicyAccess::engine(policy), Request {op, offer, black_holes, &outputs, nullptr}, staged);
}
#endif

void FlowPolicy::remove() {
  auto& pe = *engine_;
  (void)LedgerAccess::detach(pe.fg, pe.owner);
  // A new epoch (nothing below throws): the same number of cohorts, every one empty.
  for (auto& c : pe.epoch.cohorts) c.clear();
  pe.epoch.program.reset();
  pe.epoch.support.clear();
  pe.epoch.curve.reset();
  pe.epoch.head_lost.clear();
  pe.epoch.by_tunnel.clear();
  pe.epoch.requested = Counter128 {};
  pe.epoch.offered = Counter128 {};
  pe.epoch.infinite = false;
  pe.epoch.prepared = false;
}

Amount FlowPolicy::placed_demand() const {
  const auto& pe = *engine_;
  Amount total = 0;
  for (const auto& c : pe.epoch.cohorts) {
    for (const auto& s : c) total += LedgerAccess::live(pe.fg, s.id).delivered;
  }
  return total;
}

std::vector<FlowRecord> FlowPolicy::records() const {
  const auto& pe = *engine_;
  std::vector<std::int64_t> ids;
  for (const auto& c : pe.epoch.cohorts) {
    for (const auto& s : c) ids.push_back(s.id);
  }
  std::sort(ids.begin(), ids.end());
  std::vector<FlowRecord> out;
  out.reserve(ids.size());
  for (std::int64_t id : ids) out.push_back(LedgerAccess::live(pe.fg, id));
  return out;
}

FlowRecord FlowPolicy::record(const FlowId& flow_id) const {
  const auto& pe = *engine_;
  const FlowLedger& fg = pe.fg;
  // The ledger's checks first (another ledger, a retired or unknown id), then the owner.
  FlowRecord rec = fg.record(flow_id);
  if (rec.owner != pe.owner) {
    throw std::invalid_argument("FlowId (" + std::to_string(flow_id.first) + ", " + std::to_string(flow_id.second) +
                                ") is not a record of this policy");
  }
  return rec;
}

std::shared_ptr<const ForwardingProgram> FlowPolicy::program() const noexcept { return engine_->epoch.program; }

const FlowPolicyConfig& FlowPolicy::config() const noexcept { return engine_->config; }

} // namespace netgraph::core
