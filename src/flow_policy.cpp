/*
  The policy engine FlowPolicy: one placement transaction for every routing
  mode, committed through the ledger as one atomic exchange per call.

  * CostOnly EqualBalanced and Weighted (place_pooled): the epoch's first positive call prepares
    every roster source's route on the working residual at its entry and joins them into one
    program, a plan rooted at the whole roster (a router reached by several sources splits their
    added traffic once). Lossless admits the largest total of the equal origin whose no-drop
    forwarding fits V = W + own; Lossy forwards the round-robin increments of the offer on W.
    The amounts are attributed back to the sources (kernels' forward()).
  * The controllers (place_controller): Global over every eligible edge, CostOnly MaxFlow over
    the union of the roster's routes (prepared as above). A minimum-cost flow of the cumulative
    value on V = W + own, from one virtual source (the optimized origin) or with source i's arc
    capped at its round-robin share of the largest total whose shares all fit (the equal origin,
    integer bisection). Split by source in the fixed order described by decompose().
  * CapacityAware and Pinned (place_flows) visit the roster in ascending NodeId. Round robin
    gives each source an equal integer share of the cumulative request, with leftover units
    going to earlier sources. Offer each source the increase in its share, then run stages 1
    creation, 2 reoptimisation, 3 amounts and 4 growth on shared unused capacity W (the working
    residual). One `created` counter covers the call. Pinned has a single source.

  Work on copies and commit all changes together; an exception leaves ledger and policy intact.
  Update policy epoch state (accumulated since construction or remove()) after commit, without
  throwing. Store only changes to ledger residual in per-thread slots. Pinned and split CostOnly
  gather their program members and attributed records. Global reads a dense residual and solves
  over every edge (its view and flows per EdgeId); CostOnly MaxFlow gathers its view over its
  support and keeps its flows aligned with it, expanding them over all edges only to split the flow
  among several sources (decompose()) or for the trace. CapacityAware reads a dense residual copy;
  CostOnly preparation reads the working residual edge by edge.

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
#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
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
  std::int32_t bundle_index {-1};
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
  PolicyEngine& pe;
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

// A flow's plan on its route (the record's bundle) with its membership captured on the creation
// view (Weighted Residual: `dense_view` when the search read one, else the working residual W);
// with `cost`, the route's minimum walk cost.
std::shared_ptr<const KernelPlan> plan_on(const Call& c, NodeId src, const PredDAG& route,
                                          const std::vector<Amount>* dense_view, Cost* cost) {
  const auto& pe = c.pe;
  std::vector<Amount> creation;
  if (*pe.config.placement == FlowPlacement::Weighted && pe.selection.priority == EdgePriority::Residual) {
    const auto edges = PredDAGBuilder::ascending(route);
    creation.resize(edges.size());
    for (std::size_t i = 0; i < edges.size(); ++i) {
      creation[i] = dense_view != nullptr ? (*dense_view)[static_cast<std::size_t>(edges[i])] : c.W.at(edges[i]);
    }
  }
  return std::make_shared<const KernelPlan>(make_plan(c.g, *pe.config.placement, pe.config.admission,
                                                      pe.selection.priority, src, pe.dests_ok, route, creation, cost));
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

// CostOnly's creation products in one pass (kernels' cost_only_creation, identical to
// cost_only_route, bundle_of and plan_on where it applies); nullopt elsewhere.
std::optional<CreationProducts> creation_products(const Call& c, NodeId src) {
  const auto& pe = c.pe;
  const auto& sel = pe.selection;
  const WorkingPriority view {&c.W, c.g.capacity_units(), sel.priority};
  return cost_only_creation(c.g, pe.support, src, pe.dests_ok, sel, view.on_w(), *pe.config.placement,
                            pe.config.admission);
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
    for (const auto& cohort : pe.cohorts) {
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
  for (const auto& cohort : pe.cohorts) {
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

// The demand's state after a call: its records after the call (`after`: one per flow, or one
// per roster source; the record's source names its row) and the epoch's counters after it. Lossy:
// a row's source loss is what its source was offered and neither delivered nor dropped; the
// demand's is its offered total minus delivery and drops, which with an empty roster is the
// whole offer while every row stays zero.
void fill_state(const PolicyEngine& pe, std::span<const FlowRecord* const> after, const Counter128& requested,
                bool infinite, const Counter128& offered, PlacementResult& out) {
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
  for (const auto& cohort : pe.cohorts) {
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

// Global, or CostOnly MaxFlow: a controller.
bool controller(const PolicyEngine& pe) noexcept {
  return is<Global>(pe) || (is<CostOnly>(pe) && *pe.config.placement == FlowPlacement::MaxFlow);
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
  auto& pe = c.pe;
  Trace tr = Trace::map();
  if (is<Pinned>(pe)) {
    const bool captured = *pe.config.placement == FlowPlacement::Weighted &&
                          pe.selection.priority == EdgePriority::Residual;
    Trace usable = Trace::list();
    for (std::int32_t i : pe.usable) {
      const auto& b = pe.bundles[static_cast<std::size_t>(i)];
      WorkFlow f = new_flow(c, st.src, b.route_dag, b.cost, nullptr, captured ? nullptr : b.route_plan);
      f.bundle_index = i;
      f.seq = ++c.seq;
      st.flows.push_back(std::move(f));
      if (kTraced && c.trace) usable.push(Trace::integer(i));
    }
    if (kTraced && c.trace) tr.set("usable", std::move(usable));
  } else {
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

void stage3_lossy(Call& c, SourceRun& st) {
  auto& flows = st.flows;
  std::vector<Amount> offers(flows.size(), 0);
  if (st.offer.has_value()) round_robin(*st.offer, offers);
  Trace rows = Trace::list();
  std::vector<Amount> w_m;
  for (std::size_t i = 0; i < flows.size(); ++i) {
    auto& f = flows[i];
    const std::optional<Amount> offer = st.offer.has_value() ? std::optional<Amount>(offers[i]) : std::nullopt;
    const auto members = f.members();
    c.W.gather(members, w_m);
    Injection in;
    if (offer.has_value()) {
      in.per_root = std::span<const Amount>(&offers[i], 1);
    } else {
      in.first_hop = true;
    }
    const auto step = forward(*f.plan, in, w_m);
    // The 2^127 checks below can fail: finite offers bring an offered total up to the request
    // total (below 2^127), and a +inf call then adds its injection. A drop counter (at most its
    // flow's offered counter) can reach 2^127 only in such a call, and its check runs first.
    c.W.consume(members, step.carried);
    for (std::size_t m = 0; m < members.size(); ++m) {
      f.alloc[m] += step.carried[m];
      if (step.drops[m] != 0) {
        if (!f.drops[m].add(step.drops[m])) throw std::invalid_argument("drops limit exceeded");
        c.drops.emplace_back(members[m], step.drops[m]);
      }
    }
    if (!f.offered->add(step.offered)) throw std::invalid_argument("offered limit exceeded");
    f.delivered += step.delivered;
    st.injected += step.injected;
    c.offered += step.offered;
    c.source_loss += step.source_loss;
    if (kTraced && c.trace) {
      Trace row = Trace::map();
      row.set("offer", amount_or_inf(offer));
      row.set("offered", Trace::integer(step.offered));
      row.set("delivered", Trace::integer(step.delivered));
      std::vector<std::pair<EdgeId, Amount>> carried, drops;
      for (std::size_t m = 0; m < members.size(); ++m) {
        if (step.carried[m] != 0) carried.emplace_back(members[m], step.carried[m]);
        if (step.drops[m] != 0) drops.emplace_back(members[m], step.drops[m]);
      }
      row.set("carried", Trace::pairs(carried));
      row.set("drops", Trace::pairs(drops));
      rows.push(std::move(row));
    }
  }
  if (kTraced && st.tr) {
    st.tr->set("stage3", std::move(rows));
    st.tr->set("injected", Trace::integer(st.injected));
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
  auto& pe = c.pe;
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
  auto& pe = c.pe;
  if (st.flows.empty()) stage1(c, st);
  if (st.flows.empty()) return;
  if (is<CapacityAware>(pe) && std::get<CapacityAware>(pe.config.routing).reoptimize) stage2(c, st);
  if (pe.config.admission == Admission::Lossy) {
    stage3_lossy(c, st);
  } else if (*pe.config.placement == FlowPlacement::EqualBalanced) {
    stage3_eb(c, st);
  } else {
    stage3_seq(c, st);
  }
  if (is<CapacityAware>(pe)) stage4(c, st);
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
  f.bundle_index = s.bundle_index;
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
    for (const auto& f : st.flows) c.push_back(LiveFlow {f.id, f.src, f.plan, f.bundle, f.cost, f.bundle_index});
  }
  return out;
}

std::vector<std::pair<EdgeId, Amount>> merge_drops(std::vector<std::pair<EdgeId, Amount>> drops) {
  std::sort(drops.begin(), drops.end());
  std::vector<std::pair<EdgeId, Amount>> out;
  for (const auto& [e, x] : drops) {
    if (!out.empty() && out.back().first == e) {
      out.back().second += x;
    } else {
      out.emplace_back(e, x);
    }
  }
  return out;
}

std::string key(NodeId v) { return std::to_string(v); }

PlacementResult place_flows(PolicyEngine& pe, std::optional<Amount> v, const OutputOptions& outputs, Trace* trace,
                            const std::function<void()>* staged) {
  const auto& g = *pe.graph;
  Working W(pe.fg);
  Call c {pe, g, W};
#if NETGRAPH_CORE_PRIVATE_HOOKS
  c.cap = max_flows_per_call().load(std::memory_order_relaxed);
#endif
  c.trace = trace;
  const std::size_t k = pe.roster.size();
  std::vector<SourceRun> runs(k);
  if (is<Pinned>(pe)) {
    runs[0].offer = v;
  } else if (v.has_value()) {
    cumulative_offers(pe.requested, *v, k, [&](std::size_t i, Amount x) { runs[i].offer = x; });
  }
  Trace sources_tr = Trace::map();
  for (std::size_t i = 0; i < k; ++i) {
    auto& st = runs[i];
    st.ri = i;
    st.src = pe.roster[i];
    st.B = st.offer;
    st.flows.reserve(pe.cohorts[i].size());
    for (const auto& s : pe.cohorts[i]) st.flows.push_back(working_copy(pe, s));
    for (const auto& f : st.flows) {
      st.P += pe.config.admission == Admission::Lossy ? f.delivered : f.T;
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
  std::vector<Draft> drafts;
  // (run, flow) of the new flows, per thread: cleared here, read until the epoch state is set.
  thread_local std::vector<std::pair<std::size_t, std::size_t>> fresh_scratch;
  auto& fresh = fresh_scratch;
  fresh.clear();
  for (std::size_t i = 0; i < k; ++i) {
    for (std::size_t j = 0; j < runs[i].flows.size(); ++j) {
      const auto& f = runs[i].flows[j];
      if (f.is_new()) {
        fresh.emplace_back(i, j);
      } else {
        drafts.push_back(draft_of(pe, f));
      }
    }
  }
  std::sort(fresh.begin(), fresh.end(), [&](const auto& a, const auto& b) {
    return runs[a.first].flows[a.second].seq < runs[b.first].flows[b.second].seq;
  });
  for (const auto& [i, j] : fresh) drafts.push_back(draft_of(pe, runs[i].flows[j]));
  // Everything the publication of the epoch state needs is allocated before the commit.
  const bool lossy = pe.config.admission == Admission::Lossy;
  Counter128 requested = pe.requested;
  if (v.has_value()) (void)requested.add(*v);  // checked in place_impl
  Counter128 offered = pe.offered;
  if (!offered.add(c.offered)) throw std::invalid_argument("offered limit exceeded");
  auto next = next_cohorts(runs, k);
  Amount placed = 0;
  std::vector<Termination> outcomes;  // CapacityAware's per-source outcomes (and the trace's)
  if (is<CapacityAware>(pe)) outcomes.reserve(k);
  Termination term = Termination::Complete;
  bool any_member = false;
  for (const auto& st : runs) {
    Amount after = 0;
    for (const auto& f : st.flows) {
      after += lossy ? f.delivered : f.T;
      any_member = any_member || !f.plan->members.empty();
    }
    placed += after - st.P;
    if (is<CapacityAware>(pe)) {
      outcomes.push_back(local_outcome(st));
      if (rank(outcomes.back()) > rank(term)) term = outcomes.back();
    }
  }
  if (is<Pinned>(pe)) {
    // Termination's rule (the one source has its usable bundles).
    if (!any_member) {
      term = Termination::NoActiveWeight;
    } else if (lossy) {
      term = v.has_value() ? Termination::Complete : Termination::Exhausted;
    } else {
      term = v.has_value() && placed == *v ? Termination::Complete : Termination::Exhausted;
    }
    if (kTraced && trace) outcomes.assign(k, term);  // the trace's per-source outcomes
  }
  PlacementResult out;
  out.placed = placed;
  if (v.has_value()) out.remaining = *v - placed;
  out.delivery_delta = placed;
  out.termination = term;
  out.routes_changed = c.routes_changed;
  out.weights_changed = c.weights_changed;
  out.call_drops = merge_drops(std::move(c.drops));
  out.call_source_loss = c.source_loss;
  fill_state(pe, draft_records(drafts), requested, pe.infinite || !v.has_value(), offered, out);
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
  if (any_output(outputs)) prepare_outputs(pe, outputs, draft_edges(drafts), out);
  if (staged != nullptr && *staged) (*staged)();
  std::vector<FlowId> ids;
  if (!drafts.empty()) ids = LedgerAccess::commit(pe.fg, pe.owner, drafts);
  // Committed: the epoch state follows (nothing below throws).
  for (std::size_t n = 0; n < fresh.size(); ++n) next[fresh[n].first][fresh[n].second].id = ids[n].second;
  pe.cohorts.swap(next);
  pe.requested = requested;
  pe.offered = offered;
  pe.prepared = true;
  if (!v.has_value()) pe.infinite = true;
  return out;
}

// ---------------------------------------------------------------- Global

// Split flow by source in ascending NodeId. While a source has supply left, take the smallest
// EdgeId sequence with positive remaining flow, stopping at the first destination. Assign
// min(supply left, smallest remaining edge flow on that walk), then subtract it from the walk.
// Each node's out-edge cursor only advances: remaining flow never increases.
std::vector<std::vector<std::pair<EdgeId, Amount>>> decompose(const Graph& g, std::vector<Amount> rem,
                                                              std::span<const NodeId> roster,
                                                              std::span<const NodeId> dests) {
  const auto row = g.row_offsets();
  const auto edge_dst = g.edge_dst();
  const auto in_row = g.in_row_offsets();
  const auto in_adj = g.in_adj_edge_index();
  std::vector<std::uint8_t> is_dst(static_cast<std::size_t>(g.num_nodes()), 0);
  for (NodeId d : dests) is_dst[static_cast<std::size_t>(d)] = 1;
  std::vector<std::int32_t> cur(row.begin(), row.end() - 1);
  std::vector<std::vector<std::pair<EdgeId, Amount>>> out(roster.size());
  std::vector<EdgeId> walk;
  for (std::size_t i = 0; i < roster.size(); ++i) {
    const NodeId s = roster[i];
    const auto si = static_cast<std::size_t>(s);
    // The supply of s: its net outflow (transit through s from another source is not supply).
    Amount left = 0;
    for (auto e = row[si]; e < row[si + 1]; ++e) left += rem[static_cast<std::size_t>(e)];
    for (auto k = in_row[si]; k < in_row[si + 1]; ++k) left -= rem[static_cast<std::size_t>(in_adj[static_cast<std::size_t>(k)])];
    std::vector<std::pair<EdgeId, Amount>> acc;
    while (left > 0) {
      walk.clear();
      NodeId u = s;
      while (walk.empty() || is_dst[static_cast<std::size_t>(u)] == 0) {
        auto& p = cur[static_cast<std::size_t>(u)];
        while (p < row[static_cast<std::size_t>(u) + 1] && rem[static_cast<std::size_t>(p)] == 0) ++p;
        if (p == row[static_cast<std::size_t>(u) + 1]) {
          throw std::logic_error("Global decomposition: flow conservation broken (no continuation)");
        }
        walk.push_back(p);
        u = edge_dst[static_cast<std::size_t>(p)];
        if (walk.size() > static_cast<std::size_t>(g.num_nodes())) {
          throw std::logic_error("Global decomposition: a cycle with positive flow");
        }
      }
      Amount amount = left;
      for (EdgeId e : walk) amount = std::min(amount, rem[static_cast<std::size_t>(e)]);
      for (EdgeId e : walk) {
        rem[static_cast<std::size_t>(e)] -= amount;
        acc.emplace_back(e, amount);
      }
      left -= amount;
    }
    std::sort(acc.begin(), acc.end());
    auto& res = out[i];
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

// A CostOnly demand's routes for the epoch, prepared at its first positive call on the working
// residual at that call's entry: per roster source its route (the record's bundle and route
// cost), and the pooled program (EqualBalanced, Weighted: one plan rooted at the roster) or the
// routes' union (MaxFlow, the controller's support).
struct Prepared {
  std::vector<LiveFlow> flows {};  // per roster source, new (id -1)
  std::shared_ptr<const KernelPlan> program {};
  std::vector<EdgeId> edges {};    // ascending: the union of the routes
};

Prepared prepare_cost_only(Call& c) {
  auto& pe = c.pe;
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
  // One split-routing source: its one-pass creation plan is the program (its route is the union).
  if (single) {
    out.program = std::move(single);
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
  std::vector<Amount> creation;
  if (placement == FlowPlacement::Weighted && pe.selection.priority == EdgePriority::Residual) {
    creation.resize(out.edges.size());
    for (std::size_t i = 0; i < out.edges.size(); ++i) creation[i] = c.W.at(out.edges[i]);
  }
  out.program = std::make_shared<const KernelPlan>(make_plan(g, placement, pe.config.admission, pe.selection.priority,
                                                             pe.roster, pe.dests_ok, *dag, creation));
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

// The commit of a demand with one flow per roster source (CostOnly, Global): the drafts in
// roster order (all retained, or all new at the epoch's first call), then the epoch's flows.
std::vector<std::vector<LiveFlow>> commit_sources(PolicyEngine& pe, std::vector<Draft>& drafts,
                                                  const std::vector<LiveFlow>& flows) {
  std::vector<std::vector<LiveFlow>> next(flows.size());
  for (std::size_t i = 0; i < flows.size(); ++i) next[i].push_back(flows[i]);
  const auto ids = LedgerAccess::commit(pe.fg, pe.owner, drafts);
  // Committed: nothing below throws. New records take their ids in roster order.
  for (std::size_t i = 0, n = 0; i < next.size(); ++i) {
    if (next[i].front().id < 0) next[i].front().id = ids[n++].second;
  }
  return next;
}

// ---------------------------------------------------------------- CostOnly EqualBalanced and Weighted

PlacementResult place_pooled(PolicyEngine& pe, std::optional<Amount> v, const OutputOptions& outputs, Trace* trace,
                             const std::function<void()>* staged) {
  const auto& g = *pe.graph;
  Working W(pe.fg);
  Call c {pe, g, W};
  c.trace = trace;
  const std::size_t k = pe.roster.size();
  const bool lossy = pe.config.admission == Admission::Lossy;
  const bool created = pe.cohorts.front().empty();
  Prepared prep;
  std::vector<LiveFlow> flows(k);
  if (created) {
    prep = prepare_cost_only(c);
    flows = prep.flows;
  } else {
    for (std::size_t i = 0; i < k; ++i) flows[i] = pe.cohorts[i].front();
  }
  const auto program = created ? prep.program : pe.program;
  const KernelPlan& plan = *program;
  const auto members = std::span<const EdgeId>(plan.members);
  std::vector<const FlowRecord*> before(k, nullptr);
  if (!created) {
    for (std::size_t i = 0; i < k; ++i) before[i] = &LedgerAccess::live(pe.fg, flows[i].id);
  }
  std::vector<Draft> drafts(k);
  for (std::size_t i = 0; i < k; ++i) {
    if (!created) drafts[i].id = flows[i].id;
    drafts[i].record = source_record(pe, flows[i]);
  }
  // A one-source program is its own attribution.
  Attribution attr;
  Attribution* per_source = k > 1 ? &attr : nullptr;
  std::vector<Amount> inject(k, 0);
  Forwarding fw;
  Amount placed = 0;
  Counter128 offered = pe.offered;
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
      cumulative_offers(pe.offered, *v, k, [&](std::size_t i, Amount x) { inject[i] = x; });
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
  Counter128 requested = pe.requested;
  if (v.has_value()) (void)requested.add(*v);  // checked in place_impl
  PlacementResult out;
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
  out.routes_changed = created;
  fill_state(pe, draft_records(drafts), requested, pe.infinite || !v.has_value(), offered, out);
  if (kTraced && trace) {
    trace->set("prepare", created ? routes_trace(prep) : Trace::null());
    trace->set(lossy ? "offers" : "T", std::move(step_tr));
  }
  if (any_output(outputs)) prepare_outputs(pe, outputs, draft_edges(drafts), out);
  if (staged != nullptr && *staged) (*staged)();
  auto next = commit_sources(pe, drafts, flows);
  pe.cohorts.swap(next);
  pe.program = program;
  pe.requested = requested;
  pe.offered = offered;
  pe.prepared = true;
  if (!v.has_value()) pe.infinite = true;
  return out;
}

// ---------------------------------------------------------------- the controllers

// One solve of a controller from zero on the view V with value cap `cap` (nullopt: uncapped);
// `caps` (the equal origin) caps source i's arc at caps[i]. Global solves over every eligible
// edge, its V and flows per EdgeId; CostOnly MaxFlow over `support`, its V and flows aligned with
// the support.
struct Solve {
  Amount placed {0};
  std::vector<Amount> flows {};
  FlowCoreResult payload {};  // Global: the flow core's payload (the curve, the cut, the trace's values)
};

Solve solve(const PolicyEngine& pe, std::span<const EdgeId> support, std::span<const Amount> V,
            std::optional<Amount> cap, std::span<const Amount> caps) {
  const auto& g = *pe.graph;
  Solve out;
  if (is<Global>(pe)) {
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
// (Global: the eligible ones, V per EdgeId; CostOnly MaxFlow: the support's, V aligned with it);
// the supplies leave through them.
Amount roster_out(const PolicyEngine& pe, std::span<const EdgeId> support, std::span<const Amount> V) {
  const auto row = pe.graph->row_offsets();
  Amount sum = 0;
  for (NodeId s : pe.roster) {
    const auto si = static_cast<std::size_t>(s);
    if (is<Global>(pe)) {
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

PlacementResult place_controller(PolicyEngine& pe, std::optional<Amount> v, const OutputOptions& outputs, Trace* trace,
                                 const std::function<void()>* staged) {
  const auto& g = *pe.graph;
  const bool global = is<Global>(pe);
  const bool equal = pe.origin == Origin::Equal;
  const std::size_t k = pe.roster.size();
  const bool created = pe.cohorts.front().empty();
  // CostOnly MaxFlow: the support of the epoch (its routes, prepared at its first call).
  Prepared prep;
  if (!global && created) {
    Working W(pe.fg);
    Call c {pe, g, W};
    prep = prepare_cost_only(c);
  }
  const std::span<const EdgeId> support = global    ? std::span<const EdgeId> {}
                                          : created ? std::span<const EdgeId>(prep.edges)
                                                    : std::span<const EdgeId>(pe.program_edges);
  std::vector<LiveFlow> flows(k);
  for (std::size_t i = 0; i < k; ++i) {
    if (!created) {
      flows[i] = pe.cohorts[i].front();
    } else if (!global) {
      flows[i] = prep.flows[i];
    } else {
      flows[i].src = pe.roster[i];
    }
  }
  // V = W + own: the residual plus the policy's own allocation, per EdgeId (Global) or aligned
  // with the support (CostOnly MaxFlow: a record's edges lie in its source's route, inside the
  // support, both ascending).
  std::vector<Amount> V;
  if (global) {
    LedgerAccess::residual_into(pe.fg, V);
  } else {
    Working(pe.fg).gather(support, V);
  }
  Amount P = 0;
  for (const auto& f : flows) {
    if (f.id < 0) continue;
    const FlowRecord& rec = LedgerAccess::live(pe.fg, f.id);
    P += rec.delivered;
    std::size_t p = 0;
    for (const auto& [e, x] : rec.edges) {
      if (global) {
        V[static_cast<std::size_t>(e)] += x;
        continue;
      }
      while (p < support.size() && support[p] < e) ++p;
      if (p == support.size() || support[p] != e) throw std::logic_error("CostOnly MaxFlow: a record edge outside the support");
      V[p] += x;
    }
  }
  const std::optional<Amount> cap = v.has_value() ? std::optional<Amount>(P + *v) : std::nullopt;
  if (pe.capture != nullptr) {
    if (equal || !global) throw std::logic_error("a deletion capture of a controller other than optimized Global");
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
  // The records: each source's share of the flow (decompose(); one source's is the whole flow).
  std::vector<std::vector<std::pair<EdgeId, Amount>>> shares;
  // CostOnly MaxFlow's flows per EdgeId where decompose() or the trace reads them.
  std::vector<Amount> dense;
  if (!global && (k > 1 || (kTraced && trace))) {
    dense.assign(static_cast<std::size_t>(g.num_edges()), 0);
    for (std::size_t p = 0; p < support.size(); ++p) dense[static_cast<std::size_t>(support[p])] = sol.flows[p];
  }
  if (k == 1) {
    if (global) {
      shares.emplace_back();
      for (std::size_t e = 0; e < sol.flows.size(); ++e) {
        if (sol.flows[e] != 0) shares.front().emplace_back(static_cast<EdgeId>(e), sol.flows[e]);
      }
    } else {
      shares.push_back(positive(support, sol.flows));
    }
  } else {
    shares = decompose(g, global ? sol.flows : dense, pe.roster, pe.dests_ok);
  }
  bool changed = created;
  std::vector<Amount> supplies(k, 0);
  std::vector<Draft> drafts(k);
  const auto row = g.row_offsets();
  for (std::size_t i = 0; i < k; ++i) {
    const NodeId s = pe.roster[i];
    auto& f = flows[i];
    // The supply of s: what its share delivers (its out-edges' share).
    Amount supply = 0;
    for (const auto& [e, x] : shares[i]) {
      if (e >= row[static_cast<std::size_t>(s)] && e < row[static_cast<std::size_t>(s) + 1]) supply += x;
    }
    supplies[i] = supply;
    if (global) {
      // Global's record bundle is the edges its share uses.
      std::vector<EdgeId> route;
      route.reserve(shares[i].size());
      for (const auto& [e, x] : shares[i]) route.push_back(e);
      const auto old = f.bundle ? PredDAGBuilder::ascending(*f.bundle) : std::span<const EdgeId> {};
      if (!f.bundle || !std::equal(route.begin(), route.end(), old.begin(), old.end())) {
        if (f.bundle) changed = true;
        f.bundle = std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, route));
        f.cost = route.empty() ? kUnreachable : route_cost(g, *f.bundle, s, pe.dests_ok);
      }
    }
    if (f.id >= 0) drafts[i].id = f.id;
    drafts[i].record = source_record(pe, f);
    drafts[i].record.delivered = supply;
    drafts[i].record.edges = shares[i];
  }
  Counter128 requested = pe.requested;
  if (v.has_value()) (void)requested.add(*v);  // checked in place_impl
  PlacementResult out;
  out.placed = sol.placed - P;
  if (v.has_value()) out.remaining = *v - out.placed;
  out.delivery_delta = out.placed;
  out.termination = v.has_value() && out.placed == *v ? Termination::Complete : Termination::Exhausted;
  out.routes_changed = changed;
  const auto after = draft_records(drafts);
  fill_state(pe, after, requested, pe.infinite || !v.has_value(), pe.offered, out);
  out.total_cost = total_cost_of(g, after);
  const bool curve = global && !equal;
  if (curve) out.cost_distribution = sol.payload.cost_distribution;
  if (curve && outputs.min_cut) out.min_cut = sol.payload.min_cut;
  if (kTraced && trace) {
    if (!global) trace->set("prepare", created ? routes_trace(prep) : Trace::null());
    if (equal) trace->set("probes", std::move(probes));
    if (global) trace->set("V", Trace::ints(V));
    trace->set("cap", cap.has_value() ? Trace::integer(*cap) : Trace::null());
    Trace fc = Trace::map();
    fc.set("placed", Trace::integer(sol.placed));
    if (global) {
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
      fc.set("edge_flows", Trace::ints(dense));
    }
    Trace supply = Trace::map();
    Trace shares_tr = Trace::map();
    for (std::size_t i = 0; i < k; ++i) {
      supply.set(key(pe.roster[i]), Trace::integer(supplies[i]));
      shares_tr.set(key(pe.roster[i]), Trace::pairs(shares[i]));
    }
    fc.set("supply", std::move(supply));
    trace->set("flow_core", std::move(fc));
    trace->set("shares", std::move(shares_tr));
  }
  if (any_output(outputs)) prepare_outputs(pe, outputs, draft_edges(drafts), out);
  std::vector<std::pair<Cost, Amount>> next_curve;
  if (curve) next_curve = *out.cost_distribution;  // the cumulative curve a later zero call reports
  std::vector<EdgeId> next_edges;
  if (!global && created) next_edges = prep.edges;
  if (staged != nullptr && *staged) (*staged)();
  auto next = commit_sources(pe, drafts, flows);
  pe.cohorts.swap(next);
  pe.curve.swap(next_curve);
  if (!global && created) pe.program_edges.swap(next_edges);
  pe.requested = requested;
  pe.prepared = true;
  if (!v.has_value()) pe.infinite = true;
  return out;
}

// ---------------------------------------------------------------- construction

} // namespace

const char* model_name(const Routing& r) noexcept {
  if (std::holds_alternative<CostOnly>(r)) return "CostOnly";
  if (std::holds_alternative<Global>(r)) return "Global";
  if (std::holds_alternative<CapacityAware>(r)) return "T1";
  return "Pinned";
}

ResolvedAxes resolve_config(const FlowPolicyConfig& cfg, bool cache) {
  const auto& r = cfg.routing;
  if (r.valueless_by_exception()) invalid("routing: the variant holds no routing mode");
  // The enumeration values.
  if (cfg.placement.has_value()) check_member(*cfg.placement, "placement");
  check_member(cfg.admission, "admission");
  if (cfg.origin.has_value()) check_member(*cfg.origin, "origin");
  if (cfg.selection.has_value()) check_selection(*cfg.selection);
  // The model's row, its columns left to right: placement, origin, admission, selection.
  ResolvedAxes out {cfg.selection, cfg.origin};
  const bool cost_only = std::holds_alternative<CostOnly>(r);
  const bool global = std::holds_alternative<Global>(r);
  const bool capacity_aware = std::holds_alternative<CapacityAware>(r);
  if (global) {
    if (cfg.placement.has_value()) invalid("Global does not accept placement");
  } else if (!cfg.placement.has_value()) {
    invalid(std::string(cost_only ? "CostOnly" : capacity_aware ? "CapacityAware" : "Pinned") + " requires placement");
  }
  const bool max_flow = global || *cfg.placement == FlowPlacement::MaxFlow;
  if (cost_only && !max_flow && cfg.origin == Origin::Optimized) {
    invalid("CostOnly with EQUAL_BALANCED or WEIGHTED requires origin=EQUAL");
  }
  if (capacity_aware && cfg.origin == Origin::Optimized) invalid("T1 requires origin=EQUAL");
  if (std::holds_alternative<Pinned>(r) && cfg.origin.has_value()) invalid("Pinned does not accept origin");
  if (!cfg.origin.has_value() && !std::holds_alternative<Pinned>(r)) {
    out.origin = cost_only && !max_flow ? Origin::Equal : capacity_aware ? Origin::Equal : Origin::Optimized;
  }
  if (cfg.admission == Admission::Lossy) {
    if (max_flow) invalid("MAX_FLOW requires admission=LOSSLESS");
    if (capacity_aware) invalid("CapacityAware requires admission=LOSSLESS");
  }
  if (global) {
    if (cfg.selection.has_value()) invalid("Global does not accept selection");
  } else {
    if (!cfg.selection.has_value()) out.selection = PathSelection();
    check_kernel_combination(*cfg.placement, cfg.admission, out.selection->priority);
  }
  if (cache && !cost_only) invalid(std::string(model_name(r)) + " does not use DagCache");
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
  return out;
}

void check_source_count(const FlowPolicyConfig& cfg, std::size_t num_sources, std::size_t num_destinations) {
  const auto& r = cfg.routing;
  if (std::holds_alternative<Pinned>(r) && (num_sources != 1 || num_destinations != 1)) {
    invalid("Pinned routing requires single endpoints (bundles bind one source and one destination)");
  }
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

// The engine of a new policy, bound to its ledger as the last step. A supplied frozen cache's hit
// is recorded in `tally`: the public constructor publishes it once construction has
// succeeded, an evaluator once its whole call has.
std::unique_ptr<PolicyEngine> make_engine(FlowLedger& ledger, std::span<const NodeId> sources,
                                          std::span<const NodeId> destinations, FlowPolicyConfig config,
                                          const Eligibility* elig, DagCache* cache, detail::FrozenTally& tally) {
  auto pe = std::make_unique<PolicyEngine>(ledger, std::move(config));
  pe->graph = ledger.graph();
  const auto& g = *pe->graph;
  // The endpoints, then the ownership of the mask and the cache, then the configuration.
  pe->sources = detail::node_set(sources, g.num_nodes(), "sources");
  pe->destinations = detail::node_set(destinations, g.num_nodes(), "destinations");
  detail::check_disjoint(pe->sources, pe->destinations, "sources", "destinations");
  if (elig != nullptr && elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  if (cache != nullptr && cache->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign DagCache: bound to another graph");
  }
  pe->elig = elig;
  const auto axes = detail::resolve_config(pe->config, cache != nullptr);
  pe->config.selection = axes.selection;
  pe->config.origin = axes.origin;
  detail::check_source_count(pe->config, pe->sources.size(), pe->destinations.size());
  if (pe->config.selection.has_value()) pe->selection = *pe->config.selection;
  if (pe->config.origin.has_value()) pe->origin = *pe->config.origin;
  auto node_ok = [&](NodeId v) { return elig == nullptr || elig->node_ok(v); };
  for (NodeId d : pe->destinations) {
    if (node_ok(d)) pe->dests_ok.push_back(d);
  }
  if (detail::is<Pinned>(*pe)) {
    const auto& bundles = std::get<Pinned>(pe->config.routing).bundles;
    if (bundles.empty()) throw std::invalid_argument("Pinned: an empty bundle list is rejected");
    const NodeId src = pe->sources.front();
    const NodeId dst[1] = {pe->destinations.front()};
    for (std::size_t i = 0; i < bundles.size(); ++i) {
      const auto& b = bundles[i];
      if (!b) throw std::invalid_argument("Pinned: bundle " + std::to_string(i) + " is null");
      if (b->graph_uid() != g.uid()) throw std::invalid_argument("foreign PredDAG: bundle built for another graph");
      for (std::size_t j = 0; j < i; ++j) {
        if (*bundles[j] == *b) {
          throw std::invalid_argument("Pinned: duplicate bundle " + std::to_string(i) + " (equal to bundle " +
                                      std::to_string(j) + " before pruning)");
        }
      }
      Cost cost = kUnreachable;
      auto plan = detail::make_route_plan(g, *pe->config.placement, pe->config.admission, pe->selection.priority, src,
                                          dst, *b, elig, &cost);
      if (plan.route.empty() && detail::prune(g, *b, src, dst).empty()) {
        throw std::invalid_argument("Pinned: malformed bundle " + std::to_string(i) + ", no src -> dst walk (src " +
                                    std::to_string(src) + ", dst " + std::to_string(dst[0]) + ")");
      }
      detail::PreparedPinnedBundle pb;
      pb.cost = cost;
      pb.route_dag = detail::bundle_of(g, plan.route, b);
      pb.route_plan = std::make_shared<const detail::KernelPlan>(std::move(plan));
      if (!pb.route_plan->route.empty()) pe->usable.push_back(static_cast<std::int32_t>(i));
      pe->bundles.push_back(std::move(pb));
    }
    if (!pe->usable.empty()) pe->roster.push_back(src);
  } else if (detail::is<CostOnly>(*pe)) {
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
  }
  pe->cohorts.resize(pe->roster.size());
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

namespace {

// `staged` (optional) runs once the call's result and trace are final and before anything is
// published (the ledger commit, the epoch state): an exception from it aborts the call with
// nothing changed (the Python binding builds its result objects there).
PlacementResult place_impl(PolicyEngine& pe, std::optional<Amount> volume, const OutputOptions& outputs,
                           detail::Trace* trace, const std::function<void()>* staged) {
  // Validation before any mutation, also before the zero shortcut.
  if (volume.has_value() && (*volume < 0 || *volume >= kDomainLimit)) {
    throw std::invalid_argument("offer must be an integer in [0, 2^62) or permitted +inf");
  }
  if (pe.infinite && volume.has_value()) {
    throw std::invalid_argument("finite placement after infinite placement requires a new epoch");
  }
  // The only check of the cumulative request: nothing changes the requested total between here
  // and the call's publication, so the later additions (in the engines) repeat this sum and
  // cannot fail.
  if (volume.has_value()) {
    Counter128 probe = pe.requested;
    if (!probe.add(*volume)) throw std::invalid_argument("requested limit exceeded");
  }
  const bool lossy = pe.config.admission == Admission::Lossy;
  const bool optimized_global = detail::is<Global>(pe) && pe.origin == Origin::Optimized;
  PlacementResult out;
  if (volume.has_value() && *volume == 0) {
    // The current state with every flag false; before the epoch's first positive call it is all
    // zeros (no records, no curve). A zero call solves nothing, so it claims no cut.
    out.remaining = 0;
    const auto records = detail::current_records(pe);
    detail::fill_state(pe, records, pe.requested, pe.infinite, pe.offered, out);
    if (detail::controller(pe)) out.total_cost = detail::total_cost_of(*pe.graph, records);
    if (optimized_global && pe.prepared) out.cost_distribution = pe.curve;
    detail::prepare_outputs(pe, outputs, detail::current_edges(pe), out);
    if (staged != nullptr && *staged) (*staged)();
    return out;
  }
  if (pe.roster.empty()) {
    // No source reaches an allowed destination. A Lossy offer is lost before it enters (the
    // demand's source loss; no source's row has it). Global's flow is 0 and maximal: its curve
    // is empty and its cut is empty (no roster source, so the closure is empty).
    Counter128 requested = pe.requested;
    Counter128 offered = pe.offered;
    if (volume.has_value()) {
      (void)requested.add(*volume);              // checked above
      if (lossy) (void)offered.add(*volume);     // offered <= requested
    }
    out.remaining = volume;
    out.termination = Termination::NoRoute;
    if (lossy && volume.has_value()) out.call_source_loss = *volume;
    detail::fill_state(pe, {}, requested, pe.infinite || !volume.has_value(), offered, out);
    if (detail::controller(pe)) out.total_cost = Counter128 {};
    if (optimized_global) {
      out.cost_distribution.emplace();
      if (outputs.min_cut) out.min_cut.emplace();
    }
    if (detail::kTraced && trace) trace->set("roster", detail::Trace::list());
    detail::prepare_outputs(pe, outputs, detail::current_edges(pe), out);
    if (staged != nullptr && *staged) (*staged)();
    // Nothing below throws.
    pe.requested = requested;
    pe.offered = offered;
    pe.prepared = true;
    if (!volume.has_value()) pe.infinite = true;
    return out;
  }
  if (detail::controller(pe)) return detail::place_controller(pe, volume, outputs, trace, staged);
  if (detail::is<CostOnly>(pe)) return detail::place_pooled(pe, volume, outputs, trace, staged);
  return detail::place_flows(pe, volume, outputs, trace, staged);
}

} // namespace

PlacementResult FlowPolicy::place(std::optional<Amount> volume, const OutputOptions& outputs) {
  return place_impl(*engine_, volume, outputs, nullptr, nullptr);
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
PlacementResult detail::place_traced(FlowPolicy& policy, std::optional<Amount> volume, const OutputOptions& outputs,
                                     Trace& trace, const std::function<void()>* staged) {
  trace = Trace::map();
  return place_impl(PolicyAccess::engine(policy), volume, outputs, &trace, staged);
}
#endif

void FlowPolicy::remove() {
  auto& pe = *engine_;
  (void)LedgerAccess::detach(pe.fg, pe.owner);
  // A new epoch (nothing below throws).
  for (auto& c : pe.cohorts) c.clear();
  pe.program.reset();
  pe.program_edges.clear();
  pe.curve.clear();
  pe.requested = Counter128 {};
  pe.offered = Counter128 {};
  pe.infinite = false;
  pe.prepared = false;
}

Amount FlowPolicy::placed_demand() const {
  const auto& pe = *engine_;
  Amount total = 0;
  for (const auto& c : pe.cohorts) {
    for (const auto& s : c) total += LedgerAccess::live(pe.fg, s.id).delivered;
  }
  return total;
}

std::vector<FlowRecord> FlowPolicy::records() const {
  const auto& pe = *engine_;
  std::vector<std::int64_t> ids;
  for (const auto& c : pe.cohorts) {
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

std::int32_t FlowPolicy::usable_bundles() const noexcept {
  return detail::is<Pinned>(*engine_) ? static_cast<std::int32_t>(engine_->usable.size()) : 0;
}

const FlowPolicyConfig& FlowPolicy::config() const noexcept { return engine_->config; }

} // namespace netgraph::core
