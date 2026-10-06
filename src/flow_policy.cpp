/*
  The policy engine FlowPolicy (APPROACH 8 and 9): one placement transaction for every routing
  mode, committed through the ledger as one atomic exchange per call.

  * CostOnly, CapacityAware and Pinned run the per-source transaction of 9.2 and 9.3: the
    roster (fixed at construction) in ascending NodeId, each source offered the increment of its
    round-robin share of the cumulative request, stages 1 (creation), 2 (reoptimisation), 3
    (amounts) and 4 (growth) on the shared working residual W; one `created` counter per call.
    Pinned has one actor, the policy, with its one source.
  * Global re-solves the minimum-cost flow of the cumulative value on V = W + own (the flow core
    of 2.2 over the whole eligible graph) and decomposes it per source canonically (9.1).

  Every stage works on copies; nothing is published before the commit, which either succeeds or
  throws leaving the ledger and the policy unchanged; the epoch state follows the commit with
  operations that cannot throw. The working residual is a per-thread overlay over the ledger's
  residual (only the edges a call touches cost anything), so a Pinned call, and a CostOnly call
  after its cohort's creation, costs O(its flows' members), never O(E); CostOnly creation of a
  single path, CapacityAware searches and Global solves read a dense copy.

  The transaction trace (prototype 2's keys) exists in a build with the private hooks only: every
  trace block is guarded by kTraced and is compile-time dead without them.
*/
#include "netgraph/core/flow_policy.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/flow_core.hpp"
#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/node_map.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/detail/search.hpp"
#include "netgraph/core/detail/trace.hpp"

namespace netgraph::core {

namespace detail {

#if NETGRAPH_CORE_PRIVATE_HOOKS
std::atomic<std::int64_t>& max_flows_per_call() noexcept {
  static std::atomic<std::int64_t> cap {kMaxFlowsPerCall};
  return cap;
}

std::atomic<std::int64_t>& sensitivity_fail_edge() noexcept {
  static std::atomic<std::int64_t> edge {-1};
  return edge;
}
#endif

namespace {

// ================================================================ working residual

// W = baseline - allocated + delta, the delta being this call's working changes (releases of the
// policy's own allocations and new working allocations). Per-thread slots stamped with the
// call's epoch: a slot is valid only under the current epoch, so nothing is cleared between calls
// and a call that throws leaves stale slots behind an old epoch. One array of slots: growing it
// either succeeds whole or leaves the previous array (strong guarantee).
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
  // The member-aligned working residual of `members`.
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
  }
  return "?";
}

// The reduction rank of 9.4 (higher wins among unserved sources).
int rank(Termination t) noexcept {
  switch (t) {
    case Termination::LimitFlows: return 4;
    case Termination::LimitDiminishing: return 3;
    case Termination::Exhausted: return 2;
    case Termination::NoRoute: return 1;
    case Termination::Complete: return 0;
  }
  return 0;
}

Trace amount_or_inf(std::optional<Amount> v) { return v.has_value() ? Trace::integer(*v) : Trace::str("inf"); }

// cost <= f * base, exactly: f is the dyadic rational it is (f >= 1, finite), base >= 1.
// f = m * 2^x with m odd below 2^53; products stay within 128 bits (APPROACH 6).
bool stretch_le(Cost cost, Cost base, double f) {
  const auto [m, x] = dyadic(f);
  const U128 rhs = mul_u64(static_cast<std::uint64_t>(base), m);  // below 2^115
  if (x >= 0) {
    if (bits(rhs) + x > 64) return true;  // rhs * 2^x >= 2^64 > cost
    return !(shl(rhs, x) < U128 {0, static_cast<std::uint64_t>(cost)});
  }
  // f >= 1 with m < 2^53: -x < 53, so cost * 2^-x < 2^116.
  const U128 lhs = shl(U128 {0, static_cast<std::uint64_t>(cost)}, -x);
  return !(rhs < lhs);
}

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

// Round robin of a cumulative 128-bit request: the per-source increment of the share when the
// cumulative request goes from R to R + v over k sources (APPROACH 9.3). It depends only on
// r = R mod k: floor((r + v) / k) + [i < (r + v) mod k] - [i < r].
std::vector<Amount> cumulative_offers(const Counter128& requested, Amount v, std::size_t k) {
  const auto kk = static_cast<std::uint64_t>(k);
  const std::uint64_t two64 = ((std::numeric_limits<std::uint64_t>::max() % kk) + 1) % kk;
  const std::uint64_t r = ((requested.hi % kk) * two64 % kk + requested.lo % kk) % kk;
  const auto total = static_cast<std::uint64_t>(v) + r;
  const auto base = static_cast<Amount>(total / kk);
  const std::uint64_t extra = total % kk;
  std::vector<Amount> out(k);
  for (std::size_t i = 0; i < k; ++i) {
    out[i] = base + (i < extra ? 1 : 0) - (i < r ? 1 : 0);
  }
  return out;
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
  const auto& cfg = c.pe.config;
  std::vector<Amount> creation;
  if (cfg.placement == FlowPlacement::Weighted && cfg.selection.priority == EdgePriority::Residual) {
    const auto edges = PredDAGBuilder::ascending(route);
    creation.resize(edges.size());
    for (std::size_t i = 0; i < edges.size(); ++i) {
      creation[i] = dense_view != nullptr ? (*dense_view)[static_cast<std::size_t>(edges[i])] : c.W.at(edges[i]);
    }
  }
  return std::make_shared<const KernelPlan>(
      make_plan(c.g, cfg.placement, cfg.admission, cfg.selection.priority, src, c.pe.dests_ok, route, creation, cost));
}

// The record bundle of a route: `share` itself when it is the whole route, else a new DAG.
std::shared_ptr<const PredDAG> bundle_of(const Graph& g, std::span<const EdgeId> route,
                                         const std::shared_ptr<const PredDAG>& share) {
  if (share && PredDAGBuilder::ascending(*share).size() == route.size()) return share;
  return std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, route));
}

// CostOnly: the static support restricted to the walks from src (a forward reachability over the
// support's edges: B' ids rank edges by source, so a node's support out-edges are one range of
// the ascending edge list), then the selection projected on the creation view W (APPROACH 9.1,
// 10): Collapse by (cost, priority under PreferHigherPriority, lowest EdgeId), and the staged
// canonical single path for multipath = false.
Route cost_only_route(const Call& c, NodeId src) {
  const auto& g = c.g;
  const auto& pe = c.pe;
  const auto& sel = pe.config.selection;
  const auto support = PredDAGBuilder::ascending(*pe.support);
  const auto row = g.row_offsets();
  const auto edge_dst = g.edge_dst();
  const auto edge_src = g.edge_src();
  std::vector<EdgeId> route;
  {
    // Visited nodes in the per-thread node map (released before any routine that uses it).
    NodeMap& map = NodeMap::local();
    map.reset(static_cast<std::size_t>(g.num_nodes()));
    std::vector<NodeId> stack {src};
    (void)map.add(src);
    while (!stack.empty()) {
      const NodeId u = stack.back();
      stack.pop_back();
      auto it = std::lower_bound(support.begin(), support.end(), row[static_cast<std::size_t>(u)]);
      for (; it != support.end() && *it < row[static_cast<std::size_t>(u) + 1]; ++it) {
        route.push_back(*it);
        const NodeId v = edge_dst[static_cast<std::size_t>(*it)];
        if (map.of(v) < 0) {
          (void)map.add(v);
          stack.push_back(v);
        }
      }
    }
    std::sort(route.begin(), route.end());
  }
  Route out;
  if (!sel.multi_edge) {
    // One member per (u, v) by select_adjacency on the creation view W (the support has no
    // residual gating, so the projection reads W's priority values, never gates on them): the
    // members of a pair are consecutive in id order (B' ranks by (src, cost, dst, input)).
    const bool prefer = sel.tie_break == EdgeTieBreak::PreferHigherPriority;
    const auto capacity = g.capacity_units();
    auto value = [&](EdgeId e) { return priority_value(sel.priority, c.W.at(e), capacity[static_cast<std::size_t>(e)]); };
    std::vector<EdgeId> kept;
    for (std::size_t i = 0; i < route.size();) {
      const auto ui = static_cast<std::size_t>(route[i]);
      EdgeId best = route[i];
      std::size_t j = i + 1;
      for (; j < route.size() && edge_src[static_cast<std::size_t>(route[j])] == edge_src[ui] &&
             edge_dst[static_cast<std::size_t>(route[j])] == edge_dst[ui];
           ++j) {
        if (replaces(g, route[j], best, prefer, value)) best = route[j];
      }
      kept.push_back(best);
      i = j;
    }
    route = std::move(kept);
  }
  if (sel.multipath) {
    out.edges = std::move(route);
    return out;
  }
  // The staged canonical single path over the (collapsed) support of the source.
  std::vector<std::uint8_t> member(static_cast<std::size_t>(g.num_edges()), 0);
  for (EdgeId e : route) member[static_cast<std::size_t>(e)] = 1;
  std::vector<Amount> wd;
  std::optional<PriorityValues> prio;
  if (sel.tie_break == EdgeTieBreak::PreferHigherPriority) {
    c.W.dense(wd);
    prio = PriorityValues {&g, sel.priority, wd};
  }
  const NodeId s[1] = {src};
  auto path = canonical_path(g, member, s, pe.dests_ok, prio ? &*prio : nullptr);
  if (path.has_value()) {
    out.edges = path->walk;
    std::sort(out.edges.begin(), out.edges.end());
    out.cost = path->cost;
  }
  return out;
}

// CapacityAware: the canonical route of section 6 on the selection view `view` (dense): static
// eligibility, residual >= max(1, min_flow_bandwidth), Collapse on the view, the tight DAG or the
// staged single path toward the allowed destinations, then the max_path_cost ceiling.
std::optional<Route> search(const Call& c, std::size_t ri, std::span<const Amount> view) {
  const auto& pe = c.pe;
  const auto& g = c.g;
  const auto& ca = std::get<CapacityAware>(pe.config.routing);
  const auto& sel = pe.config.selection;
  const Amount threshold = std::max<Amount>(1, ca.min_flow_bandwidth);
  const auto member = member_mask(g, pe.elig, sel, view, threshold);
  const NodeId s[1] = {pe.roster[ri]};
  Route out;
  if (sel.multipath) {
    auto tight = tight_dag(g, member, s, pe.dests_ok);
    if (tight.cost == kUnreachable) return std::nullopt;
    out.edges = std::move(tight.edges);
    out.cost = tight.cost;
  } else {
    const PriorityValues prio {&g, sel.priority, view};
    const bool prefer = sel.tie_break == EdgeTieBreak::PreferHigherPriority;
    auto path = canonical_path(g, member, s, pe.dests_ok, prefer ? &prio : nullptr);
    if (!path.has_value()) return std::nullopt;
    out.edges = std::move(path->walk);
    std::sort(out.edges.begin(), out.edges.end());
    out.cost = path->cost;
  }
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

// The cohort allocation of an EqualBalanced cohort (APPROACH 16.5): T split by round robin over
// the flows in order, each flow's own plan, summed per edge over the union of members.
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
    // The search's upper bound: the view summed over the source members of every flow.
    std::vector<std::uint8_t> src_member(uni_.size(), 0);
    for (std::size_t i = 0; i < plans_.size(); ++i) {
      for (auto k = plans_[i]->src_out_begin; k < plans_[i]->src_out_end; ++k) {
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
// allocation of its new share (APPROACH 9.2 stage 3 EB).
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
  if (is<CostOnly>(pe)) {
    const Route route = cost_only_route(c, st.src);
    if (kTraced && c.trace) tr.set("route", route_trace(route.edges));
    WorkFlow f = new_flow(c, st.src, bundle_of(c.g, route.edges, pe.support), route.cost, nullptr);
    f.seq = ++c.seq;
    st.flows.push_back(std::move(f));
  } else if (is<Pinned>(pe)) {
    const bool captured = pe.config.placement == FlowPlacement::Weighted &&
                          pe.config.selection.priority == EdgePriority::Residual;
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
    if (pe.config.placement == FlowPlacement::EqualBalanced && st.offer.has_value()) {
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
      std::optional<Amount> request = pe.config.placement == FlowPlacement::EqualBalanced ? target : rem;
      auto step = lossless_step(c.g, *f.plan, 0, view_m, zeros, request);
      const Amount t = step.total;
      if (pe.config.placement != FlowPlacement::EqualBalanced && rem.has_value()) *rem -= t;
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
    if (c.pe.config.placement == FlowPlacement::Weighted) c.weights_changed = true;
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
    const auto step = lossy_step(*f.plan, offer, w_m);
    c.W.consume(members, step.carried);
    for (std::size_t m = 0; m < members.size(); ++m) {
      f.alloc[m] += step.carried[m];
      if (step.drops[m] != 0) {
        if (!f.drops[m].add(step.drops[m])) throw std::invalid_argument("a cumulative drop counter would reach 2^127");
        c.drops.emplace_back(members[m], step.drops[m]);
      }
    }
    if (!f.offered->add(step.offered)) throw std::invalid_argument("the cumulative offered counter would reach 2^127");
    f.delivered += step.delivered;
    st.injected += step.injected;
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

// fraction * v as the reduced rational "num/den" (or "num"), prototype 2's trace text.
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

// gains.sum < fraction * v, exactly (fraction the dyadic rational it is, v < 2^62).
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
    if (pe.config.selection.multipath) {
      // Revision 24 (owner): a multipath DAG already spans every equal-cost path, so a second
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
    if (pe.config.placement == FlowPlacement::EqualBalanced) {
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
  } else if (pe.config.placement == FlowPlacement::EqualBalanced) {
    stage3_eb(c, st);
  } else {
    stage3_seq(c, st);
  }
  if (is<CapacityAware>(pe)) stage4(c, st);
}

Termination local_outcome(const PolicyEngine& pe, const SourceRun& st) {
  Amount after = 0;
  for (const auto& f : st.flows) after += pe.config.admission == Admission::Lossy ? f.delivered : f.T;
  const Amount delivered = after - st.P;
  if (st.offer.has_value() && *st.offer == 0) return Termination::Complete;
  if (pe.config.admission == Admission::Lossy) {
    if (st.offer.has_value() && st.injected == *st.offer) return Termination::Complete;
  } else if (st.offer.has_value() && delivered == *st.offer) {
    return Termination::Complete;
  }
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
  r.placement = pe.config.placement;
  r.admission = pe.config.admission;
  r.bundle = f.bundle;
  r.route_cost = f.cost;
  const auto members = f.members();
  r.edges.reserve(static_cast<std::size_t>(std::count_if(f.alloc.begin(), f.alloc.end(), [](Amount x) { return x != 0; })));
  for (std::size_t m = 0; m < members.size(); ++m) {
    if (f.alloc[m] != 0) r.edges.emplace_back(members[m], f.alloc[m]);
  }
  if (pe.config.admission == Admission::Lossy) {
    r.offered = f.offered;
    r.delivered = f.delivered;
    r.drops.reserve(static_cast<std::size_t>(
        std::count_if(f.drops.begin(), f.drops.end(), [](const Counter128& x) { return !x.is_zero(); })));
    for (std::size_t m = 0; m < members.size(); ++m) {
      if (!f.drops[m].is_zero()) r.drops.emplace_back(members[m], f.drops[m]);
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
    const auto inc = cumulative_offers(pe.requested, *v, k);
    for (std::size_t i = 0; i < k; ++i) runs[i].offer = inc[i];
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
  Counter128 requested = pe.requested;
  if (v.has_value() && !requested.add(*v)) throw std::invalid_argument("the cumulative request counter would reach 2^127");
  auto next = next_cohorts(runs, k);
  Amount placed = 0;
  std::vector<std::pair<NodeId, Termination>> outcomes;
  outcomes.reserve(k);
  Termination term = Termination::Complete;
  for (const auto& st : runs) {
    Amount after = 0;
    for (const auto& f : st.flows) after += pe.config.admission == Admission::Lossy ? f.delivered : f.T;
    placed += after - st.P;
    const Termination o = local_outcome(pe, st);
    outcomes.emplace_back(st.src, o);
    if (rank(o) > rank(term)) term = o;
  }
  PlacementResult out;
  out.placed = placed;
  if (v.has_value()) out.remaining = *v - placed;
  out.termination = term;
  out.local_outcomes = std::move(outcomes);
  out.routes_changed = c.routes_changed;
  out.weights_changed = c.weights_changed;
  out.drops = merge_drops(std::move(c.drops));
  if (kTraced && trace) {
    Trace offers_tr = Trace::map();
    Trace outcomes_tr = Trace::map();
    for (std::size_t i = 0; i < k; ++i) {
      offers_tr.set(key(runs[i].src), amount_or_inf(runs[i].offer));
      outcomes_tr.set(key(runs[i].src), Trace::str(outcome_name(out.local_outcomes[i].second)));
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
  if (!v.has_value()) pe.infinite = true;
  return out;
}

// ---------------------------------------------------------------- Global

// The canonical per-source decomposition (9.1): sources ascending; for each, while it has
// supply, the lexicographically smallest EdgeId walk with positive remaining flow from it to the
// first destination it reaches, saturated by min(the walk's bottleneck, the remaining supply).
// Per node the first out-edge with remaining flow only moves forward (remaining flow never grows).
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
  std::vector<Amount> mine;
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

PlacementResult place_global(PolicyEngine& pe, std::optional<Amount> v, const OutputOptions& outputs, Trace* trace,
                             const std::function<void()>* staged) {
  const auto& g = *pe.graph;
  std::vector<Amount> V;
  LedgerAccess::residual_into(pe.fg, V);
  Amount P = 0;
  const bool created = std::all_of(pe.cohorts.begin(), pe.cohorts.end(), [](const auto& c) { return c.empty(); });
  for (std::size_t i = 0; i < pe.roster.size(); ++i) {
    for (const auto& s : pe.cohorts[i]) {
      const FlowRecord& rec = LedgerAccess::live(pe.fg, s.id);
      P += rec.delivered;
      for (const auto& [e, x] : rec.edges) V[static_cast<std::size_t>(e)] += x;
    }
  }
  const std::optional<Amount> cap = v.has_value() ? std::optional<Amount>(P + *v) : std::nullopt;
  if (pe.capture != nullptr) pe.capture->graph = pe.graph;
  auto sol = place_global_units(g, V, pe.roster, pe.dests_ok, pe.elig, cap, pe.capture);
  if (sol.placed < P) throw std::logic_error("Global: a call reduced placed_demand()");
  auto shares = decompose(g, sol.edge_flows, pe.roster, pe.dests_ok);
  bool changed = created;
  std::vector<Amount> supplies(pe.roster.size(), 0);
  std::vector<Draft> drafts;
  std::vector<std::size_t> fresh;
  std::vector<std::vector<LiveFlow>> next(pe.roster.size());
  for (std::size_t i = 0; i < pe.roster.size(); ++i) {
    const NodeId s = pe.roster[i];
    std::vector<EdgeId> route;
    route.reserve(shares[i].size());
    Amount supply = 0;
    for (const auto& [e, x] : shares[i]) route.push_back(e);
    {
      // The supply of s: what its share delivers (its out-edges' share).
      const auto row = g.row_offsets();
      for (const auto& [e, x] : shares[i]) {
        if (e >= row[static_cast<std::size_t>(s)] && e < row[static_cast<std::size_t>(s) + 1]) supply += x;
      }
    }
    std::shared_ptr<const PredDAG> bundle;
    Cost cost = kUnreachable;
    const LiveFlow* old = pe.cohorts[i].empty() ? nullptr : &pe.cohorts[i].front();
    if (old != nullptr && PredDAGBuilder::ascending(*old->bundle).size() == route.size() &&
        std::equal(route.begin(), route.end(), PredDAGBuilder::ascending(*old->bundle).begin())) {
      bundle = old->bundle;
      cost = old->cost;
    } else {
      if (old != nullptr) changed = true;
      bundle = std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, route));
      if (!route.empty()) cost = route_cost(g, *bundle, s, pe.dests_ok);
    }
    Draft d;
    if (old != nullptr) d.id = old->id;
    auto& r = d.record;
    r.sources = {s};
    r.destinations = pe.dests_ok;
    r.placement = FlowPlacement::MaxFlow;
    r.admission = Admission::Lossless;
    r.bundle = bundle;
    r.route_cost = cost;
    r.delivered = supply;
    supplies[i] = supply;
    r.edges = shares[i];
    if (old == nullptr) fresh.push_back(i);
    next[i].push_back(LiveFlow {old != nullptr ? old->id : -1, s, nullptr, bundle, cost, -1});
    drafts.push_back(std::move(d));
  }
  // New records in source order (ascending NodeId), after the retained ones.
  std::stable_partition(drafts.begin(), drafts.end(), [](const Draft& d) { return d.id.has_value(); });
  Counter128 requested = pe.requested;
  if (v.has_value() && !requested.add(*v)) throw std::invalid_argument("the cumulative request counter would reach 2^127");
  PlacementResult out;
  out.placed = sol.placed - P;
  if (v.has_value()) out.remaining = *v - out.placed;
  if (v.has_value() && out.placed == *v) {
    out.termination = Termination::Complete;
  } else if (sol.placed == 0) {
    out.termination = Termination::NoRoute;
  } else {
    out.termination = Termination::Exhausted;
  }
  out.local_outcomes.reserve(pe.roster.size());
  for (NodeId s : pe.roster) out.local_outcomes.emplace_back(s, out.termination);
  out.routes_changed = changed;
  out.cost_distribution = sol.cost_distribution;
  if (kTraced && trace) {
    trace->set("V", Trace::ints(V));
    trace->set("cap", cap.has_value() ? Trace::integer(*cap) : Trace::null());
    Trace fc = Trace::map();
    fc.set("placed", Trace::integer(sol.placed));
    fc.set("cost_distribution", Trace::pairs(sol.cost_distribution));
    fc.set("edge_flows", Trace::ints(sol.edge_flows));
    fc.set("min_cut", sol.min_cut.has_value() ? Trace::ints(*sol.min_cut) : Trace::null());
    fc.set("potentials", Trace::ints(sol.potentials));
    Trace counters = Trace::map();
    counters.set("dijkstra_runs", Trace::integer(sol.counters.dijkstra_runs));
    counters.set("dijkstra_pops", Trace::integer(sol.counters.dijkstra_pops));
    counters.set("settled", Trace::integer(sol.counters.settled));
    counters.set("arc_scans", Trace::integer(sol.counters.arc_scans));
    counters.set("dinic_phases", Trace::integer(sol.counters.dinic_phases));
    counters.set("augmentations", Trace::integer(sol.counters.augmentations));
    fc.set("counters", std::move(counters));
    Trace supply = Trace::map();
    Trace shares_tr = Trace::map();
    for (std::size_t i = 0; i < pe.roster.size(); ++i) {
      supply.set(key(pe.roster[i]), Trace::integer(supplies[i]));
      shares_tr.set(key(pe.roster[i]), Trace::pairs(shares[i]));
    }
    fc.set("supply", std::move(supply));
    trace->set("flow_core", std::move(fc));
    trace->set("shares", std::move(shares_tr));
  }
  if (outputs.bottleneck) out.min_cut = std::move(sol.min_cut);
  if (any_output(outputs)) prepare_outputs(pe, outputs, draft_edges(drafts), out);
  auto curve = *out.cost_distribution;  // the cumulative curve a later zero call reports (9.1)
  if (staged != nullptr && *staged) (*staged)();
  const auto ids = LedgerAccess::commit(pe.fg, pe.owner, drafts);
  // Committed: the epoch state follows (nothing below throws).
  for (std::size_t n = 0; n < fresh.size(); ++n) next[fresh[n]].front().id = ids[n].second;
  pe.cohorts.swap(next);
  pe.curve.swap(curve);
  pe.requested = requested;
  if (!v.has_value()) pe.infinite = true;
  return out;
}

// ---------------------------------------------------------------- construction

// The static reverse distances toward the allowed destinations over the statically eligible
// member graph (capacity >= 1 unit, no residual).
std::vector<Cost> static_distances(const Graph& g, const Eligibility* elig, std::span<const NodeId> dests) {
  const PathSelection full(EdgeTieBreak::Deterministic, EdgePriority::Capacity, true, true);
  const auto member = member_mask(g, elig, full, {}, 1);
  return distances(g, member, dests, Direction::Reverse);
}

} // namespace

void validate_config(const FlowPolicyConfig& cfg, std::size_t num_sources, std::size_t num_destinations) {
  const auto& r = cfg.routing;
  if (r.valueless_by_exception()) invalid("routing: the variant holds no routing mode");
  check_member(cfg.placement, "placement");
  check_member(cfg.admission, "admission");
  check_selection(cfg.selection);
  check_kernel(cfg.placement, cfg.admission, cfg.selection.priority);
  if (cfg.admission == Admission::Lossy &&
      (std::holds_alternative<Global>(r) || std::holds_alternative<CapacityAware>(r))) {
    invalid("Lossy admission is rejected with Global and CapacityAware routing (hop-by-hop forwarding without "
            "admission control neither grows nor reroutes)");
  }
  if (std::holds_alternative<Global>(r) && cfg.placement != FlowPlacement::MaxFlow) {
    invalid("Global routing requires placement = MaxFlow");
  }
  if (const auto* ca = std::get_if<CapacityAware>(&r)) {
    if (ca->flow_count < 1) invalid("CapacityAware flow_count must be >= 1 (got " + str(ca->flow_count) + ")");
    if (ca->max_flow_count.has_value() && *ca->max_flow_count < ca->flow_count) {
      invalid("CapacityAware max_flow_count (" + str(*ca->max_flow_count) + ") must be >= flow_count (" +
              str(ca->flow_count) + ")");
    }
    if (ca->min_flow_bandwidth < 0 || ca->min_flow_bandwidth >= kDomainLimit) {
      invalid("CapacityAware min_flow_bandwidth must be in [0, 2^62) units (got " + str(ca->min_flow_bandwidth) + ")");
    }
    if (ca->max_path_cost.has_value()) {
      if (const auto* a = std::get_if<Absolute>(&*ca->max_path_cost)) {
        if (a->cost < 0) invalid("max_path_cost Absolute cost must be >= 0 (got " + str(a->cost) + ")");
      } else {
        const double f = std::get<Stretch>(*ca->max_path_cost).factor;
        if (!(f >= 1.0) || !std::isfinite(f)) {
          invalid("max_path_cost Stretch factor must be finite and in [1, inf) (got " + std::to_string(f) + ")");
        }
      }
    }
    if (ca->diminishing.has_value()) {
      const auto& d = *ca->diminishing;
      if (d.window < 1) invalid("diminishing window must be >= 1 (got " + str(d.window) + ")");
      if (!(d.fraction > 0.0 && d.fraction < 1.0)) {
        invalid("diminishing fraction must be in (0, 1) (got " + std::to_string(d.fraction) + ")");
      }
    }
  }
  if (std::holds_alternative<Pinned>(r) && (num_sources != 1 || num_destinations != 1)) {
    invalid("Pinned routing requires single endpoints (bundles bind one source and one destination)");
  }
}

} // namespace detail

// ================================================================ FlowPolicyConfig

// The enum-typed fields inside their enumerations, in validate_config's order.
FlowPolicyConfig::FlowPolicyConfig(Routing routing_, FlowPlacement placement_, Admission admission_,
                                   PathSelection selection_)
    : routing(std::move(routing_)), placement(placement_), admission(admission_), selection(selection_) {
  detail::check_member(placement, "placement");
  detail::check_member(admission, "admission");
  detail::check_selection(selection);
}

// ================================================================ FlowPolicy

using detail::LedgerAccess;
using detail::PolicyEngine;

FlowPolicy::FlowPolicy(FlowLedger& ledger, std::span<const NodeId> sources, std::span<const NodeId> destinations,
                       FlowPolicyConfig config, const Eligibility* elig, DagCache* cache) {
  auto pe = std::make_unique<PolicyEngine>(ledger, std::move(config));
  pe->graph = ledger.graph();
  const auto& g = *pe->graph;
  if (elig != nullptr && elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  if (cache != nullptr && cache->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign DagCache: bound to another graph");
  }
  pe->elig = elig;
  pe->sources = detail::node_set(sources, g.num_nodes(), "sources");
  pe->destinations = detail::node_set(destinations, g.num_nodes(), "destinations");
  detail::check_disjoint(pe->sources, pe->destinations, "sources", "destinations");
  detail::validate_config(pe->config, pe->sources.size(), pe->destinations.size());
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
      auto plan = detail::make_route_plan(g, pe->config.placement, pe->config.admission,
                                          pe->config.selection.priority, src, dst, *b, elig, &cost);
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
                      ? cache->lookup(g, pe->destinations, Direction::Reverse, elig)
                      : std::make_shared<const PredDAG>(detail::static_support(g, pe->destinations, Direction::Reverse, elig));
    const auto edges = detail::PredDAGBuilder::ascending(*pe->support);
    const auto row = g.row_offsets();
    for (NodeId s : pe->sources) {
      if (!node_ok(s)) continue;
      const auto it = std::lower_bound(edges.begin(), edges.end(), row[static_cast<std::size_t>(s)]);
      if (it != edges.end() && *it < row[static_cast<std::size_t>(s) + 1]) pe->roster.push_back(s);
    }
  } else {
    const auto dist = detail::static_distances(g, elig, pe->dests_ok);
    for (NodeId s : pe->sources) {
      if (!node_ok(s) || dist[static_cast<std::size_t>(s)] == kUnreachable) continue;
      pe->roster.push_back(s);
      pe->baseline.push_back(dist[static_cast<std::size_t>(s)]);
    }
  }
  pe->cohorts.resize(pe->roster.size());
  pe->owner = LedgerAccess::bind(ledger);  // last: nothing is bound when construction fails
  engine_ = std::move(pe);
}

// The destruction law (section 8): the records stay in the ledger as orphans.
FlowPolicy::~FlowPolicy() { LedgerAccess::release(engine_->fg, engine_->owner); }

namespace {

// `staged` (optional) runs once the call's result and trace are final and before anything is
// published (the ledger commit, the epoch state): an exception from it aborts the call with
// nothing changed (the Python binding builds its result objects there).
PlacementResult place_impl(PolicyEngine& pe, std::optional<Amount> volume, const OutputOptions& outputs,
                           detail::Trace* trace, const std::function<void()>* staged) {
  // Validation before any mutation, also before the zero shortcut (section 8).
  if (volume.has_value() && (*volume < 0 || *volume >= kDomainLimit)) {
    throw std::invalid_argument("request " + std::to_string(*volume) +
                                " is outside [0, 2^62) units (requests are integer units or +inf)");
  }
  if (pe.infinite && volume.has_value()) {
    throw std::invalid_argument(
        "a finite request after +inf is rejected within an epoch (the cumulative finite request is undefined; "
        "rule 13)");
  }
  if (volume.has_value()) {
    Counter128 probe = pe.requested;
    if (!probe.add(*volume)) throw std::invalid_argument("the cumulative request counter would reach 2^127");
  }
  PlacementResult out;
  if (volume.has_value() && *volume == 0) {
    out.remaining = 0;
    out.local_outcomes.reserve(pe.roster.size());
    for (NodeId s : pe.roster) out.local_outcomes.emplace_back(s, Termination::Complete);
    // Global: the cumulative curve (9.1: its amounts sum to placed_demand(); empty before the
    // first placement); a zero call solves nothing, so it claims no cut (2.2's cap-zero payload).
    if (detail::is<Global>(pe)) out.cost_distribution = pe.curve;
    detail::prepare_outputs(pe, outputs, detail::current_edges(pe), out);
    if (staged != nullptr && *staged) (*staged)();
    return out;
  }
  if (pe.roster.empty()) {
    out.remaining = volume;
    out.termination = Termination::NoRoute;
    if (detail::is<Global>(pe)) {
      // No source reaches an allowed destination: the cumulative flow is 0 and maximal, so its
      // marginal curve is empty and the full-residual cut of 2.2 is empty (the closure of the
      // sources has no arc of positive view leaving it); None would claim a residual path.
      out.cost_distribution.emplace();
      if (outputs.bottleneck) out.min_cut.emplace();
    }
    if (detail::kTraced && trace) trace->set("roster", detail::Trace::list());
    detail::prepare_outputs(pe, outputs, detail::current_edges(pe), out);
    if (staged != nullptr && *staged) (*staged)();
    // The epoch records every positive request, also one that finds no roster (section 8);
    // nothing below throws (the counter was checked above).
    if (volume.has_value()) {
      (void)pe.requested.add(*volume);
    } else {
      pe.infinite = true;
    }
    return out;
  }
  if (detail::is<Global>(pe)) return detail::place_global(pe, volume, outputs, trace, staged);
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
  pe.curve.clear();
  pe.requested = Counter128 {};
  pe.infinite = false;
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
