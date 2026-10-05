/*
  A Pinned placement with one bundle (the single-flow transaction of sub-PR 1b): the
  per-flow kernel steps of detail/kernels.hpp committed through the ledger's internal entry
  points, one atomic exchange per call (APPROACH 8, 9.2).
*/
#include "netgraph/core/detail/pinned_flow.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"

namespace netgraph::core::detail {

PinnedFlow::PinnedFlow(FlowGraph& fg, NodeId src, NodeId dst, std::shared_ptr<const PredDAG> bundle,
                       FlowPlacement placement, Admission admission, EdgePriority priority,
                       std::shared_ptr<const Eligibility> elig)
    : fg_(fg), src_(src), dst_(dst), bundle_(std::move(bundle)), placement_(placement), admission_(admission),
      priority_(priority), elig_(std::move(elig)) {
  check_kernel(placement, admission, priority);
  const auto& g = *fg_.graph();
  for (NodeId v : {src, dst}) {
    if (v < 0 || v >= g.num_nodes()) {
      throw std::invalid_argument("Pinned: node " + std::to_string(v) + " out of node range [0, " +
                                  std::to_string(g.num_nodes()) + ")");
    }
  }
  if (src == dst) {
    throw std::invalid_argument("node " + std::to_string(src) +
                                " is in both endpoint sets: sources and destinations must not overlap");
  }
  if (!bundle_) throw std::invalid_argument("Pinned: the bundle is null");
  if (bundle_->graph_uid() != g.uid()) throw std::invalid_argument("foreign PredDAG: built for another graph");
  if (elig_ && elig_->graph_uid() != g.uid()) throw std::invalid_argument("foreign Eligibility: built for another graph");
  // The usable route: statically eligible routable members on src -> dst walks. An empty one
  // is either all-down (no_route at placement) or malformed (no walk at all: rejected).
  const NodeId dsts[1] = {dst};
  const auto edges = PredDAGBuilder::ascending(*bundle_);
  route_plan_ = make_route_plan(g, placement, admission, priority, src, dsts, *bundle_, elig_.get(), &route_cost_);
  if (route_plan_.route.empty() && prune(g, edges, src, dsts).empty()) {
    throw std::invalid_argument("Pinned: malformed bundle, no src -> dst walk (src " + std::to_string(src) + ", dst " +
                                std::to_string(dst) + ")");
  }
  owner_ = LedgerAccess::bind(fg_);
  bound_ = true;
}

PinnedFlow::~PinnedFlow() { release(); }

void PinnedFlow::release() noexcept {
  if (!bound_) return;
  LedgerAccess::release(fg_, owner_);
  bound_ = false;
}

void PinnedFlow::check_bound() const {
  if (!bound_) throw std::logic_error("Pinned: the policy was released (destroyed)");
}

Amount PinnedFlow::placed_demand() const {
  check_bound();
  return flow_id_.has_value() ? LedgerAccess::live(fg_, flow_id_->second).delivered : 0;
}

void PinnedFlow::remove() {
  check_bound();
  (void)LedgerAccess::detach(fg_, owner_);
  flow_id_.reset();
  plan_ = nullptr;
  captured_.reset();
  total_ = 0;
  requested_ = Counter128 {};
  infinite_ = false;
}

namespace {

// Member index of e in the ascending member list (members contain every allocated edge).
std::size_t member_index(const std::vector<EdgeId>& members, EdgeId e) {
  const auto it = std::lower_bound(members.begin(), members.end(), e);
  if (it == members.end() || *it != e) throw std::logic_error("Pinned: an allocation outside the split graph");
  return static_cast<std::size_t>(it - members.begin());
}

// Per-thread member-aligned buffers of one call (nothing survives a call).
struct CallScratch {
  std::vector<Amount> own;
  std::vector<Amount> working;
  std::vector<Amount> view;
  std::vector<Counter128> drops;
};

CallScratch& call_scratch() {
  thread_local CallScratch s;
  return s;
}

} // namespace

PinnedResult PinnedFlow::place(std::optional<Amount> volume) {
  check_bound();
  if (volume.has_value() && (*volume < 0 || *volume >= kDomainLimit)) {
    throw std::invalid_argument("request " + std::to_string(*volume) +
                                " is outside [0, 2^62) units (requests are integer units or +inf)");
  }
  if (infinite_ && volume.has_value()) {
    throw std::invalid_argument(
        "a finite request after +inf is rejected within an epoch (the cumulative finite request is undefined; "
        "rule 13)");
  }
  PinnedResult out;
  if (volume.has_value() && *volume == 0) {
    out.remaining = 0;
    return out;
  }
  const auto& g = *fg_.graph();
  const auto base = LedgerAccess::baseline(fg_);
  const auto held = LedgerAccess::allocated(fg_);
  const bool created = !flow_id_.has_value();
  std::optional<KernelPlan> fresh;
  const KernelPlan* plan = plan_;
  if (created) {
    if (route_plan_.route.empty()) {
      // Nothing is created, but the epoch's call state records every positive request
      // (the cumulative request and the finite-or-infinite state; APPROACH 8, 9.2).
      if (volume.has_value()) {
        Counter128 requested = requested_;
        if (!requested.add(*volume)) {
          throw std::invalid_argument("the cumulative requested counter would reach 2^127");
        }
        requested_ = requested;
      } else {
        infinite_ = true;
      }
      out.remaining = volume;
      out.termination = Termination::NoRoute;
      return out;
    }
    if (placement_ == FlowPlacement::Weighted && priority_ == EdgePriority::Residual) {
      // Membership and weights captured on the creation view W (route-aligned).
      auto& view = call_scratch().view;
      view.clear();
      for (EdgeId e : route_plan_.route) view.push_back(base[static_cast<std::size_t>(e)] - held[static_cast<std::size_t>(e)]);
      const NodeId dsts[1] = {dst_};
      fresh = make_plan(g, placement_, admission_, priority_, src_, dsts, route_plan_.route, view);
      plan = &*fresh;
    } else {
      plan = &route_plan_;
    }
  }
  const auto& members = plan->members;
  const FlowRecord* old = flow_id_.has_value() ? &LedgerAccess::live(fg_, flow_id_->second) : nullptr;
  CallScratch& scratch = call_scratch();
  auto& own_m = scratch.own;
  own_m.assign(members.size(), 0);
  if (old != nullptr) {
    for (const auto& [e, x] : old->edges) own_m[member_index(members, e)] = x;
  }
  auto& w_m = scratch.working;
  w_m.resize(members.size());
  for (std::size_t m = 0; m < members.size(); ++m) {
    const auto e = static_cast<std::size_t>(members[m]);
    w_m[m] = base[e] - held[e];
  }

  Draft d;
  if (flow_id_.has_value()) d.id = flow_id_->second;
  d.sources = {src_};
  d.destinations = {dst_};
  d.placement = placement_;
  d.admission = admission_;
  const Amount delivered0 = old != nullptr ? old->delivered : 0;
  bool served = false;
  Amount total_after = total_;
  if (admission_ == Admission::Lossy) {
    const auto step = lossy_step(*plan, volume, w_m);
    Counter128 offered = (old != nullptr && old->offered.has_value()) ? *old->offered : Counter128 {};
    if (!offered.add(step.offered)) {
      throw std::invalid_argument("the cumulative offered counter would reach 2^127");
    }
    d.offered = offered;
    auto& drops_m = scratch.drops;
    drops_m.assign(members.size(), Counter128 {});
    if (old != nullptr) {
      for (const auto& [e, x] : old->drops) drops_m[member_index(members, e)] = x;
    }
    for (std::size_t m = 0; m < members.size(); ++m) {
      if (step.drops[m] == 0) continue;
      if (!drops_m[m].add(step.drops[m])) {
        throw std::invalid_argument("a cumulative drop counter would reach 2^127");
      }
      out.drops.emplace_back(members[m], step.drops[m]);
    }
    for (std::size_t m = 0; m < members.size(); ++m) {
      if (own_m[m] + step.carried[m] != 0) d.edges.emplace_back(members[m], own_m[m] + step.carried[m]);
      if (!drops_m[m].is_zero()) d.drops.emplace_back(members[m], drops_m[m]);
    }
    d.delivered = delivered0 + step.delivered;
    served = volume.has_value() && step.injected == *volume;
  } else {
    auto& v_m = scratch.view;
    v_m.resize(members.size());
    for (std::size_t m = 0; m < members.size(); ++m) v_m[m] = w_m[m] + own_m[m];
    const auto step = lossless_step(g, *plan, total_, v_m, own_m, volume);
    for (std::size_t m = 0; m < members.size(); ++m) {
      if (step.alloc[m] != 0) d.edges.emplace_back(members[m], step.alloc[m]);
    }
    d.delivered = step.total;
    total_after = step.total;
    out.iterations = step.iterations;
    served = volume.has_value() && step.total - delivered0 == *volume;
  }
  // The record's bundle is the usable route: the supplied DAG itself when nothing was pruned,
  // else one DAG of the route made at the first creation and kept (the route is fixed).
  std::shared_ptr<const PredDAG> route_dag = route_dag_;
  if (!route_dag) {
    route_dag = route_plan_.route.size() == bundle_->via_edges().size()
                    ? bundle_
                    : std::make_shared<const PredDAG>(PredDAGBuilder::from_edges(g, route_plan_.route));
  }
  d.bundle = route_dag;
  d.route_cost = route_cost_;
  const Amount delivered = d.delivered;
  Counter128 requested = requested_;
  if (volume.has_value() && !requested.add(*volume)) {
    throw std::invalid_argument("the cumulative requested counter would reach 2^127");
  }
  const auto ids = LedgerAccess::commit(fg_, owner_, std::span<Draft>(&d, 1));
  // Committed: the epoch state follows (nothing below throws).
  route_dag_ = std::move(route_dag);
  if (created) {
    if (fresh.has_value()) {
      captured_ = std::move(fresh);
      plan_ = &*captured_;
    } else {
      plan_ = &route_plan_;
    }
    flow_id_ = ids.front();
  }
  total_ = total_after;
  requested_ = requested;
  if (!volume.has_value()) infinite_ = true;
  out.placed = delivered - delivered0;
  if (volume.has_value()) out.remaining = *volume - out.placed;
  out.termination = served ? Termination::Complete : Termination::Exhausted;
  out.routes_changed = created;
  return out;
}

} // namespace netgraph::core::detail
