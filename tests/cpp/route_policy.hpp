// A policy on one fixed route, in the program API: CentralizedTE restricted to the route (on_support) or
// Pinned over a program of the route as one supplied split route (split_route: Equal, or
// Weighted with its WeightBasis; RESIDUAL_SNAPSHOT reads the ledger's residual when the program
// is made). The ledger tests (flow_ledger_tests.cpp, flow_ledger_fault_tests.cpp) drive the ledger
// through RoutePolicy and compare the real PlacementResult.
#pragma once

#include <memory>
#include <optional>
#include <utility>

#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/forwarding_program.hpp"

namespace netgraph::core::test {

inline FlowPolicyConfig on_support(std::shared_ptr<const PredecessorDAG> dag, Admission admission = Admission::Lossless) {
  FlowPolicyConfig cfg;
  CentralizedTE controller(SourceAllocation::Optimized);
  controller.support = std::move(dag);
  cfg.routing = std::move(controller);
  cfg.admission = admission;
  return cfg;
}

inline FlowPolicyConfig split_route(const FlowLedger& fg, NodeId src, NodeId dst, std::shared_ptr<const PredecessorDAG> dag,
                                    RouteAllocation placement, Admission admission,
                                    WeightBasis weights = WeightBasis::Capacity) {
  const NodeId s[1] = {src};
  const NodeId d[1] = {dst};
  RouteOptions options;
  options.allocation = placement;
  if (placement == RouteAllocation::Weighted && weights != WeightBasis::Capacity) {
    options.weight_basis = weights;
    if (weights == WeightBasis::ResidualSnapshot) options.residual = fg.residual();
  }
  const RouteDefinition routes[1] = {RouteDefinition {src, std::move(dag), 1}};
  FlowPolicyConfig cfg;
  cfg.routing = Pinned {ForwardingProgram::from_routes(*fg.graph(), s, d, routes, options), nullptr};
  cfg.admission = admission;
  return cfg;
}

class RoutePolicy {
public:
  RoutePolicy(FlowLedger& fg, NodeId src, NodeId dst, FlowPolicyConfig cfg,
              std::shared_ptr<const Eligibility> elig = nullptr)
      : elig_(std::move(elig)) {
    const NodeId s[1] = {src};
    const NodeId d[1] = {dst};
    policy_ = std::make_unique<FlowPolicy>(fg, s, d, std::move(cfg), elig_.get());
  }

  PlacementResult place(std::optional<Amount> offer) { return policy_->place(offer); }
  void remove() { policy_->remove(); }
  // The FlowId of the policy's one record (its first), if it has one.
  [[nodiscard]] std::optional<FlowId> flow_id() const {
    const auto& pe = engine();
    if (pe.epoch.record_groups.empty() || pe.epoch.record_groups.front().empty()) return std::nullopt;
    return FlowId {pe.ledger.uid(), pe.epoch.record_groups.front().front().record_counter};
  }
  // The split route's plan as forwarded (Pinned), or null (CentralizedTE's support has none).
  [[nodiscard]] const detail::PreparedRoute* plan() const {
    const auto& pe = engine();
    if (!std::holds_alternative<Pinned>(pe.config.routing) || pe.view.data->split_contexts.empty()) return nullptr;
    return pe.view.prepared_route(0).get();
  }
  [[nodiscard]] const Counter128& requested() const { return engine().epoch.requested; }
  [[nodiscard]] bool infinite() const { return engine().epoch.infinite; }
  [[nodiscard]] bool lossy() const { return engine().config.admission == Admission::Lossy; }
  // The wrapped policy (for the shared comparators of state_comparators.hpp, such as epoch()).
  [[nodiscard]] FlowPolicy& policy() { return *policy_; }

private:
  [[nodiscard]] const detail::PolicyEngine& engine() const { return detail::PolicyAccess::engine(*policy_); }
  std::shared_ptr<const Eligibility> elig_;
  std::unique_ptr<FlowPolicy> policy_;
};

} // namespace netgraph::core::test
