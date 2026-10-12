// A manual placement: a FlowPolicy with Pinned routing and one bundle. The ledger tests
// (flow_ledger_tests.cpp, flow_ledger_fault_tests.cpp) drive the ledger through it and compare the
// real PlacementResult.
#pragma once

#include <memory>
#include <optional>
#include <utility>

#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"

namespace netgraph::core::test {

class PinnedOne {
public:
  PinnedOne(FlowLedger& fg, NodeId src, NodeId dst, std::shared_ptr<const PredDAG> bundle, FlowPlacement placement,
            Admission admission, EdgePriority priority, std::shared_ptr<const Eligibility> elig = nullptr)
      : elig_(std::move(elig)) {
    Pinned pinned;
    pinned.bundles.push_back(std::move(bundle));
    FlowPolicyConfig cfg {pinned, placement, admission, PathSelection(EdgeTieBreak::Deterministic, priority)};
    const NodeId s[1] = {src};
    const NodeId d[1] = {dst};
    policy_ = std::make_unique<FlowPolicy>(fg, s, d, std::move(cfg), elig_.get());
  }

  PlacementResult place(std::optional<Amount> volume) { return policy_->place(volume); }
  void remove() { policy_->remove(); }
  [[nodiscard]] std::optional<FlowId> flow_id() const {
    const auto& pe = engine();
    if (pe.cohorts.empty() || pe.cohorts.front().empty()) return std::nullopt;
    return FlowId {pe.fg.uid(), pe.cohorts.front().front().id};
  }
  [[nodiscard]] int usable_bundles() const { return policy_->usable_bundles(); }
  [[nodiscard]] const detail::KernelPlan* plan() const {
    const auto& pe = engine();
    if (pe.cohorts.empty() || pe.cohorts.front().empty()) return nullptr;
    return pe.cohorts.front().front().plan.get();
  }
  [[nodiscard]] const Counter128& requested() const { return engine().requested; }
  [[nodiscard]] bool infinite() const { return engine().infinite; }
  // The wrapped policy (for the shared comparators of state_comparators.hpp, such as epoch()).
  [[nodiscard]] FlowPolicy& policy() { return *policy_; }

private:
  [[nodiscard]] const detail::PolicyEngine& engine() const { return detail::PolicyAccess::engine(*policy_); }
  std::shared_ptr<const Eligibility> elig_;
  std::unique_ptr<FlowPolicy> policy_;
};

} // namespace netgraph::core::test
