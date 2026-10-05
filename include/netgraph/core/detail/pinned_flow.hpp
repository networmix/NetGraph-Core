/* A Pinned placement with one bundle: the single-flow transaction of sub-PR 1b (internal).

   APPROACH 8: a fixed-structure ("manual") placement is a Pinned policy with one bundle. Until
   the policy engine of sub-PR 1c exists, this driver runs that transaction (APPROACH 9.2
   restricted to one flow on one bundle) for the tests, the A/B/A kernel rows and the private
   recorder entry, composing the per-flow kernel steps with the ledger's internal entry
   points; 1c's FlowPolicy with a one-bundle Pinned routing runs the same steps and replaces it.

   * place(0) is a no-op that creates nothing; a positive request creates the flow when the
     bundle is usable (a src -> dst walk survives static eligibility and the capacity >= 1
     unit rule), else nothing is created and the call ends no_route;
   * MaxFlow (replacement): a fresh flow-core solve on V = W + own over the split graph with
     cap T + v (+inf: uncapped);
   * Weighted and EqualBalanced (incremental): T_new = max(T, min(T + v, T_max(V)));
   * Lossy: the offer split by the kernel's rule on W, per call;
   * a finite request after +inf is rejected within an epoch (also 0); remove() reverts the
     allocation, retires the id and opens a new epoch; destruction without remove() leaves the
     record as an orphan.
*/
#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/flow.hpp"
#include "netgraph/core/flow_graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"

namespace netgraph::core::detail {

// APPROACH 9.4 local outcomes (the public PlacementResult arrives with the policies in 1c).
enum class Termination { Complete, LimitFlows, LimitDiminishing, NoRoute, Exhausted };

struct PinnedResult {
  Amount placed {0};
  std::optional<Amount> remaining {};  // nullopt: +inf
  Termination termination {Termination::Complete};
  bool routes_changed {false};
  bool weights_changed {false};
  std::vector<std::pair<EdgeId, Amount>> drops {};  // this call's transit drops
  std::int32_t iterations {0};  // diagnostic: bisection steps of this call's T_max search
};

class PinnedFlow {
public:
  // Validates (std::invalid_argument): the kernel combination (Lossy with MaxFlow, Weighted
  // with Headroom), the endpoints (in range, distinct), a bundle or Eligibility of another
  // graph, a bundle without any src -> dst walk before pruning (malformed).
  PinnedFlow(FlowGraph& fg, NodeId src, NodeId dst, std::shared_ptr<const PredDAG> bundle, FlowPlacement placement,
             Admission admission, EdgePriority priority, std::shared_ptr<const Eligibility> elig = nullptr);
  PinnedFlow(const PinnedFlow&) = delete;
  PinnedFlow& operator=(const PinnedFlow&) = delete;
  ~PinnedFlow();

  // `volume` nullopt is +inf; a finite request is an integer in [0, 2^62).
  PinnedResult place(std::optional<Amount> volume);
  void remove();
  // The policy is destroyed without remove(): the record stays as an orphan; later calls throw.
  void release() noexcept;

  [[nodiscard]] std::optional<FlowId> flow_id() const noexcept { return flow_id_; }
  [[nodiscard]] Amount placed_demand() const;
  [[nodiscard]] int usable_bundles() const noexcept { return route_plan_.route.empty() ? 0 : 1; }
  // The live flow's plan, frozen at its creation (null without a flow).
  [[nodiscard]] const KernelPlan* plan() const noexcept { return plan_; }
  [[nodiscard]] std::int64_t owner() const noexcept { return owner_; }
  // The epoch's call state: the cumulative finite request (positive calls that reached the
  // ledger) and whether +inf was requested.
  [[nodiscard]] const Counter128& requested() const noexcept { return requested_; }
  [[nodiscard]] bool infinite() const noexcept { return infinite_; }

private:
  void check_bound() const;

  FlowGraph& fg_;
  NodeId src_;
  NodeId dst_;
  std::shared_ptr<const PredDAG> bundle_;
  FlowPlacement placement_;
  Admission admission_;
  EdgePriority priority_;
  std::shared_ptr<const Eligibility> elig_;
  // The usable route (pruned) as a plan with members == route: the plan of every kernel with a
  // static membership. Its cost and its DAG (the record's bundle) do not change either.
  KernelPlan route_plan_ {};
  Cost route_cost_ {kUnreachable};
  std::shared_ptr<const PredDAG> route_dag_ {};  // made at the first creation when pruned
  std::int64_t owner_ {0};
  bool bound_ {false};
  // Epoch state.
  std::optional<FlowId> flow_id_ {};
  std::optional<KernelPlan> captured_ {};  // a Weighted Residual plan, captured at creation
  const KernelPlan* plan_ {nullptr};       // &route_plan_ or &*captured_ while a flow lives
  Amount total_ {0};
  Counter128 requested_ {};
  bool infinite_ {false};
};

} // namespace netgraph::core::detail
