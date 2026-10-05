/* The integer ledger FlowGraph (APPROACH 8). */
#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <vector>

#include "netgraph/core/flow.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

namespace detail {
struct LedgerAccess;
} // namespace detail

// The ledger of one graph: an immutable baseline residual in units and the records of the
// flows placed on it. The ledger is mutated only through policies (bound to it by the
// internal entry points of detail::LedgerAccess); there is no public placement, removal or
// reservation. Aggregates are exact: allocated() is the sum of every record's edges and
// residual() = baseline - allocated(), both int64[E] snapshots in EdgeId order.
//
// Lifetime: the graph outlives the ledger (shared pointer), the ledger outlives its policies.
// A ledger with bound policies cannot be copied, moved or reset. A copy of an unbound ledger
// is a new ledger (new uid) holding the same records and ids under the new uid; a moved-from
// ledger is empty and unusable (every call throws std::logic_error). A policy destroyed
// without remove() leaves its records and allocations as orphans (readable, never changed
// again) until reset().
class FlowGraph {
public:
  // `residual`: the baseline in units, one entry per edge in [0, capacity_units[e]]
  // (std::invalid_argument otherwise); absent: the capacities.
  explicit FlowGraph(std::shared_ptr<const StrictMultiDiGraph> graph);
  FlowGraph(std::shared_ptr<const StrictMultiDiGraph> graph, std::span<const Amount> residual);
  FlowGraph(const FlowGraph& other);
  FlowGraph(FlowGraph&& other);
  FlowGraph& operator=(const FlowGraph&) = delete;
  FlowGraph& operator=(FlowGraph&&) = delete;
  ~FlowGraph() = default;

  [[nodiscard]] std::uint64_t uid() const;
  [[nodiscard]] const std::shared_ptr<const StrictMultiDiGraph>& graph() const;
  // The graph's static capacities (the ledger's view of them).
  [[nodiscard]] std::span<const Amount> capacity_units() const;
  [[nodiscard]] std::vector<Amount> residual() const;
  [[nodiscard]] std::vector<Amount> allocated() const;
  // Live records (bound and orphaned) as snapshots in FlowId order.
  [[nodiscard]] std::vector<FlowRecord> records() const;
  // One record; a FlowId of another ledger, a retired id or an unknown id throws
  // std::invalid_argument.
  [[nodiscard]] FlowRecord record(const FlowId& flow_id) const;
  // Saturation diagnostic: the edges (statically eligible, capacity >= 1 unit) leaving the set
  // reachable from the allowed sources over eligible forward arcs with residual >= 1 unit, in
  // ascending EdgeId. Not a cut; empty when the sources reach everything they can.
  [[nodiscard]] std::vector<EdgeId> bottleneck(std::span<const NodeId> sources, const Eligibility* elig = nullptr) const;
  // A new baseline (the capacities when absent); orphaned records and their allocations are
  // cleared and their ids retired. Rejected (std::invalid_argument) while a policy is bound.
  void reset();
  void reset(std::span<const Amount> residual);

private:
  friend struct detail::LedgerAccess;
  void check_usable() const;
  std::vector<Amount> checked_baseline(std::span<const Amount> residual) const;

  std::shared_ptr<const StrictMultiDiGraph> graph_ {};
  std::uint64_t uid_ {0};
  bool usable_ {true};
  std::vector<Amount> baseline_ {};
  std::vector<Amount> allocated_ {};
  std::map<std::int64_t, FlowRecord> records_ {};
  std::map<std::int64_t, std::vector<std::int64_t>> owner_ids_ {};  // owner -> live ids, ascending
  std::int64_t next_id_ {0};
  std::int64_t next_owner_ {0};
  // Every owner since the last reset: true while bound, false once orphaned (release() flips
  // the flag, so it never allocates).
  std::map<std::int64_t, bool> owners_ {};
  std::size_t num_bound_ {0};
};

} // namespace netgraph::core
