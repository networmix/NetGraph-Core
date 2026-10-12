/* The integer ledger FlowLedger. */
#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <vector>

#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

namespace detail {
struct LedgerAccess;
} // namespace detail

// A ledger tracks the flows placed on one graph. Its baseline is the available capacity in
// units before these flows; its residual is the capacity still unused. Policies alone place
// or remove flows through detail::LedgerAccess; there is no public placement, removal or
// reservation. allocated() is the exact sum of all records' edge amounts, and
// residual() = baseline - allocated(); both return int64[E] snapshots in EdgeId order.
//
// Lifetime: the graph outlives the ledger (shared pointer), the ledger outlives its policies.
// A ledger with bound policies cannot be copied, moved or reset. A copy of an unbound ledger
// is a new ledger (new uid) holding the same records and ids under the new uid; a moved-from
// ledger is empty and unusable (every call throws std::logic_error). A policy destroyed
// without remove() leaves its records and allocations as orphans (readable, never changed
// again) until reset().
class FlowLedger {
public:
  // `residual`: the baseline in units, one entry per edge in [0, capacity_units[e]]
  // (std::invalid_argument otherwise); absent: the capacities.
  explicit FlowLedger(std::shared_ptr<const Graph> graph);
  FlowLedger(std::shared_ptr<const Graph> graph, std::span<const Amount> residual);
  FlowLedger(const FlowLedger& other);
  FlowLedger(FlowLedger&& other);
  FlowLedger& operator=(const FlowLedger&) = delete;
  FlowLedger& operator=(FlowLedger&&) = delete;
  ~FlowLedger() = default;

  [[nodiscard]] std::uint64_t uid() const;
  [[nodiscard]] const std::shared_ptr<const Graph>& graph() const;
  // The graph's static capacities (the ledger's view of them).
  [[nodiscard]] std::span<const Amount> capacity_units() const;
  [[nodiscard]] std::vector<Amount> residual() const;
  [[nodiscard]] std::vector<Amount> allocated() const;
  // Live records (bound and orphaned) as snapshots in FlowId order.
  [[nodiscard]] std::vector<FlowRecord> records() const;
  // One record; a FlowId of another ledger, a retired id or an unknown id throws
  // std::invalid_argument.
  [[nodiscard]] FlowRecord record(const FlowId& flow_id) const;
  // Saturation diagnostic: follow eligible forward edges with residual >= 1 unit from the
  // allowed sources. Return statically eligible edges with capacity >= 1 unit leaving that
  // reachable set, in ascending EdgeId; empty if no such edge remains. This is not a minimum-cut
  // certificate: it does not follow reverse arcs that could cancel existing flow.
  [[nodiscard]] std::vector<EdgeId> bottleneck(std::span<const NodeId> sources, const Eligibility* elig = nullptr) const;
  // A new baseline (the capacities when absent); orphaned records and their allocations are
  // cleared and their ids retired. Rejected (std::invalid_argument) while a policy is bound.
  void reset();
  void reset(std::span<const Amount> residual);

private:
  friend struct detail::LedgerAccess;
  void check_usable() const;
  std::vector<Amount> checked_baseline(std::span<const Amount> residual) const;

  std::shared_ptr<const Graph> graph_ {};
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
