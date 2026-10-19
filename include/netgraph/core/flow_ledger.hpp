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

// A ledger tracks the flows placed on one graph in two projections that are never added: the
// carried traffic (the service baseline B, the available capacity in units before these flows;
// allocated() L, the exact sum of all records' edge amounts; residual() = B - L) and the strict
// holds of T1 tunnels (the reservation baseline R0; reserved() S, the exact sum of all records'
// reserved_edges; unreserved() = R0 - S). Each projection keeps its own bound (0 <= L <= B,
// 0 <= S <= R0 per edge); max(L, S) is no occupancy, and traffic may use bandwidth that is held.
// Policies alone place, hold or remove through detail::LedgerAccess; there is no public
// placement, removal or reservation. All four return int64[E] snapshots in EdgeId order.
//
// Lifetime: the graph outlives the ledger (shared pointer), the ledger outlives its policies.
// A ledger with bound policies cannot be copied, moved or reset. A copy of an unbound ledger
// is a new ledger (new uid) holding the same baselines, records, ids and holds under the new uid;
// a moved-from ledger is empty and unusable (every call throws std::logic_error). A policy
// destroyed without remove() leaves its records, allocations and holds as orphans (readable,
// never changed again) until reset().
class FlowLedger {
public:
  // `residual` and `unreserved`: the service and reservation baselines in units, each one entry
  // per edge in [0, capacity_units[e]] (std::invalid_argument naming the array, residual checked
  // first); absent: the capacities.
  explicit FlowLedger(std::shared_ptr<const Graph> graph);
  FlowLedger(std::shared_ptr<const Graph> graph, std::span<const Amount> residual);
  FlowLedger(std::shared_ptr<const Graph> graph, std::span<const Amount> residual, std::span<const Amount> unreserved);
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
  [[nodiscard]] std::vector<Amount> reserved() const;
  [[nodiscard]] std::vector<Amount> unreserved() const;
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
  // New baselines (the capacities when absent); orphaned records with their allocations and holds
  // are cleared and their ids retired. Rejected (std::invalid_argument) while a policy is bound.
  void reset();
  void reset(std::span<const Amount> residual);
  void reset(std::span<const Amount> residual, std::span<const Amount> unreserved);

private:
  friend struct detail::LedgerAccess;
  void check_usable() const;
  std::vector<Amount> checked_baseline(std::span<const Amount> amounts, const char* label) const;
  void check_resettable() const;
  void install(std::vector<Amount> baseline, std::vector<Amount> reservation);
  // The reservation projection costs nothing until something holds: an empty reserved_ is all
  // zero, an empty reservation_baseline_ is the static capacity.
  [[nodiscard]] Amount reservation_at(std::size_t e) const noexcept {
    return reservation_baseline_.empty() ? graph_->capacity_units()[e] : reservation_baseline_[e];
  }
  [[nodiscard]] Amount reserved_at(std::size_t e) const noexcept { return reserved_.empty() ? 0 : reserved_[e]; }

  std::shared_ptr<const Graph> graph_ {};
  std::uint64_t uid_ {0};
  bool usable_ {true};
  std::vector<Amount> baseline_ {};
  std::vector<Amount> allocated_ {};
  std::vector<Amount> reservation_baseline_ {};
  std::vector<Amount> reserved_ {};
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
