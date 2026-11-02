/* The ledger's internal entry points: what a policy uses to place on and remove from it.

   Policies place and remove flows through these entry points; an unbound ledger can also be
   reset, through its public reset(), which is not one of them. A policy binds (an owner id),
   then every call ends with one atomic exchange of all of its records (`commit`); remove()
   reverts its allocations and retires its ids (`detach`); a policy destroyed without remove()
   leaves its records as orphans (`release`). None of this is public API.
*/
#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/detail/concurrency.hpp"
#include "netgraph/core/detail/ledger_entry.hpp"

namespace netgraph::core::detail {

// One owner's complete entry in an atomic batch. A continuing epoch keeps its inventory
// and first counter; a full-offer call retires that range and receives new identities.
struct OwnerEntry {
  std::int64_t owner {0};
  LedgerEntry* entry {nullptr};
  bool retire {false};
};

// The saturation diagnostic of FlowLedger::saturated_boundary_edges on a given residual (the sources already a
// validated set): the edges (statically eligible, capacity >= 1 unit) leaving the set reachable
// from the allowed sources over eligible forward arcs with residual >= 1 unit, ascending EdgeId.
[[nodiscard]] std::vector<EdgeId> saturated_boundary_edges_of(const Graph& graph, std::span<const Amount> residual,
                                                std::span<const NodeId> sources, const Eligibility* eligibility);

struct LedgerAccess {
  [[nodiscard]] static std::unique_lock<std::mutex> lock(const FlowLedger& ledger) {
    concurrent_event(ConcurrentStage::LedgerWait, &ledger);
    std::unique_lock guard(ledger.gate_);
    concurrent_event(ConcurrentStage::LedgerAcquired, &ledger);
    return guard;
  }
  // Immutable metadata; no mutable-state check or lock. All other internal access below
  // requires the ledger gate, except release(), which acquires it for destruction.
  [[nodiscard]] static const std::shared_ptr<const Graph>& graph(const FlowLedger& ledger) noexcept { return ledger.graph_; }
  [[nodiscard]] static FlowRecord record(const FlowLedger& ledger, const FlowId& id);
  // A new owner id (0, 1, ... per ledger, below the largest int64: std::runtime_error `policy
  // identity limit exceeded` there, nothing changed). Throws on a moved-from ledger.
  [[nodiscard]] static std::int64_t bind(FlowLedger& ledger);
  // The owner's policy is destroyed without remove(): its records stay as orphans (an owner
  // without records is forgotten).
  static void release(FlowLedger& ledger, std::int64_t owner) noexcept;
  // remove(): the owner's records are deleted (allocations and holds reverted), their ids retired.
  // Throws only before any change (a moved-from ledger, an unbound owner).
  static void detach(FlowLedger& ledger, std::int64_t owner);
  // The owner's live entry, or null before its first records (requires the gate).
  [[nodiscard]] static const LedgerEntry* entry(const FlowLedger& ledger, std::int64_t owner);
  // An internal record snapshot by counter; public lookup also checks the ledger uid first.
  [[nodiscard]] static FlowRecord materialize(const FlowLedger& ledger, std::int64_t counter);
  [[nodiscard]] static std::vector<FlowRecord> records(const FlowLedger& ledger, std::int64_t owner);
  // Every buffer and map node is staged before the nonthrowing publication. No identity or
  // aggregate changes on failure. Existing nodes and values are reused only after staging.
  static void commit(FlowLedger& ledger, std::int64_t owner, LedgerEntry& entry);
  static void exchange(FlowLedger& ledger, std::int64_t owner, LedgerEntry& entry);
  static void commit_batch(FlowLedger& ledger, std::span<OwnerEntry> owners);
  // The working residual service baseline - carried() into `out` (resized to E).
  static void residual_into(const FlowLedger& ledger, std::vector<Amount>& out);
  [[nodiscard]] static std::span<const Amount> baseline(const FlowLedger& ledger);
  // The aggregate allocation (borrowed; valid until the next commit, detach or reset).
  [[nodiscard]] static std::span<const Amount> carried(const FlowLedger& ledger);
  // The reservation baseline R0 and the aggregate holds S (borrowed, as above). S is empty until a
  // commit with holds stages it as zeros (a commit that fails after that leaves it staged); empty
  // and staged S read as all zero.
  [[nodiscard]] static std::span<const Amount> reservation_baseline(const FlowLedger& fg);
  [[nodiscard]] static std::span<const Amount> reserved(const FlowLedger& fg);
  // S as stored, made E zeros first when nothing has been held (the tests' state checks write
  // through it).
  static std::vector<Amount>& reserved_storage(FlowLedger& fg);
  // Diagnostics (the tests and the private bindings).
  [[nodiscard]] static std::int64_t next_record_counter(const FlowLedger& ledger);
  // The next record and owner ids, set (the tests' checks of the identity limits).
  static void set_next_ids(FlowLedger& ledger, std::int64_t next_record_counter, std::int64_t next_policy_id);
  [[nodiscard]] static std::vector<std::int64_t> bound_policy_ids(const FlowLedger& ledger);
  [[nodiscard]] static std::vector<std::int64_t> orphan_policy_ids(const FlowLedger& ledger);
};

} // namespace netgraph::core::detail
