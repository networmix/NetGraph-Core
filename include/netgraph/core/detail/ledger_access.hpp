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

namespace netgraph::core::detail {

// A record proposed by a policy for the commit: `id` nullopt is a new flow (its id is assigned
// at the commit, in draft order); the commit sets the record's flow_id and policy_id. The record's
// `edge_flows`, `edge_drops` and `edge_reservations` ascend by EdgeId; zero entries are dropped.
struct RecordDraft {
  std::optional<std::int64_t> record_counter {};
  FlowRecord record {};
};

// One owner's part of a batch exchange (LedgerAccess::commit_batch).
struct OwnerDrafts {
  std::int64_t owner {0};
  std::span<RecordDraft> drafts {};
  bool retire {false};
};

// The saturation diagnostic of FlowLedger::saturated_boundary_edges on a given residual (the sources already a
// validated set): the edges (statically eligible, capacity >= 1 unit) leaving the set reachable
// from the allowed sources over eligible forward arcs with residual >= 1 unit, ascending EdgeId.
[[nodiscard]] std::vector<EdgeId> saturated_boundary_edges_of(const Graph& graph, std::span<const Amount> residual,
                                                std::span<const NodeId> sources, const Eligibility* eligibility);

struct LedgerAccess {
  // A new owner id (0, 1, ... per ledger, below the largest int64: std::runtime_error `policy
  // identity limit exceeded` there, nothing changed). Throws on a moved-from ledger.
  [[nodiscard]] static std::int64_t bind(FlowLedger& ledger);
  // The owner's policy is destroyed without remove(): its records stay as orphans (an owner
  // without records is forgotten).
  static void release(FlowLedger& ledger, std::int64_t owner) noexcept;
  // remove(): the owner's records are deleted (allocations and holds reverted), their ids retired.
  // Throws only before any change (a moved-from ledger, an unbound owner).
  static std::vector<std::int64_t> detach(FlowLedger& ledger, std::int64_t owner);
  // One atomic exchange: `drafts` must carry every live record of the owner (by id) plus the
  // new ones; the new aggregates must stay within [0, baseline] on every edge (carried loads
  // within the service baseline, holds within the reservation baseline). Throws (and
  // changes nothing) otherwise, an allocation failure included: every allocation happens
  // before a publication step that cannot throw, and the new ids running past the largest int64
  // (std::runtime_error `record identity limit exceeded`). Returns the FlowIds of the new records in
  // draft order. The drafts' vectors are moved from.
  static std::vector<FlowId> commit(FlowLedger& ledger, std::int64_t owner, std::span<RecordDraft> drafts);
  // The exchange of replace() and forward(): every live record of the owner is retired (its
  // allocation reverted, its id never live again) and `drafts` (new records only) are created, in
  // one atomic step with commit's checks and guarantee. Returns the new FlowIds in draft order.
  static std::vector<FlowId> exchange(FlowLedger& ledger, std::int64_t owner, std::span<RecordDraft> drafts);
  // commit (retire false) and exchange (retire true).
  static std::vector<FlowId> commit_drafts(FlowLedger& ledger, std::int64_t owner, std::span<RecordDraft> drafts, bool retire);
  // One atomic exchange for several distinct owners (apply_actions): each owner's drafts with
  // commit's (retire false) or exchange's (retire true) rule, the aggregates checked once over
  // every owner, the new ids counted for the whole batch first and assigned in owner order, then
  // draft order. Throws (and changes nothing) as commit does. `created` (one per owner) receives
  // each owner's new FlowIds; the drafts' vectors are moved from.
  static void commit_batch(FlowLedger& ledger, std::span<OwnerDrafts> owners, std::span<std::vector<FlowId>> created);
  // A live record by id (throws std::invalid_argument for an unknown or retired id).
  [[nodiscard]] static const FlowRecord& live(const FlowLedger& ledger, std::int64_t id);
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
