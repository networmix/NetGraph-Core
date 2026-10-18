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
// at the commit, in draft order); the commit sets the record's flow_id and owner. The record's
// `edges` and `drops` ascend by EdgeId; zero entries are dropped.
struct Draft {
  std::optional<std::int64_t> id {};
  FlowRecord record {};
};

// The saturation diagnostic of FlowLedger::bottleneck on a given residual (the sources already a
// validated set): the edges (statically eligible, capacity >= 1 unit) leaving the set reachable
// from the allowed sources over eligible forward arcs with residual >= 1 unit, ascending EdgeId.
[[nodiscard]] std::vector<EdgeId> bottleneck_of(const Graph& g, std::span<const Amount> residual,
                                                std::span<const NodeId> sources, const Eligibility* elig);

struct LedgerAccess {
  // A new owner id (0, 1, ... per ledger). Throws on a moved-from ledger.
  [[nodiscard]] static std::int64_t bind(FlowLedger& fg);
  // The owner's policy is destroyed without remove(): its records stay as orphans.
  static void release(FlowLedger& fg, std::int64_t owner) noexcept;
  // remove(): the owner's records are deleted (allocations reverted), their ids retired.
  // Throws only before any change (a moved-from ledger, an unbound owner).
  static std::vector<std::int64_t> detach(FlowLedger& fg, std::int64_t owner);
  // One atomic exchange: `drafts` must carry every live record of the owner (by id) plus the
  // new ones; the new aggregate must stay within [0, baseline] on every edge. Throws (and
  // changes nothing) otherwise, an allocation failure included: every allocation happens
  // before a publication step that cannot throw. Returns the FlowIds of the new records in
  // draft order. The drafts' vectors are moved from.
  static std::vector<FlowId> commit(FlowLedger& fg, std::int64_t owner, std::span<Draft> drafts);
  // The exchange of replace() and forward(): every live record of the owner is retired (its
  // allocation reverted, its id never live again) and `drafts` (new records only) are created, in
  // one atomic step with commit's checks and guarantee. Returns the new FlowIds in draft order.
  static std::vector<FlowId> exchange(FlowLedger& fg, std::int64_t owner, std::span<Draft> drafts);
  // commit (retire false) and exchange (retire true).
  static std::vector<FlowId> commit_drafts(FlowLedger& fg, std::int64_t owner, std::span<Draft> drafts, bool retire);
  // A live record by id (throws std::invalid_argument for an unknown or retired id).
  [[nodiscard]] static const FlowRecord& live(const FlowLedger& fg, std::int64_t id);
  // The working residual baseline - allocated into `out` (resized to E).
  static void residual_into(const FlowLedger& fg, std::vector<Amount>& out);
  [[nodiscard]] static std::span<const Amount> baseline(const FlowLedger& fg);
  // The aggregate allocation (borrowed; valid until the next commit, detach or reset).
  [[nodiscard]] static std::span<const Amount> allocated(const FlowLedger& fg);
  // Diagnostics (the tests and the private bindings).
  [[nodiscard]] static std::int64_t next_flow_id(const FlowLedger& fg);
  [[nodiscard]] static std::vector<std::int64_t> bound_owners(const FlowLedger& fg);
  [[nodiscard]] static std::vector<std::int64_t> orphan_owners(const FlowLedger& fg);
};

} // namespace netgraph::core::detail
