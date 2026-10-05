/* The ledger's internal entry points: what a policy uses to place on and remove from it.

   APPROACH 8: the ledger is mutated only through policies. A policy binds (an owner id), then
   every call ends with one atomic exchange of all of its records (`commit`); remove() reverts
   its allocations and retires its ids (`detach`); a policy destroyed without remove() leaves
   its records as orphans (`release`). None of this is public API.
*/
#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/flow.hpp"
#include "netgraph/core/flow_graph.hpp"

namespace netgraph::core::detail {

// A record proposed by a policy for the commit: `id` nullopt is a new flow (its id is assigned
// at the commit, in draft order). `edges` and `drops` ascending EdgeId; zero entries dropped.
struct Draft {
  std::optional<std::int64_t> id {};
  std::vector<NodeId> sources {};
  std::vector<NodeId> destinations {};
  FlowPlacement placement {FlowPlacement::MaxFlow};
  Admission admission {Admission::Lossless};
  std::shared_ptr<const PredDAG> bundle {};
  Cost route_cost {kUnreachable};
  std::optional<Counter128> offered {};
  Amount delivered {0};
  std::vector<std::pair<EdgeId, Amount>> edges {};
  std::vector<std::pair<EdgeId, Counter128>> drops {};
};

struct LedgerAccess {
  // A new owner id (0, 1, ... per ledger). Throws on a moved-from ledger.
  [[nodiscard]] static std::int64_t bind(FlowGraph& fg);
  // The owner's policy is destroyed without remove(): its records stay as orphans.
  static void release(FlowGraph& fg, std::int64_t owner) noexcept;
  // remove(): the owner's records are deleted (allocations reverted), their ids retired.
  // Throws only before any change (a moved-from ledger, an unbound owner).
  static std::vector<std::int64_t> detach(FlowGraph& fg, std::int64_t owner);
  // One atomic exchange: `drafts` must carry every live record of the owner (by id) plus the
  // new ones; the new aggregate must stay within [0, baseline] on every edge. Throws (and
  // changes nothing) otherwise, an allocation failure included: every allocation happens
  // before a publication step that cannot throw. Returns the FlowIds of the new records in
  // draft order. The drafts' vectors are moved from.
  static std::vector<FlowId> commit(FlowGraph& fg, std::int64_t owner, std::span<Draft> drafts);
  // A live record by id (throws std::invalid_argument for an unknown or retired id).
  [[nodiscard]] static const FlowRecord& live(const FlowGraph& fg, std::int64_t id);
  // The working residual baseline - allocated into `out` (resized to E).
  static void residual_into(const FlowGraph& fg, std::vector<Amount>& out);
  [[nodiscard]] static std::span<const Amount> baseline(const FlowGraph& fg);
  // The aggregate allocation (borrowed; valid until the next commit, detach or reset).
  [[nodiscard]] static std::span<const Amount> allocated(const FlowGraph& fg);
  // Diagnostics (tests, the private recorder binding).
  [[nodiscard]] static std::int64_t next_flow_id(const FlowGraph& fg);
  [[nodiscard]] static std::vector<std::int64_t> bound_owners(const FlowGraph& fg);
  [[nodiscard]] static std::vector<std::int64_t> orphan_owners(const FlowGraph& fg);
};

} // namespace netgraph::core::detail
