/* Flow identities, the 128-bit counter and the ledger records. */
#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "netgraph/core/options.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// (ledger uid, id): ids are assigned by the ledger in creation order and never reused.
using FlowId = std::pair<std::uint64_t, std::int64_t>;

// A cumulative counter (offered, per-edge drops): a sum over calls, the 128-bit value
// hi * 2^64 + lo, always below 2^127. Such sums over many calls may exceed 64 bits; the
// Python binding returns them as exact int.
struct Counter128 {
  std::uint64_t hi {0};
  std::uint64_t lo {0};
  friend bool operator==(const Counter128&, const Counter128&) = default;
  [[nodiscard]] bool is_zero() const noexcept { return hi == 0 && lo == 0; }
  // The sum with a non-negative amount; false (and no change) when the sum reaches 2^127.
  bool add(Amount x) noexcept {
    const auto ux = static_cast<std::uint64_t>(x);
    const std::uint64_t lo2 = lo + ux;
    const std::uint64_t hi2 = hi + (lo2 < lo ? 1U : 0U);
    if (hi2 >= (std::uint64_t{1} << 63)) return false;
    hi = hi2;
    lo = lo2;
    return true;
  }
};

// A snapshot of one ledger flow, returned by records() or record(id). `route` is its route
// structure; `min_route_cost` is the minimum walk cost from its sources to its destinations over
// that route, or kUnreachable if none exists. `entry_id` is the program entry the record carries
// (T1, T2 and Pinned: one record per entry, the entry identifier of its program); absent for the
// records kept per roster source (ShortestPaths, CentralizedTE). Lossy records alone have `offered`: the
// cumulative amount submitted to this flow before any drops. `edge_flows` (carried load), `edge_drops`
// and `edge_reservations` (the strict holds of an up T1 tunnel: its request on every link of its path)
// contain positive values in ascending EdgeId; per-call and cumulative drop totals are derived.
struct FlowRecord {
  FlowId flow_id {};
  std::int64_t policy_id {0};
  std::vector<NodeId> sources {};
  std::vector<NodeId> destinations {};
  Admission admission {Admission::Lossless};
  std::optional<std::int32_t> entry_id {};
  std::shared_ptr<const PredecessorDAG> route {};
  Cost min_route_cost {kUnreachable};
  std::optional<Counter128> offered {};
  Amount delivered {0};
  std::vector<std::pair<EdgeId, Amount>> edge_flows {};
  std::vector<std::pair<EdgeId, Counter128>> edge_drops {};
  std::vector<std::pair<EdgeId, Amount>> edge_reservations {};
};

} // namespace netgraph::core
