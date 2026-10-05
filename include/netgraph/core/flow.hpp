/* Placement kinds, flow identities and ledger records (APPROACH 7 and 8). */
#pragma once

#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// The kernel that computes one flow's pushes on its DAG: MaxFlow (the flow core of the
// minimum-cost flow restricted to the flow's split graph, a replacement on V = W + own),
// Weighted (hop-by-hop WCMP: highest-averages apportionment over frozen weights at every node)
// and EqualBalanced (ECMP: round robin at every node).
enum class FlowPlacement { MaxFlow, Weighted, EqualBalanced };

// Lossless admits only what fits; Lossy forwards hop by hop without admission control and
// drops what does not fit on each edge (EqualBalanced and Weighted only).
enum class Admission { Lossless, Lossy };

// The path-cost sentinel of an unreachable node or an empty route (None in Python).
inline constexpr Cost kUnreachable = std::numeric_limits<Cost>::max();

// (ledger uid, id): ids are assigned by the ledger in creation order and never reused.
using FlowId = std::pair<std::uint64_t, std::int64_t>;

// A cumulative counter (offered, per-edge drops): a sum over calls, kept in 128 bits as
// hi * 2^64 + lo, always below 2^127 (APPROACH 2.1, 2.3).
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

// One flow of the ledger, as a snapshot (records(), record(id)). `bundle` is the flow's route
// structure (its structural identity); `route_cost` the cached minimum walk cost from the
// sources to the destinations over that bundle (kUnreachable when it has no such walk);
// `offered` exists on Lossy records only (the injection basis); `edges` and `drops` list
// positive values in ascending EdgeId; per-call and cumulative drop totals are derived.
struct FlowRecord {
  FlowId flow_id {};
  std::int64_t owner {0};
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

} // namespace netgraph::core
