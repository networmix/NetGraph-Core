/* The shared small option types: the tie rule of the searches and the policies, the route
   allocations of the kernels, their weight bases and the admissions. Component configurations and result statuses stay with their
   component. */
#pragma once

#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <variant>

#include "netgraph/core/types.hpp"

namespace netgraph::core {

// Forward follows edges from the roots; Reverse follows them backwards toward the roots.
// The returned DAG (a directed graph with no cycles) keeps the edges' original orientation.
enum class Direction { Forward, Reverse };

// Inclusive cost ceiling: AbsoluteCost gives a budget >= 0; CostStretch multiplies the static
// shortest cost under the given Eligibility by an exact factor n/d >= 1: a path of cost c
// qualifies when c * d <= n * base, compared exactly. Both are checked where they are used.
struct AbsoluteCost {
  Cost cost {0};
  friend bool operator==(const AbsoluteCost&, const AbsoluteCost&) = default;
};

// The factor as one reduced fraction numerator / denominator with a denominator in [1, 2^63);
// any numerator is a valid value (below 1 is refused where the ceiling is used).
class CostStretch {
 public:
  // Reduces numerator / denominator; a denominator <= 0 throws std::invalid_argument.
  explicit CostStretch(std::int64_t numerator, std::int64_t denominator = 1) : numerator_(numerator), denominator_(denominator) {
    if (denominator < 1) throw std::invalid_argument("CostStretch fraction denominator must be in [1, 2^63)");
    // The magnitude as unsigned: -2^63 has no int64 negation.
    const std::uint64_t magnitude =
        numerator < 0 ? std::uint64_t {0} - static_cast<std::uint64_t>(numerator) : static_cast<std::uint64_t>(numerator);
    const std::uint64_t g = std::gcd(magnitude, static_cast<std::uint64_t>(denominator));
    if (g > 1) {
      // g divides the denominator (< 2^63), so g < 2^63 and the quotients fit int64.
      numerator_ = numerator < 0 ? -static_cast<std::int64_t>(magnitude / g) : static_cast<std::int64_t>(magnitude / g);
      denominator_ = static_cast<std::int64_t>(static_cast<std::uint64_t>(denominator) / g);
    }
  }
  [[nodiscard]] std::int64_t numerator() const noexcept { return numerator_; }
  [[nodiscard]] std::int64_t denominator() const noexcept { return denominator_; }
  friend bool operator==(const CostStretch&, const CostStretch&) = default;

 private:
  std::int64_t numerator_;
  std::int64_t denominator_;
};
using CostCeiling = std::variant<AbsoluteCost, CostStretch>;

// How a route choice breaks a cost tie; cost always comes first. FewerHops: then fewer hops, then the
// lexicographically smallest EdgeId sequence (compared at the first differing edge). The other three
// first prefer the larger bottleneck (the smallest value along the path) of an edge score, then
// fewer hops, then the EdgeId sequence: HigherCapacity scores static capacity, HigherAvailableCapacity
// the availability of the search's view, HigherHeadroom availability / static capacity as an exact
// ratio (0 on an edge of static capacity 0). Availability is the unreserved capacity for
// DistributedTE's CSPF, the service view for ShortestPaths, the supplied residual for spf and ksp,
// and the static capacity where no view is supplied (HigherHeadroom then reads 1 on a positive edge).
// The same order chooses among parallel edges of one (u, v) pair when parallel_edges is false (the
// cheapest member, then the higher score, then the lowest EdgeId). The tie rule never supplies
// weights and never enables a controller objective (CentralizedTE::prefer_headroom is T2's).
enum class RouteTieBreak { FewerHops, HigherCapacity, HigherAvailableCapacity, HigherHeadroom };

// A route allocation (a kernel) assigns one flow's amounts to its usable route edges (its split graph).
// Equal splits equally at each node, giving leftover units to edges in ascending EdgeId order.
// Weighted splits at each node using weights fixed at creation: give each next unit to the edge
// with largest w / (a + 1), where w is its weight and a its share so far; ties use lower EdgeId.
// Optimized replaces its allocation with a minimum-cost flow on the split graph over
// V = W + own, where W is unused capacity and own is this flow's current allocation.
enum class RouteAllocation { Equal, Weighted, Optimized };

// Where the Weighted kernel's member weights come from (ShortestPaths Weighted and supplied split
// routes): Capacity, each member's static capacity; Explicit, one integer per edge given with the
// configuration (edge_weights; a positive weight on an edge of capacity 0 becomes 0);
// ResidualSnapshot, each member's free capacity read once when the routes are prepared, against
// one common view, and frozen for the epoch. A member of weight 0 is pruned, together with the
// branches it leaves without a destination, inside the route only (no longer bypass is created).
// The tie rule's edge score never chooses weights.
enum class WeightBasis { Capacity, Explicit, ResidualSnapshot };

// Lossless admits only the largest total that the routes carry with no loss at all; Lossy
// forwards the whole finite offer hop by hop and records where it is lost (every model except the
// optimizing ones: ShortestPaths Optimized and CentralizedTE without a tunnel count).
enum class Admission { Lossless, Lossy };

// How a demand's traffic enters at its sources, separately from how it is routed. Equal: a total
// x enters as round robin over the roster (ascending NodeId, the first x mod H sources one unit
// more); a source whose share cannot go anywhere keeps it and loses it, it never moves to another
// source. Optimized: the controller chooses the integer injections (one virtual source).
// Weighted: CentralizedTE without tunnels prescribes H(total, source_weights) instead of RR.
enum class SourceAllocation { Equal, Optimized, Weighted };

} // namespace netgraph::core
