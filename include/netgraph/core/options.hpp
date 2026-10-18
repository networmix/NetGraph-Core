/* The shared small option types: the selection axes of the searches and the policies and the
   placement kinds of the kernels. Component configurations and result statuses stay with their
   component. */
#pragma once

#include <variant>

#include "netgraph/core/types.hpp"

namespace netgraph::core {

// Forward follows edges from the roots; Reverse follows them backwards toward the roots.
// The returned DAG (a directed graph with no cycles) keeps the edges' original orientation.
enum class Direction { Forward, Reverse };

// Inclusive cost ceiling: Absolute gives a budget >= 0; Stretch multiplies the static shortest
// cost under the given Eligibility by a factor in [1, inf). Use the double's exact stored
// value, a dyadic rational (an integer divided by a power of two), with no comparison rounding.
struct Absolute {
  Cost cost {0};
  friend bool operator==(const Absolute&, const Absolute&) = default;
};
struct Stretch {
  double factor {1.0};
  friend bool operator==(const Stretch&, const Stretch&) = default;
};
using CostCeiling = std::variant<Absolute, Stretch>;

// Break equal-cost ties first between parallel edges, then between walks.
// For one (u, v) pair, Deterministic picks the lowest EdgeId; PreferHigherPriority picks the
// higher priority, then lowest EdgeId. For a single path (multipath = false),
// PreferHigherPriority first maximises the smallest edge priority along the walk (bottleneck);
// Deterministic ignores capacities. Both then prefer fewer hops, then the lexicographically
// smallest EdgeId sequence: compare at the first differing edge. Choose parallel members
// first (multi_edge = false), then choose the walk among the remaining members.
enum class EdgeTieBreak { Deterministic, PreferHigherPriority };

// Per-edge priority for PreferHigherPriority (single-path and parallel-edge choices), and for
// CapacityAware's Weighted membership and weights:
// Residual is unused capacity in units on the search's view; Capacity is static capacity in
// units; Headroom is the exact ratio residual / capacity, 0 on an edge of capacity 0 (rejected
// with CapacityAware's Weighted). Without a supplied residual, a search reads the static capacity as the unused
// capacity: Residual and Capacity read the capacity, Headroom 1 (0 on an edge of capacity 0).
enum class EdgePriority { Residual, Capacity, Headroom };

// `multipath` keeps all equal-cost walks as a shortest-path DAG; false selects the canonical
// single path. `multi_edge` keeps all eligible parallel edges between each (u, v) pair; false
// keeps one by (minimum cost, higher priority under PreferHigherPriority, lowest EdgeId).
// Every field has a default: Deterministic, Capacity, and both flags true. Invalid enum values
// throw std::invalid_argument naming the field. Fields remain assignable, so every entry point
// that uses a selection validates it again.
struct PathSelection {
  PathSelection(EdgeTieBreak tie_break_ = EdgeTieBreak::Deterministic, EdgePriority priority_ = EdgePriority::Capacity,
                bool multipath_ = true, bool multi_edge_ = true);
  EdgeTieBreak tie_break;
  EdgePriority priority;
  bool multipath;
  bool multi_edge;
  friend bool operator==(const PathSelection&, const PathSelection&) = default;
};

// A placement kernel assigns one flow's amounts to its usable route edges (its split graph).
// MaxFlow replaces its allocation with a minimum-cost flow on the split graph over
// V = W + own, where W is unused capacity and own is this flow's current allocation. Weighted
// splits at each node using weights fixed at creation: give each next unit to the edge with
// largest w / (a + 1), where w is its weight and a its share so far; ties use lower EdgeId.
// EqualBalanced splits equally at each node, giving leftover units to edges in ascending
// EdgeId order.
enum class FlowPlacement { MaxFlow, Weighted, EqualBalanced };

// Where the Weighted kernel's member weights come from (CostOnly Weighted and supplied split
// routes): Capacity, each member's static capacity; Explicit, one integer per edge given with the
// configuration (edge_weights; a positive weight on an edge of capacity 0 becomes 0);
// ResidualSnapshot, each member's free capacity read once when the routes are prepared, against
// one common view, and frozen for the epoch. A member of weight 0 is pruned, together with the
// branches it leaves without a destination, inside the route only (no longer bypass is created).
// The selection's priority never chooses weights.
enum class WeightBasis { Capacity, Explicit, ResidualSnapshot };

// Lossless admits only the largest total that the routes carry with no loss at all; Lossy
// forwards the whole finite offer hop by hop and records where it is lost (EqualBalanced and
// Weighted only).
enum class Admission { Lossless, Lossy };

// How a demand's traffic enters at its sources, separately from how it is routed. Equal: a total
// x enters as round robin over the roster (ascending NodeId, the first x mod H sources one unit
// more); a source whose share cannot go anywhere keeps it and loses it, it never moves to another
// source. Optimized: the controller chooses the integer injections (one virtual source).
enum class Origin { Equal, Optimized };

} // namespace netgraph::core
