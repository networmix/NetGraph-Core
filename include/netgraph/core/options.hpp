/* The shared small option types: the selection axes of the searches and the policies (APPROACH 6)
   and the placement kinds of the kernels (APPROACH 7). Component configurations and result
   statuses stay with their component. */
#pragma once

#include <variant>

#include "netgraph/core/types.hpp"

namespace netgraph::core {

// The orientation of a search: Forward from the roots, or Reverse (rooted at the roots over
// reversed arcs; the DAG it returns is forward-oriented, toward the roots).
enum class Direction { Forward, Reverse };

// One cost ceiling knob (APPROACH 6): an Absolute cost budget (inclusive, >= 0), or a Stretch
// factor (in [1, inf)) of the static shortest cost under the given Eligibility, read as the exact
// dyadic rational it is and compared without rounding.
struct Absolute {
  Cost cost {0};
  friend bool operator==(const Absolute&, const Absolute&) = default;
};
struct Stretch {
  double factor {1.0};
  friend bool operator==(const Stretch&, const Stretch&) = default;
};
using CostCeiling = std::variant<Absolute, Stretch>;

// How ties between equal-cost candidates are broken, at both tie levels: the member chosen
// among parallel edges of one (u, v) pair and the single path chosen among equal-cost walks.
// Deterministic ignores capacities entirely (the lowest EdgeId wins); PreferHigherPriority
// prefers the higher `priority` value first (then the lowest EdgeId).
enum class EdgeTieBreak { Deterministic, PreferHigherPriority };

// The per-edge value read by PreferHigherPriority and, for the Weighted kernel, its membership
// and weights: the residual units on the view of the search (Residual), the static capacity
// units (Capacity), or the exact rational residual / capacity (Headroom; rejected with
// Weighted). A search given no residual reads static capacities for every priority.
enum class EdgePriority { Residual, Capacity, Headroom };

// The selection of a search: `multipath` keeps every equal-cost walk (the shortest-path DAG)
// instead of the one canonical single path; `multi_edge` keeps every eligible parallel member
// of a (u, v) pair instead of the one chosen by (minimum cost, priority under
// PreferHigherPriority, lowest EdgeId). tie_break and priority have no default (APPROACH 9: an
// outer axis is explicit); the two flags default to the plain form of their mechanism. The
// constructor rejects a tie_break or priority outside its enumeration (std::invalid_argument
// naming the field; rule 13); the fields stay assignable, so every entry that reads a selection
// checks it again.
struct PathSelection {
  PathSelection(EdgeTieBreak tie_break_, EdgePriority priority_, bool multipath_ = true, bool multi_edge_ = true);
  EdgeTieBreak tie_break;
  EdgePriority priority;
  bool multipath;
  bool multi_edge;
  friend bool operator==(const PathSelection&, const PathSelection&) = default;
};

// The kernel that computes one flow's pushes on its DAG: MaxFlow (the flow core of the
// minimum-cost flow restricted to the flow's split graph, a replacement on V = W + own),
// Weighted (hop-by-hop WCMP: highest-averages apportionment over frozen weights at every node)
// and EqualBalanced (ECMP: round robin at every node).
enum class FlowPlacement { MaxFlow, Weighted, EqualBalanced };

// Lossless admits only what fits; Lossy forwards hop by hop without admission control and
// drops what does not fit on each edge (EqualBalanced and Weighted only).
enum class Admission { Lossless, Lossy };

} // namespace netgraph::core
