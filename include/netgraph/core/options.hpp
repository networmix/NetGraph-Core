/* Route selection axes shared by the searches and the policies (APPROACH 6). */
#pragma once

namespace netgraph::core {

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
// outer axis is explicit); the two flags default to the plain form of their mechanism.
struct PathSelection {
  PathSelection(EdgeTieBreak tie_break_, EdgePriority priority_, bool multipath_ = true, bool multi_edge_ = true)
      : tie_break(tie_break_), priority(priority_), multipath(multipath_), multi_edge(multi_edge_) {}
  EdgeTieBreak tie_break;
  EdgePriority priority;
  bool multipath;
  bool multi_edge;
  friend bool operator==(const PathSelection&, const PathSelection&) = default;
};

} // namespace netgraph::core
