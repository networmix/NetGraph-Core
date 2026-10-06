/* The search semantics of APPROACH 6 as internal routines (the DagCache's supports, the
   policies' routes and the private _spf entry of the bindings).

   Positive costs make the shortest-path DAG from a root set unique (2.4): the tight edges
   d(u) + cost(e) = d(v) of the member graph. With targets the DAG is pruned to the walks that
   end in a nearest target (anycast), the same edge set for a forward search from the roots
   and a reverse search rooted at the targets. The single path is the staged canonical walk:
   (1) the tight DAG of the query, (2) under PreferHigherPriority the best attainable
   bottleneck B* (a max-min DP on the DAG) and the DAG restricted to priority >= B*, (3) one
   reverse hop count h from the targets and a greedy walk taking the smallest tight edge whose
   head has h(v) = h(u) - 1.
*/
#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

// An exact priority value: num / den with den >= 1 (den == 1 for Residual and Capacity, the
// capacity for Headroom); `infinite` above every finite value (the bottleneck of an empty walk).
struct Priority {
  Amount num {0};
  Amount den {1};
  bool infinite {false};
};
[[nodiscard]] int compare(const Priority& a, const Priority& b) noexcept;  // -1, 0, 1

// The priority value of an edge with `residual` units left of `capacity`: the residual units
// (Residual), the capacity units (Capacity) or the rational residual / capacity (Headroom).
[[nodiscard]] Priority priority_value(EdgePriority priority, Amount residual, Amount capacity) noexcept;

// The per-edge priority values a search reads: `residual` empty (no residual) reads the static
// capacity for every priority; else priority_value on the residual.
struct PriorityValues {
  const Graph* g {nullptr};
  EdgePriority priority {EdgePriority::Capacity};
  std::span<const Amount> residual {};
  [[nodiscard]] Priority of(EdgeId e) const noexcept;
};

// select_adjacency (APPROACH 6): whether member `e` replaces `cur`, the choice so far among the
// eligible members of one (u, v) pair seen in ascending EdgeId (cur < e): the lower cost, then
// under PreferHigherPriority (`prefer`) the higher priority value; on a tie the lower EdgeId
// (cur) stays. `value(e)` reads the priority value of an edge on the view of the selection.
template <class Value>
[[nodiscard]] bool replaces(const Graph& g, EdgeId e, EdgeId cur, bool prefer, const Value& value) {
  const auto cost = g.cost();
  const Cost ce = cost[static_cast<std::size_t>(e)];
  const Cost cc = cost[static_cast<std::size_t>(cur)];
  if (ce != cc) return ce < cc;
  return prefer && compare(value(e), value(cur)) > 0;
}

// The member graph (flags per edge): statically eligible, capacity >= 1 unit, residual >=
// min_residual when a residual is given; with multi_edge = false one member per (u, v) by
// select_adjacency.
[[nodiscard]] std::vector<std::uint8_t> member_mask(const Graph& g, const Eligibility* elig,
                                                    const PathSelection& selection,
                                                    std::span<const Amount> residual, Amount min_residual);

// Dijkstra over the members from (Forward) or toward (Reverse) the root set; kUnreachable for
// the nodes it does not reach.
[[nodiscard]] std::vector<Cost> distances(const Graph& g, std::span<const std::uint8_t> member,
                                          std::span<const NodeId> roots, Direction direction);

// The tight edges among the reached nodes of a full search (forward-oriented), ascending.
[[nodiscard]] std::vector<EdgeId> full_dag(const Graph& g, std::span<const std::uint8_t> member,
                                           std::span<const Cost> dist, Direction direction);

// The edges on the minimum-cost walks from the source set to the destination set (the walks
// of cost d* = min over destinations of the distance), ascending; cost kUnreachable and no
// edge when no destination is reachable.
struct TightDag {
  Cost cost {0};
  std::vector<EdgeId> edges {};
};
[[nodiscard]] TightDag tight_dag(const Graph& g, std::span<const std::uint8_t> member,
                                 std::span<const NodeId> sources, std::span<const NodeId> destinations);

// The staged canonical single path from the source set to the destination set over the
// members; `priority` null for Deterministic.
struct SinglePath {
  std::vector<EdgeId> walk {};  // in walk order
  Cost cost {0};
};
[[nodiscard]] std::optional<SinglePath> canonical_path(const Graph& g,
                                                       std::span<const std::uint8_t> member,
                                                       std::span<const NodeId> sources,
                                                       std::span<const NodeId> destinations,
                                                       const PriorityValues* priority);

// The search of APPROACH 6: validates (std::invalid_argument) the direction and the selection
// (enumeration values), the root set (non-empty, ids in range; duplicates removed), the targets
// when given (non-empty, ids in range, disjoint from the roots), multipath = false without
// targets, the residual (length E, entries in [0, capacity_units]), min_residual >= 1 and a
// foreign Eligibility, before any work. Excluded roots and targets are skipped (all excluded:
// an empty result). Without targets: `costs` for every node and the full DAG; with targets: the
// anycast-pruned DAG (or the canonical single path for multipath = false) and `target_cost`
// (nullopt when no target is reachable).
struct SpfOut {
  std::optional<std::vector<Cost>> costs {};
  std::vector<EdgeId> edges {};
  std::optional<Cost> target_cost {};
  std::optional<SinglePath> path {};
};
[[nodiscard]] SpfOut spf(const Graph& g, std::span<const NodeId> roots, Direction direction,
                         std::optional<std::span<const NodeId>> targets, const PathSelection& selection,
                         std::optional<std::span<const Amount>> residual, Amount min_residual,
                         const Eligibility* elig);

// The full static support of a search from (Forward) or toward (Reverse) the root set: every
// tight edge of the statically eligible member graph (multipath, multi_edge, no residual) over
// the nodes it reaches, forward-oriented (APPROACH 10: the DagCache's value, and a CostOnly
// policy's DAG when it has no cache).
[[nodiscard]] PredDAG static_support(const Graph& g, std::span<const NodeId> roots, Direction direction,
                                     const Eligibility* elig);

} // namespace netgraph::core::detail
