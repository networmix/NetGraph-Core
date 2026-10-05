/* The search semantics of APPROACH 6 as internal routines (the public spf and KSP are 1d).

   Positive costs make the shortest-path DAG from a root set unique (2.4): the tight edges
   d(u) + cost(e) = d(v) of the member graph. With targets the DAG is pruned to the walks that
   end in a nearest target (anycast), the same edge set for a forward search from the roots
   and a reverse search rooted at the targets. The single path is the staged canonical walk:
   (1) the tight DAG of the query, (2) under PreferHigherPriority the best attainable
   bottleneck B* (a max-min DP on the DAG, capped by a Yen prefix's bottleneck) and the DAG
   restricted to priority >= B*, (3) one reverse hop count h from the targets and a greedy
   walk taking the smallest tight edge whose head has h(v) = h(u) - 1.
*/
#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

enum class Direction { Forward, Reverse };

// An exact priority value: num / den with den >= 1 (den == 1 for Residual and Capacity, the
// capacity for Headroom); `infinite` above every finite value (the bottleneck of an empty walk).
struct Priority {
  Amount num {0};
  Amount den {1};
  bool infinite {false};
};
[[nodiscard]] int compare(const Priority& a, const Priority& b) noexcept;  // -1, 0, 1

// The per-edge priority values a search reads: `residual` empty (no residual) reads the static
// capacity for every priority; else Residual = residual units, Capacity = capacity units,
// Headroom = residual / capacity.
struct PriorityValues {
  const StrictMultiDiGraph* g {nullptr};
  EdgePriority priority {EdgePriority::Capacity};
  std::span<const Amount> residual {};
  [[nodiscard]] Priority of(EdgeId e) const noexcept;
};

// The member graph (flags per edge): statically eligible, capacity >= 1 unit, residual >=
// min_residual when a residual is given; with multi_edge = false one member per (u, v): the
// minimum cost, then under PreferHigherPriority the highest priority value, then the lowest
// EdgeId (select_adjacency).
[[nodiscard]] std::vector<std::uint8_t> member_mask(const StrictMultiDiGraph& g, const Eligibility* elig,
                                                    const PathSelection& selection,
                                                    std::span<const Amount> residual, Amount min_residual);

// Dijkstra over the members from (Forward) or toward (Reverse) the root set; kUnreachable for
// the nodes it does not reach.
[[nodiscard]] std::vector<Cost> distances(const StrictMultiDiGraph& g, std::span<const std::uint8_t> member,
                                          std::span<const NodeId> roots, Direction direction);

// The tight edges among the reached nodes of a full search (forward-oriented), ascending.
[[nodiscard]] std::vector<EdgeId> full_dag(const StrictMultiDiGraph& g, std::span<const std::uint8_t> member,
                                           std::span<const Cost> dist, Direction direction);

// The edges on the minimum-cost walks from the source set to the destination set (the walks
// of cost d* = min over destinations of the distance), ascending; cost kUnreachable and no
// edge when no destination is reachable.
struct TightDag {
  Cost cost {0};
  std::vector<EdgeId> edges {};
};
[[nodiscard]] TightDag tight_dag(const StrictMultiDiGraph& g, std::span<const std::uint8_t> member,
                                 std::span<const NodeId> sources, std::span<const NodeId> destinations);

// The staged canonical single path from the source set to the destination set over the
// members; `priority` null for Deterministic; `prefix_cap` caps B* (a Yen root prefix).
struct SinglePath {
  std::vector<EdgeId> walk {};  // in walk order
  Cost cost {0};
  std::optional<Priority> bottleneck {};  // under PreferHigherPriority
  std::optional<Priority> b_star {};
};
[[nodiscard]] std::optional<SinglePath> canonical_path(const StrictMultiDiGraph& g,
                                                       std::span<const std::uint8_t> member,
                                                       std::span<const NodeId> sources,
                                                       std::span<const NodeId> destinations,
                                                       const PriorityValues* priority,
                                                       std::optional<Priority> prefix_cap = std::nullopt);

// The search of APPROACH 6 (sub-PR 1d exposes it as spf()): validates (std::invalid_argument)
// the root set (non-empty, ids in range; duplicates removed), the targets when given
// (non-empty, ids in range, disjoint from the roots), multipath = false without targets, the
// residual (length E, entries in [0, capacity_units]), min_residual >= 1 and a foreign
// Eligibility, before any work. Excluded roots and targets are skipped (all excluded: an
// empty result). Without targets: `costs` for every node and the full DAG; with targets: the
// anycast-pruned DAG (or the canonical single path for multipath = false) and `target_cost`
// (nullopt when no target is reachable).
struct SpfOut {
  std::optional<std::vector<Cost>> costs {};
  std::vector<EdgeId> edges {};
  std::optional<Cost> target_cost {};
  std::optional<SinglePath> path {};
};
[[nodiscard]] SpfOut spf(const StrictMultiDiGraph& g, std::span<const NodeId> roots, Direction direction,
                         std::optional<std::span<const NodeId>> targets, const PathSelection& selection,
                         std::optional<std::span<const Amount>> residual, Amount min_residual,
                         const Eligibility* elig);

// The sorted, duplicate-free node set; throws std::invalid_argument for an empty set or an id
// out of range (`name` labels the message).
[[nodiscard]] std::vector<NodeId> node_set(std::span<const NodeId> nodes, std::int32_t num_nodes, const char* name);

// Throws std::invalid_argument when the two sorted sets share a vertex (rule 13).
void check_disjoint(std::span<const NodeId> a, std::span<const NodeId> b);

} // namespace netgraph::core::detail
