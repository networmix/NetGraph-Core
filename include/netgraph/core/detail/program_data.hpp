/* The inside of a ForwardingProgram (internal): its entries, split contexts, injection and replay
   basis, the builders the models use, and the comparisons of the change flags. The public value
   (forwarding_program.hpp) shows only the endpoints, the tunnels and the source weights. */
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/detail/route_kernels.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/forwarding_program.hpp"

namespace netgraph::core::detail {

// A tunnel's path in walk order (empty for a tunnel that is down).
[[nodiscard]] inline std::span<const EdgeId> path_of(const Tunnel& t) noexcept {
  return t.edges.has_value() ? std::span<const EdgeId>(*t.edges) : std::span<const EdgeId> {};
}

// One entry: a tunnel (tunnels()[tunnel_position]) or a split route (a root of
// split_contexts[split_context_position]). `roster_index` is the head's position in the roster;
// `route` the entry's route as a DAG (a split route's edges before weight pruning; null for a
// tunnel, whose path the record's route is made from when a Pinned policy needs it) and
// `min_route_cost` its minimum walk cost.
struct ProgramEntry {
  std::int32_t entry_id {0};
  NodeId head {0};
  std::int32_t roster_index {0};
  std::int32_t local_index {0};
  Amount weight {0};
  std::int32_t tunnel_position {-1};
  std::int32_t split_context_position {-1};
  std::shared_ptr<const PredecessorDAG> route {};
  Cost min_route_cost {kUnreachable};
};

// A split context: one plan whose roots are the heads it serves (a ShortestPaths program: the whole
// roster, pooled; a supplied split route: its head alone), and per root the entry it carries.
// `route` (ascending) is the context's edges before weight pruning, the replay basis of its
// members.
struct SplitContext {
  std::shared_ptr<const PreparedRoute> prepared_route {};  // its route: the context's edges before weight pruning
  // Per plan root: its entry's position in entries (every root has one), in the program's block
  // (initialize_storage).
  std::span<const std::int32_t> source_entry_positions {};
};

struct ProgramData {
  ProgramData() = default;
  // The lists below are views of the data's own block (initialize_storage): a copy would read the
  // original's block, so there is none; a move keeps the block and with it the views.
  ProgramData(const ProgramData&) = delete;
  ProgramData& operator=(const ProgramData&) = delete;
  ProgramData(ProgramData&&) noexcept = default;
  ProgramData& operator=(ProgramData&&) noexcept = default;
  ~ProgramData() = default;

  std::uint64_t graph_uid {0};
  std::span<const NodeId> sources {};       // ascending
  std::span<const NodeId> destinations {};  // ascending
  std::vector<Tunnel> tunnels {};
  std::optional<std::vector<std::pair<NodeId, Amount>>> source_weights {};
  // The injection: the roster (the heads with an entry, ascending; T1, T2: the policy's roster), its
  // captured weights (empty: equal sources), every entry in ascending index, per roster slot its
  // entries (ascending index).
  std::span<const NodeId> reachable_sources {};
  std::vector<Amount> source_weight_values {};
  std::vector<ProgramEntry> entries {};
  // Per roster slot its entries (entries_for_source): slot s's are source_entry_positions[source_entry_offsets[s]..source_entry_offsets[s+1]),
  // ascending index; both empty for a pooled program.
  std::span<const std::int32_t> source_entry_offsets {};
  std::span<const std::int32_t> source_entry_positions {};
  std::span<const std::int32_t> tunnel_entry_positions {};  // per tunnel: its entry
  std::vector<SplitContext> split_contexts {};
  // The replay basis: the topology under the producer's exclusions (topology: unrestricted
  // CentralizedTE), or per roster slot the edges of its original routes (basis_parts).
  bool uses_topology_basis {false};
  std::span<const NodeId> excluded_nodes {};  // ascending
  std::span<const EdgeId> excluded_edges {};  // ascending (an edge whose endpoint is excluded included)
  std::shared_ptr<const std::vector<EdgeId>> common_replay_support {};  // every slot's (a supported CentralizedTE)
  std::vector<std::shared_ptr<const PredecessorDAG>> source_supports {};    // per slot its route (ShortestPaths Optimized)
  // A slot's replay basis as edge sets, appended to `out` (each ascending, except a tunnel's path,
  // in walk order: `second` false): the shared route or edge set above, else (supplied and ShortestPaths
  // split programs) the slot's own entries' routes: a tunnel's edges, a split route's route DAG (its
  // edges before weight pruning). Nothing is copied.
  void basis_parts(std::size_t slot, std::vector<std::pair<std::span<const EdgeId>, bool>>& out) const;
  // One pooled context carries every entry: a single context whose roots are the roster, one entry
  // of positive weight per head, equal sources (a ShortestPaths program; a supplied program of one
  // split route of positive weight). Its forwarding is the plan's own kernels.
  bool pools_source_inflows {false};
  // A slot's entries, ascending index: a pooled program's slot carries the one entry its context
  // root carries (no per-slot lists are kept).
  [[nodiscard]] std::span<const std::int32_t> entries_for_source(std::size_t slot) const noexcept {
    if (pools_source_inflows) return {&split_contexts.front().source_entry_positions[slot], 1};
    const auto first = static_cast<std::size_t>(source_entry_offsets[slot]);
    return {source_entry_positions.data() + first, static_cast<std::size_t>(source_entry_offsets[slot + 1]) - first};
  }
  // The last step of every producer, after its entries, tunnels, contexts and `pooled` are set:
  // the lists above copied into one block (one allocation; a failure leaves the data without
  // lists, to be discarded), the roster (given, or the entries' heads, ascending and distinct),
  // each entry's slot and local index (its position among its slot's entries), each tunnel's
  // local index, the per-slot lists (a stable counting sort by slot), the per-tunnel entries and
  // each context's root entries (the entry whose head is the root).
  void initialize_storage(std::span<const NodeId> sources, std::span<const NodeId> destinations,
               std::optional<std::span<const NodeId>> reachable_sources, std::span<const NodeId> excluded_nodes = {},
               std::span<const EdgeId> excluded_edges = {});
  // The block the lists view (the ownership tests).
  [[nodiscard]] std::span<const std::int32_t> block() const noexcept { return store_; }
  // Some entry has a path: an up tunnel or a split route (the condition of termination rule (3);
  // a program of down tunnels only, or without entries, has none).
  [[nodiscard]] bool has_path() const noexcept {
    for (const auto& entry : entries) {
      if (entry.tunnel_position < 0 || tunnels[static_cast<std::size_t>(entry.tunnel_position)].edges.has_value()) return true;
    }
    return false;
  }

private:
  std::vector<std::int32_t> store_ {};
};

// A program and its data in one allocation (ProgramAccess::make).
struct ProgramBlock {
  explicit ProgramBlock(ProgramData&& d) : data(std::move(d)), program(&data) {}
  ProgramData data;
  ForwardingProgram program;
};

struct ProgramAccess {
  [[nodiscard]] static const ProgramData& data(const ForwardingProgram& p) noexcept { return *p.data_; }
  [[nodiscard]] static std::shared_ptr<const ForwardingProgram> make(ProgramData&& data);
  // Endpoints/options/entries have all passed their ordered checks. Construction
  // does not repeat route validation; native and Python entries share this builder.
  [[nodiscard]] static std::shared_ptr<const ForwardingProgram> from_checked_routes(
      std::span<const NodeId> sources, std::span<const NodeId> destinations,
      std::span<const RouteDefinition> routes, ProgramData&& data);
};

// Check each supplied entry in order, preparing its private storage in the same
// pass. A later refusal discards it all. Endpoints/options have already passed
// their checks; no program is published until from_checked_routes finishes it.
[[nodiscard]] ProgramData check_route_definitions(const Graph& graph, std::span<const NodeId> sources,
    std::span<const NodeId> destinations, std::span<const RouteDefinition> routes, const RouteOptions& options);

// A controller's path (CentralizedTE, ShortestPaths Optimized): its source's roster position, its edges in walk
// order, its amount.
struct OptimizedFlowPath {
  std::int32_t roster_index {0};
  std::vector<EdgeId> edges {};
  Amount weight {0};
};

// A flow cut into paths, source by source in ascending NodeId: while a source has supply
// left (its net outflow), take the smallest EdgeId sequence with positive remaining flow, stopping
// at the first destination; assign min(supply left, smallest remaining edge flow on that walk),
// subtract it, record the walk as one path (`slot`: the source's position in `roster`). `edges`
// (ascending) are the flow's edges and `rem` is aligned with them, or `edges` is empty and `rem`
// has one entry per EdgeId. Per source also its share (the edges and amounts its paths use,
// ascending) and its supply. A unit flow cut this way gives the unit walks of T2's program (each
// path repeated `amount` times). std::logic_error when the flow is not conserved, has a positive
// cycle or is not covered by the walks.
struct Decomposition {
  std::vector<OptimizedFlowPath> paths {};
  std::vector<std::vector<std::pair<EdgeId, Amount>>> shares {};
  std::vector<Amount> supplies {};
};
[[nodiscard]] Decomposition decompose(const Graph& g, std::span<const EdgeId> edges, std::vector<Amount> rem,
                                      std::span<const NodeId> roster, std::span<const NodeId> dests);

// A controller program's replay basis: the topology under the producer's exclusions (its excluded
// nodes and edges, ascending), per roster source its route, or one edge set for every source
// (both shared, not copied).
struct ReplayBasis {
  bool uses_topology_basis {false};
  std::span<const NodeId> excluded_nodes {};
  std::span<const EdgeId> excluded_edges {};
  std::vector<std::shared_ptr<const PredecessorDAG>> source_supports {};
  std::shared_ptr<const std::vector<EdgeId>> common_support {};
};

// The program of a controller: one tunnel per path in decomposition order (index = position),
// weight = amount; captured weights `supplies` (per roster source, optimized) or equal sources
// (empty).
[[nodiscard]] std::shared_ptr<const ForwardingProgram> optimized_flow_program(
    const Graph& graph, std::span<const NodeId> sources, std::span<const NodeId> destinations,
    std::span<const NodeId> reachable_sources, std::vector<OptimizedFlowPath> paths, std::span<const Amount> supplies,
    ReplayBasis basis);

// The program of a pooled split demand: one split route per roster source in the one pooled
// context `plan` (rooted at the roster); `entries` one per roster source in roster order (its
// record route, the source's route before weight pruning and its replay basis, and its route cost
// set; index, slot, weight 1 and context are set here).
[[nodiscard]] std::shared_ptr<const ForwardingProgram> shortest_path_program(
    const Graph& graph, std::span<const NodeId> sources, std::span<const NodeId> destinations,
    std::span<const NodeId> reachable_sources, std::shared_ptr<const PreparedRoute> plan, std::vector<ProgramEntry> entries);

// The one pooled DAG builder, shared by ShortestPaths Equal/Weighted and from_dag: the plan of
// `dag` rooted at the roster (the allocation, Equal or Weighted, and the weight source: explicit
// weights edge-aligned, a snapshot aligned with the DAG's ascending edges), then the program of
// shortest_path_program. A plan does not depend on the admission. ShortestPaths with one roster
// source, multipath and a static weight basis builds its plan in one pass instead
// (prepare_shortest_paths, the same plan as this builder's of that source's route) and passes it to
// shortest_path_program directly.
[[nodiscard]] std::shared_ptr<const ForwardingProgram> pooled_dag_program(
    const Graph& graph, std::span<const NodeId> sources, std::span<const NodeId> destinations,
    std::span<const NodeId> roster, const PredecessorDAG& dag, RouteAllocation allocation,
    const WeightSource& weights, std::vector<ProgramEntry> entries);

// A T1 or T2 tunnel as its computation leaves it: its roster slot and round, T1's request, its
// status and (up) its path in walk order with its cost, T2's planned bandwidth (0 for T1).
struct TunnelSetupResult {
  std::int32_t roster_index {0};
  std::int32_t local_index {0};
  Amount requested_bandwidth {0};
  TunnelStatus status {TunnelStatus::Up};
  std::vector<EdgeId> edges {};
  Cost path_cost {0};
  Amount planned {0};
};

// The program of T1 or T2: the policy's roster (every head, an entry or not: T2's dead heads have
// none), its tunnels in round order (`tunnels`: round-major, heads in roster order, so index =
// round * H + slot, a dead head's positions left as holes), weight 1 up and 0 down, reserved =
// requested up and 0 down, planned as given; equal sources; the replay basis the topology under
// the policy's exclusions (ascending).
[[nodiscard]] std::shared_ptr<const ForwardingProgram> distributed_te_program(
    const Graph& graph, std::span<const NodeId> sources, std::span<const NodeId> destinations,
    std::span<const NodeId> reachable_sources, std::vector<TunnelSetupResult> tunnels, std::span<const NodeId> excluded_nodes,
    std::span<const EdgeId> excluded_edges);

// The empty program of an empty roster (no entry).
[[nodiscard]] std::shared_ptr<const ForwardingProgram> empty_program(std::uint64_t graph_uid,
                                                                    std::span<const NodeId> sources,
                                                                    std::span<const NodeId> destinations,
                                                                    bool captured);

// The change flags (routes_changed, weights_changed, reservations_changed) between an old program (null: none) and a new one: routes (the ordered
// heads, the entries' kinds and routes, split tables, duplicates), weights (the injection, every
// entry and split-table weight, an absent entry counting as weight 0) and reservations (every
// tunnel's requested and reserved amounts by index, an absent tunnel counting as (0, 0)). Two
// programs without an entry are equal. With `reservations` false (a call that books no holds) the
// reservation comparison is skipped and reported as no change.
struct ProgramChange {
  bool routes {false};
  bool weights {false};
  bool reservations {false};
};
[[nodiscard]] ProgramChange compare_programs(const ForwardingProgram* before, const ForwardingProgram* after,
                                             bool reservations = true) noexcept;

// Validation of the options of from_routes and from_dag (`function`) before any route is read
// (in order: options on a list of paths only, the allocation, the weight basis's combination), shared with the binding, which reads the entries of edge_weights and residual itself
// after these checks.
void check_route_options(const char* function, bool any_dag, std::optional<RouteAllocation> allocation,
                         std::optional<WeightBasis> weight_basis, bool edge_weights, bool residual);

} // namespace netgraph::core::detail
