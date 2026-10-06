/* The one integer minimum-cost flow core (internal; APPROACH 2.2): the solve on a graph's
   eligible edges behind the recorder's private entry and Global routing, the solve on an edge
   subset behind the MaxFlow kernel, and the deletion re-solve behind the exact sensitivity
   model. */
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

// Potential of a vertex that no Dijkstra run of the call settled. Potentials are otherwise
// non-negative path costs, so the sentinel never collides.
inline constexpr Cost kNeverReached = -1;
// Upper bound of the internal value cap; also the capacity of the virtual source and sink arcs,
// strictly above every attainable value (below 2^62), so a virtual arc never saturates.
inline constexpr Amount kMaxFlowCap = INT64_MAX - 1;

// Operation counters of one call: diagnostics for the performance gate's drift check, identical
// across identical calls, never an oracle. Virtual vertices and arcs count like real ones.
struct FlowCoreCounters {
  std::int64_t dijkstra_runs {0};  // every run, the final failed one included
  std::int64_t dijkstra_pops {0};  // heap pops, stale ones included
  std::int64_t settled {0};        // vertices finalised, summed over runs
  std::int64_t arc_scans {0};      // arcs examined by Dijkstra, Dinic BFS and DFS advance
  std::int64_t dinic_phases {0};   // BFS level builds, the last unsuccessful one of a tier included
  std::int64_t augmentations {0};  // DFS paths pushed
};

// The payload of a solve on a graph (the recorder's _place_global_units, Global routing).
struct FlowCoreResult {
  Amount placed {0};
  // The marginal cost curve: (cost, units) with costs strictly increasing and units positive.
  std::vector<std::pair<Cost, Amount>> cost_distribution {};
  std::vector<Amount> edge_flows {};  // per edge, id order
  // Edges with positive view from the full residual closure of the real sources to its
  // complement, ascending id; nullopt while the destinations stay connected (a cap stopped the
  // solve) and on an input cap of 0.
  std::optional<std::vector<EdgeId>> min_cut {};
  std::vector<Cost> potentials {};  // per node; kNeverReached when never settled
  FlowCoreCounters counters {};
};

// A solve's network in local coordinates: real vertices 0..num_nodes-1, s* = num_nodes,
// t* = num_nodes + 1. The input is the view's edges in ascending EdgeId (tail, head, cost and
// view capacity per edge) and the endpoint sets in ascending NodeId; the arcs are slots in
// adjacency order (the slots of x are [off[x], off[x + 1])), each with its head, its cost and
// the slot of its partner, in the arc order of APPROACH 2.2.
struct Network {
  std::int32_t num_nodes {0};
  std::vector<EdgeId> edges {};
  std::vector<std::int32_t> tail {};
  std::vector<std::int32_t> head {};
  std::vector<Cost> cost {};
  std::vector<Amount> cap {};
  std::vector<std::int32_t> sources {};
  std::vector<std::int32_t> destinations {};
  std::vector<std::int32_t> off {};
  std::vector<std::int32_t> arc_head {};
  std::vector<std::int32_t> rev {};
  std::vector<Cost> arc_cost {};
  std::vector<std::int32_t> edge_slot {};  // per edge: its forward slot (at the tail)
  std::vector<std::int32_t> src_rev {};    // per source: the reverse slot of its virtual arc (at the source)
  std::vector<std::int32_t> dst_fwd {};    // per destination: the forward slot of its virtual arc
};

// The final state of a solve, captured for the deletion re-solve (APPROACH 2.2): bound to the
// solve's graph (whose owner it holds), effective view and endpoint sets (copied into its
// network, shared immutable between copies); the residual capacities, the attained value, the
// working potentials of all N + 2 vertices, the maximum certificate and the solve's counters.
// A copyable value: copies share the network only. Every field belongs to the one solve, so an
// assignment is all or nothing (APPROACH 13): the copy is built aside, then exchanged by moves,
// which never throw; a failed copy assignment leaves the target as it was.
class DeletionState {
public:
  DeletionState() = default;
  DeletionState(const DeletionState&) = default;
  DeletionState(DeletionState&&) = default;
  DeletionState& operator=(const DeletionState& other);
  DeletionState& operator=(DeletionState&&) = default;
  ~DeletionState() = default;

  [[nodiscard]] const Graph& graph() const noexcept { return *graph_; }
  [[nodiscard]] const Network& network() const noexcept { return *net_; }
  [[nodiscard]] std::span<const Amount> residual() const noexcept { return res_; }  // per slot
  [[nodiscard]] Amount value() const noexcept { return value_; }
  // The working potentials per local vertex: the published potential plus the tier maxima of
  // the later tiers; numeric for every vertex (never the published sentinel).
  [[nodiscard]] std::span<const I128> potentials() const noexcept { return pi_; }
  // The maximum certificate: the destinations are outside the full residual closure of the
  // real sources.
  [[nodiscard]] bool maximum() const noexcept { return maximum_; }
  [[nodiscard]] const FlowCoreCounters& counters() const noexcept { return counters_; }
  [[nodiscard]] std::vector<Amount> edge_flows() const;  // per EdgeId of the graph

private:
  friend class FlowSolver;
  friend class DeletionResolver;
  std::shared_ptr<const Graph> graph_ {};
  std::shared_ptr<const Network> net_ {};
  std::vector<Amount> res_ {};
  std::vector<I128> pi_ {};
  Amount value_ {0};
  bool maximum_ {false};
  FlowCoreCounters counters_ {};
};

// The exchange of the copy assignment, and every move, must not throw.
static_assert(std::is_nothrow_move_constructible_v<DeletionState> && std::is_nothrow_move_assignable_v<DeletionState>);

// The internal capture request of a solve: place_global_units stores its final state in
// `state`, bound to `graph`, the shared owner of the solve's graph. The ordinary payload does
// not change. A request, not a value: movable, never copied.
struct DeletionCapture {
  DeletionCapture() = default;
  DeletionCapture(const DeletionCapture&) = delete;
  DeletionCapture(DeletionCapture&&) = default;
  DeletionCapture& operator=(const DeletionCapture&) = delete;
  DeletionCapture& operator=(DeletionCapture&&) = default;
  ~DeletionCapture() = default;

  std::shared_ptr<const Graph> graph {};
  std::optional<DeletionState> state {};
#if NETGRAPH_CORE_PRIVATE_HOOKS
  // The tests' view of the solve: every Dijkstra's distances per local vertex (INT64_MAX when
  // not reached), before its potential update.
  std::function<void(std::span<const Cost>)> on_search {};
#endif
};
static_assert(std::is_nothrow_move_constructible_v<DeletionCapture> &&
              std::is_nothrow_move_assignable_v<DeletionCapture>);

// The result of one deletion re-solve: the repaired value and its certificates (APPROACH 2.2).
struct DeletionResult {
  Amount value {0};
  // The working potentials of all real and virtual vertices (local coordinates, s* = N and
  // t* = N + 1), valid on every positive residual arc.
  std::vector<I128> potentials {};
  // The full residual closure of the real sources to its complement on the deleted view, edges
  // of positive view capacity, ascending EdgeId: value == sum(view capacity[cut]).
  std::vector<EdgeId> cut {};
  std::vector<Amount> edge_flows {};  // the final flow per EdgeId (the deleted edge 0)
  FlowCoreCounters counters {};       // this call's, from zero
};

#if NETGRAPH_CORE_PRIVATE_HOOKS
// The tests' observation of a re-solve: the state right after the deletion (repair 0), after
// every potential update and every push of repairs 1 to 3, and at completion (repair 4).
struct DeletionEvent {
  int repair {0};
  bool push {false};
  Amount moved {0};  // by the current repair so far
  std::span<const Amount> residual {};
  std::span<const I128> potentials {};
  const FlowCoreCounters* counters {nullptr};
};
struct DeletionProbe {
  std::function<void(const DeletionEvent&)> observe {};
  // The control of the search reuse: repair 2 reruns repair 1's failed search.
  bool rerun_failed_search {false};
};
#endif

class FlowSolver;

// The deletion re-solve on one worker's copy of a captured state: each call resets the
// residuals and potentials from that copy, deletes a real arc e = (u, v) with flow f_e > 0 and
// repairs the state by three fixed-origin repairs (u -> v for f_e; u -> s*, then t* -> v, for
// the remainder). The working storage is the resolver's own, allocated at construction.
class DeletionResolver {
public:
  explicit DeletionResolver(const DeletionState& state);
  ~DeletionResolver();
  DeletionResolver(DeletionResolver&&) noexcept;
  DeletionResolver& operator=(DeletionResolver&&) noexcept;
  DeletionResolver(const DeletionResolver&) = delete;
  DeletionResolver& operator=(const DeletionResolver&) = delete;

  // Throws std::invalid_argument before any change when the state is not certified maximum or
  // e is not an edge of the view with positive flow.
  [[nodiscard]] DeletionResult resolve(EdgeId e);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  [[nodiscard]] DeletionResult resolve(EdgeId e, const DeletionProbe& probe);
#endif

private:
  // The input index of e after the precondition checks above.
  [[nodiscard]] std::size_t index_of(EdgeId e) const;
  DeletionState state_;
  std::unique_ptr<FlowSolver> solver_;
};

// Minimum-cost flow of the maximum value attainable under `cap` from the source set to the
// destination set over the eligible edges, on the capacities `view` (units, id order). Every
// solve starts from zero flow. `cap` nullopt is uncapped; an input cap of 0 returns the zero
// payload without running a search. Throws std::invalid_argument before any work, in this
// order, for a view of the wrong length or with an entry outside [0, capacity_units[e]], an
// empty or out-of-range endpoint set, a vertex in both sets, a cap outside [0, 2^63 - 2], or a
// foreign Eligibility. Excluded endpoints are skipped. With `capture`, the final state is
// stored there as well.
[[nodiscard]] FlowCoreResult place_global_units(const Graph& g, std::span<const Amount> view,
                                                std::span<const NodeId> sources,
                                                std::span<const NodeId> destinations, const Eligibility* elig,
                                                std::optional<Amount> cap, DeletionCapture* capture = nullptr);

// The solve on an edge subset (the MaxFlow kernel's replacement on V_i = W + own_i over a
// flow's split graph): flows aligned with the input edges, the value and the marginal curve.
struct SubsetFlow {
  Amount placed {0};
  std::vector<Amount> flows {};
  std::vector<std::pair<Cost, Amount>> cost_distribution {};
};

// The minimum-cost flow of the maximum value under `cap` (nullopt: uncapped) from `sources` to
// `destinations` over the edges `edges` (ascending EdgeId, distinct) with capacities `caps`
// (aligned with `edges`, each in [0, 2^62)). The same core and arc order as
// place_global_units (ascending EdgeId at every node, the virtual arc first, virtual arcs by
// ascending endpoint NodeId), from zero flow; the vertices are compressed to the ones the edges
// and endpoints touch, which changes no flow, curve or augmenting path. The endpoint sets are
// ascending, distinct and disjoint; the caller validated every input.
[[nodiscard]] SubsetFlow solve_subset(const Graph& g, std::span<const EdgeId> edges,
                                      std::span<const Amount> caps, std::span<const NodeId> sources,
                                      std::span<const NodeId> destinations, std::optional<Amount> cap);

} // namespace netgraph::core::detail
