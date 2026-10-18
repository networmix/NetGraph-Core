/* The inside of a ForwardingProgram (internal): its entries, split contexts, injection and replay
   basis, the builders the models use, and the comparisons of the change flags. The public value
   (program.hpp) shows only the endpoints, the tunnels and the source weights. */
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/program.hpp"

namespace netgraph::core::detail {

// One entry: a tunnel (tunnels()[tunnel]) or a split route (a root of contexts[context]). `slot`
// is the head's position in the roster; `bundle` the entry's route as a DAG (a split route's
// edges before weight pruning; null for a tunnel, whose path the record's bundle is made from when a
// Pinned policy needs it) and `cost` its minimum walk cost.
struct ProgramEntry {
  std::int32_t index {0};
  NodeId head {0};
  std::int32_t slot {0};
  std::int32_t local_index {0};
  Amount weight {0};
  std::int32_t tunnel {-1};
  std::int32_t context {-1};
  std::shared_ptr<const PredDAG> bundle {};
  Cost cost {kUnreachable};
};

// A split context: one plan whose roots are the heads it serves (a CostOnly program: the whole
// roster, pooled; a supplied split route: its head alone), and per root the entry it carries.
// `route` (ascending) is the context's edges before weight pruning, the replay basis of its
// members.
struct SplitContext {
  std::shared_ptr<const KernelPlan> plan {};  // its route: the context's edges before weight pruning
  // Per plan root: its entry's position in entries (every root has one), in the program's block
  // (lay_out).
  std::span<const std::int32_t> root_entry {};
};

struct ProgramData {
  ProgramData() = default;
  // The lists below are views of the data's own block (lay_out): a copy would read the
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
  // The injection: the roster (the heads with an entry, ascending), its captured weights (empty:
  // equal sources), every entry in ascending index, per roster slot its entries (ascending index).
  std::span<const NodeId> roster {};
  std::vector<Amount> head_weight {};
  std::vector<ProgramEntry> entries {};
  // Per roster slot its entries (entries_of): slot s's are slot_entries[slot_first[s]..slot_first[s+1]),
  // ascending index; both empty for a pooled program.
  std::span<const std::int32_t> slot_first {};
  std::span<const std::int32_t> slot_entries {};
  std::span<const std::int32_t> tunnel_entry {};  // per tunnel: its entry
  std::vector<SplitContext> contexts {};
  // The replay basis (8.4): the topology under the producer's exclusions (topology: unrestricted
  // Global), or per roster slot the edges of its original routes (basis_parts).
  bool topology {false};
  std::span<const NodeId> excluded_nodes {};  // ascending
  std::span<const EdgeId> excluded_edges {};  // ascending (an edge whose endpoint is excluded included)
  std::shared_ptr<const std::vector<EdgeId>> common_basis {};  // every slot's (a supported Global)
  std::vector<std::shared_ptr<const PredDAG>> slot_routes {};    // per slot its route (CostOnly MaxFlow)
  // A slot's replay basis as edge sets, appended to `out` (each ascending, except a tunnel's path,
  // in walk order: `second` false): the shared route or edge set above, else (supplied and CostOnly
  // split programs) the slot's own entries' routes: a tunnel's edges, a split route's bundle (its
  // edges before weight pruning). Nothing is copied.
  void basis_parts(std::size_t slot, std::vector<std::pair<std::span<const EdgeId>, bool>>& out) const;
  // One pooled context carries every entry: a single context whose roots are the roster, one entry
  // of positive weight per head, equal sources (a CostOnly program; a supplied program of one
  // split route of positive weight). Its forwarding is the plan's own kernels.
  bool pooled {false};
  // A slot's entries, ascending index: a pooled program's slot carries the one entry its context
  // root carries (no per-slot lists are kept).
  [[nodiscard]] std::span<const std::int32_t> entries_of(std::size_t slot) const noexcept {
    if (pooled) return {&contexts.front().root_entry[slot], 1};
    const auto first = static_cast<std::size_t>(slot_first[slot]);
    return {slot_entries.data() + first, static_cast<std::size_t>(slot_first[slot + 1]) - first};
  }
  // The last step of every producer, after its entries, tunnels, contexts and `pooled` are set:
  // the lists above copied into one block (one allocation; a failure leaves the data without
  // lists, to be discarded), the roster (given, or the entries' heads, ascending and distinct),
  // each entry's slot and local index (its position among its slot's entries), each tunnel's
  // local index, the per-slot lists (a stable counting sort by slot), the per-tunnel entries and
  // each context's root entries (the entry whose head is the root).
  void lay_out(std::span<const NodeId> sources, std::span<const NodeId> destinations,
               std::optional<std::span<const NodeId>> roster, std::span<const NodeId> excluded_nodes = {},
               std::span<const EdgeId> excluded_edges = {});
  // The block the lists view (the ownership tests).
  [[nodiscard]] std::span<const std::int32_t> block() const noexcept { return store_; }

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
};

// A controller's path (Global, CostOnly MaxFlow): its source's roster position, its edges in walk
// order, its amount.
struct ControllerPath {
  std::int32_t slot {0};
  std::vector<EdgeId> edges {};
  Amount amount {0};
};

// A controller program's replay basis: the topology under the producer's exclusions (its excluded
// nodes and edges, ascending), per roster source its route, or one edge set for every source
// (both shared, not copied).
struct ControllerBasis {
  bool topology {false};
  std::span<const NodeId> excluded_nodes {};
  std::span<const EdgeId> excluded_edges {};
  std::vector<std::shared_ptr<const PredDAG>> per_source {};
  std::shared_ptr<const std::vector<EdgeId>> common {};
};

// The program of a controller: one tunnel per path in decomposition order (index = position),
// weight = amount; captured weights `supplies` (per roster source, optimized) or equal sources
// (empty).
[[nodiscard]] std::shared_ptr<const ForwardingProgram> controller_program(
    const Graph& g, std::span<const NodeId> sources, std::span<const NodeId> destinations,
    std::span<const NodeId> roster, std::vector<ControllerPath> paths, std::span<const Amount> supplies,
    ControllerBasis basis);

// The program of a CostOnly split demand: one split route per roster source in the one pooled
// context `plan` (rooted at the roster); `entries` one per roster source in roster order (its
// record bundle, the source's route before weight pruning and its replay basis, and its route cost
// set; index, slot, weight 1 and context are set here).
[[nodiscard]] std::shared_ptr<const ForwardingProgram> cost_only_program(
    const Graph& g, std::span<const NodeId> sources, std::span<const NodeId> destinations,
    std::span<const NodeId> roster, std::shared_ptr<const KernelPlan> plan, std::vector<ProgramEntry> entries);

// The empty program of an empty roster (no entry).
[[nodiscard]] std::shared_ptr<const ForwardingProgram> empty_program(std::uint64_t graph_uid,
                                                                    std::span<const NodeId> sources,
                                                                    std::span<const NodeId> destinations,
                                                                    bool captured);

// The change flags of 10.4 between an old program (null: none) and a new one: routes (the ordered
// heads, the entries' kinds and routes, split tables, duplicates) and weights (the injection, every
// entry and split-table weight). Two programs without an entry are equal.
struct ProgramChange {
  bool routes {false};
  bool weights {false};
};
[[nodiscard]] ProgramChange compare_programs(const ForwardingProgram* before, const ForwardingProgram* after) noexcept;

// Validation of from_routes' options before any route is read (steps 3 to 5 of the order: options
// on a list of paths only, the placement, the weights' combination), shared with the binding,
// which reads the entries of edge_weights and residual itself after these checks.
void check_route_options(bool any_dag, std::optional<FlowPlacement> placement, std::optional<WeightBasis> weights,
                         bool edge_weights, bool residual);

} // namespace netgraph::core::detail
