/* The forwarding of a program (internal): one engine for every model that forwards a program
   (ShortestPaths's own pooled program, a controller's paths through forward(), and Pinned).

   A view is the program as one forwarding sees it: its replay under the forwarding policy's mask
   (docs/semantics.md, "Replay under a failure": failed tunnels, split routes without their excluded
   members, the roster rebuilt from the program's basis) and the black holes laid over it ("Stale
   forwarding": they change no split; traffic reaching one is dropped there, a black-holed head
   loses its own share at the source).

   Injection: a total enters as round robin over the view's heads (equal sources) or by highest
   averages over their captured weights; each head splits its share by highest averages over its
   surviving entries' weights. A head without a surviving entry of positive weight, or a
   black-holed head, keeps its share and loses it at the source; an entry whose route has no
   member left (a Weighted route of weights 0) loses what it receives at the source too.

   Lossless: the largest total in an interval whose no-drop loads fit the free capacity with
   no share lost at a source or on a black hole; the loads are monotone in the total (round robin
   and highest averages are house-monotone), so integer bisection is exact.
   Lossy ("Losses inside one program"): the one evaluator of the loss rule. A link that receives
   more than its free capacity V shares T = min(V, arrivals) by highest averages over the sources' arrivals, then within a
   source over its entries'; A+ forwards the no-drop loads U0 first, makes every link with
   U0 <= V a pass-through and evaluates the other links (barriers) and the split routes' routers
   once each in dependency order; barriers that wait on each other in a cycle are iterated from
   zero (serial normalization first, then rounds, 64 b rounds for b barriers, the period-two
   and round-limit errors, the independent check of the whole result before anything is
   returned). A pooled program (one split context: a ShortestPaths program, one supplied split route)
   is forwarded by its plan's kernels, which are this rule for one acyclic context.
*/
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/detail/route_kernels.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/forwarding_program.hpp"

namespace netgraph::core::detail {

struct ForwardingView {
  std::shared_ptr<const ForwardingProgram> program {};
  const ProgramData* data {nullptr};
  // A borrowed view is the program's own state (no mask, no black holes, no record route to
  // make): every head of the roster injects, every entry is up with the program's route, cost and
  // plan, and the vectors below stay empty (the accessors read the program). Otherwise they hold
  // the replay:
  bool borrowed {false};
  // The heads that inject (the replay roster, ascending NodeId) as roster slots; their captured
  // weights (empty: equal sources); black-holed heads.
  std::vector<std::int32_t> source_indices_ {};
  std::vector<Amount> source_weights_ {};
  std::vector<std::uint8_t> blackholed_sources_ {};
  // Per program entry: Failed (excluded by the mask, or a head outside the replay roster),
  // Blackholed (a tunnel whose path reaches a black hole), a down tunnel's setup status
  // (InsufficientBandwidth, CostCeilingExceeded) or Up; its route under the mask (the record's route; an empty DAG
  // when it failed or is down) and that route's cost.
  std::vector<TunnelStatus> status_ {};
  std::vector<std::shared_ptr<const PredecessorDAG>> routes_ {};
  std::vector<Cost> min_route_costs_ {};
  // Per view head: its surviving entries (positions in data->entries, ascending index).
  std::vector<std::int32_t> source_entry_offsets_ {};
  std::vector<std::int32_t> source_entry_positions_ {};
  // The split contexts as forwarded (restricted to the mask; the program's plan when nothing is
  // removed) and per entry of a context its root's index in that plan's roots (-1 otherwise).
  std::vector<std::shared_ptr<const PreparedRoute>> prepared_routes_ {};
  std::vector<std::int32_t> context_source_indices_ {};
  // Black holes: the nodes and edges this Eligibility excludes (null: none).
  const Eligibility* black_holes {nullptr};
  // The pooled form: one context, one root per view head (equal sources, every entry weight
  // positive). Its plan is plan(0), its roots the view's heads in order.
  bool pooled {false};

  [[nodiscard]] std::size_t num_sources() const noexcept { return borrowed ? data->reachable_sources.size() : source_indices_.size(); }
  [[nodiscard]] std::int32_t source_index(std::size_t h) const noexcept {
    return borrowed ? static_cast<std::int32_t>(h) : source_indices_[h];
  }
  [[nodiscard]] NodeId source(std::size_t h) const noexcept { return data->reachable_sources[static_cast<std::size_t>(source_index(h))]; }
  [[nodiscard]] std::span<const Amount> source_weights() const noexcept {
    return borrowed ? std::span<const Amount>(data->source_weight_values) : std::span<const Amount>(source_weights_);
  }
  [[nodiscard]] bool is_blackholed_source(std::size_t h) const noexcept { return !blackholed_sources_.empty() && blackholed_sources_[h] != 0; }
  // Entry e's status: a borrowed view's are the program's own (Up, or a down tunnel's setup status).
  [[nodiscard]] TunnelStatus status(std::size_t e) const noexcept {
    if (!borrowed) return status_[e];
    const auto t = data->entries[e].tunnel_position;
    return t < 0 ? TunnelStatus::Up : data->tunnels[static_cast<std::size_t>(t)].status;
  }
  [[nodiscard]] const std::shared_ptr<const PredecessorDAG>& route(std::size_t e) const noexcept {
    return borrowed ? data->entries[e].route : routes_[e];
  }
  [[nodiscard]] Cost min_route_cost(std::size_t e) const noexcept { return borrowed ? data->entries[e].min_route_cost : min_route_costs_[e]; }
  [[nodiscard]] const std::shared_ptr<const PreparedRoute>& prepared_route(std::size_t c) const noexcept {
    return borrowed ? data->split_contexts[c].prepared_route : prepared_routes_[c];
  }
  // Entry e's root in its context's plan (-1 for a tunnel or a root the mask removed).
  [[nodiscard]] std::int32_t context_source_index(std::size_t e) const noexcept {
    if (!borrowed) return context_source_indices_[e];
    const auto c = data->entries[e].split_context_position;
    if (c < 0) return -1;
    const auto& roots = data->split_contexts[static_cast<std::size_t>(c)].source_entry_positions;
    if (e < roots.size() && static_cast<std::size_t>(roots[e]) == e) return static_cast<std::int32_t>(e);  // pooled
    for (std::size_t r = 0; r < roots.size(); ++r) {
      if (static_cast<std::size_t>(roots[r]) == e) return static_cast<std::int32_t>(r);
    }
    return -1;
  }
  // View head h's surviving entries, ascending index.
  [[nodiscard]] std::span<const std::int32_t> entries_for_source(std::size_t h) const noexcept {
    if (borrowed) return data->entries_for_source(h);
    const auto b = static_cast<std::size_t>(source_entry_offsets_[h]);
    return {source_entry_positions_.data() + b, static_cast<std::size_t>(source_entry_offsets_[h + 1]) - b};
  }
};

// The view of `program` under the forwarding policy's mask `mask` (null: the program's own
// state, nothing replayed) with `dests` the allowed destinations (ascending), and the black holes
// `dark` (null or an all-allowed value: none; it must outlive the view). With `records`, every
// entry's route under the mask is made (a Pinned policy's record routes); without, a tunnel's
// route stays null.
[[nodiscard]] ForwardingView make_forwarding_view(const Graph& graph, std::shared_ptr<const ForwardingProgram> program,
                                       const Eligibility* mask, std::span<const NodeId> destinations,
                                       const Eligibility* dark, bool records = false);
// The same view with other black holes (the program, its replay and its structure shared).
[[nodiscard]] ForwardingView with_black_holes(const ForwardingView& view, const Eligibility* dark);

// One entry's part of a forwarding: its share (assigned), what it delivered, what it lost at the
// source (a route without a member), what it carried and dropped per edge (ascending EdgeId,
// positive entries).
struct EntryFlow {
  Amount assigned {0};
  Amount delivered {0};
  Amount source_drops {0};
  std::vector<std::pair<EdgeId, Amount>> edge_flows {};
  std::vector<std::pair<EdgeId, Amount>> edge_drops {};
};

// A forwarding's outcome: the heads' injections (offered; with captured weights all 0 a total
// can exceed it: the rest is lost at the sources, outside every head), lost at the sources,
// delivered; per view head its injection and (Lossy only: empty for Lossless, which loses
// nothing) the part it lost itself (black-holed, or no surviving entry of positive weight); per
// program entry its EntryFlow (zero for entries that receive nothing).
struct ProgramFlow {
  Amount admitted {0};
  Amount offered {0};
  Amount source_drops {0};
  Amount delivered {0};
  std::vector<Amount> offered_by_source {};
  std::vector<Amount> source_drops_by_source {};
  std::vector<EntryFlow> entries {};
};

// A forwarding's network numbers its positions in 32 bits. Its counts, in the order they are
// checked: the forwarded entries, the occurrences (every tunnel hop and split-route member), the
// routers of the split routes, the lossy rule's events (the links that can drop plus the routers)
// and their dependency arcs. Each is computed exactly in 64 bits: at most kForwardingLimit, else
// std::invalid_argument `forwarding <counter> limit exceeded: <n> > 2147483647`. The entry,
// occurrence and router counts are checked after validation and before the network is allocated;
// the event and dependency counts, known only after the no-drop loads, before the event
// structures are allocated. Either refusal changes nothing.
inline constexpr std::int64_t kForwardingLimit = std::numeric_limits<std::int32_t>::max();
enum class ForwardingCounter { Entry, Occurrence, Router, Event, Dependency };
void check_forwarding_count(ForwardingCounter counter, std::int64_t n);

// The free capacity of each of `edges` in the forwarding, into `out` (aligned; before black
// holes, which the engine applies): one call per evaluation, not per edge.
using AvailableCapacityFn = std::function<void(std::span<const EdgeId> edges, std::span<Amount> out)>;

// The heads' injections of a total x: round robin over the heads, or highest averages over the
// captured weights.
void source_shares(const ForwardingView& view, Amount x, std::span<Amount> out);
// The heads' increments when a cumulative 128-bit offer X grows by v (a Lossy place): RR(X + v) -
// RR(X) or H(X + v) - H(X), each below 2^62.
void source_offer_increments(const ForwardingView& view, const Counter128& X, Amount v, std::span<Amount> out);

// Lossless: the largest total in [lo, hi] (hi nullopt: unbounded) whose no-drop forwarding fits
// `free`, with its loads; `lo` must fit.
[[nodiscard]] ProgramFlow lossless(const ForwardingView& view, const AvailableCapacityFn& free, Amount lo, std::optional<Amount> hi);
// Lossy: the loss rule over per-head injections (aligned with the view's heads) on `free`.
[[nodiscard]] ProgramFlow lossy(const ForwardingView& view, const AvailableCapacityFn& free, std::span<const Amount> injections);

#if NETGRAPH_CORE_PRIVATE_HOOKS
// Test controls (hooks builds only): how many cyclic groups the last lossy() on this thread
// iterated and its rounds; a corruption of the iterated result before its check (the check's
// planted fault: RuntimeError expected).
struct LossyCounters {
  std::int64_t groups {0};
  std::int64_t rounds {0};
  bool corrupt {false};
};
LossyCounters& thread_lossy_counters() noexcept;
// Lower limits for this thread's forwardings, per ForwardingCounter (0: kForwardingLimit), so a
// test reaches each refusal with a small program.
std::array<std::int64_t, 5>& thread_forwarding_limits() noexcept;
#endif

} // namespace netgraph::core::detail
