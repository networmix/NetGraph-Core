/* The forwarding of a program (internal): one engine for every model that forwards a program
   (CostOnly's own pooled program, a controller's paths through forward(), and Pinned).

   A view is the program as one forwarding sees it: its replay under the forwarding policy's mask
   (8.4: failed tunnels, split routes without their excluded members, the roster rebuilt from the
   program's basis) and the black holes laid over it (8.6: they change no split; traffic reaching
   one is dropped there, a black-holed head loses its own share at the source).

   Injection: a total enters as round robin over the view's heads (equal sources) or by highest
   averages over their captured weights; each head splits its share by highest averages over its
   surviving entries' weights. A head without a surviving entry of positive weight, or a
   black-holed head, keeps its share and loses it at the source; an entry whose route has no
   member left (a Weighted route of weights 0) loses what it receives at the source too.

   Lossless (8.3): the largest total in an interval whose no-drop loads fit the free capacity with
   no share lost at a source or on a black hole; the loads are monotone in the total (round robin
   and highest averages are house-monotone), so integer bisection is exact.
   Lossy (9): the one evaluator of the loss rule. A link that receives more than its free capacity
   V shares T = min(V, arrivals) by highest averages over the sources' arrivals, then within a
   source over its entries'; A+ (9.2) forwards the no-drop loads U0 first, makes every link with
   U0 <= V a pass-through and evaluates the other links (barriers) and the split routes' routers
   once each in dependency order; barriers that wait on each other in a cycle are iterated from
   zero (9.4: serial normalization first, then rounds, 64 b rounds for b barriers, the period-two
   and round-limit errors, the independent check of the whole result before anything is
   returned). A pooled program (one split context: a CostOnly program, one supplied split route)
   is forwarded by its plan's kernels, which are this rule for one acyclic context.
*/
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/program.hpp"

namespace netgraph::core::detail {

struct ForwardingView {
  std::shared_ptr<const ForwardingProgram> program {};
  const ProgramData* data {nullptr};
  // A borrowed view is the program's own state (no mask, no black holes, no record bundle to
  // make): every head of the roster injects, every entry is up with the program's route, cost and
  // plan, and the vectors below stay empty (the accessors read the program). Otherwise they hold
  // the replay:
  bool borrowed {false};
  // The heads that inject (the replay roster, ascending NodeId) as roster slots; their captured
  // weights (empty: equal sources); black-holed heads.
  std::vector<std::int32_t> slots_ {};
  std::vector<Amount> weights_ {};
  std::vector<std::uint8_t> dark_head_ {};
  // Per program entry: Failed (excluded by the mask, or a head outside the replay roster),
  // Blackholed (a tunnel whose path reaches a black hole), a down tunnel's setup status
  // (NoBandwidth, CostCeiling) or Up; its route under the mask (the record's bundle; an empty DAG
  // when it failed or is down) and that route's cost.
  std::vector<TunnelStatus> status_ {};
  std::vector<std::shared_ptr<const PredDAG>> bundles_ {};
  std::vector<Cost> costs_ {};
  // Per view head: its surviving entries (positions in data->entries, ascending index).
  std::vector<std::int32_t> head_begin_ {};
  std::vector<std::int32_t> head_entries_ {};
  // The split contexts as forwarded (restricted to the mask; the program's plan when nothing is
  // removed) and per entry of a context its root's index in that plan's roots (-1 otherwise).
  std::vector<std::shared_ptr<const KernelPlan>> plans_ {};
  std::vector<std::int32_t> root_of_ {};
  // Black holes: the nodes and edges this Eligibility excludes (null: none).
  const Eligibility* dark {nullptr};
  // The pooled form: one context, one root per view head (equal sources, every entry weight
  // positive). Its plan is plan(0), its roots the view's heads in order.
  bool pooled {false};

  [[nodiscard]] std::size_t num_heads() const noexcept { return borrowed ? data->roster.size() : slots_.size(); }
  [[nodiscard]] std::int32_t slot(std::size_t h) const noexcept {
    return borrowed ? static_cast<std::int32_t>(h) : slots_[h];
  }
  [[nodiscard]] NodeId head(std::size_t h) const noexcept { return data->roster[static_cast<std::size_t>(slot(h))]; }
  [[nodiscard]] std::span<const Amount> weights() const noexcept {
    return borrowed ? std::span<const Amount>(data->head_weight) : std::span<const Amount>(weights_);
  }
  [[nodiscard]] bool dark_head(std::size_t h) const noexcept { return !dark_head_.empty() && dark_head_[h] != 0; }
  // Entry e's status: a borrowed view's are the program's own (Up, or a down tunnel's setup status).
  [[nodiscard]] TunnelStatus status(std::size_t e) const noexcept {
    if (!borrowed) return status_[e];
    const auto t = data->entries[e].tunnel;
    return t < 0 ? TunnelStatus::Up : data->tunnels[static_cast<std::size_t>(t)].status;
  }
  [[nodiscard]] const std::shared_ptr<const PredDAG>& bundle(std::size_t e) const noexcept {
    return borrowed ? data->entries[e].bundle : bundles_[e];
  }
  [[nodiscard]] Cost cost(std::size_t e) const noexcept { return borrowed ? data->entries[e].cost : costs_[e]; }
  [[nodiscard]] const std::shared_ptr<const KernelPlan>& plan(std::size_t c) const noexcept {
    return borrowed ? data->contexts[c].plan : plans_[c];
  }
  // Entry e's root in its context's plan (-1 for a tunnel or a root the mask removed).
  [[nodiscard]] std::int32_t root_of(std::size_t e) const noexcept {
    if (!borrowed) return root_of_[e];
    const auto c = data->entries[e].context;
    if (c < 0) return -1;
    const auto& roots = data->contexts[static_cast<std::size_t>(c)].root_entry;
    if (e < roots.size() && static_cast<std::size_t>(roots[e]) == e) return static_cast<std::int32_t>(e);  // pooled
    for (std::size_t r = 0; r < roots.size(); ++r) {
      if (static_cast<std::size_t>(roots[r]) == e) return static_cast<std::int32_t>(r);
    }
    return -1;
  }
  // View head h's surviving entries, ascending index.
  [[nodiscard]] std::span<const std::int32_t> head_entries(std::size_t h) const noexcept {
    if (borrowed) return data->entries_of(h);
    const auto b = static_cast<std::size_t>(head_begin_[h]);
    return {head_entries_.data() + b, static_cast<std::size_t>(head_begin_[h + 1]) - b};
  }
};

// The view of `program` under the forwarding policy's mask `mask` (null: the program's own
// state, nothing replayed) with `dests` the allowed destinations (ascending), and the black holes
// `dark` (null or an all-allowed value: none; it must outlive the view). With `records`, every
// entry's route under the mask is made (a Pinned policy's record bundles); without, a tunnel's
// bundle stays null.
[[nodiscard]] ForwardingView make_view(const Graph& g, std::shared_ptr<const ForwardingProgram> program,
                                       const Eligibility* mask, std::span<const NodeId> dests,
                                       const Eligibility* dark, bool records = false);
// The same view with other black holes (the program, its replay and its structure shared).
[[nodiscard]] ForwardingView with_black_holes(const ForwardingView& view, const Eligibility* dark);

// One entry's part of a forwarding: its share (assigned), what it delivered, what it lost at the
// source (a route without a member), what it carried and dropped per edge (ascending EdgeId,
// positive entries).
struct EntryFlow {
  Amount assigned {0};
  Amount delivered {0};
  Amount lost {0};
  std::vector<std::pair<EdgeId, Amount>> carried {};
  std::vector<std::pair<EdgeId, Amount>> dropped {};
};

// A forwarding's outcome: the heads' injections (offered; with captured weights all 0 a total
// can exceed it: the rest is lost at the sources, outside every head), lost at the sources,
// delivered; per view head its injection and (Lossy only: empty for Lossless, which loses
// nothing) the part it lost itself (black-holed, or no surviving entry of positive weight); per
// program entry its EntryFlow (zero for entries that receive nothing).
struct ProgramFlow {
  Amount total {0};
  Amount offered {0};
  Amount source_loss {0};
  Amount delivered {0};
  std::vector<Amount> head_offered {};
  std::vector<Amount> head_lost {};
  std::vector<EntryFlow> entries {};
};

// The free capacity of each of `edges` in the forwarding, into `out` (aligned; before black
// holes, which the engine applies): one call per evaluation, not per edge.
using FreeFn = std::function<void(std::span<const EdgeId> edges, std::span<Amount> out)>;

// The heads' injections of a total x: round robin over the heads, or highest averages over the
// captured weights.
void head_shares(const ForwardingView& view, Amount x, std::span<Amount> out);
// The heads' increments when a cumulative 128-bit offer X grows by v (a Lossy place): RR(X + v) -
// RR(X) or H(X + v) - H(X), each below 2^62.
void head_increments(const ForwardingView& view, const Counter128& X, Amount v, std::span<Amount> out);

// Lossless: the largest total in [lo, hi] (hi nullopt: unbounded) whose no-drop forwarding fits
// `free`, with its loads; `lo` must fit.
[[nodiscard]] ProgramFlow lossless(const Graph& g, const ForwardingView& view, const FreeFn& free, Amount lo,
                                   std::optional<Amount> hi);
// Lossy: the loss rule over per-head injections (aligned with the view's heads) on `free`.
[[nodiscard]] ProgramFlow lossy(const Graph& g, const ForwardingView& view, const FreeFn& free,
                                std::span<const Amount> injections);

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
#endif

} // namespace netgraph::core::detail
