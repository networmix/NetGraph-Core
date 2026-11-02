/*
  The integer ledger FlowLedger.

  The aggregate allocation and the aggregate holds are kept exactly (integers) and changed only by
  an atomic commit, a detach or a reset (which zeroes them); residual() and unreserved() are
  derived from them (each baseline minus its aggregate). The two projections are never added.
  Each policy epoch is one entry keyed by its first FlowId; its records are stored or derived
  when read. Ids come from a per-ledger counter and are never reused, also across reset().
*/
#include "netgraph/core/flow_ledger.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "netgraph/core/detail/edge_amounts.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/uid.hpp"
#include "netgraph/core/detail/forwarding.hpp"
#include "netgraph/core/detail/route_builders.hpp"

namespace netgraph::core {

FlowLedger::~FlowLedger() = default;

FlowLedger::FlowLedger(std::shared_ptr<const Graph> graph) {
  if (!graph) throw std::invalid_argument("FlowLedger: graph is null");
  graph_ = std::move(graph);
  service_baseline_.assign(graph_->edge_capacity_units().begin(), graph_->edge_capacity_units().end());
  carried_.assign(service_baseline_.size(), 0);
  uid_ = detail::next_uid();
}

FlowLedger::FlowLedger(std::shared_ptr<const Graph> graph, std::span<const Amount> service_baseline) {
  if (!graph) throw std::invalid_argument("FlowLedger: graph is null");
  graph_ = std::move(graph);
  service_baseline_ = checked_baseline(service_baseline, "service_baseline");
  carried_.assign(service_baseline_.size(), 0);
  uid_ = detail::next_uid();
}

FlowLedger::FlowLedger(std::shared_ptr<const Graph> graph, std::span<const Amount> service_baseline,
                       std::span<const Amount> reservation_baseline) {
  if (!graph) throw std::invalid_argument("FlowLedger: graph is null");
  graph_ = std::move(graph);
  service_baseline_ = checked_baseline(service_baseline, "service_baseline");
  reservation_baseline_ = checked_baseline(reservation_baseline, "reservation_baseline");
  carried_.assign(service_baseline_.size(), 0);
  uid_ = detail::next_uid();
}

FlowLedger::FlowLedger(const FlowLedger& other) {
  const auto guard = detail::LedgerAccess::lock(other);
  other.check_usable();
  if (other.num_bound_policies_ != 0) {
    throw std::invalid_argument("FlowLedger: a ledger with bound policies cannot be copied");
  }
  graph_ = other.graph_;
  uid_ = detail::next_uid();
  service_baseline_ = other.service_baseline_;
  carried_ = other.carried_;
  reservation_baseline_ = other.reservation_baseline_;
  reserved_ = other.reserved_;
  for (const auto& [id, entry] : other.entries_) entries_.emplace(id, std::make_unique<detail::LedgerEntry>(*entry));
  next_record_counter_ = other.next_record_counter_;
  next_policy_id_ = other.next_policy_id_;
  policy_bindings_ = other.policy_bindings_;  // orphans only
}

FlowLedger::FlowLedger(FlowLedger&& other) {
  const auto guard = detail::LedgerAccess::lock(other);
  other.check_usable();
  if (other.num_bound_policies_ != 0) {
    throw std::invalid_argument("FlowLedger: a ledger with bound policies cannot be moved");
  }
  graph_ = other.graph_;
  uid_ = other.uid_;
  service_baseline_ = std::move(other.service_baseline_);
  carried_ = std::move(other.carried_);
  reservation_baseline_ = std::move(other.reservation_baseline_);
  reserved_ = std::move(other.reserved_);
  entries_ = std::move(other.entries_);
  next_record_counter_ = other.next_record_counter_;
  next_policy_id_ = other.next_policy_id_;
  policy_bindings_ = std::move(other.policy_bindings_);
  other.usable_ = false;
  other.service_baseline_.clear();
  other.carried_.clear();
  other.reservation_baseline_.clear();
  other.reserved_.clear();
  other.entries_.clear();
  other.policy_bindings_.clear();
}

void FlowLedger::check_usable() const {
  if (!usable_) throw std::logic_error("FlowLedger: a moved-from ledger is empty and unusable");
}

std::vector<Amount> FlowLedger::checked_baseline(std::span<const Amount> amounts, const char* label) const {
  detail::check_edge_amounts(*graph_, amounts, label);
  return std::vector<Amount>(amounts.begin(), amounts.end());
}

std::uint64_t FlowLedger::uid() const {
  const auto guard = detail::LedgerAccess::lock(*this);
  check_usable();
  return uid_;
}

const std::shared_ptr<const Graph>& FlowLedger::graph() const {
  const auto guard = detail::LedgerAccess::lock(*this);
  check_usable();
  return graph_;
}

std::span<const Amount> FlowLedger::edge_capacity_units() const {
  const auto guard = detail::LedgerAccess::lock(*this);
  check_usable();
  return graph_->edge_capacity_units();
}

std::vector<Amount> FlowLedger::residual() const {
  const auto guard = detail::LedgerAccess::lock(*this);
  check_usable();
  std::vector<Amount> out(service_baseline_.size());
  for (std::size_t e = 0; e < out.size(); ++e) out[e] = service_baseline_[e] - carried_[e];
  return out;
}

std::vector<Amount> FlowLedger::carried() const {
  const auto guard = detail::LedgerAccess::lock(*this);
  check_usable();
  return carried_;
}

std::vector<Amount> FlowLedger::reserved() const {
  const auto guard = detail::LedgerAccess::lock(*this);
  check_usable();
  return reserved_.empty() ? std::vector<Amount>(service_baseline_.size(), 0) : reserved_;
}

std::vector<Amount> FlowLedger::unreserved() const {
  const auto guard = detail::LedgerAccess::lock(*this);
  check_usable();
  std::vector<Amount> out(service_baseline_.size());
  for (std::size_t e = 0; e < out.size(); ++e) out[e] = reservation_at(e) - reserved_at(e);
  return out;
}

std::vector<FlowRecord> FlowLedger::records() const {
  const auto guard = detail::LedgerAccess::lock(*this);
  check_usable();
  std::vector<FlowRecord> out;
  std::size_t count = 0;
  for (const auto& [id, entry] : entries_) count += entry->record_inventory().size();
  out.reserve(count);
  for (const auto& [id, entry] : entries_) {
    auto records = entry->records(*graph_, uid_, 0, entry->record_inventory().size());
    out.insert(out.end(), std::make_move_iterator(records.begin()), std::make_move_iterator(records.end()));
  }
  return out;
}

FlowRecord FlowLedger::record(const FlowId& flow_id) const {
  const auto guard = detail::LedgerAccess::lock(*this);
  return detail::LedgerAccess::record(*this, flow_id);
}

FlowRecord detail::LedgerAccess::record(const FlowLedger& ledger, const FlowId& flow_id) {
  ledger.check_usable();
  if (flow_id.first != ledger.uid_) {
    throw std::invalid_argument("FlowId (" + std::to_string(flow_id.first) + ", " + std::to_string(flow_id.second) +
                                ") belongs to another ledger (uid " + std::to_string(ledger.uid_) + ")");
  }
  return materialize(ledger, flow_id.second);
}

std::vector<EdgeId> FlowLedger::saturated_boundary_edges(std::span<const NodeId> sources, const Eligibility* eligibility) const {
  const auto& g = *graph_;
  const auto srcs = detail::node_set(sources, g.num_nodes(), "sources");
  if (eligibility != nullptr && eligibility->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  const auto guard = detail::LedgerAccess::lock(*this);
  std::vector<Amount> residual;
  detail::LedgerAccess::residual_into(*this, residual);
  return detail::saturated_boundary_edges_of(g, residual, srcs, eligibility);
}

namespace detail {

std::vector<EdgeId> saturated_boundary_edges_of(const Graph& graph, std::span<const Amount> residual,
                                  std::span<const NodeId> sources, const Eligibility* eligibility) {
  const auto n = static_cast<std::size_t>(graph.num_nodes());
  const auto m = static_cast<std::size_t>(graph.num_edges());
  const auto edge_sources = graph.edge_sources();
  const auto edge_destinations = graph.edge_destinations();
  const auto capacity = graph.edge_capacity_units();
  const auto row = graph.outgoing_offsets();
  std::vector<std::uint8_t> seen(n, 0);
  std::vector<NodeId> stack;
  for (NodeId s : sources) {
    if (eligibility != nullptr && !eligibility->allows_node(s)) continue;
    if (seen[static_cast<std::size_t>(s)] != 0) continue;
    seen[static_cast<std::size_t>(s)] = 1;
    stack.push_back(s);
  }
  while (!stack.empty()) {
    const auto u = static_cast<std::size_t>(stack.back());
    stack.pop_back();
    for (auto e = row[u]; e < row[u + 1]; ++e) {
      const auto ei = static_cast<std::size_t>(e);
      if (eligibility != nullptr && !eligibility->allows_edge(e)) continue;
      if (residual[ei] < 1) continue;
      const auto v = static_cast<std::size_t>(edge_destinations[ei]);
      if (seen[v] == 0) {
        seen[v] = 1;
        stack.push_back(static_cast<NodeId>(v));
      }
    }
  }
  std::vector<EdgeId> out;
  for (std::size_t e = 0; e < m; ++e) {
    if (seen[static_cast<std::size_t>(edge_sources[e])] == 0 || seen[static_cast<std::size_t>(edge_destinations[e])] != 0) continue;
    if (eligibility != nullptr && !eligibility->allows_edge(static_cast<EdgeId>(e))) continue;
    if (capacity[e] < 1) continue;
    out.push_back(static_cast<EdgeId>(e));
  }
  return out;
}

} // namespace detail

void FlowLedger::reset() {
  reset(graph_->edge_capacity_units());
}

void FlowLedger::reset(std::span<const Amount> service_baseline) {
  auto baseline = checked_baseline(service_baseline, "service_baseline");
  const auto guard = detail::LedgerAccess::lock(*this);
  check_usable();
  check_resettable();
  install(std::move(baseline), {});
}

void FlowLedger::reset(std::span<const Amount> service_baseline, std::span<const Amount> reservation_baseline) {
  auto baseline = checked_baseline(service_baseline, "service_baseline");
  auto reservation = checked_baseline(reservation_baseline, "reservation_baseline");
  const auto guard = detail::LedgerAccess::lock(*this);
  check_usable();
  check_resettable();
  install(std::move(baseline), std::move(reservation));
}

void FlowLedger::check_resettable() const {
  if (num_bound_policies_ != 0) {
    throw std::invalid_argument("FlowLedger.reset is rejected while a policy is bound (" +
                                std::to_string(num_bound_policies_) + " bound)");
  }
}

// The new baselines (an empty reservation: the capacities) and no records; carried_ keeps its
// size, so nothing here allocates.
void FlowLedger::install(std::vector<Amount> baseline, std::vector<Amount> reservation) {
  carried_.assign(baseline.size(), 0);
  reserved_.clear();
  service_baseline_ = std::move(baseline);
  reservation_baseline_ = std::move(reservation);
  entries_.clear();
  policy_bindings_.clear();
}

namespace detail {

namespace {

bool is_bound(const auto& owners, std::int64_t owner) {
  const auto it = owners.find(owner);
  return it != owners.end() && it->second.bound;
}

// Per-thread deltas indexed by EdgeId. Each commit advances an epoch (a numeric stamp).
// A matching slot stamp means the edge is already in the touched list, even if its delta
// returned to zero (-1, then +1), so each edge is published once. Reset a slot on first touch;
// later commits ignore old stamps, including after an exception. Only stamp wraparound (after
// 2^32 commits on one thread, which a long-running process reaches) needs a full stamp clear.
// Delta and stamp share one array: failed growth leaves the old array intact (the strong
// guarantee). One scratch per projection (Projection 0: carried loads, 1: holds).
template <int Projection>
class DeltaScratch {
public:
  explicit DeltaScratch(std::size_t num_edges) : storage_(storage()) {
    storage_.touched.clear();
    if (storage_.slot.size() < num_edges) storage_.slot.resize(num_edges);  // strong guarantee; new stamps 0
    if (++storage_.generation == 0) {
      for (auto& x : storage_.slot) x.stamp = 0;
      storage_.generation = 1;
    }
  }
  DeltaScratch(const DeltaScratch&) = delete;
  DeltaScratch& operator=(const DeltaScratch&) = delete;
  ~DeltaScratch() { storage_.touched.clear(); }
  void add(EdgeId e, Amount x) {
    auto& slot = storage_.slot[static_cast<std::size_t>(e)];
    if (slot.stamp != storage_.generation) {
      storage_.touched.push_back(e);  // may throw: the slot is left untouched in this commit
      slot.stamp = storage_.generation;
      slot.delta = 0;
    }
    slot.delta += x;
  }
  // Each touched edge once, in first-touch order.
  [[nodiscard]] const std::vector<EdgeId>& touched() const noexcept { return storage_.touched; }
  [[nodiscard]] Amount at(EdgeId e) const noexcept {
    const auto& slot = storage_.slot[static_cast<std::size_t>(e)];
    return slot.stamp == storage_.generation ? slot.delta : 0;
  }

private:
  struct Slot {
    Amount delta {0};
    std::uint32_t stamp {0};
  };
  struct Storage {
    std::vector<Slot> slot;
    std::vector<EdgeId> touched;
    std::uint32_t generation {0};
  };
  static Storage& storage() {
    thread_local Storage st;
    return st;
  }
  Storage& storage_;
};

// Entries publish by swap, validated scalar insertion or staged batch-node transfer.
static_assert(std::is_nothrow_move_assignable_v<LedgerEntry>);
static_assert(std::is_nothrow_swappable_v<LedgerEntry>);

} // namespace

#if NETGRAPH_CORE_PRIVATE_HOOKS
EntryWork& thread_entry_work() noexcept { thread_local EntryWork work; return work; }
#endif

std::vector<FlowRecord> LedgerEntry::records(const Graph& graph, std::uint64_t uid, std::size_t first, std::size_t count) const {
  const auto& identities = record_inventory();
  const auto& targets = destination_nodes();
  std::vector<FlowRecord> out;
  out.reserve(count);
  for (std::size_t i = first; i < first + count; ++i) {
    const auto& identity = identities[i];
    FlowRecord r;
    r.flow_id = {uid, first_counter + static_cast<std::int64_t>(i)};
    r.policy_id = owner;
    r.sources = {identity.source};
    r.destinations = targets;
    r.admission = admission;
    r.entry_id = identity.entry_id;
    if (!derived) {
      const auto& d = stored[i];
      r.route = d.route;
      r.min_route_cost = d.min_route_cost;
      r.offered = d.offered;
      r.delivered = d.delivered;
      r.edge_flows = d.edge_flows;
      r.edge_drops = d.edge_drops;
      r.edge_reservations = d.edge_reservations;
    } else {
      const auto& inputs = *derived->inputs;
      std::vector<EdgeId> route;
      if (bit_at(inputs.roster, i)) {
        const auto& basis = *ProgramAccess::data(*inputs.program).common_replay_support;
        route = reachable_route(graph, basis, identity.source, targets, inputs.basis);
        if (!inputs.basis.empty()) route = prune(graph, PredecessorDAGBuilder::of_edges(graph, route), identity.source, targets);
      }
      r.route = std::make_shared<const PredecessorDAG>(PredecessorDAGBuilder::of_edges(graph, route));
      r.min_route_cost = route.empty() ? kUnreachable : inputs.basis.empty() ?
          ProgramAccess::data(*inputs.program).entries[i].min_route_cost : min_route_cost(graph, *r.route, identity.source, targets);
      if (admission == Admission::Lossy) r.offered = derived->inputs->plan->sources.size() <= 1 ? derived->rows[i].offered : Counter128 {};
#if NETGRAPH_CORE_PRIVATE_HOOKS
      ++thread_entry_work().routes;
#endif
    }
    out.push_back(std::move(r));
  }
  if (!derived) return out;
  const auto& plan = *derived->inputs->plan;
  for (auto& r : out) {
    if (plan.sources.size() != 1 || r.sources.front() != plan.sources.front()) continue;
    r.edge_flows = carried;
    r.edge_drops = drops;
    r.delivered = delivered;
  }
  if (plan.sources.size() <= 1) return out;
  // Each forwarding is derived once for a bulk read; only the requested record DAGs were built.
  for (auto history = derived->history; history; history = history->previous) {
    std::vector<Amount> inject(plan.sources.size(), 0);
    if (history->offer) round_robin_increments(history->offset, *history->offer, inject);
    for (std::size_t h = 0; h < inject.size(); ++h) if (bit_at(history->dark_heads, h)) inject[h] = 0;
    Attribution attr;
    (void)forward(plan, Injection {inject, !history->offer.has_value()}, history->carried, &attr);
#if NETGRAPH_CORE_PRIVATE_HOOKS
    ++thread_entry_work().derivations;
#endif
    if (history->offer) round_robin_increments(history->offset, *history->offer, inject);
    for (auto& r : out) {
      const auto it = std::lower_bound(plan.sources.begin(), plan.sources.end(), r.sources.front());
      if (it == plan.sources.end() || *it != r.sources.front()) continue;
      const auto h = static_cast<std::size_t>(it - plan.sources.begin());
      r.delivered += attr.delivered[h];
      if (r.offered) {
        accumulate_counter(*r.offered, Counter128 {0, static_cast<std::uint64_t>(!history->offer ? attr.offered[h] :
            !history->raw_head_offers && bit_at(history->dark_heads, h) ? 0 : inject[h])}, "offered");
      }
      for (const auto& x : attr.entries[h]) {
        const EdgeId e = plan.forwarding_edges[static_cast<std::size_t>(x.edge_position)];
        if (x.carried != 0) r.edge_flows.emplace_back(e, x.carried);
        if (x.drops != 0) r.edge_drops.emplace_back(e, Counter128 {0, static_cast<std::uint64_t>(x.drops)});
      }
    }
  }
  // A single attribution already has unique entries sorted by member position.
  if (!derived->history || !derived->history->previous) return out;
  for (auto& r : out) {
    auto merge = [](auto& pairs) {
      std::sort(pairs.begin(), pairs.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
      std::size_t n = 0;
      for (const auto& x : pairs) {
        if (n != 0 && pairs[n - 1].first == x.first) {
          if constexpr (std::is_same_v<std::decay_t<decltype(x.second)>, Counter128>) accumulate_counter(pairs[n - 1].second, x.second, "drops");
          else pairs[n - 1].second += x.second;
        } else pairs[n++] = x;
      }
      pairs.resize(n);
    };
    merge(r.edge_flows); merge(r.edge_drops);
  }
  return out;
}

LedgerEntry derived_entry(const Graph&, std::int64_t owner, Admission admission, std::span<const NodeId> destinations,
    const ForwardingView& view, bool entry_records, const LedgerEntry* previous, std::optional<Amount> offer,
    Counter128 offset, ProgramFlow& flow) {
  LedgerEntry out;
  out.owner = owner;
  out.admission = admission;
  const auto& p = *view.data;
  const auto& plan = *view.prepared_route(0);
  out.derived = previous ? *previous->derived : DerivedRecords {};
  auto& d = *out.derived;
  if (!previous) {
    auto inputs = std::make_shared<DerivedInputs>();
    inputs->program = view.program;
    inputs->plan = view.prepared_route(0);
    inputs->basis = view.basis_;
    inputs->destinations.assign(destinations.begin(), destinations.end());
    inputs->inventory.reserve(p.entries.size());
    for (const auto& e : p.entries) inputs->inventory.push_back({e.head, entry_records ? std::optional(e.entry_id) : std::nullopt});
    inputs->roster.assign((inputs->inventory.size() + 63) / 64, 0);
    for (std::size_t i = 0; i < inputs->inventory.size(); ++i) if (view.status(i) != TunnelStatus::Failed) set_bit(inputs->roster, i);
    d.inputs = std::move(inputs);
    if (admission == Admission::Lossy) d.rows.resize(out.record_inventory().size());
  }
  const auto& members = plan.forwarding_edges;
  // Joint entry results become member-aligned once; every later consumer uses that basis.
  if (flow.carried.empty()) {
    flow.carried.assign(members.size(), 0);
    if (admission == Admission::Lossy) flow.drops.assign(members.size(), 0);
    auto position = [&](EdgeId e) { return static_cast<std::size_t>(std::lower_bound(members.begin(), members.end(), e) - members.begin()); };
    for (const auto& e : flow.entries) {
      for (const auto& [edge, x] : e.edge_flows) flow.carried[position(edge)] += x;
      for (const auto& [edge, x] : e.edge_drops) flow.drops[position(edge)] += x;
    }
  }
  const auto* before = admission == Admission::Lossy ? previous : nullptr;
  out.delivered = flow.delivered + (before ? before->delivered : 0);
  const auto nonzero = [](const auto& values) { return static_cast<std::size_t>(std::count_if(values.begin(), values.end(), [](Amount x) { return x != 0; })); };
  out.carried.reserve(std::min(members.size(), nonzero(flow.carried) + (before ? before->carried.size() : 0)));
  out.drops.reserve(nonzero(flow.drops) + (before ? before->drops.size() : 0));
  std::size_t cpos = 0, dpos = 0;
  for (std::size_t m = 0; m < members.size(); ++m) {
    const EdgeId e = members[m];
    Amount c = flow.carried[m];
    Counter128 drop {0, static_cast<std::uint64_t>(flow.drops.empty() ? 0 : flow.drops[m])};
    if (before) {
      while (cpos < before->carried.size() && before->carried[cpos].first < e) ++cpos;
      if (cpos < before->carried.size() && before->carried[cpos].first == e) c += before->carried[cpos++].second;
      while (dpos < before->drops.size() && before->drops[dpos].first < e) out.drops.push_back(before->drops[dpos++]);
      if (dpos < before->drops.size() && before->drops[dpos].first == e) accumulate_counter(drop, before->drops[dpos++].second, "drops");
    }
    if (c != 0) out.carried.emplace_back(e, c);
    if (!drop.is_zero()) out.drops.emplace_back(e, drop);
  }
  if (before) out.drops.insert(out.drops.end(), before->drops.begin() + static_cast<std::ptrdiff_t>(dpos), before->drops.end());
  std::shared_ptr<ForwardingHistory> history;
  if (plan.sources.size() > 1) {
    history = std::make_shared<ForwardingHistory>();
    history->offer = admission == Admission::Lossy ? offer : std::optional(flow.admitted);
    history->offset = offset;
    history->raw_head_offers = !entry_records;
    history->dark_heads.assign((view.num_sources() + 63) / 64, 0);
    for (std::size_t h = 0; h < view.num_sources(); ++h) if (view.is_blackholed_source(h)) set_bit(history->dark_heads, h);
    history->carried = flow.carried;
    if (admission == Admission::Lossy) {
      history->drops = flow.drops;
      if (history->drops.empty()) history->drops.assign(members.size(), 0);
    }
  }
  if (admission == Admission::Lossy) {
    for (std::size_t i = 0; i < flow.entries.size(); ++i) {
      const auto& f = flow.entries[i];
      const auto it = std::lower_bound(plan.sources.begin(), plan.sources.end(), out.record_inventory()[i].source);
      const auto h = static_cast<std::size_t>(it - plan.sources.begin());
      const Amount raw = it != plan.sources.end() && *it == out.record_inventory()[i].source ? flow.offered_by_source[h] : 0;
      const Amount lost = !entry_records ? raw - f.assigned : 0;
      accumulate_counter(d.rows[i].offered, Counter128 {0, static_cast<std::uint64_t>(f.assigned + lost)}, "offered");
      accumulate_counter(d.rows[i].source_drops, Counter128 {0, static_cast<std::uint64_t>(f.source_drops + lost)}, "drops");
      d.rows[i].delivered += f.delivered;
    }
    if (history) history->previous = d.history;
  }
  d.history = std::move(history);
  return out;
}

LedgerEntry stored_entry(const Graph& graph, std::int64_t owner, std::span<FlowRecord> records) {
  LedgerEntry out;
  out.owner = owner;
  if (records.empty()) return out;
  out.admission = records.front().admission;
  out.destinations = records.front().destinations;
  out.inventory.reserve(records.size());
  out.stored.reserve(records.size());
  std::optional<DeltaScratch<0>> carried;
  if (records.size() > 1) carried.emplace(static_cast<std::size_t>(graph.num_edges()));
  std::optional<DeltaScratch<1>> held;
  std::map<EdgeId, Counter128> drops;
  auto compact = [](auto& pairs, const char* label, auto add) {
    EdgeId previous = -1;
    std::size_t n = 0;
    for (const auto& [e, x] : pairs) {
      bool nonzero;
      if constexpr (std::is_same_v<std::remove_cv_t<decltype(x)>, Amount>) {
        if (x < 0) throw std::logic_error(std::string("ledger: a negative ") + label + " in a commit");
        nonzero = x != 0;
      } else nonzero = !x.is_zero();
      if (e <= previous) throw std::logic_error(std::string("ledger: noncanonical ") + label + " in a commit");
      previous = e;
      add(e, x);
      if (nonzero) pairs[n++] = {e, x};
    }
    pairs.resize(n);
  };
  for (auto& r : records) {
    if (r.sources.size() != 1 || r.destinations != out.destinations || r.admission != out.admission) {
      throw std::logic_error("ledger: inconsistent record inventory");
    }
    compact(r.edge_flows, "allocation", [&](EdgeId e, Amount x) { if (carried) carried->add(e, x); });
    if (carried && !r.edge_reservations.empty() && !held) held.emplace(static_cast<std::size_t>(graph.num_edges()));
    compact(r.edge_reservations, "hold", [&](EdgeId e, Amount x) { if (held) held->add(e, x); });
    compact(r.edge_drops, "drops", [&](EdgeId e, Counter128 x) { if (carried) accumulate_counter(drops[e], x, "drops"); });
    out.delivered += r.delivered;
    out.inventory.push_back({r.sources.front(), r.entry_id});
    if (!carried) {
      out.carried = r.edge_flows;
      out.holds = r.edge_reservations;
      out.drops = r.edge_drops;
    }
    out.stored.push_back({std::move(r.route), r.min_route_cost, r.offered, r.delivered,
                          std::move(r.edge_flows), std::move(r.edge_drops), std::move(r.edge_reservations)});
  }
  auto finish = [](const auto& scratch, auto& result) {
    auto edges = scratch.touched();
    std::sort(edges.begin(), edges.end());
    result.reserve(edges.size());
    for (EdgeId e : edges) if (const auto x = scratch.at(e); x != 0) result.emplace_back(e, x);
  };
  if (carried) finish(*carried, out.carried);
  if (held) finish(*held, out.holds);
  for (const auto& [e, x] : drops) if (!x.is_zero()) out.drops.emplace_back(e, x);
  return out;
}

std::int64_t LedgerAccess::bind(FlowLedger& ledger) {
  ledger.check_usable();
  if (ledger.next_policy_id_ == std::numeric_limits<std::int64_t>::max()) throw std::runtime_error("policy identity limit exceeded");
  const std::int64_t owner = ledger.next_policy_id_;
  ledger.policy_bindings_.try_emplace(owner);  // the only allocation: nothing changes if it throws
  ++ledger.next_policy_id_;
  ++ledger.num_bound_policies_;
  return owner;
}

void LedgerAccess::release(FlowLedger& ledger, std::int64_t owner) noexcept {
  const auto guard = lock(ledger);
  const auto it = ledger.policy_bindings_.find(owner);
  if (it == ledger.policy_bindings_.end() || !it->second.bound) return;
  --ledger.num_bound_policies_;
  if (it->second.first_counter >= 0) it->second.bound = false;
  else ledger.policy_bindings_.erase(it);
}

const LedgerEntry* LedgerAccess::entry(const FlowLedger& ledger, std::int64_t owner) {
  ledger.check_usable();
  const auto it = ledger.policy_bindings_.find(owner);
  return it == ledger.policy_bindings_.end() || it->second.first_counter < 0 ? nullptr : ledger.entries_.at(it->second.first_counter).get();
}

FlowRecord LedgerAccess::materialize(const FlowLedger& ledger, std::int64_t counter) {
  ledger.check_usable();
  auto it = ledger.entries_.upper_bound(counter);
  if (it != ledger.entries_.begin()) {
    --it;
    const auto offset = static_cast<std::uint64_t>(counter - it->first);
    if (offset < it->second->record_inventory().size()) return it->second->materialize(*ledger.graph_, ledger.uid_, static_cast<std::size_t>(offset));
  }
  throw std::invalid_argument("FlowId (" + std::to_string(ledger.uid_) + ", " + std::to_string(counter) + ") is retired or unknown");
}

std::vector<FlowRecord> LedgerAccess::records(const FlowLedger& ledger, std::int64_t owner) {
  std::vector<FlowRecord> out;
  if (const auto* own = entry(ledger, owner)) {
    out = own->records(*ledger.graph_, ledger.uid_, 0, own->record_inventory().size());
  }
  return out;
}

void LedgerAccess::detach(FlowLedger& ledger, std::int64_t owner) {
  ledger.check_usable();
  if (!is_bound(ledger.policy_bindings_, owner)) throw std::logic_error("ledger: detach by an unbound owner");
  const auto it = ledger.policy_bindings_.find(owner);
  if (it->second.first_counter < 0) return;
  const auto current = ledger.entries_.find(it->second.first_counter);
  for (const auto& [e, x] : current->second->carried) ledger.carried_[static_cast<std::size_t>(e)] -= x;
  for (const auto& [e, x] : current->second->holds) ledger.reserved_[static_cast<std::size_t>(e)] -= x;
  ledger.entries_.erase(current);
  it->second.first_counter = -1;
}

void LedgerAccess::commit(FlowLedger& ledger, std::int64_t owner, LedgerEntry& value) {
  OwnerEntry one {owner, &value, false};
  commit_batch(ledger, std::span<OwnerEntry>(&one, 1));
}

void LedgerAccess::exchange(FlowLedger& ledger, std::int64_t owner, LedgerEntry& value) {
  OwnerEntry one {owner, &value, true};
  commit_batch(ledger, std::span<OwnerEntry>(&one, 1));
}

void LedgerAccess::commit_batch(FlowLedger& ledger, std::span<OwnerEntry> owners) {
  ledger.check_usable();
  for (const auto& o : owners) {
    if (!is_bound(ledger.policy_bindings_, o.owner)) throw std::logic_error("ledger: commit by an unbound owner");
  }
  if (owners.size() > 1) {
    std::vector<std::int64_t> ids;
    for (const auto& o : owners) ids.push_back(o.owner);
    std::sort(ids.begin(), ids.end());
    if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) throw std::logic_error("ledger: an owner twice in a batch");
  }
  struct Prepared {
    OwnerEntry* update {nullptr};
    const LedgerEntry* value {nullptr};
    decltype(ledger.policy_bindings_)::iterator owner;
    decltype(ledger.entries_)::iterator old;
    std::int64_t key {-1};
  };
  Prepared single;
  std::optional<std::vector<Prepared>> multiple;
  if (owners.size() > 1) multiple.emplace(owners.size());
  const std::span<Prepared> prepared = owners.size() == 1 ? std::span(&single, 1) :
      multiple ? std::span(*multiple) : std::span<Prepared> {};
  std::uint64_t count = 0;
  for (std::size_t i = 0; i < owners.size(); ++i) {
    const auto& o = owners[i];
    auto& s = prepared[i];
    s.update = &owners[i];
    s.value = o.entry;
    s.owner = ledger.policy_bindings_.find(o.owner);
    s.old = ledger.entries_.find(s.owner->second.first_counter);
    const auto* old = s.old == ledger.entries_.end() ? nullptr : s.old->second.get();
    if (o.retire || old == nullptr) count += o.entry->record_inventory().size();
    else if (old->record_inventory().size() != o.entry->record_inventory().size()) throw std::logic_error("ledger: a commit must keep its entry inventory");
  }
  if (count > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max() - ledger.next_record_counter_)) {
    throw std::runtime_error("record identity limit exceeded");
  }
  bool holds = false;
  for (const auto& s : prepared) {
    holds = holds || !s.value->holds.empty() || (s.old != ledger.entries_.end() && !s.old->second->holds.empty());
  }
  const bool direct = prepared.size() == 1 && !holds;
  std::optional<DeltaScratch<0>> delta;
  if (!direct) delta.emplace(ledger.carried_.size());
  bool negative = false;
  for (const auto& s : prepared) {
    if (delta && s.old != ledger.entries_.end()) for (const auto& [e, x] : s.old->second->carried) delta->add(e, -x);
    for (const auto& [e, x] : s.value->carried) { negative = negative || x < 0; if (delta) delta->add(e, x); }
  }
  if (negative) throw std::logic_error("ledger: a negative allocation in a commit");
  // Canonical C is strictly sorted and nonzero. A scalar without holds streams its
  // old/new merge directly; batches and holds keep the shared first-touch scratch.
  auto changes = [&](auto visit) {
    if (delta) {
      for (EdgeId e : delta->touched()) visit(e, delta->at(e), false);
      return;
    }
    const auto& s = prepared.front();
    const std::span<const std::pair<EdgeId, Amount>> old = s.old == ledger.entries_.end() ?
        std::span<const std::pair<EdgeId, Amount>> {} : s.old->second->carried;
    const auto& now = s.value->carried;
    std::size_t i = 0, j = 0;
    while (i < old.size() || j < now.size()) {
      const bool was = i < old.size() && (j == now.size() || old[i].first <= now[j].first);
      const EdgeId e = was ? old[i].first : now[j].first;
      Amount x = was ? -old[i++].second : 0;
      if (j < now.size() && now[j].first == e) x += now[j++].second;
      visit(e, x, was);
    }
  };
  auto infeasible = [&](EdgeId e, Amount next) {
    throw std::logic_error("ledger: infeasible commit on edge " + std::to_string(e) + " (" + std::to_string(next) + " of " + std::to_string(ledger.service_baseline_[static_cast<std::size_t>(e)]) + ")");
  };
  std::optional<std::pair<EdgeId, Amount>> first_new_failure;
  changes([&](EdgeId e, Amount x, bool was) {
    const auto ei = static_cast<std::size_t>(e);
    const Amount next = ledger.carried_[ei] + x;
    if (next < 0 || next > ledger.service_baseline_[ei]) {
      // First-touch validation checks every old edge before any new-only edge.
      if (!direct || was) infeasible(e, next);
      if (!first_new_failure) first_new_failure.emplace(e, next);
    }
  });
  if (first_new_failure) infeasible(first_new_failure->first, first_new_failure->second);
  std::optional<DeltaScratch<1>> held;
  if (holds) {
    held.emplace(ledger.carried_.size());
    bool negative_hold = false;
    for (const auto& s : prepared) {
      const auto& o = *s.update;
      if (s.old != ledger.entries_.end()) for (const auto& [e, x] : s.old->second->holds) held->add(e, -x);
      for (const auto& [e, x] : o.entry->holds) { negative_hold = negative_hold || x < 0; held->add(e, x); }
    }
    if (negative_hold) throw std::logic_error("ledger: a negative hold in a commit");
    for (EdgeId e : held->touched()) {
      const auto ei = static_cast<std::size_t>(e);
      const Amount next = ledger.reserved_at(ei) + held->at(e);
      if (next < 0 || next > ledger.reservation_at(ei)) {
        throw std::logic_error("ledger: infeasible reservation on edge " + std::to_string(e) + " (" + std::to_string(next) + " of " + std::to_string(ledger.reservation_at(ei)) + ")");
      }
    }
  }
  // All fallible preparation precedes publication; scalar insertion allocates before linking its node.
  if (held && ledger.reserved_.empty()) ledger.reserved_.assign(ledger.carried_.size(), 0);
  std::optional<decltype(ledger.entries_)> fresh;
  std::int64_t next = ledger.next_record_counter_;
  for (auto& s : prepared) {
    auto& o = *s.update;
    if (o.entry->record_inventory().empty()) continue;
    const bool created = s.old == ledger.entries_.end();
    s.key = !created && !o.retire ? s.old->first : next;
    if (created || o.retire) next += static_cast<std::int64_t>(o.entry->record_inventory().size());
    o.entry->owner = o.owner;
    o.entry->first_counter = s.key;
    if (created) {
      if (prepared.size() > 1 && !fresh) fresh.emplace();
      // A scalar inserts last: allocation failure changes nothing; success has no throwing work left.
      auto& destination = fresh ? *fresh : ledger.entries_;
      s.value = destination.emplace(s.key, std::make_unique<LedgerEntry>(std::move(*o.entry))).first->second.get();
    }
  }
  // No allocation or validation remains. Existing entry nodes and values are reused;
  // a full replacement rekeys its node to the newly reserved identity range.
  changes([&](EdgeId e, Amount x, bool) { ledger.carried_[static_cast<std::size_t>(e)] += x; });
  if (held) for (EdgeId e : held->touched()) ledger.reserved_[static_cast<std::size_t>(e)] += held->at(e);
  for (auto& s : prepared) {
    if (s.old != ledger.entries_.end()) {
      if (s.key < 0) ledger.entries_.erase(s.old);
      else {
        std::swap(*s.old->second, *s.update->entry);
        if (s.old->first != s.key) {
          auto node = ledger.entries_.extract(s.old);
          node.key() = s.key;
          ledger.entries_.insert(std::move(node));
        }
      }
    }
    s.owner->second.first_counter = s.key;
  }
  if (fresh) ledger.entries_.merge(*fresh);
  ledger.next_record_counter_ = next;
}

void LedgerAccess::residual_into(const FlowLedger& ledger, std::vector<Amount>& out) {
  ledger.check_usable();
  out.resize(ledger.service_baseline_.size());
  for (std::size_t e = 0; e < out.size(); ++e) out[e] = ledger.service_baseline_[e] - ledger.carried_[e];
}

std::span<const Amount> LedgerAccess::baseline(const FlowLedger& ledger) {
  ledger.check_usable();
  return ledger.service_baseline_;
}

std::span<const Amount> LedgerAccess::carried(const FlowLedger& ledger) {
  ledger.check_usable();
  return ledger.carried_;
}

std::span<const Amount> LedgerAccess::reservation_baseline(const FlowLedger& fg) {
  fg.check_usable();
  if (fg.reservation_baseline_.empty()) return fg.graph_->edge_capacity_units();
  return fg.reservation_baseline_;
}

std::span<const Amount> LedgerAccess::reserved(const FlowLedger& fg) {
  fg.check_usable();
  return fg.reserved_;
}

std::vector<Amount>& LedgerAccess::reserved_storage(FlowLedger& fg) {
  fg.check_usable();
  if (fg.reserved_.empty()) fg.reserved_.assign(fg.carried_.size(), 0);
  return fg.reserved_;
}

void LedgerAccess::set_next_ids(FlowLedger& ledger, std::int64_t next_record_counter, std::int64_t next_policy_id) {
  ledger.next_record_counter_ = next_record_counter;
  ledger.next_policy_id_ = next_policy_id;
}

std::int64_t LedgerAccess::next_record_counter(const FlowLedger& ledger) {
  ledger.check_usable();
  return ledger.next_record_counter_;
}

std::vector<std::int64_t> LedgerAccess::bound_policy_ids(const FlowLedger& ledger) {
  ledger.check_usable();
  std::vector<std::int64_t> out;
  for (const auto& [owner, state] : ledger.policy_bindings_) {
    if (state.bound) out.push_back(owner);
  }
  return out;
}

std::vector<std::int64_t> LedgerAccess::orphan_policy_ids(const FlowLedger& ledger) {
  ledger.check_usable();
  std::vector<std::int64_t> out;
  for (const auto& [owner, state] : ledger.policy_bindings_) {
    if (!state.bound) out.push_back(owner);
  }
  return out;
}

} // namespace detail

} // namespace netgraph::core
