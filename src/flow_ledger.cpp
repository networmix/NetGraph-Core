/*
  The integer ledger FlowLedger.

  The aggregate allocation and the aggregate holds are kept exactly (integers) and changed only by
  an atomic commit, a detach or a reset (which zeroes them); residual() and unreserved() are
  derived from them (each baseline minus its aggregate). The two projections are never added.
  Records are kept by id (FlowId order); ids come from a per-ledger counter and are never
  reused, also across reset().
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

namespace netgraph::core {

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
  records_ = other.records_;
  policy_record_counters_ = other.policy_record_counters_;
  for (auto& [id, rec] : records_) rec.flow_id = FlowId {uid_, id};
  next_record_counter_ = other.next_record_counter_;
  next_policy_id_ = other.next_policy_id_;
  policy_bindings_ = other.policy_bindings_;  // orphans only
}

FlowLedger::FlowLedger(FlowLedger&& other) {
  other.check_usable();
  if (other.num_bound_policies_ != 0) {
    throw std::invalid_argument("FlowLedger: a ledger with bound policies cannot be moved");
  }
  graph_ = std::move(other.graph_);
  uid_ = other.uid_;
  service_baseline_ = std::move(other.service_baseline_);
  carried_ = std::move(other.carried_);
  reservation_baseline_ = std::move(other.reservation_baseline_);
  reserved_ = std::move(other.reserved_);
  records_ = std::move(other.records_);
  policy_record_counters_ = std::move(other.policy_record_counters_);
  next_record_counter_ = other.next_record_counter_;
  next_policy_id_ = other.next_policy_id_;
  policy_bindings_ = std::move(other.policy_bindings_);
  other.usable_ = false;
  other.service_baseline_.clear();
  other.carried_.clear();
  other.reservation_baseline_.clear();
  other.reserved_.clear();
  other.records_.clear();
  other.policy_record_counters_.clear();
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
  check_usable();
  return uid_;
}

const std::shared_ptr<const Graph>& FlowLedger::graph() const {
  check_usable();
  return graph_;
}

std::span<const Amount> FlowLedger::edge_capacity_units() const {
  check_usable();
  return graph_->edge_capacity_units();
}

std::vector<Amount> FlowLedger::residual() const {
  check_usable();
  std::vector<Amount> out(service_baseline_.size());
  for (std::size_t e = 0; e < out.size(); ++e) out[e] = service_baseline_[e] - carried_[e];
  return out;
}

std::vector<Amount> FlowLedger::carried() const {
  check_usable();
  return carried_;
}

std::vector<Amount> FlowLedger::reserved() const {
  check_usable();
  return reserved_.empty() ? std::vector<Amount>(service_baseline_.size(), 0) : reserved_;
}

std::vector<Amount> FlowLedger::unreserved() const {
  check_usable();
  std::vector<Amount> out(service_baseline_.size());
  for (std::size_t e = 0; e < out.size(); ++e) out[e] = reservation_at(e) - reserved_at(e);
  return out;
}

std::vector<FlowRecord> FlowLedger::records() const {
  check_usable();
  std::vector<FlowRecord> out;
  out.reserve(records_.size());
  for (const auto& [id, rec] : records_) out.push_back(rec);
  return out;
}

FlowRecord FlowLedger::record(const FlowId& flow_id) const {
  check_usable();
  if (flow_id.first != uid_) {
    throw std::invalid_argument("FlowId (" + std::to_string(flow_id.first) + ", " + std::to_string(flow_id.second) +
                                ") belongs to another ledger (uid " + std::to_string(uid_) + ")");
  }
  const auto it = records_.find(flow_id.second);
  if (it == records_.end()) {
    throw std::invalid_argument("FlowId (" + std::to_string(flow_id.first) + ", " + std::to_string(flow_id.second) +
                                ") is retired or unknown");
  }
  return it->second;
}

std::vector<EdgeId> FlowLedger::saturated_boundary_edges(std::span<const NodeId> sources, const Eligibility* eligibility) const {
  check_usable();
  const auto& g = *graph_;
  const auto srcs = detail::node_set(sources, g.num_nodes(), "sources");
  if (eligibility != nullptr && eligibility->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  return detail::saturated_boundary_edges_of(g, residual(), srcs, eligibility);
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
  check_usable();
  reset(graph_->edge_capacity_units());
}

void FlowLedger::reset(std::span<const Amount> service_baseline) {
  check_usable();
  check_resettable();
  install(checked_baseline(service_baseline, "service_baseline"), {});
}

void FlowLedger::reset(std::span<const Amount> service_baseline, std::span<const Amount> reservation_baseline) {
  check_usable();
  check_resettable();
  auto baseline = checked_baseline(service_baseline, "service_baseline");
  install(std::move(baseline), checked_baseline(reservation_baseline, "reservation_baseline"));
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
  records_.clear();
  policy_record_counters_.clear();
  policy_bindings_.clear();
}

namespace detail {

namespace {

bool is_bound(const std::map<std::int64_t, bool>& owners, std::int64_t owner) {
  const auto it = owners.find(owner);
  return it != owners.end() && it->second;
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
  [[nodiscard]] Amount at(EdgeId e) const noexcept { return storage_.slot[static_cast<std::size_t>(e)].delta; }

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

FlowRecord record_of(std::uint64_t uid, std::int64_t id, std::int64_t owner, RecordDraft& d) {
  FlowRecord rec = std::move(d.record);
  rec.flow_id = FlowId {uid, id};
  rec.policy_id = owner;
  std::erase_if(rec.edge_flows, [](const auto& p) { return p.second == 0; });
  std::erase_if(rec.edge_drops, [](const auto& p) { return p.second.is_zero(); });
  std::erase_if(rec.edge_reservations, [](const auto& p) { return p.second == 0; });
  return rec;
}

// The publication step moves staged records into place; it must not throw.
static_assert(std::is_nothrow_move_assignable_v<FlowRecord>);
static_assert(std::is_nothrow_swappable_v<std::vector<std::int64_t>>);

} // namespace

std::int64_t LedgerAccess::bind(FlowLedger& ledger) {
  ledger.check_usable();
  if (ledger.next_policy_id_ == std::numeric_limits<std::int64_t>::max()) throw std::runtime_error("policy identity limit exceeded");
  const std::int64_t owner = ledger.next_policy_id_;
  ledger.policy_bindings_.emplace(owner, true);  // the only allocation: nothing changes if it throws
  ++ledger.next_policy_id_;
  ++ledger.num_bound_policies_;
  return owner;
}

void LedgerAccess::release(FlowLedger& ledger, std::int64_t owner) noexcept {
  const auto it = ledger.policy_bindings_.find(owner);
  if (it == ledger.policy_bindings_.end() || !it->second) return;
  --ledger.num_bound_policies_;
  // An owner with live records stays as their orphaned owner; one without leaves nothing to keep
  // (its id is never issued again either way). Neither step allocates.
  const auto records = ledger.policy_record_counters_.find(owner);
  if (records != ledger.policy_record_counters_.end() && !records->second.empty()) {
    it->second = false;
  } else {
    ledger.policy_bindings_.erase(it);
  }
}

std::vector<std::int64_t> LedgerAccess::detach(FlowLedger& ledger, std::int64_t owner) {
  ledger.check_usable();
  if (!is_bound(ledger.policy_bindings_, owner)) throw std::logic_error("ledger: detach by an unbound owner");
  std::vector<std::int64_t> ids;
  const auto it = ledger.policy_record_counters_.find(owner);
  if (it == ledger.policy_record_counters_.end()) return ids;
  // Nothing below allocates: the index vector moves out, its node and the records are erased.
  ids = std::move(it->second);
  ledger.policy_record_counters_.erase(it);
  for (std::int64_t id : ids) {
    const auto rec = ledger.records_.find(id);
    for (const auto& [e, x] : rec->second.edge_flows) ledger.carried_[static_cast<std::size_t>(e)] -= x;
    for (const auto& [e, x] : rec->second.edge_reservations) ledger.reserved_[static_cast<std::size_t>(e)] -= x;
    ledger.records_.erase(rec);
  }
  return ids;
}

std::vector<FlowId> LedgerAccess::commit(FlowLedger& ledger, std::int64_t owner, std::span<RecordDraft> drafts) {
  return commit_drafts(ledger, owner, drafts, false);
}

std::vector<FlowId> LedgerAccess::exchange(FlowLedger& ledger, std::int64_t owner, std::span<RecordDraft> drafts) {
  return commit_drafts(ledger, owner, drafts, true);
}

std::vector<FlowId> LedgerAccess::commit_drafts(FlowLedger& ledger, std::int64_t owner, std::span<RecordDraft> drafts,
                                                bool retire) {
  ledger.check_usable();
  if (!is_bound(ledger.policy_bindings_, owner)) throw std::logic_error("ledger: commit by an unbound owner");
  // Every new record's id, checked as one batch before anything changes: ids stay below the
  // largest int64.
  const auto num_new = static_cast<std::size_t>(
      std::count_if(drafts.begin(), drafts.end(), [](const RecordDraft& d) { return !d.record_counter.has_value(); }));
  if (num_new > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max() - ledger.next_record_counter_)) {
    throw std::runtime_error("record identity limit exceeded");
  }
  const auto owned = ledger.policy_record_counters_.find(owner);
  if (owned == ledger.policy_record_counters_.end() && drafts.size() == 1 && !drafts[0].record_counter.has_value() &&
      drafts[0].record.edge_reservations.empty() &&
      std::adjacent_find(drafts[0].record.edge_flows.begin(), drafts[0].record.edge_flows.end(),
                         [](const auto& a, const auto& b) { return a.first >= b.first; }) ==
          drafts[0].record.edge_flows.end()) {
    // One new record without holds for an owner without live records (a fresh placement), its
    // edges strictly ascending: each touched edge's delta is the draft's own amount, touched in the draft's order,
    // so the checks below run on the draft itself, in the general path's order and with its
    // messages (negative amounts first, then feasibility edge by edge); the staging and the
    // publication are the general path's for one new record.
    RecordDraft& d = drafts[0];
    for (const auto& [e, x] : d.record.edge_flows) {
      if (x < 0) throw std::logic_error("ledger: a negative allocation in a commit");
    }
    for (const auto& [e, x] : d.record.edge_flows) {
      const auto ei = static_cast<std::size_t>(e);
      const Amount next = ledger.carried_[ei] + x;
      if (next < 0 || next > ledger.service_baseline_[ei]) {
        throw std::logic_error("ledger: infeasible commit on edge " + std::to_string(e) + " (" +
                               std::to_string(next) + " of " + std::to_string(ledger.service_baseline_[ei]) + ")");
      }
    }
    std::vector<FlowId> created;
    created.reserve(1);
    std::vector<std::int64_t> live;
    live.reserve(1);
    const std::int64_t id = ledger.next_record_counter_;
    std::map<std::int64_t, FlowRecord> fresh;
    const auto node = fresh.emplace_hint(fresh.end(), id, record_of(ledger.uid_, id, owner, d));
    created.emplace_back(ledger.uid_, id);
    live.push_back(id);
    std::map<std::int64_t, std::vector<std::int64_t>> index_node;
    index_node.emplace(owner, std::move(live));
    // Publication: nothing below allocates or throws (the record's dropped zero entries add
    // nothing).
    for (const auto& [e, x] : node->second.edge_flows) ledger.carried_[static_cast<std::size_t>(e)] += x;
    ledger.records_.merge(fresh);
    ledger.policy_record_counters_.merge(index_node);
    ledger.next_record_counter_ = id + 1;
    return created;
  }
  std::vector<FlowId> created;
  OwnerDrafts one {owner, drafts, retire};
  commit_batch(ledger, std::span<OwnerDrafts>(&one, 1), std::span<std::vector<FlowId>>(&created, 1));
  return created;
}

void LedgerAccess::commit_batch(FlowLedger& ledger, std::span<OwnerDrafts> owners, std::span<std::vector<FlowId>> created) {
  ledger.check_usable();
  if (created.size() != owners.size()) throw std::logic_error("ledger: a batch needs one id list per owner");
  for (const auto& o : owners) {
    if (!is_bound(ledger.policy_bindings_, o.owner)) throw std::logic_error("ledger: commit by an unbound owner");
  }
  if (owners.size() > 1) {
    std::vector<std::int64_t> ids;
    ids.reserve(owners.size());
    for (const auto& o : owners) ids.push_back(o.owner);
    std::sort(ids.begin(), ids.end());
    if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) throw std::logic_error("ledger: an owner twice in a batch");
  }
  // Every new record's id, counted over the whole batch before anything changes: ids stay below
  // the largest int64.
  std::uint64_t num_new = 0;
  for (const auto& o : owners) {
    num_new += static_cast<std::uint64_t>(std::count_if(o.drafts.begin(), o.drafts.end(),
                                                        [](const RecordDraft& d) { return !d.record_counter.has_value(); }));
  }
  if (num_new > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max() - ledger.next_record_counter_)) {
    throw std::runtime_error("record identity limit exceeded");
  }
  // Per owner: its index entry and the live records its drafts keep. The drafts carry exactly the
  // owner's live records plus new ones (retire: only new ones, every live record is retired).
  struct OwnerState {
    std::map<std::int64_t, std::vector<std::int64_t>>::iterator owned;
    std::size_t num_current {0};
    std::size_t num_new {0};
    std::vector<std::int64_t> kept;
  };
  std::vector<OwnerState> state(owners.size());
  for (std::size_t i = 0; i < owners.size(); ++i) {
    const auto& o = owners[i];
    auto& st = state[i];
    st.owned = ledger.policy_record_counters_.find(o.owner);
    const bool indexed = st.owned != ledger.policy_record_counters_.end();
    st.num_current = indexed ? st.owned->second.size() : 0;
    for (const auto& d : o.drafts) {
      if (d.record_counter.has_value()) {
        st.kept.push_back(*d.record_counter);
      } else {
        ++st.num_new;
      }
    }
    std::sort(st.kept.begin(), st.kept.end());
    if (o.retire ? !st.kept.empty()
                 : (st.kept.size() != st.num_current || (st.num_current != 0 && st.kept != st.owned->second))) {
      throw std::logic_error(o.retire ? "ledger: an exchange carries only new records"
                                      : "ledger: a commit must carry every live record of its policy");
    }
  }

  // 1. The new aggregates on every touched edge (every owner's old allocation and holds out, its
  // new drafts in), checked before anything changes: carried loads against the service baseline,
  // then holds against the reservation baseline (each below 2^62, so a sum of two stays below
  // 2^63: every owner's records fit the baseline before and after). Holds are gathered only when
  // some record has one.
  DeltaScratch<0> delta(ledger.carried_.size());
  bool holds = false;
  bool negative = false;
  for (std::size_t i = 0; i < owners.size(); ++i) {
    if (state[i].num_current != 0) {
      for (std::int64_t id : state[i].owned->second) {
        const FlowRecord& r = ledger.records_.at(id);
        for (const auto& [e, x] : r.edge_flows) delta.add(e, -x);
        holds = holds || !r.edge_reservations.empty();
      }
    }
    for (const auto& d : owners[i].drafts) {
      for (const auto& [e, x] : d.record.edge_flows) {
        negative = negative || x < 0;
        delta.add(e, x);
      }
      holds = holds || !d.record.edge_reservations.empty();
    }
  }
  if (negative) throw std::logic_error("ledger: a negative allocation in a commit");
  for (EdgeId e : delta.touched()) {
    const auto ei = static_cast<std::size_t>(e);
    const Amount next = ledger.carried_[ei] + delta.at(e);
    if (next < 0 || next > ledger.service_baseline_[ei]) {
      throw std::logic_error("ledger: infeasible commit on edge " + std::to_string(e) + " (" + std::to_string(next) +
                             " of " + std::to_string(ledger.service_baseline_[ei]) + ")");
    }
  }
  std::optional<DeltaScratch<1>> held;
  if (holds) {
    held.emplace(ledger.carried_.size());
    bool negative_hold = false;
    for (std::size_t i = 0; i < owners.size(); ++i) {
      if (state[i].num_current != 0) {
        for (std::int64_t id : state[i].owned->second) {
          for (const auto& [e, x] : ledger.records_.at(id).edge_reservations) held->add(e, -x);
        }
      }
      for (const auto& d : owners[i].drafts) {
        for (const auto& [e, x] : d.record.edge_reservations) {
          negative_hold = negative_hold || x < 0;
          held->add(e, x);
        }
      }
    }
    if (negative_hold) throw std::logic_error("ledger: a negative hold in a commit");
    for (EdgeId e : held->touched()) {
      const auto ei = static_cast<std::size_t>(e);
      const Amount next = ledger.reserved_at(ei) + held->at(e);
      if (next < 0 || next > ledger.reservation_at(ei)) {
        throw std::logic_error("ledger: infeasible reservation on edge " + std::to_string(e) + " (" +
                               std::to_string(next) + " of " + std::to_string(ledger.reservation_at(ei)) + ")");
      }
    }
  }

  // 2. Every allocation of the publication, staged: the holds' storage on the first hold (all
  // zero, as the empty storage reads), the records (new ones as map nodes), the owners' index
  // nodes, the returned ids and the next id. A throw here changes nothing.
  if (held.has_value() && ledger.reserved_.empty()) ledger.reserved_.assign(ledger.carried_.size(), 0);
  std::vector<std::vector<std::int64_t>> live(owners.size());  // an owner's new index, when it gains ids
  std::size_t num_retired = 0;
  std::size_t num_kept = 0;
  for (std::size_t i = 0; i < owners.size(); ++i) {
    auto& st = state[i];
    created[i].clear();
    created[i].reserve(st.num_new);
    if (st.num_new != 0) {
      live[i].reserve(st.kept.size() + st.num_new);
      live[i].assign(st.kept.begin(), st.kept.end());  // ascending; new ids are larger than every kept id
    }
    if (owners[i].retire) num_retired += st.num_current;
    num_kept += st.kept.size();
  }
  // An exchange retires the owner's live records: their map nodes leave the map before the new
  // ones are spliced in (an erase neither allocates nor throws).
  std::vector<std::map<std::int64_t, FlowRecord>::iterator> retired;
  retired.reserve(num_retired);
  std::vector<std::map<std::int64_t, FlowRecord>::iterator> targets;
  std::vector<FlowRecord> updated;
  targets.reserve(num_kept);
  updated.reserve(num_kept);
  std::map<std::int64_t, FlowRecord> fresh;
  std::map<std::int64_t, std::vector<std::int64_t>> index_nodes;
  std::int64_t next_id = ledger.next_record_counter_;
  for (std::size_t i = 0; i < owners.size(); ++i) {
    const auto& o = owners[i];
    auto& st = state[i];
    if (o.retire && st.num_current != 0) {
      for (std::int64_t id : st.owned->second) retired.push_back(ledger.records_.find(id));
    }
    for (auto& d : o.drafts) {
      if (d.record_counter.has_value()) {
        targets.push_back(ledger.records_.find(*d.record_counter));
        updated.push_back(record_of(ledger.uid_, *d.record_counter, o.owner, d));
      } else {
        const std::int64_t id = next_id++;
        fresh.emplace_hint(fresh.end(), id, record_of(ledger.uid_, id, o.owner, d));
        created[i].emplace_back(ledger.uid_, id);
        live[i].push_back(id);
      }
    }
    if (st.owned == ledger.policy_record_counters_.end() && !live[i].empty()) {
      index_nodes.emplace_hint(index_nodes.end(), o.owner, std::move(live[i]));
    }
  }

  // 3. Publication: nothing below allocates or throws (map nodes are spliced, not copied).
  for (EdgeId e : delta.touched()) ledger.carried_[static_cast<std::size_t>(e)] += delta.at(e);
  if (held.has_value()) {
    for (EdgeId e : held->touched()) ledger.reserved_[static_cast<std::size_t>(e)] += held->at(e);
  }
  for (std::size_t i = 0; i < targets.size(); ++i) targets[i]->second = std::move(updated[i]);
  for (const auto& it : retired) ledger.records_.erase(it);
  ledger.records_.merge(fresh);
  for (std::size_t i = 0; i < owners.size(); ++i) {
    auto& st = state[i];
    if (st.owned == ledger.policy_record_counters_.end()) continue;
    if (st.num_new != 0) {
      st.owned->second.swap(live[i]);
    } else if (owners[i].retire) {
      ledger.policy_record_counters_.erase(st.owned);
    }
  }
  ledger.policy_record_counters_.merge(index_nodes);
  ledger.next_record_counter_ = next_id;
}

const FlowRecord& LedgerAccess::live(const FlowLedger& ledger, std::int64_t id) {
  const auto it = ledger.records_.find(id);
  if (it == ledger.records_.end()) throw std::invalid_argument("ledger: retired or unknown flow id " + std::to_string(id));
  return it->second;
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
  for (const auto& [owner, bound] : ledger.policy_bindings_) {
    if (bound) out.push_back(owner);
  }
  return out;
}

std::vector<std::int64_t> LedgerAccess::orphan_policy_ids(const FlowLedger& ledger) {
  ledger.check_usable();
  std::vector<std::int64_t> out;
  for (const auto& [owner, bound] : ledger.policy_bindings_) {
    if (!bound) out.push_back(owner);
  }
  return out;
}

} // namespace detail

} // namespace netgraph::core
