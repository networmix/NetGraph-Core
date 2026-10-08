/*
  The integer ledger FlowLedger (APPROACH 8).

  The aggregate allocation is kept exactly (integers) and changed only by an atomic commit or a
  detach; residual() is derived from it (baseline - allocated). Records are kept by id (FlowId
  order); ids come from a per-ledger counter and are never reused, also across reset().
*/
#include "netgraph/core/flow_ledger.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/uid.hpp"

namespace netgraph::core {

FlowLedger::FlowLedger(std::shared_ptr<const Graph> graph) {
  if (!graph) throw std::invalid_argument("FlowLedger: graph is null");
  graph_ = std::move(graph);
  baseline_.assign(graph_->capacity_units().begin(), graph_->capacity_units().end());
  allocated_.assign(baseline_.size(), 0);
  uid_ = detail::next_uid();
}

FlowLedger::FlowLedger(std::shared_ptr<const Graph> graph, std::span<const Amount> residual) {
  if (!graph) throw std::invalid_argument("FlowLedger: graph is null");
  graph_ = std::move(graph);
  baseline_ = checked_baseline(residual);
  allocated_.assign(baseline_.size(), 0);
  uid_ = detail::next_uid();
}

FlowLedger::FlowLedger(const FlowLedger& other) {
  other.check_usable();
  if (other.num_bound_ != 0) {
    throw std::invalid_argument("FlowLedger: a ledger with bound policies cannot be copied");
  }
  graph_ = other.graph_;
  uid_ = detail::next_uid();
  baseline_ = other.baseline_;
  allocated_ = other.allocated_;
  records_ = other.records_;
  owner_ids_ = other.owner_ids_;
  for (auto& [id, rec] : records_) rec.flow_id = FlowId {uid_, id};
  next_id_ = other.next_id_;
  next_owner_ = other.next_owner_;
  owners_ = other.owners_;  // orphans only
}

FlowLedger::FlowLedger(FlowLedger&& other) {
  other.check_usable();
  if (other.num_bound_ != 0) {
    throw std::invalid_argument("FlowLedger: a ledger with bound policies cannot be moved");
  }
  graph_ = std::move(other.graph_);
  uid_ = other.uid_;
  baseline_ = std::move(other.baseline_);
  allocated_ = std::move(other.allocated_);
  records_ = std::move(other.records_);
  owner_ids_ = std::move(other.owner_ids_);
  next_id_ = other.next_id_;
  next_owner_ = other.next_owner_;
  owners_ = std::move(other.owners_);
  other.usable_ = false;
  other.baseline_.clear();
  other.allocated_.clear();
  other.records_.clear();
  other.owner_ids_.clear();
  other.owners_.clear();
}

void FlowLedger::check_usable() const {
  if (!usable_) throw std::logic_error("FlowLedger: a moved-from ledger is empty and unusable");
}

std::vector<Amount> FlowLedger::checked_baseline(std::span<const Amount> residual) const {
  const auto capacity = graph_->capacity_units();
  if (residual.size() != capacity.size()) {
    throw std::invalid_argument("residual has length " + std::to_string(residual.size()) + ", expected " +
                                std::to_string(capacity.size()) + " (one entry per edge)");
  }
  for (std::size_t e = 0; e < residual.size(); ++e) {
    if (residual[e] < 0 || residual[e] > capacity[e]) {
      throw std::invalid_argument("residual[" + std::to_string(e) + "] = " + std::to_string(residual[e]) +
                                  " is outside [0, capacity_units[" + std::to_string(e) + "] = " +
                                  std::to_string(capacity[e]) + "]");
    }
  }
  return std::vector<Amount>(residual.begin(), residual.end());
}

std::uint64_t FlowLedger::uid() const {
  check_usable();
  return uid_;
}

const std::shared_ptr<const Graph>& FlowLedger::graph() const {
  check_usable();
  return graph_;
}

std::span<const Amount> FlowLedger::capacity_units() const {
  check_usable();
  return graph_->capacity_units();
}

std::vector<Amount> FlowLedger::residual() const {
  check_usable();
  std::vector<Amount> out(baseline_.size());
  for (std::size_t e = 0; e < out.size(); ++e) out[e] = baseline_[e] - allocated_[e];
  return out;
}

std::vector<Amount> FlowLedger::allocated() const {
  check_usable();
  return allocated_;
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

std::vector<EdgeId> FlowLedger::bottleneck(std::span<const NodeId> sources, const Eligibility* elig) const {
  check_usable();
  const auto& g = *graph_;
  const auto srcs = detail::node_set(sources, g.num_nodes(), "sources");
  if (elig != nullptr && elig->graph_uid() != g.uid()) {
    throw std::invalid_argument("foreign Eligibility: built for another graph");
  }
  return detail::bottleneck_of(g, residual(), srcs, elig);
}

namespace detail {

std::vector<EdgeId> bottleneck_of(const Graph& g, std::span<const Amount> residual,
                                  std::span<const NodeId> sources, const Eligibility* elig) {
  const auto n = static_cast<std::size_t>(g.num_nodes());
  const auto m = static_cast<std::size_t>(g.num_edges());
  const auto edge_src = g.edge_src();
  const auto edge_dst = g.edge_dst();
  const auto capacity = g.capacity_units();
  const auto row = g.row_offsets();
  std::vector<std::uint8_t> seen(n, 0);
  std::vector<NodeId> stack;
  for (NodeId s : sources) {
    if (elig != nullptr && !elig->node_ok(s)) continue;
    if (seen[static_cast<std::size_t>(s)] != 0) continue;
    seen[static_cast<std::size_t>(s)] = 1;
    stack.push_back(s);
  }
  while (!stack.empty()) {
    const auto u = static_cast<std::size_t>(stack.back());
    stack.pop_back();
    for (auto e = row[u]; e < row[u + 1]; ++e) {
      const auto ei = static_cast<std::size_t>(e);
      if (elig != nullptr && !elig->edge_ok(e)) continue;
      if (residual[ei] < 1) continue;
      const auto v = static_cast<std::size_t>(edge_dst[ei]);
      if (seen[v] == 0) {
        seen[v] = 1;
        stack.push_back(static_cast<NodeId>(v));
      }
    }
  }
  std::vector<EdgeId> out;
  for (std::size_t e = 0; e < m; ++e) {
    if (seen[static_cast<std::size_t>(edge_src[e])] == 0 || seen[static_cast<std::size_t>(edge_dst[e])] != 0) continue;
    if (elig != nullptr && !elig->edge_ok(static_cast<EdgeId>(e))) continue;
    if (capacity[e] < 1) continue;
    out.push_back(static_cast<EdgeId>(e));
  }
  return out;
}

} // namespace detail

void FlowLedger::reset() {
  check_usable();
  reset(graph_->capacity_units());
}

void FlowLedger::reset(std::span<const Amount> residual) {
  check_usable();
  if (num_bound_ != 0) {
    throw std::invalid_argument("FlowLedger.reset is rejected while a policy is bound (" +
                                std::to_string(num_bound_) + " bound)");
  }
  auto baseline = checked_baseline(residual);
  baseline_ = std::move(baseline);
  allocated_.assign(baseline_.size(), 0);
  records_.clear();
  owner_ids_.clear();
  owners_.clear();
}

namespace detail {

namespace {

bool is_bound(const std::map<std::int64_t, bool>& owners, std::int64_t owner) {
  const auto it = owners.find(owner);
  return it != owners.end() && it->second;
}

// The dense per-edge delta of one commit, in a per-thread array of slots. An edge is touched in
// this commit when its slot carries the commit's epoch stamp: that test does not depend on the
// delta's value, so a delta that returns to zero during accumulation (-1, then +1) keeps its one
// entry in the touched list and is published once (review-1b-fixes N2). A slot is reset on its
// first touch in a commit, so nothing needs clearing between commits: a commit that throws
// leaves stale slots behind an old epoch. Delta and stamp share one array, so growing it either
// succeeds whole or leaves the previous array (review-1b-fixes N1).
class DeltaScratch {
public:
  explicit DeltaScratch(std::size_t num_edges) : s_(storage()) {
    s_.touched.clear();
    if (s_.slot.size() < num_edges) s_.slot.resize(num_edges);  // strong guarantee; new stamps 0
    if (++s_.epoch == 0) {
      for (auto& x : s_.slot) x.stamp = 0;
      s_.epoch = 1;
    }
  }
  DeltaScratch(const DeltaScratch&) = delete;
  DeltaScratch& operator=(const DeltaScratch&) = delete;
  ~DeltaScratch() { s_.touched.clear(); }
  void add(EdgeId e, Amount x) {
    auto& slot = s_.slot[static_cast<std::size_t>(e)];
    if (slot.stamp != s_.epoch) {
      s_.touched.push_back(e);  // may throw: the slot is left untouched in this commit
      slot.stamp = s_.epoch;
      slot.delta = 0;
    }
    slot.delta += x;
  }
  // Each touched edge once, in first-touch order.
  [[nodiscard]] const std::vector<EdgeId>& touched() const noexcept { return s_.touched; }
  [[nodiscard]] Amount at(EdgeId e) const noexcept { return s_.slot[static_cast<std::size_t>(e)].delta; }

private:
  struct Slot {
    Amount delta {0};
    std::uint32_t stamp {0};
  };
  struct Storage {
    std::vector<Slot> slot;
    std::vector<EdgeId> touched;
    std::uint32_t epoch {0};
  };
  static Storage& storage() {
    thread_local Storage st;
    return st;
  }
  Storage& s_;
};

FlowRecord record_of(std::uint64_t uid, std::int64_t id, std::int64_t owner, Draft& d) {
  FlowRecord rec = std::move(d.record);
  rec.flow_id = FlowId {uid, id};
  rec.owner = owner;
  std::erase_if(rec.edges, [](const auto& p) { return p.second == 0; });
  std::erase_if(rec.drops, [](const auto& p) { return p.second.is_zero(); });
  return rec;
}

// The publication step moves staged records into place; it must not throw.
static_assert(std::is_nothrow_move_assignable_v<FlowRecord>);
static_assert(std::is_nothrow_swappable_v<std::vector<std::int64_t>>);

} // namespace

std::int64_t LedgerAccess::bind(FlowLedger& fg) {
  fg.check_usable();
  const std::int64_t owner = fg.next_owner_;
  fg.owners_.emplace(owner, true);  // the only allocation: nothing changes if it throws
  ++fg.next_owner_;
  ++fg.num_bound_;
  return owner;
}

void LedgerAccess::release(FlowLedger& fg, std::int64_t owner) noexcept {
  const auto it = fg.owners_.find(owner);
  if (it == fg.owners_.end() || !it->second) return;
  it->second = false;
  --fg.num_bound_;
}

std::vector<std::int64_t> LedgerAccess::detach(FlowLedger& fg, std::int64_t owner) {
  fg.check_usable();
  if (!is_bound(fg.owners_, owner)) throw std::logic_error("ledger: detach by an unbound owner");
  std::vector<std::int64_t> ids;
  const auto it = fg.owner_ids_.find(owner);
  if (it == fg.owner_ids_.end()) return ids;
  // Nothing below allocates: the index vector moves out, its node and the records are erased.
  ids = std::move(it->second);
  fg.owner_ids_.erase(it);
  for (std::int64_t id : ids) {
    const auto rec = fg.records_.find(id);
    for (const auto& [e, x] : rec->second.edges) fg.allocated_[static_cast<std::size_t>(e)] -= x;
    fg.records_.erase(rec);
  }
  return ids;
}

std::vector<FlowId> LedgerAccess::commit(FlowLedger& fg, std::int64_t owner, std::span<Draft> drafts) {
  fg.check_usable();
  if (!is_bound(fg.owners_, owner)) throw std::logic_error("ledger: commit by an unbound owner");
  // The drafts carry exactly the owner's live records plus new ones.
  const auto owned = fg.owner_ids_.find(owner);
  const std::size_t num_current = owned == fg.owner_ids_.end() ? 0 : owned->second.size();
  if (owned == fg.owner_ids_.end() && drafts.size() == 1 && !drafts[0].id.has_value() &&
      std::adjacent_find(drafts[0].record.edges.begin(), drafts[0].record.edges.end(),
                         [](const auto& a, const auto& b) { return a.first >= b.first; }) ==
          drafts[0].record.edges.end()) {
    // One new record for an owner without live records (a fresh placement), its edges strictly
    // ascending: each touched edge's delta is the draft's own amount, touched in the draft's order,
    // so the checks below run on the draft itself, in the general path's order and with its
    // messages (negative amounts first, then feasibility edge by edge); the staging and the
    // publication are the general path's for one new record.
    Draft& d = drafts[0];
    for (const auto& [e, x] : d.record.edges) {
      if (x < 0) throw std::logic_error("ledger: a negative allocation in a commit");
    }
    for (const auto& [e, x] : d.record.edges) {
      const auto ei = static_cast<std::size_t>(e);
      const Amount next = fg.allocated_[ei] + x;
      if (next < 0 || next > fg.baseline_[ei]) {
        throw std::logic_error("ledger: infeasible commit on edge " + std::to_string(e) + " (" +
                               std::to_string(next) + " of " + std::to_string(fg.baseline_[ei]) + ")");
      }
    }
    std::vector<FlowId> created;
    created.reserve(1);
    std::vector<std::int64_t> live;
    live.reserve(1);
    const std::int64_t id = fg.next_id_;
    std::map<std::int64_t, FlowRecord> fresh;
    const auto node = fresh.emplace_hint(fresh.end(), id, record_of(fg.uid_, id, owner, d));
    created.emplace_back(fg.uid_, id);
    live.push_back(id);
    std::map<std::int64_t, std::vector<std::int64_t>> index_node;
    index_node.emplace(owner, std::move(live));
    // Publication: nothing below allocates or throws (the record's dropped zero entries add
    // nothing).
    for (const auto& [e, x] : node->second.edges) fg.allocated_[static_cast<std::size_t>(e)] += x;
    fg.records_.merge(fresh);
    fg.owner_ids_.merge(index_node);
    fg.next_id_ = id + 1;
    return created;
  }
  std::vector<std::int64_t> kept;
  std::size_t num_new = 0;
  for (const auto& d : drafts) {
    if (d.id.has_value()) {
      kept.push_back(*d.id);
    } else {
      ++num_new;
    }
  }
  std::sort(kept.begin(), kept.end());
  if (kept.size() != num_current || (num_current != 0 && kept != owned->second)) {
    throw std::logic_error("ledger: a commit must carry every live record of its policy");
  }

  // 1. The new aggregate on every touched edge (old owner allocation out, new drafts in),
  // checked before anything changes.
  DeltaScratch delta(fg.allocated_.size());
  if (num_current != 0) {
    for (std::int64_t id : owned->second) {
      for (const auto& [e, x] : fg.records_.at(id).edges) delta.add(e, -x);
    }
  }
  bool negative = false;
  for (const auto& d : drafts) {
    const FlowRecord& r = d.record;
    for (const auto& [e, x] : r.edges) {
      negative = negative || x < 0;
      delta.add(e, x);
    }
  }
  if (negative) throw std::logic_error("ledger: a negative allocation in a commit");
  for (EdgeId e : delta.touched()) {
    const auto ei = static_cast<std::size_t>(e);
    const Amount next = fg.allocated_[ei] + delta.at(e);
    if (next < 0 || next > fg.baseline_[ei]) {
      throw std::logic_error("ledger: infeasible commit on edge " + std::to_string(e) + " (" + std::to_string(next) +
                             " of " + std::to_string(fg.baseline_[ei]) + ")");
    }
  }

  // 2. Every allocation of the publication, staged: the records (new ones as map nodes), the
  // owner's index node, the returned ids and the next id. A throw here changes nothing.
  std::vector<FlowId> created;
  created.reserve(num_new);
  std::vector<std::int64_t> live;  // the owner's new index, when it gains ids
  if (num_new != 0) {
    live.reserve(kept.size() + num_new);
    live.assign(kept.begin(), kept.end());  // ascending; new ids are larger than every kept id
  }
  std::vector<std::map<std::int64_t, FlowRecord>::iterator> targets;
  std::vector<FlowRecord> updated;
  targets.reserve(kept.size());
  updated.reserve(kept.size());
  std::map<std::int64_t, FlowRecord> fresh;
  std::int64_t next_id = fg.next_id_;
  for (auto& d : drafts) {
    if (d.id.has_value()) {
      targets.push_back(fg.records_.find(*d.id));
      updated.push_back(record_of(fg.uid_, *d.id, owner, d));
    } else {
      const std::int64_t id = next_id++;
      fresh.emplace_hint(fresh.end(), id, record_of(fg.uid_, id, owner, d));
      created.emplace_back(fg.uid_, id);
      live.push_back(id);
    }
  }
  std::map<std::int64_t, std::vector<std::int64_t>> index_node;
  if (owned == fg.owner_ids_.end() && !live.empty()) index_node.emplace(owner, std::move(live));

  // 3. Publication: nothing below allocates or throws (map nodes are spliced, not copied).
  for (EdgeId e : delta.touched()) fg.allocated_[static_cast<std::size_t>(e)] += delta.at(e);
  for (std::size_t i = 0; i < targets.size(); ++i) targets[i]->second = std::move(updated[i]);
  fg.records_.merge(fresh);
  if (num_new != 0) {
    if (owned != fg.owner_ids_.end()) {
      owned->second.swap(live);
    } else {
      fg.owner_ids_.merge(index_node);
    }
  }
  fg.next_id_ = next_id;
  return created;
}

const FlowRecord& LedgerAccess::live(const FlowLedger& fg, std::int64_t id) {
  const auto it = fg.records_.find(id);
  if (it == fg.records_.end()) throw std::invalid_argument("ledger: retired or unknown flow id " + std::to_string(id));
  return it->second;
}

void LedgerAccess::residual_into(const FlowLedger& fg, std::vector<Amount>& out) {
  fg.check_usable();
  out.resize(fg.baseline_.size());
  for (std::size_t e = 0; e < out.size(); ++e) out[e] = fg.baseline_[e] - fg.allocated_[e];
}

std::span<const Amount> LedgerAccess::baseline(const FlowLedger& fg) {
  fg.check_usable();
  return fg.baseline_;
}

std::span<const Amount> LedgerAccess::allocated(const FlowLedger& fg) {
  fg.check_usable();
  return fg.allocated_;
}

std::int64_t LedgerAccess::next_flow_id(const FlowLedger& fg) {
  fg.check_usable();
  return fg.next_id_;
}

std::vector<std::int64_t> LedgerAccess::bound_owners(const FlowLedger& fg) {
  fg.check_usable();
  std::vector<std::int64_t> out;
  for (const auto& [owner, bound] : fg.owners_) {
    if (bound) out.push_back(owner);
  }
  return out;
}

std::vector<std::int64_t> LedgerAccess::orphan_owners(const FlowLedger& fg) {
  fg.check_usable();
  std::vector<std::int64_t> out;
  for (const auto& [owner, bound] : fg.owners_) {
    if (!bound) out.push_back(owner);
  }
  return out;
}

} // namespace detail

} // namespace netgraph::core
