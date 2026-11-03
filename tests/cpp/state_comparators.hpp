// What a call leaves, as the tests compare it: the ledger (LedgerState, every record field), its
// invariant (expect_balanced), a policy's epoch state (Epoch) and a placement result; the
// comparators' negative controls intercept their failures with test_utils.hpp's failures_of.
// Shared by both test executables; the fault-injection sweeps (fault_injection.hpp) compare
// through these.
#pragma once

#include <gtest/gtest.h>

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "test_utils.hpp"

namespace netgraph::core::test {

// The ledger as the tests compare it: both projections (the aggregate and the residual, the
// holds and the unreserved capacity), every record field (FlowIds keep their local ids: a
// reference sequence runs on another ledger), the next FlowId and the bound and orphaned owners.
struct EntryProjection {
  std::int64_t owner;
  std::int64_t first;
  std::size_t count;
  std::vector<std::pair<EdgeId, Amount>> carried;
  std::vector<std::pair<EdgeId, Amount>> holds;
  Amount delivered;
  std::vector<std::pair<EdgeId, Counter128>> drops;
  std::vector<std::uint64_t> derived_values {};
  bool operator==(const EntryProjection&) const = default;
};

struct LedgerState {
  std::vector<Amount> carried;
  std::vector<Amount> residual;
  std::vector<Amount> reserved;
  std::vector<Amount> unreserved;
  std::vector<FlowRecord> records;
  std::int64_t next_id {0};
  std::vector<std::int64_t> bound;
  std::vector<std::int64_t> orphans;
  std::vector<EntryProjection> entries {};
};

// Every record's FlowId carries the ledger's uid; then the uid is cleared, so records of two
// ledgers compare by their local ids.
inline void expect_uid_then_clear(std::vector<FlowRecord>& records, std::uint64_t uid) {
  for (auto& r : records) {
    EXPECT_EQ(r.flow_id.first, uid);
    r.flow_id.first = 0;
  }
}

inline LedgerState ledger_state(const FlowLedger& fg) {
  using detail::LedgerAccess;
  LedgerState s {fg.carried(), fg.residual(), fg.reserved(), fg.unreserved(), fg.records(),
                 LedgerAccess::next_record_counter(fg), LedgerAccess::bound_policy_ids(fg), LedgerAccess::orphan_policy_ids(fg)};
  expect_uid_then_clear(s.records, fg.uid());
  auto owners = s.bound;
  owners.insert(owners.end(), s.orphans.begin(), s.orphans.end());
  std::sort(owners.begin(), owners.end());
  const auto gate = LedgerAccess::lock(fg);
  for (const auto owner : owners) {
    if (const auto* entry = LedgerAccess::entry(fg, owner)) {
      s.entries.push_back({owner, entry->first_counter, entry->record_inventory().size(), entry->aggregate_carried(),
                           entry->aggregate_holds(), entry->delivered, entry->aggregate_drops()});
      auto& values = s.entries.back().derived_values;
      values.push_back(entry->derived.has_value());
      if (entry->derived) {
        const auto& derived = *entry->derived;
        auto scalar = [&](auto x) { values.push_back(static_cast<std::uint64_t>(x)); };
        auto array = [&](const auto& a) { scalar(a.size()); for (auto x : a) scalar(x); };
        const auto& inputs = *derived.inputs;
        const auto& program = detail::ProgramAccess::data(*inputs.program);
        array(inputs.destinations);
        scalar(inputs.inventory.size());
        for (const auto& identity : inputs.inventory) { scalar(identity.source); scalar(identity.entry_id.has_value()); scalar(identity.entry_id.value_or(0)); }
        array(inputs.roster); array(inputs.basis); array(*program.common_replay_support);
        for (const auto& e : program.entries) { scalar(e.head); scalar(e.entry_id); scalar(e.min_route_cost); }
        const auto& plan = *inputs.plan;
        scalar(plan.allocation); array(plan.sources); array(plan.source_positions); array(plan.destinations);
        array(plan.support_edges); array(plan.forwarding_edges); array(plan.weights); array(plan.topological_order);
        array(plan.outgoing_begin); array(plan.outgoing_end); array(plan.edge_destination_positions); array(plan.is_destination);
        scalar(derived.rows.size());
        for (const auto& row : derived.rows) {
          scalar(row.offered.hi); scalar(row.offered.lo); scalar(row.delivered);
          scalar(row.source_drops.hi); scalar(row.source_drops.lo);
        }
        for (auto history = derived.history; history; history = history->previous) {
          scalar(1); scalar(history->offer.has_value()); scalar(history->offer.value_or(0));
          scalar(history->offset.hi); scalar(history->offset.lo); scalar(history->raw_head_offers);
          array(history->dark_heads); array(history->carried); array(history->drops);
        }
        scalar(0);
      }

    }
  }
  return s;
}

inline bool same_record(const FlowRecord& a, const FlowRecord& b) {
  return a.flow_id == b.flow_id && a.policy_id == b.policy_id && a.sources == b.sources && a.destinations == b.destinations &&
         a.admission == b.admission && a.entry_id == b.entry_id && a.route && b.route && *a.route == *b.route &&
         a.min_route_cost == b.min_route_cost && a.offered == b.offered && a.delivered == b.delivered && a.edge_flows == b.edge_flows &&
         a.edge_drops == b.edge_drops && a.edge_reservations == b.edge_reservations;
}

inline void expect_same(const LedgerState& a, const LedgerState& b, const std::string& where) {
  EXPECT_EQ(a.carried, b.carried) << where;
  EXPECT_EQ(a.residual, b.residual) << where;
  EXPECT_EQ(a.reserved, b.reserved) << where;
  EXPECT_EQ(a.unreserved, b.unreserved) << where;
  EXPECT_EQ(a.next_id, b.next_id) << where;
  EXPECT_EQ(a.bound, b.bound) << where;
  EXPECT_EQ(a.orphans, b.orphans) << where;
  EXPECT_TRUE(a.entries == b.entries) << where;
  ASSERT_EQ(a.records.size(), b.records.size()) << where;
  for (std::size_t i = 0; i < a.records.size(); ++i) {
    EXPECT_TRUE(same_record(a.records[i], b.records[i])) << where << " record " << i;
  }
}

// The ledger invariant, per projection: the aggregate is the sum of the records' edges and the
// holds the sum of their edge_reservations, each within [0, its baseline] on every edge.
inline void expect_balanced(const FlowLedger& fg, const std::string& where) {
  const std::vector<Amount> carried = fg.carried();
  const std::vector<Amount> reserved = fg.reserved();
  const auto baseline = detail::LedgerAccess::baseline(fg);
  const auto reservation = detail::LedgerAccess::reservation_baseline(fg);
  std::vector<Amount> sum(carried.size(), 0);
  std::vector<Amount> held(carried.size(), 0);
  for (const auto& r : fg.records()) {
    for (const auto& [e, x] : r.edge_flows) sum[static_cast<std::size_t>(e)] += x;
    for (const auto& [e, x] : r.edge_reservations) held[static_cast<std::size_t>(e)] += x;
  }
  const auto state = ledger_state(fg);
  for (const auto& entry : state.entries) {
    std::vector<Amount> own(carried.size(), 0), own_holds(carried.size(), 0);
    Amount delivered = 0;
    for (const auto& r : state.records) if (r.policy_id == entry.owner) {
      delivered += r.delivered;
      for (const auto& [e, x] : r.edge_flows) own[static_cast<std::size_t>(e)] += x;
      for (const auto& [e, x] : r.edge_reservations) own_holds[static_cast<std::size_t>(e)] += x;
    }
    std::vector<Amount> cached(carried.size(), 0), cached_holds(carried.size(), 0);
    for (const auto& [e, x] : entry.carried) cached[static_cast<std::size_t>(e)] = x;
    for (const auto& [e, x] : entry.holds) cached_holds[static_cast<std::size_t>(e)] = x;
    EXPECT_EQ(cached, own) << where << " owner " << entry.owner;
    EXPECT_EQ(cached_holds, own_holds) << where << " owner " << entry.owner;
    EXPECT_EQ(entry.delivered, delivered) << where << " owner " << entry.owner;
  }
  EXPECT_EQ(sum, carried) << where;
  EXPECT_EQ(held, reserved) << where;
  ASSERT_EQ(baseline.size(), carried.size()) << where;
  ASSERT_EQ(reservation.size(), carried.size()) << where;
  for (std::size_t e = 0; e < carried.size(); ++e) {
    EXPECT_GE(carried[e], 0) << where << " edge " << e;
    EXPECT_LE(carried[e], baseline[e]) << where << " edge " << e;
    EXPECT_GE(reserved[e], 0) << where << " edge " << e;
    EXPECT_LE(reserved[e], reservation[e]) << where << " edge " << e;
  }
}

// A program as the tests compare it: its endpoints, its tunnels, its source weights, and per split
// context its plan's roots, members and weights (a recomputed program is a new value).
struct ProgramSnapshot {
  bool present {false};
  std::vector<NodeId> sources;
  std::vector<NodeId> destinations;
  std::vector<Tunnel> tunnels;
  std::optional<std::vector<std::pair<NodeId, Amount>>> source_weights;
  std::vector<std::vector<EdgeId>> members;
  std::vector<std::vector<Amount>> weights;
  std::vector<std::vector<NodeId>> roots;
  std::vector<std::int64_t> structure;
  bool operator==(const ProgramSnapshot&) const = default;
};

inline ProgramSnapshot program_snapshot(const ForwardingProgram* p) {
  ProgramSnapshot s;
  if (p == nullptr) return s;
  s.present = true;
  s.sources.assign(p->sources().begin(), p->sources().end());
  s.destinations.assign(p->destinations().begin(), p->destinations().end());
  s.tunnels.assign(p->tunnels().begin(), p->tunnels().end());
  s.source_weights = p->source_weights();
  const auto& data = detail::ProgramAccess::data(*p);
  auto scalar = [&](auto value) { s.structure.push_back(static_cast<std::int64_t>(value)); };
  auto array = [&](const auto& values) {
    scalar(values.size());
    for (auto value : values) scalar(value);
  };
  scalar(data.pools_source_inflows);
  scalar(data.uses_topology_basis);
  scalar(data.common_replay_support != nullptr);
  if (data.common_replay_support) array(*data.common_replay_support);
  array(data.block());
  array(data.source_weight_values);
  scalar(data.entries.size());
  for (const auto& e : data.entries) {
    for (auto value : {e.entry_id, e.head, e.roster_index, e.local_index, e.tunnel_position, e.split_context_position}) scalar(value);
    scalar(e.weight);
    scalar(e.min_route_cost);
    scalar(e.route != nullptr);
    if (e.route) {
      array(e.route->parent_offsets());
      array(e.route->parents());
      array(e.route->via_edges());
    }
  }
  for (const auto& c : data.split_contexts) {
    s.members.push_back(c.prepared_route->forwarding_edges);
    s.weights.push_back(c.prepared_route->weights);
    s.roots.emplace_back(c.prepared_route->sources.begin(), c.prepared_route->sources.end());
    const auto& plan = *c.prepared_route;
    scalar(plan.allocation);
    array(plan.support_edges);
    array(plan.sources);
    array(plan.source_positions);
    array(plan.destinations);
    array(plan.topological_order);
    array(plan.outgoing_begin);
    array(plan.outgoing_end);
    array(plan.edge_destination_positions);
    array(plan.is_destination);
    array(c.source_entry_positions);
  }
  return s;
}

// Every field of a placement result.
inline void expect_same(const PlacementResult& a, const PlacementResult& b, const std::string& where) {
  EXPECT_EQ(a.call_delivered, b.call_delivered) << where;
  EXPECT_EQ(a.call_undelivered, b.call_undelivered) << where;
  EXPECT_EQ(a.delivery_delta, b.delivery_delta) << where;
  EXPECT_EQ(a.call_edge_drops, b.call_edge_drops) << where;
  EXPECT_EQ(a.call_source_drops, b.call_source_drops) << where;
  EXPECT_EQ(a.termination, b.termination) << where;
  EXPECT_EQ(a.routes_changed, b.routes_changed) << where;
  EXPECT_EQ(a.weights_changed, b.weights_changed) << where;
  EXPECT_EQ(a.reservations_changed, b.reservations_changed) << where;
  EXPECT_TRUE(program_snapshot(a.program.get()) == program_snapshot(b.program.get())) << where;
  EXPECT_TRUE(a.plan == b.plan) << where;
  EXPECT_EQ(a.requested, b.requested) << where;
  EXPECT_EQ(a.offered, b.offered) << where;
  EXPECT_EQ(a.admitted, b.admitted) << where;
  EXPECT_EQ(a.delivered, b.delivered) << where;
  EXPECT_EQ(a.unadmitted, b.unadmitted) << where;
  EXPECT_EQ(a.source_drops, b.source_drops) << where;
  EXPECT_EQ(a.edge_drops, b.edge_drops) << where;
  EXPECT_TRUE(a.by_source == b.by_source) << where;
  EXPECT_TRUE(a.by_tunnel == b.by_tunnel) << where;
  EXPECT_EQ(a.total_cost, b.total_cost) << where;
  EXPECT_EQ(a.marginal_cost_curve, b.marginal_cost_curve) << where;
  EXPECT_EQ(a.edge_flows, b.edge_flows) << where;
  EXPECT_EQ(a.residual, b.residual) << where;
  EXPECT_EQ(a.saturated_boundary_edges, b.saturated_boundary_edges) << where;
  EXPECT_EQ(a.edge_reservations, b.edge_reservations) << where;
  EXPECT_EQ(a.min_cut, b.min_cut) << where;
}

struct Epoch {
  Counter128 requested;
  Counter128 offered;
  bool infinite {false};
  bool prepared {false};
  Amount delivered {0};
  std::optional<std::vector<std::pair<Cost, Amount>>> curve;
  const ForwardingProgram* program_id {nullptr};
  ProgramSnapshot program;
  std::optional<TunnelPlan> plan;
  std::vector<EdgeId> support;
  std::vector<Counter128> source_drops_by_source;
  std::vector<TunnelFlow> by_tunnel;
  ProgramSnapshot kept_program;
  const ForwardingProgram* kept_id {nullptr};
};

inline Epoch epoch(const FlowPolicy& p) {
  const auto& pe = detail::PolicyAccess::engine(p);
  const auto& ep = pe.epoch;
  std::vector<EdgeId> support;
  if (ep.program && std::holds_alternative<ShortestPaths>(pe.config.routing)
      && std::get<ShortestPaths>(pe.config.routing).allocation == RouteAllocation::Optimized) {
    const auto& common = detail::ProgramAccess::data(*ep.program).common_replay_support;
    if (common) support = *common;
  }
  Epoch e {ep.requested, ep.offered, ep.infinite, ep.prepared, p.delivered(), ep.marginal_cost_curve, ep.program.get(),
           program_snapshot(ep.program.get()), ep.plan ? std::optional<TunnelPlan>(*ep.plan) : std::nullopt, support,
           ep.source_drops_by_source, ep.by_tunnel, program_snapshot(pe.kept_program.get()), pe.kept_program.get()};
  return e;
}

inline void expect_same(const Epoch& a, const Epoch& b, const std::string& where) {
  EXPECT_EQ(a.requested, b.requested) << where;
  EXPECT_EQ(a.offered, b.offered) << where;
  EXPECT_EQ(a.infinite, b.infinite) << where;
  EXPECT_EQ(a.prepared, b.prepared) << where;
  EXPECT_EQ(a.delivered, b.delivered) << where;
  EXPECT_EQ(a.curve, b.curve) << where;
  EXPECT_TRUE(a.program == b.program) << where;
  EXPECT_TRUE(a.plan == b.plan) << where;
  EXPECT_EQ(a.support, b.support) << where;
  EXPECT_EQ(a.source_drops_by_source, b.source_drops_by_source) << where;
  EXPECT_TRUE(a.by_tunnel == b.by_tunnel) << where;
  EXPECT_TRUE(a.kept_program == b.kept_program) << where;
}

// The same state on the same policy: also the same program object (a failed call installs none).
inline void expect_unchanged(const Epoch& a, const Epoch& b, const std::string& where) {
  expect_same(a, b, where);
  EXPECT_EQ(a.program_id, b.program_id) << where;
}

} // namespace netgraph::core::test
