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
struct LedgerState {
  std::vector<Amount> allocated;
  std::vector<Amount> residual;
  std::vector<Amount> reserved;
  std::vector<Amount> unreserved;
  std::vector<FlowRecord> records;
  std::int64_t next_id {0};
  std::vector<std::int64_t> bound;
  std::vector<std::int64_t> orphans;
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
  LedgerState s {fg.allocated(), fg.residual(), fg.reserved(), fg.unreserved(), fg.records(),
                 LedgerAccess::next_flow_id(fg), LedgerAccess::bound_owners(fg), LedgerAccess::orphan_owners(fg)};
  expect_uid_then_clear(s.records, fg.uid());
  return s;
}

inline bool same_record(const FlowRecord& a, const FlowRecord& b) {
  return a.flow_id == b.flow_id && a.owner == b.owner && a.sources == b.sources && a.destinations == b.destinations &&
         a.admission == b.admission && a.entry == b.entry && a.bundle && b.bundle && *a.bundle == *b.bundle &&
         a.route_cost == b.route_cost && a.offered == b.offered && a.delivered == b.delivered && a.edges == b.edges &&
         a.drops == b.drops && a.reserved_edges == b.reserved_edges;
}

inline void expect_same(const LedgerState& a, const LedgerState& b, const std::string& where) {
  EXPECT_EQ(a.allocated, b.allocated) << where;
  EXPECT_EQ(a.residual, b.residual) << where;
  EXPECT_EQ(a.reserved, b.reserved) << where;
  EXPECT_EQ(a.unreserved, b.unreserved) << where;
  EXPECT_EQ(a.next_id, b.next_id) << where;
  EXPECT_EQ(a.bound, b.bound) << where;
  EXPECT_EQ(a.orphans, b.orphans) << where;
  ASSERT_EQ(a.records.size(), b.records.size()) << where;
  for (std::size_t i = 0; i < a.records.size(); ++i) {
    EXPECT_TRUE(same_record(a.records[i], b.records[i])) << where << " record " << i;
  }
}

// The ledger invariant, per projection: the aggregate is the sum of the records' edges and the
// holds the sum of their reserved_edges, each within [0, its baseline] on every edge.
inline void expect_balanced(const FlowLedger& fg, const std::string& where) {
  const std::vector<Amount> allocated = fg.allocated();
  const std::vector<Amount> reserved = fg.reserved();
  const auto baseline = detail::LedgerAccess::baseline(fg);
  const auto reservation = detail::LedgerAccess::reservation_baseline(fg);
  std::vector<Amount> sum(allocated.size(), 0);
  std::vector<Amount> held(allocated.size(), 0);
  for (const auto& r : fg.records()) {
    for (const auto& [e, x] : r.edges) sum[static_cast<std::size_t>(e)] += x;
    for (const auto& [e, x] : r.reserved_edges) held[static_cast<std::size_t>(e)] += x;
  }
  EXPECT_EQ(sum, allocated) << where;
  EXPECT_EQ(held, reserved) << where;
  ASSERT_EQ(baseline.size(), allocated.size()) << where;
  ASSERT_EQ(reservation.size(), allocated.size()) << where;
  for (std::size_t e = 0; e < allocated.size(); ++e) {
    EXPECT_GE(allocated[e], 0) << where << " edge " << e;
    EXPECT_LE(allocated[e], baseline[e]) << where << " edge " << e;
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
  for (const auto& c : detail::ProgramAccess::data(*p).contexts) {
    s.members.push_back(c.plan->members);
    s.weights.push_back(c.plan->weights);
    s.roots.emplace_back(c.plan->roots.begin(), c.plan->roots.end());
  }
  return s;
}

// Every field of a placement result.
inline void expect_same(const PlacementResult& a, const PlacementResult& b, const std::string& where) {
  EXPECT_EQ(a.placed, b.placed) << where;
  EXPECT_EQ(a.remaining, b.remaining) << where;
  EXPECT_EQ(a.delivery_delta, b.delivery_delta) << where;
  EXPECT_EQ(a.call_drops, b.call_drops) << where;
  EXPECT_EQ(a.call_source_loss, b.call_source_loss) << where;
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
  EXPECT_EQ(a.source_loss, b.source_loss) << where;
  EXPECT_EQ(a.drops, b.drops) << where;
  EXPECT_TRUE(a.by_source == b.by_source) << where;
  EXPECT_TRUE(a.by_tunnel == b.by_tunnel) << where;
  EXPECT_EQ(a.total_cost, b.total_cost) << where;
  EXPECT_EQ(a.cost_distribution, b.cost_distribution) << where;
  EXPECT_EQ(a.edge_flows, b.edge_flows) << where;
  EXPECT_EQ(a.residual, b.residual) << where;
  EXPECT_EQ(a.bottleneck, b.bottleneck) << where;
  EXPECT_EQ(a.reservation_edge_flows, b.reservation_edge_flows) << where;
  EXPECT_EQ(a.min_cut, b.min_cut) << where;
}

// One live flow as the tests compare it: its local id, its route cost, its bundle's edges and its
// program entry.
struct FlowSnapshot {
  std::int64_t id {0};
  Cost cost {0};
  std::vector<EdgeId> route;
  std::int32_t entry {-1};
  bool operator==(const FlowSnapshot&) const = default;
};

// A policy's epoch state: the call state (requested, offered, +inf, prepared, the placed
// demand), Global's cumulative curve, the installed program (its identity and its value) and T1's
// plan, CostOnly MaxFlow's support, a Pinned head's own losses, the by_tunnel rows and per cohort
// its flows.
struct Epoch {
  Counter128 requested;
  Counter128 offered;
  bool infinite {false};
  bool prepared {false};
  Amount placed {0};
  std::optional<std::vector<std::pair<Cost, Amount>>> curve;
  const ForwardingProgram* program_id {nullptr};
  ProgramSnapshot program;
  std::optional<TunnelPlan> plan;
  std::vector<EdgeId> support;
  std::vector<Counter128> head_lost;
  std::vector<TunnelFlow> by_tunnel;
  std::vector<std::vector<FlowSnapshot>> cohorts;
};

inline Epoch epoch(const FlowPolicy& p) {
  const auto& pe = detail::PolicyAccess::engine(p);
  const auto& ep = pe.epoch;
  Epoch e {ep.requested, ep.offered, ep.infinite, ep.prepared, p.placed_demand(), ep.curve, ep.program.get(),
           program_snapshot(ep.program.get()), ep.plan, ep.support, ep.head_lost, ep.by_tunnel, {}};
  for (const auto& c : ep.cohorts) {
    auto& out = e.cohorts.emplace_back();
    for (const auto& f : c) out.push_back({f.id, f.cost, f.bundle ? ascending_edges(*f.bundle) : std::vector<EdgeId> {}, f.entry});
  }
  return e;
}

inline void expect_same(const Epoch& a, const Epoch& b, const std::string& where) {
  EXPECT_EQ(a.requested, b.requested) << where;
  EXPECT_EQ(a.offered, b.offered) << where;
  EXPECT_EQ(a.infinite, b.infinite) << where;
  EXPECT_EQ(a.prepared, b.prepared) << where;
  EXPECT_EQ(a.placed, b.placed) << where;
  EXPECT_EQ(a.curve, b.curve) << where;
  EXPECT_TRUE(a.program == b.program) << where;
  EXPECT_TRUE(a.plan == b.plan) << where;
  EXPECT_EQ(a.support, b.support) << where;
  EXPECT_EQ(a.head_lost, b.head_lost) << where;
  EXPECT_TRUE(a.by_tunnel == b.by_tunnel) << where;
  EXPECT_TRUE(a.cohorts == b.cohorts) << where;
}

// The same state on the same policy: also the same program object (a failed call installs none).
inline void expect_unchanged(const Epoch& a, const Epoch& b, const std::string& where) {
  expect_same(a, b, where);
  EXPECT_EQ(a.program_id, b.program_id) << where;
}

} // namespace netgraph::core::test
