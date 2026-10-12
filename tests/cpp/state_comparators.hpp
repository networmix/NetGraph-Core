// What a call leaves, as the tests compare it: the ledger (LedgerState, every record field), its
// invariant (expect_balanced), a policy's epoch state (Epoch) and a placement result; the
// comparators' negative controls intercept their failures with test_utils.hpp's failures_of.
// Shared by both test executables; the fault-injection sweeps (fault_injection.hpp) compare
// through these.
#pragma once

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "test_utils.hpp"

namespace netgraph::core::test {

// The ledger as the tests compare it: the aggregate, the residual, every record field (FlowIds
// keep their local ids: a reference sequence runs on another ledger), the next FlowId and the
// bound and orphaned owners.
struct LedgerState {
  std::vector<Amount> allocated;
  std::vector<Amount> residual;
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
  LedgerState s {fg.allocated(), fg.residual(), fg.records(), LedgerAccess::next_flow_id(fg),
                 LedgerAccess::bound_owners(fg), LedgerAccess::orphan_owners(fg)};
  expect_uid_then_clear(s.records, fg.uid());
  return s;
}

inline bool same_record(const FlowRecord& a, const FlowRecord& b) {
  return a.flow_id == b.flow_id && a.owner == b.owner && a.sources == b.sources && a.destinations == b.destinations &&
         a.placement == b.placement && a.admission == b.admission && a.bundle && b.bundle && *a.bundle == *b.bundle &&
         a.route_cost == b.route_cost && a.offered == b.offered && a.delivered == b.delivered && a.edges == b.edges &&
         a.drops == b.drops;
}

inline void expect_same(const LedgerState& a, const LedgerState& b, const std::string& where) {
  EXPECT_EQ(a.allocated, b.allocated) << where;
  EXPECT_EQ(a.residual, b.residual) << where;
  EXPECT_EQ(a.next_id, b.next_id) << where;
  EXPECT_EQ(a.bound, b.bound) << where;
  EXPECT_EQ(a.orphans, b.orphans) << where;
  ASSERT_EQ(a.records.size(), b.records.size()) << where;
  for (std::size_t i = 0; i < a.records.size(); ++i) {
    EXPECT_TRUE(same_record(a.records[i], b.records[i])) << where << " record " << i;
  }
}

// The ledger invariant: the aggregate is the sum of the records, within [0, baseline] on every
// edge.
inline void expect_balanced(const FlowLedger& fg, const std::string& where) {
  const std::vector<Amount> allocated = fg.allocated();
  const auto baseline = detail::LedgerAccess::baseline(fg);
  std::vector<Amount> sum(allocated.size(), 0);
  for (const auto& r : fg.records()) {
    for (const auto& [e, x] : r.edges) sum[static_cast<std::size_t>(e)] += x;
  }
  EXPECT_EQ(sum, allocated) << where;
  ASSERT_EQ(baseline.size(), allocated.size()) << where;
  for (std::size_t e = 0; e < allocated.size(); ++e) {
    EXPECT_GE(allocated[e], 0) << where << " edge " << e;
    EXPECT_LE(allocated[e], baseline[e]) << where << " edge " << e;
  }
}

// Every field of a placement result.
inline void expect_same(const PlacementResult& a, const PlacementResult& b, const std::string& where) {
  EXPECT_EQ(a.placed, b.placed) << where;
  EXPECT_EQ(a.remaining, b.remaining) << where;
  EXPECT_EQ(a.termination, b.termination) << where;
  EXPECT_EQ(a.local_outcomes, b.local_outcomes) << where;
  EXPECT_EQ(a.routes_changed, b.routes_changed) << where;
  EXPECT_EQ(a.weights_changed, b.weights_changed) << where;
  EXPECT_EQ(a.drops, b.drops) << where;
  EXPECT_EQ(a.cost_distribution, b.cost_distribution) << where;
  EXPECT_EQ(a.edge_flows, b.edge_flows) << where;
  EXPECT_EQ(a.residual, b.residual) << where;
  EXPECT_EQ(a.bottleneck, b.bottleneck) << where;
  EXPECT_EQ(a.min_cut, b.min_cut) << where;
}

// One live flow as the tests compare it: its local id, its route cost and its plan's route,
// members and Weighted weights (a Global flow has no plan: its bundle's edges as route and
// members).
struct FlowSnapshot {
  std::int64_t id {0};
  Cost cost {0};
  std::vector<EdgeId> route;
  std::vector<EdgeId> members;
  std::vector<Amount> weights;
  bool operator==(const FlowSnapshot&) const = default;
};

// A policy's epoch state: the call state (requested, +inf, the placed demand), Global's
// cumulative curve and per roster source its flows.
struct Epoch {
  Counter128 requested;
  bool infinite {false};
  Amount placed {0};
  std::vector<std::pair<Cost, Amount>> curve;
  std::vector<std::vector<FlowSnapshot>> cohorts;
};

inline Epoch epoch(const FlowPolicy& p) {
  const auto& pe = detail::PolicyAccess::engine(p);
  Epoch e {pe.requested, pe.infinite, p.placed_demand(), pe.curve, {}};
  for (const auto& c : pe.cohorts) {
    auto& out = e.cohorts.emplace_back();
    for (const auto& f : c) {
      if (f.plan) {
        out.push_back({f.id, f.cost, f.plan->route, f.plan->members, f.plan->weights});
      } else {
        const auto ids = ascending_edges(*f.bundle);
        out.push_back({f.id, f.cost, ids, ids, {}});
      }
    }
  }
  return e;
}

inline void expect_same(const Epoch& a, const Epoch& b, const std::string& where) {
  EXPECT_EQ(a.requested, b.requested) << where;
  EXPECT_EQ(a.infinite, b.infinite) << where;
  EXPECT_EQ(a.placed, b.placed) << where;
  EXPECT_EQ(a.curve, b.curve) << where;
  EXPECT_TRUE(a.cohorts == b.cohorts) << where;
}

} // namespace netgraph::core::test
