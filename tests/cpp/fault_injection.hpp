// The allocation-failure hook of the fault-injection executable (netgraph_core_fault_tests) and
// the state its sweeps compare: flow_ledger_fault_tests.cpp replaces the global operator new so that
// the chosen allocation throws std::bad_alloc; every test file of the executable sweeps its
// operations through fails_at(). The counter is shared by every thread (the evaluators' workers
// allocate too): the k-th allocation in the process fails, whichever thread makes it.
#pragma once

#include <gtest/gtest.h>

#include <atomic>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"

// Allocations left before the injected failure; negative: no injection.
extern std::atomic<long> g_countdown;

// Runs fn with the k-th allocation (from 0) failing; returns whether that failure happened
// (and expects that a reached failure propagated as std::bad_alloc).
bool fails_at(long k, const std::function<void()>& fn);

// As fails_at, counting only the allocations of threads other than the caller's (the
// evaluators' worker threads): the k-th worker allocation fails.
bool fails_in_workers_at(long k, const std::function<void()>& fn);

namespace netgraph::core::test {

// The ledger as a sweep compares it: the aggregate, the residual, every record field (FlowIds
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

inline LedgerState ledger_state(const FlowLedger& fg) {
  using detail::LedgerAccess;
  LedgerState s {fg.allocated(), fg.residual(), fg.records(), LedgerAccess::next_flow_id(fg),
                 LedgerAccess::bound_owners(fg), LedgerAccess::orphan_owners(fg)};
  for (auto& r : s.records) {
    EXPECT_EQ(r.flow_id.first, fg.uid());
    r.flow_id.first = 0;
  }
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

// The ledger invariant: the aggregate is the sum of the records, within [0, baseline].
inline void expect_balanced(const FlowLedger& fg, const std::string& where) {
  std::vector<Amount> sum(fg.allocated().size(), 0);
  for (const auto& r : fg.records()) {
    for (const auto& [e, x] : r.edges) sum[static_cast<std::size_t>(e)] += x;
  }
  EXPECT_EQ(sum, fg.allocated()) << where;
  for (Amount x : fg.residual()) EXPECT_GE(x, 0) << where;
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

} // namespace netgraph::core::test
