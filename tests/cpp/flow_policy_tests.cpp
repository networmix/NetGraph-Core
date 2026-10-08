// The policy engine beyond prototype 2's fixtures (those replay in flow_policy_fixture_tests.cpp):
// configuration validation, the outputs, Global's empty roster, record bundles (sparse storage,
// sharing), the DagCache in building mode and the evaluators (capacity, batch_capacity,
// sensitivity) with workers (APPROACH 8 to 11).
#include <gtest/gtest.h>

#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/search.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

std::shared_ptr<const Graph> shared(Graph g) {
  return std::make_shared<const Graph>(std::move(g));
}

std::shared_ptr<const PredDAG> dag_of(const Graph& g, std::vector<EdgeId> edges) {
  return std::make_shared<const PredDAG>(detail::PredDAGBuilder::of_edges(g, edges));
}

PathSelection sel(EdgeTieBreak t = EdgeTieBreak::Deterministic, EdgePriority p = EdgePriority::Capacity,
                  bool multipath = true, bool multi_edge = true) {
  return PathSelection(t, p, multipath, multi_edge);
}

FlowPolicyConfig cfg(Routing r, FlowPlacement k, Admission a = Admission::Lossless, PathSelection s = sel()) {
  return FlowPolicyConfig {std::move(r), k, a, s};
}

FlowPolicyConfig global_cfg() { return cfg(Global {}, FlowPlacement::MaxFlow); }

std::string error_of(const std::function<void()>& fn) {
  try {
    fn();
  } catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
}

bool contains(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

std::vector<EdgeId> edges_of(const PredDAG& dag) {
  const auto e = detail::PredDAGBuilder::ascending(dag);
  return {e.begin(), e.end()};
}

std::vector<NodeId> v(std::initializer_list<NodeId> x) { return x; }

// Two disjoint routes 0 -> 1 -> 3 (caps a, cost 1) and 0 -> 2 -> 3 (caps b, cost c). B' ids:
// 0 = 0->1, 1 = 0->2, 2 = 1->3, 3 = 2->3.
std::shared_ptr<const Graph> diamond(Amount a, Amount b, Cost c = 1) {
  return shared(make_graph(4, {{0, 1, a, 1}, {1, 3, a, 1}, {0, 2, b, c}, {2, 3, b, c}}));
}

} // namespace

TEST(FlowPolicyConfig, EveryRejectedCombinationBeforeBinding) {
  auto g = diamond(5, 5);
  FlowLedger fg(g);
  auto rejects = [&](const FlowPolicyConfig& c, const std::string& part, std::vector<NodeId> s = {0},
                     std::vector<NodeId> d = {3}) {
    const auto msg = error_of([&] { FlowPolicy p(fg, s, d, c); });
    EXPECT_TRUE(contains(msg, part)) << "expected '" << part << "', got '" << msg << "'";
    EXPECT_TRUE(detail::LedgerAccess::bound_owners(fg).empty());
  };
  rejects(cfg(CostOnly {}, FlowPlacement::MaxFlow, Admission::Lossy), "Lossy");
  rejects(cfg(CostOnly {}, FlowPlacement::Weighted, Admission::Lossless, sel(EdgeTieBreak::Deterministic, EdgePriority::Headroom)),
          "Headroom");
  rejects(cfg(Global {}, FlowPlacement::EqualBalanced), "Global routing requires placement = MaxFlow");
  rejects(cfg(Global {}, FlowPlacement::EqualBalanced, Admission::Lossy), "Lossy");
  rejects(cfg(CapacityAware {}, FlowPlacement::EqualBalanced, Admission::Lossy), "Lossy");
  CapacityAware ca;
  ca.flow_count = 0;
  rejects(cfg(ca, FlowPlacement::MaxFlow), "flow_count");
  ca = {};
  ca.flow_count = 3;
  ca.max_flow_count = 2;
  rejects(cfg(ca, FlowPlacement::MaxFlow), "max_flow_count");
  ca = {};
  ca.min_flow_bandwidth = -1;
  rejects(cfg(ca, FlowPlacement::MaxFlow), "min_flow_bandwidth");
  ca = {};
  ca.max_path_cost = Absolute {-1};
  rejects(cfg(ca, FlowPlacement::MaxFlow), "Absolute");
  for (double f : {0.5, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
    ca = {};
    ca.max_path_cost = Stretch {f};
    rejects(cfg(ca, FlowPlacement::MaxFlow), "Stretch");
  }
  ca = {};
  ca.diminishing = DiminishingReturns {0, 0.5};
  rejects(cfg(ca, FlowPlacement::MaxFlow), "window");
  for (double f : {0.0, 1.0, std::numeric_limits<double>::quiet_NaN()}) {
    ca = {};
    ca.diminishing = DiminishingReturns {1, f};
    rejects(cfg(ca, FlowPlacement::MaxFlow), "fraction");
  }
  Pinned none;
  rejects(cfg(none, FlowPlacement::MaxFlow), "empty bundle list");
  Pinned dup;
  dup.bundles = {dag_of(*g, {0, 2}), dag_of(*g, {0, 2})};
  rejects(cfg(dup, FlowPlacement::MaxFlow), "duplicate bundle");
  Pinned malformed;
  malformed.bundles = {dag_of(*g, {0, 3})};
  rejects(cfg(malformed, FlowPlacement::MaxFlow), "malformed bundle");
  Pinned ok;
  ok.bundles = {dag_of(*g, {0, 2})};
  rejects(cfg(ok, FlowPlacement::MaxFlow), "single endpoints", {0, 1}, {3});
  auto other = diamond(5, 5);
  Pinned foreign;
  foreign.bundles = {dag_of(*other, {0, 2})};
  rejects(cfg(foreign, FlowPlacement::MaxFlow), "foreign PredDAG");
  rejects(cfg(CostOnly {}, FlowPlacement::MaxFlow), "both endpoint sets", {0, 3}, {3});
  rejects(cfg(CostOnly {}, FlowPlacement::MaxFlow), "empty", {}, {3});
  rejects(cfg(CostOnly {}, FlowPlacement::MaxFlow), "out of node range", {9}, {3});
  const auto el = other->eligibility({}, std::vector<EdgeId> {0});
  EXPECT_TRUE(contains(error_of([&] { FlowPolicy p(fg, v({0}), v({3}), cfg(CostOnly {}, FlowPlacement::MaxFlow), &el); }),
                       "foreign Eligibility"));
  DagCache foreign_cache(*other);
  EXPECT_TRUE(contains(error_of([&] {
                         FlowPolicy p(fg, v({0}), v({3}), cfg(CostOnly {}, FlowPlacement::MaxFlow), nullptr, &foreign_cache);
                       }),
                       "foreign DagCache"));
  EXPECT_TRUE(detail::LedgerAccess::bound_owners(fg).empty());
}

// Rule 13 at construction: PathSelection and FlowPolicyConfig reject a value outside its
// enumeration (here 99, reachable by a cast) with std::invalid_argument naming the field, in
// validate_config's order; a selection assigned out of its domain after its own construction is
// rejected by the configuration. Combinations are the consumers' (the next test).
TEST(FlowPolicyConfig, ConstructorsRejectEnumValuesOutsideTheirEnumeration) {
  EXPECT_EQ(error_of([] { (void)PathSelection(static_cast<EdgeTieBreak>(99), EdgePriority::Capacity); }),
            "tie_break: 99 is not a value of EdgeTieBreak");
  EXPECT_EQ(error_of([] { (void)PathSelection(EdgeTieBreak::Deterministic, static_cast<EdgePriority>(99)); }),
            "priority: 99 is not a value of EdgePriority");
  EXPECT_EQ(error_of([] { (void)PathSelection(static_cast<EdgeTieBreak>(98), static_cast<EdgePriority>(99)); }),
            "tie_break: 98 is not a value of EdgeTieBreak");
  const auto s = sel();
  EXPECT_EQ(error_of([&] {
              (void)FlowPolicyConfig {CostOnly {}, static_cast<FlowPlacement>(99), Admission::Lossless, s};
            }),
            "placement: 99 is not a value of FlowPlacement");
  EXPECT_EQ(error_of([&] {
              (void)FlowPolicyConfig {CostOnly {}, FlowPlacement::MaxFlow, static_cast<Admission>(99), s};
            }),
            "admission: 99 is not a value of Admission");
  auto tie = s;
  tie.tie_break = static_cast<EdgeTieBreak>(99);
  auto priority = s;
  priority.priority = static_cast<EdgePriority>(99);
  EXPECT_EQ(error_of([&] { (void)FlowPolicyConfig {CostOnly {}, FlowPlacement::MaxFlow, Admission::Lossless, tie}; }),
            "selection.tie_break: 99 is not a value of EdgeTieBreak");
  EXPECT_EQ(error_of([&] {
              (void)FlowPolicyConfig {CostOnly {}, FlowPlacement::MaxFlow, Admission::Lossless, priority};
            }),
            "selection.priority: 99 is not a value of EdgePriority");
  EXPECT_EQ(error_of([&] {
              (void)FlowPolicyConfig {Global {}, static_cast<FlowPlacement>(97), static_cast<Admission>(98), tie};
            }),
            "placement: 97 is not a value of FlowPlacement");
  EXPECT_EQ(error_of([&] {
              (void)FlowPolicyConfig {Global {}, FlowPlacement::MaxFlow, static_cast<Admission>(98), tie};
            }),
            "admission: 98 is not a value of Admission");
  // Every value inside its enumeration constructs, whatever the combination (Global with
  // Weighted and Lossy constructs; the policy rejects it).
  for (auto t : {EdgeTieBreak::Deterministic, EdgeTieBreak::PreferHigherPriority}) {
    for (auto p : {EdgePriority::Residual, EdgePriority::Capacity, EdgePriority::Headroom}) {
      EXPECT_NO_THROW((void)PathSelection(t, p, false, false));
      for (auto k : {FlowPlacement::MaxFlow, FlowPlacement::Weighted, FlowPlacement::EqualBalanced}) {
        for (auto a : {Admission::Lossless, Admission::Lossy}) {
          EXPECT_NO_THROW((void)(FlowPolicyConfig {Global {}, k, a, PathSelection(t, p)}));
        }
      }
    }
  }
}

// Rule 13 at every consuming entry: the fields stay assignable after construction, so a value
// outside its enumeration is a domain error named by its field, before any combination check,
// at every entry that takes the field; it never reaches a record.
TEST(FlowPolicyConfig, EnumValuesOutsideTheirEnumerationAreRejected) {
  auto g = diamond(5, 5);
  FlowLedger fg(g);
  const auto bad_placement = static_cast<FlowPlacement>(99);
  const auto bad_admission = static_cast<Admission>(99);
  const auto bad_tie = static_cast<EdgeTieBreak>(99);
  const auto bad_priority = static_cast<EdgePriority>(99);
  const auto bad_direction = static_cast<Direction>(99);
  auto tie_sel = sel();
  tie_sel.tie_break = bad_tie;
  auto priority_sel = sel();
  priority_sel.priority = bad_priority;
  Pinned pinned;
  pinned.bundles = {dag_of(*g, {0, 2})};
  const std::vector<Routing> modes = {CostOnly {}, CapacityAware {}, Global {}, pinned};
  struct Bad {
    std::function<FlowPolicyConfig(const Routing&)> make;
    const char* field;
  };
  // A valid configuration with one field assigned outside its enumeration.
  const auto assigned = [](const Routing& r, const std::function<void(FlowPolicyConfig&)>& set) {
    auto c = cfg(r, FlowPlacement::MaxFlow);
    set(c);
    return c;
  };
  const std::vector<Bad> bads = {
      {[&](const Routing& r) { return assigned(r, [&](FlowPolicyConfig& c) { c.placement = bad_placement; }); },
       "placement: 99 is not a value of FlowPlacement"},
      {[&](const Routing& r) { return assigned(r, [&](FlowPolicyConfig& c) { c.admission = bad_admission; }); },
       "admission: 99 is not a value of Admission"},
      {[&](const Routing& r) { return assigned(r, [&](FlowPolicyConfig& c) { c.selection = tie_sel; }); },
       "selection.tie_break: 99 is not a value of EdgeTieBreak"},
      {[&](const Routing& r) { return assigned(r, [&](FlowPolicyConfig& c) { c.selection = priority_sel; }); },
       "selection.priority: 99 is not a value of EdgePriority"},
  };
  for (const auto& mode : modes) {
    for (const auto& b : bads) {
      const auto c = b.make(mode);
      EXPECT_TRUE(contains(error_of([&] { FlowPolicy p(fg, v({0}), v({3}), c); }), b.field)) << b.field;
      EXPECT_TRUE(detail::LedgerAccess::bound_owners(fg).empty());
      EXPECT_TRUE(fg.records().empty());
      CapacityOptions co {c};
      EXPECT_TRUE(contains(error_of([&] { (void)capacity(*g, v({0}), v({3}), co); }), b.field)) << b.field;
      EXPECT_TRUE(contains(error_of([&] { (void)sensitivity(*g, v({0}), v({3}), SensitivityOptions {c}, 1); }), b.field))
          << b.field;
      if (!std::holds_alternative<Pinned>(mode)) {
        const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{v({0}), v({3})}};
        EXPECT_TRUE(contains(error_of([&] { (void)batch_capacity(*g, pairs, co, {}, 1); }), b.field)) << b.field;
      }
    }
  }
  // The kernel check, the search and the cache's key.
  EXPECT_TRUE(contains(error_of([&] { detail::check_kernel(bad_placement, Admission::Lossless, EdgePriority::Capacity); }),
                       "placement: 99"));
  EXPECT_TRUE(contains(error_of([&] { detail::check_kernel(FlowPlacement::MaxFlow, bad_admission, EdgePriority::Capacity); }),
                       "admission: 99"));
  EXPECT_TRUE(contains(error_of([&] { detail::check_kernel(FlowPlacement::MaxFlow, Admission::Lossless, bad_priority); }),
                       "priority: 99"));
  {
    SpfOptions o(Direction::Forward, sel());
    o.direction = bad_direction;
    EXPECT_TRUE(contains(error_of([&] { (void)spf(*g, v({0}), o); }), "direction: 99 is not a value of Direction"));
    EXPECT_TRUE(contains(error_of([&] { SpfOptions bad(bad_direction, sel()); }), "direction: 99 is not a value of Direction"));
    SpfOptions t(Direction::Forward, tie_sel);
    EXPECT_TRUE(contains(error_of([&] { (void)spf(*g, v({0}), t); }), "selection.tie_break: 99"));
    SpfOptions p(Direction::Forward, priority_sel);
    EXPECT_TRUE(contains(error_of([&] { (void)spf(*g, v({0}), p); }), "selection.priority: 99"));
  }
  DagCache cache(*g);
  EXPECT_TRUE(contains(error_of([&] { (void)cache.lookup(*g, v({3}), bad_direction, nullptr); }), "direction: 99"));
  const auto st = cache.stats();
  EXPECT_EQ(st.hits + st.misses + st.entries, 0);
}

TEST(FlowPolicyRequests, ValidationPrecedesTheZeroShortcut) {
  auto g = diamond(5, 5);
  FlowLedger fg(g);
  FlowPolicy p(fg, v({0}), v({3}), cfg(CostOnly {}, FlowPlacement::EqualBalanced));
  EXPECT_TRUE(contains(error_of([&] { (void)p.place(-1); }), "request"));
  EXPECT_TRUE(contains(error_of([&] { (void)p.place(kDomainLimit); }), "request"));
  EXPECT_EQ(p.place(std::nullopt).placed, 10);
  EXPECT_TRUE(contains(error_of([&] { (void)p.place(1); }), "after +inf"));
  EXPECT_TRUE(contains(error_of([&] { (void)p.place(0); }), "after +inf"));
  p.remove();
  EXPECT_EQ(p.place(1).placed, 1);
}

TEST(FlowPolicyOutputs, OwnAllocationLedgerResidualBottleneckAndGlobalCut) {
  auto g = diamond(5, 3);
  FlowLedger fg(g);
  // A third party holds one unit of the first route.
  Pinned third;
  third.bundles = {dag_of(*g, {0, 2})};
  FlowPolicy other(fg, v({0}), v({3}), cfg(third, FlowPlacement::MaxFlow));
  ASSERT_EQ(other.place(1).placed, 1);
  FlowPolicy p(fg, v({0}), v({3}), cfg(CostOnly {}, FlowPlacement::EqualBalanced));
  const OutputOptions all {true, true, true};
  auto r = p.place(4, all);
  ASSERT_TRUE(r.edge_flows && r.residual && r.bottleneck);
  EXPECT_FALSE(r.min_cut.has_value());  // not Global
  EXPECT_FALSE(r.cost_distribution.has_value());
  std::vector<Amount> own(4, 0);
  for (const auto& rec : p.records()) {
    for (const auto& [e, x] : rec.edges) own[static_cast<std::size_t>(e)] += x;
  }
  EXPECT_EQ(*r.edge_flows, own);
  EXPECT_EQ(*r.residual, fg.residual());
  EXPECT_EQ(*r.bottleneck, fg.bottleneck(v({0})));
  // Nothing requested: absent.
  r = p.place(1);
  EXPECT_FALSE(r.edge_flows || r.residual || r.bottleneck || r.min_cut);
  // Global: the curve always, min_cut with bottleneck when the solve is maximal (at +inf; none
  // under a cap of 2 below the maximum 8, which leaves a residual path).
  FlowLedger fg2(g);
  FlowPolicy gl(fg2, v({0}), v({3}), global_cfg());
  r = gl.place(2, OutputOptions {false, false, true});
  ASSERT_TRUE(r.cost_distribution.has_value());
  EXPECT_EQ(*r.cost_distribution, curve({{2, 2}}));
  EXPECT_LT(gl.placed_demand(), 8);
  EXPECT_FALSE(r.min_cut.has_value());
  r = gl.place(std::nullopt, OutputOptions {false, false, true});
  EXPECT_EQ(r.placed, 6);
  ASSERT_TRUE(r.min_cut.has_value());
  EXPECT_EQ(r.min_cut->size(), 2U);
}

TEST(FlowPolicyGlobal, EmptyRosterReportsAnEmptyCurveAndAnEmptyCut) {
  // No source reaches the destination: the maximum flow is 0 and maximal (APPROACH 2.2).
  auto g = shared(make_graph(3, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  FlowPolicy p(fg, v({0}), v({2}), global_cfg());
  auto r = p.place(std::nullopt, OutputOptions {false, false, true});
  EXPECT_EQ(r.termination, Termination::NoRoute);
  ASSERT_TRUE(r.cost_distribution.has_value());
  EXPECT_TRUE(r.cost_distribution->empty());
  ASSERT_TRUE(r.min_cut.has_value());
  EXPECT_TRUE(r.min_cut->empty());
  EXPECT_TRUE(p.records().empty());
  // The zero request is a no-op: no solve, no cut; the cumulative curve is empty.
  FlowPolicy q(fg, v({1}), v({2}), global_cfg());
  r = q.place(0, OutputOptions {false, false, true});
  EXPECT_FALSE(r.min_cut.has_value());
  ASSERT_TRUE(r.cost_distribution.has_value());
  EXPECT_TRUE(r.cost_distribution->empty());
}

// A Global zero call reports the cumulative curve (9.1: its amounts sum to placed_demand();
// empty before the first placement and after remove()) and changes nothing; it solves nothing,
// so it claims no cut (review-1c F5).
TEST(FlowPolicyGlobal, ZeroCallReportsTheCumulativeCurve) {
  using Curve = std::vector<std::pair<Cost, Amount>>;
  auto g = shared(make_graph(3, {{0, 1, 3, 7}, {0, 2, 4, 1}, {2, 1, 4, 1}}));
  FlowLedger fg(g);
  FlowPolicy p(fg, v({0}), v({1}), global_cfg());
  auto z = p.place(0);
  ASSERT_TRUE(z.cost_distribution.has_value());
  EXPECT_TRUE(z.cost_distribution->empty());
  auto r = p.place(5);
  EXPECT_EQ(*r.cost_distribution, (Curve {{2, 4}, {7, 1}}));
  const auto allocated = fg.allocated();
  const auto requested = detail::PolicyAccess::engine(p).requested;
  z = p.place(0, OutputOptions {true, true, true});
  EXPECT_EQ(z.placed, 0);
  EXPECT_EQ(*z.cost_distribution, *r.cost_distribution);
  EXPECT_FALSE(z.min_cut.has_value());
  EXPECT_EQ(fg.allocated(), allocated);
  EXPECT_EQ(detail::PolicyAccess::engine(p).requested, requested);
  r = p.place(1);  // cumulative 6: the solve's curve on V = W + own
  EXPECT_EQ(*r.cost_distribution, (Curve {{2, 4}, {7, 2}}));
  z = p.place(0);
  EXPECT_EQ(*z.cost_distribution, *r.cost_distribution);
  Amount sum = 0;
  for (const auto& [c, x] : *z.cost_distribution) sum += x;
  EXPECT_EQ(sum, p.placed_demand());
  p.remove();
  z = p.place(0);
  ASSERT_TRUE(z.cost_distribution.has_value());
  EXPECT_TRUE(z.cost_distribution->empty());
}

TEST(FlowPolicyRecords, BundlesAreSharedOrSparse) {
  auto g = diamond(5, 5);
  FlowLedger fg(g);
  const auto whole = dag_of(*g, {0, 1, 2, 3});
  Pinned pin;
  pin.bundles = {whole};
  FlowPolicy p(fg, v({0}), v({3}), cfg(pin, FlowPlacement::EqualBalanced));
  (void)p.place(2);
  // An unchanged supplied bundle is shared by the record (no copy).
  EXPECT_EQ(p.records().front().bundle.get(), whole.get());
  // A pruned bundle is a new DAG of the usable route.
  const auto el = g->eligibility({}, std::vector<EdgeId> {2});
  FlowLedger fg2(g);
  FlowPolicy q(fg2, v({0}), v({3}), cfg(pin, FlowPlacement::EqualBalanced), &el);
  (void)q.place(2);
  const auto rec = q.records().front();
  EXPECT_NE(rec.bundle.get(), whole.get());
  EXPECT_EQ(edges_of(*rec.bundle), (std::vector<EdgeId> {1, 3}));
  EXPECT_EQ(q.usable_bundles(), 1);
  // CostOnly: a source whose route is the whole cached support shares it.
  DagCache cache(*g);
  FlowLedger fg3(g);
  FlowPolicy c(fg3, v({0}), v({3}), cfg(CostOnly {}, FlowPlacement::EqualBalanced), nullptr, &cache);
  (void)c.place(2);
  const auto support = cache.lookup(*g, v({3}), Direction::Reverse, nullptr);
  EXPECT_EQ(c.records().front().bundle.get(), support.get());
}

TEST(PredDAGStorage, ABundleCostsItsEntriesNotTheGraph) {
  const std::int32_t n = 1'000'000;
  auto g = shared(make_graph(n, {{n - 3, n - 2, 4, 1}, {n - 2, n - 1, 4, 1}}));
  const auto dag = dag_of(*g, {0, 1});
  EXPECT_LT(detail::PredDAGBuilder::payload_bytes(*dag), 128U);  // two entries, not N + 1 offsets
  const auto off = dag->parent_offsets();  // the from_csr shape, on demand
  ASSERT_EQ(off.size(), static_cast<std::size_t>(n) + 1);
  EXPECT_EQ(off[static_cast<std::size_t>(n) - 2], 0);
  EXPECT_EQ(off[static_cast<std::size_t>(n) - 1], 1);
  EXPECT_EQ(off[static_cast<std::size_t>(n)], 2);
  EXPECT_TRUE(PredDAG::from_csr(*g, off, dag->parents(), dag->via_edges()) == *dag);
  FlowLedger fg(g);
  Pinned pin;
  pin.bundles = {dag};
  FlowPolicy p(fg, v({n - 3}), v({n - 1}), cfg(pin, FlowPlacement::MaxFlow));
  EXPECT_EQ(p.place(3).placed, 3);
  EXPECT_LT(detail::PredDAGBuilder::payload_bytes(*p.records().front().bundle), 128U);
}

TEST(DagCache, KeyStatisticsAndTheBuildingTransitions) {
  auto g = diamond(5, 5);
  DagCache cache(*g);
  EXPECT_EQ(cache.stats().entries, 0);
  // Two CostOnly policies differing only in selection share one entry (selection is a
  // projection at creation, never part of the key).
  FlowLedger fg(g);
  FlowPolicy a(fg, v({0}), v({3}), cfg(CostOnly {}, FlowPlacement::EqualBalanced), nullptr, &cache);
  FlowPolicy b(fg, v({0}), v({3}),
               cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless,
                   sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, false, false)),
               nullptr, &cache);
  auto s = cache.stats();
  EXPECT_EQ(s.misses, 1);
  EXPECT_EQ(s.hits, 1);
  EXPECT_EQ(s.entries, 1);
  EXPECT_GT(s.bytes, 0);
  // The key: roots as a set, the direction, the Eligibility uid.
  const auto x = cache.lookup(*g, v({3, 3}), Direction::Reverse, nullptr);
  EXPECT_EQ(cache.stats().hits, 2);
  (void)cache.lookup(*g, v({3}), Direction::Forward, nullptr);
  const auto el = g->eligibility({}, std::vector<EdgeId> {0});
  const auto y = cache.lookup(*g, v({3}), Direction::Reverse, &el);
  s = cache.stats();
  EXPECT_EQ(s.misses, 3);
  EXPECT_EQ(s.entries, 3);
  EXPECT_EQ(edges_of(*y), (std::vector<EdgeId> {1, 2, 3}));  // the static support under el (no 0 -> 1)
  EXPECT_EQ(edges_of(*x), (std::vector<EdgeId> {0, 1, 2, 3}));
  // Invalid queries throw and count nothing.
  auto other = diamond(5, 5);
  EXPECT_THROW((void)cache.lookup(*other, v({3}), Direction::Reverse, nullptr), std::invalid_argument);
  const auto foreign = other->eligibility({}, std::vector<EdgeId> {0});
  EXPECT_THROW((void)cache.lookup(*g, v({3}), Direction::Reverse, &foreign), std::invalid_argument);
  EXPECT_THROW((void)cache.lookup(*g, v({9}), Direction::Reverse, nullptr), std::invalid_argument);
  EXPECT_THROW((void)cache.lookup(*g, std::vector<NodeId> {}, Direction::Reverse, nullptr), std::invalid_argument);
  EXPECT_EQ(cache.stats().misses, 3);
  EXPECT_EQ(cache.stats().hits, 2);
  // clear(): empties, counters persist, held values stay valid and alive.
  const auto held = x;
  EXPECT_EQ(cache.stats().live_values, 3);
  cache.clear();
  s = cache.stats();
  EXPECT_EQ(s.entries, 0);
  EXPECT_EQ(s.bytes, 0);
  EXPECT_EQ(s.hits, 2);
  EXPECT_EQ(s.misses, 3);
  EXPECT_EQ(s.live_values, 2);  // x (= held) and y; the Forward entry died with the clear
  EXPECT_EQ(edges_of(*held), (std::vector<EdgeId> {0, 1, 2, 3}));
  // After clear() a lookup misses and inserts again.
  (void)cache.lookup(*g, v({3}), Direction::Reverse, nullptr);
  EXPECT_EQ(cache.stats().misses, 4);
}

TEST(DagCache, ValuesOutliveTheCache) {
  auto g = diamond(5, 5);
  std::shared_ptr<const PredDAG> v1;
  {
    DagCache cache(*g);
    v1 = cache.lookup(*g, v({3}), Direction::Reverse, nullptr);
  }
  EXPECT_EQ(edges_of(*v1), (std::vector<EdgeId> {0, 1, 2, 3}));
}

TEST(Evaluators, CapacityIsPlaceInfOnAFreshLedger) {
  auto g = diamond(5, 3, 2);
  for (const auto& c : {global_cfg(), cfg(CostOnly {}, FlowPlacement::EqualBalanced),
                        cfg(CapacityAware {.flow_count = 2, .max_flow_count = 4}, FlowPlacement::MaxFlow)}) {
    const auto r = capacity(*g, v({0}), v({3}), CapacityOptions {c, nullptr, nullptr, OutputOptions {true, true, true}});
    FlowLedger fg(g);
    FlowPolicy p(fg, v({0}), v({3}), c);
    const auto want = p.place(std::nullopt, OutputOptions {true, true, true});
    EXPECT_EQ(r.placed, want.placed);
    EXPECT_EQ(r.termination, want.termination);
    EXPECT_EQ(r.edge_flows, want.edge_flows);
    EXPECT_EQ(r.residual, want.residual);
    EXPECT_EQ(r.bottleneck, want.bottleneck);
    EXPECT_EQ(r.min_cut, want.min_cut);
    EXPECT_EQ(r.cost_distribution, want.cost_distribution);
  }
  // The exact maximum flow: Global's curve (cheapest walks first).
  const auto exact = capacity(*g, v({0}), v({3}), CapacityOptions {global_cfg()});
  EXPECT_EQ(exact.placed, 8);
  EXPECT_EQ(*exact.cost_distribution, curve({{2, 5}, {4, 3}}));
}

TEST(Evaluators, BatchOrderIntersectionsExceptionsAndWorkers) {
  auto g = diamond(5, 3);
  const auto common = g->eligibility({}, std::vector<EdgeId> {3});       // route 2's last edge
  const auto first_route = g->eligibility({}, std::vector<EdgeId> {0});  // route 1's first edge
  std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {{0}, {3}}, {{0}, {3}}, {{1}, {3}}, {{0}, {1, 2}}};
  const std::vector<const Eligibility*> per = {nullptr, &first_route, nullptr, nullptr};
  const CapacityOptions opts {global_cfg(), &common};
  const auto r1 = batch_capacity(*g, pairs, opts, per, 1);
  ASSERT_EQ(r1.size(), 4U);
  EXPECT_EQ(r1[0].placed, 5);  // common only: route 1
  EXPECT_EQ(r1[1].placed, 0);  // common and the pair's: both routes cut
  EXPECT_EQ(r1[2].placed, 5);
  EXPECT_EQ(r1[3].placed, 8);  // destination set {1, 2}: both first edges (edge 3 excluded is past 2)
  for (std::int32_t w : {2, 3, 8}) {
    const auto rw = batch_capacity(*g, pairs, opts, per, w);
    ASSERT_EQ(rw.size(), 4U);
    for (std::size_t i = 0; i < 4; ++i) EXPECT_EQ(rw[i].placed, r1[i].placed) << "workers " << w << " pair " << i;
  }
  EXPECT_TRUE(batch_capacity(*g, {}, opts, {}, 4).empty());
  // Validation before any worker: endpoints, lengths, foreign values, Pinned, workers, the
  // building cache's one writer.
  auto bad = pairs;
  bad[2].second = {9};
  EXPECT_THROW((void)batch_capacity(*g, bad, opts, per, 4), std::invalid_argument);
  bad = pairs;
  bad[1].first = {3};
  EXPECT_TRUE(contains(error_of([&] { (void)batch_capacity(*g, bad, opts, per, 4); }), "both endpoint sets"));
  EXPECT_THROW((void)batch_capacity(*g, pairs, opts, std::vector<const Eligibility*> {nullptr}, 1), std::invalid_argument);
  auto other = diamond(5, 3);
  const auto foreign = other->eligibility({}, std::vector<EdgeId> {0});
  const std::vector<const Eligibility*> per_foreign = {nullptr, &foreign, nullptr, nullptr};
  EXPECT_TRUE(contains(error_of([&] { (void)batch_capacity(*g, pairs, opts, per_foreign, 1); }), "foreign"));
  Pinned pin;
  pin.bundles = {dag_of(*g, {0, 2})};
  EXPECT_TRUE(contains(error_of([&] { (void)batch_capacity(*g, pairs, CapacityOptions {cfg(pin, FlowPlacement::MaxFlow)}, {}, 1); }),
                       "Pinned"));
  EXPECT_TRUE(contains(error_of([&] { (void)batch_capacity(*g, pairs, opts, per, 0); }), "workers"));
  DagCache cache(*g);
  const CapacityOptions with_cache {cfg(CostOnly {}, FlowPlacement::EqualBalanced), nullptr, &cache};
  EXPECT_TRUE(contains(error_of([&] { (void)batch_capacity(*g, pairs, with_cache, {}, 2); }), "workers == 1"));
  const auto cached = batch_capacity(*g, pairs, with_cache, {}, 1);
  EXPECT_EQ(cache.stats().misses, 2);  // destination sets {3} and {1, 2}
  EXPECT_EQ(cache.stats().hits, 2);
  EXPECT_EQ(cached[0].placed, 7);  // CostOnly EB over routes of 5 and 3: [4, 3] (odd unit to the lower id)
}

TEST(Evaluators, BatchValidatesInPairOrderBeforeAnyWorker) {
  // Two invalid pairs (1 and 3): the batch raises pair 1's error before any worker starts.
  auto g = diamond(5, 3);
  std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {{0}, {3}}, {{0}, {7}}, {{0}, {3}}, {{0}, {8}}};
  const auto msg = error_of([&] { (void)batch_capacity(*g, pairs, CapacityOptions {global_cfg()}, {}, 4); });
  EXPECT_TRUE(contains(msg, "node 7")) << msg;
}

// The sensitivity shortcuts (Global from the baseline's maximum flow, CostOnly without a
// re-evaluation of the edges outside the baseline's routes) give the losses of the definition,
// the baseline's placed amount minus capacity() with the edge excluded, on random graphs with
// parallel edges, self-loops, dispersed costs, zero capacities, endpoint sets, eligibility
// masks and every CostOnly kernel, admission and selection.
TEST(Evaluators, SensitivityEqualsTheDefinitionOnRandomGraphs) {
  std::mt19937 rng(20261005);
  auto pick = [&](int n) { return static_cast<int>(rng() % static_cast<unsigned>(n)); };
  std::vector<FlowPolicyConfig> configs = {global_cfg()};
  for (auto k : {FlowPlacement::MaxFlow, FlowPlacement::EqualBalanced, FlowPlacement::Weighted}) {
    for (bool multipath : {true, false}) {
      for (bool multi_edge : {true, false}) {
        for (auto t : {EdgeTieBreak::Deterministic, EdgeTieBreak::PreferHigherPriority}) {
          configs.push_back(cfg(CostOnly {}, k, Admission::Lossless, sel(t, EdgePriority::Capacity, multipath, multi_edge)));
        }
      }
    }
  }
  configs.push_back(cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossy));
  configs.push_back(cfg(CostOnly {}, FlowPlacement::Weighted, Admission::Lossy, sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, false, false)));
  int entries = 0;
  int global_losses = 0;
  int policy_losses = 0;
  int policy_zero = 0;
  for (int trial = 0; trial < 60; ++trial) {
    const int n = 4 + pick(7);
    std::vector<InEdge> es;
    for (int k = 0, m = n + pick(3 * n); k < m; ++k) {
      es.push_back({pick(n), pick(n), pick(5) == 0 ? 0 : 1 + pick(9), static_cast<Cost>(1 + pick(4))});
    }
    auto g = shared(make_graph(n, es));
    std::vector<NodeId> s = {0};
    std::vector<NodeId> t = {n - 1};
    if (pick(2) == 0) s.push_back(1);
    if (pick(2) == 0) t.push_back(n - 2);
    std::vector<NodeId> xn;
    std::vector<EdgeId> xe;
    if (pick(3) == 0) xn.push_back(static_cast<NodeId>(2 + pick(n - 4 > 0 ? n - 4 : 1)));
    if (pick(2) == 0) xe.push_back(static_cast<EdgeId>(pick(static_cast<int>(es.size()))));
    const auto mask = g->eligibility(xn, xe);
    const Eligibility* elig = pick(2) == 0 ? &mask : nullptr;
    for (const auto& c : configs) {
      const bool global = std::holds_alternative<Global>(c.routing);
      const auto r = sensitivity(*g, s, t, SensitivityOptions {c, elig}, 1 + pick(3));
      for (const auto& e : r.entries) {
        ASSERT_EQ(e.status, SensitivityStatus::Evaluated);
        const EdgeId one[1] = {e.edge};
        auto deleted = g->eligibility({}, one);
        if (elig != nullptr) deleted = Eligibility::intersect(*elig, deleted);
        const Amount want = r.baseline.placed - capacity(*g, s, t, CapacityOptions {c, &deleted}).placed;
        ASSERT_EQ(*e.loss, want) << "trial " << trial << " edge " << e.edge << (global ? " Global" : " CostOnly");
        ++entries;
        if (want != 0) ++(global ? global_losses : policy_losses);
        if (!global && want == 0) ++policy_zero;
      }
    }
  }
  // Not vacuous: many entries, losses under both models and zero-loss CostOnly deletions.
  EXPECT_GT(entries, 3000);
  EXPECT_GT(global_losses, 100);
  EXPECT_GT(policy_losses, 500);
  EXPECT_GT(policy_zero, 500);
}

TEST(Evaluators, SensitivityCandidatesSignedLossesAndUnevaluable) {
  // Global: candidates are the baseline's positive-flow edges; an unsaturated critical edge
  // appears (rule 1): 0 -> 1 of capacity 9 carries 4 (its tail 1 -> 2 holds 4) and its deletion
  // loses 4; the saturated-only candidate set of 0.12 missed it.
  auto g = shared(make_graph(3, {{0, 1, 9, 1}, {1, 2, 4, 1}, {0, 2, 2, 5}}));
  auto s = sensitivity(*g, v({0}), v({2}), SensitivityOptions {global_cfg()}, 1);
  EXPECT_EQ(s.baseline.placed, 6);
  EXPECT_FALSE(s.baseline.edge_flows.has_value());  // requested internally, not returned
  std::map<EdgeId, Amount> losses;
  for (const auto& e : s.entries) {
    EXPECT_EQ(e.status, SensitivityStatus::Evaluated);
    losses[e.edge] = *e.loss;
  }
  const auto inv = id_of_input(*g);
  EXPECT_EQ(losses.at(inv[0]), 4);  // unsaturated (4 of 9) and critical
  EXPECT_EQ(losses.at(inv[1]), 4);
  EXPECT_EQ(losses.at(inv[2]), 2);
  EXPECT_EQ(losses.size(), 3U);
  // A policy model: every allowed edge of >= 1 unit; a loss can be negative (EqualBalanced over
  // parallel links of 1 and 10 delivers 2; without the thin link it delivers 10).
  auto p = shared(make_graph(2, {{0, 1, 1, 1}, {0, 1, 10, 1}, {0, 1, 0, 1}}));
  const auto opts = SensitivityOptions {cfg(CostOnly {}, FlowPlacement::EqualBalanced)};
  s = sensitivity(*p, v({0}), v({1}), opts, 1);
  EXPECT_EQ(s.baseline.placed, 2);
  ASSERT_EQ(s.entries.size(), 2U);  // the zero-capacity link is absent (loss zero)
  const auto pin = id_of_input(*p);
  for (const auto& e : s.entries) {
    if (e.edge == pin[0]) EXPECT_EQ(*e.loss, -8);
    if (e.edge == pin[1]) EXPECT_EQ(*e.loss, 1);
  }
  // Workers give the same entries; an excluded edge is no candidate.
  const auto s4 = sensitivity(*p, v({0}), v({1}), opts, 4);
  ASSERT_EQ(s4.entries.size(), s.entries.size());
  for (std::size_t i = 0; i < s.entries.size(); ++i) EXPECT_EQ(s4.entries[i].loss, s.entries[i].loss);
  const auto el = p->eligibility({}, std::vector<EdgeId> {pin[0]});
  const auto se = sensitivity(*p, v({0}), v({1}), SensitivityOptions {opts.cfg, &el}, 1);
  ASSERT_EQ(se.entries.size(), 1U);
  EXPECT_EQ(se.entries[0].edge, pin[1]);
  // Invalid inputs propagate (before any worker).
  EXPECT_THROW((void)sensitivity(*p, v({0}), v({0}), opts, 1), std::invalid_argument);
  DagCache cache(*p);
  EXPECT_THROW((void)sensitivity(*p, v({0}), v({1}), SensitivityOptions {opts.cfg, nullptr, &cache}, 2),
               std::invalid_argument);
  EXPECT_EQ(sensitivity(*p, v({0}), v({1}), SensitivityOptions {opts.cfg, nullptr, &cache}, 1).entries.size(), 2U);
}

// A deletion raising a model exception is unevaluable with its message, not raised (rule 1). No
// valid configuration raises under a deletion, so the private hook injects one.
TEST(Evaluators, SensitivityReportsAModelExceptionAsUnevaluable) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  auto p = shared(make_graph(2, {{0, 1, 1, 1}, {0, 1, 10, 1}, {0, 1, 0, 1}}));
  const auto opts = SensitivityOptions {cfg(CostOnly {}, FlowPlacement::EqualBalanced)};
  const auto pin = id_of_input(*p);
  detail::sensitivity_fail_edge().store(pin[1]);
  const auto su = sensitivity(*p, v({0}), v({1}), opts, 2);
  detail::sensitivity_fail_edge().store(-1);
  ASSERT_EQ(su.entries.size(), 2U);
  for (const auto& e : su.entries) {
    if (e.edge == pin[1]) {
      EXPECT_EQ(e.status, SensitivityStatus::Unevaluable);
      EXPECT_FALSE(e.loss.has_value());
      EXPECT_TRUE(contains(e.reason, "injected model exception"));
    } else {
      EXPECT_EQ(e.status, SensitivityStatus::Evaluated);
    }
  }
#else
  GTEST_SKIP() << "needs the private hooks (NETGRAPH_CORE_PRIVATE_HOOKS)";
#endif
}
