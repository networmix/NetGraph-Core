// The policy engine beyond the reference's fixtures (those replay in flow_policy_fixture_tests.cpp):
// configuration validation, the configuration a policy keeps and the selection fields Global and
// Pinned read, request validation (a rejected request changes nothing) and the per-source round
// robin of the cumulative request, the outputs (the bottleneck under the policy's Eligibility),
// Global's empty roster, its roster without flow and its zero call's cumulative curve, record
// bundles (sparse storage, sharing, pruning), records in FlowId order, one flow to several nearest
// destinations, the reduction of the sources' outcomes, CapacityAware's diminishing guard, path-cost
// ceiling (exact at large costs and factors), stretch baseline and bandwidth floor, and the
// evaluators (capacity, batch_capacity, sensitivity) with workers and their default from the
// environment. The DagCache's own modes and key are tested in dag_cache_tests.cpp.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/detail/evaluator_checks.hpp"
#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/traversal.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/search.hpp"
#include "state_comparators.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

PathSelection sel(EdgeTieBreak t = EdgeTieBreak::Deterministic, EdgePriority p = EdgePriority::Capacity,
                  bool multipath = true, bool multi_edge = true) {
  return PathSelection(t, p, multipath, multi_edge);
}

FlowPolicyConfig cfg(Routing r, FlowPlacement k, Admission a = Admission::Lossless, PathSelection s = sel()) {
  return FlowPolicyConfig {std::move(r), k, a, s};
}

// Two disjoint routes 0 -> 1 -> 3 (caps a, cost 1) and 0 -> 2 -> 3 (caps b, cost c). EdgeIds:
// 0 = 0->1, 1 = 0->2, 2 = 1->3, 3 = 2->3.
std::shared_ptr<const Graph> diamond(Amount a, Amount b, Cost c = 1) {
  return shared(make_graph(4, {{0, 1, a, 1}, {1, 3, a, 1}, {0, 2, b, c}, {2, 3, b, c}}));
}

} // namespace

TEST(FlowPolicyConfig, EveryRejectedCombinationBeforeBinding) {
  auto g = diamond(5, 5);
  FlowLedger fg(g);
  // Each combination is refused with its whole message, and no policy is bound afterwards.
  auto rejects = [&](const FlowPolicyConfig& c, const std::string& message, std::vector<NodeId> s = {0},
                     std::vector<NodeId> d = {3}) {
    EXPECT_EQ(error_of([&] { FlowPolicy p(fg, s, d, c); }), message);
    EXPECT_TRUE(detail::LedgerAccess::bound_owners(fg).empty());
  };
  const std::string lossy_hop_by_hop =
      "Lossy admission is rejected with Global and CapacityAware routing (hop-by-hop forwarding without admission "
      "control neither grows nor reroutes)";
  rejects(cfg(CostOnly {}, FlowPlacement::MaxFlow, Admission::Lossy),
          "Lossy admission is rejected with the MaxFlow kernel: a controller never drops (EqualBalanced or Weighted only)");
  rejects(cfg(CostOnly {}, FlowPlacement::Weighted, Admission::Lossless, sel(EdgeTieBreak::Deterministic, EdgePriority::Headroom)),
          "the Weighted kernel is rejected with priority Headroom (weights are units)");
  rejects(cfg(Global {}, FlowPlacement::EqualBalanced), "Global routing requires placement = MaxFlow");
  rejects(cfg(Global {}, FlowPlacement::EqualBalanced, Admission::Lossy), lossy_hop_by_hop);
  rejects(cfg(CapacityAware {}, FlowPlacement::EqualBalanced, Admission::Lossy), lossy_hop_by_hop);
  CapacityAware ca;
  ca.flow_count = 0;
  rejects(cfg(ca, FlowPlacement::MaxFlow), "CapacityAware flow_count must be >= 1 (got 0)");
  ca = {};
  ca.flow_count = 3;
  ca.max_flow_count = 2;
  rejects(cfg(ca, FlowPlacement::MaxFlow), "CapacityAware max_flow_count (2) must be >= flow_count (3)");
  ca = {};
  ca.min_flow_bandwidth = -1;
  rejects(cfg(ca, FlowPlacement::MaxFlow), "CapacityAware min_flow_bandwidth must be in [0, 2^62) units (got -1)");
  // The upper bound: 2^62 - 1 units is a valid floor, 2^62 is not.
  ca.min_flow_bandwidth = kDomainLimit;
  rejects(cfg(ca, FlowPlacement::MaxFlow),
          "CapacityAware min_flow_bandwidth must be in [0, 2^62) units (got 4611686018427387904)");
  ca.min_flow_bandwidth = kDomainLimit - 1;
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({3}), cfg(ca, FlowPlacement::MaxFlow)); }), "");
  ca = {};
  ca.max_path_cost = Absolute {-1};
  rejects(cfg(ca, FlowPlacement::MaxFlow), "max_path_cost Absolute cost must be >= 0 (got -1)");
  for (const auto& [f, shown] : {std::pair {0.5, "0.500000"}, std::pair {std::numeric_limits<double>::quiet_NaN(), "nan"},
                                 std::pair {std::numeric_limits<double>::infinity(), "inf"}}) {
    ca = {};
    ca.max_path_cost = Stretch {f};
    rejects(cfg(ca, FlowPlacement::MaxFlow),
            std::string("max_path_cost Stretch factor must be finite and in [1, inf) (got ") + shown + ")");
  }
  ca = {};
  ca.diminishing = DiminishingReturns {0, 0.5};
  rejects(cfg(ca, FlowPlacement::MaxFlow), "diminishing window must be >= 1 (got 0)");
  for (const auto& [f, shown] : {std::pair {0.0, "0.000000"}, std::pair {1.0, "1.000000"},
                                 std::pair {std::numeric_limits<double>::quiet_NaN(), "nan"}}) {
    ca = {};
    ca.diminishing = DiminishingReturns {1, f};
    rejects(cfg(ca, FlowPlacement::MaxFlow), std::string("diminishing fraction must be in (0, 1) (got ") + shown + ")");
  }
  Pinned none;
  rejects(cfg(none, FlowPlacement::MaxFlow), "Pinned: an empty bundle list is rejected");
  // A null bundle, first or after a valid one: refused before it is read.
  Pinned null_first;
  null_first.bundles = {nullptr};
  rejects(cfg(null_first, FlowPlacement::MaxFlow), "Pinned: bundle 0 is null");
  Pinned null_second;
  null_second.bundles = {dag_of(*g, {0, 2}), nullptr};
  rejects(cfg(null_second, FlowPlacement::MaxFlow), "Pinned: bundle 1 is null");
  Pinned dup;
  dup.bundles = {dag_of(*g, {0, 2}), dag_of(*g, {0, 2})};
  rejects(cfg(dup, FlowPlacement::MaxFlow), "Pinned: duplicate bundle 1 (equal to bundle 0 before pruning)");
  Pinned malformed;
  malformed.bundles = {dag_of(*g, {0, 3})};
  rejects(cfg(malformed, FlowPlacement::MaxFlow), "Pinned: malformed bundle 0, no src -> dst walk (src 0, dst 3)");
  Pinned ok;
  ok.bundles = {dag_of(*g, {0, 2})};
  rejects(cfg(ok, FlowPlacement::MaxFlow),
          "Pinned routing requires single endpoints (bundles bind one source and one destination)", {0, 1}, {3});
  auto other = diamond(5, 5);
  Pinned foreign;
  foreign.bundles = {dag_of(*other, {0, 2})};
  rejects(cfg(foreign, FlowPlacement::MaxFlow), "foreign PredDAG: bundle built for another graph");
  rejects(cfg(CostOnly {}, FlowPlacement::MaxFlow),
          "node 3 is in both endpoint sets: sources and destinations must not overlap", {0, 3}, {3});
  rejects(cfg(CostOnly {}, FlowPlacement::MaxFlow), "sources: the node set is empty", {}, {3});
  rejects(cfg(CostOnly {}, FlowPlacement::MaxFlow), "sources: node 9 out of node range [0, 4)", {9}, {3});
  // A routing variant left without a mode (valueless) by an exception while a mode was emplaced.
  struct ThrowingConversion {
    operator Pinned() const { throw std::runtime_error("conversion failed"); }
  };
  auto valueless = cfg(CostOnly {}, FlowPlacement::MaxFlow);
  EXPECT_THROW(valueless.routing.emplace<Pinned>(ThrowingConversion {}), std::runtime_error);
  ASSERT_TRUE(valueless.routing.valueless_by_exception());
  rejects(valueless, "routing: the variant holds no routing mode");
  const auto el = other->eligibility({}, std::vector<EdgeId> {0});
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({3}), cfg(CostOnly {}, FlowPlacement::MaxFlow), &el); }),
            "foreign Eligibility: built for another graph");
  DagCache foreign_cache(*other);
  EXPECT_EQ(error_of([&] {
              FlowPolicy p(fg, v({0}), v({3}), cfg(CostOnly {}, FlowPlacement::MaxFlow), nullptr, &foreign_cache);
            }),
            "foreign DagCache: bound to another graph");
  EXPECT_TRUE(detail::LedgerAccess::bound_owners(fg).empty());
}

// Enumeration checks at construction: PathSelection and FlowPolicyConfig reject a value outside its
// enumeration (here 97 to 99, reachable by a cast) with std::invalid_argument naming the field, in
// validate_config's order; a selection assigned out of its domain after its own construction is
// rejected by the configuration. Combinations are the consumers' (EveryRejectedCombinationBefore-
// Binding above and SeveralDefectsReportTheFirstInTheContractOrder below).
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

// Enumeration checks at every consuming entry: the fields stay assignable after construction, so a value
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

// validate_config with several defects at once: one error, the first in the order
// placement, admission, the selection's tie_break and priority (the fields assigned out of their
// enumerations after construction), then the kernel combinations (Lossy with MaxFlow, Weighted
// with Headroom), then the routing rules (Lossy with Global or CapacityAware, Global without
// MaxFlow, CapacityAware's counts, Pinned's bundles). Exact texts: the selection's priority is
// named "selection.priority", never the kernel check's "priority". The same first error at
// construction and at the evaluators; nothing is bound.
TEST(FlowPolicyConfig, SeveralDefectsReportTheFirstInTheContractOrder) {
  auto g = diamond(5, 5);
  FlowLedger fg(g);
  const auto bad_placement = static_cast<FlowPlacement>(99);
  const auto bad_admission = static_cast<Admission>(98);
  const auto bad_tie = static_cast<EdgeTieBreak>(97);
  const auto bad_priority = static_cast<EdgePriority>(96);
  CapacityAware zero;
  zero.flow_count = 0;
  struct Case {
    FlowPolicyConfig c;
    std::string error;
  };
  const auto with = [](FlowPolicyConfig c, std::optional<FlowPlacement> k, std::optional<Admission> a,
                       std::optional<EdgeTieBreak> t, std::optional<EdgePriority> p) {
    if (k) c.placement = *k;
    if (a) c.admission = *a;
    auto& selection = c.selection;
    if (t) selection.tie_break = *t;
    if (p) selection.priority = *p;
    return c;
  };
  const std::string lossy_routing =
      "Lossy admission is rejected with Global and CapacityAware routing (hop-by-hop forwarding without admission "
      "control neither grows nor reroutes)";
  const std::vector<Case> cases = {
      {with(cfg(Global {}, FlowPlacement::MaxFlow, Admission::Lossy), bad_placement, bad_admission, bad_tie, bad_priority),
       "placement: 99 is not a value of FlowPlacement"},
      {with(cfg(Global {}, FlowPlacement::Weighted), std::nullopt, bad_admission, bad_tie, bad_priority),
       "admission: 98 is not a value of Admission"},
      {with(cfg(zero, FlowPlacement::MaxFlow, Admission::Lossy), std::nullopt, std::nullopt, bad_tie, bad_priority),
       "selection.tie_break: 97 is not a value of EdgeTieBreak"},
      {with(cfg(zero, FlowPlacement::MaxFlow, Admission::Lossy), std::nullopt, std::nullopt, std::nullopt, bad_priority),
       "selection.priority: 96 is not a value of EdgePriority"},
      {with(cfg(Global {}, FlowPlacement::Weighted), std::nullopt, std::nullopt, std::nullopt, bad_priority),
       "selection.priority: 96 is not a value of EdgePriority"},
      {cfg(Global {}, FlowPlacement::MaxFlow, Admission::Lossy), kLossyKernelError},
      {cfg(zero, FlowPlacement::MaxFlow, Admission::Lossy), kLossyKernelError},
      {cfg(Global {}, FlowPlacement::Weighted, Admission::Lossless, sel(EdgeTieBreak::Deterministic, EdgePriority::Headroom)),
       kHeadroomKernelError},
      {cfg(zero, FlowPlacement::Weighted, Admission::Lossy, sel(EdgeTieBreak::Deterministic, EdgePriority::Headroom)),
       kHeadroomKernelError},
      {cfg(Pinned {}, FlowPlacement::MaxFlow, Admission::Lossy), kLossyKernelError},
      {cfg(zero, FlowPlacement::EqualBalanced, Admission::Lossy), lossy_routing},
      {cfg(Global {}, FlowPlacement::Weighted, Admission::Lossy), lossy_routing},
  };
  for (const auto& k : cases) {
    EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({3}), k.c); }), k.error);
    EXPECT_TRUE(detail::LedgerAccess::bound_owners(fg).empty());
    EXPECT_TRUE(fg.records().empty());
    CapacityOptions co {k.c};
    EXPECT_EQ(error_of([&] { (void)capacity(*g, v({0}), v({3}), co); }), k.error);
    EXPECT_EQ(error_of([&] { (void)sensitivity(*g, v({0}), v({3}), SensitivityOptions {k.c}, 1); }), k.error);
    if (std::holds_alternative<Pinned>(k.c.routing)) continue;  // batch_capacity rejects Pinned first
    const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{v({0}), v({3})}};
    EXPECT_EQ(error_of([&] { (void)batch_capacity(*g, pairs, co, {}, 1); }), k.error);
  }
}

// The four axes are required: a configuration has no default.
static_assert(!std::is_default_constructible_v<FlowPolicyConfig>);

// config() is the copy taken at construction: every field as given, and a later change of the
// caller's object changes nothing.
TEST(FlowPolicyConfig, ConfigIsTheCopyTakenAtConstruction) {
  auto g = diamond(5, 5);
  FlowLedger fg(g);
  CapacityAware ca;
  ca.flow_count = 2;
  ca.max_flow_count = 3;
  ca.max_path_cost = Absolute {9};
  ca.min_flow_bandwidth = 1;
  ca.reoptimize = true;
  ca.diminishing = DiminishingReturns {2, 0.25};
  const auto selection = sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, false, false);
  auto c = cfg(ca, FlowPlacement::EqualBalanced, Admission::Lossless, selection);
  FlowPolicy p(fg, v({0}), v({3}), c);
  c.placement = FlowPlacement::MaxFlow;
  c.selection = sel();
  std::get<CapacityAware>(c.routing).flow_count = 7;
  const auto& got = p.config();
  EXPECT_NE(&got, &c);
  EXPECT_EQ(got.placement, FlowPlacement::EqualBalanced);
  EXPECT_EQ(got.admission, Admission::Lossless);
  EXPECT_TRUE(got.selection == selection);
  ASSERT_TRUE(std::holds_alternative<CapacityAware>(got.routing));
  const auto& r = std::get<CapacityAware>(got.routing);
  EXPECT_EQ(r.flow_count, 2);
  EXPECT_EQ(r.max_flow_count, std::optional<std::int32_t> {3});
  EXPECT_TRUE(r.max_path_cost == std::optional<CostCeiling> {Absolute {9}});
  EXPECT_EQ(r.min_flow_bandwidth, 1);
  EXPECT_TRUE(r.reoptimize);
  ASSERT_TRUE(r.diminishing.has_value());
  EXPECT_EQ(r.diminishing->window, 2);
  EXPECT_EQ(r.diminishing->fraction, 0.25);
  // Pinned: the supplied bundles themselves, in supply order.
  Pinned pin;
  pin.bundles = {dag_of(*g, {1, 3}), dag_of(*g, {0, 2})};
  FlowLedger fg2(g);
  FlowPolicy q(fg2, v({0}), v({3}), cfg(pin, FlowPlacement::MaxFlow));
  ASSERT_TRUE(std::holds_alternative<Pinned>(q.config().routing));
  const auto& kept = std::get<Pinned>(q.config().routing).bundles;
  ASSERT_EQ(kept.size(), 2U);
  EXPECT_EQ(kept[0].get(), pin.bundles[0].get());
  EXPECT_EQ(kept[1].get(), pin.bundles[1].get());
}

// Global and Pinned read only the selection's priority, and only with the Weighted kernel (its
// membership and weights): every other selection gives the same placement. Two parallel edges
// 0 -> 1 of 4 units (EdgeIds 0 and 1), a third party holding 3 units of edge 0. Pinned Weighted
// on both edges, request 5: with priority Capacity the weights are 4 : 4, so the flow is at most
// 1 + 1 = 2 (edge 0 has 1 unit left); with priority Residual they are the residuals 1 : 4 at
// creation, so it carries 1 + 4 = 5. Global (MaxFlow only), Pinned MaxFlow and Pinned
// EqualBalanced place the same under every selection. Control: CostOnly reads the selection (a
// single path, multipath = false, carries less than the two-edge DAG).
TEST(FlowPolicyConfig, GlobalAndPinnedReadOnlyTheWeightedPriority) {
  auto g = shared(make_graph(2, {{0, 1, 4, 1}, {0, 1, 4, 1}}));
  std::vector<PathSelection> all;
  for (auto t : {EdgeTieBreak::Deterministic, EdgeTieBreak::PreferHigherPriority}) {
    for (auto p : {EdgePriority::Residual, EdgePriority::Capacity, EdgePriority::Headroom}) {
      for (bool mp : {true, false}) {
        for (bool me : {true, false}) all.push_back(sel(t, p, mp, me));
      }
    }
  }
  using Edges = std::vector<std::pair<EdgeId, Amount>>;
  // The placement of request 5 after the third party: (placed, the policy's records' edges).
  auto run = [&](const Routing& routing, FlowPlacement k, const PathSelection& s) {
    FlowLedger fg(g);
    Pinned third;
    third.bundles = {dag_of(*g, {0})};
    FlowPolicy other(fg, v({0}), v({1}), cfg(third, FlowPlacement::MaxFlow));
    EXPECT_EQ(other.place(3).placed, 3);
    FlowPolicy p(fg, v({0}), v({1}), cfg(routing, k, Admission::Lossless, s));
    const auto r = p.place(5);
    std::vector<Edges> edges;
    for (const auto& rec : p.records()) edges.push_back(rec.edges);
    return std::pair {r.placed, edges};
  };
  Pinned both;
  both.bundles = {dag_of(*g, {0, 1})};
  for (const auto& [name, routing, k] : std::vector<std::tuple<std::string, Routing, FlowPlacement>> {
           {"Global", Global {}, FlowPlacement::MaxFlow},
           {"Pinned MaxFlow", both, FlowPlacement::MaxFlow},
           {"Pinned EqualBalanced", both, FlowPlacement::EqualBalanced}}) {
    const auto want = run(routing, k, all.front());
    EXPECT_EQ(want.first, k == FlowPlacement::EqualBalanced ? 2 : 5) << name;
    for (const auto& s : all) EXPECT_EQ(run(routing, k, s), want) << name;
  }
  for (const auto& s : all) {
    if (s.priority == EdgePriority::Headroom) continue;  // rejected with Weighted
    const bool residual = s.priority == EdgePriority::Residual;
    const auto got = run(both, FlowPlacement::Weighted, s);
    EXPECT_EQ(got.first, residual ? 5 : 2);
    EXPECT_EQ(got.second, (std::vector<Edges> {residual ? Edges {{0, 1}, {1, 4}} : Edges {{0, 1}, {1, 1}}}));
  }
  EXPECT_EQ(run(CostOnly {}, FlowPlacement::MaxFlow, sel()).first, 5);
  EXPECT_EQ(run(CostOnly {}, FlowPlacement::MaxFlow, sel(EdgeTieBreak::Deterministic, EdgePriority::Capacity, false)).first,
            1);
}

namespace {

// What a rejected call must leave as it was: the policy's epoch state (the requested total, the
// +inf flag, the placed demand, the Global curve and each roster source's flows with their routes)
// and the ledger (every record field, the allocations, the residual, the next FlowId and the
// owners). A flow's source is its roster source's, and its record's sources are compared.
struct PolicyState {
  Epoch epoch;
  LedgerState ledger;
};

PolicyState state_of(const FlowPolicy& p, const FlowLedger& fg) { return {test::epoch(p), ledger_state(fg)}; }

void expect_same(const PolicyState& a, const PolicyState& b, const std::string& where) {
  expect_same(a.epoch, b.epoch, where);
  expect_same(a.ledger, b.ledger, where);
}

// Each source's delivered total, by the records' single source (a source with a flow and no
// delivery included).
std::map<NodeId, Amount> delivered_by_source(const FlowPolicy& p) {
  std::map<NodeId, Amount> out;
  for (const auto& r : p.records()) out[r.sources.at(0)] += r.delivered;
  return out;
}

} // namespace

// A request rejected by validation (out of [0, 2^62), finite after +inf, a cumulative total
// that would reach 2^127) changes neither the policy's epoch state nor the ledger, under every
// routing model; validation precedes the shortcut of a zero request (place(0) after +inf is
// rejected too), and remove() ends the epoch, after which a finite request places again. The
// 2^127 case starts from a requested total set just below the limit.
TEST(FlowPolicyRequests, RejectedRequestsLeaveThePolicyAndTheLedger) {
  auto g = diamond(5, 3, 2);
  Pinned pinned;
  pinned.bundles = {dag_of(*g, {0, 2}), dag_of(*g, {1, 3})};
  CapacityAware ca;
  ca.max_path_cost = Absolute {4};
  const std::vector<std::pair<const char*, FlowPolicyConfig>> configs = {
      {"CostOnly", cost_only_config()},
      {"CapacityAware", cfg(ca, FlowPlacement::EqualBalanced)},
      {"Global", global_config()},
      {"Pinned", cfg(pinned, FlowPlacement::MaxFlow)},
      {"CostOnly Lossy", cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossy)}};
  const std::string range_error = " is outside [0, 2^62) units (requests are integer units or +inf)";
  const std::string after_inf =
      "a finite request after +inf is rejected within an epoch (the cumulative finite request is undefined)";
  for (const auto& [name, config] : configs) {
    FlowLedger fg(g);
    FlowPolicy p(fg, v({0}), v({3}), config);
    ASSERT_GT(p.place(4).placed, 0) << name;
    const auto finite = state_of(p, fg);
    EXPECT_EQ(error_of([&] { (void)p.place(-1); }), "request -1" + range_error) << name;
    EXPECT_EQ(error_of([&] { (void)p.place(kDomainLimit); }), "request 4611686018427387904" + range_error) << name;
    expect_same(state_of(p, fg), finite, name);
    // Just below 2^127: 5 more would reach it, 4 would not.
    auto& requested = detail::PolicyAccess::engine(p).requested;
    const Counter128 saved = requested;
    requested = Counter128 {(std::uint64_t {1} << 63) - 1, std::numeric_limits<std::uint64_t>::max() - 4};
    const auto near_limit = state_of(p, fg);
    EXPECT_EQ(error_of([&] { (void)p.place(5); }), "the cumulative request counter would reach 2^127") << name;
    expect_same(state_of(p, fg), near_limit, name);
    requested = saved;
    (void)p.place(std::nullopt);
    const auto infinite = state_of(p, fg);
    ASSERT_TRUE(infinite.epoch.infinite) << name;
    EXPECT_EQ(error_of([&] { (void)p.place(1); }), after_inf) << name;
    EXPECT_EQ(error_of([&] { (void)p.place(0); }), after_inf) << name;
    EXPECT_EQ(error_of([&] { (void)p.place(-1); }), "request -1" + range_error) << name;
    expect_same(state_of(p, fg), infinite, name);
    p.remove();
    EXPECT_EQ(p.place(1).placed, 1) << name;
  }
}

// The Lossy counters' 2^127 limit. Finite calls bring a flow's offered total up to the request
// total (below 2^127) and a +inf call then adds its injection, so both checks are reachable; the
// state about 2^65 calls would build is set here through the ledger's commit (a drop counter
// never above its flow's offered counter). On the path 0 -> 1 -> 2 with capacities 10 and 4,
// place(1) leaves residuals 9 and 3; place(+inf) then offers 9 (its injection on edge 0), edge 1
// carries 3 and drops 6. From L = 2^127: an offered total of L - 10 reaches L - 1 (accepted), L - 9
// would reach L (rejected); a drop counter of L - 7 reaches L - 1 and passes, so the offered
// check rejects the call; L - 6 would reach L (rejected by the drop check, which runs first).
// A rejected call leaves the policy and the ledger as they were.
TEST(FlowPolicyRequests, LossyCountersStopBelow2To127) {
  auto g = shared(make_graph(3, {{0, 1, 10, 1}, {1, 2, 4, 1}}));
  const auto below = [](std::uint64_t x) {
    return Counter128 {(std::uint64_t {1} << 63) - 1, std::numeric_limits<std::uint64_t>::max() - (x - 1)};
  };
  const std::string offered_error = "the cumulative offered counter would reach 2^127";
  const std::string drop_error = "a cumulative drop counter would reach 2^127";
  struct Case {
    std::uint64_t offered_below;
    std::optional<std::uint64_t> drop_below;
    std::optional<std::string> error;
  };
  for (const auto& c : {Case {10, std::nullopt, std::nullopt}, Case {9, std::nullopt, offered_error},
                        Case {7, 7, offered_error}, Case {6, 6, drop_error}}) {
    const std::string name = "L - " + std::to_string(c.offered_below);
    FlowLedger fg(g);
    FlowPolicy p(fg, v({0}), v({2}), cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossy));
    ASSERT_EQ(p.place(1).placed, 1) << name;
    std::vector<detail::Draft> drafts(1);
    drafts[0].record = p.records().at(0);
    drafts[0].id = drafts[0].record.flow_id.second;
    drafts[0].record.offered = below(c.offered_below);
    if (c.drop_below.has_value()) drafts[0].record.drops = {{1, below(*c.drop_below)}};
    (void)detail::LedgerAccess::commit(fg, detail::PolicyAccess::engine(p).owner, drafts);
    const auto before = state_of(p, fg);
    if (c.error.has_value()) {
      EXPECT_EQ(error_of([&] { (void)p.place(std::nullopt); }), *c.error) << name;
      expect_same(state_of(p, fg), before, name);
      continue;
    }
    const auto r = p.place(std::nullopt);
    EXPECT_EQ(r.placed, 3) << name;
    EXPECT_EQ(r.drops, (std::vector<std::pair<EdgeId, Amount>> {{1, 6}})) << name;
    const auto records = p.records();
    ASSERT_EQ(records.size(), 1U) << name;
    EXPECT_EQ(records[0].offered, below(1)) << name;
    EXPECT_EQ(records[0].drops, (std::vector<std::pair<EdgeId, Counter128>> {{1, Counter128 {0, 6}}})) << name;
    EXPECT_EQ(records[0].delivered, 4) << name;
  }
}

// The per-source shares of the cumulative request R (round robin from the lowest source): after
// R units source i of k has floor(R / k) + [i < R mod k]. Sources 0 .. k - 1 with their own
// 100-unit edge to k. k = 3: place(2) from R = 0 offers [1, 1, 0]; the next place(2) (R = 2)
// offers [1, 0, 1], which a rejected call that changed R would move. k = 7, where 2^64 = 2 mod 7:
// from R = 2^64 place(1) offers its unit to source 2, from R = 2^64 + 5 (= 0 mod 7) to source 0,
// so the high word of R counts with its exact weight.
TEST(FlowPolicyRequests, RoundRobinSharesFollowTheCumulativeRequest) {
  auto star = [](std::int32_t k) {
    std::vector<InEdge> rows;
    for (std::int32_t s = 0; s < k; ++s) rows.push_back({s, k, 100, 1});
    return shared(make_graph(k + 1, rows));
  };
  auto sources = [](std::int32_t k) {
    std::vector<NodeId> out(static_cast<std::size_t>(k));
    std::iota(out.begin(), out.end(), 0);
    return out;
  };
  using Shares = std::map<NodeId, Amount>;
  {
    auto g = star(3);
    FlowLedger fg(g);
    FlowPolicy p(fg, sources(3), v({3}), cost_only_config());
    EXPECT_EQ(p.place(2).placed, 2);
    EXPECT_EQ(delivered_by_source(p), (Shares {{0, 1}, {1, 1}, {2, 0}}));
    EXPECT_EQ(error_of([&] { (void)p.place(-1); }),
              "request -1 is outside [0, 2^62) units (requests are integer units or +inf)");
    EXPECT_EQ(error_of([&] { (void)p.place(kDomainLimit); }),
              "request 4611686018427387904 is outside [0, 2^62) units (requests are integer units or +inf)");
    EXPECT_EQ(p.place(2).placed, 2);
    EXPECT_EQ(delivered_by_source(p), (Shares {{0, 2}, {1, 1}, {2, 1}}));
  }
  auto g = star(7);
  for (const auto& [requested, to] : {std::pair {Counter128 {1, 0}, NodeId {2}}, std::pair {Counter128 {1, 5}, NodeId {0}}}) {
    FlowLedger fg(g);
    FlowPolicy p(fg, sources(7), v({7}), cost_only_config());
    detail::PolicyAccess::engine(p).requested = requested;
    EXPECT_EQ(p.place(1).placed, 1) << requested.lo;
    for (const auto& [s, x] : delivered_by_source(p)) EXPECT_EQ(x, s == to ? 1 : 0) << requested.lo << " source " << s;
    EXPECT_EQ(delivered_by_source(p).at(to), 1) << requested.lo;
  }
}

TEST(FlowPolicyOutputs, OwnAllocationLedgerResidualBottleneckAndGlobalCut) {
  auto g = diamond(5, 3);
  FlowLedger fg(g);
  // A third party holds one unit of the first route.
  Pinned third;
  third.bundles = {dag_of(*g, {0, 2})};
  FlowPolicy other(fg, v({0}), v({3}), cfg(third, FlowPlacement::MaxFlow));
  ASSERT_EQ(other.place(1).placed, 1);
  FlowPolicy p(fg, v({0}), v({3}), cost_only_config());
  const OutputOptions all {true, true, true};
  auto r = p.place(4, all);
  // Both routes cost 2: EqualBalanced puts 2 units on each, so every edge carries 2 of this
  // policy's units and the residuals are 5 - 1 - 2 on route 0 -> 1 -> 3 and 3 - 2 on the other.
  // Every edge keeps a residual unit, so the sources reach everything: no saturated edge.
  EXPECT_EQ(r.placed, 4);
  ASSERT_TRUE(r.edge_flows && r.residual && r.bottleneck);
  EXPECT_FALSE(r.min_cut.has_value());  // not Global
  EXPECT_FALSE(r.cost_distribution.has_value());
  std::vector<Amount> own(4, 0);
  for (const auto& rec : p.records()) {
    for (const auto& [e, x] : rec.edges) own[static_cast<std::size_t>(e)] += x;
  }
  EXPECT_EQ(*r.edge_flows, (std::vector<Amount> {2, 2, 2, 2}));
  EXPECT_EQ(*r.edge_flows, own);
  EXPECT_EQ(*r.residual, (std::vector<Amount> {2, 1, 2, 1}));
  EXPECT_EQ(*r.residual, fg.residual());
  EXPECT_EQ(*r.bottleneck, std::vector<EdgeId> {});
  EXPECT_EQ(*r.bottleneck, fg.bottleneck(v({0})));
  // Nothing requested: absent.
  r = p.place(1);
  EXPECT_FALSE(r.edge_flows || r.residual || r.bottleneck || r.min_cut);
  // Global: the curve always, min_cut with bottleneck when the solve is maximal (at +inf; none
  // under a cap of 2 below the maximum 8, which leaves a residual path).
  FlowLedger fg2(g);
  FlowPolicy gl(fg2, v({0}), v({3}), global_config());
  r = gl.place(2, OutputOptions {false, false, true});
  EXPECT_EQ(r.placed, 2);
  ASSERT_TRUE(r.cost_distribution.has_value());
  EXPECT_EQ(*r.cost_distribution, curve({{2, 2}}));
  EXPECT_EQ(gl.placed_demand(), 2);
  EXPECT_FALSE(r.min_cut.has_value());
  // At +inf the remaining 6 of the maximum 5 + 3: both edges out of node 0 saturate, so the
  // closure is {0} and the cut is 0 -> 1 and 0 -> 2 (EdgeIds 0 and 1).
  r = gl.place(std::nullopt, OutputOptions {false, false, true});
  EXPECT_EQ(r.placed, 6);
  EXPECT_EQ(gl.placed_demand(), 8);
  ASSERT_TRUE(r.min_cut.has_value());
  EXPECT_EQ(*r.min_cut, (std::vector<EdgeId> {0, 1}));
}

// The bottleneck output is the ledger's saturation report for the policy's sources under the
// policy's Eligibility. Edges (input order): a = 0 -> 3 (2 units), x = 0 -> 2 (5), y = 2 -> 3 (1). A
// third party carries 1 unit on x and y, which saturates y. The policy excludes x and carries 2
// units on a, which saturates a. Under the Eligibility the sources reach only node 0: the report is
// [a]. Without it node 2 is reachable too and y is reported as well.
TEST(FlowPolicyOutputs, BottleneckIsTheLedgersForThePolicysSourcesAndEligibility) {
  auto g = shared(make_graph(4, {{0, 3, 2, 1}, {0, 2, 5, 1}, {2, 3, 1, 1}}));
  const auto id = id_of_input(*g);
  const EdgeId a = id[0], x = id[1], y = id[2];
  FlowLedger fg(g);
  Pinned third;
  third.bundles = {dag_of(*g, {x, y})};
  FlowPolicy other(fg, v({0}), v({3}), cfg(third, FlowPlacement::MaxFlow));
  ASSERT_EQ(other.place(1).placed, 1);
  const auto el = g->eligibility({}, std::vector<EdgeId> {x});
  FlowPolicy p(fg, v({0}), v({3}), cost_only_config(), &el);
  const auto r = p.place(2, OutputOptions {false, false, true});
  EXPECT_EQ(r.placed, 2);
  ASSERT_TRUE(r.bottleneck.has_value());
  EXPECT_EQ(*r.bottleneck, std::vector<EdgeId> {a});
  EXPECT_EQ(*r.bottleneck, fg.bottleneck(v({0}), &el));
  auto without = std::vector<EdgeId> {a, y};
  std::sort(without.begin(), without.end());
  EXPECT_EQ(fg.bottleneck(v({0})), without);  // the Eligibility matters here
}

TEST(FlowPolicyGlobal, EmptyRosterReportsAnEmptyCurveAndAnEmptyCut) {
  // No source reaches the destination: the maximum flow is 0 and maximal.
  auto g = shared(make_graph(3, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  FlowPolicy p(fg, v({0}), v({2}), global_config());
  auto r = p.place(std::nullopt, OutputOptions {false, false, true});
  EXPECT_EQ(r.termination, Termination::NoRoute);
  ASSERT_TRUE(r.cost_distribution.has_value());
  EXPECT_TRUE(r.cost_distribution->empty());
  ASSERT_TRUE(r.min_cut.has_value());
  EXPECT_TRUE(r.min_cut->empty());
  EXPECT_TRUE(p.records().empty());
  // The zero request is a no-op: no solve, no cut; the cumulative curve is empty.
  FlowPolicy q(fg, v({1}), v({2}), global_config());
  r = q.place(0, OutputOptions {false, false, true});
  EXPECT_FALSE(r.min_cut.has_value());
  ASSERT_TRUE(r.cost_distribution.has_value());
  EXPECT_TRUE(r.cost_distribution->empty());
}

// Global with a roster whose maximum flow is 0: NoRoute while nothing was ever placed, Exhausted
// once something was. One edge 0 -> 1 of 2 units that a third party fills: the source statically
// reaches the destination (it is on the roster), and place(1) finds no flow: NoRoute for the call
// and the source, its record created empty (no bundle, no route cost, nothing delivered). After the
// third party leaves, place(1) delivers 1; a third party then takes the other unit, and the next
// place(1) re-solves 1 + 1 on a view of 1: nothing more, Exhausted.
TEST(FlowPolicyGlobal, ARosterWithoutFlowIsNoRouteUntilSomethingWasPlaced) {
  auto g = shared(make_graph(2, {{0, 1, 2, 1}}));
  FlowLedger fg(g);
  Pinned edge;
  edge.bundles = {dag_of(*g, {0})};
  auto third = std::make_unique<FlowPolicy>(fg, v({0}), v({1}), cfg(edge, FlowPlacement::MaxFlow));
  ASSERT_EQ(third->place(2).placed, 2);
  FlowPolicy p(fg, v({0}), v({1}), global_config());
  auto r = p.place(1);
  EXPECT_EQ(r.placed, 0);
  EXPECT_EQ(r.termination, Termination::NoRoute);
  EXPECT_EQ(r.local_outcomes, (std::vector<std::pair<NodeId, Termination>> {{0, Termination::NoRoute}}));
  ASSERT_TRUE(r.cost_distribution.has_value());
  EXPECT_TRUE(r.cost_distribution->empty());
  ASSERT_EQ(p.records().size(), 1U);
  EXPECT_TRUE(ascending_edges(*p.records()[0].bundle).empty());
  EXPECT_EQ(p.records()[0].route_cost, kUnreachable);
  EXPECT_EQ(p.records()[0].delivered, 0);
  third->remove();
  r = p.place(1);
  EXPECT_EQ(r.placed, 1);
  EXPECT_EQ(r.termination, Termination::Complete);
  FlowPolicy fourth(fg, v({0}), v({1}), cfg(edge, FlowPlacement::MaxFlow));
  ASSERT_EQ(fourth.place(1).placed, 1);
  r = p.place(1);
  EXPECT_EQ(r.placed, 0);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  EXPECT_EQ(r.local_outcomes, (std::vector<std::pair<NodeId, Termination>> {{0, Termination::Exhausted}}));
}

// A Global zero call reports the cumulative curve (its amounts sum to placed_demand();
// empty before the first placement and after remove()) and changes nothing; it solves nothing,
// so it claims no cut.
TEST(FlowPolicyGlobal, ZeroCallReportsTheCumulativeCurve) {
  using Curve = std::vector<std::pair<Cost, Amount>>;
  auto g = shared(make_graph(3, {{0, 1, 3, 7}, {0, 2, 4, 1}, {2, 1, 4, 1}}));
  FlowLedger fg(g);
  FlowPolicy p(fg, v({0}), v({1}), global_config());
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

// The Global layer on random graphs (nonuniform costs 1 to 9 times a factor of 1, 1,000,003 or
// 2^35, zero capacities, a node or edge mask, another owner's occupancy W): each call re-solves the policy's cumulative amount on the view
// W + own. After each of up to two capped calls and after place(+inf), the policy's own flow is
// feasible on the view, conserved, carries the cumulative amount at the minimum cost of that
// amount (the independent Bellman-Ford curve, not the flow core), its records sum to it, and its
// cumulative curve is the optimal curve; at +inf the reported cut separates the sources from the
// destinations with a capacity equal to the flow. capacity() on a fresh ledger reaches the same
// maximum and curve, and batch_capacity() agrees with capacity() per pair.
TEST(FlowPolicyGlobal, RandomCumulativeSolvesAreOptimalWithACertifiedCut) {
  std::mt19937_64 rng(20261010);
  const OutputOptions all {true, true, true};
  auto below = [&](std::int64_t k) { return static_cast<std::int32_t>(rng() % static_cast<std::uint64_t>(k)); };
  int capped = 0, masked = 0, occupied = 0, with_flow = 0;
  int evaluated = 0;  // trials whose evaluators are compared (no other owner)
  for (int trial = 0; trial < 400; ++trial) {
    const std::int32_t n = 4 + static_cast<std::int32_t>(rng() % 5);
    std::vector<InEdge> rows;
    const int m = 6 + static_cast<int>(rng() % 11);
    // The cost factor of the trial, outside the random stream: 1, 1,000,003 or 2^35.
    const Cost factor = trial % 3 == 0 ? 1 : trial % 3 == 1 ? 1000003 : Cost {1} << 35;
    for (int i = 0; i < m; ++i) {
      const auto u = static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n));
      auto w = static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n - 1));
      if (w >= u) ++w;
      rows.push_back({u, w, static_cast<Amount>(rng() % 7), static_cast<Cost>(1 + rng() % 9) * factor});
    }
    auto g = shared(make_graph(n, rows));
    const std::string where = "trial " + std::to_string(trial);
    // Sources 0 and maybe 1, destinations n - 1 and maybe n - 2.
    std::vector<NodeId> srcs = {0}, dsts = {n - 1};
    if (rng() % 2 == 0) srcs.push_back(1);
    if (rng() % 2 == 0) dsts.push_back(n - 2);
    std::optional<Eligibility> el;
    switch (rng() % 3) {
      case 0: el = g->eligibility(std::vector<NodeId> {below(n)}, {}); break;
      case 1: el = g->eligibility({}, std::vector<EdgeId> {below(m)}); break;
      default: break;
    }
    const Eligibility* elig = el ? &*el : nullptr;
    masked += elig != nullptr ? 1 : 0;
    FlowLedger fg(g);
    // Another owner: a Global policy between other endpoints (when the graph allows them) takes
    // part of the capacity first.
    std::optional<FlowPolicy> other;
    if (n >= 6 && rng() % 2 == 0) {
      other.emplace(fg, v({2}), v({3}), global_config());
      occupied += other->place(static_cast<Amount>(rng() % 5)).placed > 0 ? 1 : 0;
    }
    const std::vector<Amount> view = fg.residual();  // W: the policy's own flow is 0 so far
    const auto curve_pts = independent_curve(*g, view, srcs, dsts, elig);
    const Amount max_flow = curve_pts.back().first;
    FlowPolicy p(fg, srcs, dsts, global_config(), elig);
    // The policy's own per-edge flow, its cost, and the checks every call ends with.
    auto check_own = [&](Amount total, const PlacementResult& r, const std::string& at) {
      std::vector<Amount> own(view.size(), 0);
      Amount recorded = 0;
      for (const auto& rec : p.records()) {
        recorded += rec.delivered;
        for (const auto& [e, x] : rec.edges) own[static_cast<std::size_t>(e)] += x;
      }
      EXPECT_EQ(recorded, total) << at;
      EXPECT_EQ(p.placed_demand(), total) << at;
      ASSERT_TRUE(r.edge_flows.has_value() && r.residual.has_value()) << at;
      EXPECT_EQ(*r.edge_flows, own) << at;
      EXPECT_EQ(*r.residual, fg.residual()) << at;
      // The flow certificates' flow clauses (test_utils.hpp): within the view, nothing on a masked
      // edge, conserved, the cumulative amount out of the sources and into the destinations.
      std::vector<std::string> fails;
      (void)flow_clauses(*g, view, own, srcs, dsts, elig, total, fails);
      const Int128 cost = flow_cost(*g, own);
      EXPECT_TRUE(cost == curve_at(curve_pts, total)) << at;
      // The cumulative curve, by the flow certificates' clauses: positive volumes at strictly
      // increasing costs, its total the cumulative amount and its moment the flow's cost, and equal
      // to the optimal curve at every breakpoint of either up to the total.
      ASSERT_TRUE(r.cost_distribution.has_value()) << at;
      const auto pts = curve_shape_clauses(*r.cost_distribution, total, cost, fails);
      curve_optimum_clause(pts, curve_pts, total, "the optimal curve", fails);
      for (const auto& f : fails) ADD_FAILURE() << at << ": " << f;
    };
    // Up to two capped calls below the maximum (each re-solves the cumulative amount), then +inf.
    Amount total = 0;
    for (int call = 0; call < 2; ++call) {
      const Amount x = max_flow - total > 1 ? below(max_flow - total) : 0;
      if (x == 0) continue;
      ++capped;
      total += x;
      const auto r = p.place(x, all);
      EXPECT_EQ(r.placed, x) << where;
      EXPECT_FALSE(r.min_cut.has_value()) << where;  // stopped by the cap
      check_own(total, r, where + " capped call " + std::to_string(call));
    }
    const auto r = p.place(std::nullopt, all);
    EXPECT_EQ(r.placed, max_flow - total) << where;
    with_flow += max_flow > 0 ? 1 : 0;
    check_own(max_flow, r, where + " +inf");
    // The cut, as the flow certificates define it (residual_closure and closure_cut,
    // test_utils.hpp): the residual closure of the allowed sources is the nodes they still reach
    // over allowed edges, forward where the policy's own flow leaves room in the view and backward
    // where it can be undone; no destination is in it, and the cut is every allowed edge of
    // positive view that leaves it, in ascending EdgeId, with a capacity equal to the maximum flow.
    ASSERT_TRUE(r.min_cut.has_value() && r.edge_flows.has_value()) << where;
    Amount cut_capacity = 0;
    for (const EdgeId e : *r.min_cut) cut_capacity += view[static_cast<std::size_t>(e)];
    EXPECT_EQ(cut_capacity, max_flow) << where;
    const auto closure = residual_closure(*g, view, *r.edge_flows, srcs, elig);
    for (const NodeId d : dsts) EXPECT_EQ(closure[static_cast<std::size_t>(d)], 0) << where << " destination " << d;
    EXPECT_EQ(*r.min_cut, closure_cut(*g, view, closure, elig)) << where;
    // The evaluators on a fresh ledger (static capacities): the same maximum and cost.
    if (other.has_value()) continue;
    ++evaluated;
    const CapacityOptions co {global_config(), elig, nullptr, all};
    const auto cr = capacity(*g, srcs, dsts, co);
    EXPECT_EQ(cr.placed, max_flow) << where;
    ASSERT_TRUE(cr.cost_distribution.has_value()) << where;
    EXPECT_EQ(*cr.cost_distribution, *r.cost_distribution) << where;
    const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{srcs, dsts}, {v({0}), v({n - 1})}};
    const auto batch = batch_capacity(*g, pairs, co, {}, 1);
    ASSERT_EQ(batch.size(), 2U) << where;
    for (std::size_t i = 0; i < pairs.size(); ++i) {
      const auto one = capacity(*g, pairs[i].first, pairs[i].second, co);
      EXPECT_EQ(batch[i].placed, one.placed) << where << " pair " << i;
      EXPECT_EQ(batch[i].cost_distribution, one.cost_distribution) << where << " pair " << i;
      EXPECT_EQ(batch[i].min_cut, one.min_cut) << where << " pair " << i;
    }
  }
  // The population reaches every case (this seed: 289 capped calls, 257 masks, 44 occupied
  // ledgers, 241 positive maxima of 400 trials).
  EXPECT_GT(capped, 220);
  EXPECT_GT(masked, 200);
  EXPECT_GT(occupied, 33);
  EXPECT_GT(with_flow, 180);
  EXPECT_GT(evaluated, 100);
}

// CostOnly and CapacityAware placements on random graphs (costs 1 to 9 times a factor of 1,
// 1,000,003 or 2^35, zero capacities, a node or edge mask, one or two sources and destinations,
// another owner's flow first), in the kernels and admissions these models take: after each of two
// capped calls and a call of +inf, each of the policy's records is a flow of its own from the
// policy's sources to its destinations (expect_conserved, test_utils.hpp: positive amounts in
// ascending EdgeId, every transit node passing on what it receives less its drops, delivered the
// destinations' net inflow and, Lossless, the sources' net outflow), carrying and dropping nothing
// on an edge the Eligibility excludes; a Lossless policy's records deliver its placed demand, and
// the ledger's aggregate is the sum of its records (expect_balanced). These models choose routes
// by cost, not by the minimum-cost flow, so no optimality certificate applies here.
TEST(FlowPolicyRecords, RandomCostOnlyAndCapacityAwareRecordsAreFlows) {
  const auto tunnel = sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, false, false);
  const auto res = sel(EdgeTieBreak::Deterministic, EdgePriority::Residual);
  CapacityAware grow;
  grow.flow_count = 1;
  grow.max_flow_count = 4;
  CapacityAware reopt;
  reopt.flow_count = 2;
  reopt.reoptimize = true;
  reopt.max_path_cost = Stretch {2.0};
  CapacityAware dim;
  dim.flow_count = 1;
  dim.diminishing = DiminishingReturns {1, 0.5};
  dim.min_flow_bandwidth = 1;
  const std::vector<FlowPolicyConfig> configs = {
      cost_only_config(),
      cfg(CostOnly {}, FlowPlacement::Weighted),
      cfg(CostOnly {}, FlowPlacement::MaxFlow),
      cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossless, tunnel),
      cfg(CostOnly {}, FlowPlacement::EqualBalanced, Admission::Lossy),
      cfg(CostOnly {}, FlowPlacement::Weighted, Admission::Lossy, res),
      cfg(grow, FlowPlacement::EqualBalanced, Admission::Lossless, tunnel),
      cfg(reopt, FlowPlacement::MaxFlow),
      cfg(dim, FlowPlacement::Weighted, Admission::Lossless, res),
  };
  std::mt19937_64 rng(20261011);
  auto below = [&](std::int64_t k) { return static_cast<std::int32_t>(rng() % static_cast<std::uint64_t>(k)); };
  int carrying_cost_only = 0;       // record checks of a CostOnly record that carries something
  int carrying_capacity_aware = 0;  // the same for CapacityAware
  int dropping = 0;                 // record checks of a Lossy record with drops
  int masked = 0;                   // record checks of a carrying record under a mask
  int scaled = 0;                   // record checks of a carrying record with costs times 1,000,003 or 2^35
  int several = 0;                  // calls that delivered from two sources or into two destinations
  for (int trial = 0; trial < 600; ++trial) {
    const FlowPolicyConfig& config = configs[static_cast<std::size_t>(trial) % configs.size()];
    const bool cost_only = std::holds_alternative<CostOnly>(config.routing);
    const std::int32_t n = 4 + static_cast<std::int32_t>(rng() % 5);
    const int m = 6 + static_cast<int>(rng() % 11);
    const Cost factor = trial % 3 == 0 ? 1 : trial % 3 == 1 ? 1000003 : Cost {1} << 35;
    std::vector<InEdge> rows;
    for (int i = 0; i < m; ++i) {
      const auto u = static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n));
      auto w = static_cast<NodeId>(rng() % static_cast<std::uint64_t>(n - 1));
      if (w >= u) ++w;
      rows.push_back({u, w, static_cast<Amount>(rng() % 7), static_cast<Cost>(1 + rng() % 9) * factor});
    }
    auto g = shared(make_graph(n, rows));
    std::vector<NodeId> srcs = {0}, dsts = {n - 1};
    if (rng() % 2 == 0) srcs.push_back(1);
    if (rng() % 2 == 0) dsts.push_back(n - 2);
    std::optional<Eligibility> el;
    switch (rng() % 3) {
      case 0: el = g->eligibility(std::vector<NodeId> {below(n)}, {}); break;
      case 1: el = g->eligibility({}, std::vector<EdgeId> {below(m)}); break;
      default: break;
    }
    const Eligibility* elig = el ? &*el : nullptr;
    FlowLedger fg(g);
    std::optional<FlowPolicy> other;
    if (n >= 6 && rng() % 2 == 0) {
      other.emplace(fg, v({2}), v({3}), cost_only_config());
      (void)other->place(static_cast<Amount>(rng() % 5));
    }
    FlowPolicy p(fg, srcs, dsts, config, elig);
    const std::vector<std::optional<Amount>> calls = {static_cast<Amount>(rng() % 6), static_cast<Amount>(rng() % 6),
                                                      std::nullopt};
    for (std::size_t call = 0; call < calls.size(); ++call) {
      const std::string where = "trial " + std::to_string(trial) + " call " + std::to_string(call);
      (void)p.place(calls[call]);
      expect_balanced(fg, where);
      Amount delivered = 0;
      std::set<NodeId> from, into;
      for (const auto& r : p.records()) {
        const std::string at = where + " record " + std::to_string(r.flow_id.second);
        for (NodeId x : r.sources) {
          EXPECT_TRUE(std::find(srcs.begin(), srcs.end(), x) != srcs.end()) << at << " source " << x;
        }
        for (NodeId x : r.destinations) {
          EXPECT_TRUE(std::find(dsts.begin(), dsts.end(), x) != dsts.end()) << at << " destination " << x;
        }
        const bool carries = expect_conserved(*g, r, at);
        for (const auto& [e, x] : r.edges) {
          EXPECT_TRUE(elig == nullptr || elig->edge_ok(e)) << at << " carries " << x << " on excluded edge " << e;
        }
        for (const auto& [e, x] : r.drops) {
          EXPECT_TRUE(elig == nullptr || elig->edge_ok(e)) << at << " drops on excluded edge " << e;
        }
        delivered += r.delivered;
        if (r.delivered > 0) {
          from.insert(r.sources.begin(), r.sources.end());
          for (const auto& [e, x] : r.edges) {
            const NodeId head = g->edge_dst()[static_cast<std::size_t>(e)];
            if (std::find(dsts.begin(), dsts.end(), head) != dsts.end()) into.insert(head);
          }
        }
        (cost_only ? carrying_cost_only : carrying_capacity_aware) += carries ? 1 : 0;
        dropping += r.drops.empty() ? 0 : 1;
        masked += carries && elig != nullptr ? 1 : 0;
        scaled += carries && factor > 1 ? 1 : 0;
      }
      if (config.admission == Admission::Lossless) EXPECT_EQ(delivered, p.placed_demand()) << where;
      several += from.size() > 1 || into.size() > 1 ? 1 : 0;
      ASSERT_FALSE(::testing::Test::HasFailure()) << where;
    }
  }
  EXPECT_GT(carrying_cost_only, 0);
  EXPECT_GT(carrying_capacity_aware, 0);
  EXPECT_GT(dropping, 0);
  EXPECT_GT(masked, 0);
  EXPECT_GT(scaled, 0);
  EXPECT_GT(several, 0);
}

TEST(FlowPolicyRecords, BundlesAreSharedOrPruned) {
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
  EXPECT_EQ(ascending_edges(*rec.bundle), (std::vector<EdgeId> {1, 3}));
  EXPECT_EQ(q.usable_bundles(), 1);
  // CostOnly: a source whose route is the whole cached support shares it.
  DagCache cache(*g);
  FlowLedger fg3(g);
  FlowPolicy c(fg3, v({0}), v({3}), cost_only_config(), nullptr, &cache);
  (void)c.place(2);
  const auto support = cache.lookup(*g, v({3}), Direction::Reverse, nullptr);
  EXPECT_EQ(c.records().front().bundle.get(), support.get());
}

// records() lists the live records in FlowId order, also when the sources' flows interleave.
// Sources 0 and 1, each with two parallel edges of 3 and 5 units to node 2; CapacityAware single
// paths, up to two flows per source. place(2) offers 1 to each source: flows 0 (source 0) and 1
// (source 1). place(14) offers 7 to each: each source fills its flow to 3 and grows a second flow
// on its other edge (5 units), FlowIds 2 (source 0) and 3 (source 1). Per source the flows are
// [0, 2] and [1, 3]; records() gives 0, 1, 2, 3.
TEST(FlowPolicyRecords, RecordsAreInFlowIdOrderAcrossSources) {
  auto g = shared(make_graph(3, {{0, 2, 3, 1}, {0, 2, 5, 1}, {1, 2, 3, 1}, {1, 2, 5, 1}}));
  FlowLedger fg(g);
  CapacityAware ca;
  ca.max_flow_count = 2;
  FlowPolicy p(fg, v({0, 1}), v({2}), cfg(ca, FlowPlacement::MaxFlow, Admission::Lossless,
                                            sel(EdgeTieBreak::Deterministic, EdgePriority::Capacity, false)));
  ASSERT_EQ(p.place(2).placed, 2);
  ASSERT_EQ(p.place(14).placed, 14);
  std::vector<std::int64_t> ids;
  std::vector<NodeId> sources;
  for (const auto& r : p.records()) {
    ids.push_back(r.flow_id.second);
    sources.push_back(r.sources.at(0));
  }
  EXPECT_EQ(ids, (std::vector<std::int64_t> {0, 1, 2, 3}));
  EXPECT_EQ(sources, (std::vector<NodeId> {0, 1, 0, 1}));
}

// CostOnly and CapacityAware route to the nearest destinations, and one flow may deliver to
// several of them: edges (input order) 0 -> 1 (10 units, cost 1), 1 -> 2 (2, cost 1), 1 -> 3 (3,
// cost 1) and 0 -> 3 (5, cost 5); destinations 2 and 3, both at cost 2 through node 1. The one
// flow carries 2 to node 2 and 3 to node 3 (the edges entering each destination), 5 in all; the
// costlier edge 0 -> 3 is no part of the route.
TEST(FlowPolicyRecords, OneFlowDeliversToTheNearestDestinations) {
  auto g = shared(make_graph(4, {{0, 1, 10, 1}, {1, 2, 2, 1}, {1, 3, 3, 1}, {0, 3, 5, 5}}));
  const auto id = id_of_input(*g);
  CapacityAware fixed;
  fixed.max_flow_count = 1;
  for (const auto& [name, routing] : std::vector<std::pair<std::string, Routing>> {{"CostOnly", CostOnly {}},
                                                                                   {"CapacityAware", fixed}}) {
    FlowLedger fg(g);
    FlowPolicy p(fg, v({0}), v({2, 3}), cfg(routing, FlowPlacement::MaxFlow));
    const auto r = p.place(std::nullopt);
    EXPECT_EQ(r.placed, 5) << name;
    const auto records = p.records();
    ASSERT_EQ(records.size(), 1U) << name;
    const auto& rec = records[0];
    EXPECT_EQ(rec.delivered, 5) << name;
    EXPECT_EQ(rec.route_cost, 2) << name;
    EXPECT_EQ(rec.destinations, (std::vector<NodeId> {2, 3})) << name;
    std::map<NodeId, Amount> into;
    for (const auto& [e, x] : rec.edges) into[g->edge_dst()[static_cast<std::size_t>(e)]] += x;
    EXPECT_EQ(into, (std::map<NodeId, Amount> {{1, 5}, {2, 2}, {3, 3}})) << name;
    EXPECT_EQ(ascending_edges(*rec.bundle), (std::vector<EdgeId> {id[0], id[1], id[2]})) << name;
  }
}

// The call's termination is the highest local outcome among the unserved sources, LimitFlows >
// LimitDiminishing > Exhausted > NoRoute, here around LimitDiminishing. Sources 0 and 1 with
// parallel edges to node 2; CapacityAware single paths (the lowest EdgeId first), MaxFlow, one flow
// at creation, the guard DiminishingReturns(1, 0.5); request 40 offers 20 to each source, so a
// growth step that gains below 10 stops growth. Source 1 has edges of 5, 3 and 10 units: it creates
// 5, grows 3 and stops (LimitDiminishing, 8 placed). Source 0:
//   edges 5, 10, 10 and at most 2 flows: 5, then 10 (not below 10), then the ceiling of 2 flows
//   above flow_count 1 with the offer unmet: LimitFlows (15), and the call is LimitFlows;
//   edges 5 and 10, no ceiling: 5, then 10, then no route is left: Exhausted (15), and the call
//   is LimitDiminishing;
//   one edge of 5 that a third party fills: no flow (NoRoute), and the call is LimitDiminishing.
TEST(FlowPolicyOutcomes, TheCallTakesTheHighestUnservedOutcome) {
  using Outcomes = std::vector<std::pair<NodeId, Termination>>;
  struct Case {
    const char* name;
    std::vector<Amount> source0;
    std::optional<std::int32_t> max_flow_count;
    bool occupied;
    Amount placed;
    Termination first, call;
  };
  const std::vector<Case> cases = {
      {"LimitFlows over LimitDiminishing", {5, 10, 10}, 2, false, 23, Termination::LimitFlows, Termination::LimitFlows},
      {"LimitDiminishing over Exhausted", {5, 10}, std::nullopt, false, 23, Termination::Exhausted,
       Termination::LimitDiminishing},
      {"LimitDiminishing over NoRoute", {5}, std::nullopt, true, 8, Termination::NoRoute, Termination::LimitDiminishing},
  };
  for (const auto& c : cases) {
    std::vector<InEdge> rows;
    for (const Amount x : c.source0) rows.push_back({0, 2, x, 1});
    for (const Amount x : {5, 3, 10}) rows.push_back({1, 2, x, 1});
    auto g = shared(make_graph(3, rows));
    FlowLedger fg(g);
    std::optional<FlowPolicy> third;
    if (c.occupied) {
      Pinned edge;
      edge.bundles = {dag_of(*g, {0})};
      third.emplace(fg, v({0}), v({2}), cfg(edge, FlowPlacement::MaxFlow));
      ASSERT_EQ(third->place(5).placed, 5) << c.name;
    }
    CapacityAware ca;
    ca.max_flow_count = c.max_flow_count;
    ca.diminishing = DiminishingReturns {1, 0.5};
    FlowPolicy p(fg, v({0, 1}), v({2}), cfg(ca, FlowPlacement::MaxFlow, Admission::Lossless,
                                              sel(EdgeTieBreak::Deterministic, EdgePriority::Capacity, false)));
    const auto r = p.place(40);
    EXPECT_EQ(r.placed, c.placed) << c.name;
    EXPECT_EQ(r.local_outcomes, (Outcomes {{0, c.first}, {1, Termination::LimitDiminishing}})) << c.name;
    EXPECT_EQ(r.termination, c.call) << c.name;
  }
}

namespace {

// Three disjoint routes from 0 to 4 (input order): A = 0 -> 1 -> 4 of cost 1 + 1 = 2, B = 0 -> 2
// -> 4 of cost 2 + 1 = 3, C = 0 -> 3 -> 4 of cost 2 + 2 = 4, capacities a, b and c on both
// edges of each route.
std::shared_ptr<const Graph> three_routes(Amount a, Amount b, Amount c) {
  return shared(make_graph(5, {{0, 1, a, 1}, {1, 4, a, 1}, {0, 2, b, 2}, {2, 4, b, 1}, {0, 3, c, 2}, {3, 4, c, 2}}));
}

// A policy's records as (route cost, delivered, ascending bundle edges, per-edge allocation).
using RouteRow = std::tuple<Cost, Amount, std::vector<EdgeId>, std::vector<std::pair<EdgeId, Amount>>>;

std::vector<RouteRow> routes_of(const FlowPolicy& p) {
  std::vector<RouteRow> out;
  for (const auto& r : p.records()) out.emplace_back(r.route_cost, r.delivered, ascending_edges(*r.bundle), r.edges);
  return out;
}

// The route of the input edges i and j carrying x units on each.
RouteRow route(const Graph& g, Cost cost, Amount x, std::size_t i, std::size_t j) {
  const auto inv = id_of_input(g);
  std::vector<EdgeId> e = {inv[i], inv[j]};
  std::sort(e.begin(), e.end());
  return {cost, x, e, {{e[0], x}, {e[1], x}}};
}

CapacityAware ceiling(std::optional<CostCeiling> c, Amount floor = 0) {
  CapacityAware ca;
  ca.max_path_cost = c;
  ca.min_flow_bandwidth = floor;
  return ca;
}

} // namespace

// The CapacityAware path-cost ceiling (inclusive, Absolute(c) a budget, Stretch(f) f times the
// static cost-only shortest cost), MaxFlow and EqualBalanced, one source, place(10) with
// flow_count 1 and no max_flow_count. Derivation (the stages of flow_policy.cpp: 1 creation,
// 3 amounts, 4 growth) on three_routes(2, 2, 2), A, B, C of cost 2, 3, 4 and 2 units each:
//   MaxFlow: stage 1 creates a flow on the shortest admitted route A; stage 3 fills it to
//   T_max = 2; growth searches W (A now has residual 0) and adds B (cost 3) filled to 2, then
//   searches again: only C is left. Under Absolute(3) or Stretch(1.5) (1.5 * 2 = 3) C's cost 4
//   is above the ceiling: no candidate, stop, placed 4, `exhausted` (growth found no
//   candidate). Under Absolute(4) or Stretch(2) (ceiling 4, inclusive) C is added: placed 6.
//   EqualBalanced: stage 3 sets T = T_max_cohort({A}) = 2; growth adds B with T' = min(0 + 10,
//   T_max_cohort({A, B})) = 4 split [2, 2] (accepted, gain 2); the next search finds only C:
//   the same ceilings give placed 4, or 6 with C (cohort [2, 2, 2]).
//   Absolute(1): no route is admitted at all, stage 1 creates nothing: placed 0, `no_route`.
TEST(CapacityAwareGates, PathCostCeilingExcludesACostlierRoute) {
  auto g = three_routes(2, 2, 2);
  const auto A = route(*g, 2, 2, 0, 1);
  const auto B = route(*g, 3, 2, 2, 3);
  const auto C = route(*g, 4, 2, 4, 5);
  for (const auto kernel : {FlowPlacement::MaxFlow, FlowPlacement::EqualBalanced}) {
    struct Case {
      const char* name;
      std::optional<CostCeiling> c;
      Amount placed;
      Termination term;
      std::vector<RouteRow> routes;
    };
    const std::vector<Case> cases = {
        {"Absolute(3)", Absolute {3}, 4, Termination::Exhausted, {A, B}},
        {"Stretch(1.5)", Stretch {1.5}, 4, Termination::Exhausted, {A, B}},
        {"Absolute(4)", Absolute {4}, 6, Termination::Exhausted, {A, B, C}},
        {"Stretch(2)", Stretch {2.0}, 6, Termination::Exhausted, {A, B, C}},
        {"none", std::nullopt, 6, Termination::Exhausted, {A, B, C}},
        {"Absolute(1)", Absolute {1}, 0, Termination::NoRoute, {}},
    };
    for (const auto& c : cases) {
      const std::string where =
          std::string(kernel == FlowPlacement::MaxFlow ? "MaxFlow" : "EqualBalanced") + " ceiling " + c.name;
      FlowLedger fg(g);
      FlowPolicy p(fg, v({0}), v({4}), cfg(ceiling(c.c), kernel));
      const auto r = p.place(10);
      EXPECT_EQ(r.placed, c.placed) << where;
      EXPECT_EQ(r.remaining, 10 - c.placed) << where;
      EXPECT_EQ(r.termination, c.term) << where;
      EXPECT_EQ(r.local_outcomes, (std::vector<std::pair<NodeId, Termination>> {{0, c.term}})) << where;
      EXPECT_EQ(r.routes_changed, !c.routes.empty()) << where;
      EXPECT_FALSE(r.weights_changed) << where;
      EXPECT_EQ(routes_of(p), c.routes) << where;
    }
  }
}

// The Stretch baseline is the static shortest cost, residual ignored. A third party holds
// both units of A before the policy is built: the static shortest cost stays 2, so Stretch(1.5)
// admits cost 3 and not 4. Stage 1 creates the flow on B (A has residual 0), stage 3 fills it to
// 2, growth finds no admitted route: placed 2, `exhausted`. A baseline taken on the residual
// (3, from B) would admit C under 1.5 * 3 = 4.5 and place 4; Stretch(2) admits C (2 * 2 = 4)
// and places 4.
TEST(CapacityAwareGates, StretchBaselineIsTheStaticShortestCost) {
  auto g = three_routes(2, 2, 2);
  const auto inv = id_of_input(*g);
  const auto B = route(*g, 3, 2, 2, 3);
  const auto C = route(*g, 4, 2, 4, 5);
  for (const auto kernel : {FlowPlacement::MaxFlow, FlowPlacement::EqualBalanced}) {
    for (const double f : {1.5, 2.0}) {
      const std::string where = std::string(kernel == FlowPlacement::MaxFlow ? "MaxFlow" : "EqualBalanced") +
                                " Stretch " + std::to_string(f);
      FlowLedger fg(g);
      Pinned third;
      third.bundles = {dag_of(*g, {inv[0], inv[1]})};
      FlowPolicy other(fg, v({0}), v({4}), cfg(third, FlowPlacement::MaxFlow));
      ASSERT_EQ(other.place(2).placed, 2) << where;
      FlowPolicy p(fg, v({0}), v({4}), cfg(ceiling(Stretch {f}), kernel));
      const auto r = p.place(10);
      const bool with_c = f == 2.0;
      EXPECT_EQ(r.placed, with_c ? 4 : 2) << where;
      EXPECT_EQ(r.termination, Termination::Exhausted) << where;
      EXPECT_EQ(routes_of(p), with_c ? (std::vector<RouteRow> {B, C}) : (std::vector<RouteRow> {B})) << where;
    }
  }
  // Under an Eligibility the baseline is the shortest eligible cost: with A's first edge excluded
  // it is B's 3, so Stretch(1.5) admits up to 4.5, C included: B and C, 4 placed (a baseline that
  // ignored the Eligibility, 2, would admit up to 3: B alone).
  const auto el = g->eligibility({}, std::vector<EdgeId> {inv[0]});
  for (const auto kernel : {FlowPlacement::MaxFlow, FlowPlacement::EqualBalanced}) {
    FlowLedger fg(g);
    FlowPolicy p(fg, v({0}), v({4}), cfg(ceiling(Stretch {1.5}), kernel), &el);
    const auto r = p.place(10);
    EXPECT_EQ(r.placed, 4);
    EXPECT_EQ(r.termination, Termination::Exhausted);
    EXPECT_EQ(routes_of(p), (std::vector<RouteRow> {B, C}));
  }
}

// The policy's Stretch ceiling is exact at large costs and factors: route A 0-1-4 of cost
// 2^60 and route B 0-2-4 of cost 2^60 + 1, 2 units each, MaxFlow, place(10). Stretch(1) admits
// A only (in a double both costs round to 2^60): placed 2. Stretch(1 + 2^-52) admits B
// (2^60 + 2^8 >= 2^60 + 1), and so does Stretch(2^70), whose product with 2^60 needs more than
// 128 bits unless the bit lengths are compared first: placed 4.
TEST(CapacityAwareGates, StretchCeilingIsExactForLargeCostsAndFactors) {
  const Cost half = Cost {1} << 59;
  auto g = shared(make_graph(5, {{0, 1, 2, half}, {1, 4, 2, half}, {0, 2, 2, half}, {2, 4, 2, half + 1}}));
  const auto A = route(*g, 2 * half, 2, 0, 1);
  const auto B = route(*g, 2 * half + 1, 2, 2, 3);
  for (const double f : {1.0, std::nextafter(1.0, 2.0), std::ldexp(1.0, 70)}) {
    const std::string where = "Stretch " + std::to_string(f);
    FlowLedger fg(g);
    FlowPolicy p(fg, v({0}), v({4}), cfg(ceiling(Stretch {f}), FlowPlacement::MaxFlow));
    const auto r = p.place(10);
    const bool with_b = f != 1.0;
    EXPECT_EQ(r.placed, with_b ? 4 : 2) << where;
    EXPECT_EQ(r.termination, Termination::Exhausted) << where;
    EXPECT_EQ(routes_of(p), with_b ? (std::vector<RouteRow> {A, B}) : (std::vector<RouteRow> {A})) << where;
  }
}

// min_flow_bandwidth gates every CapacityAware search at residual >= max(1, floor), on
// three_routes(4, 2, 0): A of 4 units at cost 2, B of 2 units at cost 3, C of 0 units, place(10).
//   Floor 3, MaxFlow: stage 1 creates on A (residual 4 >= 3), stage 3 fills 4; growth: A has
//   residual 0 and B 2 < 3, no route: stop, placed 4, `exhausted`. Floor 2: growth adds B (2):
//   placed 6. EqualBalanced: T = T_max_cohort({A}) = 4; floor 3 stops growth there (placed 4);
//   floor 2 adds B with T' = T_max_cohort({A, B}) = 5 split [3, 2] by round robin in FlowId
//   order (6 would need 3 on B), gain 1; the next search sees A at residual 1 < 2: placed 5.
//   Floor 5 is above every residual: stage 1 finds no route, placed 0, `no_route`.
TEST(CapacityAwareGates, MinFlowBandwidthStopsGrowthBelowTheFloor) {
  auto g = three_routes(4, 2, 0);
  struct Case {
    FlowPlacement kernel;
    Amount floor;
    Amount placed;
    Termination term;
    std::vector<RouteRow> routes;
  };
  const std::vector<Case> cases = {
      {FlowPlacement::MaxFlow, 3, 4, Termination::Exhausted, {route(*g, 2, 4, 0, 1)}},
      {FlowPlacement::MaxFlow, 2, 6, Termination::Exhausted, {route(*g, 2, 4, 0, 1), route(*g, 3, 2, 2, 3)}},
      {FlowPlacement::MaxFlow, 5, 0, Termination::NoRoute, {}},
      {FlowPlacement::EqualBalanced, 3, 4, Termination::Exhausted, {route(*g, 2, 4, 0, 1)}},
      {FlowPlacement::EqualBalanced, 2, 5, Termination::Exhausted, {route(*g, 2, 3, 0, 1), route(*g, 3, 2, 2, 3)}},
      {FlowPlacement::EqualBalanced, 5, 0, Termination::NoRoute, {}},
  };
  for (const auto& c : cases) {
    const std::string where = std::string(c.kernel == FlowPlacement::MaxFlow ? "MaxFlow" : "EqualBalanced") +
                              " floor " + std::to_string(c.floor);
    FlowLedger fg(g);
    FlowPolicy p(fg, v({0}), v({4}), cfg(ceiling(std::nullopt, c.floor), c.kernel));
    const auto r = p.place(10);
    EXPECT_EQ(r.placed, c.placed) << where;
    EXPECT_EQ(r.termination, c.term) << where;
    EXPECT_EQ(r.routes_changed, !c.routes.empty()) << where;
    EXPECT_EQ(routes_of(p), c.routes) << where;
  }
}

// The diminishing guard stops growth once at least `window` growth steps were accepted in the
// call and the last `window` gains sum strictly below fraction times the source's offer. One source,
// parallel edges 0 -> 1 taken in EdgeId order (CapacityAware single paths, MaxFlow, one flow at
// creation), DiminishingReturns(window, 0.5):
//   edges 2, 6, 1, 4, offer 12, window 1 (threshold 6): creation 2; growth gains 6 (equal to the
//   threshold, not below: continue), then 1 (below, with 3 of the offer left): stop,
//   LimitDiminishing, 9 placed. Only the last gain counts: 6 + 1 would not stop it.
//   edges 2, 1, 4, offer 8, window 1 (threshold 4): creation 2; the first growth step gains 1:
//   stop at once, 3 placed.
//   edges 2, 1, 4, offer 8, window 2: the first step is not checked alone; 1 + 4 = 5 is not
//   below 4; no edge is left: Exhausted, 7 placed.
TEST(CapacityAwareGates, DiminishingStopsWhenTheLastWindowGainsFallBelowTheFraction) {
  struct Case {
    std::vector<Amount> caps;
    Amount offer;
    std::int32_t window;
    Amount placed;
    Termination term;
    std::vector<Amount> flows;
  };
  const std::vector<Case> cases = {
      {{2, 6, 1, 4}, 12, 1, 9, Termination::LimitDiminishing, {2, 6, 1}},
      {{2, 1, 4}, 8, 1, 3, Termination::LimitDiminishing, {2, 1}},
      {{2, 1, 4}, 8, 2, 7, Termination::Exhausted, {2, 1, 4}},
  };
  for (const auto& c : cases) {
    const std::string where = "window " + std::to_string(c.window) + " offer " + std::to_string(c.offer);
    std::vector<InEdge> rows;
    for (const Amount x : c.caps) rows.push_back({0, 1, x, 1});
    auto g = shared(make_graph(2, rows));
    FlowLedger fg(g);
    CapacityAware ca;
    ca.diminishing = DiminishingReturns {c.window, 0.5};
    FlowPolicy p(fg, v({0}), v({1}), cfg(ca, FlowPlacement::MaxFlow, Admission::Lossless,
                                         sel(EdgeTieBreak::Deterministic, EdgePriority::Capacity, false)));
    const auto r = p.place(c.offer);
    EXPECT_EQ(r.placed, c.placed) << where;
    EXPECT_EQ(r.termination, c.term) << where;
    std::vector<Amount> flows;
    for (const auto& rec : p.records()) flows.push_back(rec.delivered);
    EXPECT_EQ(flows, c.flows) << where;
  }
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
  // The record bundles the engine builds itself are sparse too: a supplied bundle pruned by an
  // Eligibility (the bypass n - 3 -> n - 1 excluded), and the CostOnly and Global routes.
  auto h = shared(make_graph(n, {{n - 3, n - 2, 4, 1}, {n - 2, n - 1, 4, 1}, {n - 3, n - 1, 4, 5}}));
  const auto inv = id_of_input(*h);
  Pinned all;
  all.bundles = {dag_of(*h, {inv[0], inv[1], inv[2]})};
  const auto el = h->eligibility({}, std::vector<EdgeId> {inv[2]});
  const std::vector<std::pair<const char*, FlowPolicyConfig>> built = {
      {"pruned Pinned", cfg(all, FlowPlacement::MaxFlow)}, {"CostOnly", cost_only_config()}, {"Global", global_config()}};
  for (const auto& [name, c] : built) {
    FlowLedger fh(h);
    FlowPolicy q(fh, v({n - 3}), v({n - 1}), c, std::holds_alternative<Pinned>(c.routing) ? &el : nullptr);
    EXPECT_EQ(q.place(3).placed, 3) << name;
    const auto bundle = q.records().front().bundle;
    EXPECT_NE(bundle.get(), all.bundles[0].get()) << name;
    EXPECT_EQ(ascending_edges(*bundle), (std::vector<EdgeId> {inv[0], inv[1]})) << name;
    EXPECT_LT(detail::PredDAGBuilder::payload_bytes(*bundle), 128U) << name;
  }
}

TEST(Evaluators, CapacityIsPlaceInfOnAFreshLedger) {
  // Every field of the result, in every routing model (Pinned's bundle joins the endpoints).
  auto g = diamond(5, 3, 2);
  Pinned pin;
  pin.bundles = {dag_of(*g, {0, 2}), dag_of(*g, {1, 3})};
  for (const auto& c : {global_config(), cost_only_config(),
                        cfg(CapacityAware {.flow_count = 2, .max_flow_count = 4}, FlowPlacement::MaxFlow),
                        cfg(pin, FlowPlacement::EqualBalanced)}) {
    const auto r = capacity(*g, v({0}), v({3}), CapacityOptions {c, nullptr, nullptr, OutputOptions {true, true, true}});
    FlowLedger fg(g);
    FlowPolicy p(fg, v({0}), v({3}), c);
    const auto want = p.place(std::nullopt, OutputOptions {true, true, true});
    EXPECT_GT(want.placed, 0);  // two empty results would agree too
    expect_same(r, want, "capacity against place(+inf)");
  }
  // The exact maximum flow: Global's curve (cheapest walks first).
  const auto exact = capacity(*g, v({0}), v({3}), CapacityOptions {global_config()});
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
  const CapacityOptions opts {global_config(), &common, nullptr, OutputOptions {true, true, true}};
  const auto r1 = batch_capacity(*g, pairs, opts, per, 1);
  ASSERT_EQ(r1.size(), 4U);
  EXPECT_EQ(r1[0].placed, 5);  // common only: route 1
  EXPECT_EQ(r1[1].placed, 0);  // common and the pair's: both routes cut
  EXPECT_EQ(r1[2].placed, 5);
  EXPECT_EQ(r1[3].placed, 8);  // destination set {1, 2}: both first edges (edge 3 excluded is past 2)
  EXPECT_EQ(*r1[3].cost_distribution, curve({{1, 8}}));
  EXPECT_EQ(*r1[3].min_cut, (std::vector<EdgeId> {0, 1}));
  // Every field of every result is the same for any worker count.
  for (std::int32_t w : {2, 3, 8}) {
    const auto rw = batch_capacity(*g, pairs, opts, per, w);
    ASSERT_EQ(rw.size(), 4U);
    for (std::size_t i = 0; i < 4; ++i) {
      expect_same(rw[i], r1[i], "workers " + std::to_string(w) + " pair " + std::to_string(i));
    }
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
  const CapacityOptions with_cache {cost_only_config(), nullptr, &cache};
  EXPECT_TRUE(contains(error_of([&] { (void)batch_capacity(*g, pairs, with_cache, {}, 2); }), "workers == 1"));
  const auto cached = batch_capacity(*g, pairs, with_cache, {}, 1);
  EXPECT_EQ(cache.stats().misses, 2);  // destination sets {3} and {1, 2}
  EXPECT_EQ(cache.stats().hits, 2);
  EXPECT_EQ(cached[0].placed, 7);  // CostOnly EB over routes of 5 and 3: [4, 3] (odd unit to the lower id)
  // The argument checks in their documented order, all before the pairs: with every defect
  // present the first is reported, and with that one removed the next. Pinned also carries a
  // configuration error (Lossy admission with the MaxFlow kernel); pair 2 has a destination out of
  // range. The building cache is untouched: no lookup ran.
  DagCache foreign_cache(*other);
  DagCache building(*g);
  Pinned lossy_pin;
  lossy_pin.bundles = {dag_of(*g, {0, 2})};
  const FlowPolicyConfig pinned_bad = cfg(lossy_pin, FlowPlacement::MaxFlow, Admission::Lossy);
  const FlowPolicyConfig config_bad = cfg(CostOnly {}, FlowPlacement::MaxFlow, Admission::Lossy);
  bad = pairs;
  bad[2].second = {9};
  const std::vector<const Eligibility*> short_per = {nullptr};
  auto first_error = [&](const FlowPolicyConfig& c, const Eligibility* elig, DagCache* with,
                         const std::vector<const Eligibility*>& per_pair, std::int32_t workers) {
    return error_of([&] { (void)batch_capacity(*g, bad, CapacityOptions {c, elig, with}, per_pair, workers); });
  };
  EXPECT_EQ(first_error(pinned_bad, &foreign, &foreign_cache, short_per, 0), "workers must be >= 1 (got 0)");
  EXPECT_EQ(first_error(pinned_bad, &foreign, &foreign_cache, short_per, 2),
            "foreign Eligibility: built for another graph");
  EXPECT_EQ(first_error(pinned_bad, nullptr, &foreign_cache, short_per, 2), "foreign DagCache: bound to another graph");
  EXPECT_EQ(first_error(pinned_bad, nullptr, &building, short_per, 2),
            "a supplied building DagCache has one writer: it is rejected unless workers == 1 (got workers = 2; a "
            "frozen cache may be shared by every worker)");
  EXPECT_EQ(first_error(pinned_bad, nullptr, &building, short_per, 1),
            "Pinned routing is rejected in batch_capacity (bundles bind their endpoints)");
  EXPECT_EQ(first_error(config_bad, nullptr, &building, short_per, 1), kLossyKernelError);
  EXPECT_EQ(first_error(global_config(), nullptr, &building, short_per, 1), "per_pair_elig has 1 entries for 4 pairs");
  EXPECT_EQ(first_error(global_config(), nullptr, &building, per, 1), "destinations: node 9 out of node range [0, 4)");
  EXPECT_EQ(stats_of(building.stats()), stats_of(DagCache(*g).stats()));
}

TEST(Evaluators, BatchValidatesInPairOrderBeforeAnyWorker) {
  // Two invalid pairs (1 and 3): the batch raises pair 1's error before any worker starts.
  auto g = diamond(5, 3);
  std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {{0}, {3}}, {{0}, {7}}, {{0}, {3}}, {{0}, {8}}};
  const auto msg = error_of([&] { (void)batch_capacity(*g, pairs, CapacityOptions {global_config()}, {}, 4); });
  EXPECT_TRUE(contains(msg, "node 7")) << msg;
  // Within one pair: the sources, then the destinations, then their overlap, then the pair's
  // Eligibility's graph.
  const auto elsewhere = diamond(5, 3);
  const auto foreign = elsewhere->eligibility({}, std::vector<EdgeId> {0});
  const std::vector<const Eligibility*> per = {&foreign};
  const auto pair_error = [&](std::vector<NodeId> s, std::vector<NodeId> d) {
    const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> one = {{s, d}};
    return error_of([&] { (void)batch_capacity(*g, one, CapacityOptions {global_config()}, per, 1); });
  };
  EXPECT_EQ(pair_error({}, {9}), "sources: the node set is empty");
  EXPECT_EQ(pair_error({0, 3}, {3, 9}), "destinations: node 9 out of node range [0, 4)");
  EXPECT_EQ(pair_error({0, 3}, {3}), "node 3 is in both endpoint sets: sources and destinations must not overlap");
  EXPECT_EQ(pair_error({0}, {3}), "foreign Eligibility: per_pair_elig[0] built for another graph");
}

// The sensitivity shortcuts (Global from the baseline's maximum flow, CostOnly without a
// re-evaluation of the edges outside the static shortest-path support reachable from the sources)
// and the plain re-evaluations (CapacityAware, Pinned) give the losses of the definition, the
// baseline's placed amount minus capacity() with the edge excluded, on random graphs with parallel
// edges, self-loops, dispersed costs, zero capacities, endpoint sets and eligibility masks; every
// CostOnly kernel with every tie break, multipath and multi_edge setting under the Capacity
// priority, both admissions, and one Residual single-path selection (no Headroom); two
// CapacityAware configurations (growth, and reoptimization under a Stretch ceiling); and Pinned
// with two bundles (the static support toward the destination and the second shortest path) on
// single endpoints. An edge without an entry has loss 0 by the same definition.
TEST(Evaluators, SensitivityEqualsTheDefinitionOnRandomGraphs) {
  std::mt19937 rng(20261005);
  auto pick = [&](int n) { return static_cast<int>(rng() % static_cast<unsigned>(n)); };
  std::vector<FlowPolicyConfig> configs = {global_config()};
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
  configs.push_back(cfg(CapacityAware {.flow_count = 1, .max_flow_count = 3}, FlowPlacement::EqualBalanced,
                        Admission::Lossless, sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, false, false)));
  configs.push_back(cfg(CapacityAware {.flow_count = 2, .max_path_cost = Stretch {2.0}, .reoptimize = true},
                        FlowPlacement::MaxFlow));
  int entries = 0;
  int global_losses = 0;
  int policy_losses = 0;
  int policy_zero = 0;
  int aware_losses = 0;
  int pinned_losses = 0;
  int pinned_trials = 0;
  int absent = 0;
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
    // Pinned on the single endpoints 0 and n - 1, when 0 reaches n - 1 statically.
    auto trial_configs = configs;
    const NodeId t0 = n - 1;
    const auto support = std::make_shared<const PredDAG>(
        detail::static_support(*g, std::vector<NodeId> {t0}, Direction::Reverse, nullptr));
    const auto paths = ksp(*g, 0, t0, KspOptions(2, sel()));
    if (paths.size() == 2) {
      Pinned pin;
      pin.bundles = {support, std::make_shared<const PredDAG>(paths[1].second.to_dag())};
      trial_configs.push_back(cfg(pin, FlowPlacement::Weighted, Admission::Lossless, sel()));
      ++pinned_trials;
    }
    for (const auto& c : trial_configs) {
      const bool global = std::holds_alternative<Global>(c.routing);
      const bool cost_only = std::holds_alternative<CostOnly>(c.routing);
      const bool pinned = std::holds_alternative<Pinned>(c.routing);
      const std::vector<NodeId> cs = pinned ? std::vector<NodeId> {0} : s;
      const std::vector<NodeId> ct = pinned ? std::vector<NodeId> {t0} : t;
      const auto r = sensitivity(*g, cs, ct, SensitivityOptions {c, elig}, 1 + pick(3));
      // The baseline is capacity() itself, and the entries are the candidates in ascending order:
      // Global, the baseline's positive-flow edges; CostOnly, every allowed edge of at least one
      // unit of capacity.
      const auto base = capacity(*g, cs, ct, CapacityOptions {c, elig, nullptr, OutputOptions {true, false, false}});
      ASSERT_EQ(r.baseline.placed, base.placed) << "trial " << trial;
      std::vector<EdgeId> candidates;
      for (EdgeId e = 0; e < g->num_edges(); ++e) {
        const auto ei = static_cast<std::size_t>(e);
        if (global ? (*base.edge_flows)[ei] > 0 : (elig == nullptr || elig->edge_ok(e)) && g->capacity_units()[ei] >= 1) {
          candidates.push_back(e);
        }
      }
      const std::string model = global ? " Global" : cost_only ? " CostOnly" : pinned ? " Pinned" : " CapacityAware";
      std::vector<EdgeId> listed;
      for (const auto& e : r.entries) listed.push_back(e.edge);
      ASSERT_EQ(listed, candidates) << "trial " << trial << model;
      // The loss of deleting `e`, by the definition.
      const auto loss = [&](EdgeId e) {
        const EdgeId one[1] = {e};
        auto deleted = g->eligibility({}, one);
        if (elig != nullptr) deleted = Eligibility::intersect(*elig, deleted);
        return base.placed - capacity(*g, cs, ct, CapacityOptions {c, &deleted}).placed;
      };
      for (const auto& e : r.entries) {
        const Amount want = loss(e.edge);
        ASSERT_EQ(e.loss, want) << "trial " << trial << " edge " << e.edge << model;
        ++entries;
        if (want != 0) ++(global ? global_losses : cost_only ? policy_losses : pinned ? pinned_losses : aware_losses);
        if (cost_only && want == 0) ++policy_zero;
      }
      for (EdgeId e = 0; e < g->num_edges(); ++e) {
        if (std::find(listed.begin(), listed.end(), e) != listed.end()) continue;
        EXPECT_EQ(loss(e), 0) << "trial " << trial << " edge " << e << model << ": an edge without an entry";
        ++absent;
      }
    }
  }
  // Not vacuous: many entries, losses under every model, zero-loss CostOnly deletions, Pinned
  // trials, and edges without an entry.
  EXPECT_GT(entries, 3000);
  EXPECT_GT(global_losses, 100);
  EXPECT_GT(policy_losses, 500);
  EXPECT_GT(policy_zero, 500);
  EXPECT_GT(aware_losses, 100);
  EXPECT_GT(pinned_losses, 30);
  EXPECT_GT(pinned_trials, 20);
  EXPECT_GT(absent, 500);
}

TEST(Evaluators, SensitivityCandidatesAndSignedLosses) {
  // Global: candidates are the baseline's positive-flow edges; an unsaturated critical edge
  // appears: 0 -> 1 of capacity 9 carries 4 (its tail 1 -> 2 holds 4) and its deletion
  // loses 4; a candidate set of saturated edges only would miss it.
  auto g = shared(make_graph(3, {{0, 1, 9, 1}, {1, 2, 4, 1}, {0, 2, 2, 5}}));
  auto s = sensitivity(*g, v({0}), v({2}), SensitivityOptions {global_config()}, 1);
  EXPECT_EQ(s.baseline.placed, 6);
  EXPECT_FALSE(s.baseline.edge_flows.has_value());  // requested internally, not returned
  std::map<EdgeId, Amount> losses;
  for (const auto& e : s.entries) losses[e.edge] = e.loss;
  const auto inv = id_of_input(*g);
  EXPECT_EQ(losses.at(inv[0]), 4);  // unsaturated (4 of 9) and critical
  EXPECT_EQ(losses.at(inv[1]), 4);
  EXPECT_EQ(losses.at(inv[2]), 2);
  EXPECT_EQ(losses.size(), 3U);
  // A policy model: every allowed edge of >= 1 unit; a loss can be negative (EqualBalanced over
  // parallel links of 1 and 10 delivers 2; without the thin link it delivers 10).
  auto p = shared(make_graph(2, {{0, 1, 1, 1}, {0, 1, 10, 1}, {0, 1, 0, 1}}));
  const auto opts = SensitivityOptions {cost_only_config()};
  s = sensitivity(*p, v({0}), v({1}), opts, 1);
  EXPECT_EQ(s.baseline.placed, 2);
  // The exact (edge, loss) list in ascending EdgeId: the zero-capacity link is absent (loss zero).
  const auto pin = id_of_input(*p);
  std::vector<std::pair<EdgeId, Amount>> got;
  for (const auto& e : s.entries) got.emplace_back(e.edge, e.loss);
  std::vector<std::pair<EdgeId, Amount>> want = {{pin[0], -8}, {pin[1], 1}};
  std::sort(want.begin(), want.end());
  EXPECT_EQ(got, want);
  // Workers give the same entries; an excluded edge is no candidate.
  const auto s4 = sensitivity(*p, v({0}), v({1}), opts, 4);
  std::vector<std::pair<EdgeId, Amount>> got4;
  for (const auto& e : s4.entries) got4.emplace_back(e.edge, e.loss);
  EXPECT_EQ(got4, want);
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

// The baseline is capacity() with the requested outputs, in every routing model; without a
// request it carries none.
TEST(Evaluators, SensitivityBaselineIsCapacityWithTheRequestedOutputs) {
  auto g = diamond(5, 3, 2);
  Pinned pin;
  pin.bundles = {dag_of(*g, {0, 2}), dag_of(*g, {1, 3})};
  const OutputOptions all {true, true, true};
  for (const auto& c : {global_config(), cost_only_config(),
                        cfg(CapacityAware {.flow_count = 2, .max_flow_count = 4}, FlowPlacement::MaxFlow),
                        cfg(pin, FlowPlacement::EqualBalanced)}) {
    const auto s = sensitivity(*g, v({0}), v({3}), SensitivityOptions {c, nullptr, nullptr, all}, 1);
    const auto want = capacity(*g, v({0}), v({3}), CapacityOptions {c, nullptr, nullptr, all});
    EXPECT_GT(want.placed, 0);
    ASSERT_TRUE(want.edge_flows.has_value() && want.residual.has_value() && want.bottleneck.has_value());
    expect_same(s.baseline, want, "sensitivity baseline with every output");
    const auto plain = sensitivity(*g, v({0}), v({3}), SensitivityOptions {c}, 1);
    expect_same(plain.baseline, capacity(*g, v({0}), v({3}), CapacityOptions {c}), "sensitivity baseline, no output");
  }
}

// CostOnly re-evaluates a deletion exactly when the edge lies in the static shortest-path support
// (every edge on a shortest path to the destinations) reachable from a roster source; every other
// candidate gets loss 0 without a placement. A supplied building cache shows which deletions were
// placed: each placed deletion looks up its own fresh Eligibility uid (one miss and one entry),
// after the baseline's one. Edges 0->2, 1->2, 2->3, 2->4, 3->5, 4->5 of 10 units at cost 1: all six
// form the support toward 5. From source 0 the single path is 0->2->3->5, but 2->4 and 4->5, on no
// baseline route, are in its reachable support and are re-evaluated (loss 0); 1->2 is not and is
// skipped. From sources 0 and 1 the second source reaches 1->2 too.
TEST(Evaluators, CostOnlySensitivityReevaluatesExactlyTheReachableSupport) {
  auto g = shared(make_graph(6, {{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 10, 1}, {2, 4, 10, 1}, {3, 5, 10, 1}, {4, 5, 10, 1}}));
  const auto in = id_of_input(*g);
  const auto c = cfg(CostOnly {}, FlowPlacement::MaxFlow, Admission::Lossless,
                     sel(EdgeTieBreak::PreferHigherPriority, EdgePriority::Residual, false, true));
  const auto losses_of = [&](const SensitivityResult& r) {
    std::map<std::size_t, Amount> out;  // by input position
    for (const auto& e : r.entries) {
      out[static_cast<std::size_t>(std::find(in.begin(), in.end(), e.edge) - in.begin())] = e.loss;
    }
    return out;
  };
  DagCache one(*g);
  const auto r1 = sensitivity(*g, v({0}), v({5}), SensitivityOptions {c, nullptr, &one}, 1);
  EXPECT_EQ(r1.baseline.placed, 10);
  EXPECT_EQ(losses_of(r1), (std::map<std::size_t, Amount> {{0, 10}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}}));
  EXPECT_EQ(one.stats().hits, 0);
  EXPECT_EQ(one.stats().misses, 1 + 5);  // the baseline, then every candidate but 1->2
  EXPECT_EQ(one.stats().entries, 1 + 5);
  DagCache two(*g);
  const auto r2 = sensitivity(*g, v({0, 1}), v({5}), SensitivityOptions {c, nullptr, &two}, 1);
  EXPECT_EQ(r2.baseline.placed, 20);
  EXPECT_EQ(losses_of(r2), (std::map<std::size_t, Amount> {{0, 10}, {1, 10}, {2, 10}, {3, 10}, {4, 10}, {5, 10}}));
  EXPECT_EQ(two.stats().misses, 1 + 6);  // every candidate
  EXPECT_EQ(two.stats().entries, 1 + 6);
}

// sensitivity's worker and cache rules are batch_capacity's, checked in its order before any work:
// workers >= 1, the Eligibility's graph, the cache's graph, a supplied building cache only with
// one worker; a frozen cache is shared by every worker.
TEST(Evaluators, SensitivityChecksWorkersAndCachesAsTheBatchDoes) {
  auto g = diamond(5, 3);
  auto other = diamond(5, 3);
  const auto foreign = other->eligibility({}, std::vector<EdgeId> {0});
  DagCache foreign_cache(*other);
  DagCache building(*g);
  const auto first_error = [&](const Eligibility* elig, DagCache* cache, std::int32_t workers) {
    return error_of([&] { (void)sensitivity(*g, v({0}), v({3}), SensitivityOptions {global_config(), elig, cache}, workers); });
  };
  EXPECT_EQ(first_error(&foreign, &foreign_cache, 0), "workers must be >= 1 (got 0)");
  EXPECT_EQ(first_error(&foreign, &foreign_cache, 2), "foreign Eligibility: built for another graph");
  EXPECT_EQ(first_error(nullptr, &foreign_cache, 2), "foreign DagCache: bound to another graph");
  EXPECT_EQ(first_error(nullptr, &building, 2),
            "a supplied building DagCache has one writer: it is rejected unless workers == 1 (got workers = 2; a "
            "frozen cache may be shared by every worker)");
  EXPECT_EQ(stats_of(building.stats()), stats_of(DagCache(*g).stats()));  // no lookup ran
  const auto want = sensitivity(*g, v({0}), v({3}), SensitivityOptions {global_config()}, 1);
  building.freeze();
  for (const std::int32_t workers : {1, 2, 4}) {
    const auto r = sensitivity(*g, v({0}), v({3}), SensitivityOptions {global_config(), nullptr, &building}, workers);
    ASSERT_EQ(r.entries.size(), want.entries.size());
    for (std::size_t i = 0; i < r.entries.size(); ++i) {
      EXPECT_EQ(r.entries[i].edge, want.entries[i].edge);
      EXPECT_EQ(r.entries[i].loss, want.entries[i].loss);
    }
  }
}

// Only CostOnly reads the cache: through an empty frozen cache Global, CapacityAware and Pinned
// sensitivity give the entries of a call without a cache and count nothing; CostOnly's baseline
// raises the miss.
TEST(Evaluators, SensitivityReadsTheCacheOnlyUnderCostOnly) {
  auto g = diamond(5, 3, 2);
  Pinned pin;
  pin.bundles = {dag_of(*g, {0, 2}), dag_of(*g, {1, 3})};
  DagCache empty(*g);
  empty.freeze();
  for (const auto& c : {global_config(), cfg(CapacityAware {.flow_count = 2, .max_flow_count = 4}, FlowPlacement::MaxFlow),
                        cfg(pin, FlowPlacement::EqualBalanced)}) {
    const auto want = sensitivity(*g, v({0}), v({3}), SensitivityOptions {c}, 1);
    ASSERT_FALSE(want.entries.empty());
    for (const std::int32_t workers : {1, 3}) {
      const auto r = sensitivity(*g, v({0}), v({3}), SensitivityOptions {c, nullptr, &empty}, workers);
      expect_same(r.baseline, want.baseline, "through an empty frozen cache");
      ASSERT_EQ(r.entries.size(), want.entries.size());
      for (std::size_t i = 0; i < r.entries.size(); ++i) {
        EXPECT_EQ(r.entries[i].edge, want.entries[i].edge);
        EXPECT_EQ(r.entries[i].loss, want.entries[i].loss);
      }
    }
  }
  EXPECT_EQ(stats_of(empty.stats()), stats_of(DagCache(*g).stats()));
  EXPECT_THROW((void)sensitivity(*g, v({0}), v({3}), SensitivityOptions {cost_only_config(), nullptr, &empty}, 1),
               DagCacheMiss);
  EXPECT_EQ(empty.stats().misses, 1);
}

namespace {

// Sets (value) or removes (nullptr) an environment variable of this process.
void set_env(const char* name, const char* value) {
#ifdef _WIN32
  (void)_putenv_s(name, value == nullptr ? "" : value);  // an empty value removes it on Windows
#else
  if (value == nullptr) {
    (void)unsetenv(name);
  } else {
    (void)setenv(name, value, 1);
  }
#endif
}

} // namespace

// The evaluators' default worker counts (README "Environment variables", evaluator_checks.hpp):
// NGRAPH_CORE_BATCH_THREADS for batch_capacity and NGRAPH_CORE_SENSITIVITY_THREADS for
// sensitivity, each an integer in [1, 4096] read independently of the other; 1 when unset or
// invalid. " 4" reads as 4 (leading blanks are part of the integer syntax, as for strtol and
// Python's int()); trailing characters, a sign alone and values outside the range
// give 1.
TEST(Evaluators, DefaultWorkersParseTheEnvironment) {
  const char* names[] = {"NGRAPH_CORE_BATCH_THREADS", "NGRAPH_CORE_SENSITIVITY_THREADS"};
  std::optional<std::string> saved[2];
  for (int i = 0; i < 2; ++i) {
    if (const char* x = std::getenv(names[i])) saved[i] = x;
  }
  struct Restore {
    const char* const* names;
    std::optional<std::string>* saved;
    ~Restore() {
      for (int i = 0; i < 2; ++i) set_env(names[i], saved[i] ? saved[i]->c_str() : nullptr);
    }
  } restore {names, saved};
  const std::vector<std::pair<const char*, std::int32_t>> cases = {
      {nullptr, 1}, {"1", 1},   {"4", 4},     {"4096", 4096}, {"4097", 1}, {"0", 1},  {"-4", 1},
      {"abc", 1},   {"4x", 1},  {"4 ", 1},    {" 4", 4},      {"+", 1},    {"99999999999999999999", 1},
  };
  for (const bool sensitivity : {false, true}) {
    const int self = sensitivity ? 1 : 0;
    for (const auto& [value, want] : cases) {
      const std::string where =
          std::string(names[self]) + "=" + (value == nullptr ? "(unset)" : "'" + std::string(value) + "'");
      set_env(names[self], value);
      set_env(names[1 - self], "7");  // the other variable never applies
      EXPECT_EQ(detail::default_workers(sensitivity), want) << where;
      EXPECT_EQ(detail::default_workers(!sensitivity), 7) << where;
    }
  }
}
