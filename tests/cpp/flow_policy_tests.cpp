// The policy engine beyond the reference's fixtures (those replay in flow_policy_fixture_tests.cpp):
// configuration validation, the configuration a policy keeps and the route choice ShortestPaths and
// T1 read, request validation (a rejected request changes nothing) and the per-source round
// robin of the cumulative request, the outputs (the bottleneck under the policy's Eligibility),
// CentralizedTE's empty roster, its roster without flow and its zero call's cumulative curve, record
// routes (sparse storage, sharing, pruning), records in FlowId order, one flow to several nearest
// destinations, and the evaluators (evaluate_demand, evaluate_pairs, evaluate_edge_failures) with workers and their
// default from the environment. T1 (DistributedTE) has its own file, distributed_te_tests.cpp; the
// ShortestPathCache's own modes and key are tested in shortest_path_cache_tests.cpp.
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

#include "netgraph/core/shortest_path_cache.hpp"
#include "netgraph/core/detail/evaluator_checks.hpp"
#include "netgraph/core/detail/route_kernels.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/route_builders.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/detail/trace.hpp"
#include "netgraph/core/detail/search_detail.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/search.hpp"
#include "state_comparators.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

detail::RouteChoice sel(RouteTieBreak t = RouteTieBreak::FewerHops, bool multipath = true, bool parallel_edges = true) {
  return detail::RouteChoice {t, multipath, parallel_edges};
}

// A configuration: ShortestPaths takes the allocation `k` and the route choice `s`, T1 the tie_break
// and parallel_edges of `s`; the other models take neither.
FlowPolicyConfig cfg(Routing r, RouteAllocation k, Admission a = Admission::Lossless, detail::RouteChoice s = sel()) {
  if (auto* sp = std::get_if<ShortestPaths>(&r)) {
    sp->allocation = k;
    sp->tie_break = s.tie_break;
    sp->multipath = s.multipath;
    sp->parallel_edges = s.parallel_edges;
  }
  if (auto* ca = std::get_if<DistributedTE>(&r)) {
    ca->tie_break = s.tie_break;
    ca->parallel_edges = s.parallel_edges;
  }
  return FlowPolicyConfig {std::move(r), a};
}

// A T1 configuration: N tunnels per head, no placement.
FlowPolicyConfig t1(std::int64_t n, Admission a = Admission::Lossless) {
  return FlowPolicyConfig {.routing = DistributedTE {TunnelCount {n}}, .admission = a};
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
    EXPECT_TRUE(detail::LedgerAccess::bound_policy_ids(fg).empty());
  };
  rejects(cfg(ShortestPaths {}, RouteAllocation::Optimized, Admission::Lossy), "ShortestPaths OPTIMIZED requires admission=LOSSLESS");
  // ShortestPaths's weights come from its weight basis: HIGHER_HEADROOM is a tie rule like the others.
  EXPECT_EQ(error_of([&] {
              FlowPolicy p(fg, v({0}), v({3}),
                           cfg(ShortestPaths {}, RouteAllocation::Weighted, Admission::Lossless,
                               sel(RouteTieBreak::HigherHeadroom, true, false)));
            }),
            "");
  fg.reset();
  rejects({.routing = centralized(), .admission = Admission::Lossy}, "CentralizedTE without tunnel_count requires admission=LOSSLESS");
  // T1: one exact tunnel count in [1, 10000], the ceiling's and the search minimum's own texts (the
  // T1 rows of the matrix are distributed_te_tests.cpp's).
  rejects(t1(0), "tunnel count must be in [1, 10000]");
  rejects({.routing = DistributedTE {TunnelCount {2, true}}}, "DistributedTE requires an exact tunnel count");
  DistributedTE ca {TunnelCount {1}};
  ca.min_unreserved_bandwidth = -1;
  rejects({.routing = ca}, "min_unreserved_bandwidth must be an integer in [0, 2^62)");
  // The upper bound: 2^62 - 1 units is a valid minimum, 2^62 is not.
  ca.min_unreserved_bandwidth = kDomainLimit;
  rejects({.routing = ca}, "min_unreserved_bandwidth must be an integer in [0, 2^62)");
  ca.min_unreserved_bandwidth = kDomainLimit - 1;
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({3}), {.routing = ca}); }), "");
  fg.reset();
  ca = DistributedTE {TunnelCount {1}};
  ca.cost_ceiling = AbsoluteCost {-1};
  rejects({.routing = ca}, "cost_ceiling AbsoluteCost cost must be >= 0 (got -1)");
  // A factor below 1 is a valid value, refused where the ceiling is used, with its reduced
  // components (the numerator at its lowest included).
  for (const auto& [f, shown] : {std::pair {CostStretch(1, 2), "1/2"}, std::pair {CostStretch(4, 6), "2/3"},
                                 std::pair {CostStretch(std::numeric_limits<std::int64_t>::min()), "-9223372036854775808/1"}}) {
    ca = DistributedTE {TunnelCount {1}};
    ca.cost_ceiling = f;
    rejects({.routing = ca}, std::string("cost_ceiling CostStretch factor must be >= 1 (got ") + shown + ")");
  }
  // Pinned: a program of this graph whose endpoints are the policy's.
  rejects({.routing = Pinned {}}, "Pinned requires a program");
  const Pinned ok = pinned_dags(*g, 0, 3, {dag_of(*g, {0, 2})});
  rejects({.routing = ok}, "program endpoints do not match policy endpoints", {0, 1}, {3});
  auto other = diamond(5, 5);
  rejects({.routing = pinned_dags(*other, 0, 3, {dag_of(*other, {0, 2})})}, "program belongs to another graph");
  rejects({.routing = Pinned {ok.program, std::make_shared<const Eligibility>(other->eligibility({}, std::vector<EdgeId> {0}))}},
          "black_holes belongs to another graph");
  rejects({.routing = centralized(SourceAllocation::Optimized, dag_of(*other, {0, 2}))}, "support belongs to another graph");
  rejects(cfg(ShortestPaths {}, RouteAllocation::Optimized),
          "node 3 is in both endpoint sets: sources and destinations must not overlap", {0, 3}, {3});
  rejects(cfg(ShortestPaths {}, RouteAllocation::Optimized), "sources: the node set is empty", {}, {3});
  rejects(cfg(ShortestPaths {}, RouteAllocation::Optimized), "sources: node 9 out of node range [0, 4)", {9}, {3});
  // A routing variant left without a mode (valueless) by an exception while a mode was emplaced.
  struct ThrowingConversion {
    operator Pinned() const { throw std::runtime_error("conversion failed"); }
  };
  auto valueless = cfg(ShortestPaths {}, RouteAllocation::Optimized);
  EXPECT_THROW(valueless.routing.emplace<Pinned>(ThrowingConversion {}), std::runtime_error);
  ASSERT_TRUE(valueless.routing.valueless_by_exception());
  rejects(valueless, "routing: the variant holds no routing mode");
  const auto el = other->eligibility({}, std::vector<EdgeId> {0});
  EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({3}), cfg(ShortestPaths {}, RouteAllocation::Optimized), &el); }),
            "foreign Eligibility: built for another graph");
  // The endpoints come before the ownership of the mask and the cache.
  ShortestPathCache foreign_first(*other);
  EXPECT_EQ(error_of([&] {
              FlowPolicy p(fg, std::vector<NodeId> {}, v({3}), cfg(ShortestPaths {}, RouteAllocation::Optimized), &el,
                           &foreign_first);
            }),
            "sources: the node set is empty");
  ShortestPathCache foreign_cache(*other);
  EXPECT_EQ(error_of([&] {
              FlowPolicy p(fg, v({0}), v({3}), cfg(ShortestPaths {}, RouteAllocation::Optimized), nullptr, &foreign_cache);
            }),
            "foreign ShortestPathCache: bound to another graph");
  EXPECT_TRUE(detail::LedgerAccess::bound_policy_ids(fg).empty());
}

// FlowPolicyConfig and the routing models are plain values: they hold any field, a value outside
// its enumeration included (here 96 to 99, reachable by a cast), and their consumers check it
// (EnumValuesOutsideTheirEnumerationAreRejected and SeveralDefectsReportTheFirstInTheContractOrder
// below). CentralizedTE requires its source allocation: it has no default constructor.
TEST(FlowPolicyConfig, PlainValuesHoldEnumValuesOutsideTheirEnumeration) {
  static_assert(!std::is_default_constructible_v<CentralizedTE>);
  EXPECT_NO_THROW((void)(FlowPolicyConfig {centralized(static_cast<SourceAllocation>(96)), static_cast<Admission>(98)}));
  ShortestPaths sp;
  sp.allocation = static_cast<RouteAllocation>(97);
  sp.source_allocation = static_cast<SourceAllocation>(96);
  sp.weight_basis = static_cast<WeightBasis>(95);
  sp.tie_break = static_cast<RouteTieBreak>(99);
  EXPECT_NO_THROW((void)(FlowPolicyConfig {.routing = sp}));
  DistributedTE ca {TunnelCount {1}};
  ca.tie_break = static_cast<RouteTieBreak>(99);
  EXPECT_NO_THROW((void)(FlowPolicyConfig {.routing = ca}));
}

// Enumeration checks at every consuming entry: the fields stay assignable after construction, so
// a value outside its enumeration is a domain error named by its field (allocation,
// source_allocation, weight_basis, tie_break, admission), before any matrix check, at every entry
// that takes the configuration; it never reaches a record.
TEST(FlowPolicyConfig, EnumValuesOutsideTheirEnumerationAreRejected) {
  auto g = diamond(5, 5);
  FlowLedger fg(g);
  const Pinned pinned = pinned_dags(*g, 0, 3, {dag_of(*g, {0, 2})});
  std::vector<std::pair<FlowPolicyConfig, std::string>> cases;
  const auto sp_with = [](const std::function<void(ShortestPaths&)>& set) {
    ShortestPaths sp;
    sp.allocation = RouteAllocation::Optimized;
    sp.source_allocation = SourceAllocation::Optimized;  // valid, then a later defect
    set(sp);
    return FlowPolicyConfig {.routing = sp, .admission = Admission::Lossy};  // Lossy: a matrix defect after
  };
  cases.push_back({sp_with([](ShortestPaths& sp) { sp.allocation = static_cast<RouteAllocation>(99); }),
                   "allocation: 99 is not a value of RouteAllocation"});
  cases.push_back({sp_with([](ShortestPaths& sp) { sp.source_allocation = static_cast<SourceAllocation>(99); }),
                   "source_allocation: 99 is not a value of SourceAllocation"});
  cases.push_back({sp_with([](ShortestPaths& sp) { sp.weight_basis = static_cast<WeightBasis>(99); }),
                   "weight_basis: 99 is not a value of WeightBasis"});
  cases.push_back({sp_with([](ShortestPaths& sp) { sp.tie_break = static_cast<RouteTieBreak>(99); }),
                   "tie_break: 99 is not a value of RouteTieBreak"});
  cases.push_back({{.routing = centralized(static_cast<SourceAllocation>(99)), .admission = Admission::Lossy},
                   "source_allocation: 99 is not a value of SourceAllocation"});
  DistributedTE ca {TunnelCount {0}};
  ca.tie_break = static_cast<RouteTieBreak>(99);
  cases.push_back({{.routing = ca}, "tie_break: 99 is not a value of RouteTieBreak"});
  for (const Routing& r : {Routing {ShortestPaths {}}, Routing {DistributedTE {TunnelCount {0}}}, Routing {centralized()},
                           Routing {pinned}}) {
    cases.push_back({{.routing = r, .admission = static_cast<Admission>(99)}, "admission: 99 is not a value of Admission"});
  }
  for (const auto& [c, text] : cases) {
    EXPECT_EQ(error_of([&] { FlowPolicy p(fg, v({0}), v({3}), c); }), text);
    EXPECT_TRUE(detail::LedgerAccess::bound_policy_ids(fg).empty());
    EXPECT_TRUE(fg.records().empty());
    EvaluationOptions co {c};
    EXPECT_EQ(error_of([&] { (void)evaluate_demand(*g, v({0}), v({3}), co); }), text);
    EXPECT_EQ(error_of([&] { (void)evaluate_edge_failures(*g, v({0}), v({3}), EvaluationOptions {c}, 1); }), text);
    if (!std::holds_alternative<Pinned>(c.routing)) {
      const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{v({0}), v({3})}};
      EXPECT_EQ(error_of([&] { (void)evaluate_pairs(*g, pairs, co, {}, 1); }), text);
    }
  }
  // The kernel check, the search and the cache's key.
  const auto bad_allocation = static_cast<RouteAllocation>(99);
  const auto bad_direction = static_cast<Direction>(99);
  EXPECT_TRUE(contains(error_of([&] { detail::check_kernel(bad_allocation, Admission::Lossless); }), "allocation: 99"));
  EXPECT_TRUE(contains(error_of([&] { detail::check_kernel(RouteAllocation::Optimized, static_cast<Admission>(99)); }),
                       "admission: 99"));
  {
    SpfOptions o(Direction::Forward);
    o.direction = bad_direction;
    EXPECT_TRUE(contains(error_of([&] { (void)spf(*g, v({0}), o); }), "direction: 99 is not a value of Direction"));
    EXPECT_TRUE(contains(error_of([&] { SpfOptions bad(bad_direction); }), "direction: 99 is not a value of Direction"));
    SpfOptions t(Direction::Forward);
    t.tie_break = static_cast<RouteTieBreak>(99);
    EXPECT_TRUE(contains(error_of([&] { (void)spf(*g, v({0}), t); }), "tie_break: 99 is not a value of RouteTieBreak"));
  }
  ShortestPathCache cache(*g);
  EXPECT_TRUE(contains(error_of([&] { (void)cache.lookup(*g, v({3}), bad_direction, nullptr); }), "direction: 99"));
  const auto st = cache.stats();
  EXPECT_EQ(st.hits + st.misses + st.entries, 0);
}

// A configuration with two defects at once: one error, the first in the documented order:
// the enumeration values; the matrix (route and source allocation, admission, the ineffective tie
// rule, the pooled single-path rule, prefer_headroom, the weight combinations, then a cache on a
// model that does not read it); the tunnel count, an "up to" count for T1, the ceilings, the search
// minimum and edge_weights' length and entries; the tunnel limit; a Pinned program's endpoints.
// Exact texts; the same first error at construction and at the evaluators; nothing is bound.
TEST(FlowPolicyConfig, SeveralDefectsReportTheFirstInTheContractOrder) {
  auto g = diamond(5, 5);
  FlowLedger fg(g);
  ShortestPathCache cache(*g);
  const auto weights = [](std::vector<Amount> w) { return std::make_shared<const std::vector<Amount>>(std::move(w)); };
  const auto sp = [](RouteAllocation k, std::optional<SourceAllocation> source = std::nullopt) {
    ShortestPaths out;
    out.allocation = k;
    out.source_allocation = source;
    return out;
  };
  struct Case {
    const char* why;
    FlowPolicyConfig c;
    std::vector<NodeId> sources;
    bool with_cache;
    std::string error;
  };
  std::vector<Case> cases;
  {
    ShortestPaths x = sp(RouteAllocation::Equal, SourceAllocation::Optimized);
    x.tie_break = RouteTieBreak::HigherCapacity;  // also an ineffective tie
    cases.push_back({"source allocation before the tie rule", {x, Admission::Lossless}, {0}, false,
                     "ShortestPaths with EQUAL or WEIGHTED requires source_allocation=EQUAL"});
  }
  {
    ShortestPaths x = sp(RouteAllocation::Optimized);
    x.tie_break = RouteTieBreak::HigherCapacity;
    cases.push_back({"admission before the tie rule", {x, Admission::Lossy}, {0}, false,
                     "ShortestPaths OPTIMIZED requires admission=LOSSLESS"});
    x.multipath = false;
    cases.push_back({"the pooled rule before the weights", {x, Admission::Lossless}, {0, 1}, false,
                     "ShortestPaths with multiple input sources requires tie_break=FEWER_HOPS when multipath=False"});
  }
  {
    ShortestPaths x = sp(RouteAllocation::Equal);
    x.tie_break = RouteTieBreak::HigherCapacity;
    x.weight_basis = WeightBasis::Capacity;
    cases.push_back({"the tie rule before the weights", {x, Admission::Lossless}, {0}, false,
                     "ShortestPaths tie_break requires multipath=False or parallel_edges=False"});
  }
  {
    // V3: two sources, a single path under HIGHER_CAPACITY and explicit weights of the wrong length.
    ShortestPaths x = sp(RouteAllocation::Weighted);
    x.weight_basis = WeightBasis::Explicit;
    x.edge_weights = weights({1});
    x.multipath = false;
    x.tie_break = RouteTieBreak::HigherCapacity;
    cases.push_back({"the pooled rule before the vectors", {x, Admission::Lossless}, {0, 1}, false,
                     "ShortestPaths with multiple input sources requires tie_break=FEWER_HOPS when multipath=False"});
    x.tie_break = RouteTieBreak::FewerHops;
    cases.push_back({"the vectors last", {x, Admission::Lossless}, {0, 1}, false, "edge_weights: expected 4 entries"});
  }
  {
    ShortestPaths x = sp(RouteAllocation::Equal);
    x.weight_basis = WeightBasis::Explicit;
    cases.push_back({"weights on another allocation", {x, Admission::Lossless}, {0}, true,
                     "weight_basis requires allocation=WEIGHTED"});
    x.allocation = RouteAllocation::Weighted;
    cases.push_back({"explicit weights without the array", {x, Admission::Lossless}, {0}, false,
                     "weight_basis=EXPLICIT requires edge_weights"});
    x.weight_basis = WeightBasis::Capacity;
    x.edge_weights = weights({1, 1, 1, 1});
    cases.push_back({"an array without explicit weights", {x, Admission::Lossless}, {0}, false,
                     "edge_weights require weight_basis=EXPLICIT"});
    x.weight_basis = WeightBasis::Explicit;
    x.edge_weights = weights({1, 1, kDomainLimit, -1});
    cases.push_back({"the first entry out of range", {x, Admission::Lossless}, {0}, false,
                     "edge_weights[2] must be an integer in [0, 2^62)"});
  }
  // CentralizedTE: T2's several-sources rule before its tunnel count (today it came after); a
  // live controller's admission before prefer_headroom; prefer_headroom before the cache.
  cases.push_back({"T2's sources before the count", {centralized(SourceAllocation::Optimized, {}, TunnelCount {0}), Admission::Lossy},
                   {0, 1}, false, "CentralizedTE with tunnel_count and several sources requires source_allocation=EQUAL"});
  cases.push_back({"admission before prefer_headroom", {centralized(SourceAllocation::Equal, {}, std::nullopt, true), Admission::Lossy},
                   {0}, false, "CentralizedTE without tunnel_count requires admission=LOSSLESS"});
  cases.push_back({"prefer_headroom before the cache", {centralized(SourceAllocation::Equal, {}, std::nullopt, true)}, {0}, true,
                   "prefer_headroom requires tunnel_count"});
  cases.push_back({"the cache before T2's count", {centralized(SourceAllocation::Equal, {}, TunnelCount {0})}, {0}, true,
                   "CentralizedTE does not use ShortestPathCache"});
  cases.push_back({"T2's count", {centralized(SourceAllocation::Equal, {}, TunnelCount {0})}, {0}, false,
                   "tunnel count must be in [1, 10000]"});
  // T1: the cache, then the count before up_to, up_to before the ceiling, the ceiling before the
  // search minimum.
  cases.push_back({"the cache before T1's count", {DistributedTE {TunnelCount {0}}}, {0}, true,
                   "DistributedTE does not use ShortestPathCache"});
  cases.push_back({"the count before up_to", {DistributedTE {TunnelCount {0, true}}}, {0}, false,
                   "tunnel count must be in [1, 10000]"});
  {
    DistributedTE x {TunnelCount {2, true}, AbsoluteCost {-1}};
    cases.push_back({"up_to before the ceiling", {x}, {0}, false, "DistributedTE requires an exact tunnel count"});
    x.tunnel_count.up_to = false;
    x.min_unreserved_bandwidth = kDomainLimit;
    cases.push_back({"the ceiling before the minimum", {x}, {0}, false, "cost_ceiling AbsoluteCost cost must be >= 0 (got -1)"});
    x.cost_ceiling.reset();
    cases.push_back({"the minimum", {x}, {0}, false, "min_unreserved_bandwidth must be an integer in [0, 2^62)"});
  }
  // Pinned: the cache before the program's endpoints.
  const Pinned other_endpoints = pinned_dags(*g, 0, 3, {dag_of(*g, {0, 2})});
  cases.push_back({"the cache before the endpoints", {other_endpoints}, {0, 1}, true, "Pinned does not use ShortestPathCache"});
  cases.push_back({"the endpoints", {other_endpoints}, {0, 1}, false, "program endpoints do not match policy endpoints"});
  for (const auto& k : cases) {
    ShortestPathCache* c = k.with_cache ? &cache : nullptr;
    EXPECT_EQ(error_of([&] { FlowPolicy p(fg, k.sources, v({3}), k.c, nullptr, c); }), k.error) << k.why;
    EXPECT_TRUE(detail::LedgerAccess::bound_policy_ids(fg).empty()) << k.why;
    EXPECT_TRUE(fg.records().empty()) << k.why;
    const EvaluationOptions co {k.c, nullptr, c};
    EXPECT_EQ(error_of([&] { (void)evaluate_demand(*g, k.sources, v({3}), co); }), k.error) << k.why;
    EXPECT_EQ(error_of([&] { (void)evaluate_edge_failures(*g, k.sources, v({3}), co, 1); }), k.error) << k.why;
  }
  EXPECT_EQ(cache.stats().misses + cache.stats().hits, 0);
}

// Every field has a default but CentralizedTE's source allocation: ShortestPaths() is ECMP (EQUAL),
// and the policy's snapshot resolves the absent source allocation (the allocation's default) and
// the absent weight basis of WEIGHTED (CAPACITY); SP() and SP(allocation=EQUAL) give one snapshot.
TEST(FlowPolicyConfig, ShortestPathsDefaults) {
  auto g = diamond(5, 5);
  const FlowPolicyConfig empty;
  ASSERT_TRUE(std::holds_alternative<ShortestPaths>(empty.routing));
  const auto& d = std::get<ShortestPaths>(empty.routing);
  EXPECT_EQ(d.allocation, RouteAllocation::Equal);
  EXPECT_FALSE(d.source_allocation.has_value() || d.weight_basis.has_value() || d.edge_weights);
  EXPECT_EQ(d.tie_break, RouteTieBreak::FewerHops);
  EXPECT_TRUE(d.multipath && d.parallel_edges);
  EXPECT_EQ(empty.admission, Admission::Lossless);
  const auto snapshot = [&](FlowPolicyConfig c) {
    FlowLedger own(g);
    FlowPolicy p(own, v({0}), v({3}), std::move(c));
    return p.config();
  };
  const auto resolved = [&](ShortestPaths x) { return std::get<ShortestPaths>(snapshot({.routing = x}).routing); };
  const auto ecmp = resolved(ShortestPaths {});
  EXPECT_EQ(ecmp.source_allocation, SourceAllocation::Equal);
  EXPECT_FALSE(ecmp.weight_basis.has_value());
  EXPECT_EQ(ecmp.allocation, resolved(shortest_paths(RouteAllocation::Equal)).allocation);
  EXPECT_EQ(ecmp.source_allocation, resolved(shortest_paths(RouteAllocation::Equal)).source_allocation);
  const auto weighted = resolved(shortest_paths(RouteAllocation::Weighted));
  EXPECT_EQ(weighted.source_allocation, SourceAllocation::Equal);
  EXPECT_EQ(weighted.weight_basis, WeightBasis::Capacity);
  const auto optimized = resolved(shortest_paths(RouteAllocation::Optimized));
  EXPECT_EQ(optimized.source_allocation, SourceAllocation::Optimized);
  EXPECT_FALSE(optimized.weight_basis.has_value());
  ShortestPaths equal_controller = shortest_paths(RouteAllocation::Optimized);
  equal_controller.source_allocation = SourceAllocation::Equal;
  EXPECT_EQ(resolved(equal_controller).source_allocation, SourceAllocation::Equal);
  // The other models keep what they were given.
  EXPECT_EQ(std::get<CentralizedTE>(snapshot({.routing = centralized(SourceAllocation::Equal)}).routing).source_allocation,
            SourceAllocation::Equal);
  const auto tunnels = std::get<DistributedTE>(snapshot({.routing = DistributedTE {TunnelCount {2}}}).routing);
  EXPECT_EQ(tunnels.tie_break, RouteTieBreak::FewerHops);
  EXPECT_TRUE(tunnels.parallel_edges);
}

// config() is the copy taken at construction: every field as given, and a later change of the
// caller's object changes nothing.
TEST(FlowPolicyConfig, ConfigIsTheCopyTakenAtConstruction) {
  auto g = diamond(5, 5);
  FlowLedger fg(g);
  DistributedTE ca {TunnelCount {2}};
  ca.cost_ceiling = AbsoluteCost {9};
  ca.min_unreserved_bandwidth = 1;
  ca.tie_break = RouteTieBreak::HigherAvailableCapacity;
  ca.parallel_edges = false;
  FlowPolicyConfig c {.routing = ca, .admission = Admission::Lossy};
  FlowPolicy p(fg, v({0}), v({3}), c);
  c.admission = Admission::Lossless;
  std::get<DistributedTE>(c.routing).tunnel_count.per_head = 7;
  std::get<DistributedTE>(c.routing).tie_break = RouteTieBreak::FewerHops;
  const auto& got = p.config();
  EXPECT_NE(&got, &c);
  EXPECT_EQ(got.admission, Admission::Lossy);
  ASSERT_TRUE(std::holds_alternative<DistributedTE>(got.routing));
  const auto& r = std::get<DistributedTE>(got.routing);
  EXPECT_EQ(r.tunnel_count, (TunnelCount {2, false}));
  EXPECT_TRUE(r.cost_ceiling == std::optional<CostCeiling> {AbsoluteCost {9}});
  EXPECT_EQ(r.min_unreserved_bandwidth, 1);
  EXPECT_EQ(r.tie_break, RouteTieBreak::HigherAvailableCapacity);
  EXPECT_FALSE(r.parallel_edges);
  // Pinned: the supplied program and black holes themselves (shared, never copied); explicit
  // weights: the configuration's one array, shared by every policy built from it.
  const auto dark = std::make_shared<const Eligibility>(g->eligibility({}, std::vector<EdgeId> {1}));
  const Pinned pin {pinned_dags(*g, 0, 3, {dag_of(*g, {1, 3}), dag_of(*g, {0, 2})}).program, dark};
  FlowLedger fg2(g);
  FlowPolicy q(fg2, v({0}), v({3}), {.routing = pin});
  ASSERT_TRUE(std::holds_alternative<Pinned>(q.config().routing));
  EXPECT_EQ(std::get<Pinned>(q.config().routing).program.get(), pin.program.get());
  EXPECT_EQ(std::get<Pinned>(q.config().routing).black_holes.get(), dark.get());
  ShortestPaths wcmp = shortest_paths(RouteAllocation::Weighted);
  wcmp.weight_basis = WeightBasis::Explicit;
  wcmp.edge_weights = std::make_shared<const std::vector<Amount>>(std::vector<Amount> {1, 2, 3, 4});
  const FlowPolicyConfig explicit_weights {.routing = wcmp};
  FlowPolicy a(fg2, v({0}), v({3}), explicit_weights);
  FlowPolicy b(fg2, v({0}), v({3}), explicit_weights);
  EXPECT_EQ(std::get<ShortestPaths>(a.config().routing).edge_weights.get(), wcmp.edge_weights.get());
  EXPECT_EQ(std::get<ShortestPaths>(b.config().routing).edge_weights.get(), wcmp.edge_weights.get());
}

// The tie rule is no weight. Two parallel edges 0 -> 1 of 4 units (EdgeIds 0 and 1), a third party
// holding 3 units of edge 0, request 5. ShortestPaths Weighted over both edges: its weight basis
// decides (Capacity 4 : 4 carries at most 1 + 1 = 2, edge 0 having 1 unit left; a residual
// snapshot 1 : 4 at preparation carries 1 + 4 = 5) under the one tie rule that keeps both edges,
// FEWER_HOPS (a HIGHER_* tie_break needs multipath=False or parallel_edges=False, either of which
// keeps one edge). CentralizedTE, CentralizedTE on the two edges as its support, and Pinned over one
// split route of both edges (Equal 2; Weighted with capacity weights 2, with the snapshot 5) take
// no route choice. Control: ShortestPaths reads the route choice (a single path, multipath = false,
// carries less than the two-edge DAG).
TEST(FlowPolicyConfig, TheTieRuleIsNoWeight) {
  auto g = shared(make_graph(2, {{0, 1, 4, 1}, {0, 1, 4, 1}}));
  using Edges = std::vector<std::pair<EdgeId, Amount>>;
  // The placement of request 5 after the third party: (call_delivered, the policy's records' edge_flows).
  auto run = [&](const FlowPolicyConfig& c) {
    FlowLedger fg(g);
    FlowPolicy other(fg, v({0}), v({1}), FlowPolicyConfig {.routing = centralized(SourceAllocation::Optimized, dag_of(*g, {0}))});
    EXPECT_EQ(other.place(3).call_delivered, 3);
    FlowPolicy p(fg, v({0}), v({1}), c);
    const auto r = p.place(5);
    std::vector<Edges> edges;
    for (const auto& rec : p.records()) edges.push_back(rec.edge_flows);
    return std::pair {r.call_delivered, edges};
  };
  for (auto basis : {WeightBasis::Capacity, WeightBasis::ResidualSnapshot}) {
    const bool snapshot = basis == WeightBasis::ResidualSnapshot;
    ShortestPaths wcmp = shortest_paths(RouteAllocation::Weighted);
    wcmp.weight_basis = basis;
    const auto got = run({.routing = wcmp});
    EXPECT_EQ(got.first, snapshot ? 5 : 2);
    EXPECT_EQ(got.second, (std::vector<Edges> {snapshot ? Edges {{0, 1}, {1, 4}} : Edges {{0, 1}, {1, 1}}}));
  }
  const auto both = dag_of(*g, {0, 1});
  EXPECT_EQ(run({.routing = centralized()}).first, 5);
  EXPECT_EQ(run({.routing = centralized(SourceAllocation::Optimized, both)}).first, 5);
  EXPECT_EQ(run({.routing = pinned_dags(*g, 0, 1, {both})}).first, 2);
  EXPECT_EQ(run({.routing = pinned_dags(*g, 0, 1, {both}, RouteAllocation::Weighted)}).first, 2);
  RouteOptions snap;
  snap.allocation = RouteAllocation::Weighted;
  snap.weight_basis = WeightBasis::ResidualSnapshot;
  snap.residual = std::vector<Amount> {1, 4};
  const NodeId s0[1] = {0};
  const NodeId d1[1] = {1};
  const RouteDefinition route[1] = {RouteDefinition {0, both, 1}};
  EXPECT_EQ(run({.routing = Pinned {ForwardingProgram::from_routes(*g, s0, d1, route, snap), nullptr}}).first, 5);
  EXPECT_EQ(run(cfg(ShortestPaths {}, RouteAllocation::Optimized, Admission::Lossless, sel())).first, 5);
  EXPECT_EQ(run(cfg(ShortestPaths {}, RouteAllocation::Optimized, Admission::Lossless,
                    sel(RouteTieBreak::FewerHops, false)))
                .first,
            1);
}

namespace {

// What a rejected call must leave as it was: the policy's epoch state (the requested total, the
// +inf flag, delivered(), the CentralizedTE curve and each roster source's flows with their routes)
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
// routing model with incremental placement (T1's calls: distributed_te_tests.cpp); validation
// precedes the shortcut of a zero request (place(0) after +inf is rejected too), and remove() ends
// the epoch, after which a finite request places again. The 2^127 case starts from a requested
// total set just below the limit.
TEST(FlowPolicyRequests, RejectedRequestsLeaveThePolicyAndTheLedger) {
  auto g = diamond(5, 3, 2);
  const Pinned pinned = pinned_dags(*g, 0, 3, {dag_of(*g, {0, 2}), dag_of(*g, {1, 3})});
  const std::vector<std::pair<const char*, FlowPolicyConfig>> configs = {
      {"ShortestPaths", shortest_paths_config()},
      {"CentralizedTE", global_config()},
      {"Pinned", {.routing = pinned}},
      {"CentralizedTE support", {.routing = centralized(SourceAllocation::Optimized, dag_of(*g, {0, 2}))}},
      {"ShortestPaths Lossy", cfg(ShortestPaths {}, RouteAllocation::Equal, Admission::Lossy)}};
  const std::string range_error = "offer must be an integer in [0, 2^62) or permitted +inf";
  const std::string after_inf = "finite placement after infinite placement requires a new epoch";
  for (const auto& [name, config] : configs) {
    FlowLedger fg(g);
    FlowPolicy p(fg, v({0}), v({3}), config);
    ASSERT_GT(p.place(4).call_delivered, 0) << name;
    const auto finite = state_of(p, fg);
    EXPECT_EQ(error_of([&] { (void)p.place(-1); }), range_error) << name;
    EXPECT_EQ(error_of([&] { (void)p.place(kDomainLimit); }), range_error) << name;
    expect_same(state_of(p, fg), finite, name);
    // Just below 2^127: 5 more would reach it, 4 would not.
    auto& requested = detail::PolicyAccess::engine(p).epoch.requested;
    const Counter128 saved = requested;
    requested = Counter128 {(std::uint64_t {1} << 63) - 1, std::numeric_limits<std::uint64_t>::max() - 4};
    const auto near_limit = state_of(p, fg);
    EXPECT_EQ(error_of([&] { (void)p.place(5); }), "requested limit exceeded") << name;
    expect_same(state_of(p, fg), near_limit, name);
    requested = saved;
    (void)p.place(std::nullopt);
    const auto infinite = state_of(p, fg);
    ASSERT_TRUE(infinite.epoch.infinite) << name;
    EXPECT_EQ(error_of([&] { (void)p.place(1); }), after_inf) << name;
    EXPECT_EQ(error_of([&] { (void)p.place(0); }), after_inf) << name;
    EXPECT_EQ(error_of([&] { (void)p.place(-1); }), range_error) << name;
    expect_same(state_of(p, fg), infinite, name);
    p.remove();
    EXPECT_EQ(p.place(1).call_delivered, 1) << name;
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
  const std::string offered_error = "offered limit exceeded";
  const std::string drop_error = "drops limit exceeded";
  struct Case {
    std::uint64_t offered_below;
    std::optional<std::uint64_t> drop_below;
    std::optional<std::string> error;
  };
  for (const auto& c : {Case {10, std::nullopt, std::nullopt}, Case {9, std::nullopt, offered_error},
                        Case {7, 7, offered_error}, Case {6, 6, drop_error}}) {
    const std::string name = "L - " + std::to_string(c.offered_below);
    FlowLedger fg(g);
    FlowPolicy p(fg, v({0}), v({2}), cfg(ShortestPaths {}, RouteAllocation::Equal, Admission::Lossy));
    ASSERT_EQ(p.place(1).call_delivered, 1) << name;
    std::vector<detail::RecordDraft> drafts(1);
    drafts[0].record = p.records().at(0);
    drafts[0].record_counter = drafts[0].record.flow_id.second;
    drafts[0].record.offered = below(c.offered_below);
    if (c.drop_below.has_value()) drafts[0].record.edge_drops = {{1, below(*c.drop_below)}};
    (void)detail::LedgerAccess::commit(fg, detail::PolicyAccess::engine(p).policy_id, drafts);
    const auto before = state_of(p, fg);
    if (c.error.has_value()) {
      EXPECT_EQ(error_of([&] { (void)p.place(std::nullopt); }), *c.error) << name;
      expect_same(state_of(p, fg), before, name);
      continue;
    }
    const auto r = p.place(std::nullopt);
    EXPECT_EQ(r.call_delivered, 3) << name;
    EXPECT_EQ(r.call_edge_drops, (std::vector<std::pair<EdgeId, Amount>> {{1, 6}})) << name;
    const auto records = p.records();
    ASSERT_EQ(records.size(), 1U) << name;
    EXPECT_EQ(records[0].offered, below(1)) << name;
    EXPECT_EQ(records[0].edge_drops, (std::vector<std::pair<EdgeId, Counter128>> {{1, Counter128 {0, 6}}})) << name;
    EXPECT_EQ(records[0].delivered, 4) << name;
  }
}

// The equal source allocation's shares (round robin from the lowest source): a total T gives source i of k
// floor(T / k) + [i < T mod k]. Sources 0 .. k - 1 with their own 100-unit edge to k. k = 3,
// Lossless: place(2) admits 2, shares [1, 1, 0]; the next place(2) admits 4, shares [2, 1, 1] (the
// admitted total split from zero), which a rejected call would move. Lossy injects the increments
// of the cumulative offer O; k = 7, where 2^64 = 2 mod 7: from O = 2^64 place(1) injects its unit
// at source 2, from O = 2^64 + 5 (= 0 mod 7) at source 0, so the high word of O counts with its
// exact weight.
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
    FlowPolicy p(fg, sources(3), v({3}), shortest_paths_config());
    EXPECT_EQ(p.place(2).call_delivered, 2);
    EXPECT_EQ(delivered_by_source(p), (Shares {{0, 1}, {1, 1}, {2, 0}}));
    EXPECT_EQ(error_of([&] { (void)p.place(-1); }), "offer must be an integer in [0, 2^62) or permitted +inf");
    EXPECT_EQ(error_of([&] { (void)p.place(kDomainLimit); }), "offer must be an integer in [0, 2^62) or permitted +inf");
    EXPECT_EQ(p.place(2).call_delivered, 2);
    EXPECT_EQ(delivered_by_source(p), (Shares {{0, 2}, {1, 1}, {2, 1}}));
  }
  // Lossy injects the round-robin increments of the cumulative offer O (Lossless re-splits its
  // admitted total from zero, above): the next unit goes to source O mod k.
  auto g = star(7);
  for (const auto& [offered, to] : {std::pair {Counter128 {1, 0}, NodeId {2}}, std::pair {Counter128 {1, 5}, NodeId {0}}}) {
    FlowLedger fg(g);
    FlowPolicy p(fg, sources(7), v({7}), cfg(ShortestPaths {}, RouteAllocation::Equal, Admission::Lossy));
    detail::PolicyAccess::engine(p).epoch.offered = offered;
    EXPECT_EQ(p.place(1).call_delivered, 1) << offered.lo;
    for (const auto& [s, x] : delivered_by_source(p)) EXPECT_EQ(x, s == to ? 1 : 0) << offered.lo << " source " << s;
    EXPECT_EQ(delivered_by_source(p).at(to), 1) << offered.lo;
  }
}

TEST(FlowPolicyOutputs, OwnAllocationLedgerResidualBottleneckAndGlobalCut) {
  auto g = diamond(5, 3);
  FlowLedger fg(g);
  // A third party holds one unit of the first route.
  const CentralizedTE third = centralized(SourceAllocation::Optimized, dag_of(*g, {0, 2}));
  FlowPolicy other(fg, v({0}), v({3}), cfg(third, RouteAllocation::Optimized));
  ASSERT_EQ(other.place(1).call_delivered, 1);
  FlowPolicy p(fg, v({0}), v({3}), shortest_paths_config());
  const OutputOptions all {true, true, true};
  auto r = p.place(4, all);
  // Both routes cost 2: Equal puts 2 units on each, so every edge carries 2 of this
  // policy's units and the residuals are 5 - 1 - 2 on route 0 -> 1 -> 3 and 3 - 2 on the other.
  // Every edge keeps a residual unit, so the sources reach everything: no saturated edge.
  EXPECT_EQ(r.call_delivered, 4);
  ASSERT_TRUE(r.edge_flows && r.residual && r.saturated_boundary_edges);
  EXPECT_FALSE(r.min_cut.has_value());  // not CentralizedTE
  EXPECT_FALSE(r.marginal_cost_curve.has_value());
  std::vector<Amount> own(4, 0);
  for (const auto& rec : p.records()) {
    for (const auto& [e, x] : rec.edge_flows) own[static_cast<std::size_t>(e)] += x;
  }
  EXPECT_EQ(*r.edge_flows, (std::vector<Amount> {2, 2, 2, 2}));
  EXPECT_EQ(*r.edge_flows, own);
  EXPECT_EQ(*r.residual, (std::vector<Amount> {2, 1, 2, 1}));
  EXPECT_EQ(*r.residual, fg.residual());
  EXPECT_EQ(*r.saturated_boundary_edges, std::vector<EdgeId> {});
  EXPECT_EQ(*r.saturated_boundary_edges, fg.saturated_boundary_edges(v({0})));
  // Nothing requested: absent.
  r = p.place(1);
  EXPECT_FALSE(r.edge_flows || r.residual || r.saturated_boundary_edges || r.min_cut);
  // CentralizedTE: the curve always, min_cut when requested and the solve is maximal (at +inf; none
  // under a cap of 2 below the maximum 8, which leaves a residual path); bottleneck alone does
  // not request it.
  FlowLedger fg2(g);
  FlowPolicy gl(fg2, v({0}), v({3}), global_config());
  r = gl.place(2, OutputOptions {false, false, true, true});
  EXPECT_EQ(r.call_delivered, 2);
  ASSERT_TRUE(r.marginal_cost_curve.has_value());
  EXPECT_EQ(*r.marginal_cost_curve, curve({{2, 2}}));
  EXPECT_EQ(gl.delivered(), 2);
  EXPECT_FALSE(r.min_cut.has_value());
  // At +inf the remaining 6 of the maximum 5 + 3: both edges out of node 0 saturate, so the
  // closure is {0} and the cut is 0 -> 1 and 0 -> 2 (EdgeIds 0 and 1).
  r = gl.place(std::nullopt, OutputOptions {false, false, true, true});
  EXPECT_EQ(r.call_delivered, 6);
  EXPECT_EQ(gl.delivered(), 8);
  ASSERT_TRUE(r.min_cut.has_value());
  EXPECT_EQ(*r.min_cut, (std::vector<EdgeId> {0, 1}));
  FlowLedger fg3(g);
  FlowPolicy gb(fg3, v({0}), v({3}), global_config());
  r = gb.place(std::nullopt, OutputOptions {false, false, true});
  EXPECT_TRUE(r.saturated_boundary_edges.has_value());
  EXPECT_FALSE(r.min_cut.has_value());
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
  const CentralizedTE third = centralized(SourceAllocation::Optimized, dag_of(*g, {x, y}));
  FlowPolicy other(fg, v({0}), v({3}), cfg(third, RouteAllocation::Optimized));
  ASSERT_EQ(other.place(1).call_delivered, 1);
  const auto el = g->eligibility({}, std::vector<EdgeId> {x});
  FlowPolicy p(fg, v({0}), v({3}), shortest_paths_config(), &el);
  const auto r = p.place(2, OutputOptions {false, false, true});
  EXPECT_EQ(r.call_delivered, 2);
  ASSERT_TRUE(r.saturated_boundary_edges.has_value());
  EXPECT_EQ(*r.saturated_boundary_edges, std::vector<EdgeId> {a});
  EXPECT_EQ(*r.saturated_boundary_edges, fg.saturated_boundary_edges(v({0}), &el));
  auto without = std::vector<EdgeId> {a, y};
  std::sort(without.begin(), without.end());
  EXPECT_EQ(fg.saturated_boundary_edges(v({0})), without);  // the Eligibility matters here
}

TEST(FlowPolicyGlobal, EmptyRosterReportsAnEmptyCurveAndAnEmptyCut) {
  // No source reaches the destination: the maximum flow is 0 and maximal.
  auto g = shared(make_graph(3, {{0, 1, 5, 1}}));
  FlowLedger fg(g);
  FlowPolicy p(fg, v({0}), v({2}), global_config());
  auto r = p.place(std::nullopt, OutputOptions {false, false, true, true});
  EXPECT_EQ(r.termination, Termination::NoRoute);
  ASSERT_TRUE(r.marginal_cost_curve.has_value());
  EXPECT_TRUE(r.marginal_cost_curve->empty());
  ASSERT_TRUE(r.min_cut.has_value());
  EXPECT_TRUE(r.min_cut->empty());
  EXPECT_TRUE(p.records().empty());
  // The zero request is a no-op: no solve, no cut; nothing is prepared, so no curve.
  FlowPolicy q(fg, v({1}), v({2}), global_config());
  r = q.place(0, OutputOptions {false, false, true, true});
  EXPECT_FALSE(r.min_cut.has_value());
  EXPECT_FALSE(r.marginal_cost_curve.has_value());
}

// CentralizedTE with a roster whose maximum flow is 0: Exhausted (a controller whose maximum is 0 over
// existing routes), whether or not something was placed before. One edge 0 -> 1 of 2 units that a
// third party fills: the source statically reaches the destination (it is on the roster), and
// place(1) finds no flow: Exhausted, its record created empty (no route, no min_route_cost, nothing
// delivered). After the third party leaves, place(1) delivers 1; a third party then takes the other
// unit, and the next place(1) re-solves 1 + 1 on a view of 1: nothing more, Exhausted.
TEST(FlowPolicyGlobal, ARosterWithoutFlowIsExhausted) {
  auto g = shared(make_graph(2, {{0, 1, 2, 1}}));
  FlowLedger fg(g);
  const CentralizedTE edge = centralized(SourceAllocation::Optimized, dag_of(*g, {0}));
  auto third = std::make_unique<FlowPolicy>(fg, v({0}), v({1}), cfg(edge, RouteAllocation::Optimized));
  ASSERT_EQ(third->place(2).call_delivered, 2);
  FlowPolicy p(fg, v({0}), v({1}), global_config());
  auto r = p.place(1);
  EXPECT_EQ(r.call_delivered, 0);
  EXPECT_EQ(r.termination, Termination::Exhausted);
  ASSERT_TRUE(r.marginal_cost_curve.has_value());
  EXPECT_TRUE(r.marginal_cost_curve->empty());
  ASSERT_EQ(p.records().size(), 1U);
  EXPECT_TRUE(ascending_edges(*p.records()[0].route).empty());
  EXPECT_EQ(p.records()[0].min_route_cost, kUnreachable);
  EXPECT_EQ(p.records()[0].delivered, 0);
  third->remove();
  r = p.place(1);
  EXPECT_EQ(r.call_delivered, 1);
  EXPECT_EQ(r.termination, Termination::Complete);
  FlowPolicy fourth(fg, v({0}), v({1}), cfg(edge, RouteAllocation::Optimized));
  ASSERT_EQ(fourth.place(1).call_delivered, 1);
  r = p.place(1);
  EXPECT_EQ(r.call_delivered, 0);
  EXPECT_EQ(r.termination, Termination::Exhausted);
}

// A CentralizedTE zero call reports the cumulative curve (its amounts sum to delivered(); None
// before the epoch's first positive call, also after remove()) and changes nothing; it solves
// nothing, so it claims no cut.
TEST(FlowPolicyGlobal, ZeroCallReportsTheCumulativeCurve) {
  using Curve = std::vector<std::pair<Cost, Amount>>;
  auto g = shared(make_graph(3, {{0, 1, 3, 7}, {0, 2, 4, 1}, {2, 1, 4, 1}}));
  FlowLedger fg(g);
  FlowPolicy p(fg, v({0}), v({1}), global_config());
  auto z = p.place(0);
  EXPECT_FALSE(z.marginal_cost_curve.has_value());
  auto r = p.place(5);
  EXPECT_EQ(*r.marginal_cost_curve, (Curve {{2, 4}, {7, 1}}));
  const auto carried = fg.carried();
  const auto requested = detail::PolicyAccess::engine(p).epoch.requested;
  z = p.place(0, OutputOptions {true, true, true, true});
  EXPECT_EQ(z.call_delivered, 0);
  EXPECT_EQ(*z.marginal_cost_curve, *r.marginal_cost_curve);
  EXPECT_FALSE(z.min_cut.has_value());
  EXPECT_EQ(fg.carried(), carried);
  EXPECT_EQ(detail::PolicyAccess::engine(p).epoch.requested, requested);
  r = p.place(1);  // cumulative 6: the solve's curve on V = W + own
  EXPECT_EQ(*r.marginal_cost_curve, (Curve {{2, 4}, {7, 2}}));
  z = p.place(0);
  EXPECT_EQ(*z.marginal_cost_curve, *r.marginal_cost_curve);
  Amount sum = 0;
  for (const auto& [c, x] : *z.marginal_cost_curve) sum += x;
  EXPECT_EQ(sum, p.delivered());
  p.remove();
  z = p.place(0);
  EXPECT_FALSE(z.marginal_cost_curve.has_value());
}

// The CentralizedTE layer on random graphs (nonuniform costs 1 to 9 times a factor of 1, 1,000,003 or
// 2^35, zero capacities, a node or edge mask, another owner's occupancy W): each call re-solves the policy's cumulative amount on the view
// W + own. After each of up to two capped calls and after place(+inf), the policy's own flow is
// feasible on the view, conserved, carries the cumulative amount at the minimum cost of that
// amount (the independent Bellman-Ford curve, not the flow core), its records sum to it, and its
// cumulative curve is the optimal curve; at +inf the reported cut separates the sources from the
// destinations with a capacity equal to the flow. evaluate_demand() on a fresh ledger reaches the same
// maximum and curve, and evaluate_pairs() agrees with evaluate_demand() per pair.
TEST(FlowPolicyGlobal, RandomCumulativeSolvesAreOptimalWithACertifiedCut) {
  std::mt19937_64 rng(20261010);
  const OutputOptions all {true, true, true, true};
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
    // Another owner: a CentralizedTE policy between other endpoints (when the graph allows them) takes
    // part of the capacity first.
    std::optional<FlowPolicy> other;
    if (n >= 6 && rng() % 2 == 0) {
      other.emplace(fg, v({2}), v({3}), global_config());
      occupied += other->place(static_cast<Amount>(rng() % 5)).call_delivered > 0 ? 1 : 0;
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
        for (const auto& [e, x] : rec.edge_flows) own[static_cast<std::size_t>(e)] += x;
      }
      EXPECT_EQ(recorded, total) << at;
      EXPECT_EQ(p.delivered(), total) << at;
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
      ASSERT_TRUE(r.marginal_cost_curve.has_value()) << at;
      const auto pts = curve_shape_clauses(*r.marginal_cost_curve, total, cost, fails);
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
      EXPECT_EQ(r.call_delivered, x) << where;
      EXPECT_FALSE(r.min_cut.has_value()) << where;  // stopped by the cap
      check_own(total, r, where + " capped call " + std::to_string(call));
    }
    const auto r = p.place(std::nullopt, all);
    EXPECT_EQ(r.call_delivered, max_flow - total) << where;
    with_flow += max_flow > 0 ? 1 : 0;
    check_own(max_flow, r, where + " +inf");
    // The cut, as the flow certificates define it (residual_closure and closure_cut,
    // test_utils.hpp): the residual closure of the roster sources is the nodes they still reach
    // over allowed edges, forward where the policy's own flow leaves room in the view and backward
    // where it can be undone; no destination is in it, and the cut is every allowed edge that
    // leaves it, whatever its view, in ascending EdgeId, with a capacity equal to the maximum flow.
    ASSERT_TRUE(r.min_cut.has_value() && r.edge_flows.has_value()) << where;
    Amount cut_capacity = 0;
    for (const EdgeId e : *r.min_cut) cut_capacity += view[static_cast<std::size_t>(e)];
    EXPECT_EQ(cut_capacity, max_flow) << where;
    // The closure is seeded at the roster: the allowed sources with a walk over allowed edges
    // (whatever their capacity) to an allowed destination. Another source can reach no
    // destination, so it needs no cut edge; the cut of an empty roster is [].
    std::vector<NodeId> roster;
    for (const NodeId src : srcs) {
      if (elig != nullptr && !elig->allows_node(src)) continue;
      std::vector<char> reach(static_cast<std::size_t>(g->num_nodes()), 0);
      std::vector<NodeId> stack = {src};
      reach[static_cast<std::size_t>(src)] = 1;
      while (!stack.empty()) {
        const NodeId u = stack.back();
        stack.pop_back();
        for (EdgeId e = 0; e < g->num_edges(); ++e) {
          const NodeId w = g->edge_destinations()[static_cast<std::size_t>(e)];
          if (g->edge_sources()[static_cast<std::size_t>(e)] != u || reach[static_cast<std::size_t>(w)] != 0) continue;
          if (elig != nullptr && !elig->allows_edge(e)) continue;
          reach[static_cast<std::size_t>(w)] = 1;
          stack.push_back(w);
        }
      }
      if (std::any_of(dsts.begin(), dsts.end(), [&](NodeId d) {
            return (elig == nullptr || elig->allows_node(d)) && reach[static_cast<std::size_t>(d)] != 0;
          })) {
        roster.push_back(src);
      }
    }
    EXPECT_EQ(roster, detail::PolicyAccess::engine(p).reachable_sources) << where;
    const auto closure = residual_closure(*g, view, *r.edge_flows, roster, elig);
    for (const NodeId d : dsts) EXPECT_EQ(closure[static_cast<std::size_t>(d)], 0) << where << " destination " << d;
    EXPECT_EQ(*r.min_cut, closure_cut(*g, closure, elig)) << where;
    // The evaluators on a fresh ledger (static capacities): the same maximum and cost.
    if (other.has_value()) continue;
    ++evaluated;
    const EvaluationOptions co {global_config(), elig, nullptr, all};
    const auto cr = evaluate_demand(*g, srcs, dsts, co);
    EXPECT_EQ(cr.call_delivered, max_flow) << where;
    ASSERT_TRUE(cr.marginal_cost_curve.has_value()) << where;
    EXPECT_EQ(*cr.marginal_cost_curve, *r.marginal_cost_curve) << where;
    const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {{srcs, dsts}, {v({0}), v({n - 1})}};
    const auto batch = evaluate_pairs(*g, pairs, co, {}, 1);
    ASSERT_EQ(batch.size(), 2U) << where;
    for (std::size_t i = 0; i < pairs.size(); ++i) {
      const auto one = evaluate_demand(*g, pairs[i].first, pairs[i].second, co);
      EXPECT_EQ(batch[i].call_delivered, one.call_delivered) << where << " pair " << i;
      EXPECT_EQ(batch[i].marginal_cost_curve, one.marginal_cost_curve) << where << " pair " << i;
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

// ShortestPaths placements on random graphs (costs 1 to 9 times a factor of 1, 1,000,003 or 2^35, zero
// capacities, a node or edge mask, one or two sources and destinations, another owner's flow
// first), in the kernels and admissions ShortestPaths takes: after each of two capped calls and a call
// of +inf, each of the policy's records is a flow of its own from the policy's sources to its
// destinations (expect_conserved, test_utils.hpp: positive amounts in ascending EdgeId, every
// transit node passing on what it receives less its drops, delivered the destinations' net inflow
// and, Lossless, the sources' net outflow), carrying and dropping nothing on an edge the Eligibility
// excludes; a Lossless policy's records deliver its delivered() total, and the ledger's aggregate is the
// sum of its records (expect_balanced). ShortestPaths chooses routes by cost, not by the minimum-cost
// flow, so no optimality certificate applies here (T1's records: T1Property, distributed_te_tests.cpp).
TEST(FlowPolicyRecords, RandomShortestPathsRecordsAreFlows) {
  const auto res = sel(RouteTieBreak::FewerHops);
  const std::vector<FlowPolicyConfig> configs = {
      shortest_paths_config(),
      cfg(ShortestPaths {}, RouteAllocation::Weighted),
      cfg(ShortestPaths {}, RouteAllocation::Optimized),
      // The FEWER_HOPS single path (several sources reject a HIGHER_* tie_break with it).
      cfg(ShortestPaths {}, RouteAllocation::Equal, Admission::Lossless,
          sel(RouteTieBreak::FewerHops, false, false)),
      cfg(ShortestPaths {}, RouteAllocation::Equal, Admission::Lossy),
      cfg(ShortestPaths {}, RouteAllocation::Weighted, Admission::Lossy, res),
  };
  std::mt19937_64 rng(20261011);
  auto below = [&](std::int64_t k) { return static_cast<std::int32_t>(rng() % static_cast<std::uint64_t>(k)); };
  int carrying = 0;                 // record checks of a record that carries something
  int dropping = 0;                 // record checks of a Lossy record with drops
  int masked = 0;                   // record checks of a carrying record under a mask
  int scaled = 0;                   // record checks of a carrying record with costs times 1,000,003 or 2^35
  int several = 0;                  // calls that delivered from two sources or into two destinations
  for (int trial = 0; trial < 600; ++trial) {
    const FlowPolicyConfig& config = configs[static_cast<std::size_t>(trial) % configs.size()];
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
      other.emplace(fg, v({2}), v({3}), shortest_paths_config());
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
        for (const auto& [e, x] : r.edge_flows) {
          EXPECT_TRUE(elig == nullptr || elig->allows_edge(e)) << at << " carries " << x << " on excluded edge " << e;
        }
        for (const auto& [e, x] : r.edge_drops) {
          EXPECT_TRUE(elig == nullptr || elig->allows_edge(e)) << at << " drops on excluded edge " << e;
        }
        delivered += r.delivered;
        if (r.delivered > 0) {
          from.insert(r.sources.begin(), r.sources.end());
          for (const auto& [e, x] : r.edge_flows) {
            const NodeId head = g->edge_destinations()[static_cast<std::size_t>(e)];
            if (std::find(dsts.begin(), dsts.end(), head) != dsts.end()) into.insert(head);
          }
        }
        carrying += carries ? 1 : 0;
        dropping += r.edge_drops.empty() ? 0 : 1;
        masked += carries && elig != nullptr ? 1 : 0;
        scaled += carries && factor > 1 ? 1 : 0;
      }
      if (config.admission == Admission::Lossless) EXPECT_EQ(delivered, p.delivered()) << where;
      several += from.size() > 1 || into.size() > 1 ? 1 : 0;
      ASSERT_FALSE(::testing::Test::HasFailure()) << where;
    }
  }
  EXPECT_GT(carrying, 0);
  EXPECT_GT(dropping, 0);
  EXPECT_GT(masked, 0);
  EXPECT_GT(scaled, 0);
  EXPECT_GT(several, 0);
}

TEST(FlowPolicyRecords, BundlesAreSharedOrPruned) {
  auto g = diamond(5, 5);
  FlowLedger fg(g);
  const auto whole = dag_of(*g, {0, 1, 2, 3});
  const Pinned pin = pinned_dags(*g, 0, 3, {whole});
  FlowPolicy p(fg, v({0}), v({3}), {.routing = pin});
  (void)p.place(2);
  // An unchanged supplied split route is shared by the record (no copy).
  EXPECT_EQ(p.records().front().route.get(), whole.get());
  // A route pruned by the replay's mask is a new DAG of what remains.
  const auto el = g->eligibility({}, std::vector<EdgeId> {2});
  FlowLedger fg2(g);
  FlowPolicy q(fg2, v({0}), v({3}), {.routing = pin}, &el);
  (void)q.place(2);
  const auto rec = q.records().front();
  EXPECT_NE(rec.route.get(), whole.get());
  EXPECT_EQ(ascending_edges(*rec.route), (std::vector<EdgeId> {1, 3}));
  // ShortestPaths: a source whose route is the whole cached support shares it.
  ShortestPathCache cache(*g);
  FlowLedger fg3(g);
  FlowPolicy c(fg3, v({0}), v({3}), shortest_paths_config(), nullptr, &cache);
  (void)c.place(2);
  const auto support = cache.lookup(*g, v({3}), Direction::Reverse, nullptr);
  EXPECT_EQ(c.records().front().route.get(), support.get());
}

// records() lists the live records in FlowId order, also when the sources' entries interleave.
// Sources 0 and 1, each with two parallel edges of 3 and 5 units to node 2; T1 with 2 tunnels per
// head: one record per program entry, created in index order j * H + i (round j, head i), so the
// entries 0 (source 0), 1 (source 1), 2 (source 0) and 3 (source 1) get FlowIds 0 to 3 and records()
// gives their sources as 0, 1, 0, 1 (a source's records are not contiguous).
TEST(FlowPolicyRecords, RecordsAreInFlowIdOrderAcrossSources) {
  auto g = shared(make_graph(3, {{0, 2, 3, 1}, {0, 2, 5, 1}, {1, 2, 3, 1}, {1, 2, 5, 1}}));
  FlowLedger fg(g);
  FlowPolicy p(fg, v({0, 1}), v({2}), t1(2));
  ASSERT_EQ(p.replace(8, 8).call_delivered, 8);
  std::vector<std::int64_t> ids;
  std::vector<NodeId> sources;
  std::vector<std::int32_t> entries;
  for (const auto& r : p.records()) {
    ids.push_back(r.flow_id.second);
    sources.push_back(r.sources.at(0));
    entries.push_back(r.entry_id.value_or(-1));
  }
  EXPECT_EQ(ids, (std::vector<std::int64_t> {0, 1, 2, 3}));
  EXPECT_EQ(sources, (std::vector<NodeId> {0, 1, 0, 1}));
  EXPECT_EQ(entries, (std::vector<std::int32_t> {0, 1, 2, 3}));
}

// ShortestPaths routes to the nearest destinations, and one flow may deliver to several of them: edges (input order) 0 -> 1 (10 units, cost 1), 1 -> 2 (2, cost 1), 1 -> 3 (3,
// cost 1) and 0 -> 3 (5, cost 5); destinations 2 and 3, both at cost 2 through node 1. The one
// flow carries 2 to node 2 and 3 to node 3 (the edges entering each destination), 5 in all; the
// costlier edge 0 -> 3 is no part of the route.
TEST(FlowPolicyRecords, OneFlowDeliversToTheNearestDestinations) {
  auto g = shared(make_graph(4, {{0, 1, 10, 1}, {1, 2, 2, 1}, {1, 3, 3, 1}, {0, 3, 5, 5}}));
  const auto id = id_of_input(*g);
  // Optimized fills both; Weighted's capacity weights 2 and 3 split 5 at node 1 as (2, 3).
  for (const auto& [name, kernel] : std::vector<std::pair<std::string, RouteAllocation>> {
           {"Optimized", RouteAllocation::Optimized}, {"Weighted", RouteAllocation::Weighted}}) {
    FlowLedger fg(g);
    FlowPolicy p(fg, v({0}), v({2, 3}), cfg(ShortestPaths {}, kernel));
    const auto r = p.place(std::nullopt);
    EXPECT_EQ(r.call_delivered, 5) << name;
    const auto records = p.records();
    ASSERT_EQ(records.size(), 1U) << name;
    const auto& rec = records[0];
    EXPECT_EQ(rec.delivered, 5) << name;
    EXPECT_EQ(rec.min_route_cost, 2) << name;
    EXPECT_EQ(rec.destinations, (std::vector<NodeId> {2, 3})) << name;
    std::map<NodeId, Amount> into;
    for (const auto& [e, x] : rec.edge_flows) into[g->edge_destinations()[static_cast<std::size_t>(e)]] += x;
    EXPECT_EQ(into, (std::map<NodeId, Amount> {{1, 5}, {2, 2}, {3, 3}})) << name;
    EXPECT_EQ(ascending_edges(*rec.route), (std::vector<EdgeId> {id[0], id[1], id[2]})) << name;
  }
}

TEST(PredDAGStorage, ABundleCostsItsEntriesNotTheGraph) {
  const std::int32_t n = 1'000'000;
  auto g = shared(make_graph(n, {{n - 3, n - 2, 4, 1}, {n - 2, n - 1, 4, 1}}));
  const auto dag = dag_of(*g, {0, 1});
  EXPECT_LT(detail::PredecessorDAGBuilder::payload_bytes(*dag), 128U);  // two entries, not N + 1 offsets
  const auto off = dag->parent_offsets();  // the from_csr shape, on demand
  ASSERT_EQ(off.size(), static_cast<std::size_t>(n) + 1);
  EXPECT_EQ(off[static_cast<std::size_t>(n) - 2], 0);
  EXPECT_EQ(off[static_cast<std::size_t>(n) - 1], 1);
  EXPECT_EQ(off[static_cast<std::size_t>(n)], 2);
  EXPECT_TRUE(PredecessorDAG::from_csr(*g, off, dag->parents(), dag->via_edges()) == *dag);
  FlowLedger fg(g);
  FlowPolicy p(fg, v({n - 3}), v({n - 1}), {.routing = pinned_dags(*g, n - 3, n - 1, {dag})});
  EXPECT_EQ(p.place(3).call_delivered, 3);
  EXPECT_LT(detail::PredecessorDAGBuilder::payload_bytes(*p.records().front().route), 128U);
  // The record routes the engine builds itself are sparse too: a supplied split route pruned by
  // the replay's mask (the bypass n - 3 -> n - 1 excluded), CentralizedTE's used edges on a support, and
  // the ShortestPaths and CentralizedTE routes.
  auto h = shared(make_graph(n, {{n - 3, n - 2, 4, 1}, {n - 2, n - 1, 4, 1}, {n - 3, n - 1, 4, 5}}));
  const auto inv = id_of_input(*h);
  const auto three = dag_of(*h, {inv[0], inv[1], inv[2]});
  const auto el = h->eligibility({}, std::vector<EdgeId> {inv[2]});
  const std::vector<std::pair<const char*, FlowPolicyConfig>> built = {
      {"pruned Pinned", {.routing = pinned_dags(*h, n - 3, n - 1, {three})}},
      {"CentralizedTE support", {.routing = centralized(SourceAllocation::Optimized, dag_of(*h, {inv[0], inv[1]}))}},
      {"ShortestPaths", shortest_paths_config()},
      {"CentralizedTE", global_config()}};
  for (const auto& [name, c] : built) {
    FlowLedger fh(h);
    FlowPolicy q(fh, v({n - 3}), v({n - 1}), c, std::holds_alternative<Pinned>(c.routing) ? &el : nullptr);
    EXPECT_EQ(q.place(3).call_delivered, 3) << name;
    const auto route = q.records().front().route;
    EXPECT_NE(route.get(), three.get()) << name;
    EXPECT_EQ(ascending_edges(*route), (std::vector<EdgeId> {inv[0], inv[1]})) << name;
    EXPECT_LT(detail::PredecessorDAGBuilder::payload_bytes(*route), 128U) << name;
  }
}

TEST(Evaluators, CapacityIsPlaceInfOnAFreshLedger) {
  // Every field of the result, in every routing model with incremental placement (Pinned's
  // program has the endpoints; T1's evaluate of its design: TunnelEvaluators.T1Capacity).
  auto g = diamond(5, 3, 2);
  const Pinned pin = pinned_dags(*g, 0, 3, {dag_of(*g, {0, 2}), dag_of(*g, {1, 3})});
  for (const auto& c : {global_config(), shortest_paths_config(), FlowPolicyConfig {.routing = pin},
                        FlowPolicyConfig {.routing = centralized(SourceAllocation::Optimized, dag_of(*g, {0, 2}))}}) {
    const auto r = evaluate_demand(*g, v({0}), v({3}), EvaluationOptions {c, nullptr, nullptr, OutputOptions {true, true, true}});
    FlowLedger fg(g);
    FlowPolicy p(fg, v({0}), v({3}), c);
    const auto want = p.place(std::nullopt, OutputOptions {true, true, true});
    EXPECT_GT(want.call_delivered, 0);  // two empty results would agree too
    expect_same(r, want, "capacity against place(+inf)");
  }
  // The exact maximum flow: CentralizedTE's curve (cheapest walks first).
  const auto exact = evaluate_demand(*g, v({0}), v({3}), EvaluationOptions {global_config()});
  EXPECT_EQ(exact.call_delivered, 8);
  EXPECT_EQ(*exact.marginal_cost_curve, curve({{2, 5}, {4, 3}}));
}

TEST(Evaluators, BatchOrderIntersectionsExceptionsAndWorkers) {
  auto g = diamond(5, 3);
  const auto common = g->eligibility({}, std::vector<EdgeId> {3});       // route 2's last edge
  const auto first_route = g->eligibility({}, std::vector<EdgeId> {0});  // route 1's first edge
  std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {{0}, {3}}, {{0}, {3}}, {{1}, {3}}, {{0}, {1, 2}}};
  const std::vector<const Eligibility*> per = {nullptr, &first_route, nullptr, nullptr};
  const EvaluationOptions opts {global_config(), &common, nullptr, OutputOptions {true, true, true, true}};
  const auto r1 = evaluate_pairs(*g, pairs, opts, per, 1);
  ASSERT_EQ(r1.size(), 4U);
  EXPECT_EQ(r1[0].call_delivered, 5);  // common only: route 1
  EXPECT_EQ(r1[1].call_delivered, 0);  // common and the pair's: both routes cut
  EXPECT_EQ(r1[2].call_delivered, 5);
  EXPECT_EQ(r1[3].call_delivered, 8);  // destination set {1, 2}: both first edges (edge 3 excluded is past 2)
  EXPECT_EQ(*r1[3].marginal_cost_curve, curve({{1, 8}}));
  EXPECT_EQ(*r1[3].min_cut, (std::vector<EdgeId> {0, 1}));
  // Every field of every result is the same for any worker count.
  for (std::int32_t w : {2, 3, 8}) {
    const auto rw = evaluate_pairs(*g, pairs, opts, per, w);
    ASSERT_EQ(rw.size(), 4U);
    for (std::size_t i = 0; i < 4; ++i) {
      expect_same(rw[i], r1[i], "workers " + std::to_string(w) + " pair " + std::to_string(i));
    }
  }
  EXPECT_TRUE(evaluate_pairs(*g, {}, opts, {}, 4).empty());
  // Validation before any worker: endpoints, lengths, foreign values, Pinned, workers, the
  // building cache's one writer.
  auto bad = pairs;
  bad[2].second = {9};
  EXPECT_THROW((void)evaluate_pairs(*g, bad, opts, per, 4), std::invalid_argument);
  bad = pairs;
  bad[1].first = {3};
  EXPECT_TRUE(contains(error_of([&] { (void)evaluate_pairs(*g, bad, opts, per, 4); }), "both endpoint sets"));
  EXPECT_THROW((void)evaluate_pairs(*g, pairs, opts, std::vector<const Eligibility*> {nullptr}, 1), std::invalid_argument);
  auto other = diamond(5, 3);
  const auto foreign = other->eligibility({}, std::vector<EdgeId> {0});
  const std::vector<const Eligibility*> per_foreign = {nullptr, &foreign, nullptr, nullptr};
  EXPECT_TRUE(contains(error_of([&] { (void)evaluate_pairs(*g, pairs, opts, per_foreign, 1); }), "foreign"));
  const Pinned pin = pinned_dags(*g, 0, 3, {dag_of(*g, {0, 2})});
  EXPECT_EQ(error_of([&] { (void)evaluate_pairs(*g, pairs, EvaluationOptions {{.routing = pin}}, {}, 1); }),
            "evaluate_pairs does not accept Pinned");
  EXPECT_TRUE(contains(error_of([&] { (void)evaluate_pairs(*g, pairs, opts, per, 0); }), "workers"));
  ShortestPathCache cache(*g);
  const EvaluationOptions with_cache {shortest_paths_config(), nullptr, &cache};
  EXPECT_TRUE(contains(error_of([&] { (void)evaluate_pairs(*g, pairs, with_cache, {}, 2); }), "workers == 1"));
  const auto cached = evaluate_pairs(*g, pairs, with_cache, {}, 1);
  EXPECT_EQ(cache.stats().misses, 2);  // destination sets {3} and {1, 2}
  EXPECT_EQ(cache.stats().hits, 2);
  EXPECT_EQ(cached[0].call_delivered, 7);  // ShortestPaths EB over routes of 5 and 3: [4, 3] (odd unit to the lower id)
  // The argument checks in their documented order, all before the pairs: with every defect
  // present the first is reported, and with that one removed the next. Pinned also carries a
  // configuration error (an admission outside its enumeration); pair 2 has a destination out of range. The building cache
  // is untouched: no lookup ran.
  ShortestPathCache foreign_cache(*other);
  ShortestPathCache building(*g);
  const FlowPolicyConfig pinned_bad {.routing = pin, .admission = static_cast<Admission>(97)};
  const FlowPolicyConfig config_bad = cfg(ShortestPaths {}, RouteAllocation::Optimized, Admission::Lossy);
  bad = pairs;
  bad[2].second = {9};
  const std::vector<const Eligibility*> short_per = {nullptr};
  auto first_error = [&](const FlowPolicyConfig& c, const Eligibility* elig, ShortestPathCache* with,
                         const std::vector<const Eligibility*>& per_pair, std::int32_t workers) {
    return error_of([&] { (void)evaluate_pairs(*g, bad, EvaluationOptions {c, elig, with}, per_pair, workers); });
  };
  EXPECT_EQ(first_error(pinned_bad, &foreign, &foreign_cache, short_per, 0), "workers must be >= 1 (got 0)");
  EXPECT_EQ(first_error(pinned_bad, &foreign, &foreign_cache, short_per, 2),
            "foreign Eligibility: built for another graph");
  EXPECT_EQ(first_error(pinned_bad, nullptr, &foreign_cache, short_per, 2), "foreign ShortestPathCache: bound to another graph");
  EXPECT_EQ(first_error(pinned_bad, nullptr, &building, short_per, 2),
            "evaluate_pairs does not accept Pinned");  // arguments before cache mode
  EXPECT_EQ(first_error(pinned_bad, nullptr, &building, short_per, 1), "evaluate_pairs does not accept Pinned");
  EXPECT_EQ(first_error(config_bad, nullptr, &building, short_per, 1), "ShortestPaths OPTIMIZED requires admission=LOSSLESS");
  EXPECT_EQ(first_error(global_config(), nullptr, &building, short_per, 1), "CentralizedTE does not use ShortestPathCache");
  EXPECT_EQ(first_error(shortest_paths_config(), nullptr, &building, short_per, 1), "per_pair_eligibility has 1 entries for 4 pairs");
  EXPECT_EQ(first_error(shortest_paths_config(), nullptr, &building, per, 1), "destinations: node 9 out of node range [0, 4)");
  EXPECT_EQ(stats_of(building.stats()), stats_of(ShortestPathCache(*g).stats()));
}

TEST(Evaluators, BatchValidatesInPairOrderBeforeAnyWorker) {
  // Two invalid pairs (1 and 3): the batch raises pair 1's error before any worker starts.
  auto g = diamond(5, 3);
  std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs = {
      {{0}, {3}}, {{0}, {7}}, {{0}, {3}}, {{0}, {8}}};
  const auto msg = error_of([&] { (void)evaluate_pairs(*g, pairs, EvaluationOptions {global_config()}, {}, 4); });
  EXPECT_TRUE(contains(msg, "node 7")) << msg;
  // Within one pair: the sources, then the destinations, then their overlap, then the pair's
  // Eligibility's graph.
  const auto elsewhere = diamond(5, 3);
  const auto foreign = elsewhere->eligibility({}, std::vector<EdgeId> {0});
  const std::vector<const Eligibility*> per = {&foreign};
  const auto pair_error = [&](std::vector<NodeId> s, std::vector<NodeId> d) {
    const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> one = {{s, d}};
    return error_of([&] { (void)evaluate_pairs(*g, one, EvaluationOptions {global_config()}, per, 1); });
  };
  EXPECT_EQ(pair_error({}, {9}), "sources: the node set is empty");
  EXPECT_EQ(pair_error({0, 3}, {3, 9}), "destinations: node 9 out of node range [0, 4)");
  EXPECT_EQ(pair_error({0, 3}, {3}), "node 3 is in both endpoint sets: sources and destinations must not overlap");
  EXPECT_EQ(pair_error({0}, {3}), "foreign Eligibility: per_pair_eligibility[0] built for another graph");
}

// The failure study's shortcuts (optimized CentralizedTE from the baseline's maximum flow, ShortestPaths without a
// re-evaluation of the edges outside the static shortest-path support reachable from the sources)
// and the plain re-evaluations (equal CentralizedTE, T1, Pinned) give the losses of the
// definition, the
// baseline's call_delivered minus evaluate_demand() with the edge excluded, on random graphs with parallel
// edges, self-loops, dispersed costs, zero capacities, endpoint sets and eligibility masks; every
// ShortestPaths kernel under FEWER_HOPS and HIGHER_CAPACITY with every multipath and parallel_edges
// setting the tie rule allows (the pooled single-path rule skips HIGHER_* single paths with two
// sources), both admissions, one HIGHER_AVAILABLE_CAPACITY single path (no HIGHER_HEADROOM) and the
// equal source allocation of Optimized; CentralizedTE under both source allocations (equal: a deletion that removes a
// source from the roster can raise the total, a negative loss); two T1 configurations with the
// design 6 (HIGHER_AVAILABLE_CAPACITY without parallel links, and three tunnels under a CostStretch ceiling: setup redone under
// every deletion); and Pinned
// with two route DAGs (the static support toward the destination and the second shortest path) on
// single endpoints. An edge without an entry has loss 0 by the same definition.
TEST(Evaluators, SensitivityEqualsTheDefinitionOnRandomGraphs) {
  std::mt19937 rng(20261005);
  auto pick = [&](int n) { return static_cast<int>(rng() % static_cast<unsigned>(n)); };
  std::vector<FlowPolicyConfig> configs = {global_config(), {.routing = centralized(SourceAllocation::Equal)}};
  for (auto k : {RouteAllocation::Optimized, RouteAllocation::Equal, RouteAllocation::Weighted}) {
    for (bool multipath : {true, false}) {
      for (bool parallel_edges : {true, false}) {
        for (auto t : {RouteTieBreak::FewerHops, RouteTieBreak::HigherCapacity}) {
          if (t != RouteTieBreak::FewerHops && multipath && parallel_edges) continue;  // the ineffective tie rule
          configs.push_back(cfg(ShortestPaths {}, k, Admission::Lossless, sel(t, multipath, parallel_edges)));
        }
      }
    }
  }
  configs.push_back(cfg(ShortestPaths {}, RouteAllocation::Equal, Admission::Lossy));
  configs.push_back(cfg(ShortestPaths {}, RouteAllocation::Weighted, Admission::Lossy, sel(RouteTieBreak::HigherAvailableCapacity, false, false)));
  {
    ShortestPaths equal_controller = shortest_paths(RouteAllocation::Optimized);
    equal_controller.source_allocation = SourceAllocation::Equal;
    configs.push_back({.routing = equal_controller});
  }
  configs.push_back({.routing = DistributedTE {.tunnel_count = TunnelCount {2},
                                               .tie_break = RouteTieBreak::HigherAvailableCapacity,
                                               .parallel_edges = false}});
  configs.push_back({.routing = DistributedTE {TunnelCount {3}, CostStretch(2)}});
  // T1 names its design; the other models take none.
  const auto query = [](const FlowPolicyConfig& c, const Eligibility* e, OutputOptions o = {}) {
    EvaluationOptions q {c, e, nullptr, o};
    if (std::holds_alternative<DistributedTE>(c.routing)) q.design_demand = Amount {6};
    return q;
  };
  int entries = 0;
  int global_losses = 0;
  int policy_losses = 0;
  int policy_zero = 0;
  int aware_losses = 0;
  int pinned_losses = 0;
  int pinned_trials = 0;
  int zero_edge_losses = 0;  // entries of zero-capacity edges whose deletion changes the result
  int negative = 0;          // equal-source-allocation deletions that raise the total
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
    const auto support = std::make_shared<const PredecessorDAG>(
        detail::static_support(*g, std::vector<NodeId> {t0}, Direction::Reverse, nullptr));
    const auto paths = ksp(*g, 0, t0, KspOptions(2));
    if (paths.size() == 2) {
      trial_configs.push_back({.routing = pinned_dags(*g, 0, t0,
                                                      {support, std::make_shared<const PredecessorDAG>(paths[1].second.to_dag())},
                                                      RouteAllocation::Weighted)});
      ++pinned_trials;
    }
    for (const auto& c : trial_configs) {
      const bool shortest_paths = std::holds_alternative<ShortestPaths>(c.routing);
      if (shortest_paths && s.size() > 1 && !std::get<ShortestPaths>(c.routing).multipath &&
          std::get<ShortestPaths>(c.routing).tie_break != RouteTieBreak::FewerHops) {
        continue;  // the pooled single-path rule
      }
      const bool global = std::holds_alternative<CentralizedTE>(c.routing) &&
                          std::get<CentralizedTE>(c.routing).source_allocation != SourceAllocation::Equal;
      const bool pinned = std::holds_alternative<Pinned>(c.routing);
      const std::vector<NodeId> cs = pinned ? std::vector<NodeId> {0} : s;
      const std::vector<NodeId> ct = pinned ? std::vector<NodeId> {t0} : t;
      const auto r = evaluate_edge_failures(*g, cs, ct, query(c, elig), 1 + pick(3));
      // The baseline is evaluate_demand() itself, and the entries are the candidates in ascending order:
      // CentralizedTE, the baseline's positive-flow edges; every other model, every allowed edge,
      // whatever its capacity.
      const auto base = evaluate_demand(*g, cs, ct, query(c, elig, OutputOptions {true, false, false}));
      ASSERT_EQ(r.baseline.call_delivered, base.call_delivered) << "trial " << trial;
      std::vector<EdgeId> candidates;
      for (EdgeId e = 0; e < g->num_edges(); ++e) {
        const auto ei = static_cast<std::size_t>(e);
        if (global ? (*base.edge_flows)[ei] > 0 : elig == nullptr || elig->allows_edge(e)) {
          candidates.push_back(e);
        }
      }
      const bool equal_global = std::holds_alternative<CentralizedTE>(c.routing) && !global;
      const std::string model = global          ? " CentralizedTE"
                                : equal_global  ? " CentralizedTE equal"
                                : shortest_paths     ? " ShortestPaths"
                                : pinned        ? " Pinned"
                                                : " T1";
      std::vector<EdgeId> listed;
      for (const auto& e : r.entries) listed.push_back(e.edge);
      ASSERT_EQ(listed, candidates) << "trial " << trial << model;
      // The loss of deleting `e`, by the definition.
      const auto loss = [&](EdgeId e) {
        const EdgeId one[1] = {e};
        auto deleted = g->eligibility({}, one);
        if (elig != nullptr) deleted = Eligibility::intersect(*elig, deleted);
        return base.call_delivered - evaluate_demand(*g, cs, ct, query(c, &deleted)).call_delivered;
      };
      for (const auto& e : r.entries) {
        const Amount want = loss(e.edge);
        ASSERT_EQ(e.delivery_reduction, want) << "trial " << trial << " edge " << e.edge << model;
        ++entries;
        if (want != 0) ++(global ? global_losses : shortest_paths ? policy_losses : pinned ? pinned_losses : aware_losses);
        if (shortest_paths && want == 0) ++policy_zero;
        if (g->edge_capacity_units()[static_cast<std::size_t>(e.edge)] == 0 && want != 0) ++zero_edge_losses;
        if (want < 0) ++negative;
      }
      for (EdgeId e = 0; e < g->num_edges(); ++e) {
        if (std::find(listed.begin(), listed.end(), e) != listed.end()) continue;
        EXPECT_EQ(loss(e), 0) << "trial " << trial << " edge " << e << model << ": an edge without an entry";
        ++absent;
      }
    }
  }
  // Not vacuous: many entries, losses under every model, zero-loss ShortestPaths deletions, Pinned
  // trials, deletions of zero-capacity edges that change the result, and edges without an entry.
  EXPECT_GT(entries, 3000);
  EXPECT_GT(global_losses, 100);
  EXPECT_GT(policy_losses, 500);
  EXPECT_GT(policy_zero, 500);
  EXPECT_GT(aware_losses, 100);
  EXPECT_GT(pinned_losses, 30);
  EXPECT_GT(pinned_trials, 20);
  EXPECT_GT(zero_edge_losses, 50);
  EXPECT_GT(negative, 10);
  EXPECT_GT(absent, 500);
}

TEST(Evaluators, SensitivityCandidatesAndSignedLosses) {
  // CentralizedTE: candidates are the baseline's positive-flow edges; an unsaturated critical edge
  // appears: 0 -> 1 of capacity 9 carries 4 (its tail 1 -> 2 holds 4) and its deletion
  // loses 4; a candidate set of saturated edges only would miss it.
  auto g = shared(make_graph(3, {{0, 1, 9, 1}, {1, 2, 4, 1}, {0, 2, 2, 5}}));
  auto s = evaluate_edge_failures(*g, v({0}), v({2}), EvaluationOptions {global_config()}, 1);
  EXPECT_EQ(s.baseline.call_delivered, 6);
  EXPECT_FALSE(s.baseline.edge_flows.has_value());  // requested internally, not returned
  std::map<EdgeId, Amount> losses;
  for (const auto& e : s.entries) losses[e.edge] = e.delivery_reduction;
  const auto inv = id_of_input(*g);
  EXPECT_EQ(losses.at(inv[0]), 4);  // unsaturated (4 of 9) and critical
  EXPECT_EQ(losses.at(inv[1]), 4);
  EXPECT_EQ(losses.at(inv[2]), 2);
  EXPECT_EQ(losses.size(), 3U);
  // A policy model: every allowed edge, whatever its capacity. Equal over parallel
  // links of 1, 10 and 0 (ids 0, 1, 2): round robin over the three members fits 2 units and the
  // third would cross the zero link. Without the thin link (members 1, 2): 1, loss 1; without the
  // wide link (0, 2): 1, loss 1; without the zero link (0, 1): 2, loss 0.
  auto p = shared(make_graph(2, {{0, 1, 1, 1}, {0, 1, 10, 1}, {0, 1, 0, 1}}));
  const auto opts = EvaluationOptions {shortest_paths_config()};
  s = evaluate_edge_failures(*p, v({0}), v({1}), opts, 1);
  EXPECT_EQ(s.baseline.call_delivered, 2);
  // The exact (edge, loss) list in ascending EdgeId.
  const auto pin = id_of_input(*p);
  std::vector<std::pair<EdgeId, Amount>> got;
  for (const auto& e : s.entries) got.emplace_back(e.edge, e.delivery_reduction);
  std::vector<std::pair<EdgeId, Amount>> want = {{pin[0], 1}, {pin[1], 1}, {pin[2], 0}};
  std::sort(want.begin(), want.end());
  EXPECT_EQ(got, want);
  // Workers give the same entries; an excluded edge is no candidate.
  const auto s4 = evaluate_edge_failures(*p, v({0}), v({1}), opts, 4);
  std::vector<std::pair<EdgeId, Amount>> got4;
  for (const auto& e : s4.entries) got4.emplace_back(e.edge, e.delivery_reduction);
  EXPECT_EQ(got4, want);
  const auto el = p->eligibility({}, std::vector<EdgeId> {pin[0]});
  const auto se = evaluate_edge_failures(*p, v({0}), v({1}), EvaluationOptions {opts.config, &el}, 1);
  ASSERT_EQ(se.entries.size(), 2U);
  EXPECT_EQ(se.entries[0].edge, std::min(pin[1], pin[2]));
  EXPECT_EQ(se.entries[1].edge, std::max(pin[1], pin[2]));
  // A loss can be negative, also for a zero-capacity link (the search graph S: ids 0: 0 -> 1 of
  // capacity 0, 1: 0 -> 2, 2: 2 -> 1 of 5, costs 1): the only shortest route crosses the zero
  // link, so Equal carries 0; without it the route 0 -> 2 -> 1 carries 5, a loss of -5.
  // The other two edges lie outside the reachable support (loss 0 without a placement).
  auto g8 = shared(make_graph(3, {{0, 1, 0, 1}, {0, 2, 5, 1}, {2, 1, 5, 1}}));
  const auto s8 = evaluate_edge_failures(*g8, v({0}), v({1}), opts, 1);
  EXPECT_EQ(s8.baseline.call_delivered, 0);
  std::vector<std::pair<EdgeId, Amount>> got8;
  for (const auto& e : s8.entries) got8.emplace_back(e.edge, e.delivery_reduction);
  EXPECT_EQ(got8, (std::vector<std::pair<EdgeId, Amount>> {{0, -5}, {1, 0}, {2, 0}}));
  // Invalid inputs propagate (before any worker).
  EXPECT_THROW((void)evaluate_edge_failures(*p, v({0}), v({0}), opts, 1), std::invalid_argument);
  ShortestPathCache cache(*p);
  EXPECT_THROW((void)evaluate_edge_failures(*p, v({0}), v({1}), EvaluationOptions {opts.config, nullptr, &cache}, 2),
               std::invalid_argument);
  EXPECT_EQ(evaluate_edge_failures(*p, v({0}), v({1}), EvaluationOptions {opts.config, nullptr, &cache}, 1).entries.size(), 3U);
}

// The baseline is evaluate_demand() with the requested outputs, in every routing model (T1 with its
// design); without a request it carries none.
TEST(Evaluators, SensitivityBaselineIsCapacityWithTheRequestedOutputs) {
  auto g = diamond(5, 3, 2);
  const Pinned pin = pinned_dags(*g, 0, 3, {dag_of(*g, {0, 2}), dag_of(*g, {1, 3})});
  const OutputOptions all {.edge_flows = true, .residual = true, .saturated_boundary_edges = true, .edge_reservations = true};
  for (const auto& c : {global_config(), shortest_paths_config(), t1(2), FlowPolicyConfig {.routing = pin},
                        FlowPolicyConfig {.routing = centralized(SourceAllocation::Optimized, dag_of(*g, {0, 2}))}}) {
    const auto with = [&](EvaluationOptions q) {
      if (std::holds_alternative<DistributedTE>(c.routing)) q.design_demand = Amount {6};
      return q;
    };
    const auto s = evaluate_edge_failures(*g, v({0}), v({3}), with({c, nullptr, nullptr, all}), 1);
    const auto want = evaluate_demand(*g, v({0}), v({3}), with({c, nullptr, nullptr, all}));
    EXPECT_GT(want.call_delivered, 0);
    ASSERT_TRUE(want.edge_flows.has_value() && want.residual.has_value() && want.saturated_boundary_edges.has_value() &&
                want.edge_reservations.has_value());
    expect_same(s.baseline, want, "sensitivity baseline with every output");
    const auto plain = evaluate_edge_failures(*g, v({0}), v({3}), with({c}), 1);
    expect_same(plain.baseline, evaluate_demand(*g, v({0}), v({3}), with({c})), "sensitivity baseline, no output");
  }
}

// ShortestPaths re-evaluates a deletion exactly when the edge lies in the static shortest-path support
// (every edge on a shortest path to the destinations) reachable from a roster source; every other
// candidate gets loss 0 without a placement. A supplied building cache shows which deletions were
// placed: each placed deletion looks up its own fresh Eligibility uid (one miss and one entry),
// after the baseline's one. Edges 0->2, 1->2, 2->3, 2->4, 3->5, 4->5 of 10 units at cost 1: all six
// form the support toward 5. From source 0 the Optimized controller sends 10 over 0->2->3->5 (the
// lowest EdgeIds first), but 2->4 and 4->5, carrying nothing, are in its reachable support and are
// re-evaluated (loss 0); 1->2 is not and is skipped. From sources 0 and 1 (pooled, 20) the second
// source reaches 1->2 too, and every deletion loses 10.
TEST(Evaluators, ShortestPathsSensitivityReevaluatesExactlyTheReachableSupport) {
  auto g = shared(make_graph(6, {{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 10, 1}, {2, 4, 10, 1}, {3, 5, 10, 1}, {4, 5, 10, 1}}));
  const auto in = id_of_input(*g);
  const auto c = cfg(ShortestPaths {}, RouteAllocation::Optimized);
  const auto losses_of = [&](const SensitivityResult& r) {
    std::map<std::size_t, Amount> out;  // by input position
    for (const auto& e : r.entries) {
      out[static_cast<std::size_t>(std::find(in.begin(), in.end(), e.edge) - in.begin())] = e.delivery_reduction;
    }
    return out;
  };
  ShortestPathCache one(*g);
  const auto r1 = evaluate_edge_failures(*g, v({0}), v({5}), EvaluationOptions {c, nullptr, &one}, 1);
  EXPECT_EQ(r1.baseline.call_delivered, 10);
  EXPECT_EQ(losses_of(r1), (std::map<std::size_t, Amount> {{0, 10}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}}));
  EXPECT_EQ(one.stats().hits, 0);
  EXPECT_EQ(one.stats().misses, 1 + 5);  // the baseline, then every candidate but 1->2
  EXPECT_EQ(one.stats().entries, 1 + 5);
  ShortestPathCache two(*g);
  const auto r2 = evaluate_edge_failures(*g, v({0, 1}), v({5}), EvaluationOptions {c, nullptr, &two}, 1);
  EXPECT_EQ(r2.baseline.call_delivered, 20);
  EXPECT_EQ(losses_of(r2), (std::map<std::size_t, Amount> {{0, 10}, {1, 10}, {2, 10}, {3, 10}, {4, 10}, {5, 10}}));
  EXPECT_EQ(two.stats().misses, 1 + 6);  // every candidate
  EXPECT_EQ(two.stats().entries, 1 + 6);
}

// evaluate_edge_failures's worker and cache rules are evaluate_pairs's, checked in its order before any work:
// workers >= 1, the Eligibility's graph, the cache's graph, a supplied building cache only with
// one worker; a frozen cache is shared by every worker. A model that reads no cache (CentralizedTE)
// rejects one, frozen or not, after those checks.
TEST(Evaluators, SensitivityChecksWorkersAndCachesAsTheBatchDoes) {
  auto g = diamond(5, 3);
  auto other = diamond(5, 3);
  const auto foreign = other->eligibility({}, std::vector<EdgeId> {0});
  ShortestPathCache foreign_cache(*other);
  ShortestPathCache building(*g);
  const auto first_error = [&](const Eligibility* elig, ShortestPathCache* cache, std::int32_t workers) {
    return error_of([&] { (void)evaluate_edge_failures(*g, v({0}), v({3}), EvaluationOptions {shortest_paths_config(), elig, cache}, workers); });
  };
  EXPECT_EQ(first_error(&foreign, &foreign_cache, 0), "workers must be >= 1 (got 0)");
  EXPECT_EQ(first_error(&foreign, &foreign_cache, 2), "foreign Eligibility: built for another graph");
  EXPECT_EQ(first_error(nullptr, &foreign_cache, 2), "foreign ShortestPathCache: bound to another graph");
  EXPECT_EQ(first_error(nullptr, &building, 2),
            "a supplied building ShortestPathCache has one writer: it is rejected unless workers == 1 (got workers = 2; a "
            "frozen cache may be shared by every worker)");
  EXPECT_EQ(stats_of(building.stats()), stats_of(ShortestPathCache(*g).stats()));  // no lookup ran
  building.freeze();
  for (const std::int32_t workers : {1, 2, 4}) {
    EXPECT_EQ(error_of([&] {
                (void)evaluate_edge_failures(*g, v({0}), v({3}), EvaluationOptions {global_config(), nullptr, &building}, workers);
              }),
              "CentralizedTE does not use ShortestPathCache");
  }
  EXPECT_EQ(stats_of(building.stats()), stats_of(ShortestPathCache(*g).stats()));
}

// Only ShortestPaths reads the cache: CentralizedTE, T1 and Pinned reject one (frozen too) before
// any lookup and count nothing; ShortestPaths's baseline raises the miss of an empty frozen cache.
TEST(Evaluators, SensitivityReadsTheCacheOnlyUnderShortestPaths) {
  auto g = diamond(5, 3, 2);
  const Pinned pin = pinned_dags(*g, 0, 3, {dag_of(*g, {0, 2}), dag_of(*g, {1, 3})});
  ShortestPathCache empty(*g);
  empty.freeze();
  for (const auto& [c, model] :
       {std::pair {global_config(), "CentralizedTE"},
        std::pair {t1(2), "DistributedTE"},
        std::pair {FlowPolicyConfig {.routing = pin}, "Pinned"}}) {
    for (const std::int32_t workers : {1, 3}) {
      EXPECT_EQ(error_of([&] {
                  (void)evaluate_edge_failures(*g, v({0}), v({3}), EvaluationOptions {c, nullptr, &empty}, workers);
                }),
                std::string(model) + " does not use ShortestPathCache");
    }
  }
  EXPECT_EQ(stats_of(empty.stats()), stats_of(ShortestPathCache(*g).stats()));
  EXPECT_THROW((void)evaluate_edge_failures(*g, v({0}), v({3}), EvaluationOptions {shortest_paths_config(), nullptr, &empty}, 1),
               ShortestPathCacheMiss);
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
// NGRAPH_CORE_BATCH_THREADS for evaluate_pairs and NGRAPH_CORE_SENSITIVITY_THREADS for
// evaluate_edge_failures, each an integer in [1, 4096] read independently of the other; 1 when unset or
// invalid. " 4" reads as 4 (leading blanks are part of the integer syntax, as for strtol and
// Python's int()); trailing characters, a sign alone and values outside the range
// give 1.
// ---------------------------------------------------------------- the offer of every model

// PAR(a, b): two parallel links 0 -> 1 of capacities a and b, cost 1.
std::shared_ptr<const Graph> par(Amount a, Amount b) { return shared(make_graph(2, {{0, 1, a, 1}, {0, 1, b, 1}})); }
EvaluationOptions offered(FlowPolicyConfig c, std::optional<std::optional<Amount>> offer,
                          std::optional<std::optional<Amount>> design = std::nullopt) {
  EvaluationOptions o {std::move(c)};
  o.offer = offer;
  o.design_demand = design;
  return o;
}
const std::optional<std::optional<Amount>> kNoOffer;
std::optional<std::optional<Amount>> units(Amount a) { return std::optional<std::optional<Amount>> {std::in_place, a}; }

// ECMP Lossless on PAR(4, 10): uncapped 8; an offer below the bound in full, above it up to 8.
TEST(Evaluators, OfferShortestPaths) {
  const auto g = par(4, 10);
  const auto ecmp = shortest_paths_config();
  auto r = evaluate_demand(*g, v({0}), v({1}), offered(ecmp, kNoOffer));
  EXPECT_EQ(std::tuple(r.call_delivered, r.termination), std::tuple(Amount {8}, Termination::Exhausted));
  r = evaluate_demand(*g, v({0}), v({1}), offered(ecmp, units(5)));
  EXPECT_EQ(std::tuple(r.call_delivered, r.termination), std::tuple(Amount {5}, Termination::Complete));
  r = evaluate_demand(*g, v({0}), v({1}), offered(ecmp, units(12)));
  EXPECT_EQ(std::tuple(r.call_delivered, r.termination), std::tuple(Amount {8}, Termination::Exhausted));
  r = evaluate_demand(*g, v({0}), v({1}), offered(ecmp, units(0)));
  EXPECT_EQ(std::tuple(r.call_delivered, r.termination), std::tuple(Amount {0}, Termination::Complete));
}

// T1, design 8 (one tunnel on each link): Lossless 8, or 5 at the offer 5; Lossy 12 drops 2 on
// link 0; the missing offer and the missing design.
TEST(Evaluators, OfferDistributed) {
  const auto g = par(4, 10);
  EXPECT_EQ(evaluate_demand(*g, v({0}), v({1}), offered(t1(2), kNoOffer, units(8))).call_delivered, 8);
  EXPECT_EQ(evaluate_demand(*g, v({0}), v({1}), offered(t1(2), units(5), units(8))).call_delivered, 5);
  const auto r = evaluate_demand(*g, v({0}), v({1}), offered(t1(2, Admission::Lossy), units(12), units(8)));
  EXPECT_EQ(r.call_delivered, 10);
  EXPECT_EQ(r.call_edge_drops, (std::vector<std::pair<EdgeId, Amount>> {{0, 2}}));
  EXPECT_EQ(error_of([&] { (void)evaluate_demand(*g, v({0}), v({1}), offered(t1(2, Admission::Lossy), kNoOffer, units(8))); }),
            "DistributedTE LOSSY evaluation requires a finite offer");
  EXPECT_EQ(error_of([&] { (void)evaluate_demand(*g, v({0}), v({1}), offered(t1(2), kNoOffer)); }),
            "DistributedTE evaluation requires design_demand");
}

// evaluate_edge_failures at a fixed offer: every deletion answers the same offer (EVAL-fixed-offer
// and PAR(4, 10)); without one, the first-hop experiment injects each remaining link's capacity.
TEST(Evaluators, FailuresFixedOffer) {
  const auto ly = FlowPolicyConfig {.routing = ShortestPaths {}, .admission = Admission::Lossy};
  const auto run = [&](const Graph& g, std::optional<std::optional<Amount>> offer, std::int32_t workers) {
    const auto s = evaluate_edge_failures(g, v({0}), v({1}), offered(ly, offer), workers);
    std::vector<Amount> out {s.baseline.call_delivered};
    for (const auto& e : s.entries) out.push_back(e.delivery_reduction);
    return out;
  };
  for (const std::int32_t w : {1, 2}) {
    EXPECT_EQ(run(*par(10, 10), units(4), w), (std::vector<Amount> {4, 0, 0}));
    EXPECT_EQ(run(*par(10, 10), kNoOffer, w), (std::vector<Amount> {20, 10, 10}));
    EXPECT_EQ(run(*par(4, 10), units(12), w), (std::vector<Amount> {10, 0, 6}));
    EXPECT_EQ(run(*par(4, 10), kNoOffer, w), (std::vector<Amount> {14, 4, 10}));
  }
}

// The C++ offer outside its domain (2^62, a negative value) and an explicit +inf: the one text,
// after the configuration and before the design.
TEST(Evaluators, OfferDomain) {
  const auto g = par(4, 10);
  const std::string domain = "offer must be an integer in [0, 2^62)";
  for (const auto& offer : {units(kL62), units(-5), std::optional<std::optional<Amount>> {std::in_place, std::nullopt}}) {
    EXPECT_EQ(error_of([&] { (void)evaluate_demand(*g, v({0}), v({1}), offered(shortest_paths_config(), offer)); }), domain);
    EXPECT_EQ(error_of([&] { (void)evaluate_demand(*g, v({0}), v({1}), offered(t1(2), offer, units(-1))); }), domain);
  }
  EXPECT_EQ(error_of([&] {
              (void)evaluate_demand(*g, v({0}), v({1}),
                                    offered({.routing = centralized(SourceAllocation::Equal, {}, std::nullopt, true)}, units(-5)));
            }),
            "prefer_headroom requires tunnel_count");
  EXPECT_EQ(error_of([&] { (void)evaluate_demand(*g, v({0}), v({1}), offered(shortest_paths_config(), kNoOffer, units(-1))); }),
            "ShortestPaths does not accept design_demand");
  EXPECT_EQ(error_of([&] { (void)evaluate_demand(*g, v({0}), v({1}), offered(shortest_paths_config(), units(kL62 - 1))); }), "");
}

// The failure study of every controller and ShortestPaths at finite offers below, at and above
// the maximum: the baseline is evaluate_demand's, and each allowed edge's reduction is the baseline
// minus evaluate_demand with that edge excluded too (a capped CentralizedTE baseline is no maximum
// flow: its deletions are solved afresh), workers 1 and 2.
TEST(Evaluators, FailuresAtEveryOffer) {
  const auto g = shared(make_graph(4, {{0, 1, 4, 1}, {0, 2, 6, 1}, {1, 3, 3, 2}, {2, 3, 9, 1}}));
  const std::vector<FlowPolicyConfig> configs = {
      {.routing = centralized(SourceAllocation::Optimized)}, {.routing = centralized(SourceAllocation::Equal)},
      {.routing = shortest_paths(RouteAllocation::Optimized)}, {.routing = shortest_paths(RouteAllocation::Equal)},
      {.routing = shortest_paths(RouteAllocation::Weighted)}};
  int checked = 0;
  for (std::size_t c = 0; c < configs.size(); ++c) {
    const Amount top = evaluate_demand(*g, v({0}), v({3}), offered(configs[c], kNoOffer)).call_delivered;
    ASSERT_GT(top, 2);
    for (const Amount offer : {Amount {1}, top / 2, top, top + 5}) {
      const auto o = offered(configs[c], units(offer));
      const Amount base = evaluate_demand(*g, v({0}), v({3}), o).call_delivered;
      for (const std::int32_t workers : {1, 2}) {
        const std::string where = "config " + std::to_string(c) + " offer " + std::to_string(offer) + " workers " +
                                  std::to_string(workers);
        SensitivityResult s;
        ASSERT_EQ(error_of([&] { s = evaluate_edge_failures(*g, v({0}), v({3}), o, workers); }), "") << where;
        EXPECT_EQ(s.baseline.call_delivered, base) << where;
        std::vector<Amount> got(static_cast<std::size_t>(g->num_edges()), 0);
        for (const auto& e : s.entries) got[static_cast<std::size_t>(e.edge)] = e.delivery_reduction;
        for (EdgeId e = 0; e < g->num_edges(); ++e) {
          const EdgeId one[1] = {e};
          const auto mask = g->eligibility({}, one);
          auto deleted = o;
          deleted.eligibility = &mask;
          EXPECT_EQ(got[static_cast<std::size_t>(e)], base - evaluate_demand(*g, v({0}), v({3}), deleted).call_delivered)
              << where << " edge " << e;
        }
        ++checked;
      }
    }
  }
  EXPECT_EQ(checked, 5 * 4 * 2);
}

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
  for (const bool edge_failures : {false, true}) {
    const int self = edge_failures ? 1 : 0;
    for (const auto& [value, want] : cases) {
      const std::string where =
          std::string(names[self]) + "=" + (value == nullptr ? "(unset)" : "'" + std::string(value) + "'");
      set_env(names[self], value);
      set_env(names[1 - self], "7");  // the other variable never applies
      EXPECT_EQ(detail::default_workers(edge_failures), want) << where;
      EXPECT_EQ(detail::default_workers(!edge_failures), 7) << where;
    }
  }
}

// The route choice of ShortestPaths: the single-path order read on the service view
// at preparation. TRI with the service view 4, 8, 5, 50: the record's route is e0, e2, e1, e0 for
// FEWER_HOPS, HIGHER_CAPACITY, HIGHER_AVAILABLE_CAPACITY, HIGHER_HEADROOM.
TEST(RouteChoice, ShortestPathsSinglePath) {
  const auto g = shared(make_graph(2, {{0, 1, 4, 1}, {0, 1, 10, 1}, {0, 1, 100, 1}, {0, 1, 50, 2}}));
  for (const auto& [tie, want] : {std::pair {RouteTieBreak::FewerHops, EdgeId {0}},
                                  std::pair {RouteTieBreak::HigherCapacity, EdgeId {2}},
                                  std::pair {RouteTieBreak::HigherAvailableCapacity, EdgeId {1}},
                                  std::pair {RouteTieBreak::HigherHeadroom, EdgeId {0}}}) {
    FlowLedger fg(g, std::vector<Amount> {4, 8, 5, 50});
    FlowPolicy p(fg, v({0}), v({1}), cfg(ShortestPaths {}, RouteAllocation::Equal, Admission::Lossless, sel(tie, false)));
    EXPECT_EQ(p.place(1).call_delivered, 1) << static_cast<int>(tie);
    ASSERT_EQ(p.records().size(), 1U);
    EXPECT_EQ(p.records()[0].edge_flows, (std::vector<std::pair<EdgeId, Amount>> {{want, 1}})) << static_cast<int>(tie);
  }
}
