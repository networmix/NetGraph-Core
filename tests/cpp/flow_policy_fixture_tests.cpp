// The 31 fixtures of the independent reference replayed through FlowPolicy.
//
// tests/data/policy_fixtures.json was exported from the independent reference (its
// frozen copy, the exporter, the regeneration check and the meaning of the data's historical
// "source" string: tests/reference/policy/README.md): every ledger, policy and call of every
// fixture in order, with the reference's values after each operation. The replay checks
// every one of them on the engine: each call's result (this call's values, termination and flags,
// the demand's state, the rows per source, the total cost, the cost curve, the requested outputs),
// its trace (the reference's keys),
// the ledger after it (allocations, residuals, every record field, next FlowId, owners) and the
// policy state (roster, requested, +inf, placed demand, per source the flows'
// ids, routes, T, deliveries and allocations). An operation the reference rejects must throw
// std::invalid_argument here too (the messages differ); any other exception is a failure.
// Normalised on both sides: the flow core's operation counters and ledger uids. The trace needs
// the private hooks (NETGRAPH_CORE_PRIVATE_HOOKS): without them, in the shipped configuration,
// the replay compares everything but the trace, and the corruption test corrupts what that
// replay compares: every section but the trace.
//
// The data keeps the reference's own names; the replay maps them to the library's here and nowhere
// else: the routing kinds "CostOnly" and "Global" build ShortestPaths and CentralizedTE (with
// source allocation OPTIMIZED, Global's default there), a configuration's "placement" is
// ShortestPaths' allocation ("MaxFlow" and "EqualBalanced" are RouteAllocation::Optimized and
// Equal) and its "weights" ShortestPaths' weight_basis; the configuration's "selection" is
// ShortestPaths' route choice: "Deterministic" is FEWER_HOPS
// whatever the "priority", "PreferHigherPriority" with "Residual", "Capacity" or "Headroom" is
// HIGHER_AVAILABLE_CAPACITY, HIGHER_CAPACITY or HIGHER_HEADROOM, "multipath" is multipath and
// "multi_edge" parallel_edges; a result's
// "placed", "remaining", "call_drops", "call_source_loss", "source_loss", "drops",
// "cost_distribution" and "bottleneck" are call_delivered, call_undelivered, call_edge_drops,
// call_source_drops, source_drops, edge_drops, marginal_cost_curve and saturated_boundary_edges; the
// ledger's "allocated" is carried(); a record's "owner", "bundle", "route_cost" and "edges" are
// policy_id, route, min_route_cost and edge_flows; the policy's "placed_demand" is delivered().
#include <gtest/gtest.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/route_builders.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/detail/trace.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/graph.hpp"
#include "test_utils.hpp"

#ifndef NETGRAPH_CORE_TEST_DATA
#error "NETGRAPH_CORE_TEST_DATA (the tests/data directory) must be defined by the build"
#endif

using namespace netgraph::core;
using detail::LedgerAccess;
using test::Int128;
using test::UInt128;

namespace {

// ---------------------------------------------------------------- a small JSON value

struct J {
  enum class K { Null, Bool, Int, Dbl, Str, Arr, Obj };
  K k {K::Null};
  bool b {false};
  Int128 i {0};
  double d {0.0};
  std::string s {};
  std::vector<J> a {};
  std::vector<std::pair<std::string, J>> o {};

  static J null() { return {}; }
  static J boolean(bool v) {
    J j;
    j.k = K::Bool;
    j.b = v;
    return j;
  }
  static J integer(Int128 v) {
    J j;
    j.k = K::Int;
    j.i = v;
    return j;
  }
  static J str(std::string v) {
    J j;
    j.k = K::Str;
    j.s = std::move(v);
    return j;
  }
  static J arr() {
    J j;
    j.k = K::Arr;
    return j;
  }
  static J obj() {
    J j;
    j.k = K::Obj;
    return j;
  }
  // Builders for the replay.
  J& set(const std::string& key, J v) {
    k = K::Obj;
    o.emplace_back(key, std::move(v));
    return o.back().second;
  }
  J& push(J v) {
    k = K::Arr;
    a.push_back(std::move(v));
    return a.back();
  }
  [[nodiscard]] const J* get(const std::string& key) const {
    for (const auto& [kk, v] : o) {
      if (kk == key) return &v;
    }
    return nullptr;
  }
  [[nodiscard]] const J& at(const std::string& key) const {
    const J* p = get(key);
    if (p == nullptr) throw std::runtime_error("missing key " + key);
    return *p;
  }
  // Readers for the replay.
  [[nodiscard]] bool is_null() const { return k == K::Null; }
  [[nodiscard]] std::int64_t int64() const {
    if (k != K::Int) throw std::runtime_error("not an integer");
    return static_cast<std::int64_t>(i);
  }
};

class Parser {
public:
  explicit Parser(const std::string& text) : t_(text) {}
  J parse() {
    J v = value();
    ws();
    if (p_ != t_.size()) fail("trailing data");
    return v;
  }

private:
  [[noreturn]] void fail(const std::string& why) const {
    throw std::runtime_error("JSON: " + why + " at offset " + std::to_string(p_));
  }
  void ws() {
    while (p_ < t_.size() && (t_[p_] == ' ' || t_[p_] == '\n' || t_[p_] == '\r' || t_[p_] == '\t')) ++p_;
  }
  bool lit(const char* word) {
    const std::string w(word);
    if (t_.compare(p_, w.size(), w) == 0) {
      p_ += w.size();
      return true;
    }
    return false;
  }
  std::string string() {
    if (t_[p_] != '"') fail("expected a string");
    ++p_;
    std::string out;
    while (p_ < t_.size() && t_[p_] != '"') {
      char c = t_[p_++];
      if (c == '\\') {
        const char e = t_[p_++];
        switch (e) {
          case 'n': c = '\n'; break;
          case 't': c = '\t'; break;
          case 'r': c = '\r'; break;
          case 'u': {
            const auto code = std::strtol(t_.substr(p_, 4).c_str(), nullptr, 16);
            p_ += 4;
            if (code > 0x7f) fail("non-ASCII escape");
            c = static_cast<char>(code);
            break;
          }
          default: c = e;
        }
      }
      out.push_back(c);
    }
    if (p_ >= t_.size()) fail("unterminated string");
    ++p_;
    return out;
  }
  J value() {
    ws();
    if (p_ >= t_.size()) fail("unexpected end");
    const char c = t_[p_];
    if (c == '{') {
      ++p_;
      J j = J::obj();
      ws();
      if (t_[p_] == '}') {
        ++p_;
        return j;
      }
      while (true) {
        ws();
        std::string key = string();
        ws();
        if (t_[p_++] != ':') fail("expected ':'");
        j.o.emplace_back(std::move(key), value());
        ws();
        if (t_[p_] == ',') {
          ++p_;
          continue;
        }
        if (t_[p_] == '}') {
          ++p_;
          return j;
        }
        fail("expected ',' or '}'");
      }
    }
    if (c == '[') {
      ++p_;
      J j = J::arr();
      ws();
      if (t_[p_] == ']') {
        ++p_;
        return j;
      }
      while (true) {
        j.a.push_back(value());
        ws();
        if (t_[p_] == ',') {
          ++p_;
          continue;
        }
        if (t_[p_] == ']') {
          ++p_;
          return j;
        }
        fail("expected ',' or ']'");
      }
    }
    if (c == '"') return J::str(string());
    if (lit("null")) return J::null();
    if (lit("true")) return J::boolean(true);
    if (lit("false")) return J::boolean(false);
    // A number: an integer of any size below 2^127, or a double.
    const std::size_t start = p_;
    if (t_[p_] == '-') ++p_;
    while (p_ < t_.size() && std::isdigit(static_cast<unsigned char>(t_[p_]))) ++p_;
    bool real = false;
    while (p_ < t_.size() && (t_[p_] == '.' || t_[p_] == 'e' || t_[p_] == 'E' || t_[p_] == '+' || t_[p_] == '-' ||
                              std::isdigit(static_cast<unsigned char>(t_[p_])))) {
      real = true;
      ++p_;
    }
    const std::string num = t_.substr(start, p_ - start);
    if (num.empty() || num == "-") fail("bad value");
    if (real) {
      J j;
      j.k = J::K::Dbl;
      j.d = std::strtod(num.c_str(), nullptr);
      return j;
    }
    Int128 v = 0;
    const bool neg = num[0] == '-';
    for (std::size_t n = neg ? 1 : 0; n < num.size(); ++n) v = v * 10 + (num[n] - '0');
    return J::integer(neg ? -v : v);
  }

  const std::string& t_;
  std::size_t p_ {0};
};

const J& data() {
  static const J d = [] {
    const std::string path = std::string(NETGRAPH_CORE_TEST_DATA) + "/policy_fixtures.json";
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    const std::string text = ss.str();
    return Parser(text).parse();
  }();
  return d;
}

} // namespace

TEST(PolicyFixtures, TheDataHoldsEveryReferenceFixture) {
  const J& d = data();
  EXPECT_EQ(d.at("format").s, "p2-fixtures/1");
  EXPECT_EQ(d.at("fixtures").a.size(), 31U);
}

namespace {

std::string int128_text(Int128 v) {
  if (v == 0) return "0";
  const bool neg = v < 0;
  UInt128 u = neg ? static_cast<UInt128>(-(v + 1)) + 1 : static_cast<UInt128>(v);
  std::string out;
  while (u != 0) {
    out.push_back(static_cast<char>('0' + static_cast<int>(u % 10)));
    u /= 10;
  }
  if (neg) out.push_back('-');
  return {out.rbegin(), out.rend()};
}

std::string dump(const J& j) {
  switch (j.k) {
    case J::K::Null: return "null";
    case J::K::Bool: return j.b ? "true" : "false";
    case J::K::Int: return int128_text(j.i);
    case J::K::Dbl: {
      std::ostringstream os;
      os << j.d;
      return os.str();
    }
    case J::K::Str: return "\"" + j.s + "\"";
    case J::K::Arr: {
      std::string out = "[";
      for (std::size_t n = 0; n < j.a.size(); ++n) out += (n ? "," : "") + dump(j.a[n]);
      return out + "]";
    }
    case J::K::Obj: {
      std::string out = "{";
      for (std::size_t n = 0; n < j.o.size(); ++n) out += (n ? "," : "") + ("\"" + j.o[n].first + "\":") + dump(j.o[n].second);
      return out + "}";
    }
  }
  return "?";
}

// Structural equality; objects by key set, arrays in order. `why` names the first difference.
bool same(const J& a, const J& b, const std::string& path, std::string& why) {
  if (a.k != b.k) {
    why = path + ": " + dump(a) + " (expected) != " + dump(b) + " (engine)";
    return false;
  }
  switch (a.k) {
    case J::K::Null: return true;
    case J::K::Bool:
      if (a.b == b.b) return true;
      break;
    case J::K::Int:
      if (a.i == b.i) return true;
      break;
    case J::K::Dbl:
      if (a.d == b.d) return true;
      break;
    case J::K::Str:
      if (a.s == b.s) return true;
      break;
    case J::K::Arr:
      if (a.a.size() != b.a.size()) {
        why = path + ": length " + std::to_string(a.a.size()) + " (expected) != " + std::to_string(b.a.size()) +
              " (engine): " + dump(a) + " vs " + dump(b);
        return false;
      }
      for (std::size_t n = 0; n < a.a.size(); ++n) {
        if (!same(a.a[n], b.a[n], path + "[" + std::to_string(n) + "]", why)) return false;
      }
      return true;
    case J::K::Obj:
      for (const auto& [key, v] : a.o) {
        const J* w = b.get(key);
        if (w == nullptr) {
          why = path + "." + key + ": missing in the engine's value";
          return false;
        }
        if (!same(v, *w, path + "." + key, why)) return false;
      }
      for (const auto& [key, v] : b.o) {
        if (a.get(key) == nullptr) {
          why = path + "." + key + ": only in the engine's value (" + dump(v) + ")";
          return false;
        }
      }
      return true;
  }
  why = path + ": " + dump(a) + " (expected) != " + dump(b) + " (engine)";
  return false;
}

// ---------------------------------------------------------------- the engine as JSON values

template <class T>
J ints(const T& v) {
  J j = J::arr();
  for (const auto& x : v) j.push(J::integer(static_cast<Int128>(x)));
  return j;
}

template <class A, class B>
J pairs(const std::vector<std::pair<A, B>>& v) {
  J j = J::arr();
  for (const auto& [x, y] : v) {
    J p = J::arr();
    p.push(J::integer(static_cast<Int128>(x)));
    p.push(J::integer(static_cast<Int128>(y)));
    j.push(std::move(p));
  }
  return j;
}

Int128 wide(const Counter128& c) {
  return test::int128_of(static_cast<std::int64_t>(c.hi), c.lo);
}

const char* outcome_name(Termination t) {
  switch (t) {
    case Termination::Complete: return "complete";
    case Termination::NoRoute: return "no_route";
    case Termination::Exhausted: return "exhausted";
    case Termination::NoActiveWeight: return "no_active_weight";
  }
  return "?";
}

J counter_pairs(const std::vector<std::pair<EdgeId, Counter128>>& v) {
  J j = J::arr();
  for (const auto& [e, c] : v) {
    J p = J::arr();
    p.push(J::integer(e));
    p.push(J::integer(wide(c)));
    j.push(std::move(p));
  }
  return j;
}

J result_json(const PlacementResult& r) {
  J j = J::obj();
  j.set("placed", J::integer(r.call_delivered));
  j.set("remaining", r.call_undelivered.has_value() ? J::integer(*r.call_undelivered) : J::str("inf"));
  j.set("delivery_delta", J::integer(r.delivery_delta));
  j.set("call_drops", pairs(r.call_edge_drops));
  j.set("call_source_loss", J::integer(r.call_source_drops));
  j.set("termination", J::str(outcome_name(r.termination)));
  j.set("routes_changed", J::boolean(r.routes_changed));
  j.set("weights_changed", J::boolean(r.weights_changed));
  j.set("requested", r.requested.has_value() ? J::integer(wide(*r.requested)) : J::str("inf"));
  j.set("admitted", J::integer(wide(r.admitted)));
  j.set("delivered", J::integer(r.delivered));
  j.set("source_loss", J::integer(wide(r.source_drops)));
  j.set("drops", counter_pairs(r.edge_drops));
  J rows = J::arr();
  for (const auto& x : r.by_source) {
    J row = J::arr();
    row.push(J::integer(x.source));
    row.push(x.offered.has_value() ? J::integer(wide(*x.offered)) : J::null());
    row.push(J::integer(wide(x.admitted)));
    row.push(J::integer(x.delivered));
    row.push(J::integer(wide(x.source_drops)));
    row.push(J::integer(wide(x.total_edge_drops)));
    rows.push(std::move(row));
  }
  j.set("by_source", std::move(rows));
  if (r.offered.has_value()) j.set("offered", J::integer(wide(*r.offered)));
  if (r.unadmitted.has_value()) j.set("unadmitted", J::integer(wide(*r.unadmitted)));
  if (r.total_cost.has_value()) j.set("total_cost", J::integer(wide(*r.total_cost)));
  if (r.marginal_cost_curve.has_value()) j.set("cost_distribution", pairs(*r.marginal_cost_curve));
  if (r.edge_flows.has_value()) j.set("edge_flows", ints(*r.edge_flows));
  if (r.residual.has_value()) j.set("residual", ints(*r.residual));
  if (r.saturated_boundary_edges.has_value()) j.set("bottleneck", ints(*r.saturated_boundary_edges));
  if (r.min_cut.has_value()) j.set("min_cut", ints(*r.min_cut));
  return j;
}

J ledger_json(const FlowLedger& fg) {
  J j = J::obj();
  j.set("allocated", ints(fg.carried()));
  j.set("residual", ints(fg.residual()));
  J recs = J::arr();
  for (const auto& r : fg.records()) {
    J x = J::obj();
    x.set("flow_id", J::integer(r.flow_id.second));
    x.set("owner", J::integer(r.policy_id));
    x.set("sources", ints(r.sources));
    x.set("destinations", ints(r.destinations));
    x.set("admission", J::str(r.admission == Admission::Lossy ? "Lossy" : "Lossless"));
    x.set("bundle", ints(detail::PredecessorDAGBuilder::ascending(*r.route)));
    x.set("route_cost", r.min_route_cost == kUnreachable ? J::null() : J::integer(r.min_route_cost));
    x.set("offered", r.offered.has_value() ? J::integer(wide(*r.offered)) : J::null());
    x.set("delivered", J::integer(r.delivered));
    x.set("edges", pairs(r.edge_flows));
    J drops = J::arr();
    for (const auto& [e, c] : r.edge_drops) {
      J p = J::arr();
      p.push(J::integer(e));
      p.push(J::integer(wide(c)));
      drops.push(std::move(p));
    }
    x.set("drops", std::move(drops));
    recs.push(std::move(x));
  }
  j.set("records", std::move(recs));
  j.set("next_flow_id", J::integer(LedgerAccess::next_record_counter(fg)));
  j.set("bound", ints(LedgerAccess::bound_policy_ids(fg)));
  j.set("orphans", ints(LedgerAccess::orphan_policy_ids(fg)));
  return j;
}

J policy_json(const FlowPolicy& policy) {
  const auto& pe = detail::PolicyAccess::engine(policy);
  const bool lossy = pe.config.admission == Admission::Lossy;
  J j = J::obj();
  j.set("roster", ints(pe.reachable_sources));
  j.set("requested", J::integer(wide(pe.epoch.requested)));
  j.set("infinite", J::boolean(pe.epoch.infinite));
  j.set("placed_demand", J::integer(policy.delivered()));
  J record_groups = J::arr();
  for (std::size_t i = 0; i < pe.reachable_sources.size(); ++i) {
    if (pe.epoch.record_groups[i].empty()) continue;
    J flows = J::arr();
    for (const auto& s : pe.epoch.record_groups[i]) {
      const FlowRecord& rec = LedgerAccess::live(pe.ledger, s.record_counter);
      const auto e = detail::PredecessorDAGBuilder::ascending(*s.route);
      const std::vector<EdgeId> route(e.begin(), e.end());
      J f = J::obj();
      f.set("id", J::integer(s.record_counter));
      f.set("route", ints(route));
      f.set("T", lossy ? J::null() : J::integer(rec.delivered));
      f.set("delivered", J::integer(rec.delivered));
      f.set("alloc", pairs(rec.edge_flows));
      flows.push(std::move(f));
    }
    J entry = J::arr();
    entry.push(J::integer(pe.reachable_sources[i]));
    entry.push(std::move(flows));
    record_groups.push(std::move(entry));
  }
  j.set("cohorts", std::move(record_groups));
  return j;
}

// The trace with the normalisation of the exported data: probe counts, counters and ledger uids
// dropped, FlowIds as ids.
J trace_json(const detail::Trace& t) {
  using K = detail::Trace::Kind;
  switch (t.kind) {
    case K::Null: return J::null();
    case K::Int: return J::integer(t.i);
    case K::Str: return J::str(t.s);
    case K::List: {
      J j = J::arr();
      for (const auto& x : t.items) j.push(trace_json(x));
      return j;
    }
    case K::Map: {
      J j = J::obj();
      for (const auto& [key, x] : t.fields) {
        if (key == "ledger_uid" || key == "counters") continue;
        if ((key == "flow" || key == "flow_id") && x.kind == K::List && x.items.size() == 2) {
          j.set(key, J::integer(x.items[1].i));
          continue;
        }
        j.set(key, trace_json(x));
      }
      return j;
    }
  }
  return J::null();
}

// ---------------------------------------------------------------- the replay

std::vector<std::int32_t> int32s(const J& j) {
  std::vector<std::int32_t> out;
  for (const auto& x : j.a) out.push_back(static_cast<std::int32_t>(x.int64()));
  return out;
}

std::vector<Amount> amounts(const J& j) {
  std::vector<Amount> out;
  for (const auto& x : j.a) out.push_back(x.int64());
  return out;
}

std::vector<NodeId> endpoints(const J& j) {
  if (j.k == J::K::Int) return {static_cast<NodeId>(j.int64())};
  return int32s(j);
}

// A name in the data that the replay does not know: the data and the replay disagree, which is a
// failure of the replay, never a rejection by the engine.
struct UnknownName : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// The value the data names `name` among `values`; any other name throws UnknownName.
template <typename T>
T named(const std::string& what, const std::string& name, std::initializer_list<std::pair<std::string, T>> values) {
  for (const auto& [n, v] : values) {
    if (n == name) return v;
  }
  throw UnknownName("the data names an unknown " + what + " (" + name + ")");
}

FlowPolicyConfig config_of(const Graph& g, const J& c) {
  const J& r = c.at("routing");
  const std::string kind = r.at("kind").s;
  Routing routing = ShortestPaths {};
  if (kind == "Global") {
    CentralizedTE gl(SourceAllocation::Optimized);  // the frozen data's Global: the optimized source allocation
    if (const J* support = r.get("support")) gl.support = test::dag_of(g, int32s(*support));
    routing = gl;
  } else if (kind != "CostOnly") {
    throw UnknownName("the data names an unknown routing (" + kind + ")");
  }
  FlowPolicyConfig out {.routing = std::move(routing)};
  auto* sp = std::get_if<ShortestPaths>(&out.routing);
  if (!c.at("placement").is_null()) {
    const auto allocation = named<RouteAllocation>("placement", c.at("placement").s,
                                                   {{"MaxFlow", RouteAllocation::Optimized},
                                                    {"Weighted", RouteAllocation::Weighted},
                                                    {"EqualBalanced", RouteAllocation::Equal}});
    if (sp == nullptr) throw UnknownName("the data names a placement for a routing without one (" + kind + ")");
    sp->allocation = allocation;
  }
  out.admission = named<Admission>("admission", c.at("admission").s,
                                   {{"Lossless", Admission::Lossless}, {"Lossy", Admission::Lossy}});
  if (const J& sel = c.at("selection"); !sel.is_null()) {
    const bool prefer = named<bool>("tie break", sel.at("tie_break").s,
                                    {{"Deterministic", false}, {"PreferHigherPriority", true}});
    const RouteTieBreak higher = named<RouteTieBreak>("priority", sel.at("priority").s,
                                                      {{"Residual", RouteTieBreak::HigherAvailableCapacity},
                                                       {"Capacity", RouteTieBreak::HigherCapacity},
                                                       {"Headroom", RouteTieBreak::HigherHeadroom}});
    if (sp == nullptr) throw UnknownName("the data names a selection for a routing without a route choice (" + kind + ")");
    sp->tie_break = prefer ? higher : RouteTieBreak::FewerHops;
    sp->multipath = sel.at("multipath").b;
    sp->parallel_edges = sel.at("multi_edge").b;
  }
  if (const J* w = c.get("weights")) {
    const auto basis = named<WeightBasis>("weight basis", w->s,
                                          {{"Capacity", WeightBasis::Capacity},
                                           {"Explicit", WeightBasis::Explicit},
                                           {"ResidualSnapshot", WeightBasis::ResidualSnapshot}});
    if (sp == nullptr) throw UnknownName("the data names weights for a routing without them (" + kind + ")");
    sp->weight_basis = basis;
  }
  return out;
}

// Replays one fixture; returns the failures (empty when every value matches). `compared`, when
// given, counts the expectations compared. `after_rejected_call`, when given, runs with the
// operation's index and its policy right after each call the reference rejects, before the policy
// is compared with its state before the call (the negative control of that comparison changes the
// policy there). `inside_call`, when given, runs with the operation's index inside the guarded
// part of every ledger, policy and placement operation, before the engine's call (the
// negative control of the rejection rule throws there). The engine rejects an operation with
// std::invalid_argument only; any other exception (an internal invariant's std::logic_error,
// std::bad_alloc, a reader of the data) is a failure, also where the reference rejects.
std::vector<std::string> replay(const J& fx, std::size_t* compared = nullptr,
                                const std::function<void(std::size_t, FlowPolicy&)>& after_rejected_call = {},
                                const std::function<void(std::size_t)>& inside_call = {}) {
  std::vector<std::string> fails;
  const std::string name = fx.at("name").s;
  std::vector<std::shared_ptr<const Graph>> graphs;
  for (const auto& spec : fx.at("graphs").a) {
    std::vector<test::InEdge> edges;
    for (const auto& e : spec.at("edges").a) {
      edges.push_back({static_cast<NodeId>(e.a[0].int64()), static_cast<NodeId>(e.a[1].int64()), e.a[2].int64(),
                       e.a[3].int64()});
    }
    const double q = spec.at("q").k == J::K::Dbl ? spec.at("q").d : static_cast<double>(spec.at("q").int64());
    auto g = std::make_shared<const Graph>(
        test::make_graph(static_cast<std::int32_t>(spec.at("n").int64()), edges, q));
    const auto ii = g->input_positions();
    if (std::vector<std::int32_t>(ii.begin(), ii.end()) != int32s(spec.at("input_index"))) {
      fails.push_back(name + ": the engine ranks the edges differently (input_index)");
      return fails;
    }
    graphs.push_back(std::move(g));
  }
  // Lifetimes: graphs outlive the ledgers and the Eligibility values, which outlive the policies.
  std::deque<FlowLedger> ledgers;
  std::vector<FlowLedger*> ledger_of_id;
  std::deque<Eligibility> eligs;
  std::vector<std::unique_ptr<FlowPolicy>> policies;
  std::vector<FlowLedger*> policy_ledger;
  auto check = [&](const std::string& where, const char* what, const J& expected, const J& actual) {
    if (compared != nullptr) ++*compared;
    std::string why;
    if (!same(expected, actual, what, why)) fails.push_back(name + " " + where + ": " + why);
  };
  auto not_a_rejection = [&](const std::string& where, const std::exception& e) {
    fails.push_back(name + " " + where + ": throws an exception other than a rejection: " + e.what());
    return fails;
  };
  auto planted = [&](std::size_t n) {
    if (inside_call) inside_call(n);
  };
  const auto& ops = fx.at("ops").a;
  for (std::size_t n = 0; n < ops.size(); ++n) {
    const J& op = ops[n];
    const std::string kind = op.at("op").s;
    const std::string where = "op " + std::to_string(n) + " (" + kind + ")";
    const bool error = op.at("error").b;
    const J* expect = op.get("expect");
    if (kind == "ledger") {
      const auto& g = graphs[static_cast<std::size_t>(op.at("graph").int64())];
      try {
        planted(n);
        if (op.at("baseline").is_null()) {
          ledgers.emplace_back(g);
        } else {
          ledgers.emplace_back(g, amounts(op.at("baseline")));
        }
        if (error) fails.push_back(name + " " + where + ": accepted, the reference rejects it");
        ledger_of_id.push_back(&ledgers.back());
      } catch (const std::invalid_argument& e) {
        if (!error) {
          fails.push_back(name + " " + where + ": throws " + e.what());
          return fails;  // the later operations name a ledger that does not exist
        }
      } catch (const std::exception& e) {
        return not_a_rejection(where, e);
      }
      continue;
    }
    if (kind == "policy") {
      FlowLedger& fg = *ledger_of_id[static_cast<std::size_t>(op.at("ledger").int64())];
      const auto& g = *fg.graph();
      std::unique_ptr<FlowPolicy> pol;
      try {
        planted(n);
        const Eligibility* elig = nullptr;
        if (!op.at("elig").is_null()) {
          eligs.push_back(g.eligibility(int32s(op.at("elig").at("nodes")), int32s(op.at("elig").at("edges"))));
          elig = &eligs.back();
        }
        auto cfg = config_of(g, op.at("config"));
        pol = std::make_unique<FlowPolicy>(fg, endpoints(op.at("sources")), endpoints(op.at("destinations")),
                                           std::move(cfg), elig);
      } catch (const UnknownName& e) {
        fails.push_back(name + " " + where + ": " + e.what());
        return fails;
      } catch (const std::invalid_argument& e) {
        if (!error) {
          fails.push_back(name + " " + where + ": throws " + e.what());
          return fails;  // the later operations name a policy that does not exist
        }
        continue;
      } catch (const std::exception& e) {
        return not_a_rejection(where, e);
      }
      if (error) {
        fails.push_back(name + " " + where + ": accepted, the reference rejects it");
        return fails;
      }
      policies.push_back(std::move(pol));
      policy_ledger.push_back(&fg);
      check(where, "policy", expect->at("policy"), policy_json(*policies.back()));
      continue;
    }
    const auto pi = static_cast<std::size_t>(op.at("policy").int64());
    const FlowLedger& fg = *policy_ledger[pi];
    if (kind == "remove" && error) {
      // It cannot fail here, so a reference that rejects one is a difference.
      fails.push_back(name + " " + where + ": accepted, the reference rejects it");
      return fails;
    }
    FlowPolicy& pol = *policies[pi];
    if (kind == "remove") {
      pol.remove();
      check(where, "ledger", expect->at("ledger"), ledger_json(fg));
      check(where, "policy", expect->at("policy"), policy_json(pol));
      continue;
    }
    if (kind != "place") {
      fails.push_back(name + " " + where + ": an unknown operation");
      return fails;
    }
    std::optional<Amount> volume;
    if (op.at("volume").k == J::K::Int) volume = op.at("volume").int64();
    if (op.at("volume").k == J::K::Str && op.at("volume").s != "inf") {
      fails.push_back(name + " " + where + ": unexpected volume " + dump(op.at("volume")));
      continue;
    }
    if (op.at("volume").k == J::K::Dbl) {
      // A non-integer request: the reference rejects it, as the bindings do; C++ has no such value.
      if (!error) fails.push_back(name + " " + where + ": a non-integer request accepted by the reference");
      continue;
    }
    OutputOptions outputs;
    if (!op.at("outputs").is_null()) {
      const auto& o = op.at("outputs").a;
      outputs = OutputOptions {o[0].b, o[1].b, o[2].b, o[3].b};
    }
    const J policy_before = error ? policy_json(pol) : J::null();
    detail::Trace trace;
    PlacementResult r;
    bool threw = false;
    try {
      planted(n);
#if NETGRAPH_CORE_PRIVATE_HOOKS
      r = detail::place_traced(pol, volume, outputs, trace);
#else
      r = pol.place(volume, outputs);
#endif
    } catch (const std::invalid_argument&) {
      threw = true;
    } catch (const std::exception& e) {
      return not_a_rejection(where, e);
    }
    if (threw != error) {
      fails.push_back(name + " " + where + (error ? ": accepted, the reference rejects it" : ": throws"));
      return fails;
    }
    if (error) {
      if (after_rejected_call) after_rejected_call(n, pol);
      // A rejected call leaves the ledger as the reference has it and the policy as it was (the
      // data holds no policy section for it, so this comparison is not an expectation).
      check(where, "ledger", expect->at("ledger"), ledger_json(fg));
      std::string why;
      if (!same(policy_before, policy_json(pol), "policy", why)) {
        fails.push_back(name + " " + where + ": the rejected call changed the policy: " + why);
      }
      continue;
    }
    check(where, "result", expect->at("result"), result_json(r));
    if (detail::kTraced) check(where, "trace", expect->at("trace"), trace_json(trace));
    check(where, "ledger", expect->at("ledger"), ledger_json(fg));
    check(where, "policy", expect->at("policy"), policy_json(pol));
  }
  return fails;
}

} // namespace

TEST(PolicyFixtures, EveryFixtureReplaysIdentically) {
  const J& d = data();
  std::size_t sections = 0;  // the expectations the data holds (without the hooks: no trace)
  std::size_t compared = 0;  // the expectations the replay compared
  std::size_t replayed = 0;
  for (const auto& fx : d.at("fixtures").a) {
    ++replayed;
    for (const auto& op : fx.at("ops").a) {
      if (const J* expect = op.get("expect"); expect != nullptr) {
        sections += expect->o.size() - (!detail::kTraced && expect->get("trace") != nullptr ? 1 : 0);
      }
    }
    const auto fails = replay(fx, &compared);
    for (const auto& f : fails) ADD_FAILURE() << f;
  }
  EXPECT_EQ(replayed, 31U);
  EXPECT_EQ(sections, detail::kTraced ? 693U : 548U);  // the data's count, not vacuous
  EXPECT_EQ(compared, sections);  // every expectation was compared, none skipped
}

namespace {

// One leaf changed: a null becomes 0, a boolean flips, a number grows by one, a string gains a
// character. Returns its kind: 0 integer, 1 boolean, 2 string, 3 other.
int corrupt(J& leaf) {
  switch (leaf.k) {
    case J::K::Null: leaf = J::integer(0); return 3;
    case J::K::Bool: leaf.b = !leaf.b; return 1;
    case J::K::Int: leaf.i += 1; return 0;
    case J::K::Dbl: leaf.d += 1.0; return 3;
    case J::K::Str: leaf.s += "#"; return 2;
    default: return 3;
  }
}

// The member `key` of an object, for writing.
J& member(J& x, const std::string& key) {
  for (auto& [k, v] : x.o) {
    if (k == key) return v;
  }
  throw std::out_of_range("no member " + key);
}

// The first leaf of x in document order (an empty container has none); nullptr when x has none.
J* first_leaf(J& x) {
  if (x.k == J::K::Arr) {
    for (auto& y : x.a) {
      if (J* l = first_leaf(y); l != nullptr) return l;
    }
    return nullptr;
  }
  if (x.k == J::K::Obj) {
    for (auto& [k, y] : x.o) {
      if (J* l = first_leaf(y); l != nullptr) return l;
    }
    return nullptr;
  }
  return &x;
}

// Whether the replay of fx reports a failure that starts with `prefix`.
bool reported(const J& fx, const std::string& prefix) {
  return test::has_failure(replay(fx), prefix);
}

// The replay's other failure clauses, each reported at the operation it concerns: the first leaf
// of every section of every other operation's expectations (failed placements, policies,
// removals); the error flag of the first operation of every kind and
// outcome in each fixture; the input order of the first graph; an unknown operation; a ledger the
// data accepts and the engine rejects; an unknown and a non-integer volume; the curve a CentralizedTE
// zero call reports, against the remembered one; a policy changed by a call the data rejects;
// an exception other than a rejection; and a configuration name the replay does not know.
void expect_every_other_clause_reported() {
  const J& d = data();
  std::set<std::string> sections;  // "kind error key" combinations corrupted
  std::set<std::string> flags;     // "kind error" combinations flipped
  std::size_t exercised = 0;       // fixtures whose clauses were exercised
  for (const auto& fx : d.at("fixtures").a) {
    ++exercised;
    const std::string name = fx.at("name").s;
    ASSERT_TRUE(replay(fx).empty()) << name;
    std::set<std::string> flipped;
    for (std::size_t n = 0; n < fx.at("ops").a.size(); ++n) {
      const J& op = fx.at("ops").a[n];
      const std::string kind = op.at("op").s;
      const bool error = op.at("error").b;
      const std::string combination = kind + (error ? " failed" : " ok");
      const std::string at = name + " op " + std::to_string(n) + " (" + kind + ")";
      if (const J* expect = op.get("expect"); expect != nullptr && (kind != "place" || error)) {
        for (const auto& [key, section] : expect->o) {
          J copy = fx;
          J* leaf = first_leaf(member(member(member(copy, "ops").a[n], "expect"), key));
          ASSERT_NE(leaf, nullptr) << at << " " << key;
          corrupt(*leaf);
          EXPECT_TRUE(reported(copy, at + ": " + key)) << "undetected: " << at << " " << key;
          sections.insert(combination + " " + key);
        }
      }
      if (flipped.insert(combination).second) {
        J copy = fx;
        auto& flag = member(member(copy, "ops").a[n], "error");
        flag.b = !flag.b;
        EXPECT_TRUE(reported(copy, at + ": ")) << "undetected error flag: " << at;
        flags.insert(combination);
      }
    }
    J reordered = fx;
    J* first = first_leaf(member(member(reordered, "graphs").a[0], "input_index"));
    ASSERT_NE(first, nullptr) << name;
    corrupt(*first);
    EXPECT_TRUE(reported(reordered, name + ": the engine ranks the edges differently")) << name;
  }
  EXPECT_EQ(exercised, 31U);
  // Every section of every other operation kind and outcome in the data, and every flag.
  EXPECT_EQ(sections, (std::set<std::string> {"place failed ledger", "policy ok policy", "remove ok ledger",
                                              "remove ok policy"}));
  EXPECT_EQ(flags, (std::set<std::string> {"ledger ok", "place failed", "place ok", "policy failed", "policy ok",
                                           "remove ok"}));
  // Op 3 of global_epochs renamed: an unknown operation.
  const auto& fixtures = d.at("fixtures").a;
  const auto& epochs = *std::find_if(fixtures.begin(), fixtures.end(),
                                     [](const J& fx) { return fx.at("name").s == "global_epochs"; });
  ASSERT_EQ(epochs.at("ops").a[3].at("op").s, "place");
  J renamed = epochs;
  member(member(renamed, "ops").a[3], "op") = J::str("move");
  EXPECT_TRUE(reported(renamed, "global_epochs op 3 (move): an unknown operation"));
  // A ledger the data accepts and the engine rejects: op 0's baseline, one entry for two edges.
  ASSERT_EQ(epochs.at("ops").a[0].at("op").s, "ledger");
  J rejected = epochs;
  auto& baseline = member(member(rejected, "ops").a[0], "baseline");
  baseline = J::arr();
  baseline.a.push_back(J::integer(-1));
  EXPECT_TRUE(reported(rejected, "global_epochs op 0 (ledger): throws "));
  // The volume of op 3 of global_epochs: neither "inf" nor an integer, then a non-integer number
  // the data says is accepted.
  J unknown = epochs;
  member(member(unknown, "ops").a[3], "volume") = J::str("many");
  EXPECT_TRUE(reported(unknown, "global_epochs op 3 (place): unexpected volume"));
  J fractional = epochs;
  auto& volume = member(member(fractional, "ops").a[3], "volume");
  volume.k = J::K::Dbl;
  volume.d = 0.5;
  EXPECT_TRUE(reported(fractional, "global_epochs op 3 (place): a non-integer request accepted"));
  // Op 2 of global_epochs is a CentralizedTE zero call on a fresh policy (no curve: nothing prepared),
  // op 3 a positive call with the curve [[1, 1]]. With op 3 moved before op 2, the zero call
  // reports the remembered curve [[1, 1]] where the data, recorded before any placement, has
  // none: a result difference.
  ASSERT_EQ(dump(epochs.at("ops").a[2].at("volume")), "0");
  ASSERT_EQ(epochs.at("ops").a[2].at("expect").at("result").get("cost_distribution"), nullptr);
  J swapped = epochs;
  auto& ops = member(swapped, "ops").a;
  std::swap(ops[2], ops[3]);
  EXPECT_TRUE(reported(swapped, "global_epochs op 3 (place): result"));
  // Every call the data rejects, with the policy changed after it (its +inf flag flipped; after
  // that call only): the comparison with the policy before the call reports it.
  std::size_t rejected_calls = 0;
  for (const auto& fx : fixtures) {
    for (std::size_t n = 0; n < fx.at("ops").a.size(); ++n) {
      const J& op = fx.at("ops").a[n];
      if (op.at("op").s != "place" || !op.at("error").b) continue;
      ++rejected_calls;
      const auto fails = replay(fx, nullptr, [n](std::size_t at, FlowPolicy& p) {
        if (at != n) return;
        auto& engine = detail::PolicyAccess::engine(p);
        engine.epoch.infinite = !engine.epoch.infinite;
      });
      const std::string at = fx.at("name").s + " op " + std::to_string(n) + " (place)";
      EXPECT_TRUE(test::has_failure(fails, at + ": the rejected call changed the policy")) << at;
    }
  }
  EXPECT_EQ(rejected_calls, 3U);
  // An exception other than std::invalid_argument (an internal invariant's std::logic_error)
  // thrown inside the first operation of each kind, the first one the data rejects where there
  // is one: a failure at that operation, never taken for the expected rejection.
  std::set<std::string> planted;
  for (const std::string kind : {"ledger", "policy", "place"}) {
    const J* chosen = nullptr;
    std::size_t at = 0;
    for (const auto& fx : fixtures) {
      for (std::size_t n = 0; n < fx.at("ops").a.size(); ++n) {
        const J& op = fx.at("ops").a[n];
        if (op.at("op").s != kind) continue;
        if (chosen == nullptr || (op.at("error").b && !chosen->at("ops").a[at].at("error").b)) {
          chosen = &fx;
          at = n;
        }
      }
    }
    ASSERT_NE(chosen, nullptr) << kind;
    const bool error = chosen->at("ops").a[at].at("error").b;
    planted.insert(kind + (error ? " failed" : " ok"));
    const auto fails = replay(*chosen, nullptr, {}, [at](std::size_t n) {
      if (n == at) throw std::logic_error("planted");
    });
    const std::string where = chosen->at("name").s + " op " + std::to_string(at) + " (" + kind + ")";
    EXPECT_TRUE(test::has_failure(fails, where + ": throws an exception other than a rejection: planted")) << where;
  }
  EXPECT_EQ(planted, (std::set<std::string> {"ledger ok", "place failed", "policy failed"}));
  // A name the replay does not know, in the configuration of eb_unequal_1_10's ShortestPaths policy
  // (op 1): a failure of the replay, never taken for the engine rejecting the policy.
  const auto& unequal = *std::find_if(fixtures.begin(), fixtures.end(),
                                      [](const J& fx) { return fx.at("name").s == "eb_unequal_1_10"; });
  ASSERT_EQ(unequal.at("ops").a[1].at("op").s, "policy");
  for (const std::string path :
       {"routing kind", "placement", "admission", "selection tie_break", "selection priority"}) {
    J copy = unequal;
    J* config = &member(member(copy, "ops").a[1], "config");
    std::string key = path;
    if (const auto space = path.find(' '); space != std::string::npos) {
      config = &member(*config, path.substr(0, space));
      key = path.substr(space + 1);
    }
    member(*config, key) = J::str("Unknown");
    EXPECT_TRUE(reported(copy, "eb_unequal_1_10 op 1 (policy): the data names an unknown ")) << path;
  }
}

// The comparison's structural clauses, which the leaf corruptions above never reach: equal
// values pass, and an array of another length, a key the engine lacks, a key only the engine has
// and a value of another kind each fail with their own message.
void expect_the_comparison_reports_shape_differences() {
  auto diff = [](const std::string& expected, const std::string& engine) {
    std::string why;
    return same(Parser(expected).parse(), Parser(engine).parse(), "x", why) ? std::string("same") : why;
  };
  EXPECT_EQ(diff(R"({"a": [1, 2], "b": null})", R"({"b": null, "a": [1, 2]})"), "same");
  EXPECT_EQ(diff("[1, 2]", "[1]"), "x: length 2 (expected) != 1 (engine): [1,2] vs [1]");
  EXPECT_EQ(diff("[1]", "[1, 2]"), "x: length 1 (expected) != 2 (engine): [1] vs [1,2]");
  EXPECT_EQ(diff(R"({"a": 1, "b": 2})", R"({"a": 1})"), "x.b: missing in the engine's value");
  EXPECT_EQ(diff(R"({"a": 1})", R"({"a": 1, "c": 3})"), "x.c: only in the engine's value (3)");
  EXPECT_EQ(diff(R"({"a": [1]})", R"({"a": 1})"), "x.a: [1] (expected) != 1 (engine)");
  EXPECT_EQ(diff("[[1, 2]]", "[[1, 3]]"), "x[0][1]: 2 (expected) != 3 (engine)");
}

} // namespace

// The comparison is not vacuous: a corrupted expectation is reported at the operation and by the
// section it changes. In the first successful placement of every fixture, each section (result,
// trace, ledger, policy) gets several single corruptions: its first and last integer, boolean
// and string leaves, and one leaf drawn at random. Then the replay's other clauses and the
// comparison's structural clauses (above). Without the hooks the replay compares no trace, so it
// is not corrupted; the counts below are those of each configuration (the random draws come from std::mt19937_64, whose sequence the
// standard fixes, so they are the same with every standard library).
TEST(PolicyFixtures, ACorruptedExpectationFails) {
  const J& d = data();
  std::mt19937_64 rng(837);
  std::size_t corrupted = 0;
  std::size_t detected = 0;
  std::size_t fixtures = 0;  // fixtures whose first successful placement was corrupted
  std::size_t kinds[4] = {0, 0, 0, 0};  // integers, booleans, strings, others
  for (const auto& fx : d.at("fixtures").a) {
    for (std::size_t n = 0; n < fx.at("ops").a.size(); ++n) {
      const J& op = fx.at("ops").a[n];
      if (op.at("op").s != "place" || op.at("error").b) continue;
      ++fixtures;
      for (const std::string key : {"result", "trace", "ledger", "policy"}) {
        if (!detail::kTraced && key == "trace") continue;
        // The leaves of the section in document order, as paths of child indices.
        std::vector<std::vector<std::size_t>> leaves;
        std::vector<std::size_t> path;
        std::function<void(const J&)> walk = [&](const J& x) {
          if (x.k == J::K::Arr || x.k == J::K::Obj) {
            const std::size_t count = x.k == J::K::Arr ? x.a.size() : x.o.size();
            for (std::size_t i = 0; i < count; ++i) {
              path.push_back(i);
              walk(x.k == J::K::Arr ? x.a[i] : x.o[i].second);
              path.pop_back();
            }
            return;
          }
          leaves.push_back(path);
        };
        const J& section = op.at("expect").at(key);
        walk(section);
        auto leaf = [](auto& root, const std::vector<std::size_t>& at) -> auto& {
          auto* x = &root;
          for (std::size_t i : at) x = x->k == J::K::Arr ? &x->a[i] : &x->o[i].second;
          return *x;
        };
        std::set<std::vector<std::size_t>> picks;
        for (const J::K kind : {J::K::Int, J::K::Bool, J::K::Str}) {
          std::vector<std::vector<std::size_t>> of_kind;
          for (const auto& at : leaves) {
            if (leaf(section, at).k == kind) of_kind.push_back(at);
          }
          if (!of_kind.empty()) {
            picks.insert(of_kind.front());
            picks.insert(of_kind.back());
          }
        }
        if (!leaves.empty()) picks.insert(leaves[rng() % leaves.size()]);
        for (const auto& at : picks) {
          J copy = fx;
          J* ops = nullptr;
          for (auto& [k, v] : copy.o) {
            if (k == "ops") ops = &v;
          }
          J* exp = nullptr;
          for (auto& [k, v] : ops->a[n].o) {
            if (k == "expect") exp = &v;
          }
          J* target = nullptr;
          for (auto& [k, v] : exp->o) {
            if (k == key) target = &leaf(v, at);
          }
          ASSERT_NE(target->k, J::K::Arr) << "a container is not a leaf";
          ASSERT_NE(target->k, J::K::Obj) << "a container is not a leaf";
          ++kinds[corrupt(*target)];
          ++corrupted;
          // Reported at this op, by the comparison of this section.
          const std::string where = fx.at("name").s + " op " + std::to_string(n) + " (place): " + key;
          if (reported(copy, where)) {
            ++detected;
          } else {
            ADD_FAILURE() << "undetected: " << where;
          }
        }
      }
      break;  // one placement per fixture
    }
  }
  RecordProperty("corruptions", static_cast<int>(corrupted));
  EXPECT_EQ(fixtures, 31U);
  EXPECT_EQ(corrupted, detail::kTraced ? 518U : 429U);
  EXPECT_GT(kinds[0], 200U);
  EXPECT_GT(kinds[1], 50U);
  EXPECT_GT(kinds[2], 50U);
  EXPECT_EQ(detected, corrupted);
  expect_every_other_clause_reported();
  expect_the_comparison_reports_shape_differences();
}
