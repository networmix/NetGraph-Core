// Prototype 2's 59 fixtures replayed through FlowPolicy (APPROACH 8, 9).
//
// tests/data/policy_fixtures.json was exported from the reference, prototype 2 of the design tree
// (its "source" names it; the design tree is untracked by design; its maintained exporter is
// equivalence/export_policy_fixtures.py there): every ledger, policy and call
// of every fixture in order, with the prototype's values after each operation. The replay checks
// every one of them on the engine: each call's result (placed, remaining, termination, local
// outcomes, flags, drops, the cost curve, the requested outputs), its trace (prototype 2's keys),
// the ledger after it (allocations, residuals, every record field, next FlowId, owners) and the
// policy state (roster, requested, +inf, placed demand, usable bundles, per source the flows'
// ids, routes, T, deliveries, allocations, members and weights). An operation the prototype
// rejects must throw here too (the messages differ). Normalised on both sides: the T_max probe
// counts, the flow core's operation counters and ledger uids. The replay needs the private hooks
// (the trace and the creation cap; NETGRAPH_CORE_PRIVATE_HOOKS) and is skipped without them.
#include <gtest/gtest.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <functional>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/detail/trace.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/graph.hpp"
#include "test_utils.hpp"

#ifndef NETGRAPH_CORE_TEST_DATA
#error "NETGRAPH_CORE_TEST_DATA (the tests/data directory) must be defined by the build"
#endif

#if NETGRAPH_CORE_PRIVATE_HOOKS

using namespace netgraph::core;
using detail::LedgerAccess;

namespace {

// ---------------------------------------------------------------- a small JSON value

struct J {
  enum class K { Null, Bool, Int, Dbl, Str, Arr, Obj };
  K k {K::Null};
  bool b {false};
  __int128 i {0};
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
  static J integer(__int128 v) {
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
  [[nodiscard]] bool is_null() const { return k == K::Null; }
  [[nodiscard]] std::int64_t int64() const {
    if (k != K::Int) throw std::runtime_error("not an integer");
    return static_cast<std::int64_t>(i);
  }
};

std::string int128_text(__int128 v) {
  if (v == 0) return "0";
  const bool neg = v < 0;
  unsigned __int128 u = neg ? static_cast<unsigned __int128>(-(v + 1)) + 1 : static_cast<unsigned __int128>(v);
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
    __int128 v = 0;
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
  for (const auto& x : v) j.push(J::integer(static_cast<__int128>(x)));
  return j;
}

template <class A, class B>
J pairs(const std::vector<std::pair<A, B>>& v) {
  J j = J::arr();
  for (const auto& [x, y] : v) {
    J p = J::arr();
    p.push(J::integer(static_cast<__int128>(x)));
    p.push(J::integer(static_cast<__int128>(y)));
    j.push(std::move(p));
  }
  return j;
}

__int128 wide(const Counter128& c) {
  return (static_cast<__int128>(c.hi) << 64) | static_cast<__int128>(c.lo);
}

const char* placement_name(FlowPlacement p) {
  switch (p) {
    case FlowPlacement::MaxFlow: return "MaxFlow";
    case FlowPlacement::Weighted: return "Weighted";
    case FlowPlacement::EqualBalanced: return "EqualBalanced";
  }
  return "?";
}

const char* outcome_name(Termination t) {
  switch (t) {
    case Termination::Complete: return "complete";
    case Termination::LimitFlows: return "limit_flows";
    case Termination::LimitDiminishing: return "limit_diminishing";
    case Termination::NoRoute: return "no_route";
    case Termination::Exhausted: return "exhausted";
  }
  return "?";
}

J result_json(const PlacementResult& r) {
  J j = J::obj();
  j.set("placed", J::integer(r.placed));
  j.set("remaining", r.remaining.has_value() ? J::integer(*r.remaining) : J::str("inf"));
  j.set("termination", J::str(outcome_name(r.termination)));
  j.set("routes_changed", J::boolean(r.routes_changed));
  j.set("weights_changed", J::boolean(r.weights_changed));
  j.set("drops", pairs(r.drops));
  J lo = J::arr();
  for (const auto& [s, o] : r.local_outcomes) {
    J p = J::arr();
    p.push(J::integer(s));
    p.push(J::str(outcome_name(o)));
    lo.push(std::move(p));
  }
  j.set("local_outcomes", std::move(lo));
  if (r.cost_distribution.has_value()) j.set("cost_distribution", pairs(*r.cost_distribution));
  if (r.edge_flows.has_value()) j.set("edge_flows", ints(*r.edge_flows));
  if (r.residual.has_value()) j.set("residual", ints(*r.residual));
  if (r.bottleneck.has_value()) j.set("bottleneck", ints(*r.bottleneck));
  if (r.min_cut.has_value()) j.set("min_cut", ints(*r.min_cut));
  return j;
}

J ledger_json(const FlowLedger& fg) {
  J j = J::obj();
  j.set("allocated", ints(fg.allocated()));
  j.set("residual", ints(fg.residual()));
  J recs = J::arr();
  for (const auto& r : fg.records()) {
    J x = J::obj();
    x.set("flow_id", J::integer(r.flow_id.second));
    x.set("owner", J::integer(r.owner));
    x.set("sources", ints(r.sources));
    x.set("destinations", ints(r.destinations));
    x.set("placement", J::str(placement_name(r.placement)));
    x.set("admission", J::str(r.admission == Admission::Lossy ? "Lossy" : "Lossless"));
    x.set("bundle", ints(detail::PredDAGBuilder::ascending(*r.bundle)));
    x.set("route_cost", r.route_cost == kUnreachable ? J::null() : J::integer(r.route_cost));
    x.set("offered", r.offered.has_value() ? J::integer(wide(*r.offered)) : J::null());
    x.set("delivered", J::integer(r.delivered));
    x.set("edges", pairs(r.edges));
    J drops = J::arr();
    for (const auto& [e, c] : r.drops) {
      J p = J::arr();
      p.push(J::integer(e));
      p.push(J::integer(wide(c)));
      drops.push(std::move(p));
    }
    x.set("drops", std::move(drops));
    recs.push(std::move(x));
  }
  j.set("records", std::move(recs));
  j.set("next_flow_id", J::integer(LedgerAccess::next_flow_id(fg)));
  j.set("bound", ints(LedgerAccess::bound_owners(fg)));
  j.set("orphans", ints(LedgerAccess::orphan_owners(fg)));
  return j;
}

J policy_json(const FlowPolicy& policy) {
  const auto& pe = detail::PolicyAccess::engine(policy);
  const bool lossy = pe.config.admission == Admission::Lossy;
  J j = J::obj();
  j.set("roster", ints(pe.roster));
  j.set("requested", J::integer(wide(pe.requested)));
  j.set("infinite", J::boolean(pe.infinite));
  j.set("placed_demand", J::integer(policy.placed_demand()));
  j.set("usable", J::integer(policy.usable_bundles()));
  J cohorts = J::arr();
  for (std::size_t i = 0; i < pe.roster.size(); ++i) {
    if (pe.cohorts[i].empty()) continue;
    J flows = J::arr();
    for (const auto& s : pe.cohorts[i]) {
      const FlowRecord& rec = LedgerAccess::live(pe.fg, s.id);
      std::vector<EdgeId> route;
      std::vector<EdgeId> members;
      if (s.plan) {
        route = s.plan->route;
        members = s.plan->members;
      } else {
        const auto e = detail::PredDAGBuilder::ascending(*s.bundle);
        route.assign(e.begin(), e.end());
        members = route;
      }
      J f = J::obj();
      f.set("id", J::integer(s.id));
      f.set("route", ints(route));
      f.set("T", lossy ? J::null() : J::integer(rec.delivered));
      f.set("delivered", J::integer(rec.delivered));
      f.set("alloc", pairs(rec.edges));
      f.set("members", ints(members));
      if (s.plan && s.plan->placement == FlowPlacement::Weighted) {
        J w = J::arr();
        for (std::size_t m = 0; m < s.plan->members.size(); ++m) {
          J p = J::arr();
          p.push(J::integer(s.plan->members[m]));
          p.push(J::integer(s.plan->weights[m]));
          w.push(std::move(p));
        }
        f.set("weights", std::move(w));
      } else {
        f.set("weights", J::null());
      }
      flows.push(std::move(f));
    }
    J entry = J::arr();
    entry.push(J::integer(pe.roster[i]));
    entry.push(std::move(flows));
    cohorts.push(std::move(entry));
  }
  j.set("cohorts", std::move(cohorts));
  return j;
}

// The trace with the normalisation of the exported data: probe counts, counters and ledger uids
// dropped, FlowIds as ids.
J trace_json(const detail::Trace& t) {
  using K = detail::Trace::Kind;
  switch (t.kind) {
    case K::Null: return J::null();
    case K::Bool: return J::boolean(t.b);
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
        if (key == "evals" || key == "ledger_uid" || key == "counters" || key == "evaluations") continue;
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

std::shared_ptr<const PredDAG> dag_of(const Graph& g, std::vector<EdgeId> edges) {
  std::sort(edges.begin(), edges.end());
  std::vector<std::vector<EdgeId>> per(static_cast<std::size_t>(g.num_nodes()));
  for (EdgeId e : edges) per[static_cast<std::size_t>(g.edge_dst()[static_cast<std::size_t>(e)])].push_back(e);
  std::vector<std::int32_t> offsets {0};
  std::vector<NodeId> parents;
  std::vector<EdgeId> via;
  for (const auto& entries : per) {
    for (EdgeId e : entries) {
      parents.push_back(g.edge_src()[static_cast<std::size_t>(e)]);
      via.push_back(e);
    }
    offsets.push_back(static_cast<std::int32_t>(via.size()));
  }
  return std::make_shared<const PredDAG>(PredDAG::from_csr(g, offsets, parents, via));
}

FlowPolicyConfig config_of(const Graph& g, const J& c) {
  const J& r = c.at("routing");
  const std::string kind = r.at("kind").s;
  Routing routing = CostOnly {};
  if (kind == "Global") {
    routing = Global {};
  } else if (kind == "Pinned") {
    Pinned p;
    for (const auto& b : r.at("bundles").a) p.bundles.push_back(dag_of(g, int32s(b)));
    routing = std::move(p);
  } else if (kind == "CapacityAware") {
    CapacityAware ca;
    ca.flow_count = static_cast<std::int32_t>(r.at("flow_count").int64());
    if (!r.at("max_flow_count").is_null()) ca.max_flow_count = static_cast<std::int32_t>(r.at("max_flow_count").int64());
    if (!r.at("max_path_cost").is_null()) {
      const J& mc = r.at("max_path_cost");
      if (mc.a[0].s == "Absolute") {
        ca.max_path_cost = Absolute {mc.a[1].int64()};
      } else {
        ca.max_path_cost = Stretch {mc.a[1].k == J::K::Dbl ? mc.a[1].d : static_cast<double>(mc.a[1].int64())};
      }
    }
    ca.min_flow_bandwidth = r.at("min_flow_bandwidth").int64();
    ca.reoptimize = r.at("reoptimize").b;
    if (!r.at("diminishing").is_null()) {
      const J& d = r.at("diminishing");
      ca.diminishing = DiminishingReturns {static_cast<std::int32_t>(d.a[0].int64()),
                                    d.a[1].k == J::K::Dbl ? d.a[1].d : static_cast<double>(d.a[1].int64())};
    }
    routing = ca;
  }
  const std::string placement = c.at("placement").s;
  const FlowPlacement fp = placement == "MaxFlow"    ? FlowPlacement::MaxFlow
                           : placement == "Weighted" ? FlowPlacement::Weighted
                                                     : FlowPlacement::EqualBalanced;
  const Admission adm = c.at("admission").s == "Lossy" ? Admission::Lossy : Admission::Lossless;
  const J& sel = c.at("selection");
  const EdgeTieBreak tie =
      sel.at("tie_break").s == "Deterministic" ? EdgeTieBreak::Deterministic : EdgeTieBreak::PreferHigherPriority;
  const std::string pr = sel.at("priority").s;
  const EdgePriority prio = pr == "Residual"   ? EdgePriority::Residual
                            : pr == "Capacity" ? EdgePriority::Capacity
                                               : EdgePriority::Headroom;
  return FlowPolicyConfig {std::move(routing), fp, adm,
                           PathSelection(tie, prio, sel.at("multipath").b, sel.at("multi_edge").b)};
}

// Restores the creation cap whatever happens.
struct CapGuard {
  std::int64_t old {detail::max_flows_per_call().load()};
  ~CapGuard() { detail::max_flows_per_call().store(old); }
};

// Replays one fixture; returns the failures (empty when every value matches).
std::vector<std::string> replay(const J& fx) {
  std::vector<std::string> fails;
  const std::string name = fx.at("name").s;
  CapGuard guard;
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
    const auto ii = g->input_index();
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
  std::vector<J> curves;  // per policy: the cumulative Global curve of the epoch (F5 normalisation)
  auto check = [&](const std::string& where, const char* what, const J& expected, const J& actual) {
    std::string why;
    if (!same(expected, actual, what, why)) fails.push_back(name + " " + where + ": " + why);
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
        if (op.at("baseline").is_null()) {
          ledgers.emplace_back(g);
        } else {
          ledgers.emplace_back(g, amounts(op.at("baseline")));
        }
        if (error) fails.push_back(name + " " + where + ": accepted, the prototype rejects it");
        ledger_of_id.push_back(&ledgers.back());
      } catch (const std::exception& e) {
        if (!error) fails.push_back(name + " " + where + ": throws " + e.what());
      }
      continue;
    }
    if (kind == "copy" || kind == "reset") {
      FlowLedger& fg = *ledger_of_id[static_cast<std::size_t>(op.at("ledger").int64())];
      bool threw = false;
      try {
        if (kind == "copy") {
          ledgers.emplace_back(fg);
          ledger_of_id.push_back(&ledgers.back());
        } else if (op.at("baseline").is_null()) {
          fg.reset();
        } else {
          fg.reset(amounts(op.at("baseline")));
        }
      } catch (const std::exception&) {
        threw = true;
      }
      if (threw != error) {
        fails.push_back(name + " " + where + (error ? ": accepted, the prototype rejects it" : ": throws"));
        return fails;
      }
      const FlowLedger& shown = (kind == "copy" && !error) ? ledgers.back() : fg;
      check(where, "ledger", expect->at("ledger"), ledger_json(shown));
      continue;
    }
    if (kind == "policy") {
      FlowLedger& fg = *ledger_of_id[static_cast<std::size_t>(op.at("ledger").int64())];
      const auto& g = *fg.graph();
      std::unique_ptr<FlowPolicy> pol;
      try {
        const Eligibility* elig = nullptr;
        if (!op.at("elig").is_null()) {
          eligs.push_back(g.eligibility(int32s(op.at("elig").at("nodes")), int32s(op.at("elig").at("edges"))));
          elig = &eligs.back();
        }
        auto cfg = config_of(g, op.at("config"));
        pol = std::make_unique<FlowPolicy>(fg, endpoints(op.at("sources")), endpoints(op.at("destinations")),
                                           std::move(cfg), elig);
      } catch (const std::exception& e) {
        if (!error) fails.push_back(name + " " + where + ": throws " + e.what());
        continue;
      }
      if (error) {
        fails.push_back(name + " " + where + ": accepted, the prototype rejects it");
        return fails;
      }
      policies.push_back(std::move(pol));
      policy_ledger.push_back(&fg);
      curves.push_back(J::arr());
      check(where, "policy", expect->at("policy"), policy_json(*policies.back()));
      continue;
    }
    const auto pi = static_cast<std::size_t>(op.at("policy").int64());
    const FlowLedger& fg = *policy_ledger[pi];
    if (kind == "destroy") {
      policies[pi] = nullptr;  // the destruction law: the records stay as orphans
      check(where, "ledger", expect->at("ledger"), ledger_json(fg));
      continue;
    }
    FlowPolicy& pol = *policies[pi];
    if (kind == "remove") {
      pol.remove();
      curves[pi] = J::arr();
      check(where, "ledger", expect->at("ledger"), ledger_json(fg));
      check(where, "policy", expect->at("policy"), policy_json(pol));
      continue;
    }
    // place
    std::optional<Amount> volume;
    if (op.at("volume").k == J::K::Int) volume = op.at("volume").int64();
    if (op.at("volume").k == J::K::Str && op.at("volume").s != "inf") {
      fails.push_back(name + " " + where + ": unexpected volume " + dump(op.at("volume")));
      continue;
    }
    if (op.at("volume").k == J::K::Dbl) {
      // A non-integer request: the prototype rejects it, as the bindings do; C++ has no such value.
      if (!error) fails.push_back(name + " " + where + ": a non-integer request accepted by the prototype");
      continue;
    }
    OutputOptions outputs;
    if (!op.at("outputs").is_null()) {
      const auto& o = op.at("outputs").a;
      outputs = OutputOptions {o[0].b, o[1].b, o[2].b};
    }
    detail::max_flows_per_call().store(op.at("cap").int64());
    detail::Trace trace;
    PlacementResult r;
    bool threw = false;
    try {
      r = detail::place_traced(pol, volume, outputs, trace);
    } catch (const std::exception&) {
      threw = true;
    }
    if (threw != error) {
      fails.push_back(name + " " + where + (error ? ": accepted, the prototype rejects it" : ": throws"));
      return fails;
    }
    if (error) {
      check(where, "ledger", expect->at("ledger"), ledger_json(fg));
      continue;
    }
    J actual = result_json(r);
    const J& expected_result = expect->at("result");
    if (std::holds_alternative<Global>(pol.config().routing)) {
      // Contract over prototype (1c NOTES P3, review-1c F5): a Global zero call reports the
      // policy's cumulative curve (9.1), which prototype 2 omits; it must equal the curve of the
      // epoch's last positive call ([] before any), and only then is it set aside.
      const J* got = actual.get("cost_distribution");
      if (volume.has_value() && *volume == 0 && expected_result.get("cost_distribution") == nullptr && got != nullptr) {
        std::string why;
        if (!same(curves[pi], *got, "cost_distribution", why)) {
          fails.push_back(name + " " + where + ": the zero call's curve is not the cumulative one: " + why);
        }
        std::erase_if(actual.o, [](const auto& kv) { return kv.first == "cost_distribution"; });
      }
      if (const J* c = expected_result.get("cost_distribution"); c != nullptr) curves[pi] = *c;
    }
    check(where, "result", expected_result, actual);
    check(where, "trace", expect->at("trace"), trace_json(trace));
    check(where, "ledger", expect->at("ledger"), ledger_json(fg));
    check(where, "policy", expect->at("policy"), policy_json(pol));
  }
  return fails;
}

} // namespace

TEST(PolicyFixtures, TheDataHoldsEveryPrototypeFixture) {
  const J& d = data();
  EXPECT_EQ(d.at("format").s, "p2-fixtures/1");
  EXPECT_EQ(d.at("fixtures").a.size(), 59U);
}

TEST(PolicyFixtures, EveryFixtureReplaysIdentically) {
  const J& d = data();
  std::size_t ops = 0;
  for (const auto& fx : d.at("fixtures").a) {
    ops += fx.at("ops").a.size();
    const auto fails = replay(fx);
    for (const auto& f : fails) ADD_FAILURE() << f;
  }
  EXPECT_GT(ops, 600U);
}

// The comparison is not vacuous: a corrupted expectation is reported.
TEST(PolicyFixtures, ACorruptedExpectationFails) {
  const J& d = data();
  std::size_t corrupted = 0;
  std::size_t detected = 0;
  for (const auto& fx : d.at("fixtures").a) {
    for (std::size_t n = 0; n < fx.at("ops").a.size() && corrupted < 200; ++n) {
      const J& op = fx.at("ops").a[n];
      if (op.at("op").s != "place" || op.at("error").b) continue;
      for (const char* key : {"result", "trace", "ledger", "policy"}) {
        J copy = fx;
        J* ops = nullptr;
        for (auto& [k, v] : copy.o) {
          if (k == "ops") ops = &v;
        }
        J* exp = nullptr;
        for (auto& [k, v] : ops->a[n].o) {
          if (k == "expect") exp = &v;
        }
        for (auto& [k, v] : exp->o) {
          if (k != key) continue;
          // Corrupt the first integer field of this value.
          std::function<bool(J&)> corrupt = [&](J& x) {
            if (x.k == J::K::Int) {
              x.i += 1;
              return true;
            }
            for (auto& y : x.a) {
              if (corrupt(y)) return true;
            }
            for (auto& [kk, y] : x.o) {
              if (corrupt(y)) return true;
            }
            return false;
          };
          if (!corrupt(v)) continue;
          ++corrupted;
          if (!replay(copy).empty()) ++detected;
        }
      }
      break;  // one placement per fixture
    }
  }
  EXPECT_GT(corrupted, 100U);
  EXPECT_EQ(detected, corrupted);
}

#else

// The replay reads the trace and sets the creation cap: without the private hooks it is skipped.
TEST(PolicyFixtures, TheDataHoldsEveryPrototypeFixture) { GTEST_SKIP() << "needs the private hooks"; }
TEST(PolicyFixtures, EveryFixtureReplaysIdentically) { GTEST_SKIP() << "needs the private hooks"; }
TEST(PolicyFixtures, ACorruptedExpectationFails) { GTEST_SKIP() << "needs the private hooks"; }

#endif
