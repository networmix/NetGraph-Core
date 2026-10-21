/* The transaction trace of a placement call (internal diagnostics for the tests, which compare
   it key by key with the same trace of their independent model of the policy): a small
   JSON-like value of the call's intermediate values: the roster, the routes prepared at the
   epoch's first positive call, ShortestPaths's offers or bound, and for the controllers the
   bisection probes, the solve's view and capacities, the flow core's values and the shares.
   Filled only in a build with the private hooks, when a call asks for it (detail::place_traced);
   other calls build only empty values. */
#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace netgraph::core::detail {

struct Trace {
  enum class Kind : std::uint8_t { Null, Int, Str, List, Map };
  Kind kind {Kind::Null};
  std::int64_t i {0};
  std::string s {};
  std::vector<Trace> items {};                           // List
  std::vector<std::pair<std::string, Trace>> fields {};  // Map, in insertion order

  static Trace null() { return {}; }
  static Trace integer(std::int64_t v) {
    Trace t;
    t.kind = Kind::Int;
    t.i = v;
    return t;
  }
  static Trace str(std::string v) {
    Trace t;
    t.kind = Kind::Str;
    t.s = std::move(v);
    return t;
  }
  static Trace list() {
    Trace t;
    t.kind = Kind::List;
    return t;
  }
  static Trace map() {
    Trace t;
    t.kind = Kind::Map;
    return t;
  }
  template <class T>
  static Trace ints(std::span<const T> v) {
    Trace t = list();
    t.items.reserve(v.size());
    for (const auto& x : v) t.items.push_back(integer(static_cast<std::int64_t>(x)));
    return t;
  }
  template <class T>
  static Trace ints(const std::vector<T>& v) {
    return ints(std::span<const T>(v));
  }
  // [[a, b], ...] from pairs.
  template <class A, class B>
  static Trace pairs(const std::vector<std::pair<A, B>>& v) {
    Trace t = list();
    for (const auto& [a, x] : v) {
      Trace p = list();
      p.items.push_back(integer(static_cast<std::int64_t>(a)));
      p.items.push_back(integer(static_cast<std::int64_t>(x)));
      t.items.push_back(std::move(p));
    }
    return t;
  }

  // Map: set (replacing an existing key).
  Trace& set(std::string_view key, Trace v) {
    kind = Kind::Map;
    for (auto& [k, x] : fields) {
      if (k == key) {
        x = std::move(v);
        return x;
      }
    }
    fields.emplace_back(std::string(key), std::move(v));
    return fields.back().second;
  }
  // List: append.
  Trace& push(Trace v) {
    kind = Kind::List;
    items.push_back(std::move(v));
    return items.back();
  }
};

} // namespace netgraph::core::detail
