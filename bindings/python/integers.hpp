/*
  Integer arguments from Python: Python ints, NumPy integer scalars and integer arrays, read
  without narrowing (WideInts) and validated in the C++ constructors' order before any
  conversion to the C++ types, so the reported error does not depend on how the input is
  represented; bools and floats are rejected, never truncated. Also the shape step of sequence
  arguments, the per-edge amount arrays and the node sets.
*/
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/types.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

inline bool is_bool(const py::handle& obj) {
  return PyBool_Check(obj.ptr()) || py::isinstance(obj, py::module_::import("numpy").attr("bool_"));
}

// The Python int of an integer-like object (a Python int or a NumPy integer scalar); bools and
// everything else are rejected, and an exception of the object's own __index__ propagates as
// raised. The result is an exact int, so PyLong_AsLongLongAndOverflow cannot fail on it.
inline py::object python_int(const py::handle& item, const std::string& what) {
  if (is_bool(item) || !PyIndex_Check(item.ptr())) {
    throw py::value_error(what + " must be an integer (got " +
                          std::string(py::str(py::type::of(item).attr("__name__"))) + ")");
  }
  py::object idx = py::reinterpret_steal<py::object>(PyNumber_Index(item.ptr()));
  if (!idx) throw py::error_already_set();
  return idx;
}

// One integer of [lo, hi] (scalar arguments).
inline std::int64_t integer_value(const py::handle& item, std::int64_t lo, std::int64_t hi, const std::string& what) {
  py::object idx = python_int(item, what);
  int overflow = 0;
  const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
  if (overflow != 0 || v < lo || v > hi) {
    throw py::value_error(what + " = " + std::string(py::str(idx)) + " is out of range [" +
                          std::to_string(lo) + ", " + std::to_string(hi) + "]");
  }
  return static_cast<std::int64_t>(v);
}

// A 1-D integer sequence read without narrowing: every entry's int64 value where it fits, and
// the sign of the overflow where it does not (+1 above INT64_MAX, -1 below INT64_MIN). Only type
// defects raise here (a non-integer, a bool, a float array, a non-sequence); every value bound
// is checked by the caller in the C++ constructor's order, so the reported error does not depend
// on the representation of the input.
struct WideInts {
  py::object source;
  std::vector<std::int64_t> value;
  std::vector<std::int8_t> overflow;  // empty when no entry overflows
  [[nodiscard]] std::size_t size() const noexcept { return value.size(); }
  [[nodiscard]] int sign(std::size_t i) const noexcept { return overflow.empty() ? 0 : overflow[i]; }
  // The entry as Python prints it (exact also beyond int64).
  [[nodiscard]] std::string text(std::size_t i) const {
    if (sign(i) == 0) return std::to_string(value[i]);
    return std::string(py::str(source[py::int_(i)]));
  }
  [[nodiscard]] bool in(std::size_t i, std::int64_t lo, std::int64_t hi) const noexcept {
    return sign(i) == 0 && value[i] >= lo && value[i] <= hi;
  }
  // The entry as an exact Python int (an overflowing entry compared with another). The entry is
  // converted again here, so its __index__ may now fail: that exception propagates as raised.
  [[nodiscard]] py::object exact(std::size_t i) const {
    py::object item = source[py::int_(i)];
    py::object idx = py::reinterpret_steal<py::object>(PyNumber_Index(item.ptr()));
    if (!idx) throw py::error_already_set();
    return idx;
  }
  template <typename T>
  [[nodiscard]] std::vector<T> narrow() const {  // after every bound was checked
    return std::vector<T>(value.begin(), value.end());
  }
};

inline bool not_a_sequence(const py::handle& obj) {
  return py::isinstance<py::str>(obj) || py::isinstance<py::bytes>(obj) || !py::isinstance<py::sequence>(obj);
}

inline WideInts wide_integers(const py::handle& obj, const char* name) {
  WideInts out;
  out.source = py::reinterpret_borrow<py::object>(obj);
  if (py::isinstance<py::array>(obj)) {
    auto arr = py::reinterpret_borrow<py::array>(obj);
    if (arr.ndim() != 1) {
      throw py::value_error(std::string(name) + " must be a 1-D array (got " + std::to_string(arr.ndim()) +
                            " dimensions)");
    }
    const char kind = arr.dtype().kind();
    if (kind != 'i' && kind != 'u') {
      throw py::value_error(std::string(name) + " entries must be integers (got dtype " +
                            std::string(py::str(arr.dtype())) + ")");
    }
    out.value.reserve(static_cast<std::size_t>(arr.size()));
    if (kind == 'u' && arr.itemsize() == 8) {
      auto u = py::array_t<std::uint64_t, py::array::forcecast>::ensure(arr);
      auto r = u.unchecked<1>();
      constexpr auto kMax = static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
      for (py::ssize_t i = 0; i < r.shape(0); ++i) {
        if (r(i) > kMax) {
          if (out.overflow.empty()) out.overflow.assign(static_cast<std::size_t>(r.shape(0)), 0);
          out.overflow[static_cast<std::size_t>(i)] = 1;
          out.value.push_back(std::numeric_limits<std::int64_t>::max());
        } else {
          out.value.push_back(static_cast<std::int64_t>(r(i)));
        }
      }
      return out;
    }
    auto a = py::array_t<std::int64_t, py::array::forcecast>::ensure(arr);
    auto r = a.unchecked<1>();
    for (py::ssize_t i = 0; i < r.shape(0); ++i) out.value.push_back(r(i));
    return out;
  }
  if (not_a_sequence(obj)) throw py::value_error(std::string(name) + " must be a sequence of integers");
  auto seq = py::reinterpret_borrow<py::sequence>(obj);
  const std::size_t n = seq.size();
  out.value.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    py::object idx = python_int(seq[i], std::string(name) + "[" + std::to_string(i) + "]");
    int overflow = 0;
    const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
    if (overflow != 0) {
      if (out.overflow.empty()) out.overflow.assign(n, 0);
      out.overflow[i] = static_cast<std::int8_t>(overflow);
      out.value.push_back(overflow > 0 ? std::numeric_limits<std::int64_t>::max()
                                       : std::numeric_limits<std::int64_t>::min());
    } else {
      out.value.push_back(static_cast<std::int64_t>(v));
    }
  }
  return out;
}

// An entry that is itself a sequence or an array of rank >= 1 (a shape defect of its argument).
inline bool is_nested(const py::handle& item) {
  if (py::isinstance<py::array>(item)) return py::reinterpret_borrow<py::array>(item).ndim() >= 1;
  return !py::isinstance<py::str>(item) && !py::isinstance<py::bytes>(item) && PySequence_Check(item.ptr()) != 0;
}

// The shape step: the length of a rank-1 integer sequence argument. A NumPy array of another
// rank (0-D included), a nested sequence or a non-sequence is a shape defect. Entry types and
// values are not read here; their steps come later.
inline std::size_t sequence_length(const py::handle& obj, const char* name) {
  if (py::isinstance<py::array>(obj)) {
    auto arr = py::reinterpret_borrow<py::array>(obj);
    if (arr.ndim() != 1) {
      throw py::value_error(std::string(name) + " must be a 1-D array (got " + std::to_string(arr.ndim()) +
                            " dimensions)");
    }
    return static_cast<std::size_t>(arr.size());
  }
  if (not_a_sequence(obj)) {
    throw py::value_error(std::string(name) + " must be a 1-D sequence of integers (got " +
                          std::string(py::str(py::type::of(obj).attr("__name__"))) + ")");
  }
  auto seq = py::reinterpret_borrow<py::sequence>(obj);
  const std::size_t n = seq.size();
  for (std::size_t i = 0; i < n; ++i) {
    if (is_nested(seq[i])) {
      throw py::value_error(std::string(name) + " must be a 1-D sequence of integers (entry " + std::to_string(i) +
                            " is a sequence)");
    }
  }
  return n;
}

// A per-edge integer array in [0, capacity_units[e]] (a ledger baseline, a search residual,
// the capacities `view` of _place_global_units), validated on the unnarrowed values: the length
// (shape step), then every entry in EdgeId order. `label` names the entries in a type error
// ("residual (units)", "view (capacity units)"). The checked values come back unnarrowed: a
// caller narrows them later (_place_global_units, inside its call).
inline WideInts edge_amounts(const py::object& obj, const Graph& g, const char* name, const char* label) {
  const auto num_edges = static_cast<std::size_t>(g.num_edges());
  const std::size_t len = sequence_length(obj, name);
  if (len != num_edges) {
    throw py::value_error(std::string(name) + " has length " + std::to_string(len) + ", expected " +
                          std::to_string(num_edges) + " (one entry per edge)");
  }
  auto v = wide_integers(obj, label);
  const auto capacity = g.capacity_units();
  for (std::size_t e = 0; e < num_edges; ++e) {
    if (!v.in(e, 0, capacity[e])) {
      throw py::value_error(std::string(name) + "[" + std::to_string(e) + "] = " + v.text(e) + " is outside [0, capacity_units[" +
                            std::to_string(e) + "] = " + std::to_string(capacity[e]) + "]");
    }
  }
  return v;
}

// A per-edge integer array in [0, capacity_units[e]] (a ledger baseline, a search residual),
// validated on the unnarrowed values: the length (shape step), then every entry.
inline std::vector<Amount> edge_units(const py::object& obj, const Graph& g, const char* name) {
  return edge_amounts(obj, g, name, (std::string(name) + " (units)").c_str()).narrow<Amount>();
}

// A node set's ids, non-empty and range-checked before narrowing (detail::node_set's messages).
inline std::vector<NodeId> checked_ids(const WideInts& ids, std::int32_t num_nodes, const char* name) {
  if (ids.size() == 0) throw py::value_error(detail::empty_set_error(name));
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (!ids.in(i, 0, static_cast<std::int64_t>(num_nodes) - 1)) {
      throw py::value_error(detail::node_range_error(name, ids.text(i), num_nodes));
    }
  }
  return ids.narrow<NodeId>();
}

// A node set: one NodeId or a sequence of them.
inline std::vector<NodeId> node_ids(const py::object& obj, std::int32_t num_nodes, const char* name) {
  if (!is_bool(obj) && PyIndex_Check(obj.ptr()) && !py::isinstance<py::array>(obj)) {
    py::list one;
    one.append(obj);
    return checked_ids(wide_integers(one, name), num_nodes, name);
  }
  return checked_ids(wide_integers(obj, name), num_nodes, name);
}

} // namespace netgraph::core::bindings
