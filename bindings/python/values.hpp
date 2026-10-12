/*
  Small conversions shared by the registration files: read-only NumPy views that
  keep their owner alive, the resolution's check, exact counters, int64 arrays, FlowIds,
  (EdgeId, units) pairs and the type checks of option arguments.
*/
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/types.hpp"
#include "integers.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

// Read-only NumPy view of a span owned by `owner`.
template <typename T>
py::array readonly_view(std::span<const T> s, const py::handle& owner, const py::dtype& dtype) {
  static constexpr T empty_element {};
  py::array arr(dtype, {static_cast<py::ssize_t>(s.size())}, {static_cast<py::ssize_t>(sizeof(T))},
                s.empty() ? &empty_element : s.data(), owner);
  py::detail::array_proxy(arr.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
  return arr;
}

template <typename T>
py::array readonly_view(std::span<const T> s, const py::handle& owner) {
  return readonly_view(s, owner, py::dtype::of<T>());
}

// Validates a resolution (a power of two in [2^-32, 2^32]) and returns its exponent.
inline int resolution_exponent(double resolution) {
  (void)as_float(0, resolution);  // throws the resolution message
  return std::ilogb(resolution);
}

inline py::object counter_int(const Counter128& c) {
  py::int_ hi(c.hi);
  py::int_ lo(c.lo);
  return hi.attr("__lshift__")(py::int_(64)).attr("__or__")(lo);
}

inline py::array_t<std::int64_t> int64_array(const std::vector<std::int64_t>& v) {
  py::array_t<std::int64_t> out(static_cast<py::ssize_t>(v.size()));
  std::copy(v.begin(), v.end(), out.mutable_data());
  return out;
}

inline FlowId flow_id_of(const py::object& obj) {
  if (!py::isinstance<py::tuple>(obj) || py::len(obj) != 2) {
    throw py::value_error("a FlowId is a tuple (ledger uid, id)");
  }
  auto t = py::reinterpret_borrow<py::tuple>(obj);
  py::object uid = python_int(t[0], "FlowId ledger uid");
  py::object id = python_int(t[1], "FlowId id");
  int o1 = 0, o2 = 0;
  const unsigned long long u = PyLong_AsUnsignedLongLong(uid.ptr());
  if (PyErr_Occurred()) {
    // Only the range error means no such ledger; another failure (memory running out while that
    // error is raised) propagates as raised.
    if (!PyErr_ExceptionMatches(PyExc_OverflowError)) throw py::error_already_set();
    PyErr_Clear();
    o1 = 1;
  }
  const long long i = PyLong_AsLongLongAndOverflow(id.ptr(), &o2);
  if (o1 != 0 || o2 != 0 || i < 0) {
    throw py::value_error("FlowId (" + std::string(py::str(uid)) + ", " + std::string(py::str(id)) + ") is retired or unknown");
  }
  return FlowId {static_cast<std::uint64_t>(u), static_cast<std::int64_t>(i)};
}

inline py::list pairs_list(const std::vector<std::pair<EdgeId, Amount>>& v) {
  py::list out;
  for (const auto& [e, x] : v) out.append(py::make_tuple(e, x));
  return out;
}

inline OutputOptions outputs_of(const py::object& obj) {
  if (obj.is_none()) return OutputOptions {};
  if (!py::isinstance<OutputOptions>(obj)) {
    throw py::type_error("outputs must be an OutputOptions (got " + std::string(py::str(py::type::of(obj).attr("__name__"))) +
                         ")");
  }
  return obj.cast<OutputOptions>();
}

template <class Options>
Options options_of(const py::object& obj, const char* name) {
  if (!py::isinstance<Options>(obj)) {
    throw py::type_error(std::string(name) + " expected (got " + std::string(py::str(py::type::of(obj).attr("__name__"))) +
                         ")");
  }
  return obj.cast<Options>();
}

} // namespace netgraph::core::bindings
