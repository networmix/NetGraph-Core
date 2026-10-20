/*
  Small conversions shared by the registration files: read-only NumPy views that
  keep their owner alive, the resolution's check, exact counters, int64 arrays, FlowIds,
  (EdgeId, units) pairs and the type checks of option arguments.
*/
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/typing.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <optional>
#include <span>
#include <string>
#include <type_traits>
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

// A new reference from the C API, or the pending Python error raised (pybind11's own int and
// tuple constructors fail on a failed allocation with an unchecked internal error instead).
inline py::object checked(PyObject* o) {
  if (o == nullptr) throw py::error_already_set();
  return py::reinterpret_steal<py::object>(o);
}

// hi * 2^64 + lo as an exact Python int.
inline py::object counter_int(const Counter128& c) {
  py::object lo = checked(PyLong_FromUnsignedLongLong(c.lo));
  if (c.hi == 0) return lo;
  py::object hi = checked(PyLong_FromUnsignedLongLong(c.hi));
  py::object shift = checked(PyLong_FromLong(64));
  py::object high = checked(PyNumber_Lshift(hi.ptr(), shift.ptr()));
  return checked(PyNumber_Or(high.ptr(), lo.ptr()));
}

// The readers of fields and properties build every Python value below through the checked C API
// (a new reference, or the pending MemoryError raised). pybind11's own conversions report a failed
// allocation as another class: its return conversion of a C++ value as TypeError, its int, float,
// list and tuple constructors as RuntimeError.

// An integer (any integer type but bool) or a 128-bit counter as a Python int.
template <class T>
py::int_ int_object(T v) {
  static_assert(std::is_integral_v<T> && !std::is_same_v<T, bool>);
  PyObject* o = nullptr;
  if constexpr (std::is_signed_v<T>) {
    o = PyLong_FromLongLong(static_cast<long long>(v));
  } else {
    o = PyLong_FromUnsignedLongLong(static_cast<unsigned long long>(v));
  }
  return py::reinterpret_steal<py::int_>(checked(o).release());
}
inline py::int_ int_object(const Counter128& c) { return py::reinterpret_steal<py::int_>(counter_int(c).release()); }

inline py::float_ float_object(double v) { return py::reinterpret_steal<py::float_>(checked(PyFloat_FromDouble(v)).release()); }

// None, or the value as a Python int.
template <class T>
py::typing::Optional<py::int_> optional_int(const std::optional<T>& v) {
  return py::reinterpret_steal<py::typing::Optional<py::int_>>(v.has_value() ? int_object(*v).release() : py::none().release());
}

// The reader of an integer field, and of an optional one (None when absent).
template <class C, class T>
auto int_reader(T C::*field) {
  return [field](const C& c) { return int_object(c.*field); };
}
template <class C, class T>
auto optional_int_reader(std::optional<T> C::*field) {
  return [field](const C& c) { return optional_int(c.*field); };
}

// (a, b) as a tuple of two Python ints.
template <class A, class B>
py::typing::Tuple<py::int_, py::int_> int_pair(const A& a, const B& b) {
  py::object x = int_object(a);
  py::object y = int_object(b);
  PyObject* out = checked(PyTuple_New(2)).release().ptr();
  PyTuple_SET_ITEM(out, 0, x.release().ptr());
  PyTuple_SET_ITEM(out, 1, y.release().ptr());
  return py::reinterpret_steal<py::typing::Tuple<py::int_, py::int_>>(out);
}

// Integers as a list of Python ints.
template <class T>
py::typing::List<py::int_> ints_list(std::span<const T> v) {
  auto out = py::reinterpret_steal<py::typing::List<py::int_>>(checked(PyList_New(static_cast<Py_ssize_t>(v.size()))).release());
  for (std::size_t i = 0; i < v.size(); ++i) PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), int_object(v[i]).release().ptr());
  return out;
}

// (a, b) pairs as a list (Out: its pybind11 type, which names the stub's type) or a tuple of
// two-int tuples.
template <class Out = py::list, class A, class B>
Out pairs_list(const std::vector<std::pair<A, B>>& v) {
  auto out = py::reinterpret_steal<Out>(checked(PyList_New(static_cast<Py_ssize_t>(v.size()))).release());
  for (std::size_t i = 0; i < v.size(); ++i) {
    PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), int_pair(v[i].first, v[i].second).release().ptr());
  }
  return out;
}
template <class A, class B>
py::tuple pairs_tuple(const std::vector<std::pair<A, B>>& v) {
  auto out = py::reinterpret_steal<py::tuple>(checked(PyTuple_New(static_cast<Py_ssize_t>(v.size()))).release());
  for (std::size_t i = 0; i < v.size(); ++i) {
    PyTuple_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), int_pair(v[i].first, v[i].second).release().ptr());
  }
  return out;
}

// A tuple of the given objects (made by the caller), through a checked allocation.
inline py::tuple objects_tuple(std::initializer_list<py::object> items) {
  auto out = py::reinterpret_steal<py::tuple>(checked(PyTuple_New(static_cast<Py_ssize_t>(items.size()))).release());
  Py_ssize_t i = 0;
  for (const auto& item : items) PyTuple_SET_ITEM(out.ptr(), i++, item.inc_ref().ptr());
  return out;
}

// fraction_type(n, d) (fractions.Fraction): the two ints through checked allocations and the call
// through the C API, so a failure of the call raises as raised.
inline py::object fraction_object(const py::handle& fraction_type, std::int64_t n, std::int64_t d) {
  const py::object num = int_object(n);
  const py::object den = int_object(d);
  return checked(PyObject_CallFunctionObjArgs(fraction_type.ptr(), num.ptr(), den.ptr(), nullptr));
}

// A str through a checked allocation.
inline py::str str_object(const char* text) { return py::reinterpret_steal<py::str>(checked(PyUnicode_FromString(text)).release()); }

// The hash of `key` as a Python int: a failure of the hash raises as raised, and the int is made
// through a checked allocation (pybind11's return conversion of a hash reports a failed
// allocation as TypeError).
inline py::int_ hash_int(const py::handle& key) {
  const Py_hash_t h = PyObject_Hash(key.ptr());
  if (h == -1 && PyErr_Occurred() != nullptr) throw py::error_already_set();
  return int_object(static_cast<long long>(h));
}

// Integers as a tuple of Python ints.
inline py::tuple ints_tuple(std::span<const std::int32_t> v) {
  auto out = py::reinterpret_steal<py::tuple>(checked(PyTuple_New(static_cast<Py_ssize_t>(v.size()))).release());
  for (std::size_t i = 0; i < v.size(); ++i) {
    PyTuple_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), checked(PyLong_FromLong(v[i])).release().ptr());
  }
  return out;
}

// A tuple of Python ints, each through a checked allocation.
inline py::tuple amounts_tuple(std::span<const Amount> v) {
  auto out = py::reinterpret_steal<py::tuple>(checked(PyTuple_New(static_cast<Py_ssize_t>(v.size()))).release());
  for (std::size_t i = 0; i < v.size(); ++i) {
    PyTuple_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), checked(PyLong_FromLongLong(v[i])).release().ptr());
  }
  return out;
}

inline py::array_t<std::int64_t> int64_array(const std::vector<std::int64_t>& v) {
  py::array_t<std::int64_t> out(static_cast<py::ssize_t>(v.size()));
  std::copy(v.begin(), v.end(), out.mutable_data());
  return out;
}

inline FlowId flow_id_of(const py::object& obj) {
  // Not a tuple is a wrong type, a tuple of another length a wrong value: one text.
  if (!py::isinstance<py::tuple>(obj)) throw py::type_error("a FlowId is a tuple (ledger uid, id)");
  if (py::len(obj) != 2) throw py::value_error("a FlowId is a tuple (ledger uid, id)");
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


// An amount argument: an integer number of units in [0, 2^62) or +inf (float('inf')); nullopt is
// +inf. Every other argument raises the one `text`, with the class of its cause: a type the position
// never takes (a bool, a string, None, another non-number) TypeError; a float other than +inf (NaN,
// -inf, a finite one) or an integer out of range ValueError. An exception of the object's own
// __index__ propagates as raised.
inline std::optional<Amount> units_or_inf(const py::object& value, const char* text) {
  if (PyFloat_Check(value.ptr()) || has_type(value.ptr(), numpy_number_types().second)) {
    const double x = PyFloat_AsDouble(value.ptr());
    if (x == -1.0 && PyErr_Occurred() != nullptr) throw py::error_already_set();
    if (std::isinf(x) && x > 0) return std::nullopt;
    throw py::value_error(text);
  }
  if (is_bool(value) || !PyIndex_Check(value.ptr())) throw py::type_error(text);
  const py::object idx = py::reinterpret_steal<py::object>(PyNumber_Index(value.ptr()));
  if (!idx) throw py::error_already_set();
  int overflow = 0;
  const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
  if (overflow != 0 || v < 0 || v >= kDomainLimit) throw py::value_error(text);
  return static_cast<Amount>(v);
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

// A cost ceiling argument (DistributedTE's and KspOptions'): None, an AbsoluteCost or a CostStretch;
// any other object raises one TypeError.
inline std::optional<CostCeiling> cost_ceiling_of(const py::object& obj) {
  if (obj.is_none()) return std::nullopt;
  if (py::isinstance<AbsoluteCost>(obj)) return CostCeiling {obj.cast<AbsoluteCost>()};
  if (py::isinstance<CostStretch>(obj)) return CostCeiling {obj.cast<CostStretch>()};
  throw py::type_error("cost_ceiling must be AbsoluteCost, CostStretch or None");
}

} // namespace netgraph::core::bindings
