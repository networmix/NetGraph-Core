/*
  The float converter of the resolution parameters. It takes exactly a number: a float (a Python
  float, read by its stored value, or a NumPy floating value, converted to double) or an integer (a Python int or a
  NumPy integer, read exactly). Every other type is refused, and the call fails with pybind11's
  TypeError for incompatible arguments: a bool (Python's or numpy.bool_, by type identity), a
  fractions.Fraction, a decimal.Decimal, any object known only by __float__ or __index__. An
  integer beyond the float range is a value: it becomes +inf or -inf, which the resolution's own
  check refuses with its ValueError. An exception of the value's own conversion propagates as
  raised.
*/
#pragma once

#include <pybind11/pybind11.h>

#include <limits>

#include "strict_bool.hpp"

namespace netgraph::core::bindings {

// The double value of a parameter that takes concrete Python and NumPy numbers.
struct StrictFloat {
  double value {0.0};
  constexpr operator double() const noexcept { return value; }  // NOLINT(google-explicit-constructor)
};

} // namespace netgraph::core::bindings

namespace pybind11::detail {

template <>
struct type_caster<netgraph::core::bindings::StrictFloat> {
  // The input names exactly the concrete number types accepted by load().
  PYBIND11_TYPE_CASTER(netgraph::core::bindings::StrictFloat, io_name("int | float | numpy.integer | numpy.floating", "float"));

  bool load(handle src, bool convert) {
    PyObject* o = src.ptr();
    if (PyBool_Check(o) || reinterpret_cast<PyObject*>(Py_TYPE(o)) == netgraph::core::bindings::numpy_bool_type()) return false;
    const auto& [np_integer, np_floating] = netgraph::core::bindings::numpy_number_types();
    if (PyFloat_Check(o) || netgraph::core::bindings::has_type(o, np_floating)) {
      const double x = PyFloat_AsDouble(o);
      if (x == -1.0 && PyErr_Occurred() != nullptr) throw error_already_set();
      value.value = x;
      return true;
    }
    if (!convert || (PyLong_Check(o) == 0 && !netgraph::core::bindings::has_type(o, np_integer))) return false;
    const object idx = reinterpret_steal<object>(PyNumber_Index(o));
    if (!idx) throw error_already_set();
    const double x = PyLong_AsDouble(idx.ptr());
    if (x == -1.0 && PyErr_Occurred() != nullptr) {
      if (PyErr_ExceptionMatches(PyExc_OverflowError) == 0) throw error_already_set();
      PyErr_Clear();
      const int sign = PyObject_RichCompareBool(idx.ptr(), int_(0).ptr(), Py_GT);
      if (sign < 0) throw error_already_set();
      value.value = sign == 1 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
      return true;
    }
    value.value = x;
    return true;
  }

  static handle cast(netgraph::core::bindings::StrictFloat f, return_value_policy policy, handle parent) {
    return type_caster<double>::cast(f.value, policy, parent);
  }
};

} // namespace pybind11::detail
