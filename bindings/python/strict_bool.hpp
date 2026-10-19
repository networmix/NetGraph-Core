/*
  The one bool converter of every bool parameter: exactly a Python bool or a numpy.bool_, both by
  type identity. pybind11's own bool caster in no-convert mode recognizes NumPy's bool by the name
  of its type, so a class named "numpy.bool_" passes it; this caster compares the type object.
*/
#pragma once

#include <pybind11/gil_safe_call_once.h>
#include <pybind11/pybind11.h>

namespace netgraph::core::bindings {

// A bool parameter's value (its signature says bool).
struct StrictBool {
  bool value {false};
  constexpr operator bool() const noexcept { return value; }  // NOLINT(google-explicit-constructor)
};

// numpy.bool_'s type object, looked up once.
inline PyObject* numpy_bool_type() {
  PYBIND11_CONSTINIT static pybind11::gil_safe_call_once_and_store<pybind11::object> storage;
  return storage.call_once_and_store_result([] { return pybind11::module_::import("numpy").attr("bool_"); }).get_stored().ptr();
}

} // namespace netgraph::core::bindings

namespace pybind11::detail {

template <>
struct type_caster<netgraph::core::bindings::StrictBool> {
  PYBIND11_TYPE_CASTER(netgraph::core::bindings::StrictBool, const_name("bool"));

  bool load(handle src, bool /*convert*/) {
    if (src.ptr() == Py_True || src.ptr() == Py_False) {
      value.value = src.ptr() == Py_True;
      return true;
    }
    if (reinterpret_cast<PyObject*>(Py_TYPE(src.ptr())) != netgraph::core::bindings::numpy_bool_type()) return false;
    const int truth = PyObject_IsTrue(src.ptr());
    if (truth < 0) throw error_already_set();
    value.value = truth == 1;
    return true;
  }

  static handle cast(netgraph::core::bindings::StrictBool b, return_value_policy /*policy*/, handle /*parent*/) {
    return handle(b.value ? Py_True : Py_False).inc_ref();
  }
};

} // namespace pybind11::detail
