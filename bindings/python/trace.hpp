/*
  The trace of a traced placement as Python values (the private hooks only).
*/
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "netgraph/core/detail/hooks_flag.hpp"

#if NETGRAPH_CORE_PRIVATE_HOOKS
#include <algorithm>
#include <string>

#include "netgraph/core/detail/trace.hpp"
#include "values.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

// A trace value as Python objects (map keys that are decimal integers become ints: the maps
// keyed by source NodeId), each made through a checked C API call: a failed allocation raises
// MemoryError (pybind11's int, str, list and dict constructors report it as RuntimeError).
inline py::object trace_obj(const detail::Trace& t) {
  using K = detail::Trace::Kind;
  switch (t.kind) {
    case K::Null: return py::none();
    case K::Int: return checked(PyLong_FromLongLong(t.i));
    case K::Str: return checked(PyUnicode_FromStringAndSize(t.s.data(), static_cast<Py_ssize_t>(t.s.size())));
    case K::List: {
      py::object out = checked(PyList_New(0));
      for (const auto& x : t.items) {
        if (PyList_Append(out.ptr(), trace_obj(x).ptr()) != 0) throw py::error_already_set();
      }
      return out;
    }
    case K::Map: {
      py::object out = checked(PyDict_New());
      for (const auto& [k, x] : t.fields) {
        const bool digits = !k.empty() && std::all_of(k.begin(), k.end(), [](char ch) { return ch >= '0' && ch <= '9'; });
        const py::object key = digits ? checked(PyLong_FromLongLong(std::stoll(k)))
                                      : checked(PyUnicode_FromStringAndSize(k.data(), static_cast<Py_ssize_t>(k.size())));
        if (PyDict_SetItem(out.ptr(), key.ptr(), trace_obj(x).ptr()) != 0) throw py::error_already_set();
      }
      return out;
    }
  }
  return py::none();
}

} // namespace netgraph::core::bindings

#endif
