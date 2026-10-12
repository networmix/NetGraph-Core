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

namespace py = pybind11;

namespace netgraph::core::bindings {

// A trace value as Python objects (map keys that are decimal integers become ints: the maps
// keyed by source NodeId).
inline py::object trace_obj(const detail::Trace& t) {
  using K = detail::Trace::Kind;
  switch (t.kind) {
    case K::Null: return py::none();
    case K::Bool: return py::bool_(t.b);
    case K::Int: return py::int_(t.i);
    case K::Str: return py::str(t.s);
    case K::List: {
      py::list out;
      for (const auto& x : t.items) out.append(trace_obj(x));
      return out;
    }
    case K::Map: {
      py::dict out;
      for (const auto& [k, x] : t.fields) {
        const bool digits = !k.empty() && std::all_of(k.begin(), k.end(), [](char ch) { return ch >= '0' && ch <= '9'; });
        if (digits) {
          out[py::int_(std::stoll(k))] = trace_obj(x);
        } else {
          out[py::str(k)] = trace_obj(x);
        }
      }
      return out;
    }
  }
  return py::none();
}

} // namespace netgraph::core::bindings

#endif
