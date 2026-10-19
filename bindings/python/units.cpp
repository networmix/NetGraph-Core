// to_units and as_float: physical values to integer units and back.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <string>
#include <utility>

#include "netgraph/core/graph.hpp"
#include "integers.hpp"
#include "values.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// Physical numbers as a float64 array (or a 0-d array for a scalar); bools and strings rejected.
py::array_t<double> physical_array(const py::handle& x) {
  if (is_bool(x) || py::isinstance<py::str>(x) || py::isinstance<py::bytes>(x)) {
    throw py::value_error("to_units takes physical numbers (a float or an array of floats)");
  }
  auto np = py::module_::import("numpy");
  py::array raw = np.attr("asarray")(x);
  const char kind = raw.dtype().kind();
  if (kind != 'f' && kind != 'i' && kind != 'u') {
    throw py::value_error("to_units takes physical numbers (got dtype " + std::string(py::str(raw.dtype())) + ")");
  }
  return py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(raw);
}

// x * 2^k for an unbounded Python int, correctly rounded: Python's int shift and int true
// division are exact or correctly rounded, so only a result beyond the finite doubles fails.
double scaled_python_int(const py::object& x, int k) {
  try {
    py::object r;
    if (k >= 0) {
      r = py::float_(x.attr("__lshift__")(py::int_(k)));
    } else {
      r = x.attr("__truediv__")(py::int_(1).attr("__lshift__")(py::int_(-k)));
    }
    return r.cast<double>();
  } catch (py::error_already_set& e) {
    if (e.matches(PyExc_OverflowError)) {
      throw py::value_error("as_float: " + std::string(py::str(x)) +
                            " units times the resolution exceed the largest finite double");
    }
    throw;
  }
}

// One Amount to a physical double; the resolution's exponent k is already validated.
double amount_as_float(const py::handle& item, int k, double resolution) {
  py::object idx = python_int(item, "as_float: x");
  int overflow = 0;
  const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
  if (overflow == 0) return as_float(static_cast<std::int64_t>(v), resolution);
  return scaled_python_int(idx, k);
}

} // namespace

void register_units(py::module_& m) {
  m.def(
      "to_units",
      [](const py::object& x, double resolution) -> py::object {
        (void)resolution_exponent(resolution);  // before the shape: an empty array is no exemption
        auto arr = physical_array(x);
        if (arr.ndim() == 0) return py::int_(to_units(*arr.data(), resolution));
        if (arr.ndim() != 1) throw py::value_error("to_units takes a float or a 1-D array of floats");
        py::array_t<std::int64_t> out(arr.size());
        auto in = arr.unchecked<1>();
        auto o = out.mutable_unchecked<1>();
        for (py::ssize_t i = 0; i < in.shape(0); ++i) o(i) = to_units(in(i), resolution);
        return std::move(out);
      },
      py::arg("x"), py::arg("resolution"),
      "Physical value(s) to units of the resolution q, rounding half up: an int for a scalar, an\n"
      "int64 array for an array. NaN, infinity, negative values and results at or above 2^62\n"
      "units raise ValueError; q must be a power of two in [2^-32, 2^32].");

  m.def(
      "as_float",
      [](const py::object& x, double resolution) -> py::object {
        const int k = resolution_exponent(resolution);  // the resolution before the shape
        // A sequence, a NumPy array of any rank included: wide_integers refuses a rank other than 1.
        if (!is_invalid_sequence(x)) {
          const auto v = wide_integers(x, "as_float: x");
          py::array_t<double> out(static_cast<py::ssize_t>(v.size()));
          auto o = out.mutable_unchecked<1>();
          for (std::size_t i = 0; i < v.size(); ++i) {
            o(static_cast<py::ssize_t>(i)) = v.sign(i) == 0 ? as_float(v.value[i], resolution)
                                                            : amount_as_float(v.source[py::int_(i)], k, resolution);
          }
          return std::move(out);
        }
        return py::float_(amount_as_float(x, k, resolution));
      },
      py::arg("x"), py::arg("resolution"),
      "Units to physical value(s): x * q, correctly rounded (exact below 2^53 units). An int, of\n"
      "any size, gives a float, an integer array a float64 array; a result beyond the finite\n"
      "doubles raises ValueError; q must be a power of two in [2^-32, 2^32].");
}

} // namespace netgraph::core::bindings
