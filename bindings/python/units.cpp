// to_units and as_float: physical values to integer units and back.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "integers.hpp"
#include "strict_float.hpp"
#include "values.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

constexpr const char* kNotPhysical = "to_units takes physical numbers (a float or an array of floats)";
constexpr const char* kNotOneDimensional = "to_units takes a float or a 1-D array of floats";

std::string dtype_text(const py::array& a) {
  return "to_units takes physical numbers (got dtype " + std::string(py::str(a.dtype())) + ")";
}

enum class Kind { Number, Wide, Bool, Other };

// An int's value as NumPy reads it, read from the int itself (a subclass's own methods are never
// called): through int64 or uint64, rounded to the nearest double (Number); beyond 64 bits NumPy
// holds it as an object (Wide: a value outside the domain).
Kind int_value(PyObject* o, double& value) {
  int overflow = 0;
  const long long v = PyLong_AsLongLongAndOverflow(o, &overflow);
  if (v == -1 && PyErr_Occurred() != nullptr) throw py::error_already_set();
  if (overflow == 0) {
    value = static_cast<double>(v);
    return Kind::Number;
  }
  if (overflow > 0) {
    const unsigned long long u = PyLong_AsUnsignedLongLong(o);
    if (u != static_cast<unsigned long long>(-1) || PyErr_Occurred() == nullptr) {
      value = static_cast<double>(u);
      return Kind::Number;
    }
    if (PyErr_ExceptionMatches(PyExc_OverflowError) == 0) throw py::error_already_set();
    PyErr_Clear();
  }
  return Kind::Wide;
}

// One scalar under the rule, in this order: a bool (Python's or numpy.bool_) is Bool; an int or a
// float, a subclass included, is a Number (Wide beyond 64 bits) read by its stored value (PyLong_*
// and PyFloat_AS_DOUBLE read the object itself, never its __int__, __float__ or sequence methods);
// a NumPy integer or floating scalar is a Number; anything else is Other. Every type is the
// object's real one (never its __class__).
Kind scalar_kind(const py::handle& x, double& value) {
  PyObject* o = x.ptr();
  if (is_bool(x)) return Kind::Bool;
  if (PyLong_Check(o) != 0) return int_value(o, value);
  if (PyFloat_Check(o) != 0) {
    value = PyFloat_AS_DOUBLE(o);
    return Kind::Number;
  }
  const auto& [np_integer, np_floating] = numpy_number_types();
  if (has_type(o, np_integer)) {
    const auto exact = checked(PyNumber_Index(o));
    return int_value(exact.ptr(), value);
  }
  if (has_type(o, np_floating)) {
    value = PyFloat_AsDouble(o);
    if (value == -1.0 && PyErr_Occurred() != nullptr) throw py::error_already_set();
    return Kind::Number;
  }
  return Kind::Other;
}

// Whether NumPy reads `o` without running code of the caller's: a str, bytes, bool, None, an exact
// int or float, an exact NumPy scalar or array (its type named "numpy.*"), or a list or tuple of
// those (nested at most 32 deep).
bool read_by_numpy_alone(PyObject* o, int depth = 0) {
  if (PyUnicode_CheckExact(o) || PyBytes_CheckExact(o) || PyBool_Check(o) || o == Py_None || PyLong_CheckExact(o) ||
      PyFloat_CheckExact(o)) {
    return true;
  }
  if (std::string_view(Py_TYPE(o)->tp_name).starts_with("numpy.")) return true;
  if (depth >= 32 || (PyList_CheckExact(o) == 0 && PyTuple_CheckExact(o) == 0)) return false;
  const Py_ssize_t n = PySequence_Fast_GET_SIZE(o);
  PyObject** items = PySequence_Fast_ITEMS(o);
  for (Py_ssize_t i = 0; i < n; ++i) {
    if (!read_by_numpy_alone(items[i], depth + 1)) return false;
  }
  return true;
}

// The refusal of an input that holds something other than numbers, always a TypeError, worded as
// NumPy reads it (the texts these inputs always had): a non-numeric dtype names it; a numeric
// reading of more than one dimension is the shape text; any other numeric reading the scalar text.
// An input NumPy would read through the caller's code (an object's own __array__, a subclass) is
// not given to NumPy: it holds an object, NumPy's `got dtype object`.
[[noreturn]] void refuse_physical(const py::handle& x) {
  if (!read_by_numpy_alone(x.ptr())) throw py::type_error("to_units takes physical numbers (got dtype object)");
  const auto np = py::module_::import("numpy");
  const auto raw = py::reinterpret_steal<py::array>(checked(PyObject_CallOneArg(np.attr("asarray").ptr(), x.ptr())).release());
  const char kind = raw.dtype().kind();
  if (kind != 'f' && kind != 'i' && kind != 'u') throw py::type_error(dtype_text(raw));
  if (raw.ndim() > 1) throw py::type_error(kNotOneDimensional);
  throw py::type_error(kNotPhysical);
}

// to_units' x as doubles: one scalar, or the entries of a 1-D input.
struct PhysicalValues {
  bool scalar {false};
  std::vector<double> values;
};

// The one classifier of to_units' x, in this order:
//   1. a bool, a str or a bytes: TypeError (the scalar text);
//   2. an int or a float, a subclass included, by its stored value; a NumPy integer or floating
//      scalar (scalar_kind): one number (an int beyond 64 bits: ValueError);
//   3. a NumPy array: an integer or floating dtype (else TypeError, the dtype text), of rank 0 (one
//      number) or 1 (its entries; another rank TypeError, the shape text);
//   4. a Python sequence (a memoryview of rank 1; a str or bytes never reaches here), each entry
//      through scalar_kind: any entry of another type refuses it as NumPy words it; then NumPy's
//      dtype: bools alone TypeError with `got dtype bool`, an int beyond 64 bits `got dtype
//      object` (ValueError among numbers alone, TypeError beside a bool); then bools among
//      numbers TypeError (the scalar text: NumPy would read them as numbers); numbers alone are
//      its values;
//   5. anything else: TypeError, as NumPy words it (refuse_physical).
// An int beyond 64 bits among numbers is a value outside the domain (ValueError, NumPy's `got dtype
// object`), decided after every entry is classified. No Python object is made but the checked
// copies of NumPy integers, NumPy's array of a given array and the refusal's reading.
PhysicalValues physical_values(const py::handle& x) {
  PhysicalValues out;
  if (py::isinstance<py::str>(x) || py::isinstance<py::bytes>(x)) throw py::type_error(kNotPhysical);
  constexpr const char* wide = "to_units takes physical numbers (got dtype object)";
  double value = 0.0;
  switch (scalar_kind(x, value)) {
    case Kind::Bool:
      throw py::type_error(kNotPhysical);
    case Kind::Wide:
      throw py::value_error(wide);
    case Kind::Number:
      out.scalar = true;
      out.values.push_back(value);
      return out;
    case Kind::Other:
      break;
  }
  if (py::isinstance<py::array>(x)) {
    const auto given = py::reinterpret_borrow<py::array>(x);
    const char kind = given.dtype().kind();
    if (kind != 'f' && kind != 'i' && kind != 'u') throw py::type_error(dtype_text(given));
    if (given.ndim() > 1) throw py::type_error(kNotOneDimensional);
    const py::array_t<double, py::array::c_style | py::array::forcecast> doubles(given);
    out.scalar = given.ndim() == 0;
    out.values.assign(doubles.data(), doubles.data() + doubles.size());
    return out;
  }
  if (PySequence_Check(x.ptr()) == 0) refuse_physical(x);
  if (PyMemoryView_Check(x.ptr()) != 0 && PyMemoryView_GET_BUFFER(x.ptr())->ndim != 1) {
    if (PyMemoryView_GET_BUFFER(x.ptr())->ndim > 1) throw py::type_error(kNotOneDimensional);
    refuse_physical(x);
  }
  // A list or a tuple is read in place; any other sequence by index, as NumPy reads it (its
  // own __len__ and __getitem__; no iterator, whose machinery may absorb a failure).
  const bool in_place = PyList_Check(x.ptr()) != 0 || PyTuple_Check(x.ptr()) != 0;
  const Py_ssize_t size = PySequence_Size(x.ptr());
  if (size < 0) throw py::error_already_set();
  const auto n = static_cast<std::size_t>(size);
  const auto entry = [&](std::size_t i) {
    if (in_place) return py::reinterpret_borrow<py::object>(PySequence_Fast_ITEMS(x.ptr())[i]);
    return checked(PySequence_GetItem(x.ptr(), static_cast<Py_ssize_t>(i)));
  };
  // A refusal is worded on the entries as a list: the sequence itself, or a checked copy.
  const auto refuse = [&]() {
    if (in_place) refuse_physical(x);
    const py::object copy = checked(PyList_New(size));
    for (std::size_t i = 0; i < n; ++i) PyList_SET_ITEM(copy.ptr(), static_cast<Py_ssize_t>(i), entry(i).release().ptr());
    refuse_physical(copy);
  };
  out.values.reserve(n);
  std::size_t bools = 0;
  std::size_t wides = 0;
  for (std::size_t i = 0; i < n; ++i) {
    switch (scalar_kind(entry(i), value)) {
      case Kind::Number:
        out.values.push_back(value);
        break;
      case Kind::Wide:
        ++wides;
        break;
      case Kind::Bool:
        ++bools;
        break;
      case Kind::Other:
        refuse();
    }
  }
  // NumPy's order: its dtype first (bools alone: bool; an int beyond 64 bits: object, a value
  // error only among numbers alone), then a bool among numbers.
  if (bools == n && n > 0) throw py::type_error("to_units takes physical numbers (got dtype bool)");
  if (wides > 0) {
    if (bools > 0) throw py::type_error(wide);
    throw py::value_error(wide);
  }
  if (bools > 0) throw py::type_error(kNotPhysical);
  return out;
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
      [](const py::object& x, StrictFloat resolution) -> py::object {
        (void)resolution_exponent(resolution);  // before the shape: an empty array is no exemption
        const PhysicalValues v = physical_values(x);
        if (v.scalar) return int_object(to_units(v.values[0], resolution));
        py::array_t<std::int64_t> out(static_cast<py::ssize_t>(v.values.size()));
        auto o = out.mutable_unchecked<1>();
        for (std::size_t i = 0; i < v.values.size(); ++i) o(static_cast<py::ssize_t>(i)) = to_units(v.values[i], resolution);
        return std::move(out);
      },
      py::arg("x"), py::arg("resolution"),
      "Physical value(s) to units of the resolution q, rounding half up: an int for a scalar, an\n"
      "int64 array for an array. NaN, infinity, negative values, an integer too large for NumPy and\n"
      "results at or above 2^62 units raise ValueError; a bool (also as an entry), a string or an array\n"
      "of another dtype or rank, or any type but int, float and the NumPy numbers (a Decimal, a\n"
      "Fraction, an object known by __float__, __index__ or __array__) TypeError; an int or float\n"
      "subclass is read by its value. q must be a power of two in [2^-32, 2^32]: a float, a NumPy\n"
      "float or an integer (any other type is a TypeError; an integer beyond the float range a\n"
      "ValueError).");

  m.def(
      "as_float",
      [](const py::object& x, StrictFloat resolution) -> py::object {
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
      "doubles raises ValueError. x is an integer position (any object with __index__, never a bool;\n"
      "anything else is a TypeError); q must be a power of two in [2^-32, 2^32] (a float, a NumPy float\n"
      "or an integer).");
}

} // namespace netgraph::core::bindings
