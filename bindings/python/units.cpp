// to_units and as_float: physical values to integer units and back.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/detail/unit_conversion.hpp"
#include "integers.hpp"
#include "strict_float.hpp"
#include "values.hpp"
#include "native_call.hpp"
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
    const auto doubles = numpy_conversion<double, py::array::c_style | py::array::forcecast>(given);
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
        const int k = resolution_exponent(resolution);  // before the shape: an empty array is no exemption
        const PhysicalValues v = physical_values(x);
        if (v.scalar) return int_object(to_units(v.values[0], resolution));
        // At this scale the nearest double below 2^62 is hundreds of units
        // away: rounding cannot cross the boundary. Validate with the physical
        // bound, then leave the scaling/rounding kernel to the detached phase.
        const double limit = 4611686018427387904.0 * static_cast<double>(resolution);
        for (double x : v.values) {
          if (!(x >= 0.0 && x < limit)) (void)to_units(x, resolution);
        }
        auto result = native_call([&] {
          std::vector<std::int64_t> values(v.values.size());
          for (std::size_t i = 0; i < v.values.size(); ++i) {
            values[i] = static_cast<Amount>(detail::rounded_units(v.values[i], k));
          }
          return values;
        });
        return int64_array(result);
      },
      py::arg("x"), py::arg("resolution"),
      "Physical value(s) to units of resolution q. Each input becomes binary64 b: float32 widens\n"
      "exactly; integers become the nearest double, ties to even. Return floor(b/q + 1/2),\n"
      "evaluated exactly, never by floating addition of 1/2 or by reading the printed decimal.\n"
      "Ties are b=(2n+1)*q/2, 0 <= n < 2^52; each gives n+1, never the even neighbour. Scaling\n"
      "by q=2^k is exact while normal; underflow only affects values rounding to zero. Below\n"
      "2^52 scaled units the fractional comparison is exact; above it doubles are integral.\n"
      "A scalar or 0-d array gives an int; a sequence or 1-d array gives an int64 array. Negative\n"
      "zero gives 0. NaN, infinity, negative values, an integer too large for NumPy and results\n"
      "at or above 2^62 units raise ValueError. A bool, string, wrong dtype/rank or type other\n"
      "than int, float or NumPy numbers raises TypeError; int/float subclasses use stored values.\n"
      "q must be a power of two in [2^-32, 2^32]: a float, NumPy float or integer (another type\n"
      "raises TypeError, an integer beyond the float range ValueError).\n"
      "Every result u is exactly representable, so as_float(u,q) adds no error: against b,\n"
      "-q/2 < as_float(u,q)-b <= q/2, with +q/2 exactly at ties. Against an original argument a\n"
      "changed by promotion, add b-a (2^53+1 at q=1 loses 1 before quantization). N separately\n"
      "rounded values have total error in (-N*q/2,N*q/2]; no sign is promised. With symmetric\n"
      "non-tie fractional parts, mean error is (q/2)*P(tie); continuous uniform parts have none.\n"
      "Round a demand once, then split whole units. 2^53 is no decimal precision threshold: a\n"
      "decimal is exact only if its reduced denominator is a power of two and it fits binary64.\n"
      "A decimal and its parsed double agree only in the same half-up cell; 5501195.1796 at\n"
      "q=2^-21 gives 11536842473289, whereas its decimal reading would give 11536842473288.\n"
      "CostStretch deliberately reads printed decimals for ratios; physical amounts use b.");

  m.def(
      "as_float",
      [](const py::object& x, StrictFloat resolution) -> py::object {
        const int k = resolution_exponent(resolution);  // the resolution before the shape
        if (!is_bool(x) &&
            (PyLong_Check(x.ptr()) != 0 || has_type(x.ptr(), numpy_number_types().first))) {
          return py::float_(amount_as_float(x, k, resolution));
        }
        // A sequence, a NumPy array of any rank included: wide_integers refuses a rank other than 1.
        if (!is_invalid_sequence(x)) {
          const auto v = wide_integers(x, "as_float: x");
          std::vector<double> values(v.size());
          // Arbitrary Python integer arithmetic remains attached; the ordinary
          // fixed-width conversion loop works only on the copied native values.
          for (std::size_t i = 0; i < v.size(); ++i) {
            if (v.sign(i) != 0) values[i] = amount_as_float(v.source[py::int_(i)], k, resolution);
          }
          native_call([&] {
            for (std::size_t i = 0; i < v.size(); ++i) {
              if (v.sign(i) == 0) values[i] = detail::scaled_units(v.value[i], k);
            }
          });
          py::array_t<double> out(static_cast<py::ssize_t>(values.size()));
          std::copy(values.begin(), values.end(), out.mutable_data());
          return std::move(out);
        }
        return py::float_(amount_as_float(x, k, resolution));
      },
      py::arg("x"), py::arg("resolution"),
      "Units to physical value(s): the exact product x*q correctly rounded to binary64, nearest\n"
      "with ties to even. x is any signed, unbounded integer: an int or NumPy integer gives a\n"
      "float, an integer vector a float64 array. An int or NumPy integer, a subclass included,\n"
      "is never read as a sequence. An integer position accepts __index__, never\n"
      "a bool; other types raise TypeError. q is a power of two in [2^-32, 2^32] (float, NumPy\n"
      "float or integer). The result is finite iff abs(x*q) < 2^1024 - 2^970; otherwise raise\n"
      "ValueError('as_float: <x> units times the resolution exceed the largest finite double').\n"
      "Exact for abs(x) <= 2^53 inclusive. Above 2^53 an integer of binary length L is exact\n"
      "iff divisible by 2^(L-53); 2^53+1 rounds to 2^53, 2^53+3 to 2^53+4, times q.\n"
      "to_units(as_float(x,q),q)==x holds for x in [0,2^53] at every q, and for larger x in\n"
      "[0,2^62) exactly when x is a double. The first failure is 2^53+1; 2^62-1 rounds to\n"
      "2^62*q and reconversion raises the to_units domain ValueError. A successful to_units\n"
      "result is always exact here; its error against the promoted binary64 input is in\n"
      "(-q/2,q/2]. Choose q before Graph construction, then use Graph.resolution().");
}

} // namespace netgraph::core::bindings
