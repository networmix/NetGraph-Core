/*
  Pybind11 module `_netgraph_core`: the graph on integer units and the private flow-core entry.

  - Graph accessors return read-only NumPy views that keep the graph alive (no per-call copy).
  - Integer inputs (capacities, costs, ids, views, caps) must be integers: Python ints, NumPy
    integer scalars or integer arrays; floats and bools are rejected with ValueError, never
    truncated. Values are read without narrowing (WideInts) and validated in the contract's
    order before any conversion to the C++ types (APPROACH 3, revision 25). Physical values
    enter only through to_units; the resolution is validated before anything else.
  - std::invalid_argument from Core maps to ValueError.
  - _place_global_units is the recorder's private entry point (APPROACH 2.2): it returns the
    payload as a dict of Python ints and lists.
*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

namespace py = pybind11;
using namespace netgraph::core;

namespace {

// Read-only NumPy view of a span owned by `owner`.
template <typename T>
py::array readonly_view(std::span<const T> s, const py::handle& owner, const py::dtype& dtype) {
  py::array arr(dtype, {static_cast<py::ssize_t>(s.size())}, {static_cast<py::ssize_t>(sizeof(T))},
                s.data(), owner);
  py::detail::array_proxy(arr.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
  return arr;
}

template <typename T>
py::array readonly_view(std::span<const T> s, const py::handle& owner) {
  return readonly_view(s, owner, py::dtype::of<T>());
}

bool is_bool(const py::handle& obj) {
  return PyBool_Check(obj.ptr()) || py::isinstance(obj, py::module_::import("numpy").attr("bool_"));
}

// The Python int of an integer-like object (a Python int or a NumPy integer scalar); bools and
// everything else are rejected.
py::object python_int(const py::handle& item, const std::string& what) {
  if (is_bool(item) || !PyIndex_Check(item.ptr())) {
    throw py::value_error(what + " must be an integer (got " +
                          std::string(py::str(py::type::of(item).attr("__name__"))) + ")");
  }
  py::object idx = py::reinterpret_steal<py::object>(PyNumber_Index(item.ptr()));
  if (!idx) throw py::error_already_set();
  return idx;
}

// One integer of [lo, hi] (scalar arguments).
std::int64_t integer_value(const py::handle& item, std::int64_t lo, std::int64_t hi, const std::string& what) {
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
// is checked by the caller in the contract's order, so the reported error does not depend on
// the representation of the input (APPROACH 3, revision 25).
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
  template <typename T>
  [[nodiscard]] std::vector<T> narrow() const {  // after every bound was checked
    return std::vector<T>(value.begin(), value.end());
  }
};

bool not_a_sequence(const py::handle& obj) {
  return py::isinstance<py::str>(obj) || py::isinstance<py::bytes>(obj) || !py::isinstance<py::sequence>(obj);
}

WideInts wide_integers(const py::handle& obj, const char* name) {
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
bool is_nested(const py::handle& item) {
  if (py::isinstance<py::array>(item)) return py::reinterpret_borrow<py::array>(item).ndim() >= 1;
  return !py::isinstance<py::str>(item) && !py::isinstance<py::bytes>(item) && PySequence_Check(item.ptr()) != 0;
}

// The shape step: the length of a rank-1 integer sequence argument. A NumPy array of another
// rank (0-D included), a nested sequence or a non-sequence is a shape defect. Entry types and
// values are not read here; their steps come later.
std::size_t sequence_length(const py::handle& obj, const char* name) {
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

// The unit-value step on the capacity argument whatever its shape (the shape step comes next):
// every integer entry it holds must lie in [0, 2^62). A NumPy array is read flattened; nested
// entries of a Python sequence and a non-integer non-sequence are left to the shape step.
void check_unit_values(const py::handle& capacity_units) {
  const auto out_of_domain = [](const std::string& what, const std::string& value) {
    return py::value_error(what + " = " + value + " is outside the units domain [0, 2^62)");
  };
  if (py::isinstance<py::array>(capacity_units)) {
    py::object flat = capacity_units.attr("ravel")();
    const auto w = wide_integers(flat, "capacity_units");
    for (std::size_t i = 0; i < w.size(); ++i) {
      if (!w.in(i, 0, kDomainLimit - 1)) throw out_of_domain("capacity_units[" + std::to_string(i) + "]", w.text(i));
    }
    return;
  }
  if (!not_a_sequence(capacity_units)) {
    auto seq = py::reinterpret_borrow<py::sequence>(capacity_units);
    for (std::size_t i = 0; i < seq.size(); ++i) {
      py::object item = seq[i];
      if (is_nested(item)) continue;
      const std::string what = "capacity_units[" + std::to_string(i) + "]";
      py::object idx = python_int(item, what);
      int overflow = 0;
      const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
      if (overflow != 0 || v < 0 || v >= kDomainLimit) throw out_of_domain(what, std::string(py::str(idx)));
    }
    return;
  }
  if (!is_bool(capacity_units) && PyIndex_Check(capacity_units.ptr())) {
    py::object idx = python_int(capacity_units, "capacity_units");
    int overflow = 0;
    const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
    if (overflow != 0 || v < 0 || v >= kDomainLimit) throw out_of_domain("capacity_units", std::string(py::str(idx)));
  }
}

// Validates a resolution (a power of two in [2^-32, 2^32]) and returns its exponent.
int resolution_exponent(double resolution) {
  (void)as_float(0, resolution);  // throws the resolution message
  return std::ilogb(resolution);
}

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

std::shared_ptr<StrictMultiDiGraph> make_graph(const py::object& num_nodes, const py::object& src,
                                               const py::object& dst, const py::object& capacity_units,
                                               const py::object& cost, double resolution) {
  // APPROACH 3 (revision 25): the resolution, the unit values, the shapes, zero or negative
  // costs, the remaining bounds (ids in range, then the sums); every value is read without
  // narrowing, and the messages are the constructor's, which re-checks the same order.
  (void)resolution_exponent(resolution);
  check_unit_values(capacity_units);
  // The shapes: num_nodes, then every array argument rank 1 (in signature order), then equal
  // lengths; no entry of src, dst or cost is read before the cost step.
  const auto n = integer_value(num_nodes, 0, std::numeric_limits<std::int32_t>::max(), "num_nodes");
  const std::size_t m = sequence_length(src, "src");
  const std::size_t md = sequence_length(dst, "dst");
  const std::size_t mu = sequence_length(capacity_units, "capacity_units");
  const std::size_t mc = sequence_length(cost, "cost");
  if (md != m || mu != m || mc != m) {
    throw py::value_error("src, dst, capacity_units and cost must have the same length (got " +
                          std::to_string(m) + ", " + std::to_string(md) + ", " + std::to_string(mu) + ", " +
                          std::to_string(mc) + ")");
  }
  const auto cap = wide_integers(capacity_units, "capacity_units");  // rank 1, values checked
  if (m > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    throw py::value_error("the number of edges exceeds INT32_MAX (array length)");
  }
  const auto w = wide_integers(cost, "cost");
  for (std::size_t i = 0; i < m; ++i) {
    if (w.sign(i) < 0 || (w.sign(i) == 0 && w.value[i] < 1)) {
      throw py::value_error("cost[" + std::to_string(i) + "] = " + w.text(i) + ": costs must be positive (cost >= 1)");
    }
  }
  const auto s = wide_integers(src, "src");
  const auto d = wide_integers(dst, "dst");
  for (std::size_t i = 0; i < m; ++i) {
    if (!s.in(i, 0, n - 1) || !d.in(i, 0, n - 1)) {
      throw py::value_error("edge " + std::to_string(i) + ": endpoint (" + s.text(i) + ", " + d.text(i) +
                            ") out of node range [0, " + std::to_string(n) + ")");
    }
  }
  constexpr auto kLimit = static_cast<std::uint64_t>(kDomainLimit);
  std::uint64_t cap_sum = 0;
  std::uint64_t cost_sum = 0;
  for (std::size_t i = 0; i < m; ++i) {
    cap_sum += static_cast<std::uint64_t>(cap.value[i]);  // each below 2^62: no wrap before the check
    if (cap_sum >= kLimit) {
      throw py::value_error("the sum of all capacities reaches the bound of the units domain (2^62 units)");
    }
    if (w.sign(i) > 0 || static_cast<std::uint64_t>(w.value[i]) >= kLimit - cost_sum) {
      throw py::value_error("the sum of all edge costs reaches the bound 2^62: path costs would overflow int64");
    }
    cost_sum += static_cast<std::uint64_t>(w.value[i]);
  }
  auto cap_v = cap.narrow<Amount>();
  auto cost_v = w.narrow<Cost>();
  auto src_v = s.narrow<NodeId>();
  auto dst_v = d.narrow<NodeId>();
  py::gil_scoped_release release;
  return std::make_shared<StrictMultiDiGraph>(static_cast<std::int32_t>(n), src_v, dst_v, cap_v, cost_v, resolution);
}

py::list int_list(const std::vector<std::int64_t>& v) {
  py::list out(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) out[i] = py::int_(v[i]);
  return out;
}

} // namespace

PYBIND11_MODULE(_netgraph_core, m, py::mod_gil_not_used()) {
  m.doc() = "NetGraph-Core: the graph on integer units and the minimum-cost flow core";

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
        const int k = resolution_exponent(resolution);  // before the shape (review-1a D3)
        const bool array = py::isinstance<py::array>(x) && py::reinterpret_borrow<py::array>(x).ndim() > 0;
        if (array || !not_a_sequence(x)) {
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

  py::class_<Eligibility, std::shared_ptr<Eligibility>>(
      m, "Eligibility",
      "Immutable, graph-bound static eligibility: dense node and edge flags. Built by\n"
      "StrictMultiDiGraph.eligibility(excluded_nodes, excluded_edges).")
      .def("graph_uid", &Eligibility::graph_uid, "uid of the graph this value belongs to")
      .def("uid", &Eligibility::uid, "0 for the all-allowed value, a fresh uid otherwise")
      .def(
          "node_allowed",
          [](const py::object& self) {
            const auto& el = self.cast<const Eligibility&>();
            return readonly_view(el.node_allowed(), self, py::dtype("bool"));
          },
          "bool[N] read-only view")
      .def(
          "edge_allowed",
          [](const py::object& self) {
            const auto& el = self.cast<const Eligibility&>();
            return readonly_view(el.edge_allowed(), self, py::dtype("bool"));
          },
          "bool[E] read-only view (not excluded and both endpoints allowed)");

  py::class_<StrictMultiDiGraph, std::shared_ptr<StrictMultiDiGraph>>(
      m, "StrictMultiDiGraph",
      "Immutable directed multigraph on integer units. Edges are given in input order; the\n"
      "public EdgeId is the rank in the stable order (src, cost, dst, input), input_index()[e]\n"
      "is the input position of e, and every other per-edge array is in id order. Costs are\n"
      "positive integers; capacities are int64 units of the explicit resolution q (a power of\n"
      "two in [2^-32, 2^32]); physical capacities enter through to_units.")
      .def(py::init(&make_graph), py::arg("num_nodes"), py::arg("src"), py::arg("dst"),
           py::arg("capacity_units"), py::arg("cost"), py::arg("resolution"),
           "Validates, in order, and reports the first defect as ValueError: the resolution (a power of\n"
           "two in [2^-32, 2^32]), the unit values (every capacity in [0, 2^62)), the array shapes,\n"
           "zero or negative costs (rule 28), the remaining bounds (endpoints in range, then the\n"
           "capacity and cost sums below 2^62). Integers are read without narrowing, so the error\n"
           "does not depend on how the input is represented (APPROACH 3).")
      .def("uid", &StrictMultiDiGraph::uid, "process-unique id (never reused)")
      .def("resolution", &StrictMultiDiGraph::resolution, "the resolution q (units of capacity)")
      .def("num_nodes", &StrictMultiDiGraph::num_nodes)
      .def("num_edges", &StrictMultiDiGraph::num_edges)
      .def("max_cost", &StrictMultiDiGraph::max_cost, "largest edge cost (0 without edges)")
      .def("cost_gcd", &StrictMultiDiGraph::cost_gcd, "gcd of the edge costs (0 without edges)")
      .def(
          "capacity_units", [](const py::object& self) {
            const auto& g = self.cast<const StrictMultiDiGraph&>();
            return readonly_view(g.capacity_units(), self);
          },
          "int64[E] capacities in units, id order (read-only view)")
      .def(
          "cost", [](const py::object& self) {
            const auto& g = self.cast<const StrictMultiDiGraph&>();
            return readonly_view(g.cost(), self);
          },
          "int64[E] costs, id order (read-only view)")
      .def(
          "edge_src", [](const py::object& self) {
            const auto& g = self.cast<const StrictMultiDiGraph&>();
            return readonly_view(g.edge_src(), self);
          },
          "int32[E] source node per edge (read-only view)")
      .def(
          "edge_dst", [](const py::object& self) {
            const auto& g = self.cast<const StrictMultiDiGraph&>();
            return readonly_view(g.edge_dst(), self);
          },
          "int32[E] destination node per edge (read-only view)")
      .def(
          "input_index", [](const py::object& self) {
            const auto& g = self.cast<const StrictMultiDiGraph&>();
            return readonly_view(g.input_index(), self);
          },
          "int32[E]: id -> input position (read-only view); the inverse is the caller's")
      .def(
          "row_offsets", [](const py::object& self) {
            const auto& g = self.cast<const StrictMultiDiGraph&>();
            return readonly_view(g.row_offsets(), self);
          },
          "int32[N + 1]: the out-edges of u are the ids [row_offsets[u], row_offsets[u + 1])")
      .def(
          "in_row_offsets", [](const py::object& self) {
            const auto& g = self.cast<const StrictMultiDiGraph&>();
            return readonly_view(g.in_row_offsets(), self);
          },
          "int32[N + 1] reverse CSR offsets")
      .def(
          "in_adj_edge_index", [](const py::object& self) {
            const auto& g = self.cast<const StrictMultiDiGraph&>();
            return readonly_view(g.in_adj_edge_index(), self);
          },
          "int32[E]: the in-edges of v in (cost, src, input) order, as EdgeIds")
      .def(
          "eligibility",
          [](const StrictMultiDiGraph& g, const py::object& excluded_nodes, const py::object& excluded_edges) {
            // Nodes, then edges, each range-checked before narrowing (the constructor's order).
            const auto xn = wide_integers(excluded_nodes, "excluded_nodes");
            for (std::size_t i = 0; i < xn.size(); ++i) {
              if (!xn.in(i, 0, static_cast<std::int64_t>(g.num_nodes()) - 1)) {
                throw py::value_error("excluded node " + xn.text(i) + " out of node range [0, " +
                                      std::to_string(g.num_nodes()) + ")");
              }
            }
            const auto xe = wide_integers(excluded_edges, "excluded_edges");
            for (std::size_t i = 0; i < xe.size(); ++i) {
              if (!xe.in(i, 0, static_cast<std::int64_t>(g.num_edges()) - 1)) {
                throw py::value_error("excluded edge " + xe.text(i) + " out of edge range [0, " +
                                      std::to_string(g.num_edges()) + ")");
              }
            }
            return std::make_shared<Eligibility>(g.eligibility(xn.narrow<NodeId>(), xe.narrow<EdgeId>()));
          },
          py::arg("excluded_nodes") = py::tuple(), py::arg("excluded_edges") = py::tuple(),
          "Static eligibility from excluded node ids and EdgeIds (duplicates and unsorted lists\n"
          "accepted, out-of-range ids rejected); an edge is allowed when it is not excluded and\n"
          "both endpoints are allowed. Empty lists mean all allowed (uid 0).");

  m.def(
      "_place_global_units",
      [](const std::shared_ptr<StrictMultiDiGraph>& graph, const py::object& view, const py::object& sources,
         const py::object& destinations, const std::shared_ptr<Eligibility>& elig,
         const py::object& cap) -> py::dict {
        if (!graph) throw py::value_error("graph must be a StrictMultiDiGraph");
        // The core's order (view length, view entries, sources, destinations, overlap, cap,
        // Eligibility), every value checked before narrowing; the core re-checks the same order.
        const auto num_edges = static_cast<std::size_t>(graph->num_edges());
        const std::size_t view_len = sequence_length(view, "view");
        if (view_len != num_edges) {
          throw py::value_error("view has length " + std::to_string(view_len) + ", expected " +
                                std::to_string(num_edges) + " (one entry per edge)");
        }
        const auto v = wide_integers(view, "view (capacity units)");
        const auto capacity = graph->capacity_units();
        for (std::size_t e = 0; e < num_edges; ++e) {
          if (!v.in(e, 0, capacity[e])) {
            throw py::value_error("view[" + std::to_string(e) + "] = " + v.text(e) + " is outside [0, capacity_units[" +
                                  std::to_string(e) + "] = " + std::to_string(capacity[e]) + "]");
          }
        }
        const auto endpoints = [&](const py::object& obj, const char* name) {
          const auto ids = wide_integers(obj, (std::string(name) + " node ids").c_str());
          if (ids.size() == 0) throw py::value_error(std::string(name) + ": the endpoint set is empty");
          for (std::size_t i = 0; i < ids.size(); ++i) {
            if (!ids.in(i, 0, static_cast<std::int64_t>(graph->num_nodes()) - 1)) {
              throw py::value_error(std::string(name) + ": node " + ids.text(i) + " out of node range [0, " +
                                    std::to_string(graph->num_nodes()) + ")");
            }
          }
          return ids.narrow<NodeId>();
        };
        const auto s = endpoints(sources, "sources");
        const auto d = endpoints(destinations, "destinations");
        std::vector<NodeId> ss(s);
        std::vector<NodeId> dd(d);
        std::sort(ss.begin(), ss.end());
        std::sort(dd.begin(), dd.end());
        for (std::size_t i = 0, k = 0; i < ss.size() && k < dd.size();) {
          if (ss[i] == dd[k]) {
            throw py::value_error("node " + std::to_string(ss[i]) +
                                  " is in both endpoint sets: sources and destinations must not overlap");
          }
          if (ss[i] < dd[k]) ++i; else ++k;
        }
        std::optional<Amount> c;
        if (!cap.is_none()) {
          if (is_bool(cap) || !PyIndex_Check(cap.ptr())) {
            throw py::value_error("cap must be an integer number of units in [0, 2^63 - 2] or None");
          }
          py::object idx = python_int(cap, "cap");
          int overflow = 0;
          const long long value = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
          if (overflow != 0 || value < 0 || value > kMaxFlowCap) {
            throw py::value_error("cap " + std::string(py::str(idx)) + " is outside [0, 2^63 - 2]");
          }
          c = static_cast<Amount>(value);
        }
        FlowCoreResult r;
        {
          py::gil_scoped_release release;
          const auto view_units = v.narrow<Amount>();
          r = place_global_units(*graph, view_units, s, d, elig.get(), c);
        }
        py::dict out;
        out["placed"] = py::int_(r.placed);
        py::list curve;
        for (const auto& [cost, units] : r.cost_distribution) {
          py::list entry;
          entry.append(py::int_(cost));
          entry.append(py::int_(units));
          curve.append(entry);
        }
        out["cost_distribution"] = curve;
        out["edge_flows"] = int_list(r.edge_flows);
        if (r.min_cut.has_value()) {
          py::list cut;
          for (EdgeId e : *r.min_cut) cut.append(py::int_(e));
          out["min_cut"] = cut;
        } else {
          out["min_cut"] = py::none();
        }
        out["potentials"] = int_list(r.potentials);
        py::dict counters;
        counters["dijkstra_runs"] = py::int_(r.counters.dijkstra_runs);
        counters["dijkstra_pops"] = py::int_(r.counters.dijkstra_pops);
        counters["settled"] = py::int_(r.counters.settled);
        counters["arc_scans"] = py::int_(r.counters.arc_scans);
        counters["dinic_phases"] = py::int_(r.counters.dinic_phases);
        counters["augmentations"] = py::int_(r.counters.augmentations);
        out["counters"] = counters;
        return out;
      },
      py::arg("graph"), py::arg("view"), py::arg("sources"), py::arg("destinations"), py::arg("elig"),
      py::arg("cap"),
      "Private recorder entry (APPROACH 2.2): the minimum-cost flow of the maximum value under\n"
      "cap (None: uncapped) from the source set to the destination set over the eligible edges,\n"
      "on the capacities `view` (int64 units, id order), from zero flow. Returns {placed,\n"
      "cost_distribution: [[cost, units]], edge_flows, min_cut (None while the destinations stay\n"
      "connected, and on cap 0), potentials (-1 = never reached), counters}. Invalid inputs raise\n"
      "ValueError before any work.");
}
