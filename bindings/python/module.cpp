/*
  Pybind11 module `_netgraph_core`: the graph on integer units, the flow core's private entry,
  the ledger with its records, the route structures (PredDAG, Path), Eligibility with
  intersect, the selection and placement enums, and the private entries of the 1b recorder
  population (a Pinned placement with one bundle, the internal search, the apportionment).

  - Graph accessors return read-only NumPy views that keep the graph alive (no per-call copy).
  - Integer inputs (capacities, costs, ids, views, caps) must be integers: Python ints, NumPy
    integer scalars or integer arrays; floats and bools are rejected with ValueError, never
    truncated. Values are read without narrowing (WideInts) and validated in the contract's
    order before any conversion to the C++ types (APPROACH 3, revision 25). Physical values
    enter only through to_units; the resolution is validated before anything else.
  - std::invalid_argument from Core maps to ValueError.
  - _place_global_units is the recorder's private entry point (APPROACH 2.2): it returns the
    payload as a dict of Python ints and lists.
  - Amounts and counters cross as Python ints (exact, also beyond 64 bits); per-edge arrays as
    int64 NumPy arrays in EdgeId order; +inf is accepted only as a request (float('inf')).
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

#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/pinned_flow.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/search.hpp"
#include "netgraph/core/flow.hpp"
#include "netgraph/core/flow_graph.hpp"
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"
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
  // The entry as an exact Python int (an overflowing entry compared with another).
  [[nodiscard]] py::object exact(std::size_t i) const {
    py::object item = source[py::int_(i)];
    return py::reinterpret_steal<py::object>(PyNumber_Index(item.ptr()));
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


// ---------------------------------------------------------------- 1b helpers

// The exact order of wide integers: a positive overflow above every int64, a negative one
// below, and two overflows on the same side by their exact values (review-1b D2).
bool wide_less(const WideInts& w, std::size_t i, std::size_t j) {
  if (w.sign(i) != w.sign(j)) return w.sign(i) < w.sign(j);
  if (w.sign(i) == 0) return w.value[i] < w.value[j];
  const int less = PyObject_RichCompareBool(w.exact(i).ptr(), w.exact(j).ptr(), Py_LT);
  if (less < 0) throw py::error_already_set();
  return less == 1;
}

// A per-edge integer array in [0, capacity_units[e]] (a ledger baseline, a search residual),
// validated on the unnarrowed values: the length (shape step), then every entry.
std::vector<Amount> edge_units(const py::object& obj, const StrictMultiDiGraph& g, const char* name) {
  const auto num_edges = static_cast<std::size_t>(g.num_edges());
  const std::size_t len = sequence_length(obj, name);
  if (len != num_edges) {
    throw py::value_error(std::string(name) + " has length " + std::to_string(len) + ", expected " +
                          std::to_string(num_edges) + " (one entry per edge)");
  }
  const auto v = wide_integers(obj, (std::string(name) + " (units)").c_str());
  const auto capacity = g.capacity_units();
  for (std::size_t e = 0; e < num_edges; ++e) {
    if (!v.in(e, 0, capacity[e])) {
      throw py::value_error(std::string(name) + "[" + std::to_string(e) + "] = " + v.text(e) + " is outside [0, capacity_units[" +
                            std::to_string(e) + "] = " + std::to_string(capacity[e]) + "]");
    }
  }
  return v.narrow<Amount>();
}

// A node set: one NodeId or a sequence of them, ids range-checked before narrowing.
std::vector<NodeId> node_ids(const py::object& obj, std::int32_t num_nodes, const char* name) {
  WideInts ids;
  if (!is_bool(obj) && PyIndex_Check(obj.ptr()) && !py::isinstance<py::array>(obj)) {
    py::list one;
    one.append(obj);
    ids = wide_integers(one, name);
  } else {
    ids = wide_integers(obj, name);
  }
  if (ids.size() == 0) throw py::value_error(std::string(name) + ": the node set is empty");
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (!ids.in(i, 0, static_cast<std::int64_t>(num_nodes) - 1)) {
      throw py::value_error(std::string(name) + ": node " + ids.text(i) + " out of node range [0, " +
                            std::to_string(num_nodes) + ")");
    }
  }
  return ids.narrow<NodeId>();
}

// A request: an integer number of units in [0, 2^62) or +inf (float('inf')); nullopt is +inf.
std::optional<Amount> request_units(const py::object& volume) {
  if (py::isinstance<py::float_>(volume) || py::isinstance(volume, py::module_::import("numpy").attr("floating"))) {
    const double x = volume.cast<double>();
    if (std::isinf(x) && x > 0) return std::nullopt;
    throw py::value_error("requests are integer units or +inf (a physical request goes through to_units first; got " +
                          std::string(py::str(volume)) + ")");
  }
  const auto v = integer_value(volume, 0, kDomainLimit - 1, "request (units)");
  return static_cast<Amount>(v);
}

py::object counter_int(const Counter128& c) {
  py::int_ hi(c.hi);
  py::int_ lo(c.lo);
  return hi.attr("__lshift__")(py::int_(64)).attr("__or__")(lo);
}

py::object cost_or_none(Cost c) { return c == kUnreachable ? py::object(py::none()) : py::object(py::int_(c)); }

py::array_t<std::int64_t> int64_array(const std::vector<std::int64_t>& v) {
  py::array_t<std::int64_t> out(static_cast<py::ssize_t>(v.size()));
  std::copy(v.begin(), v.end(), out.mutable_data());
  return out;
}

const char* termination_name(detail::Termination t) {
  switch (t) {
    case detail::Termination::Complete: return "complete";
    case detail::Termination::LimitFlows: return "limit_flows";
    case detail::Termination::LimitDiminishing: return "limit_diminishing";
    case detail::Termination::NoRoute: return "no_route";
    case detail::Termination::Exhausted: return "exhausted";
  }
  return "?";
}

FlowId flow_id_of(const py::object& obj) {
  if (!py::isinstance<py::tuple>(obj) || py::len(obj) != 2) {
    throw py::value_error("a FlowId is a tuple (ledger uid, id)");
  }
  auto t = py::reinterpret_borrow<py::tuple>(obj);
  py::object uid = python_int(t[0], "FlowId ledger uid");
  py::object id = python_int(t[1], "FlowId id");
  int o1 = 0, o2 = 0;
  const unsigned long long u = PyLong_AsUnsignedLongLong(uid.ptr());
  if (PyErr_Occurred()) {
    PyErr_Clear();
    o1 = 1;
  }
  const long long i = PyLong_AsLongLongAndOverflow(id.ptr(), &o2);
  if (o1 != 0 || o2 != 0 || i < 0) {
    throw py::value_error("FlowId (" + std::string(py::str(uid)) + ", " + std::string(py::str(id)) + ") is retired or unknown");
  }
  return FlowId {static_cast<std::uint64_t>(u), static_cast<std::int64_t>(i)};
}

py::list pairs_list(const std::vector<std::pair<EdgeId, Amount>>& v) {
  py::list out;
  for (const auto& [e, x] : v) out.append(py::make_tuple(e, x));
  return out;
}

py::list int_list(const std::vector<std::int64_t>& v) {
  py::list out(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) out[i] = py::int_(v[i]);
  return out;
}

} // namespace

PYBIND11_MODULE(_netgraph_core, m, py::mod_gil_not_used()) {
  m.doc() = "NetGraph-Core: the graph on integer units, the minimum-cost flow core, the ledger and the kernels";

  py::enum_<FlowPlacement>(m, "FlowPlacement",
                           "The kernel of one flow: MAX_FLOW (the minimum-cost flow core on the flow's split graph,\n"
                           "a replacement on V = W + own), WEIGHTED (hop-by-hop WCMP: highest averages over frozen\n"
                           "weights at every node), EQUAL_BALANCED (ECMP: round robin at every node).")
      .value("MAX_FLOW", FlowPlacement::MaxFlow)
      .value("WEIGHTED", FlowPlacement::Weighted)
      .value("EQUAL_BALANCED", FlowPlacement::EqualBalanced);

  py::enum_<Admission>(m, "Admission",
                       "LOSSLESS admits only what fits; LOSSY forwards hop by hop without admission control and\n"
                       "drops what does not fit on each edge (EQUAL_BALANCED and WEIGHTED only).")
      .value("LOSSLESS", Admission::Lossless)
      .value("LOSSY", Admission::Lossy);

  py::enum_<EdgeTieBreak>(m, "EdgeTieBreak",
                          "Ties between equal-cost candidates (parallel members and equal-cost walks):\n"
                          "DETERMINISTIC (lowest EdgeId; capacities ignored) or PREFER_HIGHER_PRIORITY (the higher\n"
                          "priority value first, then the lowest EdgeId).")
      .value("DETERMINISTIC", EdgeTieBreak::Deterministic)
      .value("PREFER_HIGHER_PRIORITY", EdgeTieBreak::PreferHigherPriority);

  py::enum_<EdgePriority>(m, "EdgePriority",
                          "The per-edge value PREFER_HIGHER_PRIORITY reads, and the Weighted kernel's membership and\n"
                          "weights: RESIDUAL (residual units), CAPACITY (static capacity units), HEADROOM (the exact\n"
                          "rational residual / capacity; rejected with WEIGHTED).")
      .value("RESIDUAL", EdgePriority::Residual)
      .value("CAPACITY", EdgePriority::Capacity)
      .value("HEADROOM", EdgePriority::Headroom);

  py::class_<PathSelection>(m, "PathSelection",
                            "The selection of a search: tie_break and priority are explicit; multipath = True keeps\n"
                            "every equal-cost walk (the DAG), False the one canonical single path; multi_edge = True\n"
                            "keeps every eligible parallel member, False one per (u, v) pair.")
      .def(py::init<EdgeTieBreak, EdgePriority, bool, bool>(), py::arg("tie_break"), py::arg("priority"),
           py::arg("multipath") = true, py::arg("multi_edge") = true)
      .def_readonly("tie_break", &PathSelection::tie_break)
      .def_readonly("priority", &PathSelection::priority)
      .def_readonly("multipath", &PathSelection::multipath)
      .def_readonly("multi_edge", &PathSelection::multi_edge)
      .def("__eq__", [](const PathSelection& a, const py::object& b) {
        return py::isinstance<PathSelection>(b) && a == b.cast<const PathSelection&>();
      })
      .def("__hash__", [](const PathSelection& a) {
        return py::hash(py::make_tuple(static_cast<int>(a.tie_break), static_cast<int>(a.priority), a.multipath,
                                       a.multi_edge));
      })
      .def("__repr__", [](const PathSelection& a) {
        return std::string("PathSelection(tie_break=") +
               (a.tie_break == EdgeTieBreak::Deterministic ? "DETERMINISTIC" : "PREFER_HIGHER_PRIORITY") +
               ", priority=" +
               (a.priority == EdgePriority::Residual ? "RESIDUAL"
                                                     : (a.priority == EdgePriority::Capacity ? "CAPACITY" : "HEADROOM")) +
               ", multipath=" + (a.multipath ? "True" : "False") + ", multi_edge=" + (a.multi_edge ? "True" : "False") +
               ")";
      });

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
          "bool[E] read-only view (not excluded and both endpoints allowed)")
      .def_static(
          "intersect",
          [](const std::shared_ptr<Eligibility>& a, const std::shared_ptr<Eligibility>& b) -> std::shared_ptr<Eligibility> {
            if (!a || !b) throw py::value_error("Eligibility.intersect takes two Eligibility values");
            auto r = Eligibility::intersect(*a, *b);
            if (r.uid() == a->uid()) return a;
            if (r.uid() == b->uid()) return b;
            return std::make_shared<Eligibility>(std::move(r));
          },
          py::arg("a"), py::arg("b"),
          "The intersection of two values of one graph (ValueError for values of different graphs):\n"
          "an operand itself when the result equals it (a first, then b), otherwise a new value with a\n"
          "fresh uid; equal intersections built separately get distinct uids.");

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

  // ------------------------------------------------------------------ route structures

  py::class_<PredDAG, std::shared_ptr<PredDAG>>(
      m, "PredDAG",
      "Immutable, graph-bound predecessor DAG: the entries of node v are parents()[k], via_edges()[k]\n"
      "for k in [parent_offsets()[v], parent_offsets()[v + 1]), within a node in strictly ascending\n"
      "via EdgeId (one canonical form). Equality is structural. Built by PredDAG.from_csr,\n"
      "Path.to_dag() and the searches.")
      .def_static(
          "from_csr",
          [](const std::shared_ptr<StrictMultiDiGraph>& graph, const py::object& offsets, const py::object& parents,
             const py::object& via) {
            if (!graph) throw py::value_error("graph must be a StrictMultiDiGraph");
            // APPROACH 5's order on unnarrowed values; the C++ constructor re-checks it and adds the
            // structural checks (endpoints, duplicates, acyclic, ascending).
            const auto n = static_cast<std::size_t>(graph->num_nodes());
            const auto m = static_cast<std::int64_t>(graph->num_edges());
            const std::size_t lo = sequence_length(offsets, "offsets");
            const std::size_t lp = sequence_length(parents, "parents");
            const std::size_t lv = sequence_length(via, "via");
            if (lo != n + 1) {
              throw py::value_error("PredDAG.from_csr: offsets has length " + std::to_string(lo) +
                                    ", expected N + 1 = " + std::to_string(n + 1));
            }
            const auto o = wide_integers(offsets, "offsets");
            if (!o.in(0, 0, 0)) throw py::value_error("PredDAG.from_csr: offsets[0] = " + o.text(0) + ", expected 0");
            for (std::size_t v = 0; v < n; ++v) {
              if (wide_less(o, v + 1, v)) {
                throw py::value_error("PredDAG.from_csr: offsets are not monotone at " + std::to_string(v + 1) + " (" +
                                      o.text(v) + " > " + o.text(v + 1) + ")");
              }
            }
            if (!o.in(n, static_cast<std::int64_t>(lp), static_cast<std::int64_t>(lp)) || lp != lv) {
              throw py::value_error("PredDAG.from_csr: offsets[N] = " + o.text(n) + " but parents has " +
                                    std::to_string(lp) + " and via has " + std::to_string(lv) + " entries");
            }
            const auto p = wide_integers(parents, "parents");
            const auto w = wide_integers(via, "via");
            for (std::size_t k = 0; k < lp; ++k) {
              if (!p.in(k, 0, static_cast<std::int64_t>(n) - 1)) {
                throw py::value_error("PredDAG.from_csr: parents[" + std::to_string(k) + "] = " + p.text(k) +
                                      " out of node range [0, " + std::to_string(n) + ")");
              }
              if (!w.in(k, 0, m - 1)) {
                throw py::value_error("PredDAG.from_csr: via[" + std::to_string(k) + "] = " + w.text(k) +
                                      " out of edge range [0, " + std::to_string(m) + ")");
              }
            }
            const auto ov = o.narrow<std::int32_t>();
            const auto pv = p.narrow<NodeId>();
            const auto wv = w.narrow<EdgeId>();
            return std::make_shared<PredDAG>(PredDAG::from_csr(*graph, ov, pv, wv));
          },
          py::arg("graph"), py::arg("offsets"), py::arg("parents"), py::arg("via"),
          "Validates, in order, before any traversal (ValueError on the first defect): offsets of\n"
          "length N + 1, offsets[0] == 0, monotone, offsets[N] == len(parents) == len(via); then every\n"
          "parent and via in range, every via joining (parent, child), no duplicate (child, via),\n"
          "acyclic; within a node the entries in ascending via EdgeId.")
      .def("graph_uid", &PredDAG::graph_uid, "uid of the graph this DAG belongs to")
      .def(
          "parent_offsets",
          [](const py::object& self) { return readonly_view(self.cast<const PredDAG&>().parent_offsets(), self); },
          "int32[N + 1] CSR offsets (read-only view)")
      .def(
          "parents", [](const py::object& self) { return readonly_view(self.cast<const PredDAG&>().parents(), self); },
          "int32 parent NodeId per entry (read-only view)")
      .def(
          "via_edges",
          [](const py::object& self) { return readonly_view(self.cast<const PredDAG&>().via_edges(), self); },
          "int32 via EdgeId per entry (read-only view)")
      .def("__eq__", [](const PredDAG& a, const py::object& b) {
        return py::isinstance<PredDAG>(b) && a == b.cast<const PredDAG&>();
      })
      .def("__hash__", [](const PredDAG& a) {
        return py::hash(py::make_tuple(a.graph_uid(), py::tuple(py::cast(std::vector<EdgeId>(a.via_edges().begin(),
                                                                                              a.via_edges().end())))));
      })
      .def("__repr__", [](const PredDAG& a) {
        return "PredDAG(graph_uid=" + std::to_string(a.graph_uid()) + ", entries=" +
               std::to_string(a.via_edges().size()) + ")";
      });

  py::class_<Path, std::shared_ptr<Path>>(
      m, "Path",
      "Immutable, graph-bound explicit walk of hops; a hop lists one or more parallel members of one\n"
      "(u, v) pair in the caller's order. cost = sum over hops of the minimum member cost.")
      .def(py::init([](const std::shared_ptr<StrictMultiDiGraph>& graph, const py::object& hops) {
             if (!graph) throw py::value_error("graph must be a StrictMultiDiGraph");
             if (not_a_sequence(hops)) throw py::value_error("Path: hops must be a sequence of hops (sequences of EdgeIds)");
             auto seq = py::reinterpret_borrow<py::sequence>(hops);
             std::vector<WideInts> raw;
             for (std::size_t h = 0; h < seq.size(); ++h) {
               py::object hop = seq[h];
               raw.push_back(wide_integers(hop, ("hops[" + std::to_string(h) + "]").c_str()));
             }
             // Every member in range before anything else (APPROACH 5), then C++'s order.
             const auto m = static_cast<std::int64_t>(graph->num_edges());
             for (std::size_t h = 0; h < raw.size(); ++h) {
               for (std::size_t i = 0; i < raw[h].size(); ++i) {
                 if (!raw[h].in(i, 0, m - 1)) {
                   throw py::value_error("Path: hop " + std::to_string(h) + " member " + raw[h].text(i) +
                                         " out of edge range [0, " + std::to_string(m) + ")");
                 }
               }
             }
             std::vector<std::vector<EdgeId>> narrowed;
             narrowed.reserve(raw.size());
             for (const auto& r : raw) narrowed.push_back(r.narrow<EdgeId>());
             return std::make_shared<Path>(*graph, std::move(narrowed));
           }),
           py::arg("graph"), py::arg("hops"),
           "Validates, in order (ValueError on the first defect): every member id in range, no duplicate\n"
           "member, a non-empty hop list, every hop non-empty, the members of a hop sharing (u, v),\n"
           "contiguous hops, no repeated node.")
      .def("graph_uid", &Path::graph_uid, "uid of the graph this path belongs to")
      .def("hops", [](const Path& p) { return p.hops(); }, "the hops as lists of EdgeIds (the caller's member order)")
      .def("src", &Path::src)
      .def("dst", &Path::dst)
      .def("nodes", [](const Path& p) { return p.nodes(); }, "the visited nodes, src first")
      .def("cost", &Path::cost, "sum over hops of the minimum member cost")
      .def(
          "to_dag", [](const Path& p) { return std::make_shared<PredDAG>(p.to_dag()); },
          "the equivalent PredDAG (every member of every hop)")
      .def("__eq__", [](const Path& a, const py::object& b) {
        return py::isinstance<Path>(b) && a == b.cast<const Path&>();
      })
      .def("__hash__",
           [](const Path& a) {
             py::list hops;
             for (const auto& hop : a.hops()) hops.append(py::tuple(py::cast(hop)));
             return py::hash(py::make_tuple(a.graph_uid(), py::tuple(hops)));
           })
      .def("__repr__", [](const Path& a) {
        return "Path(src=" + std::to_string(a.src()) + ", dst=" + std::to_string(a.dst()) +
               ", hops=" + std::to_string(a.hops().size()) + ", cost=" + std::to_string(a.cost()) + ")";
      });

  // ------------------------------------------------------------------ the ledger

  py::class_<FlowRecord>(
      m, "FlowRecord",
      "One flow of a ledger, as a snapshot: flow_id (ledger uid, id), owner (the policy's id),\n"
      "sources, destinations, placement, admission, bundle (PredDAG: the flow's route),\n"
      "route_cost (its minimum walk cost; None when empty), offered (Lossy only, else None),\n"
      "delivered, edges [(EdgeId, units)] and drops [(EdgeId, units)] (cumulative, ascending\n"
      "EdgeId). Counters are exact Python ints (they may exceed 64 bits).")
      .def_property_readonly("flow_id", [](const FlowRecord& r) { return py::make_tuple(r.flow_id.first, r.flow_id.second); })
      .def_readonly("owner", &FlowRecord::owner)
      .def_readonly("sources", &FlowRecord::sources)
      .def_readonly("destinations", &FlowRecord::destinations)
      .def_readonly("placement", &FlowRecord::placement)
      .def_readonly("admission", &FlowRecord::admission)
      .def_property_readonly("bundle", [](const FlowRecord& r) { return std::const_pointer_cast<PredDAG>(r.bundle); })
      .def_property_readonly("route_cost", [](const FlowRecord& r) { return cost_or_none(r.route_cost); })
      .def_property_readonly("offered", [](const FlowRecord& r) -> py::object {
        return r.offered.has_value() ? counter_int(*r.offered) : py::object(py::none());
      })
      .def_readonly("delivered", &FlowRecord::delivered)
      .def_property_readonly("edges", [](const FlowRecord& r) { return pairs_list(r.edges); })
      .def_property_readonly("drops", [](const FlowRecord& r) {
        py::list out;
        for (const auto& [e, x] : r.drops) out.append(py::make_tuple(e, counter_int(x)));
        return out;
      })
      .def("__repr__", [](const FlowRecord& r) {
        return "FlowRecord(flow_id=(" + std::to_string(r.flow_id.first) + ", " + std::to_string(r.flow_id.second) +
               "), owner=" + std::to_string(r.owner) + ", delivered=" + std::to_string(r.delivered) + ")";
      });

  py::class_<FlowGraph, std::shared_ptr<FlowGraph>>(
      m, "FlowGraph",
      "The integer ledger of one graph: an immutable baseline residual in units and the records of\n"
      "the flows placed on it by policies (the ledger has no public placement or removal).\n"
      "allocated() is the sum of every record's edges and residual() = baseline - allocated(), both\n"
      "int64[E] snapshots. reset() and copy() are rejected while a policy is bound.")
      .def(py::init([](const std::shared_ptr<StrictMultiDiGraph>& graph, const py::object& residual) {
             if (!graph) throw py::value_error("graph must be a StrictMultiDiGraph");
             if (residual.is_none()) return std::make_shared<FlowGraph>(graph);
             const auto r = edge_units(residual, *graph, "residual");
             return std::make_shared<FlowGraph>(graph, r);
           }),
           py::arg("graph"), py::arg("residual") = py::none(),
           "residual: the baseline as integer units, one entry per edge in [0, capacity_units[e]] (a\n"
           "physical baseline goes through to_units first); None: the capacities. ValueError otherwise.")
      .def("uid", &FlowGraph::uid, "process-unique ledger uid (the first half of every FlowId)")
      .def(
          "graph", [](const FlowGraph& fg) { return std::const_pointer_cast<StrictMultiDiGraph>(fg.graph()); },
          "the ledger's graph")
      .def(
          "capacity_units",
          [](const py::object& self) { return readonly_view(self.cast<const FlowGraph&>().capacity_units(), self); },
          "int64[E] the graph's capacities in units (read-only view)")
      .def(
          "residual", [](const FlowGraph& fg) { return int64_array(fg.residual()); },
          "int64[E] snapshot: baseline - allocated")
      .def(
          "allocated", [](const FlowGraph& fg) { return int64_array(fg.allocated()); },
          "int64[E] snapshot: the sum of every record's edges")
      .def("records", &FlowGraph::records, "live records (bound and orphaned) as snapshots in FlowId order")
      .def(
          "record", [](const FlowGraph& fg, const py::object& flow_id) { return fg.record(flow_id_of(flow_id)); },
          py::arg("flow_id"),
          "one record by FlowId (ValueError for a FlowId of another ledger, a retired or an unknown id)")
      .def(
          "bottleneck",
          [](const FlowGraph& fg, const py::object& sources, const std::shared_ptr<Eligibility>& elig) {
            const auto s = node_ids(sources, fg.graph()->num_nodes(), "sources");
            return fg.bottleneck(s, elig.get());
          },
          py::arg("sources"), py::arg("elig") = nullptr,
          "Saturation diagnostic: the edges (statically eligible, capacity >= 1 unit) leaving the set\n"
          "reachable from the allowed sources over eligible forward arcs with residual >= 1 unit, in\n"
          "ascending EdgeId (not a cut).")
      .def(
          "reset",
          [](FlowGraph& fg, const py::object& residual) {
            if (residual.is_none()) {
              fg.reset();
              return;
            }
            const auto r = edge_units(residual, *fg.graph(), "residual");
            fg.reset(r);
          },
          py::arg("residual") = py::none(),
          "A new baseline (None: the capacities); orphaned records and their allocations are cleared\n"
          "and their ids retired (ids are never reused). ValueError while a policy is bound.")
      .def(
          "copy", [](const FlowGraph& fg) { return std::make_shared<FlowGraph>(fg); },
          "a new ledger (new uid) with the same baseline and records, ids kept under the new uid;\n"
          "ValueError while a policy is bound")
      .def("__copy__", [](const FlowGraph& fg) { return std::make_shared<FlowGraph>(fg); })
      .def("__deepcopy__", [](const FlowGraph& fg, const py::object&) { return std::make_shared<FlowGraph>(fg); },
           py::arg("memo"))
      .def(
          "_state",
          [](const FlowGraph& fg) {
            py::dict out;
            out["next_flow_id"] = detail::LedgerAccess::next_flow_id(fg);
            out["bound"] = detail::LedgerAccess::bound_owners(fg);
            out["orphans"] = detail::LedgerAccess::orphan_owners(fg);
            return out;
          },
          "private diagnostics: {next_flow_id, bound, orphans} (owner ids)");

  // ------------------------------------------------------------------ private 1b recorder entries

  py::class_<detail::PinnedFlow>(
      m, "_PinnedFlow",
      "Private (sub-PR 1b): a Pinned placement with one bundle, the single-flow transaction of the\n"
      "1b recorder population (APPROACH 8: a manual flow is a Pinned policy with one bundle). Sub-PR\n"
      "1c's FlowPolicy with Pinned routing replaces it.")
      .def(py::init([](FlowGraph& fg, const py::object& src, const py::object& dst, const std::shared_ptr<PredDAG>& bundle,
                       FlowPlacement placement, Admission admission, EdgePriority priority,
                       const std::shared_ptr<Eligibility>& elig) {
             const auto n = fg.graph()->num_nodes();
             const auto s = integer_value(src, 0, static_cast<std::int64_t>(n) - 1, "src");
             const auto d = integer_value(dst, 0, static_cast<std::int64_t>(n) - 1, "dst");
             if (!bundle) throw py::value_error("bundle must be a PredDAG");
             return std::make_unique<detail::PinnedFlow>(fg, static_cast<NodeId>(s), static_cast<NodeId>(d), bundle,
                                                         placement, admission, priority, elig);
           }),
           py::keep_alive<1, 2>(), py::arg("flow_graph"), py::arg("src"), py::arg("dst"), py::arg("bundle"),
           py::arg("placement"), py::arg("admission"), py::arg("priority"), py::arg("elig") = nullptr)
      .def(
          "place",
          [](detail::PinnedFlow& pf, const py::object& volume) {
            const auto v = request_units(volume);
            // The GIL stays held: the call mutates the ledger, which other Python threads may
            // read (one writer per ledger, APPROACH 13).
            const detail::PinnedResult r = pf.place(v);
            py::dict out;
            out["placed"] = py::int_(r.placed);
            out["remaining"] = r.remaining.has_value() ? py::object(py::int_(*r.remaining))
                                                       : py::object(py::float_(std::numeric_limits<double>::infinity()));
            out["termination"] = termination_name(r.termination);
            out["routes_changed"] = r.routes_changed;
            out["weights_changed"] = r.weights_changed;
            out["drops"] = pairs_list(r.drops);
            out["iterations"] = r.iterations;
            return out;
          },
          py::arg("volume"),
          "place(volume: int units or float('inf')) -> {placed, remaining, termination, routes_changed,\n"
          "weights_changed, drops, iterations (diagnostic: bisection steps of the T_max search)}")
      .def("remove", &detail::PinnedFlow::remove)
      .def("release", &detail::PinnedFlow::release, "destroy without remove(): the record stays as an orphan")
      .def("flow_id",
           [](const detail::PinnedFlow& pf) -> py::object {
             const auto id = pf.flow_id();
             return id.has_value() ? py::object(py::make_tuple(id->first, id->second)) : py::object(py::none());
           })
      .def("owner", &detail::PinnedFlow::owner)
      .def("requested", [](const detail::PinnedFlow& pf) { return counter_int(pf.requested()); },
           "the epoch's cumulative finite request (units)")
      .def("infinite", &detail::PinnedFlow::infinite, "whether +inf was requested in this epoch")
      .def("placed_demand", &detail::PinnedFlow::placed_demand)
      .def("usable_bundles", &detail::PinnedFlow::usable_bundles)
      .def("plan", [](const detail::PinnedFlow& pf) -> py::object {
        const auto* plan = pf.plan();
        if (plan == nullptr) return py::none();
        py::dict out;
        out["placement"] = plan->placement;
        out["admission"] = plan->admission;
        out["route"] = plan->route;
        out["members"] = plan->members;
        out["weights"] = plan->placement == FlowPlacement::Weighted ? py::object(py::cast(plan->weights))
                                                                     : py::object(py::none());
        return out;
      });

  m.def(
      "_spf",
      [](const std::shared_ptr<StrictMultiDiGraph>& graph, const py::object& roots, const std::string& direction,
         const py::object& targets, const PathSelection& selection, const py::object& residual,
         const py::object& min_residual, const std::shared_ptr<Eligibility>& elig) {
        if (!graph) throw py::value_error("graph must be a StrictMultiDiGraph");
        if (direction != "forward" && direction != "reverse") {
          throw py::value_error("direction must be 'forward' or 'reverse'");
        }
        const auto r = node_ids(roots, graph->num_nodes(), "roots");
        std::optional<std::vector<NodeId>> t;
        if (!targets.is_none()) t = node_ids(targets, graph->num_nodes(), "targets");
        // Section 6's order: the endpoint checks (overlap, a single path without targets) come
        // before the residual is read (review-1b R15); detail::spf repeats them.
        if (t.has_value()) {
          detail::check_disjoint(detail::node_set(r, graph->num_nodes(), "roots"),
                                 detail::node_set(*t, graph->num_nodes(), "targets"));
        } else if (!selection.multipath) {
          throw py::value_error("multipath = false requires targets (the single path joins a root to a target)");
        }
        std::optional<std::vector<Amount>> res;
        if (!residual.is_none()) res = edge_units(residual, *graph, "residual");
        const auto minr = integer_value(min_residual, 1, kDomainLimit - 1, "min_residual");
        detail::SpfOut o;
        {
          py::gil_scoped_release release;
          std::optional<std::span<const NodeId>> ts;
          if (t.has_value()) ts = std::span<const NodeId>(*t);
          std::optional<std::span<const Amount>> rs;
          if (res.has_value()) rs = std::span<const Amount>(*res);
          o = detail::spf(*graph, r, direction == "forward" ? detail::Direction::Forward : detail::Direction::Reverse, ts,
                          selection, rs, minr, elig.get());
        }
        py::dict out;
        if (o.costs.has_value()) {
          py::list costs;
          for (Cost c : *o.costs) costs.append(cost_or_none(c));
          out["costs"] = costs;
        }
        out["dag"] = std::make_shared<PredDAG>(detail::PredDAGBuilder::from_edges(*graph, o.edges));
        if (t.has_value()) {
          out["target_cost"] = o.target_cost.has_value() ? py::object(py::int_(*o.target_cost)) : py::object(py::none());
        }
        if (o.path.has_value()) out["path"] = o.path->walk;
        return out;
      },
      py::arg("graph"), py::arg("roots"), py::arg("direction"), py::arg("targets"), py::arg("selection"),
      py::arg("residual") = py::none(), py::arg("min_residual") = 1, py::arg("elig") = nullptr,
      "Private (sub-PR 1b; public spf in 1d): the search of APPROACH 6. Returns {costs (full search),\n"
      "dag, target_cost (with targets; None when unreachable), path (the single path's walk)}.");

  m.def(
      "_round_robin",
      [](const py::object& total, const py::object& k) {
        const auto t = integer_value(total, 0, kDomainLimit - 1, "total");
        const auto n = integer_value(k, 1, 1 << 20, "k");
        std::vector<Amount> out(static_cast<std::size_t>(n));
        detail::round_robin(t, out);
        return out;
      },
      py::arg("total"), py::arg("k"), "Private: round robin of total units over k members (APPROACH 2.1).");

  m.def(
      "_highest_averages",
      [](const py::object& total, const py::object& weights) {
        const auto t = integer_value(total, 0, kDomainLimit - 1, "total");
        const auto w = wide_integers(weights, "weights");
        if (w.size() == 0) throw py::value_error("highest averages over no member");
        std::uint64_t sum = 0;
        for (std::size_t i = 0; i < w.size(); ++i) {
          if (!w.in(i, 1, kDomainLimit - 1)) {
            throw py::value_error("weights[" + std::to_string(i) + "] = " + w.text(i) + " is outside [1, 2^62)");
          }
          sum += static_cast<std::uint64_t>(w.value[i]);
          if (sum >= static_cast<std::uint64_t>(kDomainLimit)) throw py::value_error("the weights sum to 2^62 or more");
        }
        const auto wv = w.narrow<Amount>();
        std::vector<Amount> out(wv.size());
        detail::highest_averages(t, wv, out);
        return out;
      },
      py::arg("total"), py::arg("weights"),
      "Private: highest-averages (D'Hondt) apportionment of total units over positive weights, ties by\n"
      "the lower index (APPROACH 2.1).");
}
