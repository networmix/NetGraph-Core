/*
  Pybind11 module `_netgraph_core`: the graph on integer units with Eligibility and the unit
  conversions, the route structures (PredDAG, Path), the selection and placement axes, the
  ledger with its records, the policy engine with its configuration and results, the DagCache
  , the evaluators (capacity, batch_capacity, sensitivity) and the public search (spf); and one
  private entry: the flow core's _place_global_units (the recorder's entry, APPROACH 2.2).
  A build with the private hooks (CMake NETGRAPH_CORE_PRIVATE_HOOKS: the editable development
  build) adds the test controls and the parity diagnostics.

  - Graph accessors return read-only NumPy views that keep the graph alive (no per-call copy).
  - Integer inputs (capacities, costs, ids, views, caps) must be integers: Python ints, NumPy
    integer scalars or integer arrays; floats and bools are rejected with ValueError, never
    truncated. Values are read without narrowing (WideInts) and validated in the contract's
    order before any conversion to the C++ types (APPROACH 3, revision 25). Physical values
    enter only through to_units; the resolution is validated before anything else.
  - std::invalid_argument from Core maps to ValueError; enumeration values outside their
    enumeration are rejected at construction (rule 13).
  - _place_global_units is the recorder's private entry point (APPROACH 2.2): it returns the
    payload as a dict of Python ints and lists.
  - Amounts and counters cross as Python ints (exact, also beyond 64 bits); per-edge arrays as
    int64 NumPy arrays in EdgeId order; +inf is accepted only as a request (float('inf')).
  - Python spells constants in UPPER_SNAKE_CASE: MAX_FLOWS_PER_CALL is the C++ kMaxFlowsPerCall.
*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/evaluator_checks.hpp"
#include "netgraph/core/detail/flow_core.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/detail/traversal.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/search.hpp"
#include "netgraph/core/types.hpp"
#if NETGRAPH_CORE_PRIVATE_HOOKS
#include <atomic>

#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/trace.hpp"
#endif

namespace py = pybind11;
using namespace netgraph::core;

namespace {

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

std::shared_ptr<Graph> make_graph(const py::object& num_nodes, const py::object& src,
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
  return std::make_shared<Graph>(static_cast<std::int32_t>(n), src_v, dst_v, cap_v, cost_v, resolution);
}


// ---------------------------------------------------------------- ledger, records, requests

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
std::vector<Amount> edge_units(const py::object& obj, const Graph& g, const char* name) {
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

// A node set's ids, non-empty and range-checked before narrowing (detail::node_set's messages).
std::vector<NodeId> checked_ids(const WideInts& ids, std::int32_t num_nodes, const char* name) {
  if (ids.size() == 0) throw py::value_error(detail::empty_set_error(name));
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (!ids.in(i, 0, static_cast<std::int64_t>(num_nodes) - 1)) {
      throw py::value_error(detail::node_range_error(name, ids.text(i), num_nodes));
    }
  }
  return ids.narrow<NodeId>();
}

// A node set: one NodeId or a sequence of them.
std::vector<NodeId> node_ids(const py::object& obj, std::int32_t num_nodes, const char* name) {
  if (!is_bool(obj) && PyIndex_Check(obj.ptr()) && !py::isinstance<py::array>(obj)) {
    py::list one;
    one.append(obj);
    return checked_ids(wide_integers(one, name), num_nodes, name);
  }
  return checked_ids(wide_integers(obj, name), num_nodes, name);
}

// One NodeId of [0, num_nodes) (a single endpoint); the range message of detail/endpoints.hpp.
NodeId endpoint(const py::object& obj, std::int32_t num_nodes, const char* name) {
  const py::object idx = python_int(obj, name);
  int overflow = 0;
  const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
  if (overflow != 0 || v < 0 || v >= num_nodes) {
    throw py::value_error(detail::node_range_error(name, std::string(py::str(idx)), num_nodes));
  }
  return static_cast<NodeId>(v);
}

// A count of at least 1 (KSP's k, the resolver's max_paths): the rule-13 message below 1; a count
// above INT64_MAX means the same as INT64_MAX (no list reaches it).
std::int64_t count_of(const py::object& obj, const char* name) {
  const py::object idx = python_int(obj, name);
  if (PyObject_RichCompareBool(idx.ptr(), py::int_(1).ptr(), Py_LT) == 1) {
    throw py::value_error(std::string(name) + " must be >= 1 (got " + std::string(py::str(idx)) + ")");
  }
  int overflow = 0;
  const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
  return overflow != 0 ? std::numeric_limits<std::int64_t>::max() : static_cast<std::int64_t>(v);
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


// The Python form of CapacityOptions and SensitivityOptions: the Eligibility and the cache are
// held by reference count for the call (the C++ structs borrow raw pointers).
struct PyCapacityOptions {
  FlowPolicyConfig config;
  std::shared_ptr<Eligibility> elig;
  std::shared_ptr<DagCache> cache;
  OutputOptions outputs;
};
struct PySensitivityOptions {
  FlowPolicyConfig config;
  std::shared_ptr<Eligibility> elig;
  std::shared_ptr<DagCache> cache;
  OutputOptions outputs;
};

// The Python form of SpfOptions: the arguments as given (targets, residual and min_residual are
// validated at the call, in the contract's order, on their unnarrowed values); the Eligibility is
// held by reference count.
struct PySpfOptions {
  Direction direction;
  PathSelection selection;
  py::object targets;
  py::object residual;
  py::object min_residual;
  std::shared_ptr<Eligibility> elig;
};

// The Python form of KspOptions: the arguments as given (k, max_path_cost and residual are
// validated at the call, in the contract's order, on their unnarrowed values).
struct PyKspOptions {
  py::object k;
  PathSelection selection;
  py::object max_path_cost;
  py::object residual;
  std::shared_ptr<Eligibility> elig;
};

// The Python form of SpfResult: every field already a Python object (made by the call through
// the checked allocations), so reading one allocates nothing.
struct PySpfResult {
  py::object costs;        // a read-only int64 array (full searches) or None
  py::object dag;          // PredDAG
  py::object target_cost;  // an int, or None (no target reached, or a full search)
};

// ---------------------------------------------------------------- policies and evaluators

// An optional int field of [lo, hi] (None allowed when `optional`); bools and floats rejected.
std::optional<std::int64_t> int_field(const py::object& obj, std::int64_t lo, std::int64_t hi, const std::string& what,
                                      bool optional) {
  if (obj.is_none()) {
    if (optional) return std::nullopt;
    throw py::value_error(what + " must be an integer (got None)");
  }
  return integer_value(obj, lo, hi, what);
}

double float_field(const py::object& obj, const std::string& what) {
  if (is_bool(obj) || !(py::isinstance<py::float_>(obj) || PyIndex_Check(obj.ptr()) ||
                        py::isinstance(obj, py::module_::import("numpy").attr("floating")))) {
    throw py::value_error(what + " must be a number (got " + std::string(py::str(py::type::of(obj).attr("__name__"))) +
                          ")");
  }
  return obj.cast<double>();
}

py::object remaining_obj(const std::optional<Amount>& r) {
  return r.has_value() ? py::object(py::int_(*r)) : py::object(py::float_(std::numeric_limits<double>::infinity()));
}

Routing routing_of(const py::object& obj) {
  if (py::isinstance<CostOnly>(obj)) return obj.cast<CostOnly>();
  if (py::isinstance<CapacityAware>(obj)) return obj.cast<CapacityAware>();
  if (py::isinstance<Global>(obj)) return obj.cast<Global>();
  if (py::isinstance<Pinned>(obj)) return obj.cast<Pinned>();
  throw py::type_error("routing must be CostOnly, CapacityAware, Global or Pinned (got " +
                       std::string(py::str(py::type::of(obj).attr("__name__"))) + ")");
}

py::object routing_obj(const Routing& r) {
  return std::visit([](const auto& x) { return py::cast(x); }, r);
}

// workers: an explicit integer is passed through (validated by the evaluator); None takes the
// environment default (detail::default_workers), and 1 whenever a cache is supplied (the
// building cache's one writer; APPROACH 10, 11).
std::int32_t workers_of(const py::object& workers, bool sensitivity, bool cache_given) {
  if (!workers.is_none()) {
    return static_cast<std::int32_t>(integer_value(workers, std::numeric_limits<std::int32_t>::min(),
                                                   std::numeric_limits<std::int32_t>::max(), "workers"));
  }
  if (cache_given) return 1;
  return detail::default_workers(sensitivity);
}

OutputOptions outputs_of(const py::object& obj) {
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

// A new Python object of the bound class T owning a default-constructed T, made as pybind11's
// generic move-cast makes one (type_caster_generic::cast) but with the object allocation
// checked: pybind11's make_new_instance uses the result of tp_alloc unchecked, so a failed
// Python allocation there is a crash rather than a MemoryError. Every failure here raises, and
// the half-built object is released on unwinding (its value is deleted, nothing registered).
// The policy calls make their result object with this before anything is published.
template <class T>
std::pair<py::object, T*> new_object() {
  const py::detail::type_info* tinfo = py::detail::get_type_info(typeid(T), /*throw_if_missing=*/true);
  PyTypeObject* type = tinfo->type;
  PyObject* raw = type->tp_alloc(type, 0);
  if (raw == nullptr) throw py::error_already_set();
  auto obj = py::reinterpret_steal<py::object>(raw);
  auto* inst = reinterpret_cast<py::detail::instance*>(raw);
  inst->allocate_layout();
  inst->owned = false;
  auto* value = new T();
  py::detail::values_and_holders(inst).begin()->value_ptr() = value;
  inst->owned = true;
  tinfo->init_instance(inst, nullptr);
  return {std::move(obj), value};
}

// The value moved into its result object after the commit: nothing may allocate there.
static_assert(std::is_nothrow_move_assignable_v<PlacementResult>);
static_assert(std::is_nothrow_move_assignable_v<SensitivityResult>);
static_assert(std::is_nothrow_move_assignable_v<py::object>);

template <class T>
py::object object_of(T value) {
  auto [obj, slot] = new_object<T>();
  *slot = std::move(value);
  return std::move(obj);
}

// A new Python object of the bound class T (held by std::shared_ptr) owning `value`, through the
// same checked allocation as new_object: for the values without a default constructor (PredDAG,
// Path) that a search returns. A failure raises and releases everything made so far.
template <class T>
py::object object_from(T value) {
  const py::detail::type_info* tinfo = py::detail::get_type_info(typeid(T), /*throw_if_missing=*/true);
  PyTypeObject* type = tinfo->type;
  PyObject* raw = type->tp_alloc(type, 0);
  if (raw == nullptr) throw py::error_already_set();
  auto obj = py::reinterpret_steal<py::object>(raw);
  auto* inst = reinterpret_cast<py::detail::instance*>(raw);
  inst->allocate_layout();
  inst->owned = false;
  auto* held = new T(std::move(value));
  py::detail::values_and_holders(inst).begin()->value_ptr() = held;
  inst->owned = true;
  tinfo->init_instance(inst, nullptr);
  return obj;
}


// ---------------------------------------------------------------- private hooks

#if NETGRAPH_CORE_PRIVATE_HOOKS

// A trace value as Python objects (map keys that are decimal integers become ints: the maps
// keyed by source NodeId).
py::object trace_obj(const detail::Trace& t) {
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

// The epoch state of a policy (roster, requested, +inf, per roster source the flows).
py::dict state_obj(const FlowPolicy& policy) {
  const auto& pe = detail::PolicyAccess::engine(policy);
  py::dict out;
  const auto& r = pe.config.routing;
  out["mode"] = std::holds_alternative<CostOnly>(r)        ? "CostOnly"
                : std::holds_alternative<CapacityAware>(r) ? "CapacityAware"
                : std::holds_alternative<Global>(r)        ? "Global"
                                                           : "Pinned";
  out["roster"] = pe.roster;
  out["requested"] = counter_int(pe.requested);
  out["infinite"] = pe.infinite;
  out["placed_demand"] = policy.placed_demand();
  py::dict cohorts;
  for (std::size_t i = 0; i < pe.roster.size(); ++i) {
    py::list flows;
    for (const auto& s : pe.cohorts[i]) {
      const FlowRecord& rec = detail::LedgerAccess::live(pe.fg, s.id);
      py::dict f;
      f["flow_id"] = py::make_tuple(rec.flow_id.first, rec.flow_id.second);
      f["src"] = s.src;
      if (s.plan) {
        f["route"] = s.plan->route;
      } else {
        const auto edges = detail::PredDAGBuilder::ascending(*s.bundle);
        f["route"] = std::vector<EdgeId>(edges.begin(), edges.end());
      }
      f["members"] = s.plan ? py::cast(s.plan->members) : py::object(py::none());
      f["weights"] = s.plan && s.plan->placement == FlowPlacement::Weighted ? py::cast(s.plan->weights)
                                                                            : py::object(py::none());
      f["T"] = rec.delivered;
      f["delivered"] = rec.delivered;
      f["alloc"] = pairs_list(rec.edges);
      f["bundle_index"] = s.bundle_index < 0 ? py::object(py::none()) : py::object(py::int_(s.bundle_index));
      flows.append(f);
    }
    cohorts[py::int_(pe.roster[i])] = flows;
  }
  out["cohorts"] = cohorts;
  return out;
}

// The allocation-failure hook: the k-th object-domain allocation fails (Python's allocator
// wrapped for one call, then restored). The allocator is process-wide while the call runs, so
// the hook serves one caller at a time (the tests); its counters are atomics because the module
// does not hold the GIL for every thread (mod_gil_not_used).
PyMemAllocatorEx g_saved_allocator {};
std::atomic<long> g_fail_countdown {-1};
std::atomic<long> g_fail_seen {0};
std::atomic<bool> g_fail_hit {false};

bool fail_now() {
  g_fail_seen.fetch_add(1);
  if (g_fail_countdown.load() < 0) return false;
  if (g_fail_countdown.fetch_sub(1) == 0) {
    g_fail_hit.store(true);
    return true;
  }
  return false;
}
void* failing_malloc(void*, std::size_t n) { return fail_now() ? nullptr : g_saved_allocator.malloc(g_saved_allocator.ctx, n); }
void* failing_calloc(void*, std::size_t n, std::size_t size) {
  return fail_now() ? nullptr : g_saved_allocator.calloc(g_saved_allocator.ctx, n, size);
}
void* failing_realloc(void*, void* ptr, std::size_t n) {
  return fail_now() ? nullptr : g_saved_allocator.realloc(g_saved_allocator.ctx, ptr, n);
}
void forwarding_free(void*, void* ptr) { g_saved_allocator.free(g_saved_allocator.ctx, ptr); }

#endif

} // namespace

PYBIND11_MODULE(_netgraph_core, m, py::mod_gil_not_used()) {
  m.doc() =
      "NetGraph-Core: the graph on integer units, the route structures, the integer ledger, the policy\n"
      "engine FlowPolicy, the DagCache and the evaluators (capacity, batch_capacity, sensitivity).\n"
      "Every public call is all-or-nothing up to its caller (APPROACH 13), with two limits: (1) the pybind11\n"
      "object machinery is outside the guarantee: constructing a bound object from Python, or reading an\n"
      "attribute or method value converted into a new bound object or container outside the staged calls,\n"
      "allocates through pybind11's make_new_instance, which does not check the result of tp_alloc, so a\n"
      "Python allocation failure at that point crashes the process, with Core state unchanged; (2) a supplied\n"
      "building DagCache is a memo: it keeps the complete entries and statistics of the lookups that\n"
      "succeeded before a later failure of the call.";

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

  py::enum_<Direction>(m, "Direction",
                       "The orientation of a search: FORWARD from the roots, or REVERSE (rooted at the roots over\n"
                       "reversed arcs; the DAG it returns is forward-oriented, toward the roots).")
      .value("FORWARD", Direction::Forward)
      .value("REVERSE", Direction::Reverse);

  py::class_<PathSelection>(m, "PathSelection",
                            "The selection of a search: tie_break and priority are explicit; multipath = True keeps\n"
                            "every equal-cost walk (the DAG), False the one canonical single path; multi_edge = True\n"
                            "keeps every eligible parallel member, False one per (u, v) pair. A tie_break or priority\n"
                            "outside its enumeration raises ValueError.")
      // The constructor checks both enumerations (ValueError naming the field).
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
      "Graph.eligibility(excluded_nodes, excluded_edges).")
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

  py::class_<Graph, std::shared_ptr<Graph>>(
      m, "Graph",
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
      .def("uid", &Graph::uid, "process-unique id (never reused)")
      .def("resolution", &Graph::resolution, "the resolution q (units of capacity)")
      .def("num_nodes", &Graph::num_nodes)
      .def("num_edges", &Graph::num_edges)
      .def("max_cost", &Graph::max_cost, "largest edge cost (0 without edges)")
      .def("cost_gcd", &Graph::cost_gcd, "gcd of the edge costs (0 without edges)")
      .def(
          "capacity_units", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.capacity_units(), self);
          },
          "int64[E] capacities in units, id order (read-only view)")
      .def(
          "cost", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.cost(), self);
          },
          "int64[E] costs, id order (read-only view)")
      .def(
          "edge_src", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.edge_src(), self);
          },
          "int32[E] source node per edge (read-only view)")
      .def(
          "edge_dst", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.edge_dst(), self);
          },
          "int32[E] destination node per edge (read-only view)")
      .def(
          "input_index", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.input_index(), self);
          },
          "int32[E]: id -> input position (read-only view); the inverse is the caller's")
      .def(
          "row_offsets", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.row_offsets(), self);
          },
          "int32[N + 1]: the out-edges of u are the ids [row_offsets[u], row_offsets[u + 1])")
      .def(
          "in_row_offsets", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.in_row_offsets(), self);
          },
          "int32[N + 1] reverse CSR offsets")
      .def(
          "in_adj_edge_index", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.in_adj_edge_index(), self);
          },
          "int32[E]: the in-edges of v in (cost, src, input) order, as EdgeIds")
      .def(
          "eligibility",
          [](const Graph& g, const py::object& excluded_nodes, const py::object& excluded_edges) {
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
      [](const std::shared_ptr<Graph>& graph, const py::object& view, const py::object& sources,
         const py::object& destinations, const std::shared_ptr<Eligibility>& elig,
         const py::object& cap) -> py::dict {
        if (!graph) throw py::value_error("graph must be a Graph");
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
        const auto s = checked_ids(wide_integers(sources, "sources node ids"), graph->num_nodes(), "sources");
        const auto d = checked_ids(wide_integers(destinations, "destinations node ids"), graph->num_nodes(), "destinations");
        detail::check_disjoint(detail::node_set(s, graph->num_nodes(), "sources"),
                               detail::node_set(d, graph->num_nodes(), "destinations"), "sources", "destinations");
        std::optional<Amount> c;
        if (!cap.is_none()) {
          if (is_bool(cap) || !PyIndex_Check(cap.ptr())) {
            throw py::value_error("cap must be an integer number of units in [0, 2^63 - 2] or None");
          }
          py::object idx = python_int(cap, "cap");
          int overflow = 0;
          const long long value = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
          if (overflow != 0 || value < 0 || value > detail::kMaxFlowCap) {
            throw py::value_error("cap " + std::string(py::str(idx)) + " is outside [0, 2^63 - 2]");
          }
          c = static_cast<Amount>(value);
        }
        detail::FlowCoreResult r;
        {
          py::gil_scoped_release release;
          const auto view_units = v.narrow<Amount>();
          r = detail::place_global_units(*graph, view_units, s, d, elig.get(), c);
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
      "via EdgeId (one canonical form). Equality is structural. parents() and via_edges() are\n"
      "read-only views that keep the DAG alive; parent_offsets() is a fresh read-only array built on\n"
      "demand (the DAG is stored over its member nodes). Built by PredDAG.from_csr, Path.to_dag()\n"
      "and the searches.")
      .def_static(
          "from_csr",
          [](const std::shared_ptr<Graph>& graph, const py::object& offsets, const py::object& parents,
             const py::object& via) {
            if (!graph) throw py::value_error("graph must be a Graph");
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
          [](const PredDAG& d) {
            const auto off = d.parent_offsets();
            py::array_t<std::int32_t> arr(static_cast<py::ssize_t>(off.size()));
            std::copy(off.begin(), off.end(), arr.mutable_data());
            py::detail::array_proxy(arr.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
            return arr;
          },
          "int32[N + 1] CSR offsets, presented on demand from the sparse storage (a read-only array)")
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
      .def(py::init([](const std::shared_ptr<Graph>& graph, const py::object& hops) {
             if (!graph) throw py::value_error("graph must be a Graph");
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

  auto ledger_class = py::class_<FlowLedger, std::shared_ptr<FlowLedger>>(
      m, "FlowLedger",
      "The integer ledger of one graph: an immutable baseline residual in units and the records of\n"
      "the flows placed on it by policies (the ledger has no public placement or removal).\n"
      "allocated() is the sum of every record's edges and residual() = baseline - allocated(), both\n"
      "int64[E] snapshots. reset() and copy() are rejected while a policy is bound.")
      .def(py::init([](const std::shared_ptr<Graph>& graph, const py::object& residual) {
             if (!graph) throw py::value_error("graph must be a Graph");
             if (residual.is_none()) return std::make_shared<FlowLedger>(graph);
             const auto r = edge_units(residual, *graph, "residual");
             return std::make_shared<FlowLedger>(graph, r);
           }),
           py::arg("graph"), py::arg("residual") = py::none(),
           "residual: the baseline as integer units, one entry per edge in [0, capacity_units[e]] (a\n"
           "physical baseline goes through to_units first); None: the capacities. ValueError otherwise.")
      .def("uid", &FlowLedger::uid, "process-unique ledger uid (the first half of every FlowId)")
      .def(
          "graph", [](const FlowLedger& fg) { return std::const_pointer_cast<Graph>(fg.graph()); },
          "the ledger's graph")
      .def(
          "capacity_units",
          [](const py::object& self) { return readonly_view(self.cast<const FlowLedger&>().capacity_units(), self); },
          "int64[E] the graph's capacities in units (read-only view)")
      .def(
          "residual", [](const FlowLedger& fg) { return int64_array(fg.residual()); },
          "int64[E] snapshot: baseline - allocated")
      .def(
          "allocated", [](const FlowLedger& fg) { return int64_array(fg.allocated()); },
          "int64[E] snapshot: the sum of every record's edges")
      .def("records", &FlowLedger::records, "live records (bound and orphaned) as snapshots in FlowId order")
      .def(
          "record", [](const FlowLedger& fg, const py::object& flow_id) { return fg.record(flow_id_of(flow_id)); },
          py::arg("flow_id"),
          "one record by FlowId (ValueError for a FlowId of another ledger, a retired or an unknown id)")
      .def(
          "bottleneck",
          [](const FlowLedger& fg, const py::object& sources, const std::shared_ptr<Eligibility>& elig) {
            const auto s = node_ids(sources, fg.graph()->num_nodes(), "sources");
            return fg.bottleneck(s, elig.get());
          },
          py::arg("sources"), py::arg("elig") = nullptr,
          "Saturation diagnostic: the edges (statically eligible, capacity >= 1 unit) leaving the set\n"
          "reachable from the allowed sources over eligible forward arcs with residual >= 1 unit, in\n"
          "ascending EdgeId (not a cut).")
      .def(
          "reset",
          [](FlowLedger& fg, const py::object& residual) {
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
          "copy", [](const FlowLedger& fg) { return std::make_shared<FlowLedger>(fg); },
          "a new ledger (new uid) with the same baseline and records, ids kept under the new uid;\n"
          "ValueError while a policy is bound")
      .def("__copy__", [](const FlowLedger& fg) { return std::make_shared<FlowLedger>(fg); })
      .def("__deepcopy__", [](const FlowLedger& fg, const py::object&) { return std::make_shared<FlowLedger>(fg); },
           py::arg("memo"));
#if NETGRAPH_CORE_PRIVATE_HOOKS
  ledger_class
      .def(
          "_state",
          [](const FlowLedger& fg) {
            py::dict out;
            out["next_flow_id"] = detail::LedgerAccess::next_flow_id(fg);
            out["bound"] = detail::LedgerAccess::bound_owners(fg);
            out["orphans"] = detail::LedgerAccess::orphan_owners(fg);
            return out;
          },
          "private: {next_flow_id, bound, orphans} (owner ids)");
#else
  (void)ledger_class;
#endif

  // ------------------------------------------------------------------ the policy engine

  py::class_<Absolute>(m, "Absolute", "CostCeiling: an absolute, inclusive path-cost budget (cost >= 0).")
      .def(py::init([](const py::object& cost) {
             return Absolute {*int_field(cost, std::numeric_limits<Cost>::min(), std::numeric_limits<Cost>::max(),
                                         "Absolute cost", false)};
           }),
           py::arg("cost"))
      .def_readonly("cost", &Absolute::cost)
      .def("__eq__", [](const Absolute& a, const py::object& b) { return py::isinstance<Absolute>(b) && a == b.cast<Absolute>(); })
      .def("__hash__", [](const Absolute& a) { return py::hash(py::make_tuple("Absolute", a.cost)); })
      .def("__repr__", [](const Absolute& a) { return "Absolute(" + std::to_string(a.cost) + ")"; });

  py::class_<Stretch>(m, "Stretch",
                      "CostCeiling: a factor in [1, inf) of the static shortest cost under the policy's\n"
                      "Eligibility, read as the exact dyadic rational it is (compared without rounding).")
      .def(py::init([](const py::object& factor) { return Stretch {float_field(factor, "Stretch factor")}; }),
           py::arg("factor"))
      .def_readonly("factor", &Stretch::factor)
      .def("__eq__", [](const Stretch& a, const py::object& b) { return py::isinstance<Stretch>(b) && a == b.cast<Stretch>(); })
      .def("__hash__", [](const Stretch& a) { return py::hash(py::make_tuple("Stretch", a.factor)); })
      .def("__repr__", [](const Stretch& a) { return "Stretch(" + std::string(py::str(py::float_(a.factor))) + ")"; });

  py::class_<DiminishingReturns>(m, "DiminishingReturns",
                                 "The diminishing-returns guard of CapacityAware growth: stop once at least `window`\n"
                                 "growth steps were accepted in the call and the last `window` gains sum below\n"
                                 "`fraction` times the source's (finite) offer. window >= 1, fraction in (0, 1).")
      .def(py::init([](const py::object& window, const py::object& fraction) {
             DiminishingReturns d;
             d.window = static_cast<std::int32_t>(*int_field(window, std::numeric_limits<std::int32_t>::min(),
                                                             std::numeric_limits<std::int32_t>::max(),
                                                             "diminishing window", false));
             d.fraction = float_field(fraction, "diminishing fraction");
             return d;
           }),
           py::arg("window"), py::arg("fraction"))
      .def_readonly("window", &DiminishingReturns::window)
      .def_readonly("fraction", &DiminishingReturns::fraction)
      .def("__repr__", [](const DiminishingReturns& d) {
        return "DiminishingReturns(window=" + std::to_string(d.window) + ", fraction=" + std::string(py::str(py::float_(d.fraction))) + ")";
      });

  py::class_<CostOnly>(m, "CostOnly",
                       "Routing: IGP forwarding. One flow per roster source on the static shortest-path DAG\n"
                       "(the DagCache value; selection projected at creation), never refreshed, grown or\n"
                       "reoptimised.")
      .def(py::init<>())
      .def("__repr__", [](const CostOnly&) { return std::string("CostOnly()"); });

  py::class_<CapacityAware>(
      m, "CapacityAware",
      "Routing: distributed CSPF tunnels. flow_count flows created by the canonical search on the\n"
      "selection view (residual >= max(1, min_flow_bandwidth) units, the max_path_cost ceiling),\n"
      "grown up to max_flow_count (None: no ceiling), optionally reoptimised (reroute with retention\n"
      "fallback) and guarded by `diminishing` (None: off).")
      .def(py::init([](const py::object& flow_count, const py::object& max_flow_count, const py::object& max_path_cost,
                       const py::object& min_flow_bandwidth, bool reoptimize, const py::object& diminishing) {
             CapacityAware ca;
             constexpr auto lo32 = std::numeric_limits<std::int32_t>::min();
             constexpr auto hi32 = std::numeric_limits<std::int32_t>::max();
             ca.flow_count = static_cast<std::int32_t>(*int_field(flow_count, lo32, hi32, "flow_count", false));
             if (auto mx = int_field(max_flow_count, lo32, hi32, "max_flow_count", true)) {
               ca.max_flow_count = static_cast<std::int32_t>(*mx);
             }
             if (!max_path_cost.is_none()) {
               if (py::isinstance<Absolute>(max_path_cost)) {
                 ca.max_path_cost = max_path_cost.cast<Absolute>();
               } else if (py::isinstance<Stretch>(max_path_cost)) {
                 ca.max_path_cost = max_path_cost.cast<Stretch>();
               } else {
                 throw py::type_error("max_path_cost must be Absolute, Stretch or None");
               }
             }
             ca.min_flow_bandwidth = *int_field(min_flow_bandwidth, std::numeric_limits<Amount>::min(),
                                                std::numeric_limits<Amount>::max(), "min_flow_bandwidth (units)", false);
             ca.reoptimize = reoptimize;
             if (!diminishing.is_none()) {
               if (!py::isinstance<DiminishingReturns>(diminishing)) throw py::type_error("diminishing must be DiminishingReturns or None");
               ca.diminishing = diminishing.cast<DiminishingReturns>();
             }
             return ca;
           }),
           py::kw_only(), py::arg("flow_count") = 1, py::arg("max_flow_count") = py::none(),
           py::arg("max_path_cost") = py::none(), py::arg("min_flow_bandwidth") = 0, py::arg("reoptimize") = false,
           py::arg("diminishing") = py::none())
      .def_readonly("flow_count", &CapacityAware::flow_count)
      .def_readonly("max_flow_count", &CapacityAware::max_flow_count)
      .def_property_readonly("max_path_cost",
                             [](const CapacityAware& ca) -> py::object {
                               if (!ca.max_path_cost.has_value()) return py::none();
                               return std::visit([](const auto& x) { return py::cast(x); }, *ca.max_path_cost);
                             })
      .def_readonly("min_flow_bandwidth", &CapacityAware::min_flow_bandwidth)
      .def_readonly("reoptimize", &CapacityAware::reoptimize)
      .def_readonly("diminishing", &CapacityAware::diminishing)
      .def("__repr__", [](const CapacityAware& ca) {
        return "CapacityAware(flow_count=" + std::to_string(ca.flow_count) + ", max_flow_count=" +
               (ca.max_flow_count ? std::to_string(*ca.max_flow_count) : std::string("None")) + ")";
      });

  py::class_<Global>(m, "Global",
                     "Routing: the omniscient controller. Every call re-solves the minimum-cost flow of the\n"
                     "cumulative value over the whole eligible graph (placement MaxFlow only).")
      .def(py::init<>())
      .def("__repr__", [](const Global&) { return std::string("Global()"); });

  py::class_<Pinned>(m, "Pinned",
                     "Routing: explicit route structures, one flow per usable bundle in supply order (single\n"
                     "endpoints; a Path converts with to_dag()). The list is copied at construction.")
      .def(py::init([](const py::object& bundles) {
             if (not_a_sequence(bundles)) throw py::type_error("Pinned bundles must be a sequence of PredDAG");
             Pinned p;
             for (const auto& b : py::reinterpret_borrow<py::sequence>(bundles)) {
               if (!py::isinstance<PredDAG>(b)) {
                 throw py::type_error("Pinned bundles must be PredDAG values (a Path converts with to_dag())");
               }
               p.bundles.push_back(b.cast<std::shared_ptr<PredDAG>>());
             }
             return p;
           }),
           py::arg("bundles"))
      .def_property_readonly("bundles",
                             [](const Pinned& p) {
                               py::list out;
                               for (const auto& b : p.bundles) out.append(std::const_pointer_cast<PredDAG>(b));
                               return out;
                             })
      .def("__repr__", [](const Pinned& p) { return "Pinned(bundles=" + std::to_string(p.bundles.size()) + ")"; });

  py::class_<FlowPolicyConfig>(
      m, "FlowPolicyConfig",
      "The four explicit axes of a policy: routing (CostOnly, CapacityAware, Global or Pinned), placement\n"
      "(FlowPlacement), admission (Admission) and selection (PathSelection). A placement or admission\n"
      "outside its enumeration raises ValueError.")
      .def(py::init([](const py::object& routing, FlowPlacement placement, Admission admission,
                       const PathSelection& selection) {
             // The routing's structure first (TypeError), then the constructor's enumeration checks.
             auto r = routing_of(routing);
             return FlowPolicyConfig {std::move(r), placement, admission, selection};
           }),
           py::kw_only(), py::arg("routing"), py::arg("placement"), py::arg("admission"), py::arg("selection"))
      .def_property_readonly("routing", [](const FlowPolicyConfig& c) { return routing_obj(c.routing); })
      .def_readonly("placement", &FlowPolicyConfig::placement)
      .def_readonly("admission", &FlowPolicyConfig::admission)
      .def_readonly("selection", &FlowPolicyConfig::selection)
      .def("__repr__", [](const FlowPolicyConfig& c) {
        return "FlowPolicyConfig(routing=" + std::string(py::repr(routing_obj(c.routing))) + ")";
      });

  py::class_<OutputOptions>(m, "OutputOptions",
                            "Optional outputs of a placement (all false by default; absent outputs are None):\n"
                            "edge_flows (the policy's own allocation per edge after the call), residual and\n"
                            "bottleneck (the ledger's after the call); on Global, bottleneck adds min_cut.")
      .def(py::init<bool, bool, bool>(), py::kw_only(), py::arg("edge_flows") = false, py::arg("residual") = false,
           py::arg("bottleneck") = false)
      .def_readonly("edge_flows", &OutputOptions::edge_flows)
      .def_readonly("residual", &OutputOptions::residual)
      .def_readonly("bottleneck", &OutputOptions::bottleneck);

  py::enum_<Termination>(m, "Termination",
                         "The local outcome of a source and a call's termination (APPROACH 9.4): COMPLETE (the\n"
                         "offer is served), LIMIT_FLOWS, LIMIT_DIMINISHING, NO_ROUTE, EXHAUSTED.")
      .value("COMPLETE", Termination::Complete)
      .value("LIMIT_FLOWS", Termination::LimitFlows)
      .value("LIMIT_DIMINISHING", Termination::LimitDiminishing)
      .value("NO_ROUTE", Termination::NoRoute)
      .value("EXHAUSTED", Termination::Exhausted);

  py::class_<PlacementResult>(
      m, "PlacementResult",
      "One placement call: placed (units), remaining (request - placed; float('inf') for +inf),\n"
      "termination, local_outcomes [(source, Termination)] per roster source, routes_changed,\n"
      "weights_changed, drops [(EdgeId, units)] of this call; cost_distribution [(cost, units)] on\n"
      "Global (the marginal curve of the policy's cumulative flow; a zero call reports it unchanged), else\n"
      "None; the optional outputs\n"
      "edge_flows, residual (int64[E]), bottleneck, min_cut (Global with bottleneck; None when the\n"
      "solve was not maximal), None unless requested.")
      .def_readonly("placed", &PlacementResult::placed)
      .def_property_readonly("remaining", [](const PlacementResult& r) { return remaining_obj(r.remaining); })
      .def_readonly("termination", &PlacementResult::termination)
      .def_readonly("local_outcomes", &PlacementResult::local_outcomes)
      .def_readonly("routes_changed", &PlacementResult::routes_changed)
      .def_readonly("weights_changed", &PlacementResult::weights_changed)
      .def_property_readonly("drops", [](const PlacementResult& r) { return pairs_list(r.drops); })
      .def_readonly("cost_distribution", &PlacementResult::cost_distribution)
      .def_property_readonly("edge_flows",
                             [](const PlacementResult& r) -> py::object {
                               return r.edge_flows ? py::object(int64_array(*r.edge_flows)) : py::object(py::none());
                             })
      .def_property_readonly("residual",
                             [](const PlacementResult& r) -> py::object {
                               return r.residual ? py::object(int64_array(*r.residual)) : py::object(py::none());
                             })
      .def_readonly("bottleneck", &PlacementResult::bottleneck)
      .def_readonly("min_cut", &PlacementResult::min_cut)
      .def("__repr__", [](const PlacementResult& r) {
        return "PlacementResult(placed=" + std::to_string(r.placed) + ", termination=" +
               std::string(py::str(py::cast(r.termination).attr("name"))) + ")";
      });

  py::class_<DagCache::Stats>(m, "DagCacheStats",
                              "DagCache statistics (atomic per field, not a consistent snapshot): hits, misses,\n"
                              "entries, bytes (the entries' retained payload), live_values (values alive anywhere).")
      .def_readonly("hits", &DagCache::Stats::hits)
      .def_readonly("misses", &DagCache::Stats::misses)
      .def_readonly("entries", &DagCache::Stats::entries)
      .def_readonly("bytes", &DagCache::Stats::bytes)
      .def_readonly("live_values", &DagCache::Stats::live_values);

  py::register_exception<DagCacheMiss>(m, "DagCacheMiss", PyExc_Exception);
  m.attr("DagCacheMiss").attr("__doc__") =
      "A lookup in a frozen DagCache whose key is absent (APPROACH 10): a dedicated exception (not a\n"
      "RuntimeError or ValueError); the evaluators propagate it unchanged, never as a status.";

  py::class_<DagCache, std::shared_ptr<DagCache>>(
      m, "DagCache",
      "A cache of static shortest-path supports bound to one graph, keyed (roots as a set, direction,\n"
      "Eligibility uid). Building mode (a new cache): one writer, no concurrent readers; a miss searches\n"
      "and inserts. Frozen mode (freeze()): read-only and shareable by every worker; a miss counts one miss\n"
      "and raises DagCacheMiss; clear() raises ValueError; thaw() returns to building once the caller has\n"
      "joined every reader. CostOnly policies and the evaluators read it; selection is projected at\n"
      "creation, never cached. clear() empties a building cache (the hit and miss counters persist); held\n"
      "values stay valid. A supplied building cache is a memo (APPROACH 13, limit 2): a public call that\n"
      "fails keeps the complete entries and statistics of the lookups that succeeded before; never a\n"
      "partial entry, and a retained entry never changes a result. A frozen cache's statistics follow the\n"
      "public call that reads it: its hits count when the call succeeds, its Python result included; a\n"
      "failed call changes no statistic, except that a call raising DagCacheMiss counts one miss.")
      .def(py::init([](const std::shared_ptr<Graph>& graph) {
             if (!graph) throw py::value_error("graph must be a Graph");
             return std::make_shared<DagCache>(*graph);
           }),
           py::arg("graph"))
      .def("graph_uid", &DagCache::graph_uid, "uid of the graph this cache is bound to")
      .def("freeze", &DagCache::freeze, "enter frozen mode (no effect when frozen)")
      .def("thaw", &DagCache::thaw, "return to building mode (no effect when building); join every reader first")
      .def("clear", &DagCache::clear, "empty a building cache; ValueError when frozen")
      .def("stats", &DagCache::stats);

  auto flow_policy = py::class_<FlowPolicy>(
      m, "FlowPolicy",
      "The one placement engine (APPROACH 9), bound to its ledger at construction: place(volume,\n"
      "outputs) runs the transaction and commits one atomic exchange (a failed call raises with the ledger\n"
      "and the policy as they were, its Python result included; limits in APPROACH 13: pybind11's\n"
      "unchecked tp_alloc in object construction and conversions outside the staged calls, and a supplied\n"
      "building DagCache, a memo that keeps its successful lookups); remove() reverts every flow and\n"
      "opens a new epoch; dropping the policy without remove() leaves its records as orphans.\n"
      "sources and destinations are one NodeId or a sequence (disjoint sets). ValueError for invalid\n"
      "endpoints, configurations, foreign objects and requests; TypeError for wrong argument types.")
      .def(py::init([](const std::shared_ptr<FlowLedger>& ledger, const py::object& sources,
                       const py::object& destinations, const FlowPolicyConfig& config,
                       const std::shared_ptr<Eligibility>& elig, const std::shared_ptr<DagCache>& cache) {
             if (!ledger) throw py::value_error("ledger must be a FlowLedger");
             const auto n = ledger->graph()->num_nodes();
             const auto s = node_ids(sources, n, "sources");
             const auto d = node_ids(destinations, n, "destinations");
             return std::make_unique<FlowPolicy>(*ledger, s, d, config, elig.get(), cache.get());
           }),
           py::keep_alive<1, 2>(), py::keep_alive<1, 6>(), py::keep_alive<1, 7>(), py::arg("ledger"),
           py::arg("sources"), py::arg("destinations"), py::arg("config"), py::arg("elig") = nullptr,
           py::arg("cache") = nullptr)
      .def(
          "place",
          [](FlowPolicy& p, const py::object& volume, const py::object& outputs) {
            const auto v = request_units(volume);
            const auto o = outputs_of(outputs);
            // All or nothing across the boundary: the result's Python object exists before the
            // call, and after the commit nothing allocates (the value moves into it). The GIL
            // stays held: the call mutates the ledger (one writer per ledger).
            auto [obj, slot] = new_object<PlacementResult>();
            *slot = p.place(v, o);
            return std::move(obj);
          },
          py::arg("volume"), py::arg("outputs") = py::none(),
          "place(volume: int units or float('inf'), outputs: OutputOptions | None) -> PlacementResult")
      .def("remove", &FlowPolicy::remove, "revert every flow of the policy, retire the ids, open a new epoch")
      .def("placed_demand", &FlowPolicy::placed_demand, "the sum of delivered over the policy's live records")
      .def("records", &FlowPolicy::records, "the policy's live records as snapshots in FlowId order")
      .def(
          "record", [](const FlowPolicy& p, const py::object& flow_id) { return p.record(flow_id_of(flow_id)); },
          py::arg("flow_id"),
          "one of the policy's live records by FlowId (ValueError for a FlowId of another ledger or policy, a\n"
          "retired or an unknown id)")
      .def("usable_bundles", &FlowPolicy::usable_bundles,
           "Pinned: the bundles with a src -> dst walk after static pruning (0 for the other modes)")
      .def_property_readonly("config", &FlowPolicy::config);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  flow_policy
      .def(
          "_place_traced",
          [](FlowPolicy& p, const py::object& volume, const py::object& outputs) {
            const auto v = request_units(volume);
            const auto o = outputs_of(outputs);
            // As place(): every Python object exists before anything is published; the trace is
            // converted at the engine's staged point, after the call is complete and before the
            // commit, so a failed conversion changes nothing.
            auto [obj, slot] = new_object<PlacementResult>();
            py::tuple pair(2);
            py::object trace_py;
            detail::Trace trace;
            const std::function<void()> staged = [&] { trace_py = trace_obj(trace); };
            *slot = detail::place_traced(p, v, o, trace, &staged);
            PyTuple_SET_ITEM(pair.ptr(), 0, obj.release().ptr());
            PyTuple_SET_ITEM(pair.ptr(), 1, trace_py.release().ptr());
            return pair;
          },
          py::arg("volume"), py::arg("outputs") = py::none(),
          "private: place() with the transaction's trace (prototype 2's keys)")
      .def("_state", &state_obj, "private: the epoch state (roster, requested, +inf, cohorts)")
      .def(
          "_owner", [](const FlowPolicy& p) { return detail::PolicyAccess::engine(p).owner; },
          "private: the ledger-local owner id (FlowRecord.owner)");
#else
  (void)flow_policy;
#endif


  py::class_<PyCapacityOptions>(m, "CapacityOptions",
                                "The configuration, the static Eligibility (None: all allowed), the DagCache (None:\n"
                                "private per worker) and the optional outputs of capacity() and batch_capacity().")
      .def(py::init([](const FlowPolicyConfig& config, const std::shared_ptr<Eligibility>& elig,
                       const std::shared_ptr<DagCache>& cache, const py::object& outputs) {
             return PyCapacityOptions {config, elig, cache, outputs_of(outputs)};
           }),
           py::arg("config"), py::kw_only(), py::arg("elig") = nullptr, py::arg("cache") = nullptr,
           py::arg("outputs") = py::none())
      .def_readonly("config", &PyCapacityOptions::config)
      .def_readonly("elig", &PyCapacityOptions::elig)
      .def_readonly("cache", &PyCapacityOptions::cache)
      .def_readonly("outputs", &PyCapacityOptions::outputs);

  py::class_<PySensitivityOptions>(m, "SensitivityOptions",
                                   "The configuration, Eligibility, DagCache and baseline outputs of sensitivity().")
      .def(py::init([](const FlowPolicyConfig& config, const std::shared_ptr<Eligibility>& elig,
                       const std::shared_ptr<DagCache>& cache, const py::object& outputs) {
             return PySensitivityOptions {config, elig, cache, outputs_of(outputs)};
           }),
           py::arg("config"), py::kw_only(), py::arg("elig") = nullptr, py::arg("cache") = nullptr,
           py::arg("outputs") = py::none())
      .def_readonly("config", &PySensitivityOptions::config)
      .def_readonly("elig", &PySensitivityOptions::elig)
      .def_readonly("cache", &PySensitivityOptions::cache)
      .def_readonly("outputs", &PySensitivityOptions::outputs);

  py::enum_<SensitivityStatus>(m, "SensitivityStatus",
                               "EVALUATED (the loss is the baseline minus the capacity without the edge) or\n"
                               "UNEVALUABLE (the deletion raised a model exception: loss None, reason the message).")
      .value("EVALUATED", SensitivityStatus::Evaluated)
      .value("UNEVALUABLE", SensitivityStatus::Unevaluable);

  py::class_<SensitivityEntry>(m, "SensitivityEntry", "One candidate edge: edge, loss (int units or None), status, reason.")
      .def_readonly("edge", &SensitivityEntry::edge)
      .def_readonly("loss", &SensitivityEntry::loss)
      .def_readonly("status", &SensitivityEntry::status)
      .def_property_readonly("reason", [](const SensitivityEntry& e) -> py::object {
        return e.status == SensitivityStatus::Unevaluable ? py::object(py::str(e.reason)) : py::object(py::none());
      });

  py::class_<SensitivityResult>(m, "SensitivityResult",
                                "baseline (the PlacementResult of capacity()) and entries (ascending EdgeId; an\n"
                                "absent edge has loss 0).")
      .def_readonly("baseline", &SensitivityResult::baseline)
      .def_readonly("entries", &SensitivityResult::entries);

  m.attr("MAX_FLOWS_PER_CALL") = py::int_(kMaxFlowsPerCall);

  m.def(
      "capacity",
      [](const std::shared_ptr<Graph>& graph, const py::object& sources, const py::object& destinations,
         const py::object& opts) {
        if (!graph) throw py::value_error("graph must be a Graph");
        const auto o = options_of<PyCapacityOptions>(opts, "CapacityOptions");
        const auto s = node_ids(sources, graph->num_nodes(), "sources");
        const auto d = node_ids(destinations, graph->num_nodes(), "destinations");
        const CapacityOptions c {o.config, o.elig.get(), o.cache.get(), o.outputs};
        // One public call on a frozen cache's statistics, its result object included (NOTES R27).
        return detail::frozen_call(c.cache, [&](detail::FrozenTally& tally) {
          PlacementResult r;
          {
            py::gil_scoped_release release;
            r = detail::capacity_in_call(*graph, s, d, c, tally);
          }
          return object_of(std::move(r));
        });
      },
      py::arg("graph"), py::arg("sources"), py::arg("destinations"), py::arg("opts"),
      "place(+inf) of a FlowPolicy with opts.config on a fresh ledger of the graph -> PlacementResult.\n"
      "The exact maximum flow is the Global + MaxFlow configuration. A supplied frozen cache raises\n"
      "DagCacheMiss for an absent key, unchanged.\n"
      "A failed call raises; a supplied building cache keeps the entries and statistics of the lookups\n"
      "that succeeded before the failure (a memo: APPROACH 13, limit 2); a supplied frozen cache's\n"
      "statistics are unchanged, except one miss when the call raises DagCacheMiss.");

  m.def(
      "batch_capacity",
      [](const std::shared_ptr<Graph>& graph, const py::object& pairs, const py::object& opts,
         const py::object& per_pair_elig, const py::object& workers) {
        if (!graph) throw py::value_error("graph must be a Graph");
        // The arguments first, then each pair completely and in pair order (the same two parts
        // as the C++ validation), so the lowest failing pair index is reported.
        const auto o = options_of<PyCapacityOptions>(opts, "CapacityOptions");
        const auto w = workers_of(workers, false, static_cast<bool>(o.cache));
        const CapacityOptions c {o.config, o.elig.get(), o.cache.get(), o.outputs};
        detail::check_batch_arguments(*graph, c, w);
        if (not_a_sequence(pairs)) throw py::value_error("pairs must be a sequence of (sources, destinations)");
        const auto items = py::reinterpret_borrow<py::sequence>(pairs);
        const auto num_pairs = py::len(items);
        std::optional<py::sequence> eligs;
        if (!per_pair_elig.is_none()) {
          // An explicit sequence has one entry per pair (None entries allowed), also when empty.
          if (not_a_sequence(per_pair_elig)) throw py::value_error("per_pair_elig must be a sequence of Eligibility or None");
          eligs = py::reinterpret_borrow<py::sequence>(per_pair_elig);
          if (py::len(*eligs) != num_pairs) {
            throw py::value_error("per_pair_elig has " + std::to_string(py::len(*eligs)) + " entries for " +
                                  std::to_string(num_pairs) + " pairs");
          }
        }
        std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> ps;
        std::vector<std::shared_ptr<Eligibility>> keep;
        std::vector<const Eligibility*> per;
        const auto n = graph->num_nodes();
        for (std::size_t i = 0; i < num_pairs; ++i) {
          const py::object item = items[i];
          if (not_a_sequence(item) || py::len(item) != 2) {
            throw py::value_error("each pair must be (sources, destinations)");
          }
          auto t = py::reinterpret_borrow<py::sequence>(item);
          auto src = node_ids(t[0], n, "sources");
          auto dst = node_ids(t[1], n, "destinations");
          const Eligibility* e = nullptr;
          if (eligs) {
            const py::object entry = (*eligs)[i];
            if (!entry.is_none()) {
              if (!py::isinstance<Eligibility>(entry)) throw py::type_error("per_pair_elig entries must be Eligibility or None");
              keep.push_back(entry.cast<std::shared_ptr<Eligibility>>());
              e = keep.back().get();
            }
            per.push_back(e);
          }
          detail::check_batch_pair(*graph, src, dst, e, i);
          ps.emplace_back(std::move(src), std::move(dst));
        }
        return detail::frozen_call(c.cache, [&](detail::FrozenTally& tally) {
          std::vector<PlacementResult> rs;
          {
            py::gil_scoped_release release;
            rs = detail::batch_capacity_in_call(*graph, ps, c, per, w, tally);
          }
          py::list out(rs.size());
          for (std::size_t i = 0; i < rs.size(); ++i) {
            PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), object_of(std::move(rs[i])).release().ptr());
          }
          return py::object(std::move(out));
        });
      },
      py::arg("graph"), py::arg("pairs"), py::arg("opts"), py::arg("per_pair_elig") = py::none(),
      py::arg("workers") = py::none(),
      "capacity() per (sources, destinations) pair in pair order, opts.elig intersected with\n"
      "per_pair_elig[i] (None entries allowed; a given sequence has one entry per pair, also when\n"
      "empty); workers >= 1 (None: NGRAPH_CORE_BATCH_THREADS or 1; 1 with a supplied cache).\n"
      "Everything is validated before a worker starts, the arguments first, then each pair\n"
      "completely in pair order, so the error of the lowest failing pair index is raised; Pinned is\n"
      "rejected; an empty batch returns []. A supplied building cache needs workers == 1; a frozen cache\n"
      "is shared by every worker (pass workers explicitly) and its DagCacheMiss reaches the caller\n"
      "unchanged (the lowest failing pair's); a per-pair intersection that differs from both operands\n"
      "has a new uid, absent from a frozen cache (pass effective per-pair values with elig None, or a\n"
      "common value containing them).\n"
      "A failed call raises; a supplied building cache keeps the entries and statistics of the lookups\n"
      "that succeeded before the failure (a memo: APPROACH 13, limit 2); a supplied frozen cache's\n"
      "statistics are unchanged, except one miss when the call raises DagCacheMiss.");

  m.def(
      "sensitivity",
      [](const std::shared_ptr<Graph>& graph, const py::object& sources, const py::object& destinations,
         const py::object& opts, const py::object& workers) {
        if (!graph) throw py::value_error("graph must be a Graph");
        const auto o = options_of<PySensitivityOptions>(opts, "SensitivityOptions");
        const auto s = node_ids(sources, graph->num_nodes(), "sources");
        const auto d = node_ids(destinations, graph->num_nodes(), "destinations");
        const auto w = workers_of(workers, true, static_cast<bool>(o.cache));
        const SensitivityOptions c {o.config, o.elig.get(), o.cache.get(), o.outputs};
        return detail::frozen_call(c.cache, [&](detail::FrozenTally& tally) {
          SensitivityResult r;
          {
            py::gil_scoped_release release;
            r = detail::sensitivity_in_call(*graph, s, d, c, w, tally);
          }
          return object_of(std::move(r));
        });
      },
      py::arg("graph"), py::arg("sources"), py::arg("destinations"), py::arg("opts"), py::arg("workers") = py::none(),
      "The baseline capacity() and one entry per candidate edge (Global: the baseline's positive-flow\n"
      "edges; other models: every allowed edge of at least 1 unit): the loss is the baseline's placed\n"
      "minus capacity() without the edge, in signed units, never thresholded (Global evaluates it by the\n"
      "flow core's deletion re-solve from the baseline's final state, and CostOnly gives edges no route\n"
      "uses loss 0 without a placement, result-equivalent forms of that definition); a placed deletion\n"
      "raising a model exception is UNEVALUABLE with the message as reason.\n"
      "workers as batch_capacity (NGRAPH_CORE_SENSITIVITY_THREADS); a CostOnly deletion (the one\n"
      "model reading the cache) looks up a new Eligibility uid, so a frozen cache raises DagCacheMiss\n"
      "there (unchanged, never UNEVALUABLE); other models read no cache.\n"
      "A failed call raises; a supplied building cache keeps the entries and statistics of the lookups\n"
      "that succeeded before the failure (a memo: APPROACH 13, limit 2); a supplied frozen cache's\n"
      "statistics are unchanged, except one miss when the call raises DagCacheMiss.");

  py::class_<PyKspOptions>(
      m, "KspOptions",
      "The options of ksp (APPROACH 6, rule 10): k (>= 1) and selection are explicit (its multi_edge is\n"
      "the member rule, Expand or Collapse; multipath is ignored); max_path_cost (None, Absolute(c) with\n"
      "c >= 0, or Stretch(f) with f finite in [1, inf) times the static shortest cost under the\n"
      "Eligibility, residual ignored; inclusive); residual (int64 units per EdgeId, each in\n"
      "[0, capacity_units]: gating at 1 unit and the priority values); elig (None: everything allowed).\n"
      "The fields are validated by ksp, in the contract's order.")
      .def(py::init([](const py::object& k, const PathSelection& selection, const py::object& max_path_cost,
                       const py::object& residual, const std::shared_ptr<Eligibility>& elig) {
             return PyKspOptions {k, selection, max_path_cost, residual, elig};
           }),
           py::arg("k"), py::arg("selection"), py::kw_only(), py::arg("max_path_cost") = py::none(),
           py::arg("residual") = py::none(), py::arg("elig") = nullptr)
      .def_readonly("k", &PyKspOptions::k)
      .def_readonly("selection", &PyKspOptions::selection)
      .def_readonly("max_path_cost", &PyKspOptions::max_path_cost)
      .def_readonly("residual", &PyKspOptions::residual)
      .def_readonly("elig", &PyKspOptions::elig);

  m.def(
      "ksp",
      [](const std::shared_ptr<Graph>& graph, const py::object& src, const py::object& dst, const py::object& opts) {
        if (!graph) throw py::value_error("graph must be a Graph");
        const auto o = options_of<PyKspOptions>(opts, "KspOptions");
        const auto n = graph->num_nodes();
        // APPROACH 6's order on unnarrowed values (Core repeats every check): the selection, src,
        // dst, src == dst, k, the ceiling, the residual, then a foreign Eligibility.
        detail::check_selection(o.selection);
        const NodeId s = endpoint(src, n, "src");
        const NodeId t = endpoint(dst, n, "dst");
        if (s == t) throw py::value_error(detail::overlap_error(s, "src", "dst"));
        KspOptions c(count_of(o.k, "k"), o.selection);
        if (!o.max_path_cost.is_none()) {
          if (py::isinstance<Absolute>(o.max_path_cost)) {
            c.max_path_cost = o.max_path_cost.cast<Absolute>();
          } else if (py::isinstance<Stretch>(o.max_path_cost)) {
            c.max_path_cost = o.max_path_cost.cast<Stretch>();
          } else {
            throw py::type_error("max_path_cost must be Absolute, Stretch or None (got " +
                                 std::string(py::str(py::type::of(o.max_path_cost).attr("__name__"))) + ")");
          }
          detail::check_ceiling(c.max_path_cost);
        }
        std::optional<std::vector<Amount>> res;
        if (!o.residual.is_none()) res = edge_units(o.residual, *graph, "residual");
        if (res.has_value()) c.residual = std::span<const Amount>(*res);
        c.elig = o.elig.get();
        std::vector<std::pair<Cost, Path>> paths;
        {
          py::gil_scoped_release release;
          paths = ksp(*graph, s, t, c);
        }
        // Every object through a checked allocation: a failure raises.
        py::list out(paths.size());
        for (std::size_t i = 0; i < paths.size(); ++i) {
          py::tuple item(2);
          py::int_ cost(paths[i].first);
          py::object path = object_from<Path>(std::move(paths[i].second));
          PyTuple_SET_ITEM(item.ptr(), 0, cost.release().ptr());
          PyTuple_SET_ITEM(item.ptr(), 1, path.release().ptr());
          PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), item.release().ptr());
        }
        return out;
      },
      py::arg("graph"), py::arg("src"), py::arg("dst"), py::arg("opts"),
      "The K shortest simple src -> dst paths of the member graph (APPROACH 6, rule 10) -> [(cost, Path)],\n"
      "in increasing (cost, -bottleneck under PREFER_HIGHER_PRIORITY, hops, edge-id sequence), the\n"
      "bottleneck being the minimum priority value along a path: the minimal prefix of that order,\n"
      "truncated at k and at the inclusive ceiling; each Path has single-member hops. Validated in this\n"
      "order before any work (ValueError): src and dst in range, src == dst, k >= 1, the ceiling, the\n"
      "residual (length, then every entry), a foreign Eligibility. An excluded endpoint, no static route\n"
      "under Stretch, or no path gives [].");

  m.def(
      "resolve_to_paths",
      [](const std::shared_ptr<Graph>& graph, const std::shared_ptr<PredDAG>& dag, const py::object& src,
         const py::object& dst, bool split_parallel_edges, const py::object& max_paths) {
        if (!graph) throw py::value_error("graph must be a Graph");
        if (!dag) throw py::value_error("dag must be a PredDAG");
        if (dag->graph_uid() != graph->uid()) throw py::value_error("foreign PredDAG: built for another graph");
        const auto n = graph->num_nodes();
        const NodeId s = endpoint(src, n, "src");
        const NodeId t = endpoint(dst, n, "dst");
        if (s == t) throw py::value_error(detail::overlap_error(s, "src", "dst"));
        std::optional<std::int64_t> mp;
        if (!max_paths.is_none()) mp = count_of(max_paths, "max_paths");
        std::vector<Path> paths;
        {
          py::gil_scoped_release release;
          paths = resolve_to_paths(*graph, *dag, s, t, split_parallel_edges, mp);
        }
        py::list out(paths.size());
        for (std::size_t i = 0; i < paths.size(); ++i) {
          PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), object_from<Path>(std::move(paths[i])).release().ptr());
        }
        return out;
      },
      py::arg("graph"), py::arg("dag"), py::arg("src"), py::arg("dst"), py::arg("split_parallel_edges") = false,
      py::arg("max_paths") = py::none(),
      "The src -> dst walks of a DAG (APPROACH 5) -> [Path], in increasing (cost, hops, edge-id sequence):\n"
      "a hop costs its cheapest member and is keyed by its lowest member id. split_parallel_edges = False\n"
      "keeps every parallel member of a (u, v) pair in one hop (members ascending); True gives one path\n"
      "per member combination in the same order. max_paths (>= 1) truncates; None lists every walk.\n"
      "ValueError: a DAG of another graph, src or dst out of range, src == dst, max_paths < 1. An\n"
      "endpoint outside the DAG, or no walk, gives [].");

  // ------------------------------------------------------------------ private recorder entries

  // ------------------------------------------------------------------ the public search

  m.attr("UNREACHABLE") = py::int_(kUnreachable);

  py::class_<PySpfOptions>(
      m, "SpfOptions",
      "The options of spf (APPROACH 6): direction and selection are explicit; targets (one NodeId or a\n"
      "sequence; None: a full search) stops the search at the nearest targets (anycast); residual (int64\n"
      "units per EdgeId, each in [0, capacity_units]) adds residual gating at min_residual units (an\n"
      "integer >= 1; it gates a supplied residual only) and supplies the priority values; elig (None:\n"
      "everything allowed). A direction outside its enumeration raises ValueError here; the other fields\n"
      "are validated by spf, in the contract's order.")
      .def(py::init([](Direction direction, const PathSelection& selection, const py::object& targets,
                       const py::object& residual, const py::object& min_residual,
                       const std::shared_ptr<Eligibility>& elig) {
             detail::check_member(direction, "direction");
             return PySpfOptions {direction, selection, targets, residual, min_residual, elig};
           }),
           py::arg("direction"), py::arg("selection"), py::kw_only(), py::arg("targets") = py::none(),
           py::arg("residual") = py::none(), py::arg("min_residual") = 1, py::arg("elig") = nullptr)
      .def_readonly("direction", &PySpfOptions::direction)
      .def_readonly("selection", &PySpfOptions::selection)
      .def_readonly("targets", &PySpfOptions::targets)
      .def_readonly("residual", &PySpfOptions::residual)
      .def_readonly("min_residual", &PySpfOptions::min_residual)
      .def_readonly("elig", &PySpfOptions::elig);

  py::class_<PySpfResult>(m, "SpfResult",
                          "The result of spf: costs (full searches: a read-only int64 array, UNREACHABLE for the\n"
                          "nodes not reached; None with targets), dag (PredDAG: the unique DAG of the edges tight in\n"
                          "the traversal direction, in their original orientation; with targets pruned to the walks\n"
                          "between the roots and a nearest target, from the target to a root for REVERSE, or the\n"
                          "canonical single path for multipath = False) and target_cost (with targets: the nearest\n"
                          "targets' cost, None when no target is reached; None for a full search).")
      .def_property_readonly("costs", [](const PySpfResult& r) { return r.costs ? r.costs : py::object(py::none()); })
      .def_property_readonly("dag", [](const PySpfResult& r) { return r.dag; })
      .def_property_readonly("target_cost",
                             [](const PySpfResult& r) { return r.target_cost ? r.target_cost : py::object(py::none()); })
      .def("__repr__", [](const PySpfResult& r) {
        std::string tc = r.target_cost && !r.target_cost.is_none() ? std::string(py::str(r.target_cost)) : "None";
        return "SpfResult(costs=" + std::string(r.costs ? "array" : "None") +
               ", dag_entries=" + std::to_string(r.dag.cast<const PredDAG&>().via_edges().size()) +
               ", target_cost=" + tc + ")";
      });

  m.def(
      "spf",
      [](const std::shared_ptr<Graph>& graph, const py::object& roots, const py::object& opts) {
        if (!graph) throw py::value_error("graph must be a Graph");
        const auto o = options_of<PySpfOptions>(opts, "SpfOptions");
        const auto n = graph->num_nodes();
        // APPROACH 6's order on unnarrowed values (Core repeats every check): the enumerations,
        // the roots, the targets, their overlap, a single path without targets, the residual, then
        // min_residual and a foreign Eligibility.
        detail::check_member(o.direction, "direction");
        detail::check_selection(o.selection);
        const auto r = node_ids(roots, n, "roots");
        std::optional<std::vector<NodeId>> t;
        if (!o.targets.is_none()) {
          t = node_ids(o.targets, n, "targets");
          detail::check_disjoint(detail::node_set(r, n, "roots"), detail::node_set(*t, n, "targets"), "roots",
                                 "targets");
        } else if (!o.selection.multipath) {
          throw py::value_error("multipath = false requires targets (the single path joins a root to a target)");
        }
        std::optional<std::vector<Amount>> res;
        if (!o.residual.is_none()) res = edge_units(o.residual, *graph, "residual");
        const auto minr = python_int(o.min_residual, "min_residual");
        if (PyObject_RichCompareBool(minr.ptr(), py::int_(1).ptr(), Py_LT) == 1) {
          throw py::value_error("min_residual must be at least 1 unit (got " + std::string(py::str(minr)) + ")");
        }
        int overflow = 0;
        const long long mv = PyLong_AsLongLongAndOverflow(minr.ptr(), &overflow);
        SpfOptions c(o.direction, o.selection);
        c.targets = std::move(t);
        if (res.has_value()) c.residual = std::span<const Amount>(*res);
        // Above INT64_MAX it gates every residual edge exactly as INT64_MAX does (residuals < 2^62).
        c.min_residual = overflow != 0 ? std::numeric_limits<Amount>::max() : static_cast<Amount>(mv);
        c.elig = o.elig.get();
        std::optional<SpfResult> result;
        {
          py::gil_scoped_release release;
          result.emplace(spf(*graph, r, c));
        }
        // The result's Python objects, each through a checked allocation: a failure raises.
        auto [obj, slot] = new_object<PySpfResult>();
        if (result->costs.has_value()) {
          py::array_t<std::int64_t> costs(static_cast<py::ssize_t>(result->costs->size()));
          std::copy(result->costs->begin(), result->costs->end(), costs.mutable_data());
          py::detail::array_proxy(costs.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
          slot->costs = std::move(costs);
        }
        slot->dag = object_from<PredDAG>(std::move(result->dag));
        if (result->target_cost.has_value() && *result->target_cost != kUnreachable) {
          slot->target_cost = py::int_(*result->target_cost);
        }
        return std::move(obj);
      },
      py::arg("graph"), py::arg("roots"), py::arg("opts"),
      "Shortest paths from the roots (one NodeId or a sequence; every root at distance 0), FORWARD or\n"
      "REVERSE (APPROACH 6) -> SpfResult. Without targets every reachable node is settled: costs and the\n"
      "full DAG of tight edges. With targets the search stops at the nearest targets' cost: target_cost\n"
      "and the DAG pruned to the walks between the roots and a nearest target (anycast; for REVERSE, from\n"
      "the target to a root in the edges' orientation), or for multipath = False\n"
      "the canonical single path, the walk minimising (cost, -bottleneck under PREFER_HIGHER_PRIORITY,\n"
      "hops, edge-id sequence). Validated in this order before any work (ValueError): the roots (non-\n"
      "empty, ids in range), the targets (the same), a vertex in both sets, multipath = False without\n"
      "targets, the residual (length, then every entry in [0, capacity_units]), min_residual >= 1, a\n"
      "foreign Eligibility. Excluded roots and targets are skipped; all roots excluded gives an empty\n"
      "result.");

  // ------------------------------------------------------------------ private hooks

#if NETGRAPH_CORE_PRIVATE_HOOKS

  m.def(
      "_max_flows_per_call",
      [](const py::object& value) {
        auto& cap = detail::max_flows_per_call();
        const auto old = cap.load();
        if (!value.is_none()) cap.store(integer_value(value, 0, kMaxFlowsPerCall, "max_flows_per_call"));
        return old;
      },
      py::arg("value") = py::none(),
      "private test hook: the call-wide creation cap (MAX_FLOWS_PER_CALL by default); returns the old value");

  m.def(
      "_sensitivity_fail_edge",
      [](const py::object& edge) {
        auto& hook = detail::sensitivity_fail_edge();
        const auto old = hook.load();
        hook.store(edge.is_none() ? -1 : integer_value(edge, -1, std::numeric_limits<std::int32_t>::max(), "edge"));
        return old;
      },
      py::arg("edge") = py::none(),
      "private test hook: the candidate whose deletion raises a model exception (None: off)");

  m.def(
      "_fail_python_allocation",
      [](const py::object& fn, const py::tuple& args, long k) {
        // Calls fn(*args) with the k-th Python object allocation (from 0) failing, then
        // restores the allocator. Returns (hit, allocations seen, result or None, exception or
        // None). One caller at a time: the allocator is process-wide while the call runs.
        PyMem_GetAllocator(PYMEM_DOMAIN_OBJ, &g_saved_allocator);
        PyMemAllocatorEx hook {nullptr, failing_malloc, failing_calloc, failing_realloc, forwarding_free};
        g_fail_countdown.store(k);
        g_fail_seen.store(0);
        g_fail_hit.store(false);
        PyMem_SetAllocator(PYMEM_DOMAIN_OBJ, &hook);
        PyObject* r = PyObject_CallObject(fn.ptr(), args.ptr());
        PyMem_SetAllocator(PYMEM_DOMAIN_OBJ, &g_saved_allocator);
        g_fail_countdown.store(-1);
        py::object result = r != nullptr ? py::reinterpret_steal<py::object>(r) : py::object(py::none());
        py::object error = py::none();
        if (r == nullptr) error = py::reinterpret_steal<py::object>(PyErr_GetRaisedException());
        return py::make_tuple(g_fail_hit.load(), g_fail_seen.load(), result, error);
      },
      py::arg("fn"), py::arg("args"), py::arg("k"),
      "private test hook: fn(*args) with the k-th Python object allocation failing -> (hit, seen, result, error)");

  m.def(
      "_round_robin",
      [](const py::object& total, const py::object& k) {
        const auto t = integer_value(total, 0, kDomainLimit - 1, "total");
        const auto n = integer_value(k, 1, 1 << 20, "k");
        std::vector<Amount> out(static_cast<std::size_t>(n));
        detail::round_robin(t, out);
        return out;
      },
      py::arg("total"), py::arg("k"), "private: round robin of total units over k members (APPROACH 2.1)");

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
      "private: highest-averages (D'Hondt) apportionment of total units over positive weights, ties by\n"
      "the lower index (APPROACH 2.1)");
#endif
}
