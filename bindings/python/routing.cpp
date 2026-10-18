// The routing modes of a policy (CostOnly, CapacityAware, Global, Pinned) and their parameters:
// the cost ceilings Absolute and Stretch, the DiminishingReturns guard, Global's support and
// Pinned's program and black holes.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <variant>

#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/program.hpp"
#include "integers.hpp"
#include "publication.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// An optional int field of [lo, hi] (None allowed when `optional`); bools and floats rejected.
std::optional<std::int64_t> int_field(const py::object& obj, std::int64_t lo, std::int64_t hi, const std::string& what,
                                      bool optional) {
  if (obj.is_none()) {
    if (optional) return std::nullopt;
    throw py::value_error(what + " must be an integer (got None)");
  }
  return integer_value(obj, lo, hi, what);
}

// A float field: a Python float, an integer or a NumPy floating scalar; bools and other types
// rejected. The conversion order is Python's float(): a float as it is, then the object's own
// __float__, then its __index__. A Python integer beyond the finite doubles (the int itself, or the
// one its __index__ returned) is outside the field's domain (ValueError); every exception raised by
// the object's own __float__ or __index__, OverflowError included, propagates as raised.
double float_field(const py::object& obj, const std::string& what) {
  if (is_bool(obj) || !(py::isinstance<py::float_>(obj) || PyIndex_Check(obj.ptr()) ||
                        py::isinstance(obj, py::module_::import("numpy").attr("floating")))) {
    throw py::value_error(what + " must be a number (got " + std::string(py::str(py::type::of(obj).attr("__name__"))) +
                          ")");
  }
  if (PyFloat_Check(obj.ptr())) return PyFloat_AS_DOUBLE(obj.ptr());
  const PyNumberMethods* nb = Py_TYPE(obj.ptr())->tp_as_number;
  const bool own_float = nb != nullptr && nb->nb_float != nullptr && nb->nb_float != PyLong_Type.tp_as_number->nb_float;
  if (own_float) {  // the object's own __float__ (a NumPy scalar, a class defining it)
    const double x = PyFloat_AsDouble(obj.ptr());
    if (x == -1.0 && PyErr_Occurred() != nullptr) throw py::error_already_set();
    return x;
  }
  // A Python int (int's own __float__ only narrows it), or an object converted by its __index__.
  py::object idx = py::reinterpret_borrow<py::object>(obj);
  if (!PyLong_Check(obj.ptr())) {
    idx = py::reinterpret_steal<py::object>(PyNumber_Index(obj.ptr()));
    if (!idx) throw py::error_already_set();
  }
  const double x = PyLong_AsDouble(idx.ptr());
  if (x == -1.0 && PyErr_Occurred() != nullptr) {
    if (!PyErr_ExceptionMatches(PyExc_OverflowError)) throw py::error_already_set();
    PyErr_Clear();
    throw py::value_error(what + " must be a number in the float range (got an integer beyond the largest finite float)");
  }
  return x;
}

} // namespace

void register_routing(py::module_& m) {
  py::class_<Absolute>(m, "Absolute", "Cost ceiling: an absolute, inclusive path-cost budget (cost >= 0).")
      .def(py::init([](const py::object& cost) {
             return Absolute {*int_field(cost, std::numeric_limits<Cost>::min(), std::numeric_limits<Cost>::max(),
                                         "Absolute cost", false)};
           }),
           py::arg("cost"), "Absolute(cost): a 64-bit integer budget in cost units (checked >= 0 where it is used)")
      .def_readonly("cost", &Absolute::cost, "the budget: a path qualifies when its cost is at most this")
      .def("__eq__", [](const Absolute& a, const py::object& b) { return py::isinstance<Absolute>(b) && a == b.cast<Absolute>(); })
      .def("__hash__", [](const Absolute& a) { return py::hash(py::make_tuple("Absolute", a.cost)); })
      .def("__repr__", [](const Absolute& a) { return "Absolute(" + std::to_string(a.cost) + ")"; });

  py::class_<Stretch>(m, "Stretch",
                      "Cost ceiling: a route qualifies when its cost is at most factor times the static shortest cost\n"
                      "under the Eligibility; factor is in [1, inf), and the float is compared exactly, without rounding.")
      .def(py::init([](const py::object& factor) { return Stretch {float_field(factor, "Stretch factor")}; }),
           py::arg("factor"),
           "Stretch(factor): a number, checked finite and >= 1 where it is used; a path qualifies when its\n"
           "cost is at most factor times the static shortest cost. Here, an integer beyond the largest finite\n"
           "float raises ValueError.")
      .def_readonly("factor", &Stretch::factor, "the factor as a float")
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
           py::arg("window"), py::arg("fraction"),
           "DiminishingReturns(window, fraction): a 32-bit integer and a number, checked (window >= 1,\n"
           "fraction in (0, 1)) when a policy is built with them. Here, a window outside 32 bits or a fraction\n"
           "that is an integer beyond the largest finite float raises ValueError.")
      .def_readonly("window", &DiminishingReturns::window, "the number of recent accepted growth steps compared")
      .def_readonly("fraction", &DiminishingReturns::fraction,
                    "the share of the source's offer those steps must gain together to continue")
      .def("__repr__", [](const DiminishingReturns& d) {
        return "DiminishingReturns(window=" + std::to_string(d.window) + ", fraction=" + std::string(py::str(py::float_(d.fraction))) + ")";
      });

  py::class_<CostOnly>(m, "CostOnly",
                       "Routing like IGP forwarding (OSPF, IS-IS routers) on each roster source's static shortest-path\n"
                       "DAG (read from the DagCache; the selection is applied at the epoch's first positive call),\n"
                       "never refreshed, grown or reoptimised; one record per roster source. EQUAL_BALANCED and\n"
                       "WEIGHTED pool the sources (origin EQUAL): where their routes meet, a router splits their traffic\n"
                       "once. MAX_FLOW is a controller on the union of those routes (origin OPTIMIZED or EQUAL). An edge\n"
                       "of capacity 0 stays in the DAG: EQUAL_BALANCED still sends it its share, which it cannot\n"
                       "carry, and WEIGHTED gives it weight 0.")
      .def(py::init<>(), "CostOnly routing (no parameters)")
      .def("__repr__", [](const CostOnly&) { return std::string("CostOnly()"); });

  py::class_<CapacityAware>(
      m, "CapacityAware",
      "Routing like distributed CSPF tunnels (each tunnel head picks its own shortest path among edges\n"
      "with enough free capacity). Each source creates flow_count flows, one search at a time, on edges\n"
      "with residual >= max(1, min_flow_bandwidth) units and within max_path_cost; while its offer is\n"
      "unmet it adds flows up to max_flow_count (None: no ceiling). With reoptimize, a flow moves to a\n"
      "new route only when it keeps its whole amount there. `diminishing` (None: off) stops adding flows\n"
      "once the recent ones gain too little.")
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
           py::arg("diminishing") = py::none(),
           "The fields are checked when a policy is built with them (flow_count >= 1, max_flow_count >=\n"
           "flow_count, min_flow_bandwidth in [0, 2^62) units, the ceiling, the guard). Here, an integer field\n"
           "of the wrong type or outside its integer type (32-bit flow_count and max_flow_count, 64-bit\n"
           "min_flow_bandwidth) raises ValueError; a ceiling or guard of the wrong class, a reoptimize that is\n"
           "not a bool, a number or None (a number or None counts by its truth value), and an unknown keyword\n"
           "raise TypeError.")
      .def_readonly("flow_count", &CapacityAware::flow_count,
                    "the flows a source creates on its first positive request")
      .def_readonly("max_flow_count", &CapacityAware::max_flow_count,
                    "the most flows growth may reach per source; None: no ceiling")
      .def_property_readonly(
          "max_path_cost",
          [](const CapacityAware& ca) -> py::object {
            if (!ca.max_path_cost.has_value()) return py::none();
            return std::visit([](const auto& x) { return py::cast(x); }, *ca.max_path_cost);
          },
          "the inclusive route cost ceiling (Absolute or Stretch); None: no ceiling")
      .def_readonly("min_flow_bandwidth", &CapacityAware::min_flow_bandwidth,
                    "units in [0, 2^62): every search keeps only edges with at least max(1, this) residual units")
      .def_readonly("reoptimize", &CapacityAware::reoptimize,
                    "True: a retained flow moves to a newly found route when it keeps its whole amount there")
      .def_readonly("diminishing", &CapacityAware::diminishing, "the growth guard; None: off")
      .def("__repr__", [](const CapacityAware& ca) {
        return "CapacityAware(flow_count=" + std::to_string(ca.flow_count) + ", max_flow_count=" +
               (ca.max_flow_count ? std::to_string(*ca.max_flow_count) : std::string("None")) + ")";
      });

  py::class_<Global>(m, "Global",
                     "Routing: an ideal central controller that sees the whole network. Every positive call re-solves\n"
                     "the minimum-cost flow of what the policy already delivers plus this call's request\n"
                     "(placed_demand() + volume; what an earlier call could not place is not requested again) over\n"
                     "the eligible graph, or exactly over `support` (a PredDAG of the policy's graph, trimmed to the\n"
                     "walks from the sources to the allowed destinations under the Eligibility; every other edge\n"
                     "carries nothing), on the free capacity plus its own allocation (it may move its own traffic,\n"
                     "never another's). It takes no placement and no selection. Origin OPTIMIZED (the default): the\n"
                     "controller chooses how much enters at each source; EQUAL: every roster source sends its\n"
                     "round-robin share of the largest total whose shares all fit. Its program is the solve cut into\n"
                     "paths per source, each a tunnel weighted by its amount. A positive call on a network beyond the\n"
                     "flow core's int32 numbering raises ValueError (see FlowPolicy.place).")
      .def(py::init([](const py::object& support) {
             Global gl;
             if (!support.is_none()) {
               if (!py::isinstance<PredDAG>(support)) throw py::type_error("Global support must be a PredDAG or None");
               gl.support = support.cast<std::shared_ptr<PredDAG>>();
             }
             return gl;
           }),
           py::kw_only(), py::arg("support") = py::none(),
           "Global(*, support=None): support a PredDAG of the policy's graph (checked when a policy is built\n"
           "with it: `support belongs to another graph`) or None for every eligible edge")
      .def_property_readonly(
          "support",
          [](const Global& gl) {
            return optional_shared(gl.support);
          },
          "the supplied support (a PredDAG), or None")
      .def("__repr__", [](const Global& gl) { return std::string(gl.support ? "Global(support=...)" : "Global()"); });

  py::class_<Pinned>(m, "Pinned",
                     "Routing: forwarding over a given program, computing nothing (the program's endpoints must be\n"
                     "the policy's). The policy's Eligibility repairs failures: a tunnel touching an excluded element is\n"
                     "FAILED, a split route loses its excluded members, and the roster is rebuilt from the program's\n"
                     "basis (a head that no longer reaches a destination leaves the equal pool). black_holes: the nodes\n"
                     "and edges this Eligibility excludes are silent failures; the program is forwarded unchanged and\n"
                     "traffic reaching them is dropped there (stale forwarding; None: no black holes). One record per\n"
                     "program entry.")
      .def(py::init([](const std::shared_ptr<ForwardingProgram>& program, const std::shared_ptr<Eligibility>& black_holes) {
             if (!program) throw py::type_error("Pinned program must be a ForwardingProgram");
             Pinned p;
             p.program = program;
             p.black_holes = black_holes;
             return p;
           }),
           py::arg("program"), py::kw_only(), py::arg("black_holes") = nullptr,
           "Pinned(program, *, black_holes=None): a ForwardingProgram and an Eligibility of the same graph (both\n"
           "kept by reference; checked when a policy is built with them)")
      .def_property_readonly(
          "program", [](const Pinned& p) { return typed_shared(p.program); },
          "the forwarded program (shared, immutable)")
      .def_property_readonly(
          "black_holes",
          [](const Pinned& p) {
            return optional_shared(p.black_holes);
          },
          "the black holes (an Eligibility), or None")
      .def("__repr__", [](const Pinned& p) {
        return "Pinned(tunnels=" + std::to_string(p.program->tunnels().size()) + ")";
      });
}

} // namespace netgraph::core::bindings
