// The routing modes of a policy (CostOnly, CapacityAware, Global, Pinned) and their parameters:
// the cost ceilings Absolute and Stretch, the tunnel count Tunnels, Global's support and Pinned's
// program and black holes.
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

  py::class_<Tunnels>(m, "Tunnels",
                      "The tunnel count of T1 (CapacityAware) and T2, per head: count in [1, 10000] (checked when a\n"
                      "policy or an evaluator uses it). up_to=True (T2 only; T1 rejects it): one common count k in\n"
                      "1..count for every head, the smallest that delivers the design.")
      .def(py::init([](const py::object& count, bool up_to) {
             return Tunnels {*int_field(count, std::numeric_limits<std::int64_t>::min(),
                                        std::numeric_limits<std::int64_t>::max(), "tunnel count", false),
                             up_to};
           }),
           py::arg("count"), py::kw_only(), py::arg("up_to") = false,
           "Tunnels(count, *, up_to=False): a 64-bit integer count (a bool, a float or an integer outside 64 bits\n"
           "raises ValueError here; [1, 10000] is checked where it is used) and a bool")
      .def_readonly("count", &Tunnels::count, "the tunnels per head (the largest count for up_to)")
      .def_readonly("up_to", &Tunnels::up_to, "True: the smallest sufficient count up to `count` (T2 only)")
      .def("__eq__", [](const Tunnels& a, const py::object& b) { return py::isinstance<Tunnels>(b) && a == b.cast<Tunnels>(); })
      .def("__hash__", [](const Tunnels& a) { return py::hash(py::make_tuple("Tunnels", a.count, a.up_to)); })
      .def("__repr__", [](const Tunnels& a) {
        return "Tunnels(" + std::to_string(a.count) + (a.up_to ? ", up_to=True)" : ")");
      });

  py::class_<CostOnly>(m, "CostOnly",
                       "Routing like IGP forwarding (OSPF, IS-IS routers) on each roster source's static shortest-path\n"
                       "DAG (read from the DagCache; the selection is applied at the epoch's first positive call),\n"
                       "fixed for the epoch; one record per roster source. EQUAL_BALANCED and WEIGHTED pool the sources\n"
                       "(origin EQUAL): where their routes meet, a router splits their traffic once. MAX_FLOW is a\n"
                       "controller on the union of those routes (origin OPTIMIZED or EQUAL). An edge of capacity 0\n"
                       "stays in the DAG: EQUAL_BALANCED still sends it its share, which it cannot carry, and WEIGHTED\n"
                       "gives it weight 0.")
      .def(py::init<>(), "CostOnly routing (no parameters)")
      .def("__repr__", [](const CostOnly&) { return std::string("CostOnly()"); });

  py::class_<CapacityAware>(
      m, "CapacityAware",
      "T1: every head sets up its own strict tunnels by CSPF, no controller (like RSVP-TE head ends).\n"
      "The design volume D is split equally over the roster heads and each head's share equally over its\n"
      "tunnels.count tunnels; in round order (round j, heads ascending: tunnel index j*H + i) each tunnel\n"
      "keeps the allowed links whose unreserved capacity (the ledger's, minus the holds of every other\n"
      "owner and of this setup's earlier tunnels) is at least max(request, min_unreserved_bandwidth) (a\n"
      "request of 0 with a minimum of 0 has no filter) and takes the shortest path in the selection's\n"
      "single-path order; it comes up holding exactly its request on every link (UP), or stays down\n"
      "holding nothing: NO_BANDWIDTH (no path passes the filter) or COST_CEILING (the cheapest exceeds\n"
      "max_path_cost: Absolute is inclusive; Stretch(f) is f times the smallest static shortest cost over\n"
      "all roster heads). Traffic enters equally at the heads and splits equally over each head's up\n"
      "tunnels; a head with no tunnel up loses its share. Holds never police traffic.")
      .def(py::init([](const Tunnels& tunnels, const py::object& max_path_cost, const py::object& min_unreserved) {
             CapacityAware ca;
             ca.tunnels = tunnels;
             if (!max_path_cost.is_none()) {
               if (py::isinstance<Absolute>(max_path_cost)) {
                 ca.max_path_cost = max_path_cost.cast<Absolute>();
               } else if (py::isinstance<Stretch>(max_path_cost)) {
                 ca.max_path_cost = max_path_cost.cast<Stretch>();
               } else {
                 throw py::type_error("max_path_cost must be Absolute, Stretch or None");
               }
             }
             ca.min_unreserved_bandwidth = *int_field(min_unreserved, std::numeric_limits<Amount>::min(),
                                                      std::numeric_limits<Amount>::max(), "min_unreserved_bandwidth",
                                                      false);
             return ca;
           }),
           py::arg("tunnels"), py::kw_only(), py::arg("max_path_cost") = py::none(),
           py::arg("min_unreserved_bandwidth") = 0,
           "CapacityAware(tunnels, *, max_path_cost=None, min_unreserved_bandwidth=0): the fields are checked\n"
           "when a policy is built with them (the count in [1, 10000] and exact, the ceiling,\n"
           "min_unreserved_bandwidth in [0, 2^62) units). Here, a min_unreserved_bandwidth of the wrong type or\n"
           "outside 64 bits raises ValueError; tunnels that is not a Tunnels, a ceiling of the wrong class and\n"
           "an unknown keyword raise TypeError.")
      .def_property_readonly(
          "tunnels", [](const CapacityAware& ca) { return typed_object(ca.tunnels); },
          "the tunnel count per head (a Tunnels)")
      .def_property_readonly(
          "max_path_cost",
          [](const CapacityAware& ca) -> py::object {
            if (!ca.max_path_cost.has_value()) return py::none();
            return std::visit([](const auto& x) { return py::cast(x); }, *ca.max_path_cost);
          },
          "the inclusive route cost ceiling (Absolute or Stretch); None: no ceiling")
      .def_readonly("min_unreserved_bandwidth", &CapacityAware::min_unreserved_bandwidth,
                    "units in [0, 2^62): a search minimum only (each attempt keeps the links with at least\n"
                    "max(request, this) unreserved units), never a reservation")
      .def("__repr__", [](const CapacityAware& ca) {
        std::string ceiling = "None";
        if (ca.max_path_cost.has_value()) {
          ceiling = std::visit([](const auto& x) { return std::string(py::repr(typed_object(x))); }, *ca.max_path_cost);
        }
        return "CapacityAware(Tunnels(" + std::to_string(ca.tunnels.count) + (ca.tunnels.up_to ? ", up_to=True)" : ")") +
               ", max_path_cost=" + ceiling + ", min_unreserved_bandwidth=" + std::to_string(ca.min_unreserved_bandwidth) +
               ")";
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
