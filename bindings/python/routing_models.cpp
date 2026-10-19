// The routing modes of a policy (ShortestPaths, DistributedTE, CentralizedTE, Pinned) and their parameters:
// the cost ceilings AbsoluteCost and CostStretch, the tunnel count TunnelCount, CentralizedTE's support and Pinned's
// program and black holes.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <variant>

#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/forwarding_program.hpp"
#include "integers.hpp"
#include "strict_bool.hpp"
#include "publication.hpp"
#include "registration.hpp"
#include "values.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// A RouteTieBreak as its repr spells it, without a Python object (a repr must not fail on its own
// allocations of an enum object).
const char* tie_break_repr(RouteTieBreak t) {
  switch (t) {
    case RouteTieBreak::FewerHops: return "RouteTieBreak.FEWER_HOPS";
    case RouteTieBreak::HigherCapacity: return "RouteTieBreak.HIGHER_CAPACITY";
    case RouteTieBreak::HigherAvailableCapacity: return "RouteTieBreak.HIGHER_AVAILABLE_CAPACITY";
    case RouteTieBreak::HigherHeadroom: return "RouteTieBreak.HIGHER_HEADROOM";
  }
  return "RouteTieBreak(?)";
}

// The other enumerations of the models as their reprs spell them (C++ tables, as tie_break_repr).
const char* allocation_repr(RouteAllocation a) {
  switch (a) {
    case RouteAllocation::Equal: return "RouteAllocation.EQUAL";
    case RouteAllocation::Weighted: return "RouteAllocation.WEIGHTED";
    case RouteAllocation::Optimized: return "RouteAllocation.OPTIMIZED";
  }
  return "RouteAllocation(?)";
}
const char* source_allocation_repr(SourceAllocation a) {
  switch (a) {
    case SourceAllocation::Equal: return "SourceAllocation.EQUAL";
    case SourceAllocation::Optimized: return "SourceAllocation.OPTIMIZED";
  }
  return "SourceAllocation(?)";
}
const char* weight_basis_repr(WeightBasis b) {
  switch (b) {
    case WeightBasis::Capacity: return "WeightBasis.CAPACITY";
    case WeightBasis::Explicit: return "WeightBasis.EXPLICIT";
    case WeightBasis::ResidualSnapshot: return "WeightBasis.RESIDUAL_SNAPSHOT";
  }
  return "WeightBasis(?)";
}

// Explicit WCMP weights, one integer per edge, converted once and shared by every policy built from
// the model: a length or type defect raises here; an entry beyond int64 is outside the domain
// (ValueError); the policy checks the length against its graph and the range [0, 2^62).
std::shared_ptr<const std::vector<Amount>> edge_weights_of(const py::object& edge_weights) {
  if (edge_weights.is_none()) return {};
  (void)sequence_length(edge_weights, "edge_weights", IntegerArgument::Configuration);
  const auto w = wide_integers(edge_weights, "edge_weights", IntegerArgument::Configuration);
  for (std::size_t e = 0; e < w.size(); ++e) {
    if (w.sign(e) != 0) throw py::value_error("edge_weights[" + std::to_string(e) + "] must be an integer in [0, 2^62)");
  }
  return std::make_shared<const std::vector<Amount>>(w.narrow<Amount>());
}

// An integer field of a configuration value, in [lo, hi]: a bool, a float, None or another type is
// the configuration's wrong-type error (integers.hpp); an integer outside [lo, hi] is a ValueError.
std::int64_t int_field(const py::object& obj, std::int64_t lo, std::int64_t hi, const std::string& what) {
  return integer_value(obj, lo, hi, what, IntegerArgument::Configuration);
}

// A CostStretch factor as its one exact fraction, converted in this order: a bool (Python or
// numpy.bool_) raises TypeError; a fractions.Fraction is taken as it is; a Python float or a NumPy
// floating value is Fraction(str(f)), the shortest decimal its own type prints (a NaN or an infinity
// raises ValueError); an integer (an object with __index__) is n/1; anything else raises TypeError.
// Then the reduced components: a numerator in [-2^63, 2^63), a denominator in [1, 2^63).
CostStretch cost_stretch_of(const py::object& factor) {
  const auto type_error = [&] {
    return py::type_error("CostStretch factor must be a Fraction, a float or an integer (got " +
                          std::string(py::str(py::type::of(factor).attr("__name__"))) + ")");
  };
  if (is_bool(factor)) throw type_error();
  const py::object fraction_type = py::module_::import("fractions").attr("Fraction");
  py::object fraction;
  if (py::isinstance(factor, fraction_type)) {
    fraction = factor;
  } else if (PyFloat_Check(factor.ptr()) || py::isinstance(factor, py::module_::import("numpy").attr("floating"))) {
    const double x = PyFloat_AsDouble(factor.ptr());
    if (x == -1.0 && PyErr_Occurred() != nullptr) throw py::error_already_set();
    if (!std::isfinite(x)) throw py::value_error("CostStretch factor must be finite (got " + std::string(py::str(factor)) + ")");
    fraction = fraction_type(py::str(factor));
  } else if (PyIndex_Check(factor.ptr())) {
    const auto n = py::reinterpret_steal<py::object>(PyNumber_Index(factor.ptr()));
    if (!n) throw py::error_already_set();
    fraction = fraction_type(n);
  } else {
    throw type_error();
  }
  const auto component = [&](const char* name, std::int64_t lo, const char* text) {
    const py::object v = fraction.attr(name);
    int overflow = 0;
    const long long x = PyLong_AsLongLongAndOverflow(v.ptr(), &overflow);
    if (x == -1 && PyErr_Occurred() != nullptr) throw py::error_already_set();
    if (overflow != 0 || x < lo) throw py::value_error(text);
    return static_cast<std::int64_t>(x);
  };
  const auto n = component("numerator", std::numeric_limits<std::int64_t>::min(),
                           "CostStretch fraction numerator must be in [-2^63, 2^63)");
  const auto d = component("denominator", 1, "CostStretch fraction denominator must be in [1, 2^63)");
  return CostStretch(n, d);
}

} // namespace

void register_routing(py::module_& m) {
  py::class_<AbsoluteCost>(m, "AbsoluteCost", "Cost ceiling: an absolute, inclusive path-cost budget (cost >= 0).")
      .def(py::init([](const py::object& cost) {
             return AbsoluteCost {int_field(cost, std::numeric_limits<Cost>::min(), std::numeric_limits<Cost>::max(),
                                         "AbsoluteCost cost")};
           }),
           py::arg("cost"), "AbsoluteCost(cost): a 64-bit integer budget in cost units (checked >= 0 where it is used)")
      .def_property_readonly("cost", int_reader(&AbsoluteCost::cost), "the budget: a path qualifies when its cost is at most this")
      .def("__eq__", [](const AbsoluteCost& a, const py::object& b) { return py::isinstance<AbsoluteCost>(b) && a == b.cast<AbsoluteCost>(); })
      .def("__hash__", [](const AbsoluteCost& a) { return hash_int(objects_tuple({str_object("AbsoluteCost"), int_object(a.cost)})); })
      .def("__repr__", [](const AbsoluteCost& a) { return "AbsoluteCost(" + std::to_string(a.cost) + ")"; });

  py::class_<CostStretch>(m, "CostStretch",
                      "Cost ceiling: a path qualifies when cost * d <= n * base, exactly, where n/d is the factor as\n"
                      "one reduced fraction and base the static shortest cost under the Eligibility.")
      .def(py::init([](const py::object& factor) { return cost_stretch_of(factor); }), py::arg("factor"),
           "CostStretch(factor): in this order, a bool (or numpy.bool_) raises TypeError; a fractions.Fraction\n"
           "is taken as it is; a Python float or any NumPy floating value is Fraction(str(f)), the shortest\n"
           "decimal its own type prints (1.2, np.float64(1.2) and np.float32(1.2) are all 6/5), and a NaN or\n"
           "an infinity raises ValueError; an integer, NumPy integer or any __index__ object is n/1 exactly;\n"
           "anything else raises TypeError. Then a reduced numerator outside [-2^63, 2^63) or denominator\n"
           "outside [1, 2^63) raises ValueError. A factor below 1 is valid here and refused where the\n"
           "ceiling is used.")
      .def_property_readonly(
          "factor",
          [](const CostStretch& a) {
            return fraction_object(py::module_::import("fractions").attr("Fraction"), a.numerator(), a.denominator());
          },
          "the factor as a fractions.Fraction (built on each access)")
      .def("__eq__", [](const CostStretch& a, const py::object& b) { return py::isinstance<CostStretch>(b) && a == b.cast<CostStretch>(); })
      .def("__hash__",
           [](const CostStretch& a) {
             return hash_int(fraction_object(py::module_::import("fractions").attr("Fraction"), a.numerator(), a.denominator()));
           })
      .def("__repr__", [](const CostStretch& a) {
        return "CostStretch(Fraction(" + std::to_string(a.numerator()) + ", " + std::to_string(a.denominator()) + "))";
      });

  py::class_<TunnelCount>(m, "TunnelCount",
                      "The tunnel count of T1 (DistributedTE) and T2, per head: per_head in [1, 10000] (checked when a\n"
                      "policy or an evaluator uses it). up_to=True (T2 only; T1 rejects it): one common count k in\n"
                      "1..per_head for every head: the smallest common count whose program carries the live heads'\n"
                      "design shares; if none does, the largest live-share score, then the smallest count. It does not\n"
                      "promise whole-demand lossless delivery.")
      .def(py::init([](const py::object& count, StrictBool up_to) {
             return TunnelCount {int_field(count, std::numeric_limits<std::int64_t>::min(),
                                        std::numeric_limits<std::int64_t>::max(), "TunnelCount per_head"),
                             up_to};
           }),
           py::arg("per_head"), py::kw_only(), py::arg("up_to") = StrictBool {false},
           "TunnelCount(per_head, *, up_to=False): a 64-bit integer per_head (a bool, a float or None raises\n"
           "TypeError and an integer outside 64 bits ValueError here; [1, 10000] is checked where it is\n"
           "used) and a bool")
      .def_property_readonly("per_head", int_reader(&TunnelCount::per_head), "the tunnels per head (the largest count for up_to)")
      .def_readonly("up_to", &TunnelCount::up_to, "True: the smallest count up to per_head whose program carries the live heads' design shares (else the\n"
                    "best live-share score, then the smallest); T2 only")
      .def("__eq__", [](const TunnelCount& a, const py::object& b) { return py::isinstance<TunnelCount>(b) && a == b.cast<TunnelCount>(); })
      .def("__hash__", [](const TunnelCount& a) {
        return hash_int(objects_tuple({str_object("TunnelCount"), int_object(a.per_head), py::bool_(a.up_to)}));
      })
      .def("__repr__", [](const TunnelCount& a) {
        return "TunnelCount(" + std::to_string(a.per_head) + (a.up_to ? ", up_to=True)" : ")");
      });

  py::class_<ShortestPaths>(m, "ShortestPaths",
                       "Routing like IGP forwarding (OSPF, IS-IS routers) on each roster source's static shortest-path\n"
                       "DAG (read from the ShortestPathCache; the route choice is applied at the epoch's first positive call),\n"
                       "fixed for the epoch; one record per roster source. EQUAL (ECMP, the default) and WEIGHTED\n"
                       "(WCMP) pool the sources (source allocation EQUAL): where their routes meet, a router splits\n"
                       "their traffic once. OPTIMIZED is a controller on the union of those routes (source allocation\n"
                       "OPTIMIZED by default, or EQUAL). An edge of capacity 0 stays in the DAG: EQUAL still sends it\n"
                       "its share, which it cannot carry, and WEIGHTED gives it weight 0.")
      .def(py::init([](RouteAllocation allocation, std::optional<SourceAllocation> source_allocation,
                       std::optional<WeightBasis> weight_basis, const py::object& edge_weights, RouteTieBreak tie_break,
                       StrictBool multipath, StrictBool parallel_edges) {
             ShortestPaths sp;
             sp.allocation = allocation;
             sp.source_allocation = source_allocation;
             sp.weight_basis = weight_basis;
             sp.edge_weights = edge_weights_of(edge_weights);
             sp.tie_break = tie_break;
             sp.multipath = multipath;
             sp.parallel_edges = parallel_edges;
             return sp;
           }),
           py::kw_only(), py::arg_v("allocation", RouteAllocation::Equal, "RouteAllocation.EQUAL"),
           py::arg("source_allocation") = py::none(), py::arg("weight_basis") = py::none(),
           py::arg("edge_weights") = py::none(), py::arg_v("tie_break", RouteTieBreak::FewerHops, "RouteTieBreak.FEWER_HOPS"),
           py::arg("multipath") = StrictBool {true}, py::arg("parallel_edges") = StrictBool {true},
           "ShortestPaths(*, allocation=EQUAL, source_allocation=None, weight_basis=None, edge_weights=None,\n"
           "tie_break=FEWER_HOPS, multipath=True, parallel_edges=True): ShortestPaths() is ECMP.\n"
           "source_allocation None resolves to EQUAL for EQUAL and WEIGHTED (the only allowed value) and to\n"
           "OPTIMIZED for OPTIMIZED; weight_basis (WEIGHTED only; None: CAPACITY) and edge_weights (one\n"
           "integer per edge, given exactly when weight_basis=EXPLICIT; one array shared by every policy\n"
           "built from this model). The combination is checked when a policy or an evaluator uses the model\n"
           "(ValueError, in the order of the model matrix). Here, a value of the wrong type (a multipath or parallel_edges that\n"
           "is not a bool, an enumeration of another class, an edge_weights that is not a 1-D sequence of\n"
           "integers or has a bool or a float entry) raises TypeError, and an edge_weights entry beyond\n"
           "int64 ValueError.")
      .def_property_readonly("allocation", field_object(&ShortestPaths::allocation), "the RouteAllocation: EQUAL, WEIGHTED or OPTIMIZED")
      .def_property_readonly("source_allocation", optional_field(&ShortestPaths::source_allocation),
                    "the SourceAllocation, or None (a policy's config: resolved to the allocation's default)")
      .def_property_readonly("weight_basis", optional_field(&ShortestPaths::weight_basis),
                    "the WeightBasis, or None (a policy's config: CAPACITY resolved for WEIGHTED)")
      .def_property_readonly(
          "edge_weights",
          [](const ShortestPaths& sp) -> py::object {
            if (!sp.edge_weights) return py::none();
            return amounts_tuple(*sp.edge_weights);
          },
          "the explicit weights (one int per edge, a tuple), or None")
      .def_property_readonly("tie_break", field_object(&ShortestPaths::tie_break), "the RouteTieBreak of equal-cost choices")
      .def_readonly("multipath", &ShortestPaths::multipath,
                    "True: every equal-cost route; False: one canonical path per roster source toward the destination set")
      .def_readonly("parallel_edges", &ShortestPaths::parallel_edges,
                    "True: every eligible parallel edge; False: one per (u, v) pair by the tie rule")
      .def("__repr__", [](const ShortestPaths& sp) {
        return "ShortestPaths(allocation=" + std::string(allocation_repr(sp.allocation)) + ", source_allocation=" +
               (sp.source_allocation ? source_allocation_repr(*sp.source_allocation) : "None") +
               ", weight_basis=" + (sp.weight_basis ? weight_basis_repr(*sp.weight_basis) : "None") +
               ", edge_weights=" + (sp.edge_weights ? "(...)" : "None") +
               ", tie_break=" + tie_break_repr(sp.tie_break) + ", multipath=" + (sp.multipath ? "True" : "False") +
               ", parallel_edges=" + (sp.parallel_edges ? "True" : "False") + ")";
      });

  py::class_<DistributedTE>(
      m, "DistributedTE",
      "T1: every head sets up its own strict tunnels by CSPF, no controller (like RSVP-TE head ends).\n"
      "The design volume D is split equally over the roster heads and each head's share equally over its\n"
      "tunnel_count.per_head tunnels; in round order (round j, heads ascending: tunnel index j*H + i) each tunnel\n"
      "keeps the allowed links whose unreserved capacity (the ledger's, minus the holds of every other\n"
      "owner and of this setup's earlier tunnels) is at least max(request, min_unreserved_bandwidth) (a\n"
      "request of 0 with a minimum of 0 has no filter) and takes the shortest path in the tie_break's\n"
      "single-path order; it comes up holding exactly its request on every link (UP), or stays down\n"
      "holding nothing: INSUFFICIENT_BANDWIDTH (no path passes the filter) or COST_CEILING_EXCEEDED (the cheapest exceeds\n"
      "cost_ceiling: AbsoluteCost is inclusive; CostStretch(f) is f times the smallest static shortest cost over\n"
      "all roster heads). Traffic enters equally at the heads and splits equally over each head's up\n"
      "tunnels; a head with no tunnel up loses its share. Holds never police traffic.")
      .def(py::init([](const TunnelCount& tunnels, const py::object& cost_ceiling, const py::object& min_unreserved,
                       RouteTieBreak tie_break, StrictBool parallel_edges) {
             DistributedTE ca;
             ca.tunnel_count = tunnels;
             ca.tie_break = tie_break;
             ca.parallel_edges = parallel_edges;
             ca.cost_ceiling = cost_ceiling_of(cost_ceiling);
             ca.min_unreserved_bandwidth = int_field(min_unreserved, std::numeric_limits<Amount>::min(),
                                                      std::numeric_limits<Amount>::max(), "min_unreserved_bandwidth");
             return ca;
           }),
           py::arg("tunnel_count"), py::kw_only(), py::arg("cost_ceiling") = py::none(),
           py::arg("min_unreserved_bandwidth") = 0,
           py::arg_v("tie_break", RouteTieBreak::FewerHops, "RouteTieBreak.FEWER_HOPS"),
           py::arg("parallel_edges") = StrictBool {true},
           "DistributedTE(tunnel_count, *, cost_ceiling=None, min_unreserved_bandwidth=0, tie_break=FEWER_HOPS,\n"
           "parallel_edges=True): the fields are checked when a policy is built with them (the count in\n"
           "[1, 10000] and exact, the ceiling, min_unreserved_bandwidth in [0, 2^62) units). tie_break and\n"
           "parallel_edges are the CSPF choice, scored on unreserved bandwidth. Here, a\n"
           "min_unreserved_bandwidth outside 64 bits raises ValueError; one of the wrong type, a tunnel_count\n"
           "that is not a TunnelCount, a ceiling of the wrong class, a parallel_edges that is not a bool and an\n"
           "unknown keyword raise TypeError.")
      .def_property_readonly(
          "tunnel_count", [](const DistributedTE& ca) { return typed_object(ca.tunnel_count); },
          "the tunnel count per head (a TunnelCount)")
      .def_property_readonly(
          "cost_ceiling",
          [](const DistributedTE& ca) -> py::object {
            if (!ca.cost_ceiling.has_value()) return py::none();
            return std::visit([](const auto& x) { return object_of(x); }, *ca.cost_ceiling);
          },
          "the inclusive route cost ceiling (AbsoluteCost or CostStretch); None: no ceiling")
      .def_property_readonly("min_unreserved_bandwidth", int_reader(&DistributedTE::min_unreserved_bandwidth),
                    "units in [0, 2^62): a search minimum only (each attempt keeps the links with at least\n"
                    "max(request, this) unreserved units), never a reservation")
      .def_property_readonly("tie_break", field_object(&DistributedTE::tie_break), "the RouteTieBreak of the CSPF choice")
      .def_readonly("parallel_edges", &DistributedTE::parallel_edges,
                    "True: every eligible parallel edge; False: one per (u, v) pair by the tie rule")
      .def("__repr__", [](const DistributedTE& ca) {
        std::string ceiling = "None";
        if (ca.cost_ceiling.has_value()) {
          ceiling = std::visit([](const auto& x) { return std::string(py::repr(typed_object(x))); }, *ca.cost_ceiling);
        }
        return "DistributedTE(TunnelCount(" + std::to_string(ca.tunnel_count.per_head) + (ca.tunnel_count.up_to ? ", up_to=True)" : ")") +
               ", cost_ceiling=" + ceiling + ", min_unreserved_bandwidth=" + std::to_string(ca.min_unreserved_bandwidth) +
               ", tie_break=" + std::string(tie_break_repr(ca.tie_break)) +
               ", parallel_edges=" + (ca.parallel_edges ? "True" : "False") + ")";
      });

  py::class_<CentralizedTE>(m, "CentralizedTE",
                     "Routing: an ideal central controller that sees the whole network (with `tunnel_count`: T2, below).\n"
                     "Every positive call re-solves the minimum-cost flow of what the policy already delivers plus this\n"
                     "call's request (delivered() + offer; what an earlier call could not place is not requested\n"
                     "again) over the eligible graph, or exactly over `support` (a PredecessorDAG of the policy's graph, trimmed to the\n"
                     "walks from the sources to the allowed destinations under the Eligibility; every other edge\n"
                     "carries nothing), on the free capacity plus its own allocation (it may move its own traffic,\n"
                     "never another's). It takes no route choice. Its source_allocation is required: OPTIMIZED, the\n"
                     "controller chooses how much enters at each source; EQUAL: every roster source sends its\n"
                     "round-robin share of the largest total whose shares all fit. Its program is the solve cut into\n"
                     "paths per source, each a tunnel weighted by its amount. A positive call on a network beyond the\n"
                     "flow core's int32 numbering raises ValueError (see FlowPolicy.place).\n"
                     "T2, CentralizedTE(source_allocation=..., tunnel_count=TunnelCount(N)) (no route choice either; with\n"
                     "support, every step uses only the support's edges allowed by the policy's mask): a controller (like\n"
                     "a PCE) lays out, for every live head (one reaching a destination over links with free capacity\n"
                     ">= 1; the other heads are dead and lose their share), N single-path tunnels with one common\n"
                     "whole-unit planned bandwidth: first as large as possible (U), then for a finite design no larger\n"
                     "than needed (u), then the cheapest routes, then with prefer_headroom=True the routes that leave\n"
                     "the most spare capacity, then one fixed unit-by-unit tie. TunnelCount(L, up_to=True): the smallest\n"
                     "common count k <= L whose program carries the live heads' design shares (else the largest\n"
                     "live-share score, then the smallest count); no promise of whole-demand lossless delivery.\n"
                     "Traffic splits equally over the heads and over each head's tunnels; planned bandwidth is never\n"
                     "held.")
      .def(py::init([](SourceAllocation source_allocation, const py::object& support, const py::object& tunnels,
                       StrictBool prefer_headroom) {
             CentralizedTE gl(source_allocation);
             gl.prefer_headroom = prefer_headroom;
             if (!support.is_none()) {
               if (!py::isinstance<PredecessorDAG>(support)) throw py::type_error("CentralizedTE support must be a PredecessorDAG or None");
               gl.support = support.cast<std::shared_ptr<PredecessorDAG>>();
             }
             if (!tunnels.is_none()) {
               if (!py::isinstance<TunnelCount>(tunnels)) throw py::type_error("CentralizedTE tunnel_count must be a TunnelCount or None");
               gl.tunnel_count = tunnels.cast<TunnelCount>();
             }
             return gl;
           }),
           py::kw_only(), py::arg("source_allocation"), py::arg("support") = py::none(),
           py::arg("tunnel_count") = py::none(), py::arg("prefer_headroom") = StrictBool {false},
           "CentralizedTE(*, source_allocation, support=None, tunnel_count=None, prefer_headroom=False):\n"
           "source_allocation is required (EQUAL or OPTIMIZED; None or a missing argument raises TypeError);\n"
           "support a PredecessorDAG of the policy's graph (checked when a policy is built with it: `support\n"
           "belongs to another graph`) or None for every eligible edge; tunnel_count a TunnelCount (T2;\n"
           "the count in [1, 10000] is checked when a policy is built with it) or None; prefer_headroom (T2\n"
           "only; True without tunnel_count raises `prefer_headroom requires tunnel_count` when a policy is\n"
           "built with it): the headroom penalty as T2's objective after route cost")
      .def_property_readonly("source_allocation", field_object(&CentralizedTE::source_allocation), "the SourceAllocation: EQUAL or OPTIMIZED")
      .def_property_readonly(
          "support",
          [](const CentralizedTE& gl) {
            return optional_shared(gl.support);
          },
          "the supplied support (a PredecessorDAG), or None")
      .def_property_readonly(
          "tunnel_count", [](const CentralizedTE& gl) { return optional_of(gl.tunnel_count); },
          "T2's tunnel count per head (a TunnelCount), or None")
      .def_readonly("prefer_headroom", &CentralizedTE::prefer_headroom,
                    "True: T2 minimizes the headroom penalty after route cost (T2 only)")
      .def("__repr__", [](const CentralizedTE& gl) {
        // The source allocation and every other field that is set.
        std::string out;
        const auto add = [&](const std::string& part) { out += (out.empty() ? "" : ", ") + part; };
        add("source_allocation=" + std::string(source_allocation_repr(gl.source_allocation)));
        if (gl.support) add("support=...");
        if (gl.tunnel_count.has_value()) {
          add("tunnel_count=TunnelCount(" + std::to_string(gl.tunnel_count->per_head) + (gl.tunnel_count->up_to ? ", up_to=True)" : ")"));
        }
        if (gl.prefer_headroom) add("prefer_headroom=True");
        return "CentralizedTE(" + out + ")";
      });

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
