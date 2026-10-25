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
    case SourceAllocation::Weighted: return "SourceAllocation.WEIGHTED";
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

// Explicit edge or source weights, converted once and shared by policies built from the model.
// Scan all shape defects, then integer types, then int64 storage overflow (ValueError).
// The policy checks graph-sized length and the range [0, 2^62); only the diagnostic name differs.
std::shared_ptr<const std::vector<Amount>> weights_of(const py::object& weights, const char* name) {
  if (weights.is_none()) return {};
  (void)sequence_length(weights, name);
  const auto w = wide_integers(weights, name);
  for (std::size_t e = 0; e < w.size(); ++e) {
    if (w.sign(e) != 0) throw py::value_error(std::string(name) + "[" + std::to_string(e) + "] must be an integer in [0, 2^62)");
  }
  return std::make_shared<const std::vector<Amount>>(w.narrow<Amount>());
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
  } else if (PyFloat_Check(factor.ptr()) || has_type(factor.ptr(), numpy_number_types().second)) {
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
             return AbsoluteCost {integer_value(cost, std::numeric_limits<Cost>::min(), std::numeric_limits<Cost>::max(),
                                             "AbsoluteCost cost")};
           }),
           py::arg("cost"), "AbsoluteCost(cost): a 64-bit integer budget in cost units (checked >= 0 where it is used)")
      .def_property_readonly("cost", int_reader(&AbsoluteCost::cost), "the budget: a path qualifies when its cost is at most this")
      .def("__eq__", [](const AbsoluteCost& a, const py::object& b) { return py::isinstance<AbsoluteCost>(b) && a == b.cast<AbsoluteCost>(); })
      .def("__hash__", [](const AbsoluteCost& a) { return hash_int(objects_tuple({str_object("AbsoluteCost"), int_object(a.cost)})); })
      .def("__repr__", [](const AbsoluteCost& a) { return "AbsoluteCost(" + std::to_string(a.cost) + ")"; });

  py::class_<CostStretch>(m, "CostStretch",
                      "Cost ceiling: a path qualifies when cost * d <= n * base, exactly, where n/d is the factor as\n"
                      "one reduced fraction and base the static shortest cost under the Eligibility (T1: the smallest\n"
                      "over all roster heads).")
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
             return TunnelCount {integer_value(count, std::numeric_limits<std::int64_t>::min(),
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
                       "DAG (read from the ShortestPathCache; the route choice is applied when the program is prepared:\n"
                       "a positive place of an epoch without one, or any replace or evaluate), fixed for the epoch;\n"
                       "one record per roster source. EQUAL (ECMP, the default) and WEIGHTED\n"
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
             sp.edge_weights = weights_of(edge_weights, "edge_weights");
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
      "tunnel_count.per_head tunnels; in round order (round j, heads ascending: entry_id j*H + i) each tunnel\n"
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
             ca.min_unreserved_bandwidth = integer_value(min_unreserved, std::numeric_limits<Amount>::min(),
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
                     "Routing: an ideal central controller that sees the whole network, or only a support DAG. It takes\n"
                     "no route choice, and its source_allocation is required.\n"
                     "Without tunnel_count: a positive place re-solves the minimum-cost flow of what the policy\n"
                     "already delivers plus this call's request (delivered() + offer; what an earlier call could not\n"
                     "place is not requested again); evaluate and replace solve their whole offer afresh, and forward\n"
                     "sends its offer over the installed paths without a solve. A solve runs over the eligible graph,\n"
                     "or exactly over `support` (a\n"
                     "PredecessorDAG of the policy's graph, trimmed to the walks from the sources to the allowed\n"
                     "destinations under the Eligibility; every other edge carries nothing), on the free capacity\n"
                     "plus its own allocation: it may move its own traffic, never another's. OPTIMIZED: the controller\n"
                     "chooses how much enters at each source. EQUAL: every roster source sends its round-robin share\n"
                     "of the largest total whose shares all fit. WEIGHTED: highest-averages shares of the declared\n"
                     "source_weights; zero weights receive nothing, and all-zero roster weights end NO_ACTIVE_WEIGHT\n"
                     "on a positive live call. The program is the solve cut into paths per source, each a tunnel\n"
                     "weighted by its amount; it keeps the declared weights for weighted sources and the solved\n"
                     "supplies for optimized ones. A positive call on a network beyond the flow core's int32\n"
                     "numbering raises ValueError (see FlowPolicy.place).\n"
                     "With tunnel_count=TunnelCount(N) (T2; with support, every step uses only the support's edges\n"
                     "allowed by the policy's mask): a controller (like a PCE) lays out, for every live head (one\n"
                     "reaching a destination over links with free capacity >= 1; the other heads are dead and lose\n"
                     "their share), N single-path tunnels with one common whole-unit planned bandwidth: first as\n"
                     "large as possible (U), then for a finite design no larger than needed (u), then the cheapest\n"
                     "routes, then with prefer_headroom=True the routes that leave the most spare capacity, then one\n"
                     "fixed unit-by-unit tie. TunnelCount(L, up_to=True): the smallest common count k <= L whose\n"
                     "program carries the live heads' design shares (else the largest live-share score, then the\n"
                     "smallest count); no promise of whole-demand lossless delivery. Traffic splits equally over the\n"
                     "heads and over each head's tunnels; planned bandwidth is never held.")
      .def(py::init([](SourceAllocation source_allocation, const py::object& source_weights, const py::object& support, const py::object& tunnels,
                       StrictBool prefer_headroom) {
             CentralizedTE gl(source_allocation);
             gl.source_weights = weights_of(source_weights, "source_weights");
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
           py::kw_only(), py::arg("source_allocation"), py::arg("source_weights") = py::none(), py::arg("support") = py::none(),
           py::arg("tunnel_count") = py::none(), py::arg("prefer_headroom") = StrictBool {false},
           "CentralizedTE(*, source_allocation, source_weights=None, support=None, tunnel_count=None, prefer_headroom=False):\n"
           "source_allocation is required (EQUAL, OPTIMIZED or WEIGHTED; None or a missing argument raises TypeError).\n"
           "WEIGHTED requires source_weights and no tunnel_count: one integer per graph node in [0, 2^62),\n"
           "copied once and exposed as an immutable tuple. Other source allocations refuse source_weights.\n"
           "support a PredecessorDAG of the policy's graph (checked when a policy is built with it: `support\n"
           "belongs to another graph`) or None for every eligible edge; tunnel_count a TunnelCount (T2;\n"
           "the count in [1, 10000] is checked when a policy is built with it) or None; prefer_headroom (T2\n"
           "only; True without tunnel_count raises `prefer_headroom requires tunnel_count` when a policy is\n"
           "built with it): the headroom penalty as T2's objective after route cost")
      .def_property_readonly("source_allocation", field_object(&CentralizedTE::source_allocation), "the SourceAllocation: EQUAL, OPTIMIZED or WEIGHTED")
      .def_property_readonly("source_weights", [](const CentralizedTE& gl) -> py::object {
        if (!gl.source_weights) return py::none();
        return amounts_tuple(*gl.source_weights);
      }, "declared node-indexed source weights as an immutable tuple, or None")
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
        if (gl.source_weights) add("source_weights=...");
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
