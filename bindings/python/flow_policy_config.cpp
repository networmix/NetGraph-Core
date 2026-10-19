// FlowPolicyConfig and OutputOptions: the configuration of a placement.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/options.hpp"
#include "integers.hpp"
#include "strict_bool.hpp"
#include "publication.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

Routing routing_of(const py::object& obj) {
  if (py::isinstance<ShortestPaths>(obj)) return obj.cast<ShortestPaths>();
  if (py::isinstance<DistributedTE>(obj)) return obj.cast<DistributedTE>();
  if (py::isinstance<CentralizedTE>(obj)) return obj.cast<CentralizedTE>();
  if (py::isinstance<Pinned>(obj)) return obj.cast<Pinned>();
  throw py::type_error("routing must be ShortestPaths, DistributedTE, CentralizedTE or Pinned (got " +
                       std::string(py::str(py::type::of(obj).attr("__name__"))) + ")");
}

// The routing as a new Python object through the checked allocation (publication.hpp): a failed
// allocation raises instead of crashing.
py::object routing_obj(const Routing& r) {
  return std::visit([](const auto& x) { return object_of(x); }, r);
}

} // namespace

void register_policy_config(py::module_& m) {
  py::class_<FlowPolicyConfig>(
      m, "FlowPolicyConfig",
      "One demand's settings: the routing model, which holds every option it reads (ShortestPaths,\n"
      "DistributedTE (T1), CentralizedTE, CentralizedTE(tunnel_count=...) (T2) or Pinned), and the admission\n"
      "(LOSSLESS by default; LOSSY for ShortestPaths EQUAL and WEIGHTED, T1, T2 and Pinned). The\n"
      "combination is checked when a policy or an evaluator uses the configuration (ValueError), in\n"
      "this order: the enumeration values; the model matrix (route and source allocation, admission,\n"
      "the ineffective tie rule, the pooled single-path rule, prefer_headroom, the weight combinations,\n"
      "and last a ShortestPathCache on a model that does not read it); the tunnel count, an up_to count\n"
      "for DistributedTE, the ceilings, min_unreserved_bandwidth and edge_weights' length and entries;\n"
      "the tunnel limit (roster heads times the tunnel count at most 2^30 - 1); a Pinned program's\n"
      "endpoints.")
      .def(py::init([](const py::object& routing, Admission admission) {
             return FlowPolicyConfig {routing_of(routing), admission};
           }),
           py::kw_only(), py::arg("routing"), py::arg_v("admission", Admission::Lossless, "Admission.LOSSLESS"))
      .def_property_readonly(
          "routing", [](const FlowPolicyConfig& c) { return routing_obj(c.routing); },
          "the routing model: ShortestPaths, DistributedTE, CentralizedTE or Pinned")
      .def_property_readonly("admission", field_object(&FlowPolicyConfig::admission), "LOSSLESS or LOSSY")
      .def("__repr__", [](const FlowPolicyConfig& c) {
        return "FlowPolicyConfig(routing=" + std::string(py::repr(routing_obj(c.routing))) + ", admission=" +
               (c.admission == Admission::Lossy ? "Admission.LOSSY" : "Admission.LOSSLESS") + ")";
      });

  py::class_<OutputOptions>(m, "OutputOptions",
                            "Optional outputs of a placement (all false by default; absent outputs are None):\n"
                            "edge_flows (the policy's own allocation per edge after the call), residual and\n"
                            "saturated_boundary_edges (the ledger's after the call), edge_reservations (the policy's own\n"
                            "strict holds per edge after the call; zeros for a model without holds), and min_cut (CentralizedTE with\n"
                            "the optimized source allocation, when its solve reached the maximum flow; None for every other\n"
                            "model, source allocation or call).")
      .def(py::init([](StrictBool edge_flows, StrictBool residual, StrictBool saturated_boundary_edges, StrictBool edge_reservations,
                       StrictBool min_cut) {
             OutputOptions o;
             o.edge_flows = edge_flows;
             o.residual = residual;
             o.saturated_boundary_edges = saturated_boundary_edges;
             o.edge_reservations = edge_reservations;
             o.min_cut = min_cut;
             return o;
           }),
           py::kw_only(), py::arg("edge_flows") = StrictBool {false}, py::arg("residual") = StrictBool {false},
           py::arg("saturated_boundary_edges") = StrictBool {false}, py::arg("edge_reservations") = StrictBool {false},
           py::arg("min_cut") = StrictBool {false})
      .def_readonly("edge_flows", &OutputOptions::edge_flows, "request PlacementResult.edge_flows")
      .def_readonly("residual", &OutputOptions::residual, "request PlacementResult.residual")
      .def_readonly("saturated_boundary_edges", &OutputOptions::saturated_boundary_edges, "request PlacementResult.saturated_boundary_edges")
      .def_readonly("edge_reservations", &OutputOptions::edge_reservations,
                    "request PlacementResult.edge_reservations")
      .def_readonly("min_cut", &OutputOptions::min_cut, "request PlacementResult.min_cut");
}

} // namespace netgraph::core::bindings
