// FlowPolicyConfig and OutputOptions: the configuration of a placement.
#include <pybind11/pybind11.h>

#include <string>
#include <utility>
#include <variant>

#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/options.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

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

} // namespace

void register_policy_config(py::module_& m) {
  py::class_<FlowPolicyConfig>(
      m, "FlowPolicyConfig",
      "The four explicit axes (independent settings) of a policy: routing (CostOnly, CapacityAware,\n"
      "Global or Pinned), placement (FlowPlacement), admission (Admission) and selection\n"
      "(PathSelection). A placement or admission outside its enumeration raises ValueError.")
      .def(py::init([](const py::object& routing, FlowPlacement placement, Admission admission,
                       const PathSelection& selection) {
             // The routing's structure first (TypeError), then the constructor's enumeration checks.
             auto r = routing_of(routing);
             return FlowPolicyConfig {std::move(r), placement, admission, selection};
           }),
           py::kw_only(), py::arg("routing"), py::arg("placement"), py::arg("admission"), py::arg("selection"))
      .def_property_readonly(
          "routing", [](const FlowPolicyConfig& c) { return routing_obj(c.routing); },
          "the routing mode: CostOnly, CapacityAware, Global or Pinned")
      .def_readonly("placement", &FlowPolicyConfig::placement, "how every flow spreads over its route (FlowPlacement)")
      .def_readonly("admission", &FlowPolicyConfig::admission, "LOSSLESS or LOSSY")
      .def_readonly("selection", &FlowPolicyConfig::selection, "the PathSelection the routing's searches read")
      .def("__repr__", [](const FlowPolicyConfig& c) {
        return "FlowPolicyConfig(routing=" + std::string(py::repr(routing_obj(c.routing))) + ")";
      });

  py::class_<OutputOptions>(m, "OutputOptions",
                            "Optional outputs of a placement (all false by default; absent outputs are None):\n"
                            "edge_flows (the policy's own allocation per edge after the call), residual and\n"
                            "bottleneck (the ledger's after the call); on Global, bottleneck adds min_cut.")
      .def(py::init<bool, bool, bool>(), py::kw_only(), py::arg("edge_flows") = false, py::arg("residual") = false,
           py::arg("bottleneck") = false)
      .def_readonly("edge_flows", &OutputOptions::edge_flows, "request PlacementResult.edge_flows")
      .def_readonly("residual", &OutputOptions::residual, "request PlacementResult.residual")
      .def_readonly("bottleneck", &OutputOptions::bottleneck,
                    "request PlacementResult.bottleneck (and min_cut on Global)");
}

} // namespace netgraph::core::bindings
