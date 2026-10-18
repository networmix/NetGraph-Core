// FlowPolicyConfig and OutputOptions: the configuration of a placement.
#include <pybind11/pybind11.h>

#include <optional>
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
      "One demand's settings: routing (CostOnly, CapacityAware, Global or Pinned), placement\n"
      "(FlowPlacement: required for CostOnly, CapacityAware and Pinned; Global takes none), admission\n"
      "(LOSSLESS by default), selection (PathSelection; None: PathSelection(); Global takes none) and\n"
      "origin (Origin; None: the model's default, EQUAL for CostOnly EQUAL_BALANCED and WEIGHTED and for\n"
      "CapacityAware, OPTIMIZED for CostOnly MAX_FLOW and Global; Pinned takes none). A field that does\n"
      "not apply to the model must stay None. The combination is checked when a policy or an evaluator\n"
      "uses the configuration (ValueError): the enumeration values, then the model's row (placement,\n"
      "origin, admission, selection), then a DagCache (only CostOnly reads one), then CapacityAware's\n"
      "parameters, then the pooled single-path rule (CostOnly with several input sources rejects\n"
      "multipath=False with PREFER_HIGHER_PRIORITY).")
      .def(py::init([](const py::object& routing, std::optional<FlowPlacement> placement, Admission admission,
                       std::optional<PathSelection> selection, std::optional<Origin> origin) {
             return FlowPolicyConfig {routing_of(routing), placement, admission, selection, origin};
           }),
           py::kw_only(), py::arg("routing"), py::arg("placement") = py::none(),
           py::arg_v("admission", Admission::Lossless, "Admission.LOSSLESS"), py::arg("selection") = py::none(),
           py::arg("origin") = py::none())
      .def_property_readonly(
          "routing", [](const FlowPolicyConfig& c) { return routing_obj(c.routing); },
          "the routing mode: CostOnly, CapacityAware, Global or Pinned")
      .def_readonly("placement", &FlowPolicyConfig::placement, "the FlowPlacement, or None (Global)")
      .def_readonly("admission", &FlowPolicyConfig::admission, "LOSSLESS or LOSSY")
      .def_readonly("selection", &FlowPolicyConfig::selection,
                    "the PathSelection, or None (a policy's config: PathSelection() resolved; Global: None)")
      .def_readonly("origin", &FlowPolicyConfig::origin,
                    "the Origin, or None (a policy's config: the model's default resolved; Pinned: None)")
      .def("__repr__", [](const FlowPolicyConfig& c) {
        return "FlowPolicyConfig(routing=" + std::string(py::repr(routing_obj(c.routing))) + ")";
      });

  py::class_<OutputOptions>(m, "OutputOptions",
                            "Optional outputs of a placement (all false by default; absent outputs are None):\n"
                            "edge_flows (the policy's own allocation per edge after the call), residual and\n"
                            "bottleneck (the ledger's after the call), and min_cut (Global with the optimized origin,\n"
                            "when its solve reached the maximum flow; None for every other model, origin or call).")
      .def(py::init<bool, bool, bool, bool>(), py::kw_only(), py::arg("edge_flows") = false,
           py::arg("residual") = false, py::arg("bottleneck") = false, py::arg("min_cut") = false)
      .def_readonly("edge_flows", &OutputOptions::edge_flows, "request PlacementResult.edge_flows")
      .def_readonly("residual", &OutputOptions::residual, "request PlacementResult.residual")
      .def_readonly("bottleneck", &OutputOptions::bottleneck, "request PlacementResult.bottleneck")
      .def_readonly("min_cut", &OutputOptions::min_cut, "request PlacementResult.min_cut");
}

} // namespace netgraph::core::bindings
