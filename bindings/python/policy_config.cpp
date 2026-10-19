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
      "One demand's settings: routing (CostOnly, CapacityAware (T1), Global or Pinned), placement\n"
      "(FlowPlacement: required for CostOnly; the other models take none), admission (LOSSLESS by\n"
      "default; LOSSY for CostOnly EQUAL_BALANCED and WEIGHTED, T1 and Pinned), selection (PathSelection;\n"
      "None: PathSelection(); Global and Pinned take none; T1 reads tie_break, priority and multi_edge and\n"
      "rejects multipath=True), origin (Origin; None: the model's default, EQUAL for CostOnly EQUAL_BALANCED\n"
      "and WEIGHTED and for T1, OPTIMIZED for CostOnly MAX_FLOW and Global; Pinned takes none, its program\n"
      "has it), weights (WeightBasis, CostOnly WEIGHTED only; None: CAPACITY) and edge_weights (one integer\n"
      "per edge, given exactly when weights=EXPLICIT; one array shared by every policy built from this\n"
      "configuration). A field that does not apply to the model must stay None. The combination is checked\n"
      "when a policy or an evaluator uses the configuration (ValueError): the enumeration values, then the\n"
      "model's row (T1: an up_to count first; then placement, origin, admission, selection, weights with\n"
      "edge_weights), then a DagCache (only CostOnly reads one), then T1's parameters (the tunnel count in\n"
      "[1, 10000], max_path_cost, min_unreserved_bandwidth), then edge_weights' length and entries, then the\n"
      "roster and count rules (CostOnly with several input sources rejects multipath=False with\n"
      "PREFER_HIGHER_PRIORITY; T1's roster heads times its tunnel count at most 2^30 - 1).")
      .def(py::init([](const py::object& routing, std::optional<FlowPlacement> placement, Admission admission,
                       std::optional<PathSelection> selection, std::optional<Origin> origin,
                       std::optional<WeightBasis> weights, const py::object& edge_weights) {
             FlowPolicyConfig cfg {routing_of(routing), placement, admission, selection, origin, weights, {}};
             if (!edge_weights.is_none()) {
               // Converted once and shared: every policy copied from this configuration reads the same
               // array. An entry beyond int64 is outside the domain here; the policy checks the rest.
               (void)sequence_length(edge_weights, "edge_weights");
               const auto w = wide_integers(edge_weights, "edge_weights");
               for (std::size_t e = 0; e < w.size(); ++e) {
                 if (w.sign(e) != 0) throw py::value_error("edge_weights[" + std::to_string(e) + "] must be an integer in [0, 2^62)");
               }
               cfg.edge_weights = std::make_shared<const std::vector<Amount>>(w.narrow<Amount>());
             }
             return cfg;
           }),
           py::kw_only(), py::arg("routing"), py::arg("placement") = py::none(),
           py::arg_v("admission", Admission::Lossless, "Admission.LOSSLESS"), py::arg("selection") = py::none(),
           py::arg("origin") = py::none(), py::arg("weights") = py::none(), py::arg("edge_weights") = py::none())
      .def_property_readonly(
          "routing", [](const FlowPolicyConfig& c) { return routing_obj(c.routing); },
          "the routing mode: CostOnly, CapacityAware, Global or Pinned")
      .def_readonly("placement", &FlowPolicyConfig::placement, "the FlowPlacement, or None (Global, Pinned)")
      .def_readonly("admission", &FlowPolicyConfig::admission, "LOSSLESS or LOSSY")
      .def_readonly("selection", &FlowPolicyConfig::selection,
                    "the PathSelection, or None (a policy's config: PathSelection() resolved; Global, Pinned: None)")
      .def_readonly("origin", &FlowPolicyConfig::origin,
                    "the Origin, or None (a policy's config: the model's default resolved; Pinned: None)")
      .def_readonly("weights", &FlowPolicyConfig::weights,
                    "the WeightBasis, or None (a policy's config: CAPACITY resolved for CostOnly WEIGHTED)")
      .def_property_readonly(
          "edge_weights",
          [](const FlowPolicyConfig& c) -> py::object {
            if (!c.edge_weights) return py::none();
            return py::tuple(py::cast(*c.edge_weights));
          },
          "the explicit weights (one per edge, a tuple), or None")
      .def("__repr__", [](const FlowPolicyConfig& c) {
        return "FlowPolicyConfig(routing=" + std::string(py::repr(routing_obj(c.routing))) + ")";
      });

  py::class_<OutputOptions>(m, "OutputOptions",
                            "Optional outputs of a placement (all false by default; absent outputs are None):\n"
                            "edge_flows (the policy's own allocation per edge after the call), residual and\n"
                            "bottleneck (the ledger's after the call), reservation_edge_flows (the policy's own strict\n"
                            "holds per edge after the call; zeros for a model without holds), and min_cut (Global with\n"
                            "the optimized origin, when its solve reached the maximum flow; None for every other\n"
                            "model, origin or call).")
      .def(py::init([](bool edge_flows, bool residual, bool bottleneck, bool reservation_edge_flows, bool min_cut) {
             OutputOptions o;
             o.edge_flows = edge_flows;
             o.residual = residual;
             o.bottleneck = bottleneck;
             o.reservation_edge_flows = reservation_edge_flows;
             o.min_cut = min_cut;
             return o;
           }),
           py::kw_only(), py::arg("edge_flows") = false, py::arg("residual") = false, py::arg("bottleneck") = false,
           py::arg("reservation_edge_flows") = false, py::arg("min_cut") = false)
      .def_readonly("edge_flows", &OutputOptions::edge_flows, "request PlacementResult.edge_flows")
      .def_readonly("residual", &OutputOptions::residual, "request PlacementResult.residual")
      .def_readonly("bottleneck", &OutputOptions::bottleneck, "request PlacementResult.bottleneck")
      .def_readonly("reservation_edge_flows", &OutputOptions::reservation_edge_flows,
                    "request PlacementResult.reservation_edge_flows")
      .def_readonly("min_cut", &OutputOptions::min_cut, "request PlacementResult.min_cut");
}

} // namespace netgraph::core::bindings
