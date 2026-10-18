// Termination and PlacementResult: the outcome of one placement call.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <limits>
#include <optional>
#include <string>

#include "netgraph/core/flow_policy.hpp"
#include "values.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

py::object remaining_obj(const std::optional<Amount>& r) {
  return r.has_value() ? py::object(py::int_(*r)) : py::object(py::float_(std::numeric_limits<double>::infinity()));
}

} // namespace

void register_placement_result(py::module_& m) {
  py::enum_<Termination>(m, "Termination",
                         "How a call ended, per source (local_outcomes) and for the whole call (termination).\n"
                         "A source's offer is served when it is 0, when (LOSSLESS) the delivery equals it, or when (LOSSY,\n"
                         "finite offer) all of it was forwarded, drops included. A source's outcome is the first that\n"
                         "applies:\n"
                         "  COMPLETE: served.\n"
                         "  LIMIT_FLOWS: creating or adding flows was stopped by MAX_FLOWS_PER_CALL, or adding flows by a\n"
                         "    CapacityAware max_flow_count above flow_count.\n"
                         "  LIMIT_DIMINISHING: the diminishing guard stopped adding flows.\n"
                         "  NO_ROUTE: the source has no flows after the call (Global: the maximum flow is 0 and nothing\n"
                         "    was ever placed).\n"
                         "  EXHAUSTED: any other unserved case (every flow full and no new route helps, a full fixed-size\n"
                         "    set of flows, Global at the maximum flow, a LOSSY infinite offer).\n"
                         "Global and Pinned act as one source. The call's termination is COMPLETE for a request of 0 or\n"
                         "when every source is served, NO_ROUTE for a positive request with an empty roster, else the\n"
                         "highest outcome of the unserved sources in the order LIMIT_FLOWS > LIMIT_DIMINISHING >\n"
                         "EXHAUSTED > NO_ROUTE.")
      .value("COMPLETE", Termination::Complete)
      .value("LIMIT_FLOWS", Termination::LimitFlows)
      .value("LIMIT_DIMINISHING", Termination::LimitDiminishing)
      .value("NO_ROUTE", Termination::NoRoute)
      .value("EXHAUSTED", Termination::Exhausted);

  py::class_<PlacementResult>(
      m, "PlacementResult",
      "One placement call: placed (units), remaining (request - placed; float('inf') for an infinite\n"
      "request), termination, local_outcomes [(source, Termination)] per roster source, routes_changed,\n"
      "weights_changed, and this call's drops [(EdgeId, units)]. cost_distribution: Global only, the\n"
      "marginal cost curve [(cost, units)] of the policy's flow so far, else None: in ascending cost,\n"
      "delivering `units` more raised the minimum total flow cost by `cost` per unit (more delivery can\n"
      "reroute earlier flow, so these are not the costs of the final paths); a zero call leaves it\n"
      "unchanged. Optional outputs, None unless requested:\n"
      "edge_flows and residual (int64[E]), bottleneck, and min_cut (Global with bottleneck; None when the\n"
      "solve did not reach the maximum flow).")
      .def_readonly("placed", &PlacementResult::placed, "the units placed by this call")
      .def_property_readonly(
          "remaining", [](const PlacementResult& r) { return remaining_obj(r.remaining); },
          "the request minus placed, in units (an amount, not a status); float('inf') for an infinite request")
      .def_readonly("termination", &PlacementResult::termination, "the call's outcome (Termination)")
      .def_readonly("local_outcomes", &PlacementResult::local_outcomes,
                    "[(source, Termination)] for every roster source, ascending node id")
      .def_readonly("routes_changed", &PlacementResult::routes_changed,
                    "True when this call committed a route change: a flow created (also at zero delivery), rerouted\n"
                    "or added, or a change in the edges a Global source's share uses; a candidate equal to the\n"
                    "current route is not a change")
      .def_readonly("weights_changed", &PlacementResult::weights_changed,
                    "True when this call committed new weights for a rerouted WEIGHTED flow")
      .def_property_readonly(
          "drops", [](const PlacementResult& r) { return pairs_list(r.drops); },
          "this call's transit drops [(EdgeId, units)], ascending EdgeId (empty unless Lossy; records\n"
          "keep the cumulative drops)")
      .def_readonly("cost_distribution", &PlacementResult::cost_distribution,
                    "Global: the marginal cost curve [(cost, units)] of the policy's cumulative flow (a zero\n"
                    "call reports it unchanged); None for the other routing modes")
      .def_property_readonly(
          "edge_flows",
          [](const PlacementResult& r) -> py::object {
            return r.edge_flows ? py::object(int64_array(*r.edge_flows)) : py::object(py::none());
          },
          "if requested: the policy's own allocation per edge after the call (int64[E], EdgeId order),\n"
          "else None")
      .def_property_readonly(
          "residual",
          [](const PlacementResult& r) -> py::object {
            return r.residual ? py::object(int64_array(*r.residual)) : py::object(py::none());
          },
          "if requested: the ledger's residual after the call (int64[E], EdgeId order), else None")
      .def_readonly("bottleneck", &PlacementResult::bottleneck,
                    "if requested: the ledger's saturation report after the call (FlowLedger.bottleneck of the\n"
                    "policy's sources and Eligibility), else None")
      .def_readonly("min_cut", &PlacementResult::min_cut,
                    "Global with bottleneck requested: a minimum cut (edges whose removal separates the sources from\n"
                    "the destinations, with the smallest total capacity) of the policy's view (free capacity plus its\n"
                    "own allocation), when the solve reached the maximum flow: every allowed edge leaving the nodes the\n"
                    "sources still reach, an edge with nothing free (capacity 0, or held by others) included; None\n"
                    "otherwise (also on a zero call)")
      .def("__repr__", [](const PlacementResult& r) {
        return "PlacementResult(placed=" + std::to_string(r.placed) + ", termination=" +
               std::string(py::str(py::cast(r.termination).attr("name"))) + ")";
      });
}

} // namespace netgraph::core::bindings
