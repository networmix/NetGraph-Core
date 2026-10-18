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

py::object inf_or(const std::optional<Amount>& r) {
  return r.has_value() ? py::object(py::int_(*r)) : py::object(py::float_(std::numeric_limits<double>::infinity()));
}

py::object counter_or_none(const std::optional<Counter128>& c) {
  return c.has_value() ? counter_int(*c) : py::object(py::none());
}

py::list counter_pairs(const std::vector<std::pair<EdgeId, Counter128>>& v) {
  py::list out;
  for (const auto& [e, x] : v) out.append(py::make_tuple(e, counter_int(x)));
  return out;
}

} // namespace

void register_placement_result(py::module_& m) {
  py::enum_<Termination>(m, "Termination",
                         "How a call ended, the first that applies: COMPLETE for a zero request; NO_ROUTE for a positive\n"
                         "request with an empty roster (no eligible source statically reaches an allowed destination);\n"
                         "NO_ACTIVE_WEIGHT when routes exist but none has a positive weight (WEIGHTED over edges of\n"
                         "capacity 0 only); a LOSSY finite request: COMPLETE (delivery and losses account for the whole\n"
                         "offer, which does not mean nothing was lost); otherwise COMPLETE when the whole request was\n"
                         "delivered, else EXHAUSTED (an infinite request is EXHAUSTED). A controller whose maximum is 0\n"
                         "over existing routes is EXHAUSTED, not NO_ROUTE.\n"
                         "CapacityAware alone reduces per-source outcomes: a source is served when its offer is 0 or\n"
                         "delivered; an unserved source's outcome is LIMIT_FLOWS (creating or adding flows was stopped by\n"
                         "MAX_FLOWS_PER_CALL, or adding flows by a max_flow_count above flow_count), else\n"
                         "LIMIT_DIMINISHING (the diminishing guard), else NO_ROUTE (no flows), else EXHAUSTED; the call\n"
                         "takes the highest unserved outcome in that order, COMPLETE when every source is served.")
      .value("COMPLETE", Termination::Complete)
      .value("LIMIT_FLOWS", Termination::LimitFlows)
      .value("LIMIT_DIMINISHING", Termination::LimitDiminishing)
      .value("NO_ROUTE", Termination::NoRoute)
      .value("EXHAUSTED", Termination::Exhausted)
      .value("NO_ACTIVE_WEIGHT", Termination::NoActiveWeight);

  py::class_<SourceFlow>(m, "SourceFlow",
                         "One normalized input source's forwarding in the epoch (read-only; zero rows included: a\n"
                         "masked or unreachable source): source, offered (LOSSY: its injections; None when\n"
                         "LOSSLESS), admitted (LOSSLESS: what it delivers; LOSSY: offered), delivered, source_loss\n"
                         "(offered units lost at the source for lack of a usable route) and dropped (its losses on\n"
                         "links). LOSSY: offered = delivered + source_loss + dropped. Exact ints.")
      .def_readonly("source", &SourceFlow::source, "the source node id")
      .def_property_readonly(
          "offered", [](const SourceFlow& r) { return counter_or_none(r.offered); },
          "LOSSY: the units this source injected; None when LOSSLESS")
      .def_property_readonly(
          "admitted", [](const SourceFlow& r) { return counter_int(r.admitted); },
          "LOSSLESS: the units it delivers; LOSSY: offered")
      .def_readonly("delivered", &SourceFlow::delivered, "the units it delivers to the destinations")
      .def_property_readonly(
          "source_loss", [](const SourceFlow& r) { return counter_int(r.source_loss); },
          "LOSSY: its offered units lost at the source (no usable route); 0 when LOSSLESS")
      .def_property_readonly(
          "dropped", [](const SourceFlow& r) { return counter_int(r.dropped); },
          "LOSSY: its units dropped on links; 0 when LOSSLESS")
      .def("__repr__", [](const SourceFlow& r) {
        return "SourceFlow(source=" + std::to_string(r.source) + ", delivered=" + std::to_string(r.delivered) + ")";
      });

  py::class_<PlacementResult>(
      m, "PlacementResult",
      "The result of one call. This call: placed (units delivered by the call), remaining (request -\n"
      "placed; float('inf') for an infinite request), delivery_delta (new minus old delivery: placed for a\n"
      "placement), call_drops [(EdgeId, units)], call_source_loss, termination, routes_changed and\n"
      "weights_changed. The demand's state after the call (totals of the epoch): requested (float('inf')\n"
      "after an infinite request), offered (LOSSY: the traffic sent; None when LOSSLESS), admitted,\n"
      "delivered, unadmitted (LOSSLESS: requested - delivered, None after an infinite request; LOSSY: 0),\n"
      "source_loss, drops [(EdgeId, units)], by_source (one SourceFlow per normalized input source,\n"
      "ascending NodeId), total_cost (controllers: Global and CostOnly MAX_FLOW; else None) and\n"
      "cost_distribution (Global with the optimized origin only). LOSSY: offered = source_loss +\n"
      "sum(drops) + delivered. The rows sum to the totals, except for LOSSY with an empty roster, where the\n"
      "whole offer is the demand's source loss and every row is zero. Optional outputs, None unless\n"
      "requested: edge_flows and residual (int64[E]), bottleneck, min_cut. Amounts that accumulate are\n"
      "exact ints (they may exceed 64 bits).")
      .def_readonly("placed", &PlacementResult::placed, "the units placed (delivered) by this call")
      .def_property_readonly(
          "remaining", [](const PlacementResult& r) { return inf_or(r.remaining); },
          "the request minus placed, in units (an amount, not a status); float('inf') for an infinite request")
      .def_readonly("delivery_delta", &PlacementResult::delivery_delta,
                    "the delivery after the call minus before it (placed, for a placement)")
      .def_property_readonly(
          "call_drops", [](const PlacementResult& r) { return pairs_list(r.call_drops); },
          "this call's drops on links [(EdgeId, units)], ascending EdgeId (empty unless LOSSY)")
      .def_readonly("call_source_loss", &PlacementResult::call_source_loss,
                    "this call's offered units lost at sources without a usable route (LOSSY)")
      .def_readonly("termination", &PlacementResult::termination, "the call's outcome (Termination)")
      .def_readonly("routes_changed", &PlacementResult::routes_changed,
                    "True when this call committed a route change: routes prepared or a flow created (also at\n"
                    "zero delivery), rerouted or added, or a change in the edges a Global source's share uses; a\n"
                    "candidate equal to the current route is not a change")
      .def_readonly("weights_changed", &PlacementResult::weights_changed,
                    "True when this call committed new weights for a rerouted WEIGHTED flow")
      .def_property_readonly(
          "requested",
          [](const PlacementResult& r) -> py::object {
            return r.requested.has_value() ? counter_int(*r.requested)
                                           : py::object(py::float_(std::numeric_limits<double>::infinity()));
          },
          "the sum of the epoch's requests; float('inf') after an infinite one")
      .def_property_readonly(
          "offered", [](const PlacementResult& r) { return counter_or_none(r.offered); },
          "LOSSY: the traffic sent in the epoch (the first-hop injection of an infinite request); None\n"
          "when LOSSLESS")
      .def_property_readonly(
          "admitted", [](const PlacementResult& r) { return counter_int(r.admitted); },
          "LOSSLESS: the total carried; LOSSY: offered")
      .def_readonly("delivered", &PlacementResult::delivered, "the units the demand delivers")
      .def_property_readonly(
          "unadmitted", [](const PlacementResult& r) { return counter_or_none(r.unadmitted); },
          "LOSSLESS: requested - delivered, None after an infinite request; LOSSY: 0")
      .def_property_readonly(
          "source_loss", [](const PlacementResult& r) { return counter_int(r.source_loss); },
          "LOSSY: the offered units lost at sources (no usable route, or an empty roster); 0 when LOSSLESS")
      .def_property_readonly(
          "drops", [](const PlacementResult& r) { return counter_pairs(r.drops); },
          "the epoch's drops on links [(EdgeId, units)], ascending EdgeId (empty unless LOSSY)")
      .def_property_readonly(
          "by_source", [](const PlacementResult& r) { return py::tuple(py::cast(r.by_source)); },
          "one SourceFlow per normalized input source, ascending NodeId (zero rows included)")
      .def_property_readonly(
          "total_cost", [](const PlacementResult& r) { return counter_or_none(r.total_cost); },
          "controllers (Global, CostOnly MAX_FLOW): the exact sum of cost * flow over the policy's\n"
          "allocation; None for the other models")
      .def_readonly("cost_distribution", &PlacementResult::cost_distribution,
                    "Global with the optimized origin: the marginal cost curve [(cost, units)] of the policy's\n"
                    "cumulative flow (a zero call reports it unchanged; None before the epoch's first positive\n"
                    "call); in ascending cost, delivering `units` more raised the minimum total flow cost by `cost`\n"
                    "per unit (more delivery can reroute earlier flow, so these are not the costs of the final\n"
                    "paths). None for the other models and the equal origin")
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
                    "if requested (OutputOptions(min_cut=True)), Global with the optimized origin, when the solve\n"
                    "reached the maximum flow: the minimum cut closest to the sources, every allowed edge leaving\n"
                    "the nodes the roster sources reach in the residual network (free capacity plus the policy's\n"
                    "own allocation), an edge with nothing free (capacity 0, or held by others) included,\n"
                    "ascending EdgeId; [] for an empty roster. None otherwise (another model, the equal origin, a\n"
                    "call below the maximum, a zero call)")
      .def("__repr__", [](const PlacementResult& r) {
        return "PlacementResult(placed=" + std::to_string(r.placed) + ", termination=" +
               std::string(py::str(py::cast(r.termination).attr("name"))) + ")";
      });
}

} // namespace netgraph::core::bindings
