// Termination, SourceFlow, TunnelFlow, TunnelPlan and PlacementResult: the outcome of one policy call.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <limits>
#include <memory>
#include <optional>
#include <string>

#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/program.hpp"
#include "values.hpp"
#include "publication.hpp"
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
                         "NO_ACTIVE_WEIGHT when the program has an entry with a path (an up tunnel, a controller path or\n"
                         "a split route) but no head can use one (no positive weight: WEIGHTED over edges of capacity 0\n"
                         "only, routes of weight 0, every up tunnel failed); a program whose tunnels are all down is not\n"
                         "this case; a LOSSY finite request: COMPLETE (delivery and losses account for the whole offer,\n"
                         "which does not mean nothing was lost); otherwise COMPLETE when the whole request was delivered,\n"
                         "else EXHAUSTED (an infinite request is EXHAUSTED). A controller whose maximum is 0 over existing\n"
                         "routes is EXHAUSTED, not NO_ROUTE.")
      .value("COMPLETE", Termination::Complete)
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

  py::class_<TunnelFlow>(m, "TunnelFlow",
                         "One program tunnel's forwarding in the epoch (read-only): index (its entry identifier) and\n"
                         "head, status (TunnelStatus: UP; NO_BANDWIDTH or COST_CEILING for a tunnel T1's setup left\n"
                         "down; FAILED, BLACKHOLED in this forwarding), assigned (its share of the injections: LOSSY the\n"
                         "shares sent to it, LOSSLESS what it carries), delivered and dropped (its losses on links).\n"
                         "Exact ints (assigned and dropped may exceed 64 bits).")
      .def_readonly("index", &TunnelFlow::index, "the tunnel's entry identifier")
      .def_readonly("head", &TunnelFlow::head, "the tunnel's head")
      .def_readonly("status", &TunnelFlow::status,
                    "UP, a down tunnel's setup status (NO_BANDWIDTH, COST_CEILING), FAILED or BLACKHOLED in this\n"
                    "forwarding")
      .def_property_readonly(
          "assigned", [](const TunnelFlow& r) { return counter_int(r.assigned); },
          "its share of the injections in the epoch (LOSSY: sent to it; LOSSLESS: what it carries)")
      .def_readonly("delivered", &TunnelFlow::delivered, "the units it delivers")
      .def_property_readonly(
          "dropped", [](const TunnelFlow& r) { return counter_int(r.dropped); }, "its units dropped on links")
      .def("__repr__", [](const TunnelFlow& r) {
        return "TunnelFlow(index=" + std::to_string(r.index) + ", delivered=" + std::to_string(r.delivered) + ")";
      });

  py::class_<TunnelPlan>(m, "TunnelPlan",
                         "The tunnel models' setup summary (read-only): count (the tunnels per head), design_demand\n"
                         "(the design volume D in units), head_design ((head, RR(D, H, i)) for every roster head in\n"
                         "ascending NodeId) and reserved_total (T1's strict holds, counted once per up tunnel: the sum\n"
                         "of their requests). Exact ints.")
      .def_readonly("count", &TunnelPlan::count, "the tunnels per head")
      .def_property_readonly(
          "design_demand", [](const TunnelPlan& p) { return inf_or(p.design_demand); },
          "the design volume D in units; float('inf') when unlimited")
      .def_property_readonly(
          "head_design", [](const TunnelPlan& p) { return pairs_tuple(p.head_design); },
          "((head, share), ...): every roster head's share RR(D, H, i), ascending NodeId")
      .def_readonly("reserved_total", &TunnelPlan::reserved_total,
                    "the strict holds of the up tunnels, each tunnel's request counted once")
      .def("__repr__", [](const TunnelPlan& p) {
        return "TunnelPlan(count=" + std::to_string(p.count) + ", reserved_total=" + std::to_string(p.reserved_total) + ")";
      });

  py::class_<PlacementResult>(
      m, "PlacementResult",
      "The result of one call. This call: placed (units delivered by the call; replace, forward and\n"
      "evaluate: the full new delivery), remaining (request - placed; float('inf') for an infinite\n"
      "request), delivery_delta (new minus old delivery, may be negative), call_drops [(EdgeId, units)],\n"
      "call_source_loss, termination, routes_changed, weights_changed and reservations_changed. The\n"
      "demand's state after the call (totals of the epoch): program (the installed ForwardingProgram;\n"
      "evaluate: the computed one; None before the first preparation), plan (a TunnelPlan for T1; None for\n"
      "every other model), requested (float('inf') after an infinite request), offered (LOSSY: the traffic sent; None when LOSSLESS), admitted, delivered, unadmitted\n"
      "(LOSSLESS: requested - delivered, None after an infinite request; LOSSY: 0), source_loss, drops\n"
      "[(EdgeId, units)], by_source (one SourceFlow per normalized input source, ascending NodeId),\n"
      "by_tunnel (one TunnelFlow per program tunnel, ascending index; an alternative view, never added to\n"
      "by_source), total_cost (controllers: Global and CostOnly MAX_FLOW; else None) and cost_distribution\n"
      "(unrestricted Global with the optimized origin only). LOSSY: offered = source_loss + sum(drops) +\n"
      "delivered. The rows sum to the totals, except for LOSSY with an empty roster, where the whole offer\n"
      "is the demand's source loss and every row is zero. Optional outputs, None unless requested:\n"
      "edge_flows, residual and reservation_edge_flows (int64[E]), bottleneck, min_cut (they describe the same\n"
      "state; evaluate: the evaluated one). Amounts that accumulate are exact ints (they may exceed 64 bits).")
      .def_readonly("placed", &PlacementResult::placed, "the units placed (delivered) by this call")
      .def_property_readonly(
          "remaining", [](const PlacementResult& r) { return inf_or(r.remaining); },
          "the request minus placed, in units (an amount, not a status); float('inf') for an infinite request")
      .def_readonly("delivery_delta", &PlacementResult::delivery_delta,
                    "the delivery after the call minus before it (placed for a placement; may be negative for a\n"
                    "full-offer call)")
      .def_property_readonly(
          "call_drops", [](const PlacementResult& r) { return pairs_list(r.call_drops); },
          "this call's drops on links [(EdgeId, units)], ascending EdgeId (empty unless LOSSY)")
      .def_readonly("call_source_loss", &PlacementResult::call_source_loss,
                    "this call's offered units lost at sources without a usable route (LOSSY)")
      .def_readonly("termination", &PlacementResult::termination, "the call's outcome (Termination)")
      .def_readonly("routes_changed", &PlacementResult::routes_changed,
                    "True when the call's program differs from the installed one in its ordered heads, paths, split\n"
                    "tables or duplicates (no program counts as the empty program; two programs without an entry\n"
                    "are equal). An offer-only change (forward, Pinned, CostOnly increments) and an identical\n"
                    "replace set nothing")
      .def_readonly("weights_changed", &PlacementResult::weights_changed,
                    "True when the call's program differs from the installed one in its injection (source weights)\n"
                    "or any entry or split-table weight (a controller's re-solve that changes a path's amount; an\n"
                    "absent entry weighs 0)")
      .def_readonly("reservations_changed", &PlacementResult::reservations_changed,
                    "True when a tunnel's requested or reserved amount differs from the installed program's, by\n"
                    "index (an absent tunnel holds 0): only T1 sets it; forward and Pinned never do")
      .def_property_readonly(
          "program",
          [](const PlacementResult& r) { return optional_shared(r.program); },
          "the installed ForwardingProgram after the call (evaluate: the computed one); None before the first\n"
          "preparation")
      .def_property_readonly(
          "plan", [](const PlacementResult& r) { return optional_of(r.plan); },
          "T1: the TunnelPlan of the program after the call; None for every other model (Pinned over a T1\n"
          "program too)")
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
          "by_source", [](const PlacementResult& r) { return tuple_of(r.by_source); },
          "one SourceFlow per normalized input source, ascending NodeId (zero rows included)")
      .def_property_readonly(
          "by_tunnel", [](const PlacementResult& r) { return tuple_of(r.by_tunnel); },
          "one TunnelFlow per program tunnel, ascending index (empty without tunnels)")
      .def_property_readonly(
          "total_cost", [](const PlacementResult& r) { return counter_or_none(r.total_cost); },
          "controllers (Global, CostOnly MAX_FLOW): the exact sum of cost * flow over the policy's\n"
          "allocation; None for the other models (T1, Pinned)")
      .def_readonly("cost_distribution", &PlacementResult::cost_distribution,
                    "unrestricted Global with the optimized origin: the marginal cost curve [(cost, units)] of the policy's\n"
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
      .def_property_readonly(
          "reservation_edge_flows",
          [](const PlacementResult& r) -> py::object {
            return r.reservation_edge_flows ? py::object(int64_array(*r.reservation_edge_flows)) : py::object(py::none());
          },
          "if requested: the policy's own strict holds per edge after the call (int64[E], EdgeId order; zeros\n"
          "for a model without holds), else None")
      .def_readonly("bottleneck", &PlacementResult::bottleneck,
                    "if requested: the ledger's saturation report after the call (FlowLedger.bottleneck of the\n"
                    "policy's sources and Eligibility), else None")
      .def_readonly("min_cut", &PlacementResult::min_cut,
                    "if requested (OutputOptions(min_cut=True)), unrestricted Global with the optimized origin, when the solve\n"
                    "reached the maximum flow: the minimum cut closest to the sources, every allowed edge leaving\n"
                    "the nodes the roster sources reach in the residual network (free capacity plus the policy's\n"
                    "own allocation), an edge with nothing free (capacity 0, or held by others) included,\n"
                    "ascending EdgeId; [] for an empty roster. None otherwise (another model, the equal origin, a\n"
                    "support, a call below the maximum, a zero call, forward)")
      .def("__repr__", [](const PlacementResult& r) {
        return "PlacementResult(placed=" + std::to_string(r.placed) + ", termination=" +
               std::string(py::str(py::cast(r.termination).attr("name"))) + ")";
      });
}

} // namespace netgraph::core::bindings
