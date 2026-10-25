// Termination, SourceFlow, TunnelFlow, TunnelPlan and PlacementResult: the outcome of one policy call.
#include <pybind11/pybind11.h>
#include <pybind11/gil_safe_call_once.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <limits>
#include <memory>
#include <optional>
#include <string>

#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/forwarding_program.hpp"
#include "values.hpp"
#include "publication.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// A Termination's name, without a Python object (a repr must not fail on its own allocation of an
// enum object).
const char* termination_name(Termination t) {
  switch (t) {
    case Termination::Complete: return "COMPLETE";
    case Termination::NoRoute: return "NO_ROUTE";
    case Termination::Exhausted: return "EXHAUSTED";
    case Termination::NoActiveWeight: return "NO_ACTIVE_WEIGHT";
  }
  return "?";
}

py::object inf_or(const std::optional<Amount>& r) {
  return r.has_value() ? py::object(int_object(*r)) : py::object(float_object(std::numeric_limits<double>::infinity()));
}

py::object counter_or_none(const std::optional<Counter128>& c) {
  return c.has_value() ? counter_int(*c) : py::object(py::none());
}

// fractions.Fraction, imported once when the module is registered: TunnelPlan.max_fractional_bandwidth_per_tunnel only
// calls it (an import inside a getter can absorb a failed allocation and still succeed).
py::gil_safe_call_once_and_store<py::object>& fraction_class() {
  PYBIND11_CONSTINIT static py::gil_safe_call_once_and_store<py::object> storage;
  return storage;
}

} // namespace

void register_placement_result(py::module_& m) {
  fraction_class().call_once_and_store_result([] { return py::module_::import("fractions").attr("Fraction"); });
  py::enum_<Termination>(m, "Termination",
                         "How a call ended, the first that applies: COMPLETE for a zero request; NO_ROUTE for a positive\n"
                         "request with an empty roster (no eligible source statically reaches an allowed destination);\n"
                         "NO_ACTIVE_WEIGHT when the program has an entry with a path (an up tunnel, a controller path or\n"
                         "a split route) but no head can use one (no positive weight: WEIGHTED over edges of capacity 0\n"
                         "only, routes of weight 0, every up tunnel failed); a program whose tunnels are all down is not\n"
                         "this case. A live CentralizedTE WEIGHTED call without tunnel_count with a positive offer\n"
                         "and a nonempty all-zero source-weight roster is also NO_ACTIVE_WEIGHT, with no tunnels;\n"
                         "Pinned replay of that empty program follows the remaining rules. A LOSSY finite request: COMPLETE (delivery and losses account for the whole offer,\n"
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
                         "LOSSLESS), admitted (LOSSLESS: what it delivers; LOSSY: offered), delivered, source_drops\n"
                         "(offered units lost at the source for lack of a usable route) and total_edge_drops (its losses\n"
                         "on links). LOSSY: offered = delivered + source_drops + total_edge_drops. Exact ints.")
      .def_property_readonly("source", int_reader(&SourceFlow::source), "the source node id")
      .def_property_readonly(

          "offered", [](const SourceFlow& r) { return counter_or_none(r.offered); },
          "LOSSY: the units this source injected; None when LOSSLESS")
      .def_property_readonly(

          "admitted", [](const SourceFlow& r) { return counter_int(r.admitted); },
          "LOSSLESS: the units it delivers; LOSSY: offered")
      .def_property_readonly("delivered", int_reader(&SourceFlow::delivered), "the units it delivers to the destinations")
      .def_property_readonly(
          "source_drops", [](const SourceFlow& r) { return counter_int(r.source_drops); },
          "LOSSY: its offered units lost at the source (no usable route); 0 when LOSSLESS")
      .def_property_readonly(
          "total_edge_drops", [](const SourceFlow& r) { return counter_int(r.total_edge_drops); },
          "LOSSY: its units dropped on links; 0 when LOSSLESS")
      .def("__repr__", [](const SourceFlow& r) {
        return "SourceFlow(source=" + std::to_string(r.source) + ", delivered=" + std::to_string(r.delivered) + ")";
      });

  py::class_<TunnelFlow>(m, "TunnelFlow",
                         "One program tunnel's forwarding in the epoch (read-only): entry_id (its entry identifier) and\n"
                         "head, status (TunnelStatus: UP; INSUFFICIENT_BANDWIDTH or COST_CEILING_EXCEEDED for a tunnel\n"
                         "T1's setup left down; FAILED, BLACKHOLED in this forwarding), assigned (its share of the\n"
                         "injections: LOSSY the shares sent to it, LOSSLESS what it carries), delivered and\n"
                         "total_edge_drops (its losses on links). Exact ints (assigned and total_edge_drops may exceed 64\n"
                         "bits).")
      .def_property_readonly("entry_id", int_reader(&TunnelFlow::entry_id), "the tunnel's entry identifier")
      .def_property_readonly("head", int_reader(&TunnelFlow::head), "the tunnel's head")
      .def_property_readonly("status", field_object(&TunnelFlow::status),
                    "UP, a down tunnel's setup status (INSUFFICIENT_BANDWIDTH, COST_CEILING_EXCEEDED), FAILED or\n"
                    "BLACKHOLED in this forwarding")
      .def_property_readonly(

          "assigned", [](const TunnelFlow& r) { return counter_int(r.assigned); },
          "its share of the injections in the epoch (LOSSY: sent to it; LOSSLESS: what it carries)")
      .def_property_readonly("delivered", int_reader(&TunnelFlow::delivered), "the units it delivers")
      .def_property_readonly(
          "total_edge_drops", [](const TunnelFlow& r) { return counter_int(r.total_edge_drops); }, "its units dropped on links")
      .def("__repr__", [](const TunnelFlow& r) {
        return "TunnelFlow(entry_id=" + std::to_string(r.entry_id) + ", delivered=" + std::to_string(r.delivered) + ")";
      });

  py::class_<TunnelPlan>(m, "TunnelPlan",
                         "The tunnel models' setup summary (read-only): tunnels_per_head (T2 up to: the chosen count),\n"
                         "design_demand (the design demand D in units), design_by_head ((head, RR(D, H, i)) for every\n"
                         "roster head in ascending NodeId) and total_reserved_bandwidth (T1's strict holds, counted once\n"
                         "per up tunnel: the sum of their requests; 0 for T2). T2 only (None for T1): heads_with_capacity,\n"
                         "insufficient_resolution (live heads exist but no whole unit per tunnel fits, U == 0: the\n"
                         "program has no tunnels), max_fractional_bandwidth_per_tunnel (the exact fluid optimum s*, a\n"
                         "Fraction), max_bandwidth_per_tunnel (U), planned_bandwidth_per_tunnel (u),\n"
                         "total_planned_bandwidth (K*u), total_tunnel_cost (the route-cost sum, can exceed 2^64),\n"
                         "headroom_penalty (with prefer_headroom only) and max_lossless_offer (the program's lossless\n"
                         "maximum on the view it was planned on; a dead head's share blocks it). Exact ints. A plan\n"
                         "describes the computation that installed the tunnels: forward returns it unchanged, its\n"
                         "bounds not recomputed for the new offer, load or black holes.")
      .def_property_readonly("tunnels_per_head", int_reader(&TunnelPlan::tunnels_per_head), "the tunnels per head (T2 up to: the chosen count)")
      .def_property_readonly(

          "design_demand", [](const TunnelPlan& p) { return inf_or(p.design_demand); },
          "the design volume D in units; float('inf') when unlimited")
      .def_property_readonly(
          "design_by_head", [](const TunnelPlan& p) { return pairs_tuple(p.design_by_head); },
          "((head, share), ...): every roster head's share RR(D, H, i), ascending NodeId")
      .def_property_readonly("total_reserved_bandwidth", int_reader(&TunnelPlan::total_reserved_bandwidth),
                    "the strict holds of the up tunnels, each tunnel's request counted once")
      .def_property_readonly(
          "heads_with_capacity",
          [](const TunnelPlan& p) -> py::object {
            return p.heads_with_capacity.has_value() ? py::object(ints_tuple(*p.heads_with_capacity)) : py::object(py::none());
          },
          "T2: the heads that reach a destination over links with free capacity >= 1 (ascending); None for T1")
      .def_readonly("insufficient_resolution", &TunnelPlan::insufficient_resolution,
                    "T2: True when live heads exist but no whole unit per tunnel fits (max_bandwidth_per_tunnel 0,\n"
                    "no tunnels)")
      .def_property_readonly(
          "max_fractional_bandwidth_per_tunnel",
          [](const TunnelPlan& p) -> py::object {
            if (!p.max_fractional_bandwidth_per_tunnel.has_value()) return py::none();
            const auto& [n, d] = *p.max_fractional_bandwidth_per_tunnel;
            return fraction_object(fraction_class().get_stored(), n, d);
          },
          "T2: the exact fluid optimum s* (a fractions.Fraction; 0 without a live head); None for T1")
      .def_property_readonly("max_bandwidth_per_tunnel", optional_int_reader(&TunnelPlan::max_bandwidth_per_tunnel), "T2: U, the largest common whole-unit size")
      .def_property_readonly("planned_bandwidth_per_tunnel", optional_int_reader(&TunnelPlan::planned_bandwidth_per_tunnel), "T2: u, every tunnel's planned bandwidth")
      .def_property_readonly("total_planned_bandwidth", optional_int_reader(&TunnelPlan::total_planned_bandwidth), "T2: K*u over the K tunnels")
      .def_property_readonly(
          "total_tunnel_cost",
          [](const TunnelPlan& p) -> py::object { return p.total_tunnel_cost.has_value() ? counter_int(*p.total_tunnel_cost) : py::none(); },
          "T2: the sum of cost * tunnel count over the links (exact; can exceed 2^64)")
      .def_property_readonly(

          "headroom_penalty",
          [](const TunnelPlan& p) -> py::object {
            return p.headroom_penalty.has_value() ? counter_int(*p.headroom_penalty) : py::none();
          },
          "T2 with prefer_headroom: the sum over the links' tunnel slots j of floor(2^32 (B - V + j*u) / B)\n"
          "(B: static capacity; V: free capacity); None without prefer_headroom and for T1")
      .def_property_readonly("max_lossless_offer", optional_int_reader(&TunnelPlan::max_lossless_offer),
                    "T2: F, the program's lossless maximum on the view it was planned on (a dead head's share blocks\n"
                    "it); forward keeps the plan unchanged")
      .def("__repr__", [](const TunnelPlan& p) {
        return "TunnelPlan(tunnels_per_head=" + std::to_string(p.tunnels_per_head) + ", total_reserved_bandwidth=" + std::to_string(p.total_reserved_bandwidth) + ")";
      });

  py::class_<PlacementResult>(
      m, "PlacementResult",
      "The result of one call. This call: call_delivered (units delivered by the call; replace, forward\n"
      "and evaluate: the full new delivery), call_undelivered (request - call_delivered; float('inf') for\n"
      "an infinite request), delivery_delta (new minus old delivery, may be negative), call_edge_drops\n"
      "[(EdgeId, units)], call_source_drops, termination, routes_changed, weights_changed and\n"
      "reservations_changed. The\n"
      "demand's state after the call (totals of the epoch): program (the installed ForwardingProgram;\n"
      "evaluate: the computed one; None before the first preparation), plan (a TunnelPlan for T1 and T2;\n"
      "None for every other model), requested (float('inf') after an infinite request), offered (LOSSY: the traffic sent; None when LOSSLESS), admitted, delivered, unadmitted\n"
      "(LOSSLESS: requested - delivered, None after an infinite request; LOSSY: 0), source_drops,\n"
      "edge_drops [(EdgeId, units)], by_source (one SourceFlow per normalized input source, ascending\n"
      "NodeId), by_tunnel (one TunnelFlow per program tunnel, ascending entry_id; an alternative view, never\n"
      "added to by_source), total_cost (controllers: CentralizedTE and ShortestPaths OPTIMIZED; else None)\n"
      "and marginal_cost_curve (unrestricted CentralizedTE with the optimized source allocation only).\n"
      "LOSSY: offered = source_drops + sum(edge_drops) + delivered. The rows sum to the totals, except, under\n"
      "LOSSY, for an empty roster and for a forwarding whose remaining captured source weights are all 0\n"
      "(a controller program that sent nothing, or a mask that left only zero-weight heads): the whole\n"
      "offer is then the demand's source drops and every row is zero. Optional outputs, None unless requested: edge_flows, residual and edge_reservations\n"
      "(int64[E]), saturated_boundary_edges, min_cut (they describe the same\n"
      "state; evaluate: the evaluated one). Amounts that accumulate are exact ints (they may exceed 64 bits).")
      .def_property_readonly("call_delivered", int_reader(&PlacementResult::call_delivered), "the units placed (delivered) by this call")
      .def_property_readonly(
          "call_undelivered", [](const PlacementResult& r) { return inf_or(r.call_undelivered); },
          "the request minus call_delivered, in units (an amount, not a status); float('inf') for an infinite\n"
          "request")
      .def_property_readonly("delivery_delta", int_reader(&PlacementResult::delivery_delta),
                    "the delivery after the call minus before it (call_delivered for a placement; may be negative\n"
                    "for a full-offer call)")
      .def_property_readonly(
          "call_edge_drops", [](const PlacementResult& r) { return pairs_list(r.call_edge_drops); },
          "this call's drops on links [(EdgeId, units)], ascending EdgeId (empty unless LOSSY)")
      .def_property_readonly("call_source_drops", int_reader(&PlacementResult::call_source_drops),
                    "this call's offered units lost at sources without a usable route (LOSSY)")
      .def_property_readonly("termination", field_object(&PlacementResult::termination), "the call's outcome (Termination)")
      .def_readonly("routes_changed", &PlacementResult::routes_changed,
                    "True when the call's program differs from the installed one in its ordered heads, paths, split\n"
                    "tables or duplicates (no program counts as the empty program; two programs without an entry\n"
                    "are equal). An offer-only change (forward, Pinned, ShortestPaths increments) and an identical\n"
                    "replace set nothing")
      .def_readonly("weights_changed", &PlacementResult::weights_changed,
                    "True when the call's program differs from the installed one in its injection (source weights)\n"
                    "or any entry or split-table weight (a controller's re-solve that changes a path's amount; an\n"
                    "absent entry weighs 0)")
      .def_readonly("reservations_changed", &PlacementResult::reservations_changed,
                    "True when a tunnel's requested or reserved amount differs from the installed program's, by\n"
                    "entry_id (an absent tunnel holds 0): only T1 sets it; forward and Pinned never do")
      .def_property_readonly(

          "program",
          [](const PlacementResult& r) { return optional_shared(r.program); },
          "the installed ForwardingProgram after the call (evaluate: the computed one); None before the first\n"
          "preparation")
      .def_property_readonly(

          "plan", [](const PlacementResult& r) { return optional_of(r.plan); },
          "T1 and T2: the TunnelPlan of the program after the call; None for every other model (Pinned over\n"
          "a tunnel program too)")
      .def_property_readonly(

          "requested",
          [](const PlacementResult& r) -> py::object {
            return r.requested.has_value() ? counter_int(*r.requested)
                                           : py::object(float_object(std::numeric_limits<double>::infinity()));
          },
          "the sum of the epoch's requests; float('inf') after an infinite one")
      .def_property_readonly(

          "offered", [](const PlacementResult& r) { return counter_or_none(r.offered); },
          "LOSSY: the traffic sent in the epoch (the first-hop injection of an infinite request); None\n"
          "when LOSSLESS")
      .def_property_readonly(

          "admitted", [](const PlacementResult& r) { return counter_int(r.admitted); },
          "LOSSLESS: the total carried; LOSSY: offered")
      .def_property_readonly("delivered", int_reader(&PlacementResult::delivered), "the units the demand delivers")
      .def_property_readonly(

          "unadmitted", [](const PlacementResult& r) { return counter_or_none(r.unadmitted); },
          "LOSSLESS: requested - delivered, None after an infinite request; LOSSY: 0")
      .def_property_readonly(
          "source_drops", [](const PlacementResult& r) { return counter_int(r.source_drops); },
          "LOSSY: the offered units lost at sources (no usable route, or an empty roster); 0 when LOSSLESS")
      .def_property_readonly(
          "edge_drops", [](const PlacementResult& r) { return pairs_list(r.edge_drops); },
          "the epoch's drops on links [(EdgeId, units)], ascending EdgeId (empty unless LOSSY)")
      .def_property_readonly(

          "by_source", [](const PlacementResult& r) { return tuple_of(r.by_source); },
          "one SourceFlow per normalized input source, ascending NodeId (zero rows included)")
      .def_property_readonly(

          "by_tunnel", [](const PlacementResult& r) { return tuple_of(r.by_tunnel); },
          "one TunnelFlow per program tunnel, ascending entry_id (empty without tunnels)")
      .def_property_readonly(

          "total_cost", [](const PlacementResult& r) { return counter_or_none(r.total_cost); },
          "controllers (CentralizedTE, ShortestPaths OPTIMIZED): the exact sum of cost * flow over the policy's\n"
          "allocation; None for the other models (T1, T2, Pinned)")
      .def_property_readonly(
          "marginal_cost_curve",
          [](const PlacementResult& r) {
            using Curve = py::typing::List<py::typing::Tuple<py::int_, py::int_>>;
            return py::reinterpret_steal<py::typing::Optional<Curve>>(
                r.marginal_cost_curve ? pairs_list<Curve>(*r.marginal_cost_curve).release() : py::none().release());
          },
                    "unrestricted CentralizedTE with the optimized source allocation: the marginal cost curve [(cost, units)] of the policy's\n"
                    "flow as the call's solve computed it (a positive place: the cumulative flow; evaluate and replace,\n"
                    "a zero offer included: the whole offer); place(0) keeps the stored curve (None before the epoch's\n"
                    "first solve or forward; after a forward, None for a nonempty roster and [] for an empty one);\n"
                    "a forward with a nonempty roster solves nothing and reports None; with an empty roster every\n"
                    "call but place(0) reports []; in ascending cost, delivering `units` more raised the minimum total flow cost by `cost`\n"
                    "per unit (more delivery can reroute earlier flow, so these are not the costs of the final\n"
                    "paths). None for the other models and prescribed (EQUAL or WEIGHTED) sources")
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
          "edge_reservations",
          [](const PlacementResult& r) -> py::object {
            return r.edge_reservations ? py::object(int64_array(*r.edge_reservations)) : py::object(py::none());
          },
          "if requested: the policy's own strict holds per edge after the call (int64[E], EdgeId order; zeros\n"
          "for a model without holds), else None")
      .def_property_readonly(
          "saturated_boundary_edges",
          [](const PlacementResult& r) {
            return py::reinterpret_steal<py::typing::Optional<py::typing::List<py::int_>>>(
                r.saturated_boundary_edges ? ints_list(std::span<const EdgeId>(*r.saturated_boundary_edges)).release() : py::none().release());
          },
                    "if requested: the ledger's saturation report after the call\n"
                    "(FlowLedger.saturated_boundary_edges of the "
                    "policy's sources and Eligibility), else None")
      .def_property_readonly(
          "min_cut",
          [](const PlacementResult& r) {
            return py::reinterpret_steal<py::typing::Optional<py::typing::List<py::int_>>>(
                r.min_cut ? ints_list(std::span<const EdgeId>(*r.min_cut)).release() : py::none().release());
          },
                    "if requested (OutputOptions(min_cut=True)), unrestricted CentralizedTE with the optimized source allocation, when the solve\n"
                    "reached the maximum flow: the minimum cut closest to the sources, every allowed edge leaving\n"
                    "the nodes the roster sources reach in the residual network (free capacity plus the policy's\n"
                    "own allocation), an edge with nothing free (capacity 0, or held by others) included,\n"
                    "ascending EdgeId; [] for an empty roster, place(0) excepted. None otherwise (another model,\n"
                    "prescribed EQUAL or WEIGHTED sources, a support, a solve stopped below the maximum, place(0),\n"
                    "a forward with a nonempty roster)")
      .def("__repr__", [](const PlacementResult& r) {
        return "PlacementResult(call_delivered=" + std::to_string(r.call_delivered) + ", termination=" +
               termination_name(r.termination) + ")";
      });
}

} // namespace netgraph::core::bindings
