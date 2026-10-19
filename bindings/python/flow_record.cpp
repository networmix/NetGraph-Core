// FlowRecord: one flow of a ledger, as a snapshot.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>

#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/route_structures.hpp"
#include "values.hpp"
#include "publication.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

py::object cost_or_none(Cost c) { return c == kUnreachable ? py::object(py::none()) : py::object(py::int_(c)); }

} // namespace

void register_flow_record(py::module_& m) {
  py::class_<FlowRecord>(
      m, "FlowRecord",
      "One flow of a ledger, as a snapshot: flow_id (ledger uid, id), policy_id (the policy's id),\n"
      "sources, destinations, admission, entry_id (a T1, T2 or Pinned policy's program entry, else None),\n"
      "route (PredecessorDAG: the flow's route), min_route_cost (its minimum walk cost; None when empty),\n"
      "offered (Lossy only, else None), delivered, edge_flows [(EdgeId, units)] (carried load) and\n"
      "edge_drops [(EdgeId, units)] (cumulative, ascending EdgeId), edge_reservations [(EdgeId, units)]\n"
      "(T1's strict holds, ascending EdgeId). offered and the edge_drops units are exact Python ints (they\n"
      "may exceed 64 bits).")
      .def_property_readonly(

          "flow_id", [](const FlowRecord& r) { return py::make_tuple(r.flow_id.first, r.flow_id.second); },
          "(ledger uid, id): ids are assigned in creation order and never reused in the ledger")
      .def_readonly("policy_id", &FlowRecord::policy_id, "the ledger-local id of the policy that created the flow")
      .def_readonly("sources", &FlowRecord::sources, "the flow's source node ids, ascending")
      .def_readonly("destinations", &FlowRecord::destinations, "the flow's destination node ids, ascending")
      .def_property_readonly("admission", field_object(&FlowRecord::admission), "LOSSLESS or LOSSY")
      .def_readonly("entry_id", &FlowRecord::entry_id,
                    "the program entry the record carries (T1, T2 and Pinned: the entry identifier); None for the\n"
                    "records kept per roster source (ShortestPaths, CentralizedTE)")
      .def_property_readonly(
          "route", [](const FlowRecord& r) { return typed_shared(r.route); },
          "the flow's route (an immutable PredecessorDAG shared with the ledger)")
      .def_property_readonly(
          "min_route_cost", [](const FlowRecord& r) { return cost_or_none(r.min_route_cost); },
          "the cost of the cheapest walk of the route from the sources to the destinations; None when\n"
          "the route has no such walk")
      .def_property_readonly(

          "offered",
          [](const FlowRecord& r) -> py::object {
            return r.offered.has_value() ? counter_int(*r.offered) : py::object(py::none());
          },
          "Lossy: the units offered to the flow, summed over calls (an exact int, may exceed 64 bits);\n"
          "None for a Lossless flow")
      .def_readonly("delivered", &FlowRecord::delivered, "the units the flow delivers to its destinations (its total)")
      .def_property_readonly(
          "edge_flows", [](const FlowRecord& r) { return pairs_list(r.edge_flows); },
          "the flow's allocation [(EdgeId, units)], positive entries in ascending EdgeId")
      .def_property_readonly(
          "edge_drops",
          [](const FlowRecord& r) {
            py::list out;
            for (const auto& [e, x] : r.edge_drops) out.append(py::make_tuple(e, counter_int(x)));
            return out;
          },
          "Lossy: the units dropped per edge [(EdgeId, units)], summed over calls, positive entries in\n"
          "ascending EdgeId (exact ints); empty for a Lossless flow")
      .def_property_readonly(
          "edge_reservations", [](const FlowRecord& r) { return pairs_list(r.edge_reservations); },
          "T1: the strict holds of an up tunnel [(EdgeId, units)], its request on every link of its path,\n"
          "ascending EdgeId; empty for every other record (holds never police traffic)")
      .def("__repr__", [](const FlowRecord& r) {
        return "FlowRecord(flow_id=(" + std::to_string(r.flow_id.first) + ", " + std::to_string(r.flow_id.second) +
               "), policy_id=" + std::to_string(r.policy_id) + ", delivered=" + std::to_string(r.delivered) + ")";
      });
}

} // namespace netgraph::core::bindings
