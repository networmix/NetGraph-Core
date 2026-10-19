// FlowRecord: one flow of a ledger, as a snapshot.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>

#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/pred_dag.hpp"
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
      "One flow of a ledger, as a snapshot: flow_id (ledger uid, id), owner (the policy's id),\n"
      "sources, destinations, admission, entry (a T1, T2 or Pinned policy's program entry, else None), bundle\n"
      "(PredDAG: the flow's route),\n"
      "route_cost (its minimum walk cost; None when empty), offered (Lossy only, else None),\n"
      "delivered, edges [(EdgeId, units)] (carried load) and drops [(EdgeId, units)] (cumulative,\n"
      "ascending EdgeId), reserved_edges [(EdgeId, units)] (T1's strict holds, ascending EdgeId). offered\n"
      "and the drops units are exact Python ints (they may exceed 64 bits).")
      .def_property_readonly(
          "flow_id", [](const FlowRecord& r) { return py::make_tuple(r.flow_id.first, r.flow_id.second); },
          "(ledger uid, id): ids are assigned in creation order and never reused in the ledger")
      .def_readonly("owner", &FlowRecord::owner, "the ledger-local id of the policy that created the flow")
      .def_readonly("sources", &FlowRecord::sources, "the flow's source node ids, ascending")
      .def_readonly("destinations", &FlowRecord::destinations, "the flow's destination node ids, ascending")
      .def_readonly("admission", &FlowRecord::admission, "LOSSLESS or LOSSY")
      .def_readonly("entry", &FlowRecord::entry,
                    "the program entry the record carries (T1, T2 and Pinned: the entry identifier); None for the\n"
                    "records kept per roster source (CostOnly, Global)")
      .def_property_readonly(
          "bundle", [](const FlowRecord& r) { return typed_shared(r.bundle); },
          "the flow's route (an immutable PredDAG shared with the ledger)")
      .def_property_readonly(
          "route_cost", [](const FlowRecord& r) { return cost_or_none(r.route_cost); },
          "the cost of the cheapest walk of the bundle from the sources to the destinations; None when\n"
          "the bundle has no such walk")
      .def_property_readonly(
          "offered",
          [](const FlowRecord& r) -> py::object {
            return r.offered.has_value() ? counter_int(*r.offered) : py::object(py::none());
          },
          "Lossy: the units offered to the flow, summed over calls (an exact int, may exceed 64 bits);\n"
          "None for a Lossless flow")
      .def_readonly("delivered", &FlowRecord::delivered, "the units the flow delivers to its destinations (its total)")
      .def_property_readonly(
          "edges", [](const FlowRecord& r) { return pairs_list(r.edges); },
          "the flow's allocation [(EdgeId, units)], positive entries in ascending EdgeId")
      .def_property_readonly(
          "drops",
          [](const FlowRecord& r) {
            py::list out;
            for (const auto& [e, x] : r.drops) out.append(py::make_tuple(e, counter_int(x)));
            return out;
          },
          "Lossy: the units dropped per edge [(EdgeId, units)], summed over calls, positive entries in\n"
          "ascending EdgeId (exact ints); empty for a Lossless flow")
      .def_property_readonly(
          "reserved_edges", [](const FlowRecord& r) { return pairs_list(r.reserved_edges); },
          "T1: the strict holds of an up tunnel [(EdgeId, units)], its request on every link of its path,\n"
          "ascending EdgeId; empty for every other record (holds never police traffic)")
      .def("__repr__", [](const FlowRecord& r) {
        return "FlowRecord(flow_id=(" + std::to_string(r.flow_id.first) + ", " + std::to_string(r.flow_id.second) +
               "), owner=" + std::to_string(r.owner) + ", delivered=" + std::to_string(r.delivered) + ")";
      });
}

} // namespace netgraph::core::bindings
