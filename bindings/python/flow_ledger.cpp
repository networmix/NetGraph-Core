// FlowLedger: the integer ledger of one graph.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>

#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/graph.hpp"
#include "integers.hpp"
#include "values.hpp"
#if NETGRAPH_CORE_PRIVATE_HOOKS
#include "netgraph/core/detail/ledger_access.hpp"
#endif
#include "publication.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

void register_flow_ledger(py::module_& m) {
  [[maybe_unused]] auto ledger_class = py::class_<FlowLedger, std::shared_ptr<FlowLedger>>(
      m, "FlowLedger",
      "The integer ledger of one graph, in two projections that are never added: carried traffic (the\n"
      "service baseline; allocated() is the sum of every record's edges and residual() = baseline -\n"
      "allocated()) and the strict holds of T1 tunnels (the reservation baseline; reserved() is the sum of\n"
      "every record's reserved_edges and unreserved() = reservation baseline - reserved()), all int64[E]\n"
      "snapshots. Traffic may use bandwidth that is held: holds never police. The baselines are replaced\n"
      "only by reset(); the records are the flows placed on it by policies (the ledger has no public\n"
      "placement, removal or reservation). reset() and copy() are rejected while a policy is bound (a\n"
      "FlowPolicy built on this ledger still exists).")
      .def(py::init([](const std::shared_ptr<Graph>& graph, const py::object& residual, const py::object& unreserved) {
             if (!graph) throw py::value_error("graph must be a Graph");
             const auto caps = graph->capacity_units();
             const auto r = residual.is_none() ? std::vector<Amount>(caps.begin(), caps.end())
                                               : edge_units(residual, *graph, "residual");
             const auto u = unreserved.is_none() ? std::vector<Amount>(caps.begin(), caps.end())
                                                 : edge_units(unreserved, *graph, "unreserved");
             return std::make_shared<FlowLedger>(graph, r, u);
           }),
           py::arg("graph"), py::arg("residual") = py::none(), py::kw_only(), py::arg("unreserved") = py::none(),
           "residual and unreserved: the service and reservation baselines as integer units, one entry per\n"
           "edge in [0, capacity_units[e]] (a physical baseline goes through to_units first; residual checked\n"
           "first); None: the capacities. ValueError otherwise.")
      .def("uid", &FlowLedger::uid, "process-unique ledger uid (the first half of every FlowId)")
      .def(
          "graph", [](const FlowLedger& fg) { return typed_shared(fg.graph()); },
          "the ledger's graph")
      .def(
          "capacity_units",
          [](const py::object& self) { return readonly_view(self.cast<const FlowLedger&>().capacity_units(), self); },
          "int64[E] the graph's capacities in units (read-only view)")
      .def(
          "residual", [](const FlowLedger& fg) { return int64_array(fg.residual()); },
          "int64[E] snapshot: baseline - allocated")
      .def(
          "allocated", [](const FlowLedger& fg) { return int64_array(fg.allocated()); },
          "int64[E] snapshot: the sum of every record's edges")
      .def(
          "reserved", [](const FlowLedger& fg) { return int64_array(fg.reserved()); },
          "int64[E] snapshot: strict holds per edge (T1), the sum of every record's reserved_edges")
      .def(
          "unreserved", [](const FlowLedger& fg) { return int64_array(fg.unreserved()); },
          "int64[E] snapshot: reservation baseline minus strict holds, per edge")
      .def(
          "records", [](const FlowLedger& fg) { return list_of(fg.records()); },
          "live records as snapshots in FlowId order: those of bound policies and orphans (records left by a\n"
          "policy destroyed without remove())")
      .def(
          "record",
          [](const FlowLedger& fg, const py::object& flow_id) { return typed_object(fg.record(flow_id_of(flow_id))); },
          py::arg("flow_id"),
          "one record by FlowId (ValueError for a FlowId of another ledger, a retired or an unknown id)")
      .def(
          "bottleneck",
          [](const FlowLedger& fg, const py::object& sources, const std::shared_ptr<Eligibility>& elig) {
            const auto s = node_ids(sources, fg.graph()->num_nodes(), "sources");
            return fg.bottleneck(s, elig.get());
          },
          py::arg("sources"), py::arg("elig") = nullptr,
          "Saturation report: the eligible edges with capacity >= 1 unit that leave the set of nodes the\n"
          "allowed sources can reach (following eligible edges forward with residual >= 1 unit), in\n"
          "ascending EdgeId. It ignores the destinations, so it is not a source-destination cut.")
      .def(
          "reset",
          [](FlowLedger& fg, const py::object& residual, const py::object& unreserved) {
            const auto caps = fg.capacity_units();
            const auto r = residual.is_none() ? std::vector<Amount>(caps.begin(), caps.end())
                                              : edge_units(residual, *fg.graph(), "residual");
            const auto u = unreserved.is_none() ? std::vector<Amount>(caps.begin(), caps.end())
                                                : edge_units(unreserved, *fg.graph(), "unreserved");
            fg.reset(r, u);
          },
          py::arg("residual") = py::none(), py::kw_only(), py::arg("unreserved") = py::none(),
          "New baselines (None: the capacities); orphaned records with their allocations and holds are\n"
          "cleared and their ids retired (ids are never reused). ValueError while a policy is bound.")
      .def(
          "copy", [](const FlowLedger& fg) { return std::make_shared<FlowLedger>(fg); },
          "a new ledger (new uid) with the same baselines, records and holds, ids kept under the new uid;\n"
          "ValueError while a policy is bound")
      .def("__copy__", [](const FlowLedger& fg) { return std::make_shared<FlowLedger>(fg); })
      .def("__deepcopy__", [](const FlowLedger& fg, const py::object&) { return std::make_shared<FlowLedger>(fg); },
           py::arg("memo"));
#if NETGRAPH_CORE_PRIVATE_HOOKS
  ledger_class
      .def(
          "_state",
          [](const FlowLedger& fg) {
            py::dict out;
            out["next_flow_id"] = detail::LedgerAccess::next_flow_id(fg);
            out["bound"] = detail::LedgerAccess::bound_owners(fg);
            out["orphans"] = detail::LedgerAccess::orphan_owners(fg);
            return out;
          },
          "private: {next_flow_id, bound, orphans} (owner ids)");
#endif
}

} // namespace netgraph::core::bindings
