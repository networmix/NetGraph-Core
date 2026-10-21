// FlowLedger: the integer ledger of one graph.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>

#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/graph.hpp"
#include "integers.hpp"
#include "values.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "native_call.hpp"
#include "publication.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

void register_flow_ledger(py::module_& m) {
  [[maybe_unused]] auto ledger_class = py::class_<FlowLedger, std::shared_ptr<FlowLedger>>(
      m, "FlowLedger",
      "The integer ledger of one graph, in two projections that are never added: carried traffic (the\n"
      "service baseline; carried() is the sum of every record's edge_flows and residual() = baseline -\n"
      "carried()) and the strict holds of T1 tunnels (the reservation baseline; reserved() is the sum of\n"
      "every record's edge_reservations and unreserved() = reservation baseline - reserved()), all int64[E]\n"
      "snapshots. Traffic may use bandwidth that is held: holds never police. The baselines are replaced\n"
      "only by reset(); the records are the flows placed on it by policies (the ledger has no public\n"
      "placement, removal or reservation). reset() and copy() are rejected while a policy is bound (a\n"
      "FlowPolicy built on this ledger still exists).")
      .def(py::init([](const std::shared_ptr<Graph>& graph, const py::object& service_baseline, const py::object& reservation_baseline) {
             if (!graph) throw py::type_error("graph must be a Graph");
             const auto caps = graph->edge_capacity_units();
             const auto r = service_baseline.is_none() ? std::vector<Amount>(caps.begin(), caps.end())
                                               : edge_units(service_baseline, *graph, "service_baseline");
             const auto u = reservation_baseline.is_none() ? std::vector<Amount>(caps.begin(), caps.end())
                                                 : edge_units(reservation_baseline, *graph, "reservation_baseline");
             return native_call([&] { return std::make_shared<FlowLedger>(graph, r, u); });
           }),
           py::arg("graph"), py::arg("service_baseline") = py::none(), py::kw_only(), py::arg("reservation_baseline") = py::none(),
           "service_baseline and reservation_baseline: the baselines as integer units, one entry per edge in\n"
           "[0, the edge's capacity] (a physical baseline goes through to_units first; service_baseline\n"
           "checked first); None: the capacities. A wrong type (not a 1-D sequence of integers, a non-integer\n"
           "entry) raises TypeError, a wrong length or an entry outside its range ValueError.")
      .def("uid", [](const FlowLedger& l) { return native_call([&] { return l.uid(); }); }, "process-unique ledger uid (the first half of every FlowId)")
      .def(
          "graph", [](const FlowLedger& fg) { auto graph = native_call([&] { return std::shared_ptr<const Graph>(fg.graph()); }); return typed_shared(std::move(graph)); },
          "the ledger's graph")
      .def(
          "edge_capacity_units",
          [](const py::object& self) {
            const auto& ledger = self.cast<const FlowLedger&>();
            const auto capacities = native_call([&] { return ledger.edge_capacity_units(); });
            return readonly_view(capacities, self);
          },
          "int64[E] the graph's capacities in units (read-only view)")
      .def(
          "residual", [](const FlowLedger& fg) { return int64_array(native_call([&] { return fg.residual(); })); },
          "int64[E] snapshot: service_baseline - carried()")
      .def(
          "carried", [](const FlowLedger& fg) { return int64_array(native_call([&] { return fg.carried(); })); },
          "int64[E] snapshot: the sum of every record's edge_flows")
      .def(
          "reserved", [](const FlowLedger& fg) { return int64_array(native_call([&] { return fg.reserved(); })); },
          "int64[E] snapshot: strict holds per edge (T1), the sum of every record's edge_reservations")
      .def(
          "unreserved", [](const FlowLedger& fg) { return int64_array(native_call([&] { return fg.unreserved(); })); },
          "int64[E] snapshot: reservation baseline minus strict holds, per edge")
      .def(
          "records", [](const FlowLedger& fg) { return list_of(native_call([&] { return fg.records(); })); },
          "live records as snapshots in FlowId order: those of bound policies and orphans (records left by a\n"
          "policy destroyed without remove())")
      .def(
          "record",
          [](const FlowLedger& fg, const py::object& flow_id) { const auto id = flow_id_of(flow_id); return typed_object(native_call([&] { return fg.record(id); })); },
          py::arg("flow_id"),
          "one record by FlowId (ValueError for a FlowId of another ledger, a retired or an unknown id)")
      .def(
          "saturated_boundary_edges",
          [](const FlowLedger& fg, const py::object& sources, const std::shared_ptr<Eligibility>& elig) {
            const auto s = node_ids(sources, detail::LedgerAccess::graph(fg)->num_nodes(), "sources");
            if (elig && elig->graph_uid() != detail::LedgerAccess::graph(fg)->uid()) throw py::value_error("foreign Eligibility: built for another graph");
            return native_call([&] { return fg.saturated_boundary_edges(s, elig.get()); });
          },
          py::arg("sources"), py::arg("eligibility") = nullptr,
          "Saturation report: the eligible edges with capacity >= 1 unit that leave the set of nodes the\n"
          "allowed sources can reach (following eligible edges forward with residual >= 1 unit), in\n"
          "ascending EdgeId. It ignores the destinations, so it is not a source-destination cut.")
      .def(
          "reset",
          [](FlowLedger& fg, const py::object& service_baseline, const py::object& reservation_baseline) {
            const auto& graph = detail::LedgerAccess::graph(fg);
            const auto caps = graph->edge_capacity_units();
            const auto r = service_baseline.is_none() ? std::vector<Amount>(caps.begin(), caps.end())
                                              : edge_units(service_baseline, *graph, "service_baseline");
            const auto u = reservation_baseline.is_none() ? std::vector<Amount>(caps.begin(), caps.end())
                                                : edge_units(reservation_baseline, *graph, "reservation_baseline");
            native_call([&] { fg.reset(r, u); });
          },
          py::arg("service_baseline") = py::none(), py::kw_only(), py::arg("reservation_baseline") = py::none(),
          "New baselines (None: the capacities); orphaned records with their allocations and holds are\n"
          "cleared and their ids retired (ids are never reused). ValueError while a policy is bound.")
      .def(
          "copy", [](const FlowLedger& fg) { return native_call([&] { return std::make_shared<FlowLedger>(fg); }); },
          "a new ledger (new uid) with the same baselines, records and holds, ids kept under the new uid;\n"
          "ValueError while a policy is bound")
      .def("__copy__", [](const FlowLedger& fg) { return native_call([&] { return std::make_shared<FlowLedger>(fg); }); })
      .def("__deepcopy__", [](const FlowLedger& fg, const py::object&) { return native_call([&] { return std::make_shared<FlowLedger>(fg); }); },
           py::arg("memo"));
#if NETGRAPH_CORE_PRIVATE_HOOKS
  ledger_class
      .def(
          "_state",
          [](const FlowLedger& fg) {
            struct Snapshot {
              std::int64_t next;
              std::vector<std::int64_t> bound, orphans;
            };
            const auto state = native_call([&] {
              const auto gate = detail::LedgerAccess::lock(fg);
              return Snapshot {detail::LedgerAccess::next_record_counter(fg),
                  detail::LedgerAccess::bound_policy_ids(fg), detail::LedgerAccess::orphan_policy_ids(fg)};
            });
            py::dict out;
            out["next_flow_id"] = state.next;
            out["bound"] = state.bound;
            out["orphans"] = state.orphans;
            return out;
          },
          "private: {next_flow_id, bound, orphans} (owner ids)");
#endif
}

} // namespace netgraph::core::bindings
