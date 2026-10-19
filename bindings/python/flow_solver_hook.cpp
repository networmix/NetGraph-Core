// The private entry _solve_min_cost_flow, registered in both builds: the flow core's
// minimum-cost flow from zero on a given view, with its full payload for the flow checks.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/min_cost_flow.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/graph.hpp"
#include "integers.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

py::list int_list(const std::vector<std::int64_t>& v) {
  py::list out(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) out[i] = py::int_(v[i]);
  return out;
}

} // namespace

void register_flow_solver_hook(py::module_& m) {
  m.def(
      "_solve_min_cost_flow",
      [](const std::shared_ptr<Graph>& graph, const py::object& view, const py::object& sources,
         const py::object& destinations, const std::shared_ptr<Eligibility>& elig,
         const py::object& cap) -> py::dict {
        if (!graph) throw py::value_error("graph must be a Graph");
        // The core's order (view length, view entries, sources, destinations, overlap, cap,
        // Eligibility), every value checked before narrowing; the core re-checks the same order.
        const auto v = edge_amounts(view, *graph, "view", "view (capacity units)");
        const auto s = checked_ids(wide_integers(sources, "sources node ids"), graph->num_nodes(), "sources");
        const auto d = checked_ids(wide_integers(destinations, "destinations node ids"), graph->num_nodes(), "destinations");
        detail::check_disjoint(detail::node_set(s, graph->num_nodes(), "sources"),
                               detail::node_set(d, graph->num_nodes(), "destinations"), "sources", "destinations");
        std::optional<Amount> c;
        if (!cap.is_none()) {
          if (is_bool(cap) || !PyIndex_Check(cap.ptr())) {
            throw py::value_error("cap must be an integer number of units in [0, 2^63 - 2] or None");
          }
          py::object idx = python_int(cap, "cap");
          int overflow = 0;
          const long long value = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
          if (overflow != 0 || value < 0 || value > detail::kMaxFlowCap) {
            throw py::value_error("cap " + std::string(py::str(idx)) + " is outside [0, 2^63 - 2]");
          }
          c = static_cast<Amount>(value);
        }
        detail::MinCostFlowResult r;
        {
          py::gil_scoped_release release;
          const auto view_units = v.narrow<Amount>();
          r = detail::solve_min_cost_flow(*graph, view_units, s, d, elig.get(), c);
        }
        py::dict out;
        out["flow_value"] = py::int_(r.flow_value);
        py::list curve;
        for (const auto& [cost, units] : r.marginal_cost_curve) {
          py::list entry;
          entry.append(py::int_(cost));
          entry.append(py::int_(units));
          curve.append(entry);
        }
        out["marginal_cost_curve"] = curve;
        out["edge_flows"] = int_list(r.edge_flows);
        if (r.min_cut.has_value()) {
          py::list cut;
          for (EdgeId e : *r.min_cut) cut.append(py::int_(e));
          out["min_cut"] = cut;
        } else {
          out["min_cut"] = py::none();
        }
        out["potentials"] = int_list(r.potentials);
        py::dict counters;
        counters["dijkstra_runs"] = py::int_(r.counters.dijkstra_runs);
        counters["dijkstra_pops"] = py::int_(r.counters.dijkstra_pops);
        counters["settled"] = py::int_(r.counters.settled);
        counters["arc_scans"] = py::int_(r.counters.arc_scans);
        counters["dinic_phases"] = py::int_(r.counters.dinic_phases);
        counters["augmentations"] = py::int_(r.counters.augmentations);
        out["counters"] = counters;
        return out;
      },
      py::arg("graph"), py::arg("view"), py::arg("sources"), py::arg("destinations"), py::arg("eligibility"),
      py::arg("cap"),
      "Private flow-core entry: the largest flow up to cap (None: no cap) from the source set to the\n"
      "destination set over the eligible edges, sent at the lowest total cost, on the capacities `view`\n"
      "(int64 units, EdgeId order), starting from zero flow. Returns {flow_value, marginal_cost_curve:\n"
      "[[cost, units]], edge_flows, min_cut (None while the destinations are still reachable, because\n"
      "the cap stopped the solve, and on cap 0), potentials (per-node prices from the shortest-path runs\n"
      "that prove the flow is cheapest; -1 = never reached), counters}. Invalid inputs raise ValueError\n"
      "before any work. A solve (cap not 0) then raises ValueError, before the solver sizes anything, on\n"
      "a network beyond the flow core's int32 numbering: more than INT32_MAX - 2 nodes, or more than\n"
      "2^30 - 1 eligible edges, eligible sources and eligible destinations together.");
}

} // namespace netgraph::core::bindings
