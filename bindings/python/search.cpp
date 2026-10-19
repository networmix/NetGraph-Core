// spf with SpfOptions and SpfResult, the public search, and UNREACHABLE.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/search.hpp"
#include "netgraph/core/types.hpp"
#include "integers.hpp"
#include "publication.hpp"
#include "values.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// The Python form of SpfOptions: the arguments as given (targets, residual and min_residual are
// validated at the call, in spf's validation order, on their unnarrowed values); the Eligibility is
// held by reference count.
struct PySpfOptions {
  Direction direction;
  PathSelection selection;
  py::object targets;
  py::object residual;
  py::object min_residual;
  std::shared_ptr<Eligibility> elig;
};

// The Python form of SpfResult: every field already a Python object (made by the call through
// the checked allocations), so reading one allocates nothing.
struct PySpfResult {
  py::object costs;        // a read-only int64 array (full searches) or None
  py::object dag;          // PredDAG
  py::object target_cost;  // an int, or None (no target reached, or a full search)
};

} // namespace

void register_search(py::module_& m) {
  m.attr("UNREACHABLE") = py::int_(kUnreachable);

  py::class_<PySpfOptions>(
      m, "SpfOptions",
      "The options of spf: direction and selection are explicit; targets (one NodeId or a sequence;\n"
      "None: a full search) stops the search at the nearest targets (anycast: any nearest target will\n"
      "do); residual (int64 units per EdgeId, each in [0, capacity_units]) makes the search skip edges\n"
      "with fewer than min_residual residual units (an integer >= 1; without a residual it has no\n"
      "effect) and supplies the priority values: it is the only capacity filter, and without it an edge\n"
      "of capacity 0 is searched like any other; elig (None: everything allowed). A direction outside its\n"
      "enumeration raises ValueError here; the other fields are validated by spf, in the order its\n"
      "docstring gives.")
      .def(py::init([](Direction direction, const PathSelection& selection, const py::object& targets,
                       const py::object& residual, const py::object& min_residual,
                       const std::shared_ptr<Eligibility>& elig) {
             detail::check_member(direction, "direction");
             return PySpfOptions {direction, selection, targets, residual, min_residual, elig};
           }),
           py::arg("direction"), py::arg("selection"), py::kw_only(), py::arg("targets") = py::none(),
           py::arg("residual") = py::none(), py::arg("min_residual") = 1, py::arg("elig") = nullptr)
      .def_readonly("direction", &PySpfOptions::direction, "FORWARD or REVERSE")
      .def_readonly("selection", &PySpfOptions::selection, "the PathSelection the search reads")
      .def_readonly("targets", &PySpfOptions::targets, "the targets as given (validated by spf); None: a full search")
      .def_readonly("residual", &PySpfOptions::residual,
                    "the residual units per EdgeId as given (validated by spf); None: no capacity filter;\n"
                    "priorities read static capacity")
      .def_readonly("min_residual", &PySpfOptions::min_residual,
                    "the gating threshold in units as given (validated by spf: an integer >= 1)")
      .def_readonly("elig", &PySpfOptions::elig, "the static Eligibility; None: everything allowed");

  py::class_<PySpfResult>(m, "SpfResult",
                          "The result of spf: costs (full searches: a read-only int64 array, UNREACHABLE for the nodes not\n"
                          "reached; None with targets), dag (PredDAG of the tight edges, the edges on some shortest walk,\n"
                          "one per (u, v) pair with multi_edge = False, each in its own direction; with targets only the\n"
                          "walks between the roots and a nearest target, for REVERSE from the target to a root; for\n"
                          "multipath = False the canonical single path) and\n"
                          "target_cost (with targets: the nearest targets' cost, None when no target is reached; None for a\n"
                          "full search).")
      .def_property_readonly(
          "costs", [](const PySpfResult& r) { return r.costs ? r.costs : py::object(py::none()); },
          "full search: each node's distance from the roots for FORWARD, or to the roots for REVERSE\n"
          "(read-only int64[N], UNREACHABLE for nodes not reached); None with targets")
      .def_property_readonly(
          "dag", [](const PySpfResult& r) { return r.dag; },
          "the PredDAG of the result, in the edges' original orientation")
      .def_property_readonly(
          "target_cost", [](const PySpfResult& r) { return r.target_cost ? r.target_cost : py::object(py::none()); },
          "with targets: the cost of the nearest targets, None when no target is reached; None for a\n"
          "full search")
      .def("__repr__", [](const PySpfResult& r) {
        std::string tc = r.target_cost && !r.target_cost.is_none() ? std::string(py::str(r.target_cost)) : "None";
        return "SpfResult(costs=" + std::string(r.costs ? "array" : "None") +
               ", dag_entries=" + std::to_string(r.dag.cast<const PredDAG&>().via_edges().size()) +
               ", target_cost=" + tc + ")";
      });

  m.def(
      "spf",
      [](const std::shared_ptr<Graph>& graph, const py::object& roots, const py::object& opts) {
        if (!graph) throw py::value_error("graph must be a Graph");
        const auto o = options_of<PySpfOptions>(opts, "SpfOptions");
        const auto n = graph->num_nodes();
        // spf's validation order on unnarrowed values (Core repeats every check): the enumerations,
        // the roots, the targets, their overlap, a single path without targets, the residual, then
        // min_residual and a foreign Eligibility.
        detail::check_member(o.direction, "direction");
        detail::check_selection(o.selection);
        const auto r = node_ids(roots, n, "roots");
        std::optional<std::vector<NodeId>> t;
        if (!o.targets.is_none()) {
          t = node_ids(o.targets, n, "targets");
          detail::check_disjoint(detail::node_set(r, n, "roots"), detail::node_set(*t, n, "targets"), "roots",
                                 "targets");
        } else if (o.selection.single_path()) {
          throw py::value_error("multipath = false requires targets (the single path joins a root to a target)");
        }
        std::optional<std::vector<Amount>> res;
        if (!o.residual.is_none()) res = edge_units(o.residual, *graph, "residual");
        const auto minr = python_int(o.min_residual, "min_residual");
        int overflow = 0;
        const long long mv = PyLong_AsLongLongAndOverflow(minr.ptr(), &overflow);
        if (overflow < 0 || (overflow == 0 && mv < 1)) {
          throw py::value_error("min_residual must be at least 1 unit (got " + std::string(py::str(minr)) + ")");
        }
        SpfOptions c(o.direction, o.selection);
        c.targets = std::move(t);
        if (res.has_value()) c.residual = std::span<const Amount>(*res);
        // Above INT64_MAX it gates every residual edge exactly as INT64_MAX does (residuals < 2^62).
        c.min_residual = overflow != 0 ? std::numeric_limits<Amount>::max() : static_cast<Amount>(mv);
        c.elig = o.elig.get();
        std::optional<SpfResult> result;
        {
          py::gil_scoped_release release;
          result.emplace(spf(*graph, r, c));
        }
        // The result's Python objects, each through a checked allocation: a failure raises.
        auto [obj, slot] = new_object<PySpfResult>();
        if (result->costs.has_value()) {
          py::array_t<std::int64_t> costs(static_cast<py::ssize_t>(result->costs->size()));
          std::copy(result->costs->begin(), result->costs->end(), costs.mutable_data());
          py::detail::array_proxy(costs.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
          slot->costs = std::move(costs);
        }
        slot->dag = object_from<PredDAG>(std::move(result->dag));
        if (result->target_cost.has_value() && *result->target_cost != kUnreachable) {
          slot->target_cost = py::int_(*result->target_cost);
        }
        return std::move(obj);
      },
      py::arg("graph"), py::arg("roots"), py::arg("opts"),
      "Shortest paths from the roots for FORWARD, or to the roots for REVERSE (one NodeId or a sequence;\n"
      "every root at distance 0), over the eligible edges (only those with residual >= min_residual when\n"
      "a residual is given; an edge of capacity 0 is searched like any other) -> SpfResult. Without\n"
      "targets every reachable node gets its distance: costs and the full DAG of tight edges (the edges\n"
      "on some shortest walk; one per (u, v) pair with multi_edge = False). With targets the search\n"
      "stops at the nearest targets' cost:\n"
      "target_cost and the DAG of the walks between the roots and a nearest target (for REVERSE, from\n"
      "the target to a root, in the edges' own direction);\n"
      "with multipath = False, instead the canonical single path, the walk minimising (cost, -bottleneck\n"
      "under PREFER_HIGHER_PRIORITY, hops, EdgeId sequence). Validated in this order before any work\n"
      "(ValueError): the roots (non-empty, ids in range), the targets (the same), a node in both sets,\n"
      "multipath = False without targets, the residual (length, then every entry in [0, capacity_units]),\n"
      "min_residual >= 1, a foreign Eligibility. Excluded roots and targets are skipped; all roots\n"
      "excluded gives an empty result.");
}

} // namespace netgraph::core::bindings
