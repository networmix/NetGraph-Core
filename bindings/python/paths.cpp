// The path listings: ksp with KspOptions, and resolve_to_paths.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/traversal.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/search.hpp"
#include "integers.hpp"
#include "publication.hpp"
#include "values.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// One NodeId of [0, num_nodes) (a single endpoint); the range message of detail/endpoints.hpp.
NodeId endpoint(const py::object& obj, std::int32_t num_nodes, const char* name) {
  const py::object idx = python_int(obj, name);
  int overflow = 0;
  const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
  if (overflow != 0 || v < 0 || v >= num_nodes) {
    throw py::value_error(detail::node_range_error(name, std::string(py::str(idx)), num_nodes));
  }
  return static_cast<NodeId>(v);
}

// A count of at least 1 (KSP's k, the resolver's max_paths): the validation message below 1; a count
// above INT64_MAX means the same as INT64_MAX (no list reaches it).
std::int64_t count_of(const py::object& obj, const char* name) {
  const py::object idx = python_int(obj, name);
  int overflow = 0;
  const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
  if (overflow < 0 || (overflow == 0 && v < 1)) {
    throw py::value_error(std::string(name) + " must be >= 1 (got " + std::string(py::str(idx)) + ")");
  }
  return overflow != 0 ? std::numeric_limits<std::int64_t>::max() : static_cast<std::int64_t>(v);
}

// The Python form of KspOptions: the arguments as given (k, max_path_cost and residual are
// validated at the call, in ksp's validation order, on their unnarrowed values).
struct PyKspOptions {
  py::object k;
  PathSelection selection;
  py::object max_path_cost;
  py::object residual;
  std::shared_ptr<Eligibility> elig;
};

} // namespace

void register_paths(py::module_& m) {
  py::class_<PyKspOptions>(
      m, "KspOptions",
      "The options of ksp: k (>= 1) and selection are explicit (multi_edge = True lists each parallel\n"
      "edge as a separate path, False keeps one per (u, v) pair as EdgeTieBreak describes; multipath is\n"
      "ignored); max_path_cost (None, Absolute(c) with c >= 0, or Stretch(f) with f finite in [1, inf)\n"
      "times the static shortest cost under the Eligibility, residual ignored; inclusive); residual\n"
      "(int64 units per EdgeId, each in [0, capacity_units]: edges with less than 1 unit are skipped, and\n"
      "the values feed the priority); elig (None: everything allowed). The fields are validated by ksp,\n"
      "in the order its docstring gives.")
      .def(py::init([](const py::object& k, const PathSelection& selection, const py::object& max_path_cost,
                       const py::object& residual, const std::shared_ptr<Eligibility>& elig) {
             return PyKspOptions {k, selection, max_path_cost, residual, elig};
           }),
           py::arg("k"), py::arg("selection"), py::kw_only(), py::arg("max_path_cost") = py::none(),
           py::arg("residual") = py::none(), py::arg("elig") = nullptr)
      .def_readonly("k", &PyKspOptions::k, "the number of paths as given (validated by ksp: an integer >= 1)")
      .def_readonly("selection", &PyKspOptions::selection, "the PathSelection (multipath is ignored)")
      .def_readonly("max_path_cost", &PyKspOptions::max_path_cost,
                    "the inclusive cost ceiling as given: None, Absolute or Stretch (validated by ksp)")
      .def_readonly("residual", &PyKspOptions::residual,
                    "the residual units per EdgeId as given (validated by ksp); None: static capacities, no\n"
                    "residual gating (edges below 1 unit of capacity are still skipped)")
      .def_readonly("elig", &PyKspOptions::elig, "the static Eligibility; None: everything allowed");

  m.def(
      "ksp",
      [](const std::shared_ptr<Graph>& graph, const py::object& src, const py::object& dst, const py::object& opts) {
        if (!graph) throw py::value_error("graph must be a Graph");
        const auto o = options_of<PyKspOptions>(opts, "KspOptions");
        const auto n = graph->num_nodes();
        // ksp's validation order on unnarrowed values (Core repeats every check): the selection, src,
        // dst, src == dst, k, the ceiling, the residual, then a foreign Eligibility.
        detail::check_selection(o.selection);
        const NodeId s = endpoint(src, n, "src");
        const NodeId t = endpoint(dst, n, "dst");
        if (s == t) throw py::value_error(detail::overlap_error(s, "src", "dst"));
        KspOptions c(count_of(o.k, "k"), o.selection);
        if (!o.max_path_cost.is_none()) {
          if (py::isinstance<Absolute>(o.max_path_cost)) {
            c.max_path_cost = o.max_path_cost.cast<Absolute>();
          } else if (py::isinstance<Stretch>(o.max_path_cost)) {
            c.max_path_cost = o.max_path_cost.cast<Stretch>();
          } else {
            throw py::type_error("max_path_cost must be Absolute, Stretch or None (got " +
                                 std::string(py::str(py::type::of(o.max_path_cost).attr("__name__"))) + ")");
          }
          detail::check_ceiling(c.max_path_cost);
        }
        std::optional<std::vector<Amount>> res;
        if (!o.residual.is_none()) res = edge_units(o.residual, *graph, "residual");
        if (res.has_value()) c.residual = std::span<const Amount>(*res);
        c.elig = o.elig.get();
        std::vector<std::pair<Cost, Path>> paths;
        {
          py::gil_scoped_release release;
          paths = ksp(*graph, s, t, c);
        }
        // Every object through a checked allocation: a failure raises.
        py::list out(paths.size());
        for (std::size_t i = 0; i < paths.size(); ++i) {
          py::tuple item(2);
          py::int_ cost(paths[i].first);
          py::object path = object_from<Path>(std::move(paths[i].second));
          PyTuple_SET_ITEM(item.ptr(), 0, cost.release().ptr());
          PyTuple_SET_ITEM(item.ptr(), 1, path.release().ptr());
          PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), item.release().ptr());
        }
        return out;
      },
      py::arg("graph"), py::arg("src"), py::arg("dst"), py::arg("opts"),
      "The k shortest simple src -> dst paths over the eligible edges of at least 1 unit, as the options\n"
      "further limit them -> [(cost, Path)], each hop a single edge. Order: lowest cost; then, under\n"
      "PREFER_HIGHER_PRIORITY only, highest bottleneck (the smallest priority value along the path); then\n"
      "fewest hops; then the smallest EdgeId sequence. The result is the first k paths of that order\n"
      "within the inclusive cost ceiling. Validated in this order before any work (ValueError): src and\n"
      "dst in range, src == dst, k >= 1, the ceiling, the residual (length, then every entry), a foreign\n"
      "Eligibility. An excluded endpoint, no static route under Stretch, or no path gives [].");

  m.def(
      "resolve_to_paths",
      [](const std::shared_ptr<Graph>& graph, const std::shared_ptr<PredDAG>& dag, const py::object& src,
         const py::object& dst, bool split_parallel_edges, const py::object& max_paths) {
        if (!graph) throw py::value_error("graph must be a Graph");
        if (!dag) throw py::value_error("dag must be a PredDAG");
        if (dag->graph_uid() != graph->uid()) throw py::value_error("foreign PredDAG: built for another graph");
        const auto n = graph->num_nodes();
        const NodeId s = endpoint(src, n, "src");
        const NodeId t = endpoint(dst, n, "dst");
        if (s == t) throw py::value_error(detail::overlap_error(s, "src", "dst"));
        std::optional<std::int64_t> mp;
        if (!max_paths.is_none()) mp = count_of(max_paths, "max_paths");
        std::vector<Path> paths;
        {
          py::gil_scoped_release release;
          paths = resolve_to_paths(*graph, *dag, s, t, split_parallel_edges, mp);
        }
        py::list out(paths.size());
        for (std::size_t i = 0; i < paths.size(); ++i) {
          PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), object_from<Path>(std::move(paths[i])).release().ptr());
        }
        return out;
      },
      py::arg("graph"), py::arg("dag"), py::arg("src"), py::arg("dst"), py::arg("split_parallel_edges") = false,
      py::arg("max_paths") = py::none(),
      "The src -> dst walks of a DAG -> [Path], in increasing (cost, hops, edge-id sequence):\n"
      "a hop costs its cheapest member and is keyed by its lowest member id. split_parallel_edges = False\n"
      "keeps every parallel member of a (u, v) pair in one hop (members ascending); True gives one path\n"
      "per member combination, sorted by that path's own (cost, hops, edge-id sequence). max_paths\n"
      "(>= 1) truncates; None lists every walk.\n"
      "ValueError: a DAG of another graph, src or dst out of range, src == dst, max_paths < 1. An\n"
      "endpoint outside the DAG, or no walk, gives [].");
}

} // namespace netgraph::core::bindings
