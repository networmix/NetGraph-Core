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
#include "netgraph/core/detail/search_detail.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/search.hpp"
#include "integers.hpp"
#include "strict_bool.hpp"
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

// The Python form of KspOptions: the arguments as given (k, cost_ceiling and residual are
// validated at the call, in ksp's validation order, on their unnarrowed values).
struct PyKspOptions {
  py::object k;
  RouteTieBreak tie_break;
  bool parallel_edges;
  py::object cost_ceiling;
  py::object residual;
  std::shared_ptr<Eligibility> eligibility;
};

} // namespace

void register_paths(py::module_& m) {
  py::class_<PyKspOptions>(
      m, "KspOptions",
      "The options of ksp: k (>= 1); tie_break (RouteTieBreak) orders paths of equal cost;\n"
      "parallel_edges = True lists each parallel edge as a separate path, False keeps one per (u, v)\n"
      "pair by the tie rule (no multipath field: k sets the number of paths); cost_ceiling (None,\n"
      "AbsoluteCost(c) with c >= 0, or CostStretch(n/d) with n/d >= 1 times the static shortest\n"
      "cost under the eligibility, compared exactly, residual ignored; inclusive); residual (int64 units per\n"
      "EdgeId, each in [0, edge_capacity_units]: edges with less than 1 unit are skipped, and the values\n"
      "feed the availability; it is the only capacity filter, and without it an edge of capacity 0 is\n"
      "searched like any other); eligibility (None: everything allowed). A parallel_edges that is not a\n"
      "bool raises TypeError here; the other fields are validated by ksp, in the order its docstring\n"
      "gives.")
      .def(py::init([](const py::object& k, RouteTieBreak tie_break, StrictBool parallel_edges,
                       const py::object& cost_ceiling, const py::object& residual,
                       const std::shared_ptr<Eligibility>& eligibility) {
             return PyKspOptions {k, tie_break, parallel_edges, cost_ceiling, residual, eligibility};
           }),
           py::arg("k"), py::kw_only(), py::arg_v("tie_break", RouteTieBreak::FewerHops, "RouteTieBreak.FEWER_HOPS"),
           py::arg("parallel_edges") = StrictBool {true}, py::arg("cost_ceiling") = py::none(),
           py::arg("residual") = py::none(), py::arg("eligibility") = nullptr)
      .def_readonly("k", &PyKspOptions::k, "the number of paths as given (validated by ksp: an integer >= 1)")
      .def_property_readonly("tie_break", field_object(&PyKspOptions::tie_break), "the RouteTieBreak of paths of equal cost")
      .def_readonly("parallel_edges", &PyKspOptions::parallel_edges,
                    "True: each parallel edge a separate path; False: one per (u, v) pair by the tie rule")
      .def_readonly("cost_ceiling", &PyKspOptions::cost_ceiling,
                    "the inclusive cost ceiling as given: None, AbsoluteCost or CostStretch (validated by ksp)")
      .def_readonly("residual", &PyKspOptions::residual,
                    "the residual units per EdgeId as given (validated by ksp); None: no capacity filter;\n"
                    "priorities read static capacity")
      .def_property_readonly("eligibility", shared_field(&PyKspOptions::eligibility), "the static Eligibility; None: everything allowed");

  m.def(
      "ksp",
      [](const std::shared_ptr<Graph>& graph, const py::object& src, const py::object& dst, const py::object& opts) {
        if (!graph) throw py::value_error("graph must be a Graph");
        const auto o = options_of<PyKspOptions>(opts, "KspOptions");
        const auto n = graph->num_nodes();
        // ksp's validation order on unnarrowed values (Core repeats every check): the tie_break, the
        // source, the destination, source == destination, k, the ceiling, the residual, then a foreign
        // Eligibility.
        detail::check_member(o.tie_break, "tie_break");
        const NodeId s = endpoint(src, n, "source");
        const NodeId t = endpoint(dst, n, "destination");
        if (s == t) throw py::value_error(detail::overlap_error(s, "source", "destination"));
        KspOptions c(count_of(o.k, "k"));
        c.tie_break = o.tie_break;
        c.parallel_edges = o.parallel_edges;
        c.cost_ceiling = cost_ceiling_of(o.cost_ceiling);
        detail::check_ceiling(c.cost_ceiling);
        std::optional<std::vector<Amount>> res;
        if (!o.residual.is_none()) res = edge_units(o.residual, *graph, "residual");
        if (res.has_value()) c.residual = std::span<const Amount>(*res);
        c.eligibility = o.eligibility.get();
        std::vector<std::pair<Cost, Path>> paths;
        {
          py::gil_scoped_release release;
          paths = ksp(*graph, s, t, c);
        }
        // Every object through a checked allocation: a failure raises MemoryError.
        auto out = py::reinterpret_steal<py::list>(checked(PyList_New(static_cast<Py_ssize_t>(paths.size()))).release());
        for (std::size_t i = 0; i < paths.size(); ++i) {
          py::object cost = int_object(paths[i].first);
          py::object path = object_from<Path>(std::move(paths[i].second));
          PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), objects_tuple({cost, path}).release().ptr());
        }
        return out;
      },
      py::arg("graph"), py::arg("source"), py::arg("destination"), py::arg("options"),
      "The k shortest simple source -> destination paths over the eligible edges (only those with residual >= 1 when\n"
      "a residual is given; an edge of capacity 0 is searched like any other) -> [(cost, Path)], each hop\n"
      "a single edge. Order: lowest cost; then, under a HIGHER_* tie_break only, the highest bottleneck\n"
      "(the smallest tie value along the path); then fewest hops; then the smallest EdgeId\n"
      "sequence. The result is the first k paths of that order\n"
      "within the inclusive cost ceiling. Validated in this order before any work (ValueError): source and\n"
      "destination in range, source == destination, k >= 1, the ceiling, the residual (length, then every\n"
      "entry), a foreign Eligibility. An excluded endpoint, no static route under a CostStretch ceiling, or\n"
      "no path gives [].");

  m.def(
      "resolve_to_paths",
      [](const std::shared_ptr<Graph>& graph, const std::shared_ptr<PredecessorDAG>& dag, const py::object& src,
         const py::object& dst, StrictBool split_parallel_edges, const py::object& max_paths) {
        if (!graph) throw py::value_error("graph must be a Graph");
        if (!dag) throw py::value_error("dag must be a PredecessorDAG");
        if (dag->graph_uid() != graph->uid()) throw py::value_error("foreign PredecessorDAG: built for another graph");
        const auto n = graph->num_nodes();
        const NodeId s = endpoint(src, n, "source");
        const NodeId t = endpoint(dst, n, "destination");
        if (s == t) throw py::value_error(detail::overlap_error(s, "source", "destination"));
        std::optional<std::int64_t> mp;
        if (!max_paths.is_none()) mp = count_of(max_paths, "max_paths");
        std::vector<Path> paths;
        {
          py::gil_scoped_release release;
          paths = resolve_to_paths(*graph, *dag, s, t, split_parallel_edges, mp);
        }
        auto out = py::reinterpret_steal<py::list>(checked(PyList_New(static_cast<Py_ssize_t>(paths.size()))).release());
        for (std::size_t i = 0; i < paths.size(); ++i) {
          PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), object_from<Path>(std::move(paths[i])).release().ptr());
        }
        return out;
      },
      py::arg("graph"), py::arg("dag"), py::arg("source"), py::arg("destination"), py::arg("split_parallel_edges") = StrictBool {false},
      py::arg("max_paths") = py::none(),
      "The source -> destination walks of a DAG -> [Path], in increasing (cost, hops, edge-id sequence):\n"
      "a hop costs its cheapest member and is keyed by its lowest member id. split_parallel_edges = False\n"
      "keeps every parallel member of a (u, v) pair in one hop (members ascending); True gives one path\n"
      "per member combination, sorted by that path's own (cost, hops, edge-id sequence). max_paths\n"
      "(>= 1) truncates; None lists every walk.\n"
      "ValueError: a DAG of another graph, source or destination out of range, source == destination,\n"
      "max_paths < 1. An endpoint outside the DAG, or no walk, gives [].");
}

} // namespace netgraph::core::bindings
