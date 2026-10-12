// DagCacheStats, DagCacheMiss and DagCache: the cache of static shortest-path DAGs.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/graph.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

void register_dag_cache(py::module_& m) {
  py::class_<DagCache::Stats>(m, "DagCacheStats",
                              "DagCache statistics (atomic per field, not a consistent snapshot): hits, misses,\n"
                              "entries, bytes (the entries' retained payload), live_values (values alive anywhere).")
      .def_readonly("hits", &DagCache::Stats::hits, "lookups answered from an entry (they persist across clear())")
      .def_readonly("misses", &DagCache::Stats::misses,
                    "lookups without an entry: a building cache searched, a frozen one raised DagCacheMiss")
      .def_readonly("entries", &DagCache::Stats::entries, "the entries held")
      .def_readonly("bytes", &DagCache::Stats::bytes, "the retained payload of the entries (their DAGs' stored arrays)")
      .def_readonly("live_values", &DagCache::Stats::live_values, "the cached values alive anywhere, entries included");

  py::register_exception<DagCacheMiss>(m, "DagCacheMiss", PyExc_Exception);
  m.attr("DagCacheMiss").attr("__doc__") =
      "A lookup in a frozen DagCache whose key is absent: a dedicated exception (not a\n"
      "RuntimeError or ValueError); the evaluators propagate it unchanged.";

  py::class_<DagCache, std::shared_ptr<DagCache>>(
      m, "DagCache",
      "A cache of static shortest-path DAGs (from costs and the Eligibility only) for one graph, keyed by\n"
      "(the root set, the direction, the Eligibility uid). It stores the full DAG and applies each\n"
      "caller's selection afterwards, so the selection is never cached. CostOnly policies and the\n"
      "evaluators read it.\n"
      "Building mode (a new cache): one writer, no concurrent readers; a miss searches and inserts;\n"
      "clear() empties it (the hit and miss counters persist), and values already handed out stay valid.\n"
      "Frozen mode (freeze()): read-only and shareable by every worker; a miss counts one miss and raises\n"
      "DagCacheMiss; clear() raises ValueError; thaw() returns to building once the caller has joined\n"
      "every reader.\n"
      "When a public call fails: a building cache passed to it keeps the complete entries and statistics\n"
      "of the lookups that succeeded before the failure (never a partial entry; a kept entry never\n"
      "changes a result); a frozen cache's statistics stay unchanged, except one miss when the call\n"
      "raised DagCacheMiss. A frozen cache counts its hits only when the call succeeds, Python result\n"
      "included.")
      .def(py::init([](const std::shared_ptr<Graph>& graph) {
             if (!graph) throw py::value_error("graph must be a Graph");
             return std::make_shared<DagCache>(*graph);
           }),
           py::arg("graph"), "An empty cache in building mode, bound to graph (by its uid).")
      .def("graph_uid", &DagCache::graph_uid, "uid of the graph this cache is bound to")
      .def("freeze", &DagCache::freeze, "enter frozen mode (no effect when frozen)")
      .def("thaw", &DagCache::thaw, "return to building mode (no effect when building); join every reader first")
      .def("clear", &DagCache::clear, "empty a building cache; ValueError when frozen")
      .def("stats", &DagCache::stats, "the current DagCacheStats (each field read atomically)");
}

} // namespace netgraph::core::bindings
