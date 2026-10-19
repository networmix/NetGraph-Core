// ShortestPathCacheStats, ShortestPathCacheMiss and ShortestPathCache: the cache of static shortest-path DAGs.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>

#include "netgraph/core/shortest_path_cache.hpp"
#include "netgraph/core/graph.hpp"
#include "registration.hpp"
#include "values.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

void register_shortest_path_cache(py::module_& m) {
  py::class_<ShortestPathCache::Stats>(m, "ShortestPathCacheStats",
                              "ShortestPathCache statistics (atomic per field, not a consistent snapshot): hits, misses,\n"
                              "entries, bytes (the entries' retained payload), live_values (values alive anywhere).")
      .def_property_readonly("hits", int_reader(&ShortestPathCache::Stats::hits), "lookups answered from an entry (they persist across clear())")
      .def_property_readonly("misses", int_reader(&ShortestPathCache::Stats::misses),
                    "lookups without an entry: a building cache searched, a frozen one raised ShortestPathCacheMiss")
      .def_property_readonly("entries", int_reader(&ShortestPathCache::Stats::entries), "the entries held")
      .def_property_readonly("bytes", int_reader(&ShortestPathCache::Stats::bytes), "the retained payload of the entries (their DAGs' stored arrays)")
      .def_property_readonly("live_values", int_reader(&ShortestPathCache::Stats::live_values), "the cached values alive anywhere, entries included");

  py::register_exception<ShortestPathCacheMiss>(m, "ShortestPathCacheMiss", PyExc_Exception);
  m.attr("ShortestPathCacheMiss").attr("__doc__") =
      "A lookup in a frozen ShortestPathCache whose key is absent: a dedicated exception (not a\n"
      "RuntimeError or ValueError); the evaluators propagate it unchanged.";

  py::class_<ShortestPathCache, std::shared_ptr<ShortestPathCache>>(
      m, "ShortestPathCache",
      "A cache of static shortest-path DAGs (from costs and the Eligibility only) for one graph, keyed by\n"
      "(the root set, the direction, the Eligibility uid). It stores the full DAG and applies each\n"
      "caller's route choice afterwards, so the route choice is never cached. ShortestPaths policies and the\n"
      "evaluators read it.\n"
      "Building mode (a new cache): one writer, no concurrent readers; a miss searches and inserts;\n"
      "clear() empties it (the hit and miss counters persist), and values already handed out stay valid.\n"
      "Frozen mode (freeze()): read-only and shareable by every worker; a miss counts one miss and raises\n"
      "ShortestPathCacheMiss; clear() raises ValueError; thaw() returns to building once the caller has joined\n"
      "every reader.\n"
      "When a public call fails: a building cache passed to it keeps the complete entries and statistics\n"
      "of the lookups that succeeded before the failure (never a partial entry; a kept entry never\n"
      "changes a result); a frozen cache's statistics stay unchanged, except one miss when the call\n"
      "raised ShortestPathCacheMiss. A frozen cache counts its hits only when the call succeeds, Python result\n"
      "included.")
      .def(py::init([](const std::shared_ptr<Graph>& graph) {
             if (!graph) throw py::value_error("graph must be a Graph");
             return std::make_shared<ShortestPathCache>(*graph);
           }),
           py::arg("graph"), "An empty cache in building mode, bound to graph (by its uid).")
      .def("graph_uid", &ShortestPathCache::graph_uid, "uid of the graph this cache is bound to")
      .def("freeze", &ShortestPathCache::freeze, "enter frozen mode (no effect when frozen)")
      .def("thaw", &ShortestPathCache::thaw, "return to building mode (no effect when building); join every reader first")
      .def("clear", &ShortestPathCache::clear, "empty a building cache; ValueError when frozen")
      .def("stats", &ShortestPathCache::stats, "the current ShortestPathCacheStats (each field read atomically)");
}

} // namespace netgraph::core::bindings
