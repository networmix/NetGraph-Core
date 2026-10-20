// Path: an explicit walk of hops of one graph.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/route_structures.hpp"
#include "integers.hpp"
#include "registration.hpp"
#include "values.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

void register_path(py::module_& m) {
  py::class_<Path, std::shared_ptr<Path>>(
      m, "Path",
      "Immutable route of one graph, hop by hop: each hop lists one or more parallel edges (members)\n"
      "from the same u to the same v, in the caller's order. cost = the sum over hops of the cost of\n"
      "each hop's cheapest member.")
      .def(py::init([](const std::shared_ptr<Graph>& graph, const py::object& hops) {
             if (!graph) throw py::type_error("graph must be a Graph");
             if (is_invalid_sequence(hops)) throw py::type_error("Path: hops must be a sequence of hops (sequences of EdgeIds)");
             auto seq = py::reinterpret_borrow<py::sequence>(hops);
             std::vector<CheckedIntegerVector> raw;
             for (std::size_t h = 0; h < seq.size(); ++h) {
               py::object hop = seq[h];
               raw.push_back(wide_integers(hop, ("hops[" + std::to_string(h) + "]").c_str()));
             }
             // Every member in range before anything else (no id is narrowed unchecked), then the
             // C++ constructor's order.
             const auto edge_count = static_cast<std::int64_t>(graph->num_edges());
             for (std::size_t h = 0; h < raw.size(); ++h) {
               for (std::size_t i = 0; i < raw[h].size(); ++i) {
                 if (!raw[h].in(i, 0, edge_count - 1)) {
                   throw py::value_error("Path: hop " + std::to_string(h) + " member " + raw[h].text(i) +
                                         " out of edge range [0, " + std::to_string(edge_count) + ")");
                 }
               }
             }
             std::vector<std::vector<EdgeId>> narrowed;
             narrowed.reserve(raw.size());
             for (const auto& r : raw) narrowed.push_back(r.narrow<EdgeId>());
             return std::make_shared<Path>(*graph, std::move(narrowed));
           }),
           py::arg("graph"), py::arg("hops"),
           "Validates, in order (ValueError on the first defect): every member id in range, no duplicate\n"
           "member, a non-empty hop list, every hop non-empty, the members of a hop sharing (u, v),\n"
           "contiguous hops, no repeated node.")
      .def("graph_uid", &Path::graph_uid, "uid of the graph this path belongs to")
      .def("hops", [](const Path& p) { return p.hops(); }, "the hops as lists of EdgeIds (the caller's member order)")
      .def("source", &Path::source, "the first node of the walk")
      .def("destination", &Path::destination, "the last node of the walk")
      .def("cost", &Path::cost, "sum over hops of the minimum member cost")
      .def(
          "to_dag", [](const Path& p) { return std::make_shared<PredecessorDAG>(p.to_dag()); },
          "the equivalent PredecessorDAG (every member of every hop)")
      .def("__eq__", [](const Path& a, const py::object& b) {
        return py::isinstance<Path>(b) && a == b.cast<const Path&>();
      })
      .def("__hash__",
           [](const Path& a) {
             const auto& path_hops = a.hops();
             auto hops = py::reinterpret_steal<py::tuple>(checked(PyTuple_New(static_cast<Py_ssize_t>(path_hops.size()))).release());
             for (std::size_t i = 0; i < path_hops.size(); ++i) {
               PyTuple_SET_ITEM(hops.ptr(), static_cast<Py_ssize_t>(i), ints_tuple(path_hops[i]).release().ptr());
             }
             return hash_int(objects_tuple({int_object(a.graph_uid()), hops}));
           })
      .def("__repr__", [](const Path& a) {
        return "Path(source=" + std::to_string(a.source()) + ", destination=" + std::to_string(a.destination()) +
               ", hops=" + std::to_string(a.hops().size()) + ", cost=" + std::to_string(a.cost()) + ")";
      });
}

} // namespace netgraph::core::bindings
