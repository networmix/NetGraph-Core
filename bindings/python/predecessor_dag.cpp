// PredecessorDAG: the immutable predecessor DAG of one graph.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>
#include <vector>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/route_structures.hpp"
#include "integers.hpp"
#include "values.hpp"
#include "native_call.hpp"
#include "netgraph/core/detail/route_builders.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// The exact order of wide integers: a positive overflow above every int64, a negative one
// below, and two overflows on the same side by their exact values.
bool wide_less(const CheckedIntegerVector& w, std::size_t i, std::size_t j) {
  if (w.sign(i) != w.sign(j)) return w.sign(i) < w.sign(j);
  if (w.sign(i) == 0) return w.value[i] < w.value[j];
  const int less = PyObject_RichCompareBool(w.exact(i).ptr(), w.exact(j).ptr(), Py_LT);
  if (less < 0) throw py::error_already_set();
  return less == 1;
}

} // namespace

void register_predecessor_dag(py::module_& m) {
  py::class_<PredecessorDAG, std::shared_ptr<PredecessorDAG>>(
      m, "PredecessorDAG",
      "Immutable route structure of one graph, a predecessor DAG (directed, no cycle): for each node v,\n"
      "the edges that reach v on the route and the node each comes from (its parent). The entries of v\n"
      "are parents()[k], via_edges()[k] for k in [parent_offsets()[v], parent_offsets()[v + 1]), in\n"
      "strictly ascending EdgeId, so two DAGs of the same graph with the same entries have equal arrays\n"
      "and compare equal; DAGs of different graphs never compare equal.\n"
      "parents() and via_edges() are read-only views that keep the DAG alive; parent_offsets() is a new\n"
      "read-only array on each call (the DAG stores only the nodes that have entries). Built by\n"
      "PredecessorDAG.from_csr, Path.to_dag(), the searches and the policies (FlowRecord.route).")
      .def_static(
          "from_csr",
          [](const std::shared_ptr<Graph>& graph, const py::object& offsets, const py::object& parents,
             const py::object& via) {
            if (!graph) throw py::type_error("graph must be a Graph");
            // Validate the unnarrowed shape and ranges, then the native structural
            // rules below (endpoints, duplicates, acyclic, ascending), all before
            // detaching. The builder consumes these checked native arrays.
            const auto n = static_cast<std::size_t>(graph->num_nodes());
            const auto edge_count = static_cast<std::int64_t>(graph->num_edges());
            const std::size_t lo = sequence_length(offsets, "parent_offsets");
            const std::size_t lp = sequence_length(parents, "parents");
            const std::size_t lv = sequence_length(via, "via_edges");
            if (lo != n + 1) {
              throw py::value_error("PredecessorDAG.from_csr: parent_offsets has length " + std::to_string(lo) +
                                    ", expected N + 1 = " + std::to_string(n + 1));
            }
            const auto o = wide_integers(offsets, "parent_offsets");
            if (!o.in(0, 0, 0)) throw py::value_error("PredecessorDAG.from_csr: parent_offsets[0] = " + o.text(0) + ", expected 0");
            for (std::size_t v = 0; v < n; ++v) {
              if (wide_less(o, v + 1, v)) {
                throw py::value_error("PredecessorDAG.from_csr: parent_offsets are not monotone at " + std::to_string(v + 1) + " (" +
                                      o.text(v) + " > " + o.text(v + 1) + ")");
              }
            }
            if (!o.in(n, static_cast<std::int64_t>(lp), static_cast<std::int64_t>(lp)) || lp != lv) {
              throw py::value_error("PredecessorDAG.from_csr: parent_offsets[N] = " + o.text(n) + " but parents has " +
                                    std::to_string(lp) + " and via_edges has " + std::to_string(lv) + " entries");
            }
            const auto p = wide_integers(parents, "parents");
            const auto w = wide_integers(via, "via_edges");
            for (std::size_t k = 0; k < lp; ++k) {
              if (!p.in(k, 0, static_cast<std::int64_t>(n) - 1)) {
                throw py::value_error("PredecessorDAG.from_csr: parents[" + std::to_string(k) + "] = " + p.text(k) +
                                      " out of node range [0, " + std::to_string(n) + ")");
              }
              if (!w.in(k, 0, edge_count - 1)) {
                throw py::value_error("PredecessorDAG.from_csr: via_edges[" + std::to_string(k) + "] = " + w.text(k) +
                                      " out of edge range [0, " + std::to_string(edge_count) + ")");
              }
            }
            const auto ov = o.narrow<std::int32_t>();
            const auto pv = p.narrow<NodeId>();
            const auto wv = w.narrow<EdgeId>();
            detail::check_csr(*graph, ov, pv, wv);
            return native_call([&] { return std::make_shared<PredecessorDAG>(detail::PredecessorDAGBuilder::of_csr(*graph, ov, pv, wv)); });
          },
          py::arg("graph"), py::arg("parent_offsets"), py::arg("parents"), py::arg("via_edges"),
          "A wrong type raises TypeError; structural and range defects raise ValueError. Validates, in order,\n"
          "before any traversal: parent_offsets of\n"
          "length N + 1, parent_offsets[0] == 0, monotone, parent_offsets[N] == len(parents) ==\n"
          "len(via_edges); then every parent and via edge in range, every via edge joining (parent, child),\n"
          "no duplicate (child, via edge), acyclic; within a node the entries in ascending via EdgeId.")
      .def("graph_uid", &PredecessorDAG::graph_uid, "uid of the graph this DAG belongs to")
      .def(
          "parent_offsets",
          [](const PredecessorDAG& d) {
            const auto off = d.parent_offsets();
            py::array_t<std::int32_t> arr(static_cast<py::ssize_t>(off.size()));
            std::copy(off.begin(), off.end(), arr.mutable_data());
            py::detail::array_proxy(arr.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
            return arr;
          },
          "int32[N + 1]: the entries of node v are [parent_offsets[v], parent_offsets[v + 1]); a new read-only array per call")
      .def(
          "parents", [](const py::object& self) { return readonly_view(self.cast<const PredecessorDAG&>().parents(), self); },
          "int32 parent NodeId per entry (read-only view)")
      .def(
          "via_edges",
          [](const py::object& self) { return readonly_view(self.cast<const PredecessorDAG&>().via_edges(), self); },
          "int32 via EdgeId per entry (read-only view)")
      .def("__eq__", [](const PredecessorDAG& a, const py::object& b) {
        return py::isinstance<PredecessorDAG>(b) && a == b.cast<const PredecessorDAG&>();
      })
      .def("__hash__", [](const PredecessorDAG& a) {
        return hash_int(objects_tuple({int_object(a.graph_uid()), ints_tuple(a.via_edges())}));
      })
      .def("__repr__", [](const PredecessorDAG& a) {
        return "PredecessorDAG(graph_uid=" + std::to_string(a.graph_uid()) + ", entries=" +
               std::to_string(a.via_edges().size()) + ")";
      });
}

} // namespace netgraph::core::bindings
