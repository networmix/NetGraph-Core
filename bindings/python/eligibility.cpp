// Eligibility: the static node and edge flags of one graph.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>
#include <utility>

#include "netgraph/core/eligibility.hpp"
#include "values.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

void register_eligibility(py::module_& m) {
  py::class_<Eligibility, std::shared_ptr<Eligibility>>(
      m, "Eligibility",
      "Which nodes and edges of one graph may be used: one flag per node and per edge, fixed when built\n"
      "(it does not depend on traffic). Built by Graph.eligibility(excluded_nodes, excluded_edges) and\n"
      "Eligibility.intersect.")
      .def("graph_uid", &Eligibility::graph_uid, "uid of the graph this value belongs to")
      .def("uid", &Eligibility::uid, "0 for the all-allowed value, a fresh uid otherwise")
      .def(
          "node_allowed",
          [](const py::object& self) {
            const auto& el = self.cast<const Eligibility&>();
            return readonly_view(el.node_allowed(), self, py::dtype("bool"));
          },
          "bool[N] read-only view")
      .def(
          "edge_allowed",
          [](const py::object& self) {
            const auto& el = self.cast<const Eligibility&>();
            return readonly_view(el.edge_allowed(), self, py::dtype("bool"));
          },
          "bool[E] read-only view (not excluded and both endpoints allowed)")
      .def_static(
          "intersect",
          [](const std::shared_ptr<Eligibility>& a, const std::shared_ptr<Eligibility>& b) -> std::shared_ptr<Eligibility> {
            if (!a || !b) throw py::value_error("Eligibility.intersect takes two Eligibility values");
            auto r = Eligibility::intersect(*a, *b);
            if (r.uid() == a->uid()) return a;
            if (r.uid() == b->uid()) return b;
            return std::make_shared<Eligibility>(std::move(r));
          },
          py::arg("a"), py::arg("b"),
          "The intersection of two values of one graph (ValueError for values of different graphs):\n"
          "an operand itself when the result equals it (a first, then b), otherwise a new value with a\n"
          "fresh uid; equal intersections built separately get distinct uids.");
}

} // namespace netgraph::core::bindings
