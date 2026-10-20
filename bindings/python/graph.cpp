// Graph: the immutable graph on integer units and its validating constructor.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <limits>
#include <memory>
#include <string>

#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/graph.hpp"
#include "integers.hpp"
#include "strict_float.hpp"
#include "values.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// The unit-value step on the capacity argument whatever its shape (the shape step comes next):
// every integer entry it holds must lie in [0, 2^62). A NumPy array is read flattened; nested
// entries of a Python sequence and a non-integer non-sequence are left to the shape step.
void check_unit_values(const py::handle& edge_capacity_units) {
  const auto out_of_domain = [](const std::string& what, const std::string& value) {
    return py::value_error(what + " = " + value + " is outside the units domain [0, 2^62)");
  };
  if (py::isinstance<py::array>(edge_capacity_units)) {
    py::object flat = edge_capacity_units.attr("ravel")();
    const auto w = wide_integers(flat, "edge_capacity_units");
    for (std::size_t i = 0; i < w.size(); ++i) {
      if (!w.in(i, 0, kDomainLimit - 1)) throw out_of_domain("edge_capacity_units[" + std::to_string(i) + "]", w.text(i));
    }
    return;
  }
  if (!is_invalid_sequence(edge_capacity_units)) {
    auto seq = py::reinterpret_borrow<py::sequence>(edge_capacity_units);
    for (std::size_t i = 0; i < seq.size(); ++i) {
      py::object item = seq[i];
      if (is_nested(item)) continue;
      const std::string what = "edge_capacity_units[" + std::to_string(i) + "]";
      py::object idx = python_int(item, what);
      int overflow = 0;
      const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
      if (overflow != 0 || v < 0 || v >= kDomainLimit) throw out_of_domain(what, std::string(py::str(idx)));
    }
    return;
  }
  if (!is_bool(edge_capacity_units) && PyIndex_Check(edge_capacity_units.ptr())) {
    py::object idx = python_int(edge_capacity_units, "edge_capacity_units");
    int overflow = 0;
    const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
    if (overflow != 0 || v < 0 || v >= kDomainLimit) throw out_of_domain("edge_capacity_units", std::string(py::str(idx)));
  }
}

std::shared_ptr<Graph> make_graph(const py::object& num_nodes, const py::object& source,
                                  const py::object& destination, const py::object& edge_capacity_units,
                                  const py::object& cost, StrictFloat resolution) {
  // The constructor's validation order: the resolution, the unit values, the shapes, zero or
  // negative costs, the remaining bounds (ids in range, then the sums); every value is read
  // without narrowing, and the messages are the C++ constructor's, which re-checks the same order.
  (void)resolution_exponent(resolution);
  check_unit_values(edge_capacity_units);
  // The shapes: num_nodes, then every array argument rank 1 (in signature order), then equal
  // lengths; no entry of source, destination or cost is read before the cost step.
  const auto n = integer_value(num_nodes, 0, std::numeric_limits<std::int32_t>::max(), "num_nodes");
  const std::size_t m = sequence_length(source, "edge_sources");
  const std::size_t md = sequence_length(destination, "edge_destinations");
  const std::size_t mu = sequence_length(edge_capacity_units, "edge_capacity_units");
  const std::size_t mc = sequence_length(cost, "edge_costs");
  if (md != m || mu != m || mc != m) {
    throw py::value_error("edge_sources, edge_destinations, edge_capacity_units and edge_costs must have the same length (got " +
                          std::to_string(m) + ", " + std::to_string(md) + ", " + std::to_string(mu) + ", " +
                          std::to_string(mc) + ")");
  }
  const auto cap = wide_integers(edge_capacity_units, "edge_capacity_units");  // rank 1, values checked
  if (m > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    throw py::value_error("the number of edges exceeds INT32_MAX (array length)");
  }
  const auto w = wide_integers(cost, "edge_costs");
  for (std::size_t i = 0; i < m; ++i) {
    if (w.sign(i) < 0 || (w.sign(i) == 0 && w.value[i] < 1)) {
      throw py::value_error("edge_costs[" + std::to_string(i) + "] = " + w.text(i) + ": costs must be positive (edge_costs >= 1)");
    }
  }
  const auto s = wide_integers(source, "edge_sources");
  const auto d = wide_integers(destination, "edge_destinations");
  for (std::size_t i = 0; i < m; ++i) {
    if (!s.in(i, 0, n - 1) || !d.in(i, 0, n - 1)) {
      throw py::value_error("edge " + std::to_string(i) + ": endpoint (" + s.text(i) + ", " + d.text(i) +
                            ") out of node range [0, " + std::to_string(n) + ")");
    }
  }
  constexpr auto kLimit = static_cast<std::uint64_t>(kDomainLimit);
  std::uint64_t cap_sum = 0;
  std::uint64_t cost_sum = 0;
  for (std::size_t i = 0; i < m; ++i) {
    cap_sum += static_cast<std::uint64_t>(cap.value[i]);  // each below 2^62: no wrap before the check
    if (cap_sum >= kLimit) {
      throw py::value_error("the sum of all capacities reaches the bound of the units domain (2^62 units)");
    }
    if (w.sign(i) > 0 || static_cast<std::uint64_t>(w.value[i]) >= kLimit - cost_sum) {
      throw py::value_error("the sum of all edge costs reaches the bound 2^62: path costs would overflow int64");
    }
    cost_sum += static_cast<std::uint64_t>(w.value[i]);
  }
  auto cap_v = cap.narrow<Amount>();
  auto cost_v = w.narrow<Cost>();
  auto src_v = s.narrow<NodeId>();
  auto dst_v = d.narrow<NodeId>();
  py::gil_scoped_release release;
  return std::make_shared<Graph>(static_cast<std::int32_t>(n), src_v, dst_v, cap_v, cost_v, resolution);
}

} // namespace

void register_graph(py::module_& m) {
  py::class_<Graph, std::shared_ptr<Graph>>(
      m, "Graph",
      "Immutable directed multigraph on integer units (parallel edges and self-loops allowed). Edges are\n"
      "given in input order and numbered by sorting on (source, cost, destination), ties kept in input\n"
      "order: that number is the public EdgeId. input_positions()[e] is the input position of edge e.\n"
      "Per-edge values use EdgeId order; incoming_edges() instead lists EdgeIds in incoming-slot order.\n"
      "Costs are positive integers; capacities are int64 units of the explicit resolution q (a power of\n"
      "two in [2^-32, 2^32]); physical capacities enter through to_units.")
      .def(py::init(&make_graph), py::arg("num_nodes"), py::arg("edge_sources"), py::arg("edge_destinations"),
           py::arg("edge_capacity_units"), py::arg("edge_costs"), py::arg("resolution"),
           "Validates, in order, and reports the first defect: the resolution (a power of two in\n"
           "[2^-32, 2^32]; a bool is a TypeError), the unit values (every capacity in [0, 2^62)), the array\n"
           "shapes, zero or negative costs, the remaining bounds (endpoints in range, then the capacity\n"
           "and cost sums below 2^62); a value of a wrong type (a bool or a float for an integer, a\n"
           "non-sequence, an array of another rank or dtype) raises TypeError, a value outside its domain\n"
           "ValueError. Integers are read without narrowing, so the error does not depend on how the\n"
           "input is represented.")
      .def("uid", &Graph::uid, "process-unique id (never reused)")
      .def("resolution", &Graph::resolution, "the resolution q (units of capacity)")
      .def("num_nodes", &Graph::num_nodes, "N: node ids are 0 .. N - 1")
      .def("num_edges", &Graph::num_edges, "E: EdgeIds are 0 .. E - 1")
      .def("max_edge_cost", &Graph::max_edge_cost, "largest edge cost (0 without edges)")
      .def("edge_cost_gcd", &Graph::edge_cost_gcd, "gcd of the edge costs (0 without edges)")
      .def(
          "edge_capacity_units", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.edge_capacity_units(), self);
          },
          "int64[E] capacities in units, id order (read-only view)")
      .def(
          "edge_costs", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.edge_costs(), self);
          },
          "int64[E] costs, id order (read-only view)")
      .def(
          "edge_sources", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.edge_sources(), self);
          },
          "int32[E] source node per edge (read-only view)")
      .def(
          "edge_destinations", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.edge_destinations(), self);
          },
          "int32[E] destination node per edge (read-only view)")
      .def(
          "input_positions", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.input_positions(), self);
          },
          "int32[E]: EdgeId -> input position (read-only view); build the inverse yourself if needed")
      .def(
          "outgoing_offsets", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.outgoing_offsets(), self);
          },
          "int32[N + 1]: the out-edges of u are the ids [outgoing_offsets()[u], outgoing_offsets()[u + 1])")
      .def(
          "incoming_offsets", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.incoming_offsets(), self);
          },
          "int32[N + 1]: the in-edges of v are incoming_edges()[incoming_offsets()[v] : incoming_offsets()[v + 1]]")
      .def(
          "incoming_edges", [](const py::object& self) {
            const auto& g = self.cast<const Graph&>();
            return readonly_view(g.incoming_edges(), self);
          },
          "int32[E]: the in-edges of v in (cost, source, input) order, as EdgeIds")
      .def(
          "eligibility",
          [](const Graph& g, const py::object& excluded_nodes, const py::object& excluded_edges) {
            // Nodes, then edges, each range-checked before narrowing (the constructor's order).
            const auto xn = wide_integers(excluded_nodes, "excluded_nodes");
            for (std::size_t i = 0; i < xn.size(); ++i) {
              if (!xn.in(i, 0, static_cast<std::int64_t>(g.num_nodes()) - 1)) {
                throw py::value_error("excluded node " + xn.text(i) + " out of node range [0, " +
                                      std::to_string(g.num_nodes()) + ")");
              }
            }
            const auto xe = wide_integers(excluded_edges, "excluded_edges");
            for (std::size_t i = 0; i < xe.size(); ++i) {
              if (!xe.in(i, 0, static_cast<std::int64_t>(g.num_edges()) - 1)) {
                throw py::value_error("excluded edge " + xe.text(i) + " out of edge range [0, " +
                                      std::to_string(g.num_edges()) + ")");
              }
            }
            return std::make_shared<Eligibility>(g.eligibility(xn.narrow<NodeId>(), xe.narrow<EdgeId>()));
          },
          py::arg("excluded_nodes") = py::tuple(), py::arg("excluded_edges") = py::tuple(),
          "Static eligibility from excluded node ids and EdgeIds (duplicates and unsorted lists\n"
          "accepted, out-of-range ids rejected); an edge is allowed when it is not excluded and\n"
          "both endpoints are allowed. Empty lists mean all allowed (uid 0).");
}

} // namespace netgraph::core::bindings
