// The enumerations of the placement and selection axes, and PathSelection.
#include <pybind11/pybind11.h>

#include <string>

#include "netgraph/core/options.hpp"

#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

void register_axes(py::module_& m) {
  py::enum_<FlowPlacement>(m, "FlowPlacement",
                           "How one flow spreads over its route (its kernel). MAX_FLOW: the minimum-cost flow core, re-solved\n"
                           "on the flow's route over the free capacity plus the flow's own allocation (it may move its own\n"
                           "earlier placement, never another flow's). WEIGHTED (WCMP, weighted-cost multipath): every node\n"
                           "splits in proportion to weights fixed when the flow is created (highest averages: whole-unit\n"
                           "shares as close to proportional as possible). EQUAL_BALANCED (ECMP, equal-cost multipath): every\n"
                           "node splits equally, by round robin.")
      .value("MAX_FLOW", FlowPlacement::MaxFlow)
      .value("WEIGHTED", FlowPlacement::Weighted)
      .value("EQUAL_BALANCED", FlowPlacement::EqualBalanced);

  py::enum_<Admission>(m, "Admission",
                       "LOSSLESS admits only what fits; LOSSY forwards hop by hop without admission control and\n"
                       "drops what does not fit on each edge (EQUAL_BALANCED and WEIGHTED only).")
      .value("LOSSLESS", Admission::Lossless)
      .value("LOSSY", Admission::Lossy);

  py::enum_<EdgeTieBreak>(m, "EdgeTieBreak",
                          "How a search breaks cost ties. Between parallel members of one (u, v) pair\n"
                          "(multi_edge = False): DETERMINISTIC keeps the lowest EdgeId, PREFER_HIGHER_PRIORITY the\n"
                          "higher priority value, then the lowest EdgeId. Between equal-cost walks (multipath =\n"
                          "False): PREFER_HIGHER_PRIORITY first takes the highest bottleneck (the smallest priority\n"
                          "value along the walk), DETERMINISTIC skips that step; then the fewest hops, then the\n"
                          "lexicographically smallest EdgeId sequence.")
      .value("DETERMINISTIC", EdgeTieBreak::Deterministic)
      .value("PREFER_HIGHER_PRIORITY", EdgeTieBreak::PreferHigherPriority);

  py::enum_<EdgePriority>(m, "EdgePriority",
                          "The per-edge value PREFER_HIGHER_PRIORITY reads and, with WEIGHTED, the value that decides which\n"
                          "edges a flow uses and their weights: RESIDUAL (residual units), CAPACITY (static capacity units),\n"
                          "HEADROOM (residual / capacity, compared exactly as a fraction, 0 on an edge of capacity 0;\n"
                          "rejected with WEIGHTED). A search without a residual reads the static capacity as the residual:\n"
                          "RESIDUAL and CAPACITY read the capacity, HEADROOM 1 (0 on an edge of capacity 0).")
      .value("RESIDUAL", EdgePriority::Residual)
      .value("CAPACITY", EdgePriority::Capacity)
      .value("HEADROOM", EdgePriority::Headroom);

  py::enum_<Direction>(m, "Direction",
                       "The orientation of a search: FORWARD from the roots along the edges, or REVERSE toward the roots\n"
                       "against the edges (the DAG it returns keeps the edges' own direction, so its walks lead to the\n"
                       "roots).")
      .value("FORWARD", Direction::Forward)
      .value("REVERSE", Direction::Reverse);

  py::class_<PathSelection>(m, "PathSelection",
                            "The selection of a search: tie_break and priority are explicit; multipath = True keeps\n"
                            "every equal-cost walk (the DAG), False the one canonical single path; multi_edge = True\n"
                            "keeps every eligible parallel member, False one per (u, v) pair. A tie_break or priority\n"
                            "outside its enumeration raises ValueError.")
      // The constructor checks both enumerations (ValueError naming the field).
      .def(py::init<EdgeTieBreak, EdgePriority, bool, bool>(), py::arg("tie_break"), py::arg("priority"),
           py::arg("multipath") = true, py::arg("multi_edge") = true)
      .def_readonly("tie_break", &PathSelection::tie_break, "the EdgeTieBreak")
      .def_readonly("priority", &PathSelection::priority, "the EdgePriority")
      .def_readonly("multipath", &PathSelection::multipath,
                    "True: every equal-cost walk (a DAG); False: the canonical single path")
      .def_readonly("multi_edge", &PathSelection::multi_edge,
                    "True: every eligible parallel member; False: one member per (u, v) pair")
      .def("__eq__", [](const PathSelection& a, const py::object& b) {
        return py::isinstance<PathSelection>(b) && a == b.cast<const PathSelection&>();
      })
      .def("__hash__", [](const PathSelection& a) {
        return py::hash(py::make_tuple(static_cast<int>(a.tie_break), static_cast<int>(a.priority), a.multipath,
                                       a.multi_edge));
      })
      .def("__repr__", [](const PathSelection& a) {
        return std::string("PathSelection(tie_break=") +
               (a.tie_break == EdgeTieBreak::Deterministic ? "DETERMINISTIC" : "PREFER_HIGHER_PRIORITY") +
               ", priority=" +
               (a.priority == EdgePriority::Residual ? "RESIDUAL"
                                                     : (a.priority == EdgePriority::Capacity ? "CAPACITY" : "HEADROOM")) +
               ", multipath=" + (a.multipath ? "True" : "False") + ", multi_edge=" + (a.multi_edge ? "True" : "False") +
               ")";
      });
}

} // namespace netgraph::core::bindings
