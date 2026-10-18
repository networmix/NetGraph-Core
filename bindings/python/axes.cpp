// The enumerations of the placement, origin and selection axes, and PathSelection.
#include <pybind11/pybind11.h>

#include <string>

#include "netgraph/core/options.hpp"

#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

void register_axes(py::module_& m) {
  py::enum_<FlowPlacement>(m, "FlowPlacement",
                           "How traffic spreads over the routes (CostOnly's kernel; CapacityAware and Pinned flows too).\n"
                           "MAX_FLOW: the minimum-cost flow core over the routes (a controller restricted to them), re-solved\n"
                           "over the free capacity plus the policy's own allocation (it may move its own earlier placement,\n"
                           "never another's). WEIGHTED (WCMP, weighted-cost multipath): every node splits in proportion to\n"
                           "weights fixed when the routes are prepared (highest averages: whole-unit shares as close to\n"
                           "proportional as possible). EQUAL_BALANCED (ECMP, equal-cost multipath): every node splits\n"
                           "equally, by round robin.")
      .value("MAX_FLOW", FlowPlacement::MaxFlow)
      .value("WEIGHTED", FlowPlacement::Weighted)
      .value("EQUAL_BALANCED", FlowPlacement::EqualBalanced);

  py::enum_<Admission>(m, "Admission",
                       "LOSSLESS admits the largest total the routes carry with no loss at all; LOSSY sends the whole\n"
                       "finite offer hop by hop and records where it is lost (EQUAL_BALANCED and WEIGHTED only).")
      .value("LOSSLESS", Admission::Lossless)
      .value("LOSSY", Admission::Lossy);

  py::enum_<Origin>(m, "Origin",
                    "How a demand's traffic enters at its sources, separately from how it is routed. EQUAL: round\n"
                    "robin over the roster (ascending NodeId, the first sources one unit more); a source whose share\n"
                    "cannot go anywhere keeps it and loses it (LOSSLESS then admits only totals that give it\n"
                    "nothing). OPTIMIZED: the controller chooses the injections (one virtual source).")
      .value("EQUAL", Origin::Equal)
      .value("OPTIMIZED", Origin::Optimized);

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
                            "The selection of a search or a policy: tie_break (DETERMINISTIC by default), priority\n"
                            "(CAPACITY by default); multipath = True keeps every equal-cost walk (the DAG), False the one\n"
                            "canonical single path; multi_edge = True keeps every eligible parallel member, False one per\n"
                            "(u, v) pair. A tie_break or priority outside its enumeration raises ValueError.")
      // The constructor checks both enumerations (ValueError naming the field).
      .def(py::init<EdgeTieBreak, EdgePriority, bool, bool>(), py::arg_v("tie_break", EdgeTieBreak::Deterministic, "EdgeTieBreak.DETERMINISTIC"),
           py::arg_v("priority", EdgePriority::Capacity, "EdgePriority.CAPACITY"), py::arg("multipath") = true, py::arg("multi_edge") = true)
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
