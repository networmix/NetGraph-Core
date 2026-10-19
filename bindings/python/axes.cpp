// The enumerations of the placement, origin and selection axes, and PathSelection.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

#include "netgraph/core/options.hpp"

#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

void register_axes(py::module_& m) {
  py::enum_<FlowPlacement>(m, "FlowPlacement",
                           "How traffic spreads over the routes (CostOnly's kernel; the interior split of supplied PredDAG\n"
                           "routes, EQUAL_BALANCED or WEIGHTED).\n"
                           "MAX_FLOW: the minimum-cost flow core over the routes (a controller restricted to them), re-solved\n"
                           "over the free capacity plus the policy's own allocation (it may move its own earlier placement,\n"
                           "never another's). WEIGHTED (WCMP, weighted-cost multipath): every node splits in proportion to\n"
                           "weights fixed when the routes are prepared (highest averages: whole-unit shares as close to\n"
                           "proportional as possible; the weights come from WeightBasis). EQUAL_BALANCED (ECMP, equal-cost\n"
                           "multipath): every node splits equally, by round robin.")
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

  py::enum_<WeightBasis>(m, "WeightBasis",
                         "Where WEIGHTED member weights come from (CostOnly WEIGHTED and supplied PredDAG routes).\n"
                         "CAPACITY (the default): each member's static capacity. EXPLICIT: edge_weights, one integer per\n"
                         "edge (a positive weight on an edge of capacity 0 becomes 0). RESIDUAL_SNAPSHOT: each member's free\n"
                         "capacity read once when the routes are prepared, against one common view, and frozen for the\n"
                         "epoch. A member of weight 0 is pruned with the branches it leaves without a destination, inside\n"
                         "the route only. The selection's priority never chooses weights.")
      .value("CAPACITY", WeightBasis::Capacity)
      .value("EXPLICIT", WeightBasis::Explicit)
      .value("RESIDUAL_SNAPSHOT", WeightBasis::ResidualSnapshot);

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
                          "The per-edge value PREFER_HIGHER_PRIORITY reads (single-path and parallel-member choices):\n"
                          "RESIDUAL (residual units), CAPACITY (static capacity units), HEADROOM (residual / capacity,\n"
                          "compared exactly as a fraction, 0 on an edge of capacity 0). A search without a residual reads\n"
                          "the static capacity as the residual: RESIDUAL and CAPACITY read the capacity, HEADROOM 1 (0 on\n"
                          "an edge of capacity 0). T1's CSPF reads its unreserved view U instead (RESIDUAL U, HEADROOM\n"
                          "U / capacity).")
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
                            "(CAPACITY by default); multipath None (the default) is the model's natural shape (every\n"
                            "equal-cost walk for searches and CostOnly, one path for a tunnel model), True every equal-cost\n"
                            "walk (the DAG; tunnel models reject it), False the one canonical single path; multi_edge = True\n"
                            "keeps every eligible parallel member, False one per (u, v) pair. A tie_break or priority\n"
                            "outside its enumeration raises ValueError.")
      // The constructor checks both enumerations (ValueError naming the field).
      .def(py::init<EdgeTieBreak, EdgePriority, std::optional<bool>, bool>(), py::arg_v("tie_break", EdgeTieBreak::Deterministic, "EdgeTieBreak.DETERMINISTIC"),
           py::arg_v("priority", EdgePriority::Capacity, "EdgePriority.CAPACITY"), py::arg("multipath") = py::none(), py::arg("multi_edge") = true)
      .def_readonly("tie_break", &PathSelection::tie_break, "the EdgeTieBreak")
      .def_readonly("priority", &PathSelection::priority, "the EdgePriority")
      .def_readonly("multipath", &PathSelection::multipath,
                    "None: the model's natural shape; True: every equal-cost walk (a DAG); False: the canonical\n"
                    "single path")
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
               ", multipath=" + (!a.multipath.has_value() ? "None" : *a.multipath ? "True" : "False") +
               ", multi_edge=" + (a.multi_edge ? "True" : "False") +
               ")";
      });
}

} // namespace netgraph::core::bindings
