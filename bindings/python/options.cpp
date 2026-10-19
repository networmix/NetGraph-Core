// The enumerations: RouteAllocation, Admission, SourceAllocation, WeightBasis, RouteTieBreak and Direction.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

#include "netgraph/core/options.hpp"

#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

void register_options(py::module_& m) {
  py::enum_<RouteAllocation>(m, "RouteAllocation",
                           "How traffic spreads over the routes (ShortestPaths's kernel; the interior split of supplied PredecessorDAG\n"
                           "routes, EQUAL or WEIGHTED).\n"
                           "OPTIMIZED: the minimum-cost flow core over the routes (a controller restricted to them), re-solved\n"
                           "over the free capacity plus the policy's own allocation (it may move its own earlier placement,\n"
                           "never another's). WEIGHTED (WCMP, weighted-cost multipath): every node splits in proportion to\n"
                           "weights fixed when the routes are prepared (highest averages (D'Hondt): assign each next unit\n"
                           "to the largest weight / (assigned + 1), ties to the lower EdgeId; the weights come from\n"
                           "WeightBasis). EQUAL (ECMP, equal-cost\n"
                           "multipath): every node splits equally, by round robin.")
      .value("EQUAL", RouteAllocation::Equal)
      .value("WEIGHTED", RouteAllocation::Weighted)
      .value("OPTIMIZED", RouteAllocation::Optimized);

  py::enum_<Admission>(m, "Admission",
                       "LOSSLESS admits the largest total the routes carry with no loss at all; LOSSY sends the whole\n"
                       "finite offer hop by hop and records where it is lost (every model except the optimizing ones:\n"
                       "ShortestPaths OPTIMIZED and CentralizedTE without tunnel_count).")
      .value("LOSSLESS", Admission::Lossless)
      .value("LOSSY", Admission::Lossy);

  py::enum_<SourceAllocation>(m, "SourceAllocation",
                    "How a demand's traffic enters at its sources, separately from how it is routed. EQUAL: round\n"
                    "robin over the roster (ascending NodeId, the first sources one unit more); a source whose share\n"
                    "cannot go anywhere keeps it and loses it (LOSSLESS then admits only totals that give it\n"
                    "nothing). OPTIMIZED: the controller chooses the injections (one virtual source).")
      .value("EQUAL", SourceAllocation::Equal)
      .value("OPTIMIZED", SourceAllocation::Optimized);

  py::enum_<WeightBasis>(m, "WeightBasis",
                         "Where WEIGHTED member weights come from (ShortestPaths WEIGHTED and supplied PredecessorDAG routes).\n"
                         "CAPACITY (the default): each member's static capacity. EXPLICIT: edge_weights, one integer per\n"
                         "edge (a positive weight on an edge of capacity 0 becomes 0). RESIDUAL_SNAPSHOT: each member's free\n"
                         "capacity read once when the routes are prepared, against one common view, and frozen for the\n"
                         "epoch. A member of weight 0 is pruned with the branches it leaves without a destination, inside\n"
                         "the route only. The tie rule never chooses weights.")
      .value("CAPACITY", WeightBasis::Capacity)
      .value("EXPLICIT", WeightBasis::Explicit)
      .value("RESIDUAL_SNAPSHOT", WeightBasis::ResidualSnapshot);

  py::enum_<RouteTieBreak>(m, "RouteTieBreak",
                           "How a route choice breaks a cost tie; cost always comes first. FEWER_HOPS: then fewer hops,\n"
                           "then the smaller EdgeId sequence. HIGHER_CAPACITY / HIGHER_AVAILABLE_CAPACITY /\n"
                           "HIGHER_HEADROOM: then the larger path bottleneck of static capacity / availability /\n"
                           "availability over static capacity (exact; 0 at static capacity 0), then fewer hops, then the\n"
                           "EdgeId sequence. Availability: unreserved bandwidth for DistributedTE, the service view for\n"
                           "ShortestPaths, the residual for spf and ksp, static capacity without one. The same order\n"
                           "keeps one of each (u, v) pair's parallel edges when parallel_edges is False. It never\n"
                           "supplies weights and never enables a controller objective.")
      .value("FEWER_HOPS", RouteTieBreak::FewerHops)
      .value("HIGHER_CAPACITY", RouteTieBreak::HigherCapacity)
      .value("HIGHER_AVAILABLE_CAPACITY", RouteTieBreak::HigherAvailableCapacity)
      .value("HIGHER_HEADROOM", RouteTieBreak::HigherHeadroom);

  py::enum_<Direction>(m, "Direction",
                       "The orientation of a search: FORWARD from the roots along the edges, or REVERSE toward the roots\n"
                       "against the edges (the DAG it returns keeps the edges' own direction, so its walks lead to the\n"
                       "roots).")
      .value("FORWARD", Direction::Forward)
      .value("REVERSE", Direction::Reverse);

}

} // namespace netgraph::core::bindings
