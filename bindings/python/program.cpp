// Forwarding programs: TunnelStatus, Tunnel, Route and ForwardingProgram with from_routes.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "netgraph/core/program.hpp"
#include "integers.hpp"
#include "values.hpp"
#include "publication.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

py::typing::Union<Path, PredDAG> bundle_obj(const Route& r) {
  py::object out = std::holds_alternative<std::shared_ptr<const Path>>(r.bundle)
                       ? shared_object(std::get<std::shared_ptr<const Path>>(r.bundle))
                       : shared_object(std::get<std::shared_ptr<const PredDAG>>(r.bundle));
  return py::reinterpret_steal<py::typing::Union<Path, PredDAG>>(out.release());
}

py::tuple ids_tuple(std::span<const std::int32_t> v) {
  py::tuple out(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) out[i] = py::int_(v[i]);
  return out;
}

// A per-edge integer argument of from_routes read without narrowing (its type defects raise
// here); its values are checked later, in the order of the call's validation.
std::optional<WideInts> wide_or_none(const py::object& obj, const char* name) {
  if (obj.is_none()) return std::nullopt;
  (void)sequence_length(obj, name);
  return wide_integers(obj, name);
}

} // namespace

void register_program(py::module_& m) {
  py::enum_<TunnelStatus>(m, "TunnelStatus",
                          "A tunnel's state in one forwarding (a by_tunnel row). UP; FAILED: its path touches a node or\n"
                          "an edge the forwarding policy's Eligibility excludes, so it gets nothing and its head splits over\n"
                          "its other entries (replay); BLACKHOLED: up, but its path reaches a black hole of stale\n"
                          "forwarding, where its traffic is dropped. A program's own tunnels are all up.")
      .value("UP", TunnelStatus::Up)
      .value("FAILED", TunnelStatus::Failed)
      .value("BLACKHOLED", TunnelStatus::Blackholed);

  py::class_<Tunnel>(m, "Tunnel",
                     "One tunnel (a path entry) of a ForwardingProgram, read-only. index: the entry identifier\n"
                     "(Global and CostOnly MAX_FLOW: the path's position in the decomposition order; supplied routes:\n"
                     "the route's position in the list), never a position in program.tunnels; head; local_index: its\n"
                     "position among the head's entries of both kinds; destination: where its path ends; edges: its\n"
                     "EdgeIds in path order; route_cost: the sum of their costs; weight: its integer share among the\n"
                     "head's entries (a controller path: its solved amount).")
      .def_readonly("index", &Tunnel::index)
      .def_readonly("head", &Tunnel::head)
      .def_readonly("local_index", &Tunnel::local_index)
      .def_readonly("destination", &Tunnel::destination)
      .def_property_readonly(
          "edges", [](const Tunnel& t) { return py::tuple(py::cast(t.edges)); }, "the path's EdgeIds in path order")
      .def_readonly("route_cost", &Tunnel::route_cost)
      .def_readonly("weight", &Tunnel::weight)
      .def("__repr__", [](const Tunnel& t) {
        return "Tunnel(index=" + std::to_string(t.index) + ", head=" + std::to_string(t.head) +
               ", weight=" + std::to_string(t.weight) + ")";
      });

  py::class_<Route>(m, "Route",
                    "A supplied route for ForwardingProgram.from_routes: a Path with one edge per hop (a tunnel) or a\n"
                    "PredDAG (a split route) starting at `head`; weight is its integer share among the head's routes\n"
                    "(0 allowed: it carries nothing). from_routes validates it.")
      .def(py::init([](const py::object& head, const py::object& bundle, const py::object& weight) {
             Route r;
             r.head = static_cast<NodeId>(integer_value(head, std::numeric_limits<NodeId>::min(),
                                                        std::numeric_limits<NodeId>::max(), "head"));
             if (py::isinstance<Path>(bundle)) {
               r.bundle = std::shared_ptr<const Path>(bundle.cast<std::shared_ptr<Path>>());
             } else if (py::isinstance<PredDAG>(bundle)) {
               r.bundle = std::shared_ptr<const PredDAG>(bundle.cast<std::shared_ptr<PredDAG>>());
             } else {
               throw py::type_error("Route bundle must be a Path or a PredDAG (got " +
                                    std::string(py::str(py::type::of(bundle).attr("__name__"))) + ")");
             }
             r.weight = integer_value(weight, std::numeric_limits<Amount>::min(), std::numeric_limits<Amount>::max(),
                                      "weight");
             return r;
           }),
           py::arg("head"), py::arg("bundle"), py::kw_only(), py::arg("weight") = 1,
           "Route(head, bundle, *, weight=1): head a NodeId, bundle a Path or a PredDAG (TypeError otherwise),\n"
           "weight an integer (checked in [0, 2^62) by from_routes).")
      .def_readonly("head", &Route::head, "the route's head (a source of the program)")
      .def_property_readonly("bundle", &bundle_obj, "the Path (a tunnel) or the PredDAG (a split route)")
      .def_readonly("weight", &Route::weight, "the route's integer share among the head's routes")
      .def("__repr__", [](const Route& r) {
        return "Route(head=" + std::to_string(r.head) + ", weight=" + std::to_string(r.weight) + ")";
      });

  py::class_<ForwardingProgram, std::shared_ptr<ForwardingProgram>>(
      m, "ForwardingProgram",
      "An immutable value: one demand's routes and integer weights. sources and destinations (the normalized\n"
      "sets); tunnels: the path entries in ascending index (split routes are not listed); source_weights: the\n"
      "captured source supplies ((source, weight) for every roster source, ascending) of the optimized\n"
      "controllers, None for equal sources (round robin over the roster). Each head splits its share by\n"
      "highest averages over its entries' weights; identical routes are never merged. CostOnly and Global\n"
      "produce one (FlowPolicy.program, PlacementResult.program; CapacityAware has none); from_routes builds\n"
      "one from given routes; Pinned forwards one. No public constructor.")
      .def_property_readonly(
          "sources", [](const ForwardingProgram& p) { return ids_tuple(p.sources()); }, "the source set, ascending")
      .def_property_readonly(
          "destinations", [](const ForwardingProgram& p) { return ids_tuple(p.destinations()); },
          "the destination set, ascending")
      .def_property_readonly(
          "tunnels",
          [](const ForwardingProgram& p) {
            const auto t = p.tunnels();
            return tuple_of(std::vector<Tunnel>(t.begin(), t.end()));
          },
          "the program's tunnels (Tunnel), ascending index")
      .def_property_readonly(
          "source_weights",
          [](const ForwardingProgram& p) -> py::object {
            const auto& w = p.source_weights();
            if (!w.has_value()) return py::none();
            py::tuple out(w->size());
            for (std::size_t i = 0; i < w->size(); ++i) out[i] = py::make_tuple((*w)[i].first, (*w)[i].second);
            return out;
          },
          "the captured source weights ((source, weight), ascending), or None for equal sources")
      .def_static(
          "from_routes",
          [](const std::shared_ptr<Graph>& graph, const py::object& sources, const py::object& destinations,
             const py::object& routes, std::optional<FlowPlacement> placement, std::optional<WeightBasis> weights,
             const py::object& edge_weights, const py::object& residual) {
            if (!graph) throw py::value_error("graph must be a Graph");
            // Conversion first (TypeError and integer type defects), then the validation order of
            // from_routes; the values beyond int64 are checked here in that order.
            if (not_a_sequence(routes)) throw py::type_error("routes must be a sequence of Route");
            std::vector<Route> rs;
            for (const auto& item : py::reinterpret_borrow<py::sequence>(routes)) {
              if (!py::isinstance<Route>(item)) throw py::type_error("routes must be a sequence of Route");
              rs.push_back(item.cast<Route>());
            }
            const auto ew = wide_or_none(edge_weights, "edge_weights");
            const auto rw = wide_or_none(residual, "residual");
            const auto n = graph->num_nodes();
            const auto s = detail::node_set(node_ids(sources, n, "sources"), n, "sources");
            const auto d = detail::node_set(node_ids(destinations, n, "destinations"), n, "destinations");
            detail::check_disjoint(s, d, "sources", "destinations");
            if (rs.empty()) throw py::value_error("routes must not be empty");
            bool any_dag = false;
            for (const auto& r : rs) any_dag = any_dag || std::holds_alternative<std::shared_ptr<const PredDAG>>(r.bundle);
            detail::check_route_options(any_dag, placement, weights, ew.has_value(), rw.has_value());
            RouteOptions options;
            options.placement = placement;
            options.weights = weights;
            const auto E = static_cast<std::size_t>(graph->num_edges());
            if (ew.has_value()) {
              if (ew->size() != E) throw py::value_error("edge_weights: expected " + std::to_string(E) + " entries");
              for (std::size_t e = 0; e < E; ++e) {
                if (!ew->in(e, 0, kDomainLimit - 1)) {
                  throw py::value_error("edge_weights[" + std::to_string(e) + "] must be an integer in [0, 2^62)");
                }
              }
              options.edge_weights = ew->narrow<Amount>();
            }
            if (rw.has_value()) options.residual = edge_amounts(residual, *graph, "residual", "residual (units)").narrow<Amount>();
            return typed_shared(ForwardingProgram::from_routes(*graph, s, d, rs, options));
          },
          py::arg("graph"), py::arg("sources"), py::arg("destinations"), py::arg("routes"), py::kw_only(),
          py::arg("placement") = py::none(), py::arg("weights") = py::none(), py::arg("edge_weights") = py::none(),
          py::arg("residual") = py::none(),
          "Validate and freeze supplied routes -> ForwardingProgram (ValueError for the first failure, in this\n"
          "order): the endpoint sets; `routes must not be empty`; when no route is a PredDAG, any option given\n"
          "(`<placement, weights, edge_weights or residual> applies only to PredDAG routes`); placement MAX_FLOW\n"
          "(`from_routes placement must be EQUAL_BALANCED or WEIGHTED`); weights only with placement=WEIGHTED,\n"
          "EXPLICIT with edge_weights, RESIDUAL_SNAPSHOT with residual; edge_weights' length and entries in [0,\n"
          "2^62) and residual's (one entry per edge in [0, capacity]); then each route in list order: its\n"
          "bundle's graph, its head a source, a Path with one edge per hop that starts at the head, ends at a\n"
          "destination and meets none before, a PredDAG with a walk from the head to a destination, its weight\n"
          "in [0, 2^62). Order and duplicates are kept; a route of weight 0 is validated in full. placement\n"
          "(None: EQUAL_BALANCED) is the interior split of the PredDAG routes; each split route is its own\n"
          "split context (a router splits only that route's share). The roster is the sources with a route.");
}

} // namespace netgraph::core::bindings
