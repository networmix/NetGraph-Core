// Forwarding programs: TunnelStatus, Tunnel, RouteDefinition and ForwardingProgram with from_routes and from_dag.
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
#include "netgraph/core/route_structures.hpp"
#include "netgraph/core/forwarding_program.hpp"
#include "integers.hpp"
#include "values.hpp"
#include "native_call.hpp"
#include "publication.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

py::typing::Union<Path, PredecessorDAG> route_obj(const RouteDefinition& r) {
  py::object out = std::holds_alternative<std::shared_ptr<const Path>>(r.route)
                       ? shared_object(std::get<std::shared_ptr<const Path>>(r.route))
                       : shared_object(std::get<std::shared_ptr<const PredecessorDAG>>(r.route));
  return py::reinterpret_steal<py::typing::Union<Path, PredecessorDAG>>(out.release());
}


// A per-edge integer argument of from_routes read without narrowing (its type defects raise
// here); its values are checked later, in the order of the call's validation.
std::optional<CheckedIntegerVector> wide_or_none(const py::object& obj, const char* name) {
  if (obj.is_none()) return std::nullopt;
  (void)sequence_length(obj, name);
  return wide_integers(obj, name);
}

// Validate and narrow captured weights at the factory's weight-validation step.
std::vector<Amount> checked_edge_weights(const CheckedIntegerVector& weights, std::size_t num_edges) {
  if (weights.size() != num_edges) throw py::value_error("edge_weights: expected " + std::to_string(num_edges) + " entries");
  for (std::size_t e = 0; e < num_edges; ++e) {
    if (!weights.in(e, 0, kDomainLimit - 1)) {
      throw py::value_error("edge_weights[" + std::to_string(e) + "] must be an integer in [0, 2^62)");
    }
  }
  return weights.narrow<Amount>();
}

} // namespace

void register_program(py::module_& m) {
  py::enum_<TunnelStatus>(m, "TunnelStatus",
                          "A tunnel's state. UP; INSUFFICIENT_BANDWIDTH and COST_CEILING_EXCEEDED: T1's setup left it down (no path\n"
                          "passes the bandwidth filter, or the cheapest exceeds the cost ceiling): no path, weight 0, and it\n"
                          "stays down in every forwarding of its program; and in one forwarding only (a by_tunnel row)\n"
                          "FAILED: its path touches a node or an edge the forwarding policy's Eligibility excludes, so it\n"
                          "gets nothing and its head splits over its other entries (replay), and BLACKHOLED: up, but its\n"
                          "path reaches a black hole of stale forwarding, where its traffic is dropped.")
      .value("UP", TunnelStatus::Up)
      .value("INSUFFICIENT_BANDWIDTH", TunnelStatus::InsufficientBandwidth)
      .value("COST_CEILING_EXCEEDED", TunnelStatus::CostCeilingExceeded)
      .value("FAILED", TunnelStatus::Failed)
      .value("BLACKHOLED", TunnelStatus::Blackholed);

  py::class_<Tunnel>(m, "Tunnel",
                     "One tunnel (a path entry) of a ForwardingProgram, read-only. entry_id: the entry identifier (T1 and\n"
                     "T2: j*H + i for round j and roster position i, T2's dead heads leaving holes; CentralizedTE and ShortestPaths\n"
                     "OPTIMIZED: the path's position in the decomposition order; supplied routes: the route's position in\n"
                     "the list), never a position in program.tunnels; head; local_index: its position among the head's entries of both kinds;\n"
                     "destination: where its path ends; edges: its EdgeIds in path order; path_cost: the sum of their\n"
                     "costs (all three None for a tunnel that is down); status: UP, INSUFFICIENT_BANDWIDTH or\n"
                     "COST_CEILING_EXCEEDED; requested_bandwidth and reserved_bandwidth: T1's strict bandwidth (its request,\n"
                     "and what it holds on every link of its path: the request when UP, else 0; 0 for the other models);\n"
                     "planned_bandwidth: T2's planned bandwidth\n"
                     "u (metadata, never held; 0 for the other models); weight: its integer share among the head's\n"
                     "entries (T1: 1 up, 0 down; T2: 1; a controller path: its solved amount).")
      .def_property_readonly("entry_id", int_reader(&Tunnel::entry_id))
      .def_property_readonly("head", int_reader(&Tunnel::head))
      .def_property_readonly("local_index", int_reader(&Tunnel::local_index))
      .def_property_readonly("destination", optional_int_reader(&Tunnel::destination))
      .def_property_readonly(

          "edges",
          [](const Tunnel& t) -> py::object {
            return t.edges.has_value() ? py::object(ints_tuple(*t.edges)) : py::object(py::none());
          },
          "the path's EdgeIds in path order; None for a tunnel that is down")
      .def_property_readonly("path_cost", optional_int_reader(&Tunnel::path_cost))
      .def_property_readonly("status", [](const Tunnel& t) { return typed_object(t.status); })
      .def_property_readonly("requested_bandwidth", int_reader(&Tunnel::requested_bandwidth))
      .def_property_readonly("reserved_bandwidth", int_reader(&Tunnel::reserved_bandwidth))
      .def_property_readonly("planned_bandwidth", int_reader(&Tunnel::planned_bandwidth))
      .def_property_readonly("weight", int_reader(&Tunnel::weight))
      .def("__repr__", [](const Tunnel& t) {
        return "Tunnel(entry_id=" + std::to_string(t.entry_id) + ", head=" + std::to_string(t.head) +
               ", weight=" + std::to_string(t.weight) + ")";
      });

  py::class_<RouteDefinition>(m, "RouteDefinition",
                    "A supplied route for ForwardingProgram.from_routes: a Path with one edge per hop (a tunnel) or a\n"
                    "PredecessorDAG (a split route) starting at `head`; weight is its integer share among the head's routes\n"
                    "(0 allowed: it carries nothing). from_routes validates it.")
      .def(py::init([](const py::object& head, const py::object& route, const py::object& weight) {
             RouteDefinition r;
             r.head = static_cast<NodeId>(
                 integer_value(head, std::numeric_limits<NodeId>::min(), std::numeric_limits<NodeId>::max(), "head"));
             if (py::isinstance<Path>(route)) {
               r.route = std::shared_ptr<const Path>(route.cast<std::shared_ptr<Path>>());
             } else if (py::isinstance<PredecessorDAG>(route)) {
               r.route = std::shared_ptr<const PredecessorDAG>(route.cast<std::shared_ptr<PredecessorDAG>>());
             } else {
               throw py::type_error("RouteDefinition route must be a Path or a PredecessorDAG (got " +
                                    std::string(py::str(py::type::of(route).attr("__name__"))) + ")");
             }
             r.weight = integer_value(weight, std::numeric_limits<Amount>::min(), std::numeric_limits<Amount>::max(), "weight");
             return r;
           }),
           py::arg("head"), py::arg("route"), py::kw_only(), py::arg("weight") = 1,
           "RouteDefinition(head, route, *, weight=1): head a NodeId, route a Path or a PredecessorDAG (TypeError\n"
           "otherwise), weight an integer (checked in [0, 2^62) by from_routes).")
      .def_property_readonly("head", int_reader(&RouteDefinition::head), "the route's head (a source of the program)")
      .def_property_readonly("route", &route_obj, "the Path (a tunnel) or the PredecessorDAG (a split route)")
      .def_property_readonly("weight", int_reader(&RouteDefinition::weight), "the route's integer share among the head's routes")
      .def("__repr__", [](const RouteDefinition& r) {
        return "RouteDefinition(head=" + std::to_string(r.head) + ", weight=" + std::to_string(r.weight) + ")";
      });

  py::class_<ForwardingProgram, std::shared_ptr<ForwardingProgram>>(
      m, "ForwardingProgram",
      "An immutable value: one demand's routes and integer weights. sources and destinations (the normalized\n"
      "sets); tunnels: the path entries in ascending entry_id (split routes are not listed); source_weights: the\n"
      "captured weights ((source, weight) for every roster source, ascending): optimized controllers' solved\n"
      "supplies, weighted CentralizedTE's declared weights, or None for equal sources (round robin). Each head splits its share by\n"
      "highest averages over its entries' weights; identical routes are never merged. ShortestPaths and CentralizedTE\n"
      "produce one (FlowPolicy.program, PlacementResult.program), and so do T1 and T2 (their tunnels);\n"
      "from_routes builds one from given routes and from_dag one shared hop-by-hop program on a DAG; Pinned\n"
      "forwards one. No public constructor.")
      .def_property_readonly(

          "sources", [](const ForwardingProgram& p) { return ints_tuple(p.sources()); }, "the source set, ascending")
      .def_property_readonly(

          "destinations", [](const ForwardingProgram& p) { return ints_tuple(p.destinations()); },
          "the destination set, ascending")
      .def_property_readonly(

          "tunnels",
          [](const ForwardingProgram& p) {
            const auto t = p.tunnels();
            return tuple_of(std::vector<Tunnel>(t.begin(), t.end()));
          },
          "the program's tunnels (Tunnel), ascending entry_id")
      .def_property_readonly(

          "source_weights",
          [](const ForwardingProgram& p) -> py::object {
            const auto& w = p.source_weights();
            if (!w.has_value()) return py::none();
            return pairs_tuple(*w);
          },
          "the captured source weights ((source, weight), ascending), or None for equal sources")
      .def_static(
          "from_routes",
          [](const std::shared_ptr<Graph>& graph, const py::object& sources, const py::object& destinations,
             const py::object& routes, std::optional<RouteAllocation> allocation, std::optional<WeightBasis> weight_basis,
             const py::object& edge_weights, const py::object& residual) {
            if (!graph) throw py::type_error("graph must be a Graph");
            // Conversion first (TypeError and integer type defects), then the validation order of
            // from_routes; the values beyond int64 are checked here in that order.
            if (is_invalid_sequence(routes)) throw py::type_error("routes must be a sequence of RouteDefinition");
            std::vector<RouteDefinition> rs;
            for (const auto& item : py::reinterpret_borrow<py::sequence>(routes)) {
              if (!py::isinstance<RouteDefinition>(item)) throw py::type_error("routes must be a sequence of RouteDefinition");
              rs.push_back(item.cast<RouteDefinition>());
            }
            const auto ew = wide_or_none(edge_weights, "edge_weights");
            const auto rw = wide_or_none(residual, "residual");
            const auto n = graph->num_nodes();
            const auto s = detail::node_set(node_ids(sources, n, "sources"), n, "sources");
            const auto d = detail::node_set(node_ids(destinations, n, "destinations"), n, "destinations");
            detail::check_disjoint(s, d, "sources", "destinations");
            if (rs.empty()) throw py::value_error("routes must not be empty");
            bool any_dag = false;
            for (const auto& r : rs) any_dag = any_dag || std::holds_alternative<std::shared_ptr<const PredecessorDAG>>(r.route);
            detail::check_route_options("from_routes", any_dag, allocation, weight_basis, ew.has_value(), rw.has_value());
            RouteOptions options;
            options.allocation = allocation;
            options.weight_basis = weight_basis;
            if (ew.has_value()) options.edge_weights = checked_edge_weights(*ew, static_cast<std::size_t>(graph->num_edges()));
            if (rw.has_value()) options.residual = edge_amounts(residual, *graph, "residual", "residual (units)").narrow<Amount>();
            auto checked = native_call<false>([&] { return detail::check_route_definitions(*graph, s, d, rs, options); });
            return typed_shared(native_call([&] {
              auto prepared = std::move(checked);
              return detail::ProgramAccess::from_checked_routes(s, d, rs, std::move(prepared));
            }));
          },
          py::arg("graph"), py::arg("sources"), py::arg("destinations"), py::arg("routes"), py::kw_only(),
          py::arg("allocation") = py::none(), py::arg("weight_basis") = py::none(), py::arg("edge_weights") = py::none(),
          py::arg("residual") = py::none(),
          "Validate and freeze supplied routes -> ForwardingProgram. A wrong type raises TypeError; structural,\n"
          "range and combination defects raise ValueError, in this order: the endpoint sets; `routes must not be empty`; when no route is a PredecessorDAG, any option given\n"
          "(`<allocation, weight_basis, edge_weights or residual> applies only to PredecessorDAG routes`); allocation\n"
          "OPTIMIZED (`from_routes allocation must be EQUAL or WEIGHTED`); weight_basis only with allocation=WEIGHTED,\n"
          "EXPLICIT with edge_weights, RESIDUAL_SNAPSHOT with residual; edge_weights' length and entries in [0,\n"
          "2^62) and residual's (one entry per edge in [0, capacity]); then each route in list order: its\n"
          "route's graph, its head a source, a Path with one edge per hop that starts at the head, ends at a\n"
          "destination and meets none before, a PredecessorDAG with a walk from the head to a destination, its weight\n"
          "in [0, 2^62). Order and duplicates are kept; a route of weight 0 is validated in full. allocation\n"
          "(None: EQUAL) is the interior split of the PredecessorDAG routes; each split route is its own\n"
          "split context (a router splits only that route's share). The roster is the sources with a route.")
      .def_static(
          "from_dag",
          [](const std::shared_ptr<Graph>& graph, const py::object& sources, const py::object& destinations,
             const py::object& dag, std::optional<RouteAllocation> allocation, std::optional<WeightBasis> weight_basis,
             const py::object& edge_weights, const py::object& residual) {
            if (!graph) throw py::type_error("graph must be a Graph");
            // Conversion first (TypeError and integer type defects), then the validation order of
            // from_dag; the values beyond int64 are checked here in that order.
            if (!py::isinstance<PredecessorDAG>(dag)) {
              throw py::type_error("dag must be a PredecessorDAG (got " + std::string(py::str(py::type::of(dag).attr("__name__"))) +
                                   ")");
            }
            const auto shared_dag = dag.cast<std::shared_ptr<PredecessorDAG>>();
            const auto ew = wide_or_none(edge_weights, "edge_weights");
            const auto rw = wide_or_none(residual, "residual");
            const auto n = graph->num_nodes();
            const auto s = detail::node_set(node_ids(sources, n, "sources"), n, "sources");
            const auto d = detail::node_set(node_ids(destinations, n, "destinations"), n, "destinations");
            detail::check_disjoint(s, d, "sources", "destinations");
            if (shared_dag->graph_uid() != graph->uid()) throw py::value_error("dag belongs to another graph");
            detail::check_route_options("from_dag", true, allocation, weight_basis, ew.has_value(), rw.has_value());
            RouteOptions options;
            options.allocation = allocation;
            options.weight_basis = weight_basis;
            if (ew.has_value()) options.edge_weights = checked_edge_weights(*ew, static_cast<std::size_t>(graph->num_edges()));
            if (rw.has_value()) options.residual = edge_amounts(residual, *graph, "residual", "residual (units)").narrow<Amount>();
            return typed_shared(native_call([&] { return ForwardingProgram::from_dag(*graph, s, d, shared_dag, options); }));
          },
          py::arg("graph"), py::arg("sources"), py::arg("destinations"), py::arg("dag"), py::kw_only(),
          py::arg("allocation") = py::none(), py::arg("weight_basis") = py::none(), py::arg("edge_weights") = py::none(),
          py::arg("residual") = py::none(),
          "One shared hop-by-hop program on a supplied PredecessorDAG of the graph -> ForwardingProgram: the sources\n"
          "allocated equally and pooled at every shared router exactly as ShortestPaths EQUAL/WEIGHTED (one\n"
          "split context), the DAG trimmed to the walks from each source to its first destinations. A wrong type\n"
          "raises TypeError; structural, range and combination defects raise ValueError, in this order: the endpoint sets; `dag belongs to another graph`; allocation OPTIMIZED\n"
          "(`from_dag allocation must be EQUAL or WEIGHTED`); weight_basis only with allocation=WEIGHTED, EXPLICIT\n"
          "with edge_weights, RESIDUAL_SNAPSHOT with residual; edge_weights' length and entries in [0, 2^62) and\n"
          "residual's (one entry per edge in [0, capacity]). A dag of another type raises TypeError. The roster is\n"
          "the sources that reach a destination in the DAG (possibly none): one split route per roster source,\n"
          "identified by roster position, no tunnel; a replay rebuilds the roster on the DAG.");
}

} // namespace netgraph::core::bindings
