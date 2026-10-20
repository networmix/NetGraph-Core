// FlowPolicy: the one placement engine, bound to its ledger.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/shortest_path_cache.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/forwarding_program.hpp"
#include "integers.hpp"
#include "publication.hpp"
#include "values.hpp"
#if NETGRAPH_CORE_PRIVATE_HOOKS
#include <functional>
#include <variant>

#include "netgraph/core/detail/route_builders.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "trace.hpp"
#endif
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// A request or a full offer: units in [0, 2^62) or +inf (float('inf')); nullopt is +inf.
std::optional<Amount> request_units(const py::object& offer) {
  return units_or_inf(offer, "offer must be an integer in [0, 2^62) or permitted +inf");
}

// A design demand argument: None absent; float('inf') the inner nullopt; an integer passed on for
// the call's own checks in their order (one outside int64 as -1, outside the domain like it); any
// other value raises the domain text here (a conversion failure comes first).
std::optional<std::optional<Amount>> design_units(const py::object& design) {
  constexpr const char* text = "design_demand must be an integer in [0, 2^62) or permitted +inf";
  if (design.is_none()) return std::nullopt;
  if (is_bool(design) || !PyIndex_Check(design.ptr())) return units_or_inf(design, text);
  const py::object idx = python_int(design, "design_demand");
  int overflow = 0;
  const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
  return std::optional<Amount>(overflow != 0 ? Amount {-1} : static_cast<Amount>(v));
}

#if NETGRAPH_CORE_PRIVATE_HOOKS

// The epoch state of a policy (roster, requested, +inf, per roster source the flows).
py::dict state_obj(const FlowPolicy& policy) {
  const auto& pe = detail::PolicyAccess::engine(policy);
  py::dict out;
  const auto& r = pe.config.routing;
  out["mode"] = detail::model_name(r);
  out["roster"] = pe.reachable_sources;
  out["requested"] = counter_int(pe.epoch.requested);
  out["infinite"] = pe.epoch.infinite;
  out["placed_demand"] = policy.delivered();
  py::dict cohorts;
  // Per roster source its flows (T1, T2 and Pinned: per program entry, keyed by its position).
  const auto* gl = std::get_if<CentralizedTE>(&r);
  const bool per_entry = std::holds_alternative<Pinned>(r) || std::holds_alternative<DistributedTE>(r) ||
                         (gl != nullptr && gl->tunnel_count.has_value());
  for (std::size_t i = 0; i < pe.epoch.record_groups.size(); ++i) {
    py::list flows;
    for (const auto& s : pe.epoch.record_groups[i]) {
      const FlowRecord& rec = detail::LedgerAccess::live(pe.ledger, s.record_counter);
      py::dict f;
      f["flow_id"] = py::make_tuple(rec.flow_id.first, rec.flow_id.second);
      f["src"] = s.source;
      const auto edges = detail::PredecessorDAGBuilder::ascending(*s.route);
      f["route"] = std::vector<EdgeId>(edges.begin(), edges.end());
      f["T"] = rec.delivered;
      f["delivered"] = rec.delivered;
      f["alloc"] = pairs_list(rec.edge_flows);
      f["entry"] = s.entry_id < 0 ? py::object(py::none()) : py::object(py::int_(s.entry_id));
      flows.append(f);
    }
    cohorts[py::int_(per_entry ? static_cast<std::int64_t>(i) : static_cast<std::int64_t>(pe.reachable_sources[i]))] = flows;
  }
  out["cohorts"] = cohorts;
  return out;
}

#endif

} // namespace

void register_flow_policy(py::module_& m) {
  [[maybe_unused]] auto flow_policy = py::class_<FlowPolicy>(
      m, "FlowPolicy",
      "The placement engine for one demand, bound to its ledger at construction. Every call (place,\n"
      "evaluate, replace, forward) is all-or-nothing: a failed call raises and leaves the ledger and the\n"
      "policy as they were (with the two limits in the module docstring). remove() takes every flow of the\n"
      "policy back out, clears its program and starts a new epoch (the request total starts again from 0);\n"
      "replace and forward start one too. A policy dropped without remove()\n"
      "leaves its records in the ledger as orphans: readable, never changed again, cleared by\n"
      "FlowLedger.reset(). sources and destinations are one NodeId or a sequence; the two sets must not\n"
      "overlap. ValueError for invalid endpoints, configurations and requests and for objects of another\n"
      "graph, ledger or policy; TypeError for a ledger (None included), config, eligibility, cache or\n"
      "outputs of the wrong class, a node id or a request of a wrong type (a bool, a float node id, a\n"
      "string) and an unknown keyword.")
      .def(py::init([](const std::shared_ptr<FlowLedger>& ledger, const py::object& sources,
                       const py::object& destinations, const FlowPolicyConfig& config,
                       const std::shared_ptr<Eligibility>& elig, const std::shared_ptr<ShortestPathCache>& cache) {
             if (!ledger) throw py::type_error("ledger must be a FlowLedger");
             const auto n = ledger->graph()->num_nodes();
             const auto s = node_ids(sources, n, "sources");
             const auto d = node_ids(destinations, n, "destinations");
             return std::make_unique<FlowPolicy>(*ledger, s, d, config, elig.get(), cache.get());
           }),
           py::keep_alive<1, 2>(), py::keep_alive<1, 6>(), py::keep_alive<1, 7>(), py::arg("ledger"),
           py::arg("sources"), py::arg("destinations"), py::arg("config"), py::arg("eligibility") = nullptr,
           py::arg("cache") = nullptr,
           "Bind a policy to ledger. eligibility (None: everything allowed) and cache (None: none; only ShortestPaths\n"
           "routing reads it, another model rejects it) belong to the ledger's graph; the policy keeps the\n"
           "ledger, eligibility and cache alive (a Pinned policy shares its program and black holes). The roster\n"
           "(the sources that statically reach an allowed destination, over allowed edges whatever their\n"
           "capacity, the tunnel models' heads; Pinned: the program's heads that reach an allowed destination over its basis\n"
           "under eligibility) is fixed here. Everything is validated before the policy binds, in this order: the endpoints, the\n"
           "graph of eligibility, the program, the support, the cache and the black holes (`foreign Eligibility:\n"
           "built for another graph`, `program belongs to another graph`, `support belongs to another graph`,\n"
           "`foreign ShortestPathCache: bound to another graph`, `black_holes belongs to another graph`), the\n"
           "configuration (FlowPolicyConfig's order, the pooled single-path rule and T2's several-sources rule\n"
           "included), the tunnel limit (`tunnel limit exceeded:\n"
           "<H> heads x <N> tunnels > 1073741823`) and a Pinned program's endpoints (`program endpoints do not\n"
           "match policy endpoints`).")
      .def(
          "place",
          [](FlowPolicy& p, const py::object& offer, const py::object& outputs) {
            // The operation comes before the offer (the call's order).
            detail::check_place_operation(detail::PolicyAccess::engine(p));
            const auto v = request_units(offer);
            const auto o = outputs_of(outputs);
            // All or nothing across the boundary: the result's Python object exists before the
            // call, and after the commit nothing allocates (the value moves into it). The GIL
            // stays held: the call mutates the ledger (one writer per ledger).
            auto [obj, slot] = new_object<PlacementResult>();
            *slot = p.place(v, o);
            return std::move(obj);
          },
          py::arg("offer"), py::arg("outputs") = py::none(),
          "Place offer (int units in [0, 2^62), or float('inf')) all or nothing -> PlacementResult. T1 and T2\n"
          "raise `<DistributedTE or CentralizedTE with tunnel_count> does not support incremental placement`\n"
          "(also for 0).\n"
          "Requests add up within an epoch (from construction, the last remove(), replace() or forward());\n"
          "after float('inf') any finite request is rejected, 0 included. A request of 0 places nothing and\n"
          "reports the current state. ShortestPaths prepares its program at the epoch's first positive request;\n"
          "LOSSLESS admits the largest total in [admitted, admitted + offer] that the program carries with no\n"
          "loss (a refused unit is not offered again); LOSSY forwards the new round-robin (or captured-weight)\n"
          "injections, each head splitting its increment afresh (earlier drops are not offered again); a\n"
          "controller re-solves admitted + offer. A LOSSY float('inf') is ShortestPaths's first-hop experiment;\n"
          "Pinned raises `LOSSY forwarding requires a finite offer`. outputs (OutputOptions; None: none)\n"
          "requests edge_flows, residual, saturated_boundary_edges, edge_reservations and min_cut. Everything is\n"
          "validated before any change; a\n"
          "failed call raises with the ledger and the policy unchanged. A solve of the flow core (CentralizedTE\n"
          "or ShortestPaths OPTIMIZED) on a network beyond its int32 numbering raises ValueError:\n"
          "more than INT32_MAX - 2 nodes, or more than 2^30 - 1 edges, sources and destinations together. A\n"
          "LOSSY forwarding whose losses wait on each other in a cycle raises ValueError when its rounds\n"
          "alternate between two states or reach the round limit, RuntimeError when the result's check fails.")
      .def("remove", &FlowPolicy::remove,
           "revert every flow and hold of the policy, retire the ids, clear the program, open a new epoch")
      .def("delivered", &FlowPolicy::delivered, "the sum of delivered over the policy's live records")
      .def(
          "records", [](const FlowPolicy& p) { return list_of(p.records()); },
          "the policy's live records as snapshots in FlowId order")
      .def(
          "record",
          [](const FlowPolicy& p, const py::object& flow_id) { return typed_object(p.record(flow_id_of(flow_id))); },
          py::arg("flow_id"),
          "one of the policy's live records by FlowId (ValueError for a FlowId of another ledger or policy, a\n"
          "retired or an unknown id)")
      .def(
          "evaluate",
          [](const FlowPolicy& p, const py::object& offer, const py::object& design_demand, const py::object& outputs) {
            const auto v = request_units(offer);
            const auto d = design_units(design_demand);
            const auto o = outputs_of(outputs);
            auto [obj, slot] = new_object<PlacementResult>();
            *slot = p.evaluate(v, d, o);
            return std::move(obj);
          },
          py::arg("offer"), py::kw_only(), py::arg("design_demand") = py::none(), py::arg("outputs") = py::none(),
          "A full, independent computation and forwarding of offer (int units or float('inf')) from the\n"
          "policy's own released view (the residual plus its own carried load; T1 sets up on the unreserved\n"
          "capacity with its own holds released; T2 lays out on the released residual), as replace would do it -> PlacementResult (the computed\n"
          "program and plan, the flags replace would set, the outputs of the evaluated state); publishes\n"
          "nothing, creates no record, books no hold and leaves the epoch as it was. design_demand as for\n"
          "replace. LOSSY requires a finite offer.")
      .def(
          "replace",
          [](FlowPolicy& p, const py::object& offer, const py::object& design_demand, const py::object& outputs) {
            const auto v = request_units(offer);
            const auto d = design_units(design_demand);
            const auto o = outputs_of(outputs);
            auto [obj, slot] = new_object<PlacementResult>();
            *slot = p.replace(v, d, o);
            return std::move(obj);
          },
          py::arg("offer"), py::kw_only(), py::arg("design_demand") = py::none(), py::arg("outputs") = py::none(),
          "Recompute from the policy's own released view and forward the full offer, published atomically ->\n"
          "PlacementResult: the old records and holds are retired (their ids never live again), new ones\n"
          "created, and a new epoch starts at this offer (delivery may go down: delivery_delta). ShortestPaths\n"
          "prepares its routes again; a controller re-solves; T1 sets its tunnels up again for the design and\n"
          "books their holds; T2 lays its tunnels out again for the design (planned bandwidth, never held);\n"
          "Pinned forwards its program. design_demand (T1 and T2: int units, or float('inf'); None: equal to\n"
          "the offer): T1 needs a finite design (`DistributedTE requires a finite design demand`), T2 takes an unlimited\n"
          "one (its largest common size); every other model raises `<model> does not accept design_demand`.\n"
          "In this order: the offer, LOSSY's finite offer, then the design demand (the model, its domain,\n"
          "finiteness).")
      .def(
          "forward",
          [](FlowPolicy& p, const py::object& offer, const py::object& black_holes, const py::object& outputs) {
            // The operation and the installed program come before the offer (the call's order).
            detail::check_forward_operation(detail::PolicyAccess::engine(p));
            const auto v = request_units(offer);
            const auto bh = eligibility_or_none(black_holes, "black_holes");
            const auto o = outputs_of(outputs);
            auto [obj, slot] = new_object<PlacementResult>();
            *slot = p.forward(v, bh.get(), o);
            return std::move(obj);
          },
          py::arg("offer"), py::kw_only(), py::arg("black_holes") = py::none(), py::arg("outputs") = py::none(),
          "Forward a full new offer over the installed program, with no route computation and no setup, on the\n"
          "released view -> PlacementResult; a new epoch like replace, the program, its plan and T1's holds\n"
          "kept (no flag; T1's records are recreated with the same edge_reservations). LOSSLESS admits the\n"
          "program's lossless bound up to the offer; LOSSY sends the offer. black_holes (an Eligibility of the\n"
          "graph): the nodes and edges it excludes are black holes for this call (stale forwarding: traffic\n"
          "reaching one is dropped there; a black-holed source loses its share). In this order: Pinned raises\n"
          "`Pinned forwards with replace`; no installed program `forward requires an installed program`;\n"
          "the offer; LOSSY +inf `LOSSY forwarding requires a finite offer`; black_holes\n"
          "of another graph `black_holes belongs to another graph`.")
      .def_property_readonly(

          "program",
          [](const FlowPolicy& p) { return optional_shared(p.program()); },
          "the installed ForwardingProgram; None before the first preparation and after remove()")
      .def_property_readonly("config", [](const FlowPolicy& p) { return typed_object(p.config()); },
                             "the configuration as validated at construction, ShortestPaths' absent\n"
                             "source_allocation and weight_basis resolved to the allocation's defaults");
#if NETGRAPH_CORE_PRIVATE_HOOKS
  flow_policy
      .def(
          "_place_traced",
          [](FlowPolicy& p, const py::object& offer, const py::object& outputs) {
            detail::check_place_operation(detail::PolicyAccess::engine(p));
            const auto v = request_units(offer);
            const auto o = outputs_of(outputs);
            // As place(): every Python object exists before anything changes; the engine runs the
            // `staged` callback, which converts the trace, after the work is done and before the
            // commit, so a failed conversion changes nothing.
            auto [obj, slot] = new_object<PlacementResult>();
            auto pair = py::reinterpret_steal<py::tuple>(checked(PyTuple_New(2)).release());
            py::object trace_py;
            detail::Trace trace;
            const std::function<void()> staged = [&] { trace_py = trace_obj(trace); };
            *slot = detail::place_traced(p, v, o, trace, &staged);
            PyTuple_SET_ITEM(pair.ptr(), 0, obj.release().ptr());
            PyTuple_SET_ITEM(pair.ptr(), 1, trace_py.release().ptr());
            return pair;
          },
          py::arg("offer"), py::arg("outputs") = py::none(),
          "private: place() with the transaction's trace (the call's intermediate values, keyed as\n"
          "the tests' model of the policy keys them)")
      .def("_state", &state_obj, "private: the epoch state (roster, requested, +inf, cohorts)")
      .def(
          "_owner", [](const FlowPolicy& p) { return detail::PolicyAccess::engine(p).policy_id; },
          "private: the ledger-local owner id (FlowRecord.policy_id)");
#endif
}

} // namespace netgraph::core::bindings
