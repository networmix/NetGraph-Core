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

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "integers.hpp"
#include "publication.hpp"
#include "values.hpp"
#if NETGRAPH_CORE_PRIVATE_HOOKS
#include <functional>
#include <variant>

#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "trace.hpp"
#endif
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// A request: an integer number of units in [0, 2^62) or +inf (float('inf')); nullopt is +inf.
// Any other value (a bool, another float, a non-integer, an integer out of range) raises the one
// ValueError; an exception of the object's own __index__ propagates as raised.
std::optional<Amount> request_units(const py::object& volume) {
  const char* const text = "offer must be an integer in [0, 2^62) or permitted +inf";
  if (py::isinstance<py::float_>(volume) || py::isinstance(volume, py::module_::import("numpy").attr("floating"))) {
    const double x = PyFloat_AsDouble(volume.ptr());
    if (x == -1.0 && PyErr_Occurred() != nullptr) throw py::error_already_set();
    if (std::isinf(x) && x > 0) return std::nullopt;
    throw py::value_error(text);
  }
  if (is_bool(volume) || !PyIndex_Check(volume.ptr())) throw py::value_error(text);
  const py::object idx = py::reinterpret_steal<py::object>(PyNumber_Index(volume.ptr()));
  if (!idx) throw py::error_already_set();
  int overflow = 0;
  const long long v = PyLong_AsLongLongAndOverflow(idx.ptr(), &overflow);
  if (overflow != 0 || v < 0 || v >= kDomainLimit) throw py::value_error(text);
  return static_cast<Amount>(v);
}

#if NETGRAPH_CORE_PRIVATE_HOOKS

// The epoch state of a policy (roster, requested, +inf, per roster source the flows).
py::dict state_obj(const FlowPolicy& policy) {
  const auto& pe = detail::PolicyAccess::engine(policy);
  py::dict out;
  const auto& r = pe.config.routing;
  out["mode"] = std::holds_alternative<CostOnly>(r)        ? "CostOnly"
                : std::holds_alternative<CapacityAware>(r) ? "CapacityAware"
                : std::holds_alternative<Global>(r)        ? "Global"
                                                           : "Pinned";
  out["roster"] = pe.roster;
  out["requested"] = counter_int(pe.requested);
  out["infinite"] = pe.infinite;
  out["placed_demand"] = policy.placed_demand();
  py::dict cohorts;
  for (std::size_t i = 0; i < pe.roster.size(); ++i) {
    py::list flows;
    for (const auto& s : pe.cohorts[i]) {
      const FlowRecord& rec = detail::LedgerAccess::live(pe.fg, s.id);
      py::dict f;
      f["flow_id"] = py::make_tuple(rec.flow_id.first, rec.flow_id.second);
      f["src"] = s.src;
      if (s.plan) {
        f["route"] = s.plan->route;
      } else {
        const auto edges = detail::PredDAGBuilder::ascending(*s.bundle);
        f["route"] = std::vector<EdgeId>(edges.begin(), edges.end());
      }
      f["members"] = s.plan ? py::cast(s.plan->members) : py::object(py::none());
      f["weights"] = s.plan && s.plan->placement == FlowPlacement::Weighted ? py::cast(s.plan->weights)
                                                                            : py::object(py::none());
      f["T"] = rec.delivered;
      f["delivered"] = rec.delivered;
      f["alloc"] = pairs_list(rec.edges);
      f["bundle_index"] = s.bundle_index < 0 ? py::object(py::none()) : py::object(py::int_(s.bundle_index));
      flows.append(f);
    }
    cohorts[py::int_(pe.roster[i])] = flows;
  }
  out["cohorts"] = cohorts;
  return out;
}

#endif

} // namespace

void register_flow_policy(py::module_& m) {
  [[maybe_unused]] auto flow_policy = py::class_<FlowPolicy>(
      m, "FlowPolicy",
      "The placement engine for one demand, bound to its ledger at construction. place(volume, outputs)\n"
      "is all-or-nothing: a failed call raises and leaves the ledger and the policy as they were (with\n"
      "the two limits in the module docstring). remove() takes every flow of the policy back out and\n"
      "starts a new epoch (the request total starts again from 0). A policy dropped without remove()\n"
      "leaves its records in the ledger as orphans: readable, never changed again, cleared by\n"
      "FlowLedger.reset(). sources and destinations are one NodeId or a sequence; the two sets must not\n"
      "overlap. ValueError for invalid endpoints, configurations and requests (non-integer node ids and\n"
      "requests included) and for objects of another graph, ledger or policy; TypeError for a ledger,\n"
      "config, elig, cache or outputs of the wrong class (a ledger of None: ValueError) and for an\n"
      "unknown keyword.")
      .def(py::init([](const std::shared_ptr<FlowLedger>& ledger, const py::object& sources,
                       const py::object& destinations, const FlowPolicyConfig& config,
                       const std::shared_ptr<Eligibility>& elig, const std::shared_ptr<DagCache>& cache) {
             if (!ledger) throw py::value_error("ledger must be a FlowLedger");
             const auto n = ledger->graph()->num_nodes();
             const auto s = node_ids(sources, n, "sources");
             const auto d = node_ids(destinations, n, "destinations");
             return std::make_unique<FlowPolicy>(*ledger, s, d, config, elig.get(), cache.get());
           }),
           py::keep_alive<1, 2>(), py::keep_alive<1, 6>(), py::keep_alive<1, 7>(), py::arg("ledger"),
           py::arg("sources"), py::arg("destinations"), py::arg("config"), py::arg("elig") = nullptr,
           py::arg("cache") = nullptr,
           "Bind a policy to ledger. elig (None: everything allowed) and cache (None: none; only CostOnly\n"
           "routing reads it, another model rejects it) belong to the ledger's graph; the policy keeps the\n"
           "ledger, elig and cache alive. The roster (the sources that statically reach an allowed\n"
           "destination, over allowed edges whatever their capacity) and, for Pinned, the usable bundles\n"
           "are fixed here. Everything is validated before the policy binds, in this order: the endpoints,\n"
           "elig's and cache's graph, the configuration (FlowPolicyConfig's order), Pinned's bundles.")
      .def(
          "place",
          [](FlowPolicy& p, const py::object& volume, const py::object& outputs) {
            const auto v = request_units(volume);
            const auto o = outputs_of(outputs);
            // All or nothing across the boundary: the result's Python object exists before the
            // call, and after the commit nothing allocates (the value moves into it). The GIL
            // stays held: the call mutates the ledger (one writer per ledger).
            auto [obj, slot] = new_object<PlacementResult>();
            *slot = p.place(v, o);
            return std::move(obj);
          },
          py::arg("volume"), py::arg("outputs") = py::none(),
          "Place volume (int units in [0, 2^62), or float('inf')) all or nothing -> PlacementResult.\n"
          "Requests add up within an epoch (from construction or the last remove()); after float('inf') any\n"
          "finite request is rejected, 0 included. A request of 0 places nothing and reports the current\n"
          "state. CostOnly prepares its routes at the epoch's first positive request; LOSSLESS admits the\n"
          "largest total in [admitted, admitted + volume] that the routes carry with no loss (a refused unit\n"
          "is not offered again); LOSSY forwards the new round-robin injections (earlier drops are not\n"
          "offered again). outputs (OutputOptions; None: none) requests edge_flows, residual, bottleneck and\n"
          "min_cut. Everything is validated before any change; a failed call raises with the ledger and the\n"
          "policy unchanged. A solve of the flow core (Global, CostOnly MAX_FLOW, or a MAX_FLOW flow) on a\n"
          "network beyond its int32 numbering raises ValueError: more than INT32_MAX - 2 nodes, or more than\n"
          "2^30 - 1 edges, sources and destinations together.")
      .def("remove", &FlowPolicy::remove, "revert every flow of the policy, retire the ids, open a new epoch")
      .def("placed_demand", &FlowPolicy::placed_demand, "the sum of delivered over the policy's live records")
      .def("records", &FlowPolicy::records, "the policy's live records as snapshots in FlowId order")
      .def(
          "record", [](const FlowPolicy& p, const py::object& flow_id) { return p.record(flow_id_of(flow_id)); },
          py::arg("flow_id"),
          "one of the policy's live records by FlowId (ValueError for a FlowId of another ledger or policy, a\n"
          "retired or an unknown id)")
      .def("usable_bundles", &FlowPolicy::usable_bundles,
           "Pinned: the bundles with a src -> dst walk over eligible edges, whatever their capacity (0 for\n"
           "the other modes)")
      .def_property_readonly("config", &FlowPolicy::config,
                             "the configuration as validated at construction, its absent origin and selection\n"
                             "resolved to the model's defaults");
#if NETGRAPH_CORE_PRIVATE_HOOKS
  flow_policy
      .def(
          "_place_traced",
          [](FlowPolicy& p, const py::object& volume, const py::object& outputs) {
            const auto v = request_units(volume);
            const auto o = outputs_of(outputs);
            // As place(): every Python object exists before anything changes; the engine runs the
            // `staged` callback, which converts the trace, after the work is done and before the
            // commit, so a failed conversion changes nothing.
            auto [obj, slot] = new_object<PlacementResult>();
            py::tuple pair(2);
            py::object trace_py;
            detail::Trace trace;
            const std::function<void()> staged = [&] { trace_py = trace_obj(trace); };
            *slot = detail::place_traced(p, v, o, trace, &staged);
            PyTuple_SET_ITEM(pair.ptr(), 0, obj.release().ptr());
            PyTuple_SET_ITEM(pair.ptr(), 1, trace_py.release().ptr());
            return pair;
          },
          py::arg("volume"), py::arg("outputs") = py::none(),
          "private: place() with the transaction's trace (the call's intermediate values, keyed as\n"
          "the tests' model of the policy keys them)")
      .def("_state", &state_obj, "private: the epoch state (roster, requested, +inf, cohorts)")
      .def(
          "_owner", [](const FlowPolicy& p) { return detail::PolicyAccess::engine(p).owner; },
          "private: the ledger-local owner id (FlowRecord.owner)");
#endif
}

} // namespace netgraph::core::bindings
