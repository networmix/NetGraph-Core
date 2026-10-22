// evaluate_demand, evaluate_pairs and evaluate_edge_failures with their options and results. Each
// evaluator places on a new, empty ledger; a frozen ShortestPathCache counts its hits only when the whole Python
// call succeeds (detail::frozen_call).
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/shortest_path_cache.hpp"
#include "netgraph/core/detail/evaluator_checks.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/graph.hpp"
#include "integers.hpp"
#include "publication.hpp"
#include "native_call.hpp"
#include "values.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// The Python form of EvaluationOptions: the Eligibility and the cache are held by reference count
// for the call (the C++ structs borrow raw pointers); design_demand and offer as given.
struct PyEvaluationOptions {
  FlowPolicyConfig config;
  std::shared_ptr<Eligibility> eligibility;
  std::shared_ptr<ShortestPathCache> cache;
  OutputOptions outputs;
  py::object design_demand;
  py::object offer;
  [[nodiscard]] EvaluationOptions native() const {
    return EvaluationOptions {config, eligibility.get(), cache.get(), outputs, units_at_use(design_demand), units_at_use(offer)};
  }
};

// workers: an explicit integer is passed through (validated by the evaluator); None takes the
// environment default (detail::default_workers), and 1 whenever a cache is supplied (a building
// cache has one writer; a frozen one is shared only when the caller asks for workers).
std::int32_t workers_of(const py::object& workers, bool edge_failures, bool cache_given) {
  if (!workers.is_none()) {
    return static_cast<std::int32_t>(integer_value(workers, std::numeric_limits<std::int32_t>::min(),
                                                   std::numeric_limits<std::int32_t>::max(), "workers"));
  }
  if (cache_given) return 1;
  return detail::default_workers(edge_failures);
}

} // namespace

void register_evaluators(py::module_& m) {
  py::class_<PyEvaluationOptions>(m, "EvaluationOptions",
                             "The options of every evaluator (evaluate_demand, evaluate_pairs, evaluate_edge_failures): config\n"
                             "(the configuration), eligibility (the static Eligibility; None: all allowed), cache (the\n"
                             "ShortestPathCache; None: private per worker), outputs (the optional outputs of each result;\n"
                             "evaluate_edge_failures: of its baseline), design_demand and\n"
                             "offer (int units or a float, kept as given; another type raises TypeError). offer: a finite offer\n"
                             "for every model (None: the uncapped LOSSLESS question, or for ShortestPaths EQUAL/WEIGHTED with\n"
                             "LOSSY the first-hop experiment; required for DistributedTE, Pinned and, without a finite design,\n"
                             "CentralizedTE with tunnel_count under LOSSY). design_demand: DistributedTE always, required and\n"
                             "finite; CentralizedTE with tunnel_count: None means the offer when it is finite, else unlimited;\n"
                             "every other model leaves it None. The evaluators check both after the configuration, in this order:\n"
                             "an explicit offer (`offer must be an integer in [0, 2^62)`, +inf included), the design\n"
                             "(`<model> does not accept design_demand`, `design_demand must be an integer in [0, 2^62) or\n"
                             "permitted +inf`, `DistributedTE requires a finite design demand`), DistributedTE's missing design\n"
                             "(`DistributedTE evaluation requires design_demand`), the missing LOSSY offer (`<DistributedTE or\n"
                             "CentralizedTE or Pinned> LOSSY evaluation requires a finite offer`).")
      .def(py::init([](const FlowPolicyConfig& config, const std::shared_ptr<Eligibility>& elig,
                       const std::shared_ptr<ShortestPathCache>& cache, const py::object& outputs, const py::object& design_demand,
                       const py::object& offer) {
             auto d = units_as_given(design_demand, "design_demand");
             auto o = units_as_given(offer, "offer");
             return PyEvaluationOptions {config, elig, cache, outputs_of(outputs), std::move(d), std::move(o)};
           }),
           py::arg("config"), py::kw_only(), py::arg("eligibility") = nullptr, py::arg("cache") = nullptr,
           py::arg("outputs") = py::none(), py::arg("design_demand") = py::none(), py::arg("offer") = py::none())
      .def_property_readonly("config", field_object(&PyEvaluationOptions::config),
                    "the FlowPolicyConfig of every evaluated policy (evaluate_edge_failures: the baseline and every deletion)")
      .def_property_readonly("eligibility", shared_field(&PyEvaluationOptions::eligibility), "the static Eligibility; None: everything allowed")
      .def_property_readonly("cache", shared_field(&PyEvaluationOptions::cache), "the supplied ShortestPathCache; None: a private cache per worker")
      .def_property_readonly("outputs", field_object(&PyEvaluationOptions::outputs),
                    "the OutputOptions of every result (evaluate_edge_failures: of the baseline)")
      .def_readonly("design_demand", &PyEvaluationOptions::design_demand,
                    "the tunnel models' design volume as given (int units or a float); None: absent")
      .def_readonly("offer", &PyEvaluationOptions::offer, "the offer as given (int units or a float); None: absent");

  py::class_<SensitivityEntry>(m, "SensitivityEntry",
                                "One candidate edge: edge and delivery_reduction (int units: the baseline's delivery minus the\n"
                                "delivery without the edge; signed).")
      .def_property_readonly("edge", int_reader(&SensitivityEntry::edge))
      .def_property_readonly("delivery_reduction", int_reader(&SensitivityEntry::delivery_reduction));

  py::class_<SensitivityResult>(m, "SensitivityResult",
                                "baseline (the PlacementResult of evaluate_demand()) and entries (ascending EdgeId; an\n"
                                "absent edge has reduction 0).")
      .def_property_readonly("baseline", field_object(&SensitivityResult::baseline))
      .def_property_readonly("entries", list_field(&SensitivityResult::entries));

  m.def(
      "evaluate_demand",
      [](const std::shared_ptr<Graph>& graph, const py::object& sources, const py::object& destinations,
         const py::object& opts) {
        if (!graph) throw py::type_error("graph must be a Graph");
        const auto o = options_of<PyEvaluationOptions>(opts, "EvaluationOptions");
        const auto s = node_ids(sources, graph->num_nodes(), "sources");
        const auto d = node_ids(destinations, graph->num_nodes(), "destinations");
        const EvaluationOptions c = o.native();
        const auto roster = detail::check_demand_arguments(*graph, s, d, c, detached_argument_work);
        auto [obj, slot] = new_object<PlacementResult>();
        native_call([&] {
          return detail::frozen_call(c.cache, [&](detail::FrozenCacheHitTally& tally) {
            *slot = detail::evaluate_demand_in_call(*graph, s, d, c, tally, roster ? &*roster : nullptr);
            return true;
          });
        });
        return std::move(obj);
      },
      py::arg("graph"), py::arg("sources"), py::arg("destinations"), py::arg("options"),
      "The evaluation of a fresh FlowPolicy with options.config on a fresh ledger of graph at options.offer\n"
      "(None: +inf; DistributedTE and CentralizedTE with tunnel_count: evaluate(offer, design_demand=...))\n"
      "-> PlacementResult. LOSSLESS delivers up to a finite offer; LOSSY accounts for exactly that offer\n"
      "(ShortestPaths EQUAL/WEIGHTED LOSSY without an offer: the first-hop experiment). The options are checked as\n"
      "the EvaluationOptions docstring states.\n"
      "The exact maximum flow is CentralizedTE(source_allocation=OPTIMIZED); with source_allocation=EQUAL, the largest\n"
      "total whose round-robin source shares are all carried; with WEIGHTED (CentralizedTE without tunnel_count),\n"
      "the largest total whose highest-averages shares of source_weights are all carried.\n"
      "A supplied frozen cache raises ShortestPathCacheMiss\n"
      "for an absent key, unchanged.\n"
      "A failed call raises; a supplied ShortestPathCache then behaves as the ShortestPathCache docstring states.");

  m.def(
      "evaluate_pairs",
      [](const std::shared_ptr<Graph>& graph, const py::object& pairs, const py::object& opts,
         const py::object& per_pair_eligibility, const py::object& workers) {
        if (!graph) throw py::type_error("graph must be a Graph");
        // The arguments first, then each pair completely and in pair order (the same two parts
        // as the C++ validation), so the lowest failing pair index is reported.
        const auto o = options_of<PyEvaluationOptions>(opts, "EvaluationOptions");
        const auto w = workers_of(workers, false, static_cast<bool>(o.cache));
        const EvaluationOptions c = o.native();
        detail::check_batch_arguments(*graph, c, w);
        if (is_invalid_sequence(pairs)) throw py::type_error("pairs must be a sequence of (sources, destinations)");
        const auto items = py::reinterpret_borrow<py::sequence>(pairs);
        const auto num_pairs = py::len(items);
        std::optional<py::sequence> eligs;
        if (!per_pair_eligibility.is_none()) {
          // An explicit sequence has one entry per pair (None entries allowed), also when empty.
          if (is_invalid_sequence(per_pair_eligibility)) throw py::type_error("per_pair_eligibility must be a sequence of Eligibility or None");
          eligs = py::reinterpret_borrow<py::sequence>(per_pair_eligibility);
          if (py::len(*eligs) != num_pairs) {
            throw py::value_error("per_pair_eligibility has " + std::to_string(py::len(*eligs)) + " entries for " +
                                  std::to_string(num_pairs) + " pairs");
          }
        }
        std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> ps;
        std::vector<std::shared_ptr<Eligibility>> keep;
        std::vector<const Eligibility*> per;
        std::exception_ptr pending;
        const auto n = graph->num_nodes();
        for (std::size_t i = 0; i < num_pairs; ++i) {
          try {
            const py::object item = items[i];
            // Not a sequence is a wrong type, a sequence of another length a wrong value: one text.
            if (is_invalid_sequence(item)) throw py::type_error("each pair must be (sources, destinations)");
            if (py::len(item) != 2) throw py::value_error("each pair must be (sources, destinations)");
            auto t = py::reinterpret_borrow<py::sequence>(item);
            auto src = node_ids(t[0], n, "sources");
            auto dst = node_ids(t[1], n, "destinations");
            const Eligibility* e = nullptr;
            if (eligs) {
              const py::object entry = (*eligs)[i];
              if (!entry.is_none()) {
                if (!py::isinstance<Eligibility>(entry)) throw py::type_error("per_pair_eligibility entries must be Eligibility or None");
                keep.push_back(entry.cast<std::shared_ptr<Eligibility>>());
                e = keep.back().get();
              }
              per.push_back(e);
            }
            detail::check_batch_pair(*graph, c.config, src, dst, e, i);
            ps.emplace_back(std::move(src), std::move(dst));
          } catch (...) {
            pending = std::current_exception();
            break;
          }
        }
        const auto rosters = detail::prepare_batch_rosters(*graph, ps, c, per, w, detached_argument_work);
        if (pending) std::rethrow_exception(pending);
        auto out = py::reinterpret_steal<py::list>(checked(PyList_New(static_cast<Py_ssize_t>(num_pairs))).release());
        std::vector<PlacementResult*> slots(num_pairs);
        for (std::size_t i = 0; i < num_pairs; ++i) {
          auto [obj, slot] = new_object<PlacementResult>();
          slots[i] = slot;
          PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), obj.release().ptr());
        }
        native_call([&] {
          return detail::frozen_call(c.cache, [&](detail::FrozenCacheHitTally& tally) {
            auto rs = detail::evaluate_pairs_in_call(*graph, ps, c, per, w, tally, rosters);
            for (std::size_t i = 0; i < rs.size(); ++i) *slots[i] = std::move(rs[i]);
            return true;
          });
        });
        return out;
      },
      py::arg("graph"), py::arg("pairs"), py::arg("options"), py::arg("per_pair_eligibility") = py::none(),
      py::arg("workers") = py::none(),
      "evaluate_demand() for each (sources, destinations) pair, results in pair order. Pair i uses\n"
      "options.eligibility "
      "intersected with per_pair_eligibility[i] (None entries allowed; a given sequence has one entry per pair,\n"
      "also when empty). workers >= 1; None: NGRAPH_CORE_BATCH_THREADS or 1, and 1 when a cache is\n"
      "supplied. All checks finish before evaluation workers start: the arguments first, then each pair\n"
      "completely in pair order, so the error of the lowest failing pair is raised. First-pass tunnel\n"
      "searches run in parallel under one GIL release, before any lease or gate. Pinned is rejected\n"
      "(`evaluate_pairs does not accept Pinned`: a program has fixed endpoints);\n"
      "an empty batch returns [].\n"
      "Caches: a supplied building cache needs workers == 1; a frozen cache is shared by every worker\n"
      "(pass workers explicitly), and its ShortestPathCacheMiss reaches the caller unchanged (the lowest failing\n"
      "pair's). An intersection that differs from both operands gets a new Eligibility uid, which a\n"
      "frozen cache cannot hold: pass the combined per-pair values with eligibility None, or a common value\n"
      "that contains them.\n"
      "A failed call raises; a supplied ShortestPathCache then behaves as the ShortestPathCache docstring states.");

  m.def(
      "evaluate_edge_failures",
      [](const std::shared_ptr<Graph>& graph, const py::object& sources, const py::object& destinations,
         const py::object& opts, const py::object& workers) {
        if (!graph) throw py::type_error("graph must be a Graph");
        const auto o = options_of<PyEvaluationOptions>(opts, "EvaluationOptions");
        const auto s = node_ids(sources, graph->num_nodes(), "sources");
        const auto d = node_ids(destinations, graph->num_nodes(), "destinations");
        const auto w = workers_of(workers, true, static_cast<bool>(o.cache));
        const EvaluationOptions c = o.native();
        const auto roster = detail::check_edge_failure_arguments(*graph, s, d, c, w, detached_argument_work);
        auto [obj, slot] = new_object<SensitivityResult>();
        native_call([&] {
          return detail::frozen_call(c.cache, [&](detail::FrozenCacheHitTally& tally) {
            *slot = detail::evaluate_edge_failures_in_call(*graph, s, d, c, w, tally, roster ? &*roster : nullptr);
            return true;
          });
        });
        return std::move(obj);
      },
      py::arg("graph"), py::arg("sources"), py::arg("destinations"), py::arg("options"), py::arg("workers") = py::none(),
      "The baseline evaluate_demand() and one SensitivityEntry per candidate edge. Candidates: with unrestricted\n"
      "CentralizedTE and the optimized source allocation, the edges that carry flow in the baseline (removing any\n"
      "other edge keeps the baseline feasible, at any offer); otherwise every allowed edge, whatever its capacity (an edge of capacity 0\n"
      "can carry a route, so deleting it can change the result, and with the equal source allocation a deletion that\n"
      "removes a source from the roster can raise the total: a negative reduction). delivery_reduction =\n"
      "baseline.call_delivered minus the call_delivered of evaluate_demand() without the edge, in signed\n"
      "units, never thresholded. Two shortcuts give exactly that value. Unrestricted optimized CentralizedTE\n"
      "whose baseline is a maximum flow (no offer, or an offer above what it delivers) reads each loss from\n"
      "the baseline's residual network without a new solve: a saturated edge whose endpoints lie in different\n"
      "strongly connected components loses its whole flow, any other its flow minus what can be rerouted\n"
      "between its endpoints. At an offer the baseline delivers, each deletion is solved at the offer. ShortestPaths\n"
      "gives reduction 0 without placing to an edge outside the static shortest-path support reachable from\n"
      "the roster sources (no route can use it). CentralizedTE with a support solves every deletion afresh; Pinned\n"
      "replays the same program under each deletion. workers as evaluate_pairs\n"
      "(NGRAPH_CORE_SENSITIVITY_THREADS). Only ShortestPaths reads the cache, and each placed deletion uses a\n"
      "new Eligibility uid, so a frozen cache raises ShortestPathCacheMiss there.\n"
      "A failed deletion fails the call (the lowest failing candidate's exception, no partial result); a\n"
      "supplied ShortestPathCache then behaves as the ShortestPathCache docstring states.");
}

} // namespace netgraph::core::bindings
