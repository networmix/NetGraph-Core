// capacity, batch_capacity and sensitivity with their options and results. Each evaluator
// places on a new, empty ledger; a frozen DagCache counts its hits only when the whole Python
// call succeeds (detail::frozen_call).
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/detail/evaluator_checks.hpp"
#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/graph.hpp"
#include "integers.hpp"
#include "publication.hpp"
#include "values.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// The Python form of QueryOptions: the Eligibility and the cache are
// held by reference count for the call (the C++ structs borrow raw pointers).
struct PyQueryOptions {
  FlowPolicyConfig config;
  std::shared_ptr<Eligibility> elig;
  std::shared_ptr<DagCache> cache;
  OutputOptions outputs;
  std::optional<std::optional<Amount>> design_demand;
  std::optional<std::optional<Amount>> offered;
  [[nodiscard]] QueryOptions native() const {
    return QueryOptions {config, elig.get(), cache.get(), outputs, design_demand, offered};
  }
};

// An optional amount option: None absent, else int units or float('inf') (the inner nullopt).
std::optional<std::optional<Amount>> amount_option(const py::object& obj, const char* text) {
  if (obj.is_none()) return std::nullopt;
  return units_or_inf(obj, text);
}

py::object amount_or_none(const std::optional<std::optional<Amount>>& x) {
  if (!x.has_value()) return py::none();
  if (!x->has_value()) return py::float_(std::numeric_limits<double>::infinity());
  return py::int_(**x);
}

// workers: an explicit integer is passed through (validated by the evaluator); None takes the
// environment default (detail::default_workers), and 1 whenever a cache is supplied (a building
// cache has one writer; a frozen one is shared only when the caller asks for workers).
std::int32_t workers_of(const py::object& workers, bool sensitivity, bool cache_given) {
  if (!workers.is_none()) {
    return static_cast<std::int32_t>(integer_value(workers, std::numeric_limits<std::int32_t>::min(),
                                                   std::numeric_limits<std::int32_t>::max(), "workers"));
  }
  if (cache_given) return 1;
  return detail::default_workers(sensitivity);
}

} // namespace

void register_evaluators(py::module_& m) {
  py::class_<PyQueryOptions>(m, "QueryOptions",
                             "The options of every evaluator (capacity, batch_capacity, sensitivity): the configuration,\n"
                             "the static Eligibility (None: all allowed), the DagCache (None: private per worker), the\n"
                             "optional outputs of each result (sensitivity: of its baseline), design_demand (int units or\n"
                             "float('inf'); None: absent): T1's design, always required and finite (`T1 capacity requires\n"
                             "design_demand`, `T1 requires a finite design demand`), rejected by every other model\n"
                             "(`<model> does not accept design_demand`), and offered (int units or float('inf'); None:\n"
                             "absent), the offer of a T1 or Pinned LOSSY evaluation: required and finite there (`<T1 or\n"
                             "Pinned> LOSSY capacity requires finite offered`); a LOSSLESS model rejects it (`<model>\n"
                             "LOSSLESS capacity does not accept offered`), CostOnly LOSSY too (`CostOnly LOSSY capacity\n"
                             "uses first-hop injection`).")
      .def(py::init([](const FlowPolicyConfig& config, const std::shared_ptr<Eligibility>& elig,
                       const std::shared_ptr<DagCache>& cache, const py::object& outputs, const py::object& design_demand,
                       const py::object& offered) {
             auto d = amount_option(design_demand, "design_demand must be an integer in [0, 2^62) or permitted +inf");
             auto o = amount_option(offered, "offered must be an integer in [0, 2^62) or permitted +inf");
             return PyQueryOptions {config, elig, cache, outputs_of(outputs), d, o};
           }),
           py::arg("config"), py::kw_only(), py::arg("elig") = nullptr, py::arg("cache") = nullptr,
           py::arg("outputs") = py::none(), py::arg("design_demand") = py::none(), py::arg("offered") = py::none())
      .def_readonly("config", &PyQueryOptions::config,
                    "the FlowPolicyConfig of every evaluated policy (sensitivity: the baseline and every deletion)")
      .def_readonly("elig", &PyQueryOptions::elig, "the static Eligibility; None: everything allowed")
      .def_readonly("cache", &PyQueryOptions::cache, "the supplied DagCache; None: a private cache per worker")
      .def_readonly("outputs", &PyQueryOptions::outputs,
                    "the OutputOptions of every result (sensitivity: of the baseline)")
      .def_property_readonly(
          "design_demand", [](const PyQueryOptions& o) { return amount_or_none(o.design_demand); },
          "T1's design volume (int units, or float('inf')); None: absent")
      .def_property_readonly(
          "offered", [](const PyQueryOptions& o) { return amount_or_none(o.offered); },
          "the offer of a T1 or Pinned LOSSY evaluation (int units, or float('inf')); None: absent");

  py::class_<SensitivityEntry>(m, "SensitivityEntry",
                                "One candidate edge: edge and loss (int units: the baseline minus the capacity without the\n"
                                "edge).")
      .def_readonly("edge", &SensitivityEntry::edge)
      .def_readonly("loss", &SensitivityEntry::loss);

  py::class_<SensitivityResult>(m, "SensitivityResult",
                                "baseline (the PlacementResult of capacity()) and entries (ascending EdgeId; an\n"
                                "absent edge has loss 0).")
      .def_readonly("baseline", &SensitivityResult::baseline)
      .def_readonly("entries", &SensitivityResult::entries);

  m.def(
      "capacity",
      [](const std::shared_ptr<Graph>& graph, const py::object& sources, const py::object& destinations,
         const py::object& opts) {
        if (!graph) throw py::value_error("graph must be a Graph");
        const auto o = options_of<PyQueryOptions>(opts, "QueryOptions");
        const auto s = node_ids(sources, graph->num_nodes(), "sources");
        const auto d = node_ids(destinations, graph->num_nodes(), "destinations");
        const QueryOptions c = o.native();
        // One public call on a frozen cache's statistics, its result object included: its hits count
        // only when the Python result exists.
        return detail::frozen_call(c.cache, [&](detail::FrozenTally& tally) {
          PlacementResult r;
          {
            py::gil_scoped_release release;
            r = detail::capacity_in_call(*graph, s, d, c, tally);
          }
          return object_of(std::move(r));
        });
      },
      py::arg("graph"), py::arg("sources"), py::arg("destinations"), py::arg("opts"),
      "place(+inf) of a FlowPolicy with opts.config on a fresh ledger of the graph (Pinned LOSSY:\n"
      "place(opts.offered); T1: evaluate(+inf or opts.offered, design_demand=opts.design_demand), the\n"
      "lossless bound of the tunnels set up for the design, or what they deliver of offered) -> PlacementResult.\n"
      "The exact maximum flow is Global() (the optimized origin); with origin=EQUAL, the largest total whose\n"
      "round-robin source shares are all carried. A supplied frozen cache raises DagCacheMiss for an absent\n"
      "key, unchanged.\n"
      "A failed call raises; a supplied DagCache then behaves as the DagCache docstring states.");

  m.def(
      "batch_capacity",
      [](const std::shared_ptr<Graph>& graph, const py::object& pairs, const py::object& opts,
         const py::object& per_pair_elig, const py::object& workers) {
        if (!graph) throw py::value_error("graph must be a Graph");
        // The arguments first, then each pair completely and in pair order (the same two parts
        // as the C++ validation), so the lowest failing pair index is reported.
        const auto o = options_of<PyQueryOptions>(opts, "QueryOptions");
        const auto w = workers_of(workers, false, static_cast<bool>(o.cache));
        const QueryOptions c = o.native();
        detail::check_batch_arguments(*graph, c, w);
        if (not_a_sequence(pairs)) throw py::value_error("pairs must be a sequence of (sources, destinations)");
        const auto items = py::reinterpret_borrow<py::sequence>(pairs);
        const auto num_pairs = py::len(items);
        std::optional<py::sequence> eligs;
        if (!per_pair_elig.is_none()) {
          // An explicit sequence has one entry per pair (None entries allowed), also when empty.
          if (not_a_sequence(per_pair_elig)) throw py::value_error("per_pair_elig must be a sequence of Eligibility or None");
          eligs = py::reinterpret_borrow<py::sequence>(per_pair_elig);
          if (py::len(*eligs) != num_pairs) {
            throw py::value_error("per_pair_elig has " + std::to_string(py::len(*eligs)) + " entries for " +
                                  std::to_string(num_pairs) + " pairs");
          }
        }
        std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> ps;
        std::vector<std::shared_ptr<Eligibility>> keep;
        std::vector<const Eligibility*> per;
        const auto n = graph->num_nodes();
        for (std::size_t i = 0; i < num_pairs; ++i) {
          const py::object item = items[i];
          if (not_a_sequence(item) || py::len(item) != 2) {
            throw py::value_error("each pair must be (sources, destinations)");
          }
          auto t = py::reinterpret_borrow<py::sequence>(item);
          auto src = node_ids(t[0], n, "sources");
          auto dst = node_ids(t[1], n, "destinations");
          const Eligibility* e = nullptr;
          if (eligs) {
            const py::object entry = (*eligs)[i];
            if (!entry.is_none()) {
              if (!py::isinstance<Eligibility>(entry)) throw py::type_error("per_pair_elig entries must be Eligibility or None");
              keep.push_back(entry.cast<std::shared_ptr<Eligibility>>());
              e = keep.back().get();
            }
            per.push_back(e);
          }
          detail::check_batch_pair(*graph, c.config, src, dst, e, i);
          ps.emplace_back(std::move(src), std::move(dst));
        }
        return detail::frozen_call(c.cache, [&](detail::FrozenTally& tally) {
          std::vector<PlacementResult> rs;
          {
            py::gil_scoped_release release;
            rs = detail::batch_capacity_in_call(*graph, ps, c, per, w, tally);
          }
          py::list out(rs.size());
          for (std::size_t i = 0; i < rs.size(); ++i) {
            PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), object_of(std::move(rs[i])).release().ptr());
          }
          return py::object(std::move(out));
        });
      },
      py::arg("graph"), py::arg("pairs"), py::arg("opts"), py::arg("per_pair_elig") = py::none(),
      py::arg("workers") = py::none(),
      "capacity() for each (sources, destinations) pair, results in pair order. Pair i uses opts.elig\n"
      "intersected with per_pair_elig[i] (None entries allowed; a given sequence has one entry per pair,\n"
      "also when empty). workers >= 1; None: NGRAPH_CORE_BATCH_THREADS or 1, and 1 when a cache is\n"
      "supplied. Everything is validated before a worker starts, the arguments first, then each pair\n"
      "completely in pair order, so the error of the lowest failing pair is raised. Pinned is rejected\n"
      "(`batch_capacity does not accept Pinned`: a program has fixed endpoints);\n"
      "an empty batch returns [].\n"
      "Caches: a supplied building cache needs workers == 1; a frozen cache is shared by every worker\n"
      "(pass workers explicitly), and its DagCacheMiss reaches the caller unchanged (the lowest failing\n"
      "pair's). An intersection that differs from both operands gets a new Eligibility uid, which a\n"
      "frozen cache cannot hold: pass the combined per-pair values with elig None, or a common value\n"
      "that contains them.\n"
      "A failed call raises; a supplied DagCache then behaves as the DagCache docstring states.");

  m.def(
      "sensitivity",
      [](const std::shared_ptr<Graph>& graph, const py::object& sources, const py::object& destinations,
         const py::object& opts, const py::object& workers) {
        if (!graph) throw py::value_error("graph must be a Graph");
        const auto o = options_of<PyQueryOptions>(opts, "QueryOptions");
        const auto s = node_ids(sources, graph->num_nodes(), "sources");
        const auto d = node_ids(destinations, graph->num_nodes(), "destinations");
        const auto w = workers_of(workers, true, static_cast<bool>(o.cache));
        const QueryOptions c = o.native();
        return detail::frozen_call(c.cache, [&](detail::FrozenTally& tally) {
          SensitivityResult r;
          {
            py::gil_scoped_release release;
            r = detail::sensitivity_in_call(*graph, s, d, c, w, tally);
          }
          return object_of(std::move(r));
        });
      },
      py::arg("graph"), py::arg("sources"), py::arg("destinations"), py::arg("opts"), py::arg("workers") = py::none(),
      "The baseline capacity() and one SensitivityEntry per candidate edge. Candidates: with unrestricted\n"
      "Global and the optimized origin, the edges that carry flow in the baseline (removing any other edge cannot\n"
      "lower a maximum flow); otherwise every allowed edge, whatever its capacity (an edge of capacity 0\n"
      "can carry a route, so deleting it can change the result, and with the equal origin a deletion that\n"
      "removes a source from the roster can raise the total: a negative loss). loss = baseline.placed\n"
      "minus the placed of capacity() without the edge, in signed units, never thresholded. Two\n"
      "shortcuts give exactly that value: unrestricted optimized Global repairs a copy of the baseline's\n"
      "solve instead of solving again, and CostOnly\n"
      "gives loss 0 without placing to an edge outside the static shortest-path support reachable from\n"
      "the roster sources (no route can use it). Global with a support solves every deletion afresh; Pinned\n"
      "replays the same program under each deletion. workers as batch_capacity\n"
      "(NGRAPH_CORE_SENSITIVITY_THREADS). Only CostOnly reads the cache, and each placed deletion uses a\n"
      "new Eligibility uid, so a frozen cache raises DagCacheMiss there.\n"
      "A failed deletion fails the call (the lowest failing candidate's exception, no partial result); a\n"
      "supplied DagCache then behaves as the DagCache docstring states.");
}

} // namespace netgraph::core::bindings
