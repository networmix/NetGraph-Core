/*
  Pybind11 module `_netgraph_core`: the graph on integer units with Eligibility and the unit
  conversions, the route structures (PredDAG, Path), the selection and placement axes, the
  ledger with its records, the policy engine with its configuration and results, the DagCache,
  the evaluators (capacity, batch_capacity, sensitivity) and the public search (spf, ksp,
  resolve_to_paths); and one private entry: the flow core's _place_global_units.
  A build with the private hooks (CMake NETGRAPH_CORE_PRIVATE_HOOKS: the editable development
  build) adds the test controls and the parity diagnostics.

  - Graph accessors return read-only NumPy views that keep the graph alive (no per-call copy).
  - Integer inputs (capacities, costs, ids, views, caps) must be integers: Python ints, NumPy
    integer scalars or integer arrays; floats and bools are rejected with ValueError, never
    truncated. Values are read without narrowing (WideInts) and validated in the C++
    constructors' order before any conversion to the C++ types. Physical values
    enter only through to_units; the resolution is validated before anything else.
  - std::invalid_argument from Core maps to ValueError; enumeration values outside their
    enumeration are rejected at construction (ValueError).
  - _place_global_units is a private entry point for the flow checks: it returns the flow
    core's payload as a dict of Python ints and lists.
  - Amounts and counters cross as Python ints (exact, also beyond 64 bits); per-edge arrays as
    int64 NumPy arrays in EdgeId order; +inf is accepted only as a request (float('inf')).
  - Python spells constants in UPPER_SNAKE_CASE: MAX_FLOWS_PER_CALL is the C++ kMaxFlowsPerCall.
*/
#include <pybind11/pybind11.h>

#include "registration.hpp"

namespace py = pybind11;
using namespace netgraph::core;

PYBIND11_MODULE(_netgraph_core, m, py::mod_gil_not_used()) {
  m.doc() =
      "NetGraph-Core: the graph on integer units, the route structures, the integer ledger, the policy\n"
      "engine FlowPolicy, the DagCache, the evaluators (capacity, batch_capacity, sensitivity) and the\n"
      "searches (spf, ksp, resolve_to_paths).\n"
      "Every public call either completes or raises and leaves Core state as it was, with two limits:\n"
      "(1) pybind11 does not check every allocation when it constructs a bound object from Python, or\n"
      "converts an attribute or a returned value into a new Python object outside the library's checked\n"
      "result construction; running out of memory exactly there may crash the process, with Core state\n"
      "unchanged. (2) A DagCache passed to the call: a building one keeps the complete entries and\n"
      "statistics of the lookups that succeeded before the call failed; a frozen one counts one miss\n"
      "when the call raised DagCacheMiss.";

  // In this order: a later binding's defaults and conversions need the earlier types.
  bindings::register_axes(m);
  bindings::register_units(m);
  bindings::register_eligibility(m);
  bindings::register_graph(m);
  bindings::register_recorder(m);
  bindings::register_pred_dag(m);
  bindings::register_path(m);
  bindings::register_flow_record(m);
  bindings::register_flow_ledger(m);
  bindings::register_routing(m);
  bindings::register_policy_config(m);
  bindings::register_placement_result(m);
  bindings::register_dag_cache(m);
  bindings::register_flow_policy(m);
  bindings::register_evaluators(m);
  bindings::register_paths(m);
  bindings::register_search(m);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  bindings::register_private_hooks(m);
#endif
}
