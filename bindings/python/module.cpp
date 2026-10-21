/*
  Pybind11 module `_netgraph_core`: the graph on integer units with Eligibility and the unit
  conversions, the route structures (PredecessorDAG, Path), the route choice and placement axes, the
  ledger with its records, the forwarding programs, the policy engine with its configuration
  and results, several policies' actions in one operation (apply_actions), the ShortestPathCache,
  the evaluators (evaluate_demand, evaluate_pairs, evaluate_edge_failures) and the public search (spf, ksp,
  resolve_to_paths); and one private entry: the flow core's _solve_min_cost_flow.
  A build with the private hooks (CMake NETGRAPH_CORE_PRIVATE_HOOKS: the editable development
  build) adds the test controls and the parity diagnostics.

  - Graph accessors return read-only NumPy views that keep the graph alive (no per-call copy).
  - Integer inputs (capacities, costs, ids, views, caps) take any non-bool object with __index__, including Python ints and NumPy
    integer scalars, and sequences or integer arrays; floats and bools are rejected with TypeError, never
    truncated. Values are read without narrowing and validated in the C++
    constructors' order before any conversion to the C++ types. Physical values
    enter only through to_units; the resolution is validated before anything else.
  - One rule for the class of an input error: a type the position never takes TypeError, a value
    outside its domain ValueError. std::invalid_argument from Core maps to ValueError; enumeration
    domains are checked when a configuration or options value is consumed (ValueError).
  - _solve_min_cost_flow is a private entry point for the flow checks: it returns the flow
    core's payload as a dict of Python ints and lists.
  - Amounts and counters cross as Python ints (exact, also beyond 64 bits); per-edge arrays as
    int64 NumPy arrays in EdgeId order; +inf is accepted only as a request (float('inf')).
  - Python spells constants in UPPER_SNAKE_CASE: UNREACHABLE is the C++ kUnreachable.
*/
#include <pybind11/pybind11.h>

#include "registration.hpp"

namespace py = pybind11;
using namespace netgraph::core;

PYBIND11_MODULE(_netgraph_core, m, py::mod_gil_not_used()) {
  m.doc() =
      "NetGraph-Core: the graph on integer units, the route structures, the integer ledger, the forwarding\n"
      "programs, the policy engine FlowPolicy and its lists of actions (apply_actions), the ShortestPathCache,\n"
      "the evaluators (evaluate_demand, evaluate_pairs, evaluate_edge_failures) and the searches (spf, ksp,\n"
      "resolve_to_paths).\n"
      "One rule for the class of an input error, at every public entry: a type the position never takes\n"
      "(a bool at a numeric position, a float where only an integer is taken, a string, None where none\n"
      "is allowed, another class; for a vector a non-sequence, a nested entry, an array of another rank\n"
      "or of a non-integer dtype) raises TypeError, a value outside its domain (a range, a length, NaN,\n"
      "-inf, a finite float where floats are taken only for +inf) ValueError; the text is the same.\n"
      "Every public call either completes or raises and leaves Core state as it was, with two limits:\n"
      "(1) pybind11 does not check every allocation when it constructs a bound object from Python, parses\n"
      "the arguments of a call that passes keywords, or converts an attribute or a returned value into a\n"
      "new Python object outside the library's checked result construction; running out of memory exactly\n"
      "there may crash the process, with Core state unchanged. (2) A ShortestPathCache passed to the call: a building one keeps the complete entries and\n"
      "statistics of the lookups that succeeded before the call failed; a frozen one counts one miss\n"
      "when the call raised ShortestPathCacheMiss.";

  // In this order: a later binding's defaults and conversions need the earlier types.
  bindings::register_options(m);
  bindings::register_units(m);
  bindings::register_eligibility(m);
  bindings::register_graph(m);
  bindings::register_flow_solver_hook(m);
  bindings::register_predecessor_dag(m);
  bindings::register_path(m);
  bindings::register_flow_record(m);
  bindings::register_flow_ledger(m);
  bindings::register_program(m);
  bindings::register_routing(m);
  bindings::register_policy_config(m);
  bindings::register_placement_result(m);
  bindings::register_shortest_path_cache(m);
  bindings::register_flow_policy(m);
  bindings::register_apply_actions(m);
  bindings::register_evaluators(m);
  bindings::register_paths(m);
  bindings::register_search(m);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  bindings::register_private_hooks(m);
#endif
}
