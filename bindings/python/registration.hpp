/*
  The registration functions of the module _netgraph_core, one file each under bindings/python/.
  The module entry (module.cpp) calls them in this order, which is the order of the module's
  attributes: a default argument or a conversion of a later binding needs the types registered
  before it. The private recorder entry (_solve_min_cost_flow) is registered
  in both builds; only register_private_hooks and the private members of FlowLedger and
  FlowPolicy depend on NETGRAPH_CORE_PRIVATE_HOOKS.
*/
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "netgraph/core/detail/hooks_flag.hpp"

namespace netgraph::core::bindings {

void register_options(pybind11::module_& m);              // the axes' enumerations
void register_units(pybind11::module_& m);             // to_units, as_float
void register_eligibility(pybind11::module_& m);       // Eligibility
void register_graph(pybind11::module_& m);             // Graph
void register_flow_solver_hook(pybind11::module_& m);          // _solve_min_cost_flow (both builds)
void register_predecessor_dag(pybind11::module_& m);          // PredecessorDAG
void register_path(pybind11::module_& m);              // Path
void register_flow_record(pybind11::module_& m);       // FlowRecord
void register_flow_ledger(pybind11::module_& m);       // FlowLedger
void register_program(pybind11::module_& m);           // TunnelStatus, Tunnel, RouteDefinition, ForwardingProgram
void register_routing(pybind11::module_& m);           // AbsoluteCost, CostStretch, TunnelCount, the four routings
void register_policy_config(pybind11::module_& m);     // FlowPolicyConfig, OutputOptions
void register_placement_result(pybind11::module_& m);  // Termination, SourceFlow, TunnelFlow, PlacementResult
void register_shortest_path_cache(pybind11::module_& m);         // ShortestPathCacheStats, ShortestPathCacheMiss, ShortestPathCache
void register_flow_policy(pybind11::module_& m);       // FlowPolicy
void register_apply_actions(pybind11::module_& m);     // BatchMode, Place, Replace, Forward, apply_actions
void register_evaluators(pybind11::module_& m);        // evaluate_demand, evaluate_pairs, evaluate_edge_failures and their types
void register_paths(pybind11::module_& m);             // KspOptions, ksp, resolve_to_paths
void register_search(pybind11::module_& m);            // UNREACHABLE, SpfOptions, SpfResult, spf
#if NETGRAPH_CORE_PRIVATE_HOOKS
void register_concurrency_hooks(pybind11::module_& m);
void register_private_hooks(pybind11::module_& m);     // the test controls and diagnostics
#endif

} // namespace netgraph::core::bindings
