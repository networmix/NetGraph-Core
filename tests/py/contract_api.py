"""The promised public API (its interface listing), written out: per public name, its members
and, per callable, every overload's parameters as (name, kind, default) with kind "pos", "kw", "var"
or "varkw" and the default's source text. A member the listing states in prose rather than as a
declaration (an unchanged method or field) is "declared"; a property, field or enum member is
"attribute". test_public_surface.py compares the generated stub with it."""

CONTRACT_API = {
    "AbsoluteCost": {"__init__": [[("cost", "pos", None)]], "cost": "attribute"},
    "Admission": {"LOSSLESS": "attribute", "LOSSY": "attribute"},
    "BatchMode": {"ORDERED": "attribute", "SIMULTANEOUS": "attribute"},
    "CentralizedTE": {
        "__init__": [
            [
                ("source_allocation", "kw", None),
                ("support", "kw", "None"),
                ("tunnel_count", "kw", "None"),
                ("prefer_headroom", "kw", "False"),
            ]
        ],
        "prefer_headroom": "attribute",
        "source_allocation": "attribute",
        "support": "attribute",
        "tunnel_count": "attribute",
    },
    "CostStretch": {"__init__": [[("factor", "pos", None)]], "factor": "attribute"},
    "Direction": {"FORWARD": "declared", "REVERSE": "declared"},
    "DistributedTE": {
        "__init__": [
            [
                ("tunnel_count", "pos", None),
                ("cost_ceiling", "kw", "None"),
                ("min_unreserved_bandwidth", "kw", "0"),
                ("tie_break", "kw", "RouteTieBreak.FEWER_HOPS"),
                ("parallel_edges", "kw", "True"),
            ]
        ],
        "cost_ceiling": "attribute",
        "min_unreserved_bandwidth": "attribute",
        "parallel_edges": "attribute",
        "tie_break": "attribute",
        "tunnel_count": "attribute",
    },
    "Eligibility": {
        "edge_allowed": "declared",
        "graph_uid": "declared",
        "intersect": "declared",
        "node_allowed": "declared",
        "uid": "declared",
    },
    "EvaluationOptions": {
        "__init__": [
            [
                ("config", "pos", None),
                ("eligibility", "kw", "None"),
                ("cache", "kw", "None"),
                ("outputs", "kw", "None"),
                ("design_demand", "kw", "None"),
                ("offer", "kw", "None"),
            ]
        ],
        "cache": "attribute",
        "config": "attribute",
        "design_demand": "attribute",
        "eligibility": "attribute",
        "offer": "attribute",
        "outputs": "attribute",
    },
    "FlowLedger": {
        "__init__": [
            [
                ("graph", "pos", None),
                ("service_baseline", "pos", "None"),
                ("reservation_baseline", "kw", "None"),
            ]
        ],
        "carried": [[]],
        "copy": "declared",
        "edge_capacity_units": [[]],
        "graph": "declared",
        "record": "declared",
        "records": "declared",
        "reserved": [[]],
        "reset": [
            [
                ("service_baseline", "pos", "None"),
                ("reservation_baseline", "kw", "None"),
            ]
        ],
        "residual": "declared",
        "saturated_boundary_edges": [
            [("sources", "pos", None), ("eligibility", "pos", "None")]
        ],
        "uid": "declared",
        "unreserved": [[]],
    },
    "FlowPolicy": {
        "__init__": [
            [
                ("ledger", "pos", None),
                ("sources", "pos", None),
                ("destinations", "pos", None),
                ("config", "pos", None),
                ("eligibility", "pos", "None"),
                ("cache", "pos", "None"),
            ]
        ],
        "config": "attribute",
        "delivered": [[]],
        "evaluate": [
            [
                ("offer", "pos", None),
                ("design_demand", "kw", "None"),
                ("outputs", "kw", "None"),
            ]
        ],
        "forward": [
            [
                ("offer", "pos", None),
                ("black_holes", "kw", "None"),
                ("outputs", "kw", "None"),
            ]
        ],
        "place": [[("offer", "pos", None), ("outputs", "pos", "None")]],
        "program": "attribute",
        "record": [[("flow_id", "pos", None)]],
        "records": [[]],
        "remove": [[]],
        "replace": [
            [
                ("offer", "pos", None),
                ("design_demand", "kw", "None"),
                ("outputs", "kw", "None"),
            ]
        ],
    },
    "FlowPolicyConfig": {
        "__init__": [
            [("routing", "kw", None), ("admission", "kw", "Admission.LOSSLESS")]
        ],
        "admission": "attribute",
        "routing": "attribute",
    },
    "FlowRecord": {
        "__init__": [[("_no_public_constructor", "pos", None)]],
        "admission": "attribute",
        "delivered": "attribute",
        "destinations": "attribute",
        "edge_drops": "attribute",
        "edge_flows": "attribute",
        "edge_reservations": "attribute",
        "entry_id": "attribute",
        "flow_id": "attribute",
        "min_route_cost": "attribute",
        "offered": "attribute",
        "policy_id": "attribute",
        "route": "attribute",
        "sources": "attribute",
    },
    "Forward": {
        "__init__": [
            [
                ("policy", "pos", None),
                ("offer", "pos", None),
                ("black_holes", "kw", "None"),
            ]
        ],
        "black_holes": "attribute",
        "offer": "attribute",
        "policy": "attribute",
    },
    "ForwardingProgram": {
        "__init__": [[("_no_public_constructor", "pos", None)]],
        "destinations": "attribute",
        "from_dag": [
            [
                ("graph", "pos", None),
                ("sources", "pos", None),
                ("destinations", "pos", None),
                ("dag", "pos", None),
                ("allocation", "kw", "None"),
                ("weight_basis", "kw", "None"),
                ("edge_weights", "kw", "None"),
                ("residual", "kw", "None"),
            ]
        ],
        "from_routes": [
            [
                ("graph", "pos", None),
                ("sources", "pos", None),
                ("destinations", "pos", None),
                ("routes", "pos", None),
                ("allocation", "kw", "None"),
                ("weight_basis", "kw", "None"),
                ("edge_weights", "kw", "None"),
                ("residual", "kw", "None"),
            ]
        ],
        "source_weights": "attribute",
        "sources": "attribute",
        "tunnels": "attribute",
    },
    "Graph": {
        "__init__": [
            [
                ("num_nodes", "pos", None),
                ("edge_sources", "pos", None),
                ("edge_destinations", "pos", None),
                ("edge_capacity_units", "pos", None),
                ("edge_costs", "pos", None),
                ("resolution", "pos", None),
            ]
        ],
        "edge_capacity_units": [[]],
        "edge_cost_gcd": [[]],
        "edge_costs": [[]],
        "edge_destinations": [[]],
        "edge_sources": [[]],
        "eligibility": "declared",
        "incoming_edges": [[]],
        "incoming_offsets": [[]],
        "input_positions": [[]],
        "max_edge_cost": [[]],
        "num_edges": "declared",
        "num_nodes": "declared",
        "outgoing_offsets": [[]],
        "resolution": "declared",
        "uid": "declared",
    },
    "KspOptions": {
        "__init__": [
            [
                ("k", "pos", None),
                ("tie_break", "kw", "RouteTieBreak.FEWER_HOPS"),
                ("parallel_edges", "kw", "True"),
                ("cost_ceiling", "kw", "None"),
                ("residual", "kw", "None"),
                ("eligibility", "kw", "None"),
            ]
        ],
        "cost_ceiling": "attribute",
        "eligibility": "attribute",
        "k": "attribute",
        "parallel_edges": "attribute",
        "residual": "attribute",
        "tie_break": "attribute",
    },
    "OutputOptions": {
        "__init__": [
            [
                ("edge_flows", "kw", "False"),
                ("residual", "kw", "False"),
                ("saturated_boundary_edges", "kw", "False"),
                ("edge_reservations", "kw", "False"),
                ("min_cut", "kw", "False"),
            ]
        ],
        "edge_flows": "attribute",
        "edge_reservations": "attribute",
        "min_cut": "attribute",
        "residual": "attribute",
        "saturated_boundary_edges": "attribute",
    },
    "Path": {
        "cost": "declared",
        "destination": [[]],
        "graph_uid": "declared",
        "hops": "declared",
        "source": [[]],
        "to_dag": [[]],
    },
    "Pinned": {
        "__init__": [[("program", "pos", None), ("black_holes", "kw", "None")]],
        "black_holes": "attribute",
        "program": "attribute",
    },
    "Place": {
        "__init__": [[("policy", "pos", None), ("offer", "pos", None)]],
        "offer": "attribute",
        "policy": "attribute",
    },
    "PlacementResult": {
        "__init__": [[("_no_public_constructor", "pos", None)]],
        "admitted": "attribute",
        "by_source": "attribute",
        "by_tunnel": "attribute",
        "call_delivered": "attribute",
        "call_edge_drops": "attribute",
        "call_source_drops": "attribute",
        "call_undelivered": "attribute",
        "delivered": "attribute",
        "delivery_delta": "attribute",
        "edge_drops": "attribute",
        "edge_flows": "attribute",
        "edge_reservations": "attribute",
        "marginal_cost_curve": "attribute",
        "min_cut": "attribute",
        "offered": "attribute",
        "plan": "attribute",
        "program": "attribute",
        "requested": "attribute",
        "reservations_changed": "attribute",
        "residual": "attribute",
        "routes_changed": "attribute",
        "saturated_boundary_edges": "attribute",
        "source_drops": "attribute",
        "termination": "attribute",
        "total_cost": "attribute",
        "unadmitted": "attribute",
        "weights_changed": "attribute",
    },
    "PredecessorDAG": {
        "__init__": [[("_no_public_constructor", "pos", None)]],
        "from_csr": [
            [
                ("graph", "pos", None),
                ("parent_offsets", "pos", None),
                ("parents", "pos", None),
                ("via_edges", "pos", None),
            ]
        ],
        "graph_uid": "declared",
        "parent_offsets": "declared",
        "parents": "declared",
        "via_edges": "declared",
    },
    "Replace": {
        "__init__": [
            [
                ("policy", "pos", None),
                ("offer", "pos", None),
                ("design_demand", "kw", "None"),
            ]
        ],
        "design_demand": "attribute",
        "offer": "attribute",
        "policy": "attribute",
    },
    "RouteAllocation": {
        "EQUAL": "attribute",
        "OPTIMIZED": "attribute",
        "WEIGHTED": "attribute",
    },
    "RouteDefinition": {
        "__init__": [
            [("head", "pos", None), ("route", "pos", None), ("weight", "kw", "1")]
        ],
        "head": "attribute",
        "route": "attribute",
        "weight": "attribute",
    },
    "RouteTieBreak": {
        "FEWER_HOPS": "attribute",
        "HIGHER_AVAILABLE_CAPACITY": "attribute",
        "HIGHER_CAPACITY": "attribute",
        "HIGHER_HEADROOM": "attribute",
    },
    "SensitivityEntry": {
        "__init__": [[("_no_public_constructor", "pos", None)]],
        "delivery_reduction": "attribute",
        "edge": "attribute",
    },
    "SensitivityResult": {"baseline": "declared", "entries": "declared"},
    "ShortestPathCache": {
        "clear": "declared",
        "freeze": "declared",
        "graph_uid": "declared",
        "stats": "declared",
        "thaw": "declared",
    },
    "ShortestPathCacheMiss": {},
    "ShortestPathCacheStats": {
        "bytes": "declared",
        "entries": "declared",
        "hits": "declared",
        "live_values": "declared",
        "misses": "declared",
    },
    "ShortestPaths": {
        "__init__": [
            [
                ("allocation", "kw", "RouteAllocation.EQUAL"),
                ("source_allocation", "kw", "None"),
                ("weight_basis", "kw", "None"),
                ("edge_weights", "kw", "None"),
                ("tie_break", "kw", "RouteTieBreak.FEWER_HOPS"),
                ("multipath", "kw", "True"),
                ("parallel_edges", "kw", "True"),
            ]
        ],
        "allocation": "attribute",
        "edge_weights": "attribute",
        "multipath": "attribute",
        "parallel_edges": "attribute",
        "source_allocation": "attribute",
        "tie_break": "attribute",
        "weight_basis": "attribute",
    },
    "SourceAllocation": {"EQUAL": "attribute", "OPTIMIZED": "attribute"},
    "SourceFlow": {
        "__init__": [[("_no_public_constructor", "pos", None)]],
        "admitted": "attribute",
        "delivered": "attribute",
        "offered": "attribute",
        "source": "attribute",
        "source_drops": "attribute",
        "total_edge_drops": "attribute",
    },
    "SpfOptions": {
        "__init__": [
            [
                ("direction", "pos", None),
                ("targets", "kw", "None"),
                ("tie_break", "kw", "RouteTieBreak.FEWER_HOPS"),
                ("multipath", "kw", "True"),
                ("parallel_edges", "kw", "True"),
                ("residual", "kw", "None"),
                ("min_residual", "kw", "1"),
                ("eligibility", "kw", "None"),
            ]
        ],
        "direction": "attribute",
        "eligibility": "attribute",
        "min_residual": "attribute",
        "multipath": "attribute",
        "parallel_edges": "attribute",
        "residual": "attribute",
        "targets": "attribute",
        "tie_break": "attribute",
    },
    "SpfResult": {"costs": "declared", "dag": "declared", "target_cost": "declared"},
    "Termination": {
        "COMPLETE": "attribute",
        "EXHAUSTED": "attribute",
        "NO_ACTIVE_WEIGHT": "attribute",
        "NO_ROUTE": "attribute",
    },
    "Tunnel": {
        "__init__": [[("_no_public_constructor", "pos", None)]],
        "destination": "attribute",
        "edges": "attribute",
        "entry_id": "attribute",
        "head": "attribute",
        "local_index": "attribute",
        "path_cost": "attribute",
        "planned_bandwidth": "attribute",
        "requested_bandwidth": "attribute",
        "reserved_bandwidth": "attribute",
        "status": "attribute",
        "weight": "attribute",
    },
    "TunnelCount": {
        "__init__": [[("per_head", "pos", None), ("up_to", "kw", "False")]],
        "per_head": "attribute",
        "up_to": "attribute",
    },
    "TunnelFlow": {
        "__init__": [[("_no_public_constructor", "pos", None)]],
        "assigned": "attribute",
        "delivered": "attribute",
        "entry_id": "attribute",
        "head": "attribute",
        "status": "attribute",
        "total_edge_drops": "attribute",
    },
    "TunnelPlan": {
        "__init__": [[("_no_public_constructor", "pos", None)]],
        "design_by_head": "attribute",
        "design_demand": "attribute",
        "headroom_penalty": "attribute",
        "heads_with_capacity": "attribute",
        "insufficient_resolution": "attribute",
        "max_bandwidth_per_tunnel": "attribute",
        "max_fractional_bandwidth_per_tunnel": "attribute",
        "max_lossless_offer": "attribute",
        "planned_bandwidth_per_tunnel": "attribute",
        "total_planned_bandwidth": "attribute",
        "total_reserved_bandwidth": "attribute",
        "total_tunnel_cost": "attribute",
        "tunnels_per_head": "attribute",
    },
    "TunnelStatus": {
        "BLACKHOLED": "attribute",
        "COST_CEILING_EXCEEDED": "attribute",
        "FAILED": "attribute",
        "INSUFFICIENT_BANDWIDTH": "attribute",
        "UP": "attribute",
    },
    "WeightBasis": {
        "CAPACITY": "attribute",
        "EXPLICIT": "attribute",
        "RESIDUAL_SNAPSHOT": "attribute",
    },
    "apply_actions": {
        "__call__": [
            [
                ("ledger", "pos", None),
                ("actions", "pos", None),
                ("mode", "kw", "BatchMode.ORDERED"),
                ("outputs", "kw", "None"),
            ]
        ]
    },
    "as_float": {"__call__": [[("x", "pos", None), ("resolution", "pos", None)]]},
    "evaluate_demand": {
        "__call__": [
            [
                ("graph", "pos", None),
                ("sources", "pos", None),
                ("destinations", "pos", None),
                ("options", "pos", None),
            ]
        ]
    },
    "evaluate_edge_failures": {
        "__call__": [
            [
                ("graph", "pos", None),
                ("sources", "pos", None),
                ("destinations", "pos", None),
                ("options", "pos", None),
                ("workers", "pos", "None"),
            ]
        ]
    },
    "evaluate_pairs": {
        "__call__": [
            [
                ("graph", "pos", None),
                ("pairs", "pos", None),
                ("options", "pos", None),
                ("per_pair_eligibility", "pos", "None"),
                ("workers", "pos", "None"),
            ]
        ]
    },
    "ksp": {
        "__call__": [
            [
                ("graph", "pos", None),
                ("source", "pos", None),
                ("destination", "pos", None),
                ("options", "pos", None),
            ]
        ]
    },
    "resolve_to_paths": {
        "__call__": [
            [
                ("graph", "pos", None),
                ("dag", "pos", None),
                ("source", "pos", None),
                ("destination", "pos", None),
                ("split_parallel_edges", "pos", "False"),
                ("max_paths", "pos", "None"),
            ]
        ]
    },
    "spf": {
        "__call__": [
            [("graph", "pos", None), ("roots", "pos", None), ("options", "pos", None)]
        ]
    },
    "to_units": {"__call__": [[("x", "pos", None), ("resolution", "pos", None)]]},
}
