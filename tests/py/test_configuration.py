"""The configuration values (docs/semantics.md, "Configuration"): their construction (keywords,
required arguments, types, defaults), the model matrix in its order and its first error, the
route choice's rules (the ineffective tie rule of ShortestPaths and the pooled single-path rule),
at a policy and at an evaluator. The table cases are the worked CFG examples, run here on the
library: a valid model reads back its resolved fields from the policy's config."""

from __future__ import annotations

import ast
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
from helpers import assert_read_only, build, dag_of, exact

import netgraph_core as ngc

TIE = ngc.RouteTieBreak
A = ngc.RouteAllocation
INEFFECTIVE = "ShortestPaths tie_break requires multipath=False or parallel_edges=False"
POOLED = (
    "ShortestPaths with multiple input sources requires tie_break=FEWER_HOPS"
    " when multipath=False"
)


def config(allocation, **choice):
    return ngc.FlowPolicyConfig(
        routing=ngc.ShortestPaths(allocation=allocation, **choice)
    )


def test_ineffective_tie():
    # A HIGHER_* tie_break with every route and every parallel edge has nothing to choose:
    # rejected for every allocation (OPTIMIZED's support is the route choice too); with
    # multipath=False or parallel_edges=False it is valid; FEWER_HOPS always is.
    g = build(2, [(0, 1, 4, 1), (0, 1, 10, 1)])
    for allocation in (A.EQUAL, A.WEIGHTED, A.OPTIMIZED):
        for tie in (
            TIE.HIGHER_CAPACITY,
            TIE.HIGHER_AVAILABLE_CAPACITY,
            TIE.HIGHER_HEADROOM,
        ):
            cfg = config(allocation, tie_break=tie)
            with pytest.raises(ValueError, match=exact(INEFFECTIVE)):
                ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg)
            with pytest.raises(ValueError, match=exact(INEFFECTIVE)):
                ngc.evaluate_demand(g, 0, 1, ngc.EvaluationOptions(cfg))
            for choice in ({"multipath": False}, {"parallel_edges": False}):
                ok = config(allocation, tie_break=tie, **choice)
                assert (
                    ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, ok).place(1).call_delivered
                    == 1
                )
        assert (
            ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, config(allocation))
            .place(1)
            .call_delivered
            == 1
        )


def test_pooled_rule_counts_inputs():
    # The pooled single-path rule counts the normalized input sources before the mask:
    # [0, 0, 1] is two sources, [0, 1] with source 1 masked is still two, [0] is one.
    g = build(3, [(0, 2, 5, 1), (1, 2, 5, 1)])
    cfg = config(A.EQUAL, tie_break=TIE.HIGHER_CAPACITY, multipath=False)
    with pytest.raises(ValueError, match=exact(POOLED)):
        ngc.FlowPolicy(ngc.FlowLedger(g), [0, 0, 1], 2, cfg)
    masked = g.eligibility(excluded_nodes=[1])
    with pytest.raises(ValueError, match=exact(POOLED)):
        ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, cfg, masked)
    with pytest.raises(ValueError, match=exact(POOLED)):
        ngc.evaluate_demand(
            g, [0, 1], 2, ngc.EvaluationOptions(cfg, eligibility=masked)
        )
    assert ngc.FlowPolicy(ngc.FlowLedger(g), [0], 2, cfg).place(1).call_delivered == 1
    # FEWER_HOPS single paths pool, and HIGHER_* multipath (one member per pair) pools.
    for choice in (
        {"multipath": False},
        {"tie_break": TIE.HIGHER_CAPACITY, "parallel_edges": False},
    ):
        ok = config(A.EQUAL, **choice)
        assert (
            ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, ok).place(2).call_delivered
            == 2
        )


# -- the worked CFG examples ----------------------------------------------------------------------

# 0 -> 1 (e0), 0 -> 2 (e1), 1 -> 2 (e2): one source is [0], two are [0, 1]; the destination is 2.
EXAMPLE = [(0, 1, 5, 1), (0, 2, 5, 1), (1, 2, 5, 1)]
ENUMS = {
    "allocation": ngc.RouteAllocation,
    "source_allocation": ngc.SourceAllocation,
    "weight_basis": ngc.WeightBasis,
    "tie_break": ngc.RouteTieBreak,
}
MODELS = {
    "ShortestPaths": ngc.ShortestPaths,
    "DistributedTE": lambda tunnel_count, **kw: ngc.DistributedTE(tunnel_count, **kw),
    "CentralizedTE": ngc.CentralizedTE,
    "Pinned": lambda program, **kw: ngc.Pinned(program, **kw),
    "FlowPolicyConfig": ngc.FlowPolicyConfig,
}


CEILING = "cost_ceiling must be AbsoluteCost, CostStretch or None"


def program_of(g, spec):
    """A program of the example graph: "p" from 0 to 2, or {"sources": [0], "destinations": [d]}."""
    if spec == "p":
        spec = {"sources": [0], "destinations": [2]}
    (d,) = spec["destinations"]
    route = ngc.Path(g, [[0 if d == 1 else 1]])
    return ngc.ForwardingProgram.from_routes(
        g, spec["sources"], [d], [ngc.RouteDefinition(0, route)]
    )


def argument(g, key, value):
    """The library value of one example keyword (the examples spell values as strings)."""
    if key in ENUMS and isinstance(value, str):
        return getattr(ENUMS[key], value)
    if key == "tunnel_count":
        count, *up_to = value.split()
        return ngc.TunnelCount(int(count), up_to=up_to == ["up_to"])
    if key == "cost_ceiling":
        kind, *parts = value
        return ngc.CostStretch(Fraction(*parts)) if kind == "fraction" else parts[0]
    if key == "support":
        return dag_of(g, [1])
    if key == "program":
        return program_of(g, value)
    if key == "routing":
        return ngc.ShortestPaths()
    return value


def resolved(routing):
    """A valid model's resolved fields, as the examples print them."""

    def name(v):
        return "None" if v is None else v.name

    if isinstance(routing, ngc.ShortestPaths):
        return (
            f"allocation={name(routing.allocation)}"
            f" source_allocation={name(routing.source_allocation)}"
            f" weight_basis={name(routing.weight_basis)}"
            f" tie_break={name(routing.tie_break)}"
            f" multipath={routing.multipath} parallel_edges={routing.parallel_edges}"
        )
    if isinstance(routing, ngc.CentralizedTE):
        t = routing.tunnel_count
        count = "None" if t is None else f"{t.per_head}" + (" up_to" if t.up_to else "")
        return (
            f"source_allocation={name(routing.source_allocation)}"
            f" support={'None' if routing.support is None else 'dag'}"
            f" tunnel_count={count} prefer_headroom={routing.prefer_headroom}"
        )
    return "valid"


def outcome(
    model, kwargs, sources=1, admission="LOSSLESS", cache=False, endpoints=None
):
    """One CFG example on the library: "TypeError" (with the text where the example names one)
    at construction, the policy's first ValueError text, or the resolved fields of a valid model.
    The policy is built, and a zero call placed, only when the model is valid (a zero offer skips
    no validation)."""
    g = build(3, EXAMPLE)
    try:
        value = MODELS[model](**{k: argument(g, k, v) for k, v in kwargs.items()})
    except TypeError as e:
        return f"TypeError: {e}" if str(e) == CEILING else "TypeError"
    if model == "FlowPolicyConfig":
        return "constructed"
    cfg = ngc.FlowPolicyConfig(
        routing=value, admission=getattr(ngc.Admission, admission)
    )
    src, dst = ([0], [2]) if sources == 1 else ([0, 1], [2])
    if endpoints is not None:
        src, dst = endpoints
    tunnels = isinstance(value, ngc.DistributedTE) or (
        isinstance(value, ngc.CentralizedTE) and value.tunnel_count is not None
    )
    try:
        policy = ngc.FlowPolicy(
            ngc.FlowLedger(g),
            src,
            dst,
            cfg,
            None,
            ngc.ShortestPathCache(g) if cache else None,
        )
    except ValueError as e:
        return str(e)
    (policy.replace if tunnels else policy.place)(0)
    return resolved(policy.config.routing)


def run_cases(cases, expected):
    got = [outcome(**case) for case in cases]
    assert got == expected


def test_cfg_defaults():
    # CFG-defaults: every resolved field of ShortestPaths' allocations, and a T2 model.
    run_cases(
        [
            dict(model="ShortestPaths", kwargs={}),
            dict(model="ShortestPaths", kwargs={"allocation": "EQUAL"}),
            dict(model="ShortestPaths", kwargs={"allocation": "WEIGHTED"}),
            dict(model="ShortestPaths", kwargs={"allocation": "OPTIMIZED"}),
            dict(
                model="CentralizedTE",
                kwargs={"source_allocation": "EQUAL", "tunnel_count": "16"},
            ),
        ],
        [
            "allocation=EQUAL source_allocation=EQUAL weight_basis=None tie_break=FEWER_HOPS multipath=True parallel_edges=True",
            "allocation=EQUAL source_allocation=EQUAL weight_basis=None tie_break=FEWER_HOPS multipath=True parallel_edges=True",
            "allocation=WEIGHTED source_allocation=EQUAL weight_basis=CAPACITY tie_break=FEWER_HOPS multipath=True parallel_edges=True",
            "allocation=OPTIMIZED source_allocation=OPTIMIZED weight_basis=None tie_break=FEWER_HOPS multipath=True parallel_edges=True",
            "source_allocation=EQUAL support=None tunnel_count=16 prefer_headroom=False",
        ],
    )


def test_cfg_construction():
    # CFG-construction: every case is a TypeError when the value is built.
    run_cases(
        [
            dict(model="CentralizedTE", kwargs={}),
            dict(model="CentralizedTE", kwargs={"source_allocation": None}),
            dict(model="CentralizedTE", kwargs={"tunnel_count": "16"}),
            dict(
                model="CentralizedTE",
                kwargs={"source_allocation": "EQUAL", "selection": "default"},
            ),
            dict(
                model="DistributedTE",
                kwargs={"tunnel_count": "4", "source_allocation": "EQUAL"},
            ),
            dict(
                model="DistributedTE", kwargs={"tunnel_count": "4", "multipath": False}
            ),
            dict(model="Pinned", kwargs={"program": "p", "allocation": "EQUAL"}),
            dict(
                model="FlowPolicyConfig",
                kwargs={"routing": "ShortestPaths()", "allocation": "EQUAL"},
            ),
            dict(model="ShortestPaths", kwargs={"multipath": None}),
            dict(model="ShortestPaths", kwargs={"placement": "EQUAL"}),
        ],
        ["TypeError"] * 10,
    )


def test_cfg_matrix():
    # CFG-matrix: one result per row of the matrix (the rows checked at a policy, each valid
    # model also placing a zero call).
    run_cases(
        [
            dict(model="ShortestPaths", kwargs={"source_allocation": "OPTIMIZED"}),
            dict(
                model="ShortestPaths",
                kwargs={"allocation": "OPTIMIZED"},
                admission="LOSSY",
            ),
            dict(model="ShortestPaths", kwargs={"tie_break": "HIGHER_CAPACITY"}),
            dict(
                model="ShortestPaths",
                kwargs={"tie_break": "HIGHER_CAPACITY", "multipath": False},
                sources=2,
            ),
            dict(
                model="ShortestPaths",
                kwargs={"tie_break": "HIGHER_CAPACITY", "multipath": False},
            ),
            dict(
                model="ShortestPaths",
                kwargs={"tie_break": "HIGHER_HEADROOM", "parallel_edges": False},
                sources=2,
            ),
            dict(model="ShortestPaths", kwargs={"weight_basis": "CAPACITY"}),
            dict(
                model="ShortestPaths",
                kwargs={"allocation": "WEIGHTED", "weight_basis": "EXPLICIT"},
            ),
            dict(
                model="CentralizedTE",
                kwargs={"source_allocation": "OPTIMIZED", "tunnel_count": "2"},
                sources=2,
            ),
            dict(
                model="CentralizedTE",
                kwargs={"source_allocation": "OPTIMIZED", "tunnel_count": "2"},
            ),
            dict(
                model="CentralizedTE",
                kwargs={"source_allocation": "EQUAL"},
                admission="LOSSY",
            ),
            dict(
                model="CentralizedTE",
                kwargs={"source_allocation": "EQUAL", "prefer_headroom": True},
            ),
            dict(
                model="CentralizedTE",
                kwargs={
                    "source_allocation": "EQUAL",
                    "support": "dag",
                    "tunnel_count": "2",
                },
                admission="LOSSY",
            ),
            dict(model="DistributedTE", kwargs={"tunnel_count": "4 up_to"}),
        ],
        [
            "ShortestPaths with EQUAL or WEIGHTED requires source_allocation=EQUAL",
            "ShortestPaths OPTIMIZED requires admission=LOSSLESS",
            INEFFECTIVE,
            POOLED,
            "allocation=EQUAL source_allocation=EQUAL weight_basis=None tie_break=HIGHER_CAPACITY multipath=False parallel_edges=True",
            "allocation=EQUAL source_allocation=EQUAL weight_basis=None tie_break=HIGHER_HEADROOM multipath=True parallel_edges=False",
            "weight_basis requires allocation=WEIGHTED",
            "weight_basis=EXPLICIT requires edge_weights",
            "CentralizedTE with tunnel_count and several sources requires source_allocation=EQUAL",
            "source_allocation=OPTIMIZED support=None tunnel_count=2 prefer_headroom=False",
            "CentralizedTE without tunnel_count requires admission=LOSSLESS",
            "prefer_headroom requires tunnel_count",
            "source_allocation=EQUAL support=dag tunnel_count=2 prefer_headroom=False",
            "DistributedTE requires an exact tunnel count",
        ],
    )


def test_cfg_first_error():
    # CFG-first-error: two defects each; the first in the model matrix's order is reported.
    run_cases(
        [
            dict(
                model="CentralizedTE",
                kwargs={"source_allocation": "EQUAL", "prefer_headroom": True},
                admission="LOSSY",
            ),
            dict(
                model="ShortestPaths",
                kwargs={"source_allocation": "OPTIMIZED", "weight_basis": "CAPACITY"},
            ),
            dict(
                model="ShortestPaths",
                kwargs={"allocation": "OPTIMIZED", "tie_break": "HIGHER_CAPACITY"},
                admission="LOSSY",
            ),
            dict(
                model="ShortestPaths",
                kwargs={"tie_break": "HIGHER_CAPACITY", "weight_basis": "CAPACITY"},
            ),
            dict(
                model="ShortestPaths",
                kwargs={
                    "tie_break": "HIGHER_CAPACITY",
                    "multipath": False,
                    "weight_basis": "CAPACITY",
                },
                sources=2,
            ),
            dict(
                model="CentralizedTE",
                kwargs={
                    "source_allocation": "OPTIMIZED",
                    "tunnel_count": "0",
                    "prefer_headroom": True,
                },
                sources=2,
            ),
            dict(
                model="DistributedTE",
                kwargs={"tunnel_count": "0 up_to", "cost_ceiling": ["fraction", 1, 2]},
            ),
            dict(
                model="DistributedTE",
                kwargs={"tunnel_count": "4 up_to", "cost_ceiling": ["fraction", 1, 2]},
            ),
            dict(
                model="DistributedTE",
                kwargs={
                    "tunnel_count": "4",
                    "cost_ceiling": ["fraction", 1, 2],
                    "min_unreserved_bandwidth": -1,
                },
            ),
            dict(
                model="DistributedTE",
                kwargs={"tunnel_count": "4", "min_unreserved_bandwidth": -1},
            ),
        ],
        [
            "CentralizedTE without tunnel_count requires admission=LOSSLESS",
            "ShortestPaths with EQUAL or WEIGHTED requires source_allocation=EQUAL",
            "ShortestPaths OPTIMIZED requires admission=LOSSLESS",
            INEFFECTIVE,
            POOLED,
            "CentralizedTE with tunnel_count and several sources requires source_allocation=EQUAL",
            "tunnel count must be in [1, 10000]",
            "DistributedTE requires an exact tunnel count",
            "cost_ceiling CostStretch factor must be >= 1 (got 1/2)",
            "min_unreserved_bandwidth must be an integer in [0, 2^62)",
        ],
    )


def test_cfg_cache_order():
    # CFG-cache-order: the cache refusal is the matrix's last check (after prefer_headroom,
    # before the tunnel count), the program's endpoints come after the matrix.
    other = {"sources": [0], "destinations": [2]}
    run_cases(
        [
            dict(model="DistributedTE", kwargs={"tunnel_count": "0"}, cache=True),
            dict(
                model="Pinned",
                kwargs={"program": other},
                cache=True,
                endpoints=[[0], [1]],
            ),
            dict(
                model="CentralizedTE",
                kwargs={"source_allocation": "EQUAL", "tunnel_count": "0"},
                cache=True,
            ),
            dict(
                model="CentralizedTE",
                kwargs={"source_allocation": "EQUAL", "prefer_headroom": True},
                cache=True,
            ),
            dict(model="ShortestPaths", kwargs={"allocation": "WEIGHTED"}, cache=True),
            dict(model="DistributedTE", kwargs={"tunnel_count": "0"}),
            dict(model="Pinned", kwargs={"program": other}, endpoints=[[0], [1]]),
            dict(
                model="Pinned",
                kwargs={"program": {"sources": [0], "destinations": [1]}},
                endpoints=[[0], [1]],
            ),
        ],
        [
            "DistributedTE does not use ShortestPathCache",
            "Pinned does not use ShortestPathCache",
            "CentralizedTE does not use ShortestPathCache",
            "prefer_headroom requires tunnel_count",
            "allocation=WEIGHTED source_allocation=EQUAL weight_basis=CAPACITY tie_break=FEWER_HOPS multipath=True parallel_edges=True",
            "tunnel count must be in [1, 10000]",
            "program endpoints do not match policy endpoints",
            "valid",
        ],
    )


def test_cfg_ceiling_type():
    # CFG-ceiling-type: a ceiling of another type is refused when DistributedTE is built, before
    # any policy check (the count 0 is never reached).
    ceiling = f"TypeError: {CEILING}"
    run_cases(
        [
            dict(
                model="DistributedTE",
                kwargs={"tunnel_count": "0", "cost_ceiling": ["other", 1.5]},
            ),
            dict(
                model="DistributedTE",
                kwargs={"tunnel_count": "4", "cost_ceiling": ["other", 2]},
            ),
            dict(
                model="DistributedTE",
                kwargs={"tunnel_count": "4", "cost_ceiling": ["fraction", 3, 2]},
            ),
        ],
        [ceiling, ceiling, "valid"],
    )


# -- construction ------------------------------------------------------------------------------


def test_shortest_paths_defaults():
    # ShortestPaths() is ShortestPaths(allocation=EQUAL): the same fields as built and the same
    # resolved snapshot (the source allocation resolved to EQUAL).
    g = build(3, EXAMPLE)
    fields = (
        "allocation",
        "source_allocation",
        "weight_basis",
        "edge_weights",
        "tie_break",
        "multipath",
        "parallel_edges",
    )
    built = (A.EQUAL, None, None, None, TIE.FEWER_HOPS, True, True)
    snapshots = []
    for sp in (ngc.ShortestPaths(), ngc.ShortestPaths(allocation=A.EQUAL)):
        assert tuple(getattr(sp, f) for f in fields) == built
        cfg = ngc.FlowPolicyConfig(routing=sp)
        routing = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 2, cfg).config.routing
        snapshots.append(tuple(getattr(routing, f) for f in fields))
    assert (
        snapshots
        == [
            (
                A.EQUAL,
                ngc.SourceAllocation.EQUAL,
                None,
                None,
                TIE.FEWER_HOPS,
                True,
                True,
            )
        ]
        * 2
    )


def test_resolved_snapshot():
    # WEIGHTED resolves source EQUAL and basis CAPACITY; OPTIMIZED source OPTIMIZED and no
    # basis; OPTIMIZED with source EQUAL keeps EQUAL. The model value itself is unchanged.
    g = build(3, EXAMPLE)
    S = ngc.SourceAllocation
    for sp, want in (
        (ngc.ShortestPaths(allocation=A.WEIGHTED), (S.EQUAL, ngc.WeightBasis.CAPACITY)),
        (ngc.ShortestPaths(allocation=A.OPTIMIZED), (S.OPTIMIZED, None)),
        (
            ngc.ShortestPaths(allocation=A.OPTIMIZED, source_allocation=S.EQUAL),
            (S.EQUAL, None),
        ),
    ):
        before = repr(sp)
        cfg = ngc.FlowPolicyConfig(routing=sp)
        routing = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, cfg).config.routing
        assert (routing.source_allocation, routing.weight_basis) == want, sp
        assert repr(sp) == before and repr(cfg.routing) == before


def test_centralized_requires_source_allocation():
    with pytest.raises(TypeError):
        ngc.CentralizedTE()
    with pytest.raises(TypeError):
        ngc.CentralizedTE(source_allocation=None)
    c = ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL)
    assert (c.source_allocation, c.support, c.tunnel_count, c.prefer_headroom) == (
        ngc.SourceAllocation.EQUAL,
        None,
        None,
        False,
    )


def test_unknown_keywords():
    # A field of another model, or a retired one, is no keyword: pybind11's own TypeError.
    g = build(3, EXAMPLE)
    program = program_of(g, "p")
    sp = ngc.ShortestPaths()
    cases = [
        lambda: ngc.FlowPolicyConfig(routing=sp, placement=A.EQUAL),
        lambda: ngc.FlowPolicyConfig(routing=sp, selection=None),
        lambda: ngc.FlowPolicyConfig(routing=sp, origin=ngc.SourceAllocation.EQUAL),
        lambda: ngc.FlowPolicyConfig(routing=sp, weights=ngc.WeightBasis.CAPACITY),
        lambda: ngc.FlowPolicyConfig(routing=sp, edge_weights=[1, 1, 1]),
        lambda: ngc.FlowPolicyConfig(routing=sp, allocation=A.EQUAL),
        lambda: ngc.DistributedTE(
            ngc.TunnelCount(2), source_allocation=ngc.SourceAllocation.EQUAL
        ),
        lambda: ngc.DistributedTE(ngc.TunnelCount(2), multipath=False),
        lambda: ngc.Pinned(program, allocation=A.EQUAL),
        lambda: ngc.Pinned(program, source_allocation=ngc.SourceAllocation.EQUAL),
        lambda: ngc.KspOptions(2, multipath=False),
    ]
    for case in cases:
        with pytest.raises(
            TypeError, match="incompatible (function|constructor) arguments"
        ):
            case()


def test_missing_required():
    for case in (ngc.FlowPolicyConfig, ngc.DistributedTE, ngc.Pinned):
        with pytest.raises(TypeError):
            case()


class _Truthy:
    def __bool__(self):
        return True


# Classes that only borrow the type name of NumPy's bool (both spellings): a bool parameter
# checks the type itself.
FAKE_NUMPY_BOOLS = [
    type(name, (), {"__bool__": lambda self: True})
    for name in ("numpy.bool_", "numpy.bool")
]


def _stub_bool_parameters():
    # Every parameter the installed stub annotates `bool`, as "Owner.method:name" or
    # "function:name".
    stub = ast.parse(Path(ngc.__file__).with_name("__init__.pyi").read_text())
    out = set()

    def walk(body, owner=None):
        for node in body:
            if isinstance(node, ast.ClassDef):
                walk(node.body, node.name)
            elif isinstance(node, ast.FunctionDef):
                args = node.args
                for a in args.posonlyargs + args.args + args.kwonlyargs:
                    if a.annotation is not None and ast.unparse(a.annotation) == "bool":
                        out.add(f"{owner + '.' if owner else ''}{node.name}:{a.arg}")

    walk(stub.body)
    return out


def _bool_parameters():
    # Every bool parameter of the API (the configuration values and the search options), with
    # how to read the value back; test_bools_are_strict checks the list against the stub.
    g = build(2, [(0, 1, 1, 1), (0, 1, 1, 1)])
    dag = ngc.PredecessorDAG.from_csr(g, [0, 0, 2], [0, 0], [0, 1])
    S = ngc.SourceAllocation.EQUAL
    params = [
        (
            "ShortestPaths.__init__:multipath",
            lambda v: ngc.ShortestPaths(multipath=v).multipath,
        ),
        (
            "ShortestPaths.__init__:parallel_edges",
            lambda v: ngc.ShortestPaths(parallel_edges=v).parallel_edges,
        ),
        (
            "DistributedTE.__init__:parallel_edges",
            lambda v: ngc.DistributedTE(
                ngc.TunnelCount(1), parallel_edges=v
            ).parallel_edges,
        ),
        (
            "CentralizedTE.__init__:prefer_headroom",
            lambda v: ngc.CentralizedTE(
                source_allocation=S, prefer_headroom=v
            ).prefer_headroom,
        ),
        ("TunnelCount.__init__:up_to", lambda v: ngc.TunnelCount(2, up_to=v).up_to),
        (
            "KspOptions.__init__:parallel_edges",
            lambda v: ngc.KspOptions(2, parallel_edges=v).parallel_edges,
        ),
        (
            "SpfOptions.__init__:multipath",
            lambda v: ngc.SpfOptions(ngc.Direction.FORWARD, multipath=v).multipath,
        ),
        (
            "SpfOptions.__init__:parallel_edges",
            lambda v: ngc.SpfOptions(
                ngc.Direction.FORWARD, parallel_edges=v
            ).parallel_edges,
        ),
        (
            "resolve_to_paths:split_parallel_edges",
            lambda v: len(ngc.resolve_to_paths(g, dag, 0, 1, split_parallel_edges=v))
            == 2,
        ),
    ]
    for field in (
        "edge_flows",
        "residual",
        "saturated_boundary_edges",
        "edge_reservations",
        "min_cut",
    ):
        params.append(
            (
                f"OutputOptions.__init__:{field}",
                lambda v, field=field: getattr(ngc.OutputOptions(**{field: v}), field),
            )
        )
    return params


def test_bools_are_strict():
    # A bool parameter takes exactly a bool or a numpy.bool_, both known by their type: not a
    # class named like NumPy's bool, an object with __bool__, None, an int or a string.
    params = _bool_parameters()
    assert {name for name, _ in params} == _stub_bool_parameters()
    assert len(params) == 14
    fakes = [fake() for fake in FAKE_NUMPY_BOOLS]
    for name, make in params:
        for bad in (*fakes, _Truthy(), None, 0, 1, "x", np.int64(1), np.array([True])):
            with pytest.raises(TypeError):
                make(bad)
                pytest.fail(f"{name} took {bad!r}")
        for good in (True, False, np.bool_(True), np.bool_(False)):
            assert make(good) is bool(good), name


def test_integer_scalars():
    # A configuration value's integer field: a value of another type is a TypeError naming the
    # field; beyond 64 bits a ValueError; in 64 bits it constructs and its domain is checked
    # where the value is used.
    for bad, kind in ((True, "bool"), (1.5, "float"), (None, "NoneType")):
        with pytest.raises(
            TypeError,
            match=exact(f"TunnelCount per_head must be an integer (got {kind})"),
        ):
            ngc.TunnelCount(bad)
    with pytest.raises(
        ValueError,
        match=exact(
            "TunnelCount per_head = 18446744073709551616 is out of range "
            "[-9223372036854775808, 9223372036854775807]"
        ),
    ):
        ngc.TunnelCount(2**64)
    g = build(3, EXAMPLE)

    def policy(routing):
        ngc.FlowPolicy(ngc.FlowLedger(g), 0, 2, ngc.FlowPolicyConfig(routing=routing))

    with pytest.raises(ValueError, match=exact("tunnel count must be in [1, 10000]")):
        policy(ngc.DistributedTE(ngc.TunnelCount(0)))
    with pytest.raises(
        TypeError, match=exact("AbsoluteCost cost must be an integer (got NoneType)")
    ):
        ngc.AbsoluteCost(None)
    with pytest.raises(
        ValueError, match=exact("cost_ceiling AbsoluteCost cost must be >= 0 (got -1)")
    ):
        policy(ngc.DistributedTE(ngc.TunnelCount(1), cost_ceiling=ngc.AbsoluteCost(-1)))
    with pytest.raises(
        TypeError,
        match=exact("min_unreserved_bandwidth must be an integer (got float)"),
    ):
        ngc.DistributedTE(ngc.TunnelCount(1), min_unreserved_bandwidth=1.0)
    with pytest.raises(
        ValueError, match=r"^min_unreserved_bandwidth = .* is out of range"
    ):
        ngc.DistributedTE(ngc.TunnelCount(1), min_unreserved_bandwidth=2**63)
    with pytest.raises(
        ValueError,
        match=exact("min_unreserved_bandwidth must be an integer in [0, 2^62)"),
    ):
        policy(ngc.DistributedTE(ngc.TunnelCount(1), min_unreserved_bandwidth=2**62))


def test_wrong_types():
    g = build(3, EXAMPLE)
    S = ngc.SourceAllocation.EQUAL
    for case, text in (
        (lambda: ngc.ShortestPaths(allocation="EQUAL"), None),
        (lambda: ngc.ShortestPaths(tie_break=0), None),
        (
            lambda: ngc.CentralizedTE(source_allocation=S, support=5),
            "CentralizedTE support must be a PredecessorDAG or None",
        ),
        (
            lambda: ngc.CentralizedTE(source_allocation=S, tunnel_count=2),
            "CentralizedTE tunnel_count must be a TunnelCount or None",
        ),
        (
            lambda: ngc.FlowPolicyConfig(routing=1),
            "routing must be ShortestPaths, DistributedTE, CentralizedTE or Pinned (got int)",
        ),
        (
            lambda: ngc.DistributedTE(ngc.TunnelCount(1), cost_ceiling=1.5),
            "cost_ceiling must be AbsoluteCost, CostStretch or None",
        ),
        (lambda: ngc.Pinned(None), "Pinned program must be a ForwardingProgram"),
    ):
        with pytest.raises(TypeError, match=None if text is None else exact(text)):
            case()
    assert g.num_edges() == 3


def test_edge_weights():
    # One integer per edge, read back as a tuple; an entry beyond 64 bits is refused at
    # construction; the array is the model's, shared by every policy built from it.
    g = build(4, [(0, 1, 10, 1), (0, 2, 10, 1), (1, 3, 10, 1), (2, 3, 10, 1)])
    sp = ngc.ShortestPaths(
        allocation=A.WEIGHTED,
        weight_basis=ngc.WeightBasis.EXPLICIT,
        edge_weights=np.array([1, 3, 1, 1]),
    )
    assert sp.edge_weights == (1, 3, 1, 1)
    assert all(type(w) is int for w in sp.edge_weights)
    with pytest.raises(
        ValueError, match=exact("edge_weights[0] must be an integer in [0, 2^62)")
    ):
        ngc.ShortestPaths(
            allocation=A.WEIGHTED,
            weight_basis=ngc.WeightBasis.EXPLICIT,
            edge_weights=[2**64, 1, 1, 1],
        )
    # A vector of the wrong type or shape is a TypeError at construction; a value out of range
    # stays a ValueError (above) or is checked where the model is used.
    for bad in (
        [True, 1, 1, 1],
        [1.5, 1, 1, 1],
        [np.bool_(True), 1, 1, 1],
        [None, 1, 1, 1],
        ["1", 1, 1, 1],
        [[1], 1, 1, 1],
        3,
        "1111",
        b"1111",
        {1, 3},
        np.array([True, False, True, True]),
        np.array([1.0, 3.0, 1.0, 1.0]),
        np.array([[1, 3, 1, 1]]),
        np.array(1),
    ):
        with pytest.raises(TypeError, match=r"^edge_weights"):
            ngc.ShortestPaths(
                allocation=A.WEIGHTED,
                weight_basis=ngc.WeightBasis.EXPLICIT,
                edge_weights=bad,
            )
    with pytest.raises(
        TypeError, match=exact("edge_weights[0] must be an integer (got bool)")
    ):
        ngc.ShortestPaths(edge_weights=[True])
    cfg = ngc.FlowPolicyConfig(routing=sp)
    for _ in range(2):
        r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, cfg).place(8)
        assert r.call_delivered == 8
    assert ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, cfg).config.routing.edge_weights == (
        1,
        3,
        1,
        1,
    )


def test_immutable_values():
    g = build(3, EXAMPLE)
    program = program_of(g, "p")
    for value, names in (
        (ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()), ("admission",)),
        (
            ngc.ShortestPaths(),
            (
                "allocation",
                "source_allocation",
                "weight_basis",
                "edge_weights",
                "tie_break",
                "multipath",
                "parallel_edges",
            ),
        ),
        (
            ngc.DistributedTE(ngc.TunnelCount(1)),
            (
                "tunnel_count",
                "cost_ceiling",
                "min_unreserved_bandwidth",
                "tie_break",
                "parallel_edges",
            ),
        ),
        (
            ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL),
            ("source_allocation", "support", "tunnel_count", "prefer_headroom"),
        ),
        (ngc.Pinned(program), ("program", "black_holes")),
        (ngc.TunnelCount(2), ("per_head", "up_to")),
        (ngc.AbsoluteCost(2), ("cost",)),
        (ngc.CostStretch(2), ("factor",)),
        (ngc.RouteDefinition(0, ngc.Path(g, [[1]])), ("head", "route", "weight")),
    ):
        assert_read_only(value, *names)
    # The routing reads back as a new Python object each time (a model has no equality).
    cfg = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths())
    assert type(ngc.FlowPolicyConfig.routing).__name__ == "property"
    assert ngc.FlowPolicyConfig.routing.fset is None
    with pytest.raises(AttributeError):
        cfg.routing = ngc.ShortestPaths()


# -- the matrix's order -----------------------------------------------------------------------


def test_pooled_rule_before_vectors():
    # The pooled single-path rule is part of the matrix: before edge_weights' length.
    g = build(3, [(0, 2, 5, 1), (1, 2, 5, 1)])
    cfg = ngc.FlowPolicyConfig(
        routing=ngc.ShortestPaths(
            allocation=A.WEIGHTED,
            weight_basis=ngc.WeightBasis.EXPLICIT,
            edge_weights=[1],
            multipath=False,
            tie_break=TIE.HIGHER_CAPACITY,
        )
    )
    with pytest.raises(ValueError, match=exact(POOLED)):
        ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, cfg)
    with pytest.raises(ValueError, match=exact("edge_weights: expected 2 entries")):
        ngc.FlowPolicy(ngc.FlowLedger(g), [0], 2, cfg)


def test_cache_rule_place():
    # The cache refusal comes after the matrix's other checks and before the tunnel count.
    g = build(3, EXAMPLE)
    cache = ngc.ShortestPathCache(g)
    for routing, text in (
        (
            ngc.DistributedTE(ngc.TunnelCount(0)),
            "DistributedTE does not use ShortestPathCache",
        ),
        (
            ngc.CentralizedTE(
                source_allocation=ngc.SourceAllocation.EQUAL, prefer_headroom=True
            ),
            "prefer_headroom requires tunnel_count",
        ),
    ):
        cfg = ngc.FlowPolicyConfig(routing=routing)
        with pytest.raises(ValueError, match=exact(text)):
            ngc.FlowPolicy(ngc.FlowLedger(g), 0, 2, cfg, None, cache)
        with pytest.raises(ValueError, match=exact(text)):
            ngc.evaluate_demand(g, 0, 2, ngc.EvaluationOptions(cfg, cache=cache))


def test_ownership_first():
    # The endpoints and the ownership of the inputs come before the matrix.
    g = build(3, EXAMPLE)
    other = build(3, EXAMPLE)
    for support in (dag_of(other, [1]), dag_of(other, [])):
        routing = ngc.CentralizedTE(
            source_allocation=ngc.SourceAllocation.EQUAL,
            support=support,
            tunnel_count=ngc.TunnelCount(2),
            prefer_headroom=True,
        )
        cfg = ngc.FlowPolicyConfig(routing=routing, admission=ngc.Admission.LOSSY)
        with pytest.raises(ValueError, match=exact("support belongs to another graph")):
            ngc.FlowPolicy(ngc.FlowLedger(g), 0, 2, cfg)
    bad = ngc.FlowPolicyConfig(
        routing=ngc.ShortestPaths(source_allocation=ngc.SourceAllocation.OPTIMIZED)
    )
    with pytest.raises(ValueError, match=exact("sources: the node set is empty")):
        ngc.FlowPolicy(
            ngc.FlowLedger(g), [], 2, bad, other.eligibility(excluded_edges=[0])
        )
    with pytest.raises(
        ValueError, match=exact("foreign Eligibility: built for another graph")
    ):
        ngc.FlowPolicy(ngc.FlowLedger(g), 0, 2, bad, other.eligibility())
