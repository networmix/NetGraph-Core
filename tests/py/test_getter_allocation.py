"""Every property and repr of one value of every public class under Python allocation failure (the
private hook _fail_python_allocation): the k-th Python object allocation of the read fails, for
every k in turn, and the read either raises MemoryError or RuntimeError with no value, or returns
the value of a read without failure. A getter that casts a bound value through pybind11's unchecked
allocation crashes here instead.

Each value's reads run in a child interpreter, so a crash fails this test instead of the session:
python -B test_getter_allocation.py NAME sweeps every read of one value and prints its count."""

from __future__ import annotations

import os
import re
import subprocess
import sys

if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from helpers import build, dag_of
from private_hooks import needs_hooks

import netgraph_core as ngc


def values():
    """One value of every public class with readable fields, by name."""
    inf = float("inf")
    g = build(4, [(0, 1, 8, 1), (0, 2, 8, 1), (1, 3, 8, 1), (2, 3, 8, 1)])
    mask = g.eligibility(excluded_edges=[3])
    search = ngc.spf(g, 0, ngc.SpfOptions(ngc.Direction.FORWARD))
    path = ngc.ksp(g, 0, 3, ngc.KspOptions(1))[0][1]
    wcmp = ngc.ShortestPaths(
        allocation=ngc.RouteAllocation.WEIGHTED,
        weight_basis=ngc.WeightBasis.EXPLICIT,
        edge_weights=[1, 2, 1, 1],
    )
    dte = ngc.DistributedTE(ngc.TunnelCount(2), cost_ceiling=ngc.CostStretch(1.5))
    t2 = ngc.CentralizedTE(
        source_allocation=ngc.SourceAllocation.EQUAL,
        support=dag_of(g, [0, 1, 2, 3]),
        tunnel_count=ngc.TunnelCount(2),
        prefer_headroom=True,
    )
    controller = ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED)
    ledger = ngc.FlowLedger(g)
    cfg = ngc.FlowPolicyConfig(routing=wcmp)
    policy = ngc.FlowPolicy(ledger, 0, 3, cfg)
    result = policy.place(
        4,
        ngc.OutputOptions(
            edge_flows=True, residual=True, saturated_boundary_edges=True
        ),
    )
    t1_result = ngc.FlowPolicy(
        ngc.FlowLedger(g),
        0,
        3,
        ngc.FlowPolicyConfig(routing=dte, admission=ngc.Admission.LOSSY),
    ).replace(6, design_demand=6, outputs=ngc.OutputOptions(edge_reservations=True))
    t2_result = ngc.FlowPolicy(
        ngc.FlowLedger(g), 0, 3, ngc.FlowPolicyConfig(routing=t2)
    ).replace(inf)
    exact = ngc.FlowPolicy(
        ngc.FlowLedger(g), 0, 3, ngc.FlowPolicyConfig(routing=controller)
    ).place(inf, ngc.OutputOptions(min_cut=True))
    failures = ngc.evaluate_edge_failures(
        g, 0, 3, ngc.EvaluationOptions(ngc.FlowPolicyConfig(routing=controller)), 1
    )
    cache = ngc.ShortestPathCache(g)
    return {
        "Graph": g,
        "Eligibility": mask,
        "SpfResult": search,
        "PredecessorDAG": search.dag,
        "Path": path,
        "ShortestPaths": wcmp,
        "DistributedTE": dte,
        "CentralizedTE": t2,
        "Pinned": ngc.Pinned(t1_result.program, black_holes=mask),
        "TunnelCount": ngc.TunnelCount(2, up_to=True),
        "AbsoluteCost": ngc.AbsoluteCost(3),
        "CostStretch": ngc.CostStretch(1.2),
        "FlowPolicyConfig": cfg,
        "OutputOptions": ngc.OutputOptions(min_cut=True),
        "EvaluationOptions": ngc.EvaluationOptions(cfg, eligibility=mask, offer=5),
        "SpfOptions": ngc.SpfOptions(
            ngc.Direction.FORWARD,
            targets=[3],
            tie_break=ngc.RouteTieBreak.HIGHER_CAPACITY,
            multipath=False,
        ),
        "KspOptions": ngc.KspOptions(2, cost_ceiling=ngc.CostStretch(2)),
        "FlowLedger": ledger,
        "FlowRecord": ledger.records()[0],
        "FlowPolicy": policy,
        "PlacementResult": result,
        "SourceFlow": result.by_source[0],
        "PlacementResult.T1": t1_result,
        "TunnelPlan.T1": t1_result.plan,
        "Tunnel": t1_result.program.tunnels[0],
        "TunnelFlow": t1_result.by_tunnel[0],
        "PlacementResult.T2": t2_result,
        "TunnelPlan.T2": t2_result.plan,
        "PlacementResult.exact": exact,
        "ForwardingProgram": ngc.ForwardingProgram.from_dag(
            g, [0], [3], dag_of(g, [0, 1, 2, 3])
        ),
        "RouteDefinition": ngc.RouteDefinition(0, path, weight=2),
        "SensitivityResult": failures,
        "SensitivityEntry": failures.entries[0],
        "ShortestPathCache": cache,
        "ShortestPathCache.Stats": cache.stats(),
    }


def reads(value):
    cls = type(value)
    names = sorted(
        n
        for n in dir(cls)
        if not n.startswith("_") and isinstance(getattr(cls, n), property)
    )
    return names + ["__repr__"]


def shown(value):
    # A bound value without its own repr prints its address, which a fresh wrapper changes.
    return re.sub(r" at 0x[0-9a-f]+", "", repr(value))


def sweep(name):
    """Every read of values()[name] with each of its Python allocations failing in turn; returns
    the failures injected."""
    import _netgraph_core as ext

    value = values()[name]
    injected = 0
    for attr in reads(value):
        if attr == "__repr__":

            def read():
                return repr(value)

        else:

            def read(attr=attr):
                return getattr(value, attr)

        read()  # the interpreter's first-use allocations belong to no binding
        want = shown(read())
        for k in range(600):
            hit, _, got, error = ext._fail_python_allocation(read, (), k)
            if not hit:
                assert error is None and shown(got) == want, (name, attr, k, error)
                break
            injected += 1
            assert got is None and isinstance(error, (MemoryError, RuntimeError)), (
                name,
                attr,
                k,
                repr(error),
            )
        else:
            raise AssertionError((name, attr, "no completing read in 600 allocations"))
    return injected


if __name__ == "__main__":
    print(sweep(sys.argv[1]))
    sys.exit(0)


@needs_hooks
def test_every_read_survives_python_allocation_failures():
    failed = []
    for name in values():
        done = subprocess.run(
            [sys.executable, "-B", os.path.abspath(__file__), name],
            capture_output=True,
            text=True,
            timeout=300,
        )
        if done.returncode != 0 or int(done.stdout.split()[-1]) < 1:
            failed.append((name, done.returncode, done.stderr[-300:]))
    assert failed == []
