"""Every property, repr and own hash of one value of every public class under Python allocation
failure (the private hook _fail_python_allocation): the k-th Python object allocation of the read
fails, for every k in turn, and the read raises MemoryError with no value, or returns the value of
a read without failure. A full collection (gc.collect()) before each attempt empties CPython's
free lists, so an int or a float a read creates is a new object, never a reused one. A getter that
converts a value through pybind11's own machinery raises another class (TypeError from a return
conversion, RuntimeError from its list and tuple constructors) or crashes here instead.

The integers CPython caches (-5 to 256) are never created anew, so a read of a small int cannot
fail: the number sweep reads wide_values(), where every property that holds an int or a float
(by the stub's annotation) holds an uncached one, large node and edge ids, entry ids,
counts and amounts, and every property that holds a Fraction holds one with both components
uncached, one with only the numerator uncached and one with only the denominator uncached.

Each value's reads run in a child interpreter, so a crash fails this test instead of the session:
python -B test_getter_allocation.py NAME sweeps every read of one value and prints its count;
python -B test_getter_allocation.py --wide NAME sweeps the number reads (and the own hash) of one
wide value."""

from __future__ import annotations

import ast
import gc
import os
import re
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pytest
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


def wide_values():
    """Values whose int and float properties hold uncached objects: node ids from 296 up, edge ids
    from 260 up (260 one-unit filler edges come first), capacities of 2^40 units, costs of 2^30, 300
    tunnels per head (entry ids and local indexes up to 299), policy ids beyond 300, cache counts
    beyond 256; and whose Fraction properties hold each shape of uncached components."""
    cap, cost = 2**40, 2**30
    rows = [(i, i + 1, 1, 1) for i in range(260)]
    rows += [
        (296, 299, cap, cost),
        (297, 298, cap, cost),
        (298, 299, cap // 2, cost),
        (297, 299, cap, 3 * cost),
    ]
    g = build(300, rows)
    lossy = ngc.Admission.LOSSY
    shortest = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths())
    out = {}
    # Source 296 black-holed: its share is lost at the source; 297's share drops on 298 -> 299.
    sp = ngc.FlowPolicy(
        ngc.FlowLedger(g),
        [296, 297],
        299,
        ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(), admission=lossy),
    )
    sp.replace(2**41)
    dark = sp.forward(2**41, black_holes=g.eligibility(excluded_nodes=[296]))
    out["PlacementResult.dark"] = dark
    out["SourceFlow.dark"] = dark.by_source[0]
    out["SourceFlow.drops"] = dark.by_source[1]
    out["FlowRecord.lossy"] = sp.records()[-1]
    lossless = ngc.FlowPolicy(ngc.FlowLedger(g), 297, 299, shortest).place(2**42)
    out["PlacementResult.lossless"] = lossless
    out["SourceFlow.lossless"] = lossless.by_source[0]
    controller = ngc.FlowPolicyConfig(
        routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED)
    )
    exact = ngc.FlowPolicy(ngc.FlowLedger(g), 297, 299, controller).place(float("inf"))
    out["PlacementResult.exact"] = exact
    out["Tunnel.exact"] = exact.program.tunnels[0]
    failures = ngc.evaluate_edge_failures(
        g, 297, 299, ngc.EvaluationOptions(controller), 1
    )
    out["SensitivityEntry"] = failures.entries[-1]
    ledger = ngc.FlowLedger(g)
    for _ in range(300):
        ngc.FlowPolicy(ledger, 297, 299, shortest)
    t1 = ngc.FlowPolicy(
        ledger,
        297,
        299,
        ngc.FlowPolicyConfig(
            routing=ngc.DistributedTE(ngc.TunnelCount(300)), admission=lossy
        ),
    )
    t1_result = t1.replace(2**41, design_demand=300 * cost)
    out["PlacementResult.T1"] = t1_result
    out["TunnelPlan.T1"] = t1_result.plan
    out["Tunnel.T1"] = t1_result.program.tunnels[-1]
    out["TunnelFlow.T1"] = t1_result.by_tunnel[-1]
    out["FlowRecord.T1"] = t1.records()[-1]
    t2 = ngc.CentralizedTE(
        source_allocation=ngc.SourceAllocation.EQUAL,
        tunnel_count=ngc.TunnelCount(300),
        prefer_headroom=True,
    )
    t2_result = ngc.FlowPolicy(
        ngc.FlowLedger(g), 297, 299, ngc.FlowPolicyConfig(routing=t2)
    ).replace(2**40)
    out["TunnelPlan.T2"] = t2_result.plan  # fluid optimum 2^37/25
    out["Tunnel.T2"] = t2_result.program.tunnels[-1]
    # One link of capacity c split over 300 tunnels: the fluid optimum c/300.
    for name, c in (("TunnelPlan.T2.thin", 1), ("TunnelPlan.T2.wide", cap + 1)):
        link = build(2, [(0, 1, c, 1)])
        out[name] = (
            ngc.FlowPolicy(ngc.FlowLedger(link), 0, 1, ngc.FlowPolicyConfig(routing=t2))
            .replace(1)
            .plan
        )
    path = ngc.ksp(g, 297, 299, ngc.KspOptions(1))[0][1]
    out["RouteDefinition"] = ngc.RouteDefinition(297, path, weight=cap)
    cache = ngc.ShortestPathCache(g)
    cached = ngc.FlowLedger(g)
    for _ in range(2):
        for s in range(260):
            ngc.FlowPolicy(cached, s, s + 1, shortest, cache=cache)
    out["ShortestPathCacheStats"] = cache.stats()
    out["SpfResult"] = ngc.spf(
        g, 297, ngc.SpfOptions(ngc.Direction.FORWARD, targets=[299])
    )
    out["AbsoluteCost"] = ngc.AbsoluteCost(cap)
    for name, factor in (
        ("CostStretch", Fraction(2**62 + 1, 2**62 - 1)),
        ("CostStretch.numerator", Fraction(2**62 + 1, 3)),
        ("CostStretch.denominator", Fraction(3, 2**62 - 1)),
    ):
        out[name] = ngc.CostStretch(factor)
    out["TunnelCount"] = ngc.TunnelCount(300)
    out["DistributedTE"] = ngc.DistributedTE(
        ngc.TunnelCount(300), min_unreserved_bandwidth=cap
    )
    out["KspOptions"] = ngc.KspOptions(300)
    out["EvaluationOptions"] = ngc.EvaluationOptions(
        controller, offer=cap, design_demand=cap
    )
    out["SpfOptions"] = ngc.SpfOptions(ngc.Direction.FORWARD, min_residual=cap)
    out["Place"] = ngc.Place(sp, cap)
    out["Replace"] = ngc.Replace(sp, cap, design_demand=cap)
    out["Forward"] = ngc.Forward(sp, cap)
    return out


def number(x):
    return isinstance(x, (int, float, Fraction)) and not isinstance(x, bool)


def uncached(x):
    """An int or a float CPython does not take from its small-int cache."""
    return isinstance(x, float) or (
        isinstance(x, int) and not isinstance(x, bool) and not -5 <= x <= 256
    )


def shape(x):
    """A number's uncached components: (numerator, denominator) for a Fraction."""
    if isinstance(x, Fraction):
        return (uncached(x.numerator), uncached(x.denominator))
    return (uncached(x),)


def number_properties():
    """(class, property) of every property the stub annotates as an int, a float or a Fraction
    (with None), but an enumeration's value (a small int)."""
    tree = ast.parse((Path(ngc.__file__).parent / "__init__.pyi").read_text())
    out = set()
    for node in tree.body:
        if not isinstance(node, ast.ClassDef) or hasattr(
            getattr(ngc, node.name), "__members__"
        ):
            continue
        for f in node.body:
            if isinstance(f, ast.FunctionDef) and any(
                isinstance(d, ast.Name) and d.id == "property" for d in f.decorator_list
            ):
                parts = {p.strip() for p in ast.unparse(f.returns).split("|")}
                numbers = {"int", "float", "fractions.Fraction"}
                if parts & numbers and parts <= numbers | {"None"}:
                    out.add((node.name, f.name))
    return out


def properties(value):
    cls = type(value)
    return sorted(
        n
        for n in dir(cls)
        if not n.startswith("_") and isinstance(getattr(cls, n), property)
    )


def reads(value):
    """The properties, the repr and, for a class with its own hash (the value types), the hash."""
    own_hash = vars(type(value)).get("__hash__") is not None
    return properties(value) + ["__repr__"] + (["__hash__"] if own_hash else [])


def shown(value):
    # A bound value without its own repr prints its address, which a fresh wrapper changes.
    return re.sub(r" at 0x[0-9a-f]+", "", repr(value))


def sweep(name, wide=False):
    """Every read of values()[name] (wide: the number reads and the own hash of wide_values()[name])
    with each of its Python allocations failing in turn, the free lists emptied first; returns the
    failures injected."""
    import _netgraph_core as ext

    value = (wide_values if wide else values)()[name]
    names = reads(value)
    if wide:
        names = [a for a in properties(value) if number(getattr(value, a))]
        names += [a for a in reads(value) if a == "__hash__"]
    injected = 0
    for attr in names:
        if attr == "__repr__":

            def read():
                return repr(value)

        elif attr == "__hash__":

            def read():
                return hash(value)

        else:

            def read(attr=attr):
                return getattr(value, attr)

        read()  # the interpreter's first-use allocations belong to no binding
        want = shown(read())
        for k in range(600):
            gc.collect()
            hit, _, got, error = ext._fail_python_allocation(read, (), k)
            if not hit:
                assert error is None and shown(got) == want, (name, attr, k, error)
                break
            injected += 1
            assert got is None and type(error) is MemoryError, (
                name,
                attr,
                k,
                repr(error),
            )
        else:
            raise AssertionError((name, attr, "no completing read in 600 allocations"))
    return injected


if __name__ == "__main__":
    if sys.argv[1] == "--wide":
        print(sweep(sys.argv[2], wide=True))
    else:
        print(sweep(sys.argv[1]))
    sys.exit(0)


def test_wide_values_hold_uncached_numbers():
    # Every int, float or Fraction property of the stub is read on some wide value, and holds an
    # uncached object there (else its sweep could not fail where it converts); a Fraction property
    # holds each shape of uncached components.
    seen = {}
    for value in wide_values().values():
        for attr in properties(value):
            x = getattr(value, attr)
            if number(x):
                seen.setdefault((type(value).__name__, attr), set()).add(shape(x))
    want = number_properties()
    assert len(want) > 50
    assert sorted(want - set(seen)) == []
    fractions = {
        ("CostStretch", "factor"),
        ("TunnelPlan", "max_fractional_bandwidth_per_tunnel"),
    }
    assert fractions <= want
    every = {(True, True), (True, False), (False, True)}
    assert sorted(k for k in fractions if not every <= seen[k]) == []
    assert sorted(k for k in want - fractions if (True,) not in seen[k]) == []


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


@needs_hooks
def test_every_number_read_raises_memory_error():
    failed = []
    for name in wide_values():
        done = subprocess.run(
            [sys.executable, "-B", os.path.abspath(__file__), "--wide", name],
            capture_output=True,
            text=True,
            timeout=300,
        )
        if done.returncode != 0:
            failed.append((name, done.returncode, done.stderr[-300:]))
    assert failed == []


@needs_hooks
def test_a_failed_integer_read_raises_memory_error():
    # A valid result's int field of 2^60 units: the read's allocations fail one at a time after a
    # full collection; each failure is a MemoryError with the ledger and the policy unchanged, and
    # the read then returns the value.
    import _netgraph_core as ext

    g = build(2, [(0, 1, 2**61, 1)])
    ledger = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(ledger, 0, 1, ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()))
    r = p.place(2**60)

    def read():
        return r.call_delivered

    read()
    faults = 0
    for k in range(20):
        gc.collect()
        hit, _, got, error = ext._fail_python_allocation(read, (), k)
        assert p.delivered() == 2**60
        assert ledger.carried().tolist() == [2**60]
        if not hit:
            assert got == 2**60 and error is None
            assert faults >= 1  # the int's allocation (2^60 is never cached)
            return
        faults += 1
        assert got is None and type(error) is MemoryError, (k, repr(error))
    pytest.fail("no completing read in 20 allocation points")
