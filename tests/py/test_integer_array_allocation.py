"""Integer-array casts under Python allocation failure, each reader/dtype in its own process.

Constructors use an instance shell created before the fault window, so the sweep covers input
conversion rather than pybind11's documented unchecked allocation of a new bound instance.
The private hook counts only allocations inside NumPy casts, excluding result wrapping.
"""

from __future__ import annotations

import functools
import gc
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

if __name__ == "__main__":
    spec = importlib.util.spec_from_file_location("_netgraph_core", sys.argv[3])
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["_netgraph_core"] = module
    spec.loader.exec_module(module)

import _netgraph_core as ext
import numpy as np
from helpers import full_ledger_state, plain
from private_hooks import needs_hooks

import netgraph_core as ngc

KINDS = ("int32", "uint64", "swapped_uint64")
# Every public reader of wide_integers, with separate rows for each array argument. Deferred
# options are exercised at their consuming call; nested hops and pair endpoints are included.
READERS = (
    "Graph.sources",
    "Graph.destinations",
    "Graph.capacities",
    "Graph.costs",
    "eligibility.nodes",
    "eligibility.edges",
    "from_csr.offsets",
    "from_csr.parents",
    "from_csr.edges",
    "Path.hop",
    "ShortestPaths.weights",
    "FlowLedger.service",
    "FlowLedger.reservation",
    "reset.service",
    "reset.reservation",
    "saturated_boundary_edges.sources",
    "FlowPolicy.sources",
    "FlowPolicy.destinations",
    "spf.roots",
    "spf.targets",
    "spf.residual",
    "ksp.residual",
    "from_routes.sources",
    "from_routes.destinations",
    "from_routes.weights",
    "from_routes.residual",
    "from_dag.sources",
    "from_dag.destinations",
    "from_dag.weights",
    "from_dag.residual",
    "evaluate_demand.sources",
    "evaluate_demand.destinations",
    "evaluate_pairs.sources",
    "evaluate_pairs.destinations",
    "evaluate_edge_failures.sources",
    "evaluate_edge_failures.destinations",
    "as_float",
    "to_units",
)


def make(reader, kind):
    arrays = []

    def array(values):
        dtype = (
            np.dtype("uint64").newbyteorder("S")
            if kind == "swapped_uint64"
            else np.dtype(kind)
        )
        a = np.array(values * 1000 if len(values) == 1 else values, dtype=dtype)
        arrays.append(a)
        return a

    large = reader.startswith(
        (
            "Graph.",
            "from_csr.",
            "Path.",
            "ShortestPaths.",
            "FlowLedger.",
            "reset.",
            "from_routes.",
            "from_dag.",
        )
    ) or reader in ("spf.residual", "ksp.residual")
    n = 1001 if reader.startswith("from_csr.") else 2
    edges = 1000 if large else 1
    g = ngc.Graph(n, [0] * edges, [1] * edges, [10] * edges, [1] * edges, 1.0)
    ledger = ngc.FlowLedger(g)
    cfg = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths())
    options = ngc.EvaluationOptions(cfg, offer=4)
    dag = ngc.Path(g, [[0]]).to_dag()
    fn, args = None, ()

    def check():
        return None

    def initialize(cls, args, **kwargs):
        obj = cls.__new__(cls)
        return functools.partial(obj.__init__, **kwargs), args, obj

    if reader.startswith("Graph."):
        values = [[0] * 1000, [1] * 1000, [10] * 1000, [1] * 1000]
        position = ("sources", "destinations", "capacities", "costs").index(
            reader.split(".")[1]
        )
        values[position] = array(values[position])
        fn, args, obj = initialize(ngc.Graph, (2, *values, 1.0))

        def check():
            return (
                obj.edge_sources().tolist(),
                obj.edge_destinations().tolist(),
                obj.edge_capacity_units().tolist(),
                obj.edge_costs().tolist(),
            )
    elif reader.startswith("eligibility."):
        fn = g.eligibility
        args = (array([0]), []) if reader.endswith("nodes") else ([], array([0]))
    elif reader.startswith("from_csr."):
        values = [[0, 0] + [1000] * 1000, [0] * 1000, list(range(1000))]
        position = ("offsets", "parents", "edges").index(reader.split(".")[1])
        values[position] = array(values[position])
        fn, args = ngc.PredecessorDAG.from_csr, (g, *values)
    elif reader == "Path.hop":
        fn, args, obj = initialize(ngc.Path, (g, [array(list(range(1000)))]))

        def check():
            return plain(obj)
    elif reader.startswith("ShortestPaths."):
        weights = array([1])
        fn, args, obj = initialize(
            ngc.ShortestPaths,
            (),
            allocation=ngc.RouteAllocation.WEIGHTED,
            weight_basis=ngc.WeightBasis.EXPLICIT,
            edge_weights=weights,
        )

        def check():
            return plain(obj.edge_weights)
    elif reader.startswith("FlowLedger."):
        if reader.endswith("service"):
            fn, args, obj = initialize(ngc.FlowLedger, (g, array([8])))
        else:
            fn, args, obj = initialize(
                ngc.FlowLedger, (g,), reservation_baseline=array([8])
            )

        def check():
            return (obj.residual().tolist(), obj.unreserved().tolist())
    elif reader.startswith("reset."):
        fn = (
            ledger.reset
            if reader.endswith("service")
            else functools.partial(ledger.reset, reservation_baseline=array([8]))
        )
        args = (array([8]),) if reader.endswith("service") else ()
    elif reader == "saturated_boundary_edges.sources":
        fn, args = ledger.saturated_boundary_edges, (array([0]),)
    elif reader.startswith("FlowPolicy."):
        sources = array([0]) if reader.endswith("sources") else [0]
        destinations = array([1]) if reader.endswith("destinations") else [1]
        fn, args, obj = initialize(ngc.FlowPolicy, (ledger, sources, destinations, cfg))

        def check():
            return obj.delivered()
    elif reader.startswith("spf."):
        roots = array([0]) if reader.endswith("roots") else [0]
        opts = ngc.SpfOptions(
            ngc.Direction.FORWARD,
            targets=array([1]) if reader.endswith("targets") else None,
            residual=array([8]) if reader.endswith("residual") else None,
        )
        fn, args = ngc.spf, (g, roots, opts)
    elif reader == "ksp.residual":
        fn, args = ngc.ksp, (g, 0, 1, ngc.KspOptions(1, residual=array([8])))
    elif reader.startswith(("from_routes.", "from_dag.")):
        factory, field = reader.split(".")
        sources = array([0]) if field == "sources" else [0]
        destinations = array([1]) if field == "destinations" else [1]
        kwargs = {}
        if field == "weights":
            kwargs = dict(
                allocation=ngc.RouteAllocation.WEIGHTED,
                weight_basis=ngc.WeightBasis.EXPLICIT,
                edge_weights=array([1]),
            )
        elif field == "residual":
            kwargs = dict(
                allocation=ngc.RouteAllocation.WEIGHTED,
                weight_basis=ngc.WeightBasis.RESIDUAL_SNAPSHOT,
                residual=array([8]),
            )
        route = [ngc.RouteDefinition(0, dag)] if factory == "from_routes" else dag
        fn = functools.partial(getattr(ngc.ForwardingProgram, factory), **kwargs)
        args = (g, sources, destinations, route)
    elif reader.startswith("evaluate_"):
        entry, field = reader.split(".")
        sources = array([0]) if field == "sources" else [0]
        destinations = array([1]) if field == "destinations" else [1]
        fn = getattr(ngc, entry)
        args = (
            (g, [(sources, destinations)], options)
            if entry == "evaluate_pairs"
            else (g, sources, destinations, options)
        )
    elif reader in ("as_float", "to_units"):
        fn, args = getattr(ngc, reader), (array(list(range(1000))), 1.0)
    else:
        raise AssertionError(reader)

    def snapshot():
        return [a.tobytes() for a in arrays], full_ledger_state(ledger)

    return fn, args, snapshot, check


def sweep(reader, kind):
    fn, args, snapshot, check = make(reader, kind)
    expected = plain(fn(*args))
    expected_check = check()
    expected_after = snapshot()
    failures = 0
    recovered = 0
    for point in range(2000):
        fn, args, snapshot, check = make(reader, kind)
        before = snapshot()
        gc.collect()
        hit, seen, result, error = ext._fail_python_allocation(fn, args, point, True)
        if not hit:
            assert error is None, (point, repr(error))
            assert plain(result) == expected, point
            assert check() == expected_check, point
            assert snapshot() == expected_after, point
            assert seen == failures, (seen, failures)
            assert failures > 0 or kind == "uint64"
            print(
                json.dumps(
                    dict(
                        reader=reader,
                        dtype=kind,
                        failures=failures,
                        recovered=recovered,
                        terminal="success",
                    )
                ),
                flush=True,
            )
            return
        failures += 1
        if error is None:
            # NumPy can recover a temporary-buffer allocation internally. A recovered cast
            # must still give every original result and the same complete state.
            assert plain(result) == expected, point
            assert check() == expected_check, point
            assert snapshot() == expected_after, point
            recovered += 1
        else:
            assert result is None and type(error) is MemoryError, (point, repr(error))
            assert snapshot() == before, point
        assert plain(fn(*args)) == expected, point
        assert check() == expected_check, point
        assert snapshot() == expected_after, point
    raise AssertionError("allocation sweep did not finish")


@needs_hooks
def test_integer_array_readers_propagate_allocation_failure_in_separate_processes():
    for reader in READERS:
        for kind in KINDS:
            result = subprocess.run(
                [sys.executable, "-B", __file__, reader, kind, ext.__file__],
                capture_output=True,
                text=True,
                timeout=180,
            )
            assert result.returncode == 0, (reader, kind, result.stdout, result.stderr)
            row = json.loads(result.stdout.splitlines()[-1])
            assert row["reader"] == reader and row["dtype"] == kind
            assert row["terminal"] == "success"
            assert row["failures"] > 0 or kind == "uint64"
            if reader != "to_units" and kind != "uint64":
                assert row["failures"] > row["recovered"]


if __name__ == "__main__":
    assert Path(ext.__file__).resolve() == Path(sys.argv[3]).resolve()
    sweep(sys.argv[1], sys.argv[2])
