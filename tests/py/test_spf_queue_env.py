"""NGRAPH_CORE_SPF_QUEUE selects the SPF frontier queue; both routes must agree."""

from __future__ import annotations

import os
import subprocess
import sys

SCRIPT = r"""
import hashlib
import numpy as np
import netgraph_core as ngc

rng = np.random.default_rng(7)
n = 300
src, dst, cost, cap = [], [], [], []
for u in range(n):
    for k in range(6):
        src.append(u); dst.append(int(rng.integers(0, n)))
        cost.append(0 if rng.integers(0, 5) == 0 else int(rng.integers(1, 32)))
        cap.append(float(1 + rng.integers(0, 100) / 7.0))
g = ngc.StrictMultiDiGraph.from_arrays(
    n,
    src=np.array(src, dtype=np.int32), dst=np.array(dst, dtype=np.int32),
    capacity=np.array(cap, dtype=np.float64), cost=np.array(cost, dtype=np.int64),
    ext_edge_ids=np.arange(len(src), dtype=np.int64),
)
algs = ngc.Algorithms(ngc.Backend.cpu())
gh = algs.build_graph(g)
h = hashlib.sha256()
residual = np.array(cap, dtype=np.float64); residual[::13] = 0.0
for s in (0, 17, 150):
    for dstn in (None, 5, 299):
        for multipath in (True, False):
            for me in (True, False):
                sel = ngc.EdgeSelection(multi_edge=me, require_capacity=not me, tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL)
                dist, dag = algs.spf(gh, s, dstn, selection=sel, residual=residual, multipath=multipath)
                for arr in (np.asarray(dist), np.asarray(dag.parent_offsets), np.asarray(dag.parents), np.asarray(dag.via_edges)):
                    h.update(np.ascontiguousarray(arr).tobytes())
    dist, dag = algs.spf_to(gh, s, multipath=True)
    for arr in (np.asarray(dist), np.asarray(dag.parent_offsets), np.asarray(dag.parents), np.asarray(dag.via_edges)):
        h.update(np.ascontiguousarray(arr).tobytes())
flow, summary = algs.max_flow(gh, 0, 299, with_edge_flows=True)
h.update(repr((flow, list(summary.costs), list(summary.flows))).encode())
print(h.hexdigest())
"""


def _run(queue: str | None) -> str:
    env = dict(os.environ)
    env.pop("NGRAPH_CORE_SPF_QUEUE", None)
    if queue is not None:
        env["NGRAPH_CORE_SPF_QUEUE"] = queue
    out = subprocess.run(
        [sys.executable, "-c", SCRIPT],
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    return out.stdout.strip()


def test_queue_routes_agree_bit_for_bit() -> None:
    heap = _run("heap")
    bucket = _run("bucket")
    auto = _run(None)
    assert len(heap) == 64
    assert heap == bucket == auto
