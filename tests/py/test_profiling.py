from __future__ import annotations

import os
import re
import subprocess
import sys
import textwrap


def _run_profile_subprocess(
    code: str, *, extra_env: dict[str, str] | None = None
) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["NGRAPH_CORE_PROFILE"] = "1"
    if extra_env:
        env.update(extra_env)
    return subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )


def _calls_for(stderr: str, scope: str) -> int:
    match = re.search(rf"^{re.escape(scope)}:\n  calls: (\d+)\n", stderr, re.MULTILINE)
    assert match, stderr
    return int(match.group(1))


def test_profiling_dump_and_reset_lifecycle_via_python_api() -> None:
    code = textwrap.dedent(
        """
        import sys
        import numpy as np
        import netgraph_core as ngc

        src = np.array([0, 1], dtype=np.int32)
        dst = np.array([1, 2], dtype=np.int32)
        cap = np.array([1.0, 1.0], dtype=np.float64)
        cost = np.array([1, 1], dtype=np.int64)
        ext = np.array([0, 1], dtype=np.int64)
        graph = ngc.StrictMultiDiGraph.from_arrays(3, src, dst, cap, cost, ext)
        algs = ngc.Algorithms(ngc.Backend.cpu())
        gh = algs.build_graph(graph)

        sys.stderr.write("===FIRST===\\n")
        for _ in range(3):
            algs.spf(gh, 0, 2)
        ngc.profiling_dump()
        ngc.profiling_reset()
        sys.stderr.write("===SECOND===\\n")
        ngc.profiling_dump()
        print(ngc.profiling_enabled())
        """
    )
    proc = _run_profile_subprocess(code)
    assert proc.stdout.strip() == "True"
    first = proc.stderr.split("===FIRST===\n", 1)[1].split("===SECOND===\n", 1)[0]
    second = proc.stderr.split("===SECOND===\n", 1)[1]
    assert _calls_for(first, "shortest_paths_core") == 3
    assert "shortest_paths_core:" not in second


def test_profiling_threaded_scope_counts_match_serial() -> None:
    code = textwrap.dedent(
        """
        import concurrent.futures
        import os
        import numpy as np
        import netgraph_core as ngc

        mode = os.environ["P8_MODE"]
        src = np.array([0, 1], dtype=np.int32)
        dst = np.array([1, 2], dtype=np.int32)
        cap = np.array([1.0, 1.0], dtype=np.float64)
        cost = np.array([1, 1], dtype=np.int64)
        ext = np.array([0, 1], dtype=np.int64)
        graph = ngc.StrictMultiDiGraph.from_arrays(3, src, dst, cap, cost, ext)
        algs = ngc.Algorithms(ngc.Backend.cpu())
        gh = algs.build_graph(graph)

        def worker(iters: int) -> None:
            for _ in range(iters):
                algs.spf(gh, 0, 2)

        ngc.profiling_reset()
        if mode == "serial":
            worker(20)
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
                list(pool.map(worker, [5, 5, 5, 5]))
        ngc.profiling_dump()
        """
    )
    serial = _run_profile_subprocess(code, extra_env={"P8_MODE": "serial"})
    threaded = _run_profile_subprocess(code, extra_env={"P8_MODE": "threaded"})
    assert _calls_for(serial.stderr, "shortest_paths_core") == 20
    assert _calls_for(threaded.stderr, "shortest_paths_core") == 20
