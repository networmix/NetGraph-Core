import hashlib
import os
import subprocess
from pathlib import Path

r = Path(__file__).resolve().parent
out = r / "results"
repo = Path.cwd()
bench = repo / "dev/perf/spf_queue/bench_netgraph.py"
py = "/Users/networmix/ws/NetGraph/venv/bin/python"
workloads = "spf_full_clos_200x200,spf_full_grid_100x100_weighted,maxflow_grid_100x100_weighted_combine,placement_backbone_clos"
for label, root, heap in [
    ("old1", Path("/Users/networmix/ws/NetGraph-Core-main"), False),
    ("new", repo, False),
    ("heap", repo, True),
    ("old2", Path("/Users/networmix/ws/NetGraph-Core-main"), False),
]:
    env = os.environ.copy()
    env.pop("NGRAPH_CORE_SPF_QUEUE", None)
    env.update(
        PYTHONHASHSEED="0",
        PYTHONDONTWRITEBYTECODE="1",
        PYTHONPATH=f"{root}/venv/lib/python3.14/site-packages:{root}/python",
    )
    if heap:
        env["NGRAPH_CORE_SPF_QUEUE"] = "heap"
    with (out / "netgraph-environment.txt").open("a") as f:
        f.write(label + " " + subprocess.check_output(["uptime"], text=True))
        for p in (root / "venv/lib/python3.14/site-packages").glob(
            "_netgraph_core*.so"
        ):
            f.write(hashlib.sha256(p.read_bytes()).hexdigest() + " " + str(p) + "\n")
        f.write(
            subprocess.check_output(
                [
                    py,
                    "-c",
                    "import _netgraph_core,netgraph_core;print(_netgraph_core.__file__);print(netgraph_core.__file__)",
                ],
                env=env,
                text=True,
            )
        )
    with (
        (out / ("netgraph-" + label + ".csv")).open("w") as f,
        (out / ("netgraph-" + label + ".stderr")).open("w") as e,
    ):
        subprocess.run(
            [py, str(bench), label, "5", workloads],
            env=env,
            stdout=f,
            stderr=e,
            check=True,
        )
    print("netgraph", label, flush=True)
