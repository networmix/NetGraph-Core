import os
import subprocess
from pathlib import Path

r = Path(__file__).resolve().parent
out = r / "results"
b = Path("build/metal_pathfinding")
phases = [
    ("old-b1", "old", False),
    ("new-b1", "baseline", False),
    ("settled-b", "settled", False),
    ("inline", "inline", False),
    ("append", "append", False),
    ("frontpop", "frontpop", True),
    ("heap-b", "baseline", True),
    ("new-b2", "baseline", False),
    ("old-b2", "old", False),
]
for name, v, heap in phases:
    with (out / "load-second.txt").open("a") as f:
        f.write(name + " " + subprocess.check_output(["uptime"], text=True))
    env = os.environ.copy()
    env.pop("NGRAPH_CORE_SPF_QUEUE", None)
    if heap:
        env["NGRAPH_CORE_SPF_QUEUE"] = "heap"
    with (out / (name + ".csv")).open("w") as f:
        subprocess.run(
            [str(b / ("review-" + v)), "focus", "7"], env=env, stdout=f, check=True
        )
    print(name, flush=True)
for name, v in [
    ("old1", "old"),
    ("new1", "baseline"),
    ("gcd", "gcd"),
    ("new2", "baseline"),
    ("old2", "old"),
]:
    with (out / "load-second.txt").open("a") as f:
        f.write("build-" + name + " " + subprocess.check_output(["uptime"], text=True))
    with (out / ("build-" + name + ".csv")).open("w") as f:
        subprocess.run([str(b / ("review-build-" + v))], stdout=f, check=True)
for mode in ["reverse", "boundary"]:
    for name, v, heap in [
        ("old1", "old", False),
        ("new1", "baseline", False),
        ("heap", "baseline", True),
        ("settled", "settled", False),
        ("old2", "old", False),
    ]:
        with (out / "load-second.txt").open("a") as f:
            f.write(
                mode + "-" + name + " " + subprocess.check_output(["uptime"], text=True)
            )
        env = os.environ.copy()
        env.pop("NGRAPH_CORE_SPF_QUEUE", None)
        if heap:
            env["NGRAPH_CORE_SPF_QUEUE"] = "heap"
        with (out / (mode + "-" + name + ".csv")).open("w") as f:
            subprocess.run(
                [str(b / ("review-" + v)), mode, "7"], env=env, stdout=f, check=True
            )
        print(mode, name, flush=True)
