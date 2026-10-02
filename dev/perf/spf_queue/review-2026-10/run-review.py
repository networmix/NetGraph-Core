import hashlib
import os
import subprocess
from pathlib import Path

r = Path(__file__).resolve().parent
out = r / "results"
b = Path("build/metal_pathfinding")
with (out / "environment.txt").open("w") as f:
    for cmd in [
        ["date", "-u"],
        ["uname", "-a"],
        ["sysctl", "-n", "machdep.cpu.brand_string"],
        ["git", "rev-parse", "HEAD"],
        ["clang++", "--version"],
    ]:
        f.write(subprocess.check_output(cmd, text=True))
    for p in list(b.glob("review-*")) + [
        b / "cpu-queue",
        b / "cpu-queue-old",
        Path("src/shortest_paths.cpp"),
    ]:
        if p.is_file():
            f.write(hashlib.sha256(p.read_bytes()).hexdigest() + " " + str(p) + "\n")
phases = [
    ("old-a1", "old", False),
    ("new-a1", "baseline", False),
    ("heap-a1", "baseline", True),
    ("settled", "settled", False),
    ("settled-heap", "settled", True),
    ("singleton", "singleton", False),
    ("heapitem", "heapitem", True),
    ("notouch", "notouch", True),
    ("new-a2", "baseline", False),
    ("heap-a2", "baseline", True),
    ("old-a2", "old", False),
]
for name, v, heap in phases:
    with (out / "load.txt").open("a") as f:
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
with (out / "load.txt").open("a") as f:
    f.write("end " + subprocess.check_output(["uptime"], text=True))
