import subprocess
from pathlib import Path

r = Path(__file__).resolve().parent
out = r / "results"
b = Path("build/metal_pathfinding")
for name, v in [
    ("old-c1", "old"),
    ("new-c1", "baseline"),
    ("append-c", "append"),
    ("combined", "combined"),
    ("inline-c", "inline"),
    ("new-c2", "baseline"),
    ("old-c2", "old"),
]:
    with (out / "load-third.txt").open("a") as f:
        f.write(name + " " + subprocess.check_output(["uptime"], text=True))
    with (out / (name + ".csv")).open("w") as f:
        subprocess.run([str(b / ("review-" + v)), "focus", "7"], stdout=f, check=True)
    print(name, flush=True)
for name, v in [("new1", "baseline"), ("gcd", "gcd"), ("new2", "baseline")]:
    with (out / "load-third.txt").open("a") as f:
        f.write(
            "build-repeat-"
            + name
            + " "
            + subprocess.check_output(["uptime"], text=True)
        )
    with (out / ("build-repeat-" + name + ".csv")).open("w") as f:
        subprocess.run([str(b / ("review-build-" + v))], stdout=f, check=True)
