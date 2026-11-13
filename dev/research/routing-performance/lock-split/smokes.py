"""Sequential stage-2 controls; every workload here is explicitly shortened."""

from pathlib import Path
import sys

from measure import run
from reuse import ROOT


def main(out):
    out.mkdir()
    task = Path(__file__).resolve().parent
    ref = str(ROOT / "scratch/consumer/ref/venv/bin/python")
    scratch = str(ROOT / "scratch/split/venv/bin/python")
    run(
        [scratch, "-B", str(task / "controls.py"), str(out / "controls.json")],
        out / "controls.log",
    )
    for workload, workers in (
        ("dcbb96", 1),
        ("dcbb96", 8),
        ("nsfnet", 1),
        ("nsfnet", 8),
        ("backbone", 1),
        ("backbone", 28),
    ):
        label = f"{workload}-{workers}"
        run(
            [
                ref,
                "-B",
                str(task / "measure.py"),
                "--workload",
                workload,
                "--workers",
                str(workers),
                "--smoke",
                "--out",
                str(out / label),
            ],
            out / (label + ".log"),
        )
        print(f"PASS: shortened {label}; export and closure", flush=True)
    run(
        [
            ref,
            "-B",
            str(task / "capture.py"),
            "--workload",
            "dcbb96",
            "--workers",
            "8",
            "--seconds",
            "1",
            "--no-sample",
            "--smoke",
            "--expected",
            str(out / "dcbb96-8/plain"),
            "--out",
            str(out / "capture-no-sample"),
        ],
        out / "capture-no-sample.log",
    )
    print(
        "SUMMARY: PASS controls, six shortened export/closure rows and reference capture launcher without sampling",
        flush=True,
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]))
