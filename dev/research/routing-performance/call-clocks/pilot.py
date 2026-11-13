"""One P1 process at each pilot delta, plus its untimed plain export control."""

import json
from pathlib import Path
import subprocess
import sys
import time

from run_workload import ROOT, equal
from capture import REF

SCRATCH = ROOT / "scratch/clocks"
SHA = "ff7559d674b9a2253f6a3f8f51ca005b89ca63999dae8d4f669d1ae4e6dc4950"


def main(out):
    out.mkdir()
    rows = []
    expected = None
    for label, lane, delta, python, sha in [
        ("plain", "plain", 0, ROOT / "scratch/consumer/ref/venv/bin/python", REF),
        ("d200", "delay", 200000, SCRATCH / "venv-delay/bin/python", SHA),
        ("d0", "delay", 0, SCRATCH / "venv-delay/bin/python", SHA),
    ]:
        target = out / label
        command = [
            str(python),
            "-B",
            str(Path(__file__).with_name("run_workload.py")),
            "--workload",
            "dcbb96",
            "--workers",
            "8",
            "--lane",
            lane,
            "--delta-ns",
            str(delta),
            "--native-sha",
            sha,
            "--out",
            str(target),
        ]
        with (out / (label + ".log")).open("x") as stream:
            tick = time.monotonic()
            result = subprocess.run(
                command, stdout=stream, stderr=subprocess.STDOUT, timeout=180
            )
            wall = time.monotonic() - tick
            stream.write(f"\nSUMMARY exit={result.returncode} elapsed_s={wall}\n")
        assert result.returncode == 0, (label, result.returncode)
        actual = json.loads((target / "results.json").read_text())
        assert actual["metadata"]["core_extension_sha256"] == sha
        if expected is None:
            expected = actual
        else:
            equal(expected, actual, identity=True)
        summary = json.loads((target / "summary.json").read_text())
        row = {
            "label": label,
            "wall_s": wall,
            "phases": summary["phases"],
            "call_count": summary["parallel_call_count"],
            "helper_elapsed_ns": summary["helper_elapsed_ns"],
            "full_export_equal": True,
        }
        rows.append(row)
        print(json.dumps(row), flush=True)
    assert rows[1]["call_count"] == rows[2]["call_count"] > 0
    # 63 fresh processes; bound every delta by the slowest pilot and reserve
    # 30% for machine-state reads, export comparisons and scheduler variation.
    projection = 63 * max(rows[1]["wall_s"], rows[2]["wall_s"]) * 1.30
    report = {
        "exploratory": True,
        "rows": rows,
        "projected_window_s": projection,
        "projection_formula": "63 * max(pilot d0 wall, pilot d200 wall) * 1.30",
        "over_45_minutes": projection > 2700,
        "summary": "STOP: projected window exceeds 45 minutes"
        if projection > 2700
        else "PASS: pilot exports equal and projected timing window fits 45 minutes",
    }
    with (out / "summary.json").open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
