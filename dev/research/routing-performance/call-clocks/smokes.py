"""Sequential, bounded, untimed stage-2 workload and launcher controls."""

import json
from pathlib import Path
import subprocess
import sys

from run_workload import ROOT
from capture import REF


def main(out):
    out.mkdir()
    python = ROOT / "scratch/consumer/ref/venv/bin/python"
    task = Path(__file__).parent
    commands = [
        [
            str(ROOT / "scratch/clocks/venv-delay/bin/python"),
            "-B",
            str(task / "controls.py"),
            str(out / "controls.json"),
        ],
        [
            str(python),
            "-B",
            str(task / "calibrate_parallel.py"),
            str(out / "parallel-null.json"),
            "--calls",
            "1000",
        ],
    ]
    for workload, workers, smoke in [
        ("dcbb96", 8, True),
        ("nsfnet", 8, True),
        ("nsfnet", 8, False),
        ("backbone", 8, False),
        ("backbone", 28, False),
    ]:
        label = f"{workload}-{workers}-{'short' if smoke else 'full'}"
        plain = out / (label + "-plain")
        base = [
            str(python),
            "-B",
            str(task / "run_workload.py"),
            "--workload",
            workload,
            "--workers",
            str(workers),
            "--native-sha",
            REF,
        ]
        commands.append(
            base
            + ["--lane", "plain", "--out", str(plain)]
            + (["--smoke"] if smoke else [])
        )
        if smoke:
            clocks = out / (label + "-clocks")
            commands.append(
                base
                + [
                    "--lane",
                    "clocks",
                    "--calibration-calls",
                    "1000",
                    "--out",
                    str(clocks),
                    "--expected",
                    str(plain),
                    "--smoke",
                ]
            )
            commands.append(
                [
                    str(python),
                    "-B",
                    str(task / "analyze_clocks.py"),
                    str(clocks),
                    str(out / (label + "-analysis.json")),
                ]
            )
        else:
            commands.append(
                [
                    str(python),
                    "-B",
                    str(task / "capture.py"),
                    "--workload",
                    workload,
                    "--workers",
                    str(workers),
                    "--seconds",
                    "1",
                    "--no-sample",
                    "--out",
                    str(out / (label + "-no-sample")),
                    "--expected",
                    str(plain),
                ]
            )
    for index, command in enumerate(commands):
        with (out / f"{index:02d}.log").open("x") as log:
            result = subprocess.run(
                command, stdout=log, stderr=subprocess.STDOUT, timeout=120
            )
            log.write(f"\nSUMMARY exit={result.returncode}\n")
        print(
            json.dumps({"index": index, "exit": result.returncode, "command": command}),
            flush=True,
        )
        if result.returncode:
            raise SystemExit(result.returncode)
    (out / "summary.json").write_text(
        json.dumps(
            {
                "commands": commands,
                "summary": "PASS: controls, short plain/clock export equality and closure, three full plain/no-sample capture export comparisons",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]))
