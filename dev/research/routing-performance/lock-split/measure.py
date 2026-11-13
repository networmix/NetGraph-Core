"""One bounded reference/export and scratch clock pair; smoke mode is untimed."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from reuse import ROOT


def run(command, log):
    with log.open("x") as stream:
        result = subprocess.run(
            command, stdout=stream, stderr=subprocess.STDOUT, timeout=240
        )
        stream.write(f"\nSUMMARY exit={result.returncode}\n")
    assert result.returncode == 0, (str(log), result.returncode)


def main(args):
    args.out.mkdir()
    task = Path(__file__).resolve().parent
    ref = ROOT / "scratch/consumer/ref/venv"
    scratch = ROOT / "scratch/split/venv"
    identity = json.loads((ROOT / "scratch/split/build-identity.json").read_text())
    run(
        [
            str(scratch / "bin/python"),
            "-B",
            str(task / "calibrate.py"),
            str(args.out / "parallel-null.json"),
            "--workers",
            str(args.workers),
            "--calls",
            "200" if args.smoke else "3000",
        ],
        args.out / "parallel-null.log",
    )
    for lane, venv, folder in (("plain", ref, "plain"), ("clocks", scratch, "clocks")):
        (extension,) = (venv / "lib/python3.14/site-packages").glob(
            "_netgraph_core*.so"
        )
        sha = hashlib.sha256(extension.read_bytes()).hexdigest()
        assert (
            sha == identity["reference_sha256" if lane == "plain" else "scratch_sha256"]
        )
        command = [
            str(venv / "bin/python"),
            "-B",
            str(task / "run_workload.py"),
            "--workload",
            args.workload,
            "--workers",
            str(args.workers),
            "--lane",
            lane,
            "--native-sha",
            sha,
            "--out",
            str(args.out / folder),
        ]
        if lane == "clocks":
            command += ["--expected", str(args.out / "plain")]
        if args.smoke:
            command += ["--smoke", "--calibration-calls", "200"]
        run(command, args.out / (folder + ".log"))
        exported = json.loads((args.out / folder / "results.json").read_text())
        assert exported["metadata"]["core_extension_sha256"] == sha
    # Analyze in a separate untimed process, with the reference interpreter.
    run(
        [
            str(ref / "bin/python"),
            "-B",
            str(task / "analyze.py"),
            str(args.out / "clocks"),
            str(args.out / "analysis.json"),
            "--parallel-calibration",
            str(args.out / "parallel-null.json"),
        ],
        args.out / "analysis.log",
    )
    record = {
        "workload": args.workload,
        "workers": args.workers,
        "smoke": args.smoke,
        "summary": "PASS: reference/scratch full export equal; clock and native counter closure",
    }
    (args.out / "index.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workload", choices=("dcbb96", "nsfnet", "backbone"), required=True
    )
    parser.add_argument("--workers", choices=(1, 8, 28), type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    main(parser.parse_args())
