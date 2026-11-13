"""Coordinator-only native capture; --no-sample is an untimed launcher control."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import selectors
import subprocess
import sys
import time

from run_workload import ROOT, equal

REF = "208487c7aed158fb418791c22e23d09332533595053c0e4810923b31244f7c84"


def child(args):
    import _netgraph_core as native

    assert hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest() == REF
    assert Path(native.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
    sys.path.insert(0, str(ROOT / "worktrees/netgraph-orient"))
    import ngraph
    import yaml
    from ngraph.scenario import Scenario

    source = ROOT / f"scratch/orient-add/workloads/{args.workload}.yml"
    manifest = json.loads((source.parent / "manifest.json").read_text())
    assert (
        hashlib.sha256(source.read_bytes()).hexdigest()
        == manifest[args.workload]["sha256"]
    )
    document = yaml.safe_load(source.read_text())
    for step in document["workflow"]:
        if step["type"] in ("TrafficMatrixPlacement", "MaxFlow"):
            step["workers"] = args.workers
    prepared = ngraph.prepare(
        Scenario.from_yaml(yaml.safe_dump(document, sort_keys=False))
    )
    print(
        json.dumps({"ready": True, "workers": args.workers, "native_sha256": REF}),
        flush=True,
    )
    assert sys.stdin.readline().strip() == "GO"
    begin, durations = time.perf_counter(), []
    while time.perf_counter() - begin < args.seconds + 2:
        tick = time.perf_counter()
        result = ngraph.run(prepared)
        durations.append(time.perf_counter() - tick)
    with (args.out / "results.json").open("x") as stream:
        ngraph.write_json(result, stream)
    record = {
        "workers": args.workers,
        "workload": args.workload,
        "native_sha256": REF,
        "workflow_wall_s": durations,
        "complete_workflows": len(durations),
        "summary": "complete repeated workflows; final export after capture",
    }
    with (args.out / "worker.json").open("x") as stream:
        json.dump(record, stream, indent=2)
    print(json.dumps(record), flush=True)


def capture(args):
    args.out.mkdir()
    command = [
        sys.executable,
        "-B",
        __file__,
        "--worker",
        "--workload",
        args.workload,
        "--workers",
        str(args.workers),
        "--seconds",
        str(args.seconds),
        "--out",
        str(args.out),
    ]
    worker = subprocess.Popen(
        command,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    sampler = None
    begin = time.monotonic()
    with (
        (args.out / "worker.log").open("x") as log,
        (args.out / "sample.log").open("x") as slog,
    ):
        try:
            with selectors.DefaultSelector() as selector:
                selector.register(worker.stdout, selectors.EVENT_READ)
                if not selector.select(60):
                    raise TimeoutError("worker readiness timed out")
                line = worker.stdout.readline()
                log.write(line)
                assert json.loads(line)["ready"]
            if not args.no_sample:
                sampler = subprocess.Popen(
                    [
                        "/usr/bin/sample",
                        str(worker.pid),
                        str(args.seconds),
                        "1",
                        "-file",
                        str(args.out / "sample.txt"),
                    ],
                    stdout=slog,
                    stderr=subprocess.STDOUT,
                )
                time.sleep(0.4)
                assert sampler.poll() in (None, 0), "sample attach failed"
            worker.stdin.write("GO\n")
            worker.stdin.flush()
            output, _ = worker.communicate(
                timeout=max(1, 120 - (time.monotonic() - begin))
            )
            log.write(output)
            assert worker.returncode == 0, "profile worker failed"
            if sampler:
                assert sampler.wait(timeout=20) == 0
                assert (args.out / "sample.txt").stat().st_size > 0
            equal(
                json.loads((args.expected / "results.json").read_text()),
                json.loads((args.out / "results.json").read_text()),
            )
            record = {
                "sample_exit": None if sampler is None else sampler.returncode,
                "worker_exit": worker.returncode,
                "full_export_equal": True,
                "native_sha256": REF,
                "summary": "PASS: launcher and full export"
                if args.no_sample
                else "PASS: captured all threads; full export equal",
            }
            (args.out / "index.json").write_text(json.dumps(record, indent=2) + "\n")
            print(json.dumps(record), flush=True)
        finally:
            if worker.poll() is None:
                worker.kill()
                worker.communicate()
            if sampler is not None and sampler.poll() is None:
                sampler.kill()
                sampler.wait()
            log.write(f"\nSUMMARY worker_exit={worker.returncode}\n")
            slog.write(
                f"\nSUMMARY sample_exit={None if sampler is None else sampler.returncode}\n"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workload", choices=("nsfnet", "backbone"), required=True)
    parser.add_argument("--workers", type=int, choices=(8, 28), required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--expected", type=Path)
    parser.add_argument("--seconds", type=int, default=15)
    parser.add_argument("--no-sample", action="store_true")
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.seconds < 1 or (not args.worker and args.expected is None):
        parser.error("positive duration and a plain export reference are required")
    child(args) if args.worker else capture(args)
