"""Sample DC-BB execution after loading; run the launcher outside the sandbox.

Each mode gets a fresh process. The child loads and prepares once, signals
readiness, then runs complete fresh workflows on that prepared scenario until
the sample window is covered. Export and equality checks follow sampling.
"""

import argparse
import hashlib
import json
import os
import selectors
import subprocess
import sys
import time
from pathlib import Path

EFFORT = Path(__file__).resolve().parents[3]
SCRATCH = EFFORT / "scratch/orient-consumers"


def normalize(value, parent=None):
    if isinstance(value, dict):
        return {
            key: None
            if key == "workers" and parent == "metadata"
            else normalize(item, key)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [normalize(item) for item in value]
    return value


def worker(args):
    sys.path.insert(0, str(SCRATCH / "consumer-copy"))
    import yaml

    import ngraph
    from ngraph.scenario import Scenario

    text = args.input.read_text()
    if args.workers == "1":
        document = yaml.safe_load(text)
        for step in document["workflow"]:
            if step["type"] == "TrafficMatrixPlacement":
                step["workers"] = 1
        text = yaml.safe_dump(document, sort_keys=False)
    scenario = Scenario.from_yaml(text)
    prepared = ngraph.prepare(scenario)
    assert Path(ngraph.__file__).is_relative_to(SCRATCH / "consumer-copy")
    print(json.dumps({"ready": True, "workers": args.workers}), flush=True)
    if sys.stdin.readline().strip() != "GO":
        print(json.dumps({"summary": "workflow not started"}), flush=True)
        return 2
    start = time.perf_counter()
    durations = []
    result = None
    while time.perf_counter() - start < args.seconds + 2:
        tick = time.perf_counter()
        result = ngraph.run(prepared)
        durations.append(time.perf_counter() - tick)
    compute_s = time.perf_counter() - start
    with (args.out / "results.json").open("x") as stream:
        ngraph.write_json(result, stream)
    summary = {
        "exploratory": True,
        "workers": args.workers,
        "input": str(args.input),
        "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "complete_workflows": len(durations),
        "workflow_wall_s": durations,
        "compute_window_s": compute_s,
        "python": sys.version,
        "summary": "complete DC-BB workflows finished; last result exported after sampling window",
    }
    with (args.out / "worker.json").open("x") as stream:
        json.dump(summary, stream, indent=2)
    print(json.dumps(summary), flush=True)
    return 0


def capture(args, mode):
    out = args.out / mode
    out.mkdir()
    command = [
        sys.executable,
        "-B",
        str(Path(__file__).resolve()),
        "--worker",
        "--workers",
        mode,
        "--seconds",
        str(args.seconds),
        "--input",
        str(args.input),
        "--out",
        str(out),
    ]
    start = time.monotonic()
    process = subprocess.Popen(
        command,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    sampler = None
    sample_code = None
    with (
        (out / "worker.log").open("x") as log,
        (out / "sample.log").open("x") as sample_log,
    ):
        try:
            assert process.stdout is not None and process.stdin is not None
            with selectors.DefaultSelector() as selector:
                selector.register(process.stdout, selectors.EVENT_READ)
                ready = False
                while not ready:
                    remaining = args.timeout - (time.monotonic() - start)
                    if remaining <= 0 or not selector.select(remaining):
                        raise subprocess.TimeoutExpired(command, args.timeout)
                    line = process.stdout.readline()
                    log.write(line)
                    log.flush()
                    if not line:
                        raise RuntimeError("worker exited before readiness")
                    ready = json.loads(line).get("ready", False)
            if not args.no_sample:
                sampler = subprocess.Popen(
                    [
                        "/usr/bin/sample",
                        str(process.pid),
                        str(args.seconds),
                        "1",
                        "-file",
                        str(out / "sample.txt"),
                    ],
                    stdout=sample_log,
                    stderr=subprocess.STDOUT,
                )
                # Establish the capture before releasing the prepared workload.
                time.sleep(0.4)
                if sampler.poll() not in (None, 0):
                    raise RuntimeError("sample could not attach; see sample.log")
            process.stdin.write("GO\n")
            process.stdin.flush()
            stdout, _ = process.communicate(
                timeout=max(1, args.timeout - (time.monotonic() - start))
            )
            log.write(stdout)
            if sampler is not None:
                sample_code = sampler.wait(timeout=20)
            if process.returncode != 0 or sample_code not in (None, 0):
                raise RuntimeError(f"worker={process.returncode}, sample={sample_code}")
            expected = normalize(
                json.loads((SCRATCH / "dcbb96-1-02/results.json").read_text())
            )
            actual = normalize(json.loads((out / "results.json").read_text()))
            assert json.dumps(actual) == json.dumps(expected), (
                "export differs from pinned DC-BB result"
            )
            if not args.no_sample:
                assert (out / "sample.txt").stat().st_size > 0, "empty sample output"
            record = {
                "mode": mode,
                "pid": process.pid,
                "worker_exit": process.returncode,
                "sample_exit": sample_code,
                "full_export_equal_except_metadata_workers": True,
                "summary": "launcher control passed without sampling"
                if args.no_sample
                else "native profile captured",
            }
        except BaseException:
            if process.poll() is None:
                process.kill()
            stdout, _ = process.communicate()
            log.write(stdout)
            log.write(f"\nSUMMARY failed worker_exit={process.returncode}\n")
            if sampler is not None and sampler.poll() is None:
                sampler.kill()
                sampler.wait()
            sample_log.write(
                f"\nSUMMARY failed sample_exit={None if sampler is None else sampler.returncode}\n"
            )
            raise
        log.write("\nSUMMARY " + json.dumps(record) + "\n")
        sample_log.write(
            f"\nSUMMARY sample_exit={sample_code} no_sample={args.no_sample}\n"
        )
    with (out / "index.json").open("x") as stream:
        json.dump(record, stream, indent=2)
    print(json.dumps(record), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=SCRATCH / "workloads/dcbb96.yml")
    parser.add_argument("--seconds", type=int, default=15)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--workers", choices=("1", "default", "both"), default="both")
    parser.add_argument(
        "--no-sample",
        action="store_true",
        help="launcher control only; no native profile",
    )
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.out = args.out.resolve()
    if args.seconds < 1 or args.timeout < args.seconds + 10:
        parser.error("need seconds >=1 and timeout >= seconds+10")
    if args.worker:
        return worker(args)
    args.out.mkdir()
    for mode in ("1", "default") if args.workers == "both" else (args.workers,):
        capture(args, mode)
    print("SUMMARY native profile launcher completed", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
