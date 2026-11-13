"""Run one bounded recording or a sharded differential check campaign."""

import argparse
import json
import os
import subprocess
from collections import Counter
from pathlib import Path

from events import BASE, PYPI, REFERENCE, digest
from values import dumps

HERE = Path(__file__).resolve().parent


def storage_guard(root):
    size = sum(path.stat().st_size for path in root.rglob("*") if path.is_file())
    if size > 20 * 1024**3:
        raise ValueError(f"aggregate evidence exceeds 20 GiB: {root}: {size}")


def prune_streams(out, shard, sides, comparison_logs):
    """Retain the two plant children; audit deletion of later completed children."""
    if shard < 2:
        return
    for log in comparison_logs:
        if '"summary":"PASS exact observations"' not in log.read_text():
            raise ValueError(f"cannot prune before a clean comparison: {log}")
    files = []
    for side, _ in sides:
        directory = (out / f"s{shard:02d}-{side}").resolve()
        for path in sorted(
            [*directory.glob("*.jsonl.gz"), *(directory / "payloads").glob("*.npz")]
        ):
            files.append(
                {
                    "path": str(path),
                    "bytes": path.stat().st_size,
                    "sha256": digest(path),
                }
            )
    listed = out / f"s{shard:02d}-deletion-listed.json"
    with listed.open("x") as stream:
        json.dump(
            {"files": files, "comparisons": list(map(str, comparison_logs))},
            stream,
            indent=2,
        )
    print(dumps({"deletion_list": str(listed), "files": len(files)}), flush=True)
    for item in files:
        path = Path(item["path"])
        if path.stat().st_size != item["bytes"] or digest(path) != item["sha256"]:
            raise ValueError(f"stream changed before deletion: {path}")
    for item in files:
        Path(item["path"]).unlink()
    if any(Path(item["path"]).exists() for item in files):
        raise ValueError("listed stream deletion incomplete")
    with (out / f"s{shard:02d}-deletion-done.json").open("x") as stream:
        json.dump(
            {
                "listed": str(listed),
                "bytes": sum(v["bytes"] for v in files),
                "summary": "PASS listed streams removed after all comparisons",
            },
            stream,
            indent=2,
        )


def bounded(command, log, seconds):
    child = [
        str(REFERENCE / "bin/python"),
        "-B",
        str(BASE / "scratch/orient-add/scripts/bounded.py"),
        str(seconds),
        str(log),
        *map(str, command),
    ]
    print(dumps({"starting": list(map(str, command)), "limit_s": seconds}), flush=True)
    subprocess.run(
        ["bash", str(BASE / "rules/tools/heavy.sh"), *child],
        env={**os.environ, "HEAVY_OWNER": "replay", "PYTHONDONTWRITEBYTECODE": "1"},
        check=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    rec = sub.add_parser("record")
    rec.add_argument("workload", choices=("dcbb96", "nsfnet", "backbone"))
    rec.add_argument("out", type=Path)
    check = sub.add_parser("check")
    check.add_argument("recording", type=Path)
    check.add_argument("out", type=Path)
    check.add_argument("--workers", required=True, type=int)
    check.add_argument("--lane", required=True, choices=("a", "b"))
    check.add_argument("--shards", type=int, default=1)
    check.add_argument("--pypi", action="store_true")
    check.add_argument("--reference", type=Path)
    check.add_argument("--state-sample", action="store_true")
    check.add_argument("--prune", action="store_true")
    check.add_argument("--storage-root", type=Path)
    matrix = sub.add_parser("matrix")
    matrix.add_argument("recordings", type=Path)
    matrix.add_argument("out", type=Path)
    matrix.add_argument("--p1-reference", type=Path)
    matrix.add_argument("--p1-shards", type=int, default=2)
    args = parser.parse_args()
    if args.command == "record":
        bounded(
            [
                REFERENCE / "bin/python",
                "-B",
                HERE / "replay.py",
                "record",
                args.workload,
                args.out,
            ],
            Path(str(args.out) + ".log"),
            180,
        )
        return
    if args.command == "matrix":
        args.out.mkdir()
        rows = []
        for workload, scenario in (("p1", 8), ("p2", 8), ("p3", 28)):
            for lane, workers in (("a", 1), ("a", scenario), ("b", scenario)):
                name = f"{workload}-{lane}{workers}"
                command = [
                    str(REFERENCE / "bin/python"),
                    "-B",
                    str(Path(__file__)),
                    "check",
                    str(args.recordings / f"record-{workload}-02"),
                    str(args.out / name),
                    "--workers",
                    str(workers),
                    "--lane",
                    lane,
                ]
                command.append("--prune")
                command.extend(("--storage-root", str(args.recordings)))
                if lane == "a" and workers == 1:
                    command.append("--pypi")
                if workload == "p1":
                    command.append("--state-sample")
                    command.extend(("--shards", str(args.p1_shards)))
                    if args.p1_reference is not None and lane == "a" and workers == 8:
                        command.extend(("--reference", str(args.p1_reference)))
                subprocess.run(command, check=True)
                rows.append(
                    {"row": name, "summary": str(args.out / name / "summary.json")}
                )
        (args.out / "summary.json").write_text(
            json.dumps(
                {"rows": rows, "summary": "PASS sequential stage-1 matrix"}, indent=2
            )
            + "\n"
        )
        return
    if args.reference is not None and args.shards != 1:
        raise ValueError("an existing reference requires one observation shard")
    args.out.mkdir()
    report = []
    observed, state_coverage = Counter(), Counter()
    observed_jobs, state_jobs = set(), set()
    for shard in range(args.shards):
        storage_guard(args.storage_root or args.out)
        sides = [("ref1", REFERENCE), ("ref2", REFERENCE)]
        if args.pypi:
            sides.append(("pypi", PYPI))
        reference = (
            args.reference
            if args.reference is not None
            else args.out / f"s{shard:02d}-ref1"
        )
        if args.reference is not None:
            header = json.loads((reference / "manifest.json").read_text())
            if (
                header["workers"] != args.workers
                or header["lane"] != args.lane
                or header["shard"] != [0, 1]
            ):
                raise ValueError("existing reference settings differ")
        for side, environment in sides:
            if side == "ref1" and args.reference is not None:
                continue
            name = f"s{shard:02d}-{side}"
            bounded(
                [
                    environment / "bin/python",
                    "-B",
                    HERE / "replay.py",
                    "check",
                    args.recording,
                    args.out / name,
                    "--workers",
                    str(args.workers),
                    "--lane",
                    args.lane,
                    "--shard",
                    f"{shard}/{args.shards}",
                    *(["--state-sample"] if args.state_sample else []),
                ],
                args.out / (name + ".log"),
                600,
            )
            storage_guard(args.storage_root or args.out)
        comparison_logs = []
        for side, _ in sides[1:]:
            log = args.out / f"s{shard:02d}-compare-{side}.log"
            bounded(
                [
                    REFERENCE / "bin/python",
                    "-B",
                    HERE / "replay.py",
                    "compare",
                    reference,
                    args.out / f"s{shard:02d}-{side}",
                ],
                log,
                600,
            )
            report.append(
                {
                    "shard": shard,
                    "side": side,
                    "comparison": str(log),
                    "outcome": "PASS",
                }
            )
            comparison_logs.append(log)
        header = json.loads((reference / "manifest.json").read_text())
        if observed_jobs.intersection(header["jobs"]):
            raise ValueError("overlapping returned-value job coverage")
        observed_jobs.update(header["jobs"])
        observed.update(header["observed"])
        state_coverage.update(header["state_coverage"])
        state_jobs.update(job for job in header["state_jobs"] if job in header["jobs"])
        if args.prune:
            prune_streams(args.out, shard, sides, comparison_logs)
    recording = json.loads((args.recording / "manifest.json").read_text())
    if observed != Counter(recording["coverage"]) or observed_jobs != set(
        recording["counts"]
    ):
        raise ValueError("incomplete returned-value coverage across shards")
    if state_jobs != set(header["state_jobs"]):
        raise ValueError("incomplete sampled-state coverage across shards")
    (args.out / "summary.json").write_text(
        json.dumps(
            {
                "rows": report,
                "returned_coverage": dict(observed),
                "state_coverage": dict(state_coverage),
                "state_jobs": sorted(state_jobs),
                "summary": "PASS all differential shards",
            },
            indent=2,
        )
        + "\n"
    )
    print(
        dumps({"summary": "PASS all differential shards", "comparisons": len(report)}),
        flush=True,
    )


if __name__ == "__main__":
    main()
