"""Calibrate, run, and analyse the six-row E0/E7 reference window."""

import argparse
from collections import Counter
import contextlib
import copy
import gzip
import json
import math
import os
import re
from pathlib import Path
import signal
import statistics
import subprocess
import sys
import threading
import time

from kernel import digest
from events import BASE, NATIVE_HASHES, REFERENCE
from statecheck import exact

sys.path.insert(0, str(BASE / "worktrees/core-ref/dev/perf"))
from aba import bootstrap_interval, machine_state, peak_rss, rss_balanced, side_clusters  # noqa: E402

HERE = Path(__file__).resolve().parent
PYTHON = REFERENCE / "bin/python"
HOLDER = Path("/tmp/netgraph-core-heavy.lock.d/holder")


def save(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2)


def sources():
    return {
        str(path): digest(path)
        for path in [*HERE.glob("*.py"), *HERE.parent.glob("*.py")]
    }


def machine():
    state = machine_state()
    processes = subprocess.run(
        ["/bin/ps", "-axo", "pid=,pcpu=,comm="],
        capture_output=True,
        text=True,
        timeout=10,
    )
    if processes.returncode:
        raise ValueError("process activity unavailable")
    state["other_cpu"] = []
    for line in processes.stdout.splitlines():
        fields = line.strip().split(None, 2)
        if len(fields) == 3 and int(fields[0]) != os.getpid() and float(fields[1]) > 50:
            state["other_cpu"].append(
                {
                    "pid": int(fields[0]),
                    "cpu": float(fields[1]),
                    "command": Path(fields[2]).name,
                }
            )
    for key, command in (
        ("power", ["/usr/bin/pmset", "-g", "batt"]),
        ("swap", ["/usr/bin/vm_stat"]),
        ("sleep", ["/usr/sbin/sysctl", "kern.sleeptime"]),
    ):
        result = subprocess.run(command, capture_output=True, text=True, timeout=10)
        if result.returncode:
            raise ValueError(f"machine state unavailable: {command}: {result.stderr}")
        state[key] = result.stdout
    state["vm_stat_raw"] = state["swap"]
    state["swap"] = {
        name: int(value)
        for name, value in re.findall(
            r"^(Swapins|Swapouts):\s+(\d+)", state["vm_stat_raw"], re.MULTILINE
        )
    }
    if set(state["swap"]) != {"Swapins", "Swapouts"}:
        raise ValueError("swap-in/out counters unavailable")
    return state


def execute(command, out, timeout=600):
    out.mkdir()
    timed_out = []
    with (
        (out / "stdout.log").open("x") as stdout,
        (out / "stderr.log").open("x") as stderr,
    ):
        start = time.perf_counter()
        process = subprocess.Popen(
            list(map(str, command)),
            stdout=stdout,
            stderr=stderr,
            start_new_session=True,
            cwd=BASE / "worktrees/netgraph-orient",
            env={
                **os.environ,
                "PYTHONPATH": str(BASE / "worktrees/netgraph-orient"),
                "PYTHONDONTWRITEBYTECODE": "1",
            },
        )

        def stop():
            timed_out.append(True)
            with contextlib.suppress(ProcessLookupError):
                os.killpg(process.pid, signal.SIGKILL)

        timer = threading.Timer(timeout, stop)
        timer.start()
        try:
            _, status, usage = os.wait4(process.pid, 0)
        except BaseException:
            stop()
            os.wait4(process.pid, 0)
            raise
        finally:
            timer.cancel()
            timer.join()
        process.returncode = os.waitstatus_to_exitcode(status)
    report = {
        "command": list(map(str, command)),
        "pid": process.pid,
        "wall_s": time.perf_counter() - start,
        "rss_bytes": peak_rss(usage.ru_maxrss),
        "process_cpu_s": usage.ru_utime + usage.ru_stime,
        "exit": process.returncode,
        "timeout": bool(timed_out),
    }
    save(out / "process.json", report)
    if process.returncode or timed_out:
        raise ValueError(f"child failed: {out}: {report}")
    return report


def load_export(path):
    with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path) as stream:
        return json.load(stream)


def equal_export(expected, actual, path="$"):
    if type(expected) is dict and type(actual) is dict:
        if list(expected) != list(actual):
            raise ValueError(f"{path}: export keys/order differ")
        for key in expected:
            equal_export(expected[key], actual[key], path + "." + key)
    elif type(expected) is list and type(actual) is list:
        if len(expected) != len(actual):
            raise ValueError(f"{path}: export lengths differ")
        for index, (left, right) in enumerate(zip(expected, actual, strict=True)):
            equal_export(left, right, path + f"[{index}]")
    else:
        exact(expected, actual, path)


def archive_export(path, expected=None, keep=True):
    raw = path.read_bytes()
    actual = json.loads(raw)
    if expected is not None:
        equal_export(load_export(expected), actual)
    if not keep and expected is None:
        raise ValueError("cannot remove an export without its full comparison")
    archive = path.with_suffix(path.suffix + ".gz") if keep else None
    if archive is not None:
        with archive.open("xb") as stream:
            stream.write(gzip.compress(raw, mtime=0))
        if gzip.decompress(archive.read_bytes()) != raw:
            raise ValueError("export archive differs from original bytes")
    save(
        path.with_suffix(".archive.json"),
        {
            "listed_before_removal": str(path.resolve()),
            "bytes": len(raw),
            "sha256": __import__("hashlib").sha256(raw).hexdigest(),
            "compared_to": expected,
            "lossless_copy": str(archive.resolve()) if archive is not None else None,
        },
    )
    path.unlink()
    return str(archive) if archive is not None else None, actual["metadata"]


def one(row, mode, out, repeats=1, expected=None, timeout=600, keep_export=True):
    started = time.perf_counter()
    if digest(row["input"]) != row["input_sha256"]:
        raise ValueError("derived input changed")
    if mode == "cli":
        out.mkdir()
        command = [
            PYTHON,
            "-B",
            HERE / "cli.py",
            out / "cli-gc.json",
            "--quiet",
            "run",
            row["input"],
            "--results",
            out / "export.json",
        ]
        process = execute(command, out / "process", timeout)
        report = {
            "union_s": process["wall_s"],
            "elapsed_s": process["wall_s"],
            "repeats": 1,
            "rss_bytes": process["rss_bytes"],
            "fresh_wall_s": process["wall_s"],
            "process_cpu_s": process["process_cpu_s"],
            "gc": json.loads((out / "cli-gc.json").read_text()),
            "launcher_metadata_cost": "included in fresh-process wall; never subtracted",
        }
    else:
        command = [
            PYTHON,
            "-B",
            HERE / "worker.py",
            "replay" if mode == "replay" else "consumer",
            row["recording"] if mode == "replay" else row["input"],
            out / "data",
            "--workers",
            str(row["workers"]),
            "--repeats",
            str(repeats),
        ]
        out.mkdir()
        process = execute(command, out / "process", timeout)
        report = json.loads((out / "data/report.json").read_text())
        if report["source_sha256"] != sources():
            raise ValueError("timing source changed")
        if mode == "consumer" and report["input_sha256"] != row["input_sha256"]:
            raise ValueError("matched consumer input changed")
        if (
            report["coverage"]
            != json.loads((Path(row["recording"]) / "manifest.json").read_text())[
                "coverage"
            ]
        ):
            raise ValueError("matched call counts differ")
        if mode == "consumer" and any(
            group["workers"] != row["workers"] for group in report["groups"]
        ):
            raise ValueError("consumer Python worker count differs")
        if row["workload"] == "p2" and any(
            value["workers"] != row["workers"]
            for value in report["evaluator_workers"]
            if value["phase"] == "run"
        ):
            raise ValueError("native evaluator worker count differs")
        report["fresh_wall_s"] = process["wall_s"]
        report["fresh_process_cpu_s"] = process["process_cpu_s"]
        report["rss_bytes"] = max(report["rss_bytes"], process["rss_bytes"])
    if mode != "replay":
        path = out / ("export.json" if mode == "cli" else "data/export.json")
        report["export_raw_bytes"] = path.stat().st_size
        report["export"], report["export_metadata"] = archive_export(
            path, expected, keep_export
        )
        report["export_receipt"] = str(path.with_suffix(".archive.json"))
        if (
            report["export_metadata"]["core_extension_sha256"]
            != NATIVE_HASHES[str(REFERENCE)]
        ):
            raise ValueError("CLI loaded a different Core build")
        report["export_archive_bytes"] = (
            Path(report["export"]).stat().st_size if report["export"] is not None else 0
        )
    report["parent_overhead_s"] = time.perf_counter() - started - report["fresh_wall_s"]
    save(out / "block.json", report)
    return report


def aggregate(reports):
    result = {
        key: sum(row[key] for row in reports)
        for key in ("union_s", "elapsed_s", "repeats", "fresh_wall_s")
    }
    result["rss_bytes"] = max(row["rss_bytes"] for row in reports)
    kinds = {}
    for report in reports:
        for row in report.get("kinds", []):
            key = row["phase"] + ":" + row["op"]
            value = kinds.setdefault(
                key, {"calls": 0, "wall_s": 0.0, "thread_cpu_s": 0.0}
            )
            for field in value:
                value[field] += row[field]
    result["kinds"] = kinds
    released = Counter()
    for report in reports:
        released.update(report.get("released_wrappers", {}))
    result["released_wrappers"] = dict(released)
    result["process_cpu_s"] = sum(row["process_cpu_s"] for row in reports)
    result["fresh_process_cpu_s"] = sum(
        row.get("fresh_process_cpu_s", row["process_cpu_s"]) for row in reports
    )
    result["gc"] = [
        row.get("gc", {"before": row.get("gc_before"), "after": row.get("gc_after")})
        for row in reports
    ]
    return result


def pilot(args):
    import yaml

    inputs = json.loads(args.input.read_text())["rows"]
    frozen_sources = sources()
    start = time.monotonic()
    pilot_deadline = start + 595

    def remaining():
        value = pilot_deadline - time.monotonic()
        if value <= 0:
            raise ValueError("pilot exceeds its ten-minute bound")
        return min(600, value)

    def probe(*values, **options):
        return one(*values, **options, timeout=remaining())

    initial = None if args.smoke else machine()
    if not args.smoke and (
        not os.environ.get("HEAVY_SLOT_HELD")
        or os.environ.get("HEAVY_OWNER") != "replay-timing"
        or initial["load"] is None
        or initial["load"][0] >= 4
        or "AC Power" not in initial["power"]
    ):
        raise ValueError("quiet reference timing reservation required")
    if initial:
        allowed = {
            "mds",
            "mds_stores",
            "mdworker_shared",
            "mediaanalysisd",
            "photoanalysisd",
            "appstoreagent",
        }
        if len(initial["other_cpu"]) > 2 or any(
            p["cpu"] > 100 or p["command"] not in allowed for p in initial["other_cpu"]
        ):
            raise ValueError("background activity prevents a quiet timing window")
    rows = []
    startup = execute(
        [
            PYTHON,
            "-B",
            "-c",
            "import ngraph, netgraph_core; print(ngraph.__file__); print(netgraph_core.__version__)",
        ],
        args.out / "import-control",
        remaining(),
    )
    for index, row in enumerate(inputs):
        if digest(row["input"]) != row["input_sha256"]:
            raise ValueError("derived workload changed")
        recorded = json.loads((Path(row["recording"]) / "manifest.json").read_text())
        if row.get("original_sha256", row["input_sha256"]) != recorded["input_sha256"]:
            raise ValueError("consumer input does not match recording")
        if "original" in row:
            if digest(row["original"]) != row["original_sha256"]:
                raise ValueError("original input changed")
            original = yaml.safe_load(Path(row["original"]).read_text())
            restored = yaml.safe_load(Path(row["input"]).read_text())
            for i, step in enumerate(original["workflow"]):
                if "workers" in step:
                    restored["workflow"][i]["workers"] = step["workers"]
                else:
                    restored["workflow"][i].pop("workers", None)
            exact(original, restored)
        empty = probe(
            {
                **row,
                "input": row["empty_input"],
                "input_sha256": digest(row["empty_input"]),
            },
            "cli",
            args.out / f"pilot-{index}-empty",
        )
        expected = probe(row, "cli", args.out / f"pilot-{index}-cli")
        plain = args.out / f"pilot-{index}-unmodified-cli"
        plain.mkdir()
        plain_process = execute(
            [
                PYTHON,
                "-B",
                "-m",
                "ngraph",
                "--quiet",
                "run",
                row["input"],
                "--results",
                plain / "export.json",
            ],
            plain / "process",
            remaining(),
        )
        plain_archive, _ = archive_export(plain / "export.json", expected["export"])
        if index == 0:
            poisoned = copy.deepcopy(load_export(expected["export"]))
            poisoned["metadata"]["seed"] += 1
            poison = args.out / "export-poison.json.gz"
            with gzip.open(poison, "wt") as stream:
                json.dump(poisoned, stream)
            try:
                equal_export(load_export(expected["export"]), load_export(poison))
            except (AssertionError, ValueError):
                save(args.out / "export-poison-result.json", {"caught": True})
            else:
                raise ValueError("export poison accepted")
        replay = probe(row, "replay", args.out / f"pilot-{index}-replay")
        matched = probe(
            row,
            "consumer",
            args.out / f"pilot-{index}-consumer",
            expected=expected["export"],
        )
        target = 0.001 if args.smoke else 2
        counts = {
            "replay": max(1, math.ceil(target * 1.1 / replay["union_s"])),
            "cli": max(1, math.ceil(target * 1.1 / expected["union_s"])),
            "matched": max(
                1, math.ceil(target * 1.1 / min(replay["union_s"], matched["union_s"]))
            ),
        }

        # E7 uses identical full-sequence counts on both paths.
        def replay_cost(count):
            return (
                replay["fresh_wall_s"]
                - replay["outer_wall_s"]
                + count * replay["outer_wall_s"]
                + replay["parent_overhead_s"]
            )

        projected = (
            21
            * (
                replay_cost(counts["replay"])
                + counts["cli"] * expected["fresh_wall_s"]
                + counts["cli"] * expected["parent_overhead_s"]
            )
            + 7
            * (
                replay_cost(counts["matched"])
                + 2
                * counts["matched"]
                * (matched["fresh_wall_s"] + matched["parent_overhead_s"])
            )
            + 63 * 0.25
        )
        rows.append(
            {
                **row,
                "counts": counts,
                "target_s": target,
                "expected_export": expected["export"],
                "empty_run": empty,
                "unmodified_cli_control": {
                    **plain_process,
                    "export": plain_archive,
                    "gc_launcher_ratio_diagnostic": expected["fresh_wall_s"]
                    / plain_process["wall_s"],
                },
                "null": replay["null"],
                "projected_s": projected,
                "projected_archive_bytes": (21 * counts["cli"] + 14 * counts["matched"])
                * expected["export_archive_bytes"],
                "maximum_export_bytes": max(
                    expected["export_raw_bytes"], matched["export_raw_bytes"]
                ),
            }
        )
    plan = {
        "rows": rows,
        "source_sha256": frozen_sources,
        "started_monotonic": start,
        "holder": HOLDER.read_text() if HOLDER.exists() else None,
        "initial_machine": initial,
        "import_control": startup,
        "pilot_s": time.monotonic() - start,
        "projected_s": sum(row["projected_s"] for row in rows) * 1.2,
        "smoke": args.smoke,
        "allocator_environment": {
            key: value
            for key, value in os.environ.items()
            if key.startswith(("Malloc", "MALLOC_"))
            or key
            in (
                "PYTHONMALLOC",
                "PYTHONHASHSEED",
                "DYLD_INSERT_LIBRARIES",
                "LD_PRELOAD",
                "GLIBC_TUNABLES",
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
            )
        },
    }
    disk = subprocess.check_output(
        ["du", "-sk", str(BASE / "scratch/replay")], text=True, timeout=30
    )
    plan["existing_disk_bytes"] = int(disk.split()[0]) * 1024
    plan["projected_peak_bytes"] = (
        plan["existing_disk_bytes"]
        + 2 * sum(row["projected_archive_bytes"] for row in rows)
        + max(row["maximum_export_bytes"] for row in rows)
        + 128 * 1024**2
    )
    if sources() != frozen_sources:
        raise ValueError("source changed during pilot")
    save(args.out / "plan.json", plan)
    print(json.dumps(plan), flush=True)


def run(args):
    plan = json.loads(args.input.read_text())
    if plan["source_sha256"] != sources() or args.smoke != plan["smoke"]:
        raise ValueError("plan/source mismatch")
    deadline = (
        time.monotonic() + 600
        if args.smoke
        else plan["started_monotonic"] + args.minutes * 60
    )
    if not args.smoke and (
        not os.environ.get("HEAVY_SLOT_HELD")
        or os.environ.get("HEAVY_OWNER") != "replay-timing"
        or HOLDER.read_text() != plan["holder"]
        or time.monotonic() + plan["projected_s"] > deadline
        or plan["projected_peak_bytes"] > 20 * 1024**3
    ):
        raise ValueError("window reservation or remaining two-hour budget insufficient")
    blocks = []
    print(
        json.dumps(
            {
                "projected_s_with_reserve": plan["projected_s"],
                "remaining_s": deadline - time.monotonic(),
                "rows": len(plan["rows"]),
                "projected_peak_bytes": plan["projected_peak_bytes"],
            }
        ),
        flush=True,
    )
    for index, row in enumerate(plan["rows"]):
        for experiment in ("replay", "cli", "matched"):
            for triple in range(1 if args.smoke else 7):
                for label in ("A1", "B", "A2"):
                    if sources() != plan["source_sha256"]:
                        raise ValueError("source changed during timing window")
                    mode = (
                        ("replay" if label == "B" else "consumer")
                        if experiment == "matched"
                        else experiment
                    )
                    count = row["counts"][experiment]
                    before = None if args.smoke else machine()
                    reports = []
                    for repeat in range(1 if mode == "replay" else count):
                        remaining = deadline - time.monotonic()
                        if remaining <= 0:
                            raise ValueError("window deadline reached")
                        reports.append(
                            one(
                                row,
                                mode,
                                args.out
                                / f"r{index}-{experiment}-{triple}-{label}-{repeat}",
                                count if mode == "replay" else 1,
                                row["expected_export"],
                                min(600, remaining),
                                keep_export=False,
                            )
                        )
                    after = None if args.smoke else machine()
                    if sources() != plan["source_sha256"]:
                        raise ValueError("source changed during timing block")
                    invalid = after and (
                        "AC Power" not in after["power"]
                        or after["sleep"] != plan["initial_machine"]["sleep"]
                        or time.monotonic() > deadline
                    )
                    block = {
                        "row": index,
                        "experiment": experiment,
                        "triple": triple,
                        "slot": label,
                        "before": before,
                        "after": after,
                        **aggregate(reports),
                    }
                    save(
                        args.out / f"block-r{index}-{experiment}-{triple}-{label}.json",
                        block,
                    )
                    blocks.append(block)
                    stored = sum(
                        path.stat().st_size
                        for path in args.out.rglob("*")
                        if path.is_file()
                    )
                    if (
                        plan["existing_disk_bytes"]
                        + 2 * stored
                        + max(row["maximum_export_bytes"] for row in plan["rows"])
                        > 20 * 1024**3
                    ):
                        raise ValueError(
                            "stored window evidence exceeds its 20 GiB budget"
                        )
                    if invalid:
                        raise ValueError("power, sleep or deadline voids timing window")
                print(
                    json.dumps(
                        {
                            "row": index,
                            "experiment": experiment,
                            "triple": triple,
                            "completed": True,
                        }
                    ),
                    flush=True,
                )
    save(
        args.out / "window.json",
        {
            "plan": plan,
            "blocks": blocks,
            "smoke": args.smoke,
            "summary": "complete raw E0/E7 window; analyse before qualification",
        },
    )


def statistics_for(triples, metric):
    times = [[slot[metric] / slot["repeats"] for slot in triple] for triple in triples]
    ratios = [b / ((a + c) / 2) for a, b, c in times]
    drift = statistics.median(abs(a / c - 1) for a, _, c in times)
    ci = bootstrap_interval(ratios, 10000, 20261108)
    clusters = {
        "A": side_clusters([row[k] for row in times for k in (0, 2)], 0.03),
        "B": side_clusters([row[1] for row in times], 0.03),
    }
    return {
        "median_ratio": statistics.median(ratios),
        "ci95": ci,
        "reference_drift": drift,
        "clusters": clusters,
        "qualified_noise": len(triples) >= 7
        and drift <= 0.03
        and not any(clusters.values()),
        "beyond_noise_gain": ci[1] / 0.97,
        "native_bar_gain": 1.15 * ci[1],
    }


def analyse(args):
    window = json.loads(args.input.read_text())
    rows = []
    for index, row in enumerate(window["plan"]["rows"]):
        for experiment in ("replay", "cli", "matched"):
            blocks = [
                b
                for b in window["blocks"]
                if b["row"] == index and b["experiment"] == experiment
            ]
            triples = [
                [
                    next(b for b in blocks if b["triple"] == t and b["slot"] == slot)
                    for slot in ("A1", "B", "A2")
                ]
                for t in sorted({b["triple"] for b in blocks})
            ]
            metrics = {
                metric: statistics_for(triples, metric)
                for metric in ("union_s", "elapsed_s")
            }
            memory = rss_balanced(
                [
                    {
                        label: {"rss_tree_bytes": slot["rss_bytes"]}
                        for label, slot in zip(("A1", "B", "A2"), triple, strict=True)
                    }
                    for triple in triples
                ],
                1.10,
                8 * 1024**2,
            )
            short = any(b["union_s"] < 2 for b in blocks)
            swap = any(
                b["before"] and b["before"]["swap"] != b["after"]["swap"]
                for b in blocks
            )
            directions = {
                "faster"
                if m["ci95"][1] < 0.97
                else "slower"
                if m["ci95"][0] > 1.03
                else "unresolved"
                for m in metrics.values()
            }
            disagreement = "faster" in directions and "slower" in directions
            kinds = {}
            if experiment == "matched":
                release_match = all(
                    b["released_wrappers"] == triples[0][0]["released_wrappers"]
                    for b in blocks
                )
                for name in sorted(set.union(*(set(b["kinds"]) for b in blocks))):
                    if name.endswith(":CoreObjects.release"):
                        continue
                    if any(name not in b["kinds"] for b in blocks):
                        kinds[name] = {
                            "qualified": False,
                            "counts_match": False,
                            "reason": "phase/kind missing from a matched block",
                        }
                        continue
                    count_match = (
                        len({b["kinds"][name]["calls"] / b["repeats"] for b in blocks})
                        == 1
                    )
                    derived = [
                        [{**b, "union_s": b["kinds"][name]["wall_s"]} for b in triple]
                        for triple in triples
                    ]
                    stat = statistics_for(derived, "union_s")
                    mean_ns = min(
                        b["kinds"][name]["wall_s"] / b["kinds"][name]["calls"] * 1e9
                        for b in blocks
                    )
                    null_ns = max(r["in_span_ns"] for r in row["null"]["rows"])
                    stat.update(
                        {
                            "counts_match": count_match,
                            "release_counts_match": release_match,
                            "calls_per_iteration": blocks[0]["kinds"][name]["calls"]
                            / blocks[0]["repeats"],
                            "block_means_ns": [
                                {
                                    "slot": b["slot"],
                                    "triple": b["triple"],
                                    "wall": b["kinds"][name]["wall_s"]
                                    * 1e9
                                    / b["kinds"][name]["calls"],
                                    "thread_cpu": b["kinds"][name]["thread_cpu_s"]
                                    * 1e9
                                    / b["kinds"][name]["calls"],
                                }
                                for b in blocks
                            ],
                            "null_fraction": null_ns / mean_ns,
                            "qualified": count_match
                            and release_match
                            and stat["qualified_noise"]
                            and 0.95 <= stat["ci95"][0] <= stat["ci95"][1] <= 1.05
                            and null_ns / mean_ns < 0.05
                            and not short
                            and not swap
                            and not window["smoke"],
                        }
                    )
                    kinds[name] = stat
            rows.append(
                {
                    "workload": row["workload"],
                    "workers": row["workers"],
                    "experiment": experiment,
                    "metrics": metrics,
                    "memory": memory,
                    "short_blocks": short,
                    "swap_changed": swap,
                    "union_elapsed_disagree": disagreement,
                    "kinds": kinds,
                    "release_qualification": "not qualified; replay drops versus historical lifetime diagnostics only",
                    "denominators": "fresh CLI process wall including GC metadata launcher"
                    if experiment == "cli"
                    else "Core span union and first-to-last elapsed; matched consumer blocks sum fresh-process intervals",
                    "aggregate_fidelity": "not qualified as a whole; E7 qualifies API kinds individually and excludes release decisions"
                    if experiment == "matched"
                    else "E0 noise only",
                    "qualified_noise": not window["smoke"]
                    and not short
                    and not swap
                    and not disagreement
                    and all(m["qualified_noise"] for m in metrics.values()),
                }
            )
    save(
        args.out / "analysis.json",
        {
            "rows": rows,
            "smoke": window["smoke"],
            "summary": "raw data analysed; per-kind E7 qualification is explicit",
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("pilot", "run", "analyse"))
    parser.add_argument("input", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--minutes", type=float, default=120)
    args = parser.parse_args()
    args.input, args.out = args.input.resolve(), args.out.resolve()
    if not 0 < args.minutes <= 120:
        raise ValueError("window must be at most two hours")
    args.out.mkdir()
    {"pilot": pilot, "run": run, "analyse": analyse}[args.mode](args)


if __name__ == "__main__":
    main()
