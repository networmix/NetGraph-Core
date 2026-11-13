"""Fresh-process replay, matched consumer, and pinned worker-count inputs."""

import argparse
import copy
from collections import Counter
import gc
import json
from pathlib import Path
import resource
import sys
import time

from kernel import Clock, Plan, null_control
from events import BASE, Recorder, digest, provenance
from statecheck import exact


def consumer(args, core):
    start, cpu, before = time.perf_counter_ns(), time.process_time_ns(), gc.get_stats()
    clock = Clock(core)
    recorder = (
        Recorder(core, args.out / "recording") if args.mode == "record" else clock
    )
    recorder.install()
    sys.path.insert(0, str(BASE / "worktrees/netgraph-orient"))
    import ngraph
    from ngraph.scenario import Scenario

    if (
        Path(ngraph.__file__).resolve()
        != BASE / "worktrees/netgraph-orient/ngraph/__init__.py"
    ):
        raise ValueError("consumer import mismatch")
    recorder.scheduler()
    recorder.phase = "load"
    scenario = Scenario.from_yaml(args.input.read_text())
    recorder.phase = "prepare"
    prepared = ngraph.prepare(scenario)
    recorder.phase = "run"
    result = ngraph.run(prepared)
    recorder.phase = "export"
    with (args.out / "export.json").open("x") as stream:
        ngraph.write_json(result, stream)
    recorder.phase = "cleanup"
    del result, prepared, scenario
    gc.collect()
    if args.mode == "record":
        recorder.finish(
            {
                "workload": "controls",
                "workers": args.workers,
                "input_sha256": digest(args.input),
            }
        )
        return {"summary": "complete control recording"}
    report = clock.summarize()
    clock.active = False
    coverage = Counter()
    for row in report["kinds"]:
        coverage[row["op"]] += row["calls"]
    report.update(
        {
            "repeats": 1,
            "coverage": dict(coverage),
            "groups": clock.groups,
            "build": provenance(core),
            "input_sha256": digest(args.input),
            "rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "gc_before": before,
            "gc_after": gc.get_stats(),
            "outer_wall_s": (time.perf_counter_ns() - start) / 1e9,
            "process_cpu_s": (time.process_time_ns() - cpu) / 1e9,
        }
    )
    return report


def prepare(manifest_path, out):
    import yaml

    manifest = json.loads(manifest_path.read_text())
    rows = []
    for workload, scenario, key in (
        ("p1", 8, "dcbb96"),
        ("p2", 8, "nsfnet"),
        ("p3", 28, "backbone"),
    ):
        source = Path(manifest[key]["path"])
        if digest(source) != manifest[key]["sha256"]:
            raise ValueError("pinned source changed")
        for workers in (1, scenario):
            original = yaml.safe_load(source.read_text())
            document = copy.deepcopy(original)
            changed = []
            for index, step in enumerate(document["workflow"]):
                if "workers" in step or step["type"] in (
                    "TrafficMatrixPlacement",
                    "MaxFlow",
                    "TunnelCountAnalysis",
                ):
                    changed.append(
                        {
                            "path": f"workflow[{index}].workers",
                            "before": step.get("workers"),
                            "after": workers,
                        }
                    )
                    step["workers"] = workers
            path = out / f"{workload}-w{workers}.yml"
            with path.open("x") as stream:
                yaml.safe_dump(document, stream, sort_keys=False)
            restored = yaml.safe_load(path.read_text())
            for index, step in enumerate(original["workflow"]):
                if "workers" in step:
                    restored["workflow"][index]["workers"] = step["workers"]
                else:
                    restored["workflow"][index].pop("workers", None)
            exact(original, restored)
            document["workflow"] = []
            empty = out / f"{workload}-w{workers}-empty.yml"
            with empty.open("x") as stream:
                yaml.safe_dump(document, stream, sort_keys=False)
            rows.append(
                {
                    "workload": workload,
                    "workers": workers,
                    "input": str(path),
                    "input_sha256": digest(path),
                    "empty_input": str(empty),
                    "original": str(source),
                    "original_sha256": digest(source),
                    "changes": changed,
                    "recording": str(BASE / f"scratch/replay/record-{workload}-02"),
                }
            )
    return {
        "rows": rows,
        "summary": "six pinned inputs; only workflow worker fields changed",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("replay", "consumer", "record", "prepare"))
    parser.add_argument("input", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--dry", action="store_true")
    parser.add_argument("--expected-prefix")
    parser.add_argument("--expected-sha256")
    parser.add_argument("--expected-version")
    args = parser.parse_args()
    args.input, args.out = args.input.resolve(), args.out.resolve()
    if args.workers < 1 or args.repeats < 1:
        raise ValueError("positive workers and repeats required")
    args.out.mkdir()
    if args.mode == "prepare":
        report = prepare(args.input, args.out)
    else:
        import netgraph_core as core

        expected = None
        if any((args.expected_prefix, args.expected_sha256, args.expected_version)):
            if not all(
                (args.expected_prefix, args.expected_sha256, args.expected_version)
            ):
                raise ValueError("all expected-build fields are required")
            expected = {
                "prefix": args.expected_prefix,
                "sha256": args.expected_sha256,
                "version": args.expected_version,
            }
        provenance(core, expected)
        if args.mode == "replay":
            calibration = null_control(core)
            plan = Plan(core, args.input, args.workers, expected)
            report = plan.block(
                args.repeats, args.dry, args.out / "observations" if args.dry else None
            )
            report["null"] = calibration
        else:
            if expected is not None:
                raise ValueError(
                    "matched consumer qualification uses the frozen reference"
                )
            report = consumer(args, core)
        report["qualification"] = (
            "unqualified worker output; E0/E7 window analysis required"
        )
    report["source_sha256"] = {
        str(path): digest(path)
        for path in [
            *Path(__file__).parent.glob("*.py"),
            *Path(__file__).parent.parent.glob("*.py"),
        ]
    }
    with (args.out / "report.json").open("x") as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps({"summary": "completed", "mode": args.mode, "out": str(args.out)}))


if __name__ == "__main__":
    main()
