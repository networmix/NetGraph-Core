"""Clock, plain, or public-call-delay lane of the pinned NetGraph workflow."""

import argparse
from collections import Counter
from contextlib import contextmanager
import gc
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import threading
import time

from collector import Collector, install

ROOT = Path(__file__).resolve().parents[2]
SCRATCH = ROOT / "scratch/clocks"
sys.path.insert(0, str(SCRATCH / "t1-recorder-r1"))
from events import public_api  # noqa: E402


def equal(a, b, path="root", *, identity=False):
    if type(a) is not type(b):
        raise AssertionError((path, type(a).__name__, type(b).__name__))
    if isinstance(a, dict):
        assert list(a) == list(b), (path, "keys/order")
        for key in a:
            if identity and path == "root.metadata" and key == "core_extension_sha256":
                assert isinstance(a[key], str) and isinstance(b[key], str)
                continue
            equal(a[key], b[key], path + "." + key, identity=identity)
    elif isinstance(a, list):
        assert len(a) == len(b), (path, "length")
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            equal(x, y, f"{path}[{i}]", identity=identity)
    else:
        assert a == b, (path, a, b)


def run(args):
    args.out.mkdir()
    import netgraph_core as core
    import _netgraph_core as native

    extension = Path(native.__file__).resolve()
    assert extension.is_relative_to(Path(sys.prefix).resolve()), extension
    sha = hashlib.sha256(extension.read_bytes()).hexdigest()
    assert sha == args.native_sha, (sha, args.native_sha)
    assert sys._is_gil_enabled() and core.__version__ == "0.13.0"
    api = public_api(core)
    collector = Collector(core, api) if args.lane == "clocks" else None
    local = threading.local()
    delay_rows = []
    if collector:
        collector.registry.install()
        install(core, api, collector.wrap)
    elif args.lane == "delay":

        def factory(name, original):
            def wrapped(*a, **kw):
                if not getattr(local, "case", False) or getattr(local, "inside", False):
                    return original(*a, **kw)
                local.inside = True
                try:
                    elapsed = native._clocks_delay(args.delta_ns)
                    local.row["calls"][name] += 1
                    local.row["helper_elapsed_ns"] += elapsed
                    return original(*a, **kw)
                finally:
                    local.inside = False

            return wrapped

        install(core, api, factory)

    phases, snapshots, steps, parallels = {}, {}, [], []

    @contextmanager
    def phase(name):
        if collector:
            collector.phase = len(phases)
            collector.active = True
            snap0 = collector.registry.snapshot()
        p0, tick = time.process_time_ns(), time.perf_counter_ns()
        c0 = time.thread_time_ns()
        try:
            if collector:
                with collector.scope(name):
                    yield
            else:
                yield
        finally:
            c1 = time.thread_time_ns()
            stop_wall = time.perf_counter_ns()
            wall, process = stop_wall - tick, time.process_time_ns() - p0
            phases[name] = {
                "wall_s": wall / 1e9,
                "process_cpu_s": process / 1e9,
                "start_wall_ns": tick,
                "stop_wall_ns": stop_wall,
                "caller_start_cpu_ns": c0,
                "caller_stop_cpu_ns": c1,
            }
            if collector:
                snap1 = collector.registry.snapshot()
                snapshots[name] = [snap0[:5], snap1[:5]]
                collector.active = False

    calibration = []
    if collector:
        # No Core calls and no result encoding in null calibration.
        def null():
            return None

        name = "calibration.null"
        wrapped = collector.wrap(name, null)
        for _ in range(5):
            count = args.calibration_calls
            tick, cpu = time.perf_counter_ns(), time.thread_time_ns()
            for _ in range(count):
                null()
            direct = (time.perf_counter_ns() - tick, time.thread_time_ns() - cpu)
            collector.active = True
            tick, cpu = time.perf_counter_ns(), time.thread_time_ns()
            for _ in range(count):
                wrapped()
            full = (time.perf_counter_ns() - tick, time.thread_time_ns() - cpu)
            collector.active = False
            rows = collector.events
            calibration.append(
                {
                    "calls": count,
                    "full_wrap_wall_ns": (full[0] - direct[0]) / count,
                    "full_wrap_cpu_ns": (full[1] - direct[1]) / count,
                    "null_span_wall_ns": sum(
                        rows[i + 4] - rows[i + 3] for i in range(0, len(rows), 17)
                    )
                    / count,
                    "null_span_cpu_ns": sum(
                        rows[i + 6] - rows[i + 5] for i in range(0, len(rows), 17)
                    )
                    / count,
                }
            )
            del rows[:]
            del collector.vectors[:]

    repo = ROOT / "worktrees/netgraph-orient"
    sys.path.insert(0, str(repo))
    with phase("import"):
        import ngraph
        import yaml
        from ngraph.scenario import Scenario
    assert Path(ngraph.__file__).resolve() == repo / "ngraph/__init__.py"
    import ngraph.analysis.failures as failures
    import ngraph.analysis.placement as placement
    import ngraph.analysis.msd as msd
    import ngraph.analysis.capacity as capacity

    def scoped(original):
        def wrapped(*a, **kw):
            if not collector or not collector.active:
                return original(*a, **kw)
            with collector.scope(original.__name__):
                return original(*a, **kw)

        return wrapped

    # Replace aliases as well as defining-module references, without source edits.
    replacements = (
        {
            placement._place_state: scoped(placement._place_state),
            msd.maximum_supported_demand: scoped(msd.maximum_supported_demand),
            capacity.capacity: scoped(capacity.capacity),
        }
        if collector
        else {}
    )
    original_schedule = failures.schedule_cases

    def schedule(cases, solve, workers):
        def case(job):
            parallel = (
                workers > 1
                and threading.current_thread() is not threading.main_thread()
            )
            local.case = parallel
            if parallel and args.lane == "delay":
                local.row = {
                    "thread": threading.get_ident(),
                    "calls": Counter(),
                    "helper_elapsed_ns": 0,
                }
                delay_rows.append(local.row)
            try:
                return scoped(solve)(job)
            finally:
                local.case = False

        tick = time.perf_counter_ns()
        result = original_schedule(cases, case, workers)
        if workers > 1 and len(cases) > 1:
            parallels.append(
                {
                    "workers": workers,
                    "cases": len(cases),
                    "wall_s": (time.perf_counter_ns() - tick) / 1e9,
                }
            )
        return result

    replacements[original_schedule] = schedule
    replacements_by_id = {id(key): value for key, value in replacements.items()}
    for name, module in list(sys.modules.items()):
        if name.startswith("ngraph.") and module is not None:
            for attr, value in list(vars(module).items()):
                # Some module values are unhashable.
                if id(value) in replacements_by_id:
                    setattr(module, attr, replacements_by_id[id(value)])

    manifest = json.loads(
        (ROOT / "scratch/orient-add/workloads/manifest.json").read_text()
    )
    source = Path(manifest[args.workload]["path"])
    assert (
        hashlib.sha256(source.read_bytes()).hexdigest()
        == manifest[args.workload]["sha256"]
    )
    text = source.read_text()
    if args.workers != 8 or args.workload == "backbone" or args.smoke:
        document = yaml.safe_load(text)
        for step in document["workflow"]:
            if step["type"] in (
                "TrafficMatrixPlacement",
                "MaxFlow",
                "TunnelCountAnalysis",
            ):
                step["workers"] = args.workers
            if args.smoke and "iterations" in step:
                step["iterations"] = 2
        text = yaml.safe_dump(document, sort_keys=False)
    with phase("load"):
        scenario = Scenario.from_yaml(text)
    with phase("prepare"):
        prepared = ngraph.prepare(scenario)

    @contextmanager
    def hook(step):
        tick = time.perf_counter_ns()
        yield
        steps.append(
            {"name": step.name, "wall_s": (time.perf_counter_ns() - tick) / 1e9}
        )

    with phase("run"):
        result = ngraph.run(prepared, step_hook=hook)
    with phase("export"):
        with (args.out / "results.json").open("x") as stream:
            ngraph.write_json(result, stream)

    def cleanup():
        nonlocal result, prepared, scenario
        del result, prepared, scenario
        gc.collect()

    with phase("cleanup"):
        if collector:
            collector.measure("CoreObjects.release", cleanup)
        else:
            cleanup()
    if collector:
        collector.registry.uninstall()
    counts = sum((r["calls"] for r in delay_rows), Counter())
    summary = {
        "exploratory": True,
        "lane": args.lane,
        "workload": args.workload,
        "workers": args.workers,
        "delta_ns": args.delta_ns,
        "native_sha256": sha,
        "native": str(extension),
        "python": sys.version,
        "prefix": str(Path(sys.prefix).resolve()),
        "consumer": str(repo),
        "phases": phases,
        "phase_snapshots": snapshots,
        "steps": steps,
        "parallel_calls": dict(counts),
        "parallel_call_count": sum(counts.values()),
        "helper_elapsed_ns": sum(r["helper_elapsed_ns"] for r in delay_rows),
        "parallel_schedules": parallels,
        "calibration": calibration,
        "loadavg": os.getloadavg(),
        "smoke": args.smoke,
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    if collector:
        summary.update(collector.save(args.out))
    if args.expected:
        expected = json.loads((args.expected / "results.json").read_text())
        actual = json.loads((args.out / "results.json").read_text())
        assert actual["metadata"]["core_extension_sha256"] == sha
        equal(expected, actual, identity=True)
        summary["full_export_equal"] = True
    summary["summary"] = "workflow complete; phase clocks and full export retained"
    with (args.out / "summary.json").open("x") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")
    print(
        json.dumps(
            {
                k: v
                for k, v in summary.items()
                if k
                not in (
                    "threads",
                    "names",
                    "parallel_calls",
                    "phase_snapshots",
                    "calibration",
                )
            }
        ),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workload", choices=("dcbb96", "nsfnet", "backbone"), required=True
    )
    parser.add_argument("--workers", type=int, choices=(1, 8, 28), required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--expected", type=Path)
    parser.add_argument("--native-sha", required=True)
    parser.add_argument("--lane", choices=("plain", "clocks", "delay"), required=True)
    parser.add_argument("--delta-ns", type=int, choices=(0, 100000, 200000), default=0)
    parser.add_argument("--calibration-calls", type=int, default=30000)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.lane == "delay" and (args.workload != "dcbb96" or args.workers != 8):
        parser.error("delay lane is P1/8 only")
    if args.calibration_calls < 1:
        parser.error("calibration calls must be positive")
    run(args)


if __name__ == "__main__":
    main()
