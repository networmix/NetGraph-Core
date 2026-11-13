"""Payload-free Core clocks and a compiled full-sequence replay."""

from array import array
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import gc
import heapq
import itertools
import json
from pathlib import Path
import resource
import sys
import threading
import time
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from events import FIELDS, Blobs, Recorder, digest, provenance, public_api, read_stream  # noqa: E402
from replay import Replay  # noqa: E402
from values import Values  # noqa: E402


def interval_summary(buffers):
    streams = []
    for buffer in buffers:
        iterator = iter(buffer)
        streams.append(zip(iterator, iterator))
    first = last = edge = None
    covered = 0
    for start, stop in heapq.merge(*streams):
        if stop < start:
            raise ValueError("clock moved backwards")
        first = start if first is None else first
        covered += max(0, stop - max(start, edge if edge is not None else start))
        edge = max(edge or stop, stop)
        last = edge
    return {
        "union_s": covered / 1e9,
        "elapsed_s": (last - first) / 1e9 if first is not None else 0.0,
    }


class Clock(Recorder):
    """Recorder interception with payload recording off, identical on both paths."""

    def __init__(self, core):
        self.core, self.api = core, public_api(core)
        self.local, self.rows, self.groups = threading.local(), [], []
        self.phase, self.active = "import", True
        self.evaluators = Counter()
        self.ids, self.refs, self.released = {}, {}, Counter()
        self.counter, self.ordinals, self.lock = (
            itertools.count(),
            itertools.count(),
            threading.RLock(),
        )

    def row(self):
        if not hasattr(self.local, "row"):
            self.local.row = (
                defaultdict(lambda: [0, 0, 0]),
                defaultdict(lambda: array("q")),
            )
            self.rows.append(self.local.row)
        return self.local.row

    def emit(self, value):
        if value["kind"] == "release":
            self.released[value["type"]] += 1

    def track(self, value):
        if type(value).__name__ in FIELDS or isinstance(value, np.ndarray):
            oid = self.ids.get(id(value))
            if oid is None or self.refs[oid]() is not value:
                self.register(value)
        elif type(value) in (list, tuple):
            for child in value:
                self.track(child)
        elif type(value) is dict:
            for key, child in value.items():
                self.track(key)
                self.track(child)

    def wrap(self, name, original):
        def wrapped(*args, **kwargs):
            if not self.active or getattr(self.local, "inside", False):
                return original(*args, **kwargs)
            self.local.inside = True
            self.track(args)
            self.track(kwargs)
            stats, intervals = self.row()
            phase = getattr(self.local, "phase", self.phase)
            if name == "evaluate_pairs":
                self.evaluators[(phase, kwargs.get("workers"))] += 1
            w0, c0 = time.perf_counter_ns(), time.thread_time_ns()
            try:
                result = original(*args, **kwargs)
            finally:
                c1, w1 = time.thread_time_ns(), time.perf_counter_ns()
                row = stats[(phase, name)]
                row[0] += 1
                row[1] += w1 - w0
                row[2] += c1 - c0
                intervals[phase].extend((w0, w1))
                self.local.inside = False
            self.track(result)
            return result

        return wrapped

    def summarize(self):
        stats, phases = defaultdict(lambda: [0, 0, 0]), defaultdict(list)
        for counts, intervals in self.rows:
            for key, values in counts.items():
                stats[key] = [a + b for a, b in zip(stats[key], values, strict=True)]
            for phase, values in intervals.items():
                phases[phase].append(values)
        return {
            **interval_summary(
                [values for buffers in phases.values() for values in buffers]
            ),
            "phases": {
                phase: interval_summary(values) for phase, values in phases.items()
            },
            "released_wrappers": dict(self.released),
            "evaluator_workers": [
                {"phase": phase, "workers": workers, "calls": calls}
                for (phase, workers), calls in self.evaluators.items()
            ],
            "kinds": [
                {
                    "phase": phase,
                    "op": name,
                    "calls": n,
                    "wall_s": wall / 1e9,
                    "thread_cpu_s": cpu / 1e9,
                }
                for (phase, name), (n, wall, cpu) in sorted(stats.items())
            ],
        }


def null_control(core, count=3000):
    rows = []
    for arguments in ((), tuple(range(8))):
        clock = Clock(core)
        clock.phase = "calibration"

        def null(*args):
            return None

        wrapped = clock.wrap("calibration.null", null)
        start = time.perf_counter_ns()
        for _ in range(count):
            null(*arguments)
        direct = time.perf_counter_ns() - start
        start = time.perf_counter_ns()
        for _ in range(count):
            wrapped(*arguments)
        full = time.perf_counter_ns() - start
        measured = clock.summarize()["kinds"][0]
        rows.append(
            {
                "arguments": len(arguments),
                "calls": count,
                "in_span_ns": measured["wall_s"] * 1e9 / count,
                "full_wrapper_ns": (full - direct) / count,
            }
        )
    return {
        "rows": rows,
        "perf_counter_resolution_s": time.get_clock_info("perf_counter").resolution,
        "thread_time_resolution_s": time.get_clock_info("thread_time").resolution,
    }


class Plan(Values):
    bind = Replay.bind

    def __init__(self, core, recording, workers, expected=None):
        self.recording, self.workers = Path(recording), workers
        self.manifest = json.loads((self.recording / "manifest.json").read_text())
        if self.manifest.get("format") != 1 or self.manifest.get("fields") != FIELDS:
            raise ValueError("unknown timing recording schema")
        self.declared = expected
        build = provenance(core, expected)
        if any(
            build[key] != self.manifest["provenance"][key] for key in ("python", "gil")
        ):
            raise ValueError("timing interpreter differs from recording")
        super().__init__(core, Blobs(self.recording / "payloads"))
        self.clock = Clock(core)
        self.functions = {
            name: self.clock.wrap(name, function)
            for name, function in self.clock.api.items()
        }
        self.release = self.clock.wrap("CoreObjects.release", self.drop)
        self.jobs, self.arrays, self.payloads, self.owners = {}, {}, {}, {}
        self.observer = None
        self.expected, self.releases = Counter(), Counter()
        for job in self.manifest["counts"]:
            events = []
            for event in read_stream(self.recording, job):
                for key in ("enter", "exit", "ordinal"):
                    event.pop(key, None)
                event["phase"] = sys.intern(event["phase"])
                if event["kind"] == "call":
                    op = event["op"] = sys.intern(event["op"])
                    if op in ("FlowLedger.uid", "FlowRecord.flow_id"):
                        event["uid_object"] = event["args"]["$tuple"][0]["$object"]
                    if (
                        op == "evaluate_pairs"
                        and event["phase"] == "run"
                        and self.manifest["workload"] == "nsfnet"
                    ):
                        event["kwargs"]["$dict"] = [
                            [key, workers if key == "workers" else value]
                            for key, value in event["kwargs"]["$dict"]
                        ]
                    if op == "FlowPolicy.__init__":
                        self.owners[event["target"]] = event["args"]["$tuple"][0][
                            "$object"
                        ]
                    if op in (
                        "FlowLedger.records",
                        "FlowPolicy.records",
                        "FlowLedger.record",
                        "FlowPolicy.record",
                    ):
                        owner = event["args"]["$tuple"][0]["$object"]
                        ledger = self.owners.get(owner, owner)
                        self.record_owners(event["result"], ledger)
                    event["arguments"] = self.compile(event["args"])
                    event["keywords"] = self.compile(event["kwargs"])
                    self.expected[op] += 1
                elif event["kind"] == "release":
                    self.releases[event["type"]] += 1
                events.append(event)
            self.jobs[job] = events
        if self.expected != Counter(self.manifest["coverage"]):
            raise ValueError("incomplete timing plan")
        self.blobs = SimpleNamespace(get=lambda spec: self.payloads[spec["$array"]])
        self.pool = (
            ThreadPoolExecutor(max_workers=workers)
            if workers > 1 and self.manifest["groups"]
            else None
        )
        self.clock.row()
        if self.pool:
            barrier = threading.Barrier(workers + 1)

            def ready():
                self.clock.row()
                barrier.wait()

            futures = [self.pool.submit(ready) for _ in range(workers)]
            barrier.wait()
            for future in futures:
                future.result()

    def record_owners(self, spec, ledger):
        if isinstance(spec, dict):
            if spec.get("type") == "FlowRecord":
                self.owners[spec["$object"]] = ledger
            for value in spec.values():
                self.record_owners(value, ledger)
        elif isinstance(spec, list):
            for value in spec:
                self.record_owners(value, ledger)

    def compile(self, spec):
        if isinstance(spec, dict) and ("$object" in spec or "$numpy" in spec):
            oid = spec.get("$object", spec.get("id"))
            if "$numpy" in spec and not spec["produced"]:
                self.arrays[oid] = spec
                self.payloads[spec["$numpy"]["$array"]] = self.blobs.get(spec["$numpy"])
            return lambda: self.slots[oid]
        if isinstance(spec, dict) and "$flowid" in spec:
            ledger, number = spec["$flowid"]
            return lambda: (self.uids[ledger], self.literal(number, self.decode))
        if isinstance(spec, dict) and ("$tuple" in spec or "$list" in spec):
            key = "$tuple" if "$tuple" in spec else "$list"
            children = [self.compile(value) for value in spec[key]]
            if any(callable(value) for value in children):
                return lambda: (tuple if key == "$tuple" else list)(
                    value() if callable(value) else value for value in children
                )
            return (tuple if key == "$tuple" else list)(children)
        if isinstance(spec, dict) and "$dict" in spec:
            children = [
                (self.compile(key), self.compile(value)) for key, value in spec["$dict"]
            ]
            if any(callable(key) or callable(value) for key, value in children):
                return lambda: {
                    (key() if callable(key) else key): (
                        value() if callable(value) else value
                    )
                    for key, value in children
                }
            return dict(children)
        return self.literal(spec, self.decode)

    def decode(self, spec):
        value = self.compile(spec)
        return value() if callable(value) else value

    def put(self, oid, value, job):
        if self.observer is not None:
            self.observer.put(oid, value, job)
        else:
            self.slots[oid] = value

    def drop(self, job, event):
        if self.observer is not None:
            self.observer.release(job, event)
        else:
            del self.slots[event["oid"]]

    def run_job(self, job):
        for event in self.jobs[job]:
            self.clock.local.phase = event["phase"]
            if event["kind"] == "parallel":
                if self.observer is not None:
                    self.observer.release_boundary(job)
                if self.pool:
                    list(self.pool.map(self.run_job, event["jobs"]))
                else:
                    for child in event["jobs"]:
                        self.run_job(child)
            elif event["kind"] == "release":
                self.release(job, event)
            else:
                if self.observer is not None:
                    self.observer.release_boundary(job)
                op = event["op"]
                args, kwargs = event["arguments"], event["keywords"]
                args = args() if callable(args) else args
                kwargs = kwargs() if callable(kwargs) else kwargs
                if self.observer is not None:
                    expected = (
                        self.observer.decode(event["args"]),
                        self.observer.decode(event["kwargs"]),
                    )
                    self.observer.unchanged(
                        self.observer.guard(expected), (args, kwargs), True
                    )
                    del expected
                if event["target"] is not None:
                    shell = getattr(self.core, op.split(".")[0]).__new__(
                        getattr(self.core, op.split(".")[0])
                    )
                    self.slots[event["target"]] = shell
                    args = (shell, *args)
                    del shell
                guard = (
                    self.observer.guard((args, kwargs))
                    if self.observer is not None
                    else None
                )
                error, result = None, None
                try:
                    result = self.functions[op](*args, **kwargs)
                except Exception as exc:
                    error = {
                        "class": type(exc).__module__ + "." + type(exc).__name__,
                        "text": str(exc),
                    }
                if error != event["error"]:
                    raise ValueError(
                        f"unexpected timed exception: {job} {event['i']} {error}"
                    )
                if error is None:
                    if event["target"] is not None:
                        result = args[0]
                    Replay.bind(self, event["result"], result, job)
                    if op == "FlowLedger.uid":
                        self.uids[event["uid_object"]] = result
                    elif op == "FlowRecord.flow_id":
                        oid = event["uid_object"]
                        if oid in self.owners:
                            self.uids[self.owners[oid]] = result[0]
                if self.observer is not None:
                    self.observer.unchanged(guard)
                    self.observer.observe(job, event, result, error)
                del args, kwargs, result, guard
        if self.observer is not None:
            self.observer.release_boundary(job)

    def block(self, repeats, dry=False, observation_root=None):
        contexts = []
        observers = []
        if observation_root is not None:
            Path(observation_root).mkdir()
        for _ in range(repeats):
            context = SimpleNamespace(blobs=self.blobs, slots={})
            for spec in self.arrays.values():
                Replay.decoded_array(context, spec)
            contexts.append(context.slots)
            observers.append(
                Replay(
                    self.core,
                    self.recording,
                    Path(observation_root) / str(len(observers)),
                    self.workers,
                    "a",
                    "0/1",
                )
                if observation_root is not None
                else None
            )
        before = gc.get_stats()
        gc_control = {"enabled": gc.isenabled(), "threshold": gc.get_threshold()}
        if observation_root is None:
            for events in self.jobs.values():
                for event in events:
                    event.pop("args", None)
                    event.pop("kwargs", None)
        cpu, start = time.process_time_ns(), time.perf_counter_ns()
        for context, observer in zip(contexts, observers, strict=True):
            self.slots, self.uids = context, {}
            self.observer = observer
            if observer is not None:
                observer.slots = self.slots
            self.run_job("main")
            if dry and set(self.slots) != set(self.manifest["retained"]):
                raise ValueError("timed retained inventory differs")
            if observer is not None:
                observer.streams.close()
        outer, process = time.perf_counter_ns() - start, time.process_time_ns() - cpu
        report = self.clock.summarize()
        actual = Counter()
        for row in report["kinds"]:
            actual[row["op"]] += row["calls"]
        drops = actual.pop("CoreObjects.release", 0)
        if (
            actual
            != Counter({key: value * repeats for key, value in self.expected.items()})
            or drops != sum(self.releases.values()) * repeats
        ):
            raise ValueError("actual timed call/release coverage differs")
        report.update(
            {
                "repeats": repeats,
                "outer_wall_s": outer / 1e9,
                "process_cpu_s": process / 1e9,
                "gc_before": before,
                "gc_control": gc_control,
                "gc_after": gc.get_stats(),
                "rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                "coverage": dict(self.expected),
                "releases": dict(self.releases),
                "python_workers": self.workers if self.manifest["groups"] else 1,
                "requested_workers": self.workers,
                "recording_sha256": digest(self.recording / "manifest.json"),
                "build": provenance(self.core, self.declared),
                "gaps_s": report["elapsed_s"] - report["union_s"],
                "dry": dry,
            }
        )
        if self.pool:
            self.pool.shutdown()
        return report
