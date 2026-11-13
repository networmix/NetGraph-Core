"""Record, check and compare a complete Core boundary sequence in wheel processes."""

import argparse
import gc
import json
import itertools
import sys
import threading
import time
import weakref
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from events import (
    FIELDS,
    TOOL_HASHES,
    Blobs,
    Streams,
    digest,
    provenance,
    public_api,
    read_stream,
    record,
)
from values import Snapshots, Values, dumps, first_difference


class Replay(Values):
    def __init__(
        self,
        core,
        recording,
        out,
        workers,
        lane,
        shard,
        expected=None,
        state_sample=False,
    ):
        self.expected = expected
        self.recording, self.out = Path(recording), Path(out)
        self.manifest = json.loads((self.recording / "manifest.json").read_text())
        if (
            self.manifest.get("format") != 1
            or self.manifest.get("fields") != FIELDS
            or not self.manifest.get("coverage")
        ):
            raise ValueError("unknown or empty recording schema")
        self.out.mkdir()
        super().__init__(core, Blobs(self.recording / "payloads"))
        self.streams = Streams(self.out, 20 * 1024**3)
        self.snapshots = Snapshots(
            core, Blobs(self.out / "payloads", write=True, charge=self.streams.charge)
        )
        self.api = public_api(core)
        self.slots, self.owner_jobs = {}, {}
        self.pending_dead = {}
        self.workers, self.lane = workers, lane
        self.shard, self.shards = map(int, shard.split("/"))
        if not 0 <= self.shard < self.shards:
            raise ValueError("invalid observation shard")
        self.job_indices = {
            job: i
            for i, job in enumerate(
                job for group in self.manifest["groups"] for job in group["jobs"]
            )
        }
        sampled = {
            job for job, index in self.job_indices.items() if index % 64 in (0, 1)
        }
        for group in self.manifest["groups"]:
            if group["jobs"] and not sampled.intersection(group["jobs"]):
                sampled.add(group["jobs"][0])
        self.state_jobs = (
            {"main", *sampled} if state_sample else {"main", *self.job_indices}
        )
        self.state_sampling = (
            "mod64:0,1+uncovered-group-first" if state_sample else "all"
        )
        self.coverage, self.observed = Counter(), Counter()
        self.state_coverage = Counter()
        self.executed, self.released = set(), Counter()
        self.lock = threading.RLock()
        self.want = core.OutputOptions(**dict.fromkeys(FIELDS["OutputOptions"], True))

    def selected(self, job):
        return (
            self.shard == 0
            if job == "main"
            else self.job_indices[job] % self.shards == self.shard
        )

    def decoded_array(self, spec):
        payload = self.blobs.get(spec["$numpy"])
        oid = spec["id"]
        if oid in self.slots:
            obj = self.slots[oid]
            if (
                obj.dtype != payload.dtype
                or obj.shape != payload.shape
                or obj.tobytes() != payload.tobytes()
            ):
                raise ValueError(f"native input array changed before call: {oid}")
            return obj
        if spec["produced"]:
            raise ValueError(f"missing array producer {oid}")
        strides, shape = spec["strides"], payload.shape
        low = sum(
            min(0, (n - 1) * stride) for n, stride in zip(shape, strides, strict=True)
        )
        high = sum(
            max(0, (n - 1) * stride) for n, stride in zip(shape, strides, strict=True)
        )
        storage = bytearray(
            max(payload.dtype.itemsize, high - low + payload.dtype.itemsize)
        )
        obj = np.ndarray(
            shape,
            dtype=payload.dtype,
            buffer=storage,
            offset=-low,
            strides=tuple(strides),
        )
        obj[...] = payload
        obj.flags.writeable = not spec["readonly"]
        self.slots[oid] = obj
        return obj

    def decode(self, spec):
        if isinstance(spec, dict):
            if "$object" in spec:
                oid = spec["$object"]
                if (
                    oid not in self.slots
                    or type(self.slots[oid]).__name__ != spec["type"]
                ):
                    raise ValueError(f"missing or wrong native producer {oid}")
                return self.slots[oid]
            if "$numpy" in spec:
                return self.decoded_array(spec)
            if "$flowid" in spec:
                ledger, number = spec["$flowid"]
                matches = [
                    uid
                    for uid, owner in self.snapshots.uids["FlowLedger"].items()
                    if owner == ledger
                ]
                if len(matches) != 1:
                    raise ValueError(f"FlowId has no unique ledger {ledger}")
                return (matches[0], self.decode(number))
        return self.literal(spec, self.decode)

    def put(self, oid, obj, job):
        if oid in self.slots and self.slots[oid] is not obj:
            raise ValueError(f"object alias differs: {oid}")
        self.slots[oid] = obj
        self.snapshots.register(obj, oid)
        if type(obj).__name__ in ("FlowLedger", "FlowPolicy", "ShortestPathCache"):
            self.owner_jobs[oid] = job

    def bind(self, spec, obj, job):
        if not isinstance(spec, dict):
            return
        if "$object" in spec:
            if type(obj).__name__ != spec["type"]:
                raise ValueError(f"returned object type differs: {spec['$object']}")
            self.put(spec["$object"], obj, job)
        elif "$numpy" in spec:
            if not isinstance(obj, np.ndarray):
                raise ValueError("native array return type differs")
            owner = obj.base
            while isinstance(owner, np.ndarray):
                owner = owner.base
            if spec["owner"] is not None and owner is not self.slots.get(spec["owner"]):
                raise ValueError(f"array owner differs: {spec['id']}")
            self.put(spec["id"], obj, job)
        elif "$list" in spec or "$tuple" in spec:
            key = "$list" if "$list" in spec else "$tuple"
            if type(obj) is not (list if key == "$list" else tuple) or len(obj) != len(
                spec[key]
            ):
                raise ValueError("native return container differs")
            for expected, actual in zip(spec[key], obj, strict=True):
                self.bind(expected, actual, job)
        elif "$dict" in spec:
            for key, value in spec["$dict"]:
                self.bind(value, obj[self.decode(key)], job)

    def guard(self, obj):
        if isinstance(obj, np.ndarray):
            return (
                "array",
                obj,
                obj.copy(),
                obj.dtype.str,
                obj.shape,
                obj.strides,
                obj.flags.writeable,
            )
        if type(obj) in (tuple, list):
            return ("sequence", obj, type(obj), len(obj), [self.guard(x) for x in obj])
        if type(obj) is dict:
            return ("mapping", obj, list(obj), [self.guard(v) for v in obj.values()])
        if type(obj).__name__ in FIELDS:
            return ("identity", obj)
        return ("atom", obj, type(obj), self.atom(obj, lambda x: x))

    def unchanged(self, guard, current=None, supplied=False):
        tag, original, *data = guard
        obj = current if supplied else original
        if tag == "identity":
            if obj is not original:
                raise ValueError("Core replaced a supplied native reference")
            return
        if tag == "atom":
            kind, value = data
            if type(obj) is not kind or self.atom(obj, lambda x: x) != value:
                raise ValueError("Core mutated a supplied scalar")
            return
        if tag == "array":
            before, dtype, shape, strides, writable = data
            if (
                not isinstance(obj, np.ndarray)
                or obj.dtype.str != dtype
                or obj.shape != shape
                or obj.strides != strides
                or obj.flags.writeable != writable
                or obj.tobytes() != before.tobytes()
            ):
                raise ValueError("Core mutated a supplied input array")
        elif tag == "sequence":
            kind, length, children = data
            if type(obj) is not kind or len(obj) != length:
                raise ValueError("Core mutated a supplied input sequence")
            for actual, child in zip(obj, children, strict=True):
                self.unchanged(child, actual, True)
        else:
            keys, children = data
            if list(obj) != keys:
                raise ValueError("Core mutated a supplied input mapping")
            for actual, child in zip(obj.values(), children, strict=True):
                self.unchanged(child, actual, True)

    def observe(self, job, event, obj, error):
        if not self.selected(job):
            return
        op = event["op"]
        if error is not None:
            returned = None
        elif event["target"] is not None and type(obj).__name__ in (
            "FlowLedger",
            "FlowPolicy",
            "ShortestPathCache",
        ):
            returned = {"$owner": event["target"]}
        elif op in ("Graph.uid", "FlowLedger.uid", "Eligibility.uid"):
            returned = self.snapshots.uid(op.split(".")[0], obj)
        elif op == "FlowRecord.flow_id":
            returned = self.snapshots.snapshot(obj, "flow_id")
        elif op.endswith(".graph_uid"):
            returned = self.snapshots.uid("Graph", obj)
        else:
            returned = self.snapshots.snapshot(obj)
        owners = {}
        if job in self.state_jobs:
            for oid, owner_job in list(self.owner_jobs.items()):
                if owner_job == job:
                    owners[oid] = self.snapshots.state(oid, self.slots[oid])
            self.state_coverage[op] += 1
        self.streams.write(
            job,
            {
                "call": event["i"],
                "op": op,
                "result": self.snapshots.pack(returned),
                "error": error,
                "state": owners,
                "inputs_unchanged": True,
            },
        )
        self.observed[op] += 1

    def call(self, job, event):
        op = event["op"]
        if op not in self.api:
            raise ValueError(f"unknown operation {op}")
        args, kwargs = self.decode(event["args"]), self.decode(event["kwargs"])
        if (
            op == "evaluate_pairs"
            and event["phase"] == "run"
            and self.manifest["workload"] == "nsfnet"
        ):
            kwargs["workers"] = self.workers
        if self.lane == "b":
            if op == "OutputOptions.__init__":
                kwargs = dict.fromkeys(FIELDS["OutputOptions"], True)
            elif op == "EvaluationOptions.__init__":
                kwargs["outputs"] = self.want
            elif op in (
                "apply_actions",
                "FlowPolicy.evaluate",
                "FlowPolicy.replace",
                "FlowPolicy.forward",
            ):
                kwargs["outputs"] = self.want
            elif op == "FlowPolicy.place":
                if len(args) >= 3:
                    args = (*args[:2], self.want)
                else:
                    kwargs["outputs"] = self.want
        if event["target"] is not None:
            cls = getattr(self.core, op.split(".")[0])
            shell = cls.__new__(cls)
            self.slots[event["target"]] = shell
            args = (shell, *args)
        inputs = self.guard((args, kwargs))
        error, result = None, None
        try:
            result = self.api[op](*args, **kwargs)
        except Exception as exc:
            error = {
                "class": type(exc).__module__ + "." + type(exc).__name__,
                "text": str(exc),
            }
        self.unchanged(inputs)
        if error != event["error"]:
            raise ValueError(
                f"job={job} call={event['i']} {op}: exception mismatch: recorded={event['error']!r}, actual={error!r}"
            )
        if error is None:
            if event["target"] is not None:
                result = args[0]
            self.bind(event["result"], result, job)
        self.observe(job, event, result, error)
        self.coverage[op] += 1

    def release(self, job, event):
        oid = event["oid"]
        if oid not in self.slots:
            raise ValueError(f"release of missing object {oid}")
        ref = weakref.ref(self.slots[oid])
        self.snapshots.owners.pop(id(self.slots[oid]), None)
        self.slots.pop(oid)
        self.owner_jobs.pop(oid, None)
        had_state = self.snapshots.previous.pop(oid, None) is not None
        self.pending_dead.setdefault(job, []).append((oid, ref))
        if self.selected(job):
            self.streams.write(
                job,
                {
                    "call": event["i"],
                    "op": "release",
                    "oid": oid,
                    "type": event["type"],
                    "state_removed": had_state,
                },
            )
        self.released[event["type"]] += 1

    def release_boundary(self, job):
        # A native parent can destroy children before its own weakref callback.
        # Drop all adjacent recorded references, then check before the next call.
        for oid, ref in self.pending_dead.pop(job, []):
            if ref() is not None:
                raise ValueError(f"wrapper retained across a call boundary: {oid}")

    def run_job(self, job):
        last_phase = None
        with self.lock:
            if job in self.executed:
                raise ValueError(f"duplicate replay job {job}")
            self.executed.add(job)
        for event in read_stream(self.recording, job):
            try:
                if event["kind"] == "call":
                    self.release_boundary(job)
                    self.call(job, event)
                elif event["kind"] == "release":
                    self.release(job, event)
                elif event["kind"] == "parallel":
                    print(
                        dumps(
                            {
                                "progress": "parallel_begin",
                                "job": job,
                                "event": event["i"],
                                "elapsed_s": time.monotonic() - self.started,
                                "cases": len(event["jobs"]),
                            }
                        ),
                        flush=True,
                    )
                    self.release_boundary(job)
                    if self.workers == 1:
                        for child in event["jobs"]:
                            self.run_job(child)
                    else:
                        with ThreadPoolExecutor(max_workers=self.workers) as pool:
                            list(pool.map(self.run_job, event["jobs"]))
                else:
                    raise ValueError("unknown event kind")
            except Exception as exc:
                raise RuntimeError(
                    f"job={job} call={event['i']} kind={event.get('op', event['kind'])}: {exc}"
                ) from exc
            if job == "main" and (
                event["phase"] != last_phase
                or event["phase"] == "run"
                and event["i"] % 25 == 0
            ):
                print(
                    dumps(
                        {
                            "progress": "event_complete",
                            "job": job,
                            "event": event["i"],
                            "phase": event["phase"],
                            "elapsed_s": time.monotonic() - self.started,
                        }
                    ),
                    flush=True,
                )
            last_phase = event["phase"]

        self.release_boundary(job)

    def run(self):
        self.started = time.monotonic()
        build = provenance(self.core, self.expected)
        for field in ("python", "gil"):
            if build[field] != self.manifest["provenance"][field]:
                raise ValueError(f"interpreter setting differs from recording: {field}")
        from controls import small_calls, lookup_preflight

        small_calls()
        preflight = lookup_preflight()
        print(dumps({"preflight": preflight, "tiny_stream": "PASS"}), flush=True)
        self.run_job("main")
        if self.coverage != Counter(self.manifest["coverage"]) or self.executed != set(
            self.manifest["counts"]
        ):
            raise ValueError("incomplete call/job coverage")
        if set(self.slots) != set(self.manifest["retained"]):
            raise ValueError("retained object inventory differs")
        if TOOL_HASHES != {
            p.name: digest(p) for p in Path(__file__).parent.glob("*.py")
        }:
            raise ValueError("tool source changed during replay")
        self.streams.close()
        size = sum(p.stat().st_size for p in self.out.rglob("*") if p.is_file())
        if size > 20 * 1024**3:
            raise ValueError("observation output exceeds 20 GiB")
        report = {
            "format": 1,
            "preflight": {**preflight, "tiny_stream": "PASS"},
            "tool_sha256": TOOL_HASHES,
            "recording": digest(self.recording / "manifest.json"),
            "provenance": build,
            "workers": self.workers,
            "lane": self.lane,
            "shard": [self.shard, self.shards],
            "state_sampling": self.state_sampling,
            "state_jobs": sorted(self.state_jobs),
            "state_coverage": dict(self.state_coverage),
            "coverage": dict(self.coverage),
            "observed": dict(self.observed),
            "jobs": dict(self.streams.counts),
            "releases": dict(self.released),
            "bytes": size,
            "gc": gc.get_stats(),
            "summary": "complete checked observation stream",
        }
        (self.out / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
        print(
            dumps(
                {
                    "summary": report["summary"],
                    "calls": sum(self.coverage.values()),
                    "observed": sum(self.observed.values()),
                    "bytes": size,
                }
            ),
            flush=True,
        )


def compare(left, right):
    left, right = Path(left), Path(right)
    a, b = (json.loads((p / "manifest.json").read_text()) for p in (left, right))
    for key in (
        "format",
        "preflight",
        "tool_sha256",
        "recording",
        "workers",
        "lane",
        "shard",
        "state_sampling",
        "state_jobs",
        "state_coverage",
        "coverage",
        "observed",
        "jobs",
        "releases",
    ):
        if a[key] != b[key]:
            raise ValueError(
                f"comparison manifest mismatch: {key}: {a[key]!r} != {b[key]!r}"
            )
    if not a["observed"] or not a["jobs"]:
        raise ValueError("empty comparison")
    calls = 0
    lb, rb = Blobs(left / "payloads"), Blobs(right / "payloads")
    for job in sorted(a["jobs"]):
        for x, y in itertools.zip_longest(
            read_stream(left, job), read_stream(right, job)
        ):
            if x is None or y is None:
                raise ValueError(
                    f"job={job} call={(x or y).get('call')} field=$stream reference={x!r} candidate={y!r}"
                )
            if difference := first_difference(x, y, lb, rb):
                field, expected, actual = difference
                raise ValueError(
                    f"job={job} call={x.get('call')} field={field} reference={expected!r} candidate={actual!r}"
                )
            calls += 1
    print(
        dumps(
            {
                "summary": "PASS exact observations",
                "events": calls,
                "jobs": len(a["jobs"]),
            }
        ),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected-prefix")
    parser.add_argument("--expected-sha256")
    parser.add_argument("--expected-version")
    sub = parser.add_subparsers(dest="command", required=True)
    rec = sub.add_parser("record")
    rec.add_argument("workload", choices=("dcbb96", "nsfnet", "backbone"))
    rec.add_argument("out", type=Path)
    check = sub.add_parser("check")
    check.add_argument("recording", type=Path)
    check.add_argument("out", type=Path)
    check.add_argument("--workers", type=int, required=True)
    check.add_argument("--lane", choices=("a", "b"), required=True)
    check.add_argument("--shard", default="0/1")
    check.add_argument("--state-sample", action="store_true")
    comp = sub.add_parser("compare")
    comp.add_argument("reference", type=Path)
    comp.add_argument("candidate", type=Path)
    args = parser.parse_args()
    if args.command == "compare":
        compare(args.reference, args.candidate)
    else:
        import netgraph_core as core

        expected = None
        if any((args.expected_prefix, args.expected_sha256, args.expected_version)):
            if not all(
                (args.expected_prefix, args.expected_sha256, args.expected_version)
            ):
                raise ValueError(
                    "all three expected-build arguments are required together"
                )
            expected = {
                "prefix": args.expected_prefix,
                "sha256": args.expected_sha256,
                "version": args.expected_version,
            }
        provenance(core, expected)
        if args.command == "record":
            record(core, args.workload, args.out)
        else:
            if args.workers < 1:
                raise ValueError("workers must be positive")
            Replay(
                core,
                args.recording,
                args.out,
                args.workers,
                args.lane,
                args.shard,
                expected,
                args.state_sample,
            ).run()


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(
            dumps({"summary": "FAIL", "error": str(exc)}), file=sys.stderr, flush=True
        )
        raise
