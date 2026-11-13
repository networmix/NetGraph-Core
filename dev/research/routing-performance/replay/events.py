"""Record public Core calls and consumer case scheduling without editing NetGraph."""

import gc
import gzip
import itertools
import json
import sys
import threading
import weakref
from collections import Counter
from pathlib import Path

import numpy as np

from values import Blobs, FIELDS, Values, digest, dumps

TOOL_HASHES = {p.name: digest(p) for p in Path(__file__).parent.glob("*.py")}

BASE = Path("/Users/networmix/ws/work/2026-11-core-routing-research.noindex")
REFERENCE = BASE / "scratch/consumer/ref/venv"
PYPI = BASE / "worktrees/netgraph-orient/venv"
NATIVE_HASHES = {
    str(REFERENCE): "208487c7aed158fb418791c22e23d09332533595053c0e4810923b31244f7c84",
    str(PYPI): "6db37900d83800b0b5c227104889f47259c5761dd640c1d0d8684fac61009ff1",
}


def provenance(core, expected=None):
    import _netgraph_core

    native = Path(_netgraph_core.__file__).resolve()
    package = Path(core.__file__).resolve()
    prefix = str(Path(sys.prefix).resolve())
    if expected is None:
        if prefix not in NATIVE_HASHES:
            raise ValueError(
                "supply explicit expected prefix, native SHA-256 and version for this build"
            )
        expected = {
            "prefix": prefix,
            "sha256": NATIVE_HASHES[prefix],
            "version": "0.13.0",
        }
    if (
        prefix != str(Path(expected["prefix"]).resolve())
        or not native.is_relative_to(prefix)
        or not package.is_relative_to(prefix)
    ):
        raise ValueError("Core import is outside the declared wheel environment")
    if digest(native) != expected["sha256"] or core.__version__ != expected["version"]:
        raise ValueError("Core build identity mismatch")
    return {
        "python": sys.version,
        "prefix": prefix,
        "gil": sys._is_gil_enabled(),
        "package": str(package),
        "native": str(native),
        "sha256": digest(native),
        "version": core.__version__,
    }


def public_api(core):
    """Capture unbound originals; property readers and constructor bodies are calls."""
    api = {}
    for name, obj in list(vars(core).items()):
        if name.startswith("_"):
            continue
        if type(obj).__name__ == "builtin_function_or_method":
            api[name] = obj
        elif type(obj).__name__ == "pybind11_type":
            for attr, member in list(vars(obj).items()):
                if attr.startswith("_") and attr != "__init__":
                    continue
                if isinstance(member, property):
                    allowed = (
                        ("name", "value")
                        if hasattr(obj, "__members__")
                        else FIELDS.get(name, ())
                    )
                    if attr not in allowed:
                        raise ValueError(f"unhandled public field {name}.{attr}")
                    api[name + "." + attr] = member.fget
                elif isinstance(member, staticmethod):
                    api[name + "." + attr] = member.__func__
                elif type(member).__name__ == "instancemethod":
                    api[name + "." + attr] = getattr(obj, attr)
    return api


class Streams:
    def __init__(self, root, limit):
        self.root, self.limit = Path(root), limit
        self.files, self.locks, self.counts = {}, {}, Counter()
        self.lock = threading.RLock()
        self.bytes = 0

    def charge(self, count):
        with self.lock:
            self.bytes += count
            if self.bytes > self.limit:
                raise ValueError("stream and payload byte budget exceeded")

    def write(self, job, record):
        with self.lock:
            if job not in self.files:
                self.files[job] = gzip.open(
                    self.root / (job + ".jsonl.gz"), "xt", encoding="utf-8"
                )
                self.locks[job] = threading.RLock()
                self.charge(256)
        with self.locks[job]:
            record = {"i": self.counts[job], **record}
            line = dumps(record) + "\n"
            self.charge(len(line))
            self.files[job].write(line)
            self.counts[job] += 1
        return record["i"]

    def close(self):
        for job, stream in self.files.items():
            stream.write(dumps({"footer": self.counts[job]}) + "\n")
            stream.close()


def read_stream(root, job):
    count, ended = 0, False
    with gzip.open(Path(root) / (job + ".jsonl.gz"), "rt", encoding="utf-8") as stream:
        for line in stream:
            record = json.loads(line)
            if ended:
                raise ValueError("data follows footer")
            if "footer" in record:
                if record != {"footer": count}:
                    raise ValueError("stream footer count mismatch")
                ended = True
            else:
                if record.get("i") != count:
                    raise ValueError("stream event index mismatch")
                count += 1
                yield record
    if not ended or count == 0:
        raise ValueError("empty or incomplete event stream")


class Recorder(Values):
    def __init__(self, core, out):
        self.out = Path(out)
        self.out.mkdir()
        self.streams = Streams(self.out, 4 * 1024**3)
        super().__init__(
            core, Blobs(self.out / "payloads", write=True, charge=self.streams.charge)
        )
        self.api = public_api(core)
        self.local = threading.local()
        self.phase = "import"
        self.ids, self.refs = {}, {}
        self.ledger_uids = {}
        self.native_arrays = set()
        self.counter, self.ordinals = itertools.count(), itertools.count()
        self.groups, self.coverage = [], Counter()
        self.active = True
        self.lock = threading.RLock()

    @property
    def job(self):
        return getattr(self.local, "job", "main")

    def emit(self, value):
        return self.streams.write(self.job, {"phase": self.phase, **value})

    def register(self, obj):
        key = id(obj)
        old = self.ids.get(key)
        if old is not None and self.refs[old]() is obj:
            return old
        oid = f"{self.job}:{next(self.counter)}"
        self.ids[key] = oid

        def release(ref):
            with self.lock:
                if self.ids.get(key) == oid:
                    del self.ids[key]
                self.refs.pop(oid, None)
                if self.active:
                    self.emit(
                        {
                            "kind": "release",
                            "oid": oid,
                            "type": type_name,
                            "ordinal": next(self.ordinals),
                        }
                    )

        type_name = type(obj).__name__
        self.refs[oid] = weakref.ref(obj, release)
        return oid

    def encode(self, obj, output=False):
        kind = type(obj).__name__
        if kind in FIELDS:
            oid = self.ids.get(id(obj))
            if oid is None or self.refs[oid]() is not obj:
                if not output:
                    raise ValueError(f"native input has no producer: {kind}")
                oid = self.register(obj)
            if kind == "FlowLedger":
                self.ledger_uids[self.api["FlowLedger.uid"](obj)] = oid
            return {"$object": oid, "type": kind}
        if isinstance(obj, np.ndarray):
            oid = self.register(obj)
            if output:
                self.native_arrays.add(oid)
            owner = obj.base
            while isinstance(owner, np.ndarray):
                owner = owner.base
            return {
                "$numpy": self.blobs.put(obj),
                "id": oid,
                "strides": list(obj.strides),
                "readonly": not obj.flags.writeable,
                "owner": self.ids.get(id(owner)),
                "produced": oid in self.native_arrays,
            }
        return self.atom(obj, lambda v: self.encode(v, output))

    def wrap(self, name, original):
        def wrapped(*args, **kwargs):
            if getattr(self.local, "inside", False) or not self.active:
                return original(*args, **kwargs)
            constructor = name.endswith(".__init__")
            self.local.inside = True
            try:
                target = self.register(args[0]) if constructor else None
                inputs = self.encode(args[1:] if constructor else args)
                keywords = self.encode(kwargs)
                if name in ("FlowLedger.record", "FlowPolicy.record"):
                    flow_id = args[1]
                    if (
                        type(flow_id) is tuple
                        and len(flow_id) == 2
                        and type(flow_id[0]) is int
                        and flow_id[0] in self.ledger_uids
                    ):
                        inputs["$tuple"][1] = {
                            "$flowid": [
                                self.ledger_uids[flow_id[0]],
                                self.encode(flow_id[1]),
                            ]
                        }

            finally:
                self.local.inside = False
            enter = next(self.ordinals)
            error, result = None, None
            try:
                result = original(*args, **kwargs)
                return result
            except Exception as exc:
                error = {
                    "class": type(exc).__module__ + "." + type(exc).__name__,
                    "text": str(exc),
                }
                raise
            finally:
                stop = next(self.ordinals)
                self.local.inside = True
                try:
                    returned = (
                        self.encode(args[0] if constructor else result, output=True)
                        if error is None
                        else None
                    )
                    self.emit(
                        {
                            "kind": "call",
                            "op": name,
                            "target": target,
                            "args": inputs,
                            "kwargs": keywords,
                            "result": returned,
                            "error": error,
                            "enter": enter,
                            "exit": stop,
                        }
                    )
                    self.coverage[name] += 1
                finally:
                    self.local.inside = False

        return wrapped

    def install(self):
        for name, original in self.api.items():
            wrapper = self.wrap(name, original)
            if "." not in name:
                setattr(self.core, name, wrapper)
                continue
            cls, attr = name.split(".")
            kind = getattr(self.core, cls)
            member = vars(kind)[attr]
            if isinstance(member, property):
                wrapper = property(wrapper, member.fset, member.fdel, member.__doc__)
            elif isinstance(member, staticmethod):
                wrapper = staticmethod(wrapper)
            setattr(kind, attr, wrapper)

    def scheduler(self):
        import ngraph.analysis.failures as failures

        original = failures.schedule_cases

        def schedule(cases, solve, workers):
            group = f"g{len(self.groups):04d}"
            jobs = [f"{group}-{i:06d}" for i in range(len(cases))]
            self.groups.append({"id": group, "jobs": jobs, "workers": workers})
            self.emit({"kind": "parallel", "jobs": jobs, "workers": workers})
            positions = {id(case): job for case, job in zip(cases, jobs, strict=True)}

            def run(case):
                old = self.job
                self.local.job = positions[id(case)]
                try:
                    return solve(case)
                finally:
                    self.local.job = old

            return original(cases, run, workers)

        for name, module in list(sys.modules.items()):
            if name.startswith("ngraph.") and module is not None:
                for attr, value in list(vars(module).items()):
                    if value is original:
                        setattr(module, attr, schedule)

    def finish(self, metadata):
        if TOOL_HASHES != {
            p.name: digest(p) for p in Path(__file__).parent.glob("*.py")
        }:
            raise ValueError("tool source changed during recording")
        self.active = False
        self.streams.close()
        metadata.update(
            {
                "format": 1,
                "provenance": provenance(self.core),
                "groups": self.groups,
                "counts": dict(self.streams.counts),
                "coverage": dict(self.coverage),
                "retained": {
                    oid: type(ref()).__name__
                    for oid, ref in self.refs.items()
                    if ref() is not None
                },
                "fields": FIELDS,
                "tool_sha256": TOOL_HASHES,
                "summary": "complete full-input recording",
            }
        )
        metadata["bytes"] = sum(
            p.stat().st_size for p in self.out.rglob("*") if p.is_file()
        )
        if metadata["bytes"] > 4 * 1024**3:
            raise ValueError("recording exceeds 4 GiB")
        (self.out / "manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
        print(
            dumps(
                {
                    "summary": metadata["summary"],
                    "calls": sum(self.coverage.values()),
                    "bytes": metadata["bytes"],
                }
            ),
            flush=True,
        )


def record(core, workload, out):
    if Path(sys.prefix).resolve() != REFERENCE:
        raise ValueError("recordings must use the frozen reference wheel")
    recorder = Recorder(core, out)
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
    manifest = json.loads(
        (BASE / "scratch/orient-add/workloads/manifest.json").read_text()
    )
    source = Path(manifest[workload]["path"])
    if digest(source) != manifest[workload]["sha256"]:
        raise ValueError("workload input mismatch")
    recorder.phase = "load"
    scenario = Scenario.from_yaml(source.read_text())
    recorder.phase = "prepare"
    prepared = ngraph.prepare(scenario)
    recorder.phase = "run"
    result = ngraph.run(prepared)
    recorder.phase = "export"
    with (Path(out) / "export.json").open("x") as stream:
        ngraph.write_json(result, stream)
    recorder.phase = "cleanup"
    del result, prepared, scenario
    gc.collect()
    recorder.finish(
        {
            "workload": workload,
            "source": str(source),
            "input_sha256": digest(source),
            "consumer_commit": "8a641a1caa0b8f4b7d8c4df0aecf8625310cdf47",
            "workers": {"dcbb96": 8, "nsfnet": 8, "backbone": 28}[workload],
        }
    )
