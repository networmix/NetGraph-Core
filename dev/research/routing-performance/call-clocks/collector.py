"""Public-boundary diagnostic clocks, with retained native-object release scopes."""

from array import array
from collections import Counter
from contextlib import contextmanager
import threading
import time
import weakref

from mach_clock import Registry

WIDTH = 17


class Collector:
    def __init__(self, core, api):
        self.registry = Registry()
        self.api = api
        self.types = {
            v
            for v in vars(core).values()
            if type(v).__name__ == "pybind11_type" and not hasattr(v, "__members__")
        }
        self.events, self.vectors = array("q"), array("q")
        self.names = list(api) + ["CoreObjects.release", "calibration.null"]
        self.ids = {name: i for i, name in enumerate(self.names)}
        self.local = threading.local()
        self.phase = 0
        self.releases = Counter()
        self.active = False

    def retain(self, value):
        scopes = getattr(self.local, "holds", [])
        if not scopes:
            return
        if type(value) in self.types:
            scopes[-1][id(value)] = value
        elif isinstance(value, (tuple, list)):
            for item in value:
                if type(item) in self.types:
                    scopes[-1][id(item)] = item

    def wrap(self, name, original):
        def wrapped(*args, **kwargs):
            if not self.active or getattr(self.local, "inside", False):
                return original(*args, **kwargs)
            self.local.inside = True
            try:
                result = self.measure(name, original, args, kwargs)
                self.retain(args[0] if name.endswith(".__init__") else result)
                return result
            finally:
                self.local.inside = False

        return wrapped

    def measure(self, name, original, args=(), kwargs=None):
        start = self.registry.snapshot()
        w0, c0 = time.perf_counter_ns(), time.thread_time_ns()
        try:
            return original(*args, **(kwargs or {}))
        finally:
            c1, w1 = time.thread_time_ns(), time.perf_counter_ns()
            stop = self.registry.snapshot()
            with self.registry.lock:
                self.record(name, w0, w1, c0, c1, start, stop)

    def record(self, name, w0, w1, c0, c1, start, stop):
        offset0 = len(self.vectors)
        self.vectors.extend(start[5])
        offset1 = len(self.vectors)
        self.vectors.extend(stop[5])
        self.events.extend(
            (
                self.ids[name],
                self.phase,
                self.registry.local.row["index"],
                w0,
                w1,
                c0,
                c1,
                start[2],
                start[3],
                start[4],
                stop[2],
                stop[3],
                stop[4],
                offset0,
                len(start[5]),
                offset1,
                len(stop[5]),
            )
        )

    @contextmanager
    def scope(self, label):
        if not self.active:
            yield
            return
        holds = getattr(self.local, "holds", None)
        if holds is None:
            holds = self.local.holds = []
        hold = {}
        holds.append(hold)
        try:
            yield
        finally:
            holds.pop()
            refs = [
                (type(value).__name__, weakref.ref(value)) for value in hold.values()
            ]
            self.measure("CoreObjects.release", hold.clear)
            self.releases.update(kind for kind, ref in refs if ref() is None)

    def save(self, out):
        with (out / "events.i64").open("xb") as stream:
            self.events.tofile(stream)
        with (out / "python-vectors.i64").open("xb") as stream:
            self.vectors.tofile(stream)
        return {
            "names": self.names,
            "event_width": WIDTH,
            "events": len(self.events) // WIDTH,
            "byteorder": __import__("sys").byteorder,
            "threads": [
                {k: v for k, v in row.items() if k != "events"}
                for row in self.registry.rows
            ],
            "destroyed_wrappers": dict(self.releases),
        }


def install(core, api, factory):
    for name, original in api.items():
        wrapped = factory(name, original)
        if "." not in name:
            setattr(core, name, wrapped)
            continue
        cls, attr = name.split(".", 1)
        kind = getattr(core, cls)
        member = vars(kind)[attr]
        if isinstance(member, property):
            wrapped = property(wrapped, member.fset, member.fdel, member.__doc__)
        elif isinstance(member, staticmethod):
            wrapped = staticmethod(wrapped)
        setattr(kind, attr, wrapped)
