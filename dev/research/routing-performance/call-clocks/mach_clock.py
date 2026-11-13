"""Read registered Python threads' Darwin CPU counters without releasing the GIL."""

import ctypes as C
import threading
import time


class TimeValue(C.Structure):
    _fields_ = [("seconds", C.c_int), ("microseconds", C.c_int)]


class BasicInfo(C.Structure):
    _fields_ = [
        ("user", TimeValue),
        ("system", TimeValue),
        ("cpu_usage", C.c_int),
        ("policy", C.c_int),
        ("run_state", C.c_int),
        ("flags", C.c_int),
        ("suspend_count", C.c_int),
        ("sleep_time", C.c_int),
    ]


lib = C.PyDLL("/usr/lib/libSystem.B.dylib")
lib.pthread_self.restype = C.c_void_p
lib.pthread_mach_thread_np.argtypes = [C.c_void_p]
lib.pthread_mach_thread_np.restype = C.c_uint
lib.thread_info.argtypes = [C.c_uint, C.c_int, C.c_void_p, C.POINTER(C.c_uint)]
lib.thread_info.restype = C.c_int


def own_port():
    return lib.pthread_mach_thread_np(lib.pthread_self())


def cpu_ns(port):
    info = BasicInfo()
    count = C.c_uint(C.sizeof(info) // C.sizeof(C.c_int))
    code = lib.thread_info(port, 3, C.byref(info), C.byref(count))
    if code or count.value != C.sizeof(info) // C.sizeof(C.c_int):
        raise RuntimeError(f"thread_info failed: code={code}, count={count.value}")
    return (info.user.seconds + info.system.seconds) * 1_000_000_000 + (
        info.user.microseconds + info.system.microseconds
    ) * 1000


class Registry:
    def __init__(self):
        self.rows = []
        self.local = threading.local()
        self.original = None
        self.live = {}
        self.retired = 0
        self.lock = threading.RLock()

    def begin(self, thread=None):
        with self.lock:
            return self._begin(thread)

    def _begin(self, thread):
        port = own_port()
        row = {
            "index": len(self.rows),
            "name": (thread or threading.current_thread()).name,
            "ident": threading.get_ident(),
            "port": port,
            "start_wall": time.perf_counter_ns(),
            "start_self": time.thread_time_ns(),
            "start_mach": cpu_ns(port),
            "end_mach": None,
            "events": [],
        }
        self.rows.append(row)
        self.live[row["index"]] = row
        self.local.row = row
        return row

    def end(self, row):
        with self.lock:
            row["end_mach"] = cpu_ns(row["port"])
            row["end_self"] = time.thread_time_ns()
            row["end_wall"] = time.perf_counter_ns()
            self.retired += row["end_mach"] - row["start_mach"]
            del self.live[row["index"]]

    def snapshot(self):
        with self.lock:
            return self._snapshot()

    def _snapshot(self):
        wall0, process0 = time.perf_counter_ns(), time.process_time_ns()
        values, total = [], self.retired
        for row in self.live.values():
            value = cpu_ns(row["port"])
            values.extend((row["index"], value))
            total += value - row["start_mach"]
        process1, wall1 = time.process_time_ns(), time.perf_counter_ns()
        return [wall0, wall1, process0, process1, total, values]

    def install(self):
        existing = threading.enumerate()
        if existing != [threading.main_thread()]:
            raise RuntimeError(
                f"unregistered Python threads before installation: {existing}"
            )
        self.begin()
        self.original = threading.Thread._bootstrap_inner
        registry = self

        def bootstrap(thread):
            # Thread.name is available before bootstrap assigns its ident.
            row = registry.begin(thread)
            try:
                return registry.original(thread)
            finally:
                registry.end(row)

        threading.Thread._bootstrap_inner = bootstrap

    def uninstall(self):
        threading.Thread._bootstrap_inner = self.original
        self.end(self.rows[0])
