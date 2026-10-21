"""Parent-enforced progress checks with an explicitly selected native extension."""

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import _netgraph_core

ROOT = Path(__file__).resolve().parents[2]
PREAMBLE = r"""
import gc
import importlib.util
import os
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pytest
path = os.environ['CORE_CONCURRENCY_BINARY']
spec = importlib.util.spec_from_file_location('_netgraph_core', path)
ext = importlib.util.module_from_spec(spec)
sys.modules['_netgraph_core'] = ext
spec.loader.exec_module(ext)
assert os.path.realpath(ext.__file__) == os.path.realpath(path)
import netgraph_core as ngc

_probes, _threads = [], []
if hasattr(ext, '_ConcurrencyProbe'):
    _native_probe = ext._ConcurrencyProbe
    def new_probe(*args, **kwargs):
        probe = _native_probe(*args, **kwargs)
        _probes.append(probe)
        return probe

def graph(capacity=10):
    return ngc.Graph(2, [0], [1], [capacity], [3], 1.0)

def policy(ledger, *, cache=None, lossy=False):
    config = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(), admission=ngc.Admission.LOSSY if lossy else ngc.Admission.LOSSLESS)
    return ngc.FlowPolicy(ledger, 0, 1, config, cache=cache)

def start(function, probe=None):
    values, errors = [], []
    def call():
        try:
            values.append(probe.run(function) if probe else function())
        except BaseException as error:
            print('CHILD: worker error', type(error).__name__, str(error), flush=True)
            errors.append(error)
    thread = threading.Thread(target=call, daemon=True)
    _threads.append(thread)
    thread.start()
    return thread, values, errors

def finish(task):
    thread, values, errors = task
    thread.join(5)
    assert not thread.is_alive(), 'native thread did not finish'
    if errors:
        raise errors[0]
    assert len(values) == 1
    return values[0]

print('CHILD: selected extension', path, flush=True)
"""


def run_python(source, *, seconds=25, cleanup=True):
    env = dict(os.environ)
    env["CORE_CONCURRENCY_BINARY"] = str(Path(_netgraph_core.__file__).resolve())
    env["PYTHONPATH"] = os.pathsep.join([str(ROOT / "python"), str(ROOT / "tests/py")])
    # Assertion failures under a plant must also release native rendezvous and
    # join workers. The parent deadline remains the backstop for a real deadlock.
    if cleanup:
        program = (
            PREAMBLE + "\ntry:\n" + textwrap.indent(textwrap.dedent(source), "    ")
        )
        program += "\nfinally:\n    for probe in _probes:\n        probe.release()\n    for thread in _threads:\n        thread.join(5)\n"
    else:
        # Interpreter-finalization tests deliberately leave the daemon active.
        program = PREAMBLE + "\n" + textwrap.dedent(source)
    program += "\nprint('CHILD: COMPLETE', flush=True)\n"
    result = subprocess.run(
        [sys.executable, "-B", "-u", "-c", program],
        env=env,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=seconds,
        check=False,
    )
    assert result.returncode == 0, f"child exit={result.returncode}\n{result.stdout}"
    if cleanup:
        assert result.stdout.rstrip().endswith("CHILD: COMPLETE"), result.stdout
    else:
        assert "CHILD: COMPLETE" in result.stdout.splitlines(), result.stdout
    return result.stdout


def run_test_in_subprocess(path, function):
    run_python(
        f"""
        spec = importlib.util.spec_from_file_location('bounded_test', {str(Path(path).resolve())!r})
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        print('CHILD: entering {function}', flush=True)
        getattr(module, {function!r})()
        print('CHILD: threads joined and outcomes checked', flush=True)
        """
    )
