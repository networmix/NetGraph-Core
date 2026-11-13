"""Load read-only frozen T4 tools by their actual paths."""

import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
T4 = ROOT / "tasks/clocks/snapshot-r3/task"


def load(name):
    key = "split_t4_" + name
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, T4 / (name + ".py"))
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]
