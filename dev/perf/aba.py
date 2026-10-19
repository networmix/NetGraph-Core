#!/usr/bin/env python3
"""A/B/A comparison of two Core timing harness binaries (README.md in this directory).

    venv/bin/python dev/perf/aba.py --old-bin REF --new-bin CAND [options]

A (``--old-bin``) is the reference build and B (``--new-bin``) the candidate, each the rows of
``aba_rows.hpp`` compiled against one Core (``build.sh``). Every selected row runs one warm-up
block per side, then A1, B, A2 triples, each block a fresh harness process; the row's statistics
and verdicts follow, and a row whose blocks do not all report one and the same digest fails.
README.md, "Protocol" and "Inputs", states every rule; test_tool.py tests each one.

The run prints its plan, one line per row, and writes ``report.json``, ``summary.md``,
``blocks.csv`` and ``meta.txt``. Exit status: 0 the run completed, whatever its verdict (the
summary's first lines and ``report.json`` give the verdict and the conformance); 1 it aborted; 2 it
was refused before the first block; 130 it was interrupted (Ctrl-C) before the first block.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import datetime as dt
import fnmatch
import hashlib
import io
import json
import math
import os
import platform
import random
import re
import signal
import statistics
import subprocess
import sys
import tempfile
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PYTHON_MIN = (3, 10)
PROTOCOL_BLOCK_SECONDS = 2.0  # a timed block lasts at least 2 s
PROTOCOL_TRIPLES = 7  # at least seven A-B-A triples
RSS_FLOOR_MIB = 8.0  # memory allowance max(1.10 * reference, reference + 8 MiB)
# Child-environment variables that change the allocator: recorded, never set or cleared (an
# override such as MallocLargeCache=0 would hide the allocator-cache levels and with them a
# candidate's genuine peak, so a run with one is not conforming).
ALLOCATOR_PREFIXES = ("Malloc", "MALLOC_")
ALLOCATOR_NAMES = ("DYLD_INSERT_LIBRARIES", "LD_PRELOAD", "GLIBC_TUNABLES")
# Orchestrator and driver sources: hashed into meta.txt at the start of a run.
SOURCE_FILES = ("aba.py", "harness.cpp", "aba_rows.hpp", "aba_api.hpp")
# The names the tool accepts (README "Names"). They become parts of file names
# (blocks/<row>.<label>.<mode>.stderr) and fields of meta.txt's tab-separated lines: lower case
# only (two names that differ in case would share a log on a case-insensitive file system), and
# no dot in a label (the log name then splits back into its row, label and mode).
ROW_NAME = re.compile(r"[a-z0-9._-]{1,64}")  # a row's name or its group
LABEL = re.compile(r"[a-z0-9_-]{1,32}")
DIGEST = re.compile(r"[0-9a-f]{16}")  # as the harness prints it
ROWS_TEXT = frozenset("abcdefghijklmnopqrstuvwxyz0123456789._-*?[]!, ")


class Abort(RuntimeError):
    """The run cannot go on once its blocks have begun: the message is the abort's reason, which
    names the block (side, mode and row) or the row it comes from where there is one."""


class ProvenanceError(Abort):
    """A side's binary changed or became unreadable during the run."""


class ResponseError(Abort):
    """A harness response that is not what its command answers (``response_records``)."""


class Refused(Exception):
    """An input refused before the first block (exit 2)."""


def sha256_file(path: Path | str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


SHA256_HEX = frozenset("0123456789abcdef")


def is_sha256(text) -> bool:
    return isinstance(text, str) and len(text) == 64 and set(text) <= SHA256_HEX


def unique_keys(pairs: list[tuple[str, object]]) -> dict:
    out = {}
    for key, value in pairs:
        if key in out:
            raise ValueError(f"key {clip(key)} twice in one object")
        out[key] = value
    return out


def load_json(text: str):
    """JSON text with no key twice in an object: ``json.loads`` keeps the last value of a
    repeated key, so one record could carry two values of a field and be read by the second.
    ValueError otherwise (``json.JSONDecodeError`` is one), also for values nested too deeply
    for Python's reader (it raises RecursionError)."""
    try:
        return json.loads(text, object_pairs_hook=unique_keys)
    except RecursionError as exc:
        raise ValueError(f"JSON nested too deeply ({exc})") from exc


def clip(value, limit: int = 120) -> str:
    """An input value for a message: an array or an object by its kind only (one nested deeply
    cannot even be printed: Python's printer and JSON writer recurse, and fail where the reader
    did not), anything else by its ``repr`` cut to ``limit`` characters."""
    if isinstance(value, (list, dict)):
        return "an object" if isinstance(value, dict) else "an array"
    text = repr(value)
    return text if len(text) <= limit else text[: limit - 3] + "..."


def shown(rec: dict) -> str:
    """A record for a message: JSON cut at 200 characters, a nested value by its kind only
    (``clip``)."""
    flat = {
        k: f"<{clip(v)}>" if isinstance(v, (list, dict)) else v for k, v in rec.items()
    }
    return json.dumps(flat)[:200]


def unknown_keys(rec: dict, known) -> str | None:
    """Why ``rec`` has a key outside ``known`` (the first five named), or None: an input carries
    only the keys its rule lists, so nothing unchecked reaches the report."""
    unknown = [k for k in rec if k not in known]
    if not unknown:
        return None
    names = ", ".join(clip(k) for k in unknown[:5]) + (
        ", ..." if len(unknown) > 5 else ""
    )
    return f"unknown key{'s' if len(unknown) > 1 else ''} {names}"


def response_records(out: bytes, what: str) -> list[dict]:
    """A harness's complete standard output: UTF-8 text of whole lines, each ended by a newline
    and each one JSON object (``load_json``). Nothing else is accepted (no blank line, no text and
    no partial record before or after the records), so a usable line is never picked out of a
    broken response. ResponseError otherwise."""
    try:
        text = out.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ResponseError(f"{what}: the output is not UTF-8 text ({exc})") from exc
    if text and not text.endswith("\n"):
        raise ResponseError(
            f"{what}: malformed output line (not ended by a newline): {text[-200:]!r}"
        )
    lines = text.split("\n")[:-1]
    records = []
    for i, line in enumerate(lines, 1):
        try:
            rec = load_json(line)
        except ValueError as exc:
            rec, why = None, str(exc)
        else:
            why = "not a JSON object"
        if not isinstance(rec, dict):
            raise ResponseError(
                f"{what}: malformed output line ({why}; line {i} of {len(lines)}): "
                f"{line[:200]}"
            )
        records.append(rec)
    return records


def one_record(out: bytes, what: str, row: str) -> dict:
    """The response of a command that answers for one row: exactly one record
    (``response_records``), and it names that row. ResponseError otherwise."""
    records = response_records(out, what)
    if len(records) != 1:
        raise ResponseError(
            f"{what}: {len(records)} output records, exactly one expected"
        )
    if records[0].get("row") != row:
        raise ResponseError(
            f"{what}: the record names row {clip(records[0].get('row'))}, not {row!r}"
        )
    return records[0]


# The keys of an aba-list record: whether each is required (every one a string).
ROW_LIST_FIELDS = {"row": True, "group": True, "desc": False, "queries": False}


def row_list(out: bytes, what: str) -> dict[str, dict]:
    """``aba-list``'s complete response (``response_records``): at least one row, every record
    with a string ``row`` and a string ``group``, each a name of ``ROW_NAME`` (``desc`` and
    ``queries`` strings when present) and no other key, no row twice. ResponseError otherwise."""
    rows: dict[str, dict] = {}
    for rec in response_records(out, what):
        bad = [
            k
            for k, required in ROW_LIST_FIELDS.items()
            if (k in rec or required) and not isinstance(rec.get(k), str)
        ]
        if bad:
            raise ResponseError(
                f"{what}: a row record without a string {', '.join(bad)}: {shown(rec)}"
            )
        for k in ("row", "group"):
            if not ROW_NAME.fullmatch(rec[k]):
                raise ResponseError(
                    f"{what}: {k} {clip(rec[k])} is not 1 to 64 of the characters a-z, 0-9,"
                    f" '.', '_' and '-': {shown(rec)}"
                )
        if unknown := unknown_keys(rec, ROW_LIST_FIELDS):
            raise ResponseError(f"{what}: row {rec['row']}: {unknown}: {shown(rec)}")
        if rec["row"] in rows:
            raise ResponseError(f"{what}: row {rec['row']!r} listed twice")
        rows[rec["row"]] = rec
    if not rows:
        raise ResponseError(f"{what}: no row listed")
    return rows


@dataclass
class Side:
    """One harness binary, its hash at the start of the run and its build record."""

    key: str  # "A" or "B" in meta.txt
    label: str
    binary: str
    sha256: str
    build_record: dict = field(default_factory=dict)

    def cmd(self, mode: str, *args: str) -> list[str]:
        return [self.binary, mode, *args]

    def verify(self, when: str) -> None:
        """ProvenanceError unless the binary still has its hash."""
        try:
            got = sha256_file(self.binary)
        except OSError as exc:
            raise ProvenanceError(
                f"{self.label}: harness {self.binary} unreadable {when} ({exc})"
            ) from exc
        if got != self.sha256:
            raise ProvenanceError(
                f"{self.label}: harness changed {when}: {self.binary} "
                f"(meta.txt {self.sha256[:12]}..., now {got[:12]}...)"
            )


TEXT = (str,)
TEXT_OR_NULL = (str, type(None))
# The build record build.sh writes (README "Build records"): every key with the JSON types of its
# value, a nested object with its own keys. Beyond the types: binary_sha256 and
# core_sources_sha256 are sha256 values, dirty_core_paths is a nonnegative integer, and
# driver_sha256's content is this checkout's driver (same_driver).
BUILD_RECORD = {
    "binary": TEXT,
    "binary_sha256": TEXT,
    "rev": TEXT,
    "commit": TEXT,
    "dirty_core_paths": (int,),
    "core_sources_sha256": TEXT,
    "driver_sha256": (dict,),
    "compiler": TEXT,
    "flags": TEXT,
    "compiler_path": TEXT_OR_NULL,
    "cxx_env": TEXT_OR_NULL,
    "isysroot": TEXT_OR_NULL,
    "sdk": {
        "path": TEXT_OR_NULL,
        "version": TEXT_OR_NULL,
        "sdkroot_env": TEXT_OR_NULL,
        "source": TEXT,
    },
    "deployment_target": {
        "requested": TEXT_OR_NULL,
        "effective_minos": TEXT_OR_NULL,
        "binary_sdk": TEXT_OR_NULL,
        "platform": TEXT_OR_NULL,
        "raw": TEXT_OR_NULL,
        "source": TEXT,
    },
    "built_by": TEXT,
}


def record_problems(rec: dict, schema: dict, prefix: str = "") -> list[str]:
    """Why ``rec`` is not an object of ``schema`` (``BUILD_RECORD``): every key of the schema
    present, no other key, each value of its JSON types (a bool is no integer)."""
    problems = [f"no {prefix}{k}" for k in schema if k not in rec]
    if unknown := unknown_keys(rec, schema):
        problems.append(f"{unknown} in {prefix[:-1]}" if prefix else unknown)
    for k, kind in schema.items():
        if k not in rec:
            continue
        value = rec[k]
        if isinstance(kind, dict):
            if isinstance(value, dict):
                problems += record_problems(value, kind, f"{prefix}{k}.")
            else:
                problems.append(f"{prefix}{k} is {clip(value)}, not an object")
        elif isinstance(value, bool) or not isinstance(value, kind):
            problems.append(f"{prefix}{k} is {clip(value)}, of another type")
    return problems


def build_record_problems(rec: dict) -> list[str]:
    """Why ``rec`` is not a build record as build.sh writes it (``BUILD_RECORD``)."""
    problems = record_problems(rec, BUILD_RECORD)
    for k in ("binary_sha256", "core_sources_sha256"):
        if isinstance(rec.get(k), str) and not is_sha256(rec[k]):
            problems.append(f"{k} {clip(rec[k])} is not a sha256")
    count = rec.get("dirty_core_paths")
    if isinstance(count, int) and not isinstance(count, bool):
        if p := domain_problem("dirty_core_paths", count, (int,), NONNEGATIVE):
            problems.append(p)
    return problems


def make_side(key: str, label: str, binary: str) -> Side:
    """One side, Refused unless ``binary`` is a readable file whose resolved path is printable
    text (``meta.txt`` records it on one line) and its build record (``build.sh`` writes
    ``<binary>.json``), when there is one, is a JSON object whose ``binary_sha256`` is the
    binary's sha256 and which is build.sh's record (``build_record_problems``). The record is
    read beside the resolved path; a ``<binary>.json`` that is not that file (beside a symbolic
    link to the binary) is refused, so a harness given through a link keeps its record and the
    checks of it. A refusal of the record names the side."""
    try:
        path = Path(binary).resolve()  # RuntimeError for a symlink loop before 3.13
        if not str(path).isprintable():
            raise Refused(
                f"{label}: {binary!r}: a path with a character that is not printable"
            )
        if not path.is_file():
            raise Refused(f"{label}: {binary}: no such file")
        sha = sha256_file(path)
        record_path = Path(f"{path}.json")
        has_record = record_path.exists()
        given = Path(binary + ".json")
        if given.exists() and not (has_record and given.samefile(record_path)):
            raise Refused(
                f"{label}: {given}: a build record beside {binary}, not beside the binary it names"
                f" ({path}); keep the record beside the binary"
            )
    except (OSError, RuntimeError) as exc:
        raise Refused(f"{label}: {binary}: cannot read it ({exc})") from exc
    record = {}
    if has_record:
        try:
            record = load_json(record_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise Refused(
                f"{label}: {record_path}: not a readable JSON build record ({exc})"
            ) from exc
        if not isinstance(record, dict) or not is_sha256(record.get("binary_sha256")):
            raise Refused(
                f"{label}: {record_path}: a build record without binary_sha256"
            )
        if record["binary_sha256"] != sha:
            raise Refused(
                f"{label}: {binary}: differs from its build record {record_path}"
            )
        if problems := build_record_problems(record):
            raise Refused(
                f"{label}: {record_path}: not a build record as build.sh writes it"
                f" ({'; '.join(problems)}); rebuild the harness with build.sh"
            )
    return Side(key, label, str(path), sha, record)


def same_driver(side: Side, sources: dict[str, str]) -> None:
    """Refused when the side's binary has a build record whose ``driver_sha256`` (the hashes of
    the driver sources it was compiled from, under ``build.sh``'s names ``dev/perf/<file>``) is
    not that of this checkout's driver sources (``sources``, as ``meta.txt`` records them): a
    harness built before the rows or the adapter changed would run code ``meta.txt`` does not
    name."""
    if not side.build_record:
        return
    have = side.build_record.get("driver_sha256")
    want = {f"dev/perf/{n}": sha for n, sha in sources.items() if n != "aba.py"}
    if have != want:
        files = (
            sorted(k for k in have.keys() | want.keys() if have.get(k) != want.get(k))
            if isinstance(have, dict)
            else ["driver_sha256"]
        )
        raise Refused(
            f"{side.label}: harness {side.binary} was not built from this checkout's driver"
            f" sources (its build record differs in {', '.join(files)}); rebuild it with"
            " build.sh"
        )


def allocator_record(env: dict[str, str]) -> dict:
    """The allocator variables of the harness's environment: whether ``MallocLargeCache`` is
    set, every variable named by ``ALLOCATOR_PREFIXES`` / ``ALLOCATOR_NAMES`` with its value, and
    ``normal`` when there is none."""
    variables = {
        k: env[k]
        for k in sorted(env)
        if k.startswith(ALLOCATOR_PREFIXES) or k in ALLOCATOR_NAMES
    }
    return {
        "MallocLargeCache": {
            "set": "MallocLargeCache" in env,
            "value": env.get("MallocLargeCache"),
        },
        "variables": variables,
        "normal": not variables,
    }


def allocator_line(rec: dict) -> str:
    cache = rec["MallocLargeCache"]
    text = (
        f"MallocLargeCache={cache['value']}"
        if cache["set"]
        else "MallocLargeCache unset"
    )
    others = {k: v for k, v in rec["variables"].items() if k != "MallocLargeCache"}
    text += (
        ", " + ", ".join(f"{k}={v}" for k, v in others.items())
        if others
        else ", no other allocator variable"
    )
    verdict = "normal" if rec["normal"] else "**allocator override: not gate evidence**"
    return f"{text} ({verdict})"


def peak_rss(ru_maxrss: float) -> float:
    """``ru_maxrss`` in bytes (macOS reports bytes, Linux KiB)."""
    return float(ru_maxrss) if sys.platform == "darwin" else float(ru_maxrss) * 1024.0


# The domains of the numbers aba.py reads (block lines, options). Every number must also be
# finite: json.loads reads NaN and Infinity, and NaN compares false with every bound, so a NaN
# duration would pass any check written as a comparison.
POSITIVE, NONNEGATIVE, ANY = "positive", "nonnegative", "any"
NUMBER = (int, float)


def domain_problem(name: str, value, types: tuple, domain: str) -> str | None:
    """Why ``value`` is not a ``name`` of these JSON types (a bool is no number) in this
    domain, or None. A number is finite: a float that is not, and an integer beyond the
    largest float (JSON integers have no bound; arithmetic with a float would raise), are
    refused."""
    if isinstance(value, bool) or not isinstance(value, types):
        return f"missing or mistyped {name}"
    if isinstance(value, str):
        return None
    if isinstance(value, float) and not math.isfinite(value):
        return f"{name} {value!r} is not finite"
    if isinstance(value, int) and abs(value) > sys.float_info.max:
        return f"{name} is not finite (an integer of {len(str(abs(value)))} digits)"
    if domain == POSITIVE and value <= 0:
        return f"{name} {value!r} is not positive"
    if domain == NONNEGATIVE and value < 0:
        return f"{name} {value!r} is negative"
    return None


# The fields of a block line the rules read, with their JSON types and domains (the row a block
# answers for is checked against the requested one: ``one_record``; the digest is 16 lowercase
# hexadecimal digits: ``block_problems``). The iteration count and the
# median iteration time are positive (the ratio and the calibration divide by the median).
# block_s and setup_s may be 0: the harness prints them rounded to microseconds, so a loop or a
# setup shorter than half a microsecond reads 0 (a timed block that short fails the block-length
# rule). The harness's own peak RSS is positive (a process that ran held pages); its children's is
# 0 when it starts none, as the harness never does.
BLOCK_FIELDS = {
    "row": ((str,), ANY),
    "iterations": ((int,), POSITIVE),
    "median_s": (NUMBER, POSITIVE),
    "block_s": (NUMBER, NONNEGATIVE),
    "setup_s": (NUMBER, NONNEGATIVE),
    "digest": ((str,), ANY),
    "rss_self_bytes": (NUMBER, POSITIVE),
    "rss_children_bytes": (NUMBER, NONNEGATIVE),
}
# The harness's other iteration-time statistics: no rule reads them (blocks.csv copies the mean,
# minimum and maximum); checked when present, and 0 is accepted because an iteration shorter
# than one clock tick reads 0.
BLOCK_STATS = ("mean_s", "min_s", "p10_s", "p90_s", "max_s")
# The numbers parse_block turns into floats: all but the iteration count.
FLOAT_FIELDS = (
    *(k for k, (types, _) in BLOCK_FIELDS.items() if types is NUMBER),
    *BLOCK_STATS,
)


def block_problems(rec: dict) -> list[str]:
    """Why ``rec`` is not a block line (``BLOCK_FIELDS``, ``BLOCK_STATS``, no other key)."""
    problems = [
        p
        for k, (types, domain) in BLOCK_FIELDS.items()
        if (p := domain_problem(k, rec.get(k), types, domain))
    ]
    problems += [
        p
        for k in BLOCK_STATS
        if k in rec and (p := domain_problem(k, rec[k], NUMBER, NONNEGATIVE))
    ]
    digest = rec.get("digest")
    if isinstance(digest, str) and not DIGEST.fullmatch(digest):
        problems.append(f"digest {clip(digest)} is not 16 lowercase hexadecimal digits")
    if unknown := unknown_keys(rec, (*BLOCK_FIELDS, *BLOCK_STATS)):
        problems.append(unknown)
    count = rec.get("iterations")
    if type(count) is int and count > MAX_COUNT:
        problems.append(f"iterations above {MAX_COUNT}")
    return problems


def require_finite(what: str, value) -> None:
    """Stop the run (Abort) when a derived value, or any number inside it (lists and
    dicts are searched), is not finite. Every number a block line carries is finite
    (``BLOCK_FIELDS``), but arithmetic on huge ones is not: the mean of two medians near 1e308 is
    infinite, and a ratio against it would read 0, a false pass. Checked before the value is
    recorded or printed, so no verdict rests on it."""
    if isinstance(value, dict):
        for k, v in value.items():
            require_finite(f"{what} {k}", v)
    elif isinstance(value, list):
        for v in value:
            require_finite(what, v)
    elif isinstance(value, float) and not math.isfinite(value):
        raise Abort(
            f"{what} is {value!r}, not finite: the blocks' numbers are outside the range the "
            "statistics can use"
        )


def system_number(what: str, value, *, types=NUMBER, domain=NONNEGATIVE):
    """An OS reading of its declared numeric type and domain, checked before storage."""
    if problem := domain_problem(what, value, types, domain):
        raise Abort(problem)
    return value


def parse_block(out: bytes, what: str, row: str) -> dict:
    """A block's complete response: exactly one record for ``row`` (``one_record``) whose
    fields are of their types and domains, with no other key (``block_problems``); ResponseError
    (the run aborts) otherwise. The times and peaks become floats (the harness prints a peak as
    an integer), so every value derived from them is float arithmetic, where an overflow is an
    infinity that ``require_finite`` catches, never an exception."""
    rec = one_record(out, what, row)
    bad = block_problems(rec)
    if bad:
        raise ResponseError(
            f"{what}: malformed output line ({'; '.join(bad)}): {shown(rec)}"
        )
    rec.update({k: float(rec[k]) for k in FLOAT_FIELDS if k in rec})
    return rec


def run_block(
    side: Side, mode: str, row: str, arg: str, env: dict[str, str], log_dir: Path
) -> dict:
    """Run one block in a fresh process and return its record (``one_block``); an
    operating-system error during the block aborts the run (Abort) and an interrupt is raised
    again, each naming the block (side, mode and row)."""
    block = f"{side.label} {mode} {row}"
    try:
        return one_block(side, mode, row, arg, env, log_dir)
    except OSError as exc:
        raise Abort(f"{block}: operating-system error: {exc}") from exc
    except KeyboardInterrupt as exc:
        raise KeyboardInterrupt(block) from exc
    except Abort:
        raise
    except Exception as exc:
        raise Abort(f"{block}: {type(exc).__name__}: {exc}") from exc


def one_block(
    side: Side, mode: str, row: str, arg: str, env: dict[str, str], log_dir: Path
) -> dict:
    """One block in a fresh process and its record.

    The side's binary is verified before and after the block. Peak RSS of the harness process
    tree, two readings:
    * ``rss_self_bytes + rss_children_bytes``: reported by the harness at exit from
      ``getrusage(RUSAGE_SELF)`` and ``getrusage(RUSAGE_CHILDREN)`` (the sum is an upper
      bound on the simultaneous peak of the tree; exact without children);
    * ``rss_wait4_bytes``: ``ru_maxrss`` from ``os.wait4`` on the harness process, which
      covers the process and its waited-for descendants (largest member).
    ``rss_tree_bytes`` is the larger of the two and is what the memory rule reads.
    ``block_s`` is the harness's own wall-clock measurement of its timed loop.

    The harness's stderr is appended to ``log_dir/<row>.<label>.<mode>.stderr``, one log for
    every block of the row, side and mode. A log that is still empty after the block is
    removed, also when the harness could not start or was interrupted, and an abort names the
    log only when this block wrote to it.
    """
    side.verify(f"before {mode} {row}")
    log_dir.mkdir(parents=True, exist_ok=True)
    load_before = system_number(
        f"{side.label} {mode} {row} load before", os.getloadavg()[0]
    )
    t0 = time.perf_counter()
    err_path = log_dir / f"{row}.{side.label}.{mode}.stderr"
    try:
        with open(err_path, "ab") as err:
            logged = os.fstat(err.fileno()).st_size  # the earlier blocks' text
            try:
                proc = subprocess.Popen(
                    side.cmd(mode, row, arg),
                    stdout=subprocess.PIPE,
                    stderr=err,
                    env=env,
                    cwd=REPO,
                )
            except OSError as exc:
                raise Abort(
                    f"{side.label} {mode} {row}: cannot start {side.binary} ({exc})"
                ) from exc
            try:
                out = proc.stdout.read() if proc.stdout else b""
                _, status, ru = os.wait4(proc.pid, 0)
            except BaseException:
                # an interrupt (Ctrl-C): the block's harness does not run on
                proc.kill()
                proc.wait()
                raise
            finally:
                if proc.stdout:
                    proc.stdout.close()
            returncode = os.waitstatus_to_exitcode(status)
            proc.returncode = returncode  # reaped by os.wait4, so Popen waits no more
    except BaseException:
        # the run aborts for the first error; a log that cannot be checked is left as it is
        with contextlib.suppress(OSError):
            if err_path.stat().st_size == 0:
                err_path.unlink()
        raise
    wall = time.perf_counter() - t0
    size = err_path.stat().st_size
    if size == 0:
        err_path.unlink()
    stderr = f"stderr: {err_path}" if size > logged else "nothing on stderr"
    if returncode != 0:
        raise Abort(f"{side.label} {mode} {row} exited {returncode} ({stderr})")
    side.verify(f"after {mode} {row}")
    rec = parse_block(out, f"{side.label} {mode} {row} ({stderr})", row)
    rec["side"] = side.label
    rec["wall_s"] = system_number(f"{side.label} {mode} {row} wall time", wall)
    rec["load1_before"] = load_before
    rec["load1_after"] = system_number(
        f"{side.label} {mode} {row} load after", os.getloadavg()[0]
    )
    usage = system_number(f"{side.label} {mode} {row} wait4 peak", ru.ru_maxrss)
    rec["rss_wait4_bytes"] = system_number(
        f"{side.label} {mode} {row} wait4 bytes", peak_rss(usage)
    )
    rec["rss_tree_bytes"] = max(
        rec["rss_self_bytes"] + rec["rss_children_bytes"], rec["rss_wait4_bytes"]
    )
    require_finite(f"{side.label} {mode} {row} rss_tree_bytes", rec["rss_tree_bytes"])
    return rec


def list_rows(side: Side, env: dict[str, str]) -> dict[str, dict]:
    """The side's row list (``row_list`` of ``aba-list``); Refused when the harness does not exit
    0 or its response is not one."""
    cmd = side.cmd("aba-list")
    what = " ".join(cmd)
    try:
        proc = subprocess.run(cmd, capture_output=True, env=env, cwd=REPO)
    except OSError as exc:
        raise Refused(f"{side.label}: {what}: {exc}") from exc
    except KeyboardInterrupt as exc:
        raise KeyboardInterrupt(f"{side.label} aba-list") from exc
    if proc.returncode != 0:
        tail = proc.stderr.decode("utf-8", "replace").strip()[-300:]
        raise Refused(f"{side.label}: {what}: exit {proc.returncode}: {tail}")
    try:
        return row_list(proc.stdout, what)
    except ResponseError as exc:
        raise Refused(f"{side.label}: {exc}") from exc


def median(xs: list[float]) -> float:
    return statistics.median(xs)


def bootstrap_interval(
    x: list[float], resamples: int, seed: int, level: float = 0.95
) -> tuple[float, float]:
    """Percentile interval of the median over whole-triple resamples."""
    rng = random.Random(seed)
    n = len(x)
    stats = sorted(
        median([x[rng.randrange(n)] for _ in range(n)]) for _ in range(resamples)
    )
    lo_q, hi_q = (1 - level) / 2, 1 - (1 - level) / 2
    lo = stats[int(math.floor(lo_q * (resamples - 1)))]
    hi = stats[int(math.ceil(hi_q * (resamples - 1)))]
    return lo, hi


def side_clusters(times: list[float], tolerance: float) -> dict | None:
    """Same-side clusters: one side's timed per-iteration times (positive, ``BLOCK_FIELDS``)
    sorted; a split between consecutive blocks whose relative gap exceeds the tolerance, with at
    least two blocks on each side of it, makes the side bimodal. Returns the widest such split, or
    None."""
    s = sorted(times)
    best = None
    for i in range(2, len(s) - 1):
        gap = s[i] / s[i - 1] - 1
        if gap > tolerance and (best is None or gap > best["gap"]):
            best = {
                "gap": gap,
                "low": {"blocks": i, "min_s": s[0], "max_s": s[i - 1]},
                "high": {"blocks": len(s) - i, "min_s": s[i], "max_s": s[-1]},
            }
    return best


def machine_state() -> dict:
    """Load averages and the top CPU consumers (a machine-state line at the start and the
    end of every report), diagnostics no rule reads: the load is None when the system gives
    none, and a line of ``ps``'s output that is not a finite nonnegative CPU share and a command
    name is left out."""
    cmd = (
        ["ps", "-Ao", "pcpu=,comm=", "-r"]
        if sys.platform == "darwin"
        else ["ps", "-eo", "pcpu=,comm=", "--sort=-pcpu"]
    )
    try:
        out = subprocess.run(
            cmd, capture_output=True, text=True, errors="replace", timeout=10
        ).stdout
    except (OSError, subprocess.SubprocessError):
        out = ""
    top = []
    for line in out.splitlines()[:5]:
        cpu, _, comm = line.strip().partition(" ")
        try:
            share = float(cpu)
        except ValueError:
            continue
        if math.isfinite(share) and share >= 0 and comm.strip():
            top.append([share, Path(comm.strip()).name])
    try:
        load = [system_number("machine load", x) for x in os.getloadavg()]
    except OSError:
        load = None
    if load is not None and len(load) != 3:
        raise Abort("machine load: expected three readings")
    return {
        "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "load": load,
        "top_cpu": top,
    }


def machine_line(state: dict) -> str:
    load = "/".join(f"{x:.2f}" for x in state["load"]) if state["load"] else "n/a"
    top = ", ".join(f"{name} {cpu:.0f}%" for cpu, name in state["top_cpu"]) or "n/a"
    return f"{state['utc']}: load {load}; top CPU {top}"


def cpu_name() -> str:
    """The CPU's name, a diagnostic: ``sysctl`` on macOS, else ``platform.processor()``."""
    if sys.platform != "darwin":
        return platform.processor()
    try:
        cmd = ["sysctl", "-n", "machdep.cpu.brand_string"]
        proc = subprocess.run(cmd, capture_output=True, text=True, errors="replace")
        return proc.stdout.strip()
    except OSError:
        return ""


def cpu_count():
    value = os.cpu_count()
    return (
        None
        if value is None
        else system_number("CPU count", value, types=(int,), domain=POSITIVE)
    )


def rss_balanced(triples: list[dict], budget: float, floor_bytes: float) -> dict:
    """The memory rule: the maximum of the n A1 block peaks against the maximum of the n B block
    peaks, equal sample counts on both sides. A1 is fixed in advance for every row: a larger A2
    peak or a larger peak elsewhere never raises the reference. Only the timed triples count
    (warm-up blocks are not in ``triples``). The allowance is in bytes, so the 8 MiB floor
    applies even when the ratio exceeds 1.10."""
    peaks = {k: [t[k]["rss_tree_bytes"] for t in triples] for k in ("A1", "B", "A2")}
    reference, candidate = max(peaks["A1"]), max(peaks["B"])
    allowance = max(budget * reference, reference + floor_bytes)
    return {
        "rss_a1_max_bytes": reference,
        "rss_b_max_bytes": candidate,
        "rss_balanced_ratio": candidate / reference,
        "rss_allowance_bytes": allowance,
        "rss_a2_max_bytes": max(peaks["A2"]),
        "rss_verdict": "pass" if candidate <= allowance else "fail",
    }


def row_verdict(
    time_verdict: str,
    rss_verdict: str,
    block_verdict: str,
    cluster_verdict: str,
    digest_verdict: str,
) -> str:
    """``fail`` if any check fails, else ``inconclusive`` if timing is, else ``pass``."""
    verdicts = [
        time_verdict,
        rss_verdict,
        block_verdict,
        cluster_verdict,
        digest_verdict,
    ]
    if "fail" in verdicts:
        return "fail"
    return "inconclusive" if time_verdict == "inconclusive" else "pass"


def row_statistics(
    triples: list[dict],
    tolerance: float,
    budget: float,
    floor_bytes: float,
    min_block_s: float,
    resamples: int,
    seed: int,
) -> dict:
    """The statistics and verdicts of one measured row; the row fails unless every timed block
    of both sides reports one and the same digest. A derived value that is not finite aborts the
    run instead (``require_finite``)."""
    a1 = [t["A1"]["median_s"] for t in triples]
    b = [t["B"]["median_s"] for t in triples]
    a2 = [t["A2"]["median_s"] for t in triples]
    means = [(x + y) / 2 for x, y in zip(a1, a2, strict=True)]
    require_finite("the mean of A1 and A2", means)
    ratios = [bt / m for bt, m in zip(b, means, strict=True)]
    x = [r - 1 for r in ratios]
    stat = median(x)
    lo, hi = bootstrap_interval(x, resamples, seed)
    drift = median([abs(p / q - 1) for p, q in zip(a1, a2, strict=True)])
    if drift > tolerance:
        time_verdict = "inconclusive"
    elif hi <= tolerance:
        time_verdict = "pass"
    else:
        time_verdict = "fail"
    rss = rss_balanced(triples, budget, floor_bytes)
    # Block length: the harness's measured wall time of every timed block.
    blocks = [(i, k, t[k]) for i, t in enumerate(triples, 1) for k in ("A1", "B", "A2")]
    short = [
        {"triple": i, "slot": k, "block_s": blk["block_s"]}
        for i, k, blk in blocks
        if blk["block_s"] < min_block_s
    ]
    block_verdict = "fail" if short else "pass"
    # Same-side clusters: every timed block's per-iteration time per side.
    blocks_a = [t[k]["median_s"] for t in triples for k in ("A1", "A2")]
    blocks_b = [t["B"]["median_s"] for t in triples]
    clusters = {
        side: side_clusters(times, tolerance)
        for side, times in (("A", blocks_a), ("B", blocks_b))
    }
    bimodal = [side for side, c in clusters.items() if c is not None]
    cluster_verdict = "fail" if bimodal else "pass"
    digests_a = {t[k]["digest"] for t in triples for k in ("A1", "A2")}
    digests_b = {t["B"]["digest"] for t in triples}
    digest_match = len(digests_a) == 1 and digests_a == digests_b
    digest_verdict = "pass" if digest_match else "fail"
    durations = [blk["block_s"] for _, _, blk in blocks]
    stats = {
        "ratios": ratios,
        "statistic": stat,
        "ci95": [lo, hi],
        "drift": drift,
        "time_verdict": time_verdict,
        **rss,
        "block_s_min": min(durations),
        "block_s_median": median(durations),
        "short_blocks": short,
        "block_verdict": block_verdict,
        "blocks_a_s": blocks_a,
        "blocks_b_s": blocks_b,
        "clusters": clusters,
        "bimodal": bimodal,
        "cluster_verdict": cluster_verdict,
        "median_a_s": median(a1 + a2),
        "median_b_s": median(b),
        "load1_max": max(
            max(blk["load1_before"], blk["load1_after"]) for _, _, blk in blocks
        ),
        "digest_a": sorted(digests_a),
        "digest_b": sorted(digests_b),
        "digest_match": digest_match,
        "digest_verdict": digest_verdict,
        "verdict": row_verdict(
            time_verdict,
            rss["rss_verdict"],
            block_verdict,
            cluster_verdict,
            digest_verdict,
        ),
    }
    require_finite("row statistic", stats)
    return stats


def glob_matches(name: str, meta: dict, pattern: str) -> bool:
    """A pattern matches a row's name or its group, case-sensitively; ``all`` is ``*``."""
    pattern = "*" if pattern == "all" else pattern
    return fnmatch.fnmatchcase(name, pattern) or fnmatch.fnmatchcase(
        meta["group"], pattern
    )


def split_patterns(text: str) -> list[str]:
    """``--rows`` as globs: comma-separated, spaces around each glob ignored, empty items
    dropped; no text at all means ``all``."""
    if not text.strip():
        return ["all"]
    return [p.strip() for p in text.split(",") if p.strip()]


def groups_of(rows: dict[str, dict]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for name, rec in rows.items():
        out.setdefault(rec["group"], []).append(name)
    return out


def choose_rows(rows: dict[str, dict], rows_arg: str) -> dict:
    """Resolve ``--rows`` against the harnesses' row list.

    Errors (the run is refused before any block): a text that names no glob; a glob that matches
    no row or group on its own. The record lists, per glob, the rows it selected and the rows of
    their groups it did not, and per group the rows the whole selection did not select.
    """
    pats = split_patterns(rows_arg)
    errors = [] if pats else [f"--rows {rows_arg!r} names no glob"]
    groups = groups_of(rows)
    globs = []
    for pat in pats:
        hit = [n for n, rec in rows.items() if glob_matches(n, rec, pat)]
        touched = sorted({rows[n]["group"] for n in hit})
        globs.append(
            {
                "glob": pat,
                "matched": hit,
                "groups": touched,
                "group_rows_not_matched": [
                    n for g in touched for n in groups[g] if n not in hit
                ],
            }
        )
        if not hit:
            errors.append(f"--rows glob {pat!r} matches no row or group")
    names = [n for n in rows if any(n in g["matched"] for g in globs)]
    return {
        "rows": rows_arg,
        "selected": names,
        "globs": globs,
        "not_selected": {
            g: [n for n in members if n not in names]
            for g, members in groups.items()
            if any(n not in names for n in members)
        },
        "errors": errors,
    }


def same_rows(rows_old: dict[str, dict], rows_new: dict[str, dict]) -> None:
    """Refused unless both harnesses list the same rows, record for record, in the same order
    (``build.sh`` compiles this checkout's rows into both)."""
    if list(rows_old.items()) == list(rows_new.items()):
        return
    only_old = [n for n in rows_old if n not in rows_new]
    only_new = [n for n in rows_new if n not in rows_old]
    other = [n for n in rows_new if n in rows_old and rows_old[n] != rows_new[n]]
    raise Refused(
        "the two harnesses list different rows (build both with this checkout's"
        f" build.sh): only the reference has {only_old or 'none'}, only the candidate"
        f" {only_new or 'none'}, other records {other or 'none'}"
        + ("" if only_old or only_new or other else ", another order")
    )


def row_choice_error_text(sel: dict, rows: dict[str, dict]) -> str:
    listing = "\n".join(
        f"  {group} ({len(names)}): {', '.join(names)}"
        for group, names in groups_of(rows).items()
    )
    return (
        "; ".join(sel["errors"])
        + '.\nThe groups that exist (the harness\'s list; README "Selecting rows"; a'
        " pattern matches a row's name or its group):\n" + listing
    )


# The largest count the harness reads (a 64-bit signed integer).
MAX_COUNT = 2**63 - 1


def parse_count(text: str) -> tuple[int, int]:
    """``N`` (both sides) or ``NA:NB`` (per side), each a positive integer the harness reads (at
    most ``MAX_COUNT``); ValueError otherwise (``NA:`` too: a count is never implied)."""
    a, sep, b = text.partition(":")
    try:
        pair = (int(a), int(b if sep else a))
    except ValueError:
        raise ValueError(f"iteration count {text!r} is not N or NA:NB") from None
    if min(pair) < 1:
        raise ValueError(f"iteration count {text!r} below 1")
    if max(pair) > MAX_COUNT:
        raise ValueError(f"iteration count {text!r} above {MAX_COUNT}")
    return pair


def pinned_iterations(text: str | None, names: list[str]) -> dict[str, tuple[int, int]]:
    """Pinned per-side counts ``row -> (n_A, n_B)`` from ``--iters row=N`` (both sides) or
    ``row=NA:NB``, every item naming a selected row (``names``), none twice; ValueError for any
    item that is not valid."""
    pins: dict[str, tuple[int, int]] = {}
    for item in (text or "").split(","):
        if not item.strip():
            continue  # empty items are dropped, as in --rows
        k, sep, v = item.partition("=")
        k = k.strip()
        if not sep or not k:
            raise ValueError(f"--iters item {item!r} is not row=N or row=NA:NB")
        if k not in names:
            raise ValueError(f"--iters row {k!r} is not a selected row")
        if k in pins:
            raise ValueError(f"--iters row {k!r} given twice")
        pins[k] = parse_count(v.strip())
    return pins


def side_counts(
    warm_a: dict,
    warm_b: dict,
    target_s: float,
    margin: float,
    pin: tuple[int, int] | None,
) -> tuple[int, int, str]:
    """The iteration counts of the timed blocks: each side ``ceil(target * (1 + margin) /
    median_side)`` from its own warm-up block, so every block on either side is calibrated to
    the target; the block value is the median iteration time, so the ratio does not depend on
    the counts. A pin (``--iters``) overrides both. A count that is not finite aborts the run,
    naming the row and the side (the warm-up record's ``row`` and ``side``)."""
    if pin is not None:
        return pin[0], pin[1], "pinned"
    want = target_s * (1 + margin)
    quotients = [want / warm_a["median_s"], want / warm_b["median_s"]]
    for warm, quotient in zip((warm_a, warm_b), quotients, strict=True):
        require_finite(
            f"{warm['row']} {warm['side']}: the calibrated iteration count", quotient
        )
        if quotient > MAX_COUNT:
            raise Abort(
                f"{warm['side']} aba-calibrate {warm['row']}: calibrated iteration count above {MAX_COUNT}"
            )
    return (
        max(1, math.ceil(quotients[0])),
        max(1, math.ceil(quotients[1])),
        "calibrated per side",
    )


def source_hashes() -> dict[str, str]:
    """sha256 of every file of ``SOURCE_FILES``; Refused when one cannot be read or the
    directory's path is not printable text (``meta.txt`` records it on one line)."""
    if not str(HERE).isprintable():
        raise Refused(f"{str(HERE)!r}: a path with a character that is not printable")
    try:
        return {name: sha256_file(HERE / name) for name in SOURCE_FILES}
    except OSError as exc:
        raise Refused(f"the tool's sources: {exc}") from exc


META_HEADER = """\
# aba.py provenance: sha256 of both harness binaries (each the driver compiled with that side's
# Core; re-hashed before and after every block) and of the orchestrator and driver sources
# (hashed at the start of the run). Columns: side, role, name, sha256, path (tab separated).
"""


def meta_lines(sides: list[Side], sources: dict[str, str]) -> list[str]:
    """``meta.txt``'s hash lines: one per binary, one per source."""
    lines = [f"{s.key}\tbinary\tharness\t{s.sha256}\t{s.binary}" for s in sides]
    lines += [
        f"run\tsource\t{name}\t{sha}\t{HERE / name}" for name, sha in sources.items()
    ]
    return lines


def meta_table(lines: list[str]) -> dict[str, str]:
    """``side/role/name -> sha256`` of hash lines."""
    return {"/".join(cols[:3]): cols[3] for cols in (ln.split("\t") for ln in lines)}


def read_meta(path: Path) -> dict[str, str]:
    """A ``meta.txt`` as ``write_meta`` writes it, Refused otherwise: every line a comment, a
    ``label`` line of side A or B, or a hash line (side A, B or ``run``, role, name, sha256,
    path), no hash named twice and at least one hash."""
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise Refused(f"--expect-meta {path}: {exc}") from exc
    table = {}
    labels = set()
    if text and not text.endswith("\n"):
        raise Refused(f"--expect-meta {path}: last line is not ended by a newline")
    for i, line in enumerate(text.splitlines(), 1):
        cols = line.split("\t")
        if line.startswith("#"):
            continue
        if (
            len(cols) == 3
            and cols[0] == "label"
            and cols[1] in ("A", "B")
            and LABEL.fullmatch(cols[2])
            and cols[1] not in labels
        ):
            labels.add(cols[1])
            continue
        hash_line = (
            len(cols) == 5
            and is_sha256(cols[3])
            and cols[4].isprintable()
            and (
                (cols[0] in ("A", "B") and cols[1:3] == ["binary", "harness"])
                or (cols[0:2] == ["run", "source"] and cols[2] in SOURCE_FILES)
            )
        )
        if not hash_line:
            raise Refused(
                f"--expect-meta {path}: line {i} is not a meta.txt line: {line!r}"
            )
        key = "/".join(cols[:3])
        if key in table:
            raise Refused(f"--expect-meta {path}: {key} named twice (line {i})")
        table[key] = cols[3]
    if not table:
        raise Refused(f"--expect-meta {path}: no hash")
    return table


def write_meta(path: Path, sides: list[Side], lines: list[str]) -> None:
    labels = [f"label\t{s.key}\t{s.label}" for s in sides]
    path.write_text(META_HEADER + "\n".join(labels + lines) + "\n", encoding="utf-8")


CSV_FIELDS = (
    "row",
    "slot",
    "triple",
    "side",
    "iterations",
    "iterations_requested",
    "block_s",
    "min_block_s",
    "meets_min",
    "median_s",
    "mean_s",
    "min_s",
    "max_s",
    "setup_s",
    "process_wall_s",
    "rss_self_bytes",
    "rss_children_bytes",
    "rss_wait4_bytes",
    "rss_tree_bytes",
    "digest",
    "load1_before",
    "load1_after",
)


def csv_rows(report: dict) -> list[dict]:
    """Every block run so far, one line each; ``block_s`` is the measured wall time of the
    block's timed loop (warm-up blocks: the calibration loop, not checked). A row not run, or one
    whose warm-up did not finish, has no warm-up blocks; a row not timed yet has no triples."""
    min_block = report["protocol"]["target_block_s"]
    out = []
    for rec in report["rows"]:
        for k, blk in rec.get("warmup", {}).items():
            out.append(
                {
                    **blk,
                    "row": rec["row"],
                    "slot": f"W-{k}",
                    "process_wall_s": blk["wall_s"],
                }
            )
        for i, t in enumerate(rec.get("triples", []), 1):
            for k in ("A1", "B", "A2"):
                blk = t[k]
                out.append(
                    {
                        **blk,
                        "row": rec["row"],
                        "slot": k,
                        "triple": i,
                        "min_block_s": min_block,
                        "meets_min": blk["block_s"] >= min_block,
                        "iterations_requested": rec[
                            "iterations_b" if k == "B" else "iterations_a"
                        ],
                        "process_wall_s": blk["wall_s"],
                    }
                )
    return out


def csv_text(report: dict) -> str:
    """blocks.csv's text (``csv_rows``)."""
    f = io.StringIO()
    w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
    w.writeheader()
    w.writerows(csv_rows(report))
    return f.getvalue()


def write_csv(report: dict, path: Path) -> None:
    path.write_text(csv_text(report), encoding="utf-8", newline="")


def report_text(report: dict) -> str:
    """report.json's text: standard JSON, so a value that is not finite (NaN or an infinity,
    which Python's writer would print as nonstandard literals) raises ValueError instead of
    reaching the file."""
    return json.dumps(report, indent=1, allow_nan=False)


def write_report(report: dict, path: Path) -> None:
    path.write_text(report_text(report), encoding="utf-8")


def run_verdict(report: dict) -> str:
    """``aborted``, else the worst row verdict. A run that did not abort measured every selected
    row, and a selection has at least one (``choose_rows``)."""
    if "aborted" in report:
        return "aborted"
    verdicts = [r["verdict"] for r in report["rows"] if r["status"] == "measured"]
    if "fail" in verdicts:
        return "fail"
    return "inconclusive" if "inconclusive" in verdicts else "pass"


def conforming(report: dict) -> bool:
    """Gate evidence: not aborted, no ``--allow-short`` (without it ``main`` refuses fewer
    triples or a smaller target than the protocol), no short block and the normal allocator."""
    return (
        "aborted" not in report
        and not report["protocol"]["allow_short"]
        and not any(r.get("short_blocks") for r in report["rows"])
        and report["allocator_env"]["normal"]
    )


def start_run(args, old: Side, new: Side):
    """Everything before the first block that writes nothing: ``--out``, each build record's
    driver (``same_driver``), ``meta.txt``'s lines, ``--expect-meta``, both row lists, the
    selection, ``--iters`` (the pins) and the report's start. Refused for any input outside its
    rule; an operating-system error escapes (``run``)."""
    out = (
        Path(args.out)
        if args.out
        else HERE / "results" / dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    )
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        # A file of an earlier run would mix into this one's (the stderr logs append).
        raise Refused(f"--out {out}: must be absent or an empty directory")
    sources = source_hashes()
    for side in (old, new):
        same_driver(side, sources)
    lines = meta_lines([old, new], sources)
    if args.expect_meta:
        expected = read_meta(Path(args.expect_meta))
        have = meta_table(lines)
        diff = sorted(
            k for k in expected.keys() | have.keys() if expected.get(k) != have.get(k)
        )
        if diff:
            raise Refused(
                f"--expect-meta {args.expect_meta}: hashes differ: {', '.join(diff)}"
            )
    env = dict(os.environ)  # the harnesses run with the caller's environment, unchanged
    rows = list_rows(new, env)
    same_rows(list_rows(old, env), rows)
    chosen = choose_rows(rows, args.rows)
    if chosen["errors"]:
        raise Refused(row_choice_error_text(chosen, rows))
    names = chosen["selected"]
    try:
        pins = pinned_iterations(args.iters, names)
    except ValueError as exc:
        raise Refused(f"--iters: {exc}") from exc
    report = {
        "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "argv": list(sys.argv),
        "protocol": {
            "triples": args.triples,
            "protocol_triples": PROTOCOL_TRIPLES,
            "target_block_s": args.target_seconds,
            "protocol_block_s": PROTOCOL_BLOCK_SECONDS,
            "allow_short": args.allow_short,
            "calibration_margin": args.calibration_margin,
            "tolerance": args.tolerance,
            "budget": args.budget,
            "rss_floor_bytes": args.rss_floor_mib * 2**20,
            "bootstrap_resamples": args.resamples,
            "bootstrap_seed": args.seed,
        },
        "selection": chosen,
        "machine": {
            "platform": platform.platform(),
            "cpu": cpu_name(),
            "ncpu": cpu_count(),
            "state_at_start": machine_state(),
        },
        "sides": {
            s.key: {
                "label": s.label,
                "binary": s.binary,
                "sha256": s.sha256,
                "build_record": s.build_record,
            }
            for s in (old, new)
        },
        "meta": meta_table(lines),
        "allocator_env": allocator_record(env),
        # Every selected row from the start, so an abort at any point reports each one.
        "rows": [
            {
                "row": name,
                "group": rows[name]["group"],
                "desc": rows[name].get("desc"),
                "status": "not run",
            }
            for name in names
        ],
    }
    if args.expect_meta:
        report["expect_meta"] = args.expect_meta
    return out, report, lines, env, pins


def create_outputs(
    out: Path, sides: list[Side], lines: list[str], report: dict | None = None
) -> None:
    """``--out`` and ``meta.txt``, the last step before the first block. On an operating-system
    error (Refused) or an interrupt, what this step created is removed again, so a run refused or
    interrupted before the first block leaves no partial setup; main writes its refusal reason."""
    created: list[Path] = []
    try:
        created = [d for d in (out, *out.parents) if not d.exists()]  # deepest first
        out.mkdir(parents=True, exist_ok=True)
        write_meta(out / "meta.txt", sides, lines)
        if report is not None:
            pending = {**report, "aborted": "run has not finished"}
            atomic_text(out / "report.json", report_text(minimal_report(pending)))
            atomic_text(out / "summary.md", minimal_summary(pending))
    except (Exception, KeyboardInterrupt) as exc:
        for name in ("meta.txt", "report.json", "summary.md"):
            with contextlib.suppress(OSError):
                (out / name).unlink(missing_ok=True)
        for d in created:
            with contextlib.suppress(OSError):
                d.rmdir()
        if isinstance(exc, OSError):
            raise Refused(f"--out {out}: {exc}") from exc
        raise


def run(args, old: Side, new: Side) -> int:
    """One run; the exit status as in the module docstring (0 completed, 1 aborted, 2 refused;
    an interrupt before the first block reaches ``main``). An operating-system error refuses the
    run before the first block and aborts it from there on."""
    try:
        out, report, lines, env, pins = start_run(args, old, new)
    except OSError as exc:
        raise Refused(f"before the first block: {exc}") from exc
    report_path = out / "report.json"
    t_run = time.perf_counter()
    # One SIGINT handler serves the rest of the run: it raises KeyboardInterrupt at the first
    # interrupt only and records every later one, so a second Ctrl-C cannot cut short the
    # handling of the first. An interrupt while the outputs are created leaves nothing and
    # reaches main; from then on it aborts the run. The final outputs are written with no
    # interrupt raised; an interrupt during them, recorded or raised before they began, writes
    # them again as aborted.
    with interrupts_once() as interrupts:
        create_outputs(out, [old, new], lines, report)
        try:
            try:
                # First the warm-up blocks of every selected row (counts and the plan).
                for rec in report["rows"]:
                    prepare_row(args, rec, old, new, env, out, pins)
                    write_csv(report, out / "blocks.csv")
                report["plan"] = run_plan(report)
                write_report(report, report_path)
                print(plan_text(report), flush=True)
                # Then the timed triples, row by row.
                for rec in report["rows"]:
                    run_triples(args, report, report_path, rec, old, new, env)
                    write_csv(report, out / "blocks.csv")
                # The reference's last block (A2) was verified after it ran; the candidate's
                # binary could still change during that block.
                new.verify("at the end of the run")
            # A refused response, a failed harness, a changed binary, a value that is not
            # finite; a block's own operating-system errors among them (run_block).
            except Abort as exc:
                report["aborted"] = str(exc)
            # An output written during the run that cannot be written, a closed stdout.
            except OSError as exc:
                report["aborted"] = f"operating-system error: {exc}{during(report)}"
            except KeyboardInterrupt as exc:
                report["aborted"] = f"interrupted (SIGINT, Ctrl-C){during(report, exc)}"
            # A defect of aba.py: aborted like the others, with its traceback on stderr.
            except Exception as exc:
                with contextlib.suppress(OSError, ValueError):
                    traceback.print_exc()
                report["aborted"] = (
                    f"internal error: {type(exc).__name__}: {exc}{during(report)}"
                )
            before = len(interrupts)
            interrupts.hold()
            status = finish(report, out, t_run)
            if len(interrupts) > before:
                status = finish_interrupted(report, out, t_run)
        except KeyboardInterrupt:
            status = finish_interrupted(report, out, t_run)
    return status


def during(report: dict, exc: BaseException | None = None) -> str:
    """Where an abort reason that does not name its block came from: the block an interrupt
    stopped (``run_block`` raises it again with the block's name), else the row step in
    progress, else nothing (the plan, the end of the run)."""
    if exc is not None and exc.args:
        return f" during {exc}"
    for r in report["rows"]:
        if r["status"] in ("warm-up", "timing"):
            step = "warm-up" if r["status"] == "warm-up" else "timed"
            return f" during the {step} step of {r['row']}"
    return ""


class Interrupts(list):
    """The SIGINTs (Ctrl-C) of a run, in order, as ``interrupts_once`` records them: the first
    one raises KeyboardInterrupt unless ``hold`` came first, and every other one is only
    recorded."""

    def __init__(self):
        super().__init__()
        self._armed = {"first": True}

    def hold(self) -> None:
        """From now on every interrupt is only recorded."""
        self._armed.clear()

    def receive(self, signum: int, _frame) -> None:
        # dict.pop is one call, so two interrupts in a row cannot both take "first": at most one
        # KeyboardInterrupt is ever raised.
        first = self._armed.pop("first", False)
        self.append(signum)
        if first:
            raise KeyboardInterrupt


@contextlib.contextmanager
def interrupts_once():
    """SIGINT handled by the yielded ``Interrupts``, so that a second interrupt cannot cut short
    the handling of the first or a write of the final outputs; the previous handler is restored
    afterwards."""
    interrupts = Interrupts()
    previous = signal.signal(signal.SIGINT, interrupts.receive)
    try:
        yield interrupts
    finally:
        signal.signal(signal.SIGINT, previous)


def warn(text: str) -> None:
    """A line on stderr; a closed stderr is no error."""
    with contextlib.suppress(OSError, ValueError):
        print(text, file=sys.stderr, flush=True)


def finish(report: dict, out: Path, t_run: float) -> int:
    """The verdict, the machine state at the end and the outputs (``write_outputs``); 1 if the
    run aborted, else 0. An output that cannot be made or written aborts the run: the reason
    goes to stderr and, unless the run had aborted already, into every output, each written
    again where it can be (in its minimal form where its full text cannot be made)."""
    if "aborted" in report:
        warn(f"aborted: {report['aborted']}")
    report["finished_utc"] = dt.datetime.now(dt.timezone.utc).isoformat(
        timespec="seconds"
    )
    try:
        report["wall_s"] = system_number("run wall time", time.perf_counter() - t_run)
        report["machine"]["state_at_end"] = machine_state()
    except Exception as exc:
        why = f"while finishing the run: {type(exc).__name__}: {exc}"
        report.setdefault("aborted", why)
        warn(f"aborted: {why}")
    report["verdict"] = run_verdict(report)
    report["conforming"] = conforming(report)
    problems = write_outputs(report, out)
    if problems:
        why = "; ".join(problems)
        if any(p.startswith("operating-system error") for p in problems):
            why += f"; the outputs in {out} are incomplete"
        warn(f"aborted: {why}")
        report.setdefault("aborted", why)
        report["verdict"], report["conforming"] = "aborted", False
        write_outputs(report, out)
        return 1
    # The outputs are written: a stdout that is closed or cannot show a path is no abort.
    with contextlib.suppress(OSError, ValueError):
        print(
            f"report: {out / 'report.json'}\nsummary: {out / 'summary.md'}\n"
            f"blocks: {out / 'blocks.csv'}\nmeta: {out / 'meta.txt'}",
            flush=True,
        )
    return 1 if "aborted" in report else 0


def finish_interrupted(report: dict, out: Path, t_run: float) -> int:
    """An interrupt while the final outputs were written: the run is aborted (an earlier
    reason is kept) and the outputs are written again, each one that can be; 1."""
    why = "interrupted (SIGINT, Ctrl-C) while the outputs were written"
    report.setdefault("aborted", why)
    warn(f"aborted: {why}")
    report.setdefault(
        "finished_utc",
        dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    )
    if "wall_s" not in report:
        try:
            report["wall_s"] = system_number(
                "run wall time", time.perf_counter() - t_run
            )
        except Abort as exc:
            warn(str(exc))
    report["verdict"], report["conforming"] = "aborted", False
    if problems := write_outputs(report, out):
        warn(f"the outputs in {out} are incomplete: {'; '.join(problems)}")
    return 1


def write_outputs(report: dict, out: Path) -> list[str]:
    """Write the final outputs into ``out`` (report.json, blocks.csv, summary.md), each one that
    can be: one whose text cannot be made (a defect of aba.py, as every value of the report
    passed a rule) is written in its minimal form once the run is aborted, and one that cannot
    be written does not stop the next. Why each output not written in full failed."""
    problems = []
    for name, full, minimal in (
        ("report.json", report_text, lambda r: report_text(minimal_report(r))),
        ("blocks.csv", csv_text, lambda r: ",".join(CSV_FIELDS) + "\r\n"),
        ("summary.md", summary_text, minimal_summary),
    ):
        try:
            text = full(report)
        except Exception as exc:  # a defect: the run's abort, not a traceback
            problems.append(f"{name} could not be made ({type(exc).__name__}: {exc})")
            if "aborted" not in report:
                continue  # finish aborts the run and writes the outputs again
            text = minimal(report)
        try:
            atomic_text(out / name, text)
        except Exception as exc:
            kind = (
                "operating-system error"
                if isinstance(exc, OSError)
                else type(exc).__name__
            )
            problems.append(f"{kind} writing {name}: {exc}")
            # A directory may have become read-only while its existing output files remain
            # writable. Keep the abort reason there too when replacement is unavailable.
            if isinstance(exc, OSError) and "aborted" in report:
                try:
                    (out / name).write_text(text, encoding="utf-8", newline="")
                except Exception as retry:
                    problems.append(f"{name} retry failed: {retry}")
    return problems


def atomic_text(path: Path, text: str) -> None:
    """Replace an output only after its complete UTF-8 text is written; remove a failed draft."""
    draft = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=".aba-",
            delete=False,
        ) as stream:
            draft = Path(stream.name)
            stream.write(text)
        os.replace(draft, path)
    finally:
        if draft is not None:
            draft.unlink(missing_ok=True)


def minimal_report(report: dict) -> dict:
    """report.json in its minimal form, for an aborted run whose full report cannot be made:
    the start, the verdict and its reason, and every selected row with the step it reached
    (strings and booleans only, so it can always be written)."""
    return {
        "started_utc": report["started_utc"],
        "verdict": "aborted",
        "conforming": False,
        "aborted": report["aborted"],
        "minimal": "the full report could not be made; why is on stderr",
        "rows": [{"row": r["row"], "status": r["status"]} for r in report["rows"]],
    }


def minimal_summary(report: dict) -> str:
    """summary.md in its minimal form (``minimal_report``), its first lines as the full one's."""
    lines = [
        "# A/B/A report",
        "",
        "* **not gate evidence**: the run aborted (a minimal summary: the full one could not"
        " be made; why is on stderr)",
        f"* started {report['started_utc']}, overall verdict **aborted**",
        f"* **aborted**: {report['aborted']}",
        "",
        "Rows and the step each reached: "
        + ", ".join(f"{r['row']} ({r['status']})" for r in report["rows"])
        + ".",
    ]
    return "\n".join(lines) + "\n"


def prepare_row(args, rec, old, new, env, out, pins) -> None:
    """The warm-up step for one row: warm-up blocks, per-side counts and the expected duration;
    the row's status goes from ``warm-up`` to ``planned``."""
    name = rec["row"]
    rec["status"] = "warm-up"
    target = str(args.target_seconds)
    rec["warmup"] = {
        "A": run_block(old, "aba-calibrate", name, target, env, out / "blocks"),
        "B": run_block(new, "aba-calibrate", name, target, env, out / "blocks"),
    }
    warm_a, warm_b = rec["warmup"]["A"], rec["warmup"]["B"]
    n_a, n_b, source = side_counts(
        warm_a, warm_b, args.target_seconds, args.calibration_margin, pins.get(name)
    )
    rec["iterations_a"], rec["iterations_b"] = n_a, n_b
    rec["iterations_source"] = source
    # Per-block cost outside the timed loop (process start, setup, the untimed first
    # iteration): the warm-up block's process wall time minus its timed loop.
    blk_a = n_a * warm_a["median_s"] + max(0.0, warm_a["wall_s"] - warm_a["block_s"])
    blk_b = n_b * warm_b["median_s"] + max(0.0, warm_b["wall_s"] - warm_b["block_s"])
    row_s = args.triples * (2 * blk_a + blk_b)
    require_finite(f"{name} expected duration", [blk_a, blk_b, row_s])
    rec["expected"] = {"block_a_s": blk_a, "block_b_s": blk_b, "row_s": row_s}
    rec["status"] = "planned"


def run_plan(report: dict) -> dict:
    planned = [r for r in report["rows"] if r["status"] == "planned"]
    total = sum(r["expected"]["row_s"] for r in planned)
    require_finite("the expected duration of the run", total)
    try:
        finish = dt.datetime.now() + dt.timedelta(seconds=total)
    except OverflowError as exc:
        raise Abort(
            f"the expected duration of the run, {total:.3g} s, ends beyond any date ({exc})"
        ) from exc
    return {
        "rows": len(planned),
        "expected_s": total,
        "expected_finish_local": finish.isoformat(timespec="minutes"),
    }


def plan_text(report: dict) -> str:
    lines = [
        "plan (from the warm-up blocks, before the first timed block): "
        f"{report['protocol']['triples']} triples per row",
        f"  {'row':34s} {'iters A / B':>17s} {'block A / B s':>15s} {'expected':>10s}",
    ]
    for r in report["rows"]:
        if r["status"] != "planned":
            lines.append(f"  {r['row']:34s} {r['status']}")
            continue
        e = r["expected"]
        iters = f"{r['iterations_a']} / {r['iterations_b']}"
        blocks = f"{e['block_a_s']:.2f} / {e['block_b_s']:.2f}"
        lines.append(f"  {r['row']:34s} {iters:>17s} {blocks:>15s} {e['row_s']:9.0f}s")
    p = report["plan"]
    lines.append(
        f"  total: {p['rows']} rows, expected {p['expected_s']:.0f} s "
        f"({p['expected_s'] / 60:.1f} min), finishing about {p['expected_finish_local']}"
    )
    return "\n".join(lines)


def run_triples(args, report, report_path, rec, old, new, env) -> None:
    """The timed step for one row: the timed triples with the per-side counts (each block must
    report the count it was given), then the row statistics; the row's status goes from
    ``timing`` to ``measured``."""
    name = rec["row"]
    rec["status"] = "timing"
    t_row = time.perf_counter()
    counts = {
        "A1": rec["iterations_a"],
        "B": rec["iterations_b"],
        "A2": rec["iterations_a"],
    }
    rec["triples"] = []
    for _ in range(args.triples):
        trip = {}
        for slot, side in (("A1", old), ("B", new), ("A2", old)):
            blk = run_block(
                side,
                "aba-run",
                name,
                str(counts[slot]),
                env,
                report_path.parent / "blocks",
            )
            if blk["iterations"] != counts[slot]:
                raise Abort(
                    f"{side.label} aba-run {name}: the harness ran {blk['iterations']} "
                    f"iterations, {counts[slot]} requested"
                )
            trip[slot] = blk
        rec["triples"].append(trip)
        write_report(report, report_path)
    try:
        stats = row_statistics(
            rec["triples"],
            args.tolerance,
            args.budget,
            args.rss_floor_mib * 2**20,
            args.target_seconds,
            args.resamples,
            args.seed,
        )
    except Abort as exc:
        raise Abort(f"{name}: {exc}") from exc  # a statistic of this row's blocks
    rec.update(stats)
    rec["row_wall_s"] = system_number(
        f"{name} row wall time", time.perf_counter() - t_row
    )
    rec["status"] = "measured"
    lo, hi = rec["ci95"]
    print(
        f"{name}: {rec['verdict']:12s} stat {rec['statistic'] * 100:+6.2f}% "
        f"CI [{lo * 100:+6.2f}, {hi * 100:+6.2f}]% drift {rec['drift'] * 100:5.2f}% "
        f"rss max B/A1 x{rec['rss_balanced_ratio']:.3f} {rec['rss_verdict']} "
        f"iters {rec['iterations_a']}/{rec['iterations_b']} "
        f"block min {rec['block_s_min']:.2f}s ({len(rec['short_blocks'])} short) "
        f"digest {'same' if rec['digest_match'] else 'DIFFER'} "
        f"({rec['row_wall_s']:.0f}s, expected {rec['expected']['row_s']:.0f}s)",
        flush=True,
    )
    write_report(report, report_path)


def pct(x: float) -> str:
    return f"{x * 100:+.2f}"


def iters_text(rec: dict) -> str:
    return f"{rec['iterations_a']} / {rec['iterations_b']}"


def row_choice_lines(report: dict) -> list[str]:
    """The rows each glob selected and the rows of the groups it did not."""
    sel = report["selection"]
    status = {r["row"]: r["status"] for r in report["rows"]}
    lines = [f"Selection (`--rows {sel['rows']}`):", ""]
    for g in sel["globs"]:
        rows = ", ".join(f"{n} ({status[n]})" for n in g["matched"])
        lines.append(f"* `{g['glob']}` selected {len(g['matched'])}: {rows}")
        if g["group_rows_not_matched"]:
            lines.append(
                f"  * rows of its groups ({', '.join(g['groups'])}) it did not select: "
                + ", ".join(g["group_rows_not_matched"])
            )
    lines.append(
        "* rows the selection did not select, by group: "
        + (
            "; ".join(
                f"{grp}: {', '.join(names)}"
                for grp, names in sel["not_selected"].items()
            )
            or "none"
        )
    )
    return lines


def plan_lines(report: dict) -> list[str]:
    """The plan with each row's expected and actual duration; none when the run aborted during
    the warm-up blocks."""
    plan = report.get("plan")
    if plan is None:
        return []
    lines = [
        f"Plan printed before the first timed block: {plan['rows']} rows, expected"
        f" {plan['expected_s'] / 60:.1f} min (finishing about"
        f" {plan['expected_finish_local']}); per row expected / actual s:",
        "",
    ]
    for r in report["rows"]:
        if "expected" not in r:
            continue
        actual = r.get("row_wall_s")
        lines.append(
            f"* {r['row']}: iters {iters_text(r)}, expected"
            f" {r['expected']['row_s']:.0f} s (blocks A {r['expected']['block_a_s']:.2f}"
            f" s, B {r['expected']['block_b_s']:.2f} s), actual "
            + (
                f"{actual:.0f} s"
                if actual is not None
                else "not finished"
                if r["status"] == "timing"
                else "not run"
            )
        )
    return lines


def summary_text(report: dict) -> str:
    """summary.md: the run's conformance and verdict first, then the table of rows."""
    p = report["protocol"]
    min_block = p["target_block_s"]
    floor_mib = p["rss_floor_bytes"] / 2**20
    machine = report["machine"]
    lines = [
        "# A/B/A report",
        "",
        "* **conforming run** (the protocol's evidence: accepted when the overall verdict"
        " is pass)"
        if report["conforming"]
        else "* **not gate evidence**: the run does not conform to the protocol (see the"
        " protocol, allocator and aborted lines)",
        f"* started {report['started_utc']}, wall {report['wall_s'] / 60:.1f} min,"
        f" overall verdict **{report['verdict']}**",
        f"* A: `{report['sides']['A']['label']}`; B: `{report['sides']['B']['label']}`",
        f"* protocol: triples {p['triples']} (protocol {p['protocol_triples']}),"
        f" minimum timed block {min_block} s (protocol {p['protocol_block_s']} s)"
        + (
            ", **--allow-short given: below the protocol, not gate evidence**"
            if p["allow_short"]
            else ""
        )
        + f"; tolerance {p['tolerance'] * 100:.1f}%;"
        f" memory: max B block peak <= max(x{p['budget']:.2f}, +{floor_mib:.0f} MiB) of max"
        f" A1 block peak; bootstrap {p['bootstrap_resamples']} resamples",
        f"* machine {machine['cpu']}, {machine['ncpu']} CPUs",
        f"* machine state at start: {machine_line(machine['state_at_start'])}",
        f"* machine state at end: {machine_line(machine['state_at_end'])}",
        f"* allocator environment: {allocator_line(report['allocator_env'])}",
        "* provenance: binary and source hashes in `meta.txt`; each binary re-hashed before and"
        " after every block"
        + (
            f"; every hash equal to `{report['expect_meta']}`"
            if "expect_meta" in report
            else ""
        ),
        *([f"* **aborted**: {report['aborted']}"] if "aborted" in report else []),
        "",
        "Statistic = median over triples of `B / mean(A1, A2) - 1`; CI = 95% bootstrap"
        " over triples; drift = median `|A1/A2 - 1|`. Times are per-iteration medians."
        " Block s = shortest / median measured timed-block wall time (every block is in"
        " `blocks.csv`). RSS max A1 / B = the largest A1 block peak (A1, the first"
        " reference block of every triple, fixed in advance) and the largest B block peak;"
        f" RSS allowance = max({p['budget']:.2f} x max A1, max A1 + {floor_mib:.0f} MiB), the"
        " bound on max B; RSS max B/A1 = their ratio (warm-up blocks never count). Iters A / B"
        " = the per-side iteration counts of the timed blocks (each side calibrated on its own"
        " warm-up, unless pinned). Digest = whether every timed block of both sides reported"
        " one and the same digest.",
        "",
        "| Row | Iters A / B | Block s min / med | A ms | B ms | Stat % | 95% CI % |"
        " Drift % | RSS max A1 / B MiB | RSS allowance MiB | RSS max B/A1 | Digest |"
        " Load max | Verdict |",
        "|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---|---:|---|",
    ]
    for r in report["rows"]:
        if r["status"] != "measured":
            lines.append(f"| {r['row']} | | | | | | | | | | | | | {r['status']} |")
            continue
        lo, hi = r["ci95"]
        short = f" ({len(r['short_blocks'])} short)" if r["short_blocks"] else ""
        lines.append(
            f"| {r['row']} | {iters_text(r)} |"
            f" {r['block_s_min']:.2f} / {r['block_s_median']:.2f}{short} |"
            f" {r['median_a_s'] * 1e3:.3f} | {r['median_b_s'] * 1e3:.3f} |"
            f" {pct(r['statistic'])} | [{pct(lo)}, {pct(hi)}] | {r['drift'] * 100:.2f} |"
            f" {r['rss_a1_max_bytes'] / 2**20:.1f} / {r['rss_b_max_bytes'] / 2**20:.1f} |"
            f" {r['rss_allowance_bytes'] / 2**20:.1f} | {r['rss_balanced_ratio']:.3f} |"
            f" {'same' if r['digest_match'] else 'DIFFER'} | {r['load1_max']:.1f} |"
            f" {r['verdict']} |"
        )
    measured = [r for r in report["rows"] if r["status"] == "measured"]
    if measured:
        abs_stat = sorted(abs(r["statistic"]) for r in measured)
        widths = sorted(r["ci95"][1] - r["ci95"][0] for r in measured)
        drifts = sorted(r["drift"] for r in measured)
        counts = {
            v: sum(r["verdict"] == v for r in measured)
            for v in ("pass", "fail", "inconclusive")
        }

        def rows_where(test, text=lambda r: r["row"]) -> str:
            return ", ".join(text(r) for r in measured if test(r)) or "none"

        lines += [
            "",
            f"Rows measured: {len(measured)} ({counts['pass']} pass, {counts['fail']} fail,"
            f" {counts['inconclusive']} inconclusive); not measured:"
            f" {len(report['rows']) - len(measured)}.",
            f"Noise floor over measured rows: |statistic| median {median(abs_stat) * 100:.2f}%,"
            f" max {abs_stat[-1] * 100:.2f}%; CI width median {median(widths) * 100:.2f}%,"
            f" max {widths[-1] * 100:.2f}%; drift median {median(drifts) * 100:.2f}%,"
            f" max {drifts[-1] * 100:.2f}%.",
            f"Rows with timed blocks below {min_block} s (fail): "
            + rows_where(
                lambda r: r["short_blocks"],
                lambda r: f"{r['row']} ({len(r['short_blocks'])} blocks, shortest"
                f" {r['block_s_min']:.2f} s)",
            )
            + ".",
            "Bimodal sides (same-side clusters, fail): "
            + rows_where(
                lambda r: r["bimodal"],
                lambda r: f"{r['row']} ({', '.join(r['bimodal'])}: gap"
                f" {max(r['clusters'][x]['gap'] for x in r['bimodal']) * 100:.1f}%)",
            )
            + ".",
            "Memory failures: "
            + rows_where(lambda r: r["rss_verdict"] == "fail")
            + ".",
            "Inconclusive rows: "
            + rows_where(lambda r: r["verdict"] == "inconclusive")
            + ".",
            "Digest mismatches (fail): "
            + rows_where(lambda r: not r["digest_match"])
            + ".",
            "",
            "Timed blocks per side (ms per iteration; A in triple order A1, A2):",
            "",
        ]
        for r in measured:
            lines.append(
                f"* {r['row']}: A "
                + ", ".join(f"{x * 1e3:.4g}" for x in r["blocks_a_s"])
                + "; B "
                + ", ".join(f"{x * 1e3:.4g}" for x in r["blocks_b_s"])
            )
        lines += [
            "",
            "Peak RSS per timed block (MiB, triple order; the memory rule compares the largest"
            " A1 with the largest B; A2 is a diagnostic):",
            "",
        ]
        for r in measured:
            peaks = {
                k: ", ".join(
                    f"{t[k]['rss_tree_bytes'] / 2**20:.1f}" for t in r["triples"]
                )
                for k in ("A1", "B", "A2")
            }
            lines.append(
                f"* {r['row']}: A1 {peaks['A1']}; B {peaks['B']}; A2 {peaks['A2']}"
            )
    lines += ["", *plan_lines(report), "", *row_choice_lines(report)]
    return "\n".join(lines) + "\n"


def number_option(kind: type, domain: str):
    """An argparse ``type``: ``kind(text)`` in ``domain`` (``domain_problem``)."""

    def parse(text: str):
        value = kind(text)
        if bad := domain_problem("value", value, NUMBER, domain):
            raise argparse.ArgumentTypeError(bad)
        return value

    return parse


class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise Refused(message)


def parser() -> argparse.ArgumentParser:
    # Options by their full names only, each at most once (main): argparse would otherwise take a
    # prefix of a name for the option, and the last of an option given twice.
    ap = ArgumentParser(
        allow_abbrev=False,
        description=__doc__.split("\n\n")[0],
        epilog="Exit status: 0 the run completed, whatever its verdict (summary.md and"
        " report.json give the verdict and the conformance), 1 aborted, 2 refused before the"
        " first block, 130 interrupted (Ctrl-C) before the first block.",
    )
    ap.add_argument("--old-bin", required=True, help="the reference harness (side A)")
    ap.add_argument("--new-bin", required=True, help="the candidate harness (side B)")
    ap.add_argument("--old-label", default="old")
    ap.add_argument("--new-label", default="new")
    ap.add_argument(
        "--rows",
        default="all",
        help="comma list of row/group globs (spaces around a glob ignored, case-sensitive); "
        "a glob that matches no row or group is an error (the groups that exist are "
        "listed)",
    )
    ap.add_argument(
        "--triples",
        type=number_option(int, POSITIVE),
        default=PROTOCOL_TRIPLES,
        help=f"A-B-A triples per row (protocol: {PROTOCOL_TRIPLES}; fewer needs "
        "--allow-short)",
    )
    ap.add_argument("--tolerance", type=number_option(float, NONNEGATIVE), default=0.03)
    ap.add_argument(
        "--budget",
        type=number_option(float, POSITIVE),
        default=1.10,
        help="RSS budget factor",
    )
    ap.add_argument(
        "--rss-floor-mib",
        type=number_option(float, NONNEGATIVE),
        default=RSS_FLOOR_MIB,
        help="absolute RSS allowance: allowance = max(budget * ref, ref + floor)",
    )
    ap.add_argument(
        "--target-seconds",
        type=number_option(float, POSITIVE),
        default=PROTOCOL_BLOCK_SECONDS,
        help="minimum measured timed-block length (protocol: 2 s; shorter needs "
        "--allow-short); every shorter block fails its row",
    )
    ap.add_argument(
        "--calibration-margin",
        type=number_option(float, NONNEGATIVE),
        default=0.10,
        help="iterations are calibrated for target * (1 + margin)",
    )
    ap.add_argument(
        "--allow-short",
        action="store_true",
        help="permit fewer triples or shorter blocks than the protocol; the summary "
        "labels the run as not gate evidence",
    )
    ap.add_argument(
        "--expect-meta",
        help="a previous meta.txt: refuse to run unless every hash matches it",
    )
    ap.add_argument("--resamples", type=number_option(int, POSITIVE), default=10000)
    ap.add_argument("--seed", type=number_option(int, ANY), default=20261003)
    ap.add_argument(
        "--iters",
        help="row=N,... (both sides) or row=NA:NB,... (per side) pinned iteration counts; "
        "every row a selected one, none twice",
    )
    ap.add_argument(
        "--out",
        help="output directory, absent or empty (default dev/perf/results/<timestamp>)",
    )
    return ap


def main_run() -> int:
    if sys.version_info < PYTHON_MIN:
        raise Refused(
            f"aba.py needs Python {'.'.join(map(str, PYTHON_MIN))} or newer, not "
            f"{platform.python_version()} ({sys.executable})"
        )
    ap = parser()
    args = ap.parse_args()
    # Every option at most once: argparse keeps the last of two, so one would be ignored. After a
    # successful parse an argument that starts with "--" is an option, alone or joined to its
    # value by "=" (argparse takes a separate value that starts with "--" for an option).
    given = [a.split("=", 1)[0] for a in sys.argv[1:] if a.startswith("--")]
    if twice := sorted({a for a in given if given.count(a) > 1}):
        ap.error(f"{', '.join(twice)} given more than once")
    # The labels name the sides in every block, stderr log and meta.txt.
    for option, label in (
        ("--old-label", args.old_label),
        ("--new-label", args.new_label),
    ):
        if not LABEL.fullmatch(label):
            ap.error(
                f"{option} {label!r}: a label is 1 to 32 of the characters a-z, 0-9, '_' and"
                " '-'"
            )
    if args.old_label == args.new_label:
        ap.error("--old-label and --new-label must be two different names")
    if set(args.rows) - ROWS_TEXT:
        ap.error(
            f"--rows {args.rows!r}: only the characters of row names (a-z, 0-9, '.', '_',"
            " '-'), the glob characters '*', '?', '[', ']' and '!', commas and spaces"
        )
    for option, path in (("--out", args.out), ("--expect-meta", args.expect_meta)):
        if path is not None and not path.isprintable():
            ap.error(
                f"{option} {path!r}: a path with a character that is not printable"
            )
    short = []
    if args.triples < PROTOCOL_TRIPLES:
        short.append(f"--triples {args.triples} < {PROTOCOL_TRIPLES}")
    if args.target_seconds < PROTOCOL_BLOCK_SECONDS:
        short.append(
            f"--target-seconds {args.target_seconds} < {PROTOCOL_BLOCK_SECONDS}"
        )
    if short and not args.allow_short:
        ap.error(
            "below the protocol: "
            + "; ".join(short)
            + " (pass --allow-short for a run labelled as not gate evidence)"
        )
    # The options' own domains hold (number_option); the values derived from them must be
    # finite too (a huge --target-seconds or --rss-floor-mib overflows).
    for what, value in (
        (
            "--target-seconds * (1 + --calibration-margin)",
            args.target_seconds * (1 + args.calibration_margin),
        ),
        ("--rss-floor-mib in bytes", args.rss_floor_mib * 2**20),
    ):
        if not math.isfinite(value):
            ap.error(f"{what} is {value!r}, not finite")
    old = make_side("A", args.old_label, args.old_bin)
    new = make_side("B", args.new_label, args.new_bin)
    return run(args, old, new)


def refusal_outputs(reason: str) -> None:
    """A refusal's minimal report, when an unambiguous, printable, unused output path is writable.
    An invalid or occupied path is never used; failure to write it is reported on stderr."""
    paths = []
    argv = sys.argv[1:]
    for i, arg in enumerate(argv):
        if arg.startswith("--out="):
            paths.append(arg[6:])
        elif arg == "--out":
            paths.append(
                argv[i + 1]
                if i + 1 < len(argv) and not argv[i + 1].startswith("--")
                else ""
            )
    if len(paths) > 1 or (paths and not paths[0].isprintable()):
        warn("cannot write refusal outputs: --out is ambiguous or not a printable path")
        return
    out = (
        Path(paths[0])
        if paths
        else HERE / "results" / dt.datetime.now().strftime("%Y%m%dT%H%M%S%f")
    )
    try:
        if out.exists() and (not out.is_dir() or any(out.iterdir())):
            warn(f"cannot write refusal outputs: {out} must be absent or empty")
            return
        out.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        warn(f"cannot write refusal outputs: {exc}")
        return
    report = {
        "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "aborted": reason,
        "rows": [],
    }
    for name, text in (
        ("report.json", report_text(minimal_report(report))),
        ("summary.md", minimal_summary(report)),
    ):
        try:
            atomic_text(out / name, text)
        except Exception as exc:
            warn(f"cannot write {out / name}: {exc}")


def main() -> int:
    with interrupts_once() as interrupts:
        try:
            return main_run()
        except KeyboardInterrupt as exc:
            interrupts.hold()
            reason = "interrupted (SIGINT, Ctrl-C)" + (
                f" during {exc}" if exc.args else ""
            )
            status = 130
        except Exception as exc:
            interrupts.hold()
            reason = (
                str(exc) if isinstance(exc, Refused) else f"{type(exc).__name__}: {exc}"
            )
            status = 2
        warn(f"error: {reason}")
        refusal_outputs(reason)
        return status


if __name__ == "__main__":
    raise SystemExit(main())
