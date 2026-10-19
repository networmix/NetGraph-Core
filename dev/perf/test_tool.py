"""Tests of the timing tool itself (README.md, "The tool's own tests"); nothing is timed.

    venv/bin/python -B dev/perf/test_tool.py [-v]        # or: make perf-tool-test

The rules README.md states are tested here:

* ``aba.py`` end to end with two fake harnesses (Python scripts that answer ``aba-list``,
  ``aba-calibrate`` and ``aba-run`` like the real one; ``Driver``), and its parts directly (the
  memory rule, the verdicts, the bootstrap interval, the cluster rule, the row selection, the
  harness responses, the block line, derived values that are not finite, the SIGINT handler,
  the hashes of the tool's sources);
* ``build.sh`` with a fake compiler, ``xcrun`` and ``vtool`` in a scratch git repository
  (``BuildScript``);
* the harness's driver and the checks of ``aba_rows.hpp``, compiled unoptimised against this
  checkout's Core with a few test rows (``Harness``);
* the frozen files against the values compiled into ``aba_rows.hpp`` (``FrozenFiles``).
"""

from __future__ import annotations

import concurrent.futures
import contextlib
import csv
import glob
import hashlib
import io
import json
import math
import os
import platform
import random
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from fractions import Fraction
from pathlib import Path
from unittest import mock

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ABA = HERE / "aba.py"
sys.path.insert(0, str(HERE))

import aba  # noqa: E402

MIB = 2**20


def sha256sums(directory: Path) -> dict[str, str]:
    """``directory``'s SHA256SUMS as {file name: sha256}; ValueError for a file named twice."""
    lines = (directory / "SHA256SUMS").read_text().splitlines()
    out = dict(reversed(ln.split("  ", 1)) for ln in lines)
    if len(out) != len(lines):
        raise ValueError(f"{directory}/SHA256SUMS names a file twice")
    return out


def block(rss: float) -> dict:
    return {
        "rss_self_bytes": rss,
        "rss_children_bytes": 0,
        "rss_wait4_bytes": rss,
        "rss_tree_bytes": rss,
    }


def triples(a1: list[float], b: list[float], a2: list[float]) -> list[dict]:
    return [
        {"A1": block(x), "B": block(y), "A2": block(z)}
        for x, y, z in zip(a1, b, a2, strict=True)
    ]


def rule(a1, b, a2=None) -> dict:
    return aba.rss_balanced(triples(a1, b, a2 or a1), 1.10, 8 * MIB)


def exact_pass(reference: int, candidate: int) -> bool:
    """The allowance in exact arithmetic (1.10 = 11/10)."""
    return candidate <= max(Fraction(11, 10) * reference, reference + 8 * MIB)


class MemoryBoundary(unittest.TestCase):
    def check(self, ref: int, allowance: int) -> None:
        a1 = [ref - 4096] * 6 + [ref]
        at = rule(a1, [allowance] + [ref] * 6)
        over = rule(a1, [allowance + 1] + [ref] * 6)
        self.assertEqual((at["rss_verdict"], over["rss_verdict"]), ("pass", "fail"))
        self.assertEqual(at["rss_a1_max_bytes"], ref)
        self.assertEqual(at["rss_b_max_bytes"], allowance)
        self.assertTrue(exact_pass(ref, allowance))
        self.assertFalse(exact_pass(ref, allowance + 1))

    def test_floor_branch_equal_and_one_byte_over(self):
        # 10 MiB: 1.10 x R = 11 MiB < R + 8 MiB = 18 MiB.
        self.check(10 * MIB, 18 * MIB)
        self.assertEqual(rule([10 * MIB] * 7, [0] * 7)["rss_allowance_bytes"], 18 * MIB)

    def test_ratio_branch_equal_and_one_byte_over(self):
        # 100 MiB: 1.10 x R = 110 MiB > R + 8 MiB; the float allowance is a hair above
        # 115,343,360, which changes no integer verdict.
        self.check(100 * MIB, 110 * MIB)
        self.assertGreaterEqual(
            rule([100 * MIB] * 7, [0] * 7)["rss_allowance_bytes"], 110 * MIB
        )

    def test_branches_meet_at_80_mib(self):
        self.check(80 * MIB, 88 * MIB)

    def test_floor_applies_above_the_ratio(self):
        # +7 MiB on 40 MiB is ratio 1.175: above 1.10 but within the 8 MiB floor. A check
        # of the ratio against 1.10 alone would wrongly fail it.
        r = rule([40 * MIB] * 7, [47 * MIB] * 7)
        self.assertGreater(r["rss_balanced_ratio"], 1.10)
        self.assertEqual(r["rss_verdict"], "pass")
        self.assertEqual(r["rss_allowance_bytes"], 48 * MIB)

    def test_ratio_branch_above_the_floor(self):
        # 200 MiB: +16 MiB is within 10 percent (220 MiB), +21 MiB is not.
        self.assertEqual(rule([200 * MIB] * 7, [216 * MIB] * 7)["rss_verdict"], "pass")
        self.assertEqual(rule([200 * MIB] * 7, [221 * MIB] * 7)["rss_verdict"], "fail")

    def test_float_agrees_with_exact_arithmetic_on_page_multiples(self):
        page = 16384
        for ref in range(page, 400 * MIB, 7 * MIB + 3 * page):
            r = rule([ref] * 7, [0] * 7)
            edge = int(r["rss_allowance_bytes"]) // page * page
            for cand in (edge - page, edge, edge + page):
                got = rule([ref] * 7, [cand] * 7)["rss_verdict"] == "pass"
                self.assertEqual(got, exact_pass(ref, cand), (ref, cand))


class MemoryFixedA1(unittest.TestCase):
    def test_larger_a2_cannot_rescue(self):
        r = rule([40 * MIB] * 7, [50 * MIB] * 7, [40 * MIB] * 6 + [60 * MIB])
        self.assertEqual(r["rss_verdict"], "fail")
        self.assertEqual(r["rss_a2_max_bytes"], 60 * MIB)
        self.assertEqual(r["rss_a1_max_bytes"], 40 * MIB)

    def test_slot_not_value_selects_the_reference(self):
        # The same seven pairs of A values with the large one moved from A2 to A1: the
        # reference follows the slot.
        a_low, a_high = [40 * MIB] * 7, [40 * MIB] * 6 + [60 * MIB]
        b = [50 * MIB] * 7
        self.assertEqual(rule(a_low, b, a_high)["rss_verdict"], "fail")
        self.assertEqual(rule(a_high, b, a_low)["rss_verdict"], "pass")

    def test_one_high_b_block_fails(self):
        r = rule([40 * MIB] * 7, [40 * MIB] * 6 + [49 * MIB])
        self.assertEqual(r["rss_verdict"], "fail")


class Verdicts(unittest.TestCase):
    def test_row_verdict_combines(self):
        # any failing check fails the row; else timing's inconclusive; else pass
        for i in range(5):
            checks = ["pass"] * 5
            checks[i] = "fail"
            self.assertEqual(aba.row_verdict(*checks), "fail", i)
        self.assertEqual(
            aba.row_verdict("inconclusive", "fail", "pass", "pass", "pass"), "fail"
        )
        self.assertEqual(
            aba.row_verdict("inconclusive", "pass", "pass", "pass", "pass"),
            "inconclusive",
        )
        self.assertEqual(aba.row_verdict(*["pass"] * 5), "pass")


def expected_interval(x: list[float], resamples: int, seed: int) -> tuple[float, float]:
    """The bootstrap interval computed independently of aba.py: ``resamples`` draws of
    ``len(x)`` whole triples with ``random.Random(seed).randrange`` (the documented procedure),
    the middle element of each sorted draw (``len(x)`` is odd), and the 2.5 % and 97.5 %
    order statistics at floor(0.025 * (R - 1)) and ceil(0.975 * (R - 1))."""
    rng = random.Random(seed)
    meds = []
    for _ in range(resamples):
        draw = sorted(x[rng.randrange(len(x))] for _ in range(len(x)))
        meds.append(draw[len(draw) // 2])
    meds.sort()
    return (
        meds[math.floor(0.025 * (resamples - 1))],
        meds[math.ceil(0.975 * (resamples - 1))],
    )


def timed_block(
    median_s: float,
    digest: str = "d",
    rss: float = 40 * MIB,
    block_s: float = 2.5,
) -> dict:
    """A timed block as report.json stores it (2,000 iterations)."""
    return {
        **block(rss),
        "row": "w.one",
        "iterations": 2000,
        "median_s": median_s,
        "block_s": block_s,
        "setup_s": 0.01,
        "digest": digest,
        "load1_before": 0.0,
        "load1_after": 0.0,
        "side": "old",
        "wall_s": 2.6,
    }


def timed_triples(
    a1: list[float], b: list[float], a2: list[float], **block
) -> list[dict]:
    """Timed triples of these median iteration times (``timed_block``'s other fields from
    ``block``)."""
    return [
        {
            "A1": timed_block(x, **block),
            "B": timed_block(y, **block),
            "A2": timed_block(z, **block),
        }
        for x, y, z in zip(a1, b, a2, strict=True)
    ]


class Bootstrap(unittest.TestCase):
    SEED = 20261003
    # Seven per-triple statistics: median 0.01 (inside a 3 % tolerance), long right tail.
    X_TAIL = [-0.05, -0.04, 0.0, 0.01, 0.02, 0.10, 0.12]
    X_TIGHT = [-0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02]

    def test_endpoints_equal_an_independent_computation(self):
        for x in (self.X_TAIL, self.X_TIGHT):
            got = aba.bootstrap_interval(x, 10000, self.SEED)
            want = expected_interval(x, 10000, self.SEED)
            self.assertEqual(got, want)
            self.assertLess(got[0], got[1])  # nonconstant
        self.assertEqual(
            expected_interval(self.X_TAIL, 10000, self.SEED), (-0.04, 0.10)
        )

    def row(self, x: list[float]) -> dict:
        # A1 = A2 = 1 ms (no drift), B = 1 ms * (1 + x_t): r_t - 1 = x_t up to rounding.
        a = [1e-3] * len(x)
        triples = timed_triples(a, [1e-3 * (1 + v) for v in x], a)
        return aba.row_statistics(triples, 0.03, 1.10, 8 * MIB, 2.0, 10000, self.SEED)

    def test_upper_bound_decides_the_timing_verdict(self):
        tail = self.row(self.X_TAIL)
        self.assertAlmostEqual(tail["statistic"], 0.01)  # the median alone would pass
        self.assertAlmostEqual(tail["ci95"][1], 0.10)
        self.assertEqual(tail["time_verdict"], "fail")
        self.assertEqual(tail["verdict"], "fail")
        tight = self.row(self.X_TIGHT)
        self.assertAlmostEqual(tight["ci95"][1], 0.01)
        self.assertLessEqual(tight["ci95"][1], 0.03)
        self.assertEqual(tight["time_verdict"], "pass")


class TimingRules(unittest.TestCase):
    """The timing statistics and the block-length rule on exact binary fractions (2^-10 s and
    1.25 x 2^-10 s), so every boundary is met exactly."""

    U = 2.0**-10

    def test_drift_takes_precedence_over_a_failing_interval(self):
        data = timed_triples([1.0] * 7, [2.0] * 7, [1.1] * 7)
        rec = aba.row_statistics(data, 0.03, 1.10, 8 * MIB, 2.0, 100, 7)
        self.assertGreater(rec["ci95"][0], 0.03)
        self.assertGreater(rec["drift"], 0.03)
        self.assertEqual(rec["time_verdict"], "inconclusive")

    def stats(self, a1, b, a2, tolerance=0.25, block_s=2.5):
        a1, b, a2 = ([v * self.U for v in side] for side in (a1, b, a2))
        triples = timed_triples(a1, b, a2, block_s=block_s)
        return aba.row_statistics(triples, tolerance, 1.10, 8 * MIB, 2.0, 200, 1)

    def test_ratio_against_the_mean_of_a1_and_a2(self):
        r = self.stats([1.0] * 7, [1.25] * 7, [1.5] * 7, tolerance=1.0)
        self.assertEqual(r["statistic"], 1.25 / 1.25 - 1)
        self.assertEqual(r["ratios"], [1.0] * 7)

    def test_drift_is_the_median_of_a1_over_a2(self):
        # |A1/A2 - 1| = 0, 0, 0, 0.5, 0.5, 0.5, 0.5: median 0.5 (mean 0.29; |A2/A1 - 1| 0.33)
        r = self.stats([1.0] * 3 + [1.5] * 4, [1.0] * 7, [1.0] * 7, tolerance=1.0)
        self.assertEqual(r["drift"], 0.5)

    def test_boundaries_are_inclusive(self):
        # the interval's upper bound equal to the tolerance passes
        r = self.stats([1.0] * 7, [1.25] * 7, [1.0] * 7)
        self.assertEqual((r["ci95"][1], r["time_verdict"]), (0.25, "pass"))
        # a drift equal to the tolerance is not inconclusive (B 1.125 = mean of A1 and A2)
        r = self.stats([1.25] * 7, [1.125] * 7, [1.0] * 7)
        self.assertEqual((r["drift"], r["time_verdict"]), (0.25, "pass"))
        # a block as long as the minimum is not short
        r = self.stats([1.0] * 7, [1.0] * 7, [1.0] * 7, block_s=2.0)
        self.assertEqual((r["short_blocks"], r["block_verdict"]), ([], "pass"))


class Clusters(unittest.TestCase):
    """The same-side cluster rule: a gap above the tolerance between consecutive sorted blocks,
    with at least two blocks on each side of it."""

    def test_two_blocks_on_each_side(self):
        for times, want in (
            ([1.0] * 5 + [1.5] * 2, (5, 2)),  # two high blocks: bimodal
            ([1.0] * 2 + [1.5] * 5, (2, 5)),  # two low blocks: bimodal
            ([1.0] * 6 + [1.5], None),  # one outlier above: not a cluster
            ([0.5] + [1.0] * 6, None),  # one outlier below: not a cluster
            ([1.0] * 4 + [1.25] * 3, None),  # a gap of exactly the tolerance, not above
        ):
            got = aba.side_clusters(times, 0.25)
            self.assertEqual(
                None if got is None else (got["low"]["blocks"], got["high"]["blocks"]),
                want,
                times,
            )

    def test_the_widest_split_is_reported(self):
        got = aba.side_clusters([1.0, 1.0, 1.1, 1.1, 1.5, 1.5], 0.03)
        self.assertAlmostEqual(got["gap"], 1.5 / 1.1 - 1)
        self.assertEqual((got["low"]["blocks"], got["high"]["blocks"]), (4, 2))


class DerivedValues(unittest.TestCase):
    """Finite block numbers whose derived values are not finite abort the run (RuntimeError)
    before any statistic or verdict is returned."""

    def stats(self, a1, b, a2, rss=40 * MIB):
        triples = timed_triples(a1, b, a2, rss=rss)
        return aba.row_statistics(triples, 0.03, 1.10, 8 * MIB, 2.0, 200, 1)

    def aborts(self, why: str, *args, **kwargs) -> None:
        with self.assertRaises(RuntimeError) as ctx:
            self.stats(*args, **kwargs)
        self.assertIn(why, str(ctx.exception))

    def test_ordinary_values_pass(self):
        self.assertEqual(
            self.stats([1e-3] * 7, [1e-3] * 7, [1e-3] * 7)["verdict"], "pass"
        )

    def test_mean_overflow_in_one_triple(self):
        # Only the mean of A1 and A2 of triple 7 overflows (1.7e308 + 1e307): every value the
        # row records stays finite, so only the mean's own check stops the ratio 1e300 / inf = 0.
        a1 = [1e300] * 6 + [1.7e308]
        a2 = [1e300] * 6 + [1e307]
        self.aborts("the mean of A1 and A2 is inf, not finite", a1, [1e300] * 7, a2)

    def test_ratio_and_allowance_overflow(self):
        self.aborts(
            "row statistic ratios is inf", [1e-300] * 7, [1e10] * 7, [1e-300] * 7
        )
        # a peak near the largest double: 1.10 * peak overflows
        self.aborts(
            "row statistic rss_allowance_bytes is inf",
            [1e-3] * 7,
            [1e-3] * 7,
            [1e-3] * 7,
            rss=1.7e308,
        )

    def test_calibrated_count_overflow(self):
        # the count of the side whose warm-up median is tiny, named with its row
        ok = {"row": "r", "side": "old", "median_s": 1e-3}
        tiny = {"row": "r", "side": "new", "median_s": 5e-324}
        self.assertEqual(aba.side_counts(ok, ok, 2.0, 0.1, None)[:2], (2200, 2200))
        for warm_a, warm_b, side in ((tiny, ok, "old"), (ok, tiny, "new")):
            with self.assertRaises(RuntimeError) as ctx:
                aba.side_counts(
                    {**warm_a, "side": "old"}, {**warm_b, "side": "new"}, 2.0, 0.1, None
                )
            self.assertIn(
                f"r {side}: the calibrated iteration count is inf", str(ctx.exception)
            )


ROWS = {
    "a.one": {"group": "a"},
    "a.two": {"group": "a"},
    "b.one": {"group": "b"},
    "c.one": {"group": "c"},
}


class RowResolution(unittest.TestCase):
    def plan(self, rows: str) -> dict:
        return aba.choose_rows(ROWS, rows)

    def test_row_name_group_and_all(self):
        self.assertEqual(self.plan("a.one")["selected"], ["a.one"])
        self.assertEqual(self.plan("a")["selected"], ["a.one", "a.two"])
        self.assertEqual(self.plan("a.*")["selected"], ["a.one", "a.two"])
        self.assertEqual(self.plan("all")["selected"], list(ROWS))
        self.assertEqual(self.plan("")["selected"], list(ROWS))

    def test_glob_characters(self):
        # '?' one character, '[...]' one of a set, '[!...]' one outside it
        self.assertEqual(self.plan("?.one")["selected"], ["a.one", "b.one", "c.one"])
        self.assertEqual(self.plan("[ab].t*")["selected"], ["a.two"])
        self.assertEqual(self.plan("[!a].one")["selected"], ["b.one", "c.one"])

    def test_glob_matching_nothing_is_an_error(self):
        sel = self.plan("a.*,placement.*")
        self.assertEqual(
            sel["errors"], ["--rows glob 'placement.*' matches no row or group"]
        )
        text = aba.row_choice_error_text(sel, ROWS)
        self.assertIn("a (2): a.one, a.two", text)
        self.assertIn("b (1): b.one", text)
        self.assertEqual(self.plan(" , ")["errors"], ["--rows ' , ' names no glob"])

    def test_records(self):
        sel = self.plan("a.one,b.one")
        self.assertEqual(sel["errors"], [])
        first = sel["globs"][0]
        self.assertEqual(first["group_rows_not_matched"], ["a.two"])
        self.assertEqual(sel["not_selected"], {"a": ["a.two"], "c": ["c.one"]})

    def test_rows_parsing(self):
        cases = {
            " a.* , b.one ": ["a.one", "a.two", "b.one"],
            "a.one,,c": ["a.one", "c.one"],
            "all": list(ROWS),
        }
        for text, want in cases.items():
            self.assertEqual(self.plan(text)["selected"], want, text)
        self.assertEqual(aba.split_patterns(" a.* , b.one "), ["a.*", "b.one"])

    def test_both_sides_list_the_same_rows(self):
        aba.same_rows(ROWS, dict(ROWS))
        for old, new, why in (
            (ROWS, {**ROWS, "d.one": {"group": "d"}}, "only the candidate ['d.one']"),
            (
                {**ROWS, "d.one": {"group": "d"}},
                ROWS,
                "only the reference has ['d.one']",
            ),
            (ROWS, {**ROWS, "c.one": {"group": "e"}}, "other records ['c.one']"),
            (ROWS, dict(reversed(ROWS.items())), "another order"),
        ):
            with self.assertRaises(aba.Refused) as ctx:
                aba.same_rows(old, new)
            self.assertIn(why, str(ctx.exception))


class Responses(unittest.TestCase):
    """A harness's whole response, never a usable line picked out of it."""

    LINE = json.dumps(
        {
            "row": "r",
            "iterations": 5,
            "median_s": 0.001,
            "block_s": 0.005,
            "setup_s": 0.0,
            "digest": "00000000000000ff",
            "rss_self_bytes": 41943040,
            "rss_children_bytes": 0,
        }
    )

    def refused(self, out: bytes, why: str) -> None:
        with self.assertRaises(aba.ResponseError) as ctx:
            aba.parse_block(out, "fake run r", "r")
        self.assertIn(why, str(ctx.exception))

    def test_one_block_record_and_nothing_else(self):
        line = self.LINE
        self.assertEqual(
            aba.parse_block(f"{line}\n".encode(), "fake run r", "r")["median_s"], 0.001
        )
        negative = line.replace('"median_s": 0.001', '"median_s": -0.001')
        for out, why in (
            # a broken record, a negative median or text around the valid record
            (f'{{"row": "r", "median_s": 0.001,\n{line}\n', "; line 1 of 2): {"),
            (f"{negative}\n{line}\n", "2 output records, exactly one expected"),
            (f"{line}\nNOT JSON\n", "; line 2 of 2): NOT JSON"),
            (f"{line}\n\n", "; line 2 of 2)"),
            (f"[1]\n{line}\n", "not a JSON object; line 1 of 2"),
            (f"{line}", "not ended by a newline"),
            ("", "0 output records, exactly one expected"),
            # one record, but the key twice (the valid value last) or another row's
            ('{"median_s": -1.0, ' + line[1:] + "\n", "key 'median_s' twice"),
            (line.replace('"row": "r"', '"row": "q"') + "\n", "names row 'q', not 'r'"),
        ):
            self.refused(out.encode(), why)
        self.refused(b"\xff\n", "the output is not UTF-8 text")

    def test_row_list(self):
        def rows(*lines: str) -> dict:
            return aba.row_list(("".join(f"{ln}\n" for ln in lines)).encode(), "list")

        good = [json.dumps({"row": n, "group": "g", "desc": "d"}) for n in ("a", "b")]
        self.assertEqual(list(rows(*good)), ["a", "b"])
        for lines, why in (
            ((*good, good[0]), "row 'a' listed twice"),
            (("junk", *good), "malformed output line"),
            ((*good, json.dumps({"row": "c"})), "without a string group"),
            ((*good, json.dumps({"row": 3, "group": "g"})), "without a string row"),
            ((*good, json.dumps({"row": "c", "group": "g", "desc": 3})), "string desc"),
            (
                (*good, json.dumps({"row": "c", "group": "g", "queries": 1})),
                "string queries",
            ),
            (
                (*good, json.dumps({"row": "c", "group": "g", "x": float("nan")})),
                'list: row c: unknown key \'x\': {"row": "c", "group": "g", "x": NaN}',
            ),
            ((), "no row listed"),
        ):
            with self.assertRaises(aba.ResponseError) as ctx:
                rows(*lines)
            self.assertIn(why, str(ctx.exception))
        # a row's name and its group: 1 to 64 of a-z, 0-9, '.', '_' and '-'
        names = ("a.b_c-9", "x" * 64)
        self.assertEqual(
            list(rows(*(json.dumps({"row": n, "group": n}) for n in names))),
            list(names),
        )
        for name in (
            "",
            "probe/row",
            "Upper",
            "x" * 65,
            "tab\there",
            "\ud800",
            "\u00e9",
        ):
            for key in ("row", "group"):
                rec = {"row": "a", "group": "g", key: name}
                with self.assertRaises(aba.ResponseError, msg=(key, name)) as ctx:
                    rows(json.dumps(rec))
                self.assertIn(f"{key} {name!r} is not 1 to 64 of", str(ctx.exception))


class BlockLine(unittest.TestCase):
    """``aba.parse_block``, the harness boundary: every number of a block line of its type,
    finite and in its domain, else the run aborts."""

    LINE = {
        "row": "r",
        "iterations": 5,
        "median_s": 0.001,
        "mean_s": 0.001,
        "min_s": 0.0009,
        "p10_s": 0.0009,
        "p90_s": 0.0011,
        "max_s": 0.002,
        "block_s": 0.005,
        "setup_s": 0.01,
        "digest": "00000000000000ff",
        "rss_self_bytes": 41943040,
        "rss_children_bytes": 0,
    }

    def parse(self, **change) -> dict:
        line = json.dumps({**self.LINE, **change}) + "\n"
        return aba.parse_block(line.encode(), "fake run r", "r")

    def test_valid_line_and_the_zero_cases(self):
        # The writer refuses non-finite values even if a future caller misses an input check.
        for value in (float("nan"), float("inf"), -float("inf")):
            with self.assertRaises(ValueError):
                aba.report_text({"value": value})
        rec = self.parse()
        self.assertEqual(rec, self.LINE)
        # every number a float but the count, so arithmetic on them overflows to an infinity
        # (caught as not finite), never to an exception
        self.assertEqual(
            {k for k, v in rec.items() if type(v) is float},
            set(self.LINE) - {"row", "iterations", "digest"},
        )
        self.assertIs(type(rec["iterations"]), int)
        # 0 where a real block prints it: block_s and setup_s (rounded to microseconds), the
        # children's peak (no child), the unread statistics (an iteration under a clock tick)
        try:
            self.parse(
                block_s=0.0, setup_s=0.0, rss_children_bytes=0, min_s=0.0, p10_s=0
            )
        except RuntimeError as exc:
            self.fail(f"a real 0 refused: {exc}")

    def test_numbers_outside_their_domain_abort(self):
        for change, why in (
            ({"median_s": -0.001}, "median_s -0.001 is not positive"),
            ({"median_s": 0}, "median_s 0 is not positive"),
            ({"iterations": 0}, "iterations 0 is not positive"),
            ({"rss_self_bytes": 0}, "rss_self_bytes 0 is not positive"),
            ({"block_s": -0.5}, "block_s -0.5 is negative"),
            ({"setup_s": -1e-06}, "setup_s -1e-06 is negative"),
            ({"rss_children_bytes": -1}, "rss_children_bytes -1 is negative"),
            ({"max_s": -1.0}, "max_s -1.0 is negative"),
            ({"p90_s": "fast"}, "missing or mistyped p90_s"),
            ({"iterations": 2.0}, "missing or mistyped iterations"),
            ({"iterations": True}, "missing or mistyped iterations"),
            ({"digest": 255}, "missing or mistyped digest"),
            (
                {"digest": "00000000000000FF"},
                "digest '00000000000000FF' is not 16 lowercase hexadecimal digits",
            ),
            ({"digest": "ff"}, "digest 'ff' is not 16 lowercase hexadecimal digits"),
            (
                {"digest": "\ud800" * 16},
                f"digest {chr(0xD800) * 16!r} is not 16 lowercase hexadecimal digits",
            ),
            ({"rss_children_bytes": None}, "missing or mistyped rss_children_bytes"),
            # JSON integers have no bound: one beyond the largest float is not finite
            (
                {"rss_self_bytes": 10**400},
                "rss_self_bytes is not finite (an integer of 401 digits)",
            ),
            ({"min_s": -(10**309)}, "min_s is not finite (an integer of 310 digits)"),
        ):
            with self.assertRaises(RuntimeError) as ctx:
                self.parse(**change)
            self.assertIn(
                f"fake run r: malformed output line ({why}", str(ctx.exception)
            )
        # NaN and Infinity are JSON literals json.loads reads
        for literal, why in (
            ("NaN", "median_s nan is not finite"),
            ("Infinity", "median_s inf is not finite"),
            ("-Infinity", "median_s -inf is not finite"),
        ):
            line = json.dumps(self.LINE).replace(
                '"median_s": 0.001', f'"median_s": {literal}'
            )
            with self.assertRaises(RuntimeError) as ctx:
                aba.parse_block((line + "\n").encode(), "fake run r", "r")
            self.assertIn(why, str(ctx.exception))

    def test_keys_outside_the_rule_abort(self):
        # A block line carries the keys of its rule only: any other aborts the run whatever its
        # value (NaN; nested deeper than Python's writer goes), and no reason prints a nested
        # value (Python 3.14's reader reads this nesting; older ones refuse it as too deep).
        deep = "[" * 110000 + "]" * 110000
        line = json.dumps(self.LINE)[:-1]
        for extra, why in (
            (', "x": 1}', "unknown key 'x'"),
            (', "x": NaN}', "unknown key 'x'"),
            (', "x": ' + deep + "}", "unknown key 'x'"),
            (
                ', "a": 1, "b": 2, "c": 3, "d": 4, "e": 5, "f": 6}',
                "unknown keys 'a', 'b', 'c', 'd', 'e', ...",
            ),
            (', "' + "k" * 200 + '": 1}', "unknown key '" + "k" * 116 + "..."),
        ):
            with self.assertRaises(aba.ResponseError) as ctx:
                aba.parse_block((line + extra + "\n").encode(), "fake run r", "r")
            msg = str(ctx.exception)
            ok = f"fake run r: malformed output line ({why}" in msg
            self.assertTrue(ok or "JSON nested too deeply" in msg, msg)
        # a nested value where the rule reads a number or the row: refused by its type
        for key, why in (
            ("median_s", "malformed output line (missing or mistyped median_s): "),
            ("row", "the record names row an array, not 'r'"),
        ):
            text = json.dumps({**self.LINE, key: "@"}).replace('"@"', deep)
            with self.assertRaises(aba.ResponseError) as ctx:
                aba.parse_block((text + "\n").encode(), "fake run r", "r")
            msg = str(ctx.exception)
            self.assertTrue(why in msg or "JSON nested too deeply" in msg, msg)
        self.assertIn('"median_s": "<an array>"', aba.shown({"median_s": [[1]]}))


# The driver files build.sh compiles and records.
DRIVER = ("aba_api.hpp", "aba_rows.hpp", "harness.cpp")


def checkout_driver(perf: Path = HERE) -> dict[str, str]:
    """The ``driver_sha256`` build.sh writes for the driver in ``perf`` (this checkout's)."""
    return {
        f"dev/perf/{n}": hashlib.sha256((perf / n).read_bytes()).hexdigest()
        for n in DRIVER
    }


def fake_record(binary: Path, rev: str) -> dict:
    """A build record as build.sh writes it for ``binary`` (every field, this checkout's
    driver)."""
    return {
        "binary": str(binary),
        "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "rev": rev,
        "commit": "0" * 40,
        "dirty_core_paths": 0,
        "core_sources_sha256": "0" * 64,
        "driver_sha256": checkout_driver(),
        "compiler": "fake clang version 1.0",
        "flags": "-O3",
        "compiler_path": "/usr/bin/clang++",
        "cxx_env": None,
        "isysroot": None,
        "sdk": {
            "path": None,
            "version": None,
            "sdkroot_env": None,
            "source": "no xcrun SDK on this host: the compiler's own system headers",
        },
        "deployment_target": {
            "requested": None,
            "effective_minos": None,
            "binary_sdk": None,
            "platform": None,
            "raw": None,
            "source": "no LC_BUILD_VERSION read (no vtool, or it printed none)",
        },
        "built_by": "dev/perf/build.sh",
    }


class CreateOutputs(unittest.TestCase):
    """``aba.create_outputs``, the last step before the first block: an operating-system error
    or an interrupt there leaves nothing the run created (``meta.txt`` cut short by a full disk,
    ``--out`` and the parents it made)."""

    def test_nothing_left(self):
        def partial(path, sides, lines):
            Path(path).write_text("# aba.py provenance\n")
            raise exc

        saved = aba.write_meta
        self.addCleanup(setattr, aba, "write_meta", saved)
        aba.write_meta = partial
        for exc, caught in (
            (OSError(28, "No space left on device"), aba.Refused),
            (KeyboardInterrupt(), KeyboardInterrupt),
        ):
            with tempfile.TemporaryDirectory() as tmp:
                top = Path(tmp)
                with self.assertRaises(caught) as ctx:
                    aba.create_outputs(top / "a" / "out", [], [])
                self.assertEqual(list(top.iterdir()), [], exc)
                if caught is aba.Refused:
                    self.assertIn("No space left on device", str(ctx.exception))
                # an --out that existed, empty, stays and stays empty
                (top / "kept").mkdir()
                with self.assertRaises(caught):
                    aba.create_outputs(top / "kept", [], [])
                self.assertEqual(list((top / "kept").iterdir()), [], exc)


class BlockLogs(unittest.TestCase):
    """``aba.run_block``'s stderr log beside a log that earlier blocks of the same row, side
    and mode filled: an abort names it only when this block wrote to it, and a log still empty
    after a block is removed, also when the harness could not start or was interrupted."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.blocks = self.dir / "blocks"
        self.log = self.blocks / "r.new.aba-run.stderr"

    def tearDown(self):
        self.tmp.cleanup()

    def block(self, script: str) -> str:
        """Run one block of a harness made of ``script``; the abort's message."""
        path = self.dir / "harness"
        path.write_text(f"#!/bin/sh\n{script}\n")
        path.chmod(0o755)
        side = aba.Side("B", "new", str(path), aba.sha256_file(path))
        with self.assertRaises(RuntimeError) as ctx:
            aba.run_block(side, "aba-run", "r", "1", clean_env(), self.blocks)
        return str(ctx.exception)

    def test_named_only_when_written(self):
        self.blocks.mkdir()
        self.log.write_text("an earlier block\n")
        why = self.block("exit 3")
        self.assertEqual(why, "new aba-run r exited 3 (nothing on stderr)")
        self.assertEqual(self.log.read_text(), "an earlier block\n")
        why = self.block("echo this block >&2; exit 3")
        self.assertEqual(why, f"new aba-run r exited 3 (stderr: {self.log})")
        self.assertEqual(self.log.read_text(), "an earlier block\nthis block\n")
        why = self.block("echo this block again >&2; echo garbage")
        self.assertIn(f"new aba-run r (stderr: {self.log}): ", why)
        why = self.block("echo garbage")
        self.assertIn("new aba-run r (nothing on stderr): ", why)
        self.assertEqual(
            self.log.read_text(), "an earlier block\nthis block\nthis block again\n"
        )

    def test_empty_log_removed(self):
        self.assertEqual(
            self.block("exit 3"), "new aba-run r exited 3 (nothing on stderr)"
        )
        self.assertEqual(os.listdir(self.blocks), [])
        # a harness that cannot start (here a script without its interpreter)
        path = self.dir / "harness"
        path.write_text("#!/no/such/interpreter\n")
        side = aba.Side("B", "new", str(path), aba.sha256_file(path))
        with self.assertRaises(RuntimeError) as ctx:
            aba.run_block(side, "aba-run", "r", "1", clean_env(), self.blocks)
        self.assertIn("new aba-run r: cannot start", str(ctx.exception))
        self.assertEqual(os.listdir(self.blocks), [])
        # an interrupt while the harness runs
        saved = os.wait4
        self.addCleanup(setattr, os, "wait4", saved)

        def interrupted(*args):
            raise KeyboardInterrupt

        os.wait4 = interrupted
        path.write_text("#!/bin/sh\nexec sleep 30 >&-\n")  # stdout at its end at once
        side = aba.Side("B", "new", str(path), aba.sha256_file(path))
        with self.assertRaises(KeyboardInterrupt):
            aba.run_block(side, "aba-run", "r", "1", clean_env(), self.blocks)
        os.wait4 = saved
        self.assertEqual(os.listdir(self.blocks), [])


class InterruptsOnce(unittest.TestCase):
    """``aba.interrupts_once``, the run's SIGINT handler: the first interrupt raises
    KeyboardInterrupt, every later one, and every one after ``hold``, is only recorded, and the
    previous handler comes back."""

    def test_first_raises_the_others_are_recorded(self):
        before = signal.getsignal(signal.SIGINT)
        with aba.interrupts_once() as interrupts:
            with self.assertRaises(KeyboardInterrupt):
                signal.raise_signal(signal.SIGINT)
            signal.raise_signal(signal.SIGINT)
            signal.raise_signal(signal.SIGINT)
            self.assertEqual(interrupts, [signal.SIGINT] * 3)
        self.assertIs(signal.getsignal(signal.SIGINT), before)
        with aba.interrupts_once() as interrupts:
            interrupts.hold()
            signal.raise_signal(signal.SIGINT)
            self.assertEqual(interrupts, [signal.SIGINT])
        self.assertIs(signal.getsignal(signal.SIGINT), before)
        # the previous handler comes back after an exception too
        with self.assertRaises(KeyboardInterrupt), aba.interrupts_once():
            signal.raise_signal(signal.SIGINT)
        self.assertIs(signal.getsignal(signal.SIGINT), before)


class SourceHashes(unittest.TestCase):
    """``aba.source_hashes``: the hash of every orchestrator and driver source, refused when the
    tool's directory is not printable text or a source cannot be read."""

    def test_refusals(self):
        saved = aba.HERE
        self.addCleanup(setattr, aba, "HERE", saved)
        self.assertEqual(
            aba.source_hashes(),
            {
                n: hashlib.sha256((HERE / n).read_bytes()).hexdigest()
                for n in aba.SOURCE_FILES
            },
        )
        with tempfile.TemporaryDirectory() as tmp:
            for name in aba.SOURCE_FILES:
                shutil.copy(HERE / name, Path(tmp) / name)
            aba.HERE = Path(tmp) / "a\tb"
            with self.assertRaises(aba.Refused) as ctx:
                aba.source_hashes()
            self.assertIn(
                "a\\tb': a path with a character that is not", str(ctx.exception)
            )
            aba.HERE = Path(tmp)
            self.assertEqual(set(aba.source_hashes()), set(aba.SOURCE_FILES))
            (Path(tmp) / "aba_rows.hpp").unlink()
            with self.assertRaises(aba.Refused) as ctx:
                aba.source_hashes()
            self.assertIn("the tool's sources: ", str(ctx.exception))
            self.assertIn("aba_rows.hpp", str(ctx.exception))


# A fake harness: rows of two groups, the candidate twice as fast, special rows for every rule.
# Every call appends "side mode row argument pid" to DIR/events.log, and a block's environment
# and working directory as one JSON line to DIR/env.log (a block line carries only the keys of
# its rule); FAKE_STDERR names a side whose blocks write one line to stderr.
FAKE = r"""#!{py}
import json, os, shutil, signal, sys, time
SIDE = {side!r}
ROWS = ["n.fast", "n.digest", "n.drift", "n.bimodal", "n.rss", "n.short", "n.tamper",
        "n.cross_tamper", "n.end_tamper", "n.late_a1", "n.late_b", "n.late_a2", "n.garbled",
        "n.negative", "n.other", "n.warmfail", "n.crash", "n.fewer", "n.vanish", "n.noexec",
        "n.wander", "n.noisy", "n.tinyrss", "n.huge", "n.hugewarm", "n.longplan", "n.rsswarm",
        "n.intwarm", "n.sigint_warm", "n.sigint", "n.rmlogs", "n.lockout", "n.lockreport",
        "n.extra", "n.nan", "n.deep", "o.ok", "o.big", "o.warm"]
mode, row = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else None)
me = os.path.realpath(__file__)
top = os.path.dirname(os.path.dirname(me))
if mode == "aba-list":
    listing = [json.dumps({{"row": name, "group": name.split(".")[0]}}) for name in ROWS]
    broken = os.environ.get("FAKE_LIST")  # a broken list on request
    if broken == "junk":
        listing.insert(0, "a line that is no row record")
    if broken == "twice":
        listing.append(listing[0])
    if broken == "nogroup":
        listing.append(json.dumps({{"row": "n.extra"}}))
    if broken == "slash":
        listing.append(json.dumps({{"row": "n/extra", "group": "n"}}))
    if broken == "unknown" and SIDE == "new":  # a key no rule lists
        listing.append(json.dumps({{"row": "n.extra", "group": "n", "x": 1}}))
    if broken == "deepdesc" and SIDE == "new":  # a description nested 110,000 deep
        listing.append('{{"row": "n.extra", "group": "n", "desc": ' + "[" * 110000 + "]" * 110000 + "}}")
    if broken == "extra" and SIDE == "new":  # a row the reference does not list
        listing.append(json.dumps({{"row": "n.list_extra", "group": "n"}}))
    if broken == "fail":
        sys.exit(5)
    if broken == "sigint" and SIDE == "new":  # Ctrl-C before the first block
        os.kill(os.getppid(), signal.SIGINT)
        time.sleep(30)
    print("\n".join(listing))
    sys.exit(0)
with open(os.path.join(top, "events.log"), "a") as f:
    f.write(f"{{SIDE}} {{mode}} {{row}} {{sys.argv[3]}} {{os.getpid()}}\n")
with open(os.path.join(top, "env.log"), "a") as f:
    f.write(json.dumps({{
        "side": SIDE, "mode": mode, "row": row,
        "core": {{k: v for k, v in os.environ.items() if k.startswith("NGRAPH_CORE_")}},
        "pythonpath": os.environ.get("PYTHONPATH"), "pythonhome": os.environ.get("PYTHONHOME"),
        "malloc": os.environ.get("MallocLargeCache"), "cwd": os.getcwd()}}) + "\n")
if os.environ.get("FAKE_STDERR") == SIDE:  # one stderr line in every block of that side
    print(f"{{SIDE}} {{mode}} {{row}}: a message", file=sys.stderr)
# Ctrl-C during the candidate's warm-up block or its first timed block: SIGINT to aba.py, then
# a block that would last 30 s unless aba.py stops it.
if SIDE == "new" and (row, mode) in (("n.sigint_warm", "aba-calibrate"), ("n.sigint", "aba-run")):
    os.kill(os.getppid(), signal.SIGINT)
    time.sleep(30)
iters = 2000 if mode == "aba-calibrate" else int(sys.argv[3])
t = 0.0005 if SIDE == "new" else 0.001
digest = "00000000000000bb" if row == "n.digest" and SIDE == "new" else "00000000000000aa"
n = 0
if mode == "aba-run":
    counter = os.path.join(top, f"n-{{SIDE}}-{{row}}")
    n = int(open(counter).read()) if os.path.exists(counter) else 0
    open(counter, "w").write(str(n + 1))
    if row == "n.drift" and SIDE == "old":
        t = 0.001 * (1.10 if n % 2 else 1.0)  # A1 and A2 apart by 10 %
    if row == "n.bimodal" and SIDE == "new":
        t = 0.00060 if n in (1, 4) else 0.00050  # B split 5/2 at a 20 % gap
    if row == "n.wander" and SIDE == "old":
        # A drifts 2 % per triple and A2 is 4 % above A1: drift above 3 %, no cluster
        t = 0.001 * (1 + 0.02 * (n // 2)) * (1.04 if n % 2 else 1.0)
    if row == "n.noisy" and SIDE == "new":
        t = 0.0005 * (1 + 0.02 * n)  # a different ratio in every triple
    # A digest change only in triple 5 of 7: the reference runs A1 then A2 of triple k as its
    # calls 2k and 2k + 1, the candidate runs B of triple k as its call k.
    if (row, SIDE, n) in (("n.late_a1", "old", 8), ("n.late_a2", "old", 9), ("n.late_b", "new", 4)):
        digest = "00000000000000cc"
mib = 2**20
rss = 40 * mib  # above the fake process's own peak, so os.wait4 does not mask it
if row in ("n.rss", "o.big") and SIDE == "new":
    rss = 60 * mib
if row == "o.warm" and mode == "aba-calibrate" and SIDE == "old":
    rss = 500 * mib  # a warm-up peak: must not raise the reference
if row == "o.warm" and SIDE == "new" and mode == "aba-run":
    rss = 49 * mib  # over 40 + 8 MiB: fails unless the warm-up peak counted
if row == "n.tinyrss":
    rss = 1  # below the process's real peak: os.wait4's reading is the larger
rec = {{"row": row, "iterations": iters, "median_s": t,
       "mean_s": t, "min_s": t, "max_s": t, "setup_s": 0.0,
       "block_s": (0.5 if row == "n.short" and mode == "aba-run" else iters * t),
       "digest": digest, "rss_self_bytes": rss, "rss_children_bytes": 0}}
# Finite numbers whose derived values are not: medians near 1e308 in the timed blocks (the
# mean of A1 and A2 overflows), a warm-up median that makes the expected duration overflow or
# end beyond any date, and warm-up peaks whose sum overflows.
if row == "n.huge" and mode == "aba-run":
    rec.update(median_s=1.5e308 if SIDE == "new" else 1e308, block_s=2.5)
if row == "n.hugewarm" and mode == "aba-calibrate":
    rec.update(median_s=1e308, block_s=2.5)
if row == "n.longplan" and mode == "aba-calibrate":
    rec.update(median_s=1e15, block_s=2.5)
if row == "n.rsswarm" and mode == "aba-calibrate":
    rec.update(rss_self_bytes=1.7e308, rss_children_bytes=1.7e308)
if row == "n.intwarm" and mode == "aba-calibrate":
    rec.update(median_s=10**300, block_s=2.5)  # an integer: times a pinned count, above any float
if mode == "aba-run":
    # A binary changed or lost: by itself during a block, by the reference before the
    # candidate's first block or during the last block of the run; removed; made unexecutable.
    if row == "n.tamper" and SIDE == "new" or row == "n.cross_tamper" and SIDE == "old" and n == 0:
        with open(os.path.join(top, "new", "fake-new"), "a") as f:
            f.write("# rebuilt during the run\n")
    if row == "n.end_tamper" and SIDE == "old" and n == 13:
        with open(os.path.join(top, "new", "fake-new"), "a") as f:
            f.write("# rebuilt during the run\n")
    if row == "n.vanish" and SIDE == "new":
        os.remove(me)
    if row == "n.noexec" and SIDE == "new":
        os.chmod(me, 0o644)
    if row == "n.fewer" and SIDE == "new":
        rec["iterations"] = iters - 1
    # An operating-system error in aba.py after the first block: the logs' directory removed,
    # the output directory made read-only, or report.json (FAKE_OUT names the directory).
    if row == "n.rmlogs" and SIDE == "new":
        shutil.rmtree(os.path.join(os.environ["FAKE_OUT"], "blocks"))
    if row == "n.lockout" and SIDE == "new":
        os.chmod(os.environ["FAKE_OUT"], 0o555)
    if row == "n.lockreport" and SIDE == "new":
        os.chmod(os.path.join(os.environ["FAKE_OUT"], "report.json"), 0o444)
if (os.environ.get("FAKE_CHANGE_SIDE") == SIDE and
        os.environ.get("FAKE_CHANGE_MODE") == mode):
    rec.update(json.loads(os.environ["FAKE_CHANGE"]))
line = json.dumps(rec)
if SIDE == "new" and row == "n.warmfail" and mode == "aba-calibrate":
    print("not a block line")  # the warm-up block's response is broken
    sys.exit(0)
# A key no rule lists: a plain one in the candidate's timed block; NaN in the reference's
# warm-up block (written to report.json as nonstandard JSON when it was copied there); one
# nested 110,000 deep in the candidate's timed block (read by Python 3.14's JSON reader, which
# its writer cannot write back, so copying it crashed the run without outputs; older readers
# refuse it as nested too deeply).
if row == "n.extra" and SIDE == "new" and mode == "aba-run":
    line = line[:-1] + ', "x": 1}}'
if row == "n.nan" and SIDE == "old" and mode == "aba-calibrate":
    line = line[:-1] + ', "x": NaN}}'
if row == "n.deep" and SIDE == "new" and mode == "aba-run":
    line = line[:-1] + ', "x": ' + "[" * 110000 + "]" * 110000 + "}}"
if SIDE == "new" and mode == "aba-run":
    if row == "n.crash":
        sys.exit(3)  # a harness that fails
    if row == "n.garbled":
        line = '{{"row": "n.garbled", "median_s": 0.0005,'  # a truncated line
    if row == "n.negative":
        line = json.dumps({{**rec, "median_s": -t}})
    if row == "n.other":
        line = json.dumps({{**rec, "row": "n.fast"}})  # another row's record
print(line)
"""

# A fake ps, first on PATH: aba.py runs ps for the machine state at the start of a run and at its
# end. At its second call, while the final outputs are written, it sends SIGINT to aba.py TIMES
# times. Each call appends "x" to CALLS.
FAKE_PS = r"""#!{py}
import os, signal
with open({calls!r}, "a+") as f:
    f.seek(0)
    n = len(f.read())
    f.write("x")
if n == 1:
    for _ in range({times}):
        os.kill(os.getppid(), signal.SIGINT)
print(" 1.0 /usr/bin/fake")
"""

# Variables removed from the tests' own environment so a developer's shell cannot change a
# result (a test sets what it needs): the allocator variables aba.py records, PYTHONPATH and
# PYTHONHOME, which reach the Python of aba.py and of the fake harnesses, and NGRAPH_CORE_*, which
# the fake harness reports.
SPECIAL = ("PYTHONPATH", "PYTHONHOME", "NGRAPH_CORE_", *aba.ALLOCATOR_PREFIXES)


def clean_env(**extra: str) -> dict[str, str]:
    env = {
        k: v
        for k, v in os.environ.items()
        if not k.startswith(SPECIAL) and k not in aba.ALLOCATOR_NAMES
    }
    env.update(extra)
    return env


# n.intwarm's warm-up median (an integer near 1e300) times the largest count overflows a float.
INTWARM_PIN = ("--iters", f"n.intwarm={aba.MAX_COUNT}")


def strict_json(text: str):
    """Standard JSON only: Python's reader would also accept NaN and the infinities."""

    def refuse(literal: str):
        raise ValueError(f"{literal} is not standard JSON")

    return json.loads(text, parse_constant=refuse)


def report_of(out: Path) -> dict:
    """``out``'s report.json, read as standard JSON."""
    return strict_json((out / "report.json").read_text())


def by_row(rep: dict) -> dict[str, dict]:
    return {r["row"]: r for r in rep["rows"]}


def blocks_of(rec: dict) -> list[dict]:
    return [t[k] for t in rec["triples"] for k in ("A1", "B", "A2")]


class Driver(unittest.TestCase):
    """aba.py end to end with two fake harnesses (a few hundred short processes)."""

    QUICK = ("--resamples", "200")

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name).resolve()
        self.fakes()

    def fakes(self):
        """Both fake harnesses and their build records, written afresh."""
        for side in ("old", "new"):
            path = self.dir / side / f"fake-{side}"
            path.parent.mkdir(exist_ok=True)
            path.write_text(FAKE.format(py=sys.executable, side=side))
            path.chmod(0o755)
            # the build record build.sh writes beside a harness (the report's side info), built
            # from this checkout's driver
            path.with_name(path.name + ".json").write_text(
                json.dumps(fake_record(path, side))
            )

    def tearDown(self):
        self.tmp.cleanup()

    def aba(self, *args: str, env: dict[str, str] | None = None, script: Path = ABA):
        return subprocess.run(
            [sys.executable, "-B", str(script), *args],
            capture_output=True,
            text=True,
            env=clean_env(**(env or {})),
            timeout=600,
        )

    def sides(self) -> tuple[str, ...]:
        return (
            "--old-bin",
            str(self.dir / "old" / "fake-old"),
            "--new-bin",
            str(self.dir / "new" / "fake-new"),
        )

    def run_aba(self, out: str, *extra: str, env: dict[str, str] | None = None):
        sides = () if "--old-bin" in extra else self.sides()
        proc = self.aba(*sides, "--out", str(self.dir / out), *extra, env=env)
        path = self.dir / out / "report.json"
        return proc, (report_of(path.parent) if path.exists() else None)

    def summary(self, out: str) -> list[str]:
        return (self.dir / out / "summary.md").read_text().splitlines()

    def events(self, row: str) -> list[list[str]]:
        return [
            ln.split()
            for ln in (self.dir / "events.log").read_text().splitlines()
            if ln.split()[2] == row
        ]

    def envs(self, row: str) -> list[dict]:
        """The environment and working directory each block of ``row`` reported (env.log)."""
        lines = (self.dir / "env.log").read_text().splitlines()
        return [rec for rec in map(json.loads, lines) if rec["row"] == row]

    def blocks_run(self) -> int:
        log = self.dir / "events.log"
        return len(log.read_text().splitlines()) if log.exists() else 0

    def refused(self, out: str, *args: str, env: dict[str, str] | None = None) -> str:
        """Run aba.py; it must refuse before any block (exit 2, no traceback,
        no harness block run) and return its stderr. Both files carry the reason."""
        before = self.blocks_run()
        proc, rep = self.run_aba(out, *args, env=env)
        self.assertEqual(proc.returncode, 2, (args, proc.stderr[-500:]))
        self.assertNotIn("Traceback", proc.stderr)
        self.assertIsNotNone(rep, args)
        self.assertEqual((rep["verdict"], rep["conforming"]), ("aborted", False), args)
        self.assertIn(rep["aborted"], proc.stderr)
        self.assertIn(rep["aborted"], (self.dir / out / "summary.md").read_text())
        self.assertEqual(self.blocks_run(), before, args)
        shutil.rmtree(self.dir / out)
        return proc.stderr

    def aborted(self, out: str, *args: str, env: dict[str, str] | None = None) -> dict:
        """Run aba.py; it must abort (exit 1, no traceback, verdict aborted, not conforming,
        the reason in the summary) and return the report."""
        proc, rep = self.run_aba(out, *args, *self.QUICK, env=env)
        self.assertEqual(proc.returncode, 1, (args, proc.stderr[-500:]))
        self.assertNotIn("Traceback", proc.stderr)
        self.assertEqual(rep["verdict"], "aborted", args)
        self.assertIs(rep["conforming"], False, args)
        summary = self.summary(out)
        self.assertIn("overall verdict **aborted**", summary[3])
        self.assertIn(f"* **aborted**: {rep['aborted']}", summary)
        self.assertIn(f"aborted: {rep['aborted']}", proc.stderr)
        return rep

    def test_clean_run(self):
        proc, rep = self.run_aba("same", "--rows", " n.fast , o.ok ")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stderr, "")
        fast = by_row(rep)["n.fast"]
        self.assertEqual(rep["verdict"], "pass")
        self.assertIs(rep["conforming"], True)
        self.assertEqual(rep["selection"]["selected"], ["n.fast", "o.ok"])
        self.assertEqual(rep["protocol"]["bootstrap_resamples"], 10000)
        self.assertEqual(len(fast["triples"]), 7)
        # The candidate is twice as fast: per-side counts ceil(2 s * 1.10 / median),
        # statistic -50 %.
        self.assertEqual((fast["iterations_a"], fast["iterations_b"]), (2200, 4400))
        self.assertEqual(fast["iterations_source"], "calibrated per side")
        self.assertEqual(
            [b["iterations"] for b in blocks_of(fast)], [2200, 4400, 2200] * 7
        )
        self.assertAlmostEqual(fast["statistic"], -0.5)
        self.assertEqual(fast["digest_verdict"], "pass")
        for blk in [*fast["warmup"].values(), *blocks_of(fast)]:
            for key in ("load1_before", "load1_after"):
                self.assertTrue(math.isfinite(blk[key]) and blk[key] >= 0, blk)
        for when in ("state_at_start", "state_at_end"):
            load = rep["machine"][when]["load"]
            self.assertEqual(len(load), 3)
            self.assertTrue(all(math.isfinite(x) and x >= 0 for x in load))
        # meta.txt: both binaries and every source, as the report records them
        meta = (self.dir / "same" / "meta.txt").read_text().splitlines()
        want = {
            "A/binary/harness": hashlib.sha256(
                (self.dir / "old" / "fake-old").read_bytes()
            ).hexdigest(),
            "B/binary/harness": hashlib.sha256(
                (self.dir / "new" / "fake-new").read_bytes()
            ).hexdigest(),
            **{
                f"run/source/{n}": hashlib.sha256((HERE / n).read_bytes()).hexdigest()
                for n in aba.SOURCE_FILES
            },
        }
        self.assertEqual(rep["meta"], want)
        self.assertEqual(aba.read_meta(self.dir / "same" / "meta.txt"), want)
        self.assertIn("label\tA\told", meta)
        self.assertEqual(rep["sides"]["B"]["build_record"]["rev"], "new")
        summary = self.summary("same")
        # the summary's first lines: the conformance, then the verdict
        self.assertEqual(
            summary[2],
            "* **conforming run** (the protocol's evidence: "
            "accepted when the overall verdict is pass)",
        )
        self.assertIn("overall verdict **pass**", summary[3])
        self.assertIn("| n.fast | 2200 / 4400 |", "\n".join(summary))
        # blocks.csv: every block, warm-ups first
        csv_lines = (self.dir / "same" / "blocks.csv").read_text().splitlines()
        self.assertEqual(csv_lines[0].split(","), list(aba.CSV_FIELDS))
        self.assertEqual(len(csv_lines), 1 + 2 * (2 + 21))
        self.assertTrue(csv_lines[1].startswith("n.fast,W-A,,old,2000,,"))
        self.assertTrue(csv_lines[3].startswith("n.fast,A1,1,old,2200,2200,"))

    def test_block_order_counts_and_fresh_processes(self):
        proc, rep = self.run_aba("order", "--rows", "o.ok,n.fast", *self.QUICK)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        log = [ln.split() for ln in (self.dir / "events.log").read_text().splitlines()]
        # the warm-up blocks of every selected row first, one per side, given the target
        self.assertEqual(
            [e[:4] for e in log[:4]],
            [
                ["old", "aba-calibrate", "n.fast", "2.0"],
                ["new", "aba-calibrate", "n.fast", "2.0"],
                ["old", "aba-calibrate", "o.ok", "2.0"],
                ["new", "aba-calibrate", "o.ok", "2.0"],
            ],
        )
        # then row by row seven reference/candidate/reference repetitions with their counts
        want = [
            [side, "aba-run", row, count]
            for row in ("n.fast", "o.ok")
            for _ in range(7)
            for side, count in (("old", "2200"), ("new", "4400"), ("old", "2200"))
        ]
        self.assertEqual([e[:4] for e in log[4:]], want)
        pids = [e[4] for e in log]
        self.assertEqual(len(set(pids)), len(pids))  # a fresh process per block
        # each slot holds its side's answer (the reference's median 1 ms, the candidate's 0.5 ms)
        slots = [
            (s, t[s]["side"], t[s]["median_s"])
            for t in by_row(rep)["o.ok"]["triples"]
            for s in ("A1", "B", "A2")
        ]
        self.assertEqual(
            slots,
            [("A1", "old", 0.001), ("B", "new", 0.0005), ("A2", "old", 0.001)] * 7,
        )
        # the plan, printed before the first timed block, with each row's expected duration
        plan_at = proc.stdout.index("plan (from the warm-up blocks")
        self.assertLess(plan_at, proc.stdout.index("n.fast: pass"))
        self.assertEqual(rep["plan"]["rows"], 2)
        for rec in rep["rows"]:
            wa, wb = rec["warmup"]["A"], rec["warmup"]["B"]
            blk_a = 2200 * 0.001 + max(0.0, wa["wall_s"] - wa["block_s"])
            blk_b = 4400 * 0.0005 + max(0.0, wb["wall_s"] - wb["block_s"])
            self.assertAlmostEqual(rec["expected"]["row_s"], 7 * (2 * blk_a + blk_b))
        self.assertAlmostEqual(
            rep["plan"]["expected_s"], sum(r["expected"]["row_s"] for r in rep["rows"])
        )

    def test_options_take_effect(self):
        # counts calibrated with the margin; 0 accepted where the domain allows it
        proc, rep = self.run_aba(
            "zero",
            "--rows",
            "n.fast,n.drift,n.rss",
            "--calibration-margin",
            "0",
            "--tolerance",
            "0.2",
            "--rss-floor-mib",
            "0",
            "--budget",
            "1.6",
            "--seed",
            "-5",
            *self.QUICK,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        res = by_row(rep)
        self.assertEqual(
            (res["n.fast"]["iterations_a"], res["n.fast"]["iterations_b"]), (2000, 4000)
        )
        self.assertEqual(res["n.drift"]["time_verdict"], "pass")  # drift 9 % <= 20 %
        self.assertEqual(res["n.rss"]["rss_verdict"], "pass")  # 60 <= 1.6 x 40 MiB
        self.assertEqual(res["n.rss"]["rss_allowance_bytes"], 64 * MIB)
        # a tolerance of 0: equal blocks are no drift and no cluster
        proc, rep = self.run_aba(
            "tolerance0", "--rows", "n.fast", "--tolerance", "0", *self.QUICK
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(rep["protocol"]["tolerance"], 0.0)
        self.assertEqual(by_row(rep)["n.fast"]["time_verdict"], "pass")
        # the floor alone: 40 + 30 MiB covers 60 MiB
        proc, rep = self.run_aba(
            "floor", "--rows", "n.rss", "--rss-floor-mib", "30", *self.QUICK
        )
        self.assertEqual(by_row(rep)["n.rss"]["rss_allowance_bytes"], 70 * MIB)
        # --seed and --resamples reach the bootstrap; --triples and --target-seconds the blocks
        proc, rep = self.run_aba(
            "seed",
            "--rows",
            "n.noisy,n.short",
            "--seed",
            "7",
            "--resamples",
            "5",
            "--triples",
            "5",
            "--target-seconds",
            "0.4",
            "--allow-short",
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rec = by_row(rep)["n.noisy"]
        self.assertEqual(len(rec["triples"]), 5)
        x = [r - 1 for r in rec["ratios"]]
        self.assertEqual(rec["ci95"], list(aba.bootstrap_interval(x, 5, 7)))
        # (the interval depends on both: the default seed or count gives another one)
        self.assertNotEqual(rec["ci95"], list(aba.bootstrap_interval(x, 5, 20261003)))
        self.assertNotEqual(rec["ci95"], list(aba.bootstrap_interval(x, 10000, 7)))
        short = by_row(rep)["n.short"]
        self.assertEqual(short["block_verdict"], "pass")  # 0.5 s >= 0.4 s
        self.assertEqual(rec["iterations_a"], math.ceil(0.4 * 1.1 / 0.001))
        warmups = [e[:4] for e in self.events("n.noisy")[:2]]
        self.assertEqual(
            warmups, [[s, "aba-calibrate", "n.noisy", "0.4"] for s in ("old", "new")]
        )

    def test_iters_pin_takes_effect(self):
        proc, rep = self.run_aba(
            "pin", "--rows", "o.ok,n.fast", "--iters", "o.ok=3:5,n.fast=4", *self.QUICK
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        for name, want in (("o.ok", [3, 5, 3]), ("n.fast", [4, 4, 4])):
            rec = by_row(rep)[name]
            self.assertEqual(rec["iterations_source"], "pinned")
            self.assertEqual(
                (rec["iterations_a"], rec["iterations_b"]), tuple(want[:2])
            )
            self.assertEqual([b["iterations"] for b in blocks_of(rec)], want * 7)
            self.assertEqual(
                [e[3] for e in self.events(name)[2:]], [str(c) for c in want] * 7
            )

    def test_option_domains(self):
        for repeated in (
            ("--triples", "7", "--triples", "8"),
            ("--seed=1", "--seed=2"),
            ("--allow-short", "--allow-short"),
        ):
            self.assertIn(
                "given more than once",
                self.refused("repeated", "--rows", "o.ok", *repeated),
            )
        self.assertIn(
            "unrecognized arguments",
            self.refused("prefix", "--rows", "o.ok", "--target", "2"),
        )
        for option, value, why in (
            ("--triples", "0", "value 0 is not positive"),
            ("--triples", "2.5", "invalid parse value: '2.5'"),
            ("--tolerance", "nan", "value nan is not finite"),
            ("--tolerance", "-0.1", "value -0.1 is negative"),
            ("--budget", "0", "value 0.0 is not positive"),
            ("--rss-floor-mib", "-1", "value -1.0 is negative"),
            ("--target-seconds", "-2", "value -2.0 is not positive"),
            ("--target-seconds", "inf", "value inf is not finite"),
            ("--calibration-margin", "-0.5", "value -0.5 is negative"),
            ("--resamples", "0", "value 0 is not positive"),
            ("--seed", "x", "invalid parse value: 'x'"),
            ("--triples", "9" * 400, "value is not finite (an integer of 400 digits)"),
        ):
            stderr = self.refused("option", "--rows", "o.ok", option, value)
            self.assertIn(f"argument {option}: {why}", stderr)
        # the values derived from them must be finite
        for option, value, why in (
            (
                "--target-seconds",
                "1.7e308",
                "--target-seconds * (1 + --calibration-margin) is inf, not finite",
            ),
            ("--rss-floor-mib", "1e303", "--rss-floor-mib in bytes is inf, not finite"),
        ):
            self.assertIn(why, self.refused("derived", "--rows", "o.ok", option, value))
        # below the protocol only with --allow-short
        for args in (("--triples", "3"), ("--target-seconds", "1")):
            stderr = self.refused("short", "--rows", "o.ok", *args)
            self.assertIn("below the protocol", stderr)
            self.assertIn("--allow-short", stderr)

    def test_rows_glob_characters_are_accepted(self):
        # main's character check passes '?', '[', ']' and '!' on to the row selection
        globs = ["?.ok", "[on].ok", "[!n].ok"]
        proc, rep = self.run_aba("globs", "--rows", ",".join(globs), *self.QUICK)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(rep["selection"]["selected"], ["o.ok"])
        self.assertEqual(
            [(g["glob"], g["matched"]) for g in rep["selection"]["globs"]],
            [(g, ["o.ok"]) for g in globs],
        )

    def test_input_rules(self):
        o_ok = ("--rows", "o.ok", *self.QUICK)
        # --rows
        stderr = self.refused("nomatch", "--rows", "o.ok,placement.*")
        self.assertIn("'placement.*' matches no row or group", stderr)
        self.assertIn("o (3): o.ok, o.big, o.warm", stderr)
        self.assertIn("names no glob", self.refused("noglob", "--rows", ","))
        for rows in ("O.OK", "o.ok;x", "o/ok", "o.ok\t", "o.\u00e9"):
            why = self.refused("rowtext", "--rows", rows)
            self.assertIn(f"--rows {rows!r}: only the characters of row names", why)
        # --iters: every item row=N or row=NA:NB, naming a selected row, once
        for iters, why in (
            ("o.ok", "--iters item 'o.ok' is not row=N or row=NA:NB"),
            ("=5", "--iters item '=5' is not row=N or row=NA:NB"),
            ("o.big=5", "--iters row 'o.big' is not a selected row"),
            ("o.ok=5,o.ok=6", "--iters row 'o.ok' given twice"),
            ("o.ok=5:", "iteration count '5:' is not N or NA:NB"),
            ("o.ok=0:5", "iteration count '0:5' below 1"),
            (f"o.ok=5:{2**63}", f"iteration count '5:{2**63}' above {2**63 - 1}"),
        ):
            self.assertIn(why, self.refused("iters", *o_ok, "--iters", iters))
        # the build record: a JSON object whose binary_sha256 is the binary's (each refusal
        # names the side)
        record = self.dir / "new" / "fake-new.json"
        for text, why in (
            ("{}", "a build record without binary_sha256"),
            ("[]", "a build record without binary_sha256"),
            ("{not json", "not a readable JSON build record"),
            # deeper than Python's reader goes (it raises RecursionError)
            ('{"a":' * 10**6 + "1" + "}" * 10**6, "JSON nested too deeply"),
        ):
            record.write_text(text)
            self.assertIn(why, self.refused("record", *o_ok))
        record.write_text(json.dumps({"binary_sha256": "0" * 64}))
        self.assertIn("differs from its build record", self.refused("record", *o_ok))
        # and exactly build.sh's record: its keys only, all of them, each of its type
        full = fake_record(self.dir / "new" / "fake-new", "new")
        sdk, target = full["sdk"], full["deployment_target"]
        for name, rec, why in (
            ("unknown key", {**full, "x": 1}, "unknown key 'x'"),
            ("unknown keys", {**full, "x": 1, "y": 2}, "unknown keys 'x', 'y'"),
            ("no key", {k: v for k, v in full.items() if k != "flags"}, "no flags"),
            (
                "no driver",
                {k: v for k, v in full.items() if k != "driver_sha256"},
                "no driver_sha256",
            ),
            (
                "text null",
                {**full, "compiler": None},
                "compiler is None, of another type",
            ),
            ("text number", {**full, "rev": 7}, "rev is 7, of another type"),
            (
                "count text",
                {**full, "dirty_core_paths": "2"},
                "dirty_core_paths is '2', of another type",
            ),
            (
                "count bool",
                {**full, "dirty_core_paths": True},
                "dirty_core_paths is True, of another type",
            ),
            (
                "count negative",
                {**full, "dirty_core_paths": -1},
                "dirty_core_paths -1 is negative",
            ),
            (
                "count huge",
                {**full, "dirty_core_paths": 10**400},
                "dirty_core_paths is not finite (an integer of 401 digits)",
            ),
            (
                "not sha",
                {**full, "core_sources_sha256": "x"},
                "core_sources_sha256 'x' is not a sha256",
            ),
            (
                "driver null",
                {**full, "driver_sha256": None},
                "driver_sha256 is None, of another type",
            ),
            (
                "driver list",
                {**full, "driver_sha256": []},
                "driver_sha256 is an array, of another type",
            ),
            ("map text", {**full, "sdk": "x"}, "sdk is 'x', not an object"),
            (
                "map unknown",
                {**full, "sdk": {**sdk, "y": None}},
                "unknown key 'y' in sdk",
            ),
            (
                "map no key",
                {
                    **full,
                    "deployment_target": {
                        k: v for k, v in target.items() if k != "raw"
                    },
                },
                "no deployment_target.raw",
            ),
            (
                "map type",
                {**full, "sdk": {**sdk, "source": None}},
                "sdk.source is None, of another type",
            ),
            (
                "map nested",
                {**full, "sdk": {**sdk, "path": [["x"]]}},
                "sdk.path is an array, of another type",
            ),
        ):
            record.write_text(json.dumps(rec))
            why = f"new: {record}: not a build record as build.sh writes it ({why}"
            self.assertIn(why, self.refused("record", *o_ok), name)
        # a NaN or a deeply nested value where no key is listed, and NaN where a text is: refused,
        # never stored
        text = json.dumps(full)[:-1]
        for literal, why in (
            (', "x": NaN}', "unknown key 'x'"),
            (', "x": ' + "[" * 110000 + "]" * 110000 + "}", "unknown key 'x'"),
        ):
            record.write_text(text + literal)
            refusal = self.refused("record", *o_ok)
            self.assertTrue(
                why in refusal or "JSON nested too deeply" in refusal, refusal
            )
        record.write_text(json.dumps(full).replace('"rev": "new"', '"rev": NaN'))
        self.assertIn("rev is nan, of another type", self.refused("record", *o_ok))
        record.unlink()
        # a binary that is no file
        proc = self.aba(
            "--old-bin",
            str(self.dir / "old" / "fake-old"),
            "--new-bin",
            str(self.dir / "missing"),
            "--out",
            str(self.dir / "nobin"),
        )
        self.assertEqual(proc.returncode, 2, proc.stderr)
        self.assertIn(f"error: new: {self.dir / 'missing'}: no such file", proc.stderr)
        self.assertEqual(report_of(self.dir / "nobin")["verdict"], "aborted")
        # an unreadable binary: refused, not a traceback
        (self.dir / "new" / "fake-new").chmod(0o000)
        try:
            why = self.refused("unreadable", *o_ok)
        finally:
            (self.dir / "new" / "fake-new").chmod(0o755)
        self.assertIn(f"{self.dir / 'new' / 'fake-new'}: cannot read it", why)
        # a build record of other driver sources than this checkout's, on either side: one
        # file changed after the build (a stale harness, even when both sides agree), a file
        # missing or added, a hash of another type (no driver_sha256, or one that is not an
        # object, is no build record of build.sh's: below); a side without a build record is
        # not checked
        here = checkout_driver()
        rows = "dev/perf/aba_rows.hpp"
        for sides, driver, files in (
            (("new",), {**here, rows: "0" * 64}, rows),
            (("old",), {**here, rows: "0" * 64}, rows),
            (("old", "new"), {**here, rows: "0" * 64}, rows),
            (("new",), {k: v for k, v in here.items() if k != rows}, rows),
            (("new",), {**here, "dev/perf/extra.hpp": "0" * 64}, "dev/perf/extra.hpp"),
            (("new",), {**here, rows: [[here[rows]]]}, rows),
        ):
            self.fakes()
            for side in sides:
                record = self.dir / side / f"fake-{side}.json"
                rec = json.loads(record.read_text())
                rec["driver_sha256"] = driver
                record.write_text(json.dumps(rec))
            why = self.refused("drivers", *o_ok)
            first = self.dir / sides[0] / f"fake-{sides[0]}"
            self.assertIn(
                f"{sides[0]}: harness {first} was not built from this checkout's driver"
                f" sources (its build record differs in {files}); rebuild it with build.sh",
                why,
            )
        self.fakes()
        (self.dir / "new" / "fake-new.json").unlink()
        proc, rep = self.run_aba("no-record", *o_ok)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(rep["sides"]["B"]["build_record"], {})
        # a harness given through a symbolic link: its build record is the one beside the
        # binary it names (read and checked, here a stale driver), and a record beside the link
        # that is not that file is refused
        self.fakes()
        link = self.dir / "link" / "fake-new"
        link.parent.mkdir()
        link.symlink_to(self.dir / "new" / "fake-new")
        through_link = (*self.sides()[:3], str(link), *o_ok)
        proc, rep = self.run_aba("link-ok", *through_link)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(rep["sides"]["B"]["build_record"]["rev"], "new")
        record = self.dir / "new" / "fake-new.json"
        rec = json.loads(record.read_text())
        record.write_text(json.dumps({**rec, "driver_sha256": {}}))
        why = self.refused("link-stale", *through_link)
        self.assertIn(
            f"harness {self.dir / 'new' / 'fake-new'} was not built from", why
        )
        for target_record in (False, True):
            self.fakes()
            if not target_record:
                record.unlink()
            link.with_name("fake-new.json").write_text(json.dumps(rec))
            why = self.refused("link-record", *through_link)
            self.assertIn(
                f"{link}.json: a build record beside {link}, not beside the binary it names",
                why,
            )
        link.with_name("fake-new.json").unlink()
        self.fakes()
        # two different side labels, each 1 to 32 of a-z, 0-9, '_' and '-'
        self.assertIn(
            "two different names",
            self.refused("labels", *o_ok, "--old-label", "x", "--new-label", "x"),
        )
        for label in ("", "a/b", "a.b", "Ref", "x" * 33, "a\tb", "a\nb", "\u00e9"):
            why = self.refused("labels", *o_ok, "--new-label", label)
            self.assertIn(f"--new-label {label!r}: a label is 1 to 32 of", why)
        proc, rep = self.run_aba(
            "label-ok", *o_ok, "--old-label", "0-ref_1", "--new-label", "c" * 32
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        # paths of printable characters only (meta.txt records them on one line)
        why = self.refused("path", *o_ok, "--expect-meta", "a\nb")
        self.assertIn("--expect-meta 'a\\nb': a path with a character", why)
        proc = self.aba(*self.sides(), "--out", str(self.dir / "a\nb"), *o_ok)
        self.assertEqual(proc.returncode, 2, proc.stderr)
        self.assertIn(
            "a\\nb': a path with a character that is not printable", proc.stderr
        )
        self.assertFalse((self.dir / "a\nb").exists())
        tab = self.dir / "new" / "fake\tnew"
        shutil.copy(self.dir / "new" / "fake-new", tab)
        proc = self.aba(
            *self.sides()[:3], str(tab), "--out", str(self.dir / "tab"), *o_ok
        )
        self.assertEqual(proc.returncode, 2, proc.stderr)
        self.assertIn("a path with a character that is not printable", proc.stderr)
        self.assertEqual(report_of(self.dir / "tab")["verdict"], "aborted")
        # --out absent or empty: a file of an earlier run is never mixed in
        (self.dir / "used").mkdir()
        (self.dir / "used" / "old.stderr").write_text("an earlier run\n")
        proc, rep = self.run_aba("used", *o_ok)
        self.assertEqual(proc.returncode, 2, proc.stderr)
        self.assertIn("must be absent or an empty directory", proc.stderr)
        self.assertEqual(
            [f.name for f in (self.dir / "used").iterdir()], ["old.stderr"]
        )
        (self.dir / "file").write_text("")
        proc, rep = self.run_aba("file", *o_ok)
        self.assertEqual(proc.returncode, 2, proc.stderr)
        self.assertIn("must be absent or an empty directory", proc.stderr)
        (self.dir / "empty").mkdir()
        proc, rep = self.run_aba("empty", *o_ok)
        self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_expect_meta(self):
        o_ok = ("--rows", "o.ok", *self.QUICK)
        proc, prev = self.run_aba("prev", *o_ok)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        meta = (self.dir / "prev" / "meta.txt").read_text().splitlines()
        # a meta.txt as a run writes it is accepted, and the summary says so
        path = str(self.dir / "prev" / "meta.txt")
        proc, rep = self.run_aba("pinned", *o_ok, "--expect-meta", path)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(rep["expect_meta"], path)
        self.assertIn(
            f"every hash equal to `{path}`", "\n".join(self.summary("pinned"))
        )
        # lines as meta.txt writes them, no hash twice, at least one
        pin = next(ln for ln in meta if ln.startswith("B\tbinary\t"))
        cols = pin.split("\t")
        changed = "\t".join([*cols[:3], "0" * 64, cols[4]])
        for name, lines, why in (
            ("twice", [changed, *meta], "B/binary/harness named twice"),
            ("junk", ["junk line", *meta], "line 1 is not a meta.txt line"),
            (
                "short",
                [ln.rsplit("\t", 1)[0] if ln == pin else ln for ln in meta],
                "is not a meta.txt line",
            ),
            (
                "nosha",
                [ln.replace(cols[3], "x" * 64) for ln in meta],
                "is not a meta.txt line",
            ),
            ("empty", ["# no hash", "label\tA\told"], "no hash"),
            (
                "role",
                [pin.replace("binary", "other"), *meta[1:]],
                "not a meta.txt line",
            ),
            ("label", ["label\tA\tBAD", *meta], "not a meta.txt line"),
            (
                "label-twice",
                ["label\tA\told", "label\tA\told", pin],
                "not a meta.txt line",
            ),
            ("path", [pin + "\x01"], "not a meta.txt line"),
            # every hash must equal this run's: one changed, one missing
            (
                "changed",
                [changed if ln == pin else ln for ln in meta],
                "hashes differ: B/binary/harness",
            ),
            (
                "missing",
                [ln for ln in meta if ln != pin],
                "hashes differ: B/binary/harness",
            ),
        ):
            path = self.dir / f"meta-{name}.txt"
            path.write_text("\n".join(lines) + "\n")
            self.assertIn(
                why, self.refused("meta", *o_ok, "--expect-meta", str(path)), name
            )
        self.assertIn(
            "No such file",
            self.refused("meta", *o_ok, "--expect-meta", str(self.dir / "none")),
        )

    def test_row_list_response_is_checked_whole(self):
        for broken, why in (
            ("junk", "malformed output line (Expecting value"),
            ("twice", "row 'n.fast' listed twice"),
            ("nogroup", "a row record without a string group"),
            ("fail", "aba-list: exit 5"),
            (
                "unknown",
                f"error: new: {self.dir / 'new' / 'fake-new'} aba-list: row n.extra: unknown key 'x'",
            ),
        ):
            stderr = self.refused("list", "--rows", "o.ok", env={"FAKE_LIST": broken})
            self.assertIn(why, stderr)
        # a nested value where a text is read: refused by its type, never printed (Python 3.14's
        # reader reads this nesting; older ones refuse it as too deep)
        stderr = self.refused("list", "--rows", "o.ok", env={"FAKE_LIST": "deepdesc"})
        why = 'a row record without a string desc: {"row": "n.extra", "group": "n", "desc": "<an array>"}'
        self.assertTrue(why in stderr or "JSON nested too deeply" in stderr, stderr)

    def test_child_environment(self):
        # the caller's environment reaches every harness unchanged, in the repository root
        passed = {
            "NGRAPH_CORE_SENSITIVITY_THREADS": "8",
            "NGRAPH_CORE_FOO": "x",
            "PYTHONPATH": str(self.dir),
            "PYTHONHOME": sys.base_prefix,
        }
        proc, rep = self.run_aba("env", "--rows", "o.ok", *self.QUICK, env=passed)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        blocks = self.envs("o.ok")
        self.assertEqual(len(blocks), 2 + 3 * 7)  # the warm-up pair and seven triples
        for blk in blocks:
            self.assertEqual(
                blk["core"],
                {"NGRAPH_CORE_SENSITIVITY_THREADS": "8", "NGRAPH_CORE_FOO": "x"},
            )
            self.assertEqual(blk["pythonpath"], str(self.dir))
            self.assertEqual(blk["pythonhome"], sys.base_prefix)
            self.assertEqual(blk["cwd"], str(REPO))

    def test_allocator_variable_recorded_not_changed(self):
        proc, rep = self.run_aba(
            "malloc",
            "--rows",
            "o.ok",
            *self.QUICK,
            env={"MallocLargeCache": "0", "LD_PRELOAD": "x.so"},
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        alloc = rep["allocator_env"]
        self.assertEqual(alloc["MallocLargeCache"], {"set": True, "value": "0"})
        self.assertEqual(
            alloc["variables"], {"LD_PRELOAD": "x.so", "MallocLargeCache": "0"}
        )
        self.assertIs(alloc["normal"], False)
        self.assertIs(rep["conforming"], False)
        malloc = [blk["malloc"] for blk in self.envs("o.ok")]
        self.assertEqual(malloc, ["0"] * (2 + 3 * 7))
        summary = self.summary("malloc")
        self.assertTrue(summary[2].startswith("* **not gate evidence**"))
        self.assertIn(
            "* allocator environment: MallocLargeCache=0, LD_PRELOAD=x.so (**allocator override:"
            " not gate evidence**)",
            summary,
        )
        proc, rep = self.run_aba("normal", "--rows", "o.ok", *self.QUICK)
        self.assertEqual(rep["allocator_env"]["variables"], {})
        # every name the rule lists, by prefix or in full, beside near misses that are none
        # (the record itself, not a run: set, they would load into or retune the allocator of
        # aba.py and of every process it starts)
        named = {
            "MallocScribble": "1",
            "MALLOC_ARENA_MAX": "2",
            "DYLD_INSERT_LIBRARIES": "x.dylib",
            "LD_PRELOAD": "x.so",
            "GLIBC_TUNABLES": "glibc.malloc.tcache_count=0",
        }
        misses = {"MALLOCX": "1", "malloc_x": "1", "XMalloc": "1", "LD_PRELOAD_X": "1"}
        rec = aba.allocator_record({**named, **misses, "PATH": "/bin"})
        self.assertEqual(rec["variables"], dict(sorted(named.items())))
        self.assertEqual(
            (rec["MallocLargeCache"], rec["normal"]),
            ({"set": False, "value": None}, False),
        )
        self.assertIn(
            "* allocator environment: MallocLargeCache unset, no other allocator variable (normal)",
            self.summary("normal"),
        )

    def test_machine_state_keeps_only_readable_ps_lines(self):
        # ps output is a diagnostic: a line that is not a finite share and a name is left out
        fake = self.dir / "fakebin"
        fake.mkdir()
        (fake / "ps").write_text(
            "#!/bin/sh\nprintf 'nan evil\\n12.5 /usr/bin/real\\ninf x\\n7 \\nbusy\\n'\n"
        )
        (fake / "ps").chmod(0o755)
        path = f"{fake}{os.pathsep}{os.environ.get('PATH', '')}"
        proc, rep = self.run_aba(
            "ps", "--rows", "o.ok", *self.QUICK, env={"PATH": path}
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        for when in ("state_at_start", "state_at_end"):
            self.assertEqual(rep["machine"][when]["top_cpu"], [[12.5, "real"]])
        lines = [ln for ln in self.summary("ps") if "top CPU real 12%" in ln]
        self.assertEqual(
            [ln.split(":")[0] for ln in lines],
            ["* machine state at start", "* machine state at end"],
        )

    def test_both_sides_list_the_same_rows(self):
        for broken, why in (
            ("extra", "only the candidate ['n.list_extra']"),
            ("slash", "row 'n/extra' is not 1 to 64 of the characters"),
        ):
            stderr = self.refused("rows", "--rows", "o.ok", env={"FAKE_LIST": broken})
            self.assertIn(why, stderr)

    def test_os_errors(self):
        # before the first block: a refusal that leaves nothing in --out (here an empty,
        # read-only directory, so meta.txt cannot be written)
        locked = self.dir / "locked"
        locked.mkdir()
        locked.chmod(0o555)
        try:
            proc, rep = self.run_aba("locked", "--rows", "o.ok")
            sub_proc, sub_rep = self.run_aba("locked/sub", "--rows", "o.ok")
            sub = sub_proc.stderr
            self.assertEqual(sub_proc.returncode, 2)
            self.assertIsNone(sub_rep)
            self.assertIn("cannot write refusal outputs", sub)
        finally:
            locked.chmod(0o755)
        self.assertEqual(proc.returncode, 2, proc.stderr)
        self.assertIn(f"--out {locked}: [Errno 13] Permission denied", proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)
        self.assertEqual(list(locked.iterdir()), [])
        self.assertIn(f"--out {locked / 'sub'}: [Errno 13] Permission denied", sub)
        self.assertEqual(self.blocks_run(), 0)
        # after it: the run aborts, with both outputs
        out = self.dir / "rmlogs"
        rep = self.aborted(
            "rmlogs", "--rows", "n.rmlogs,o.ok", env={"FAKE_OUT": str(out)}
        )
        self.assertIn("operating-system error: [Errno 2] No such file", rep["aborted"])
        self.assertEqual([r["status"] for r in rep["rows"]], ["timing", "planned"])
        self.assertTrue((out / "blocks.csv").is_file())
        # and when the outputs themselves cannot be written: exit 1, the reason on stderr
        out = self.dir / "lockout"
        try:
            proc, rep = self.run_aba(
                "lockout",
                "--rows",
                "n.lockout",
                "--triples",
                "2",
                "--allow-short",
                *self.QUICK,
                env={"FAKE_OUT": str(out)},
            )
        finally:
            out.chmod(0o755)
        # Both files existed before the first block. They can still be written directly when
        # the directory stops allowing a replacement; the final reason records that failure.
        self.assertEqual(proc.returncode, 1, proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)
        self.assertIn("Permission denied", proc.stderr)
        self.assertEqual(rep["verdict"], "aborted")
        self.assertIn("Permission denied", rep["aborted"])
        self.assertIn(rep["aborted"], (out / "summary.md").read_text())

    def test_default_output_directory(self):
        # Without --out the report lands in results/<second stamp>/ beside aba.py. The run uses a
        # copy of the tool's sources in this test's own directory, so it writes nothing into the
        # checkout and two runs at once never share a results directory.
        perf = self.dir / "checkout" / "dev" / "perf"
        perf.mkdir(parents=True)
        for name in aba.SOURCE_FILES:
            shutil.copy(HERE / name, perf / name)
        proc = self.aba(
            *self.sides(),
            "--rows",
            "o.ok",
            "--triples",
            "2",
            "--allow-short",
            *self.QUICK,
            script=perf / "aba.py",
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        report = Path(re.search(r"^report: (.*)$", proc.stdout, re.M).group(1))
        self.assertEqual(report.parent.parent, perf / "results")
        self.assertRegex(report.parent.name, r"^\d{8}T\d{6}$")
        self.assertTrue(report.is_file())

    def test_peak_rss_reading(self):
        # the larger of the harness's own reading (self + children) and os.wait4's, in bytes
        proc, rep = self.run_aba("rss", "--rows", "n.tinyrss,o.ok", *self.QUICK)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        for name in ("n.tinyrss", "o.ok"):
            rec = by_row(rep)[name]
            for blk in [*rec["warmup"].values(), *blocks_of(rec)]:
                own = blk["rss_self_bytes"] + blk["rss_children_bytes"]
                self.assertEqual(
                    blk["rss_tree_bytes"], max(own, blk["rss_wait4_bytes"])
                )
                # a Python process's real peak: some MiB, never GiB (a KiB/byte mix-up)
                self.assertTrue(MIB < blk["rss_wait4_bytes"] < 40 * MIB, blk)
                self.assertEqual(
                    blk["rss_tree_bytes"],
                    blk["rss_wait4_bytes"] if own == 1 else 40 * MIB,
                )

    def test_digest_change_in_a_later_block_fails(self):
        rows = "n.late_a1,n.late_b,n.late_a2,o.ok"
        proc, rep = self.run_aba("late", "--rows", rows, *self.QUICK)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        res = by_row(rep)
        for name, slot in (("n.late_a1", "A1"), ("n.late_b", "B"), ("n.late_a2", "A2")):
            rec = res[name]
            # the change sits in triple 5 only: the first triples agree
            self.assertEqual(
                rec["triples"][0]["A1"]["digest"], rec["triples"][0]["B"]["digest"]
            )
            self.assertEqual(
                rec["triples"][4][slot]["digest"], "00000000000000cc", name
            )
            self.assertIs(rec["digest_match"], False, name)
            self.assertEqual(rec["digest_verdict"], "fail", name)
            self.assertEqual(rec["verdict"], "fail", name)
        self.assertEqual(res["o.ok"]["digest_verdict"], "pass")

    def test_rules_apply_to_every_row(self):
        rowset = "n.digest,n.drift,n.bimodal,n.rss,n.short,n.wander,o.big,o.warm"
        proc, rep = self.run_aba("rules", "--rows", rowset, *self.QUICK)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        res = by_row(rep)
        self.assertEqual(res["n.digest"]["verdict"], "fail")
        self.assertEqual(res["n.digest"]["digest_verdict"], "fail")
        self.assertEqual(res["n.digest"]["time_verdict"], "pass")
        self.assertEqual(res["n.drift"]["time_verdict"], "inconclusive")
        self.assertGreater(res["n.drift"]["drift"], 0.03)
        self.assertAlmostEqual(res["n.drift"]["statistic"], 0.0005 / 0.00105 - 1)
        self.assertEqual(
            res["n.drift"]["bimodal"], ["A"]
        )  # A1 and A2 apart: two clusters
        wander = res["n.wander"]
        self.assertGreater(wander["drift"], 0.03)
        self.assertEqual(
            [wander[k] for k in ("time_verdict", "cluster_verdict", "verdict")],
            ["inconclusive", "pass", "inconclusive"],
        )
        self.assertEqual(res["n.bimodal"]["bimodal"], ["B"])
        self.assertEqual(res["n.bimodal"]["verdict"], "fail")
        self.assertEqual(res["n.rss"]["rss_verdict"], "fail")
        self.assertEqual(res["o.big"]["verdict"], "fail")
        warm = res["o.warm"]
        self.assertEqual(warm["warmup"]["A"]["rss_tree_bytes"], 500 * MIB)
        self.assertEqual(warm["rss_a1_max_bytes"], 40 * MIB)
        self.assertEqual(warm["rss_verdict"], "fail")
        self.assertEqual(res["n.short"]["block_verdict"], "fail")
        self.assertEqual(len(res["n.short"]["short_blocks"]), 21)
        self.assertEqual(rep["verdict"], "fail")
        # every other check passed: the short blocks alone make the run non-conforming
        self.assertIs(rep["allocator_env"]["normal"], True)
        self.assertIs(rep["conforming"], False)
        summary = self.summary("rules")
        header = next(ln for ln in summary if ln.startswith("| Row |"))
        table = [ln for ln in summary if ln.startswith(("| n.", "| o."))]
        self.assertEqual(len(table), 8)
        self.assertTrue(all(ln.count("|") == header.count("|") for ln in table))
        self.assertIn(
            "Rows measured: 8 (0 pass, 7 fail, 1 inconclusive); not measured: 0.",
            summary,
        )
        self.assertIn("Inconclusive rows: n.wander.", summary)
        self.assertIn("Digest mismatches (fail): n.digest.", summary)
        self.assertIn("Memory failures: n.rss, o.big, o.warm.", summary)
        # an inconclusive row makes a run without failing rows inconclusive
        proc, rep = self.run_aba("wander", "--rows", "n.wander,o.ok", *self.QUICK)
        self.assertEqual(rep["verdict"], "inconclusive")

    def test_failing_or_short_run_exits_0_with_its_verdict(self):
        # aba.py's exit status only says the run completed; the summary's first lines and
        # report.json give the verdict and the conformance.
        proc, rep = self.run_aba("fail", "--rows", "n.bimodal,o.ok", *self.QUICK)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(rep["verdict"], "fail")
        self.assertIs(rep["conforming"], True)
        self.assertIn("overall verdict **fail**", self.summary("fail")[3])
        proc, rep = self.run_aba(
            "short", "--rows", "o.ok", "--triples", "2", "--allow-short", *self.QUICK
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(rep["verdict"], "pass")
        self.assertEqual(by_row(rep)["o.ok"]["short_blocks"], [])
        self.assertIs(rep["conforming"], False)  # --allow-short alone
        summary = self.summary("short")
        self.assertTrue(summary[2].startswith("* **not gate evidence**"), summary[2])
        self.assertIn("--allow-short given", "\n".join(summary))

    def test_python_version(self):
        out = self.dir / "py39"
        code = (
            "import runpy, sys; sys.version_info = (3, 9, 6, 'final', 0); "
            f"sys.argv = ['aba.py', *{list(self.sides())!r}, '--out', {str(out)!r}]; "
            f"runpy.run_path({str(ABA)!r}, run_name='__main__')"
        )
        proc = subprocess.run(
            [sys.executable, "-B", "-c", code],
            capture_output=True,
            text=True,
            env=clean_env(),
        )
        self.assertEqual(proc.returncode, 2, proc.stderr)
        self.assertIn("error: aba.py needs Python 3.10 or newer", proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)
        self.assertEqual(report_of(out)["verdict"], "aborted")
        self.assertIn(report_of(out)["aborted"], (out / "summary.md").read_text())

    def test_stderr_logs(self):
        # Every block of a row, side and mode appends its stderr to one log; a log still empty
        # after a block is removed, and an abort names a log only when the failing block wrote
        # to it.
        proc, rep = self.run_aba("said", "--rows", "o.ok", env={"FAKE_STDERR": "new"})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        blocks = self.dir / "said" / "blocks"
        self.assertEqual(
            sorted(os.listdir(blocks)),
            ["o.ok.new.aba-calibrate.stderr", "o.ok.new.aba-run.stderr"],
        )
        self.assertEqual(
            (blocks / "o.ok.new.aba-calibrate.stderr").read_text(),
            "new aba-calibrate o.ok: a message\n",
        )
        self.assertEqual(
            (blocks / "o.ok.new.aba-run.stderr").read_text(),
            "new aba-run o.ok: a message\n" * 7,
        )
        for row in ("n.crash", "n.garbled"):
            rep = self.aborted(f"said-{row}", "--rows", row, env={"FAKE_STDERR": "new"})
            blocks = self.dir / f"said-{row}" / "blocks"
            log = blocks / f"{row}.new.aba-run.stderr"
            self.assertIn(f"new aba-run {row} ", rep["aborted"])
            self.assertIn(f" (stderr: {log})", rep["aborted"])
            self.assertEqual(log.read_text(), f"new aba-run {row}: a message\n")
            # the reference's A1 block before it wrote nothing: no log
            self.assertEqual(
                sorted(os.listdir(blocks)),
                [f"{row}.new.aba-calibrate.stderr", log.name],
            )
            rep = self.aborted(f"quiet-{row}", "--rows", row)
            self.assertIn(f"new aba-run {row} ", rep["aborted"])
            self.assertIn(" (nothing on stderr)", rep["aborted"])
            self.assertNotIn("stderr: ", rep["aborted"])
            self.assertEqual(os.listdir(self.dir / f"quiet-{row}" / "blocks"), [])

    def test_nonfinite_derived_values_abort_before_any_verdict(self):
        for row, why in (
            # A1 = A2 = 1e308, B = 1.5e308: the mean of A1 and A2 overflows; a ratio
            # against it would read 0, statistic -1, a false pass of a 50 % slowdown
            ("n.huge", "the mean of A1 and A2 is inf, not finite"),
            ("n.hugewarm", "n.hugewarm expected duration is inf, not finite"),
            ("n.longplan", "ends beyond any date"),
            ("n.rsswarm", "aba-calibrate n.rsswarm rss_tree_bytes is inf, not finite"),
        ):
            rep = self.aborted(row, "--rows", row)
            self.assertIn(why, rep["aborted"], row)
            rec = by_row(rep)[row]
            for key in ("verdict", "statistic", "time_verdict", "rss_verdict"):
                self.assertNotIn(key, rec, (row, key))
            # nothing written is a non-finite number
            text = (self.dir / row / "report.json").read_text()
            self.assertNotIn("Infinity", text, row)
            self.assertNotIn("NaN", text, row)
        # a finite target over a short warm-up median: the calibrated count overflows
        rep = self.aborted("count", "--rows", "o.ok", "--target-seconds", "1e308")
        self.assertIn(
            "o.ok old: the calibrated iteration count is inf, not finite",
            rep["aborted"],
        )
        self.assertNotIn("iterations_a", by_row(rep)["o.ok"])
        # a warm-up median printed as an integer, times a pinned count, is beyond any float:
        # an infinity the run reports, not an exception
        rep = self.aborted("intwarm", "--rows", "n.intwarm", *INTWARM_PIN)
        self.assertIn("n.intwarm expected duration is inf, not finite", rep["aborted"])

    def test_both_block_modes_and_sides_refuse_unruled_values(self):
        for side in ("old", "new"):
            for mode in ("aba-calibrate", "aba-run"):
                for change in (
                    {"x": 1},
                    {"x": float("nan")},
                    {"x": [[1]]},
                    {"median_s": [1]},
                    {"iterations": True},
                    {"p90_s": float("inf")},
                    {"iterations": 2**63},
                ):
                    with self.subTest(side=side, mode=mode, change=change):
                        out = f"schema-{side}-{mode}-{len(list(self.dir.iterdir()))}"
                        rep = self.aborted(
                            out,
                            "--rows",
                            "o.ok",
                            "--triples",
                            "1",
                            "--allow-short",
                            env={
                                "FAKE_CHANGE_SIDE": side,
                                "FAKE_CHANGE_MODE": mode,
                                "FAKE_CHANGE": json.dumps(change),
                            },
                        )
                        self.assertIn(f"{side} {mode} o.ok", rep["aborted"])
                        self.assertNotIn(
                            '"x":', (self.dir / out / "report.json").read_text()
                        )

    def test_final_output_exceptions_keep_both_files(self):
        # A serializer or writer can fail after a successful comparison, or after a failed
        # harness. Each output is still attempted; a failed serializer gets a minimal record.
        saved_report = aba.report_text
        saved_atomic = aba.atomic_text
        for row in ("o.ok", "n.crash"):
            for fault in (
                "json",
                "summary",
                "csv",
                "write",
                "machine",
                "closed_stderr",
                "closed_traceback",
            ):
                with self.subTest(row=row, fault=fault):
                    out = self.dir / f"output-{row}-{fault}"
                    argv = [
                        "aba.py",
                        *self.sides(),
                        "--rows",
                        row,
                        "--out",
                        str(out),
                        "--triples",
                        "1",
                        "--allow-short",
                        *self.QUICK,
                    ]

                    def bad_json(report):
                        if "verdict" in report and "protocol" in report:
                            raise RecursionError("full JSON failed")
                        return saved_report(report)

                    attempts = []

                    def bad_write(path, text, attempts=attempts):
                        # The first two writes make the initial files. Fail the final JSON
                        # replacement once; the retry must overwrite both files as aborted.
                        attempts.append(path.name)
                        if len(attempts) == 3:
                            raise RuntimeError("write failed")
                        return saved_atomic(path, text)

                    state_calls = []
                    saved_state = aba.machine_state

                    def bad_state(
                        state_calls=state_calls, saved_state=saved_state, fault=fault
                    ):
                        state_calls.append(1)
                        if len(state_calls) > 1:
                            if fault == "closed_stderr":
                                sys.stderr.close()
                            raise ValueError("machine diagnostic failed")
                        return saved_state()

                    def bad_csv(report):
                        sys.stderr.close()
                        raise RuntimeError("csv failed")

                    patch = {
                        "json": mock.patch.object(aba, "report_text", bad_json),
                        "summary": mock.patch.object(
                            aba,
                            "summary_text",
                            side_effect=ValueError("summary failed"),
                        ),
                        "csv": mock.patch.object(
                            aba, "csv_text", side_effect=RuntimeError("csv failed")
                        ),
                        "closed_traceback": mock.patch.object(aba, "csv_text", bad_csv),
                        "write": mock.patch.object(aba, "atomic_text", bad_write),
                        "machine": mock.patch.object(aba, "machine_state", bad_state),
                        "closed_stderr": mock.patch.object(
                            aba, "machine_state", bad_state
                        ),
                    }[fault]
                    with (
                        mock.patch.object(sys, "argv", argv),
                        patch,
                        contextlib.redirect_stdout(io.StringIO()),
                        contextlib.redirect_stderr(io.StringIO()),
                    ):
                        status = aba.main()
                    rep = report_of(out)
                    self.assertEqual(status, 1)
                    self.assertEqual(
                        (rep["verdict"], rep["conforming"]), ("aborted", False)
                    )
                    self.assertIn(rep["aborted"], (out / "summary.md").read_text())
                    self.assertEqual([r["row"] for r in rep["rows"]], [row])
                    self.assertFalse(list(out.glob(".aba-*")))

    def test_system_numbers_are_checked_before_storage(self):
        for value in (float("nan"), float("inf"), -1, True, [1], "1"):
            with self.subTest(value=value):
                with mock.patch.object(
                    os, "getloadavg", return_value=(value, 0.0, 0.0)
                ):
                    out = self.dir / f"system-{len(list(self.dir.iterdir()))}"
                    argv = [
                        "aba.py",
                        *self.sides(),
                        "--out",
                        str(out),
                        "--rows",
                        "o.ok",
                    ]
                    with (
                        mock.patch.object(sys, "argv", argv),
                        contextlib.redirect_stderr(io.StringIO()),
                    ):
                        self.assertEqual(aba.main(), 2)
                    rep = report_of(out)
                    self.assertIn("machine load", rep["aborted"])
                    self.assertIn(rep["aborted"], (out / "summary.md").read_text())

        # An interrupt before the final elapsed-time reading must check that reading too.
        out = self.dir / "interrupted-clock"
        out.mkdir()
        report = {"started_utc": "test", "rows": []}
        with (
            mock.patch.object(time, "perf_counter", return_value=float("nan")),
            contextlib.redirect_stderr(io.StringIO()),
        ):
            self.assertEqual(aba.finish_interrupted(report, out, 0.0), 1)
        self.assertNotIn("wall_s", report)
        saved = report_of(out)
        self.assertIn(saved["aborted"], (out / "summary.md").read_text())

    def test_interrupt_aborts_and_stops_the_block(self):
        # Ctrl-C (SIGINT to aba.py) during a warm-up or a timed block: the run aborts like any
        # other, and the harness of that block is stopped at once, not left running
        for row, block in (
            ("n.sigint_warm", "new aba-calibrate n.sigint_warm"),
            ("n.sigint", "new aba-run n.sigint"),
        ):
            started = time.monotonic()
            rep = self.aborted(row, "--rows", row)
            self.assertEqual(
                rep["aborted"], f"interrupted (SIGINT, Ctrl-C) during {block}"
            )
            self.assertLess(time.monotonic() - started, 20, row)
            pid = int(self.events(row)[-1][4])
            with self.assertRaises(ProcessLookupError, msg=row):
                os.kill(pid, 0)
        # Before the first block: exit 130, both outputs record the interrupted listing.
        started = time.monotonic()
        proc, rep = self.run_aba(
            "listing", "--rows", "o.ok", env={"FAKE_LIST": "sigint"}
        )
        self.assertEqual(proc.returncode, 130, proc.stderr)
        self.assertIn(
            "interrupted (SIGINT, Ctrl-C) during new aba-list", rep["aborted"]
        )
        self.assertIn(rep["aborted"], (self.dir / "listing" / "summary.md").read_text())
        self.assertLess(time.monotonic() - started, 20)

    def test_interrupt_while_the_outputs_are_written(self):
        # Ctrl-C while the final outputs are written (once, or twice): held until they are
        # written, then the run aborts like any other, its outputs written again
        for times in (1, 2):
            tools = self.dir / f"tools{times}"
            tools.mkdir()
            ps = tools / "ps"
            calls = tools / "calls"
            ps.write_text(
                FAKE_PS.format(py=sys.executable, calls=str(calls), times=times)
            )
            ps.chmod(0o755)
            path = f"{tools}{os.pathsep}{os.environ['PATH']}"
            rep = self.aborted(f"ps{times}", "--rows", "o.ok", env={"PATH": path})
            self.assertEqual(
                rep["aborted"],
                "interrupted (SIGINT, Ctrl-C) while the outputs were written",
            )
            self.assertEqual(by_row(rep)["o.ok"]["status"], "measured")
            # ps ran twice: at the start and at the end
            self.assertEqual(calls.read_text(), "xx", times)
            self.assertEqual(rep["machine"]["state_at_end"]["top_cpu"], [[1.0, "fake"]])
            with (self.dir / f"ps{times}" / "blocks.csv").open() as f:
                self.assertEqual(len(list(csv.DictReader(f))), 2 + 3 * 7, times)
        # An interrupt raised before the hold began (summary_text raises it once): the same
        # abort; a run that aborted already keeps its reason.
        saved = aba.summary_text
        self.addCleanup(setattr, aba, "summary_text", saved)
        self.addCleanup(setattr, sys, "argv", sys.argv)
        for row, why in (
            ("o.ok", "interrupted (SIGINT, Ctrl-C) while the outputs were written"),
            ("n.crash", "new aba-run n.crash exited 3"),
        ):
            verdicts: list[str] = []

            def interrupted_once(report, verdicts=verdicts):
                verdicts.append(report["verdict"])
                if len(verdicts) == 1:
                    raise KeyboardInterrupt
                return saved(report)

            aba.summary_text = interrupted_once
            out = self.dir / f"raised-{row}"
            argv = [
                "aba.py",
                *self.sides(),
                "--rows",
                row,
                "--out",
                str(out),
                *self.QUICK,
            ]
            err = io.StringIO()
            sys.argv = argv
            with (
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(err),
            ):
                status = aba.main()
            self.assertEqual(status, 1, row)
            rep = report_of(out)
            self.assertIn(why, rep["aborted"], row)
            self.assertEqual(
                (rep["verdict"], rep["conforming"]), ("aborted", False), row
            )
            self.assertEqual(verdicts[-1], "aborted", row)
            summary = (out / "summary.md").read_text().splitlines()
            self.assertIn("overall verdict **aborted**", summary[3], row)
            self.assertIn(f"* **aborted**: {rep['aborted']}", summary, row)
            self.assertIn(
                "aborted: interrupted (SIGINT, Ctrl-C) while the outputs were written",
                err.getvalue(),
                row,
            )

    def test_second_interrupt_while_the_first_is_handled(self):
        # Ctrl-C during a block (the fake harness sends it), and again while the run records
        # that interrupt: the second is held, so the run aborts for the first, with every
        # output written (exit 1, not the 130 of an interrupt before the first block)
        class SecondInterrupt(dict):
            def __setitem__(self, key, value):
                if key == "aborted":
                    signal.raise_signal(signal.SIGINT)
                super().__setitem__(key, value)

        saved = aba.start_run

        def start_run(*args):
            out, report, *rest = saved(*args)
            return (out, SecondInterrupt(report), *rest)

        self.addCleanup(setattr, aba, "start_run", saved)
        self.addCleanup(setattr, sys, "argv", sys.argv)
        aba.start_run = start_run
        out = self.dir / "second"
        sys.argv = ["aba.py", *self.sides(), "--rows", "n.sigint", "--out", str(out)]
        sys.argv += self.QUICK
        err = io.StringIO()
        with (
            contextlib.redirect_stdout(io.StringIO()),
            contextlib.redirect_stderr(err),
        ):
            status = aba.main()
        self.assertEqual(status, 1, err.getvalue())
        rep = report_of(out)
        self.assertEqual(
            rep["aborted"], "interrupted (SIGINT, Ctrl-C) during new aba-run n.sigint"
        )
        self.assertEqual((rep["verdict"], rep["conforming"]), ("aborted", False))
        self.assertIn("overall verdict **aborted**", self.summary("second")[3])
        self.assertEqual(
            sorted(os.listdir(out)),
            ["blocks", "blocks.csv", "meta.txt", "report.json", "summary.md"],
        )
        self.assertIs(signal.getsignal(signal.SIGINT), signal.default_int_handler)

    # Every way a run aborts: the step it aborts in, the row that makes it abort, the options
    # that needs, and the reason ({new} stands for the candidate harness's path).
    ABORTS = (
        ("warm-up", "n.warmfail", (), "new aba-calibrate n.warmfail"),
        ("warm-up", "n.hugewarm", (), "n.hugewarm expected duration is inf"),
        ("warm-up", "n.rsswarm", (), "n.rsswarm rss_tree_bytes is inf"),
        (
            "warm-up",
            "n.intwarm",
            INTWARM_PIN,
            "n.intwarm expected duration is inf",
        ),
        (
            "warm-up",
            "n.sigint_warm",
            (),
            "interrupted (SIGINT, Ctrl-C) during new aba-calibrate n.sigint_warm",
        ),
        (
            "warm-up",
            "n.nan",
            (),
            "old aba-calibrate n.nan (nothing on stderr): malformed output line (unknown"
            " key 'x'): ",
        ),
        ("plan", "n.longplan", (), "ends beyond any date"),
        (
            "timing",
            "n.garbled",
            (),
            "new aba-run n.garbled (nothing on stderr): malformed output line (Expecting"
            " property name",
        ),
        ("timing", "n.negative", (), "median_s -0.0005 is not positive"),
        ("timing", "n.other", (), "the record names row 'n.fast', not 'n.other'"),
        ("timing", "n.crash", (), "new aba-run n.crash exited 3 (nothing on stderr)"),
        (
            "timing",
            "n.fewer",
            (),
            "new aba-run n.fewer: the harness ran 4399 iterations, 4400 requested",
        ),
        (
            "timing",
            "n.tamper",
            (),
            "new: harness changed after aba-run n.tamper: {new} ",
        ),
        (
            "timing",
            "n.cross_tamper",
            (),
            "new: harness changed before aba-run n.cross_tamper: {new} ",
        ),
        (
            "timing",
            "n.vanish",
            (),
            "new: harness {new} unreadable after aba-run n.vanish",
        ),
        ("timing", "n.noexec", (), "new aba-run n.noexec: cannot start {new} "),
        ("timing", "n.huge", (), "n.huge: the mean of A1 and A2 is inf"),
        (
            "timing",
            "n.sigint",
            (),
            "interrupted (SIGINT, Ctrl-C) during new aba-run n.sigint",
        ),
        (
            "timing",
            "n.extra",
            (),
            "new aba-run n.extra (nothing on stderr): malformed output line (unknown key"
            " 'x'): ",
        ),
        (
            "timing",
            "n.deep",
            (),
            # Python 3.14's reader reads it (the key is refused), older ones refuse the nesting
            "new aba-run n.deep (nothing on stderr): malformed output line (",
        ),
        (
            "end",
            "n.end_tamper",
            (),
            "new: harness changed at the end of the run: {new} ",
        ),
    )
    # The status of the rows before, at and after the aborting one, by step.
    REACHED = {
        "warm-up": ("planned", "warm-up", "not run"),
        "plan": ("planned", "planned", "planned"),
        "timing": ("measured", "timing", "planned"),
        "end": ("measured", "measured"),
    }

    def test_every_abort_with_several_rows(self):
        # Each abort with a row before the aborting one and one after it (none after the end
        # of the run): report.json and summary.md still list every selected row with the step
        # it reached, and the summary's plan says which rows did not finish.
        for step, row, extra, why in self.ABORTS:
            self.fakes()  # the previous case changed, removed or locked the candidate
            rows = ["n.fast", row, *(["o.ok"] if step != "end" else [])]
            # the end-of-run change needs all seven triples (the last A2 changes the candidate)
            short = () if step == "end" else ("--triples", "2", "--allow-short")
            out = f"several-{row}"
            rep = self.aborted(out, "--rows", ",".join(rows), *short, *extra)
            new = self.dir / "new" / "fake-new"
            self.assertIn(why.format(new=new), rep["aborted"], row)
            want = dict(zip(rows, self.REACHED[step], strict=True))
            self.assertEqual({r["row"]: r["status"] for r in rep["rows"]}, want, row)
            # blocks.csv: the blocks of the finished warm-up pairs and triples, row by row
            with open(self.dir / out / "blocks.csv", newline="") as f:
                slots = [(b["row"], b["slot"]) for b in csv.DictReader(f)]
            finished = [
                (r["row"], s)
                for r in rep["rows"]
                for s in [
                    *(["W-A", "W-B"] if "warmup" in r else []),
                    *["A1", "B", "A2"] * len(r.get("triples", [])),
                ]
            ]
            self.assertEqual(slots, finished, row)
            summary = self.summary(out)
            table = [ln for ln in summary if ln.startswith(("| n.", "| o."))]
            self.assertEqual([ln.split(" | ")[0][2:] for ln in table], rows, row)
            text = "\n".join(summary)
            for name, status in want.items():
                self.assertIn(f"{name} ({status})", text, row)
            if step == "timing":
                actual = {
                    ln.split(":")[0][2:]: ln.rsplit(", actual ", 1)[1]
                    for ln in summary
                    if ln.startswith("* ") and ", actual " in ln
                }
                self.assertEqual(
                    (actual[row], actual["o.ok"]), ("not finished", "not run"), row
                )
        # an option that aborts the first row's warm-up step: nothing after it ran
        self.fakes()
        rep = self.aborted(
            "several-count",
            "--rows",
            "n.fast,o.ok",
            "--target-seconds",
            "1e308",
        )
        self.assertEqual([r["status"] for r in rep["rows"]], ["warm-up", "not run"])


def git(cwd: Path, *args: str) -> str:
    """git without the user's hooks, signing or identity."""
    cmd = [
        "git",
        "-c",
        "core.hooksPath=/dev/null",
        "-c",
        "commit.gpgsign=false",
        "-c",
        "user.name=test",
        "-c",
        "user.email=test@example.com",
        *args,
    ]
    return subprocess.run(
        cmd, cwd=cwd, check=True, capture_output=True, text=True
    ).stdout


def tree_digest(root: Path) -> str:
    """build.sh's core_sources_sha256, computed independently (dot-files left out)."""
    lines = sorted(
        (p.relative_to(root).as_posix(), hashlib.sha256(p.read_bytes()).hexdigest())
        for sub in ("include", "src")
        for p in (root / sub).rglob("*")
        if p.is_file()
        and not any(part.startswith(".") for part in p.relative_to(root).parts)
    )
    return hashlib.sha256("".join(f"{h}  {p}\n" for p, h in lines).encode()).hexdigest()


# The fake compiler logs its arguments and writes a small binary; xcrun names $FAKE_SDK (none
# when it is empty); vtool prints one LC_BUILD_VERSION, or fails when $FAKE_VTOOL is "fail".
FAKE_CXX = """#!/bin/sh
if [ "$1" = --version ]; then echo "fake clang version 1.0"; echo second line; exit 0; fi
if [ -n "$FAKE_CXX_FAIL" ]; then echo "fake: error: does not compile" >&2; exit 1; fi
printf '%s\\n' "$@" > "$FAKE_LOG"
while [ $# -gt 0 ]; do if [ "$1" = -o ]; then printf 'binary\\n' > "$2"; fi; shift; done
"""
FAKE_XCRUN = """#!/bin/sh
[ -n "$FAKE_SDK" ] || exit 1
case "$1" in --show-sdk-path) echo "$FAKE_SDK" ;; --show-sdk-version) echo 27.0 ;; esac
"""
FAKE_VTOOL = """#!/bin/sh
[ "$FAKE_VTOOL" = fail ] && exit 1
printf '%s:\\nLoad command 9\\n      cmd LC_BUILD_VERSION\\n platform MACOS\\n    minos 15.0\\n      sdk 27.0\\n' "$2"
if [ "$FAKE_VTOOL" = several ]; then
  printf 'Load command 10\\n platform OTHER\\n    minos 99.0\\n      sdk 100.0\\n'
fi
"""


class BuildScript(unittest.TestCase):
    """build.sh in a scratch git repository with this checkout's build.sh and driver files, a
    small include/ and src/, and a fake compiler, xcrun and vtool on PATH."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        top = Path(self.tmp.name).resolve()
        self.repo = top / "repo"
        perf = self.repo / "dev" / "perf"
        perf.mkdir(parents=True)
        for name in ("build.sh", *DRIVER):
            shutil.copy2(HERE / name, perf / name)
        for path, text in (
            ("include/netgraph/core/a.hpp", "a\n"),
            ("src/a.cpp", "a\n"),
            ("src/b.cpp", "b\n"),
            # dot-files: left out of core_sources_sha256
            ("src/.hidden", "h\n"),
            ("include/.d/c.hpp", "c\n"),
        ):
            (self.repo / path).parent.mkdir(parents=True, exist_ok=True)
            (self.repo / path).write_text(text)
        git(self.repo, "init", "-q")
        git(self.repo, "add", ".")
        git(self.repo, "commit", "-q", "-m", "base")
        self.base = git(self.repo, "rev-parse", "HEAD").strip()
        # HEAD is a second revision with another src/a.cpp: a revision build of the first one
        # must take the first one's sources, not HEAD's
        (self.repo / "src" / "a.cpp").write_text("a at HEAD\n")
        git(self.repo, "commit", "-q", "-am", "head")
        self.head = git(self.repo, "rev-parse", "HEAD").strip()
        self.bin = top / "bin"
        self.bin.mkdir()
        for name, text in (
            ("clang++", FAKE_CXX),
            ("mycxx", FAKE_CXX),  # a CXX that is not the PATH's clang++
            ("xcrun", FAKE_XCRUN),
            ("vtool", FAKE_VTOOL),
        ):
            (self.bin / name).write_text(text)
            (self.bin / name).chmod(0o755)
        self.sdk = top / "sdk"
        self.sdk.mkdir()
        self.log = top / "compiler.log"
        self.top = top
        self.flags = (
            re.search(r"^flags=\((.*)\)$", (HERE / "build.sh").read_text(), re.M)
            .group(1)
            .split()
        )

    def tearDown(self):
        self.tmp.cleanup()

    def build(self, *args: str, **env: str) -> subprocess.CompletedProcess:
        base = {
            k: v
            for k, v in os.environ.items()
            if k not in ("CXX", "SDKROOT", "MACOSX_DEPLOYMENT_TARGET")
        }
        base.update(
            PATH=f"{self.bin}{os.pathsep}{os.environ['PATH']}",
            FAKE_LOG=str(self.log),
            FAKE_SDK=str(self.sdk),
            GIT_CONFIG_GLOBAL=os.devnull,
        )
        base.update(env)
        return subprocess.run(
            ["bash", str(self.repo / "dev" / "perf" / "build.sh"), *args],
            cwd=self.top,
            capture_output=True,
            text=True,
            env=base,
            timeout=120,
        )

    def record(self, out: Path) -> dict:
        return json.loads((out / "harness.json").read_text())

    def test_refusals(self):
        for args, why in (
            ((), "usage: build.sh REV|worktree OUT_DIR"),
            (("worktree",), "usage"),
            (("", "out"), "usage"),
            (("worktree", "out", "extra"), "usage"),
            (("no-such-revision", "out"), "no revision no-such-revision"),
        ):
            proc = self.build(*args)
            self.assertEqual(proc.returncode, 2, (args, proc.stderr))
            self.assertIn(why, proc.stderr)
            self.assertFalse((self.top / "out").exists())
        used = self.top / "used"
        used.mkdir()
        (used / "x").write_text("x")
        proc = self.build("worktree", str(used))
        self.assertEqual(proc.returncode, 2, proc.stderr)
        self.assertIn("must be absent or an empty directory", proc.stderr)
        self.assertEqual(os.listdir(used), ["x"])
        self.assertFalse(self.log.exists())  # the compiler never ran
        # a failing compile: the compiler's status, its message, no binary and no record
        for rev in ("worktree", self.base):
            out = self.top / f"failed-{rev}"
            proc = self.build(rev, str(out), FAKE_CXX_FAIL="1")
            self.assertEqual(proc.returncode, 1, proc.stderr)
            self.assertIn("fake: error: does not compile", proc.stderr)
            self.assertFalse((out / "harness").exists())
            self.assertFalse((out / "harness.json").exists())

    def test_failed_record_leaves_no_harness(self):
        # the compile succeeds, the record's writer fails (before writing, or half-way): a
        # harness without its record would be timed without the driver check, so none is left
        for name, script in (
            ("silent", "exit 1"),
            ("partial", 'printf "{" > "$OUT.json"; exit 1'),
        ):
            (self.bin / "python3").write_text(f"#!/bin/sh\n{script}\n")
            (self.bin / "python3").chmod(0o755)
            for rev in ("worktree", self.base):
                out = self.top / f"{name}-{rev}"
                proc = self.build(rev, str(out))
                self.assertEqual(proc.returncode, 1, (name, rev, proc.stderr))
                kept = ["src"] if rev == self.base else []
                self.assertEqual(sorted(os.listdir(out)), kept, (name, rev))
                self.assertIn(
                    f"-o {out / 'harness'}", self.log.read_text().replace("\n", " ")
                )

    def test_worktree_build_and_record(self):
        (self.repo / "src" / "a.cpp").write_text("changed\n")
        (self.repo / "src" / "new.cpp").write_text("new\n")
        (self.top / "empty").mkdir()
        # SDKROOT must name a real SDK on macOS: git there is a shim that reads it
        sdkroot = "/sdkroot"
        if platform.system() == "Darwin":
            cmd = ["/usr/bin/xcrun", "--show-sdk-path"]
            sdkroot = subprocess.run(cmd, capture_output=True, text=True).stdout.strip()
        proc = self.build("worktree", "empty", SDKROOT=sdkroot)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        out = self.top / "empty"
        rec = self.record(out)
        self.assertEqual(sorted(os.listdir(out)), ["harness", "harness.json"])
        argv = self.log.read_text().splitlines()
        src = sorted(str(p) for p in (self.repo / "src").glob("*.cpp"))
        self.assertEqual(
            argv,
            [
                *self.flags,
                f"-I{self.repo}/include",
                "dev/perf/harness.cpp",
                *src,
                "-o",
                str(out / "harness"),
            ],
        )
        perf = self.repo / "dev" / "perf"
        want = {
            "binary": str(out / "harness"),
            "binary_sha256": hashlib.sha256(b"binary\n").hexdigest(),
            "rev": "worktree",
            "commit": self.head,
            "dirty_core_paths": 2,
            "core_sources_sha256": tree_digest(self.repo),
            "driver_sha256": checkout_driver(perf),
            "compiler": "fake clang version 1.0",
            "flags": " ".join(self.flags),
            "compiler_path": str(self.bin / "clang++"),
            "cxx_env": None,
            "isysroot": None,
            "sdk": {
                "path": str(self.sdk),
                "version": "27.0",
                "sdkroot_env": sdkroot,
                "source": "xcrun --show-sdk-path / --show-sdk-version",
            },
            "deployment_target": {
                "requested": None,
                "effective_minos": "15.0",
                "binary_sdk": "27.0",
                "platform": "MACOS",
                "raw": "Load command 9; cmd LC_BUILD_VERSION; platform MACOS; minos 15.0; sdk 27.0",
                "source": "vtool -show-build",
            },
            "built_by": "dev/perf/build.sh",
        }
        self.assertEqual(rec, want)
        self.assertIn(f"worktree {self.head[:12]} core ", proc.stdout)

    def test_vtool_records_the_first_values_and_all_lines(self):
        proc = self.build("worktree", "several", FAKE_VTOOL="several")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        target = self.record(self.top / "several")["deployment_target"]
        self.assertEqual(
            [target[k] for k in ("platform", "effective_minos", "binary_sdk")],
            ["MACOS", "15.0", "27.0"],
        )
        self.assertEqual(
            target["raw"],
            "Load command 9; cmd LC_BUILD_VERSION; platform MACOS; minos 15.0; sdk 27.0; Load command 10; platform OTHER; minos 99.0; sdk 100.0",
        )

    def test_revision_build_with_cxx(self):
        (self.repo / "src" / "a.cpp").write_text("uncommitted\n")  # not in the revision
        cxx = str(self.bin / "mycxx")
        proc = self.build(
            self.base[:10],
            "rev",
            CXX=cxx,
            MACOSX_DEPLOYMENT_TARGET="15.0",
            FAKE_VTOOL="fail",
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        out = self.top / "rev"
        rec = self.record(out)
        extracted = sorted(
            p.relative_to(out / "src").as_posix()
            for p in (out / "src").rglob("*")
            if p.is_file()
        )
        self.assertEqual(
            extracted,
            [
                "include/.d/c.hpp",
                "include/netgraph/core/a.hpp",
                "src/.hidden",
                "src/a.cpp",
                "src/b.cpp",
            ],
        )
        self.assertEqual((out / "src" / "src" / "a.cpp").read_text(), "a\n")
        darwin = platform.system() == "Darwin"
        flags = [*self.flags, *(["-isysroot", str(self.sdk)] if darwin else [])]
        argv = self.log.read_text().splitlines()
        self.assertEqual(
            argv[: len(flags) + 2],
            [*flags, f"-I{out}/src/include", "dev/perf/harness.cpp"],
        )
        self.assertEqual(
            argv[len(flags) + 2 : -2],
            [str(out / "src" / "src" / n) for n in ("a.cpp", "b.cpp")],
        )
        self.assertEqual(
            (rec["rev"], rec["commit"], rec["dirty_core_paths"]),
            (self.base[:10], self.base, 0),
        )
        self.assertEqual(rec["core_sources_sha256"], tree_digest(out / "src"))
        # the first revision's tree, not HEAD's (extracted the same way for the comparison)
        head = self.top / "head"
        head.mkdir()
        tar = subprocess.run(
            ["git", "archive", self.head, "include", "src"],
            cwd=self.repo,
            check=True,
            capture_output=True,
        ).stdout
        subprocess.run(["tar", "-x", "-C", str(head)], input=tar, check=True)
        self.assertNotEqual(rec["core_sources_sha256"], tree_digest(head))
        self.assertEqual((rec["cxx_env"], rec["flags"]), (cxx, " ".join(flags)))
        self.assertEqual(rec["compiler_path"], cxx)  # CXX, not the PATH's clang++
        self.assertEqual(rec["isysroot"], str(self.sdk) if darwin else None)
        self.assertIsNone(rec["sdk"]["sdkroot_env"])
        self.assertEqual(
            rec["deployment_target"],
            {
                "requested": "15.0",
                "effective_minos": None,
                "binary_sdk": None,
                "platform": None,
                "raw": None,
                "source": "no LC_BUILD_VERSION read (no vtool, or it printed none)",
            },
        )

    def test_no_sdk(self):
        # CXX without an SDK: refused on macOS; the clang++ shim needs none
        proc = self.build("worktree", "cxx", CXX=str(self.bin / "clang++"), FAKE_SDK="")
        if platform.system() == "Darwin":
            self.assertEqual(proc.returncode, 2, proc.stderr)
            self.assertIn("names no macOS SDK", proc.stderr)
            self.assertFalse(self.log.exists())
            # refused before OUT_DIR is made or a revision extracted into it, so the
            # corrected retry can use the same OUT_DIR
            for rev in ("worktree", self.base):
                out = self.top / f"no-sdk-{rev}"
                proc = self.build(
                    rev, str(out), CXX=str(self.bin / "clang++"), FAKE_SDK=""
                )
                self.assertEqual(proc.returncode, 2, proc.stderr)
                self.assertFalse(out.exists())
                proc = self.build(rev, str(out))
                self.assertEqual(proc.returncode, 0, proc.stderr)
        else:
            self.assertEqual(proc.returncode, 0, proc.stderr)
        proc = self.build("worktree", "shim", FAKE_SDK="")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        sdk = self.record(self.top / "shim")["sdk"]
        self.assertEqual(
            (sdk["path"], sdk["version"], sdk["source"]),
            (
                None,
                None,
                "no xcrun SDK on this host: the compiler's own system headers",
            ),
        )


# The test rows appended to the harness's own: a constant digest (a slow and a fast iteration),
# a digest that changes at the first or the second timed iteration, an iteration that throws (the
# untimed first one or a timed one), and one row per setup check of aba_rows.hpp fed an input it
# must refuse. The extra command test-times compares IterationTimes with the sorted list of the
# same durations.
HARNESS_TEST = r"""
#include "aba_rows.hpp"

#include <memory>
#include <random>
#include <stdexcept>

using namespace abarows;

namespace {
// About 20,000 additions: an iteration well above the clock's resolution (a median of 0 is refused).
Iter constant() {
  return Iter([] {
    volatile std::uint64_t x = 0;
    for (int i = 0; i < 20000; ++i) x = x + 1;
    return std::uint64_t {7};
  });
}

Row check(const char* name, std::function<void()> f) {
  return Row {name, "test", "a setup check that must fail", [f]() -> Iter { f(); return constant(); }};
}

std::shared_ptr<Graph> graph(const char* input) {
  return std::make_shared<Graph>(abaapi::build(make_storage_input(input)));
}

// Durations below, at and above IterationTimes::kShort, many repeated, and more long ones than
// its preallocated slots: every statistic the block line prints must equal the sorted list's.
int test_times() {
  std::mt19937_64 rng(5);
  const std::uint64_t k = IterationTimes::kShort;
  int cases = 0;
  for (std::size_t n : {1, 2, 3, 10, 1001, 200000}) {
    IterationTimes times;
    std::vector<double> sorted;
    for (std::size_t i = 0; i < n; ++i) {
      const std::uint64_t pick = rng() % 4;
      const std::uint64_t t = pick == 0 ? k - 1 + rng() % 2 : pick == 1 ? 700 + rng() % 20 : pick == 2 ? rng() % k : k + rng() % 1000000000;
      const auto d = Clock::duration(Clock::rep(t));
      times.add(d);
      sorted.push_back(std::chrono::duration<double>(d).count());
    }
    times.finish();
    std::sort(sorted.begin(), sorted.end());
    double sum = 0;
    for (double x : sorted) sum += x;
    std::vector<std::size_t> ks = {0, n / 2, n - 1, (n - 1) / 2, std::size_t(0.1 * double(n - 1) + 0.5)};
    for (int j = 0; j < 200; ++j) ks.push_back(std::size_t(rng() % n));
    bool same = times.size() == n && times.sum() == sum;
    for (std::size_t i : ks) same = same && times.at(i) == sorted[i];
    if (!same) { std::printf("IterationTimes differs from the sorted list at n = %zu\n", n); return 1; }
    ++cases;
  }
  std::printf("IterationTimes equals the sorted list in %d cases\n", cases);
  return 0;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "test-times") return test_times();
  auto r = rows();
  r.push_back({"test.constant", "test", "every iteration the same digest", [] { return constant(); }});
  r.push_back({"test.fast", "test", "every iteration the same digest, at once", [] {
                 return Iter([] { return std::uint64_t {7}; });
               }});
  r.push_back({"test.changes", "test", "the digest changes at the second timed iteration", [] {
                 auto n = std::make_shared<int>(0);
                 return Iter([n] { return std::uint64_t(++*n > 2); });
               }});
  r.push_back({"test.changes_first", "test", "the digest changes at the first timed iteration", [] {
                 auto n = std::make_shared<int>(0);
                 return Iter([n] { return std::uint64_t(++*n > 1); });
               }});
  r.push_back({"test.throws_first", "test", "the untimed first iteration throws", []() -> Iter {
                 return Iter([]() -> std::uint64_t { throw std::runtime_error("the first iteration failed"); });
               }});
  r.push_back({"test.throws", "test", "the second timed iteration throws", [] {
                 auto n = std::make_shared<int>(0);
                 return Iter([n]() -> std::uint64_t {
                   if (++*n > 2) throw std::runtime_error("an iteration failed");
                   return 0;
                 });
               }});
  r.push_back(check("test.query_input", [] {
    auto a = gen_random(16384, 8, 8);  // spf.dag.far16k's arrays, one cost changed
    a.cost[0] += 1;
    require_frozen(a, kFar16kInput, "far16k");
  }));
  r.push_back(check("test.storage_cost", [] {
    auto a = make_storage_input("width65536");
    a.cost[0] = a.cost[0] == 65534 ? 65535 : 65534;  // the shape's counts unchanged, the digest not
    checked_storage_input("width65536", a);
  }));
  r.push_back(check("test.storage_capacity", [] {
    auto a = make_storage_input("width65536");
    a.cap[0] += 1;  // only the full-content digest sees a capacity
    checked_storage_input("width65536", a);
  }));
  r.push_back(check("test.storage_unknown", [] { checked_storage_input("nonesuch", make_storage_input("width65536")); }));
  r.push_back(check("test.roots_rule", [] { frozen_roots("spf.to.width65536", 1000); }));
  r.push_back(check("test.roots_unknown", [] { frozen_roots("nonesuch", 16384); }));
  // the library's width wrong (wgrid256's graph: W = 32, buckets), the input's right
  r.push_back(check("test.queue_library", [] { require_queue(*graph("wgrid256"), make_storage_input("width65536"), 65536, true); }));
  // the library's choice right, the input's width wrong
  r.push_back(check("test.queue_input", [] { require_queue(*graph("width65536"), make_storage_input("width65537"), 65536, true); }));
  // the width right everywhere, buckets or heap wrong
  r.push_back(check("test.queue_kind", [] { require_queue(*graph("width65536"), make_storage_input("width65536"), 65536, false); }));
  return harness_main(argc, argv, r);
}
"""


class Harness(unittest.TestCase):
    """The harness's driver (harness_main) and the setup checks of aba_rows.hpp, compiled
    unoptimised with this checkout's Core and the test rows above."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        top = Path(cls.tmp.name)
        (top / "harness_test.cpp").write_text(HARNESS_TEST)
        cls.binary = top / "harness_test"
        cmd = [
            "clang++",
            "-std=c++20",
            "-O0",
            "-Wall",
            "-Wextra",
            f"-I{REPO / 'include'}",
            f"-I{HERE}",
            str(top / "harness_test.cpp"),
            *sorted(glob.glob(str(REPO / "src" / "*.cpp"))),
            "-o",
            str(cls.binary),
        ]
        cls.compiled = subprocess.run(cmd, capture_output=True, text=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def harness(self, *args: str) -> subprocess.CompletedProcess:
        # a timeout turns a harness that never stops (a count it should refuse) into a failure
        cmd = [str(self.binary), *args]
        return subprocess.run(cmd, capture_output=True, text=True, timeout=60)

    def test_compiles_without_a_diagnostic(self):
        self.assertEqual(self.compiled.returncode, 0, self.compiled.stderr[-2000:])
        self.assertEqual(self.compiled.stderr, "")

    def test_list(self):
        proc = self.harness("aba-list")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = aba.row_list(proc.stdout.encode(), "aba-list")
        real = {k: v for k, v in rows.items() if not k.startswith("test.")}
        self.assertEqual(len(real), 65)
        sums = sha256sums(HERE / "queries")
        lists = {}
        for name in sums:
            rec = json.loads((HERE / "queries" / name).read_text())
            lists[rec["row"]] = f"queries/{name} sha256:{sums[name]}"
        self.assertEqual(
            {k: v["queries"] for k, v in real.items() if "queries" in v}, lists
        )
        self.assertTrue(all(v["group"] and v["desc"] for v in real.values()))

    def test_blocks(self):
        proc = self.harness("aba-run", "test.constant", "5")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        run = aba.parse_block(proc.stdout.encode(), "aba-run", "test.constant")
        self.assertEqual((run["iterations"], run["digest"]), (5, "0000000000000007"))
        self.assertLessEqual(run["min_s"], run["median_s"])
        self.assertLessEqual(run["median_s"], run["max_s"])
        proc = self.harness("aba-calibrate", "test.constant", "0.05")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        cal = aba.parse_block(proc.stdout.encode(), "aba-calibrate", "test.constant")
        self.assertGreaterEqual(cal["block_s"], 0.05)
        # its own peak RSS in bytes (macOS reports bytes, Linux KiB): a process that loaded the
        # driver holds some MiB, never a GiB (a KiB/byte mix-up is a factor of 1024); it starts
        # no child
        for rec in (run, cal):
            self.assertTrue(MIB < rec["rss_self_bytes"] < 1024 * MIB, rec)
            self.assertEqual(rec["rss_children_bytes"], 0)
        # every real row, one timed iteration: its setup (the frozen-input, root and queue
        # checks included) and its adapter calls run against this checkout's Core
        listed = aba.row_list(self.harness("aba-list").stdout.encode(), "aba-list")
        real = [row for row in listed if not row.startswith("test.")]
        with concurrent.futures.ThreadPoolExecutor(os.cpu_count() or 1) as pool:
            procs = list(pool.map(lambda row: self.harness("aba-run", row, "1"), real))
        self.assertEqual(len(procs), 65)
        for row, proc in zip(real, procs, strict=True):
            self.assertEqual(proc.returncode, 0, (row, proc.stderr))
            self.assertEqual(
                aba.parse_block(proc.stdout.encode(), row, row)["iterations"], 1
            )

    def test_memory_does_not_grow_with_the_count(self):
        # the per-iteration bookkeeping is allocated before the setup and does not grow: the
        # peak of a block of 2,000,000 iterations is that of a block of 1,000 (keeping every
        # time in a list would add 16 MiB, and a faster candidate would hold more of it)
        peaks = []
        for count in ("1000", "2000000"):
            proc = self.harness("aba-run", "test.fast", count)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            # json, not parse_block: an iteration this fast has a median of 0 clock ticks
            rec = json.loads(proc.stdout)
            self.assertEqual(rec["iterations"], int(count))
            peaks.append(rec["rss_self_bytes"])
        self.assertLess(peaks[1] - peaks[0], MIB / 4, peaks)

    def test_statistics_equal_the_sorted_list(self):
        proc = self.harness("test-times")
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertEqual(
            proc.stdout, "IterationTimes equals the sorted list in 6 cases\n"
        )

    def test_exit_statuses(self):
        usage = "usage: harness aba-list | aba-run <row> <iterations> | aba-calibrate <row> <seconds>"
        for args in (
            (),
            ("aba-list", "x"),
            ("aba-run",),
            ("aba-run", "test.constant"),
            ("aba-run", "test.constant", "1", "x"),
            ("aba-bogus", "test.constant", "1"),
            ("aba-run", "test.constant", "0"),
            ("aba-run", "test.constant", "-1"),
            ("aba-run", "test.constant", "1x"),
            ("aba-run", "test.constant", ""),
            ("aba-run", "test.constant", "99999999999999999999"),
            ("aba-calibrate", "test.constant", "0"),
            ("aba-calibrate", "test.constant", "nan"),
            ("aba-calibrate", "test.constant", "inf"),
            ("aba-calibrate", "test.constant", "1e999"),
            ("aba-calibrate", "test.constant", "0.1s"),
        ):
            proc = self.harness(*args)
            self.assertEqual((proc.returncode, proc.stdout), (2, ""), args)
            self.assertEqual(proc.stderr.strip(), usage, args)
        for args, status, why in (
            (("aba-run", "nonesuch", "1"), 2, "unknown row nonesuch"),
            (
                ("aba-run", "test.changes", "5"),
                3,
                "nondeterministic digest in test.changes at iteration 1",
            ),
            (
                ("aba-calibrate", "test.changes", "1"),
                3,
                "nondeterministic digest in test.changes",
            ),
            (
                ("aba-run", "test.changes_first", "5"),
                3,
                "nondeterministic digest in test.changes_first at iteration 0",
            ),
            (("aba-run", "test.throws", "5"), 1, "FAIL: an iteration failed"),
            # a library failure in the untimed first iteration is an iteration's, not setup's
            (
                ("aba-run", "test.throws_first", "5"),
                1,
                "FAIL: the first iteration failed",
            ),
        ):
            proc = self.harness(*args)
            self.assertEqual((proc.returncode, proc.stdout), (status, ""), args)
            self.assertIn(why, proc.stderr)

    def test_setup_checks_refuse(self):
        for row, why in (
            (
                "test.query_input",
                "the frozen list far16k was not frozen from this row's arrays",
            ),
            (
                "test.storage_cost",
                "storage input width65536 differs from its frozen shape",
            ),
            (
                "test.storage_capacity",
                "storage input width65536 differs from its frozen shape",
            ),
            ("test.storage_unknown", "storage input nonesuch is not frozen"),
            (
                "test.roots_rule",
                "the frozen roots of spf.to.width65536 differ from the rule",
            ),
            ("test.roots_unknown", "no frozen roots for nonesuch"),
            (
                "test.queue_library",
                "queue choice: the adapter reports W = 32 (buckets), the input gives"
                " W = 65536, the row expects W = 65536 (buckets)",
            ),
            (
                "test.queue_input",
                "queue choice: the adapter reports W = 65536 (buckets), the input gives"
                " W = 65537, the row expects W = 65536 (buckets)",
            ),
            (
                "test.queue_kind",
                "queue choice: the adapter reports W = 65536 (buckets), the input gives"
                " W = 65536, the row expects W = 65536 (heap)",
            ),
        ):
            proc = self.harness("aba-run", row, "1")
            self.assertEqual((proc.returncode, proc.stdout), (4, ""), row)
            self.assertIn(f"row {row} setup failed: {why}", proc.stderr)


class FrozenFiles(unittest.TestCase):
    """The frozen inputs: the SHA256SUMS files hash the frozen JSON files, and the values
    compiled into aba_rows.hpp are exactly theirs."""

    QUERIES = HERE / "queries"
    LISTS = {  # list -> (its pairs constant or None, its input digest constant)
        "near65k": ("kNear65k", "kNear65kInput"),
        "far16k": ("kFar16k", "kFar16kInput"),
        "canonical16k": ("kCanonical16k", "kCanonical16kInput"),
        "targetset_clos": (None, "kTargetSetInput"),
    }

    def block(self, name: str) -> str:
        text = (HERE / "aba_rows.hpp").read_text()
        begin, end = f"// BEGIN {name} (", f"// END {name}\n"
        self.assertEqual((text.count(begin), text.count(end)), (1, 1), name)
        return text[text.index(begin) : text.index(end)]

    def test_hashes(self):
        for directory, files in (
            (self.QUERIES, [f"{k}.json" for k in self.LISTS]),
            (self.QUERIES / "storage", ["inputs.json"]),
        ):
            got = sha256sums(directory)
            self.assertEqual(sorted(got), sorted(files), directory)
            self.assertEqual(
                sorted(p.name for p in directory.glob("*.json")), sorted(files)
            )
            for name, sha in got.items():
                want = hashlib.sha256((directory / name).read_bytes()).hexdigest()
                self.assertEqual(sha, want, name)

    def test_query_lists_block(self):
        block = self.block("FROZEN QUERY LISTS")
        sums = sha256sums(self.QUERIES)
        listed = re.findall(
            r'\{"([^"]+)", "queries/([^"]+) sha256:([0-9a-f]{64})"\}', block
        )
        self.assertEqual(len(listed), len(self.LISTS))
        counts = []
        for key, (pairs_name, digest_name) in self.LISTS.items():
            rec = json.loads((self.QUERIES / f"{key}.json").read_text())
            self.assertEqual(rec["list"], key)
            self.assertIn((rec["row"], f"{key}.json", sums[f"{key}.json"]), listed)
            digest = re.search(rf"\b{digest_name} = (\d+)ULL;", block)
            self.assertEqual(int(digest.group(1)), rec["input"]["input_digest"], key)
            if pairs_name is None:
                found = re.search(
                    r"std::array<int, (\d+)> kTargetSetDestinations\{([^}]*)\};", block
                )
                values = [int(x) for x in found.group(2).split(",")]
                self.assertEqual(values, rec["destinations"], key)
            else:
                found = re.search(
                    rf"std::array<std::pair<int, int>, (\d+)> {pairs_name}\{{\{{(.*?)\}}\}};",
                    block,
                    re.S,
                )
                values = [
                    [int(s), int(t)]
                    for s, t in re.findall(r"\{(\d+), (\d+)\}", found.group(2))
                ]
                self.assertEqual(values, rec["pairs"], key)
            self.assertEqual(int(found.group(1)), len(values), key)
            counts.append(len(values))
        self.assertEqual(sorted(counts), [8, 8, 50, 64])  # no list can be empty

    def test_storage_block(self):
        block = self.block("FROZEN STORAGE INPUTS")
        reg = json.loads((self.QUERIES / "storage" / "inputs.json").read_text())
        fields = list(next(iter(reg["inputs"].values()))["shape"])
        inputs, roots_part = block.split("kStorageRoots")
        compiled = re.findall(r'\{"([^"]+)", \{([^}]*)\}\},', inputs)
        self.assertEqual([name for name, _ in compiled], list(reg["inputs"]))
        for name, values in compiled:
            shape = reg["inputs"][name]["shape"]
            self.assertEqual(list(shape), fields, name)
            want = [
                ("true" if v else "false") if isinstance(v, bool) else str(int(v))
                for v in shape.values()
            ]
            got = [v.strip().removesuffix("ULL") for v in values.split(",")]
            self.assertEqual(got, want, name)
        roots = re.findall(r'\{"([^"]+)", \{([^}]*)\}\},', roots_part)
        searches = {r: s for r, s in reg["rows"].items() if s["kind"] == "search"}
        self.assertEqual([r for r, _ in roots], list(searches))
        for row, values in roots:
            n = reg["inputs"][searches[row]["input"]]["shape"]["n"]
            rule = [(i * 997 + 1) % n for i in range(8)]
            self.assertEqual([int(v) for v in values.split(",")], rule, row)
            self.assertEqual(searches[row]["roots"], rule, row)
        # every row of the registry is an input's row, and the counts are the frozen ones
        self.assertEqual(
            sorted(reg["rows"]),
            sorted(r for v in reg["inputs"].values() for r in v["rows"]),
        )
        self.assertEqual((len(compiled), len(reg["rows"]), len(roots)), (17, 20, 6))

    def test_construct_table(self):
        # the storage construction rows' table in aba_rows.hpp: each row on its registry input
        reg = json.loads((self.QUERIES / "storage" / "inputs.json").read_text())
        text = (HERE / "aba_rows.hpp").read_text()
        start = text.index("const C constructs[] = {")
        table = re.findall(
            r'\{"(construct\.[^"]+)", "([^"]+)", "',
            text[start : text.index("};", start)],
        )
        want = {
            r: s["input"] for r, s in reg["rows"].items() if s["kind"] == "construct"
        }
        self.assertEqual(dict(table), want)
        self.assertEqual(len(table), 14)

    def test_search_table(self):
        # the storage search rows' table in aba_rows.hpp: the registry's input, direction and
        # queue, and every width up to 65,536 on the buckets, every wider one on the heap
        reg = json.loads((self.QUERIES / "storage" / "inputs.json").read_text())
        searches = {r: s for r, s in reg["rows"].items() if s["kind"] == "search"}
        text = (HERE / "aba_rows.hpp").read_text()
        table = re.findall(
            r'\{"([a-z]+\.[^"]+)", "([^"]+)", (true|false), (\d+), (true|false), "',
            text,
        )
        self.assertEqual(sorted(r for r, *_ in table), sorted(searches))
        for row, inp, reverse, width, buckets in table:
            spec = searches[row]
            self.assertEqual(spec["input"], inp, row)
            self.assertEqual(
                spec["direction"], "reverse" if reverse == "true" else "forward", row
            )
            self.assertEqual(
                spec["queue"], {"width": int(width), "buckets": buckets == "true"}, row
            )
            self.assertEqual(int(width) <= 65536, buckets == "true", row)
            shape = reg["inputs"][inp]["shape"]
            self.assertEqual(shape["max_cost"] // shape["gcd"] + 1, int(width), row)


if __name__ == "__main__":
    unittest.main()
