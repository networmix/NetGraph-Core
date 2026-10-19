"""Integer units: the one numeric domain (APPROACH revision 24, sections 2.1 and 2.3).

Every amount is a Python ``int`` number of units of the graph's resolution ``q``.
There is no finer grid; physical values enter only through ``to_units`` (round half
up) and leave only through ``as_float``. Costs are positive integers and are never
converted.
"""

from __future__ import annotations

import math
from fractions import Fraction

LIMIT = 1 << 62  # exclusive bound: capacities, their sum, finite requests, cost sums
VIRTUAL_CAP = (1 << 63) - 2  # virtual arcs never saturate (2.2, revision 20)
INT64_MAX = (1 << 63) - 1
UNREACHABLE = INT64_MAX  # path-cost sentinel (None in Python accessors)
NEVER_REACHED = -1  # potential sentinel of the private flow-core payload (2.2)
INTERNAL_CAP_MAX = (1 << 63) - 2  # internal flow-core cap domain [0, 2^63 - 2]
COUNTER_LIMIT = 1 << 127  # signed 128-bit cumulative counters (offered, drops)
INF = math.inf  # the only non-integer request: +inf


class ContractError(ValueError):
    """Invalid input rejected before any mutation (rule 13)."""


def check_resolution(q) -> Fraction:
    """A power of two in [2^-32, 2^32], given explicitly (no default in Core)."""
    if isinstance(q, bool):
        raise ContractError("resolution must be a number")
    if isinstance(q, float) and not math.isfinite(q):
        raise ContractError("resolution must be finite")
    try:
        f = Fraction(q)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"resolution {q!r} is not a number") from exc
    if f <= 0:
        raise ContractError("resolution must be positive")
    num, den = f.numerator, f.denominator
    if num & (num - 1) or den & (den - 1):
        raise ContractError(f"resolution {q!r} is not a power of two")
    exponent = num.bit_length() - den.bit_length()
    if not -32 <= exponent <= 32:
        raise ContractError(f"resolution 2^{exponent} outside [2^-32, 2^32]")
    return f


def to_units(x, q) -> int:
    """Round half up of ``x / q``; NaN, infinity and negative values rejected.

    The scalar form of the conversion; arrays convert element by element.
    """
    qf = check_resolution(q)
    if isinstance(x, bool):
        raise ContractError("physical amount must be a number")
    if isinstance(x, float) and not math.isfinite(x):
        raise ContractError("physical amount must be finite")
    value = Fraction(x)
    if value < 0:
        raise ContractError("physical amount must be non-negative")
    scaled = value / qf + Fraction(1, 2)
    return scaled.numerator // scaled.denominator


def to_units_array(xs, q) -> list[int]:
    return [to_units(x, q) for x in xs]


def as_float(units: int, q) -> float:
    """Units to a physical double: exact below 2^53 units (q is a power of two)."""
    qf = check_resolution(q)
    if not isinstance(units, int) or isinstance(units, bool):
        raise ContractError("as_float takes integer units")
    return float(Fraction(units) * qf)


def check_units(name: str, x, *, upper: int = LIMIT) -> int:
    """An integer amount in [0, upper)."""
    if not isinstance(x, int) or isinstance(x, bool):
        raise ContractError(f"{name} must be an integer number of units")
    if x < 0 or x >= upper:
        raise ContractError(f"{name} = {x} outside [0, {upper})")
    return x


def check_request(volume) -> int | None:
    """A request in units: an int in [0, 2^62) or +inf (returned as None)."""
    if isinstance(volume, float):
        if volume == math.inf:
            return None
        raise ContractError("requests are integer units or +inf")
    return check_units("request", volume)


def bits(x: int) -> int:
    return abs(x).bit_length()
