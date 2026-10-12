"""to_units (round half up, the only entry of physical numbers) and as_float."""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest
from helpers import exact

import netgraph_core as ngc


@pytest.mark.parametrize("q", [2.0**-12, 1.0, 2.0**-32, 2.0**32])
def test_half_up_boundaries(q):
    assert ngc.to_units(2.5 * q, q) == 3
    assert ngc.to_units(0.5 * q, q) == 1
    assert ngc.to_units(0.49 * q, q) == 0
    assert ngc.to_units(7.5 * q, q) == 8
    assert ngc.to_units(math.nextafter(7.5 * q, 0), q) == 7
    assert ngc.to_units(0.0, q) == 0
    assert ngc.to_units(-0.0, q) == 0
    assert ngc.as_float(12345, q) == 12345 * q
    assert type(ngc.to_units(2.5 * q, q)) is int


def test_arrays():
    q = 2.0**-12
    out = ngc.to_units([0.5 * q, 1.49 * q, 2.5 * q, 1e15], q)
    assert out.dtype == np.int64
    assert out.tolist() == [1, 1, 3, 4096 * 10**15]
    assert ngc.to_units(np.array([1.5, 2.5], dtype=np.float32), 1).tolist() == [2, 3]
    assert ngc.to_units(np.float64(2.5), 1) == 3
    assert ngc.to_units(np.array([], dtype=float), 1).tolist() == []
    assert ngc.as_float([1, 2, 3], 0.5).tolist() == [0.5, 1.0, 1.5]
    assert ngc.as_float(np.array([4096], dtype=np.int64), q).tolist() == [1.0]
    # An int gives a float, an integer array a float64 array, also for whole results.
    assert type(ngc.as_float(4096, q)) is float
    assert type(ngc.as_float(2**70, 1.0)) is float
    assert ngc.as_float([1, 2, 3], 1.0).dtype == np.float64
    assert ngc.as_float(np.array([4096], dtype=np.int64), q).dtype == np.float64


def test_rounding_is_per_value_not_additive():
    # Each value is converted to whole units on its own: two values of 0.6 units give 1 + 1,
    # their sum 1.2 gives 1.
    q = 2.0**-12
    assert ngc.to_units(0.6 * q, q) + ngc.to_units(0.6 * q, q) == 2
    assert ngc.to_units(1.2 * q, q) == 1


def test_exactness_limits():
    assert ngc.to_units(2.0**62 - 512, 1) == (1 << 62) - 512
    assert ngc.as_float((1 << 53) + 1, 1) == 2.0**53
    assert ngc.as_float(1 << 70, 2.0**-12) == 2.0**58


NAN = "to_units: NaN has no unit value"
INFINITE = "to_units: an infinite amount is not finite and has no unit value"
BOUND = "to_units: the amount is at or above the bound of the units domain (2^62 units)"
NOT_PHYSICAL = "to_units takes physical numbers (a float or an array of floats)"
RESOLUTION = "resolution must be a power of two in [2^-32, 2^32]"


@pytest.mark.parametrize(
    "x, message",
    [
        (float("nan"), NAN),
        (float("inf"), INFINITE),
        (-float("inf"), INFINITE),
        (-1e-300, "to_units: a negative amount has no unit value"),
        (2.0**62, BOUND),
        ([1.0, float("nan")], NAN),
        ("1.0", NOT_PHYSICAL),
        (b"1.0", NOT_PHYSICAL),
        (True, NOT_PHYSICAL),
        (["a"], "to_units takes physical numbers (got dtype <U1)"),
        (np.array([True]), "to_units takes physical numbers (got dtype bool)"),
        (np.zeros((1, 1)), "to_units takes a float or a 1-D array of floats"),
        ([[1.0]], "to_units takes a float or a 1-D array of floats"),
    ],
)
def test_rejected(x, message):
    with pytest.raises(ValueError, match=exact(message)):
        ngc.to_units(x, 1)


def test_resolution_validated():
    for q in (3.0, 2.0**-33, 2.0**33, 0.0, -1.0, math.nan, math.inf):
        with pytest.raises(ValueError, match=exact(RESOLUTION)):
            ngc.to_units(1.0, q)
        with pytest.raises(ValueError, match=exact(RESOLUTION)):
            ngc.as_float(1, q)
    with pytest.raises(ValueError, match=exact(BOUND)):
        ngc.to_units(1.7976931348623157e308, 2.0**32)


def test_as_float_rejects_non_integers():
    with pytest.raises(
        ValueError, match=exact("as_float: x must be an integer (got float)")
    ):
        ngc.as_float(1.5, 1)
    with pytest.raises(
        ValueError,
        match=exact("as_float: x entries must be integers (got dtype float64)"),
    ):
        ngc.as_float(np.array([1.5]), 1)


@pytest.mark.parametrize("q", [0.0, 3.0, math.nan, 2.0**33, -1.0, math.inf])
@pytest.mark.parametrize(
    "empty", [[], (), np.array([], dtype=np.float64), np.array([], dtype=np.int64)]
)
def test_empty_arrays_validate_the_resolution(q, empty):
    # The resolution is validated before the shape, also for empty input.
    with pytest.raises(ValueError, match=exact(RESOLUTION)):
        ngc.to_units(empty, q)
    with pytest.raises(ValueError, match=exact(RESOLUTION)):
        ngc.as_float(empty, q)
    assert ngc.to_units(empty, 1.0).tolist() == []
    if isinstance(empty, np.ndarray) and empty.dtype.kind == "f":
        with pytest.raises(  # not units, even when empty
            ValueError,
            match=exact("as_float: x entries must be integers (got dtype float64)"),
        ):
            ngc.as_float(empty, 1.0)
    else:
        assert ngc.as_float(empty, 1.0).tolist() == []


@pytest.mark.parametrize(
    "x, q",
    [
        (2**1024, 2.0**-32),
        (-(2**1024), 2.0**-32),
        (2**1024 - 1, 2.0**-32),
        (2**1055 + 3, 2.0**-32),
        (2**64 + 1, 1.0),
        (2**63, 2.0**-12),
        (3 * 2**900 + 1, 2.0**32),
        ((1 << 63) - 1, 2.0**-32),
    ],
)
def test_as_float_unbounded_ints_are_correctly_rounded(x, q):
    # Python Amounts are unbounded; the scaled value is the correctly rounded
    # double of the exact product, and no Python error leaks.
    expect = float(Fraction(x) * Fraction(q))
    assert ngc.as_float(x, q) == expect
    assert ngc.as_float([x, 1], q).tolist() == [expect, q]
    assert ngc.as_float(np.int64(1), q) == q


def test_as_float_beyond_the_finite_doubles():
    def beyond(x):
        return exact(
            f"as_float: {x} units times the resolution exceed the largest finite double"
        )

    with pytest.raises(ValueError, match=beyond(2**1100)):
        ngc.as_float(2**1100, 2.0**32)
    with pytest.raises(ValueError, match=beyond(-(2**1100))):
        ngc.as_float([1, -(2**1100)], 1.0)
    assert ngc.as_float(np.array([2**64 - 1], dtype=np.uint64), 1.0).tolist() == [
        2.0**64
    ]
