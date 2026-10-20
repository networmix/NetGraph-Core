"""to_units (round half up, the only entry of physical numbers) and as_float."""

from __future__ import annotations

import math
from decimal import Decimal
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


# A value outside the domain is a ValueError; a type to_units never takes (a string, a bool, an
# array of another dtype or rank) a TypeError.
@pytest.mark.parametrize(
    "x, cls, message",
    [
        (float("nan"), ValueError, NAN),
        (float("inf"), ValueError, INFINITE),
        (-float("inf"), ValueError, INFINITE),
        (-1e-300, ValueError, "to_units: a negative amount has no unit value"),
        (2.0**62, ValueError, BOUND),
        ([1.0, float("nan")], ValueError, NAN),
        ("1.0", TypeError, NOT_PHYSICAL),
        (b"1.0", TypeError, NOT_PHYSICAL),
        (True, TypeError, NOT_PHYSICAL),
        (np.True_, TypeError, NOT_PHYSICAL),
        (["a"], TypeError, "to_units takes physical numbers (got dtype <U1)"),
        (
            np.array([True]),
            TypeError,
            "to_units takes physical numbers (got dtype bool)",
        ),
        (
            np.zeros((1, 1)),
            TypeError,
            "to_units takes a float or a 1-D array of floats",
        ),
        ([[1.0]], TypeError, "to_units takes a float or a 1-D array of floats"),
    ],
)
def test_rejected(x, cls, message):
    with pytest.raises((TypeError, ValueError), match=exact(message)) as error:
        ngc.to_units(x, 1)
    assert error.type is cls


@pytest.mark.parametrize("q", [True, False, np.True_, np.False_])
def test_a_bool_resolution_is_a_wrong_type(q):
    # A bool at a numeric position is a TypeError (pybind11's text for incompatible arguments),
    # never read as 1.0 or 0.0; an int resolution is still taken.
    for call in (
        lambda: ngc.to_units(1.0, q),
        lambda: ngc.as_float(1, q),
        lambda: ngc.Graph(2, [0], [1], [1], [1], q),
    ):
        with pytest.raises(
            TypeError, match="incompatible (function|constructor) arguments"
        ) as error:
            call()
        assert error.type is TypeError
    assert ngc.to_units(1.0, 1) == 1 and ngc.as_float(1, 1) == 1.0
    assert ngc.Graph(2, [0], [1], [1], [1], 1).resolution() == 1.0


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
        TypeError, match=exact("as_float: x must be an integer (got float)")
    ):
        ngc.as_float(1.5, 1)
    with pytest.raises(
        TypeError,
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
            TypeError,
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


class _IndexFails:
    def __index__(self):
        raise ZeroDivisionError("own index")


class _FloatFails(np.float32):
    # A NumPy floating value whose own conversion to a float fails.
    def __float__(self):
        raise ZeroDivisionError("own float")


class _Index:
    # An integer known only by Python's integer protocol.
    def __init__(self, n):
        self.n = n

    def __index__(self):
        return self.n


class _FloatOnly:
    # A number by duck typing alone: not a type a number position takes.
    def __float__(self):
        return 1.0


class _FloatAndIndex:
    def __float__(self):
        return 1.0

    def __index__(self):
        return 2


# The physical-number boundaries of the one rule for error classes: an accepted type with a value
# outside the domain is a ValueError, a type the position never takes a TypeError, whatever the
# position (a value to convert, a resolution).
PHYSICAL_BOUNDARIES = {
    "large Python int": (lambda: ngc.to_units(2**64, 1.0), ValueError),
    "large Python int in list": (lambda: ngc.to_units([2**64], 1.0), ValueError),
    "large negative Python int in list": (
        lambda: ngc.to_units([1, -(2**70)], 1.0),
        ValueError,
    ),
    "large int beside a float": (lambda: ngc.to_units([2**64, 1.5], 1.0), ValueError),
    "float beside a large negative int": (
        lambda: ngc.to_units([1.5, -(2**70)], 1.0),
        ValueError,
    ),
    "large int beside a NumPy float": (
        lambda: ngc.to_units([2**70, np.float32(1.5)], 1.0),
        ValueError,
    ),
    "large object-dtype ndarray": (
        lambda: ngc.to_units(np.array([2**64], dtype=object), 1.0),
        TypeError,
    ),
    "large int beside None": (lambda: ngc.to_units([2**64, None], 1.0), TypeError),
    "mixed Python bool": (lambda: ngc.to_units([True, 1.0], 1.0), TypeError),
    "mixed NumPy bool": (lambda: ngc.to_units([1, np.True_], 1.0), TypeError),
    "Graph Fraction resolution": (
        lambda: ngc.Graph(2, [0], [1], [1], [1], Fraction(1)),
        TypeError,
    ),
    "to_units Fraction resolution": (
        lambda: ngc.to_units(1.0, Fraction(1)),
        TypeError,
    ),
    "as_float Fraction resolution": (lambda: ngc.as_float(1, Fraction(1)), TypeError),
    "to_units Decimal resolution": (
        lambda: ngc.to_units(1.0, Decimal(1)),
        TypeError,
    ),
    "Graph Decimal resolution": (
        lambda: ngc.Graph(2, [0], [1], [1], [1], Decimal(1)),
        TypeError,
    ),
    "as_float Decimal resolution": (lambda: ngc.as_float(1, Decimal(1)), TypeError),
    "Graph __float__ resolution": (
        lambda: ngc.Graph(2, [0], [1], [1], [1], _FloatOnly()),
        TypeError,
    ),
    "to_units __float__ resolution": (
        lambda: ngc.to_units(2.0, _FloatOnly()),
        TypeError,
    ),
    "as_float __float__ resolution": (
        lambda: ngc.as_float(1, _FloatOnly()),
        TypeError,
    ),
    "Graph __float__ and __index__ resolution": (
        lambda: ngc.Graph(2, [0], [1], [1], [1], _FloatAndIndex()),
        TypeError,
    ),
    "to_units __float__ and __index__ resolution": (
        lambda: ngc.to_units(2.0, _FloatAndIndex()),
        TypeError,
    ),
    "as_float __float__ and __index__ resolution": (
        lambda: ngc.as_float(1, _FloatAndIndex()),
        TypeError,
    ),
    "Decimal entry": (lambda: ngc.to_units([Decimal(1), 1.0], 1.0), TypeError),
    "__index__ value": (lambda: ngc.to_units(_Index(1), 1.0), TypeError),
    "as_float Decimal entry": (lambda: ngc.as_float([Decimal(1)], 1.0), TypeError),
    "Fraction value": (lambda: ngc.to_units(Fraction(1, 2), 1.0), TypeError),
    "Graph huge integer resolution": (
        lambda: ngc.Graph(2, [0], [1], [1], [1], 2**1024),
        ValueError,
    ),
    "to_units huge integer resolution": (
        lambda: ngc.to_units(1.0, 2**1024),
        ValueError,
    ),
    "as_float huge integer resolution": (lambda: ngc.as_float(1, 2**1024), ValueError),
    "as_float huge negative integer resolution": (
        lambda: ngc.as_float(1, -(2**1024)),
        ValueError,
    ),
}


@pytest.mark.parametrize(
    ("call", "cls"), PHYSICAL_BOUNDARIES.values(), ids=PHYSICAL_BOUNDARIES.keys()
)
def test_physical_boundaries(call, cls):
    with pytest.raises((TypeError, ValueError)) as error:
        call()
    assert error.type is cls
    if cls is ValueError and "resolution" in str(error.value):
        assert str(error.value) == RESOLUTION


@pytest.mark.parametrize(
    "x", [[True, False], (np.True_, np.False_), [True]], ids=["list", "numpy", "one"]
)
def test_an_all_bool_sequence_keeps_its_dtype_text(x):
    # Refused before amendment 6 too: the text stays NumPy's dtype, the class TypeError.
    with pytest.raises(
        TypeError, match=exact("to_units takes physical numbers (got dtype bool)")
    ) as error:
        ngc.to_units(x, 1.0)
    assert error.type is TypeError


def test_physical_inputs_still_taken():
    # The accepted types at each position: NumPy numbers as resolutions, an exact int (2^63 fits
    # NumPy's uint64), floats beside ints.
    for q in (np.float32(0.5), np.float64(2.0), np.int64(4), 8):
        assert ngc.to_units(8.0, q) == int(8 / float(q))
        assert ngc.as_float(1, q) == float(q)
    assert ngc.to_units(2**40, 2.0**-12) == 2**52
    assert ngc.to_units([2**40, 1.5], 1.0).tolist() == [2**40, 2]
    assert ngc.Graph(2, [0], [1], [1], [1], np.float64(1.0)).resolution() == 1.0
    assert ngc.as_float(np.int64(3), 1.0) == 3.0
    assert ngc.as_float([np.int64(1), 2], 1.0).tolist() == [1.0, 2.0]
    # as_float's amount is an integer position: Python's integer protocol, as elsewhere.
    assert ngc.as_float(_Index(3), 1.0) == 3.0
    assert ngc.as_float([1, _Index(3)], 1.0).tolist() == [1.0, 3.0]
    g = ngc.Graph(_Index(2), [_Index(0)], [1], [_Index(5)], [1], 1.0)
    assert g.eligibility(excluded_nodes=[_Index(1)]).node_allowed().tolist() == [
        True,
        False,
    ]
    # A resolution's own conversion that fails passes through unchanged; an object known only by
    # __index__ is not a number type, refused before any conversion.
    with pytest.raises(ZeroDivisionError, match="own float"):
        ngc.to_units(1.0, _FloatFails(1.0))
    with pytest.raises(TypeError, match="incompatible function arguments"):
        ngc.to_units(1.0, _IndexFails())
