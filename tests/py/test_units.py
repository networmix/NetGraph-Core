"""to_units (round half up, the only entry of physical numbers) and as_float."""

from __future__ import annotations

import math
import random
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
    # The detached array kernel reuses the resolution checked before entry.
    # These hand-derived values distinguish half-up from floor and ties-to-even.
    values = [
        2.5 * q,
        0.5 * q,
        0.49 * q,
        7.5 * q,
        math.nextafter(7.5 * q, 0),
        0.0,
        -0.0,
    ]
    assert ngc.to_units(values, q).tolist() == [3, 1, 0, 8, 7, 0, 0]
    assert ngc.as_float(
        [0, 1, (1 << 53) + 1, (1 << 53) + 3, -((1 << 53) + 3)], q
    ).tolist() == [
        0.0,
        q,
        (1 << 53) * q,
        ((1 << 53) + 4) * q,
        -((1 << 53) + 4) * q,
    ]


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


# Exact binary values, the last half-unit tie, promotion, round trips and signed products.
# Expectations are fixed arithmetic witnesses, independent of the production converter.
UNITS_EXAMPLES = [
    pytest.param(
        [
            {"op": "to_units", "x": 5501195.1796, "q": ["pow2", -21]},
            {"op": "error", "x": 5501195.1796, "q": ["pow2", -21]},
            {
                "op": "to_units",
                "x": math.nextafter(5501195.1796, 0.0),
                "q": ["pow2", -21],
            },
        ],
        [
            11536842473289,
            dict(
                units=11536842473289,
                error="1/2",
                inside=True,
                error_against_argument="1/2",
            ),
            11536842473288,
        ],
        id="UNITS-literal-tie",
    ),
    pytest.param(
        [
            {"op": "to_units", "x": 400 + 2 ** (-19), "q": ["pow2", -18]},
            {"op": "error", "x": 400 + 2 ** (-19), "q": ["pow2", -18]},
            {"op": "to_units", "x": 400.0, "q": ["pow2", -18]},
        ],
        [
            104857601,
            dict(
                units=104857601, error="1/2", inside=True, error_against_argument="1/2"
            ),
            104857600,
        ],
        id="UNITS-computed-tie",
    ),
    pytest.param(
        [
            {"op": "to_units", "x": ["numpy", "float32", 0.1], "q": ["pow2", -32]},
            {
                "op": "to_units",
                "x": ["array", "float32", [0.1, 0.1]],
                "q": ["pow2", -32],
            },
            {"op": "error", "x": ["numpy", "float32", 0.1], "q": ["pow2", -32]},
            {"op": "to_units", "x": ["numpy", "float64", 0.1], "q": ["pow2", -32]},
            {"op": "to_units", "x": 0.1, "q": ["pow2", -32]},
        ],
        [
            429496736,
            [429496736, 429496736],
            dict(units=429496736, error="0", inside=True, error_against_argument="0"),
            429496730,
            429496730,
        ],
        id="UNITS-float32",
    ),
    pytest.param(
        [
            {"op": "to_units", "x": 2 ** (-33), "q": ["pow2", -32]},
            {
                "op": "to_units",
                "x": math.nextafter(2 ** (-33), 0.0),
                "q": ["pow2", -32],
            },
            {"op": "to_units", "x": 2.0**31, "q": ["pow2", 32]},
            {"op": "to_units", "x": 0.5, "q": 1.0},
            {"op": "to_units", "x": 1.5, "q": 1.0},
            {"op": "to_units", "x": 2.5, "q": 1.0},
            {"op": "to_units", "x": 3.5, "q": 1.0},
            {"op": "to_units", "x": 2.0**52 - 0.5, "q": 1.0},
            {"op": "to_units", "x": 2.0**52 + 1, "q": 1.0},
            {"op": "to_units", "x": 2.0**53 + 2, "q": 1.0},
            {"op": "to_units", "x": -0.0, "q": 1.0},
            {"op": "to_units", "x": 5e-324, "q": ["pow2", 32]},
        ],
        [1, 0, 1, 1, 2, 3, 4, 2**52, 2**52 + 1, 2**53 + 2, 0, 0],
        id="UNITS-ties",
    ),
    pytest.param(
        [
            {"op": "round_trip", "u": 2**53, "q": 1.0},
            {"op": "round_trip", "u": 2**53, "q": ["pow2", -32]},
            {"op": "round_trip", "u": 2**53, "q": ["pow2", 32]},
            {"op": "round_trip", "u": 2**53 + 1, "q": 1.0},
            {"op": "as_float", "u": 2**53 + 1, "q": 1.0},
            {"op": "round_trip", "u": 2**53 + 2, "q": 1.0},
            {"op": "round_trip", "u": 2**62 - 2**9, "q": ["pow2", -32]},
            {"op": "round_trip", "u": 2**62 - 1, "q": 1.0},
            {"op": "round_trip", "u": 120401795907070192, "q": ["pow2", -27]},
        ],
        [
            2**53,
            2**53,
            2**53,
            2**53,
            9007199254740992.0,
            2**53 + 2,
            2**62 - 2**9,
            "ValueError: to_units: the amount is at or above the bound of the units domain (2^62 units)",
            120401795907070192,
        ],
        id="UNITS-round-trip",
    ),
    pytest.param(
        [
            {"op": "error", "x": 2**53 + 1, "q": 1.0},
            {"op": "error", "x": 2**53 + 3, "q": 1.0},
            {"op": "to_units", "x": ["list", 2**53 + 1, 2**53 + 3], "q": 1.0},
        ],
        [
            dict(units=2**53, error="0", inside=True, error_against_argument="-1"),
            dict(units=2**53 + 4, error="0", inside=True, error_against_argument="1"),
            [2**53, 2**53 + 4],
        ],
        id="UNITS-promoted-integer",
    ),
    pytest.param(
        [
            {"op": "as_float", "u": -1, "q": 0.5},
            {"op": "as_float", "u": 2**53 + 3, "q": 1.0},
            {"op": "as_float", "u": 2**70 + 1, "q": 1.0},
            {"op": "as_float", "u": 2**1024 - 2**970 - 1, "q": 1.0},
            {"op": "as_float", "u": 2**1024 - 2**970, "q": 1.0},
            {"op": "as_float", "u": -(2**1024 - 2**970 - 1), "q": 1.0},
            {"op": "as_float", "u": (2**1024 - 2**970) * 2**32 - 1, "q": ["pow2", -32]},
            {"op": "as_float", "u": (2**1024 - 2**970) * 2**32, "q": ["pow2", -32]},
        ],
        [
            -0.5,
            9007199254740996.0,
            1.1805916207174113e21,
            1.7976931348623157e308,
            f"ValueError: as_float: {2**1024 - 2**970} units times the resolution exceed the largest finite double",
            -1.7976931348623157e308,
            1.7976931348623157e308,
            f"ValueError: as_float: {(2**1024 - 2**970) * 2**32} units times the resolution exceed the largest finite double",
        ],
        id="UNITS-as-float-range",
    ),
    pytest.param(
        [
            {"op": "sum", "xs": [0.5, 0.5], "q": 1.0},
            {"op": "sum", "xs": [0.25, 0.25], "q": 1.0},
        ],
        [
            dict(units=[1, 1], sum_of_units=2, total="1", units_of_total=1),
            dict(units=[0, 0], sum_of_units=0, total="1/2", units_of_total=1),
        ],
        id="UNITS-sums",
    ),
]


def _units_input(spec):
    if not isinstance(spec, list):
        return spec
    if spec[0] == "pow2":
        return 2.0 ** spec[1]
    if spec[0] == "numpy":
        return getattr(np, spec[1])(spec[2])
    if spec[0] == "array":
        return np.array(spec[2], dtype=spec[1])
    assert spec[0] == "list"
    return spec[1:]


def _units_to_python(value):
    if isinstance(value, np.ndarray):
        assert value.dtype == np.int64
        assert value.ndim == 1
        return value.tolist()
    assert type(value) is int
    return value


def _units_example(cases):
    results = []
    for case in cases:
        q = _units_input(case["q"])
        op = case["op"]
        try:
            if op == "to_units":
                result = _units_to_python(ngc.to_units(_units_input(case["x"]), q))
            elif op == "as_float":
                result = ngc.as_float(case["u"], q)
                assert type(result) is float
            elif op == "round_trip":
                result = ngc.to_units(ngc.as_float(case["u"], q), q)
                assert type(result) is int
            elif op == "error":
                a = _units_input(case["x"])
                b = float(a)
                u = ngc.to_units(a, q)
                value = Fraction(ngc.as_float(u, q))
                error = (value - Fraction(b)) / Fraction(q)
                exact_argument = Fraction(a) if isinstance(a, int) else Fraction(b)
                result = dict(
                    units=u,
                    error=str(error),
                    inside=-Fraction(1, 2) < error <= Fraction(1, 2),
                    error_against_argument=str((value - exact_argument) / Fraction(q)),
                )
            else:
                assert op == "sum"
                xs = case["xs"]
                each = _units_to_python(ngc.to_units(xs, q))
                total = sum(map(Fraction, xs), Fraction(0))
                result = dict(
                    units=each,
                    sum_of_units=sum(each),
                    total=str(total / Fraction(q)),
                    units_of_total=ngc.to_units(float(total), q),
                )
        except ValueError as exc:
            result = f"ValueError: {exc}"
        results.append(result)
    return results


def _assert_units_result(got, expected):
    # Equality alone lets 1, True and 1.0 pass interchangeably; pin the types and every key.
    assert type(got) is type(expected), (got, expected)
    if isinstance(expected, dict):
        assert got.keys() == expected.keys()
        for key in expected:
            _assert_units_result(got[key], expected[key])
    elif isinstance(expected, list):
        assert len(got) == len(expected)
        for value, want in zip(got, expected, strict=True):
            _assert_units_result(value, want)
    else:
        assert got == expected


@pytest.mark.parametrize("cases, expected", UNITS_EXAMPLES)
def test_binary64_examples(cases, expected):
    _assert_units_result(_units_example(cases), expected)


def _exact_binary_units(b, q):
    # An independent rational oracle: never calls a production conversion helper.
    return math.floor(Fraction(b) / Fraction(q) + Fraction(1, 2))


def _units_sweep():
    rng = random.Random(10)
    counts = dict.fromkeys(
        [
            "checked",
            "agree",
            "ties",
            "ties_up",
            "inside",
            "decimal_differs",
            "round_trips",
            "returned",
        ],
        0,
    )

    def one(b, q, tie_n=None):
        expected = _exact_binary_units(b, q)
        u = ngc.to_units(b, q)
        assert type(u) is int
        assert u == expected, (b, q, u, expected)
        # The array kernel must obey the same rule, including both ndarray ranks.
        assert ngc.to_units(np.array(b), q) == expected
        array = ngc.to_units([b, b], q)
        assert array.dtype == np.int64 and array.tolist() == [expected, expected]
        array = ngc.to_units(np.array([b, 0.0, b])[::2], q)
        assert array.dtype == np.int64 and array.tolist() == [expected, expected]
        counts["checked"] += 1
        counts["agree"] += u == expected
        if tie_n is not None:
            assert expected == tie_n + 1
            counts["ties"] += 1
            counts["ties_up"] += u == tie_n + 1
        value = ngc.as_float(u, q)
        assert type(value) is float
        assert Fraction(value) == Fraction(u) * Fraction(q)
        assert ngc.as_float([u, u], q).tolist() == [value, value]
        e = Fraction(value) - Fraction(b)
        assert -Fraction(q) / 2 < e <= Fraction(q) / 2
        counts["inside"] += 1
        counts["decimal_differs"] += (
            math.floor(Fraction(str(b)) / Fraction(q) + Fraction(1, 2)) != expected
        )

    for k in range(-32, 33):
        q = 2.0**k
        ns = [0, 1, 2**26, 2**52 - 1] + [rng.randrange(2**52) for _ in range(64)]
        for n in ns:
            tie = (2 * n + 1) * q / 2
            one(tie, q, n)
            one(math.nextafter(tie, 0.0), q)
            one(math.nextafter(tie, math.inf), q)
        for edge in (2**52, 2**53):
            one(edge * q, q)
            one(math.nextafter(edge * q, 0.0), q)
        one(5e-324, q)
        for _ in range(64):
            one(math.ldexp(rng.random(), rng.randrange(-40, 63)) * q, q)
        for u in [0, 1, 2**53 - 1, 2**53] + [
            rng.randrange(2**53 + 1) for _ in range(64)
        ]:
            value = ngc.as_float(u, q)
            assert Fraction(value) == Fraction(u) * Fraction(q)
            returned = ngc.to_units(value, q)
            assert returned == u
            counts["round_trips"] += 1
            counts["returned"] += returned == u
    return counts


def test_binary64_double_steps():
    _assert_units_result(
        _units_sweep(),
        dict(
            checked=17745,
            agree=17745,
            ties=4420,
            ties_up=4420,
            inside=17745,
            decimal_differs=2248,
            round_trips=4420,
            returned=4420,
        ),
    )


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_binary64_float32_containers(dtype):
    # float32 promotion is exact; the shortest decimal of the source float32 is not its value.
    value = dtype(0.1)
    expected = 429496736 if dtype is np.float32 else 429496730
    q = 2.0**-32
    assert ngc.to_units(value, q) == expected
    assert ngc.to_units(float(value), q) == expected
    assert ngc.to_units(np.array(value), q) == expected
    for values in (
        [value, value],
        (value, value),
        np.array([value, value], dtype=dtype),
    ):
        result = ngc.to_units(values, q)
        assert result.dtype == np.int64
        assert result.tolist() == [expected, expected]


@pytest.mark.parametrize("q", [2.0**k for k in range(-32, 33)])
def test_binary64_signed_round_trip_and_refusal(q):
    # Both directions around the consecutive-integer limit; vector and scalar paths.
    for u in [
        -(2**53 + 3),
        -(2**53 + 1),
        -(2**53),
        -1,
        0,
        1,
        2**53,
        2**53 + 1,
        2**53 + 3,
    ]:
        expected = float(Fraction(u) * Fraction(q))
        assert ngc.as_float(u, q) == expected
        assert ngc.as_float([u, u], q).tolist() == [expected, expected]
    with pytest.raises(ValueError, match=exact(BOUND)):
        ngc.to_units(ngc.as_float(2**62 - 1, q), q)
    assert ngc.to_units(ngc.as_float(2**62 - 512, q), q) == 2**62 - 512


@pytest.mark.parametrize("q", [1.0, 2.0**-32])
@pytest.mark.parametrize("sign", [-1, 1])
def test_binary64_unbounded_finite_edge(q, sign):
    first_invalid = (2**1024 - 2**970) * (1 if q == 1.0 else 2**32)
    last_valid = sign * (first_invalid - 1)
    expected = sign * float.fromhex("0x1.fffffffffffffp+1023")
    assert ngc.as_float(last_valid, q) == expected
    result = ngc.as_float([1, last_valid], q)
    assert result.dtype == np.float64 and result.tolist() == [q, expected]
    u = sign * first_invalid
    message = exact(
        f"as_float: {u} units times the resolution exceed the largest finite double"
    )
    for value in (u, [1, u]):
        with pytest.raises(ValueError, match=message) as error:
            ngc.as_float(value, q)
        assert error.type is ValueError


class _IntegerSequence(int):
    def __len__(self):
        return 2

    def __getitem__(self, index):
        if index in (0, 1):
            return (3, 4)[index]
        raise IndexError(index)


class _NumpySignedSequence(np.int64):
    __len__ = _IntegerSequence.__len__
    __getitem__ = _IntegerSequence.__getitem__


class _NumpyUnsignedSequence(np.uint64):
    __len__ = _IntegerSequence.__len__
    __getitem__ = _IntegerSequence.__getitem__


@pytest.mark.parametrize(
    "scalar_type", [_IntegerSequence, _NumpySignedSequence, _NumpyUnsignedSequence]
)
def test_as_float_integer_scalar_precedes_sequence(scalar_type):
    value = scalar_type(7)
    result = ngc.as_float(value, 0.5)
    assert type(result) is float
    assert result == 3.5
    assert result.hex() == (3.5).hex()


def test_as_float_noninteger_sequence_keeps_its_path():
    class IndexSequence(list):
        def __index__(self):
            return 7

    result = ngc.as_float(IndexSequence([3, 4]), 0.5)
    assert type(result) is np.ndarray
    assert result.dtype == np.float64 and result.shape == (2,)
    assert result.tolist() == [1.5, 2.0]


@pytest.mark.parametrize("dtype", [np.int64, np.uint64])
def test_as_float_integer_array_paths_are_unchanged(dtype):
    with pytest.raises(
        TypeError,
        match=exact("as_float: x must be a 1-D array (got 0 dimensions)"),
    ):
        ngc.as_float(np.array(7, dtype=dtype), 0.5)
    result = ngc.as_float(np.array([3, 4], dtype=dtype), 0.5)
    assert type(result) is np.ndarray
    assert result.dtype == np.float64 and result.shape == (2,)
    assert result.tolist() == [1.5, 2.0]
    with pytest.raises(
        TypeError,
        match=exact("as_float: x must be a 1-D array (got 2 dimensions)"),
    ):
        ngc.as_float(np.array([[7]], dtype=dtype), 0.5)


@pytest.mark.parametrize("value", [True, False, np.True_, np.False_])
def test_as_float_scalar_classification_still_rejects_bool(value):
    with pytest.raises(TypeError):
        ngc.as_float(value, 0.5)
    assert ngc.as_float(1, 0.5) == 0.5


@pytest.mark.parametrize("dtype", [np.int64, np.uint64])
@pytest.mark.parametrize(
    "value, expected", [(2**53 + 1, 2**53), (2**53 + 3, 2**53 + 4)]
)
def test_to_units_integer_zero_dimensional_array_promotes(dtype, value, expected):
    result = ngc.to_units(np.array(value, dtype=dtype), 1.0)
    assert type(result) is int
    assert result == expected


@pytest.mark.parametrize("dtype", [np.int64, np.uint64])
def test_to_units_integer_array_promotes(dtype):
    values = np.array([2**53 + 1, 2**53 + 3], dtype=dtype)
    result = ngc.to_units(values, 1.0)
    assert type(result) is np.ndarray
    assert result.shape == (2,)
    assert result.dtype == np.int64
    assert result.tolist() == [2**53, 2**53 + 4]
