"""The types a number position takes (the one rule for error classes, with the number types).

A non-integer number position (a resolution of Graph, to_units and as_float; to_units' value)
takes exactly an int (not a bool), a float, a NumPy integer or floating scalar, an int or float
subclass by its stored value (never through its own __int__, __float__ or __array__), and for
to_units an array of an integer or floating dtype or a sequence of those scalars; any other type is
a TypeError. An integer position (as_float's amount, an offer) takes any object with __index__,
never a bool, and an exception of that __index__ reaches the caller as raised."""

from __future__ import annotations

import array
import gc
import itertools
import math
from collections import UserList
from collections.abc import Sequence
from decimal import Decimal
from fractions import Fraction

import numpy as np
import pytest
from helpers import exact, sweep_python_allocation
from private_hooks import needs_hooks

import netgraph_core as ngc

RESOLUTION = "resolution must be a power of two in [2^-32, 2^32]"
OBJECT = "to_units takes physical numbers (got dtype object)"
BOOL = "to_units takes physical numbers (got dtype bool)"


class Index:
    def __init__(self, value):
        self.value = value

    def __index__(self):
        return self.value


class IndexFails:
    def __index__(self):
        raise ZeroDivisionError("own index")


class Both:
    def __float__(self):
        return 1.0

    def __index__(self):
        return 2


class FloatOnly:
    def __float__(self):
        return 1.0


class ArrayOnly:
    def __array__(self, dtype=None, copy=None):
        return np.array(2.0, dtype=dtype)


# Subclasses whose own conversions must not be called (each raises).
class IntValue(int):
    def __float__(self):
        raise AssertionError("int __float__ called")

    def __int__(self):
        raise AssertionError("int __int__ called")


class FloatValue(float):
    def __float__(self):
        raise AssertionError("float __float__ called")


class IntArrayValue(int):
    def __array__(self, dtype=None, copy=None):
        raise AssertionError("int __array__ called")


class FloatArrayValue(float):
    def __array__(self, dtype=None, copy=None):
        raise AssertionError("float __array__ called")


# Subclasses whose own conversions disagree with their stored value 2.
class IntDifferent(int):
    def __int__(self):
        return 9


class FloatDifferent(float):
    def __float__(self):
        return 9.0


# Numeric subclasses that also look like a sequence of one 9: a number is read as the number.
class IntSequence(int):
    def __len__(self):
        return 1

    def __getitem__(self, i):
        if i == 0:
            return 9.0
        raise IndexError(i)


class FloatSequence(float):
    def __len__(self):
        return 1

    def __getitem__(self, i):
        if i == 0:
            return 9.0
        raise IndexError(i)


class ColdResolution(np.float32):
    # A resolution whose own conversion empties CPython's free lists, so the conversion's next
    # Python allocations are real ones the hook can fail.
    def __float__(self):
        gc.collect()
        return 1.0


RESOLUTION_CALLS = {
    "Graph": (lambda q: ngc.Graph(2, [0], [1], [1], [1], q).resolution(), 1.0),
    "to_units": (lambda q: ngc.to_units(2.0, q), 2),
    "as_float": (lambda q: ngc.as_float(1, q), 1.0),
}
REFUSED_RESOLUTIONS = {
    "Decimal": Decimal(1),
    "Fraction": Fraction(1),
    "Index": Index(1),
    "Both": Both(),
    "FloatOnly": FloatOnly(),
    "array0": np.array(1.0),
    "bool": True,
    "numpy-bool": np.True_,
}
TAKEN_RESOLUTIONS = {
    "int": 1,
    "float": 1.0,
    "numpy-int": np.int64(1),
    "numpy-float": np.float32(1.0),
    "int-subclass": IntValue(1),
    "float-subclass": FloatValue(1.0),
}


@pytest.mark.parametrize("name", RESOLUTION_CALLS)
@pytest.mark.parametrize("q", REFUSED_RESOLUTIONS.values(), ids=REFUSED_RESOLUTIONS)
def test_a_resolution_of_another_type_is_refused(name, q):
    with pytest.raises(TypeError, match="incompatible") as error:
        RESOLUTION_CALLS[name][0](q)
    assert error.type is TypeError


@pytest.mark.parametrize("name", RESOLUTION_CALLS)
@pytest.mark.parametrize("q", [2**1024, -(2**1024)], ids=["positive", "negative"])
def test_a_resolution_int_beyond_the_float_range_is_a_value(name, q):
    with pytest.raises(ValueError, match=exact(RESOLUTION)) as error:
        RESOLUTION_CALLS[name][0](q)
    assert error.type is ValueError


@pytest.mark.parametrize("name", RESOLUTION_CALLS)
@pytest.mark.parametrize("q", TAKEN_RESOLUTIONS.values(), ids=TAKEN_RESOLUTIONS)
def test_a_resolution_of_a_number_type_is_taken(name, q):
    call, want = RESOLUTION_CALLS[name]
    assert call(q) == want


@pytest.mark.parametrize(
    "x",
    [[2**64, 1.5], [1.5, -(2**70)], [2**70, np.float32(1.5)]],
    ids=["large-float", "negative-large-float", "large-numpy-float"],
)
def test_ints_beyond_64_bits_beside_floats_are_values(x):
    with pytest.raises(ValueError, match=exact(OBJECT)) as error:
        ngc.to_units(x, 1.0)
    assert error.type is ValueError


@pytest.mark.parametrize(
    "x", [[True, False], (np.True_, np.False_)], ids=["bool-list", "bool-tuple"]
)
def test_bools_alone_keep_their_dtype_text(x):
    with pytest.raises(TypeError, match=exact(BOOL)):
        ngc.to_units(x, 1.0)


@pytest.mark.parametrize(
    ("x", "want"),
    [
        (Index(3), 3.0),
        (Index(2**63), float(2**63)),
        ([1, Index(3)], [1.0, 3.0]),
        (Both(), 2.0),
    ],
    ids=["Index", "Index-large", "Index-list", "Both"],
)
def test_as_float_takes_the_integer_protocol(x, want):
    got = ngc.as_float(x, 1.0)
    assert (got.tolist() if isinstance(got, np.ndarray) else got) == want


@pytest.mark.parametrize("x", [IndexFails(), [1, IndexFails()]], ids=["scalar", "list"])
def test_as_float_passes_an_own_index_failure_through(x):
    with pytest.raises(ZeroDivisionError, match="own index"):
        ngc.as_float(x, 1.0)


def test_as_float_refuses_an_object_array():
    with pytest.raises(
        TypeError,
        match=exact("as_float: x entries must be integers (got dtype object)"),
    ):
        ngc.as_float(np.array([Index(3)], dtype=object), 1.0)


@pytest.mark.parametrize(
    "x",
    [
        Decimal(1),
        Fraction(1),
        Index(1),
        Both(),
        FloatOnly(),
        ArrayOnly(),
        [ArrayOnly()],
        [np.array(1.0), 2.0],
    ],
    ids=[
        "Decimal",
        "Fraction",
        "Index",
        "Both",
        "FloatOnly",
        "array-protocol",
        "array-protocol-entry",
        "array0-entry",
    ],
)
def test_to_units_refuses_other_types(x):
    with pytest.raises(TypeError) as error:
        ngc.to_units(x, 1.0)
    assert error.type is TypeError


def test_to_units_takes_a_real_zero_dimensional_array():
    assert ngc.to_units(np.array(2.0), 1.0) == 2


@pytest.mark.parametrize(
    "x",
    [
        IntValue(2),
        FloatValue(2.0),
        IntArrayValue(2),
        FloatArrayValue(2.0),
        IntDifferent(2),
        FloatDifferent(2.0),
    ],
    ids=[
        "int-subclass",
        "float-subclass",
        "int-array-protocol",
        "float-array-protocol",
        "int-subclass-other-int",
        "float-subclass-other-float",
    ],
)
def test_to_units_reads_a_subclass_by_its_value(x):
    assert ngc.to_units(x, 1.0) == 2
    assert ngc.to_units([x], 1.0).tolist() == [2]


def test_an_offer_keeps_the_integer_protocol():
    g = ngc.Graph(2, [0], [1], [5], [1], 1.0)
    p = ngc.FlowPolicy(
        ngc.FlowLedger(g), 0, 1, ngc.FlowPolicyConfig(routing=ngc.ShortestPaths())
    )
    assert p.place(Index(3)).call_delivered == 3


@pytest.mark.parametrize("cls", [IntSequence, FloatSequence])
def test_a_numeric_subclass_with_sequence_methods_is_a_number(cls):
    got = ngc.to_units(cls(2), 1.0)
    assert got == 2 and type(got) is int
    assert ngc.to_units([cls(2)], 1.0).tolist() == [2]
    assert ngc.to_units([9.0], 1.0).tolist() == [9]


def outcome(x, q=1.0):
    """to_units(x, q) as comparable data: its units, or its refusal (class and text). A
    MemoryError is not caught."""
    try:
        got = ngc.to_units(x, q)
    except (TypeError, ValueError) as e:
        return (type(e).__name__, str(e))
    return got.tolist() if isinstance(got, np.ndarray) else got


# Every branch of the converter: numbers of each type, subclasses, NumPy scalars, arrays of each
# rank, lists and tuples, ints beyond 64 bits among numbers, bools, other types, nested lists.
CONVERTER_INPUTS = [
    2.5,
    2**40,
    IntSequence(2),
    FloatSequence(2.0),
    np.int64(3),
    np.uint64(2**63),
    np.float32(1.5),
    np.array(2.5),
    np.array([1.5, 2.5], dtype=np.float32),
    np.array([1, 2], dtype=np.int64),
    [1, 2.5, np.int64(3), np.float32(0.5), IntSequence(2)],
    (np.float64(1.0), 2),
    [],
    [2**64, 1.5],
    2**64,
    [True, False],
    [True, 1.5],
    True,
    "1",
    [Decimal(1), 1.5],
    Decimal(1),
    ArrayOnly(),
    [ArrayOnly()],
    [[1, 2], [3, 4]],
    np.array([[1.0]]),
    np.array([True]),
    np.array([1], dtype=object),
    range(3),
    UserList([1, 2.5]),
    UserList(["a"]),
    memoryview(np.array([1.0, 2.0])),
    [2**64, True],
]


@needs_hooks
def test_the_converter_raises_memory_error_on_every_python_allocation_failure():
    # Each input's conversion with each of its Python allocations failing in turn (the free lists
    # emptied first): MemoryError exactly, then the same outcome as without failure.
    for x in CONVERTER_INPUTS:
        sweep_python_allocation(lambda x=x: (outcome, (x,), lambda: None))


@needs_hooks
@pytest.mark.filterwarnings("ignore::pytest.PytestUnraisableExceptionWarning")
def test_a_cold_conversion_raises_memory_error():
    # The review's witness: a resolution whose own conversion empties the free lists, then each
    # Python allocation of to_units([1.0, 2.0], q) failing in turn. A failure inside the
    # collection itself is recovered there (the call completes); every other one is a MemoryError.
    import _netgraph_core as ext

    x, q = [1.0, 2.0], ColdResolution(1.0)
    assert ngc.to_units(x, q).tolist() == [1, 2]
    failures = 0
    for at in range(400):
        hit, _, got, error = ext._fail_python_allocation(ngc.to_units, (x, q), at)
        assert ngc.to_units(x, q).tolist() == [1, 2]
        if not hit:
            assert error is None and got.tolist() == [1, 2]
            break
        if got is not None:
            assert error is None and got.tolist() == [1, 2]
            continue
        failures += 1
        assert type(error) is MemoryError, (at, repr(error))
    else:
        pytest.fail("no completing call in 400 allocations")
    assert failures >= 1


@pytest.mark.parametrize(
    "x", [[2**64, Decimal(1)], [Decimal(1), 2**64]], ids=["wide-first", "other-first"]
)
def test_an_int_beyond_64_bits_beside_another_type_is_a_wrong_type(x):
    # NumPy's object dtype from another class is a TypeError, wherever the wide int stands.
    with pytest.raises(TypeError, match=exact(OBJECT)) as error:
        ngc.to_units(x, 1.0)
    assert error.type is TypeError


class CustomSequence(Sequence):
    def __len__(self):
        return 3

    def __getitem__(self, i):
        return [0, 1, 2][i]


@pytest.mark.parametrize(
    "x",
    [
        range(3),
        memoryview(array.array("q", [0, 1, 2])),
        array.array("q", [0, 1, 2]),
        UserList([0, 1, 2]),
        CustomSequence(),
        [0, 1, 2],
        (0, 1, 2),
    ],
    ids=["range", "memoryview", "array.array", "UserList", "custom", "list", "tuple"],
)
def test_every_python_sequence_of_numbers_is_taken(x):
    assert ngc.to_units(x, 1.0).tolist() == [0, 1, 2]


def test_a_sequence_of_other_types_keeps_numpys_text():
    # Already refused before: the text NumPy's reading gives, whatever the sequence.
    for x in (["a"], UserList(["a"]), ("a",)):
        with pytest.raises(
            TypeError, match=exact("to_units takes physical numbers (got dtype <U1)")
        ):
            ngc.to_units(x, 1.0)
    with pytest.raises(
        TypeError, match=exact("to_units takes a float or a 1-D array of floats")
    ):
        ngc.to_units(memoryview(np.zeros((2, 2))), 1.0)


@pytest.mark.parametrize(
    "x",
    [[2**64, True], [False, 2**64], [np.True_, 2**64], [-(2**64), True]],
    ids=["wide-bool", "bool-wide", "numpy-bool-wide", "negative-wide-bool"],
)
def test_numpys_dtype_comes_before_the_bool_rule(x):
    # NumPy reads an int beyond 64 bits beside a bool as an object: its text, a wrong type.
    with pytest.raises(TypeError, match=exact(OBJECT)) as error:
        ngc.to_units(x, 1.0)
    assert error.type is TypeError


class ClaimedInteger:
    @property
    def __class__(self):
        return np.int64

    def __index__(self):
        return 2


class ClaimedFloating:
    @property
    def __class__(self):
        return np.float32

    def __float__(self):
        return 2.0


@pytest.mark.parametrize("cls", [ClaimedInteger, ClaimedFloating])
def test_a_claimed_numpy_class_is_not_a_number(cls):
    # Classification reads the object's real type, not its __class__.
    x = cls()
    assert isinstance(x, np.generic) and not issubclass(cls, np.generic)
    for value in (x, [x]):
        with pytest.raises(TypeError):
            ngc.to_units(value, 1.0)
    for call in (
        lambda: ngc.to_units(1.0, x),
        lambda: ngc.Graph(2, [0], [1], [1], [1], x),
    ):
        with pytest.raises(TypeError, match="incompatible"):
            call()
    assert ngc.to_units(np.int64(2), 1.0) == 2
    assert ngc.to_units(np.float32(2), 1.0) == 2


NOT_PHYSICAL = "to_units takes physical numbers (a float or an array of floats)"
RANK = "to_units takes a float or a 1-D array of floats"


def contract_outcome(x, q):
    """The outcome 14.3 gives to_units(x, q), computed without the library: the resolution; a
    bool, str or bytes; NumPy's dtype (an object from accepted numbers alone a value error, any
    other non-numeric dtype a wrong type); the rank; a bool among numbers; then each value."""

    def number(v):
        return not isinstance(v, (bool, np.bool_)) and isinstance(
            v, (int, float, np.integer, np.floating)
        )

    if q == 0.75:
        return "ValueError", RESOLUTION
    if isinstance(x, (bool, np.bool_, str, bytes)):
        return "TypeError", NOT_PHYSICAL
    raw = np.asarray(x)
    if raw.dtype.kind not in "fiu":
        accepted = (
            all(number(v) for v in x) if isinstance(x, (list, tuple)) else number(x)
        )
        cause = "ValueError" if raw.dtype.kind == "O" and accepted else "TypeError"
        return cause, f"to_units takes physical numbers (got dtype {raw.dtype})"
    if raw.ndim > 1:
        return "TypeError", RANK
    if isinstance(x, (list, tuple)) and any(isinstance(v, (bool, np.bool_)) for v in x):
        return "TypeError", NOT_PHYSICAL
    units = []
    for v in raw.astype(np.float64).reshape(-1):
        v = float(v)
        if math.isnan(v):
            return "ValueError", "to_units: NaN has no unit value"
        if math.isinf(v):
            return (
                "ValueError",
                "to_units: an infinite amount is not finite and has no unit value",
            )
        if v < 0:
            return "ValueError", "to_units: a negative amount has no unit value"
        y = v / q
        rounded = math.floor(y) + (y - math.floor(y) >= 0.5)
        if rounded >= 2**62:
            return (
                "ValueError",
                "to_units: the amount is at or above the bound of the units domain (2^62 units)",
            )
        units.append(rounded)
    return "valid", units[0] if raw.ndim == 0 else units


def test_the_check_order_over_pairs_and_triples():
    # The review's order matrix: every pair of a pool of scalars and every triple of a smaller
    # one, as lists and tuples, at three resolutions, against the contract's order.
    pool = [0, 2.5, -1.0, math.nan, math.inf, 2**62, 2**64, True, np.True_]
    pool += [Decimal(1), "1", None]
    triples = [2.5, -1.0, 2**64, True, Decimal(1)]
    inputs = list(pool)
    for size, values in ((2, pool), (3, triples)):
        for row in itertools.product(values, repeat=size):
            inputs += [list(row), tuple(row)]
    inputs += [
        np.array([0, 1, 2]),
        np.array([0.5, 2.5], dtype=np.float32),
        np.array([True]),
        np.array([1], dtype=object),
        np.array([[1, 2]]),
        np.array([[True]]),
    ]
    bad = []
    for q in (1.0, 0.5, 0.75):
        for x in inputs:
            try:
                got = ngc.to_units(x, q)
            except (TypeError, ValueError) as e:
                actual = (type(e).__name__, str(e))
            else:
                actual = ("valid", got.tolist() if isinstance(got, np.ndarray) else got)
            if actual != contract_outcome(x, q):
                bad.append((x, q, actual))
    assert bad == []
