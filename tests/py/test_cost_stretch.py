"""CostStretch as one exact fraction: the accepted forms and their conversion in
the documented order (docs/semantics.md, "Cost ceilings and the exact stretch"), the components at their limits, and the value semantics (factor, equality,
hash, repr). The refusal below 1 where the ceiling is used is in test_ksp.py and
test_distributed_te.py."""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest
from helpers import assert_read_only, exact

import netgraph_core as ngc

NUMERATOR = "CostStretch fraction numerator must be in [-2^63, 2^63)"
DENOMINATOR = "CostStretch fraction denominator must be in [1, 2^63)"


class _Index:
    # An integer known only through __index__.
    def __init__(self, n):
        self.n = n

    def __index__(self):
        return self.n


class _OwnFloatFails(np.float32):
    # A NumPy floating value whose own conversion to a float fails.
    def __float__(self):
        raise ArithmeticError("own float")


def reduced(factor):
    f = ngc.CostStretch(factor).factor
    return f"{f.numerator}/{f.denominator}"


def test_values():
    # STRETCH-values: a float is the decimal it prints (1.2 is 6/5, 0.1 + 0.2 its own 17 digits),
    # a Fraction as it is (2^64 / 2^62 reduces to 4), an integer n/1, an __index__ object too.
    for factor, want in (
        (1.2, "6/5"),
        (1.1, "11/10"),
        (1.5, "3/2"),
        (0.1 + 0.2, "7500000000000001/25000000000000000"),
        (Fraction(6, 5), "6/5"),
        (Fraction(2**64, 2**62), "4/1"),
        (2, "2/1"),
        (_Index(3), "3/1"),
    ):
        assert reduced(factor) == want, factor
    with pytest.raises(
        ValueError, match=exact("CostStretch factor must be finite (got inf)")
    ):
        ngc.CostStretch(float("inf"))
    with pytest.raises(ValueError, match=exact(NUMERATOR)):
        ngc.CostStretch(1e19)
    with pytest.raises(ValueError, match=exact(DENOMINATOR)):
        ngc.CostStretch(Fraction(1, 2**63))
    for bad in (True, False, "1.2", None, [1]):
        with pytest.raises(TypeError, match=r"^CostStretch factor must be"):
            ngc.CostStretch(bad)
    # Equal factors are equal values: 1.5 and 3/2, 1.2 and 6/5; the binary value of 1.2 is not.
    assert ngc.CostStretch(1.5) == ngc.CostStretch(Fraction(3, 2))
    assert ngc.CostStretch(1.2) == ngc.CostStretch(Fraction(6, 5))
    assert ngc.CostStretch(1.2) != ngc.CostStretch(
        Fraction(5404319552844595, 4503599627370496)
    )


def test_numpy():
    # STRETCH-numpy: every NumPy floating type prints its own shortest decimal (float32's 1.2 is
    # 6/5, not 1.2000000476837158), a NumPy integer is n/1, numpy.bool_ is no number.
    for factor, want in (
        (np.float64(1.2), "6/5"),
        (np.float32(1.2), "6/5"),
        (np.float16(1.2), "6/5"),
        (np.float32(0.30000001192092896), "3/10"),
        (np.int64(2), "2/1"),
    ):
        assert reduced(factor) == want, factor
    for factor, shown in ((np.float64("nan"), "nan"), (np.float32("inf"), "inf")):
        with pytest.raises(
            ValueError, match=exact(f"CostStretch factor must be finite (got {shown})")
        ):
            ngc.CostStretch(factor)
    with pytest.raises(ValueError, match=exact(NUMERATOR)):
        ngc.CostStretch(np.float64(1e19))
    with pytest.raises(TypeError, match=r"^CostStretch factor must be"):
        ngc.CostStretch(np.bool_(True))
    assert ngc.CostStretch(np.float32(1.2)) == ngc.CostStretch(Fraction(6, 5))
    assert ngc.CostStretch(np.float64(1.2)) == ngc.CostStretch(1.2)
    # The object's own conversion failing reaches the caller unchanged.
    with pytest.raises(ArithmeticError, match=exact("own float")):
        ngc.CostStretch(_OwnFloatFails(1.2))


def test_component_limits():
    # The reduced components are checked as exact integers, whatever their size before reduction.
    with pytest.raises(ValueError, match=exact(NUMERATOR)):
        ngc.CostStretch(Fraction(2**63, 1))
    with pytest.raises(ValueError, match=exact(NUMERATOR)):
        ngc.CostStretch(2**63)
    with pytest.raises(ValueError, match=exact(DENOMINATOR)):
        ngc.CostStretch(Fraction(1, 2**63))
    assert (
        reduced(Fraction(-(2**63), 1)) == f"{-(2**63)}/1"
    )  # valid; refused where used
    assert reduced(Fraction(2**63 - 1, 2**63 - 1)) == "1/1"
    assert reduced(Fraction(2**64, 2**62)) == "4/1"
    assert reduced(Fraction(2**63 - 1, 2**63 - 2)) == f"{2**63 - 1}/{2**63 - 2}"
    assert (
        reduced(0.5) == "1/2"
    )  # below 1: a valid value, refused where the ceiling is used


def test_value_semantics():
    # factor is a fractions.Fraction; equality, hash and repr follow the fraction.
    f = ngc.CostStretch(1.5).factor
    assert type(f) is Fraction and f == Fraction(3, 2)
    assert ngc.CostStretch(1.5) == ngc.CostStretch(Fraction(3, 2))
    assert hash(ngc.CostStretch(1.5)) == hash(ngc.CostStretch(Fraction(3, 2)))
    assert hash(ngc.CostStretch(1.5)) == hash(Fraction(3, 2))
    assert ngc.CostStretch(2) != ngc.AbsoluteCost(2)
    assert repr(ngc.CostStretch(1.2)) == "CostStretch(Fraction(6, 5))"
    assert repr(ngc.CostStretch(-(2**63))) == f"CostStretch(Fraction({-(2**63)}, 1))"
    assert_read_only(ngc.CostStretch(Fraction(6, 5)), "factor")


def test_ceiling_products_beyond_64_bits():
    # Two parallel links 0 -> 1 of capacity 1, costs B = 2^32 (the shortest) and C. The ceiling
    # compares C * den with num * B, products of about 2^64 to 2^65: with C = 2^32 + 1 and
    # f = 2^32 / (2^32 - 1), C * den = 2^64 - 1 <= 2^64 = num * B (admitted; the low 64 bits alone
    # say 2^64 - 1 > 0); with C = 2^33 and f = (2^32 + 2) / (2^32 + 1), C * den = 2^65 + 2^33 >
    # 2^64 + 2^33 (refused; the low 64 bits alone are equal). T1's second tunnel and ksp's second
    # path follow the exact comparison.
    b = 1 << 32
    for c, factor, admitted in (
        (b + 1, Fraction(b, b - 1), True),
        (2 * b, Fraction(b + 2, b + 1), False),
    ):
        g = ngc.Graph(2, [0, 0], [1, 1], [1, 1], [b, c], 1.0)
        ceiling = ngc.CostStretch(factor)
        paths = ngc.ksp(g, 0, 1, ngc.KspOptions(2, cost_ceiling=ceiling))
        assert [cost for cost, _ in paths] == ([b, c] if admitted else [b]), c
        t1 = ngc.FlowPolicyConfig(
            routing=ngc.DistributedTE(ngc.TunnelCount(2), cost_ceiling=ceiling)
        )
        r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, t1).replace(2, design_demand=2)
        status = ngc.TunnelStatus
        assert [t.status for t in r.program.tunnels] == [
            status.UP,
            status.UP if admitted else status.COST_CEILING_EXCEEDED,
        ], c
        assert r.call_delivered == (2 if admitted else 1), c
