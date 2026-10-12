"""House-monotone integer apportionment (APPROACH revision 24, section 2.1).

* ``round_robin(I, k)``: EqualBalanced and Lossy EqualBalanced. ``floor(I / k)``
  each, one more unit to the first ``I mod k`` members (members are listed in
  ascending EdgeId).
* ``highest_averages(I, w)``: Weighted and Lossy Weighted, the Jefferson / D'Hondt
  rule with EdgeId ties. Closed-form evaluator: every member starts at its lower
  quota ``floor(I * w_m / W)`` (exactly the quotients ``w_m / j >= W / I``), then the
  remaining ``I - sum`` units (fewer than ``k``) go one at a time by a heap of the
  exact next quotients ``w_m / (a_m + 1)``, ties by ascending EdgeId. Quotients are
  compared by cross-multiplication; every product is checked to stay below 2^127
  (the 128-bit claim of the contract).
* ``highest_averages_slow``: the sequential definition, for tests only.
* ``hamilton``: largest remainders, the rejected rule (review-13 D1), controls only.
"""

from __future__ import annotations

import heapq
from dataclasses import dataclass, field

PRODUCT_LIMIT = 1 << 127


@dataclass
class ApportionStats:
    """Work counters of the closed-form evaluator (diagnostics)."""

    calls: int = 0
    heap_pops: int = 0
    max_product_bits: int = 0
    pops_per_call: list[int] = field(default_factory=list)


STATS = ApportionStats()


def round_robin(total: int, k: int) -> list[int]:
    if k <= 0:
        raise ValueError("round robin over no member")
    if total < 0:
        raise ValueError("negative amount")
    base, extra = divmod(total, k)
    return [base + (1 if i < extra else 0) for i in range(k)]


def _cross(a: int, b: int) -> int:
    p = a * b
    if p >= PRODUCT_LIMIT:
        raise OverflowError("quotient comparison exceeds 128 bits")
    bl = p.bit_length()
    if bl > STATS.max_product_bits:
        STATS.max_product_bits = bl
    return p


class _Quotient:
    """Heap entry for ``w / d`` with ascending-EdgeId ties (index = member rank)."""

    __slots__ = ("w", "d", "i")

    def __init__(self, w: int, d: int, i: int) -> None:
        self.w, self.d, self.i = w, d, i

    def __lt__(self, other: _Quotient) -> bool:
        # Larger quotient first; equal quotients by lower member rank (EdgeId).
        left = _cross(self.w, other.d)
        right = _cross(other.w, self.d)
        if left != right:
            return left > right
        return self.i < other.i


def highest_averages(total: int, weights: list[int]) -> list[int]:
    k = len(weights)
    if k == 0:
        raise ValueError("highest averages over no member")
    if total < 0:
        raise ValueError("negative amount")
    if any(w <= 0 for w in weights):
        raise ValueError("weights must be positive")
    STATS.calls += 1
    wsum = sum(weights)
    alloc = [_cross(total, w) // wsum for w in weights]
    left = total - sum(alloc)
    assert 0 <= left < k, (total, weights, alloc)
    heap = [_Quotient(w, alloc[i] + 1, i) for i, w in enumerate(weights)]
    heapq.heapify(heap)
    pops = 0
    for _ in range(left):
        q = heapq.heappop(heap)
        pops += 1
        alloc[q.i] += 1
        heapq.heappush(heap, _Quotient(q.w, alloc[q.i] + 1, q.i))
    STATS.heap_pops += pops
    STATS.pops_per_call.append(pops)
    return alloc


def highest_averages_slow(total: int, weights: list[int]) -> list[int]:
    """Sequential D'Hondt: each unit to the member maximising ``w / (a + 1)``."""
    alloc = [0] * len(weights)
    for _ in range(total):
        best = 0
        for i in range(1, len(weights)):
            lhs = weights[i] * (alloc[best] + 1)
            rhs = weights[best] * (alloc[i] + 1)
            if lhs > rhs:
                best = i
        alloc[best] += 1
    return alloc


def hamilton(total: int, weights: list[int]) -> list[int]:
    """Largest remainders (rejected; not house-monotone: [1, 3, 3] at 3 and 4)."""
    wsum = sum(weights)
    alloc = [total * w // wsum for w in weights]
    order = sorted(range(len(weights)), key=lambda i: (-(total * weights[i] % wsum), i))
    for i in order[: total - sum(alloc)]:
        alloc[i] += 1
    return alloc


def split(rule: str, total: int, weights: list[int] | None, k: int) -> list[int]:
    """The node split of a kernel: ``"rr"`` or ``"ha"`` (weights in member order)."""
    if rule == "rr":
        return round_robin(total, k)
    if rule == "ha":
        assert weights is not None and len(weights) == k
        return highest_averages(total, weights)
    raise ValueError(f"unknown split rule {rule!r}")
