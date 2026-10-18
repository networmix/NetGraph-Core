"""The apportionment through the private hooks (``_round_robin``, ``_highest_averages``;
docs/semantics.md#1-units-the-one-numeric-domain)."""

from __future__ import annotations

import _netgraph_core as ext
import pytest
from helpers import exact
from private_hooks import needs_hooks


@needs_hooks
def test_apportionment_entries():
    assert ext._round_robin(7, 3) == [3, 2, 2]
    assert ext._highest_averages(10, [3, 9, 9]) == [1, 5, 4]
    assert ext._highest_averages(11, [3, 9, 9]) == [1, 5, 5]
    assert ext._highest_averages(4, [1, 3, 3]) == [0, 2, 2]
    big = (1 << 62) - 2
    a = ext._highest_averages(big, [(1 << 60) + 1, (1 << 60) - 1, 7])
    # By hand: the lower quotas floor(big * w / W), W = 2^61 + 7, are 2^61 - 6, 2^61 - 10 and
    # 13; the one unit left goes to the largest w / (a + 1), the second member's.
    assert a == [(1 << 61) - 6, (1 << 61) - 9, 13]
    assert sum(a) == big
    with pytest.raises(ValueError, match=exact("weights[0] = 0 is outside [1, 2^62)")):
        ext._highest_averages(1, [0, 1])
    # Legal weights whose sum reaches 2^64 are summed exactly. Eight of 2^61: one each. Two of
    # 2^61 sharing 3: the lower quotas are 1 and 1; the unit left goes to the tie's first member
    # (as 9/5 above).
    assert ext._highest_averages(8, [1 << 61] * 8) == [1] * 8
    assert ext._highest_averages(3, [1 << 61, 1 << 61]) == [2, 1]
