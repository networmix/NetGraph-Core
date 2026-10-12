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
    with pytest.raises(ValueError, match=exact("the weights sum to 2^62 or more")):
        ext._highest_averages(1, [1 << 61, 1 << 61])
