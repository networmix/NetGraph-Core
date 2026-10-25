# Round robin and traffic that arrives in several calls

Status: mathematical model proofs. Native comparisons cover stated finite cases; they do not
prove that the C++ implementation agrees for every machine input.

## The question and the quantities

Suppose three sources take turns sending whole units. The first unit belongs to source 0,
the second to source 1, the third to source 2, and the fourth starts the circle again at 0.
If traffic arrives in several calls, the next call must continue that circle.

Here `total` is the number of units distributed so far, `members` is the number of sources,
and `recipient` is a source's position in the fixed order, starting at zero. The same arithmetic
can describe ordered outgoing links or tunnel entries. One unit is the graph's resolution `q`;
all calculations here count units, rather than rounded physical rates. For example, five units
at `q = 0.25` are a physical quantity of `1.25`.

`total / members` means the number of complete circles. `total % members` is the number of
units left after those circles. All theorem statements require `members > 0`; statements about
one recipient also state when its index must be below `members`. Zero recipients have no valid
index: that case is handled by the calling operation's empty-roster rules, not by dividing by
zero or claiming that a positive offer has been distributed.

The Lean definitions live in `proofs/NetGraphCore/Apportionment/RoundRobin.lean`:

- `sequential` counts how often one recipient is selected by the circle, one unit at a time.
- `roundRobin` uses division once: complete circles, plus one for each early leftover recipient.
- `increments` computes the extra share from the old remainder and the new offer, using a
  signed subtraction. Its nonnegativity is proved, rather than assumed by using a natural type.
- `amountLimit` is `2^62`, the exclusive bound on one native traffic amount.

Examples are checked by Lean in `Checks/Examples.lean` and compared with independently recorded
exact inputs and outputs in `proofs/cases.json`. The sequence example is also used by the native
bridge on both `ShortestPaths` and replay through `Pinned`.

<a id="NGC-AP-001"></a>
## NGC-AP-001: the formula describes the actual circle

Claim (`sequential_eq_roundRobin`): for any number of units and any valid recipient, distributing
units one at a time around the circle gives exactly the count from the division formula.

For ten units and three recipients, units 1, 4, 7 and 10 go to source 0; units 2, 5 and 8 go to
source 1; units 3, 6 and 9 go to source 2. The counts are `[4, 3, 3]`. Division gives three complete
circles (three units each), followed by one unit for the first source. This is case `rr-ten`.

Boundary cases: zero units give `[0, 0, 0]` (`rr-zero`); two units give `[1, 1, 0]` (`rr-short`);
exactly three units give `[1, 1, 1]` (`rr-circle`). One recipient gets every unit (`rr-one`).

The valid-index restriction matters. With three recipients and three units, nonexistent
position 5 is never selected by the circle, so its count is zero. The division formula would
give it `3 / 3 = 1` unit. The restriction excludes such phantom recipients; this contrast is
also checked in `Checks/Examples.lean`.

<a id="roundRobin_succ"></a>
Proof steps:

1. With zero units, both methods give every recipient zero.
2. Assume both methods agree after some total. The next sequential recipient is the old remainder.
3. Helper `roundRobin_succ` shows that the division formula changes the same recipient only.
   If the leftover count has not reached the end, the quotient stays fixed and just one more
   recipient enters the leftover group. At the end, the quotient rises by one and the leftover
   group becomes empty. Earlier recipients exchange their old leftover unit for a full-circle
   unit; the last recipient gains one. Therefore both methods make the same next step.
4. Repeating that argument covers every nonnegative total. This reasoning technique is called
   induction: prove the start and prove that every correct step gives a correct next step.

Why order matters: `[3, 3, 4]` also sums to ten and is almost equal, but puts the leftover unit
on the last recipient. Conservation alone would accept it; the circle-equivalence claim rejects
that rule. Reordering recipients changes who owns each turn, so membership and order are fixed.

<a id="NGC-AP-002"></a>
## NGC-AP-002: adding traffic cannot take a unit away

Claim (`roundRobin_mono`): increasing the total by a nonnegative amount never decreases a valid
recipient's share, while the recipient set and order stay fixed.

After five units, shares are `[2, 2, 1]`; after six they are `[2, 2, 2]`. Source 2 gains one and
the others keep their counts (`inc-continue`). Adding zero changes nothing (`inc-zero`).

The proof starts with an addition of zero. Helper `roundRobin_succ` says each following unit
adds either zero or one to a recipient. Combining those steps proves the claim for any addition.
This is the arithmetic monotonicity needed for later no-drop load arguments. It says nothing
about a roster change: five units over two members give `[3, 2]`, while five over three give
`[2, 2, 1]`, reducing the first member's share. It also does not by itself prove admission search.

<a id="NGC-AP-003"></a>
## NGC-AP-003: all input units are accounted for

Claim (`roundRobin_conservation`): adding the shares of all valid recipients gives the original
total exactly. The `Finset.range members` in Lean is just the finite list of valid positions
`0` through `members - 1`; the sum does not include any other indices.

For `rr-ten`, `4 + 3 + 3 = 10`; for `rr-short`, `1 + 1 + 0 = 2`; for `rr-zero`, all terms are zero.
The proof begins with the zero sum. At each new unit, `roundRobin_succ` describes every share's
change. Exactly one valid index equals the old remainder, so the sum rises by exactly one.

The imported finite-sum facts distribute a sum over additions and count the one index selected
by an equality. They do not impose a traffic model. This result rules out dropped or duplicated
units in the model, but permits a wrong allocation such as `[10, 0, 0]` if used alone; use it with
NGC-AP-001. A positive offer with no recipients is outside the theorem's stated domain.

<a id="NGC-AP-004"></a>
## NGC-AP-004: continue from the previous call using only its remainder

Claim (`increments_eq_difference`): the signed increment formula equals the new cumulative share
minus the previous cumulative share. The previous total can be arbitrarily large in this model.

With three sources, five units leave counts `[2, 2, 1]`. The previous remainder is two. Adding one
makes `2 + 1 = 3`: one complete circle and no leftovers. For each recipient the increment is
one full-circle unit, minus whether that recipient already had a leftover unit. That gives
`[1-1, 1-1, 1-0] = [0, 0, 1]`, so the final counts are `[2, 2, 2]` (`inc-continue`).

If the addition is zero, the old leftover flags are subtracted from themselves and every change
is zero (`inc-zero`). Starting from two units and adding two crosses a circle boundary: counts
change from `[1, 1, 0]` to `[2, 1, 1]`, an increment of `[1, 0, 1]` (`inc-wrap`).

The same argument works beyond 64 bits. `2^64 = 18446744073709551616` leaves remainder one when
divided by three. Adding two therefore assigns one to source 1 and one to source 2, yielding
`[0, 1, 1]` (`inc-wide`). Near the public cumulative limit, `2^127-1` leaves remainder one on
division by seven. Adding nine gives one new full circle and leftovers on the first three
positions; subtracting the old leftover on position 0 gives `[1, 2, 2, 1, 1, 1, 1]` (`inc-limit`).
That addition crosses the public cumulative-counter limit: it illustrates the unrestricted
mathematical identity, and is not an example of a public call the library must accept.

Proof steps:

1. Split the previous total into complete circles and a remainder. The imported division fact
   `Nat.mod_add_div` is the exact identity "remainder + members * quotient = original total".
2. <a id="roundRobin_periodic"></a>Helper `roundRobin_periodic` shows that adding complete circles increases every count equally.
   Division changes the quotient; its remainder and therefore leftover positions stay the same.
   For example, `16 = 10 + 3 * 2` units give `[6, 5, 5]`: the ten-unit shares `[4, 3, 3]`
   plus two complete circles, adding two to each source. This is checked in `Checks/Examples.lean`.
3. The complete circles of the old total occur in both cumulative counts and cancel in their
   difference. Only the old remainder, the extra offer and the two leftover flags remain.
4. That remaining expression is precisely `increments`, including its signed subtraction.

A wrong rule that restarts each call assigns the sixth unit to source 0, producing `[3, 2, 1]`.
The total still adds up, so this example exposes a mistake that conservation cannot detect.
The identity applies to equal round robin. It does not justify replacing weighted allocation
with equal allocation or forgetting the high word when computing a large total's remainder.

<a id="NGC-AP-005"></a>
## NGC-AP-005: an increment is a real part of the new offer

Claim (`increments_bounds`): a valid recipient's increment is between zero and the extra offer,
inclusive. This gives meaning to converting the signed expression back into a traffic amount.

Adding two units after total two gives `[1, 0, 1]` (`inc-wrap`): each entry is nonnegative and no
larger than two. Adding zero gives exactly zero (`inc-zero`). With one recipient an addition
belongs entirely to it, reaching the upper bound.

NGC-AP-004 expresses the increment as a difference. NGC-AP-002 bounds that difference below by
zero. <a id="roundRobin_add_le"></a>Helper `roundRobin_add_le` bounds it above: every sequential step can give this recipient
at most one, so an addition of `extra` cannot increase its share by more than `extra`. Combining
those facts yields both inequalities. This proves the output bound; it does not yet establish
the behaviour of a particular C++ subtraction instruction.

For example, the `inc-wrap` addition raises source 0's count from one to two. Its gain is one,
which is no more than the two added units. Source 1 gains zero and source 2 gains one; neither
can gain more than the whole offer. With a single source, its gain equals the offer.

<a id="NGC-AP-006"></a>
## NGC-AP-006: the remainder calculation has room in the native amount domain

Claim (`native_increment_bounds`): if there are fewer than `2^31` recipients and an offer smaller
than `2^62`, then "old remainder + offer" is below `2^63`. The resulting signed increment is
nonnegative and below `2^62`. The recipient must be valid; the count must be positive.

For a source split, the count bound follows from the graph's signed 32-bit node IDs: it has
fewer than `2^31` possible source members. Other uses of the arithmetic must establish their
own count bound; this theorem does not infer it from the word "recipient".

For `inc-wrap`, the remainder is two, the offer is two and their sum is four. The increment
`[1, 0, 1]` fits the claimed bounds. For `inc-limit`, the previous total is huge but its remainder
is only one. Thus the intermediate for offer nine is ten: the size of the previous total does
not make this addition large. At the upper input boundary, offer `2^62-1` and remainder at most
`2^31-2` still sum to less than `2^63`. An offer equal to `2^62` is outside the native domain,
even when this particular arithmetic calculation would fit; that is a contract restriction.

The proof uses the division rule that a remainder is smaller than its positive divisor. Add
the stated upper bounds on recipient count and offer to bound the intermediate. NGC-AP-005
then bounds each final increment by the offer. Numeric normalization computes the fixed powers
of two inside a checked proof; it does not approximate them with floating point.

The model does not prove the implementation of `mod_u128`, the representation of `Counter128`,
pointer accesses or exception safety. The native bridge exercises cumulative offers past `2^64`
through actual public calls, and records that finite coverage explicitly. Values near `2^127`
are model examples here, not a claim that the public bridge reaches that total in a few calls.

## How these results help the next stages

These proofs establish a precise arithmetic building block for source injection and ordered
equal splits. To prove a forwarding bound, we still need fixed-program load definitions and
their composition rules. To cover WCMP, we need separate highest-averages proofs. To certify
native behaviour for every input, we need a further refinement argument relating the model to
the actual program and its machine representation.
