module

public import Mathlib.Data.Nat.ModEq
public import Mathlib.Algebra.BigOperators.Group.Finset.Basic
public import Mathlib.Algebra.BigOperators.Group.Finset.Piecewise
public import Mathlib.Tactic.NormNum
import Lean.Elab.Tactic.Omega

/-!+# Whole-unit round robin

This model has ordered recipients numbered from zero. Giving the next unit to
`total % members` means continuing the same circle, rather than restarting it.
The closed formula is introduced independently and then proved to agree with
that sequence. See `docs/formal/theorems/round-robin.md` for the worked argument.
All arithmetic here is unbounded natural arithmetic; native bounds are separate.
-/

namespace NetGraphCore.Apportionment

public section

/-- The share computed in one division: a full-round share and a possible leftover unit.
Only indices below a positive `members` count describe a recipient. -/
@[expose] def roundRobin (total members recipient : Nat) : Nat :=
  total / members + if recipient < total % members then 1 else 0

/-- Count units given to one recipient by handing them out in a circle, one at a time. -/
@[expose] def sequential (total members recipient : Nat) : Nat :=
  match total with
  | 0 => 0
  | n + 1 => sequential n members recipient + if n % members = recipient then 1 else 0

/-- Increasing the total by one changes precisely the share of the next recipient.
When the remainder wraps, a complete round replaces every previous leftover unit. -/
theorem roundRobin_succ (total members recipient : Nat)
    (hm : 0 < members) (hi : recipient < members) :
    roundRobin (total + 1) members recipient = roundRobin total members recipient +
      if total % members = recipient then 1 else 0 := by
  have hr := Nat.mod_lt total hm
  -- A remainder advances by one, or returns to zero after the final recipient.
  have hs : (total + 1) % members =
      if total % members + 1 = members then 0 else total % members + 1 := by
    rw [← Nat.add_mod_eq_add_mod_right 1 (Nat.mod_mod total members)]
    by_cases h : total % members + 1 = members
    · simp [h]
    · rw [Nat.mod_eq_of_lt (by omega)]
      simp [h]
  -- The full-round count changes only on that wrap. The remaining comparisons
  -- distinguish the chosen recipient, earlier recipients and later recipients.
  by_cases hw : total % members + 1 = members
  · simp [roundRobin, Nat.succ_div, Nat.dvd_iff_mod_eq_zero, hs, hw]
    split_ifs <;> omega
  · simp [roundRobin, Nat.succ_div, Nat.dvd_iff_mod_eq_zero, hs, hw]
    split_ifs <;> omega

/-- NGC-AP-001: the direct formula gives exactly the counts of the unit-by-unit circle.
The induction follows one real assignment, not an assumed conservation equation. -/
theorem sequential_eq_roundRobin (total members recipient : Nat)
    (hm : 0 < members) (hi : recipient < members) :
    sequential total members recipient = roundRobin total members recipient := by
  induction total with
  | zero => simp [sequential, roundRobin]
  | succ n ih =>
      -- The sequential process and formula add one to the same recipient.
      simpa [sequential, ih] using (roundRobin_succ n members recipient hm hi).symm

/-- NGC-AP-002: adding traffic never takes a previously assigned unit from a recipient. -/
theorem roundRobin_mono (total extra members recipient : Nat)
    (hm : 0 < members) (hi : recipient < members) :
    roundRobin total members recipient ≤ roundRobin (total + extra) members recipient := by
  induction extra with
  | zero => simp
  | succ v ih =>
      -- Each next step adds either zero or one to this recipient's count.
      have h := roundRobin_succ (total + v) members recipient hm hi
      rw [show total + (v + 1) = total + v + 1 by omega, h]
      omega

/-- NGC-AP-003: the ordered shares account for every input unit exactly once. -/
theorem roundRobin_conservation (total members : Nat) (hm : 0 < members) :
    ∑ i ∈ Finset.range members, roundRobin total members i = total := by
  induction total with
  | zero => simp [roundRobin]
  | succ n ih =>
      -- Sum the one-step equation over every valid recipient.
      have h : (∑ i ∈ Finset.range members, roundRobin (n + 1) members i) =
          ∑ i ∈ Finset.range members,
            (roundRobin n members i + if n % members = i then 1 else 0) := by
        apply Finset.sum_congr rfl
        intro i hi
        exact roundRobin_succ n members i hm (Finset.mem_range.mp hi)
      rw [h, Finset.sum_add_distrib, ih]
      -- The remainder names exactly one member of this finite set.
      simp [Finset.mem_range, Nat.mod_lt n hm]

/-- Adding complete circles increases every share by the number of circles. -/
theorem roundRobin_periodic (total circles members recipient : Nat) (hm : 0 < members) :
    roundRobin (total + members * circles) members recipient =
      roundRobin total members recipient + circles := by
  -- Complete circles change the quotient, leaving the leftover recipients unchanged.
  simp only [roundRobin, Nat.add_mul_div_left _ _ hm, Nat.add_mul_mod_self_left]
  omega

/-- The signed arithmetic used for an increment, using only the previous remainder.
Using integers keeps the subtraction honest; nonnegativity is proved below. -/
@[expose] def increments (previous extra members recipient : Nat) : Int :=
  ((previous % members + extra) / members : Nat) +
    (if recipient < (previous % members + extra) % members then 1 else 0) -
    (if recipient < previous % members then 1 else 0)

/-- NGC-AP-004: continuing a circle is computed from the old remainder alone.
This is the mathematical identity used by Core's round_robin_increments. -/
theorem increments_eq_difference (previous extra members recipient : Nat) (hm : 0 < members) :
    increments previous extra members recipient =
      (roundRobin (previous + extra) members recipient : Int) -
      (roundRobin previous members recipient : Int) := by
  -- Write the old total as complete circles plus a remainder; the circles cancel.
  have hx := Nat.mod_add_div previous members
  have h : previous + extra = previous % members + extra + members * (previous / members) := by
    omega
  conv_rhs => lhs; rw [h, roundRobin_periodic _ _ _ _ hm]
  unfold increments roundRobin
  simp only [Nat.cast_add]
  split_ifs <;> omega

/-- A recipient cannot receive more than the number of units added in this call. -/
theorem roundRobin_add_le (total extra members recipient : Nat)
    (hm : 0 < members) (hi : recipient < members) :
    roundRobin (total + extra) members recipient ≤ roundRobin total members recipient + extra := by
  induction extra with
  | zero => simp
  | succ v ih =>
      -- Each new unit raises this recipient's count by at most one.
      have h := roundRobin_succ (total + v) members recipient hm hi
      rw [show total + (v + 1) = total + v + 1 by omega, h]
      split_ifs <;> omega

/-- NGC-AP-005: the signed increment is nonnegative and no larger than the added traffic. -/
theorem increments_bounds (previous extra members recipient : Nat)
    (hm : 0 < members) (hi : recipient < members) :
    0 ≤ increments previous extra members recipient ∧
      increments previous extra members recipient ≤ (extra : Int) := by
  -- Translate the difference identity using the lower and upper bounds on a share's growth.
  rw [increments_eq_difference _ _ _ _ hm]
  have lower := roundRobin_mono previous extra members recipient hm hi
  have upper := roundRobin_add_le previous extra members recipient hm hi
  omega

/-- The exclusive bound on one native traffic amount, in integer units. -/
@[expose] def amountLimit : Nat := 2 ^ 62

/-- NGC-AP-006: the remainder-plus-offer intermediate fits unsigned 64-bit arithmetic,
and the final signed increment fits Amount. This does not verify the C++ implementation
of 128-bit remainder; that connection is explicitly tested, not formally refined. -/
theorem native_increment_bounds (previous extra members recipient : Nat)
    (hm : 0 < members) (hi : recipient < members)
    (hk : members < 2 ^ 31) (hv : extra < amountLimit) :
    previous % members + extra < 2 ^ 63 ∧
      0 ≤ increments previous extra members recipient ∧
      increments previous extra members recipient < (amountLimit : Int) := by
  -- The remainder is smaller than the recipient count. Combining the two stated
  -- input bounds leaves room for their sum, independently of the old total's size.
  have hr := Nat.mod_lt previous hm
  have hb := increments_bounds previous extra members recipient hm hi
  norm_num [amountLimit] at hv ⊢
  norm_num at hk
  omega

end
end NetGraphCore.Apportionment
