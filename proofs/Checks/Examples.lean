module

public import NetGraphCore

/-! Executable examples of the formal model; these do not run the native implementation. -/

namespace NetGraphCore.Checks

public section

/-- Ten units give the first recipient the single leftover unit. -/
example : (List.range 3).map (Apportionment.roundRobin 10 3) = [4, 3, 3] := by decide

/-- Zero traffic, fewer units than recipients, and one complete circle are separate cases. -/
example : (List.range 3).map (Apportionment.roundRobin 0 3) = [0, 0, 0] ∧
    (List.range 3).map (Apportionment.roundRobin 2 3) = [1, 1, 0] ∧
    (List.range 3).map (Apportionment.roundRobin 3 3) = [1, 1, 1] := by decide

/-- The sixth unit continues the circle and goes to recipient 2, not recipient 0. -/
example : (List.range 3).map (Apportionment.increments 5 1 3) = [0, 0, 1] := by decide

/-- A nonexistent recipient is excluded: the direct formula would invent a full-circle share. -/
example : Apportionment.sequential 3 3 5 = 0 ∧
    Apportionment.roundRobin 3 3 5 = 1 := by decide

/-- Adding complete circles adds the same amount to each share. -/
example : (List.range 3).map (Apportionment.roundRobin (10 + 3 * 2) 3) =
    [6, 5, 5] := by decide

/-- A high word changes the remainder: 2^64 leaves remainder 1 on division by 3. -/
example : (List.range 3).map (Apportionment.increments (2^64) 2 3) = [0, 1, 1] := by decide

/-- This model-only addition crosses the public counter limit; the identity still holds. -/
example : (List.range 7).map (Apportionment.increments (2^127-1) 9 7) =
    [1, 2, 2, 1, 1, 1, 1] := by decide

/-- Inputs whose exact outputs are compared to independently recorded expected values. -/
def roundRobinCases : List (String × Nat × Nat) :=
  [("rr-ten", 10, 3), ("rr-zero", 0, 3), ("rr-short", 2, 3),
   ("rr-circle", 3, 3), ("rr-one", 7, 1)]

/-- Increment cases exercise circle boundaries and totals that do not fit 64 bits. -/
def incrementCases : List (String × Nat × Nat × Nat) :=
  [("inc-continue", 5, 1, 3), ("inc-zero", 5, 0, 3),
   ("inc-wrap", 2, 2, 3), ("inc-wide", 2^64, 2, 3),
   ("inc-limit", 2^127-1, 9, 7)]

/-- These actual offers exercise public ShortestPaths and Pinned calls beyond 2^64 in the bridge. -/
def publicOffers : List Nat := [5, 1, 2^62-1, 2^62-1, 2^62-1, 2^62-1, 2^62-1, 4, 1]

/-- Export model computations so the native bridge can compare exact values and inputs. -/
def caseJson : Lean.Json :=
  let rr := roundRobinCases.map fun (name, total, members) => Lean.Json.mkObj [
    ("id", Lean.toJson name), ("kind", Lean.toJson "rr"),
    ("total", Lean.toJson total), ("members", Lean.toJson members),
    ("expected", Lean.toJson ((List.range members).map
      (Apportionment.roundRobin total members)))]
  let inc := incrementCases.map fun (name, previous, extra, members) => Lean.Json.mkObj [
    ("id", Lean.toJson name), ("kind", Lean.toJson "increment"),
    ("previous", Lean.toJson previous), ("extra", Lean.toJson extra),
    ("members", Lean.toJson members),
    ("expected", Lean.toJson ((List.range members).map
      (Apportionment.increments previous extra members)))]
  let seq := publicOffers.foldl (fun (state : Nat × List Lean.Json) offer =>
    let total := state.1 + offer
    (total, state.2 ++ [Lean.Json.mkObj [
      ("offer", Lean.toJson offer), ("total", Lean.toJson total),
      ("expected", Lean.toJson ((List.range 3).map (Apportionment.roundRobin total 3)))]]))
    (0, [])
  Lean.Json.mkObj [("cases", Lean.toJson (rr ++ inc)), ("public_sequence", Lean.toJson seq.2)]

end
end NetGraphCore.Checks
