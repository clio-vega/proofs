/-
SCOPE ONLY (LEAN 2026-11-23 c2): what remains of `hokF` after G1–G4 (library `FwdInv.lean`).
Scratch, NOT in the library; every `sorry` is a named gap. The `#eval` sweep
(`FwdInvSweep.lean`, `g5row`) checks `RunInv` at the start of every slide of the hokF run:
310/310, `F ∪ U = B^∨` 310/310, slide stops within fuel 310/310; control (holes in FORWARD order
`U`): 156/310 pass, so the description is not vacuous.
-/
import TworowD4Kernel.FwdInv

namespace TworowD4Kernel.G5Scope

open TworowD4Kernel.SlideSwitch TworowD4Kernel.StdSlide TworowD4Kernel.PhiOut
open TworowD4Kernel.RRStd TworowD4Kernel.PhiGlue TworowD4Kernel.SkewColInj TworowD4Kernel.FwdInv

/-- SlidesOK-free description of the state before slide `i`: `T` = holes still to fill
(`U.reverse.drop i`). Support and holes are disjoint, their union is convex, no support cell lies
weakly south-west of a remaining hole, and the support is semistandard. -/
def RunInv (S : Finset Cell) (X : Cell → ℕ) (T : List Cell) : Prop :=
  (∀ t ∈ T, t ∉ S) ∧ Skew (S ∪ T.toFinset) ∧
  (∀ t ∈ T, ∀ z ∈ S, ¬ (z.1 ≤ t.1 ∧ z.2 ≤ t.2)) ∧ SSOn S X

/-- G5a: `RunInv` + the head hole maximal among the remaining holes ⇒ `FwdInv` at the start
(`HoleOK` is vacuous: the west/south neighbours of `h` are not in `S`). -/
theorem fwdInv_of_runInv (S : Finset Cell) (X : Cell → ℕ) (h : Cell) (T : List Cell)
    (hr : RunInv S X (h :: T)) (hmax : ∀ t ∈ T, ¬ (h.1 ≤ t.1 ∧ h.2 ≤ t.2)) :
    FwdInv S X h := sorry

/-- G5b: fuel. A forward slide from a `FwdInv` state stops within `S.card + 1` moves (each move
raises `h.1 + h.2` by one inside the finite `insert h S`). -/
theorem slide_stops (k : ℕ) (S : Finset Cell) (X : Cell → ℕ) (h : Cell) (hi : FwdInv S X h)
    (hk : S.card + 1 ≤ k) :
    fwdChoice (slide fwdRule k S X h).1 (slide fwdRule k S X h).2.1
      (slide fwdRule k S X h).2.2.1 = none := sorry

/-- G5c: `RunInv` passes to the next slide (uses `fwd_end`, `fwdInv_slide`, `slide_stops`,
`insert y (moveS S h y) = insert h S`: the final hole is a maximal cell of `S ∪ T`). -/
theorem runInv_next (k : ℕ) (S : Finset Cell) (X : Cell → ℕ) (h : Cell) (T : List Cell)
    (hr : RunInv S X (h :: T)) (hmax : ∀ t ∈ T, ¬ (h.1 ≤ t.1 ∧ h.2 ≤ t.2))
    (hk : S.card + 1 ≤ k) :
    RunInv (slide fwdRule k S X h).1 (slide fwdRule k S X h).2.1 T := sorry

/-- G5d: the start list `U.reverse` (U = cells of the partition `A^∨`, row-major) has each head
maximal among its tail. Pure partition combinatorics. -/
theorem rev_heads_maximal (R C : ℕ) (A : List ℕ) :
    ∀ i, ∀ h ∈ ((skewCells (compl R C A) []).reverse.drop i).head?,
      ∀ t ∈ ((skewCells (compl R C A) []).reverse.drop (i + 1)),
        ¬ (h.1 ≤ t.1 ∧ h.2 ≤ t.2) := sorry

/-- G5e = G7 (lem:rrstd(a), the SSOn half) + `F ∪ U = B^∨`: the initial `RunInv`. -/
theorem runInv_init (R C : ℕ) (A B : List ℕ) (f : Cell → ℕ)
    (hbox : InBox R C (skewCells A B).toFinset)
    (hA : ∀ i j, i ≤ j → A.getD j 0 ≤ A.getD i 0) (hB : ∀ i j, i ≤ j → B.getD j 0 ≤ B.getD i 0)
    (hBA : ∀ i, B.getD i 0 ≤ A.getD i 0)
    (hrow : ∀ Y ∈ skewCells A B, (Y.1, Y.2 + 1) ∈ skewCells A B → f Y ≤ f (Y.1, Y.2 + 1))
    (hup : ∀ Y ∈ skewCells A B, (Y.1 + 1, Y.2) ∈ skewCells A B → f Y < f (Y.1 + 1, Y.2)) :
    RunInv (holesWith tilde R C A B f).toFinset
      (rr R C (((skewCells A B).map f).foldr max 0) f) (skewCells (compl R C A) []).reverse :=
  sorry

/-- G5f: the support size is constant along the run (`|F|`), so fuel `|F| + 1` suffices for
every slide. -/
theorem card_slide (k : ℕ) (S : Finset Cell) (X : Cell → ℕ) (h : Cell) (hh : h ∉ S) :
    (slide fwdRule k S X h).1.card = S.card := sorry

/-- The composite: `hokF` from the G5 pieces, via the PROVED `slidesOK_of_starts`. -/
theorem hokF_of_ssyt (R C : ℕ) (A B : List ℕ) (f : Cell → ℕ)
    (hbox : InBox R C (skewCells A B).toFinset)
    (hA : ∀ i j, i ≤ j → A.getD j 0 ≤ A.getD i 0) (hB : ∀ i j, i ≤ j → B.getD j 0 ≤ B.getD i 0)
    (hBA : ∀ i, B.getD i 0 ≤ A.getD i 0)
    (hrow : ∀ Y ∈ skewCells A B, (Y.1, Y.2 + 1) ∈ skewCells A B → f Y ≤ f (Y.1, Y.2 + 1))
    (hup : ∀ Y ∈ skewCells A B, (Y.1 + 1, Y.2) ∈ skewCells A B → f Y < f (Y.1 + 1, Y.2)) :
    SlidesOK fwdRule ((holesWith tilde R C A B f).length + 1)
      (holesWith tilde R C A B f).toFinset
      (rr R C (((skewCells A B).map f).foldr max 0) f)
      (skewCells (compl R C A) []).reverse := by
  apply slidesOK_of_starts
  -- remaining: induction over the run, RunInv → StartsOK, using G5a–G5f
  sorry

end TworowD4Kernel.G5Scope
