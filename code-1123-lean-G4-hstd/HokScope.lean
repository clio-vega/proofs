/-
SCOPE ONLY (LEAN 2026-11-23, Task 3): what would discharge `hokR` / `hokF` of
`G4AssemblyStd.thm_G4_of_slideOK'`. Scratch file, NOT in the library; every `sorry` is a gap.
Paper fact: jdt slides preserve semistandardness (Fulton, Young Tableaux §1.2, RECALLED in the tex
`lem:A` proof; no arXiv locator). `SlidesOK` = every start hole is outside the current support and
every intermediate filling is column-injective, along the whole run.
-/
import TworowD4Kernel.G4AssemblyStd

namespace TworowD4Kernel.HokScope

open TworowD4Kernel.SlideSwitch TworowD4Kernel.StdSlide TworowD4Kernel.PhiOut
open TworowD4Kernel.RRStd TworowD4Kernel.PhiGlue

/-- Semistandard on a set of cells: weak along rows, strict up columns (adjacent cells only). -/
def SSOn (S : Finset Cell) (X : Cell → ℕ) : Prop :=
  (∀ a ∈ S, (a.1, a.2 + 1) ∈ S → X a ≤ X (a.1, a.2 + 1)) ∧
  (∀ a ∈ S, (a.1 + 1, a.2) ∈ S → X a < X (a.1 + 1, a.2))

/-- Convex in the product order (a skew shape). -/
def Skew (S : Finset Cell) : Prop :=
  ∀ a ∈ S, ∀ b ∈ S, ∀ c : Cell, a.1 ≤ c.1 → a.2 ≤ c.2 → c.1 ≤ b.1 → c.2 ≤ b.2 → c ∈ S

/-- Fulton's "semistandard with a hole": inequalities ACROSS the hole `h` (west ≤ east, south <
north), plus the two diagonal ones the move needs. Exact form is part of gap G1. -/
def HoleOK (S : Finset Cell) (X : Cell → ℕ) (h : Cell) : Prop := sorry

/-- The forward-slide invariant (gap G1). -/
def FwdInv (S : Finset Cell) (X : Cell → ℕ) (h : Cell) : Prop :=
  h ∉ S ∧ Skew (insert h S) ∧ SSOn S X ∧ HoleOK S X h

/-- G2 (Fulton 1.2 core, forward): one move keeps the invariant. Case split north/east move. -/
theorem fwd_step (S : Finset Cell) (X : Cell → ℕ) (h y : Cell) (hi : FwdInv S X h)
    (hc : fwdRule.ch S X h = some y) : FwdInv (moveS S h y) (moveX X h y) y := sorry

/-- G3: on a skew shape, semistandard ⇒ column-injective (columns are intervals). This also
makes `hcol` of `thm_G4_of_slideOK'` redundant. -/
theorem colInj_of_SS (S : Finset Cell) (X : Cell → ℕ) (hS : Skew S) (hX : SSOn S X) :
    ColInj S X := sorry

/-- G4: induction over fuel using G2 + G3 (needs `Skew S` from `Skew (insert h S)`). -/
theorem slideOK_fwd (k : ℕ) (S : Finset Cell) (X : Cell → ℕ) (h : Cell) (hi : FwdInv S X h) :
    SlideOK fwdRule k S X h := sorry

/-- G5 (the load-bearing one): after each slide of Σ's run the NEXT start hole is an inner corner
of the current support, with the invariant. Needs the support after `k` slides = the chain window
(existing: `slides_colsW_take` / `V_rowsW_of_slidesOK_fuel` give the window form only UNDER
SlidesOK, so it cannot be reused as is — circular). Needs a fuel-free description of the support. -/
theorem fwd_next_hole_ok : True := sorry

/-- G6/G7: mirror of G1–G4 for `revRule` (tie south), and lem:rrstd(a) as a standalone
`SSOn` statement for `rr f` on `rot(A/B)` (only the `std` half was formalised today). -/
theorem rev_mirror : True := sorry

/-- The composite. -/
theorem hokF_of_ssyt (R C : ℕ) (A B : List ℕ) (f : Cell → ℕ)
    (hbox : InBox R C (skewCells A B).toFinset) (hcol : ColInj (skewCells A B).toFinset f)
    (hA : ∀ i j, i ≤ j → A.getD j 0 ≤ A.getD i 0) (hB : ∀ i j, i ≤ j → B.getD j 0 ≤ B.getD i 0)
    (hBA : ∀ i, B.getD i 0 ≤ A.getD i 0)
    (hrow : ∀ Y ∈ skewCells A B, (Y.1, Y.2 + 1) ∈ skewCells A B → f Y ≤ f (Y.1, Y.2 + 1))
    (hup : ∀ Y ∈ skewCells A B, (Y.1 + 1, Y.2) ∈ skewCells A B → f Y < f (Y.1 + 1, Y.2)) :
    SlidesOK fwdRule ((holesWith tilde R C A B f).length + 1)
      (holesWith tilde R C A B f).toFinset
      (rr R C (((skewCells A B).map f).foldr max 0) f)
      (skewCells (compl R C A) []).reverse := sorry

theorem hokR_of_ssyt (R C : ℕ) (A B : List ℕ) (f : Cell → ℕ)
    (hbox : InBox R C (skewCells A B).toFinset) (hcol : ColInj (skewCells A B).toFinset f)
    (hA : ∀ i j, i ≤ j → A.getD j 0 ≤ A.getD i 0) (hB : ∀ i j, i ≤ j → B.getD j 0 ≤ B.getD i 0)
    (hBA : ∀ i, B.getD i 0 ≤ A.getD i 0)
    (hrow : ∀ Y ∈ skewCells A B, (Y.1, Y.2 + 1) ∈ skewCells A B → f Y ≤ f (Y.1, Y.2 + 1))
    (hup : ∀ Y ∈ skewCells A B, (Y.1 + 1, Y.2) ∈ skewCells A B → f Y < f (Y.1 + 1, Y.2)) :
    SlidesOK revRule (R + C) (skewCells (compl R C A) []).toFinset uVal
      (holesWith tilde R C A B f) := sorry

/-- `thm:G4` with NO named slide hypotheses, modulo the gaps above. -/
theorem thm_G4 (R C : ℕ) (A B : List ℕ) (f : Cell → ℕ)
    (hbox : InBox R C (skewCells A B).toFinset) (hcol : ColInj (skewCells A B).toFinset f)
    (hA : ∀ i j, i ≤ j → A.getD j 0 ≤ A.getD i 0) (hB : ∀ i j, i ≤ j → B.getD j 0 ≤ B.getD i 0)
    (hBA : ∀ i, B.getD i 0 ≤ A.getD i 0)
    (hrow : ∀ Y ∈ skewCells A B, (Y.1, Y.2 + 1) ∈ skewCells A B → f Y ≤ f (Y.1, Y.2 + 1))
    (hup : ∀ Y ∈ skewCells A B, (Y.1 + 1, Y.2) ∈ skewCells A B → f Y < f (Y.1 + 1, Y.2)) :
    (phi R C A B f).1 = (sigmaRR R C A B f).1 ∧
    ∀ a ∈ (phi R C A B f).1, (phi R C A B f).2 a = (sigmaRR R C A B f).2 a :=
  G4AssemblyStd.thm_G4_of_slideOK' R C A B f hbox hcol hA hB hBA hrow hup
    (hokR_of_ssyt R C A B f hbox hcol hA hB hBA hrow hup)
    (hokF_of_ssyt R C A B f hbox hcol hA hB hBA hrow hup)

end TworowD4Kernel.HokScope
