/- SCRATCH, NOT IN THE LIBRARY. Assembly signature for `thm:G4` with a `have` per link.
   Every `sorry` here is a gap that remains; the file exists to enumerate them
   ([[a-gap-list-written-before-the-assembly-undercounts]]). Not registered. -/
import TworowD4Kernel.RowsSlidesBridge
import TworowD4Kernel.Destd
import TworowD4Kernel.ChainInj
import TworowD4Kernel.RankList
open TworowD4Kernel.SlideSwitch TworowD4Kernel.StdSlide TworowD4Kernel.StdSlideBridge
open TworowD4Kernel.SlideIteration TworowD4Kernel.SlidesColsBridge TworowD4Kernel.PhiOut
open TworowD4Kernel.PhiGlue TworowD4Kernel.Destd TworowD4Kernel.RowsSlidesBridge
open TworowD4Kernel.ChainInj TworowD4Kernel.RankList


theorem thm_G4_of_slideOK (R C : ℕ) (A B : List ℕ) (f : Cell → ℕ)
    -- the paper's standing hypotheses on `f` (LR tableau on `A/B` in the box) enter only
    -- through the links below; stated here as the named hypotheses those links need
    -- standing hypotheses of `def:phi`(1): A/B in the box, f column-strict (semistandard)
    (hbox : InBox R C (skewCells A B).toFinset) (hcol : ColInj (skewCells A B).toFinset f)
    (hstd : StdFrom ∅ (skewCells (compl R C A) [] ++ holesWith tilde R C A B f))
    (hokR : SlidesOK revRule (R + C) (skewCells (compl R C A) []).toFinset uVal
      (holesWith tilde R C A B f))
    (hokF : SlidesOK fwdRule ((holesWith tilde R C A B f).length + 1)
      (holesWith tilde R C A B f).toFinset
      (rr R C (((skewCells A B).map f).foldr max 0) f)
      (skewCells (compl R C A) []).reverse) :
    (phi R C A B f).1 = (sigmaRR R C A B f).1 ∧
    ∀ a ∈ (phi R C A B f).1, (phi R C A B f).2 a = (sigmaRR R C A B f).2 a := by
  set U := skewCells (compl R C A) []
  set F := holesWith tilde R C A B f
  set X := rr R C (((skewCells A B).map f).foldr max 0) f
  -- L0 (proved): the word has no repeated cell
  have L0 : (U ++ F).Nodup := (stdFrom_nodup ∅ _ hstd).1
  -- L1 (proved 11-22, ChainInj.rowsW_eq_colsW)
  have L1 : rowsW F.length U.length (U ++ F) = colsW U.length 0 F.length (U ++ F) :=
    rowsW_eq_colsW U.length F.length (U ++ F) (by simp) hstd
  -- L2 (gap, glue only): Φ = (U'.toFinset, gU U ∘ valOf U'), U' = (colsW …).drop p
  --   from phi_values_destd + slides_colsW_take + slides_fuel_eq (fuel R+C ↔ |U|+1)
  have L2 : (phi R C A B f).1 = ((colsW U.length 0 F.length (U ++ F)).drop F.length).toFinset ∧
      ∀ a ∈ (phi R C A B f).1, (phi R C A B f).2 a =
        gU U (valOf ((colsW U.length 0 F.length (U ++ F)).drop F.length) a) := by
    have := phi_values_destd R C A B f hokR
    sorry
  -- L3 (proved 11-22, RankList.L3_holes_toFinset)
  have L3 : F.toFinset = ((skewCells A B).map (rot R C)).toFinset := L3_holes_toFinset R C A B f hbox hcol
  -- L4 (proved 11-22, RankList.L4_std_eq_valOf)
  have L4 : ∀ a ∈ F, std F.toFinset X a = valOf F a := L4_std_eq_valOf R C A B f hbox hcol
  -- L5 (proved given L3, L4): Σ's vacated cells are the window form
  have L5 : (slides fwdRule (((skewCells A B).map (rot R C)).toFinset.card + (compl R C A).sum + 1)
      ((skewCells A B).map (rot R C)).toFinset X (sigmaHoles (compl R C A))).2.2 =
      ((rowsW F.length U.length (U ++ F)).drop F.length).reverse := by
    rw [sigmaHoles_eq, ← L3]
    exact V_rowsW_of_slidesOK_fuel _ U F X L0 L4 hokF
      (by have := List.toFinset_card_of_nodup (List.nodup_append.1 L0).2.1; omega)
  -- L6 (gap, list combinatorics): with W := rowsW … = colsW … = Fpart ++ V.reverse,
  --   Σ's support V.toFinset = U'.toFinset and Σ's label of V[k] (row of U[m-1-k] + 1)
  --   = gU U (valOf U' V[k])
  have L6 : True := trivial
  sorry
