/-! ## Control (kernel): every hypothesis of `thm_G4_of_slideOK` holds on the 310 instances -/

instance decInBox (R C : ℕ) (S : Finset Cell) : Decidable (InBox R C S) := by
  unfold InBox; infer_instance

instance decColInj (S : Finset Cell) (X : Cell → ℕ) : Decidable (ColInj S X) := by
  unfold ColInj; infer_instance

/-- All five named hypotheses of `thm_G4_of_slideOK` on one data row. -/
def hypCheckG4 (d : ℕ × ℕ × List ℕ × List ℕ × List (Cell × ℕ) × List (Cell × ℕ)) : Bool :=
  let R := d.1; let C := d.2.1; let A := d.2.2.1; let B := d.2.2.2.1
  let f := ofList d.2.2.2.2.1
  let U := skewCells (compl R C A) []
  let F := holesWith tilde R C A B f
  decide (InBox R C (skewCells A B).toFinset) &&
  decide (ColInj (skewCells A B).toFinset f) &&
  decide (StdFrom ∅ (U ++ F)) &&
  decide (SlidesOK revRule (R + C) U.toFinset uVal F) &&
  decide (SlidesOK fwdRule (F.length + 1) F.toFinset
    (rr R C (((skewCells A B).map f).foldr max 0) f) U.reverse)

/-- **Control.** The hypotheses of the assembled theorem hold on all 310 `phiData33` instances
(`|A/B| ≤ 5`, `R, C ≤ 3`), so it is not vacuous there. An evaluation, not a theorem. -/
theorem g4_hyps_control_33 : PhiOutData.phiData33.all hypCheckG4 = true := by decide +kernel

end TworowD4Kernel.G4Assembly
