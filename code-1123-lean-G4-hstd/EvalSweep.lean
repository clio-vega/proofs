import TworowD4Kernel.G4AssemblyStd
open TworowD4Kernel.SlideSwitch TworowD4Kernel.StdSlide TworowD4Kernel.PhiOut TworowD4Kernel.PhiGlue
open TworowD4Kernel.G4AssemblyStd

def row (d : ℕ × ℕ × List ℕ × List ℕ × List (Cell × ℕ) × List (Cell × ℕ)) : List Bool :=
  let R := d.1; let C := d.2.1; let A := d.2.2.1; let B := d.2.2.2.1
  let f := ofList d.2.2.2.2.1
  let cs := skewCells A B
  let U := skewCells (compl R C A) []
  let F := holesWith tilde R C A B f
  [ decide (∀ Y ∈ cs, (Y.1, Y.2 + 1) ∈ cs → f Y ≤ f (Y.1, Y.2 + 1)),
    decide (∀ Y ∈ cs, (Y.1 + 1, Y.2) ∈ cs → f Y < f (Y.1 + 1, Y.2)),
    (List.range 5).all fun i => (List.range 5).all fun j => !(i ≤ j) || (A.getD j 0 ≤ A.getD i 0 && B.getD j 0 ≤ B.getD i 0 && B.getD i 0 ≤ A.getD i 0),
    decide (StdFrom ∅ (U ++ F)),
    decide (SlidesOK revRule (R + C) U.toFinset uVal F),
    decide (SlidesOK fwdRule (F.length + 1) F.toFinset (rr R C ((cs.map f).foldr max 0) f) U.reverse) ]

#eval let rs := TworowD4Kernel.PhiOutData.phiData33.map row
  (rs.length, (List.range 6).map fun k => (rs.filter fun r => r.getD k false).length)
