import TworowD4Kernel.G4Assembly
import TworowD4Kernel.PhiOutData
open TworowD4Kernel.SlideSwitch TworowD4Kernel.StdSlide TworowD4Kernel.StdSlideBridge
open TworowD4Kernel.SlideIteration TworowD4Kernel.SlidesColsBridge TworowD4Kernel.PhiOut
open TworowD4Kernel.PhiGlue TworowD4Kernel.Destd TworowD4Kernel.G4Assembly

/-- L2 window check at fuel `k(U,F,R,C)`: support of the revRule run = U'.toFinset and labels agree. -/
def l2At (fuel : ℕ → ℕ → ℕ → ℕ → ℕ)
    (d : ℕ × ℕ × List ℕ × List ℕ × List (Cell × ℕ) × List (Cell × ℕ)) : Bool :=
  let R := d.1; let C := d.2.1; let A := d.2.2.1; let B := d.2.2.2.1
  let f := ofList d.2.2.2.2.1
  let U := skewCells (compl R C A) []
  let F := holesWith tilde R C A B f
  let U' := (colsW U.length 0 F.length (U ++ F)).drop F.length
  let r := slides revRule (fuel U.length F.length R C) U.toFinset (valOf U) F
  decide (r.1 = U'.toFinset) && U'.all (fun a => r.2.1 a == valOf U' a)

-- control: Φ's own fuel R+C; lesions: |U|, |U|-1, 1
#eval (TworowD4Kernel.PhiOutData.phiData33.filter fun d => !l2At (fun _ _ R C => R + C) d).length
#eval (TworowD4Kernel.PhiOutData.phiData33.filter fun d => !l2At (fun u _ _ _ => u + 1) d).length
#eval (TworowD4Kernel.PhiOutData.phiData33.filter fun d => !l2At (fun u _ _ _ => u) d).length
#eval (TworowD4Kernel.PhiOutData.phiData33.filter fun d => !l2At (fun u _ _ _ => u - 1) d).length
#eval (TworowD4Kernel.PhiOutData.phiData33.filter fun d => !l2At (fun _ _ _ _ => 1) d).length
#eval TworowD4Kernel.PhiOutData.phiData33.length

/-- L6 lesion: U not reversed in Σ's partner list. -/
theorem l6_lesion_unreversed : ¬ ∀ (U U' : List Cell), U'.length = U.length → U'.Nodup →
    ∀ a ∈ U', ofList (U'.reverse.zip (U.map fun h => h.1 + 1)) a = gU U (valOf U' a) := by
  intro h
  have := h [(0,0), (1,0)] [(5,5), (6,6)] rfl (by decide) (5,5) (by simp)
  revert this; decide

/-- L6 control on the same instance (the true statement holds). -/
example : ofList ([(5,5), (6,6)].reverse.zip (([(0,0), (1,0)] : List Cell).reverse.map fun h => h.1 + 1)) (5,5)
    = gU [(0,0), (1,0)] (valOf [(5,5), (6,6)] (5,5)) := by decide

#print axioms l6_lesion_unreversed
theorem sorry_control : 1 = 2 := sorry
#print axioms sorry_control
