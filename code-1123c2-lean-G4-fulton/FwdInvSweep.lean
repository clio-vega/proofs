/- SCRATCH (LEAN 2026-11-23 c2): #eval FwdInv at every step of every forward slide of the hokF
run, on the 310 phiData33 instances; plus the tie-break lesion (tie → east). -/
import TworowD4Kernel.SkewColInj
open TworowD4Kernel.SlideSwitch TworowD4Kernel.StdSlide TworowD4Kernel.PhiOut TworowD4Kernel.PhiGlue
open TworowD4Kernel.SkewColInj

/-- Bool form of `HoleOK` (row-weak and column-strict ACROSS the hole). -/
def holeOKb (S : Finset Cell) (X : Cell → ℕ) (h : Cell) : Bool :=
  decide ((∀ w ∈ S, (w.1, w.2 + 1) = h → (h.1, h.2 + 1) ∈ S → X w ≤ X (h.1, h.2 + 1)) ∧
    (∀ w ∈ S, (w.1 + 1, w.2) = h → (h.1 + 1, h.2) ∈ S → X w < X (h.1 + 1, h.2)))

/-- cells of `S` inside a 7 × 7 window (all phiData33 states live in 3 × 3 boxes). -/
def cellsOf (S : Finset Cell) : List Cell :=
  ((List.range 7).flatMap fun r => (List.range 7).map fun c => (r, c)).filter (· ∈ S)

def skewb (S : Finset Cell) : Bool :=
  (cellsOf S).all fun a => (cellsOf S).all fun b =>
    (List.range (b.1 + 1)).all fun r => (List.range (b.2 + 1)).all fun c =>
      !(decide (a.1 ≤ r) && decide (a.2 ≤ c)) || decide ((r, c) ∈ S)

def invb (S : Finset Cell) (X : Cell → ℕ) (h : Cell) : Bool :=
  decide (h ∉ S) && skewb (insert h S) && decide (SSOn S X) && holeOKb S X h

/-- Check the invariant at every step of one slide, with choice function `ch`. -/
def traceb (ch : Finset Cell → (Cell → ℕ) → Cell → Option Cell) :
    ℕ → Finset Cell → (Cell → ℕ) → Cell → Bool × ℕ
  | 0, S, X, h => (invb S X h, 1)
  | k + 1, S, X, h =>
    match ch S X h with
    | none => (invb S X h, 1)
    | some y => let r := traceb ch k (moveS S h y) (moveX X h y) y
                (invb S X h && r.1, r.2 + 1)

/-- One slide (same as `slide` but with a raw choice function). -/
def slideb (ch : Finset Cell → (Cell → ℕ) → Cell → Option Cell) :
    ℕ → Finset Cell → (Cell → ℕ) → Cell → Finset Cell × (Cell → ℕ)
  | 0, S, X, _ => (S, X)
  | k + 1, S, X, h =>
    match ch S X h with
    | none => (S, X)
    | some y => slideb ch k (moveS S h y) (moveX X h y) y

def runb (ch : Finset Cell → (Cell → ℕ) → Cell → Option Cell) (k : ℕ) :
    Finset Cell → (Cell → ℕ) → List Cell → Bool × ℕ
  | _, _, [] => (true, 0)
  | S, X, h :: hs =>
    let t := traceb ch k S X h
    let s := slideb ch k S X h
    let r := runb ch k s.1 s.2 hs
    (t.1 && r.1, t.2 + r.2)

/-- Lesion: the tie goes EAST. -/
def fwdChoiceE (S : Finset Cell) (X : Cell → ℕ) (h : Cell) : Option Cell :=
  if nN h ∈ S then
    if nE h ∈ S then (if X (nN h) < X (nE h) then some (nN h) else some (nE h))
    else some (nN h)
  else if nE h ∈ S then some (nE h) else none

def rowb (ch : Finset Cell → (Cell → ℕ) → Cell → Option Cell)
    (d : ℕ × ℕ × List ℕ × List ℕ × List (Cell × ℕ) × List (Cell × ℕ)) : Bool × ℕ :=
  let R := d.1; let C := d.2.1; let A := d.2.2.1; let B := d.2.2.2.1
  let f := ofList d.2.2.2.2.1
  let cs := skewCells A B
  let U := skewCells (compl R C A) []
  let F := holesWith tilde R C A B f
  runb ch (F.length + 1) F.toFinset (rr R C ((cs.map f).foldr max 0) f) U.reverse

#eval let rs := TworowD4Kernel.PhiOutData.phiData33.map (rowb fwdChoice)
  (rs.length, (rs.filter (·.1)).length, (rs.map (·.2)).sum)
#eval let rs := TworowD4Kernel.PhiOutData.phiData33.map (rowb fwdChoiceE)
  (rs.length, (rs.filter (·.1)).length, (rs.map (·.2)).sum)
-- ablation: drop each HoleOK clause; does the invariant still hold everywhere? (informative only)

/-! ## G5 description (SlidesOK-free): run invariant at the start of slide i
`T_i := (U.reverse.drop i)` = remaining holes. Check: S_i ∩ T_i = ∅; S_i ∪ T_i convex;
no cell of S_i lies weakly SW of a cell of T_i (T_i is a lower set of S_i ∪ T_i);
SSOn S_i X_i; the slide stops within fuel (final choice = none); S_0 ∪ T_0 = B^∨. -/
def runInvb (S : Finset Cell) (X : Cell → ℕ) (T : List Cell) : Bool :=
  (T.all fun t => decide (t ∉ S)) && skewb (S ∪ T.toFinset) &&
  (T.all fun t => (cellsOf S).all fun z => !(decide (z.1 ≤ t.1) && decide (z.2 ≤ t.2))) &&
  decide (SSOn S X)

def endsb (k : ℕ) (S : Finset Cell) (X : Cell → ℕ) (h : Cell) : Bool :=
  let r := slide fwdRule k S X h
  (fwdChoice r.1 r.2.1 r.2.2.1).isNone

def g5run (k : ℕ) : Finset Cell → (Cell → ℕ) → List Cell → Bool
  | S, X, [] => runInvb S X []
  | S, X, h :: hs =>
    runInvb S X (h :: hs) && endsb k S X h &&
    g5run k (slide fwdRule k S X h).1 (slide fwdRule k S X h).2.1 hs

def g5row (d : ℕ × ℕ × List ℕ × List ℕ × List (Cell × ℕ) × List (Cell × ℕ)) : Bool × Bool :=
  let R := d.1; let C := d.2.1; let A := d.2.2.1; let B := d.2.2.2.1
  let f := ofList d.2.2.2.2.1
  let cs := skewCells A B
  let U := skewCells (compl R C A) []
  let F := holesWith tilde R C A B f
  (g5run (F.length + 1) F.toFinset (rr R C ((cs.map f).foldr max 0) f) U.reverse,
   decide (F.toFinset ∪ U.toFinset = (skewCells (compl R C B) []).toFinset))

#eval let rs := TworowD4Kernel.PhiOutData.phiData33.map g5row
  (rs.length, (rs.filter (·.1)).length, (rs.filter (·.2)).length)
-- control: with the holes in FORWARD order (U, not U.reverse) the lower-set/maximal-start fails
def g5rowFwd (d : ℕ × ℕ × List ℕ × List ℕ × List (Cell × ℕ) × List (Cell × ℕ)) : Bool :=
  let R := d.1; let C := d.2.1; let A := d.2.2.1; let B := d.2.2.2.1
  let f := ofList d.2.2.2.2.1
  let cs := skewCells A B
  let U := skewCells (compl R C A) []
  let F := holesWith tilde R C A B f
  g5run (F.length + 1) F.toFinset (rr R C ((cs.map f).foldr max 0) f) U
#eval ((TworowD4Kernel.PhiOutData.phiData33.map g5rowFwd).filter id).length
