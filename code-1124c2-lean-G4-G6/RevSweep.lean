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


/-! ## G6 (2026-11-24 c2): the hokR run, revRule. -/

def revChoiceW' (S : Finset Cell) (X : Cell → ℕ) (h : Cell) : Option Cell :=
  if nS h ∈ S then
    if nW h ∈ S then (if X (nW h) < X (nS h) then some (nS h) else some (nW h))
    else some (nS h)
  else if nW h ∈ S then some (nW h) else none

def runInvRb (S : Finset Cell) (X : Cell → ℕ) (T : List Cell) : Bool :=
  (T.all fun t => decide (t ∉ S)) && skewb (S ∪ T.toFinset) &&
  (T.all fun t => (cellsOf S).all fun z => !(decide (t.1 ≤ z.1) && decide (t.2 ≤ z.2))) &&
  decide (SSOn S X)

def headsMinb : List Cell → Bool
  | [] => true
  | h :: hs => (hs.all fun b => !(decide (b.1 ≤ h.1) && decide (b.2 ≤ h.2))) && headsMinb hs

/-- run invariant at every slide start -/
def runRb (ch : Finset Cell → (Cell → ℕ) → Cell → Option Cell) (k : ℕ) :
    Finset Cell → (Cell → ℕ) → List Cell → Bool
  | S, X, [] => runInvRb S X []
  | S, X, h :: hs =>
    let s := slideb ch k S X h
    runInvRb S X (h :: hs) && runRb ch k s.1 s.2 hs

def rowR (ch : Finset Cell → (Cell → ℕ) → Cell → Option Cell)
    (d : ℕ × ℕ × List ℕ × List ℕ × List (Cell × ℕ) × List (Cell × ℕ)) : Bool × ℕ × Bool × Bool :=
  let R := d.1; let C := d.2.1; let A := d.2.2.1; let B := d.2.2.2.1
  let f := ofList d.2.2.2.2.1
  let U := skewCells (compl R C A) []
  let F := holesWith tilde R C A B f
  let t := runb ch (R + C) U.toFinset uVal F
  (t.1, t.2, runRb ch (R + C) U.toFinset uVal F, headsMinb F)

-- (#runs, #runs with FwdInv at every step, #steps, #RunInvR at every start, #heads minimal)
#eval let rs := TworowD4Kernel.PhiOutData.phiData33.map (rowR revChoice)
  (rs.length, (rs.filter (·.1)).length, (rs.map (·.2.1)).sum,
   (rs.filter (·.2.2.1)).length, (rs.filter (·.2.2.2)).length)
-- lesion: tie → west
#eval let rs := TworowD4Kernel.PhiOutData.phiData33.map (rowR revChoiceW')
  (rs.length, (rs.filter (·.1)).length, (rs.map (·.2.1)).sum,
   (rs.filter (·.2.2.1)).length, (rs.filter (·.2.2.2)).length)

/-! Symmetry check: at every step state of the hokR run, revChoice S X h = rot (fwdChoice
(rot S) (K - X ∘ rot) (rot h)), K = 100. -/
def rotS (R C : ℕ) (S : Finset Cell) : Finset Cell := S.image (rot R C)
def symb (R C : ℕ) (S : Finset Cell) (X : Cell → ℕ) (h : Cell) : Bool :=
  revChoice S X h == (fwdChoice (rotS R C S) (fun y => 100 - X (rot R C y)) (rot R C h)).map (rot R C)
def symTrace (R C : ℕ) : ℕ → Finset Cell → (Cell → ℕ) → Cell → Bool × ℕ
  | 0, _, _, _ => (true, 0)
  | k + 1, S, X, h =>
    match revChoice S X h with
    | none => (symb R C S X h, 1)
    | some y => let r := symTrace R C k (moveS S h y) (moveX X h y) y
                (symb R C S X h && r.1, r.2 + 1)
def symRun (R C k : ℕ) : Finset Cell → (Cell → ℕ) → List Cell → Bool × ℕ
  | _, _, [] => (true, 0)
  | S, X, h :: hs =>
    let t := symTrace R C k S X h
    let s := slideb revChoice k S X h
    let r := symRun R C k s.1 s.2 hs
    (t.1 && r.1, t.2 + r.2)
#eval let rs := TworowD4Kernel.PhiOutData.phiData33.map fun d =>
    let R := d.1; let C := d.2.1; let A := d.2.2.1; let B := d.2.2.2.1
    let f := ofList d.2.2.2.2.1
    symRun R C (R + C) (skewCells (compl R C A) []).toFinset uVal (holesWith tilde R C A B f)
  (rs.length, (rs.filter (·.1)).length, (rs.map (·.2)).sum)
