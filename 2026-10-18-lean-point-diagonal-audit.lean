import TworowD4Kernel.PointDiagonalNoGo
import Lean
open Lean Elab Command

/-- The audit, reused verbatim from `2026-10-16-lean-Q473-rearrangement-audit.lean`.  ONE
copy of the code over two domains: the real module (imported) and a plant (declared
locally).  Records are counted from the ENVIRONMENT, not by grepping lines, so `@[simp]`,
wrapped signatures and multi-line proofs cannot hide a declaration -- and DECL_COUNT is
asserted against RECORD_COUNT so "no records" cannot read as "no contamination". -/
def audit (names : Array Name) (label : String) : CommandElabM Unit := do
  let env ← getEnv
  let keep := names.filter fun n => !n.isInternalDetail
  let mut records := 0
  let mut contaminated := 0
  let mut nonstandard := 0
  for n in keep do
    let axs ← collectAxioms n
    let kind := match env.find? n with
      | some (.thmInfo _)   => "thm"
      | some (.defnInfo _)  => "def"
      | some (.axiomInfo _) => "AXIOM"
      | some (.inductInfo _) => "induct"
      | some (.ctorInfo _)  => "ctor"
      | some _              => "other"
      | none                => "MISSING"
    let axStr := if axs.isEmpty then "NONE" else String.intercalate "," (axs.toList.map toString)
    let bad := axs.any fun a => a == ``sorryAx
    let nonstd := axs.any fun a =>
      a != ``propext && a != ``Classical.choice && a != ``Quot.sound
    if bad then contaminated := contaminated + 1
    if nonstd then nonstandard := nonstandard + 1
    records := records + 1
    let flag := if bad then "  <<<SORRYAX" else ""
    logInfo s!"{label}-RECORD {records} :: {kind} :: {n} :: {axStr}{flag}"
  logInfo s!"{label}-DECL_COUNT {keep.size}"
  logInfo s!"{label}-RECORD_COUNT {records}"
  logInfo s!"{label}-SORRYAX_RECORDS {contaminated}"
  logInfo s!"{label}-NONSTANDARD_AXIOM_RECORDS {nonstandard}"
  if keep.size != records then
    logError s!"{label}: DECL_COUNT != RECORD_COUNT"

run_cmd do
  let env ← getEnv
  let some idx := env.getModuleIdx? `TworowD4Kernel.PointDiagonalNoGo
    | throwError "module TworowD4Kernel.PointDiagonalNoGo not in environment"
  audit env.header.moduleData[idx]!.constNames "REAL"

/-! ## Plant control.  Same `audit` function, second domain.

Without this the audit's "SORRYAX_RECORDS 0" is a claim about an instrument that has never
been seen to fire. -/

namespace PointDiagonalPlant
open TworowD4Kernel.PointDiagonal

variable {ι K : Type*}

/-- PLANT 1: the no-go with its proof removed.  It is the kernel of every no-go theorem in
the module, so it is the right place to plant. -/
theorem no_single_cell_move_plant [CommMonoid K] {χ : ι → K} {S : Finset ι}
    (_h : ∀ m ∈ S, Psi χ {m} = 1) : Psi χ S = 1 := sorry

/-- PLANT 2: a theorem that *uses* plant 1, with no `sorry` of its own.  Tests propagation
through the proof term, not just direct detection. -/
theorem no_single_cell_move_pair_plant [DecidableEq ι] [CommMonoid K] (χ : ι → K) (a b : ι) :
    ¬(Psi χ {a, b} ≠ 1 ∧ Psi χ {a} = 1 ∧ Psi χ {b} = 1) := by
  rintro ⟨hne, ha, hb⟩
  refine hne (no_single_cell_move_plant ?_)
  intro m hm
  rcases Finset.mem_insert.mp hm with rfl | hm
  · exact ha
  · rw [Finset.mem_singleton] at hm
    subst hm
    exact hb

/-- CLEAN CONTROL inside the plant domain: if the audit flags this one too it is reporting
the domain rather than the contamination. -/
theorem clean_control [CommMonoid K] (χ : ι → K) (m : ι) : Psi χ {m} = χ m :=
  psi_singleton χ m

end PointDiagonalPlant

run_cmd do
  let env ← getEnv
  let locals := env.constants.map₂.toList.toArray.map (·.1)
  audit (locals.filter fun n => (`PointDiagonalPlant).isPrefixOf n) "PLANT"
