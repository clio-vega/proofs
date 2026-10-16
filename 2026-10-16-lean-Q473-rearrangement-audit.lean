import TworowD4Kernel.Rearrangement
import Lean
open Lean Elab Command

/-- The audit, reused verbatim from `2026-10-15-lean-Q473-audit.lean`.  ONE copy of the code
over two domains: the real module (imported) and a plant (declared locally).  Records are
counted from the ENVIRONMENT, not by grepping lines, so `@[simp]`, wrapped signatures and
multi-line proofs cannot hide a declaration — and DECL_COUNT is asserted against
RECORD_COUNT so "no records" cannot read as "no contamination". -/
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
  let some idx := env.getModuleIdx? `TworowD4Kernel.Rearrangement
    | throwError "module TworowD4Kernel.Rearrangement not in environment"
  audit env.header.moduleData[idx]!.constNames "REAL"

/-! ## Plant control.  Same `audit` function, second domain.

Without this the audit's "SORRYAX_RECORDS 0" is a claim about an instrument that has never
been seen to fire. -/

namespace Q473Plant
open Q473

/-- PLANT 1: the level-set characterisation with its proof removed.  It is the kernel of
every theorem in the module, so it is the right place to plant. -/
theorem antitone_le_iff_countAbove_le_plant {N : ℕ} {a : Fin N → ℕ} (_ha : Antitone a)
    (j : Fin N) (t : ℕ) : a j ≤ t ↔ countAbove a t ≤ (j : ℕ) := sorry

/-- PLANT 2: a theorem that *uses* plant 1, with no `sorry` of its own.  Tests propagation
through the proof term, not just direct detection. -/
theorem le_of_countAbove_le_plant {N : ℕ} {a' b' : Fin N → ℕ} (ha' : Antitone a')
    (hb' : Antitone b') (hle : ∀ t, countAbove a' t ≤ countAbove b' t) (j : Fin N) :
    a' j ≤ b' j := by
  rw [antitone_le_iff_countAbove_le_plant ha' j]
  exact le_trans (hle _) ((antitone_le_iff_countAbove_le_plant hb' j _).mp le_rfl)

/-- CLEAN CONTROL inside the plant domain: if the audit flags this one too it is reporting
the domain rather than the contamination. -/
theorem clean_control {N : ℕ} (a : Fin N → ℕ) (t : ℕ) :
    countAbove a t = (Finset.univ.filter (fun i => t < a i)).card := countAbove_filter a t

end Q473Plant

run_cmd do
  let env ← getEnv
  let locals := env.constants.map₂.toList.toArray.map (·.1)
  audit (locals.filter fun n => (`Q473Plant).isPrefixOf n) "PLANT"
