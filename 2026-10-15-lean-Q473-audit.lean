import TworowD4Kernel.HalfspaceOrder
import Lean
open Lean Elab Command

/-- The audit.  ONE copy of the code, used over two domains: the real module (imported)
and the plant (declared locally below).  Records are counted from the ENVIRONMENT, so
`@[simp]`, wrapped signatures and multi-line proofs cannot hide a declaration. -/
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
  let some idx := env.getModuleIdx? `TworowD4Kernel.HalfspaceOrder
    | throwError "module TworowD4Kernel.HalfspaceOrder not in environment"
  audit env.header.moduleData[idx]!.constNames "REAL"

/-! ## Plant control.  Same `audit` function, second domain. -/

namespace Q473Plant
open Q473 Q473.SP

/-- PLANT 1: a leaf lemma with its proof removed. -/
theorem finite_downset_plant (κ : SP) : {κ' : SP | κ' ≤ κ}.Finite := sorry

/-- PLANT 2: a theorem that *uses* plant 1.  Tests propagation, not just direct detection. -/
theorem wf_plant : WellFounded ((· < ·) : SP → SP → Prop) :=
  Q473.wellFounded_of_finite_lowerSet fun κ =>
    (finite_downset_plant κ).subset fun _ (hx : _ < κ) => hx.le

/-- CLEAN CONTROL inside the plant domain: if the audit flags this one too, it is
reporting the domain rather than the contamination. -/
theorem clean_control (κ : SP) : (⊥ : SP) ≤ κ := bot_le

end Q473Plant

run_cmd do
  let env ← getEnv
  let locals := env.constants.map₂.toList.toArray.map (·.1)
  audit (locals.filter fun n => (`Q473Plant).isPrefixOf n) "PLANT"
