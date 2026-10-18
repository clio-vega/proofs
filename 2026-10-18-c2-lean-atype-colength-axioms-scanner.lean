import TworowD4Kernel.AtypeColength

open Lean Elab Command

run_cmd do
  let env ← Lean.getEnv
  let ns := (env.constants.toList.map (·.1)).filter fun n => (`AtypeColength).isPrefixOf n
  let sorted := ns.toArray.qsort (fun a b => a.toString < b.toString)
  let mut lines : Array String := #[]
  let mut nsorry := 0
  let mut nonstd : Array Name := #[]
  for n in sorted do
    let ax ← Lean.collectAxioms n
    let axs := ax.toList
    if axs.contains ``sorryAx then nsorry := nsorry + 1
    for a in axs do
      unless a == ``propext || a == ``Classical.choice || a == ``Quot.sound do
        nonstd := nonstd.push a
    lines := lines.push s!"{n} :: {axs}"
  logInfo (String.intercalate "\n" lines.toList)
  logInfo s!"RECORDS = {sorted.size}   sorryAx-carrying = {nsorry}   nonstandard = {nonstd.toList}"
