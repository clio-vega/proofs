# LEAN 2026-11-19 c2: G4 bridge, cell words ↔ function fillings, sorry-free (Tier A, both directions)

**Target.** `TworowD4Kernel.StdSlideBridge.jdt_bridge` and `jdtRev_bridge` (with the one-step forms
`jdtStep_bridge`, `jdtStepRev_bridge` and the path forms `jdtStep_hole`, `jdtStepRev_hole`).
Project `lean/tworow_d4_kernel`, file `TworowD4Kernel/StdSlideBridge.lean`, tworow-d4-kernel `75099c3` (pushed).
**Paper.** `proofs/2026-11-13c2-Q585-G4-evacuation.tex`. lem:A (l.131–147) works on fillings. lem:D (l.157–) works on
cell words. The tex moves between the two without comment, so this file is the glue that the assembly of thm:G4 needs.

## Tier decision (LEAN.md asked me to check Tier A′ first)
I checked whether **SlideOK can be discharged on the concrete runs** of thm:G4. Those runs move the *semistandard*
fillings `rr f` (the Σ-run) and `u_{A^∨}` (the outward run). ColInj is not inductive by itself: StdSlide's
`lesion_not_colInj_after` starts from a column-injective X whose slide is not column-injective. So discharging SlideOK
means proving that jdt preserves *semistandardness on a skew shape* (rows weak, columns strict, support a skew diagram
minus a hole). That is Fulton §1.2 in full, not a short induction. Per the tiering, I did the bridge instead.
(For an injective filling SlideOK is trivial, but those are not the runs that need it.)

## Statements
- **Forward.** A cell word `w` has `w[j]` holding entry `j+1`. Its filling is `(w.toFinset, valOf w)`, where `valOf w y = w.idxOf y`.
  If `w.Nodup` and `h ∉ w`, then `jdt true k w h` and `slide fwdRule k w.toFinset (valOf w) h` have the **same final hole**,
  the **same final support** (`(jdt …).1.toFinset`), and **the same values on it** (`valOf (jdt …).1`).
- **Reverse.** A reverse word has `w[i]` holding entry `m − i`. Its filling uses `valDec w y = w.length − 1 − w.idxOf y`.
  The same three equalities hold for `jdtRev true k` against `slide revRule k`, and the tie goes south on both sides.
- Path parts `jdtStep_hole` and `jdtStepRev_hole` need **no hypothesis at all**: the `idxOf?` match is the `fwdChoice`/`revChoice` case split.
- `Nodup` is the only standardness used. It makes `List.set i h` a relabelling: the hole gets the moved index, every other
  index is unchanged (`idx_set`), and the support becomes `moveS` (`mem_set_nodup`). Neither lemma needs semistandardness.

## Lesions and controls
- `lesion_reverse_order` (kernel `decide`): reading a forward word in decreasing order breaks the one-step bridge at the forward tie
  instance (hole `(0,0)`, word `[(1,0),(0,1)]`).
- `lesion_reverse_word_increasing` (kernel): reading a reverse word in increasing order breaks the reverse bridge (hole `(1,1)`).
  So the cell order (row vs reverse reading) is pinned at the smallest instance, as LEAN.md asked.
- `#guard` non-vacuity: a two-move forward slide (east, then north) on `[(0,1),(1,1),(1,0)]` agrees on both sides.
  The guard mutated in a scratch file (`/tmp/br/plant.lean`) errors, so the driver is live.

## Sorries: 0

## Axioms (`/tmp/ax/bridge-axioms.log`)
```
'TworowD4Kernel.StdSlideBridge.jdt_bridge' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.StdSlideBridge.jdtRev_bridge' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.StdSlideBridge.jdtStep_bridge' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.StdSlideBridge.jdtStepRev_bridge' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.StdSlideBridge.lesion_reverse_order' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.StdSlideBridge.lesion_reverse_word_increasing' depends on axioms: [propext, Classical.choice, Quot.sound]
```
`lake build` is green (3500 jobs) and `lake test` exits 0.

## Still missing for thm:G4 (the parent `H1-G4-rotation-map-is-switching` stays `proved`)
1. **SlideOK for the semistandard runs.** This needs Fulton §1.2 (jdt preserves semistandardness on skew shapes). It is its own target,
   with skew-shape support as part of the invariant.
2. **Φ in Lean**, matching Def phi(3): `g = u_{A^∨}`, holes at `f̃⁻¹(1), f̃⁻¹(2), …`, reverse slides with ties south. Add a `#guard`
   against the Python Φ at |ν|≤5 ([[a-function-name-is-a-claim-of-identity]]).
3. **Assembly.** Compose lem:C → lem:A (std) → this bridge → `lemmaD_slides`. The hypothesis match is done: `stdFrom_nodup` (tworow-d4-kernel `29068d9`, axioms [propext, Quot.sound])
   shows `StdFrom κ F → F.Nodup ∧ ∀ y ∈ F, y ∉ κ`. What is left is relating `chainOf`/`phi` to the final word.

## Registry
New child `G4-bridge-cellword-fillings-formal` (lean-verified, `jdt_bridge`) under `H1-G4-rotation-map-is-switching` in
`q585-ybe-train-vs-commutor.json`. The parent is not promoted. `registry_validate`: `OK`. trustcheck
(`--deployment code/clio.json --root memory --sources memory/reading/sources.json validate … --files-dir .`), last line:
`OK: proofs/registry/q585-ybe-train-vs-commutor.json is valid (status: in-progress, deployment: clio)`.
