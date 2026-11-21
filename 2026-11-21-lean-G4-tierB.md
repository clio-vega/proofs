# LEAN 2026-11-21: G4 Tier A done (sorry-free), Tier B not reached

**Target.** The brief named `thm_G4_of_slideOK` (Tier B). What landed is **Tier A**: gaps 1 and 3 are proved, gap 2 is resolved
as a paper gap, and a capstone lemma `phi_V_colsW` connects `PhiOut.phi`'s own run to the window form.
Project `lean/tworow_d4_kernel`, new module `TworowD4Kernel/PhiGlue.lean` (in the root imports). Commit `22f094b`, pushed.
`origin/main` = `22f094b`, checked with `ls-remote`.
**Paper.** `proofs/2026-11-13c2-Q585-G4-evacuation.tex`, `def:phi` (l.98–106) and the proof of `thm:G4`. Purbhoo arXiv:0705.1184 §5.1 Cor 5.3.

## Gap 2, from reading (done first, as the brief asked)
The tex's `def:phi`(1) says that `f̃` "is a standard tableau on `B^∨/A^∨`" and gives **no proof**. Lemma D needs it, because Lemma D
takes `F` to be standard. `lem:C` gives `f̃ = std(rr f)`. That makes `f̃` standard once two further facts are known, and the tex
states neither: (a) `rr` preserves semistandardness, and (b) `std` of a semistandard tableau is standard.
**Paper gap, recorded and not patched.** In Lean the capstone below needs only `Nodup (U ++ holes)`. That is a named
hypothesis (`hnd`). `StdFrom ∅ (U ++ holes)` is not used in this file. Tier B will need it as a named hypothesis.

## What is proved (0 sorries)
- **Gap 1, `std_uVal_eq_valOf`.** The word is `U = skewCells κ []` itself (rows bottom to top, each west to east), which is exactly
  the `(value, column)` order of `uVal = row + 1`. The result holds for **every** list `κ` (no partition hypothesis). The proof uses
  `idxOf_eq_card_lt`: in a strictly lex-sorted list, a cell's index equals the number of cells below it. It also uses `skewCells_sorted`
  (via `List.pairwise_flatMap`). No `Finset.sort` is used.
- **Gap 3, `slidesOK_fuel_iff`.** `SlidesOK revRule` at fuel `k ≥ |S|` is equivalent to `SlidesOK revRule` at any `k'` covering every hole's `r+c`.
  The building blocks are `slideOK_fuel_succ`, `_add` and `_iff`, which mirror `slide_fuel_succ` (when the fuel is that small, no move is possible).
- **`holes_dsum_le`.** Every hole of Φ is a rotated cell, so `r+c ≤ R+C`. This holds unconditionally, so Φ's fuel covers it.
- **Capstone `phi_V_colsW`.** Let `U = skewCells (A^∨) []` and `F = holesWith tilde R C A B f`. Assume `Nodup (U ++ F)` and the
  named `SlidesOK revRule (R+C) U.toFinset uVal F`. Then the cells vacated by Φ's run equal `take |F| (colsW |U| 0 |F| (U ++ F))`.

## Lesions and controls (kernel `decide +kernel`)
- `capstone_control_33`: on all 310 `phiData33` instances, `Nodup`, `SlidesOK` and the conclusion all hold. So the hypotheses are met and the
  statement is not vacuous. 211/310 instances have `|U| ≥ 2`.
- `lesion_U_reversed`: `U.reverse` (still nodup, so the lesion breaks only `std = valOf`) fails on 98/310.
  The explicit instance is `A=(1,1)` in the 2×2 box, where the true `U` passes.
- **Not done:** the brief's "assembly without SlideOK must fail" lesion. It belongs to Tier B. `SlidesOK` holds on all 310 instances,
  so dropping it there cannot fire. The 11-19 counterexample would be needed.

## Build, axioms
`lake build` is green (3504 jobs) and `lake test` exits 0. `#print axioms` (`proofs/code-1121-lean-tierB/axioms.log`) gives 6/6 lines equal to
`[propext, Classical.choice, Quot.sound]`: `phi_V_colsW`, `std_uVal_eq_valOf`, `slidesOK_fuel_iff`, `holes_dsum_le`,
`capstone_control_33` and `lesion_U_reversed`.

## Tier B: what remains
1. **Gap 4, the forward mirror.** `rowsW` ↔ `slides fwdRule` via `jdt_bridge`, plus a forward fuel lemma (`sigmaRR` fuel `|S|+Σκ+1`). Not attempted.
2. **Gap 5 (new, missed by the 11-20 note): destandardising values.** `lemmaA_slides` gives the final values of the `std`-run
   as `std(final X-values)`. But Φ outputs the final *unstandardised* `uVal` values, and Σ outputs row labels. Both need
   a lemma that recovers `X' = (sorted values of X) ∘ std X'` from the multiset of values being preserved. Alternatively,
   a "same path ⇒ values transport" lemma with `X = g ∘ Y`, `g` monotone. This is bookkeeping, but it is
   real Lean work that the 4-gap list did not count.
3. `StdFrom ∅ (U ++ F)` as a named hypothesis (the gap 2 paper gap), and then the chain `lemmaD_slides` → `jdt_bridge` → `lemmaA_slides fwdRule` → `lemmaC`.

## Registry
New child `G4-tierA-phi-V-colsW` (lean-verified, `lean` = `TworowD4Kernel.PhiGlue.phi_V_colsW`) under
`H1-G4-rotation-map-is-switching`. The parent stays **proved**. `registry_validate`: `OK`. trustcheck (from `memory/`,
`--deployment ../code/clio.json --sources reading/sources.json validate … --files-dir /home/clio/projects`) ends with
`OK: ../proofs/registry/q585-ybe-train-vs-commutor.json is valid (status: in-progress, deployment: clio)`.
