# LEAN 2026-11-18 c2: G4 Lemma D in slide form, sorry-free (tier C)

**Target.** `TworowD4Kernel.SlideIteration.lemmaD_slides`, together with `sigma_eq_rows_product` and
`outward_eq_columns_product`, in project `lean/tworow_d4_kernel` (`TworowD4Kernel/SlideIteration.lean`),
tworow-d4-kernel `43eaae1`.
**Paper.** `proofs/2026-11-13c2-Q585-G4-evacuation.tex`, `lem:D`, l.172–198. This continues c1
(`2026-11-18-lean-lemmaD-rows-eq-columns.md`) and closes its "Not covered" glue.

## Statements
Chains are recorded by cell words: `lab D₀ w t = D₀ ∪ {w₀,…,w_{t−1}}`. `phi_chainOf` shows this is `phi (chainOf D₀ w)`.

- **(a) Window, forward** `window_forward`: if `|a| = |A|` and `ChainInv (D₀ ∪ A) h F`, then
  `swSeq [|A|, …, |A|+p−1] (a ++ U F^p b, lab D₀ (A ++ h::F ++ B)) = (a ++ F^p U b, lab D₀ (A ++ carryW h F ++ B))`.
  `window_forward_outside`: `η_t` is unchanged for `t ≤ |A|` and for `t ≥ |A|+p+1`.
- **(a) Window, reverse** `window_reverse`: this is the mirror statement. Positions run `|A|+k−1, …, |A|`, the window is `U ++ [h]` and its output is `carryWRev h U`.
- **(b) Propagation**: `stdFrom_carryW` (forward) and `carryRev_chainOf.2` (reverse, already in the library). After a slide the
  whole word is still `StdFrom ∅`, so the next row's or column's slide theorem applies.
- `slideW_eq_carryW` and `revSlideW_eq_carryWRev`: the carries are the jdt slides (`jdt true`, `jdtRev true`).
- **`sigma_eq_rows_product`**: for `StdFrom ∅ w` and `|w| = m+p`, applying `swSeq (rowsOrder p m)` to `(U^m F^p, chain w)` gives `(F^p U^m, chain (rowsW p m w))`.
  `rowsW` means: slide `U_m` forward through the next p cells, then `U_{m−1}`, and so on (this is Σ).
- **`outward_eq_columns_product`**: the same with `colsOrder` and `colsW m 0 p w`. `colsW` means: slide `F_1` backward through the m cells before it, then `F_2`, and so on.
- **`lemmaD_slides`**: `phi (chainOf ∅ (rowsW p m w)) = phi (chainOf ∅ (colsW m 0 p w))`. So Σ and the outward
  procedure give the same chain, hence every letter gets the same cell, and in particular the U-part is the same.

## Sorries: 0

## Axioms
```
'TworowD4Kernel.SlideIteration.lemmaD_slides' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.SlideIteration.sigma_eq_rows_product' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.SlideIteration.outward_eq_columns_product' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.SlideIteration.window_forward' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.SlideIteration.window_reverse' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.SlideIteration.window_forward_outside' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.SlideIteration.control_window' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.SlideIteration.lesion_window_shift' depends on axioms: [propext, Classical.choice, Quot.sound]
```
`lake build` green (3498 jobs, 0 warnings in the new module). `lake test` green. The driver failed on its own on a missing `open`
before I fixed it, so it is live.

## Control, lesion, evaluation
- `control_window` (kernel): on `A=[(0,0)]`, `h=(1,0)`, `F=[(0,1)]` the window `[1]` gives `D R D` and changes `η_2`.
  An `example` instantiates the general `window_forward` at the same point and agrees.
- `lesion_window_shift` (kernel): a window shifted by one position, right (`[2]`) or left (`[0]`), is a no-op there and misses
  the window lemma's output. (As c1 found, the lesion can only show up in the word. Inside the rows context a shifted switch meets `RR` or `DD`.)
- `#guard`s evaluate the definitions directly; none of them uses a theorem. On all 231 instances with `|ν| ≤ 5`, `rowsW p m w = colsW m 0 p w` **as cell words**,
  which is a stronger observable than the chain equality the theorem states. Σ moves the filling on 120 of the 231.

## Not covered
- Cell-word equality `rowsW = colsW` as a theorem. It needs "lab determines the word for fresh words" plus length and
  standardness outputs of `rowsW`/`colsW`. That is mechanical, but I did not do it. The chain equality already fixes every letter's cell.
- The parent `H1-G4-rotation-map-is-switching` stays `proved`, because lem:A, lem:C and thm:G4 are unformalised. Lem:D is now formalised in the tex's own terms.

## Registry
New child `G4-lemmaD-slide-form-formal` (lean-verified, `lean` = `TworowD4Kernel.SlideIteration.lemmaD_slides`) under
`H1-G4-rotation-map-is-switching` in `q585-ybe-train-vs-commutor.json`. The parent is not promoted. `registry_validate` OK. trustcheck
(`--root memory --sources memory/reading/sources.json validate … --files-dir .`) last line: `OK: … is valid`.
