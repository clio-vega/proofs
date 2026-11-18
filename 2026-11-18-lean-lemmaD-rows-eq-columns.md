# LEAN 2026-11-18: G4 Lemma D, switch form (C3, the GrowthConfluence dictionary), sorry-free

**Target.** `TworowD4Kernel.LemmaD.lemmaD_rows_eq_columns`, in project `lean/tworow_d4_kernel`
(`TworowD4Kernel/LemmaD.lean`), tworow-d4-kernel `9df5be9`.
**Paper.** `proofs/2026-11-13c2-Q585-G4-evacuation.tex`, `lem:D`, l.187–198.

## Statement
For all `m, p` and every chain `g : ℕ → Finset Cell`,
`swSeq (rowsOrder p m) (initWord p m, g) = swSeq (colsOrder m p) (initWord p m, g)`.
The elementary switches applied by rows ((m,1),…,(m,p),(m−1,1),…) and by columns ((m,1),(m−1,1),…,(1,1),(m,2),…)
to `(U_1⋯U_m F_1⋯F_p, η)` give the same word `F^p U^m` and the same chain, so in particular the same U-part.

## Encoding (decided, with the reason)
U ↦ `Dir.D`, F ↦ `Dir.R`, `f t = η_t`. This is forced: `Flip` rewrites `D R ↦ R D`, and both orders move U rightward past F.
So the paper's `[m]×[p]` grid is GrowthConfluence's with the parameters swapped (`U^m F^p = initWord p m`).
`lesion_encoding` (kernel) shows that under U ↦ R every switch of the rows order is illegal.

## What builds (0 sorries)
- `switch_is_flip` / `switch_is_flip'`: one switch on a UF pair = one `Flip switchL` step.
- `carry_right`, `carry_left`: one letter carried past k letters is a flip sequence with the expected word.
- `rows_reach`, `cols_reach`: both products are `ReflTransGen (Flip switchL)` from `(U^m F^p, g)` to `(F^p U^m, ·)`, for all m, p.
- `lemmaD_rows_eq_columns`: from `confluence_final`.
- `lemmaDCheck_true N` for **every** N (corollary); `lemmaDCheck_le5` (231 instances; 120 nontrivial, `#guard nontrivCount 5 = 120`).
- Control `control_exW` (kernel): the rows product really moves the chain on (3,1). Lesion `lesion_order` (kernel): `[1,2,0,0]` does not reach F²U².

## Axioms
```
'TworowD4Kernel.LemmaD.lemmaD_rows_eq_columns' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.LemmaD.lemmaDCheck_le5' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.LemmaD.lesion_order' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.LemmaD.lesion_encoding' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.LemmaD.control_exW' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.LemmaD.switch_is_flip'' depends on axioms: [propext, Classical.choice, Quot.sound]
```
`lake build` green (3497 jobs) and `lake test` green. The test driver failed on its own on a missing `open`, which shows that it fires.

## Findings
1. **A chain-level lesion cannot bite.** `confluence_final` holds for every step-blind rule L. `switchL` is also involutive (`switchL_invol`),
   so by `good_unflip` *any* switch sequence that reaches F^pU^m gives the same chain. The brief's "L ≠ switchL must be refuted" therefore
   has no step-blind instance. Only the word (`lesion_order`) or a step-indexed rule (existing `lesion_step_indexed`) can witness a lesion.
2. **The brief's labels are swapped; the tex is right.** Rows (U_i carried past F_1…F_p) are *forward* slides (`slide_eq_switch_product`),
   i.e. Σ. Columns (F_j carried past U_m…U_1) are the *reverse* slides of the outward procedure (`reverse_slide_eq_switch_product`).

## Not covered (so the parent stays `proved`)
- **Slide-iteration glue**, which is mechanical and needs no new mathematics. (a) A window lemma: switches at positions i…i+p−1 act on the full
  chain as `carry p` on `η_i ⊂ ⋯`. (b) Propagation of `IsDiag`/`StdFrom` from one row or column to the next (from `switchL_saturated`).
  Together these would make "Σ = rows product" and "outward = columns product" Lean statements about slides.
- The parent `H1-G4-rotation-map-is-switching` also rests on lem:A, lem:C and thm:G4, which are unformalised.

## Registry
New child `G4-lemmaD-rows-eq-columns-formal` (lean-verified, `lean` = `TworowD4Kernel.LemmaD.lemmaD_rows_eq_columns`). The parent stays `proved`
and its approach records why. `registry_validate` OK. trustcheck `--files-dir .` last line: `OK: … is valid`.
