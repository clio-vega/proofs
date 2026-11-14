# LEAN 2026-11-14: Q559 `thm:main` (and `cor:criterion`) is sorry-free

**Target.** `thm:main` of `proofs/2026-10-23-Q559-vanishing-converse.tex` L360–L393: for |M| ≥ r and distinct ρ-shifts
m_j = θ_j + j − 1, Ψ^C_{(1^r),θ}(M) ≠ 0. This includes Case B (the reflection) and the induction on r. These were the two items on
the 2026-11-13 c2 "still unformalised" list.
**Project.** `lean/tworow_d4_kernel` (github.com/clio-vega/tworow-d4-kernel), new module `TworowD4Kernel/ThmMain.lean`
(337 lines). Commit **`d851f73`** is pushed, and the upstream SHA was confirmed by `gh api`, with an HTTP 422 negative control on a
fake SHA. Full `lake build` rc=0 (read from `lake` itself, not from a pipe), 3488 jobs. ThmMain has **0 warnings** by `awk`;
the same `awk` caught the one ThmMain linter warning in an earlier build, so that zero is a real measurement.

## Statement (encoding unchanged from LemOrder)
```
theorem thm_main [Infinite K] (M : Finset ι) (r : ℕ) (θ : Fin r → ℤ) (hcard : r ≤ M.card)
    (hinj : Function.Injective (rhoShiftPaper θ)) :
    ∃ w : ι → K, Admissible w M ∧ psiHead w M r θ ≠ 0
```
"≠ 0 in ℚ(w)" is encoded as: some admissible point (nonzero and pairwise distinct on M) has Ψ ≠ 0. This is the **global** form that
the 11-13 notes said Case B needs. `cor_criterion` (`.tex` L395) is an **iff**:
(∀ admissible w, Ψ(w) = 0) ↔ |M| < r ∨ ∃ j<k, θ_j − θ_k = k − j.

## What builds sorry-free (sorries: **0**)
The estimate was ~150 lines and the actual count was ~170 lines of proof. No new mathematics was needed.
- `psiHead_rank_zero` gives Ψ(M,0,()) = 1.
- `exists_fresh`, `admissible_update_insert` and `exists_admissible` handle admissible points over an infinite K (Finset induction).
- `psiHead_congr`: Ψ(M) sees only w on M.
- `glueOrd`, `glueOrd_ordOK` and `glueOrd_compat` define the ordering: `snocOrd` on T ∋ n, and any `OrdOK` ordering on T ∌ n.
- `thetaHat_injective`: θ̂'s shifts are m ∘ succAbove.
- `caseA_global` (IH at rank k ⇒ some m_j ≥ 0 gives a witness) takes the IH at (M∖{n}, θ̂), extends the witness by a fresh w_n, and
  calls `lem_order`.
- `caseB_global` (all m_j < 0) runs Case A at θ∨ (m∨ = −m) and then applies `lem_reflect` at the witness w. The point **w⁻¹** has
  Ψ_θ(w⁻¹) = (−1)^r e(w)^{2r} Ψ_{θ∨}(w) ≠ 0. The paper's step "ι is injective" becomes "w⁻¹ is admissible" (`admissible_inv`).
- `mainAt_all` is the induction, and `thm_main` and `cor_criterion` follow from it.

**Deviation recorded:** n is *any* element of M. The paper's "largest index" enters only through `SnocCompat`, which `glueOrd`
supplies for every n. This makes the Lean statement more general, not weaker.

## Paper conventions (the 0-indexed ρ-shift trap)
`paper_m_00`: rhoShiftPaper (0,0) = (0,1), which is distinct. `paper_m_10`: rhoShiftPaper (1,0) = (1,1), a tie on the paper's
hyperplane θ_1 − θ_2 = 2 − 1. `inst_thm_main` instantiates `thm_main` at M = {0,1}, θ = (0,0) over ℚ. `inst_value` is the direct
value **Ψ_{(1,1),(0,0)}({2,3}) = 5/9** = (3/4)(8/9)(5/6)·det[[1,2],[1,3]]/(3−2), which is nonzero and involves no existential.

## Lesions on the final theorem (Lean derives False). Every non-lesioned hypothesis is discharged
- `thmMainWith_control`: `ThmMainWith False False` **holds**. It is `thm_main` itself.
- **L-card** (drop |M| ≥ r): **refuted** at M = ∅, r = 1, θ = (0). There is no head injection, so Ψ ≡ 0.
- **L-dist** (drop distinctness): **refuted** at M = {0,1}, r = 2, θ = (1,0) (m = (1,1), |M| = r holds). The forward law kills it.
- **Pair**: **refuted** a fortiori. Each single is refuted with the other hypothesis discharged, so joint fatality is not a
  separate question here.
- **L-infinite** (drop `Infinite K`): **refuted** at K = ZMod 2, M = {0,1}. No admissible point exists there at all. So the
  encoding needs K to be large, as ℚ(w) is.

## `#print axioms` (same run as a `sorryAx` control)
```
thm_main, mainAt_all, caseA_global, caseB_global, cor_criterion, inst_thm_main, inst_value,
thmMainWith_control, lesion_{card,dist',pair',infinite}_refutable, psiHead_congr, exists_admissible
                                   : [propext, Classical.choice, Quot.sound]   (14/14)
sorryControl1114                   : [sorryAx]
```

## Registry (`proofs/registry/rick-q488-converse.json`, backup `.bak-1114-lean`)
- **`clio-Q559-two-block-leading-order` is PROMOTED to `lean-verified`.** Every step it cites (block B, block A, order
  comparison) is now a lean-verified child. The 11-13 reason for not promoting it (the case split and the induction) belonged to
  thm:main, not to this node.
- **New node `clio-Q559-thm-main-formal`** (`lean-verified`) under the top node.
- **Top node `clio-Q553-vanishing-converse-Q559` NOT promoted.** Its child `clio-Q559-branching-rule` (cor:branch) is unformalised,
  and the bridge between "pointwise at admissible points" and "as an element of ℚ(w)" is meta-level. That bridge is the
  evaluation homomorphism in one direction and Zariski density of admissible ℚ-points in the other, and neither is in Lean.
- `registry_validate --proofs-dir .` reports **OK**. Planted control (an invalid trust value, a missing file, and a block child
  downgraded to `computed`): **4 problems**. All 3 plants fired, and the 4th line is the boundary rule cascading up from the bad
  trust value.

## Still unformalised
1. `cor:branch` (the w_n → ∞ limit statement, node `clio-Q559-branching-rule`). The pointwise ingredients are already there
   (`psiHead_order_identity`, `regAt0_orderF`, `hBlock_zero`, `orderF_zero`). What is missing is a *limit* or *valuation* statement.
2. The ℚ(w) ↔ pointwise bridge (above).
