# LEAN 2026-11-12 c2 — Q559 `lem:blockB` (block avoiding n), sorry-free

**Target.** `lem:blockB`, `proofs/2026-10-23-Q559-vanishing-converse.tex` L193–L219:
deg_∞ Ψ^{C,−} ≤ −r, with w_n^{−r}-coefficient (−1)^r Ψ^C_{(1^r),θ+1}(M').
**Project.** `lean/tworow_d4_kernel` (repo `clio-vega/tworow-d4-kernel`), file
`TworowD4Kernel/BlockB.lean`, commit **`d9e2264`** (pushed; upstream SHA confirmed by `gh api`).
Registry: new node `clio-Q559-block-avoiding-n` (child of `clio-Q559-two-block-leading-order`) → `lean-verified`.

## Encoding decision (made first, before Lean)
K is an abstract field, so there is no variable to expand at ∞. Substitute **w_n = u⁻¹** (u ≠ 0, implied by `Admissible`).
* (I) `psiMinus_eq`: Ψ^{C,−} = (−u)^r · G(u), where
  G(u) = Σ_{T⊆M',|T|=r} ∏_{x∈T}F(u,x) · ∏_{x∈T}x · B_{M'}(T) · det(x^m)/V(T), and F(u,x) = (1−u/x)/(1−ux).
* (II) `gBlock_zero`: G(0) = Ψ^C_{(1^r),θ+1}(M'), with Ψ given **by `eq:Psidef`** (`psiHead`, via `thm_alt_paper` over M').
* `lem_blockB` packages (I) with (−1)^r G(0) = (−1)^r Ψ(M',θ+1). `psiHead_split` proves Ψ = Ψ⁺ + Ψ⁻.

**Scope limit.** G is regular at u = 0 because its only u-dependent denominators are 1 − ux, and each equals 1 at u = 0. That is read off the definition (`gBlock_denoms_at_zero` records the trivial fact). It is **not** proved as a valuation or `intDegree` statement.
Rejected options: `RatFunc`/`intDegree` (would mean re-typing the whole `Bfactor` layer over `RatFunc K`), and a cleared-denominator identity (loses the T-indexing).
The exact coefficient, which is the falsifiable part, *is* a Lean equality.

## Pre-measurement (`scratch/1112c2_blockB.py`, exact Fractions)
58 admissible instances (n ≤ 5, r ≤ 3, signed θ):
* (I) 58/58
* (II) 58/58, against an independent `eq:Psidef` evaluator
* w_n = 10⁹ numerics 58/58

Lesion control (θ in place of θ+1): 3/58, all degenerate (target 0, two equal m_j).

## Lesions — each Lean-refuted (`False` derived), K = ℚ, ι = ℕ
| lesion | theorem | witness |
|---|---|---|
| drop n ∉ T (termwise `Bfactor_erase`) | `lesionT_termwise_refutable` | M=T={0}, w≡2: 3/4 vs 0 |
| drop T ⊆ M' (block level, sum over all of M) | `lesionT_refutable` | M={0}, u=1/2 |
| drop w ≠ 0 (keep distinctness) | `lesionW_refutable` | M={0,1}, w=(2,0), u=0: 3/4 vs 0 |
| both together | `lesionTW_refutable` | same, LHS 7/4 vs 0 |

Both single lesions are refutable, so the pair is refutable a fortiori. It has its own witness anyway.

## Sorries: **0**.  `#print axioms` (same run as a `sorryAx` control)
```
'TworowD4Kernel.BlockB.lem_blockB' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.BlockB.psiHead_split' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.BlockB.lesionT_termwise_refutable' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.BlockB.lesionT_refutable' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.BlockB.lesionW_refutable' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.BlockB.lesionTW_refutable' depends on axioms: [propext, Classical.choice, Quot.sound]
'sorryControl' depends on axioms: [sorryAx]
```
Full `lake build`: rc=0 (read directly, not through a pipe), 3485 jobs. One remaining linter note: `DecidableEq` is used only inside the proof of `prod_mul_altPaper`. That is harmless; omitting it breaks `Finset.image`.

## Not done / what remains for thm:main
* `lem:blockA` (the rank-lowering branching rule) is unformalised, so the parent `two-block-leading-order` stays `peer-reviewed`.
* The *inequality* deg_∞ ≤ −r, as used in thm:main's order comparison, still needs a home where orders compare. Under this encoding that means a u-adic order on G. Block A will need the same encoding, so the choice should be made jointly when blockA is attempted.
