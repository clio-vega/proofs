# LEAN 2026-11-13 — Q559 `lem:blockA` (the block containing n), sorry-free

**Target.** `lem:blockA`, `proofs/2026-10-23-Q559-vanishing-converse.tex` L224–L296.
**Project.** `lean/tworow_d4_kernel` (github.com/clio-vega/tworow-d4-kernel), module
`TworowD4Kernel/BlockA.lean`, commit **`fa6683b`** (pushed; upstream SHA confirmed by `gh api`).
`lake build` rc=0, 3486 jobs, no warnings in BlockA.

## Encoding (shared with blockB, as the 11-12 c2 report asked)
r = k+1. The indices j are 0-indexed, so the paper's sign (−1)^{r+j*} is (−1)^{k+j*} here. m = `rhoShiftPaper θ`;
θ̂_l = m_{j*.succAbove l} − l. At w_n = u⁻¹:
- (I)  `u^(m_{j*} − k) · psiPlus = hBlock u` (`psiPlus_eq`), with
  H(u) = Σ_{T'⊆M',|T'|=k} B_{M'}(T') β(u) L(u) / (V(T') ∏(1 − u x)), and
  L(u) = Σ_j (−1)^{k+j} u^{expo j} D_j, where expo j = (m_{j*} − m_j).toNat ∈ ℕ.
- (II) `hBlock 0 = (−1)^(k+j*) · psiHead M' k θ̂` (`hBlock_zero`), with psiHead given by eq:Psidef.
Regularity of H at 0 is read off its definition: the exponents are natural numbers and the denominators have the form 1 − u·x.
This is the same scope limit as blockB. The ordering of T'∪{n} is the hypothesis `SnocCompat`: T' in order, then n.
`snocOrd_compat` and `snocOrd_ordOK` show that this hypothesis is satisfiable by a genuine ordering.

## The three computations, each a separate lemma
(i) `Bfactor_insert` + `betaRaw_at_inv`. (ii) `vand_succ_last` + `vand_insert_at_inv`.
(iii) `det_laplace_last` (Mathlib `det_succ_row` at `Fin.last`) + `lapReg_zero`.
Recombination is done by `combine_core`, a scalar zpow identity. **Where the hypotheses enter, Lean-checked:** (I) is proved
under `m_j ≤ m_{j*}` only, so distinctness is not used there. Strictness enters only in `lapReg_zero`
(expo j ≠ 0 for j ≠ j*). This matches the .tex's "the one place the distinctness hypothesis enters".

## Pre-measure (`scratch/1113_blockA.py`, exact Fractions, 60 instances, n ≤ 5, r ≤ 3)
(I) 60/60. (II) 60/60 against an independent eq:Psidef evaluator. The w_n = 10⁹ numerics give 60/60.
- **Sign-dropped control:** agrees 49/60, but all 49 agreements are uninformative (sign +1 or target 0). On the 11 informative instances it agrees 0/11.
- **Non-argmax control:** agrees 0/60. Its target is nonzero in 37/60, so this control is not silent.

## Lesions (Lean derives False). Instance: M={0,1}, n=1, k=1, w=(2,3)
- L-dist: drop distinctness and keep the argmax. θ=(1,0), m=(1,1): the tied Laplace terms cancel, giving H(0)=0, while the target is −3/2. **Refuted.**
- L-arg: keep distinctness and use a non-argmax j*. θ=(0,0), m=(0,1), j*=0: H(0)=−3/4 against −3/2. **Refuted.**
- L-sign: keep every hypothesis and drop the sign. θ=(2,0), j*=0, where the sign is −1: −3/2 against +3/2. **Refuted.**
- Pair (L-dist ∧ L-arg): refutable a fortiori.
- `inst_control`: (II) *holds* on the same instance with the sign present, so the lesions are refuting the lesioned statements and not the instance.

## `#print axioms` (same run as a `sorryAx` control)
```
lem_blockA, psiPlus_eq, hBlock_zero, snocOrd_ordOK, lesion_{dist,arg,sign,pair}_refutable, inst_control
  : [propext, Classical.choice, Quot.sound]        (9/9)
sorryControl : [sorryAx]
```
Sorries: **0**.

## What is still informal in thm:main's order comparison
blockA, blockB and reflect are now all Lean-verified. Each, however, is a **pointwise identity in an abstract field at a substituted value u = 1/w_n**, together with a formula for the value at u = 0. Combined with `psiHead_split` they give
u^{m*−k}Ψ = H(u) + (−1)^{k+1}u^{m*+1}G(u) for every admissible u. Nothing formal yet turns that into "Ψ ≠ 0 in ℚ(w_1..w_n)". The missing step is the passage from "this expression equals H(0) ≠ 0 at u = 0" to "it is not identically zero". It needs one of two things:
- an instantiation over K = Frac ℚ[w], with H and G shown to lie in the local ring at u = 0 and the specialisation u ↦ 0 shown to be a ring map that computes the residue (the "regular at 0" fact that both files read off their definitions);
- or a cleared-denominator polynomial-in-u argument plus "a nonzero polynomial has a non-root" over an infinite field, avoiding the finitely many inadmissible u.

Also unformalised: the case split on max(m) ≥ 0 versus ≤ −1, the use of the reflection to cover the second case, and the induction on r. So the comparison itself is **not** formal, and the parent `clio-Q559-two-block-leading-order` stays `peer-reviewed`.
