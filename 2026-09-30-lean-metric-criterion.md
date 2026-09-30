# LEAN 2026-09-30 — `metric-criterion` and `l3-det-reduction`, both sorry-free

**Project:** `lean/tworow_d4_kernel` (`clio-vega/tworow-d4-kernel`).
Session start `git rev-parse HEAD` = `git ls-remote origin HEAD` = **`c100e75`** — printed, not
inherited from the brief. Work pushed as **`b90c3ed`**; local `git rev-parse HEAD` == 
`git ls-remote origin HEAD` == `b90c3ede12a51b42f71e2a3d2dba6e28001d9375`, both printed after the push.

**File:** `TworowD4Kernel/MetricCriterion.lean`, imported from the root `TworowD4Kernel.lean`
(line 28), so it is inside the axiom audit's import closure.

**Paper proof:** `proofs/2026-09-29-c1-cylindric-lorentzian-ell3.tex`, `lem:metric` and
`prop:reduction`. Registry `proofs/registry/cylindric-lorentzian.json`.

## Build evidence

```
Build completed successfully (3175 jobs).
⚠ [3173/3175] Built TworowD4Kernel.MetricCriterion (10s)
✔ [3174/3175] Built TworowD4Kernel (2.0s)
```
Exit 0. The module really was compiled (10s), so this is not a silent no-op build.
`grep -cE "sorry|^axiom|admit"` on the file: **0**. The only warning is a `push_neg`
deprecation notice.

## `#print axioms`

```
'…MetricCriterion.metric_criterion'      depends on axioms: [propext, Classical.choice, Quot.sound]
'…MetricCriterion.l3_det_reduction'      depends on axioms: [propext, Classical.choice, Quot.sound]
'…MetricCriterion.condA_not_sufficient'  depends on axioms: [propext, Classical.choice, Quot.sound]
'…MetricCriterion.W_not_condB'           depends on axioms: [propext, Classical.choice, Quot.sound]
'…MetricCriterion.W_det'                 depends on axioms: [propext, Classical.choice, Quot.sound]
'…MetricCriterion.atMostOnePos_iff_card' depends on axioms: [propext, Classical.choice, Quot.sound]
```
Exactly the standard three everywhere. No fourth axiom.

## The two target declarations

```lean
theorem l3_det_reduction (hM : M.IsHermitian) (hnn : ∀ i j, 0 ≤ M i j) :
    AtMostOnePos hM ↔ (e2 M ≤ 0 ∧ 0 ≤ M.det)

theorem metric_criterion (hM : M.IsHermitian) (hnn : ∀ i j, 0 ≤ M i j)
    (hA : CondA M) (hB : CondB M) :
    e2 M ≤ 0 ∧ 0 ≤ M.det ∧ AtMostOnePos hM
```

Both **re-prove** every step of the paper argument. Nothing is assumed, and in particular
nothing is *defined into* Lean: `M : Matrix (Fin 3) (Fin 3) ℝ`, the eigenvalues are Mathlib's
`Matrix.IsHermitian.eigenvalues` (spectral theorem), and (A), (B) are the paper's quantified
conditions verbatim. This is the contrast with `AffineAdditive`, where the criterion is a
definition and the type checker is green over a statement it cannot falsify.

Two guards against the naming carrying content:
- `AtMostOnePos` is stated pairwise; `atMostOnePos_iff_card` proves it is literally
  `(univ.filter fun i => 0 < eigenvalues i).card ≤ 1`. The name is not load-bearing.
- `trace_e2_det_eq_eigen` **derives** `tr`, `e₂`, `det` as the elementary symmetric functions of
  the eigenvalues, from `Matrix.IsHermitian.charpoly_eq` via `charpoly_eq_cubic`. The bridge from
  entries to spectrum is proved, not posited.

## The brief's three suspect because-clauses — what happened to each

1. **`xyz = P t²`.** True, and a pure `ring` identity in the six entries. Checked numerically
   first (0 failures / 20000 random integer matrices, entries allowed **negative** — so it is an
   identity, not a consequence of nonnegativity), then absorbed into `ring` inside
   `det_nonneg_alg`. No defect.

2. **"(B) multiplied through is exactly `x,y,z ≤ t`", and the `t = 0` branch.**
   **The brief was wrong that the `t = 0` branch is unwritten.** It *is* in the paper, in Step 2:
   "Suppose first `t=0`. Then (B1)–(B3) force `x=y=z=0`, and (A) multiplied out gives
   `t² ≥ P²`, so `P = 0`." I formalised exactly that and it is sound — the three (A) inequalities
   have nonnegative sides, so they multiply (`mul_le_mul`, four nonnegativity side goals), giving
   `(abc)² ≤ (pqr)² = 0`. **No gap.** It is *not* absorbed by `positivity`: the branch is an
   explicit `rcases ht.eq_or_lt`, and `sq_eq_zero_iff` does the final step.
   The multiplier's sign is handled honestly: each (B) instance is multiplied by the
   off-diagonal entry via `mul_le_mul_of_nonneg_left`, with the `hnn` hypothesis supplying the
   sign. Zero diagonal entries need nothing special — the whole proof divides by nothing.

3. **"multilinear hence minimised at a vertex".** This is a real lemma, and **it is not needed.**
   Substituting `P = xyz/t²` and clearing `t²` turns the claim into a polynomial inequality with
   an explicit positivity certificate:

   `xyz + 2t³ − t²(x+y+z) = x(t−y)(t−z) + t(t−x)(t−y) + t(t−x)(t−z)`

   a `ring` identity (verified in sympy before any tactic) whose three summands are products of
   nonnegatives under `0 ≤ x,y,z ≤ t`. This is `cube_certificate`, and it is strictly better
   evidence than the vertex argument: it exhibits the witness instead of quantifying over extreme
   points, and it needs no Mathlib multilinear-minimum lemma. **This is the one place the
   formalisation improves on the paper**, and the paper should be updated to it.

So the paper proof survived formalisation intact. The deviation is a simplification, not a repair.

## Negative controls

**Lean-checked, and they bite:**
- `condA_not_sufficient : ¬ AtMostOnePos W_isHermitian` for `W = !![4,6,4; 6,1,4; 4,4,4]`.
  Note the route: it is derived from the **`⇒` direction of `l3_det_reduction`** plus
  `W_det : W.det = -16`, so no eigenvalue computation is needed — and it is therefore
  simultaneously a live test of that direction of the iff, which nothing else in the file
  exercises. `W_condA` and `W_nonneg` are proved by `norm_num`, i.e. **re-derived, not copied**
  from the brief.
- `W_not_condB : ¬ CondB W`, at `(i,j,k) = (0,1,2)`: `16 < 24`. So (B) is precisely the
  hypothesis whose absence produced the recorded `rlc-implies-l3` counterexample.

**Not formalised — recorded as owed, in the file itself as well as here:**
- *Nonnegativity is load-bearing.* It enters only through `0 ≤ M.trace`, and
  `l3_det_reduction_trace` isolates exactly that, so the dependency is visible in Lean. That the
  weaker hypothesis cannot be dropped is witnessed by `M = -(1)` (`e₂ = 3 > 0`, `det = -1 < 0`,
  no positive eigenvalue) — but the clause "`-1` has no positive eigenvalue" is **not** a Lean
  theorem; it needs the eigenvalues of a diagonal matrix, which this file never computes. I
  attempted it, the `nlinarith` closing step failed, and I removed it rather than leave a sorry.
- *`3 × 3`-ness is load-bearing.* `metric_criterion` is false at `ℓ = 4`:
  `!![3,6,2,6; 6,0,1,4; 2,1,0,4; 6,4,4,1]` satisfies (A) and (B) with nonnegative entries and has
  spectrum `≈ {13.493, 0.081, −4.574, −5}` — **two** positive eigenvalues (re-verified in Python
  this session, not copied). **Python only.** Formalising it needs "a 2-dimensional positive
  subspace forces two positive eigenvalues" (Cauchy interlacing / min–max), and the `det`/`e₂`
  criterion cannot substitute because that criterion is itself `3 × 3`-only. So
  `metric-criterion-general-ell` stays refuted **by computation, not by Lean**, and stays
  `dead-end`.

The honest reading: the brief called the `ℓ=4` control "the sharpest available and the first
thing to check". I checked it numerically and could not formalise it in the time available. What
Lean *does* certify about `3 × 3`-ness is indirect — every step is a `Fin 3` computation
(`det_fin_three`, `coeff_cubic`, six entries) and none of it would typecheck at `Fin 4`.

## Sorry count

**0.** No declaration in this file has a sorry, including in docstrings.
