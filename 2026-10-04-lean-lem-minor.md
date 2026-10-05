# LEAN 2026-10-04 c4 — `lem:minor` formalised

**Target:** `lem:minor` from `proofs/2026-10-04-c3-Q-lorentzian.tex`, lines 487–510
(`proofs @ 94e6b0d`).

**Project:** `~/projects/lean/tworow_d4_kernel`, new file
`TworowD4Kernel/MinorBound.lean`, wired into the root aggregator
`TworowD4Kernel.lean` (one added import line).

**Toolchain:** `leanprover/lean4:v4.30.0`, Mathlib `v4.30.0`.
`lake` is NOT on PATH here; it lives at
`/home/clio/.elan/toolchains/leanprover--lean4---v4.30.0/bin/lake`.

## Main declaration — sorry-free

```lean
theorem TworowD4Kernel.minor_nonpos_of_sigPos_le_one
    {n : Type*} [Fintype n] [DecidableEq n]
    (M : Matrix n n ℝ) (hsymm : M.IsSymm)
    (hdiag : ∀ i, 0 ≤ M i i) (hsig : sigPos M.toQuadraticForm' ≤ 1)
    (i j : n) : M i i * M j j - M i j ^ 2 ≤ 0
```

Supporting declarations, all sorry-free:

- `toQuadraticForm'_apply` — `M.toQuadraticForm' x = x ⬝ᵥ (M *ᵥ x)`
- `sum_eq_pair` — a sum over `univ` of a function vanishing off `{i,j}` is `f i + f j`
- `toQuadraticForm'_apply_of_mem_spanSubset` — on the coordinate plane `span {eᵢ, eⱼ}`
  the form is the binary form `M i i · xᵢ² + 2·M i j · xᵢxⱼ + M j j · xⱼ²`

## Verification, read honestly

- `lake build` (whole project, 3191 jobs): **exit 0**, read from the process's own status
  via `PIPESTATUS[0]`, not from `tail`. Zero warnings, zero errors on `MinorBound`.
- **Sorries: 0.** Scan is `os.walk`-based, not `grep` (the `grep` here is a `ugrep`
  wrapper honouring `.gitignore`, and `lean/.gitignore` ignores `tworow_d4_kernel/`
  wholesale). The scan was **validated against a planted canary**: with
  `TworowD4Kernel/Canary.lean` containing `theorem canary : 1 = 1 := by sorry` present the
  scan reported that line; with it removed the scan reports only one hit, which is the word
  `sorry` inside a prose comment in `ParityObstruction.lean:18`. 51 `.lean` files walked.
- `#print axioms`: all four declarations depend on
  **`[propext, Classical.choice, Quot.sound]`** — the standard three, nothing else.

## Which rung of the fallback ladder

**Better than rung 1, and by a route the brief did not anticipate.**

The brief's named risk was: *"max dimension of a subspace on which a symmetric form is
positive definite = number of positive eigenvalues"* — is it in Mathlib? The answer is a
clean split, and it is worth stating precisely because the two halves have different status:

- **The subspace-dimension half IS in Mathlib**, in a file new enough to be dated 2026:
  `Mathlib/LinearAlgebra/QuadraticForm/Signature.lean` (author David Loeffler). It defines
  `sigPos Q` as *the maximal `finrank` of a subspace on which `Q` is positive definite* —
  the Sylvester inertia index — and supplies exactly the lemma the proof needs,
  `le_sigPos_of_posDef : (Q.restrict V).PosDef → finrank R V ≤ sigPos Q`.
  Note `sigPos` lives at the **root** namespace, not in `QuadraticForm`, even though its
  neighbours `QuadraticForm.sigPos_weightedSumSquares` etc. do not. That cost one build.
- **The eigenvalue bridge is NOT in Mathlib.** `sigPos` occurs in **no other Mathlib file** —
  I checked by `command grep` over all 8094 `.lean` files, with the instrument first
  validated against a known-present string. Nothing connects
  `Matrix.IsHermitian.eigenvalues` to `sigPos`, and no file mentions both
  `weightedSumSquares` and `eigenvalue`. There is also **no Cauchy interlacing** in Mathlib
  (`grep -i interlac` over all of Mathlib: zero hits), so the one-line interlacing proof
  (`λ₂(A) ≤ λ₂(M) ≤ 0` and `λ₁(A) ≥ M i i ≥ 0`, hence `det A = λ₁λ₂ ≤ 0`) is not available.

So I did not need the brief's rung-1 hand-rolled pigeonhole: Mathlib's
`sigPos_add_finrank_le_of_nonpos` already *is* that pigeonhole, and
`le_sigPos_of_posDef` packages the direction I wanted. The whole proof is ~45 lines.

## The one thing not formalised, stated exactly

The hypothesis is `sigPos M.toQuadraticForm' ≤ 1`, not a count of positive eigenvalues.
These are equal by the uniqueness half of Sylvester's law of inertia applied to the
spectral decomposition, and Mathlib has both ingredients separately but not the join. The
missing lemma is:

```lean
lemma sigPos_toQuadraticForm'_eq_card_pos_eigenvalues
    (M : Matrix n n ℝ) (hM : M.IsHermitian) :
    sigPos M.toQuadraticForm' = {k | 0 < hM.eigenvalues k}.ncard
```

The proof is available in principle: `Matrix.IsHermitian.spectral_theorem` gives
`M = U * diagonal (eigenvalues) * star U` with `U = hM.eigenvectorUnitary`, so
`x ↦ (star U) *ᵥ x` is a `QuadraticMap.IsometryEquiv` from `M.toQuadraticForm'` to
`weightedSumSquares ℝ hM.eigenvalues`, and then
`QuadraticForm.sigPos_of_equiv_weightedSumSquares` closes it. Constructing that
`IsometryEquiv` in matrix form is the work; I did not start it, because at ~3 minutes per
build cycle in this container it was not finishable inside the hour, and a half-built
bridge would have put a `sorry` in the file to no purpose.

**This is a gap in coverage, not a gap in the mathematics.** `sigPos ≤ 1` is the
signature-theoretic reading of "at most one positive eigenvalue", and — this matters — it
is what the paper's proof *actually uses*. The proof never counts an eigenvalue; it
produces a 2-dimensional positive-definite subspace and contradicts a dimension bound.

## Did the paper proof have to change to survive the type checker?

No gap. Two observations, both small, both honest:

1. **One of the paper's steps is redundant.** The paper writes: *"so `M_ii`, `M_jj` are
   nonzero, hence both strictly positive by the hypothesis `M_ii ≥ 0`; with `det A > 0`
   this makes `A` positive definite."* Only **`M_ii > 0`** is load-bearing. Positive
   definiteness of the binary form follows from `M_ii > 0` and `det A > 0` alone by
   completing the square, `M_ii · q(x) = (M_ii xᵢ + M_ij xⱼ)² + (det A)·xⱼ²`, which is what
   `nlinarith [sq_nonneg (M i i * x i + M i j * x j)]` discharges. The claim about `M_jj` is
   true but unused.
2. **The `i = j` case is an identity, not an inequality.** `M i i * M i i - M i i ^ 2 = 0`,
   closed by `ring` then `linarith`. The brief asked me to do this case honestly rather
   than by hand-waving; it is two lines and uses no hypothesis at all — in particular it
   does not need `hdiag`, which is worth noticing since the informal proof cites
   nonnegativity in the same breath.

## The point of the exercise

The defect that made this lemma worth machine-checking was an eigenvalue **miscount**:
sympy's `Poly.count_roots` counts *distinct* roots, so `4·I₂` scored one positive
eigenvalue, and the 600-matrix validation run saw 0 disagreements because random integer
matrices almost never have a repeated eigenvalue. That instrument produced a **false**
counterexample to `lem:minor`.

The formalisation **never counts an eigenvalue.** Multiplicity cannot be miscounted by a
proof that works with `finrank` of a subspace, because `finrank` of a span *is* the
multiplicity information, counted once and by construction. The Lean statement is
structurally immune to the class of bug that produced the false counterexample — which is
a better outcome than a second script agreeing with the first.

## Registry

**No promotion made, deliberately.** `proofs/registry/cylindric-lorentzian.json` has 111
nodes and **none of them is `lem:minor`**. The string `lem:minor` occurs only inside the
node `Q-stepC-raw-hessian-general-m` (`trust: proved`), where it names an internal step.

Per the brief's standing rule — *a Lean child does not promote its parent* —
`Q-stepC-raw-hessian-general-m` and the (Q) root **stay where they are**. Promoting a
`proved` parent because one of its internal steps is now `lean-verified` is exactly the
wrong inference.

The right follow-up is to **add** a dedicated child node `lem-minor-nonpos` with
`trust: lean-verified` and
`lean: TworowD4Kernel.minor_nonpos_of_sigPos_le_one`, carrying the caveat above that the
Lean hypothesis is the signature form. I did not make that edit at the session buzzer
without time to run `python3 code/registry_validate.py` and read its output properly —
an unvalidated registry edit is worse than none.

## Predicates, named

Using the three words the brief insists I separate:

- `lem:minor` is **proved** (paper, `2026-10-04-c3-Q-lorentzian.tex` §lines 487–510).
- `lem:minor` is now **formalised** in Lean, in the signature form of its hypothesis.
- The formalisation is **checked**: `lake build` exit 0, 0 sorries by a validated scan,
  standard three axioms.
- The eigenvalue-count form of the hypothesis is **not formalised**. The exact missing
  lemma is displayed above.

## Non-vacuity — the check I keep forgetting to run

A sorry-free theorem whose hypotheses nothing satisfies is worth nothing, and I have filed
that mistake before. So, explicitly:

- **Hypotheses satisfiable, for every `n`:** the all-ones matrix `J_n` is symmetric, has
  `diag = 1 ≥ 0`, and has eigenvalues `(n, 0, …, 0)` — **exactly one** positive. Verified
  numerically for `n = 2,3,4,5`.
- **Conclusion has content, i.e. is not just the `i = j` identity:** for
  `M = [[1,2],[2,1]]` (eigenvalues `−1, 3`, one positive, `diag ≥ 0`) the minor is
  `1·1 − 2² = −3`, strictly negative. So the theorem asserts a real inequality, not `0 ≤ 0`.

A Lean-level non-vacuity control — e.g. `sigPos (J 2).toQuadraticForm' = 1` by exhibiting
`span {(1,−1)}` as a 1-dimensional nonpositive subspace and applying
`sigPos_add_finrank_le_of_nonpos` — is the natural next addition to this file and is
*not* yet present. The numbers above are hand/numpy-verified, not machine-checked.
