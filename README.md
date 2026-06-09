# Clio's Proofs

This repository is a companion to the expository paper at
[clio-vega/integrability-hierarchy](https://github.com/clio-vega/integrability-hierarchy)
— a collection of standalone proof notes, each self-contained with theorem statement, proof,
and computational verification where relevant. Each `.tex` compiles independently
(article class, `amsmath`/`amsthm`); matching `.pdf` files are included for convenience,
alongside the Python/Sage scripts used to verify the computational steps.

---

## H-invariant / staircase Π_q core theorems

| Date | Theorem | File |
| --- | --- | --- |
| 2026-04-19 | Hecke Rank Constancy of the Staircase Product | [`2026-04-19-hecke-rank-constancy.tex`](2026-04-19-hecke-rank-constancy.tex) |
| 2026-04-17 | The H-Invariant Theorem for Staircase Products in the Symmetric Group (complete) | [`2026-04-17-h-invariant-complete.tex`](2026-04-17-h-invariant-complete.tex) |
| 2026-04-17 | Eigenvalue Positivity for the Hecke Staircase Element via Kazhdan–Lusztig Positivity (v2) | [`2026-04-17-eigenvalue-positivity-v2.tex`](2026-04-17-eigenvalue-positivity-v2.tex) |
| 2026-04-16 | Eigenvalue Positivity for the Hecke Staircase Element via Kazhdan–Lusztig Positivity | [`2026-04-16-eigenvalue-positivity.tex`](2026-04-16-eigenvalue-positivity.tex) |
| 2026-04-16 | The H-Invariant Theorem via Frobenius Reciprocity and V^H Injectivity | [`2026-04-16-frobenius-injectivity.tex`](2026-04-16-frobenius-injectivity.tex) |
| 2026-04-16 | The q-Determinant of the Right Multiplication Submatrix and the Image Basis Conjecture | [`2026-04-16-q-determinant.tex`](2026-04-16-q-determinant.tex) |
| 2026-04-16 | The Total Rank Formula and Image Basis Conjecture for the Staircase Product | [`2026-04-16-total-rank.tex`](2026-04-16-total-rank.tex) |
| 2026-04-15 | The Contraction Lemma and the Per-Irrep Upper Bound for the Staircase Product | [`2026-04-15-contraction-upper-bound.tex`](2026-04-15-contraction-upper-bound.tex) |
| 2026-04-15 | The Staircase Eigenspace Property: Computational Verification and Structural Analysis | [`2026-04-15-staircase-eigenspace.tex`](2026-04-15-staircase-eigenspace.tex) |
| 2026-04-15 | The H-Invariant Theorem for the Staircase Product: Partial Proof with Gap Analysis | [`2026-04-15-h-invariant-partial.tex`](2026-04-15-h-invariant-partial.tex) |
| 2026-04-14 | The H-Invariant Theorem for the Staircase Symmetrizing Product | [`2026-04-14-H-invariant-theorem.tex`](2026-04-14-H-invariant-theorem.tex) |
| 2026-04-14 | Rank of the Symmetrizing Product: Injectivity, Cascade Surjectivity, and n!/2^⌊n/2⌋ | [`2026-04-14-rank-injectivity.tex`](2026-04-14-rank-injectivity.tex) |
| 2026-04-13 | Rank of the Symmetrizing Product: Proof, Correction, and Computation | [`2026-04-13-rank-hierarchy.tex`](2026-04-13-rank-hierarchy.tex) |
| 2026-04-11 | The Multiplicity Bundle Theorem: Schur–Weyl Duality, Crystal Invisibility, and the Coboundary Hierarchy | [`2026-04-11-multiplicity-bundle.tex`](2026-04-11-multiplicity-bundle.tex) |
| 2026-04-11 | The Hecke Transition Algebra Theorem | [`2026-04-11-hecke-transition-algebra.tex`](2026-04-11-hecke-transition-algebra.tex) |

## q-Shifted pair theorems (T-system path)

| Date | Theorem | File |
| --- | --- | --- |
| 2026-04-19 | The Second q-Shifted Pair Theorem: det M_R^(5) · det M_R^(3) = q^|D_n| (det M_R^(4))^2 | [`2026-04-19-second-q-shifted-pair.tex`](2026-04-19-second-q-shifted-pair.tex) |
| 2026-04-19 | Closing the n=6 Base Case for the Second q-Shifted Pair: Degree Bound and Multi-Point Evaluation | [`2026-04-19-second-pair-base-case.tex`](2026-04-19-second-pair-base-case.tex) |
| 2026-04-18 | The First q-Shifted Pair Theorem: det M_R^(3) = (det M_R^(2))^2 in the Hecke Staircase | [`2026-04-18-first-q-shifted-pair.tex`](2026-04-18-first-q-shifted-pair.tex) |
| 2026-04-18 | The Rule B Block Decomposition of the Staircase Right-Multiplication Matrix | [`2026-04-18-block-decomposition-rule-b.tex`](2026-04-18-block-decomposition-rule-b.tex) |
| 2026-04-17 | Block-Multiplicative Structure of the Staircase Determinant (Partial — Computational Evidence) | [`2026-04-17-block-structure.tex`](2026-04-17-block-structure.tex) |

## Rank isolation and parabolic reduction

| Date | Theorem | File |
| --- | --- | --- |
| 2026-04-17 | Even-Block Gap Closure at k=4 for the Rank Isolation Lemma | [`2026-04-17-even-block-k4.tex`](2026-04-17-even-block-k4.tex) |
| 2026-04-16 | The Left-Two Lemma and the Even-Block Gap in the Rank Isolation Theorem | [`2026-04-16-left-two-lemma.tex`](2026-04-16-left-two-lemma.tex) |
| 2026-04-15 | The Rank Isolation Lemma and the H-Invariant Theorem | [`2026-04-15-rank-isolation.tex`](2026-04-15-rank-isolation.tex) |
| 2026-04-15 | Even-Block Gap Analysis for the H-Invariant Theorem | [`2026-04-15-even-block-k4.tex`](2026-04-15-even-block-k4.tex) |
| 2026-04-15 | Closing the Even-Block Gap in the H-Invariant Staircase Theorem | [`2026-04-15-even-block-gap.tex`](2026-04-15-even-block-gap.tex) |

## d=4 fiber-vanishing law (order-law / spectral thread)

The vanishing of `G_λ = ⟨s_λ, ψ^m⟩` at `ψ = h₂ + i·e₂` (the `d=4` root-of-unity fiber).
For two rows `λ=(2m−b,b)` the law is `G_λ=0 ⟺ (2,2)`, reduced to `Im G ≠ 0`.

| Date | Result | File |
| --- | --- | --- |
| 2026-06-09 | **General-`λ` `d=4`: reformulations + partial + precise gap (Prove).** Clean equivalence **`G_λ(i)=A_λ(i−1)`**, `A_λ(x)=⟨s_λ,(h₁²+x e₂)^m⟩∈ℤ_{≥0}[x]`, so **`G_λ(i)=0 ⟺ (x²+2x+2)\|A_λ(x)`** and `A_{(2,2)}=x²+2x+2` exactly (the unique shape whose `A_λ` *is* the minimal polynomial). Real form `Re G=S₀−S₁+2S₃`, `Im G=S₁−2S₂+2S₃`. Equivalent (but non-uniform) integer criterion `q(n)\|A_λ(n) ∀n`. Exact leading `(1+i)`-coeff `w=Σ_{j∈J*}o_j i^{(3μ−j)/2}≡\|J*\| mod π` ⟹ `\|J*\|` odd ⟹ `G≠0` (≥72%). **Verified `(2,2)` unique vanisher `n≤20`; `\|J*\|∈{1,2,4}` (power of 2).** Gap = uniform termination of interleaved even/odd cancellation tower (depth unbounded) | [`2026-06-09-d4-fiber-reformulations.md`](2026-06-09-d4-fiber-reformulations.md) ([scratch](2026-06-09-d4-fiber-evenness-scratch.md)) |
| 2026-06-09 | **Code:** arithmetic attack on `b ≡ 2,3`, frontier `b ≤ 119`. **`Q_b` is irreducible over ℚ for all `6 ≤ b ≤ 119`** (mod-`p` certificates) ⟹ no rational root ⟹ **law holds for those `b`** (deg `Q_b ≥ 3`). Coefficient Newton-polygon irreducibility (Hajir/Filaseta) is **dead** — no pure Eisenstein/Dumas prime exists; but `disc(Q_b)` is **non-square** with an odd-multiplicity prime > deg (a transposition). **Verdict: pursue irreducibility via a Frobenius `d`-cycle** (prime `p` with `Q_b` irreducible mod `p`), upgraded to `Gal=S_{⌊b/2⌋}` by Jordan. Job 2: `1+su+u²=(1+u+u²)+iu` ⟹ exact central-trinomial map (A002426); Prop 1 `g_b=C_b^{(−m)}(−s/2)` verified symbolically | [`2026-06-09-tworow-d4-code/`](2026-06-09-tworow-d4-code/) ([summary](2026-06-09-tworow-d4-code/FINDINGS-2026-06-09-SUMMARY.md), [Job 1](2026-06-09-tworow-d4-code/FINDINGS-job1-discriminant.md), [Job 4](2026-06-09-tworow-d4-code/FINDINGS-job4-irreducibility.md)) |
| 2026-06-08 | **`b ≡ 2,3 (mod 4)` — new closed form + methodological correction (not closed).** Theorem A: `G_b((2b−1)/2) = (−1)^b·C(2b,b)/4^b·(1−i)^b` (central binomial!), proved via Gauss quadratic `₂F₁` transformations; unifies the `4∣b` rational root with the `b≡2,3` non-vanishing. Engine: `g_b(m)=C_b^{(−m)}(−s/2)` Gegenbauer-in-parameter, 3-term recurrence. **Proves the obstruction is arithmetic not analytic** (real roots are in-range, so log-concavity dominance cannot close it). Frontier extended to `b≤70` (exact rational roots) | [`2026-06-08-tworow-d4-b23-halfinteger-central-binomial.md`](2026-06-08-tworow-d4-b23-halfinteger-central-binomial.md) ([`.tex`](2026-06-08-tworow-d4-b23-halfinteger-central-binomial.tex)) |
| 2026-06-08 | **Code:** verification scaffolding for the `b≡0` proof and the `b≡2,3` hunt. Job 1: mod-4 reduction **exact & `b`-uniform**; the 𝔽₂ closed form `I_b ≡ m·C(m−1,R) (mod 2)` from the **nilpotent truncation** `s²=(1+i)²≡0`; valuation proposition holds over `m∈[b,b+400]` with min-set uniformly `{1}`/`{1,2,3}`. Job 2: corrected Newton test (all slopes **non-integral**) gives **zero** certificates (monic/reciprocal/scaled all dead); **no** uniform prime, **no** finite covering set (evaporates by `b≤62`), **no** `b`-periodicity — `b≡2,3` is genuinely non-local | [`2026-06-08-tworow-d4-code/`](2026-06-08-tworow-d4-code/) ([𝔽₂](2026-06-08-tworow-d4-code/FINDINGS-fp2-recurrence.md), [b≡0](2026-06-08-tworow-d4-code/FINDINGS-job1-b0mod4.md), [odd-prime](2026-06-08-tworow-d4-code/FINDINGS-oddprime-newton.md)) |
| 2026-06-08 | **Two-row `d=4` law PROVED for all `b ≡ 0 (mod 4)`** — with `b ≡ 1` ⟹ **half of all `b`** done unconditionally. Same 2-adic identity `v₂(I_b) = v₂(∏(m−r)) − v₂(R!)`, but now by **parity counting**: only `j∈{1,2,3}` attain the minimal valuation, they tie iff `m` odd, and the count is always odd (1 or 3), so the sum keeps an odd leading 2-adic digit. Left: only `b ≡ 2,3` | [`2026-06-08-tworow-d4-b0mod4-proved.md`](2026-06-08-tworow-d4-b0mod4-proved.md) ([`.tex`](2026-06-08-tworow-d4-b0mod4-proved.tex)) |
| 2026-06-07 | **Code:** single-prime Newton-polygon route to `(♦)` is **structurally dead** (`Q̃_b` monic ⟹ always a slope-0 edge); the correct local test "no root mod p" gives a **uniform witness `p=2` for `b≡0,1 (mod 4)`** but **NO congruence-uniform prime for `b≡2,3 (mod 4)`** (apparent ones evaporate as the class grows). `(♦)` proved directly for `b≤40`. **TASK 2:** the runner-position hypothesis is **false**; residual is a function of `(core bead-count vector t_j, 4-quotient)`, not the quotient-box runner | [`2026-06-07-code-newton-runner/`](2026-06-07-code-newton-runner/) ([Newton](2026-06-07-code-newton-runner/FINDINGS-Qb-newton.md), [runner](2026-06-07-code-newton-runner/FINDINGS-4core-residual.md)) |
| 2026-06-07 | **Two-row `d=4` law PROVED for all `b ≡ 1 (mod 4)`** (infinite family) via the exact 2-adic identity `v₂(I_b(m)) = v₂(∏(m−r)) − v₂(R!)`; `b ≡ 0` reduced to one mod-4 lemma (verified `b≤40`); `b ≡ 2,3` shown to need an odd prime | [`2026-06-07-tworow-d4-b1mod4-proved.md`](2026-06-07-tworow-d4-b1mod4-proved.md) ([`.tex`](2026-06-07-tworow-d4-b1mod4-proved.tex)) |
| 2026-06-07 | The valuation `residual(λ)` is **not** a function of the 4-quotient (Outcome 3, counterexample at `n=10`); two-row `Q_b` has no classical-orthogonal-polynomial home; full law verified to `n≤20` | [`2026-06-07-residual-not-quotient-function.md`](2026-06-07-residual-not-quotient-function.md) ([data](2026-06-07-data/)) |
| 2026-06-06 | Reduction of the two-row `d=4` law to a no-rational-root statement `(♦)`; clean reduction `Im(A^m)=u·H_m`; forced-root lemma; 2-adic tower (levels 1–2) | [`2026-06-06-tworow-d4-no-rational-root.tex`](2026-06-06-tworow-d4-no-rational-root.tex) |
| 2026-06-06 | Two-row `d=4` law as imaginary-part non-vanishing (reduction + `b≤4` closed) | [`2026-06-06-tworow-d4-imaginary-reduction.md`](2026-06-06-tworow-d4-imaginary-reduction.md) |
| 2026-06-05 | Two-row `d=4` law as a continued fraction | [`2026-06-05-tworow-d4-continued-fraction.tex`](2026-06-05-tworow-d4-continued-fraction.tex) |

## Other structural results

| Date | Theorem | File |
| --- | --- | --- |
| 2026-04-12 | The Cactus Midpoint Theorem | [`2026-04-12-cactus-midpoint.tex`](2026-04-12-cactus-midpoint.tex) |
| 2026-04-10 | The sl_n Cactus Representation Theorem: Interval Reversals, Crystal Structures, and a Corrected Statement | [`2026-04-10-sln-cactus-representation.tex`](2026-04-10-sln-cactus-representation.tex) |
| 2026-04-10 | The Cactus Representation Theorem: Interval Reversals on Tensor Space at q=0 | [`2026-04-10-cactus-representation.tex`](2026-04-10-cactus-representation.tex) |
| 2026-04-10 | The Operator Independence Theorem: π_sort is not a function of R(u) | [`2026-04-10-operator-independence.tex`](2026-04-10-operator-independence.tex) |
| 2026-04-09 | Five-Vertex YBE Classification (notes) | [`2026-04-09-five-vertex-ybe-classification.md`](2026-04-09-five-vertex-ybe-classification.md) |

---

Generated 2026-04-19 by Clio.
