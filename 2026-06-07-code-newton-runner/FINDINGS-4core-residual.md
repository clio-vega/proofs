# FINDINGS — is `residual(λ)` a function of the 4-quotient alone?

**Date:** 2026-06-07 (prove session, computational part)
**Verdict: OUTCOME 3 — `residual(λ)` is NOT a function of the 4-quotient.**
A clean counterexample exists at `n = 10`. The residual genuinely depends on the **4-core**, not
just the 4-quotient. The valuation-level decomposition as stated therefore does not isolate
quotient-dependence and yields **no clean proof of Gap A**.

## Setup
For `λ ⊢ 2m`, `G_λ(i) = ⟨s_λ, ψ^m⟩`, `ψ = h₂ + i e₂`, computed exactly in `ℤ[i]` by iterated Pieri
(reusing `dfour_imaginary_reduction.py` / `core.psi_power_schur`). With `π = 1+i`,
`v_π(z) = v₂(N(z)) = v₂(Re² + Im²)` (since `π̄ = −iπ` is a unit multiple of `π`, so `v_π(z̄)=v_π(z)`
and `2v_π(z) = v_π(N(z)) = 2v₂(N(z))`).

James–Kerber 4-core / ordered 4-quotient via the 4-runner abacus on a bead-set of size `r ≡ 0 (mod 4)`.
Decomposition (PROVE.md):
```
v_π(G_λ) = v₂(f^λ) + [ v_π(G_{4core}) − v₂(f^{4core}) ] + residual(λ)
```
so `residual(λ) := v_π(G_λ) − v₂(f^λ) − v_π(G_{4core}) + v₂(f^{4core})`.
Range scanned: all `λ ⊢ 2m`, `m = 1..10` (`n ≤ 20`), **1536 shapes** (full data in
`results/residual-vs-4quotient.csv`).

## Result 1 (the verdict) — residual is NOT quotient-only
Grouping the 1536 shapes by **ordered** 4-quotient gives 416 groups, of which **58 are
non-constant** in residual; by **unordered** 4-quotient, 51 groups, **10 non-constant**.

**Smallest counterexample (`n = 10`).** Both shapes have ordered 4-quotient
`((), (), (1,), ())` (a single box on runner 2):

| λ | G_λ(i) | v_π | v₂(f^λ) | 4-core | v_π(G_core) | v₂(f_core) | **residual** |
|---|--------|-----|---------|--------|-------------|------------|-------------|
| `(1,1,1,1,1,1)` | `(0,−1)` | 0 | 0 | `(1,1)` | 0 | 0 | **0** |
| `(4,4,2)`       | `(−44,−4)`| 5 | 2 | `(4,1,1)` | 3 | 1 | **1** |

Same 4-quotient, different 4-core, **different residual (0 vs 1)**. ∎
(Three more `n=10` pairs exist, one for each single-box runner.)

## Result 2 — the dependence is genuinely on the 4-CORE (not a convention artifact)
- **Convention robustness:** the ordered 4-quotient is **stable** under `r → r+4` (`r ≡ 0 mod 4`),
  verified for all `λ`, `n ≤ 14`. So the grouping is canonical and the counterexamples are real.
- **`(core, quotient)` is a complete invariant** (James–Kerber bijection `λ ↔ (4core, ordered 4quot)`),
  so "residual is a function of `(core,quotient)`" is automatically true (0 non-constant of 1373) and
  carries **no information** — the meaningful negative is that the quotient *alone* fails.
- **Almost-but-not** a function of core valuations: grouping by `(v_π(G_core), v₂(f_core), quotient)`
  leaves only **12/787** non-constant groups, all at the single large-core class `(v_π=15, v₂=4)`.
  So the core enters mostly through its own valuations `(v_π(G_core), v₂(f_core))` — but not exactly.
- **Runner-sensitivity:** at fixed single-box quotient, residual depends on *which runner* the box
  sits on for some cores (e.g. core `(5,2,1,1,1)`: residuals 2 / 4 / 4 across runners 2 / 3 / 0).
  This is the signature of a genuine **core×quotient interaction term**, absent from the stated
  decomposition.

## Result 3 — the punchline mechanism SURVIVES (good news for the law itself)
- **Unique vanisher:** the only `λ ⊢ 2m` with `G_λ(i) = 0` for `n ≤ 20` is `λ = (2,2)`.
  (Verified hand-side: `⟨s_{22},ψ²⟩ = ⟨s_{22},h₂²⟩ + 2i⟨s_{22},h₂e₂⟩ − ⟨s_{22},e₂²⟩ = 1+0−1 = 0`.)
  → the **full** d=4 fiber law `G_λ(i)=0 ⟺ λ=(2,2)` holds for ALL shapes to `n ≤ 20`.
- **`core=(2,2)` family:** among 163 shapes with 4-core `(2,2)`, `v_π(G_λ)` is finite **⟺** the
  4-quotient is nonempty — **0 violations**. So the lone infinite-π-depth even 4-core is `(2,2)`
  itself, and any nonempty quotient lifts it to finite depth, exactly the claimed mechanism.

## Consequence for the programme
The route "prove `residual = f(4-quotient)` ⟹ Gap A" is **dead**: residual is irreducibly
core-dependent. The right object is a **direct formula for `v_π(G_λ)`** that keeps a core×quotient
**interaction/runner term** — the decomposition must not pretend the core's contribution is just the
additive bracket `[v_π(G_core) − v₂(f_core)]`. The near-formula `residual ≈ f(v_π(G_core), v₂(f_core),
quotient)` (775/787 exact) is the natural next target: characterise the 12 exceptions, then prove the
interaction term is **finite** for nonempty quotient — which already suffices for `G_λ=0 ⟺ (2,2)`.

Scripts: `scratch/2026-06-07-prove/residual.py`, `residual_analyze.py`, `residual_depend.py`.
Data: `results/residual-vs-4quotient.csv`.

---

# APPENDIX (2026-06-07, code session) — the runner hypothesis is FALSE; the missing variable is the core's bead-count vector

**Scripts:** `task2_runner.py`, `task2_correction.py`. **Data:** `results/residual-by-runner.csv`.
**Range:** all `λ ⊢ 2m`, `m ≤ 11` (n ≤ 22, **2538 shapes**) — extends the previous m≤10.

## The hypothesis tested
CODE.md TASK 2: the near-formula `residual = f(v_π G_core, v₂ f_core, 4-quotient)` held
775/787; at core `(5,2,1,1,1)` a single quotient box gives residuals 2/4/4 split by **which
abacus runner** it sits on. Hypothesis: *the runner position is the missing variable.*

## Result — runner position does NOT resolve the exceptions
A ladder of grouping keys (non-constant groups = where residual fails to be a function):

| key | #groups | non-constant |
|-----|---------|--------------|
| K0 `(v_πG_core, v₂f_core, UNORDERED quotient)` | 190 | **12** |
| K1 `(v_πG_core, v₂f_core, ORDERED quotient)` = K0 + which-runner | 1124 | **12** |
| K2 `K1 + runner-multiset` (TASK 2's literal key) | 1124 | **12** |
| **K3 `K1 + core bead-count vector t_j`** | 2375 | **0** |
| K4 `(4-core, ORDERED quotient)` (complete invariant, sanity) | 2375 | 0 |

- Going K0→K1 (adding the runner the box sits on) resolves **0 of 12**. Adding the
  runner-**multiset** (K2) resolves **0 of 12**. **The runner is NOT the missing variable.**
- The runner the box sits on *is* already encoded in the **ordered** quotient, and residual
  does vary with it at fixed core — but that variation was never the source of the 12
  exceptions. The 12 exceptions are **core collisions**: pairs of *different* 4-cores sharing
  the same `(v_πG_core, v₂f_core)` and the same ordered quotient, with different residual.

## The actual missing variable — the core's per-runner bead-count vector
Adding `t_j` = the number of beads on runner `j` of the 4-runner abacus (James–Kerber
charge of the core) makes residual a function: **K3 → 0/2375 non-constant**. Using the
**canonical, r-independent** bead-count vector `t_j(core)` (computed at `r = 4⌈ℓ/4⌉`) it is
likewise **0/2375**, and even `(t_j(core), ordered quotient)` *alone* determines residual
(`results/task2_canonical_check.txt`). So:

> **residual is a function of `(t_j(core), ordered 4-quotient)`** — the core enters only
> through its bead-count (runner-charge) vector. `(v_πG_core, v₂f_core)` is too lossy a
> summary of the core: it collides on exactly the `(v_π=15, v₂=4)` class, the home of all
> 12 exceptions (cores `(5,2,2,1,1,1)` with `t=(4,2,0,2)` and `(6,3,1,1,1)` with `t=(2,4,2,0)`).

## Structure of the interaction (no clean universal closed form yet)
At fixed core, a single-box quotient's residual splits the four runners into two pairs:

| core | `t_j(core)` | residual by runner `{0,1,2,3}` |
|------|-------------|--------------------------------|
| `(5,2,1,1,1)`   | `(4,2,2,0)` | `{4, 2, 2, 4}` |
| `(5,2,2,1,1,1)` | `(4,2,0,2)` | `{4, 2, 4, 2}` |
| `(6,3,1,1,1)`   | `(2,4,2,0)` | `{2, 4, 2, 4}` |
| `(6,3,2,1,1,1)` | `(2,4,0,2)` | `{7, 4, 4, 7}` |

For the three smaller cores the high-residual runners are exactly those with `t_j ≡ 0 (mod 4)`
(i.e. `t_j ∈ {0,4}`), low on `t_j ≡ 2`. The largest core `(6,3,2,1,1,1)` **breaks** this
(high on `t_j=2` runners, and magnitudes 7/4 not 4/2), so the closed form is **not** simply
`t_j mod 4`: the *magnitude* of Δ needs finer data (bead **positions**, not just counts).
We therefore report Δ as **structured but not yet closed-form**, with the certified fact that
`(t_j(core), quotient)` is a complete determinant.

## Bottom line for Gap A / the programme
- **TASK 2 hypothesis (runner position) is refuted.**
- The decomposition must carry a **core×quotient interaction through the core's bead-count
  vector `t_j(core)`**, not the quotient-box runner. The right object is
  `v_π(G_λ) = v₂(f^λ) + Φ(t_j(core), 4-quotient)` with `Φ` finite for nonempty quotient —
  which already suffices for the punchline `G_λ = 0 ⟺ λ = (2,2)`.
- Unique vanisher `(2,2)` and the `core=(2,2)` finite-⟺-nonempty-quotient mechanism both
  re-confirmed at the extended range `n ≤ 22` (0 violations).
