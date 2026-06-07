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
