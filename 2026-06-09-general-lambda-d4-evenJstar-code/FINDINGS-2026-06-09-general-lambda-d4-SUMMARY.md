# Code session 2026-06-09 — arming the even-|J\*| proof (general-λ, d=4)

Goal (from `state/CODE.md`): generate data to crack the **even-|J\*| crux** and close
`G_λ(i)=0 ⟺ λ=(2,2)` for all `λ ⊢ 2m`. Four jobs; all four delivered.

## Headline results

- **Job 1 — frontier m ≤ 14** (10 268 shapes). `|J*| ∈ {1,2,4}` only — *never* odd-≥3,
  never 6/8. **Zero odd-tie violations. The only `G=0` is `(2,2)`.**
  Engine cross-checked two independent ways (p-basis scalar product; master identity
  `⟨s_λ,ψ^m⟩`, 294/294). Files: `FINDINGS-job1-tie-census.md`, `results/{jstar-census,tie-shapes}.csv`.

- **Job 2 — THE PRIZE.** For every tie, **`J*` is an affine 2-adic box**
  `j₀ + {Σ_{a∈S}2ᵃ}` (generators distinct powers of two), so `|J*|=2^{|S|}`.
  Verified **4318/4318 ties to m≤14**. The smallest-generator toggle is a
  **fixed-point-free, val-preserving involution (1624/1624)** ⟹ even count.
  Crisp lever for PROVE: the **exact valuation step law**
  `v₂(C(m,j)M_j) − v₂(C(m,j+2ᵃ)M_{j+2ᵃ}) = 2^{a−1}` on every pair (**2398/2398**).
  The box is driven by `v₂(M_j)`, **not** the binomial part (binomial-only argmin: 0/1624).
  File: `FINDINGS-job2-mechanism.md`, `results/job2-box-structure.csv`.

- **Job 3 — Step 2 verdict.** Every tie ≠ (2,2) has **finite** `v_π(G)`; the surviving
  order `d=v_π(G)−μ` ranges 1..22. **The naive reduction conjecture is FALSE**: in
  918/1624 ties the survivor needs the `val>μ` tail, not J\* alone. So Step 2 cannot be
  closed by re-deflating and re-minimising on J\*. File: `FINDINGS-job3-survivor.md`.

- **Job 4 — Hidaka–Itoh.** My `G_λ` is **not** the fake degree (9/294 at d=4). The HI
  `{−1,0,1}` gate (≤2 odd prime factors) and my richness gate (`d≡2 mod4`) are **distinct
  arithmetic** — e.g. `d=210` is rich but fails HI. Orthogonal; no bearing on the crux.
  File: `FINDINGS-job4-hidaka-itoh.md`.

## Which mechanism the data favours

The evenness is **fully explained structurally**: `J*` is a 2-adic boolean subcube, and
the smallest-generator toggle is a canonical fixed-point-free involution. This reduces the
*entire* even-|J\*| phenomenon to **one 2-adic identity** — the step law
`Δv₂(C(m,j)M_j)=2^{a−1}` — which in turn reduces (on binomial-flat pairs) to a clean law
on `v₂(M_j)` alone. The recommended PROVE target is a **Kummer/Lucas formula for
`v₂(M_j)`**, `M_j=⟨s_λ,h₁^{2(m−j)}e₂ʲ⟩`: it would yield the step law, the box, the
involution, and hence leading-order cancellation for every non-(2,2) shape in one stroke.

The remaining genuine gap is **Step 2's depth** (Job 3): non-vanishing past the leading
order is empirically universal but its valuation is governed by the full tail, so it needs
an invariant beyond "re-run the Newton polygon on J\*".
