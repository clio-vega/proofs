# FINDINGS — Job 1: even-|J\*| tie census (general-λ, d=4)

**Date:** 2026-06-09 (code session)
**Scripts:** `job1_tie_census.py` (+ `verify_master.py` cross-check)
**Data:** `results/jstar-census.csv` (all shapes), `results/tie-shapes.csv` (ties only)

## Setup (exact, ℤ[i]; π = 1+i, v_π(N) = 2 v₂(N) for N∈ℤ)

`G_λ(i) = Σ_{j=0}^m C(m,j) iʲ (1+i)ʲ M_j`, `M_j = ⟨s_λ, h₁^{2(m−j)} e₂ʲ⟩ ∈ ℤ_{≥0}`,
`val(j) = j + 2 v₂(C(m,j) M_j)`, `J* = argmin val`, `μ = min val`.

**New fast engine (no Sage).** Using `h₁=p₁`, `e₂=(p₁²−p₂)/2`:
> `M_j = 2^{−j} Σ_{b=0}^j C(j,b)(−1)ᵇ χ^λ(2ᵇ 1^{2m−2b})`,

so every `M_j` follows from Murnaghan–Nakayama characters on `(2ᵇ,1^{a})`-classes.
Census to **m ≤ 14 runs in ~4 s**. Cross-checked two ways:
(i) shortcut `M_j` == p-basis scalar product `⟨s_λ, h₁^{2(m−j)}e₂ʲ⟩` (symfunc);
(ii) `G_λ(i)` from the M_j expansion == `⟨s_λ, ψ^m⟩`, `ψ = h₁²+i(1+i)e₂`
(294/294 shapes m≤7, `(2,2)=0` both ways).

## Results

**Frontier m ≤ 14** (10 268 shapes ⊢ 2..28; the tie CSV stores m≤12 detail, census m≤14):

| m | shapes | unique-min (|J\*|=1) | ties (|J\*|≥2) | odd-tie violations |
|---|--------|------|------|------|
| 7 | 135 | 98 | 37 | 0 |
| 10 | 627 | 360 | 267 | 0 |
| 12 | 1575 | 864 | 711 | 0 |
| 14 | 3718 | 1955 | 1763 | 0 |

- **even-|J\*| confirmed for every λ ≠ (2,2), all m ≤ 14.** Zero odd-tie violations.
- **|J\*| distribution (all shapes, all m≤14): {1: 5950, 2: 3529, 4: 789}.**
  Sizes **6 and 8 never appear** — ties are size 2 or 4 only, even at m=14.
- **No G=0 outlier:** the unique total-vanisher is `(2,2)` (`v_π(G)=∞`); every other
  shape — including all 4318 ties — has finite `v_π(G)`.

## Reading

`|J*|` is always a **power of 2** (1, 2, or 4), never an odd number > 1 and never
6/8. The unique-min case `|J*|=1` is handled directly by the (1+i)-adic Newton
polygon; the whole remaining problem is the even ties, whose mechanism is Job 2.
