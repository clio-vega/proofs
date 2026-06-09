# General-λ d=4: arming the even-|J\*| proof (2026-06-09 code session)

Computational support for closing `G_λ(i)=0 ⟺ λ=(2,2)` for all `λ ⊢ 2m`, where
`G_λ(i)=⟨s_λ,ψ^m⟩`, `ψ = h₁² + i(1+i)e₂`. Via the (1+i)-adic Newton polygon the
210/292 unique-min shapes are immediate; the open problem is the **even ties**.

## Machinery (exact, ℤ[i]; π=1+i)

`G_λ(i) = Σ_{j=0}^m C(m,j) iʲ (1+i)ʲ M_j`, `M_j = ⟨s_λ, h₁^{2(m−j)} e₂ʲ⟩ ∈ ℤ_{≥0}`,
`val(j)=j+2v₂(C(m,j)M_j)`, `J*=argmin val`. No Sage: with `h₁=p₁`, `e₂=(p₁²−p₂)/2`,
`M_j = 2^{−j} Σ_b C(j,b)(−1)ᵇ χ^λ(2ᵇ1^{2m−2b})` from Murnaghan–Nakayama characters.

## Scripts

- `symfunc.py` — exact power-sum symmetric-function engine (MN characters, scalar product).
- `job1_tie_census.py [MMAX]` — census; writes `results/jstar-census.csv`, `results/tie-shapes.csv`.
- `verify_master.py` — cross-check: M_j expansion == `⟨s_λ,ψ^m⟩` (294/294, m≤7).
- `job2_mechanism.py` — J\* is an affine 2-adic box; involution; `results/job2-box-structure.csv`.
- `job2c_v2Mj.py` — shows the box living in `v₂(M_j)` on representative ties.
- `job2d_steplaw.py` — the exact valuation step law on involution pairs (2398/2398).
- `job3_survivor.py` — surviving order `v_π(G)−μ`; `results/job3-survivor.csv`.
- `job4_hidaka_itoh.py` — Hidaka–Itoh (2403.10817) gate comparison.

## Results (frontier m ≤ 14, 10 268 shapes)

- `|J*| ∈ {1,2,4}` only; **zero odd-tie violations**; **only `G=0` is `(2,2)`**.
- **Mechanism:** every tie has `J* = j₀ + {Σ_{a∈S}2ᵃ}`, an affine 2-adic box (4318/4318);
  the smallest-generator toggle is a fixed-point-free, val-preserving **involution** ⟹ even.
- **Lever for PROVE:** exact step law `v₂(C(m,j)M_j)−v₂(C(m,j+2ᵃ)M_{j+2ᵃ})=2^{a−1}`.
- **Step 2:** every tie ≠ (2,2) has finite `v_π(G)`; but the survivor needs the `val>μ`
  tail in 918/1624 cases — the naive J\*-only reduction conjecture is **false**.

See `FINDINGS-2026-06-09-general-lambda-d4-SUMMARY.md` for the full writeup.
