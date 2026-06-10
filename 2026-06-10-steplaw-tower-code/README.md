# Reproducing scripts — 2026-06-10 step-law data + Step-2 deflation tower

Pure-Python (Murnaghan–Nakayama in `job1_tie_census.py`; no Sage). Run any script from this
directory. Two jobs, serving the two open fronts of the general-λ d=4 fiber order law
`ord_{x=q²} Z_λ = τ(τ+1)/2`.

Recall: `G_λ(i)=Σ_j C(m,j) iʲ(1+i)ʲ M_j`, `M_j=⟨s_λ,h₁^{2(m−j)}e₂ʲ⟩`, `M₀=f^λ`;
`val(j)=j+2v₂(C(m,j)M_j)`, `J*=argmin val`, `μ=min val`, `v_π(N)=2v₂(N)` for `N∈ℤ`.

## Job A — a 2-adic handle on `v₂(M_j)` (arms PROVE)

- `jobA_v2Mj.py [MMAX]` — builds `results/v2Mj-table.csv` (one row per `(λ,j)`, `m≤MMAX`, default 12),
  tests global closed-form candidates, isolates the **step-law decomposition** on `J*` pairs, and
  re-confirms `S⊆{1,2,3}`.
- `jobA_probe.py`, `jobA_probe2.py` — determination tests (which invariants pin `v₂(M_j)`) and the
  `v₂(f^λ)` 2-quotient factorisation.

**Result:** no clean global Kummer/Lucas/additive form for `v₂(M_j)` (the 4-core obstructs any
quotient-only form). But the step law `Δv₂(M)=−2^{a−1}−Δv₂(bin)` is **exact (6685/6685 pairs to
m=14)**, its binomial part is Kummer (provable), and `v₂(f^λ)` factors cleanly through the 2-quotient.
`S⊆{1,2,3}` re-confirmed to m=14. See `FINDINGS-2026-06-10-jobA-v2Mj.md`.

## Job B — the Step-2 deflation / depth tower

- `jobB_tower.py [MMAX]` — carry-honest deflation trace per tie; builds `results/step2-tower.csv`
  (`d=v_π(G)−μ`, within-J\* depth `d1`, surviving-order index set, 4-core/4-quotient).
- `jobB_probe.py`, `jobB_probe2.py` — predictor tests for `d`; stress cases `(2,2)`, `(3,3,1,1)`,
  largest-d tie `(6,6,4,4,2,2)` (d=22).

**Result:** the optimistic "deflate to a second-level box" picture is **refuted** — surviving-order
set is a box only 29% of the time, and `d` is governed by the full λ (not the 4-core/quotient or any
simple quotient statistic). The tail is a genuine **multi-level π-adic cascade**. One clean shard:
core `(2,2)`, odd `|quotient|` ⟹ d=2. See `FINDINGS-2026-06-10-jobB-tower.md`.

## Shared engine

`symfunc.py` (p-basis cross-check), `job1_tie_census.py` (MN → `M_j`, `G_λ(i)`, `v_π`),
`job2_mechanism.py` (`val_pieces`, `jstar_of`, `box_generators`), `core_quotient.py` (NEW —
d-core/d-quotient via the abacus; self-tests `python3 core_quotient.py`).
