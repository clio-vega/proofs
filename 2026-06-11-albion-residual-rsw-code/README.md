# 2026-06-11 — Albion z-asymmetric residual (Job A) + graded RSW (Job B)

Two decisive, never-run probes from the d=4 fiber-vanishing thread. **Both clean negatives.**

## Job A — `jobA_albion.py`
Tests whether Albion 2501.18520's **z-asymmetric 4-core/4-quotient** datum
`(charge vector κ₄(4-core), ordered 4-quotient)` supplies a closed form for the residual
`residual(λ) = v_π(G_λ) − v₂(f^λ) − v_π(G_core) + v₂(f_core)`.

- `κ₄` is a bijective invariant of the 4-core ⟹ `(κ₄, quotient)` = complete invariant of λ, so
  "residual = f(datum)" is a vacuous tautology (corrects the 06-07 `(t_j, quotient)` framing).
- Decisive ladder (n≤18): residual is **not** a function of `(charge, sizes)`, `(charge, size
  multiset)`, or `(charge, per-runner v₂f)`. The residual reads `v₂`-hook structure *inside* the
  quotient + a charge-mediated cross-box interaction. Killer witness: trivial charge, 3 boxes on
  runner 2, `(3,)→0`, `(1,1,1)→0`, `(2,1)→1`.
- Run: `PYTHONPATH=. python3 jobA_albion.py` → `results/jobA-albion-charge.csv`.

## Job B — `jobB_graded_rsw.py`
The graded `G_λ(q) = Σ_{T∈SYT(λ)} q^{s(T)}`, `s(T)=Σ_{i∈Des(T)} w_i`, `w_i=2i−1` if `n−i` odd else 0.

- **STEP 0 (the win):** `Σ_T q^{s(T)}` equals the master-identity polynomial
  `⟨s_λ, h₁^e ∏_{k∈A}(h₂+e₂t^{2k−1})⟩` **exactly, 271/271 for n≤12** — the SYT and symmetric-function
  representations are one object.
- STEP 1–2: `G_λ(q)` is **not** a q-shift of any principal specialisation `s_λ(1,q,…,q^{k−1})` for
  d≥3 (the RSW escape hatch is closed).
- STEP 3: vanishers transparent — `G_{(2,2)}=q⁶+1` dies at ζ₄, `G_{(3,1)}=q⁵+q+1` and
  `G_{(2,1,1)}=q⁶+q⁵+q` die at ζ₃.
- Run: `PYTHONPATH=. python3 jobB_graded_rsw.py 9` → `results/jobB-graded-rsw.csv`.

## Files
`jobA_albion.py`, `jobB_graded_rsw.py` (probes); `fiber_engine.py`, `core.py`, `symfunc.py`,
`core_quotient.py`, `quotient.py`, `job1_tie_census.py` (engines); `FINDINGS-2026-06-11-jobAB-albion-rsw.md`
(full writeup); `results/` (data). Reproducible with `PYTHONPATH=.`.
