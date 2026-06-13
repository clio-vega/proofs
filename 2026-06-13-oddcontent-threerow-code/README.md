# 2026-06-13 code session — odd-content probe + three-row M_j

Two decisive computational jobs from the order-law / `Z_λ` spectral thread.

## Job A — the odd-content probe (`2026-06-13-jobA-oddcontent.py`)

Tests whether the proved grade order law `ord_{q=−1} G_λ = ⌊|2-core(λ)|/2⌋` equals the
odd-content count `OC(λ) = #{cells : content odd}` that the Armon–Swanson content product
`∏(1+t q^{c})` reads off at `q=−1`.

**Verdict: MISMATCH (decisive).** Fails already at the 2-core δ₃ (`ord=3` vs `OC=2`). Closed
forms: `OC(δ_k)=⌊k/2⌋(⌊k/2⌋+1)`, `ord(δ_k)=⌊k(k+1)/4⌋`. Cross-check: SYT-sum `ord` matches the
closed form on all n≤9 (machinery sound). Kills the "smaj content-product → order law" route;
leaves the multiset identity `{d_j}={s(T)}` untouched. See `FINDINGS-...-jobA-oddcontent.md`.

## Job B — three-row M_j + sharp J* box (`2026-06-13-jobB-threerow.py`)

`M_j = ⟨s_λ, h₁^{2(m−j)} e₂^j⟩` for three-row `λ=(a,b,c)`, `m≤12`.

- **Exact closed form = signed S₃ Jacobi–Trudi determinant of trinomials** (0 failures / 2226
  pairs vs the vertical-2-strip chain engine). This is the PROVE handle.
- **Dimension conjecture `M_j = f^{(a−j,b−j,c−j)}` is FALSE** (break 1592/2226; fails even when
  the shift is a partition). The two-row shortcut does not lift.
- **Sharp three-row J* box = `{0,2,4,6}` (generators 2,4)**, not `{0,2}`. `|J*|∈{1,2,4}`, all
  single-parity XOR 2-adic boxes. First `|J*|=4` at the scaled staircase `(9,6,3)`.

See `FINDINGS-...-jobB-threerow.md`.

## Files
- `job_jstar_engine.py` — shared M_j / val / J* engine (vertical-2-strip chain model).
- `2026-06-13-jobA-oddcontent.py`, `2026-06-13-jobB-threerow.py` — the two jobs (self-checking).

Pure Python (no SageMath in container). Run: `python3 2026-06-13-jobB-threerow.py`.
