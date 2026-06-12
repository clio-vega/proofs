# 2026-06-12 code session — `|J*|`-even SYT-chain model + Izergin probe

Data-gathering for the PROVE target *"`|J*|` even on tie shapes via the vertical-2-strip chain model
of `M_j`"*, plus a cheap seed-native Izergin sanity probe.

`M_j(λ) = ⟨s_λ, h₁^{2(m−j)} e₂^j⟩ = Σ_μ N^{(j)}_{λ→μ} f^μ` (vertical-2-strip removal chains).

## Files
- `job_jstar_engine.py` — shared engine: `M_j` two ways (build-up Pieri / chain removal), partitions,
  hook-length `f^λ`, 2-adic helpers, `val`/`J*` analysis.
- `job_jstar_crosscheck.py` — the two `M_j` computations agree (4123 `(λ,j)` pairs, `m ≤ 8`).
- `job_jstar_table.py` — full `|J*|` table, `m ≤ 16`. **Run: `python3 job_jstar_table.py 16`.**
  Writes `results-jstar-ties.json` (10 727 tie shapes; ~3 MB, regenerable, not committed).
- `job_jstar_newton.py` — 2-adic Newton polygon of the half-poly; slope-`(−1)` edge data.
- `job_jstar_involution.py` — fpf-involution search over the chain model.
- `job_izergin_probe.py` — Job B: `ord_{q=-1} G_λ = ⌊|2-core|/2⌋` + Bezoutian rank-drop.
- `FINDINGS-2026-06-12-jstar-even-syt-model.md` — full write-up.

## Headlines
- `|J*|` even whenever `≥2`, **no odd `|J*|≥3`**, all `m ≤ 16` (wall holds). `|J*|=8` debuts at `m=15`.
- All 10 727 tie `J*` are XOR-closed 2-adic boxes; edge reduced-coeffs all `≡1 (mod 2)`;
  edge poly `≡ (1+y)^{Σ gₐ} (mod 2)`.
- **The "fpf chain involution" goal is mis-aimed** — `|J*|` is a run-length, not a cardinality.
- **Positive lead:** `v₂(M_j)` is carried *entirely* by the **empty-2-core (2-quotient) sector** of
  the chain ends (2064/2064, `m ≤ 11`); the rest is 2-adically negligible. PROVE should aim there.
- Job B: graded law re-confirmed `n ≤ 14`; `ord_{q=-1}` is a (generic) Bezoutian rank-drop, **no
  cheap Izergin honeycomb determinant**; the `M(q²)`-corank candidate is flagged for the future.
