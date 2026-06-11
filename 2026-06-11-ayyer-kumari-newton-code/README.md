# 2026-06-11 — Ayyer–Kumari test + sharp Newton polygon (code session)

d=4 fiber `G_λ(i)=0 ⟺ λ=(2,2)`, even-|J*| crux. See `FINDINGS-2026-06-11.md` for the full writeup.

## Headline results

- **Job A — Ayyer–Kumari hook-Schur factorisation (arXiv:2501.00275): DEAD as an H1/H2 lever.**
  The genuine rectangular-class 2-quotient factorisation `|χ^λ(2^m)|=C(m;|q⁰|,|q¹|)f^{q⁰}f^{q¹}`
  holds (433/433) but only pins the *leading* coefficient of `Φ(z)=Σ_k C(m,k)χ_k z^k`. `Φ` does
  **not** split into linear factors (0/82 tie shapes), the factorisation does not determine the
  interior `v₂(M_j)` defining the box, and vanishing needs the **irreducible quadratic** `(z²+1)`
  — carried uniquely by `(2,2)`, itself a 4-core. Confirms/sharpens the prior Route-B pruning.

- **Job B — sharp `J*` Newton polygon (m ≤ 8, box census m ≤ 10).**
  - **H1 (single parity) is TRIVIAL:** `val(j)=j+2v₂(·)≡j (mod 2)` (0 exceptions / 3832 pairs).
  - **H2 (affine 2-adic box) robust:** 172/172 (m≤8), 560/560 (m≤10); `|J*|∈{1,2,4}`; generator
    offsets even (`{2},{4},{8},{2,4},{2,8}`), refining `S⊆{0,1}` to `S⊆{1,2,3,…}`.
  - S is **not** a function of the 4-quotient (19/59 ambiguous); no determinant in {m, v₂(m)};
    `v₂(M)` jumps non-constant. H2's mechanism = the e₂-mod-2 layer wall (unchanged).

## Files

| file | what |
|---|---|
| `quotient.py` | t-core / t-quotient (abacus); verifies rectangular character formula |
| `jobA_ayyer_kumari.py` | Job A: T1 factorisation, T2 linear-split, T3 valuation lever, T4 (2,2) |
| `jobB_newton_polygon.py` | Job B: H1 triviality, H2 box, S-vs-quotient, anchor; writes the CSV |
| `jobB_probe_S.py` | probe for a determinant of the generator set S (none clean found) |
| `job1_tie_census.py` | core engine: MN characters, `M_j`, `val`, `J*` (p-basis cross-checked) |
| `job_box_detail.py`, `job_hierarchy.py`, `symfunc.py` | supporting machinery |
| `results/jobB-sharp-Jstar.csv` | 172 ties m≤8 with j₀, generators, 4-core, 4-quotient |
| `results/job*_run.txt` | captured run logs |

Run: `python3 jobA_ayyer_kumari.py`, `python3 jobB_newton_polygon.py`. Self-check:
`python3 -c "import job1_tie_census as J; J.selfcheck()"`.
