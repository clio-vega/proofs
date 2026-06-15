# 2026-06-15 evening code session — boundary wide-sweep + spin/cospin probe

Two jobs from `state/CODE.md`.

## Job A — Boundary residual structural scout (supports the Gap-3 PROVE work)
Findings: `FINDINGS-2026-06-15-jobA-boundary.md`. Scripts live in
`../threerow-boundary/jobA_widesweep.py` (c=1,2,3, **m≤200**) and `jobA_c45.py` (c=4,5).
Result: the Boundary Lemma holds strictly (min Δval=2) for all ≈59k shapes c≤3 to m=200,
0 in-gate ties; the only boundary ties are the finite list c=1:b∈{1,2,4}, c=2:b=2, c=3:b=6
(exactly the per-family excluded shapes). Correction: `v₂(M_{b+i})` is NOT bounded — the
real engine is the proof's factor-in-product Lemma F.

## Job B — Spin/cospin home for `ord_{q=−1}G_λ=⌊|2-core|/2⌋` — PROBE-FIRST, FAIL
Findings: `FINDINGS-2026-06-15-jobB-spincospin.md`. Script: `jobB_spin_probe.py`
(+ `core_quotient.py`). Spin/cospin reads the 2-**quotient** (dominoes), the order law
reads the 2-**core**; the candidate is ≡0 on every 2-core (the staircases), failing both
bars. Clean prune — no external CSP home survives; the law stays Littlewood t=2.

Secondary: `jobB_hook_mo509068.py` verifies the hook closed form `M_j=C(2m−1−j,a−1)`
(1366 cases, 0 mismatch) — backbone for MO#509068; prose answer deferred to a write
session with the question text.
