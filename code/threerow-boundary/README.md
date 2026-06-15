# threerow-boundary — Boundary Lemma (Gap 3) verification

Verifies the Boundary Lemma for the three-row d=4 family (proof:
`projects/proofs/2026-06-15-boundary-lemma-threerow.md`).

- `mn.py`          — Murnaghan–Nakayama M_j (shared, from threerow-c1)
- `boundary_forms.py` — fits boundary M_{b+i} closed forms (Lemmas T, H), c=1..5
- `general_top.py` — verifies hook formula (H), top M_{b+c} (T), clean Δ(b+c), c=1..5
- `c1_clean.py` / `c1_bdry.py` / `c1_fullproof.py` — c=1 Prop 1 + boundary lemma conclusion
- `c2_check.py` / `c2_final.py` / `c2_prod.py` — c=2 Δ(b+1),Δ(b+2) identities + full boundary check
- `c3_aeventop.py` / `c3_status.py` — c=3 a-even top identity; full c=3 boundary scan
- `lemma_only.py`  — Lemma F (factor-in-product) pure number-theory check
- `jobA_widesweep.py` — **(2026-06-15 evening)** gated val-gap sweep, c=1,2,3, **m≤200**
  (pure-Python closed forms + MN cross-check): boundary loses by min Δ=2, 0 in-gate ties;
  the only boundary ties are the finite list c=1:b∈{1,2,4}, c=2:b=2, c=3:b=6. Also the
  boundedness probe showing `v₂(M_{b+i})` is NOT bounded (the real engine is Lemma F).
- `jobA_c45.py` — c=4,5 confirmation via the exact 3-row Jacobi–Trudi determinant (min Δ≥8).

All checks: 0 mismatches / 0 violations in stated ranges.
See `../FINDINGS-2026-06-15-jobA-boundary.md` for the structural write-up.
