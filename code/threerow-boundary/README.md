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

All checks: 0 mismatches / 0 violations in stated ranges.
