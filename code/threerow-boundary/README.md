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

### 2026-06-16 — c=3 boundary CLOSED (proof: `../2026-06-16-c3-boundary-complete.md`)
- `explore_N.py` / `explore_N1.py` — the exact factorisations N₂=2[(P+1)(b+3)+b],
  N₁ around (P+1),(P+2)
- `numeric_check.py` — shows the 06-15 standalone bounds v₂(N₂)≥v₂(P+1)−1 etc. are **FALSE**
  (42 violations each) yet Δ(b+i)>−θ holds
- `direct_forms.py` — **Lemma D**: the three direct Δ(b+i) carry-formulas, vs MN (m≤40, 0 mismatch)
- `mn_crosscheck.py` — M_{b+1},M_{b+2},M_{b+3} closed forms (incl. N₁,N₂) vs MN
- `aeven_check.py` — a-even product forms (Δ(b+1),Δ(b+2)), v₂(N₂)=1, W₁≥3
- `aodd_check2.py` — a-odd product forms + the four key Lemma-F inequalities
- `lemmaF2.py` — **two-factor Lemma F2** (Q≥6) + single-factor F1 (Q≥4) + Q=5 sharpness
- `final_endtoend.py` — end-to-end Δ(b+i)>−θ, both parities, **m≤60, 3828 indices, 0 violation**

All checks: 0 mismatches / 0 violations in stated ranges.
See `../FINDINGS-2026-06-15-jobA-boundary.md` for the structural write-up.
