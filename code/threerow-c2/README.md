# threerow-c2 — even-|J*| for λ=(a,b,2) (d=4)

Companion code for `projects/proofs/2026-06-13-threerow-c2-Jstar-even.md`.
Reuses `../threerow-c1/mn.py` (Murnaghan–Nakayama M_j) and `dj.py` (S₃ trinomial determinant).

- `c2factor.py`,`c2num.py` — derive/verify Lemma 1 (closed form for M_j); Q(a,b,j).
- `c2verify2.py`,`c2bdry.py` — verify closed form b≥2 and boundary M_{b+1},M_{b+2}.
- `c2prop2.py` — verify Proposition 2 (Δ(j) Kummer formula).
- `c2census.py`,`c2full.py` — J* census; full theorem J*∈{{0},{0,2},{0,4}}.
- `c2complemma.py`,`c2Lpp.py`,`c2numlem3.py` — Compensation Lemma + reductions + Number Lemma.
- `c2bulk.py` — min T(j) pattern (=1−v₂(j)).
- `c2delta.py`,`c2zero.py` — Δ(2),Δ(4) tie conditions mod 4.
- `c2checks.py` — b=2 hand-formulas + never-both (T(4)≥0 on Δ(2)=0 classes).
- `c2resid.py` — residual boundary inequality (b≥3), m≤80.
