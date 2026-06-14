# threerow-c3 — `(a,b,3)` even-|J*| (d=4)

Companion code to `proofs/2026-06-14-threerow-c3-Jstar-even.md`. The **first family with `|J*|=4`**
(both generators 2 and 4 active), minimal witness `(9,6,3)`.

- `mn.py` — Murnaghan–Nakayama `M_j` (shared ground truth).
- `c3factor.py` — derives the closed form (Lemma 1): base `C(N,b−j)`, prefactor `(a−b+1)`, sextic
  numerator `Q_3`, denominator `6(a+4−j)(b+1−j)(b+2−j)(b+3−j)`.
- `c3verify.py` — Lemma 1 vs MN (`m≤24`, 2170 cases, 0 mismatch) + binomial-basis decomposition.
- `c3prop2.py` — Proposition 2 Kummer `Δ(j)` vs direct (`m≤21`, 0 mismatch).
- `c3decomp.py` — the structural identity `Q_3 = (a−1)(b−2)H − 720 C(j,6)`.
- `c3census.py` — box by `(a,b) mod 8` (`m≤29`): `|J*|∈{1,2,4}`, offset by parity of `a`.
- `c3detail.py` — `Δ(j)` by `(a,b) mod 4`; forced descent and relative ties.
- `c3comp.py`,`c3reduc.py`,`c3lemtest.py` — Compensation Lemmas A/B, the (L3″) reductions, the
  forced descent, and `\tildeΔ≥0`, all verified `m≤79`.
- `c3handcheck.py` — confirms the §4 hand-derived formulas (`Δ(1)=1`, `Δ(2)=2(v₂H(2)−2)`, tie
  congruences), 0 failures.
- `c3full.py` — full Theorem vs MN (`m≤39`, 630 shapes, 0 mismatch).

**Status:** closed form + structural identity + Prop-2 + offset theorem proved; even-|J*| reduced to
two explicit 2-adic inequalities (Gaps 1–2), verified `m≤79`. See proof §7 for the named gaps.

## 2026-06-14 code-session cross-check (independent engine)

- `job_jstar_engine.py` — independent `M_j` via vertical-2-strip removal chains (≠ the `D`-sum derivation).
- `2026-06-14-jobB-census-sawin.py` — **independently verifies Lemma 1** (1627 checks, 0 fail); box
  congruence rule `S` by `(a,b) mod 4`; identity `G_λ(i)=Σ_j C(m,j)M_j(i−1)^j` with `v_π=val(j)`; the
  Sawin adjacent-pair test (**refuted**; the fpf involution is the top-generator toggle `j↔j+4`).
- `2026-06-14-jobB-fit.py`, `2026-06-14-jobB-c3-explore.py` — closed-form fitting scratch (how the
  denominator `6(a+4−j)(b+1−j)(b+2−j)(b+3−j)` was re-found).
- `FINDINGS-2026-06-14-census-sawin.md` — full write-up of the three items.
