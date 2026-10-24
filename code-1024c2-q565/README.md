# code-1024c2-q565 — verification for Q576 (2026-10-24, cycle 2)

Artifact: `proofs/2026-10-24-c2-Q576-transversal-cancellation.tex`.
All arithmetic exact (`fractions.Fraction`) or symbolic (SymPy). No floating point.

`eng.py`   — this session's engine: `psi_pair` = alternant expansion with `B(T)` in the
             pair form. Cross-checked against `psi_headsum` (the head-injection sum of
             `code-1023-q559/psi.py`) — a different mechanism — in `t0_crosscheck.py`.
`tr.py`    — stable trajectory handle (absolute log_dir + pinned session tag, so every
             process appends to one file instead of restarting evt-ids).

Scripts that produced numbers quoted in the paper:

| script | claim | reported |
|---|---|---|
| `t0_crosscheck.py`  | engine vs. head sum, signs allowed, no hygiene | 350/350 |
| `t1_refute.py`      | Thm 19: conj:transversal's converse is false | 4 exact witnesses |
| `t2_grid.py`        | identical vanishing iff r >= f+2 (delta=0) | 89/89 cells |
| `t3_classes.py`     | every (Z/2)^K-orbit sum vanishes on its own | 6 cells, all orbits |
| `t4_bridge.py`      | Lemma 8 (orbit factorisation), Lemma 10 | 63/63 orbits; OK |
| `t6_scope.py`       | theta = (c,...,c+r-1) | 126/126 |
| `t7_parity.py`      | Lemma 5 (parity split); same-parity m | 216/216; 170/170 |
| `t9_corrected.py`   | Prop 23 (f=0 sign); parity-pattern scan | 140/140; 16 patterns |
| `t10_proofinputs.py`| Prop 13 + Step 1 degree bound + Step 2 divisibility | 0 violations |
| `t11_sharp.py`      | sharpened threshold r >= f+2-delta, 2 controls | 160/160, 368/368 |

Not cited in the paper, kept as session record:

- `t5_identity.py` — first attempt at the core identity symbolically via `sympy.simplify`.
  **Abandoned: too slow to finish** (never produced output; killed). Superseded by
  `t10_proofinputs.py`, which puts the sum over the common denominator in closed form and
  evaluates the numerator at exact rational points instead of simplifying.
- `t8_diag.py` — the run that refuted my first (sign-incomplete) version of Prop 23.
  Its "(G) 26/60" line is the control firing that found the missing factor
  `(-1)^{sum_j m_j}`; `t9_corrected.py` is the repair.
