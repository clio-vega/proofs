# Q147 — free height weights for the ribbon operator

Engine written fresh for `2026-09-11-c2-Q147-multi-t-bracket.tex`.  It never forms
`(-t)^N` or `t^ht`: the weight enters only as an opaque symbol, so it *derives* the
locus rather than comparing against a known answer.

| file | what it does |
|---|---|
| `engine.py` | two independent ribbon enumerations (Young-diagram brute force vs Maya bead move), cross-checked; `apply_R`, `commutator` with an arbitrary weight callable `w(e,N)` |
| `run1.py`, `run2.py` | the commuting locus, **shared** weights, 11 pairs `(e,f)` |
| `run3.py` | the commuting locus, **independent** weights, 8 pairs |
| `run4.py` | Lemma 5 (one-bead matrix element) checked against the engine, 316 configurations |
| `run5.py` | case (B) of the one-bead sector alone; and the unnormalised ($w_0$ free) locus |
| `run6.py`, `run10.py` | sector split: one-bead locus vs two-bead locus, shared and independent |
| `run7.py` | Lemma 12 (two-bead matrix element) checked, 42 configurations; negative controls N1–N3 |
| `run8.py` | the anchor identification $R_e^{((-1)^\bullet)}=p_e\cdot$ via SSYT → Kostka → monomial → inverse Kostka (no ribbon/abacus/MN input), 76/76 |
| `run9.py` | the degenerate stratum ($w_0=0$), independent weights |

Run from this directory: `python3 engine.py` first (cross-check), then any `runN.py`.
