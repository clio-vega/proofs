# code-q227 — the TASEP dictionary for Theorem 4.1's exchange move

Runs on top of the unmodified `../code-q220/` harness (`affine.py`, `emptyX.py`).

| file | what it does |
|---|---|
| `beads.py` | the instrument: Lam's nilCoxeter action `u_i` in bead coordinates (= the TASEP exclusion rule), and the cyclically decreasing word |
| `lemma1.py` | exhaustive check of Lemma 2.1 (block travel), 85860 `(S,A)` pairs, N=3..8 |
| `order.py` | settles the operator order: `u_T` acts first (1145/1145 at N=6); nilCoxeter consistency of the exchange move (1326/1326) |
| `census.py` | first pass; superseded by `verdict.py` (kept: its C1 test conflated "e=m" with "every run usable", which is how the case analysis got tightened) |
| `verdict.py` | the decisive run: realisable ⟺ 321-avoiding, `e=m`, the replacement count, `(C)` |
| `kill.py` | the two refutations of `(C)`, structural and numerical |
| `invar.py` | non-degenerate witnesses; the arc-multiset invariance (Cor 4.3) |

Paper: `../2026-09-22-c2-exchange-move-is-a-retiming.tex`.
