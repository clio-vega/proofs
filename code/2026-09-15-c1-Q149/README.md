# Q149 — free SHAPE weights for the ribbon operator

Engine for `proofs/2026-09-15-c1-Q149-shape-weights.tex`.  The weight is indexed by the
**occupancy word** `u in {0,1}^{e-1}` of the open interval `(b, b+e)` — equivalently by
the composition of `e` giving the ribbon's row lengths — so it is an opaque symbol per
ribbon SHAPE, not per height.  Nothing forms `(-1)^{ht}` or `t^{ht}`.

| file | what it does |
|---|---|
| `engine.py` | two independent ribbon enumerations (Young row-composition vs Maya occupancy word), cross-checked on 1375 ribbons; `word_to_comp` **reverses** (larger bead = higher row) |
| `solver.py` | commutator equations from the raw operator; two state enumerations (all charge-0 Maya sets in an LxL box; all window configs `Z_{<0} u T`) that agree exactly |
| `quotient.py` | the `d`-quotient weights `MNd(e,d)` and the general character weights `gamma_weight` |
| `controls.py` | realisability (C2), calibration to Q147 (C1), `e=f` (C4), exhaustive classification over `F_p` (C5) |
| `degen.py` | nonvanishing lemma at `e=3`: every proper zero-pattern of `W` forces `Wbar=0` (C6) |
| `gaugetest2.py` | cocycle/spanning-tree test: is `W^gamma / W^MN` a diagonal gauge of MN? (C7) |
| `run_controls.py` | driver for `controls.py` |

Run `python3 engine.py` first (cross-check), then any other module.

**Result.** With `d = gcd(e,f)`, the commuting locus is
`W(u) = alpha * prod_j (-gamma(j mod d))^{u_j}` for `gamma : Z/d -> F^x` with
`gamma(0)=1` and `gamma(x)gamma(-x)=1`.  One-bead sector => character + d-periodic;
two-bead sector => transpose-equivariant.  `d=1` collapses to Murnaghan-Nakayama.
