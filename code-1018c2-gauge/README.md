# PROVE 2026-10-18 c2 -- gauge theorem and (Q)-divisors

`q484_engine.py` is the canonical copy, md5 `f53c0258` (checked before use).

Scripts, in the order they were run:

| script | what it tests | status |
|---|---|---|
| `t1_smoke.py`    | symbolic N=2, d=1 (SRC) and d=2 (QC) through the generic `AtilU` route; F-, E- and Atil-ratios **separately** | green; found the engine's `n==0` branch bug (returns 1, correct only for r=0) |
| `t3_qdivisor.py` | first (Q)-divisor probe, `sqrt`-based uniformiser | **superseded** -- the radical defeats the order test |
| `t5_qvanish.py`  | second (Q)-divisor probe, also `sqrt`-based | **superseded**; disagreed with t3 on the *same object* (`ord(det Phi)` 0 vs 2) |
| `t6_qvanish2.py` | (Q)-divisor structure, **rational** parametrisation `q -> (1-t)/m(w)` | green 6/6 + control; this is the reported run |
| `t7_laurent.py`  | is `Atil` a Laurent polynomial? | 3/3, 10/10, 35/35 at N=1,2,3 |
| `t8_fvanish.py`  | the F-vanishing lemma, with MUST-FAIL plants | **had not returned** when the paper was built; gap G1 |
| `t9_fast.py`     | fast exact-rational census, **no genericity guard** | 11 "defects", every one a degenerate seed |
| `t10_fast2.py`   | same, with a genericity guard and skipped-count reported | green except the repeated-value rows, which is the real finding |
| `t2_census.py`   | full symbolic census | **incomplete** (N=3 block done at 1113 s, still in plants); nothing in the paper rests on it |

The two superseded scripts and `t9` are kept deliberately: `t3`/`t5` are the
recorded instrument fault (two measurements of one object, opposite answers), and
`t9` is the run whose 11 red lines were all seed degeneracies.
