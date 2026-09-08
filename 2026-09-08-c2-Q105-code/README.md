# Q105 verification code

Primary engines, both written 8 September 2026 (cycle 2) and code-disjoint from each other:

- `abacus.py`  — Maya sets and bead moves only.  No Young diagrams, no symmetric functions.
- `vertex.py`  — power-sum dictionaries.  No h-basis, no abacus.

Cross-checks use `clio-vega/rick-review@314f24a:2026-09-08-selfreview-code/{ribbon.py,hbasis.py}`
(border strips on Young diagrams; h-basis dictionaries).  Neither was written this session and
neither knows anything about Q105.  Run `xcheck.py` from a directory where both are importable.

| script | what it establishes | result |
|---|---|---|
| `check1.py` | abacus engine vs Q96 Thm 4.1(1), generic symbolic (t,s) | 1232/1232, 597 two-bead entries |
| `check2.py` | on s=1/t: two-bead vanishes; one-bead zero iff kappa=0; order exactly 1 | 1232 nonzero all order 1; 1356 zeros; 0 failures |
| `check3.py` | h_N[(1+t)X] by two routes; ord = N mod 2 with the full ladder | N=0..10 |
| `check4.py` | Schur matrix elements palindromic; order histograms | 513/513 palindromic |
| `check5.py` | Q99 Thm B re-derived, at u=t and at u=-t | 252/252 each |
| `check6.py` | zero/pole decomposition with four independent symbols a,b,c,d | 64+64 |
| `check7.py` | THE VARYING TEST: diagonal sum = h_N[(1+b)X] with a FREE | 35/35 |
| `check8.py` | the closed form (-1)^N (-t)^P (1-(-t)^kappa) | 2588/2588 |
| `xcheck.py` | both sides on the two code-disjoint review engines | 764/764 and m+n=0..8 |
