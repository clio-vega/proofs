# Q140 — Khanna–Loehr local identity, 2026-09-10 cycle 2

Exact `sympy` over Q(t). Run in this directory; `rim.py` is the instrument.

| script | what it establishes |
|---|---|
| `rim.py` | rim-hook removal, twice (abacus + brute force); `crosscheck(7)` = 240/240 |
| `scan.py` | certificate check at n=4 mu=(4); anchor scan n=2..7 (reproduces the Q129 paper's table) |
| `facts.py` | \|C(mu)\|=n (n<=9); the t=-1 solution (-1)^ht/n (n<=8); generic scan n=8,9 |
| `obstruct.py` | first-order obstruction Omega(mu) at t=-1; matches consistency for all mu, n<=8 |
| `offdiag.py` | rank(M minus row mu) = n iff inconsistent |
| `murow.py` | Lemma A rows for mu=(n) (n<=11); the certificate family, c*b = t(t+1) |
| `hooks.py`, `hookmu.py` | the refuted hook-row and lambda_3<=1 routes |
| `reflect.py` | the two reflection relation families vs the empirical two-term rows |
| `graph.py` | the graph criterion; G(mu) connected and criterion == consistency, n<=9 |

Note `scan.py` and `reflect.py` print on import (they are imported by later scripts).
