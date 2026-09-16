# Q150 — the nonvanishing lemma for ribbon shape weights is false

Code for `proofs/2026-09-16-c1-Q150-nonvanishing.tex`.

| file | what it does |
|---|---|
| `direct.py` | independent numeric commutator checker, written from `conv:def` of the Q149 paper; shares no algebra with `solver.py`. Calibrated against it on `(2,3)` before use (V0). |
| `reduce.py` | the structural reduction: `(1B)+(2B)` <=> `(C)` + `(2B*)`, the two-bead instance lists, `tensor()` for `tau^(x)k`, admissible-support enumeration, exhaustive `F_p` solvers |
| `verify.py` | the verification suite V1–V5 quoted in §8 of the paper. Run `python3 verify.py`. |
| `engine.py`, `solver.py` | copies of the Q149 engine, used only as the independent second code path |
| `run_2026-09-16.log` | output of the smoke test committed alongside |

**Result.** With `d = gcd(e,f)`: the one-bead sector says exactly `sigmabar = sigma (x) rho =
rho (x) sigma`, so `sigma` and `rho` commute under concatenation; commuting weights are powers
of a common weight `tau` of rank `d`; and the weights vanish nowhere iff `tau` does — which the
two-bead relation forces when `d <= 2` and does not when `d >= 3`.

**Two calibration bugs found and fixed before use**, both in `direct.py`: it dropped legal moves
originating below the frozen tail, and it compared residuals over `Z` against solutions over
`F_3` (four false failures, residuals ±3, ±6). It now tracks the whole window and takes a
modulus.

**One scope error found after the fact.** An earlier scan restricted to *singleton* supports
reported "no zero-bearing solution" at `(6,9),(9,12),(8,12),(10,15)` — all with `gcd >= 3` — and
suggested the criterion `e | f`. The composite `tau^(x)2` at `(6,9)` refutes that. Counts of
singletons do not answer questions about weights.
