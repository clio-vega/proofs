# Q484 computational verification — Motegi type B/C weight functions

Companion to `../2026-10-15-Q484-motegi-puzzle-Gm.tex`.  Pure Python + SymPy
(SageMath is not installed in this container).  Exact arithmetic over
`C(q)(w_1,...,w_N)` throughout -- no floating point.

| file | what it is |
|---|---|
| `q484_engine.py` | the engine: `F_{N,n,r}` from Motegi's symmetrisation formula by brute-force `S_n` sum; `E^C` normaliser; `Atilde_{M,I}` via the stability relation at `n=s`; the independent positive-root product for the diagonal; Motegi's order `>=` |
| `q484_plants.py`  | engine plants P1 `(u_j u_k - q) -> (u_j u_k - 1)` and P2 `(u_j^2 - 1) -> u_j^2` |
| `q484_branch.py`  | Control 3: cross-sector branching (Motegi L2886) as a consistency control with the `j>r` blocks in its domain |
| `q484_biorth.py`  | Control 4: all-cardinality biorthogonality (Motegi L4364).  Clean, but INSENSITIVE to both plants -- reported as such |
| `q484_struct.py`  | structure constants by forward substitution (Theorem B) + tests of Theorems A, C, D and Laurent polynomiality |
| `q484_N3.py`, `q484_N3.log` | the `N=3` run: Theorem A 512/512, Cor C 35/35, Laurent 176/195 |

## Reproducing the headline numbers

```
python3 -c "import sys;sys.path.insert(0,'.');from q484_engine import *;import sympy as sp
for N in (1,2,3):
    w=wvars(N)
    print(N, all(sp.simplify(sp.cancel(Atilde(N,w,I,I,'C')-diag_rootproduct(N,w,I)))==0 for I in subsets(N)))"
python3 q484_struct.py 2     # N=2: Thm A 64/64, Cor C 10/10, Thm D 64/64
python3 q484_branch.py       # Control 3 clean 27/27; plant P1 fires; plant P2 spared
python3 q484_plants.py       # Control 1 vs the plants: P2 fires, P1 does NOT (the domain fault)
python3 q484_N3.py           # ~15 min
```

## The domain fault, in one line

Control 1 runs on `M = I`, where `s = r`, hence `n = r`, hence the index range
`r+1 <= j < k <= n` is **empty**: three of the seven product blocks of Motegi's
formula cannot appear there.  Its clean 14/14 validates **4 of 7 blocks**.
Control 3 exists to cover the other three, and it fires on the plant that
Control 1 cannot see.
