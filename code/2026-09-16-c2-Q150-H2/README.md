# Q150 gap (H2) — prove cycle 2, 2026-09-16

Companion code for `proofs/2026-09-16-c2-Q150-H2-single-shape.tex`.

| file | what it does |
|---|---|
| `independent_maya.py` | the operator `R_e^W` written from Q149 Conv. 1.2 alone — **no code shared** with the Q149/Q150-c1 engines. Acts on Maya sets directly. |
| `twobead.py` | the relations (2B*_d), (2B**_d) and the support condition (S). |
| `supports.py` | DPLL enumeration of all admissible supports (`python3 supports.py 7`). |
| `analyse.py` | cylinder containment + singleton criterion vs the data. |
| `singleton.py` | the border criterion and the antipalindrome criterion vs direct evaluation, d ≤ 9. |
| `structure.py` | Theorem 6.1 (flip-invariance below the minimal j) and "1_Z is always a solution". |
| `lib.py` | side-effect-free helpers (border/criterion/tensor power/witness). |
| `witness.py`, `hard.py` | the **constructed** two-bead negative control: 112/112 hard failing words fire, d ≤ 10. |
| `liveness.py` | exhaustive window scan at (3,6) and (4,8) with a two-bead liveness count. |
| `values.py` | Smith normal form of the value lattice per support, d ≤ 6. |
| `newcounter.py`, `mayatest.py`, `rerun.py`, `nc7.py` | earlier scans, kept because they record two instrument failures (see below). |

## Two instrument failures worth keeping

1. **DPLL backtracking bug.** `supports.py` first reported 2 and 4 admissible sets at d=4,5
   where brute force gives 8 and 26: `propagate` mutated the assignment before detecting a
   conflict and the caller did not undo those. Brute force at d ≤ 5 is the calibration.
2. **A silent negative control.** Scanning small *partitions* for a nonzero commutator is
   degenerate: the two-bead witness needs a window of width `e+f`, which no partition of size
   ≤ 14 realises. The fix is to sample Maya sets directly, and better, to *construct* the
   witness (`witness.py`). An off-by-one in that construction (δ = n+1−p instead of n−p) then
   made it report zero for the positive cases too — which is how it was caught.
