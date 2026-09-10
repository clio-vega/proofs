# Q129 verification scripts (2026-09-10)

Exact `sympy` computations over Q(t) backing
`../2026-09-10-Q129-reciprocity-does-not-deform.tex`.

Run from this directory (`python3 <script>.py`).

| script | what it checks | paper |
|---|---|---|
| `ribbon.py` | core: partitions, abacus e-ribbon add/remove, e-core/e-weight | Conv. 2.1 |
| `check_indep.py` | **instrument check.** abacus heights vs direct Young-diagram border-strip enumeration (120/120); `R_e(-1)=mult(p_e)` numerically via bialternant Schur (76/76) — uses neither abacus nor Murnaghan–Nakayama | §6 |
| `gate1.py` | Z rectangular; w_e and core_e commute with conjugation | Thm 3.1 |
| `repairA.py` | Zc(I-zR)=I; conjugation reciprocity holds; inversion form fails, e=2,3,4 | Thms 4.2–4.4 |
| `repairB.py` | sorting condition: t=-1 is the unique common root | Thm 5.5 |
| `repairB2.py` | KL orthogonality defect; every entry divisible by (1+t) | §6 |
| `freeweights.py` | Thm 5.3 reciprocity; free-diagonal test via `solve`, n=3,4,5 | Thm 5.3, §6 |
| `rankfast.py` | same by rank test + t=-1 calibration against x=1/Z_beta, n=3,4 | §6 |
| `local.py` | KL local identity, free weights, all mu, n=2..6 | Thm 6.1 |
| `cert.py` | explicit 5x4 certificate at n=4, mu=(4) | Thm 6.1 proof |
| `scan0.py` | anchor scan t=-1,0,1,generic over all mu, n=4..7 | §6 table |

Note: `repairB.py` prints its own output when imported by `repairB2.py`/`freeweights.py`/`rankfast.py`; ignore the duplicated header lines.
