import sys, sympy as sp
sys.path.insert(0, '/home/clio/projects/probes/2026-09-06-Q84')
import engine
from bead import R_bead, trim, maya, open_interval_load, t

agree = 0; tot = 0; bad = []
loads = {}
for e in (2, 3, 4, 5, 6):
    for n in range(0, 10):
        for lam in (engine.parts_of(n) if n else ((),)):
            lhs = {k: sp.expand(v) for k, v in engine.op_R({trim(lam): 1}, e).items()}
            rhs = {k: sp.expand(v) for k, v in R_bead(lam, e).items()}
            tot += 1
            if lhs == rhs: agree += 1
            else: bad.append((e, lam, lhs, rhs))
            N, b = open_interval_load(lam, e)
            loads[N] = loads.get(N, 0) + 1
print(f"WIDE: e in 2..6, |lam| <= 9 :  {agree}/{tot}")
for b in bad[:3]: print("  MISMATCH", b[0], b[1], b[2], b[3])
print("distribution of max open-interval load:", dict(sorted(loads.items())))
