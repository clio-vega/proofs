"""SEPARATOR for the dictionary claim.  KMS/Uglov's B_{-m} generate a HEISENBERG
algebra: [B_m,B_m']=0 for m+m'!=0.  Do my R_e(t) commute?"""
import sys, sympy as sp
sys.path.insert(0, '/home/clio/projects/probes/2026-09-06-Q84')
import engine
from bead import trim, t

def R(state, e):   return engine.op_R(state, e)
def comm(lam, e, f):
    a = R(R({trim(lam): 1}, f), e)
    b = R(R({trim(lam): 1}, e), f)
    return engine.sub(a, b)

print("[R_e, R_f] on s_lam:")
worst = None
for (e, f) in [(1,2),(1,3),(2,3),(2,4),(3,4)]:
    nz = 0; tot = 0; ex = None
    for n in range(0, 6):
        for lam in (engine.parts_of(n) if n else ((),)):
            c = comm(lam, e, f); tot += 1
            if c:
                nz += 1
                if ex is None: ex = (lam, {k: sp.factor(v) for k, v in c.items()})
    print(f"  e={e},f={f}: nonzero on {nz}/{tot} partitions  ", end="")
    if ex: print(f" e.g. lam={ex[0]}: {ex[1]}")
    else:  print(" -- COMMUTE")
