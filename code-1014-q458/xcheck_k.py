"""Engine M1 (sympy Warnaar determinant) vs census engine, general k."""
import sys, time
sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
sys.path.insert(0,'/home/clio/projects/proofs/code-1014-q458')
import sympy as sp
from core import Engine
from warn import G_times_prodx
from censusk import c_1n
t=sp.Symbol('t')

def M1(n,k,ell):
    E=Engine(n); P=sp.expand(G_times_prodx(E,k,ell))
    for x in E.x: P=sp.expand(P).coeff(x,1)
    return sp.expand(P)

ok=fail=0
cases=[(k,ell,n) for k in (1,2,3) for ell in (0,1) for n in range(1,7)]
for k,ell,n in cases:
    if k==3 and n>4: continue
    t0=time.time()
    try: a=M1(n,k,ell)
    except Exception as ex:
        print(f"k={k} ell={ell} n={n}: M1 FAILED {type(ex).__name__}", flush=True); continue
    d=c_1n(n,k,ell)
    b=sp.expand(sum(v*t**e for e,v in d.items())) if d else sp.Integer(0)
    same = sp.expand(a-b)==0
    print(f"k={k} ell={ell} n={n}: M1={a}  census={b}  MATCH={same}  ({time.time()-t0:.1f}s)", flush=True)
    ok+= same; fail+= (not same)
print(f"XCHECK ok={ok} fail={fail}")
