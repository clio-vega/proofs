import sys, time; sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from ground import c_alpha_M1
t = sp.Symbol('t')
JOBS = [(1,1,n) for n in range(2,8)] + [(2,1,n) for n in range(2,7)] + [(3,0,n) for n in range(2,7)] + [(1,2,n) for n in range(2,7)]
for (k, ell, n) in JOBS:
    t0=time.time()
    try:
        c = c_alpha_M1(k, ell, tuple([1]*n))
    except Exception as ex:
        print(f"k={k} ell={ell} n={n}: ERROR {type(ex).__name__} {ex}", flush=True); continue
    print(f"k={k} ell={ell} h={2*k} w={2*ell+2} n={n}: c=({sp.factor(c)})  c(1)={sp.expand(c).subs(t,1)}  [{time.time()-t0:.1f}s]", flush=True)
