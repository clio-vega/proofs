import sys, time; sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from ground import c_alpha_M1
t = sp.Symbol('t')
for (k,ell,n) in [(1,1,8),(1,1,10)]:
    t0=time.time()
    c = c_alpha_M1(k, ell, tuple([1]*n))
    print(f"k={k} ell={ell} h={2*k} w={2*ell+2} n={n}: c=({sp.factor(c)})  c(1)={sp.expand(c).subs(t,1)}  [{time.time()-t0:.1f}s]", flush=True)
