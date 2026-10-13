import sys, time; sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from core import Engine
from warn import G_times_prodx
t = sp.Symbol('t')
def c_alpha_M1(k, ell, alpha):
    n = len(alpha); E = Engine(n); F = G_times_prodx(E, k, ell)
    return sp.expand(sp.Poly(F, *E.x).coeff_monomial(sp.prod([E.x[i]**alpha[i] for i in range(n)])))
if __name__ == '__main__':
    for (k, ell) in [(1,0),(2,0),(3,0),(1,1),(2,1)]:
        for n in range(1, 9):
            t0=time.time()
            try:
                c = c_alpha_M1(k, ell, tuple([1]*n))
            except Exception as ex:
                print(f"k={k} ell={ell} n={n}: ERROR {ex}", flush=True); continue
            print(f"k={k} ell={ell} h={2*k} w={2*ell+2} n={n}: c=({sp.factor(c)})  c(1)={c.subs(t,1)}  [{time.time()-t0:.1f}s]", flush=True)
