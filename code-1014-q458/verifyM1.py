"""Mechanism 1: the actual Warnaar determinant engine (proofs/code-1010-q410).
Compare the FULL polynomial c_{(1^n)}(t), not just its degree."""
import sys, time
sys.path.insert(0, '/home/clio/projects/proofs/code-1010-q410')
sys.path.insert(0, '/home/clio/projects/proofs/code-1014-q458')
import sympy as sp
from core import Engine
from warn import G_times_prodx
from census import closed_form
t = sp.Symbol('t')

def coeff_alpha_1n(n, ell):
    """c_{(1^n)}(t) via engine M1, k=1."""
    E = Engine(n)
    P = G_times_prodx(E, 1, ell)
    # extract coefficient of x_1^1 ... x_n^1
    expr = sp.expand(P)
    for x in E.x:
        expr = sp.expand(expr).coeff(x, 1)
    return sp.expand(expr)

if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    ok = fail = 0
    for ell in (0,1,2):
        for n in range(1, nmax+1):          # include ODD n -- predicted identically 0
            t0 = time.time()
            got = coeff_alpha_1n(n, ell)
            pred = closed_form(n, ell)
            same = sp.simplify(sp.expand(got - pred)) == 0
            print(f"ell={ell} n={n}  M1={got}   pred={pred}   MATCH={same}  ({time.time()-t0:.1f}s)",
                  flush=True)
            if same: ok += 1
            else: fail += 1
    print(f"RESULT ok={ok} fail={fail}")
