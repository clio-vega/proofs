import sympy as sp, itertools
from core import Engine
t = sp.Symbol('t')

def G(E, h, ell, tval=t, perturb=None):
    """Warnaar Eq_generalise, det size h, w=2*ell+2, M=h+ell+1, N=2M=2h+w.
    perturb: optional callable (i,j,y)->extra integer added to the t-exponent (control arm)."""
    M = h + ell + 1; N = 2*M; n = E.n
    ymax = (2*n + 2*h + 4)//N + 1
    tot = sp.Integer(0)
    for ys in itertools.product(range(-ymax, ymax+1), repeat=h):
        def ent(i, j):
            i += 1; j += 1; y = ys[i-1]
            ex = M*y*y - j*y
            if perturb is not None: ex += perturb(i, j, y)
            body = E.Edot(n-i+j-N*y) - E.Edot(n+i+j-N*y)
            return tval**ex * body
        Mx = sp.Matrix(h, h, ent)
        if all(all(Mx[a,b]==0 for b in range(h)) for a in range(h)): continue
        tot += Mx.det()
    return sp.expand(tot)

def G_times_prodx(E, h, ell, tval=t, perturb=None):
    return sp.expand(sp.together(G(E,h,ell,tval,perturb) * E.prodx**h))
