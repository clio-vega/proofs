import sympy as sp, itertools
from core import Engine

def cminus(la, h, w):
    """c^-_{2h,w}(la), la a tuple of length exactly 2h (zero-padded)."""
    assert len(la) == 2*h
    if all(la[2*i] == la[2*i+1] for i in range(h)):
        return 1
    if la[0]-la[2*h-1] == w and all(la[2*i+1] == la[2*i+2] for i in range(h-1)):
        return -1
    return 0

def cyl_schur_transpose(E, la, m, w):
    """s_{la[m,w]'}(x) via HKKO eq:CJT. la = tuple length m (zero-padded)."""
    N = m + w
    n = E.n
    # k in Z^m, sum=0; entries la_i-i+j+N k_i must land in [0,n] for some j
    kmax = (n + m + 2)//N + 1
    tot = sp.Integer(0)
    for ks in itertools.product(range(-kmax, kmax+1), repeat=m):
        if sum(ks) != 0: continue
        Mx = sp.Matrix(m, m, lambda i, j: E.E(la[i]-(i+1)+(j+1)+N*ks[i]))
        if all(all(Mx[i,j]==0 for j in range(m)) for i in range(m)): continue
        tot += Mx.det()
    return sp.expand(tot)

def hkko_C_minus_RHS(E, h, w, lamax=None):
    """sum_{la in Par(2h,w)} c^-(la) s_{la[2h,w]'}(x)."""
    m = 2*h; n = E.n
    if lamax is None: lamax = 2*n
    tot = sp.Integer(0); used = []
    for la in itertools.combinations_with_replacement(range(lamax, -1, -1), m):
        la = tuple(la)  # weakly decreasing
        if la[0]-la[m-1] > w: continue
        c = cminus(la, h, w)
        if c == 0: continue
        s = cyl_schur_transpose(E, la, m, w)
        if s == 0: continue
        tot += c*s; used.append((la, c))
    return sp.expand(tot), used

def hkko_C_minus_LHS(E, h, w):
    N = 2*h + w
    Mx = sp.Matrix(h, h, lambda i, j: E.F(-(i+1)+(j+1), N) - E.F((i+1)+(j+1), N))
    return sp.expand(Mx.det())
