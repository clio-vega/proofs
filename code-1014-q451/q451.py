"""Q451: are the individual k-terms of eq:CJT individually meaningful?

eq:CJT  (HKKO 2301.13117, Thm 3.3 / thm:ssyt2, line 877 of the .tex, PUBLISHED):
    s_{la[h,w]'}(x) = sum_{k in Z^h, sum k_i = 0}  det_{1<=i,j<=h}( e_{la_i-i+j+(h+w)k_i}(x) )

prop:cssyt=path (line ~860): K^alpha_{la[h,w]'} = # lattice paths 0 -> la in
    R = {x_1>=...>=x_h>=x_1-w}, i-th step a 0/1 vector with alpha_i ones.
So the DIRECT mechanism is   s = sum_alpha K^alpha x^alpha  -- the same region
machinery as the Q430 cylindric tableaux. Two mechanisms, no shared code path.
"""
import sys, itertools, functools
from fractions import Fraction
import sympy as sp

# ---------- engine 1: elementary symmetric / determinant (eq:CJT) ----------
def elem_polys(xs):
    n = len(xs)
    e = [sp.Integer(0)]*(n+1); e[0] = sp.Integer(1)
    for x in xs:
        for m in range(n, 0, -1):
            e[m] = sp.expand(e[m] + x*e[m-1])
    return e

class SymEngine:
    def __init__(self, n):
        self.n = n
        self.x = sp.symbols(f'x1:{n+1}')
        self.e = elem_polys(list(self.x))
    def E(self, m):
        return self.e[m] if 0 <= m <= self.n else sp.Integer(0)

def cjt_terms(eng, h, w, lam, kmax=2):
    """Return dict k-tuple -> expanded determinant, for nonzero terms."""
    out = {}
    for ks in itertools.product(range(-kmax, kmax+1), repeat=h):
        if sum(ks) != 0: continue
        M = sp.Matrix(h, h, lambda a, b: eng.E(lam[a]-(a+1)+(b+1)+(h+w)*ks[a]))
        if all(M[a,b] == 0 for a in range(h) for b in range(h)): continue
        d = sp.expand(M.det())
        if d != 0: out[ks] = d
    return out

# ---------- engine 2: direct lattice-path enumeration (prop:cssyt=path) ----
def in_region(x, w):
    return all(x[i] >= x[i+1] for i in range(len(x)-1)) and x[-1] >= x[0]-w

def direct_cyl_schur(eng, h, w, lam, n):
    """sum over paths 0->lam in R of n 0/1-steps, of x^alpha, alpha_i = #ones in step i."""
    x = eng.x
    cur = {tuple([0]*h): sp.Integer(1)}
    steps = [s for r in range(h+1) for s in itertools.combinations(range(h), r)]
    for a in range(n):
        nxt = {}
        for pos, wt in cur.items():
            for s in steps:
                np_ = list(pos)
                for i in s: np_[i] += 1
                np_ = tuple(np_)
                if not in_region(np_, w): continue
                # prune: cannot exceed lam coordinatewise, and need enough steps left
                if any(np_[i] > lam[i] for i in range(h)): continue
                nxt[np_] = nxt.get(np_, 0) + wt * x[a]**len(s)
        cur = nxt
    return sp.expand(cur.get(tuple(lam), sp.Integer(0)))

# ---------- Schur expansion ----------
def partitions_of(d, maxlen):
    def rec(rem, mx, cur):
        if rem == 0:
            yield tuple(cur); return
        for p in range(min(rem, mx), 0, -1):
            if len(cur) < maxlen:
                yield from rec(rem-p, p, cur+[p])
    return list(rec(d, d, []))

def conj(mu):
    if not mu: return ()
    return tuple(sum(1 for p in mu if p >= j) for j in range(1, mu[0]+1))

def schur_poly(eng, mu):
    """s_mu(x_1..x_n) via dual Jacobi-Trudi: det(e_{mu'_i - i + j}), size l(mu')."""
    mup = conj(mu); m = len(mup)
    if m == 0: return sp.Integer(1)
    M = sp.Matrix(m, m, lambda a, b: eng.E(mup[a]-(a+1)+(b+1)))
    return sp.expand(M.det())

def schur_expand(eng, P, d):
    """Expand homogeneous symmetric P (deg d) in the Schur basis {s_mu : l(mu)<=n}."""
    n = eng.n; x = eng.x
    if P == 0: return {}
    parts = partitions_of(d, n)
    if not parts: return {}
    svec = {mu: sp.Poly(schur_poly(eng, mu), *x) for mu in parts}
    Pp = sp.Poly(P, *x)
    mons = sorted(set(Pp.monoms()) | {m for s in svec.values() for m in s.monoms()})
    A = sp.Matrix([[svec[mu].coeff_monomial(sp.prod([x[i]**m[i] for i in range(n)]))
                    for mu in parts] for m in mons])
    b = sp.Matrix([[Pp.coeff_monomial(sp.prod([x[i]**m[i] for i in range(n)]))] for m in mons])
    sol, params = A.gauss_jordan_solve(A, b) if False else (None, None)
    sol = A.solve_least_squares(b) if A.rows != A.cols else A.LUsolve(b)
    # verify exactly
    assert sp.simplify(A*sol - b) == sp.zeros(A.rows, 1), "P is NOT in the Schur span"
    return {mu: sol[i] for i, mu in enumerate(parts) if sol[i] != 0}
