"""Engine C: symmetric functions as polynomials in p_1,p_2,... via Jacobi-Trudi.
Completely disjoint from the abacus / border-strip code."""
import sympy as sp
from functools import lru_cache

P = {k: sp.Symbol('p%d' % k) for k in range(1, 25)}

@lru_cache(maxsize=None)
def h(n):
    """complete homogeneous h_n in terms of power sums: n*h_n = sum_{k=1..n} p_k h_{n-k}."""
    if n < 0: return sp.Integer(0)
    if n == 0: return sp.Integer(1)
    return sp.expand(sum(P[k]*h(n-k) for k in range(1, n+1)) / n)

@lru_cache(maxsize=None)
def schur(lam):
    lam = tuple(lam)
    n = len(lam)
    if n == 0: return sp.Integer(1)
    Mx = sp.Matrix(n, n, lambda i, j: h(lam[i] - (i+1) + (j+1)))
    return sp.expand(Mx.det())

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield ()
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k):
            yield (k,) + rest

def _pmonos(deg):
    return [tuple(mu) for mu in partitions(deg)]

def _vec(expr, deg):
    """coefficient vector of a homogeneous degree-deg poly in the p's, over the
    p-monomial basis.  Coefficients may involve OTHER symbols (e.g. t); those must
    survive, so we make the p's the Poly generators and keep the full coefficient."""
    monos = _pmonos(deg)
    gens = [P[k] for k in range(1, deg+1)]
    pol = sp.Poly(sp.expand(expr), *gens)
    out = {m: sp.Integer(0) for m in monos}
    for pd, c in zip(pol.monoms(), pol.coeffs()):
        mu = []
        for k in range(1, deg+1):
            mu += [k]*pd[k-1]
        mu = tuple(sorted(mu, reverse=True))
        out[mu] = sp.expand(out.get(mu, 0) + c)
    return [out[m] for m in monos]

@lru_cache(maxsize=None)
def _schur_matrix(deg):
    monos = _pmonos(deg)
    parts = list(partitions(deg))
    M = sp.Matrix([[_vec(schur(mu), deg)[i] for mu in parts] for i in range(len(monos))])
    return M, parts

def to_schur(expr, deg):
    if deg == 0:
        return {(): sp.nsimplify(sp.expand(expr))}
    M, parts = _schur_matrix(deg)
    b = sp.Matrix(_vec(expr, deg))
    sol = M.solve(b)
    return {mu: sp.nsimplify(sol[i]) for i, mu in enumerate(parts) if sol[i] != 0}
