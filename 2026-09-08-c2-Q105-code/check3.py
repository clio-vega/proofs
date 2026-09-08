"""CHECK 3: the vertex side.  h_N[(1+t)X] and its order of vanishing at t=-1.

Two independent routes to h_N[(1+t)X]:
  (P) power-sum route: h_N = sum_{mu|-N} p_mu / z_mu, then p_k -> (1+t^k) p_k.
  (H) creation route:  sigma[Xz1]sigma[Xz2]|_{z1=t z2} => sum_k t^{N-k} h_k h_{N-k}.
They are computed with disjoint code and compared in the Schur basis.
"""
import sympy as sp
from itertools import product
from functools import lru_cache
t = sp.symbols('t')

# ---------- symmetric functions in the power-sum basis, as dicts mu(tuple sorted desc)->coeff
def pmul(a, b):
    out = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            mu = tuple(sorted(m1 + m2, reverse=True))
            out[mu] = sp.expand(out.get(mu, 0) + c1 * c2)
    return {k: v for k, v in out.items() if v != 0}

def zmu(mu):
    z = sp.Integer(1); i = 1
    from collections import Counter
    for k, m in Counter(mu).items():
        z *= k**m * sp.factorial(m)
    return z

def parts(n):
    if n == 0:
        yield (); return
    def rec(rem, mx):
        if rem == 0:
            yield (); return
        for k in range(min(rem, mx), 0, -1):
            for tail in rec(rem-k, k):
                yield (k,)+tail
    yield from rec(n, n)

def h_p(n):
    """h_n in power sums."""
    return {mu: sp.Rational(1, 1)/zmu(mu) for mu in parts(n)}

# ---------- route P
def hP(N):
    return {mu: sp.expand(c * sp.prod([(1 + t**k) for k in mu])) for mu, c in h_p(N).items()}

# ---------- route H
def hH(N):
    out = {}
    for k in range(N+1):
        term = pmul(h_p(k), h_p(N-k))
        for mu, c in term.items():
            out[mu] = sp.expand(out.get(mu, 0) + t**(N-k) * c)
    return {k: v for k, v in out.items() if v != 0}

print("route P vs route H (power-sum coefficients):")
for N in range(0, 9):
    a, b = hP(N), hH(N)
    keys = set(a) | set(b)
    assert all(sp.expand(a.get(k,0)-b.get(k,0)) == 0 for k in keys), N
print("  agree for N=0..8")

def ordvan(poly):
    p = sp.Poly(sp.expand(poly), t)
    if p.is_zero: return None
    o = 0
    while p.eval(-1) == 0:
        o += 1; p = p.diff(t)
        if p.is_zero: return None
    return o

print("\nord_{t=-1} h_N[(1+t)X]  (min over power-sum coefficients):")
for N in range(0, 11):
    a = hP(N)
    os_ = [ordvan(c) for c in a.values()]
    print(f"  N={N:2d}  ord={min(os_)}   N mod 2 = {N%2}   (per-mu orders: {sorted(set(os_))})")
