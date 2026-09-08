"""Q105 vertex-side engine.  Power-sum basis, written this session.
Code-disjoint from reviews/2026-09-08-selfreview-code/hbasis.py (h-basis) and from
scratch/q99/engine.py (which it deliberately does not import).

Htilde_{a,b}(z) f[X] = sigma[Xz] * f[X - (a-b)/z],  alphabet p_n -> p_n - (a^n-b^n) z^{-n}.
H_u = Htilde_{1,u}.
"""
import sympy as sp
from collections import Counter

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

def zmu(mu):
    z = sp.Integer(1)
    for k, m in Counter(mu).items():
        z *= k**m * sp.factorial(m)
    return z

def pmul(A, B):
    out = {}
    for m1, c1 in A.items():
        for m2, c2 in B.items():
            mu = tuple(sorted(m1+m2, reverse=True))
            out[mu] = sp.expand(out.get(mu, 0) + c1*c2)
    return {k: v for k, v in out.items() if sp.expand(v) != 0}

def h(n):
    """h_n in the p-basis; h_n = 0 for n<0, h_0 = 1."""
    if n < 0: return {}
    return {mu: sp.Rational(1)/zmu(mu) for mu in parts(n)}

def shift(f, cfun):
    """f[X - c/z] as dict j -> (element), where the coefficient of z^{-j}."""
    out = {}
    for mu, coef in f.items():
        # expand prod_i (p_{mu_i} - c_{mu_i} z^{-mu_i}) over subsets
        terms = [{(): sp.Integer(1)}]
        cur = {(0, ()): sp.Integer(1)}   # (degree of z^-1, remaining partition) -> coeff
        acc = {(0, ()): sp.Integer(1)}
        for k in mu:
            nxt = {}
            for (d, nu), c in acc.items():
                key1 = (d, tuple(sorted(nu+(k,), reverse=True)))
                nxt[key1] = sp.expand(nxt.get(key1, 0) + c)
                key2 = (d+k, nu)
                nxt[key2] = sp.expand(nxt.get(key2, 0) - c*cfun(k))
            acc = nxt
        for (d, nu), c in acc.items():
            out.setdefault(d, {})
            out[d][nu] = sp.expand(out[d].get(nu, 0) + coef*c)
    return {d: {k: v for k, v in e.items() if sp.expand(v) != 0} for d, e in out.items()}

def Hmode(f, m, cfun):
    """mode m of sigma[Xz] f[X - c/z]:  sum_j h_{m+j} * c_j(f)."""
    res = {}
    for j, elt in shift(f, cfun).items():
        prod = pmul(h(m+j), elt)
        for k, v in prod.items():
            res[k] = sp.expand(res.get(k, 0) + v)
    return {k: v for k, v in res.items() if sp.expand(v) != 0}

def sub(A, B):
    out = dict(A)
    for k, v in B.items():
        out[k] = sp.expand(out.get(k, 0) - v)
    return {k: v for k, v in out.items() if sp.expand(v) != 0}

def smul(c, A):
    return {k: sp.expand(c*v) for k, v in A.items() if sp.expand(c*v) != 0}

def h_alphabet(N, afun):
    """h_N[A X] where p_n -> afun(n) * p_n."""
    if N < 0: return {}
    return {mu: sp.expand(sp.Rational(1)/zmu(mu) * sp.prod([afun(k) for k in mu]))
            for mu in parts(N)}
