"""controls.py -- negative controls for the order harness.

[[degenerate-evidence-has-a-kernel]] / [[negative-control-must-move-a-prediction]]:
a harness that reports "never vanishes" is worthless unless it DOES vanish,
at the right depth, on operators of known finite order.  Four controls:

  C0  M_{h_2}                     ord 0   -> C^(1) = 0
  C1  M_{p_1} . p_1^perp          ord 1   -> C^(2) = 0,  C^(1) != 0
  C2  M_{p_e p_1 p_1} . p_1^perp p_1^perp
                                  ord 2   -> C^(3) = 0,  C^(2) != 0
      (degree +e and built from the same M/perp primitives as N_e: the
       structurally parallel control)
  C3  N_2 vs M_{h_2}: the MINIMAL PAIR.  Same 2-box enumeration, one predicate
      flipped (2 components vs horizontal strip).  ord 0 against ord infinity.
"""
import itertools
from math import comb
from fast import op_p, op_N, L, beta, unbeta, trim
from functools import lru_cache

@lru_cache(maxsize=None)
def coribbons(lam, a):
    """remove a connected a-ribbon: (nu, ht)."""
    B = set(beta(lam)); out = []
    for b in sorted(B):
        if b - a < 0 or b - a in B: continue
        h = sum(1 for c in B if b - a < c < b)
        out.append((unbeta((B - {b}) | {b - a}), h))
    return tuple(out)

def op_pperp(state, a):
    out = {}
    for lam, w in state.items():
        for (nu, h) in coribbons(lam, a):
            out[nu] = out.get(nu, 0) + w * (-1)**h
    return {k: v for k, v in out.items() if v}

@lru_cache(maxsize=None)
def hstrip2(lam):
    """mu/lam a horizontal 2-strip (<=1 box per column) -- gives M_{h_2}."""
    lam = trim(lam); nr = len(lam) + 2
    lamp = lam + (0,) * (nr - len(lam)); res = []
    def rec(i, prev, rem, mu):
        if rem == 0: res.append(tuple(mu) + tuple(lamp[i:])); return
        if i == nr: return
        for v in range(lamp[i], min(prev, lamp[i] + rem) + 1):
            rec(i + 1, v, rem - (v - lamp[i]), mu + [v])
    rec(0, 10**9, 2, [])
    out = []
    for mu in res:
        mu = trim(mu)
        cells = [(i + 1, j + 1) for i in range(len(mu))
                 for j in range((lamp[i] if i < len(lamp) else 0), mu[i])]
        cols = [j for (i, j) in cells]
        if len(set(cols)) == len(cols): out.append(mu)
    return tuple(out)

def op_h2(state, e=None):
    out = {}
    for lam, w in state.items():
        for mu in hstrip2(lam): out[mu] = out.get(mu, 0) + w
    return {k: v for k, v in out.items() if v}

def mk_C1():
    def X(state, e=None):
        return op_p(op_pperp(state, 1), 1)
    return X

def mk_C2(e):
    def X(state, _=None):
        s = op_pperp(op_pperp(state, 1), 1)
        s = op_p(op_p(op_p(s, e), 1), 1)
        return s
    return X

def nested(X, e, d, lam):
    out = {}
    for r in range(d + 1):
        c = (-1)**r * comb(d, r)
        st = {lam: 1}
        for _ in range(d - r): st = op_p(st, 1)
        st = X(st, e)
        for _ in range(r): st = op_p(st, 1)
        for k, v in st.items(): out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v}

def parts(n, cap=None):
    if cap is None: cap = n
    if n == 0: return [()]
    return [(p,) + r for p in range(min(n, cap), 0, -1) for r in parts(n - p, p)]

if __name__ == '__main__':
    lams = [l for n in range(0, 6) for l in parts(n)]
    tests = [('C0  M_{h_2}                  (ord 0)', op_h2, None, 3),
             ('C1  M_{p_1} p_1^perp         (ord 1)', mk_C1(), None, 4),
             ('C2  M_{p_3 p_1 p_1} (p_1^perp)^2 (ord 2)', mk_C2(3), None, 4),
             ('N   N_2                      (ord ?)', op_N, 2, 6),
             ('N   N_3                      (ord ?)', op_N, 3, 6)]
    for (name, X, e, dmax) in tests:
        print(name)
        for d in range(1, dmax + 1):
            nz = sum(1 for lam in lams if nested(X, e, d, lam))
            print('     C^(%d): %2d/%2d lam nonzero%s' % (d, nz, len(lams),
                  '   <== VANISHES' if nz == 0 else ''))
