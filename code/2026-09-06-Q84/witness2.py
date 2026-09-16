"""witness2.py -- the SECOND witness family, mu = (e,1^d), which sees t=0.
[[vary-the-source-not-only-the-target]] / [[a-function-can-have-a-ladder-of-degenerate-slices]]
"""
from math import comb
from fast import op_p, ribbons
import sympy as sp

def op_R_int(state, e, tv):
    out = {}
    for lam, w in state.items():
        for (mu, h) in ribbons(lam, e):
            out[mu] = out.get(mu, 0) + w * tv**h
    return {k: v for k, v in out.items() if v}

def ad_d(e, d, tv, mu):
    out = 0
    for r in range(d + 1):
        c = (-1)**r * comb(d, r)
        st = {(): 1}
        for _ in range(d - r): st = op_p(st, 1)
        st = op_R_int(st, e, tv)
        for _ in range(r): st = op_p(st, 1)
        out += c * st.get(mu, 0)
    return out

T = sp.Symbol('t')
print(' e  d   <s_(e,1^d), ad^d R_e(t).1>   (interpolated in t)')
for e in (2, 3, 4, 5):
    for d in (1, 2, 3, 4):
        if d + e > 11: continue
        mu = (e,) + (1,) * d
        pts = [(v, ad_d(e, d, v, mu)) for v in range(-4, 4)]
        poly = sp.expand(sp.interpolate(pts, T))
        # confirm the interpolation on two extra points
        chk = all(sp.Integer(ad_d(e, d, v, mu)) == poly.subs(T, v) for v in (5, 7))
        print('%2d %2d   %-40s %s' % (e, d, sp.factor(poly), 'ok' if chk else 'EXTRAPOLATION FAILS'))
