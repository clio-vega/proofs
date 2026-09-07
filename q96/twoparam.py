"""Verify the two-parameter commutator [R_e(t), R_f(s)] matrix-element formulas."""
import sympy as sp
from engine import beta, from_beta

t, s = sp.symbols('t s')

def Rmaya(bs, e, par, n):
    """R_e(par) on a Maya set (as frozenset of beta numbers with n beads)."""
    out = {}
    for b in bs:
        if b + e in bs: continue
        ht = sum(1 for x in bs if b < x < b + e)
        nb = frozenset((bs - {b}) | {b + e})
        out[nb] = out.get(nb, 0) + par**ht
    return out

def comm_maya(bs, e, f, n):
    """[R_e(t), R_f(s)] |bs>"""
    out = {}
    for m1, c1 in Rmaya(bs, f, s, n).items():
        for m2, c2 in Rmaya(m1, e, t, n).items():
            out[m2] = out.get(m2, 0) + c1*c2
    for m1, c1 in Rmaya(bs, e, t, n).items():
        for m2, c2 in Rmaya(m1, f, s, n).items():
            out[m2] = out.get(m2, 0) - c1*c2
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}

def predict(M, Mp, e, f):
    """Predicted <Mp|[R_e(t),R_f(s)]|M> from the two-parameter theorem."""
    D = M ^ Mp
    rem = sorted(M - Mp); add = sorted(Mp - M)
    m = lambda x: 1 if x in M else 0
    cnt = lambda lo, hi: sum(1 for x in M if lo < x < hi)
    if len(D) == 2:                                     # one-bead sector
        a = rem[0]
        if add[0] != a + e + f: return sp.Integer(0)
        N = cnt(a, a+e+f); A = cnt(a, a+f); B = cnt(a, a+e)
        me, mf = m(a+e), m(a+f)
        return sp.expand(s**A*t**(N-1-A)*(t - mf*(1+t)) - t**B*s**(N-1-B)*(s - me*(1+s)))
    if len(D) == 4:                                     # two-bead sector
        # sum over LEGAL ASSIGNMENTS of which operator moves which bead.
        # e != f : exactly one assignment.  e == f : both.
        tot = sp.Integer(0)
        for b in rem:
            c = (set(rem) - {b}).pop()
            if b + e not in add or c + f not in add: continue
            if len({b, c, b+e, c+f}) != 4: continue
            P = cnt(b, b+e); Q = cnt(c, c+f)
            k = (1 if c < b+e < c+f else 0) - (1 if c < b < c+f else 0)
            tot += t**(P-k)*s**Q*(1 - (t*s)**k)
        return sp.expand(tot)
    return sp.Integer(0)
