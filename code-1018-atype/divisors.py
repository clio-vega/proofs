"""Classify a denominator factor by which geometric divisor it cuts out.

Three families of divisors matter in w-space, named by the ROOT that cuts them
(see Remark 2 of 2026-10-17-c2, commit 7002d00 -- the old long/short names were
inverted for type C):

  A-type   {w_i = w_j},  i != j     root  e_i - e_j     <- THIS SESSION'S TARGET
  pair     {w_i w_j = 1}, i < j     root  e_i + e_j
  boundary {w_m^h = 1}              root  h . e_m       (h=2 in type C: w_m = +-1)

Anything carrying q is a q-DEFORMED divisor: at generic q these are disjoint from
all three geometric families and are NOT poles of the geometric object.  They are
reported in their own bucket, never merged into one of the three -- that merge is
how a census miscounts: {w_i = q w_j} and {w_i = w_j} differ as SETS and agree
only at q=1.
"""
import sympy as sp
from itertools import combinations

def classify(fac, w, q):
    """Return (kind, data).  kind in {'A','pair','bdy','qdef','other','unit'}.

    Decided by VANISHING, not by pattern-matching the printed form: substitute the
    defining relation of each candidate divisor and ask whether fac dies.  A
    pattern match on strings would miss (w1 - w2)*(-1) or 1 - w2/w1.
    """
    fac = sp.together(fac)
    N = len(w)
    if sp.Poly(sp.expand(sp.numer(fac)), *w).total_degree() == 0 and q not in fac.free_symbols:
        return ('unit', None)
    # q-deformed: contains q
    hasq = q in fac.free_symbols
    for i, j in combinations(range(N), 2):
        if sp.simplify(fac.subs(w[j], w[i])) == 0:
            return ('qdef' if hasq else 'A', (i+1, j+1))
    for i, j in combinations(range(N), 2):
        if sp.simplify(fac.subs(w[j], 1/w[i])) == 0:
            return ('qdef' if hasq else 'pair', (i+1, j+1))
    for i in range(N):
        for eps in (1, -1):
            if sp.simplify(fac.subs(w[i], eps)) == 0:
                return ('qdef' if hasq else 'bdy', (i+1, eps))
    return ('qdef' if hasq else 'other', sp.sstr(fac))


def den_classify(f, w, q):
    """Reduced denominator of f, factored and classified.  Monomials dropped
    (they are units in every local ring we care about: w_i is invertible there).
    Returns dict kind -> list of (data, factor, multiplicity)."""
    num, den = sp.fraction(sp.cancel(sp.together(f)))
    out = {}
    for fac, mult in sp.factor_list(sp.expand(den))[1]:
        if len(sp.Poly(fac, *w).monoms()) == 1:
            continue                       # monomial in w: a unit
        kind, data = classify(fac, w, q)
        if kind == 'unit':
            continue
        out.setdefault(kind, []).append((data, sp.factor(fac), mult))
    return out
