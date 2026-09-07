"""Engine C: pure ROUTE-COUNTING.  No Clifford algebra, no fermion signs.
Matrix elements of [R_e,R_f] obtained by enumerating the (at most two) ways a
product of an e-move and an f-move on the abacus reaches a given target.

  one-bead sector, target b -> b+e+f :   (m(b+e) - m(b+f)) * t^(N-1) * (1+t)
        N = #(M cap (b,b+e+f))
  two-bead sector, targets b->b+e, c->c+f (four sites distinct):
        t^(P+Q) * (t^(-k) - t^(k)),   P=#(M cap (b,b+e)), Q=#(M cap (c,c+f)),
        k = [c<b+e<c+f] - [c<b<c+f]
"""
import sympy as sp
t = sp.Symbol('t')
W = 24

def maya(lam):
    lam = list(lam) + [0]*(W-len(lam))
    return frozenset(lam[j-1]-j for j in range(1, W+1))

def unmaya(M):
    a = sorted(M, reverse=True)
    lam = [a[j-1]+j for j in range(1, len(a)+1)]
    while lam and lam[-1]==0: lam.pop()
    return tuple(lam)

def cnt(M, lo, hi): return sum(1 for j in range(lo+1, hi) if j in M)
def m(M, a): return 1 if a in M else 0

def routes(lam, e, f):
    M = maya(lam); out = {}
    def add(Mp, c):
        k = unmaya(Mp); out[k] = sp.expand(out.get(k, 0) + c)
    # --- one-bead sector ---
    for b in range(-W+1, W):
        if b+e+f >= W: continue
        if b not in M or (b+e+f) in M: continue
        N = cnt(M, b, b+e+f)
        coeff = (m(M, b+e) - m(M, b+f)) * t**(N-1) * (1+t)
        if coeff != 0: add(M - {b} | {b+e+f}, coeff)
    # --- two-bead sector ---
    for b in range(-W+1, W):
        for c in range(-W+1, W):
            if b+e >= W or c+f >= W: continue
            if len({b, c, b+e, c+f}) != 4: continue
            if b not in M or c not in M: continue
            if (b+e) in M or (c+f) in M: continue
            k = (1 if c < b+e < c+f else 0) - (1 if c < b < c+f else 0)
            if k == 0: continue
            P = cnt(M, b, b+e); Q = cnt(M, c, c+f)
            add(M - {b, c} | {b+e, c+f}, t**(P+Q)*(t**(-k) - t**k))
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}
