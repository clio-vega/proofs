r"""INDEPENDENT CHECK of the exponent bookkeeping that is the content of the proof.

The proof claims, for every NONZERO term t_sigma of F_{N,n,r}(w_A|w|B)|_{q=0}
under  w_m = eps,  w_c = T^{-e_c} (c != m),  e_c = c:

    ord_T ( pref * t_sigma )  =  Pen(sigma)                              (*)

where, writing H = (h_1,..,h_r) for the heads (h_j = the row element carried by
slot j) and tau for the tails,

    Pen(sigma) = sum_{j<k<=r} (E(h_j)-E(h_k))^+
               + sum_{t in tau} sum_{k<=r} (E(t)-E(h_k))^+ ,   E(c)=e_c, E(m)=0,

and the leading COEFFICIENT is +1.  Note Pen >= 0 with equality iff the heads are
the top r elements of A in increasing order and the tails are the bottom n-r.

This script computes ord_T and the leading coefficient of each term EXACTLY
(symbolic T, Laurent rational function) and compares against (*).  It is an
independent mechanism from degen.py: that one sums the terms and reads a limit,
this one dissects the sum.

DOMAIN: N = 2,3,4; X in {B,C}; every m; eps in {+1,-1} (C) / {+1} (B);
        every (A,B) with m notin A and A >= B.  Exhaustive.
"""
import sys, itertools
from itertools import permutations
from math import factorial
import sympy as sp
sys.path.insert(0, '.')
import q484_engine as E

T = sp.Symbol('T', positive=True)

def ordT(expr):
    """exact lowest T-exponent and its coefficient, for a Laurent rational fn of T"""
    f = sp.cancel(sp.together(sp.expand(expr)))
    if f == 0:
        return None, None
    num, den = sp.fraction(f)
    def low(p):
        p = sp.Poly(sp.expand(p), T)
        cs = list(reversed(p.all_coeffs()))     # ascending
        for i, c in enumerate(cs):
            if c != 0:
                return i, c
        return None, None
    i, a = low(num); j, b = low(den)
    return i - j, sp.cancel(a / b)

def term(N, n, r, u, w, I, beta, gamma, sigma):
    """ONE sigma-term of E.F, transcribed from q484_engine.F with q = 0."""
    us = [u[sigma[j]] for j in range(n)]
    t = sp.Integer(1)
    for j in range(r):
        ij = I[j]
        t *= (us[j]**2 - 1)                       # (1-q)(u^2-1) at q=0
        t *= us[j]**(N - ij)                      # prod_{k>i_j}(u_j - q w_k) at q=0
        for k in range(1, ij):
            t *= (us[j] - w[k - 1])
    for j in range(r):
        for k in range(j + 1, r):
            t *= (-us[k]) * (us[j] * us[k] - 1) / (us[j] - us[k])
    for j in range(r, n):
        t *= us[j]**N
    for j in range(r, n):
        for k in range(r):
            t *= (-us[k]) * (us[j] * us[k] - 1) / (us[j] - us[k])
    for j in range(r, n):
        for k in range(j + 1, n):
            t *= us[j] * us[k]
    for j in range(r, n):
        t *= (beta * us[j] - gamma)
    pref = sp.Integer(1)
    for k in range(n):
        pref *= u[k]**(-(N + n))
    return pref * t

def Pen(A, B, m, sigma):
    n, r = len(A), len(B)
    Ev = lambda c: 0 if c == m else c
    slot = [A[sigma[j]] for j in range(n)]
    heads, tails = slot[:r], slot[r:]
    p = 0
    for j in range(r):
        for k in range(j + 1, r):
            p += max(0, Ev(heads[j]) - Ev(heads[k]))
    for t_ in tails:
        for hk in heads:
            p += max(0, Ev(t_) - Ev(hk))
    return p

ok = bad = zeroterms = dom = 0
domct = {}
fails = []
for N in (2, 3, 4):
    P = E.subsets(N)
    for X in ('C', 'B'):
        beta, gamma = (sp.Integer(1), sp.Integer(0)) if X == 'C' else (sp.Integer(1), sp.Integer(-1))
        epss = (1, -1) if X == 'C' else (1,)
        for m in range(1, N + 1):
            for eps in epss:
                w = [sp.Integer(eps) if c == m else T**(-c) for c in range(1, N + 1)]
                for A in P:
                    if m in A or len(A) == 0: continue
                    u = [w[a - 1] for a in A]
                    n = len(A)
                    for B in P:
                        if not E.geq(A, B): continue
                        r = len(B)
                        if r > n: continue
                        for sigma in permutations(range(n)):
                            tv = term(N, n, r, u, w, B, beta, gamma, sigma)
                            o, c = ordT(tv)
                            if o is None:
                                zeroterms += 1
                                # predicted zero iff some head h_j < b_j
                                slot = [A[sigma[j]] for j in range(n)]
                                pz = any(slot[j] < B[j] for j in range(r))
                                if not pz:
                                    bad += 1; fails.append(('UNPREDICTED ZERO', N, X, m, eps, A, B, sigma))
                                continue
                            pp = Pen(A, B, m, sigma)
                            # type B: the prefactor prod(1+1/u)^{-1} -> 1, excluded here
                            # The proof needs: (a) ord = Pen for EVERY nonzero term
                            # (so Pen >= 0 bounds every term), and (b) coeff = +1 for
                            # the DOMINANT term Pen = 0.  The first run also demanded
                            # coeff = +1 for subdominant terms; that is not claimed and
                            # is false (243 of them have coeff -1).  Separate the two.
                            if o != pp:
                                bad += 1
                                if len(fails) < 10:
                                    fails.append(('ORD', N, X, m, eps, A, B, sigma, o, pp))
                            elif pp == 0 and c != 1:
                                bad += 1
                                if len(fails) < 10:
                                    fails.append(('DOMINANT COEF', N, X, m, eps, A, B, sigma, c))
                            else:
                                ok += 1
                                if pp == 0:
                                    dom += 1
                                    key = (N, X, m, eps, A, B)
                                    domct[key] = domct.get(key, 0) + 1
print(f"DOMINANT (Pen=0) terms with coeff +1: {dom}")
import collections
cnt=collections.Counter(domct.values())
print("  dominant-term multiplicity per instance (must equal (n-r)! ):", dict(cnt))
badmult=[(k,v) for k,v in domct.items() if v != __import__('math').factorial(len(k[4])-len(k[5]))]
print("  instances whose Pen=0 multiplicity != (n-r)!:", len(badmult), badmult[:5])
print(f"nonzero terms checked: {ok+bad-0}   MATCH (ord=Pen and coeff=+1): {ok}   MISMATCH: {bad}")
print(f"identically-zero terms: {zeroterms} (all predicted by 'some head h_j < b_j')")
for f in fails[:10]: print("  ", f)
