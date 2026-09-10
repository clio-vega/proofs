"""Core: partitions, e-ribbon addition/removal, R_e(t), e-core/e-weight.
Verified against Murnaghan-Nakayama (t=-1 must give p_e multiplication)."""
import sympy as sp
from sympy import Rational, Poly, symbols
from itertools import combinations
t = sp.Symbol('t')
z = sp.Symbol('z')

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k):
            yield (k,) + rest

def all_parts_upto(N):
    out = []
    for n in range(N+1): out += list(partitions(n))
    return out

def conj(lam):
    if not lam: return ()
    m = lam[0]
    return tuple(sum(1 for p in lam if p >= j) for j in range(1, m+1))

def beta_set(lam, L):
    """first-column hook lengths / beta-numbers with L beads: lam_i + L - i, i=1..L"""
    lam = list(lam) + [0]*(L - len(lam))
    return sorted(lam[i] + L - 1 - i for i in range(L))

def from_beta(bs, L):
    bs = sorted(bs)
    lam = [bs[i] - i for i in range(len(bs))]
    lam = sorted([x for x in lam if x > 0], reverse=True)
    return tuple(lam)

def add_ribbons(lam, e, L):
    """All (mu, height) with mu/lam a connected e-ribbon.  Abacus: move a bead
    from position b to b+e (must be empty).  height = #beads strictly between."""
    bs = set(beta_set(lam, L))
    out = []
    for b in sorted(bs):
        if b + e in bs: continue
        if b + e > max(bs) + e + 5: pass
        nb = (bs - {b}) | {b+e}
        ht = sum(1 for x in bs if b < x < b+e)
        out.append((from_beta(nb, L), ht))
    return out

def remove_ribbons(lam, e, L):
    bs = set(beta_set(lam, L))
    out = []
    for b in sorted(bs):
        if b - e < 0 or (b-e) in bs: continue
        nb = (bs - {b}) | {b-e}
        ht = sum(1 for x in bs if b-e < x < b)
        out.append((from_beta(nb, L), ht))
    return out

def ecore_weight(lam, e, L=None):
    if L is None: L = len(lam) + sum(lam) + 2
    cur = lam; w = 0
    while True:
        r = remove_ribbons(cur, e, L)
        if not r: return cur, w
        cur = r[0][0]; w += 1

# ---- R_e(t) as a matrix in the Schur basis, on partitions of size <= N ----
def R_matrix(e, N, tval=t):
    """dict (lam,mu) -> t^ht  where lam/mu is an e-ribbon, |lam|<=N."""
    L = N + 2
    M = {}
    for mu in all_parts_upto(N - e):
        for lam, ht in add_ribbons(mu, e, L):
            if sum(lam) <= N:
                M[(lam, mu)] = M.get((lam,mu), 0) + tval**ht
    return M
