# NOTE: cyl.py is VENDORED into this directory (copied from
# projects/scratch/q191/code/cyl.py, the code behind the 2026-09-20 cylindric
# M-convexity paper) so that everything committed here runs without reaching
# outside the repository.  projects/scratch/ is not under version control.
"""
Cylindric shapes in C^{n,m} as n-periodic bead (Maya) configurations.

A cylindric shape is a bi-infinite strictly increasing sequence (x_i)_{i in Z}
with x_{i+m} = x_i + n.  Stored as the tuple (x_1,...,x_m), x_1<...<x_m<x_1+n.

  nu subset rho            <=>  x_i(nu) <= x_i(rho)  for all i
  rho/nu horizontal strip  <=>  x_i(nu) <= x_i(rho) < x_{i+1}(nu)  for all i
  |rho/nu| (per period)    =    sum_{i=1}^m (x_i(rho) - x_i(nu))
"""
from itertools import product
from functools import lru_cache


def ext(x, i, n, m):
    """x_i for arbitrary i>=1 using x_{i+m}=x_i+n; x indexed 0-based tuple."""
    q, r = divmod(i - 1, m)
    return x[r] + q * n


def is_shape(x, n, m):
    return all(x[i] < x[i + 1] for i in range(m - 1)) and x[m - 1] < x[0] + n


def contains(rho, nu):
    return all(a >= b for a, b in zip(rho, nu))


def hstrip(nu, rho, n, m):
    """rho/nu a cylindric horizontal strip?"""
    for i in range(m):
        if not (nu[i] <= rho[i]):
            return False
        nxt = nu[i + 1] if i + 1 < m else nu[0] + n
        if not (rho[i] < nxt):
            return False
    return True


def size(nu, rho):
    return sum(a - b for a, b in zip(rho, nu))


def interval_shapes(mu, lam, n, m):
    """all cylindric shapes kappa with mu <= kappa <= lam"""
    rngs = [range(mu[i], lam[i] + 1) for i in range(m)]
    out = []
    for x in product(*rngs):
        if is_shape(x, n, m):
            out.append(x)
    return out


def weights(mu, lam, n, m, ell):
    """set of weight compositions (alpha_1..alpha_ell) of cylindric SSYT of shape lam/mu"""
    S = interval_shapes(mu, lam, n, m)
    idx = {s: j for j, s in enumerate(S)}
    # adjacency: hstrips
    adj = [[] for _ in S]
    for a, sa in enumerate(S):
        for b, sb in enumerate(S):
            if hstrip(sa, sb, n, m):
                adj[a].append((b, size(sa, sb)))
    # DP over chains of length ell from mu to lam
    cur = {(idx[mu],): ()}
    states = {idx[mu]: {()}}
    for t in range(ell):
        nxt = {}
        for a, ws in states.items():
            for b, sz in adj[a]:
                nxt.setdefault(b, set()).update(w + (sz,) for w in ws)
        states = nxt
    return states.get(idx[lam], set())


# ---------- symmetric-function bookkeeping ----------
def sort_part(alpha):
    return tuple(sorted([a for a in alpha if a > 0], reverse=True))


def dominates(lam, mu):
    """lam >= mu in dominance (same size)"""
    if sum(lam) != sum(mu):
        return False
    s1 = s2 = 0
    L = max(len(lam), len(mu))
    for i in range(L):
        s1 += lam[i] if i < len(lam) else 0
        s2 += mu[i] if i < len(mu) else 0
        if s2 > s1:
            return False
    return True


def partitions(d, maxlen=None):
    res = []
    def rec(rem, mx, cur):
        if rem == 0:
            res.append(tuple(cur)); return
        for p in range(min(rem, mx), 0, -1):
            if maxlen is not None and len(cur) + 1 > maxlen:
                break
            rec(rem - p, p, cur + [p])
    rec(d, d, [])
    return res


def is_M_convex_symmetric(P, d, ell):
    """P = set of partitions occurring; W = all rearrangements in N^ell.
       Return (is_ideal, maximal_elements)."""
    allp = [p for p in partitions(d, maxlen=ell)]
    ideal = all(q in P for p in P for q in allp if dominates(p, q))
    maxl = [p for p in P if not any(q != p and dominates(q, p) for q in P)]
    return ideal, maxl
