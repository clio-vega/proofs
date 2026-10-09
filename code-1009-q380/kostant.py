"""ENGINE B: Kostka-Foulkes via the Lusztig t-analogue of Kostant's partition
function.  Macdonald III, Section 6, Example 4 (deep-read 2026-10-09):

    K_{lam,mu}(t) = sum_{w in S_n} eps(w) * P( w^{-1}(lam+delta) - (mu+delta) ; t)

with delta = (n-1, n-2, ..., 1, 0), R^+ = {e_i - e_j : i < j}, and

    P(xi; t) = sum over families (m_alpha)_{alpha in R^+} of non-negative integers
               with xi = sum m_alpha alpha,  of  t^{sum m_alpha}.

Macdonald's Example 4 also records: P(xi;t) != 0 iff xi = sum_i eta_i (e_i - e_{i+1})
with all eta_i >= 0, i.e. iff all partial sums xi_1 + ... + xi_i >= 0 (and total 0);
and then P is monic of degree sum eta_i = <xi, delta>.

This mechanism is disjoint from charge: it is a Weyl-group alternating sum over a
root-lattice partition count.  No tableau, no reading word, no charge is used.
"""
from itertools import permutations
from functools import lru_cache

# positive roots of A_{n-1} as (i,j), i<j, meaning e_i - e_j


def pfun(xi):
    """t-analogue of Kostant's partition function, returned as dict exp->coeff.
    xi: tuple of integers summing to 0, length n."""
    n = len(xi)
    # necessary condition: all partial sums >= 0
    s = 0
    for v in xi:
        s += v
        if s < 0:
            return {}
    if s != 0:
        return {}
    # eta_i = xi_1 + ... + xi_i  for i=1..n-1 ;  xi = sum eta_i (e_i - e_{i+1})
    # Count fillings: recursive over simple-root coordinates.
    # Direct DP: process coordinates left to right, tracking the "flow".
    # Writing xi = sum_{i<j} m_{ij} (e_i - e_j):  think of m_{ij} as a unit of mass
    # sent from i to j.  Equivalently: a non-negative integer matrix (m_{ij})_{i<j}
    # with row-minus-column sums = xi.  t exponent = sum m_{ij}.
    roots = [(i, j) for i in range(n) for j in range(i + 1, n)]

    res = {}

    def rec(k, rem, used):
        if k == len(roots):
            if all(v == 0 for v in rem):
                res[used] = res.get(used, 0) + 1
            return
        i, j = roots[k]
        # prune: coordinates < i can no longer be changed
        for a in range(i):
            if rem[a] != 0:
                return
        # m_{ij} can be 0..rem[i] (if rem[i] >= 0)
        if rem[i] < 0:
            return
        maxm = rem[i]
        for m in range(maxm + 1):
            nr = list(rem)
            nr[i] -= m
            nr[j] += m
            rec(k + 1, tuple(nr), used + m)

    rec(0, tuple(xi), 0)
    return res


def KB(lam, mu):
    """Kostka-Foulkes polynomial via Engine B. Returns dict exp->coeff."""
    n = max(len(lam), len(mu))
    lam = tuple(lam) + (0,) * (n - len(lam))
    mu = tuple(mu) + (0,) * (n - len(mu))
    delta = tuple(range(n - 1, -1, -1))
    lpd = tuple(lam[i] + delta[i] for i in range(n))
    mpd = tuple(mu[i] + delta[i] for i in range(n))
    acc = {}
    for perm in permutations(range(n)):
        # w^{-1}(lam+delta): we sum over all w, and w^{-1} ranges over all of S_n too,
        # with the same sign, so just permute lpd directly.
        sign = perm_sign(perm)
        xi = tuple(lpd[perm[i]] - mpd[i] for i in range(n))
        p = pfun(xi)
        for e, c in p.items():
            acc[e] = acc.get(e, 0) + sign * c
    return {e: c for e, c in acc.items() if c != 0}


def perm_sign(perm):
    n = len(perm)
    seen = [False] * n
    sign = 1
    for i in range(n):
        if seen[i]:
            continue
        L = 0
        j = i
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            L += 1
        if L % 2 == 0:
            sign = -sign
    return sign
