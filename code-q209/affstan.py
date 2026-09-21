"""Affine Stanley symmetric polynomials straight from (def-affine-stanley), WZZ 2401.14632 l.797.

Affine symmetric group S~_n (n = k+1) in window notation: a tuple (w(1),...,w(n)) of integers,
distinct mod n, with sum = n(n+1)/2.  Right multiplication by s_i permutes POSITIONS.
"""
from itertools import combinations
from functools import lru_cache

# ---------------------------------------------------------------- group ops

def identity(n):
    return tuple(range(1, n + 1))

def rmul_s(w, i, n):
    """w * s_i  (i in 0..n-1).  s_i swaps positions i, i+1; s_0 swaps positions 0 and 1,
    where w(0) = w(n) - n and w(n+1) = w(1) + n."""
    v = list(w)
    if i == 0:
        v[0], v[n - 1] = w[n - 1] - n, w[0] + n
    else:
        v[i - 1], v[i] = w[i], w[i - 1]
    return tuple(v)

def length(w, n):
    """Affine length: sum over 1<=i<j<=n of |floor((w(j)-w(i))/n)|."""
    tot = 0
    for i in range(n):
        for j in range(i + 1, n):
            tot += abs((w[j] - w[i]) // n)
    return tot

# ------------------------------------------------- cyclically decreasing elements

@lru_cache(maxsize=None)
def cyc_dec_elements(n):
    """dict: frozenset S (proper subset of Z/n) -> (window, |S|).
    u_S = product of s_i, i in S, taken in cyclically decreasing order.
    Order: pick j not in S, read j-1, j-2, ..., j-n+1 mod n, keep those in S."""
    out = {}
    for size in range(n):                       # size <= n-1 : S proper
        for S in combinations(range(n), size):
            Sf = frozenset(S)
            j = next(a for a in range(n) if a not in Sf)
            word = [(j - t) % n for t in range(1, n) if (j - t) % n in Sf]
            w = identity(n)
            for i in word:
                w = rmul_s(w, i, n)
            out[Sf] = (w, size, tuple(word))
    return out

def check_cyc_dec_welldefined(n):
    """u_S must not depend on which j notin S we start from."""
    bad = []
    for size in range(n):
        for S in combinations(range(n), size):
            Sf = frozenset(S)
            wins = set()
            for j in range(n):
                if j in Sf:
                    continue
                word = [(j - t) % n for t in range(1, n) if (j - t) % n in Sf]
                w = identity(n)
                for i in word:
                    w = rmul_s(w, i, n)
                if length(w, n) != size:
                    bad.append(("not reduced", S, j))
                wins.add(w)
            if len(wins) > 1:
                bad.append(("ambiguous", S, wins))
    return bad

# ------------------------------------------------------------- the definition

def affine_stanley_support(w, n, r):
    """Set of exponent vectors alpha=(l(w^1),...,l(w^r)) over cyclically decreasing
    length-additive factorisations w = w^1...w^r.  Returns dict alpha -> coefficient."""
    cd = cyc_dec_elements(n)
    lw = length(w, n)
    # DP forwards over prefixes: state = (window, length)
    layer = {identity(n): {(): 1}}
    for step in range(r):
        new = {}
        for u, paths in layer.items():
            lu = length(u, n)
            for Sf, (uS, size, _word) in cd.items():
                if lu + size > lw:
                    continue
                v = u
                ok = True
                for i in _word:
                    v2 = rmul_s(v, i, n)
                    if length(v2, n) != length(v, n) + 1:
                        ok = False
                        break
                    v = v2
                if not ok:
                    continue
                tgt = new.setdefault(v, {})
                for a, c in paths.items():
                    key = a + (size,)
                    tgt[key] = tgt.get(key, 0) + c
        layer = new
    return layer.get(w, {})

def elements_of_length(n, L):
    """BFS all affine permutations of length <= L, grouped by length."""
    byl = {0: [identity(n)]}
    seen = {identity(n): 0}
    frontier = [identity(n)]
    for d in range(1, L + 1):
        nxt = []
        for w in frontier:
            for i in range(n):
                v = rmul_s(w, i, n)
                if v not in seen and length(v, n) == d:
                    seen[v] = d
                    nxt.append(v)
        byl[d] = nxt
        frontier = nxt
    return byl
