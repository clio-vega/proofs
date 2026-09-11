"""
COMPARATOR / answer key (engine E2): build M^{(mu)}(t) directly from rim-hook
combinatorics on Young diagrams -- NO abacus, NO hook lengths -- and read off the
rows with exactly two nonzero entries.  This is the source's own calibration
(prop:reflect "span exactly the same space as ... all rows of M having exactly
two nonzero entries", all mu, n<=8).

A row lambda with exactly two nonzero entries t^a (col g1), t^b (col g2) gives
   t^a w(g1) + t^b w(g2) = 0   <=>   w(g2) = -t^{a-b} w(g1).
"""
from graph import partitions

def cells(mu): return {(i+1, j+1) for i in range(len(mu)) for j in range(mu[i])}

def is_rim_hook(lam, gam):
    """lam/gam a rim hook?  connected, no 2x2.  Assumes gam subset lam."""
    S = cells(lam) - cells(gam)
    if not S: return None
    for (i, j) in S:                                   # no 2x2
        if (i,j+1) in S and (i+1,j) in S and (i+1,j+1) in S: return None
    start = next(iter(S)); seen = {start}; st = [start]
    while st:
        (i,j) = st.pop()
        for c in ((i+1,j),(i-1,j),(i,j+1),(i,j-1)):
            if c in S and c not in seen: seen.add(c); st.append(c)
    if seen != S: return None
    return len({i for (i,j) in S}) - 1                 # height

def contains(lam, gam):
    if len(gam) > len(lam): return False
    return all(gam[i] <= lam[i] for i in range(len(gam)))

def matrix(mu):
    n = sum(mu)
    cols = []
    for L in range(1, n+1):
        for g in partitions(n-L) if n-L > 0 else [()]:
            if contains(mu, g) and is_rim_hook(mu, g) is not None:
                cols.append(g)
    rows = list(partitions(n))
    M = {}
    for lam in rows:
        for g in cols:
            if contains(lam, g):
                h = is_rim_hook(lam, g)
                if h is not None: M[(lam, g)] = h     # entry t^h
    return rows, cols, M

def two_term_relations(mu):
    rows, cols, M = matrix(mu)
    rels = []
    for lam in rows:
        nz = [(g, M[(lam, g)]) for g in cols if (lam, g) in M]
        if len(nz) == 2:
            (g1, a), (g2, b) = nz
            rels.append((g1, g2, a - b))              # w(g2) = -t^{a-b} w(g1)
    return cols, rels
