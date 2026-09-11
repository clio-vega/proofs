"""
DERIVER (engine E1b): build G(mu) from the HOOK-LENGTH description proved in
the write-up (Theorem 1).  Independent of the abacus code in graph.py.

Vertices: cells (i,k) of mu, 1<=i<=ell, 1<=k<=mu_i.
Edges:  (i,k)~(i,m), k<m  iff  h[i][k]+h[i][m] is NOT a hook length of row i.
        (i,k)~(i',k), i<i' iff  h[i][k]+h[i'][k] is NOT a hook length of column k.
"""
from itertools import combinations

def conj(mu):
    if not mu: return ()
    return tuple(sum(1 for x in mu if x >= j) for j in range(1, mu[0]+1))

def hooks(mu):
    mp = conj(mu)
    return {(i+1, j+1): mu[i]-(j+1) + mp[j]-(i+1) + 1
            for i in range(len(mu)) for j in range(mu[i])}

def hook_edges(mu):
    h = hooks(mu); mp = conj(mu); E = []
    for i in range(1, len(mu)+1):                       # rows
        Hi = {h[(i, j)] for j in range(1, mu[i-1]+1)}
        for k, m in combinations(range(1, mu[i-1]+1), 2):
            if h[(i,k)] + h[(i,m)] not in Hi:
                E.append(((i,k), (i,m), 'row'))
    for k in range(1, (mu[0] if mu else 0)+1):          # columns
        Hk = {h[(i, k)] for i in range(1, mp[k-1]+1)}
        for i, ip in combinations(range(1, mp[k-1]+1), 2):
            if h[(i,k)] + h[(ip,k)] not in Hk:
                E.append(((i,k), (ip,k), 'col'))
    return E
