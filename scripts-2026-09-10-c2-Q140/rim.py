"""Rim-hook combinatorics + the Khanna-Loehr local system (Q140).

Two independent implementations of "remove an L-rim hook from lambda":
  (1) beta-number / abacus:  beta -> beta - L
  (2) brute force over Young diagrams: connectivity + no 2x2
They are cross-checked before anything else is computed.
"""
from itertools import combinations
import sympy as sp

t = sp.symbols('t')

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k):
            yield (k,) + rest

# ---------- (1) abacus ----------
def beta_set(lam, N):
    """first-column hook lengths padded to N beads: lam_i + N - i, i=1..N"""
    lam = list(lam) + [0]*(N - len(lam))
    return sorted(lam[i] + (N-1-i) for i in range(N))

def from_beta(beta):
    beta = sorted(beta)
    N = len(beta)
    lam = [beta[i] - i for i in range(N)]
    lam = sorted([x for x in lam if x > 0], reverse=True)
    return tuple(lam)

def remove_rim_abacus(lam, L, N=None):
    """all (gamma, height) from removing one L-rim hook, via beta numbers"""
    if N is None: N = sum(lam) + len(lam) + 2
    B = beta_set(lam, N)
    Bs = set(B)
    out = []
    for b in B:
        if b - L >= 0 and (b - L) not in Bs:
            newB = (Bs - {b}) | {b - L}
            ht = sum(1 for x in Bs if b - L < x < b)
            out.append((from_beta(newB), ht))
    return out

# ---------- (2) brute force ----------
def cells(lam):
    return set((i, j) for i, r in enumerate(lam) for j in range(r))

def remove_rim_brute(lam, L):
    n = sum(lam)
    out = []
    for gam in partitions(n - L):
        cg, cl = cells(gam), cells(lam)
        if not cg <= cl: continue
        skew = cl - cg
        if len(skew) != L: continue
        # connected (edge-adjacency)?
        start = next(iter(skew)); seen = {start}; stack = [start]
        while stack:
            (i, j) = stack.pop()
            for (a, b) in ((i+1,j),(i-1,j),(i,j+1),(i,j-1)):
                if (a,b) in skew and (a,b) not in seen:
                    seen.add((a,b)); stack.append((a,b))
        if seen != skew: continue
        # no 2x2
        if any((i,j) in skew and (i+1,j) in skew and (i,j+1) in skew and (i+1,j+1) in skew
               for (i,j) in skew): continue
        rows = set(i for (i,j) in skew)
        out.append((gam, len(rows)-1))
    return out

def crosscheck(maxn=7):
    bad = 0; tot = 0
    for n in range(1, maxn+1):
        for lam in partitions(n):
            for L in range(1, n+1):
                a = sorted(remove_rim_abacus(lam, L))
                b = sorted(remove_rim_brute(lam, L))
                tot += 1
                if a != b:
                    bad += 1
                    print("MISMATCH", lam, L, a, b)
    print(f"crosscheck: {tot-bad}/{tot} (lambda,L) pairs agree, |lambda|<={maxn}")
    return bad == 0

# ---------- the local system ----------
def predecessors(mu):
    """C(mu) = {gamma : mu/gamma is a rim hook}, as dict gamma -> height"""
    n = sum(mu)
    d = {}
    for L in range(1, n+1):
        for g, h in remove_rim_abacus(mu, L):
            assert g not in d
            d[g] = h
    return d

def local_system(mu):
    n = sum(mu)
    rows = list(partitions(n))
    cols = sorted(predecessors(mu).keys(), key=lambda g: (sum(g), g))
    M = sp.zeros(len(rows), len(cols))
    for i, lam in enumerate(rows):
        for j, gam in enumerate(cols):
            L = n - sum(gam)
            if L == 0: continue
            for g, h in remove_rim_abacus(lam, L):
                if g == gam:
                    M[i, j] = t**h
    b = sp.zeros(len(rows), 1)
    b[rows.index(tuple(mu)), 0] = 1
    return rows, cols, M, b
