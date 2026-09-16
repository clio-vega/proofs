"""
Independent re-implementation of the ribbon shape-weight operator R_e^W, written
from Q149 Convention 2.3 (conv:def) ALONE, to re-verify the Q150 c1 counterexample
without reusing any of the c1/Q149 engine.

Maya set M is represented as (frozenset S, L) meaning  M = {n : n < -L}  u  S,
with S a subset of [-L, inf).  All operations keep L fixed.
"""
from itertools import product
from collections import defaultdict

L = 40   # everything strictly below -L is occupied

def maya_of_partition(lam):
    """M(lam) = {lam_j - j : j>=1}, as (S,L)."""
    lam = list(lam)
    n = len(lam) + L + 1
    lam = lam + [0]*(n-len(lam))
    S = set()
    for j in range(1, n+1):
        v = lam[j-1] - j
        if v >= -L:
            S.add(v)
    return frozenset(S)

def inM(S, n):
    return n < -L or n in S

def word(S, b, e):
    """occupancy word u_i = [b+i in M], i=1..e-1"""
    return tuple(1 if inM(S, b+i) else 0 for i in range(1, e))

def apply_R(vec, e, W):
    """vec: dict S -> coeff.  Returns R_e^W vec."""
    out = defaultdict(int)
    for S, c in vec.items():
        if c == 0: continue
        # legal b: b in M, b+e not in M.  b+e not in M forces b+e >= -L.
        for b in range(-L-e, L+2):
            if not inM(S, b): continue
            if inM(S, b+e): continue
            w = W[word(S, b, e)]
            if w == 0: continue
            T = set(S)
            if b >= -L: T.discard(b)
            T.add(b+e)
            out[frozenset(T)] += c*w
    return {k:v for k,v in out.items() if v != 0}

def commutator(vec, e, W, f, Wbar):
    a = apply_R(apply_R(vec, e, W), f, Wbar)
    b = apply_R(apply_R(vec, f, Wbar), e, W)
    out = defaultdict(int)
    for k,v in a.items(): out[k] += v
    for k,v in b.items(): out[k] -= v
    return {k:v for k,v in out.items() if v != 0}

def alpha(u, e):
    """Q149 lem:dict: composition of the ribbon, top row first."""
    c = [0] + [i for i in range(1, e) if u[i-1] == 1] + [e]
    return tuple(c[j]-c[j-1] for j in range(len(c)-1, 0, -1))

def W_from_compositions(e, table, default=0):
    """table: dict composition -> value."""
    W = {}
    for u in product((0,1), repeat=e-1):
        W[u] = table.get(alpha(u, e), default)
    return W

def partitions(n):
    if n == 0:
        yield ()
        return
    def rec(n, mx):
        if n == 0:
            yield ()
            return
        for k in range(min(n, mx), 0, -1):
            for rest in rec(n-k, k):
                yield (k,) + rest
    yield from rec(n, n)

def all_partitions_upto(N):
    out = []
    for n in range(N+1):
        out.extend(partitions(n))
    return out
