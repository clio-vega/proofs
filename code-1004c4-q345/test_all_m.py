"""Exact test of Samuel's identity for ALL label words of ALL lengths at once.

Key observation: M_p has support only on pairs with l(w)-l(u) = 1, so the span
A_m of all products M_{p_1}...M_{p_m} is supported on l(w)-l(u) = m.  Samuel's
identity is LINEAR in the matrix M_p, so it holds for every word of length m iff
it holds on a BASIS of A_m.  A_m is built inductively: A_m = span( A_{m-1} . M_p ).
dim A_m is tiny (it is the degree-m part of H^*(Fl_n) acting on itself).
"""
import sys
from fractions import Fraction
from itertools import permutations
from schub import length, structure_constants
from monk import monk_matrix

def reduce_rows(vectors):
    """row-reduce a list of dict-vectors over Q; return basis (list of dicts)."""
    basis = []          # list of (pivot_key, dict)
    for v in vectors:
        v = {k: Fraction(c) for k, c in v.items() if c}
        for pk, bv in basis:
            if pk in v:
                f = v[pk]/bv[pk]
                for k, c in bv.items():
                    v[k] = v.get(k, 0) - f*c
                    if v[k] == 0:
                        del v[k]
        if v:
            pk = min(v)
            basis.append((pk, v))
    return [bv for _, bv in basis]

def run(n):
    C, _ = structure_constants(n, n+1)
    perms = list(permutations(range(1, n+1)))
    e = tuple(range(1, n+1))
    byl = {}
    for v in perms:
        byl.setdefault(length(v), []).append(v)
    idx = {w: i for i, w in enumerate(perms)}
    Ms = {}
    for p in range(1, n):
        M, _, _ = monk_matrix(n, p)
        Ms[p] = {(u, v): M[idx[u]][idx[v]]
                 for u in perms for v in perms if M[idx[u]][idx[v]]}
    maxm = max(length(w) for w in perms)
    # A_0 = identity
    cur = [{(u, u): 1 for u in perms}]
    total_bad = 0
    print(f"  n={n}: #(u,w) pairs and basis dimension per degree m")
    for m in range(0, maxm+1):
        pairs = [(u, w) for u in perms for w in perms if length(w)-length(u) == m]
        bad = 0
        for X in cur:
            for (u, w) in pairs:
                lhs = X.get((u, w), 0)
                rhs = sum(C.get((u, v, w), 0)*X.get((e, v), 0) for v in byl.get(m, []))
                if lhs != rhs:
                    bad += 1
                    if total_bad+bad == 1:
                        print("   MISMATCH m=", m, u, w, lhs, rhs)
        total_bad += bad
        print(f"    m={m:2d}:  pairs {len(pairs):5d}   dim A_m = {len(cur):3d}"
              f"   checks {len(cur)*len(pairs):7d}   bad {bad}")
        if m == maxm:
            break
        nxt = []
        for X in cur:
            for p in range(1, n):
                Y = {}
                Mp = Ms[p]
                for (a, b), c in X.items():
                    for (b2, d), c2 in Mp.items():
                        if b2 == b:
                            Y[(a, d)] = Y.get((a, d), 0) + c*c2
                nxt.append({k: v for k, v in Y.items() if v})
        cur = reduce_rows(nxt)
        if not cur:
            break
    print(f"  n={n}: TOTAL disagreements across all degrees = {total_bad}")
    return total_bad

if __name__ == "__main__":
    for n in (3, 4, 5):
        run(n)
        print()
