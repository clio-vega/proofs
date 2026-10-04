"""Three further tests.

(T1) The dual-basis identity:  [alpha_{p_1} (x) ... (x) alpha_{p_m}] T_{w/u}
     = (M_{p_1} ... M_{p_m})_{u,w}.
     This is a REIMPLEMENTATION check of Side A, not independent corroboration:
     both sides enumerate the same labelled chains.  Recorded as such.

(T2) Samuel's identity at n=5 via transfer matrices, exact, all label words.

(T3) Is the relation (a,b)+(b,c)=(a,c) load-bearing?  Same chain sum in the FREE
     vector space on transpositions.
"""
import sys, random
from itertools import permutations, product
from collections import Counter
from schub import length, structure_constants
from chains import T_tensor, T_free
from monk import monk_matrix

def matmulvec(M, v):
    out = [0]*len(v)
    for i, vi in enumerate(v):
        if vi:
            Mi = M[i]
            for j, m in enumerate(Mi):
                if m:
                    out[j] += vi*m
    return out

def T1(n):
    perms = list(permutations(range(1, n+1)))
    Ms = {p: monk_matrix(n, p)[0] for p in range(1, n)}
    idx = {w: i for i, w in enumerate(perms)}
    bad = tested = 0
    for u in perms:
        for w in perms:
            m = length(w)-length(u)
            if m < 0:
                continue
            T = T_tensor(u, w, n)
            for pp in product(range(1, n), repeat=m):
                vec = [0]*len(perms); vec[idx[u]] = 1
                for p in pp:
                    vec = matmulvec(Ms[p], vec)
                tested += 1
                if vec[idx[w]] != T.get(pp, 0):
                    bad += 1
                    if bad == 1:
                        print("   MISMATCH", u, w, pp, vec[idx[w]], T.get(pp, 0))
    print(f"(T1) n={n}: dual-basis coefficient = Monk transfer-matrix entry: "
          f"{tested-bad}/{tested} label words agree, {bad} mismatches")

def T2(n, mmax=None):
    C, _ = structure_constants(n, n+1)
    perms = list(permutations(range(1, n+1)))
    idx = {w: i for i, w in enumerate(perms)}
    e = tuple(range(1, n+1))
    Ms = {p: monk_matrix(n, p)[0] for p in range(1, n)}
    byl = {}
    for v in perms:
        byl.setdefault(length(v), []).append(v)
    maxm = max(length(w) for w in perms)
    if mmax is not None:
        maxm = min(maxm, mmax)
    bad = 0; words = 0; triples = set()
    for m in range(0, maxm+1):
        # row vectors e_u * M_{p1}...M_{pm} for every u, built by one sweep per word
        for pp in product(range(1, n), repeat=m):
            words += 1
            rows = {}
            for u in perms:
                vec = [0]*len(perms); vec[idx[u]] = 1
                for p in pp:
                    vec = matmulvec(Ms[p], vec)
                rows[u] = vec
            base = rows[e]
            for u in perms:
                for w in perms:
                    if length(w)-length(u) != m:
                        continue
                    triples.add((u, w))
                    lhs = rows[u][idx[w]]
                    rhs = sum(C.get((u, v, w), 0)*base[idx[v]] for v in byl.get(m, []))
                    if lhs != rhs:
                        bad += 1
                        if bad == 1:
                            print("   MISMATCH", u, w, pp, lhs, rhs)
    print(f"(T2) n={n}, m<={maxm}: Samuel's identity, exact, all label words: "
          f"{words} words x {len(triples)} (u,w) pairs, {bad} disagreements")

def T3(n):
    """load-bearing test for the relation (a,b)+(b,c)=(a,c)"""
    C, _ = structure_constants(n, n+1)
    perms = list(permutations(range(1, n+1)))
    e = tuple(range(1, n+1))
    byl = {}
    for v in perms:
        byl.setdefault(length(v), []).append(v)
    Tfree = {v: T_free(e, v, n) for v in perms}
    ok = bad = 0; witness = None
    for u in perms:
        for w in perms:
            m = length(w)-length(u)
            if m < 1:
                continue
            lhs = T_free(u, w, n)
            rhs = Counter()
            for v in byl.get(m, []):
                c = C.get((u, v, w), 0)
                if c:
                    for t, k in Tfree[v].items():
                        rhs[t] += c*k
            if Counter({k: v for k, v in lhs.items() if v}) == Counter({k: v for k, v in rhs.items() if v}):
                ok += 1
            else:
                bad += 1
                if witness is None and m >= 1:
                    witness = (u, w, dict(lhs), dict(rhs))
    print(f"(T3) n={n}: SAME identity in the FREE space on transpositions "
          f"(relation NOT imposed): agree {ok}, disagree {bad}")
    if witness:
        u, w, l, r = witness
        print("   smallest witness:  u =", u, " w =", w)
        print("     free LHS =", l)
        print("     free RHS =", r)

if __name__ == "__main__":
    for n in (3, 4):
        T1(n)
    print()
    T3(3); T3(4)
