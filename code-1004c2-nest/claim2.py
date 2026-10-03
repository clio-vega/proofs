"""CLAIM 2, tested in the abstract class it is stated in:
   Omega M-convex in Z^N, A subset of N  ==>  a -> #{z in Omega : sum_A z = a} log-concave?

Exhaustive over ALL M-convex subsets of the simplex {z in Z_{>=0}^N : sum z = D}
for small N, D.  (A subset S of a hyperplane is M-convex iff for all z,z' in S and
all p with z_p>z'_p there is q with z_q<z'_q and z-e_p+e_q, z'+e_p-e_q both in S.)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen import is_pf2
from itertools import product, combinations
from collections import Counter

def simplex(N, D):
    out = []
    def rec(i, rem, cur):
        if i == N-1: out.append(tuple(cur+[rem])); return
        for t in range(rem+1): rec(i+1, rem-t, cur+[t])
    rec(0, D, [])
    return out

def is_Mconvex(S):
    Sset = set(S); N = len(S[0])
    for z in S:
        for zp in S:
            for p in range(N):
                if z[p] <= zp[p]: continue
                ok = False
                for q in range(N):
                    if z[q] >= zp[q]: continue
                    a = list(z); a[p] -= 1; a[q] += 1
                    b = list(zp); b[p] += 1; b[q] -= 1
                    if tuple(a) in Sset and tuple(b) in Sset: ok = True; break
                if not ok: return False
    return True

def fibre_counts(S, A):
    c = Counter()
    for z in S: c[sum(z[i] for i in A)] += 1
    lo, hi = min(c), max(c)
    return [c.get(i, 0) for i in range(lo, hi+1)]

if __name__ == '__main__':
    for N, D in ((3,3),(3,4),(3,5),(4,3)):
        pts = simplex(N, D); n = len(pts)
        st = Counter(); wit = []
        for mask in range(1, 1 << n):
            S = [pts[i] for i in range(n) if mask >> i & 1]
            if len(S) < 2: continue
            if not is_Mconvex(S): continue
            st['Mconvex'] += 1
            for r in range(1, N):
                for A in combinations(range(N), r):
                    f = fibre_counts(S, A)
                    st['pf2' if is_pf2(f) else 'FAIL'] += 1
                    if not is_pf2(f) and len(wit) < 4: wit.append((N, D, sorted(S), A, f))
        print(f"N={N} D={D}: {dict(st)}")
        for w in wit: print("   FAIL", w)
