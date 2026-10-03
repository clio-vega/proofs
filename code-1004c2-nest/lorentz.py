"""Independent numerical check of the LORENTZIAN ROUTE to (Q).

   P~_i(X,E,W) = sum_{x,e,w>=0, x+e+w=h_i, x+e>=l_i} X^x E^e W^w
   -- homogeneous of degree h_i, coefficients all 1, support = box cap hyperplane.
   k(a) = [X^a E^{D-a} W^{H-D}] prod_i P~_i,  H = sum h_i.

Lorentzianity of N(f) for f homogeneous of degree d in 3 variables is checked by
raw-hessian-lemma (PROVED, this registry): the Hessian of d^beta N(f) has entries
M_ij = c_{beta+e_i+e_j}, the RAW coefficients.  Then l3-det-reduction
(LEAN-VERIFIED): a symmetric 3x3 with nonnegative entries has at most one positive
eigenvalue IFF e_2(M) <= 0 and det(M) >= 0.  Exact integer arithmetic throughout.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from itertools import product
from collections import Counter

def Ptilde(l, h):
    """dict (x,e,w) -> 1"""
    return {(x, e, h-x-e): 1 for x in range(h+1) for e in range(h+1-x)
            if l <= x+e <= h}

def mul(f, g):
    out = Counter()
    for a, ca in f.items():
        for b, cb in g.items():
            out[(a[0]+b[0], a[1]+b[1], a[2]+b[2])] += ca*cb
    return dict(out)

def e2_det(M):
    e2 = (M[0][0]*M[1][1]-M[0][1]**2) + (M[0][0]*M[2][2]-M[0][2]**2) + (M[1][1]*M[2][2]-M[1][2]**2)
    det = (M[0][0]*(M[1][1]*M[2][2]-M[1][2]**2)
           - M[0][1]*(M[0][1]*M[2][2]-M[1][2]*M[0][2])
           + M[0][2]*(M[0][1]*M[1][2]-M[1][1]*M[0][2]))
    return e2, det

def is_Mconvex3(S):
    Ss = set(S)
    for z in S:
        for zp in S:
            for p in range(3):
                if z[p] <= zp[p]: continue
                ok = False
                for q in range(3):
                    if z[q] >= zp[q]: continue
                    A = list(z); A[p] -= 1; A[q] += 1
                    B = list(zp); B[p] += 1; B[q] -= 1
                    if tuple(A) in Ss and tuple(B) in Ss: ok = True; break
                if not ok: return False
    return True

def N_is_lorentzian(f, d):
    """f: dict alpha->c (homogeneous degree d, 3 vars).  Returns (ok, reason)."""
    if d < 2: return (True, 'deg<2')
    if not is_Mconvex3([a for a, c in f.items() if c]): return (False, 'support not M-convex')
    for b0 in range(d-1):
        for b1 in range(d-1-b0):
            b2 = d-2-b0-b1
            beta = (b0, b1, b2)
            M = [[0]*3 for _ in range(3)]
            E = [(1,0,0),(0,1,0),(0,0,1)]
            for i in range(3):
                for j in range(3):
                    key = tuple(beta[t]+E[i][t]+E[j][t] for t in range(3))
                    M[i][j] = f.get(key, 0)
            if any(M[i][j] < 0 for i in range(3) for j in range(3)): return (False, 'neg entry')
            e2, det = e2_det(M)
            if e2 > 0 or det < 0: return (False, f'beta={beta} e2={e2} det={det}')
    return (True, 'ok')

if __name__ == '__main__':
    print("--- (i) N(P~_i) Lorentzian?  (one bead) ---", flush=True)
    st = Counter(); wit = []
    for h in range(0, 11):
        for l in range(0, h+1):
            f = Ptilde(l, h)
            ok, why = N_is_lorentzian(f, h)
            st['lorentzian' if ok else 'FAIL'] += 1
            if not ok and len(wit) < 4: wit.append((l, h, why))
    print(dict(st), flush=True)
    for w in wit: print("   FAIL", w, flush=True)

    print("--- (ii) N(prod_i P~_i) Lorentzian?  (m=2,3) ---", flush=True)
    st = Counter(); wit = []
    for m in (2, 3):
        rng = 5 if m == 2 else 3
        for h in product(range(0, rng+1), repeat=m):
            for l in product(*[range(0, h[i]+1) for i in range(m)]):
                f = {(0,0,0): 1}
                for i in range(m): f = mul(f, Ptilde(l[i], h[i]))
                f = {k: v for k, v in f.items() if v}
                if not f: continue
                ok, why = N_is_lorentzian(f, sum(h))
                st[f'm{m}_lorentzian' if ok else f'm{m}_FAIL'] += 1
                if not ok and len(wit) < 4: wit.append((m, l, h, why))
    print(dict(st), flush=True)
    for w in wit: print("   FAIL", w, flush=True)

    print("--- (iii) REFUSAL CONTROL: a NON-M-convex-support factor must break it ---", flush=True)
    st = Counter()
    for h in range(2, 8):
        f = Ptilde(0, h)
        # delete one interior support point -> support no longer M-convex
        for key in list(f):
            if key[0] >= 1 and key[1] >= 1 and key[2] >= 1:
                g = dict(f); del g[key]
                ok, why = N_is_lorentzian(g, h)
                st['REFUSED' if not ok else 'passed'] += 1
                break
    print(dict(st), flush=True)
