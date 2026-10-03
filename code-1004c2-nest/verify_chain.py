"""Verification of this session's proof chain.

(1) ONE-SIDED REDUCTION.  On a slice with zero horizontal-strip defect (all y_i>=0),
    eq:w is LINEAR: w_i = g_i+1-y_i, so W = {w} is a BOX SLICE (M-convex) and the shift
    sum_i(L_i-nu_i) = sum_i (y_i)_+ = sigma is CONSTANT.  So (A) there == (Q).

(2) (Q)  <=>  k(a) = #{(x,e) in Z_{>=0}^{2m} : l_i<=x_i+e_i<=h_i, sum x=a, sum e=D-a}
    log-concave, with supp k = [0,D] exactly (interval, free).

(3) CLAIM 1.  Omega = {(x,e)>=0 : l_i<=x_i+e_i<=h_i, sum x+sum e = D} is M-CONVEX.

(4) THE EXCHANGE MAP.  z in K_{a-1}, z' in K_{a+1}: pick p with x'_p>x_p and q with
    e'_q<e_q; then z+e_{x_p}-e_{e_q} and z'-e_{x_p}+e_{e_q} are BOTH in K_a.

(5) FREE-D VERSION is PF2 by pf2-convolution (no second hyperplane).
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen import is_pf2, conv_intervals
from itertools import product
from collections import Counter

def Omega(l, h, D):
    m = len(l); out = []
    for v in product(*[range(l[i], h[i]+1) for i in range(m)]):
        if sum(v) != D: continue
        for x in product(*[range(0, v[i]+1) for i in range(m)]):
            out.append((x, tuple(v[i]-x[i] for i in range(m))))
    return out

def is_Mconvex_pairs(S):
    """S: list of (x,e) in Z^{2m}; flatten to z in Z^{2m} and test symmetric exchange."""
    Z = [tuple(list(x)+list(e)) for (x, e) in S]; Zs = set(Z); N = len(Z[0])
    for z in Z:
        for zp in Z:
            for p in range(N):
                if z[p] <= zp[p]: continue
                ok = False
                for q in range(N):
                    if z[q] >= zp[q]: continue
                    A = list(z); A[p] -= 1; A[q] += 1
                    B = list(zp); B[p] += 1; B[q] -= 1
                    if tuple(A) in Zs and tuple(B) in Zs: ok = True; break
                if not ok: return (z, zp, p)
    return None

def Ka(l, h, D, a):
    return [(x, e) for (x, e) in Omega(l, h, D) if sum(x) == a]

def exchange_check(l, h, D, a):
    """(4): for EVERY (z,z') in K_{a-1} x K_{a+1} and the canonical p,q,
       are both images in K_a?  Also records degrees for the bipartite test."""
    m = len(l)
    Km, K0, Kp = Ka(l,h,D,a-1), Ka(l,h,D,a), Ka(l,h,D,a+1)
    K0s = set(K0)
    bad = 0; tot = 0; degP = []; degQ = Counter()
    for (x,e) in Km:
        for (xp,ep) in Kp:
            tot += 1
            # all valid transfers t = (+1 on x_p, -1 on e_q)
            good = []
            for p in range(m):
                if xp[p] <= x[p]: continue
                for q in range(m):
                    if ep[q] >= e[q]: continue
                    z3 = (tuple(x[i]+(1 if i==p else 0) for i in range(m)),
                          tuple(e[i]-(1 if i==q else 0) for i in range(m)))
                    z2 = (tuple(xp[i]-(1 if i==p else 0) for i in range(m)),
                          tuple(ep[i]+(1 if i==q else 0) for i in range(m)))
                    if z3 in K0s and z2 in K0s:
                        good.append((p,q)); degQ[(z3,z2)] += 1
            # the canonical choice: smallest p, smallest q
            P = [p for p in range(m) if xp[p] > x[p]]
            Q = [q for q in range(m) if ep[q] < e[q]]
            pq = None
            same = [i for i in P if i in Q]
            p0, q0 = (same[0], same[0]) if same else (P[0], Q[0])
            z3 = (tuple(x[i]+(1 if i==p0 else 0) for i in range(m)),
                  tuple(e[i]-(1 if i==q0 else 0) for i in range(m)))
            z2 = (tuple(xp[i]-(1 if i==p0 else 0) for i in range(m)),
                  tuple(ep[i]+(1 if i==q0 else 0) for i in range(m)))
            if not (z3 in K0s and z2 in K0s): bad += 1
            degP.append(len(good))
    return bad, tot, (min(degP) if degP else None), (max(degQ.values()) if degQ else None), len(Km), len(K0), len(Kp)

if __name__ == '__main__':
    print("--- (3) CLAIM 1: Omega M-convex ---", flush=True)
    st = Counter(); wit = []
    for m in (2,3,4):
        for l in product(range(0,3), repeat=m):
            for h in product(*[range(l[i], 4) for i in range(m)]):
                for D in range(sum(l), sum(h)+1):
                    S = Omega(list(l),list(h),D)
                    if len(S) < 2: continue
                    st['sets'] += 1; st['maxsize'] = max(st.get('maxsize',0), len(S))
                    if len(S) > 120: st['skipped_big'] += 1; continue
                    r = is_Mconvex_pairs(S)
                    st['Mconvex' if r is None else 'NOT_Mconvex'] += 1
                    if r is not None and len(wit) < 3: wit.append((m,l,h,D,r))
    print(dict(st), flush=True)
    for w in wit: print("   NOT M-convex:", w, flush=True)

    print("--- (3) negative control: drop the PAIR structure (l_i<=x_i<=h_i only) ---", flush=True)
    # replace pair sums by single-coordinate bounds on x only: should NOT be M-convex in general
    def Omega2(l,h,D):
        m=len(l); out=[]
        for x in product(*[range(l[i],h[i]+1) for i in range(m)]):
            for e in product(*[range(0,3) for i in range(m)]):
                if sum(x)+sum(e)==D: out.append((x,e))
        return out
    st=Counter()
    for m in (2,3):
        for l in product(range(0,2),repeat=m):
            for h in product(*[range(l[i],3) for i in range(m)]):
                for D in range(0,sum(h)+2*m+1):
                    S=Omega2(list(l),list(h),D)
                    if len(S)<2 or len(S)>120: continue
                    st['Mconvex' if is_Mconvex_pairs(S) is None else 'REFUSED']+=1
    print(dict(st), flush=True)

    print("--- (4) exchange map + bipartite degrees ---", flush=True)
    st = Counter(); deg_ok = 0; deg_bad = []
    for m in (2,3):
        for l in product(range(0,3), repeat=m):
            for h in product(*[range(l[i], 4) for i in range(m)]):
                for D in range(sum(l), sum(h)+1):
                    for a in range(1, D):
                        r = exchange_check(list(l),list(h),D,a)
                        bad, tot, mnP, mxQ, nm, n0, np_ = r
                        if tot == 0: continue
                        st['pairs'] += tot
                        st['canonical_ok' if bad==0 else 'canonical_BAD'] += 1
                        if mnP is not None and mxQ is not None:
                            if mnP >= mxQ: deg_ok += 1
                            else: deg_bad.append((m,l,h,D,a,mnP,mxQ,nm,n0,np_))
    print(dict(st), "degree_test_ok:", deg_ok, "degree_test_fail:", len(deg_bad), flush=True)
    for w in deg_bad[:4]: print("   mindeg(P) < maxdeg(Q):", w, flush=True)

    print("--- (5) FREE-D version: product form, PF2 by pf2-convolution ---", flush=True)
    st=Counter()
    for m in (2,3,4,5):
        for l in product(range(0,3),repeat=m):
            for h in product(*[range(l[i],5) for i in range(m)]):
                # sum over ALL D: coefficient sequence = prod_i (sum_{v=l_i}^{h_i} [v+1]_z)
                cur=[1]
                for i in range(m):
                    fac=[0]*(h[i]+1)
                    for v in range(l[i],h[i]+1):
                        for aa in range(v+1): fac[aa]+=1
                    new=[0]*(len(cur)+len(fac)-1)
                    for j,c in enumerate(cur):
                        if c:
                            for kk,d in enumerate(fac): new[j+kk]+=c*d
                    cur=new
                st['PF2' if is_pf2(cur) else 'FAIL']+=1
                # and each factor PF2?
                for i in range(m):
                    fac=[0]*(h[i]+1)
                    for v in range(l[i],h[i]+1):
                        for aa in range(v+1): fac[aa]+=1
                    st['factor_PF2' if is_pf2(fac) else 'factor_FAIL']+=1
    print(dict(st), flush=True)
