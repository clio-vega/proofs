"""Q417 STEP 2, CORRECTED (b), plus the BUDGET LEMMA test.

CORRECTION, recorded: step2.py's `wraps(P,w)` tested Q416's PER-STRIP seam predicate
(h in J_cyclic(theta), i.e. theta_h=0 & theta_1=1).  That is NOT the set PROVE.md sec 3.3
names.  Set (b) = the tableaux that are ORDINARY, i.e. whose ENDPOINT lambda satisfies
lambda_1 - lambda_h < w (the cylindric identification never bites).  At ell=0 that is the
2 SSYT of shape (2,2) -- the set the Q410 paper proved sums to 2-3t+t^3.
Two different predicates; one is on strips, one is on the endpoint.
"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1010-q410')
from korff import paths, weight_path, mvec, J_open, t
import sympy as sp

def cminus(la, k, w):
    m = 2*k
    if all(la[2*i] == la[2*i+1] for i in range(k)): return 1
    if la[0]-la[m-1] == w and all(la[2*i+1] == la[2*i+2] for i in range(k-1)): return -1
    return 0

def ordinary(la, w):
    """lambda does not reach the cylindric seam: lambda_1 - lambda_h < w."""
    return la[0] - la[-1] < w

def thetas(P):
    h = len(P[0])
    return [tuple(P[a+1][i]-P[a][i] for i in range(h)) for a in range(len(P)-1)]

def wd(th):
    return all(th[i] >= th[i+1] for i in range(len(th)-1))

def target(M):
    return sp.expand(sp.catalan(M) - (2*M-1)*t**(M-1) + t**(M+1))

def ord_t1(f):
    if sp.expand(f) == 0: return None
    u = sp.Symbol('u'); P = sp.Poly(sp.expand(sp.expand(f).subs(t, 1+u)), u)
    return min(m[0] for m, c in zip(P.monoms(), P.coeffs()) if c != 0)

def run(k, ell):
    h = 2*k; w = 2*ell+2; M = k+ell+1; n = 2*M
    alpha = (1,)*n
    allP = paths(h, w, alpha)
    sa = sp.Integer(0); sb = sp.Integer(0); A=B=0; nb=0; nwd=0
    for P in allP:
        c = cminus(P[-1], k, w)
        if all(wd(th) for th in thetas(P)): nwd += 1
        if c == 0: continue
        A += (c==1); B += (c==-1)
        W = weight_path(P, w, "Psi_open")
        sa += c*W
        if ordinary(P[-1], w):
            sb += c*W; nb += 1
    sa, sb, tg = sp.expand(sa), sp.expand(sb), target(M)
    print(f"\n===== k={k} ell={ell}  h={h} w={w} M={M} n=N={n} alpha=(1^{n}) =====")
    print(f"  tableaux {len(allP)};  A={A} B={B} A-B={A-B};  all-weakly-decreasing T: {nwd}")
    print(f"  TARGET  c_alpha = {tg}   c_alpha(1)={tg.subs(t,1)}  ord_(t=1)={ord_t1(tg)}")
    print(f"  (a) ALL TABLEAUX      sum_T c^- Psi_open = {sa}")
    print(f"      at t=1 = {sa.subs(t,1)}   ord_(t=1) = {ord_t1(sa)}")
    print(f"      (a) - target = {sp.expand(sa-tg)}   [zero? {sp.expand(sa-tg)==0}]")
    print(f"  (b) ORDINARY ONLY ({nb} T)  = {sp.factor(sb)}  = {sb}")
    print(f"      at t=1 = {sb.subs(t,1)}")
    print(f"      (b) - target = {sp.expand(sb-tg)}   [zero? {sp.expand(sb-tg)==0}]")

if __name__ == "__main__":
    for (k, ell) in [(1,0),(1,1),(1,2),(1,3)]:
        run(k, ell)
