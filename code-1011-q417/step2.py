"""Q417 STEP 2 (PROVE.md sec 3.2-3.3): compute BOTH sums, label which is which.

(a) Psi_open over ALL cylindric tableaux of the fibre  <-- THIS IS Q417
(b) Psi_open over the NON-WRAPPING tableaux only       <-- a PROVED theorem, ell=0 anchor

Report exact polynomials and the DIFFERENCE, never a boolean (PROVE.md sec 4).
"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1010-q410')
sys.path.insert(0, '/home/clio/projects/proofs/code-1011-q417')
from korff import paths, weight_path, mvec, t
import sympy as sp

def cminus(la, k, w):
    m = 2*k
    if all(la[2*i] == la[2*i+1] for i in range(k)): return 1
    if la[0]-la[m-1] == w and all(la[2*i+1] == la[2*i+2] for i in range(k-1)): return -1
    return 0

def wraps(P, w):
    """T wraps iff some strip's wrap multiplicity m_h is ACTIVE, i.e. the wrap factor
    (1 - t^{w-(x_1-x_h)}) differs from the no-wrap reading.  Operationally (Q416):
    the strip wraps iff theta_h=0, theta_1=1 -- the cyclic descent at the seam."""
    h = len(P[0])
    for a in range(len(P)-1):
        th = tuple(P[a+1][i]-P[a][i] for i in range(h))
        if th[h-1] == 0 and th[0] == 1:
            return True
    return False

def target(M):
    return sp.expand(sp.catalan(M) - (2*M-1)*t**(M-1) + t**(M+1))

def run(k, ell, kind="Psi_open"):
    h = 2*k; w = 2*ell+2; M = k+ell+1; n = 2*M
    alpha = (1,)*n
    allP = paths(h, w, alpha)
    sa = sp.Integer(0); sb = sp.Integer(0); A = 0; B = 0; nb = 0
    for P in allP:
        c = cminus(P[-1], k, w)
        if c == 0: continue
        if c == 1: A += 1
        else: B += 1
        W = weight_path(P, w, kind)
        sa += c*W
        if not wraps(P, w):
            sb += c*W; nb += 1
    sa = sp.expand(sa); sb = sp.expand(sb); tg = target(M)
    print(f"\n===== k={k} ell={ell} (h={h}, w={w}, M={M}, n=N={n}, alpha=(1^{n})) =====")
    print(f"  |all tableaux|={len(allP)}   A_alpha={A}  B_alpha={B}  (A-B={A-B})")
    print(f"  TARGET      c_alpha(t)                = {tg}")
    print(f"              c_alpha(1)                = {tg.subs(t,1)}    "
          f"ord_(t=1) = {sp.Poly(tg,t).as_expr().subs(t,1+sp.Symbol('u')) and None or ''}", end="")
    P1 = sp.Poly(sp.expand(tg.subs(t, 1+sp.Symbol('u'))), sp.Symbol('u'))
    o = min([mo[0] for mo, co in zip(P1.monoms(), P1.coeffs()) if co != 0])
    print(f"{o}")
    print(f"  (a) SUM OVER ALL TABLEAUX   sum_T c^- {kind}_T = {sa}")
    print(f"      at t=1 : {sa.subs(t,1)}")
    print(f"      DIFFERENCE (a) - target = {sp.expand(sa - tg)}")
    print(f"  (b) SUM OVER NON-WRAPPING ONLY ({nb} tableaux) = {sb}")
    print(f"      factored: {sp.factor(sb)}      at t=1 : {sb.subs(t,1)}")
    print(f"      DIFFERENCE (b) - target = {sp.expand(sb - tg)}")
    return sa, sb, tg

if __name__ == "__main__":
    for (k, ell) in [(1,0),(1,1),(1,2),(1,3)]:
        run(k, ell)
