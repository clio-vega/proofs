"""Q417 STEP 6: the EXACT leading behaviour at k=1, as a closed form (proved, then checked).

CLAIM.  At k=1, alpha=(1^{2M}), w=2M-2: the cylindric tableaux with EXACTLY ONE
J_open-nonempty strip are exactly the 2M-2 paths
   0 -> (1,0) -> ... -> (p,0) -> (p,1) -> (p+1,1) -> ... -> (2M-1,1),   p = 1..2M-2,
each with c^-(T) = -1 and Psi_open_T = 1 - t^p.  Every other T has >= 2 bad strips.
Hence
   sum_T c^-(T) Psi_open_T = -sum_{p=1}^{2M-2} (1 - t^p) + O((t-1)^2),
so  d/dt |_{t=1}  =  sum_{p=1}^{2M-2} p  =  (M-1)(2M-1),
and ord_(t=1) = 1 EXACTLY (not merely >= 1), since (M-1)(2M-1) > 0 for M >= 2.

Also tested: is {T : exactly one bad strip} = {T : c^- = -1} as SETS, not just as
cardinalities?  (memory: two-sets-of-the-same-size-are-not-the-same-set)
"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1010-q410')
from korff import paths, J_open, weight_path, t
import sympy as sp

def thetas(P):
    h = len(P[0])
    return [tuple(P[a+1][i]-P[a][i] for i in range(h)) for a in range(len(P)-1)]
def cminus(la, k, w):
    m = 2*k
    if all(la[2*i] == la[2*i+1] for i in range(k)): return 1
    if la[0]-la[m-1] == w and all(la[2*i+1] == la[2*i+2] for i in range(k-1)): return -1
    return 0

for ell in [0,1,2,3,4]:
    k = 1; h = 2; w = 2*ell+2; M = k+ell+1; n = 2*M
    one = []; S = sp.Integer(0); neg = []
    for P in paths(h, w, (1,)*n):
        b = sum(1 for th in thetas(P) if J_open(th))
        c = cminus(P[-1], k, w)
        if b == 1: one.append((P, c, weight_path(P, w, "Psi_open")))
        if c == -1: neg.append(P)
        if c: S += c*weight_path(P, w, "Psi_open")
    S = sp.expand(S)
    pset = sorted(sp.Poly(wt, t).degree() for _, _, wt in one)
    lead = sp.diff(S, t).subs(t, 1)
    pred = (M-1)*(2*M-1)
    setseq = set(P for P,_,_ in one) == set(neg)
    print(f"ell={ell} M={M} n={n} w={w}:")
    print(f"   #T with exactly one bad strip = {len(one)}   (predicted 2M-2 = {2*M-2})")
    print(f"   their c^- values = {sorted(set(c for _,c,_ in one))}   (predicted [-1])")
    print(f"   their Psi_open = 1 - t^p with p = {pset}   (predicted 1..{2*M-2})")
    print(f"   {{exactly one bad strip}} == {{c^- = -1}} as SETS? {setseq}"
          f"   (|B| = {len(neg)})")
    print(f"   d/dt sum_T c^- Psi_open at t=1 = {lead}   predicted (M-1)(2M-1) = {pred}"
          f"   {'MATCH' if lead == pred else 'MISMATCH'}")
    print(f"   ord_(t=1) of the sum = 1 exactly? {sp.expand(S).subs(t,1)==0 and lead!=0}")
