"""Q430: does a Gaussian-binomial single-box weight reproduce c_{(1^4)} = 3 - t^2?

Structure of this script:
  PART 0  the fibre, printed in full: paths, endpoints, c^-, and the MULTIPLICITY
          READING SEQUENCE of each path (the only thing a slot-blind weight sees)
  PART 1  CONTROLS 1,2 as polynomials (omega = 1 -> 2 ; omega = 1-t^m -> known)
  PART 2  the structural reduction, verified with FREE SYMBOLS f1..f4
  PART 3  the candidate, and the RESIDUAL
  PART 4  CONTROL 3 (constant fibre k=2,ell=2,n=4) printed as integers
  PART 5  CONTROL 4, the arm I expect to be redundant
  last line prints a SENTINEL
"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from core import Engine
from warn import G_times_prodx
from korff import paths, J_cyclic, mvec

t = sp.Symbol('t')

def c_alpha_M1(k, ell, alpha):
    """Warnaar determinant engine (M1): the TARGET.  Independent of HKKO."""
    n = len(alpha); E = Engine(n); F = G_times_prodx(E, k, ell)
    return sp.expand(sp.Poly(F, *E.x).coeff_monomial(
        sp.prod([E.x[i]**alpha[i] for i in range(n)])) or sp.Integer(0))

def sign_hkko(la, k, w):
    """HKKO Thm 3.3 c^-_{2h,w}.  NOTE: implementation from a reading whose source
    I do not hold (see holdings row); audited only by CONTROL 2."""
    h = 2*k
    if all(la[2*i] == la[2*i+1] for i in range(k)): return 1
    if la[0]-la[h-1] == w and all(la[2*i+1] == la[2*i+2] for i in range(k-1)): return -1
    return 0

def readings(P, w):
    """The multiplicity read at each step of path P.  J_cyclic(e_i) is a SINGLETON,
    so each single-box step reads exactly one multiplicity.  Returns (slot, m) pairs."""
    out = []
    for a in range(len(P)-1):
        th = tuple(P[a+1][i]-P[a][i] for i in range(len(P[a])))
        J = J_cyclic(th)
        assert len(J) == 1, (th, J)            # single box => singleton index set
        i = J[0]; m = mvec(P[a], w)[i]
        assert m >= 1, (P[a], th, m)           # region condition forbids reading m=0
        out.append((i, m))
    return out

def sigma(k, ell, n, f, slot_aware=None):
    """sum_T c^-(T) W_T  with W_T = prod over steps of f(m)  (or g(slot,m))."""
    h, w = 2*k, 2*ell+2
    tot = sp.Integer(0)
    for P in paths(h, w, tuple([1]*n)):
        c = sign_hkko(P[-1], k, w)
        if c == 0: continue
        W = sp.Integer(1)
        for (i, m) in readings(P, w):
            W *= (slot_aware(i, m) if slot_aware else f(m))
        tot += c*W
    return sp.expand(tot)

# ----------------------------------------------------------------- PART 0
k, ell, n = 2, 1, 4
h, w = 2*k, 2*ell+2
print("="*78)
print(f"PART 0 -- THE FIBRE  k={k} ell={ell} n={n}  =>  h=2k={h}, w=2ell+2={w}, alpha=(1^{n})")
print(f"   in Korff's domain?  h={h} > 2 : {h>2} ;  n={n} <= w={w} : {n<=w}")
P_all = paths(h, w, tuple([1]*n))
print(f"   #cylindric tableaux in the fibre = {len(P_all)}")
target = c_alpha_M1(k, ell, tuple([1]*n))
print(f"   TARGET  c_alpha(t) = {sp.expand(target)}   [Warnaar determinant engine]")
print(f"           c_alpha(1) = {sp.expand(target).subs(t,1)} ,  deg_t c_alpha = {sp.Poly(target,t).degree()}")
print(f"           squarefree over QQ? {sp.factor_list(sp.expand(target))}")
print()
print("   path (slots of the boxes, in order) | endpoint | c^- | readings (slot,m) | m-sequence")
rows = []
for P in P_all:
    rd = readings(P, w)
    slots = tuple(j[0]+1 for j in rd)          # 1-indexed slot of each box
    c = sign_hkko(P[-1], k, w)
    ms = tuple(j[1] for j in rd)
    rows.append((slots, P[-1], c, rd, ms))
    print(f"   {str(slots):<16} | {str(P[-1]):<14} | {c:+d}  | {str(rd):<34} | {ms}")
print()
nz = [r for r in rows if r[2] != 0]
print(f"   tableaux with c^- != 0 : {len(nz)} of {len(rows)}")
print(f"   sum of c^- over the fibre = {sum(r[2] for r in rows)}  (= c_alpha(1)? {sum(r[2] for r in rows) == sp.expand(target).subs(t,1)})")
print()
print("   LEMMA A (forced first step): every path's first step is at slot 1 and reads m_h(0)=w.")
first = set((r[3][0][0]+1, r[3][0][1]) for r in rows)
print(f"      set of (slot,m) over ALL {len(rows)} first steps = {first}   (singleton? {len(first)==1})")
print()
print("   CANCELLATION: m-sequences of the c^-!=0 tableaux, grouped:")
from collections import defaultdict
g = defaultdict(list)
for r in nz: g[r[4]].append((r[0], r[1], r[2]))
for ms, lst in sorted(g.items()):
    print(f"      m-seq {ms}: net c^- = {sum(x[2] for x in lst):+d}  from {[(x[0],x[2]) for x in lst]}")
