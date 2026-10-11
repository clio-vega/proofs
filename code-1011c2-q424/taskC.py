"""TASK C (Q418): the pre-registered screen on Korff 1906.02565's ribbon weight.

Weight read verbatim from 1906.02565-src/draft_02.tex eq (Hecke_wt), line 466:
    wt(b) = (t-1)^{#b-1} prod_{h in b} (-1)^{r(h)-1} t^{c(h)-1}
cylindric extension, line 1258, same formula with #b = number of DISTINCT cylindric
rim hooks and r,c counted in rows 1<=i<=k of lambda[d].

Pre-registered kill was: "if it is again prod_{i in J}(1-t^{m_i}) over a cyclic index
set, Q418 dies by Q410's own theorem."  IT IS NOT.  So instead we run the screen.

SCREEN (node weight-must-be-one-at-t1): need W_T(1) != 0 for some T using both unit
steps.  For a SINGLE BOX b: |b|=1, #b=1, r(h)=c(h)=1, so
    wt(box) = (t-1)^0 * (-1)^0 * t^0 = 1    exactly, for every t.
=> on any fibre alpha=(1^n) the ribbon weight is IDENTICALLY 1, W_T(t) = 1 for all T.
=> sum_T c^-(T) W_T(t) = sum_T c^-(T) = c_alpha(1), a CONSTANT in t.
So the ribbon weight passes the t=1 screen but can only match c_alpha(t) where
c_alpha(t) is constant.  NEW TEST: is c_{(1^n)}(t) constant throughout Korff's domain?
"""
import sys, itertools
sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from core import Engine
from warn import G_times_prodx
from walks import AB_from_walks
from korff import paths, J_cyclic, weight_path
t = sp.Symbol('t')

def c_alpha_M1(k, ell, alpha):
    n=len(alpha); E=Engine(n); F=G_times_prodx(E,k,ell)
    mon=sp.prod([E.x[i]**alpha[i] for i in range(n)])
    return sp.expand(sp.Poly(F,*E.x).coeff_monomial(mon) or sp.Integer(0))

print("wt(single box) = (t-1)^0 * (-1)^(1-1) * t^(1-1) =", sp.expand((t-1)**0*(-1)**0*t**0))
print()
print("TEST: c_{(1^n)}(t) for k>=2 (so n_K=2k>2) and n<=w=2ell+2 (so Korff's P is not 0)")
print(f"{'k':>2}{'ell':>4}{'h':>3}{'w':>3}{'n':>3}{'#T':>5}  {'c_(1^n)(t)  [Warnaar engine]':<34}{'const?':>7}{'c(1)':>6}{'A-B':>5}")
rows=[]
for k in (2,3):
    for ell in range(1,4):
        w=2*ell+2; h=2*k
        for n in range(2, min(w,5)+1):
            al=tuple([1]*n)
            nT=len(paths(h,w,al))
            c=c_alpha_M1(k,ell,al)
            A,B,TOT=AB_from_walks(k,w,al)
            const = (sp.expand(c)-sp.expand(c).subs(t,1))==0
            rows.append((k,ell,n,const,sp.expand(c)))
            print(f"{k:>2}{ell:>4}{h:>3}{w:>3}{n:>3}{nT:>5}  {str(sp.expand(c)):<34}{str(const):>7}{str(sp.expand(c).subs(t,1)):>6}{A-B:>5}")
nonconst=[r for r in rows if not r[3]]
print()
print(f"NON-CONSTANT in-domain points: {len(nonconst)} of {len(rows)}")
for r in nonconst: print(f"   k={r[0]} ell={r[1]} n={r[2]}: c = {r[4]}")
print()
print("OUT-OF-DOMAIN CONTROL (k=1, ell=0, n=4; n>w so Korff sets P=0):")
c=c_alpha_M1(1,0,(1,1,1,1)); print(f"   c_(1^4) = {sp.expand(c)}  constant? {(sp.expand(c)-sp.expand(c).subs(t,1))==0}")
