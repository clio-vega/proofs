"""Q430 PART 7: Theorem B as a FAMILY test.  For each fibre alpha=(1^n), the
Gaussian-binomial weight forces  [w+1]_t | c_alpha(t)  in QQ[t], and in particular
(w+1) | c_alpha(1) in ZZ.  Test it across a grid."""
import sys, time
sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from core import Engine
from warn import G_times_prodx
from korff import paths
t=sp.Symbol('t')
def c_alpha_M1(k,ell,alpha):
    n=len(alpha); E=Engine(n); F=G_times_prodx(E,k,ell)
    return sp.expand(sp.Poly(F,*E.x).coeff_monomial(sp.prod([E.x[i]**alpha[i] for i in range(n)])) or sp.Integer(0))
def sign_hkko(la,k,w):
    h=2*k
    if all(la[2*i]==la[2*i+1] for i in range(k)): return 1
    if la[0]-la[h-1]==w and all(la[2*i+1]==la[2*i+2] for i in range(k-1)): return -1
    return 0
def gb(m): return sp.expand(sum(t**j for j in range(m+1)))

print("PART 7 -- the family test.  domain: h=2k>2 and n<=w=2ell+2.  k=1 shown too (out of domain).")
print(f"{'k':>2}{'ell':>4}{'n':>3} | {'h':>2}{'w':>3} | in-dom | {'c_alpha(t)':<34} | deg | c(1) | w+1 | (w+1)|c(1) | [w+1]_t | c_alpha")
print("-"*126)
rows=[]
for k in (1,2,3):
    for ell in (0,1,2):
        w=2*ell+2; h=2*k
        for n in range(2, min(w,5)+1):
            t0=time.time()
            try:
                c=sp.expand(c_alpha_M1(k,ell,tuple([1]*n)))
            except Exception as ex:
                print(f"{k:>2}{ell:>4}{n:>3} | ERROR {type(ex).__name__}: {ex}"); continue
            indom = (h>2 and n<=w)
            if c==0:
                deg=-1; c1=0
            else:
                deg=sp.Poly(c,t).degree(); c1=int(c.subs(t,1))
            g=gb(w)
            if c==0: divpoly=True
            else:
                _,r=sp.div(sp.Poly(c,t), sp.Poly(g,t)); divpoly=(r.as_expr()==0)
            divint = (c1 % (w+1) == 0)
            rows.append((k,ell,n,h,w,indom,c,deg,c1,divpoly,divint))
            print(f"{k:>2}{ell:>4}{n:>3} | {h:>2}{w:>3} | {str(indom):>6} | {str(c):<34} | {deg:>3} | {c1:>4} | {w+1:>3} | {str(divint):>10}  | {str(divpoly):>7}   [{time.time()-t0:.1f}s]")
print("-"*126)
nk = [r for r in rows if not r[9]]
print(f"fibres tested: {len(rows)} ; fibres where [w+1]_t does NOT divide c_alpha(t) (=> Gaussian weight DEAD): {len(nk)}")
print(f"   of those, in Korff's domain: {len([r for r in nk if r[5]])}")
print(f"fibres where it DOES divide (Gaussian not killed by Theorem B): {len([r for r in rows if r[9]])}")
for r in rows:
    if r[9]: print(f"   SURVIVES Thm B: k={r[0]} ell={r[1]} n={r[2]} c_alpha={r[6]}")
print(f"fibres where (w+1) does NOT divide c_alpha(1) (cheap integer kill): {len([r for r in rows if not r[10]])} of {len(rows)}")
print()
print("SENTINEL: q430c.py part 7 complete.")
