"""Does the q-graded wrap (Korff 1906.02565 eq (cylchi): q^d chi^{la/d/mu}_t) rescue
the ribbon weight on the fibre alpha=(1^n)?

d is read from Korff 1110.6356 proof of Thm 5.6: the row configuration is
sigma=(sigma^1,...,sigma^n,sigma^1) with "sigma^1 := sigma^{n+1} := d" and
sigma^i = la[d]'_i - mu[0]'_i = theta_i.  Hence d_a = theta^(a)_1 for each step,
and the total degree is d(T) = sum_a theta^(a)_1 = x^(final)_1 (since x^(0)=0).
Depends only on the ENDPOINT -- consistent with d being part of the shape la/d/mu.
"""
import sys, itertools
sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from core import Engine
from warn import G_times_prodx
from korff import paths
t, q = sp.symbols('t q')

def c_alpha_M1(k, ell, alpha):
    n=len(alpha); E=Engine(n); F=G_times_prodx(E,k,ell)
    mon=sp.prod([E.x[i]**alpha[i] for i in range(n)])
    return sp.expand(sp.Poly(F,*E.x).coeff_monomial(mon) or sp.Integer(0))

def sign_hkko(la,k,w):
    h=2*k
    if all(la[2*i]==la[2*i+1] for i in range(k)): return 1
    if la[0]-la[h-1]==w and all(la[2*i+1]==la[2*i+2] for i in range(k-1)): return -1
    return 0

for (k,ell,n) in [(2,1,4),(3,1,4),(2,1,2),(2,2,4)]:
    h,w=2*k,2*ell+2; al=tuple([1]*n)
    P=paths(h,w,al)
    tot=sp.Integer(0); tbl={}
    for p in P:
        la=p[-1]; s=sign_hkko(la,k,w); d=la[0]
        tot += s*q**d
        tbl[(la,s,d)] = tbl.get((la,s,d),0)+1
    c=c_alpha_M1(k,ell,al)
    print(f"k={k} ell={ell} n={n} (h={h},w={w}): #T={len(P)}")
    for (la,s,d),mult in sorted(tbl.items()):
        print(f"    endpoint {la}  c^-={s:+d}  d=theta_1-sum={d}  (x{mult} paths)")
    print(f"    sum_T c^-(T) q^d(T)  = {sp.expand(tot)}        [ribbon wt is 1 on single boxes]")
    print(f"    target c_alpha(t)    = {sp.expand(c)}")
    sols = sp.solve(sp.Eq(sp.expand(tot), sp.expand(c)), q)
    print(f"    solve for q = f(t):  {sols}")
    print()
