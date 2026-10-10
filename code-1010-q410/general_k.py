import sys, itertools
sys.path.insert(0,'.')
from korff import *
from walks import AB_from_walks
from core import Engine
from warn import G_times_prodx
import sympy as sp
t = sp.Symbol('t')
def ordt1(f):
    f = sp.expand(f)
    if f == 0: return None
    k = 0
    while True:
        if sp.simplify(sp.diff(f,t,k).subs(t,1)) != 0: return k
        k += 1
        if k > 30: return ">30"
print("alpha = (1^n), general k: A_alpha - B_alpha = c_alpha(1); kill is immediate when != 0")
print(f"{'k':>2} {'ell':>3} {'n':>2} {'h':>2} {'w':>2} {'#T':>5} {'A':>5} {'B':>5} {'c(1)=A-B':>9} {'nu(T)':>6} {'need ord>=':>10}")
for (k,ell,n) in [(1,0,4),(1,1,6),(1,2,8),(2,0,4),(2,0,5),(2,0,6),(2,1,6),(2,1,7),(3,0,6),(3,1,8)]:
    h=2*k; w=2*ell+2; al=tuple([1]*n)
    A,B,TOT = AB_from_walks(k,w,al)
    P=paths(h,w,al)
    nus={sum(len(J_cyclic(tuple(p[a+1][i]-p[a][i] for i in range(h)))) for a in range(len(p)-1)) for p in P}
    print(f"{k:>2} {ell:>3} {n:>2} {h:>2} {w:>2} {TOT:>5} {A:>5} {B:>5} {A-B:>9} {str(sorted(nus)):>6} {n:>10}")
print()
print("exact c_alpha(t) at alpha=(1^n) for k=2 (Warnaar determinant engine), with ord_{t=1}:")
for (k,ell,n) in [(2,0,4),(2,0,5),(2,1,6)]:
    E=Engine(n); F=G_times_prodx(E,k,ell)
    ca=sp.Poly(F,*E.x).coeff_monomial(sp.prod(E.x))
    print(f"   k={k} ell={ell} n={n}: c_alpha = {sp.expand(ca)} ; c(1)={sp.expand(ca).subs(t,1)} ;"
          f" ord_{{t=1}} = {ordt1(ca)} ; required >= {n}")
