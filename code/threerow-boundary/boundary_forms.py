from mn import Mj
from fractions import Fraction
from math import comb
import itertools

# For each c, each boundary offset i (j=b+i, i=1..c), fit M_{b+i} as polynomial in (a,b).
# Collect data points (a,b) -> M_{b+i}, then solve for polynomial coeffs up to total degree D.

def collect(c, i, amax=40):
    pts=[]
    for b in range(c, 20):
        for a in range(b, amax):
            if (a+b+c)%2!=0: continue
            m=(a+b+c)//2
            j=b+i
            if j>m: continue
            M=Mj((a,b,c),j,m)
            pts.append((a,b,M))
    return pts

def fit_poly(pts, maxdeg=4):
    # try to fit M = P(a,b)/known? Just fit integer poly of total degree maxdeg in a,b.
    # monomials a^p b^q, p+q<=maxdeg
    monos=[(p,q) for p in range(maxdeg+1) for q in range(maxdeg+1) if p+q<=maxdeg]
    import sympy as sp
    A=sp.Matrix([[a**p*b**q for (p,q) in monos] for (a,b,M) in pts])
    y=sp.Matrix([M for (a,b,M) in pts])
    # least squares exact via solve (overdetermined, need consistent)
    # use first len(monos) independent rows then verify
    sol=A.solve_least_squares(y)
    # check residual
    res=A*sol-y
    ok=all(r==0 for r in res)
    return monos, sol, ok

import sympy as sp
a,b=sp.symbols('a b')
for c in [1,2,3,4,5]:
    for i in range(1,c+1):
        pts=collect(c,i)
        if len(pts)<15: 
            print(f"c={c} i={i}: too few pts ({len(pts)})"); continue
        for D in range(0,6):
            monos,sol,ok=fit_poly(pts,D)
            if ok:
                expr=sum(sol[k]*a**p*b**q for k,(p,q) in enumerate(monos))
                expr=sp.nsimplify(sp.expand(expr))
                print(f"c={c} i={i} (j=b+{i}): M = {sp.factor(expr)}   [deg {D}, {len(pts)} pts]")
                break
        else:
            print(f"c={c} i={i}: no integer poly fit up to deg 5")
