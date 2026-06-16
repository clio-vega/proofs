import sympy as sp
from math import comb, factorial
from mn import Mj
a,b = sp.symbols('a b')

def collect(c,i,amax=70,bmax=30):
    pts=[]
    for bb in range(c,bmax):
        for aa in range(bb,amax):
            if (aa+bb+c)%2!=0: continue
            m=(aa+bb+c)//2
            j=bb+i
            if j>m: continue
            pts.append((aa,bb,Mj((aa,bb,c),j,m)))
    return pts

def fit(pts,maxdeg):
    monos=[(p,q) for p in range(maxdeg+1) for q in range(maxdeg+1) if p+q<=maxdeg]
    if len(pts)<len(monos)+3: return None
    A=sp.Matrix([[aa**p*bb**q for (p,q) in monos] for (aa,bb,_) in pts])
    y=sp.Matrix([M for (_,_,M) in pts])
    sol=A.solve_least_squares(y)
    if any(A*sol-y): return None
    expr=sum(sol[k]*a**p*b**q for k,(p,q) in enumerate(monos))
    return sp.factor(sp.expand(expr))

c=4
for i in [4,3,2,1]:
    f=fit(collect(c,i), maxdeg=2*c+2)
    print("M_{b+%d} (c=4) ="%i)
    print("   ", f)
    print("   factor_list:", sp.factor_list(f)[1] if f is not None else None)
    print()
