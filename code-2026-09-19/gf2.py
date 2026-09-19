import sympy as sp, itertools
from sixvertex import Z, schur_weights
x = sp.symbols('x1:6'); eps = sp.Symbol('eps')
def ffw(xi,B):
    return {'a1':sp.Integer(1),'a2':sp.Integer(1),'b1':xi,'b2':B,'c1':sp.Integer(1),'c2':1+xi*B}
n,m,p = 3,6,2
wd=[ffw(x[i],eps) for i in range(n)]
top=[1 if j<p else 0 for j in range(m)]
for botset in [(3,4),(3,5)]:
    bot=[1 if j in botset else 0 for j in range(m)]
    ze=sp.expand(Z(n,m,top,bot,[0]*n,[0]*n,wd))
    print("bot=%s"%(botset,), flush=True)
    print("   Z_eps =", sp.factor(ze), flush=True)
    poly = sp.Poly(ze, eps)
    for d,c in zip(range(poly.degree(),-1,-1), poly.all_coeffs()):
        if c!=0: print("     eps^%d : %s"%(d, sp.factor(c)), flush=True)
    print()
