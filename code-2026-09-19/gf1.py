import sympy as sp, itertools
from sixvertex import Z, schur_weights
x = sp.symbols('x1:6'); eps = sp.Symbol('eps')
def gen_ff(xi,A,B,C):
    return {'a1':sp.Integer(1),'a2':A,'b1':xi,'b2':B,'c1':C,'c2':sp.together((A+xi*B)/C)}
print("=== k counts ROWS or PARTICLES?  n=3 rows, p=2 particles ===", flush=True)
n,m,p = 3,6,2
w0=[schur_weights(x[i]) for i in range(n)]
wd=[gen_ff(x[i],1,eps,1) for i in range(n)]
top=[1 if j<p else 0 for j in range(m)]
for botset in itertools.combinations(range(m),p):
    bot=[1 if j in botset else 0 for j in range(m)]
    z0=Z(n,m,top,bot,[0]*n,[0]*n,w0)
    if z0==0: continue
    ze=Z(n,m,top,bot,[0]*n,[0]*n,wd)
    found=None
    for k in range(0,2*m+1):
        if sp.expand(ze - z0*sp.prod([(1+eps*x[i])**k for i in range(n)]))==0:
            found=k; break
    print("   n=%d p=%d bot=%s  Z_0=%s  -> k = %s"%(n,p,botset,z0,found), flush=True)
