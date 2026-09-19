import sympy as sp, itertools
from sixvertex import Z, schur_weights, check_ff

x = sp.symbols('x1:6'); eps = sp.Symbol('eps')
A, B, C = sp.symbols('A B C')

def gen_ff(xi, A, B, C):
    """General free-fermion row weights: a1=1,a2=A,b1=xi,b2=B,c1=C,c2=(A+xi*B)/C."""
    return {'a1':sp.Integer(1),'a2':A,'b1':xi,'b2':B,'c1':C,'c2':(A+xi*B)/C}

print("FF identity holds symbolically:", check_ff(gen_ff(x[0],A,B,C)))

print("\n=== (1) does k count ROWS or PARTICLES?  n=3 rows, p=2 particles ===")
n, m, p = 3, 6, 2
w0 = [schur_weights(x[i]) for i in range(n)]
wd = [gen_ff(x[i], 1, eps, 1) for i in range(n)]
top = [1 if j < p else 0 for j in range(m)]
for botset in itertools.combinations(range(m), p):
    bot = [1 if j in botset else 0 for j in range(m)]
    z0 = Z(n, m, top, bot, [0]*n, [0]*n, w0)
    if z0 == 0: continue
    ze = Z(n, m, top, bot, [0]*n, [0]*n, wd)
    found = None
    for k in range(0, 2*m+1):
        if sp.simplify(sp.expand(ze - z0*sp.prod([(1+eps*x[i])**k for i in range(n)]))) == 0:
            found = k; break
    print("   n=%d p=%d bot=%s -> k = %s" % (n, p, botset, found))

print("\n=== (2) GENERAL free-fermion deformation: does Z still factor? ===")
n, m = 2, 5
wg = [gen_ff(x[i], A, B, C) for i in range(n)]
top = [1,1,0,0,0]
for botset in [(2,3),(2,4),(3,4)]:
    bot = [1 if j in botset else 0 for j in range(m)]
    z0 = Z(n, m, top, bot, [0]*n, [0]*n, [schur_weights(x[i]) for i in range(n)])
    zg = sp.simplify(sp.expand(Z(n, m, top, bot, [0]*n, [0]*n, wg)))
    q  = sp.simplify(sp.cancel(zg / z0))
    print("   bot=%s : Z_gen / Z_0 = %s" % (botset, sp.factor(sp.simplify(q))))
