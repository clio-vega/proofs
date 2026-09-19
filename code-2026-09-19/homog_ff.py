"""Homogeneous free-fermion six-vertex models and the Lorentzian question (Q177).

For N(.) and "Lorentzian" to be well posed, Z must be HOMOGENEOUS with nonnegative
coefficients.  Give row i two variables (x_i, y_i) and let all six weights be linear
forms with nonnegative coefficients in (x_i, y_i).  Then Z is homogeneous of degree n*m.

RIGIDITY OBSERVATION (standard convention). If a1,a2,c1,c2 are multiples of y and
b1,b2 are multiples of x, then a1a2+b1b2 = alpha*y^2 + beta*x^2 while c1c2 = gamma*y^2,
so the free-fermion condition forces beta = b1b2 = 0: the model degenerates to FIVE
vertices, Z is a Schur polynomial, and Lorentzian-ness is exactly HMMS.

So a genuine homogeneous nonneg free-fermion six-vertex model must MIX the variables.
The family below does:  a1 = A x, a2 = y, b1 = B x, b2 = y, c1 = (A+B) x, c2 = y,
for which a1a2 + b1b2 = (A+B) x y = c1 c2 exactly, with every weight nonzero.
"""
import sympy as sp, itertools
from fractions import Fraction
from sixvertex import Z, SIX
from lorentzian import is_lorentzian

def hom_ff(xi, yi, A, B):
    return {'a1':A*xi, 'a2':yi, 'b1':B*xi, 'b2':yi, 'c1':(A+B)*xi, 'c2':yi}

def ff_ok(w):
    return sp.simplify(w['a1']*w['a2'] + w['b1']*w['b2'] - w['c1']*w['c2']) == 0

def to_coeff_dict(expr, vars_):
    p = sp.Poly(sp.expand(expr), *vars_)
    return {tuple(m): Fraction(int(c.p), int(c.q)) if hasattr(c,'p') else Fraction(int(c))
            for m, c in zip(p.monoms(), p.coeffs())}

def normalize(cd):
    from lorentzian import multi_factorial
    return {a: c / multi_factorial(a) for a, c in cd.items()}

for (A,B) in [(1,1),(1,2),(2,1),(1,3),(3,2)]:
    print("=== A=%s B=%s ===" % (A,B), flush=True)
    for n,m in [(2,2),(2,3),(3,2)]:
        xs = sp.symbols('x1:%d'%(n+1)); ys = sp.symbols('y1:%d'%(n+1))
        w = [hom_ff(xs[i], ys[i], sp.Integer(A), sp.Integer(B)) for i in range(n)]
        assert all(ff_ok(wi) for wi in w), "free-fermion condition FAILED"
        vars_ = list(xs)+list(ys)
        for top in itertools.product([0,1], repeat=m):
            for bot in itertools.product([0,1], repeat=m):
                z = Z(n, m, list(top), list(bot), [0]*n, [0]*n, w)
                if z == 0: continue
                cd = to_coeff_dict(z, vars_)
                d = {sum(a) for a in cd}
                ok, why = is_lorentzian(normalize(cd), len(vars_))
                if not ok:
                    print("   n=%d m=%d top=%s bot=%s deg=%s : NOT LORENTZIAN -- %s"
                          % (n,m,top,bot,d,why), flush=True)
        print("   n=%d m=%d done" % (n,m), flush=True)
