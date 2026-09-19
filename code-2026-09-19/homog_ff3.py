"""WHY the first family was vacuous, and a family that is not.

In a row with left=right=0 the vertex counts satisfy
    n_{c1} = n_{c2} =: c,   n_{a2}+n_{b2} = p - c,   n_{a1}+n_{b1} = m - p - c,
where p = #particles (conserved).  In the family a1=Ax, a2=y, b1=Bx, b2=y,
c1=(A+B)x, c2=y the x-degree is n_{a1}+n_{b1}+n_{c1} = (m-p-c)+c = m-p: CONSTANT.
So Z is forced to be a monomial and "Lorentzian" is vacuous -- a kernel.

To break the kernel, a1 and b1 (or a2 and b2) must carry DIFFERENT variables.
Free-fermion + nonnegative + homogeneous is still satisfiable:
    a1 = alpha*x, a2 = alpha'*y, b1 = beta*y, b2 = beta'*x,
    c1*c2 = (alpha*alpha' + beta*beta') * x*y.
Now x-degree = n_{a1} + n_{b2} + (x-content of the c's), which genuinely varies.
"""
import sympy as sp, itertools
from fractions import Fraction
from sixvertex import Z
from lorentzian import is_lorentzian, multi_factorial

def mix_ff(xi, yi, al, alp, be, bep):
    s = al*alp + be*bep
    return {'a1':al*xi, 'a2':alp*yi, 'b1':be*yi, 'b2':bep*xi, 'c1':s*xi, 'c2':yi}

def ff_ok(w):
    return sp.expand(w['a1']*w['a2'] + w['b1']*w['b2'] - w['c1']*w['c2']) == 0

def cdict(expr, vars_):
    p = sp.Poly(sp.expand(expr), *vars_)
    return {tuple(m): Fraction(int(c)) for m, c in zip(p.monoms(), p.coeffs())}

def norm(cd):
    return {a: c/multi_factorial(a) for a,c in cd.items()}

I=sp.Integer
print("family FF-valid:", ff_ok(mix_ff(sp.Symbol('x'),sp.Symbol('y'),I(1),I(1),I(1),I(1))))
tot=mono=nontriv=fails=0; sizes=[]; failex=[]
for params in [(1,1,1,1),(1,2,1,1),(2,1,1,3),(1,1,2,2),(3,1,1,2)]:
    for n,m in [(2,2),(2,3),(3,2),(2,4),(3,3)]:
        xs=sp.symbols('x1:%d'%(n+1)); ys=sp.symbols('y1:%d'%(n+1))
        w=[mix_ff(xs[i],ys[i],*[I(v) for v in params]) for i in range(n)]
        assert all(ff_ok(wi) for wi in w)
        vars_=list(xs)+list(ys)
        for top in itertools.product([0,1],repeat=m):
            for bot in itertools.product([0,1],repeat=m):
                z=Z(n,m,list(top),list(bot),[0]*n,[0]*n,w)
                if z==0: continue
                cd=cdict(z,vars_); tot+=1
                if len(cd)==1: mono+=1
                else: nontriv+=1; sizes.append(len(cd))
                ok,why=is_lorentzian(norm(cd),len(vars_))
                if not ok:
                    fails+=1
                    if len(failex)<6: failex.append((params,n,m,top,bot,why,len(cd),sp.factor(z)))
print("\nAUDIT  total nonzero Z: %d | monomial(vacuous): %d | NON-MONOMIAL: %d | support sizes %s"
      %(tot,mono,nontriv,sorted(set(sizes))))
print("NOT Lorentzian: %d"%fails)
for e in failex:
    print("\n  params=%s n=%d m=%d top=%s bot=%s (%d monomials)\n    %s\n    Z = %s"%(e[0],e[1],e[2],e[3],e[4],e[6],e[5],e[7]))
