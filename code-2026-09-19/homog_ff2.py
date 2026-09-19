import sympy as sp, itertools
from fractions import Fraction
from sixvertex import Z
from lorentzian import is_lorentzian, multi_factorial
from homog_ff import hom_ff, ff_ok, to_coeff_dict, normalize

print("AUDIT: how many partition functions ran, and were any NON-TRIVIAL?")
print("(a monomial is trivially Lorentzian -- it must not dominate the sample)\n")
tot=0; mono=0; nontriv=0; fails=0; supports=[]
examples=[]
for (A,B) in [(1,1),(1,2),(3,2)]:
    for n,m in [(2,2),(2,3),(3,2),(3,3)]:
        xs=sp.symbols('x1:%d'%(n+1)); ys=sp.symbols('y1:%d'%(n+1))
        w=[hom_ff(xs[i],ys[i],sp.Integer(A),sp.Integer(B)) for i in range(n)]
        vars_=list(xs)+list(ys)
        for top in itertools.product([0,1],repeat=m):
            for bot in itertools.product([0,1],repeat=m):
                z=Z(n,m,list(top),list(bot),[0]*n,[0]*n,w)
                if z==0: continue
                cd=to_coeff_dict(z,vars_); tot+=1
                if len(cd)==1: mono+=1
                else:
                    nontriv+=1; supports.append(len(cd))
                    if len(examples)<3 and len(cd)>3:
                        examples.append(((A,B),n,m,top,bot,sp.factor(z),len(cd)))
                ok,why=is_lorentzian(normalize(cd),len(vars_))
                if not ok: fails+=1
print("total nonzero Z      :", tot)
print("  MONOMIAL (vacuous) :", mono)
print("  non-monomial       :", nontriv, " support sizes:", sorted(set(supports)))
print("  NOT Lorentzian     :", fails)
print("\nexample non-trivial partition functions:")
for e in examples:
    print("  A,B=%s n=%d m=%d top=%s bot=%s  (%d monomials)\n     Z = %s"%(e[0],e[1],e[2],e[3],e[4],e[6],e[5]))
