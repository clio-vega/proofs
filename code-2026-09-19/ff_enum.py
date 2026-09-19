"""Q177 DIRECT: enumerate homogeneous nonnegative FREE-FERMION six-vertex weight
families and test whether the normalized partition function is Lorentzian.

Weights: each of a1,a2,b1,b2,c1,c2 is p*x_i + q*y_i with p,q in {0,1,2}, not both 0.
Free-fermion: a1*a2 + b1*b2 = c1*c2 identically in (x,y).
"""
import sympy as sp, itertools
from fractions import Fraction
from sixvertex import Z
from lorentzian import is_lorentzian, is_M_convex, multi_factorial

X,Y = sp.symbols('X Y')
CO = [(p,q) for p in range(3) for q in range(3) if (p,q)!=(0,0)]

def lf(pq,xi,yi): return pq[0]*xi+pq[1]*yi

fams=[]
for pat in itertools.product(CO,repeat=6):
    a1,a2,b1,b2,c1,c2=[lf(p,X,Y) for p in pat]
    if sp.expand(a1*a2+b1*b2-c1*c2)==0:
        fams.append(pat)
print("homogeneous nonneg FREE-FERMION weight families with coeffs in {0,1,2}: %d"%len(fams))
# keep genuinely six-vertex ones (all six weights nonzero) and dedupe by x<->y symmetry
six=[p for p in fams if all(q!=(0,0) for q in p)]
print("  of these, all six weights nonzero: %d"%len(six))

def cdict(expr,vars_):
    P=sp.Poly(sp.expand(expr),*vars_)
    return {tuple(m):Fraction(int(c)) for m,c in zip(P.monoms(),P.coeffs())}

tot=mono=nontriv=fails=mconv_fail=eig_fail=0
passing=[]
for pat in six[:40]:
    for n,m in [(2,2),(2,3)]:
        xs=sp.symbols('x1:%d'%(n+1)); ys=sp.symbols('y1:%d'%(n+1))
        w=[{k:lf(p,xs[i],ys[i]) for k,p in zip(['a1','a2','b1','b2','c1','c2'],pat)}
           for i in range(n)]
        vars_=list(xs)+list(ys)
        for top in itertools.product([0,1],repeat=m):
            for bot in itertools.product([0,1],repeat=m):
                z=Z(n,m,list(top),list(bot),[0]*n,[0]*n,w)
                if z==0: continue
                cd=cdict(z,vars_); tot+=1
                if len(cd)==1: mono+=1; continue
                nontriv+=1
                nd={a:c/multi_factorial(a) for a,c in cd.items()}
                ok,why=is_lorentzian(nd,len(vars_))
                if ok: passing.append((pat,n,m,top,bot,len(cd),sp.factor(z)))
                else:
                    fails+=1
                    if '(L2)' in why: mconv_fail+=1
                    else: eig_fail+=1
print("\n--- results over %d FF families x boundaries ---"%min(40,len(six)))
print("total nonzero Z : %d | monomial: %d | NON-MONOMIAL: %d"%(tot,mono,nontriv))
print("  NON-MONOMIAL and LORENTZIAN     : %d"%(nontriv-fails))
print("  NON-MONOMIAL and NOT Lorentzian : %d  (M-convexity %d, eigenvalue %d)"%(fails,mconv_fail,eig_fail))
print("\nexamples of NON-MONOMIAL LORENTZIAN free-fermion partition functions:")
for e in passing[:8]:
    print("   weights=%s n=%d m=%d top=%s bot=%s (%d monomials)\n      Z = %s"%(e[0],e[1],e[2],e[3],e[4],e[5],e[6]))
