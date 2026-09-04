"""E-basis diagnosis of the Day 158 'log W = d Xi is FALSE' claim. Clio 2026-09-03."""
import sympy as sp
from sympy.polys.polyfuncs import symmetrize

u1,u2,u3 = sp.symbols('u1 u2 u3'); U=(u1,u2,u3)
E1,E2,E3 = sp.symbols('E1 E2 E3')
NB=6
V=sp.expand((u1-u2)*(u1-u3)*(u2-u3)); e2=u1*u2+u1*u3+u2*u3
rising=lambda x,m: sp.prod([x+j for j in range(m)]) if m>0 else sp.Integer(1)
def Tplus(p):
    P=sp.Poly(sp.expand(p),u1,u2,u3); return sp.expand(sum(c*sp.prod([rising(U[i],a[i]) for i in range(3)]) for a,c in P.terms()))
def PsiPlus(f):
    q,r=sp.div(sp.Poly(Tplus(sp.expand(f*V)),u1),sp.Poly(V,u1)); assert r.is_zero; return sp.expand(q.as_expr())
FP=[];b_=sp.Integer(1)
for b in range(NB+1):
    FP.append(sp.expand(PsiPlus(b_)/sp.factorial(b))); b_=sp.expand(b_*e2)
N=NB
mul=lambda a,b:[sp.expand(sum(a[i]*b[n-i] for i in range(n+1))) for n in range(N+1)]
def inv(a):
    r=[sp.Integer(0)]*(N+1); r[0]=sp.Integer(1)
    for n in range(1,N+1): r[n]=sp.expand(-sum(a[j]*r[n-j] for j in range(1,n+1)))
    return r
def logser(a):
    ap=[(m+1)*a[m+1] for m in range(N)]+[sp.Integer(0)]; d=mul(ap,inv(a))
    return [sp.Integer(0)]+[sp.expand(d[n-1]/n) for n in range(1,N+1)]
def homog(p,deg):
    if sp.expand(p)==0: return sp.Integer(0)
    P=sp.Poly(sp.expand(p),u1,u2,u3)
    return sp.expand(sum(c*u1**a*u2**b*u3**cc for (a,b,cc),c in P.terms() if a+b+cc==deg))
def toE(p):
    if sp.expand(p)==0: return sp.Integer(0)
    s,rem,_=symmetrize(sp.expand(p),[u1,u2,u3],formal=True)
    assert sp.expand(rem)==0, f"not symmetric: {rem}"
    s1,s2,s3=sp.symbols('s1 s2 s3')
    return sp.expand(s.subs({s1:E1,s2:E2,s3:E3}))

X=logser(FP)
Xi=[homog(X[n],n+1) for n in range(N+1)]; X0=[homog(X[n],n) for n in range(N+1)]
tauFP=[sp.expand(f.subs({u1:u1+1,u2:u2+1,u3:u3+1},simultaneous=True)) for f in FP]
H=mul(tauFP,inv(FP)); Wc=[homog(H[n],n) for n in range(N+1)]; logW=logser(Wc)

XiE=[toE(x) for x in Xi]; X0E=[toE(x) for x in X0]; logWE=[toE(x) for x in logW]
d_E = lambda f: sp.expand(3*sp.diff(f,E1)+2*E1*sp.diff(f,E2)+E2*sp.diff(f,E3))

print("="*76); print("A. d = 3 d_E1 + 2E1 d_E2 + E2 d_E3  agrees with  sum_i d/du_i  on Xi"); print("="*76)
for n in range(1,N+1):
    lhs=d_E(XiE[n]); rhs=toE(sp.expand(sum(sp.diff(Xi[n],v) for v in U)))
    print(f"   n={n}: {'OK' if sp.expand(lhs-rhs)==0 else 'MISMATCH'}")

print("\n"+"="*76); print("B. (P1) log W = d Xi   in the E-basis, E_3 FREE"); print("="*76)
for n in range(1,N+1):
    print(f"   n={n}: {'MATCH' if sp.expand(logWE[n]-d_E(XiE[n]))==0 else 'DIFFERS'}")

print("\n"+"="*76)
print("C. THE SLIP: restrict to E_3=0.  d has a term E_2 d_E3 that the slice cannot see.")
print("="*76)
z={E3:0}
for n in range(1,N+1):
    full   = sp.expand(d_E(XiE[n]).subs(z))            # differentiate THEN restrict  (correct)
    insl   = sp.expand(sp.expand(3*sp.diff(XiE[n],E1)+2*E1*sp.diff(XiE[n],E2)).subs(z))
    missing= sp.expand(E2*sp.diff(XiE[n],E3).subs(z))
    lw     = sp.expand(logWE[n].subs(z))
    print(f"   n={n}:  logW|E3=0 - (d Xi)|E3=0 = {sp.expand(lw-full)}"
          f"   |  dropped term E_2 dXi/dE_3|E3=0 = {sp.factor(missing)}"
          f"   |  accounts for gap: {sp.expand((lw-insl)-missing)==0}")

print("\n"+"="*76)
print("D. Day 156 line 130: D := X^(0) - (1/2)log W is O(E_3), leading 4E3T^3+15E1E3T^4+(36E1^2+24E2)E3T^5")
print("="*76)
for n in range(1,N+1):
    D=sp.expand(X0E[n]-logWE[n]/2)
    print(f"   [T^{n}] D = {sp.factor(D)}   (E_3=0 part: {sp.expand(D.subs(z))})")
