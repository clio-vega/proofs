"""Does Day 158 Thm 2 actually remove Day 156's obstruction? Clio 2026-09-03."""
import sympy as sp
from sympy.polys.polyfuncs import symmetrize
u1,u2,u3=sp.symbols('u1 u2 u3'); U=(u1,u2,u3); E1,E2,E3,T=sp.symbols('E1 E2 E3 T')
NB=7; V=sp.expand((u1-u2)*(u1-u3)*(u2-u3)); e2=u1*u2+u1*u3+u2*u3
rising=lambda x,m: sp.prod([x+j for j in range(m)]) if m>0 else sp.Integer(1)
def Tplus(p):
    P=sp.Poly(sp.expand(p),u1,u2,u3); return sp.expand(sum(c*sp.prod([rising(U[i],a[i]) for i in range(3)]) for a,c in P.terms()))
def PsiPlus(f):
    q,r=sp.div(sp.Poly(Tplus(sp.expand(f*V)),u1),sp.Poly(V,u1)); assert r.is_zero; return sp.expand(q.as_expr())
FP=[];b_=sp.Integer(1)
for b in range(NB+1): FP.append(sp.expand(PsiPlus(b_)/sp.factorial(b))); b_=sp.expand(b_*e2)
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
    s,rem,_=symmetrize(sp.expand(p),[u1,u2,u3],formal=True); assert sp.expand(rem)==0
    s1,s2,s3=sp.symbols('s1 s2 s3'); return sp.expand(s.subs({s1:E1,s2:E2,s3:E3}))

X=logser(FP); tauFP=[sp.expand(f.subs({u1:u1+1,u2:u2+1,u3:u3+1},simultaneous=True)) for f in FP]
H=mul(tauFP,inv(FP))
XiE=[toE(homog(X[n],n+1)) for n in range(N+1)]
X0E=[toE(homog(X[n],n))   for n in range(N+1)]
WcE=[toE(homog(H[n],n))   for n in range(N+1)]
Lm1E=[toE(homog(H[n],n-1)) for n in range(N+1)]        # ell^top_{-1}(H)
dE=lambda f: sp.expand(3*sp.diff(f,E1)+2*E1*sp.diff(f,E2)+E2*sp.diff(f,E3))
z={E3:0}

# Theorem C.5: ell_{-1}(H)|E3=0 = 6T/q^4,  q^2 = 1 - 2 T E1 + T^2(E1^2-4E2)
q2=1-2*T*E1+T**2*(E1**2-4*E2)
tgt=sp.series(6*T/q2**2,T,0,N+1).removeO()
print("="*76); print("C.5:  ell^top_{-1}(H)|_{E3=0}  ==  6T/q^4 ?"); print("="*76)
for n in range(1,N):
    lhs=sp.expand(Lm1E[n].subs(z)); rhs=sp.expand(sp.expand(tgt).coeff(T,n))
    print(f"   n={n}: {'MATCH' if sp.expand(lhs-rhs)==0 else 'DIFFERS '+str(sp.expand(lhs-rhs))}")

print("\n"+"="*76)
print("Day 156 Lemma: M^(-1) = d X^(0) + (1/2) d^2 Xi,  and ell_{-1}(H) = W * M^(-1)")
print("Compare (i) d applied in 3 vars THEN restricted   vs  (ii) d applied INSIDE the E_3=0 slice")
print("="*76)
Mfull=[sp.expand(dE(X0E[n])+sp.Rational(1,2)*dE(dE(XiE[n]))) for n in range(N+1)]
dslice=lambda f: sp.expand(3*sp.diff(f,E1)+2*E1*sp.diff(f,E2))
Mslice=[sp.expand(dslice(X0E[n].subs(z))+sp.Rational(1,2)*dslice(dslice(XiE[n].subs(z)))) for n in range(N+1)]
# W * M^(-1)
def conv(a,b,n): return sp.expand(sum(a[i]*b[n-i] for i in range(n+1)))
for n in range(1,N):
    lhs=sp.expand(Lm1E[n].subs(z))
    r_full =sp.expand(lhs-conv([w.subs(z) for w in WcE],[m.subs(z) for m in Mfull],n))
    r_slice=sp.expand(lhs-conv([w.subs(z) for w in WcE],Mslice,n))
    print(f"   n={n}:  correct (d in 3 vars): {'MATCH' if r_full==0 else 'DIFFERS'}"
          f"   |  slice-only d: {'MATCH' if r_slice==0 else 'DIFFERS by '+str(sp.factor(r_slice))}")
print("\n  => If the second column DIFFERS, then knowing X^(0)|_{E3=0} (Day 158 Thm 2) is")
print("     NOT sufficient to complete Day 156 C.5: the route needs d X^(0)|_{E3=0}, which")
print("     depends on the E_3-linear part of X^(0) and is invisible in the E_3=0 slice.")
