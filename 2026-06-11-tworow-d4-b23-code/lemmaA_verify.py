import dfour_tworow as T
import sympy as sp
from fractions import Fraction
m=T.m

def v2int(n):
    n=abs(int(n))
    if n==0:return 10**9
    c=0
    while n%2==0:n//=2;c+=1
    return c
def v2frac(x):
    x=Fraction(x)
    if x==0:return 10**9
    return v2int(x.numerator)-v2int(x.denominator)

def monic_mod2_ok(b):
    """Return True if monic Q_b is 2-integral and reduces to X*(X^2+X+1)^((d-1)/2)."""
    Q=T.Qb(b); lc=Q.LC(); d=b//2
    coeffs=[Fraction(sp.Rational(c)) for c in (Q*(1/lc)).all_coeffs()]
    for c in coeffs:
        if v2frac(c)<0: return False
    x=sp.symbols('x')
    F2=[ (c.numerator*pow(c.denominator,-1,2))%2 if c!=0 else 0 for c in coeffs]
    p2=sum(int(F2[i])*x**(len(F2)-1-i) for i in range(len(F2)))
    target=(x*(x**2+x+1)**((d-1)//2))
    return sp.Poly(p2-target,x,modulus=2).is_zero

def hensel_root(b,K):
    Qm=T.qb_monic(b); coeffs=[int(c) for c in Qm.all_coeffs()]
    def ev(x,mod):
        r=0
        for c in coeffs: r=(r*x+c)%mod
        return r
    x=0
    for k in range(1,K+1):
        mod2=2**k
        if ev(x,mod2)%mod2==0: continue
        val=(ev(x,mod2)//(2**(k-1)))%2
        x=x+val*2**(k-1)
    return x%(2**K)

def lemmaA_offset(b, span=300):
    """Return set of c_b = v2(I)-forced-v2(m-rho) over a span of m; should be singleton."""
    P=T.Ib(b); R=(b-1)//2
    K=2*span.bit_length()+200; rho=hensel_root(b,K)
    s=set()
    for m0 in range(b,b+span):
        lhs=v2frac(P.eval(m0))
        forced=sum(v2int(m0-r) for r in range(0,R+1))
        dd=(m0-rho)%(2**K)
        vrho=v2int(dd) if dd!=0 else K
        s.add(lhs-forced-vrho)
    return s

if __name__=='__main__':
    import sys
    bmax=int(sys.argv[1]) if len(sys.argv)>1 else 80
    bad_fac=[]; bad_lemma=[]
    for b in range(6,bmax+1):
        if b%4 not in (2,3): continue
        if not monic_mod2_ok(b): bad_fac.append(b)
        s=lemmaA_offset(b, span=120)
        if len(s)!=1: bad_lemma.append((b,s))
    print("monic mod2 factorization X(X^2+X+1)^((d-1)/2) FAILS at:", bad_fac if bad_fac else "NONE", "(b<=%d)"%bmax)
    print("Lemma A single-offset FAILS at:", bad_lemma if bad_lemma else "NONE")
