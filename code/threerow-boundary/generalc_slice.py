"""
generalc_slice.py -- prove session 2026-06-17

Claim to PROVE: for the deepest indices where the deficit factor is out of range,
N_i (restricted to the relevant parity slice) is divisible by that deficit factor.

  - a-c+2 = 2(P+beta+1) on the slice b=2beta (b even).  Test: (P+beta+1) | N_i?
    i.e. N_i(a=2P+2beta+c, b=2beta) == 0 when P = -(beta+1).
  - a-b+1 = 2(P+gamma), gamma=(c+1)/2 (c odd).  Test: (P+gamma) | N_i on b-even slice?
    i.e. set P = -gamma.

We compute N_i on slices and the polynomial remainder.  Also: what is the FULL
factorization of N_i on each parity slice (b=2beta and b=2beta+1)?  This reveals
the structural content (the b(b+1)-even-trick analogue).
"""
import sympy as sp
from math import factorial
from mn import Mj

a, b, P, beta = sp.symbols('a b P beta')

def collect(c, i, amax=64, bmax=28):
    pts=[]
    for bb in range(c,bmax):
        for aa in range(bb,amax):
            if (aa+bb+c)%2: continue
            m=(aa+bb+c)//2; j=bb+i
            if j>m: continue
            pts.append((aa,bb,Mj((aa,bb,c),j,m)))
    return pts

def fit_poly(pts,maxdeg):
    monos=[(p,q) for p in range(maxdeg+1) for q in range(maxdeg+1) if p+q<=maxdeg]
    if len(pts)<len(monos)+3: return None
    A=sp.Matrix([[aa**p*bb**q for (p,q) in monos] for (aa,bb,_) in pts])
    y=sp.Matrix([M for (_,_,M) in pts]); sol=A.solve_least_squares(y)
    if any(A*sol-y): return None
    return sp.factor(sp.expand(sum(sol[k]*a**p*b**q for k,(p,q) in enumerate(monos))))

def extract_N(c,i):
    M=fit_poly(collect(c,i),maxdeg=2*c+1)
    if M is None: return None
    bprod=(b-c+1)*sp.prod([b+t for t in range(2,i+1)])
    if i==1: bprod=(b-c+1)*(a-b+1)
    num,den=sp.fraction(sp.together(sp.cancel(M/bprod)))
    num=sp.expand(num); poly=sp.Poly(num,a,b); g=0
    for cc in poly.coeffs(): g=sp.igcd(g,int(cc))
    return sp.expand(num/g)

for c in [4,5,6]:
    print("="*70); print("c=%d"%c)
    gamma=(c+1)//2 if c%2 else None
    for i in range(1,c+1):
        k=c-i
        N=extract_N(c,i)
        if N is None: print("  i=%d no fit"%i); continue
        # b even slice: a=2P+2beta+c, b=2beta
        Neven=sp.expand(N.subs({a:2*P+2*beta+c, b:2*beta}))
        # factor out content 2-power and full factorization
        fe=sp.factor(Neven)
        # test (P+beta+1) | Neven
        rem_acm2 = sp.expand(Neven.subs(P,-beta-1))
        # b odd slice: a=2P+2beta+1+c, b=2beta+1
        Nodd=sp.expand(N.subs({a:2*P+2*beta+1+c, b:2*beta+1}))
        fo=sp.factor(Nodd)
        line="  i=%d(k=%d):"%(i,k)
        line+=" (P+beta+1)|N[beven]? %s"%("YES" if rem_acm2==0 else "no(rem=%s)"%rem_acm2)
        if c%2==1 and k<=c-2:
            rem_abp1=sp.expand(Neven.subs(P,-gamma))
            line+=" ; (P+gamma)|N[beven]? %s"%("YES" if rem_abp1==0 else "no")
            rem_abp1o=sp.expand(Nodd.subs(P,-gamma))
            line+=" ; (P+gamma)|N[bodd]? %s"%("YES" if rem_abp1o==0 else "no")
        print(line)
        if i<=2:  # show factorizations of the deep ones
            print("       N[b even] = %s"%fe)
            print("       N[b odd]  = %s"%fo)
