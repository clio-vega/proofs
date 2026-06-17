"""Test polynomial divisibility of N_i by deficit factors (a-b+1), (a-c+2)."""
import sympy as sp
from mn import Mj
a,b=sp.symbols('a b')

def collect(c,i,amax=64,bmax=28):
    pts=[]
    for bb in range(c,bmax):
        for aa in range(bb,amax):
            if (aa+bb+c)%2: continue
            m=(aa+bb+c)//2; j=bb+i
            if j>m: continue
            pts.append((aa,bb,Mj((aa,bb,c),j,m)))
    return pts
def fit(pts,md):
    mo=[(p,q) for p in range(md+1) for q in range(md+1) if p+q<=md]
    if len(pts)<len(mo)+3: return None
    A=sp.Matrix([[x**p*y**q for(p,q)in mo]for(x,y,_)in pts]);Y=sp.Matrix([M for(_,_,M)in pts])
    s=A.solve_least_squares(Y)
    if any(A*s-Y):return None
    return sp.expand(sum(s[k]*a**p*b**q for k,(p,q)in enumerate(mo)))
def extractN(c,i):
    M=fit(collect(c,i),2*c+1)
    if M is None:return None
    bp=(b-c+1)*sp.prod([b+t for t in range(2,i+1)])
    if i==1: bp=(b-c+1)*(a-b+1)
    num,den=sp.fraction(sp.together(sp.cancel(M/bp)));num=sp.expand(num)
    g=0
    for cc in sp.Poly(num,a,b).coeffs(): g=sp.igcd(g,int(cc))
    return sp.expand(num/g)

for c in [4,5,6]:
    print("c=%d"%c)
    for i in range(1,c+1):
        N=extractN(c,i)
        if N is None: print("  i=%d nofit"%i);continue
        div_abp1 = sp.expand(N.subs(a,b-1))==0      # (a-b+1)|N
        div_acm2 = sp.expand(N.subs(a,c-2))==0      # (a-c+2)|N
        print("  i=%d(k=%d): (a-b+1)|N=%s  (a-c+2)|N=%s"%(i,c-i,div_abp1,div_acm2))
