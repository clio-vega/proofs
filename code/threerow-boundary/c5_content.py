"""c5_content.py -- explicit c=5 slice factorizations + content bounds (for the proof)."""
import sympy as sp
from mn import Mj
a,b,P,B=sp.symbols('a b P B')   # B=beta

def collect(c,i,amax=64,bmax=28):
    pts=[]
    for bb in range(c,bmax):
        for aa in range(bb,amax):
            if (aa+bb+c)%2: continue
            m=(aa+bb+c)//2;j=bb+i
            if j>m:continue
            pts.append((aa,bb,Mj((aa,bb,c),j,m)))
    return pts
def fit(pts,md):
    mo=[(p,q) for p in range(md+1) for q in range(md+1) if p+q<=md]
    A=sp.Matrix([[x**p*y**q for(p,q)in mo]for(x,y,_)in pts]);Y=sp.Matrix([M for(_,_,M)in pts])
    s=A.solve_least_squares(Y)
    assert not any(A*s-Y)
    return sp.expand(sum(s[k]*a**p*b**q for k,(p,q)in enumerate(mo)))
def extractN(c,i):
    M=fit(collect(c,i),2*c+1)
    bp=(b-c+1)*sp.prod([b+t for t in range(2,i+1)])
    if i==1: bp=(b-c+1)*(a-b+1)
    num,den=sp.fraction(sp.together(sp.cancel(M/bp)));num=sp.expand(num)
    g=0
    for cc in sp.Poly(num,a,b).coeffs(): g=sp.igcd(g,int(cc))
    return sp.expand(num/g)

c=5
for i in range(1,c+1):
    N=extractN(c,i)
    Ne=sp.factor(sp.expand(N.subs({a:2*P+2*B+c, b:2*B})))      # b even = 2B
    No=sp.factor(sp.expand(N.subs({a:2*P+2*B+1+c, b:2*B+1}))) # b odd = 2B+1
    print("i=%d (k=%d):"%(i,c-i))
    print("   b even: %s"%Ne)
    print("   b odd : %s"%No)
    # divisibility (a-b+1)=2P+6 -> (P+3); check
    print("   (a-b+1)|N (poly, a=b-1 -> 0): %s"%(sp.expand(N.subs(a,b-1))==0))
