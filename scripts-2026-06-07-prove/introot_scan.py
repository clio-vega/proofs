"""Exact integer scan: find ALL integer (m,b), m>=b>=1, with I_b(m)=0.
I_b(m) = Im( [u^b]((1-u)(1+su+u^2)^m) ), s=1+i.
Compute (1+su+u^2)^m mod u^(b+1) with Gaussian-integer coeffs (exact), then
[u^b]((1-u)P) = P[b]-P[b-1]; take imaginary part.
"""
def Ib_value(b, m):
    # represent poly as list of (re,im) up to degree b
    N = b+1
    # base poly 1 + (1+i)u + u^2  -> coeffs
    base = [(0,0)]*N
    base[0]=(1,0)
    if N>1: base[1]=(1,1)
    if N>2: base[2]=(1,0)
    # result = base^m mod u^N via binary exponentiation
    def mul(A,B):
        C=[(0,0)]*N
        for i,(ar,ai) in enumerate(A):
            if ar==0 and ai==0: continue
            for j,(br,bi) in enumerate(B):
                if i+j>=N: break
                if br==0 and bi==0: continue
                cr,ci=C[i+j]
                C[i+j]=(cr+ar*br-ai*bi, ci+ar*bi+ai*br)
        return C
    R=[(0,0)]*N; R[0]=(1,0)
    P=base; e=m
    while e>0:
        if e&1: R=mul(R,P)
        e>>=1
        if e>0: P=mul(P,P)
    # [u^b]((1-u)P) = R[b]-R[b-1]
    rb=R[b]; rbm1=R[b-1] if b-1>=0 else (0,0)
    im = rb[1]-rbm1[1]
    return im

zeros=[]
import sys
BMAX=40
for b in range(1, BMAX+1):
    # scan m from b up to a generous bound past largest real root (~ c*b^2)
    Mmax = max(3*b*b, 50)
    for m in range(b, Mmax+1):
        if Ib_value(b,m)==0:
            zeros.append((b,m))
    # progress
print("All integer zeros (b,m) with m>=b, b<=%d, m<=3b^2:"%BMAX)
for z in zeros: print("   b=%d m=%d"%z)
print("count:", len(zeros))
