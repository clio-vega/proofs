from mn import Mj
from math import comb, factorial
from fractions import Fraction as Fr
def v2(n):
    n=abs(n)
    if n==0: return 99
    k=0
    while n%2==0: n//=2; k+=1
    return k
def N3(a,b): return a*(b+4)-(b**2+3*b+8)
def N2(a,b): return a**2*b**2+7*a**2*b+12*a**2-2*a*b**3-9*a*b**2-21*a*b-20*a+b**4+2*b**3+13*b**2+16*b-8
def N1(a,b): return a**2*b**3+9*a**2*b**2+26*a**2*b+24*a**2-2*a*b**4-7*a*b**3-13*a*b**2-14*a*b+b**5-2*b**4+17*b**3-4*b**2-84*b-96
def v2prod(lo,hi):
    s=0
    for t in range(lo,hi+1): s+=v2(t)
    return s
def fhook(a,b,m): return factorial(2*m)*(a-b+1)*(a-2)*(b-3)//(factorial(a+2)*factorial(b+1)*24)

# reduced v2(R_i) formula
def v2R_reduced(a,b,m,i):
    P=m-b-4
    if b%2==1:  # b odd
        ell=(b+7)//2; E=(b+3)//2
    else:       # b even
        ell=(b+8)//2; E=(b+2)//2
    Ndef = {4:0, 3:v2(N3(a,b)), 2:v2(N2(a,b)), 1:v2(N1(a,b))}[i]
    c = {4:0,3:0,2:1,1:1}[i]
    k=5-i
    prodv = v2prod(P+k, P+ell-1)   # v2 of prod_{t=k}^{ell-1}(P+t)
    return Ndef - c - v2(a-2) + prodv - E

# true v2(R_i)
def v2R_true(a,b,m,i):
    lam=(a,b,4); G=comb(m,b+i)*Mj(lam,b+i,m); f=fhook(a,b,m)
    R=Fr(G,f)
    return v2(R.numerator)-v2(R.denominator)

bad=[]; tested=0
for m in range(8,44):
    for b in range(6,m):
        a=2*m-4-b
        if a<b: continue
        P=m-b-4
        if P<0: continue
        for i in [1,2,3,4]:
            if b+i>m: continue
            t=v2R_true(a,b,m,i); r=v2R_reduced(a,b,m,i)
            tested+=1
            if t!=r: bad.append((a,b,m,i,t,r))
print("reduced v2(R_i) match: tested=%d mismatch=%d"%(tested,len(bad)))
print(bad[:15])
