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
def fhook(a,b,m):  # f^lambda for (a,b,4)
    return factorial(2*m)*(a-b+1)*(a-2)*(b-3)//(factorial(a+2)*factorial(b+1)*24)
def prod_range(lo,hi):  # inclusive product
    r=1
    for s in range(lo,hi+1): r*=s
    return r

def R_form(a,b,m,i):
    P=m-b-4
    den_consec = prod_range(2*P+b+7, 2*P+2*b+8)   # b+2 terms
    if i==4:
        num = prod_range(P+1,P+b+4)
        return Fr(num, (2*P+5)*(2*P+b+2)*den_consec)
    if i==3:
        num = N3(a,b)*prod_range(P+2,P+b+4)
        return Fr(num, (2*P+5)*(2*P+b+2)*den_consec)
    if i==2:
        num = N2(a,b)*prod_range(P+3,P+b+4)
        return Fr(num, 2*(2*P+5)*(2*P+b+2)*den_consec)
    if i==1:
        num = N1(a,b)*prod_range(P+4,P+b+4)
        return Fr(num, 6*(2*P+b+2)*den_consec)

bad=[]; tested=0
for m in range(8,42):
    for b in range(4,m):
        a=2*m-4-b
        if a<b: continue
        P=m-b-4
        if P<0: continue
        lam=(a,b,4); f=fhook(a,b,m)
        for i in [1,2,3,4]:
            j=b+i
            if j>m: continue
            G=comb(m,j)*Mj(lam,j,m)
            Rtrue=Fr(G,f)
            Rform=R_form(a,b,m,i)
            tested+=1
            if Rtrue!=Rform:
                bad.append((lam,m,i,Rtrue,Rform))
print("exact rational R_i match: tested=%d mismatches=%d"%(tested,len(bad)))
print(bad[:10])
