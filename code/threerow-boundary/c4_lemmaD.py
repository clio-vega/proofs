from mn import Mj
from math import comb

def v2(n):
    n=abs(n)
    if n==0: return 99
    k=0
    while n%2==0: n//=2; k+=1
    return k
def s2(n): return bin(abs(n)).count('1')
def valMN(lam,j,m):
    M=Mj(lam,j,m)
    if M==0: return None
    return j+2*v2(comb(m,j)*M)

# closed-form deficits (c=4)
def N3(a,b): return a*(b+4)-(b**2+3*b+8)
def N2(a,b): return a**2*b**2+7*a**2*b+12*a**2-2*a*b**3-9*a*b**2-21*a*b-20*a+b**4+2*b**3+13*b**2+16*b-8
def N1(a,b): return a**2*b**3+9*a**2*b**2+26*a**2*b+24*a**2-2*a*b**4-7*a*b**3-13*a*b**2-14*a*b+b**5-2*b**4+17*b**3-4*b**2-84*b-96

def v2M(a,b,i):
    if i==4: return v2(b-3)+v2(b+2)+v2(b+3)+v2(b+4)-3
    if i==3: return v2(b-3)+v2(b+2)+v2(b+3)+v2(N3(a,b))-3
    if i==2: return v2(b-3)+v2(b+2)+v2(N2(a,b))-4
    if i==1: return v2(b-3)+v2(a-b+1)+v2(N1(a,b))-4

def Delta_formula(a,b,m,i):
    P=m-b-4
    smooth = 2*(2 - s2(2*P+b+6) - s2(b+1) + s2(b+i) + s2(P+4-i) - v2(a-b+1) - v2(a-2) - v2(b-3))
    return (b+i) + 2*v2M(a,b,i) + smooth

bad=[]; tested=0
for m in range(8,46):
    for b in range(4,m):
        a=2*m-4-b
        if a<b: continue
        P=m-b-4
        if P<0: continue
        lam=(a,b,4)
        v0=valMN(lam,0,m)
        if v0 is None: continue
        for i in [1,2,3,4]:
            j=b+i
            if j>m: continue
            vj=valMN(lam,j,m)
            if vj is None: continue
            D=vj-v0
            F=Delta_formula(a,b,m,i)
            tested+=1
            if D!=F:
                bad.append((lam,m,i,D,F))
print("tested=%d  mismatches=%d"%(tested,len(bad)))
print(bad[:20])
