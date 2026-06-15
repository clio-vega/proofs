from mn import Mj, fdim
from math import comb
def v2(n):
    n=abs(int(n))
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2;k+=1
    return k
def s2(n): return bin(int(n)).count('1')

# c=2 formulas
bad=0;tested=0
for m in range(3,60):
    for b in range(2,m):
        a=2*m-2-b
        if a<b: continue
        lam=(a,b,2)
        val0=2*v2(fdim(lam))
        # b+1
        j=b+1
        if j<=m:
            M=Mj(lam,j,m)
            if M>0:
                D=j+2*v2(comb(m,j)*M)-val0
                W1=a*(b+2)-b*(b+1)
                Df=(b+3)+2*s2(m-b-1)-2*s2(a+2)+2*v2(W1)-2*v2(a)-2*v2(a-b+1)
                if D!=Df: 
                    bad+=1
                    if bad<=5:print(f"  b+1 MISMATCH {lam} m={m}: D={D} f={Df}")
                tested+=1
        # b+2 (top)
        j=b+2
        if j<=m:
            M=Mj(lam,j,m)
            if M>0:
                D=j+2*v2(comb(m,j)*M)-val0
                Df=(b+2)+4+2*s2(m-b-2)-2*s2(a+2)-2*v2(a)-2*v2(a-b+1)
                if D!=Df:
                    bad+=1
                    if bad<=8:print(f"  b+2 MISMATCH {lam} m={m}: D={D} f={Df}")
                tested+=1
print(f"c=2 formulas: tested={tested} bad={bad}")

# margins: min Delta(b+1), Delta(b+2) over valid range b>=3 (box interior b>=4 actually for {0,4})
# residual is b>=3 (b=2 handled). Let's see min Delta for b>=3.
mind1=10**9;mind2=10**9;w1=None;w2=None
for m in range(3,200):
    for b in range(3,m):
        a=2*m-2-b
        if a<b: continue
        lam=(a,b,2)
        W1=a*(b+2)-b*(b+1)
        D1=(b+3)+2*s2(m-b-1)-2*s2(a+2)+2*v2(W1)-2*v2(a)-2*v2(a-b+1) if (b+1)<=m else None
        D2=(b+2)+4+2*s2(m-b-2)-2*s2(a+2)-2*v2(a)-2*v2(a-b+1) if (b+2)<=m else None
        if D1 is not None and D1<mind1: mind1=D1;w1=(lam,m,D1)
        if D2 is not None and D2<mind2: mind2=D2;w2=(lam,m,D2)
print(f"min Delta(b+1)={mind1} at {w1}")
print(f"min Delta(b+2)={mind2} at {w2}")
