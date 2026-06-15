from mn import Mj, fdim
from math import comb, factorial
def v2(n):
    n=abs(int(n))
    if n==0: return None
    k=0
    while n%2==0: n//=2;k+=1
    return k
def s2(n): return bin(int(n)).count('1')

# Test hook formula f^{(a,b,c)} = (2m)!(a-c+2)(a-b+1)(b-c+1)/[(a+2)!(b+1)!c!]
badf=0; tf=0
for c in range(1,6):
    for m in range(c+1,30):
        for b in range(c,m):
            a=2*m-b-c
            if a<b: continue
            lam=(a,b,c)
            f=fdim(lam)
            num=factorial(2*m)*(a-c+2)*(a-b+1)*(b-c+1)
            den=factorial(a+2)*factorial(b+1)*factorial(c)
            assert num%den==0
            if num//den!=f: badf+=1; 
            tf+=1
print(f"hook formula: tested={tf} bad={badf}")

# Test top boundary M_{b+c}=(b-c+1)(b+c)!/[(b+1)!c!]
badm=0;tm=0
for c in range(1,6):
    for m in range(c+1,30):
        for b in range(c,m):
            a=2*m-b-c
            if a<b: continue
            j=b+c
            if j>m: continue
            M=Mj((a,b,c),j,m)
            pred=(b-c+1)*factorial(b+c)//(factorial(b+1)*factorial(c))
            if M!=pred: badm+=1
            tm+=1
print(f"top M_b+c formula: tested={tm} bad={badm}")

# Test top Delta(b+c)-val(0) clean formula:
# = (b+c)+4+2 s2(m-b-c) - 2 s2(a+2) - 2 v2(a-c+2) - 2 v2(a-b+1)
badd=0;td=0
for c in range(1,6):
    for m in range(c+1,40):
        for b in range(c,m):
            a=2*m-b-c
            if a<b: continue
            j=b+c
            if j>m: continue
            M=Mj((a,b,c),j,m)
            if M==0: continue
            valbc=j+2*v2(comb(m,j)*M)
            val0=2*v2(fdim((a,b,c)))
            D=valbc-val0
            Dform=(b+c)+4+2*s2(m-b-c)-2*s2(a+2)-2*v2(a-c+2)-2*v2(a-b+1)
            if D!=Dform: 
                badd+=1
                if badd<=8: print(f"  D MISMATCH c={c} lam={(a,b,c)} m={m}: D={D} form={Dform}")
            td+=1
print(f"top Delta(b+c) formula: tested={td} bad={badd}")
