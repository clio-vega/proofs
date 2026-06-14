from math import comb
from sympy import factorial
def v2(n):
    n=int(n)
    if n==0: return 10**9
    r=0
    while n%2==0: n//=2; r+=1
    return r
def falling(j,k):
    p=1
    for i in range(k): p*=(j-i)
    return p

# Step 1: subset identity  C(F+2c-1,j)*C(j,2c) == C(F+2c-1,2c)*C(F-1,j-2c)
bad1=0
for c in range(1,6):
    for F in range(2,80,2):
        for j in range(2*c,F+2*c):
            L=comb(F+2*c-1,j)*comb(j,2*c)
            R=comb(F+2*c-1,2*c)*comb(F-1,j-2*c)
            if L!=R: bad1+=1
print("subset identity failures:",bad1)

# Step: reduced bound is j-independent:  lhs >= v2 C(F+2c-1,2c)+v2((2c)!) = sum_{i=0}^{2c-1} v2(F+i)
bad2=0
for c in range(1,6):
    for F in range(2,80,2):
        anchor=v2(comb(F+2*c-1,2*c))+v2(factorial(2*c))
        ssum=sum(v2(F+i) for i in range(2*c))
        if anchor!=ssum: bad2+=1
print("anchor = sum_{i=0}^{2c-1} v2(F+i) failures:",bad2)

# Step: min over even F of sum_{i=1}^{2c-1} v2(F+i) equals (c-1)+v2((c-1)!)
for c in range(1,7):
    mn=min(sum(v2(F+i) for i in range(1,2*c)) for F in range(2,5000,2))
    print(f"c={c} min_F sum_1^(2c-1) v2(F+i) = {mn}  vs (c-1)+v2((c-1)!) = {(c-1)+v2(factorial(c-1))}")
