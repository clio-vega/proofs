from math import comb
def C(n,k):
    if k<0 or k>n or n<0: return 0
    return comb(n,k)
def v2(x):
    x=abs(x); c=0
    while x%2==0:x//=2;c+=1
    return c
def s2(n): return bin(n).count('1')
def Mfac(a,b,m,j):
    if j<=b-1:
        N=2*(m-j); return C(N,b-j-1)*(a-b+1)*(b*(a+1)-j*(j-1))//((b-j)*(b-j+1))
    if j==b: return 2*b*(m-b)
    if j==b+1: return b
    return 0
def valf(a,b,m,j):
    M=Mfac(a,b,m,j); cm=comb(m,j) if 0<=j<=m else 0
    if cm*M==0: return None
    return j+2*v2(cm*M)

# 1) standard ineq v2C(n,k)>=v2(n)-v2(k)
i1=all(v2(C(n,k))>=v2(n)-v2(k) for n in range(1,200) for k in range(1,n+1))
# 2) compensation lemma, ALL c=1, all 1<=j<=b-1
i2=True
for m in range(2,30):
    n=2*m
    for a in range(1,n+1):
        b=n-a-1
        if b<1 or b>a: continue
        K=b*(a+1)
        for j in range(1,b):
            if v2(C(a+2,j))+v2(K-j*(j-1)) < v2(a+1): i2=False
# 3) full box + tie conditions, m<=80
box_ok=True; tie_ok=True
from collections import Counter
for m in range(2,81):
    n=2*m
    for a in range(1,n+1):
        b=n-a-1
        if b<1 or b>a: continue
        vals={j:valf(a,b,m,j) for j in range(b+2) if valf(a,b,m,j) is not None}
        V=min(vals.values()); J=sorted(j for j,v in vals.items() if v==V)
        if a%2==0:   # K odd
            if not set(J)<={0,2}: box_ok=False
            tie_pred=(a%4==0 and b%4==1)
        else:        # K even
            if not set(J)<={3,5}: box_ok=False
            tie_pred=(a%4==3 and b%4==0)
        is_tie=(len(J)==2)
        if is_tie!=tie_pred: tie_ok=False
print("1) v2C(n,k)>=v2(n)-v2(k):", i1)
print("2) compensation lemma (all c=1, 1<=j<=b-1):", i2)
print("3a) box J* subset {0,2}(a even)/{3,5}(a odd), m<=80:", box_ok)
print("3b) tie<=> [a even:a=0,b=1 mod4] / [a odd:a=3,b=0 mod4], m<=80:", tie_ok)
# 4) boundary never strictly beats j0 (except legit ties already counted): implied by box_ok
print("==> all four:", i1 and i2 and box_ok and tie_ok)
