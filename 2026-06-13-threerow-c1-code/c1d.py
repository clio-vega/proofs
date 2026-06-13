from mn import Mj
from math import comb
def C(n,k):
    if k<0 or k>n or n<0: return 0
    return comb(n,k)
def Mclean(a,b,m,j):
    N=2*(m-j); beta=b-j
    return b*C(N,beta+1) - (b-j)*C(N,beta) - (j-1)*C(N,beta-1)
bad=0
for m in range(2,12):
    n=2*m
    for a in range(1,n+1):
        b=n-a-1
        if b<1 or b>a: continue
        for j in range(0,m+1):
            if Mj((a,b,1),j,m)!=Mclean(a,b,m,j):
                bad+=1
print("clean c=1 form mismatches=",bad)

def v2(x):
    if x==0: return None
    c=0
    while x%2==0:x//=2;c+=1
    return c
def val(a,b,m,j):
    M=Mclean(a,b,m,j); cm=comb(m,j)
    if cm*M==0: return None
    return j+2*v2(cm*M)
print("\nval(j) tables for c=1 (showing shapes incl. ties):")
for (a,b) in [(4,1),(8,1),(12,1),(7,4),(11,8),(5,4),(8,5),(6,3),(10,3),(9,4),(11,2)]:
    if (a+b+1)%2: continue
    m=(a+b+1)//2
    row=[val(a,b,m,j) for j in range(m+1)]
    V=min(x for x in row if x is not None)
    J=[j for j,x in enumerate(row) if x==V]
    print(f"  ({a},{b},1) m={m}: val={row}  J*={J}")
