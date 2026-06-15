import sys; sys.setrecursionlimit(100000)
from mn import Mj, fdim
from math import comb
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def v2prod(lo,hi):return sum(v2(t) for t in range(lo,hi+1))
# c=3 a-even top: Delta(b+3) = 2*[v2prod(P+1,P+Q') - v2(P+2)], Q'=(b+5)//2 ; need >=2
bad=0;t=0
for m in range(5,46):
    for b in range(7,m,2):     # b odd (a even), b>=7
        a=2*m-3-b
        if a<b or a%2!=0: continue
        j=b+3
        if j>m: continue
        M=Mj((a,b,3),j,m)
        if M==0: continue
        D=j+2*v2(comb(m,j)*M)-2*v2(fdim((a,b,3)))
        P=m-b-3; Q=(b+5)//2
        Dpred=2*(v2prod(P+1,P+Q)-v2(P+2))
        if D!=Dpred: bad+=1
        if D<2: print("D<2!",a,b,m,D)
        t+=1
print(f"c=3 a-even TOP: tested={t} identity-bad={bad}")
