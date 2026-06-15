from mn import Mj, fdim
from math import comb
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def s2(n):return bin(int(n)).count('1')
def v2prod(lo,hi): # v2 of prod_{t=lo}^{hi} t
    s=0
    for t in range(lo,hi+1): s+=v2(t)
    return s

# (1) product lemma: for Q>=4, all P>=0: v2(prod_{t=1}^Q (P+t)) - 1 >= v2(P+Q-1)
bad=0
for Q in range(4,60):
    for P in range(0,500):
        lhs=v2prod(P+1,P+Q)-1
        rhs=v2(P+Q-1)
        if lhs<rhs: bad+=1
print("product lemma Q>=4 bad=",bad)
# also check Q=3 (would be needed if b=2, but b>=3 even => Q>=4; check Q=3 fails sometimes)
b3=0
for P in range(0,200):
    if v2prod(P+1,P+3)-1 < v2(P+3-1): b3+=1
print("Q=3 failures:",b3,"(expected some, that's why b>=3 needed)")

# (2) c=2 both-even: Delta(b+2) = 2 + 2[ v2prod(P+1,P+Q) - 1 - v2(P+Q-1) ], Q=b/2+2
bad2=0;tested=0
for m in range(3,200):
    for b in range(4,m,2):  # b even >=4
        a=2*m-2-b
        if a<b or a%2!=0: continue
        lam=(a,b,2); j=b+2
        if j>m: continue
        M=Mj(lam,j,m)
        if M==0: continue
        D=j+2*v2(comb(m,j)*M)-2*v2(fdim(lam))
        P=m-b-2; Q=b//2+2
        Dpred=2+2*(v2prod(P+1,P+Q)-1-v2(P+Q-1))
        if D!=Dpred: bad2+=1; 
        tested+=1
print(f"c=2 both-even Delta(b+2) product identity: tested={tested} bad={bad2}")
