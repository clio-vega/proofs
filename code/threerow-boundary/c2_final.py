import sys; sys.setrecursionlimit(100000)
from mn import Mj, fdim
from math import comb
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def s2(n):return bin(int(n)).count('1')
def carries(x,y):return s2(x)+s2(y)-s2(x+y)
def v2prod(lo,hi):return sum(v2(t) for t in range(lo,hi+1))

bad_id=0; bad_pos=0; tested=0
# c=2 subtop both-even exact identity + positivity
for m in range(3,46):
    for b in range(4,m,2):
        a=2*m-2-b
        if a<b or a%2!=0: continue
        j=b+1
        if j>m: continue
        M=Mj((a,b,2),j,m)
        if M==0: continue
        D=j+2*v2(comb(m,j)*M)-2*v2(fdim((a,b,2)))
        beta=b//2; X=m-b-1   # P+1
        Did=2*beta+3-2*s2(beta+1)+2*carries(X,beta+1)+2*v2(2*(beta+1)*X+beta)-2*v2(X+beta)
        if D!=Did: bad_id+=1; 
        # positivity certificate: v2prod(X+1,X+beta+1) + v2(2(beta+1)X+beta) - v2(X+beta) >=0, and X+beta is a factor
        cert = v2prod(X+1,X+beta+1) + v2(2*(beta+1)*X+beta) - v2(X+beta)
        if cert<0: print("CERT FAIL",a,b,m)
        if D<1: bad_pos+=1; print("D<1!",a,b,m,D)
        tested+=1
print(f"c=2 subtop both-even: tested={tested} identity-bad={bad_id} D<1 cases={bad_pos}")

# Full c=2 boundary check: ALL boundary indices, ALL parities, Delta>0, b>=3
allbad=0;allt=0
for m in range(3,46):
    for b in range(3,m):
        a=2*m-2-b
        if a<b: continue
        val0=2*v2(fdim((a,b,2)))
        for j in (b+1,b+2):
            if j>m: continue
            M=Mj((a,b,2),j,m)
            if M==0: continue
            D=j+2*v2(comb(m,j)*M)-val0
            if D<=0: allbad+=1; print("BOUNDARY VIOL",a,b,m,j,D)
            allt+=1
print(f"c=2 ALL boundary Delta>0 (b>=3): tested={allt} viol={allbad}")
