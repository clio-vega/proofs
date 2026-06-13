import sys; sys.path.insert(0,'/home/clio/projects/code/threerow-c1')
from mn import Mj
from math import comb
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def val(a,b,j):
    m=(a+b+2)//2
    Mjv=Mj((a,b,2),j,m)
    if Mjv==0: return None
    return j+2*(v2(comb(m,j))+v2(Mjv))
# residual: for b>=3, val(b+1)>val(0) and val(b+2)>val(0)
minmargin=99; bad=0; worst=None
for A in range(3,160):
    for B in range(3,A+1):
        if (A+B)%2: continue
        m=(A+B+2)//2
        v0=val(A,B,0)
        for J in (B+1,B+2):
            if J>m: continue
            vj=val(A,B,J)
            if vj is None: continue
            d=vj-v0
            if d<=0: bad+=1
            if d<minmargin: minmargin=d; worst=(A,B,J,d)
print("residual val(b+1),val(b+2) > val(0) for b>=3: violations:",bad)
print("min margin:",minmargin,"at",worst)
# Also: which boundary point is the even one, and its margin specifically
print("\neven boundary point margin (the only tie-candidate):")
mm=99;w=None
for A in range(3,160):
    for B in range(3,A+1):
        if (A+B)%2: continue
        m=(A+B+2)//2
        v0=val(A,B,0)
        Je = B+2 if B%2==0 else B+1
        if Je>m: continue
        vj=val(A,B,Je)
        if vj is None: continue
        if vj-v0<mm: mm=vj-v0; w=(A,B,Je,vj-v0)
print("min even-boundary margin:",mm,"at",w)
