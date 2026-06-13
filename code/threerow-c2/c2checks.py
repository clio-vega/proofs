import sys; sys.path.insert(0,'/home/clio/projects/code/threerow-c1')
from mn import Mj
from math import comb
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def vC(n,k):
    if k<0 or k>n: return 10**9
    return v2(comb(n,k))
def val(a,b,j):
    m=(a+b+2)//2
    M=Mj((a,b,2),j,m)
    if M==0: return None
    return j+2*(v2(comb(m,j))+v2(M))
# b=2 hand formulas
bad=0
for A in range(4,200,2):  # a even, b=2
    m=(A+4)//2
    v0=val(A,2,0)
    # val(4)-val(0)=2[v2(m-1)+v2(m-3)]
    if val(A,2,4) is not None and val(A,2,4)-v0 != 2*(v2(m-1)+v2(m-3)): bad+=1
    # val(2)-val(0)=2 v2(a+2)
    if val(A,2,2)-v0 != 2*v2(A+2): bad+=1
    # val(3)-val(0)=1+2 v2(m-1)
    if val(A,2,3)-v0 != 1+2*v2(m-1): bad+=1
print("b=2 hand-formula mismatches:",bad)
# never-both: T(4)>=0 on (0,0),(3,1) mod4 classes
def Q(a,b,j): return a*(b-1)*((a+3)*(b+2)-2*j*j)+j*(j-1)*(j-2)*(j-3)
badT=0
for A in range(4,400):
    for B in range(4,A+1):
        if (A+B)%2: continue
        if (A%4,B%4) in [(0,0),(3,1)]:
            T4=vC(A+3,4)+vC(B+2,4)+v2(Q(A,B,4))-v2(Q(A,B,0))
            if T4<0: badT+=1
print("T(4)<0 on Delta(2)=0 classes:",badT)
