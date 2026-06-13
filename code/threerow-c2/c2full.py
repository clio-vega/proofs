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
    if k<0 or k>n:return 10**9
    return v2(comb(n,k))
def jstar_full(a,b):
    m=(a+b+2)//2
    vals={}
    for j in range(0,m+1):
        Mjv=Mj((a,b,2),j,m)
        if Mjv==0: continue
        vals[j]=j+2*(vC(m,j)+v2(Mjv))
    V=min(vals.values())
    return sorted(j for j in vals if vals[j]==V), m
# Check J* in {{0},{0,2},{0,4}} for ALL c=2 shapes, and whether any minimizer j>b
bad=0; boundary_min=0; ex=[]
for A in range(2,70):
    for B in range(2,A+1):
        if (A+B)%2: continue
        m=(A+B+2)//2
        if m>34: continue
        js,m=jstar_full(A,B)
        if js[0]!=0 or set(js)-{0,2,4} or len(js)>2:
            bad+=1
            if len(ex)<10: ex.append((A,B,js))
        if any(j>B for j in js):
            boundary_min+=1
            if len(ex)<10: ex.append(('BD',A,B,js,'b=',B))
print("c=2 full J* check (m<=34): violations of J*∈{{0},{0,2},{0,4}}:",bad)
print("cases where a minimizer j>b occurs:",boundary_min)
print("examples:",ex)
# B(j) checks
def s2(n): return bin(n).count('1')
def B(j): return j+2-2*s2(j)-2*v2(j)
b0=min(B(j) for j in range(1,100000))
bne=min(B(j) for j in range(1,100000) if j not in (2,4))
print("min B(j) over j>=1:",b0,"  min B(j) over j∉{2,4}:",bne)
