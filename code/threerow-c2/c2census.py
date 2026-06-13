import sys; sys.path.insert(0,'/home/clio/projects/code/threerow-c1')
from mn import Mj
from math import comb
def v2(n):
    n=abs(int(n));
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def vC(n,k):
    if k<0 or k>n:return 10**9
    return v2(comb(n,k))
def jstar(a,b):
    m=(a+b+2)//2
    vals={}
    for j in range(0,m+1):
        Mjv=Mj((a,b,2),j,m)
        if Mjv==0: continue
        vals[j]=j+2*(vC(m,j)+v2(Mjv))
    V=min(vals.values())
    return sorted([j for j in vals if vals[j]==V])

from collections import Counter
cnt=Counter(); examples={}
for A in range(2,40):
    for B in range(2,A+1):
        if (A+B)%2!=0: continue
        m=(A+B+2)//2
        if m>30: continue
        js=jstar(A,B)
        off=tuple(x-js[0] for x in js)
        cnt[off]+=1
        if off not in examples: examples[off]=(A,B,js)
print("c=2 offset census (m<=30):")
for off,c in sorted(cnt.items(), key=lambda x:(len(x[0]),x[0])):
    print(f"  {off}: count={c}  example a,b,J*={examples[off]}")
# parity of a within each box
print("\n--- box vs (a mod 4, b mod 4) ---")
from collections import defaultdict
patt=defaultdict(Counter)
for A in range(2,40):
    for B in range(2,A+1):
        if (A+B)%2!=0: continue
        m=(A+B+2)//2
        if m>30: continue
        js=jstar(A,B); off=tuple(x-js[0] for x in js)
        patt[off][(A%4,B%4)]+=1
for off in sorted(patt,key=lambda x:(len(x),x)):
    print(f" off={off}: {dict(patt[off])}")
