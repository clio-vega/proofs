import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from mn import Mj
from math import comb
from collections import Counter, defaultdict

def v2(n):
    n=abs(n)
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2; k+=1
    return k
def val(lam,j,m):
    M=Mj(lam,j,m)
    if M==0: return None
    return j+2*v2(comb(m,j)*M)
def Jstar_interior(a,b,m):
    lam=(a,b,5); vals={}
    for j in range(0,b+1):
        v=val(lam,j,m)
        if v is not None: vals[j]=v
    if not vals: return None
    V=min(vals.values())
    return sorted(j for j,v in vals.items() if v==V)

sizedist=Counter()
tieclass=defaultdict(Counter)  # (parity) -> Counter of (a%4,b%4)->count of ties, sizes
maxsize=0; bigexamples=[]
for m in range(5,42):
    n=2*m
    for a in range(5,n):
        b=n-a-5
        if b<10 or b>a: continue   # box-interior b>=2c=10
        J=Jstar_interior(a,b,m)
        if J is None: continue
        sizedist[len(J)]+=1
        if len(J)>maxsize:
            maxsize=len(J); bigexamples=[(a,b,m,J)]
        elif len(J)==maxsize and len(J)>2 and len(bigexamples)<5:
            bigexamples.append((a,b,m,J))
        if len(J)>1:
            par='aeven' if a%2==0 else 'aodd'
            tieclass[par][(a%4,b%4,len(J))]+=1
print("interior(b>=10) size dist:", dict(sizedist))
print("max |J*| =", maxsize, " examples:", bigexamples)
for par in ('aeven','aodd'):
    print(f"\n{par} tie classes (a%4,b%4,size)->count:")
    for k in sorted(tieclass[par]): print("  ",k, tieclass[par][k])
