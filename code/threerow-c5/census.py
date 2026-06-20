import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from mn import Mj
from math import comb

def v2(n):
    n=abs(n)
    if n==0: return 10**9
    k=0
    while n%2==0:
        n//=2; k+=1
    return k

def val(lam,j,m):
    M=Mj(lam,j,m)
    if M==0: return None
    return j+2*v2(comb(m,j)*M)

def Jstar(a,b,m):
    lam=(a,b,5)
    vals={}
    for j in range(0,b+1):  # interior indices only (0..b)
        v=val(lam,j,m)
        if v is not None:
            vals[j]=v
    if not vals: return None,None
    V=min(vals.values())
    J=sorted(j for j,v in vals.items() if v==V)
    return J,V

from collections import Counter
dist=Counter()
ties_aeven=[]
ties_aodd=[]
examples={}
for m in range(5,30):
    n=2*m
    for a in range(5,n):
        b=n-a-5
        if b<5 or b>a: continue
        # box-interior: require b>=2c=10 for clean interior; but record all
        J,V=Jstar(a,b,m)
        if J is None: continue
        key=tuple(d-J[0] for d in J)  # offset-normalized box
        dist[(len(J),)]+=1
        if len(J)>1:
            rec=(a,b,m,J,a%4,b%4)
            if a%2==0: ties_aeven.append(rec)
            else: ties_aodd.append(rec)
print("size distribution:",dict(dist))
print("\n--- a EVEN ties (j0, box, a%4,b%4) ---")
for r in ties_aeven[:40]:
    a,b,m,J,a4,b4=r
    print(f" a={a} b={b} m={m} J={J} box={[d-J[0] for d in J]} a%4={a4} b%4={b4}")
print(f" total a-even ties: {len(ties_aeven)}")
print("\n--- a ODD ties (j0, box, a%4,b%4) ---")
for r in ties_aodd[:40]:
    a,b,m,J,a4,b4=r
    print(f" a={a} b={b} m={m} J={J} box={[d-J[0] for d in J]} a%4={a4} b%4={b4}")
print(f" total a-odd ties: {len(ties_aodd)}")
