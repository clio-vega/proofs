import sys
sys.path.insert(0,'/home/clio/projects/proofs/code/threerow-c3')
from mn import Mj
from math import comb

def v2(n):
    n=abs(n)
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2; k+=1
    return k

def jstar(a,b,m):
    # val(j) = j + 2 v2(C(m,j) M_j); interior 0<=j<=b
    vals={}
    for j in range(0,b+1):
        M=Mj((a,b,4),j,m)
        if M==0: continue
        g = comb(m,j)*M
        if g==0: continue
        vals[j] = j + 2*v2(g)
    V=min(vals.values())
    return sorted(j for j,v in vals.items() if v==V), V, vals

from collections import Counter
dist=Counter()
boxes=Counter()
examples={}
for m in range(4,28):
    n=2*m
    for a in range(4,n):
        b=n-4-a
        if b<4 or b>a: continue
        J,V,vals = jstar(a,b,m)
        dist[len(J)]+=1
        # box generators relative to min
        j0=J[0]
        offsets=tuple(sorted(x-j0 for x in J))
        boxes[offsets]+=1
        if len(J)>=2 and offsets not in examples:
            examples[offsets]=(a,b,m,J)

print("m up to 27")
print("|J*| distribution:", dict(sorted(dist.items())))
print("\nbox offset-patterns (relative to minJ*):")
for off,cnt in sorted(boxes.items(), key=lambda kv:(len(kv[0]),kv[0])):
    ex = examples.get(off,"")
    print(f"  {off}: count={cnt}  example={ex}")

# Check: is min J* always 0?
print("\nAny shape with minJ* != 0?")
bad=[]
for m in range(4,24):
    n=2*m
    for a in range(4,n):
        b=n-4-a
        if b<4 or b>a: continue
        J,V,vals=jstar(a,b,m)
        if J[0]!=0: bad.append((a,b,m,J))
print("  count:",len(bad), bad[:10])
