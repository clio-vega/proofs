import sys
sys.path.insert(0,'/home/clio/projects/proofs/code/threerow-c3')
from mn import Mj
from math import comb
from collections import Counter

def v2(n):
    n=abs(n)
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2; k+=1
    return k

def jstar(a,b,m):
    vals={}
    for j in range(0,b+1):
        M=Mj((a,b,4),j,m)
        if M==0: continue
        g = comb(m,j)*M
        if g==0: continue
        vals[j] = j + 2*v2(g)
    V=min(vals.values())
    return sorted(j for j,v in vals.items() if v==V)

tie_j = Counter()    # which nonzero j ever in J*
boxes = Counter()
examples = {}
dist = Counter()
for m in range(4,42):
    n=2*m
    for a in range(4,n):
        b=n-4-a
        if b<4 or b>a: continue
        J=jstar(a,b,m)
        dist[len(J)]+=1
        for x in J[1:]:
            tie_j[x]+=1
        off=tuple(sorted(x-J[0] for x in J))
        boxes[off]+=1
        if len(J)>=2 and off not in examples:
            examples[off]=(a,b,m,J)

print("m up to 41")
print("|J*| dist:", dict(sorted(dist.items())))
print("nonzero tie j-values:", dict(sorted(tie_j.items())))
print("box offsets:")
for off,cnt in sorted(boxes.items(), key=lambda kv:(len(kv[0]),kv[0])):
    print(f"  {off}: {cnt}  ex={examples.get(off,'')}")
