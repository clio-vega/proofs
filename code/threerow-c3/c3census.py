from math import comb
from mn import Mj

def v2(n):
    n=int(n)
    if n==0: return None
    c=0
    while n%2==0: n//=2; c+=1
    return c

def Jstar(av,bv,mv):
    vals={}
    for jv in range(0,bv+4):
        if jv>mv: break
        M=Mj((av,bv,3),jv,mv)
        if M==0: continue
        cm=comb(mv,jv)
        if cm==0: continue
        vals[jv]=jv+2*(v2(cm)+v2(M))
    V=min(vals.values())
    J=sorted(j for j,v in vals.items() if v==V)
    return J,V

# census by (a,b) mod 8
from collections import defaultdict
box_by_class=defaultdict(set)
size_count=defaultdict(int)
examples={}
for mv in range(3,30):
    n=2*mv
    for av in range(1,n+1):
        for bv in range(3,av+1):
            if n-av-bv!=3: continue
            J,V=Jstar(av,bv,mv)
            j0=min(J)
            box=tuple(sorted(x-j0 for x in J))
            cls=(av%8,bv%8)
            box_by_class[cls].add((j0%8,box))
            size_count[len(J)]+=1
            if len(J)==4 and 4 not in examples:
                examples[4]=(av,bv,mv,J)
            if len(J)==3:
                examples.setdefault('THREE',[]).append((av,bv,mv,J))
            if len(J)>=5:
                examples.setdefault('BIG',[]).append((av,bv,mv,J))

print("size distribution of |J*|:", dict(size_count))
print("\nexample |J*|=4:", examples.get(4))
print("any |J*|=3:", examples.get('THREE','NONE')[:5] if 'THREE' in examples else 'NONE')
print("any |J*|>=5:", examples.get('BIG','NONE')[:5] if 'BIG' in examples else 'NONE')
print("\n(a%8,b%8) -> set of (j0%8, box):")
for cls in sorted(box_by_class):
    print(" ",cls,"->",sorted(box_by_class[cls]))
