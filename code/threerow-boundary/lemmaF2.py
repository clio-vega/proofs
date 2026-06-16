def v2(n):
    n=abs(n); k=0
    if n==0: return 99
    while n%2==0: n//=2; k+=1
    return k
import itertools
# Two-factor Lemma F: Q>=6, R>=0, any two distinct t1,t2 in [1,Q]:
#   v2(prod_{t=1}^Q (R+t)) - v2(R+t1) - v2(R+t2) >= 1
fails=[]
for Q in range(6,40):
    for R in range(0,600):
        tot=sum(v2(R+t) for t in range(1,Q+1))
        for t1,t2 in itertools.combinations(range(1,Q+1),2):
            if tot - v2(R+t1) - v2(R+t2) < 1:
                fails.append((Q,R,t1,t2))
print("two-factor F (Q>=6) failures:", fails[:5], "count", len(fails))
# show Q=5 CAN fail (justifies Q>=6 threshold):
q5fail=[]
for R in range(0,200):
    tot=sum(v2(R+t) for t in range(1,6))
    for t1,t2 in itertools.combinations(range(1,6),2):
        if tot - v2(R+t1) - v2(R+t2) < 1: q5fail.append((R,t1,t2))
print("Q=5 example failures:", q5fail[:3], "count", len(q5fail))
# single-factor F: Q>=4
sf=[]
for Q in range(4,40):
    for R in range(0,600):
        tot=sum(v2(R+t) for t in range(1,Q+1))
        for t1 in range(1,Q+1):
            if tot - v2(R+t1) < 1: sf.append((Q,R,t1))
print("single-factor F (Q>=4) failures:", sf[:5], "count", len(sf))
