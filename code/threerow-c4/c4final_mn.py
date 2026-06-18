import sys
sys.path.insert(0,'/home/clio/projects/proofs/code/threerow-c3')
from mn import Mj
from math import comb
def v2(n):
    n=abs(n)
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def jstar(a,b,m):
    vals={}
    for j in range(0,b+1):
        M=Mj((a,b,4),j,m)
        if M==0:continue
        g=comb(m,j)*M
        if g==0:continue
        vals[j]=j+2*v2(g)
    V=min(vals.values())
    return sorted(j for j,vv in vals.items() if vv==V)
# Full MN verification to m<=48: J* subset {0,2}; tie iff (a,b)=(2,2),(1,3) mod4
bad=0; tested=0; ntie=0
for m in range(4,49):
    n=2*m
    for a in range(4,n):
        b=n-4-a
        if b<4 or b>a: continue
        J=jstar(a,b,m); tested+=1
        if not set(J)<= {0,2}:
            bad+=1
            if bad<10: print("J* not subset {0,2}:",(a,b,m),J)
        # tie prediction
        tie_pred = (a%4,b%4) in {(2,2),(1,3)}
        tie_act = (J==[0,2])
        if tie_pred!=tie_act:
            bad+=1
            if bad<10: print("tie mismatch",(a,b,m),J,"pred",tie_pred)
        if tie_act: ntie+=1
print(f"tested={tested} ties={ntie} bad={bad}  (MN, m<=48)")
