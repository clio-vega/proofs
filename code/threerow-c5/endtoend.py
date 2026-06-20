import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from mn import Mj
from math import comb
def v2(n):
    n=abs(int(n))
    if n==0: return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def Jstar_interior(a,b,m):
    vals={}
    for jj in range(0,b+1):
        M=Mj((a,b,5),jj,m)
        if M!=0: vals[jj]=jj+2*v2(comb(m,jj)*M)
    V=min(vals.values())
    return sorted(jj for jj,v in vals.items() if v==V), V
def predict(a,b):
    if a%2==0: # b odd
        return [0,2] if (a%4,b%4)==(0,1) else [0]
    else: # b even
        return [3,5] if (a%4,b%4)==(3,0) else [3]
mism=0; tested=0; ties=0
for m in range(5,46):
    n=2*m
    for a in range(5,n):
        b=n-a-5
        if b<10 or b>a: continue   # box-interior
        J,V=Jstar_interior(a,b,m)
        pred=predict(a,b)
        tested+=1
        if len(J)>1: ties+=1
        if J!=pred:
            mism+=1
            if mism<10: print("MISMATCH",(a,b,m),"got",J,"pred",pred)
print(f"END-TO-END MN: tested={tested} ties={ties} mismatches={mism} (m<=45, box-interior b>=10)")
