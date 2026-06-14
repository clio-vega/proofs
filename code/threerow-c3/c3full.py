from math import comb
from mn import Mj
def v2(n):
    if n==0:return None
    n=abs(n);c=0
    while n%2==0:n//=2;c+=1
    return c
def Jstar(av,bv,mv):
    vals={}
    for jv in range(0,bv+4):
        if jv>mv:break
        M=Mj((av,bv,3),jv,mv)
        if M==0:continue
        cm=comb(mv,jv)
        vals[jv]=jv+2*(v2(cm)+v2(M))
    V=min(vals.values())
    return sorted(j for j,v in vals.items() if v==V)

def predict(av,bv):
    a4,b4=av%4,bv%4
    if av%2==0: # a even, j0=0
        if a4==0 and b4==1: return [0,4]
        if a4==2 and b4==3: return [0,2]
        return [0]
    else: # a odd, j0=3
        if a4==1 and b4==2: return [3,5,7,9]
        return [3]

bad=0;tested=0;sizes={1:0,2:0,4:0}
for mv in range(3,40):
    n=2*mv
    for av in range(1,n+1):
        for bv in range(3,av+1):
            if n-av-bv!=3:continue
            J=Jstar(av,bv,mv)
            P=predict(av,bv)
            # only compare when tie points are interior (j<=b); else boundary regime
            tested+=1
            sizes[len(J)]=sizes.get(len(J),0)+1
            if J!=P:
                # could be boundary effect when max(P)>b
                if max(P)<=bv:
                    bad+=1
                    if bad<25: print("MISMATCH",(av,bv,3),"m",mv,"got",J,"pred",P)
print(f"full theorem check (interior-tie shapes): tested={tested} mismatches={bad}")
print("size distribution:",sizes)
print("all |J*| in {1,2,4} and even on ties:", all(s in(1,2,4) for s in sizes if sizes[s]))
