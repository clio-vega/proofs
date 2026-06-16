from mn import Mj
from math import comb
def v2(n):
    n=abs(n); k=0
    if n==0: return 99
    while n%2==0: n//=2; k+=1
    return k
def val(lam,j,m):
    M=Mj(lam,j,m)
    if M==0: return None
    return j+2*v2(comb(m,j)*M)
viol=[]; cnt=0
for m in range(8,61):
    for b in range(3,m):
        a=2*m-3-b
        if a<b: continue
        c=3; P=m-b-3
        if P<0: continue
        lam=(a,b,c)
        # box-interior range only (per note): a even b>=6 ; a odd b>=9
        if a%2==0 and b<6: continue
        if a%2==1 and b<9: continue
        theta = 0 if a%2==0 else 3
        v0=val(lam,0,m)
        if v0 is None: continue
        for i in (1,2,3):
            j=b+i
            if j>m: continue
            vj=val(lam,j,m)
            if vj is None: continue
            cnt+=1
            if not (vj-v0 > -theta):
                viol.append((lam,m,i,vj-v0,theta))
print(f"checked {cnt} boundary indices, m<=60. violations:", viol[:10], "count", len(viol))
