from mn import Mj
from math import comb
def v2(n):
    n=abs(n)
    if n==0: return 99
    k=0
    while n%2==0: n//=2; k+=1
    return k
def valMN(lam,j,m):
    M=Mj(lam,j,m)
    if M==0: return None
    return j+2*v2(comb(m,j)*M)

# per-index, per-parity minimum Delta over box-interior range
# box-interior gate: interior box {0,2,4,6} must lie in [0,b], i.e. b>=6.  Use b>=6.
from collections import defaultdict
mins=defaultdict(lambda:(999,None))
viol=0
for m in range(8,70):
    for b in range(4,m):
        a=2*m-4-b
        if a<b: continue
        P=m-b-4
        if P<0: continue
        if b<6: continue
        lam=(a,b,4)
        v0=valMN(lam,0,m)
        if v0 is None: continue
        apar = 'aeven' if a%2==0 else 'aodd'
        for i in [1,2,3,4]:
            j=b+i
            if j>m: continue
            vj=valMN(lam,j,m)
            if vj is None: continue
            D=vj-v0
            if D<=0: viol+=1
            key=(apar,i)
            if D<mins[key][0]:
                mins[key]=(D,(a,b,m))
print("violations (Delta<=0):",viol)
print("per (parity,i) minimum Delta and extremal shape:")
for k in sorted(mins):
    print("  ",k, mins[k])
