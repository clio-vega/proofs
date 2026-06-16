from mn import Mj
from math import comb
def v2(n):
    n=abs(n)
    if n==0:return 99
    k=0
    while n%2==0:n//=2;k+=1
    return k
def val(lam,j,m):
    M=Mj(lam,j,m)
    if M==0:return None
    return j+2*v2(comb(m,j)*M)
# full claim: for c=4, J* subset {0,2,4,6}, |J*| in {1,2,4} (even-cardinality), j0=0
bad=[]
sizes={}
for m in range(6,60):
    for b in range(4,m):
        a=2*m-4-b
        if a<b:continue
        lam=(a,b,4)
        vals={}
        for j in range(0,min(b+5,m+1)):
            vj=val(lam,j,m)
            if vj is not None:vals[j]=vj
        if not vals:continue
        V=min(vals.values())
        Js=sorted(j for j in vals if vals[j]==V)
        sizes[len(Js)]=sizes.get(len(Js),0)+1
        if any(j not in (0,2,4,6) for j in Js): bad.append(('not-in-box',a,b,m,Js))
        if len(Js)%2==0 and 0 not in Js: pass
        # |J*| even check: |J*| in {1,2,4}? actually parity of count: should be that J* is affine box -> size power of 2
        if len(Js) not in (1,2,4): bad.append(('badsize',a,b,m,Js))
print("full-theorem violations (m<60, incl small b):",len(bad),bad[:10])
print("|J*| distribution:",sizes)

# product-bound lemmas sanity
import itertools
# E: v2(prod of Q consec) >= floor(Q/2); remove 1 factor: >= floor(Q/2)-1
ef=0
for Q in range(2,40):
    for R in range(0,300):
        prod_v=sum(v2(R+t) for t in range(1,Q+1))
        if prod_v < Q//2: ef+=1
        for drop in range(1,Q+1):
            rem=prod_v - v2(R+drop)
            if rem < Q//2 -1: ef+=1
print("product-bound (evens-count) failures:",ef)
