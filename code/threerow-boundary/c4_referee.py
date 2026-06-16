# Referee check: the boundary lemma needs val(b+i) > V (global min), NOT just > val(0).
# We prove val(b+i) > val(0); since val(0) >= V trivially, val(b+i) > V follows.
# Confirm val(0) IS the global min (j0=0, theta=0) so nothing is lost, over wide range.
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
    return None if M==0 else j+2*v2(comb(m,j)*M)
bad=0
for m in range(8,80):
    for b in range(6,m):
        a=2*m-4-b
        if a<b:continue
        lam=(a,b,4)
        vs={j:val(lam,j,m) for j in range(0,min(b+5,m+1)) if val(lam,j,m) is not None}
        if 0 not in vs:continue
        V=min(vs.values())
        # check: val(0)==V (j0=0,theta=0) AND every boundary index > val(0)
        if vs[0]!=V: bad+=1
        for i in [1,2,3,4]:
            if b+i in vs and vs[b+i]<=vs[0]: bad+=1
print("referee: val(0)==V and val(b+i)>val(0) for all b>=6, m<80 -- violations:",bad)
