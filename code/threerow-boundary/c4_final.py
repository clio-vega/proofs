from mn import Mj
from math import comb
def v2(n):
    n=abs(n)
    if n==0:return 99
    k=0
    while n%2==0:n//=2;k+=1
    return k
# end-to-end: Delta(b+i)>0 for all c=4 shapes b>=6, plus record exact min margins per (parity,i)
from collections import defaultdict
mins=defaultdict(lambda:99); viol=0; tested=0
for m in range(8,92):
    for b in range(6,m):
        a=2*m-4-b
        if a<b:continue
        lam=(a,b,4)
        M0=Mj(lam,0,m)
        if M0==0:continue
        v0=2*v2(M0)  # comb(m,0)=1
        for i in [1,2,3,4]:
            if b+i>m:continue
            Mi=Mj(lam,b+i,m)
            if Mi==0:continue
            D=(b+i)+2*v2(comb(m,b+i)*Mi)-v0
            tested+=1
            if D<=0:viol+=1
            key=('odd' if b%2 else 'even',i)
            mins[key]=min(mins[key],D)
print("end-to-end Delta(b+i)>0, b>=6, m<92: tested=%d violations=%d"%(tested,viol))
print("exact min Delta per (b-parity,i):")
for k in sorted(mins):print("   ",k,mins[k])
