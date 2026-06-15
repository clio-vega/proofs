"""Verify the full TIP proof chain end-to-end, and the heart (heart): 2 v2C(j,6)+2 s2(j) <= j-1.
Also verify v2C(j,6)=2+s2(j-6)-s2(j) and reduce heart to 2 s2(j-6) <= (j-6)+1.
"""
from math import comb
def v2(n):
    if n==0:return 10**9
    n=abs(n);c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n):return bin(int(n)).count('1')
def v2C(N,k):
    if k<0 or k>N:return 10**9
    return s2(k)+s2(N-k)-s2(N)

# (heart): 2 v2C(j,6)+2 s2(j) <= j-1 for j>=6
bad=0
for j in range(6,20000):
    if 2*v2C(j,6)+2*s2(j) > j-1: bad+=1
print("heart 2v2C(j,6)+2s2(j)<=j-1 fails:",bad)
# v2C(j,6)=2+s2(j-6)-s2(j)
bad2=sum(1 for j in range(6,5000) if v2C(j,6)!=2+s2(j-6)-s2(j))
print("v2C(j,6)=2+s2(j-6)-s2(j) fails:",bad2)
# Lemma A: 2 s2(k) <= k+1
badA=sum(1 for k in range(0,20000) if 2*s2(k)>k+1)
print("Lemma A 2s2(k)<=k+1 fails:",badA)

# Full end-to-end tip check: 2 v2 P2 >= 2 v2P(0)+2 s2(j)-(j+3), all a odd,b even,a>=b>=4,a+b=2m-3, 6<=j<=b
def Hf(av,bv,jv):
    return ((av+3)*(av+4)*(bv+2)*(bv+3) -6*(av+3)*(bv+2)*comb(jv,1)
        -6*(av*bv+av+2*bv)*comb(jv,2)+36*comb(jv,3)+72*comb(jv,4))
ftip=0
for mv in range(4,80):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            vP0=v2((av-1)*(bv-2)*Hf(av,bv,0))
            for jv in range(6,bv+1):
                P2=720*comb(av+4,jv)*comb(bv+3,jv)*comb(jv,6)
                if 2*v2(P2) < 2*vP0+2*s2(jv)-(jv+3): ftip+=1
print("end-to-end TIP (star2) fails:",ftip)
