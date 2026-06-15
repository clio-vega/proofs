"""Test the split P = P1 - P2 against the target.
P(j)=C(a+4,j)C(b+3,j)Q3(j), Q3=(a-1)(b-2)H - 720 C(j,6).
P1=(a-1)(b-2)C(a+4,j)C(b+3,j)H(j)   (heavy)
P2=720 C(a+4,j)C(b+3,j)C(j,6)        (tip)
Need: 2 v2 P_i(j) >= 2 v2P(0) + 2 s2(j) - j - 3   for 4<=j<=b   (=> tildeDelta>=0).
"""
from math import comb
import collections
def Hf(av,bv,jv):
    return ((av+3)*(av+4)*(bv+2)*(bv+3) -6*(av+3)*(bv+2)*comb(jv,1)
        -6*(av*bv+av+2*bv)*comb(jv,2)+36*comb(jv,3)+72*comb(jv,4))
def v2(n):
    if n==0:return 10**9
    n=abs(n);c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n):return bin(int(n)).count('1')

f1=f2=fboth=0
worst1={}; worst2={}
for mv in range(4,80):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            vP0 = v2((av-1)*(bv-2)*Hf(av,bv,0))
            target_x2 = 2*vP0 + 2*s2_cache if False else None
            for jv in range(4,bv+1):
                Ca=comb(av+4,jv); Cb=comb(bv+3,jv)
                P1 = (av-1)*(bv-2)*Ca*Cb*Hf(av,bv,jv)
                P2 = 720*Ca*Cb*comb(jv,6)
                tx2 = 2*vP0 + 2*s2(jv) - jv - 3   # need 2 v2 Pi >= tx2
                d1 = 2*v2(P1)-tx2
                d2 = 2*v2(P2)-tx2
                if d1<0: f1+=1
                if d2<0: f2+=1
                if d1<0 and d2<0: fboth+=1
                worst1[v2(jv)]=min(worst1.get(v2(jv),10**9),d1)
                worst2[v2(jv)]=min(worst2.get(v2(jv),10**9),d2)
print(f"P1 (heavy) misses target: {f1}")
print(f"P2 (tip)   misses target: {f2}")
print(f"BOTH miss simultaneously: {fboth}")
print("worst margin d1 (2v2P1-target) by v2(j):",dict(sorted(worst1.items())))
print("worst margin d2 (2v2P2-target) by v2(j):",dict(sorted(worst2.items())))
