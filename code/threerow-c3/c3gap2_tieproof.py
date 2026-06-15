"""Verify the HAND tie-analysis formulas and case splits for j=5,7,9.
Hand claims:
 U(5)=v2(a+1)+v2(b)+v2H(5)-6,        v2H(5)=2 unless a3p+b2p=2 (then >=2)
 U(7)=v2(a+1)+v2(b)+v2Q3(7)-8,       v2H(7)=1 always; v2Q3(7) per case
 U(9)=v2(a+1)+v2(a-3)+v2(b)+v2(b-4)+v2Q3(9)-14
 tie(j) <=> Uj = g(j):  g(5)=g(7)=-2, g(9)=-4
Also: v2C(a+4,5)=a3p+v2(a+1)-3 ; v2C(b+3,5)=b2p+v2(b)-3 ; analogues for 7,9.
"""
from math import comb
def v2(n):
    if n==0:return 10**9
    n=abs(n);c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n):return bin(int(n)).count('1')
def Hf(av,bv,jv):
    return ((av+3)*(av+4)*(bv+2)*(bv+3)-6*(av+3)*(bv+2)*comb(jv,1)
        -6*(av*bv+av+2*bv)*comb(jv,2)+36*comb(jv,3)+72*comb(jv,4))
def Qf(av,bv,jv): return (av-1)*(bv-2)*Hf(av,bv,jv)-720*comb(jv,6)
def Ureal(av,bv,jv):
    return v2(comb(av+4,jv))+v2(comb(bv+3,jv))+v2(Qf(av,bv,jv))-v2(Qf(av,bv,0))

err=dict.fromkeys(["vC5a","vC5b","vC7a","vC7b","vC9a","vC9b","U5","U7","U9",
                   "vH7=1","vH5","vH9","U5form","U7form","U9form"],0)
for mv in range(7,80):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            a1,a3,b2 = v2(av-1),v2(av+3),v2(bv+2)
            # binomial valuation formulas
            if comb(av+4,5)!=0 and v2(comb(av+4,5))!=a3+v2(av+1)-3: err["vC5a"]+=1
            if v2(comb(bv+3,5))!=b2+v2(bv)-3: err["vC5b"]+=1
            if v2(comb(av+4,7))!=a1+a3+v2(av+1)-4: err["vC7a"]+=1
            if v2(comb(bv+3,7))!=v2(bv-2)+b2+v2(bv)-4: err["vC7b"]+=1
            if v2(comb(av+4,9))!=a1+a3+v2(av+1)+v2(av-3)-7: err["vC9a"]+=1
            if v2(comb(bv+3,9))!=v2(bv-2)+b2+v2(bv)+v2(bv-4)-7: err["vC9b"]+=1
            # vH values
            if v2(Hf(av,bv,7))!=1: err["vH7=1"]+=1
            if a3+b2>2 and v2(Hf(av,bv,5))!=2: err["vH5"]+=1
            if a3+b2>3 and v2(Hf(av,bv,9))!=3: err["vH9"]+=1
            # U formulas
            if bv>=5 and Ureal(av,bv,5)!=v2(av+1)+v2(bv)+v2(Hf(av,bv,5))-6: err["U5form"]+=1
            if bv>=7 and Ureal(av,bv,7)!=v2(av+1)+v2(bv)+v2(Qf(av,bv,7))-8: err["U7form"]+=1
            if bv>=9 and Ureal(av,bv,9)!=v2(av+1)+v2(av-3)+v2(bv)+v2(bv-4)+v2(Qf(av,bv,9))-14: err["U9form"]+=1
print({k:v for k,v in err.items()})

# margin for odd j>=11 (strictness off the box): min (U(j)-g(j)) over shapes
import collections
mar=collections.defaultdict(lambda:10**9)
for mv in range(7,80):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            for jv in range(11,bv+1,2):  # odd j>=11
                g2=2*s2(jv)-(jv+3)
                mar[jv]=min(mar[jv], 2*Ureal(av,bv,jv)-g2)  # = tildeDelta(j)
print("min tildeDelta for odd j>=11 (j: minval):", dict(sorted(mar.items())))
