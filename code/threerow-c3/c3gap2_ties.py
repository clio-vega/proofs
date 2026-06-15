"""Tie classification. tildeDelta(j)=0 claimed iff j in {5,7,9} and a=1,b=2 mod4 (all-or-nothing).
Check: (1) for j not in {5,7,9}, is td>0 always?  (2) does td(5),td(7),td(9) depend only on (a%4,b%4)?
       (3) tabulate td(5),td(7),td(9) by (a%4,b%4).  (4) Also dependence on higher mod for safety.
"""
from math import comb
import collections
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
def td(av,bv,jv):
    vP0=v2(Qf(av,bv,0))
    U=v2(comb(av+4,jv)*comb(bv+3,jv)*Qf(av,bv,jv))-vP0
    return jv+3-2*s2(jv)+2*U

# (1) td>0 for j not in {5,7,9}, j>=4 interior
viol=0
for mv in range(7,60):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            for jv in range(4,bv+1):
                if jv in (5,7,9):continue
                if td(av,bv,jv)==0: viol+=1; print("UNEXPECTED tie at",av,bv,jv)
print("non-{5,7,9} ties:",viol)

# (2)+(3) tabulate td(j) by (a%4,b%4) and check determinism by that pair
for jv in (5,7,9):
    tab=collections.defaultdict(set)
    for mv in range(9,80):
        n=2*mv
        for av in range(1,n+1,2):
            for bv in range(max(jv,9-0),av+1,2):  # ensure interior (b>=9 for j=9)
                if n-av-bv!=3:continue
                if bv<jv: continue
                tab[(av%4,bv%4)].add(td(av,bv,jv))
    print(f"\nj={jv}: td values by (a%4,b%4):")
    for key in sorted(tab): print(f"   a%4={key[0]} b%4={key[1]}: td in {sorted(tab[key])}")
