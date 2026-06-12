def s2(n): return bin(n).count('1')
# Lemma A: 2 s2(j) <= j+1 for all j>=0, with equality iff j in {1,3}.
# stronger: 2 s2(j) - j: =0 iff j in {0,2}; =1 iff j in{1,3}; <=-1 for j>=4.
bad=[]
for j in range(0,200000):
    d=2*s2(j)-j
    if d>1: bad.append(("d>1",j,d))
    if d==1 and j not in (1,3): bad.append(("d=1 extra",j))
    if d==0 and j not in (0,2): bad.append(("d=0 extra",j))
    if j>=4 and d>=0: bad.append(("j>=4 nonneg",j,d))
print("Lemma A (2s2(j)<=j+1, equality j in {1,3}; j>=4 => 2s2-j<=-1):", "OK" if not bad else bad[:10])

# Lemma C(n,3) odd iff n==3 mod4 ; C(n,2) odd iff n in {2,3} mod 4
from math import comb as C
b3=all(((C(n,3)%2==1)==(n%4==3)) for n in range(3,5000))
b2=all(((C(n,2)%2==1)==(n%4 in (2,3))) for n in range(2,5000))
print("C(n,3) odd <=> n=3 mod4:", b3)
print("C(n,2) odd <=> n=2,3 mod4:", b2)
