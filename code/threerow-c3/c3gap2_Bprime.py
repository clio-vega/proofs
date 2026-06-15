"""B'  (a odd):  v2(a+1) + v2C(a-2,j-6) >= v2C(j,6)+s2(j)-(j+3)/2+3 =: R(j),  6<=j<=a.
Test B', and compute max R(j); also test the even-cleaner claim R(j)<=1 (=> B' trivial since LHS>=1).
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

# R(j) is independent of a. Compute (x2 to be integer): 2R = 2 v2C(j,6)+2 s2(j)-(j+3)+6
maxR2=-10**9; argmax=None
for jv in range(6,4000):
    R2=2*v2C(jv,6)+2*s2(jv)-(jv+3)+6
    if R2>maxR2: maxR2=R2; argmax=jv
print(f"max 2R(j) over 6<=j<4000 = {maxR2} at j={argmax}  (R={maxR2/2})")
# show first several R
print("2R(j) j=6..20:", [2*v2C(j,6)+2*s2(j)-(j+3)+6 for j in range(6,21)])

# direct B' test over odd a, 6<=j<=a
fB=0; tight=0
for av in range(5,1500,2):
    for jv in range(6,av+1):
        L2=2*(v2(av+1)+v2C(av-2,jv-6))
        R2=2*v2C(jv,6)+2*s2(jv)-(jv+3)+6
        if L2<R2: fB+=1
        if L2==R2: tight+=1
print(f"B' fails={fB}  tight={tight}  (odd a<1500)")
