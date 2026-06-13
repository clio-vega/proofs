import sys; sys.path.insert(0,'/home/clio/projects/code/threerow-c1')
from math import comb
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def s2(n): return bin(n).count('1')
def vC(n,k):
    if k<0 or k>n:return 10**9
    return v2(comb(n,k))
def Q(a,b,j): return a*(b-1)*((a+3)*(b+2)-2*j*j)+j*(j-1)*(j-2)*(j-3)
def Delta(a,b,j):
    return j-2*s2(j)+2*vC(a+3,j)+2*vC(b+2,j)+2*(v2(Q(a,b,j))-v2(Q(a,b,0)))

# 1) confirm: for j not in {0,2,4}, Delta>0 always (interior 0<=j<=b)
minpos={}
viol=0
for A in range(2,80):
    for B in range(2,A+1):
        if (A+B)%2: continue
        for j in range(1,B+1):
            d=Delta(A,B,j)
            if j not in (2,4):
                if d<=0:
                    viol+=1
                    if viol<10: print("VIOL j!in{0,2,4}",A,B,j,d)
print("interior j∉{0,2,4} with Delta<=0:",viol)

# 2) Delta(2),Delta(4) as fn of (a%8,b%8) -- find governing modulus
from collections import defaultdict
d2=defaultdict(set); d4=defaultdict(set)
for A in range(4,200):
    for B in range(4,A+1):
        if (A+B)%2: continue
        if B>=3: d2[(A%8,B%8)].add(Delta(A,B,2))
        if B>=5: d4[(A%8,B%8)].add(Delta(A,B,4))
# is Delta(2) determined by (a,b) mod 8? mod 4?
def determined(dd, mod, src):
    m=defaultdict(set)
    for A in range(4,200):
        for B in range(4,A+1):
            if (A+B)%2:continue
            if src=='d2' and B<3: continue
            if src=='d4' and B<5: continue
            val=Delta(A,B, 2 if src=='d2' else 4)
            m[(A%mod,B%mod)].add(val)
    return all(len(v)==1 for v in m.values()), m
for mod in [2,4,8,16]:
    ok2,_=determined(None,mod,'d2')
    ok4,_=determined(None,mod,'d4')
    print(f"mod {mod}: Delta(2) determined={ok2}, Delta(4) determined={ok4}")
