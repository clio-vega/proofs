import sys; sys.setrecursionlimit(100000)
from mn import Mj, fdim
from math import comb
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def s2(n):return bin(int(n)).count('1')
# c=1: verify the FULL boundary lemma proof conclusion: for b not in {1,2,4}, val(b+1)>V
viol=0;t=0
for m in range(2,50):
    for b in range(1,m):
        a=2*m-1-b
        if a<b: continue
        if b in (1,2,4): continue
        lam=(a,b,1)
        vals=[]
        for j in range(0,b+2):
            if j>m: continue
            M=Mj(lam,j,m)
            if M==0: continue
            vals.append(j+2*v2(comb(m,j)*M))
        V=min(vals)
        j=b+1
        if j<=m:
            M=Mj(lam,j,m)
            if M!=0:
                vb1=j+2*v2(comb(m,j)*M)
                if vb1<=V: viol+=1; print("c1 VIOL",lam,m,vb1,V)
                t+=1
print(f"c=1 boundary (b not in 1,2,4): tested={t} viol={viol}")
