import sys; sys.setrecursionlimit(100000)
from mn import Mj, fdim
from math import comb
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
# c=3: for b large enough that box interior, check boundary adds no minimizer
# a even -> box {0,2,4,6}? need b>=6 ; a odd -> box {3,5,7,9}, need b>=9? offset3+6=9<=b => b>=9
viol=0;t=0;mins=[]
for m in range(4,40):
    for b in range(3,m):
        a=2*m-3-b
        if a<b: continue
        bmin = 6 if a%2==0 else 9
        if b<bmin: continue
        lam=(a,b,3)
        vals={}
        for j in range(0,b+4):
            if j>m: continue
            M=Mj(lam,j,m)
            if M!=0: vals[j]=j+2*v2(comb(m,j)*M)
        if not vals: continue
        V=min(vals.values())
        for i in (1,2,3):
            j=b+i
            if j in vals:
                marg=vals[j]-V
                if marg<=0: viol+=1; print("c3 VIOL",lam,m,j,vals[j],V)
                t+=1
print(f"c=3 boundary (box interior): tested={t} viol={viol}")
