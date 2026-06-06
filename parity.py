from math import comb
from functools import lru_cache
@lru_cache(maxsize=None)
def trow(N):
    cur=(1,)
    for _ in range(N):
        nxt=[0]*(len(cur)+2)
        for j,v in enumerate(cur): nxt[j]+=v;nxt[j+1]+=v;nxt[j+2]+=v
        cur=tuple(nxt)
    return cur
def Tr(N,j):
    r=trow(N); return r[j] if 0<=j<=2*N else 0
def ImG(M,b):
    s=0
    for k in range(1,b+1,2):
        s+=((-1)**((k-1)//2))*comb(M,k)*(Tr(M-k,b-k)-Tr(M-k,b-k-1))
    return s
# parity pattern of ImG over (m,b)
print("ImG mod 2, rows m=3..14, cols b=1..m:")
for m in range(3,15):
    print(f" m={m:2d}: ", "".join(str(ImG(m,b)&1) for b in range(1,m+1)))
# When is ImG odd? count fraction
odd=0;tot=0
for m in range(3,120):
    for b in range(1,m+1):
        tot+=1; 
        if ImG(m,b)&1: odd+=1
print(f"fraction odd among m=3..119: {odd}/{tot}")
# check: ImG odd  <=> ? maybe depends on m,b mod 4. Tabulate parity vs (m%4,b%4)
from collections import Counter
patt=Counter()
for m in range(3,200):
    for b in range(1,m+1):
        patt[(m%4,b%4, ImG(m,b)&1)]+=1
# is parity determined by (m%4,b%4)? check no (a,b,0) and (a,b,1) both occur
import itertools
amb=set((x,y) for x,y,_ in patt)
det=True
for (x,y) in amb:
    if patt[(x,y,0)]>0 and patt[(x,y,1)]>0: det=False
print("parity determined by (m mod4, b mod4)?", det)
