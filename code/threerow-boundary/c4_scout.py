import sympy as sp
from math import comb, factorial
from mn import Mj

def v2(n):
    n=abs(n)
    if n==0: return 99
    k=0
    while n%2==0: n//=2; k+=1
    return k

def val(lam,j,m):
    M=Mj(lam,j,m)
    if M==0: return None
    return j+2*v2(comb(m,j)*M)

# ---- 1. Determine interior structure: j0, V, val(0), theta, by parity/mod4 ----
print("="*70)
print("c=4 interior structure: theta=val(0)-V by parity of a, (a%4,b%4)")
print("="*70)
c=4
from collections import defaultdict
table=defaultdict(lambda: defaultdict(int))
theta_by_apar=defaultdict(set)
j0_by_apar=defaultdict(set)
for m in range(8,46):
    for b in range(c,m):
        a=2*m-c-b
        if a<b: continue
        P=m-b-c
        if P<0: continue
        if b < 2*c: continue   # crude box-interior gate
        lam=(a,b,c)
        vals={}
        for j in range(0,b+c+1):
            vj=val(lam,j,m)
            if vj is not None:
                vals[j]=vj
        if 0 not in vals: continue
        V=min(vals.values())
        Jstar=[j for j in vals if vals[j]==V]
        j0=min(Jstar)
        theta=vals[0]-V
        apar = a%2
        theta_by_apar[apar].add(theta)
        j0_by_apar[apar].add(j0)
        table[(a%4,b%4)][(j0,theta)] += 1

print("theta values by parity of a:", {k:sorted(v) for k,v in theta_by_apar.items()})
print("j0 values by parity of a:", {k:sorted(v) for k,v in j0_by_apar.items()})
print()
print("(a%4,b%4) -> {(j0,theta):count}")
for k in sorted(table):
    print("  ",k, dict(table[k]))
