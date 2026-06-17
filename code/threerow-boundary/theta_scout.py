"""theta_scout.py -- interior offset theta = val(0) - V and j0=min J*, by c and parity of a."""
from math import comb
from mn import Mj

def v2(n):
    n=abs(int(n))
    if n==0: return 999
    k=0
    while n%2==0: n//=2; k+=1
    return k

def val(lam,j,m):
    if j>m: return None
    try:
        M=Mj(lam,j,m)
    except AssertionError:
        return None
    if M==0: return None
    return j+2*v2(comb(m,j)*M)

for c in range(2,8):
    data={}  # (a%2) -> set of (theta, j0)
    for m in range(c+4,34):
        for bb in range(2*c, m):   # box-interior gate b>=2c
            aa=2*m-c-bb
            if aa<bb: continue
            lam=(aa,bb,c)
            vals={}
            for j in range(0,bb+c+1):
                v=val(lam,j,m)
                if v is not None: vals[j]=v
            if 0 not in vals: continue
            V=min(vals.values())
            Jstar=sorted(j for j,v in vals.items() if v==V)
            j0=Jstar[0]
            theta=vals[0]-V
            data.setdefault(aa%2,set()).add((theta,j0))
    print("c=%d:"%c)
    for par in sorted(data):
        ths=sorted(data[par])
        print("   a%%2=%d : (theta,j0) seen = %s"%(par, ths[:12]))
    print("   c(c-1)/2=%d  c(c+1)/2=%d"%(c*(c-1)//2, c*(c+1)//2))
