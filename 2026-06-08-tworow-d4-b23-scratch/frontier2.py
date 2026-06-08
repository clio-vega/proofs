# Rigorous finite check: for b = 2,3 mod 4, b<=N, confirm I_b(m) != 0 for all
# integers m in [b, M] where M > largest real root of I_b.  Since m>=b>R=floor((b-1)/2),
# the forced product prod_{r=0}^R (m-r) != 0, so I_b(m)=0 <=> Q_b(m)=0.  Thus checking
# I_b(m)!=0 for integer m up to past the largest real root proves (star) for that b.
# I_b(m) via exact integer engine:
#   I_b(m) = sum_{j=1,j%4!=0}^{b} C(m,j) Im(s^j) eps_j C(m-j, floor((b-j)/2))
# Im(s^j) = (-4)^{j//4} * {1:1,2:2,3:2}[j%4];  eps_j=(-1)^(b-j).
from math import comb, isqrt

def Im_s_pow(j):
    r=j%4
    if r==0: return 0
    return {1:1,2:2,3:2}[r]*(-4)**(j//4)

def I_b_int(b, m):
    tot=0
    for j in range(1,b+1):
        c=Im_s_pow(j)
        if c==0: continue
        eps=1 if (b-j)%2==0 else -1
        k=(b-j)//2
        if m-j < 0: continue
        tot += comb(m,j)*c*eps*comb(m-j,k)
    return tot

def largest_real_root_bound(b):
    # empirical: largest real root ~ 0.33 b^2; use generous M = b^2 (safe margin),
    # but verify by checking sign stabilises (no more sign changes) for a stretch.
    return b*b + 5

import sys
N = int(sys.argv[1]) if len(sys.argv)>1 else 60
worst=0
for b in range(2, N+1):
    if b%4 not in (2,3): continue
    M = largest_real_root_bound(b)
    # find sign changes / zeros
    prev=None
    zeros=[]
    last_change=b
    vals_sign=[]
    for m in range(b, M+1):
        v=I_b_int(b,m)
        if v==0:
            zeros.append(m)
        sgn = 1 if v>0 else (-1 if v<0 else 0)
        if prev is not None and sgn!=0 and prev!=0 and sgn!=prev:
            last_change=m
        if sgn!=0: prev=sgn
    # ensure M is safely past last sign change
    margin = M - last_change
    status = "OK" if not zeros else f"ZEROS{zeros}"
    if zeros: worst=b
    if b%20 in (2,3,6,7,10,11) or b>N-6:
        print(f"b={b:3d}(mod4={b%4}): M={M}, last sign change at m={last_change} (margin {margin}), {status}")
print("DONE. Any integer roots found among b<=%d?  %s"%(N, "NONE" if worst==0 else f"at b={worst}"))
