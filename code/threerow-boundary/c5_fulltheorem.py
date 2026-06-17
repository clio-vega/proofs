"""c5_fulltheorem.py -- verify NO boundary index is a minimizer of val, ALL c=5 shapes."""
from math import comb
from mn import Mj
def v2(n):
    n=abs(int(n))
    if n==0: return 999
    k=0
    while n%2==0: n//=2;k+=1
    return k
c=5; bad=0; checked=0; minJ=set()
for m in range(c+1, 50):
    for bb in range(c, m+1):
        aa=2*m-c-bb
        if aa<bb: continue
        lam=(aa,bb,c)
        vals={}
        for j in range(0, bb+c+1):
            if j>m: break
            try: M=Mj(lam,j,m)
            except AssertionError: continue
            if M==0: continue
            vals[j]=j+2*v2(comb(m,j)*M)
        if not vals: continue
        V=min(vals.values())
        Jstar=sorted(j for j,v in vals.items() if v==V)
        minJ.add(len(Jstar))
        # boundary indices b+1..b+c must NOT be in Jstar
        for i in range(1,c+1):
            j=bb+i
            if j in vals:
                checked+=1
                if vals[j]==V:
                    bad+=1
                    if bad<=8: print("  BOUNDARY TIE: lam=%s m=%d j=%d (i=%d)"%(lam,m,j,i))
print("c=5 ALL shapes m<50: boundary indices checked=%d  boundary-minimizers=%d"%(checked,bad))
print("|J*| values seen:", sorted(minJ))
