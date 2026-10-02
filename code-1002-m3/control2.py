"""Refined control panel for (C), with all four outcomes PREDICTED IN ADVANCE from the
telescoping

    sum_{i=1}^m (L_i + R_i) = sum_{i=1}^m (L_i + R_{i-1})  +  R_m - R_0
                            = sum_i (nu_i + lam_{i-1})     +  (R_m - R_0).

Two things must happen: (a) each LOCAL pair L_i + R_{i-1} collapses to nu_i + lam_{i-1}
  -- needs the two clamp OFFSETS to match;
(b) the single UNMATCHED boundary term R_m - R_0 must be nu-free
  -- needs the lam-period and the nu-period to be the SAME n.

Predictions:
  P1  p=q=k, any k            -> HOLDS   (offset VALUE is free)
  P2  p != q                  -> FAILS   (local collapse destroyed)
  P3  lam-period n' != nu-period n -> FAILS (boundary term becomes nu-dependent)
  P4  no wrap at all (ordinary skew: L_1=nu_1, R_m=lam_m) -> HOLDS
          (both boundary terms degenerate separately: nu_1 is in S_nu, lam_m is a constant)
"""
import gen, sys
from collections import Counter

def LRgen(nu, lam, n, m, p=1, q=1, nprime=None, linear=False):
    nprime = n if nprime is None else nprime
    out=[]
    for i in range(m):
        if i==0:
            lp = lam[m-1]-nprime          # lam_0 = lam_m - n'
            Lo = nu[0] if linear else max(nu[0], lp+p)
        else:
            Lo = max(nu[i], lam[i-1]+p)
        if i==m-1:
            Ro = lam[m-1] if linear else min(nu[0]+n-q, lam[m-1])
        else:
            Ro = min(nu[i+1]-q, lam[i])
        out.append((Lo,Ro))
    return out

def test(m,nmax,dmax,**kw):
    st=Counter(); wit=[]
    for n in range(m,nmax+1):
        for (mu,lam) in gen.pairs(n,m,dmax):
            for b in range(0,sum(lam)-sum(mu)+1):
                nus=gen.slice_nus(mu,lam,n,m,b)
                if len(nus)<2: continue
                vals={sum(L+R for (L,R) in LRgen(nu,lam,n,m,**kw)) for nu in nus}
                st['slices']+=1
                if len(vals)==1: st['const']+=1
                else:
                    st['nonconst']+=1
                    if len(wit)<1: wit.append((n,mu,lam,b,sorted(vals)))
    return st,wit

panel=[(dict(),'P0  cylindric, p=q=1, n\'=n'),
       (dict(p=2,q=2),'P1  p=q=2        (offset value changed, match kept)'),
       (dict(p=4,q=4),'P1  p=q=4        (offset value changed, match kept)'),
       (dict(p=2,q=1),'P2  p=2,q=1      (OFFSET MISMATCH)'),
       (dict(p=1,q=3),'P2  p=1,q=3      (OFFSET MISMATCH)'),
       (dict(nprime=None),'    (placeholder)'),
       (dict(linear=True),'P4  ordinary skew, no wrap (L_1=nu_1, R_m=lam_m)')]
for m in (2,3,4):
    print(f"--- m={m}, n<={m+4}, d<=7 ---")
    for kw,label in panel:
        if label.startswith('    '): 
            # P3 handled separately: n' = n+1 and n' = n-1
            for dn in (1,-1):
                res=[]
                for n in range(m,m+5):
                    st,wit=test_one=None,None
                res=None
            continue
        st,wit=test(m,m+4,7,**kw)
        v="HOLDS" if st['nonconst']==0 else f"FAILS {st['nonconst']}/{st['slices']}"
        print(f"   {label:52s} -> {v}")
        if wit: print(f"        e.g. n={wit[0][0]} mu={wit[0][1]} lam={wit[0][2]} b={wit[0][3]} vals={wit[0][4]}")
    # P3: lam-period differs from nu-period
    for dn in (1,-1,2):
        st=Counter(); wit=[]
        for n in range(m,m+5):
            for (mu,lam) in gen.pairs(n,m,7):
                for b in range(0,sum(lam)-sum(mu)+1):
                    nus=gen.slice_nus(mu,lam,n,m,b)
                    if len(nus)<2: continue
                    vals={sum(L+R for (L,R) in LRgen(nu,lam,n,m,nprime=n+dn)) for nu in nus}
                    st['slices']+=1
                    if len(vals)>1:
                        st['nonconst']+=1
                        if len(wit)<1: wit.append((n,mu,lam,b,sorted(vals)))
        v="HOLDS" if st['nonconst']==0 else f"FAILS {st['nonconst']}/{st['slices']}"
        print(f"   {f'P3  lam-period n+({dn}) != nu-period n':52s} -> {v}")
        if wit: print(f"        e.g. n={wit[0][0]} mu={wit[0][1]} lam={wit[0][2]} b={wit[0][3]} vals={wit[0][4]}")
    sys.stdout.flush()
