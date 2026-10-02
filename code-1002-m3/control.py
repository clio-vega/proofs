"""REFUSAL CONTROLS for (C).

The brief says: break the cylindric alignment (A=D-n+1, C=B+1) and confirm (C) fails.
My derivation says the load-bearing fact is an OFFSET MATCH, not periodicity.  So I test
a two-parameter deformation and record WHICH deformations break (C):

    L_i^{p}(nu) = max(nu_i,   lam_{i-1} + p)      (p=1 is cylindric)
    R_i^{q}(nu) = min(nu_{i+1} - q, lam_i)        (q=1 is cylindric)

and separately a PERIODICITY break: replace lam_0 = lam_m - n by lam_m - n + eps.
Prediction (to be refuted if wrong):
  * p != q  ==> (C) FAILS                      [offset match is load-bearing]
  * p == q  ==> (C) HOLDS for every p          [offset VALUE is irrelevant]
  * eps != 0 ==> (C) STILL HOLDS               [periodicity is NOT load-bearing]
"""
import gen, sys
from collections import Counter

def LRdef(nu, lam, n, m, p=1, q=1, eps=0, eps_i=0):
    """eps shifts the lam-clamp seen by L_{eps_i} only (eps_i=0 is the wrap, lam_0)."""
    out=[]
    for i in range(m):
        lp = gen.lam_prev(lam,i,n,m) + (eps if i==eps_i else 0)
        nn = gen.nu_next(nu,i,n,m)
        out.append((max(nu[i], lp+p), min(nn-q, lam[i])))
    return out

def test(m, nmax, dmax, **kw):
    st=Counter(); wit=[]
    for n in range(m, nmax+1):
        for (mu,lam) in gen.pairs(n,m,dmax):
            d=sum(lam)-sum(mu)
            for b in range(0,d+1):
                nus=gen.slice_nus(mu,lam,n,m,b)
                if len(nus)<2: continue    # need >=2 summands for "constant on slice" to have content
                vals={sum(L+R for (L,R) in LRdef(nu,lam,n,m,**kw)) for nu in nus}
                st['slices']+=1
                if len(vals)==1: st['const']+=1
                else:
                    st['nonconst']+=1
                    if len(wit)<2: wit.append((n,mu,lam,b,sorted(vals)))
    return st,wit

for m in (2,3,4):
    print(f"--- m={m}, n<={m+5}, d<=8, slices with >=2 summands ---")
    for kw,label in [(dict(),'cylindric (p=q=1, eps=0)'),
                     (dict(p=2,q=2),'p=q=2  (offset matched, value changed)'),
                     (dict(p=3,q=3),'p=q=3  (offset matched, value changed)'),
                     (dict(p=2,q=1),'p=2,q=1  (OFFSET MISMATCH)'),
                     (dict(p=1,q=2),'p=1,q=2  (OFFSET MISMATCH)'),
                     (dict(eps=1),'lam_0 -> lam_m-n+1  (PERIODICITY BROKEN at the wrap)'),
                     (dict(eps=-2),'lam_0 -> lam_m-n-2  (PERIODICITY BROKEN at the wrap)'),
                     (dict(eps=2,eps_i=1),'lam_1 -> lam_1+2 seen by L_2 only (ADJACENCY BROKEN, interior)'),
                     ]:
        st,wit=test(m, m+5, 8, **kw)
        verdict = "HOLDS" if st['nonconst']==0 else f"FAILS on {st['nonconst']}/{st['slices']}"
        print(f"   {label:58s} -> (C) {verdict}")
        if wit: print(f"        e.g. n={wit[0][0]} mu={wit[0][1]} lam={wit[0][2]} b={wit[0][3]} values={wit[0][4]}")
    sys.stdout.flush()
