from regions import *
from itertools import combinations
import sys

PAIRS = list(combinations(['I','II','III','IV'],2))

def check_instance(n,A,B,u,Tm,Tp):
    """Return (pairseen dict, failures list)."""
    G,tset = region_sums(n,A,B,u,Tm,Tp)
    nz = [R for R in ['I','II','III','IV'] if G[R]]
    seen = []; fails=[]
    for R,Rp in combinations(nz,2):
        seen.append((R,Rp))
        bad = star_fails(G[R],G[Rp])
        if bad: fails.append((R,Rp,bad,as_seq(G[R]),as_seq(G[Rp])))
    return seen,fails,nz,G

# ---------------- (a) genuine cylindric sweep, n<=10, d<=14 ----------------
print("=== SWEEP (a): cylindric, 2<=n<=10, m=2, d<=14 ===")
cnt = {p:0 for p in PAIRS}; nfail=0; nslice=0; failwit=[]
both_II_IV = 0
for n in range(2,11):
  for mu1 in range(0,n+1):
    for mu2 in range(mu1, mu1+n+1):
      for lam1 in range(mu1, mu1+15):
        for lam2 in range(max(mu2,lam1), lam1+n+1):
          d = (lam1-mu1)+(lam2-mu2)
          if d>14 or d<1: continue
          if not valid_shape(n,(mu1,mu2),(lam1,lam2)): continue
          for b in range(0,d+1):
            nn,A,B,u,Tm,Tp = cyl_params(n,(mu1,mu2),(lam1,lam2),b)
            if Tm>Tp: continue
            nslice+=1
            seen,fails,nz,G = check_instance(n,A,B,u,Tm,Tp)
            if 'II' in nz and 'IV' in nz: both_II_IV+=1
            for p in seen: cnt[tuple(sorted(p,key=lambda x:['I','II','III','IV'].index(x)))]+=1
            if fails:
                nfail+=len(fails)
                if len(failwit)<5: failwit.append((n,(mu1,mu2),(lam1,lam2),b,fails))
print(f"  nonempty slices: {nslice}")
for p in PAIRS: print(f"  pair {p[0]:3s}-{p[1]:3s}: occurrences {cnt[p]}")
print(f"  total pair-instances: {sum(cnt.values())}")
print(f"  slices where BOTH II and IV nonempty: {both_II_IV}")
print(f"  (G) FAILURES: {nfail}")
for w in failwit: print("   WITNESS:",w)
