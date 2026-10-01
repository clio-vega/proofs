"""Extend the (A)-at-m=2 verification range well past anything swept before.
Direct trapezoid sum (NOT the radial shortcut), so this is independent of the proof."""
from regions import *
import sys
tot=0; fail=0; ex=None; maxn=0; maxd=0
NMAX,DMAX=16,20
for n in range(2,NMAX+1):
  for mu2 in range(0,n+1):               # mu1=0 wlog (translation invariance)
    for lam1 in range(0,DMAX+1):
      for lam2 in range(max(mu2,lam1),lam1+n+1):
        d=lam1+(lam2-mu2)
        if d<1 or d>DMAX: continue
        if not valid_shape(n,(0,mu2),(lam1,lam2)): continue
        for b in range(0,d+1):
          nn,A,B,u,Tm,Tp=cyl_params(n,(0,mu2),(lam1,lam2),b)
          if Tm>Tp: continue
          G={}
          for t in range(Tm,Tp+1):
            for s,v in conv_interval(*endpoints(t,n,A,B,u)).items(): G[s]=G.get(s,0)+v
          G={k:v for k,v in G.items() if v>0}
          if not G: continue
          tot+=1; maxn=max(maxn,n); maxd=max(maxd,d)
          if not is_pf2(G):
            fail+=1
            if ex is None: ex=(n,(0,mu2),(lam1,lam2),b,as_seq(G))
print(f"range n<={NMAX}, d<={DMAX}  (mu1=0 wlog)")
print(f"nonzero slices: {tot}   max n seen {maxn}, max d seen {maxd}")
print(f"condition (A) FAILURES: {fail}   {ex}")
# control: the instrument must be able to say no
print("control is_pf2 on a known non-log-concave sequence (1,1,2,1,1):",
      is_pf2({0:1,1:1,2:2,3:1,4:1}))
