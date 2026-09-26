import sys, itertools
sys.path.insert(0,'/home/clio/projects/proofs/code-q254'); sys.path.insert(0,'/home/clio/projects/proofs/code-2026-09-19'); sys.path.insert(0,'/home/clio/projects/proofs/code-q255')
import cyl as C
from winding import weights_counted, winds
from lorentzian import is_M_convex
from falsify import cyl_shapes, L3_exact, rlc
tot=wind=f2=f3=fr=0
fails=[]
for n in range(2,7):
  for m in range(1,n+1):
    for mu in cyl_shapes(n,m):
      for lam in itertools.product(*[range(mu[i],mu[i]+n+1) for i in range(m)]):
        if not C.is_shape(lam,n,m) or not C.contains(lam,mu): continue
        d=C.size(mu,lam)
        if d==0 or d>9: continue
        for ell in range(2,7):
          if d-2>=0 and len(list(itertools.combinations(range(d-2+ell-1),ell-1)))>3000: continue
          K=weights_counted(mu,lam,n,m,ell)
          K={a:c for a,c in K.items() if c>0}
          if not K: continue
          w,_,_=winds(mu,lam,n,m,ell); tot+=1; wind+=1 if w else 0
          if not is_M_convex(set(K)): f2+=1; fails.append(('L2',n,m,mu,lam,ell,d,w))
          ok3,wit=L3_exact(K,ell,d,want_witness=True)
          if not ok3: f3+=1; fails.append(('L3',n,m,mu,lam,ell,d,w,wit[:2]))
          b=rlc(K,ell,raw=True)
          if b: fr+=1; fails.append(('RLC',n,m,mu,lam,ell,d,w,b[0]))
print("EXTENDED: instances=%d winding=%d (%.1f%%)  L2fail=%d L3fail=%d RLCfail=%d"%(tot,wind,100*wind/max(1,tot),f2,f3,fr))
for f in fails[:15]: print("  ",f)
