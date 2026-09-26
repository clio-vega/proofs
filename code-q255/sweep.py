"""Q256 falsifier sweep: root log-concavity (raw + normalised) and exact (L3),
over cylindric skew shapes, with the WINDING fraction reported."""
import sys, itertools
from collections import defaultdict
sys.path.insert(0,'/home/clio/projects/proofs/code-q254')
sys.path.insert(0,'/home/clio/projects/proofs/code-2026-09-19')
sys.path.insert(0,'/home/clio/projects/proofs/code-q255')
import cyl as C
from winding import weights_counted, winds
from lorentzian import is_M_convex
from falsify import cyl_shapes, L3_exact, rlc

def sweep(nmax=5, ellmax=4, Nmax=7, verbose=False):
    rows=[]
    for n in range(2,nmax+1):
      for m in range(1,n+1):
        for mu in cyl_shapes(n,m):
          for lam in itertools.product(*[range(mu[i],mu[i]+n+1) for i in range(m)]):
            if not C.is_shape(lam,n,m) or not C.contains(lam,mu): continue
            d=C.size(mu,lam)
            if d==0 or d>Nmax: continue
            for ell in range(2,ellmax+1):
              K=weights_counted(mu,lam,n,m,ell)
              K={a:c for a,c in K.items() if c>0}
              if not K: continue
              w,_,_=winds(mu,lam,n,m,ell)
              mc=is_M_convex(set(K))
              ok3,wit=L3_exact(K,ell,d,want_witness=True)
              bad_raw=rlc(K,ell,raw=True)
              bad_nrm=rlc(K,ell,raw=False)
              rows.append(dict(n=n,m=m,mu=mu,lam=lam,ell=ell,d=d,winds=w,
                               mconvex=mc,L3=ok3,wit=wit,
                               nraw=len(bad_raw),nnrm=len(bad_nrm),
                               raw0=bad_raw[0] if bad_raw else None,
                               nrm0=bad_nrm[0] if bad_nrm else None,
                               nterms=len(K)))
    return rows

if __name__=="__main__":
    rows=sweep()
    W=[r for r in rows if r['winds']]
    print(f"instances: {len(rows)}   winding: {len(W)} ({100*len(W)/len(rows):.1f}%)   non-winding: {len(rows)-len(W)}")
    print(f"(L2) M-convex failures      : {sum(1 for r in rows if not r['mconvex'])}  (winding: {sum(1 for r in W if not r['mconvex'])})")
    print(f"(L3) exact failures         : {sum(1 for r in rows if not r['L3'])}  (winding: {sum(1 for r in W if not r['L3'])})")
    print(f"RAW  root-log-concav. fails : {sum(1 for r in rows if r['nraw'])}  (winding: {sum(1 for r in W if r['nraw'])})")
    print(f"NORM root-log-concav. fails : {sum(1 for r in rows if r['nnrm'])}  (winding: {sum(1 for r in W if r['nnrm'])})")
    print()
    f3=[r for r in rows if not r['L3']]
    if f3:
        f3.sort(key=lambda r:(r['d'],r['n'],r['ell']))
        print("smallest (L3) failures:")
        for r in f3[:6]:
            print("  n=%d m=%d mu=%s lam=%s ell=%d d=%d winds=%s wit=%s"%(r['n'],r['m'],r['mu'],r['lam'],r['ell'],r['d'],r['winds'],r['wit'][:2]))
    fr=[r for r in rows if r['nraw']]
    if fr:
        fr.sort(key=lambda r:(r['d'],r['n'],r['ell']))
        print("smallest RAW log-concavity failures:")
        for r in fr[:6]:
            print("  n=%d m=%d mu=%s lam=%s ell=%d d=%d winds=%s viol=%s"%(r['n'],r['m'],r['mu'],r['lam'],r['ell'],r['d'],r['winds'],r['raw0']))
    fn=[r for r in rows if r['nnrm']]
    if fn:
        fn.sort(key=lambda r:(r['d'],r['n'],r['ell']))
        print("smallest NORM log-concavity failures:")
        for r in fn[:4]:
            print("  n=%d m=%d mu=%s lam=%s ell=%d d=%d winds=%s viol=%s"%(r['n'],r['m'],r['mu'],r['lam'],r['ell'],r['d'],r['winds'],r['nrm0']))
