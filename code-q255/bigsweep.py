"""Hard stress test of (A) and (B) for cylindric K^c.
  (A) M_ij^2 >= M_ii M_jj              [RLC]
  (B) M_ik M_kj >= M_ij M_kk           [3-index 'triangle' condition]
Tested at ell=3 over a much wider range, AND at ell=4,5,6 to see whether the
criterion is ell-uniform.  Exact integer arithmetic throughout."""
import sys, itertools
sys.path.insert(0,'/home/clio/projects/proofs/code-q254')
sys.path.insert(0,'/home/clio/projects/proofs/code-2026-09-19')
sys.path.insert(0,'/home/clio/projects/proofs/code-q255')
import cyl as C
from winding import weights_counted, winds
from lorentzian import compositions, is_M_convex
from falsify import cyl_shapes, n_positive_eigs_exact
from cond import hessians, condA, condB

def run(nmax, dmax, ells, tag):
    tot=hessn=wind=0; fA=fB=fL3=0; exA=[]; exB=[]; exL=[]
    for n in range(2,nmax+1):
      for m in range(1,n+1):
        for mu in cyl_shapes(n,m):
          for lam in itertools.product(*[range(mu[i],mu[i]+n+1) for i in range(m)]):
            if not C.is_shape(lam,n,m) or not C.contains(lam,mu): continue
            d=C.size(mu,lam)
            if d<2 or d>dmax: continue
            for ell in ells:
              nb=len(list(compositions(d-2,ell)))
              if nb>4000: continue
              K=weights_counted(mu,lam,n,m,ell)
              K={a:c for a,c in K.items() if c>0}
              if not K: continue
              w,_,_=winds(mu,lam,n,m,ell)
              tot+=1; wind+= 1 if w else 0
              for beta,M in hessians(K,ell,d):
                if all(v==0 for row in M for v in row): continue
                hessn+=1
                a=condA(M,ell); b=condB(M,ell)
                if a: fA+=1; exA.append((n,m,mu,lam,ell,d,beta,a[0],M))
                if b: fB+=1; exB.append((n,m,mu,lam,ell,d,beta,b[0],M))
                if n_positive_eigs_exact(M)>1:
                    fL3+=1; exL.append((n,m,mu,lam,ell,d,beta,M))
    print("[%s] shape-instances=%d (winding %d)  nonzero Hessians=%d"%(tag,tot,wind,hessn))
    print("      (A) failures=%d   (B) failures=%d   (L3) failures=%d"%(fA,fB,fL3))
    for t,ex in (("A",exA),("B",exB),("L3",exL)):
        for e in ex[:3]: print("        %s:"%t, e)
    return fA,fB,fL3

if __name__=="__main__":
    run(8, 12, [3], "ell=3 WIDE  n<=8 d<=12")
    run(6, 9,  [4,5,6], "ell=4,5,6  n<=6 d<=9")
