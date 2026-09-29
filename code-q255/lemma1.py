"""Controls on Lemma 1 and on the interlacing reformulation of K^c at ell=3."""
import sys, itertools, random
from fractions import Fraction
sys.path.insert(0,'/home/clio/projects/proofs/code-q254')
sys.path.insert(0,'/home/clio/projects/proofs/code-2026-09-19')
sys.path.insert(0,'/home/clio/projects/proofs/code-q255')
import cyl as C
from winding import weights_counted
from lorentzian import compositions
from falsify import cyl_shapes, n_positive_eigs_exact
from cond import hessians, condA, condB, det3

# ---------- (I) Lemma 1 as pure linear algebra: random search for a counterexample
print("=== (I) Lemma 1 adversarial search: (A)&(B) nonneg symmetric 3x3 => <=1 positive eig ===")
rng=random.Random(11); tested=0; bad=0; detneg=0; e2pos=0
for trial in range(400000):
    V=[rng.randint(0,8) for _ in range(6)]
    M=[[V[0],V[1],V[2]],[V[1],V[3],V[4]],[V[2],V[4],V[5]]]
    if condA(M,3) or condB(M,3): continue
    tested+=1
    if det3(M)<0: detneg+=1
    if n_positive_eigs_exact(M)>1: bad+=1; print("   COUNTEREXAMPLE",M)
print("   %d matrices satisfy (A)&(B); %d have det<0; %d have >1 positive eigenvalue"%(tested,detneg,bad))

# negative control: the lemma's hypotheses must actually bite
rng=random.Random(12); t2=0; b2=0
for trial in range(200000):
    V=[rng.randint(0,8) for _ in range(6)]
    M=[[V[0],V[1],V[2]],[V[1],V[3],V[4]],[V[2],V[4],V[5]]]
    if not condA(M,3): continue          # (A) holds, (B) may fail
    t2+=1
    if n_positive_eigs_exact(M)>1: b2+=1
print("   NEGATIVE control: (A) alone: %d matrices, %d with >1 positive eigenvalue (must be >0)"%(t2,b2))
# the brief's own witness
W=[[4,6,4],[6,1,4],[4,4,4]]
print("   brief witness [[4,6,4],[6,1,4],[4,4,4]]: condA fails=%s condB fails=%s npos=%d det=%d"
      %(bool(condA(W,3)),bool(condB(W,3)),n_positive_eigs_exact(W),det3(W)))

# ---------- (II) prediction: det M == 0  <=>  some M_ii == 0
print()
print("=== (II) prediction from Lemma 1: all M_ii>0 => det M > 0 strictly ===")
viol=0; tot=0; z0=0; zd=0
for n in range(2,7):
  for m in range(1,n+1):
    for mu in cyl_shapes(n,m):
      for lam in itertools.product(*[range(mu[i],mu[i]+n+1) for i in range(m)]):
        if not C.is_shape(lam,n,m) or not C.contains(lam,mu): continue
        d=C.size(mu,lam)
        if d<2 or d>9: continue
        K=weights_counted(mu,lam,n,m,3); K={a:c for a,c in K.items() if c>0}
        if not K: continue
        for beta,M in hessians(K,3,d):
            if all(v==0 for r in M for v in r): continue
            tot+=1
            posdiag=all(M[i][i]>0 for i in range(3)); dt=det3(M)
            if dt==0: z0+=1; zd+= 1 if not posdiag else 0
            if posdiag and dt<=0: viol+=1; print("   VIOLATION",beta,M,dt)
print("   %d Hessians; %d have det==0, of which %d have a zero diagonal entry; violations=%d"%(tot,z0,zd,viol))

# ---------- (III) the interlacing reformulation of K^c
print()
print("=== (III) interlacing model:  kappa_i in [mu_i,mu_{i+1}-1], sigma_i in [lam_{i-1}+1,lam_i],")
print("         kappa_i <= sigma_i <= kappa_{i+1}-1   (cyclic, x_{i+m}=x_i+n) ===")
def ext(x,i,n,m):
    q,r=divmod(i,m); return x[r]+q*n
def model(mu,lam,n,m):
    """returns dict (a,b,c) -> count, built from the interlacing description"""
    from collections import defaultdict
    out=defaultdict(int)
    Ar=[range(mu[i], ext(mu,i+1,n,m)) for i in range(m)]
    for kap in itertools.product(*Ar):
        Br=[]
        ok=True
        for i in range(m):
            lo=max(kap[i], ext(lam,i-1,n,m)+1)
            hi=min(ext(kap,i+1,n,m)-1, lam[i])
            if lo>hi: ok=False; break
            Br.append(range(lo,hi+1))
        if not ok: continue
        for sig in itertools.product(*Br):
            a=sum(kap)-sum(mu); b=sum(sig)-sum(kap); c=sum(lam)-sum(sig)
            out[(a,b,c)]+=1
    return dict(out)
agree=0; disagree=0
for n in range(2,7):
  for m in range(1,n+1):
    for mu in cyl_shapes(n,m):
      for lam in itertools.product(*[range(mu[i],mu[i]+n+1) for i in range(m)]):
        if not C.is_shape(lam,n,m) or not C.contains(lam,mu): continue
        d=C.size(mu,lam)
        if d>8: continue
        K=weights_counted(mu,lam,n,m,3); K={a:c for a,c in K.items() if c>0}
        Mo=model(mu,lam,n,m); Mo={a:c for a,c in Mo.items() if c>0}
        if K==Mo: agree+=1
        else:
            disagree+=1
            if disagree<=3: print("   MISMATCH n=%d m=%d mu=%s lam=%s"%(n,m,mu,lam)); print("     dp=%s"%K); print("     model=%s"%Mo)
print("   interlacing model vs transfer-matrix DP: %d agree, %d disagree"%(agree,disagree))
