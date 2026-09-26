"""Control on Theorem A (product form) and on the m=1 characterisation."""
import sys, itertools, random
sys.path.insert(0,'/home/clio/projects/proofs/code-q254'); sys.path.insert(0,'/home/clio/projects/proofs/code-2026-09-19'); sys.path.insert(0,'.')
import cyl as C
from winding import weights_counted
from falsify import cyl_shapes, L3_exact, n_positive_eigs_exact
from lorentzian import compositions, is_M_convex

# (1) m=1 characterisation: K^c_alpha = [all alpha_i <= n-1]
bad=0; tot=0
for n in range(2,8):
  for L in range(0,n+1):
    for ell in range(1,6):
      K=weights_counted((0,),(L,),n,1,ell); K={a:c for a,c in K.items() if c>0}
      pred={a:1 for a in compositions(L,ell) if all(v<=n-1 for v in a)}
      tot+=1
      if K!=pred: bad+=1; print("  MISMATCH n=%d L=%d ell=%d"%(n,L,ell))
print("(1) m=1 characterisation K^c_alpha=[alpha_i<=n-1]: %d/%d agree"%(tot-bad,tot))

# (2) Theorem A: product form with log-concave f_t  =>  Lorentzian.
#     random log-concave f_t (indicator-of-interval times log-concave values)
def randlc(kmax,rng):
    # random log-concave nonneg sequence on a random interval
    a=rng.randint(0,kmax); b=rng.randint(a,kmax)
    vals={}
    cur=rng.randint(1,4); step=[]
    # build log-concave by decreasing increments of log
    r=rng.uniform(0.3,3.0)
    v=float(rng.randint(1,5))
    for k in range(a,b+1):
        vals[k]=max(1,int(round(v)))
        v=v*r; r=r*rng.uniform(0.3,1.0)   # ratios decreasing => log-concave-ish
    # enforce exact log-concavity by repair
    ks=sorted(vals)
    for _ in range(50):
        okk=True
        for idx in range(1,len(ks)-1):
            k=ks[idx]
            if vals[k]**2 < vals[ks[idx-1]]*vals[ks[idx+1]]:
                vals[ks[idx+1]]=vals[k]**2//max(1,vals[ks[idx-1]]); okk=False
        if okk: break
    vals={k:v for k,v in vals.items() if v>0}
    return vals

def is_lc(f):
    ks=sorted(f)
    if ks!=list(range(ks[0],ks[-1]+1)): return False
    for i in range(1,len(ks)-1):
        if f[ks[i]]**2 < f[ks[i-1]]*f[ks[i+1]]: return False
    return True

rng=random.Random(3)
tested=0; fail=0
for trial in range(4000):
    ell=rng.randint(2,4); d=rng.randint(2,6)
    fs=[randlc(4,rng) for _ in range(ell)]
    if not all(is_lc(f) for f in fs): continue
    K={}
    for g in compositions(d,ell):
        p=1
        for t in range(ell):
            p*= fs[t].get(g[t],0)
            if p==0: break
        if p>0: K[g]=p
    if not K: continue
    tested+=1
    if not is_M_convex(set(K)): fail+=1; print("  (L2) FAIL",fs,d); continue
    ok,wit=L3_exact(K,ell,d)
    if not ok: fail+=1; print("  (L3) FAIL fs=%s d=%d wit=%s"%(fs,d,wit))
print("(2) Theorem A control: %d product-form instances with log-concave f_t, %d failures"%(tested,fail))

# (3) NEGATIVE control on Theorem A: drop log-concavity of f, must be able to FAIL
tested=0; fail=0
for trial in range(6000):
    ell=rng.randint(2,4); d=rng.randint(2,6)
    fs=[{k:rng.randint(0,5) for k in range(0,4)} for _ in range(ell)]
    fs=[{k:v for k,v in f.items() if v>0} for f in fs]
    if any(not f for f in fs): continue
    if all(is_lc(f) for f in fs): continue          # want NON-log-concave
    K={}
    for g in compositions(d,ell):
        p=1
        for t in range(ell):
            p*=fs[t].get(g[t],0)
            if p==0: break
        if p>0: K[g]=p
    if not K: continue
    tested+=1
    ok,_=L3_exact(K,ell,d)
    if not ok: fail+=1
print("(3) NEGATIVE control (f_t NOT log-concave): %d instances, %d fail (L3) -- must be >0"%(tested,fail))
