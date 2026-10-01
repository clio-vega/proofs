"""ESTIMAND CHECK: does my G(s) really equal k(a,b) = K^c_{(a,b,d-a-b)}?

Instrument is INDEPENDENT of regions.py: it enumerates chains mu < nu < kappa < lam
directly with the vendored primitives from proofs/code-q254/cyl.py (the code behind the
2026-09-20 M-convexity paper), never touching a_t,b_t,c_t,d_t."""
import sys, itertools
sys.path.insert(0,'/home/clio/projects/proofs/code-q254')
from cyl import is_shape, hstrip, size
from regions import cyl_params, endpoints, conv_interval, valid_shape, as_seq, is_pf2

def chain_counts(n, mu, lam, m=2):
    """N[(a,b)] = #{(nu,kappa): mu<nu<kappa<lam hstrips, |nu/mu|=a, |kappa/nu|=b}"""
    shapes=[x for x in itertools.product(range(min(mu)-1,max(lam)+2),repeat=m)
            if is_shape(x,n,m)]
    N={}
    for nu in shapes:
        if not hstrip(mu,nu,n,m): continue
        a=size(mu,nu)
        for ka in shapes:
            if not hstrip(nu,ka,n,m): continue
            if not hstrip(ka,lam,n,m): continue
            b=size(nu,ka)
            N[(a,b)]=N.get((a,b),0)+1
    return N

def G_from_regions(n,A,B,u,Tm,Tp):
    G={}
    for t in range(Tm,Tp+1):
        for s,v in conv_interval(*endpoints(t,n,A,B,u)).items(): G[s]=G.get(s,0)+v
    return G

print("=== CONTROL: instrument must reproduce the symmetry k(a,b)=k(b,a) of s^c ===")
sym_ok=sym_bad=0
for n in range(2,8):
  for mu2 in range(0,n+1):
    for lam1 in range(0,7):
      for lam2 in range(max(mu2,lam1),lam1+n+1):
        mu,lam=(0,mu2),(lam1,lam2)
        if not valid_shape(n,mu,lam) or not is_shape(mu,n,2) or not is_shape(lam,n,2): continue
        N=chain_counts(n,mu,lam)
        for (a,b),v in N.items():
            if N.get((b,a),0)==v: sym_ok+=1
            else: sym_bad+=1
print(f"   symmetric entries {sym_ok}, asymmetric {sym_bad}")
assert sym_bad==0 and sym_ok>0, "INSTRUMENT BROKEN: s^c symmetry fails"

print()
print("=== MAIN: k(a,b) vs G(u+a) from the four-region reduction ===")
tot=0; bad=0; ex=None; nonlc=0
for n in range(2,8):
  for mu2 in range(0,n+1):
    for lam1 in range(0,8):
      for lam2 in range(max(mu2,lam1),lam1+n+1):
        mu,lam=(0,mu2),(lam1,lam2)
        if not valid_shape(n,mu,lam) or not is_shape(mu,n,2) or not is_shape(lam,n,2): continue
        d=(lam[0]-mu[0])+(lam[1]-mu[1])
        if d<1 or d>11: continue
        N=chain_counts(n,mu,lam)
        for b in range(0,d+1):
            nn,A,B,u,Tm,Tp=cyl_params(n,mu,lam,b)
            G = G_from_regions(n,A,B,u,Tm,Tp) if Tm<=Tp else {}
            # k(a,b) in the (nu,kappa) parametrisation: |nu/mu|=b, |kappa/nu|=a
            kb = {a:v for (aa,a),v in N.items() if aa==b}
            Gs = {s-u:v for s,v in G.items() if v>0}
            tot+=1
            if kb!={k:v for k,v in Gs.items() if v>0}:
                bad+=1
                if ex is None: ex=(n,mu,lam,b,dict(sorted(kb.items())),dict(sorted(Gs.items())))
            if kb and not is_pf2({k:v for k,v in kb.items()}): nonlc+=1
print(f"   (shape,b) pairs compared: {tot}")
print(f"   k(.,b) != G(u+.) : {bad}   {ex}")
print(f"   k(.,b) not PF2 (direct chain count): {nonlc}")
