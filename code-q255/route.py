"""Where would a proof of (A),(B) break?  P(u,v)=#{(kappa,sigma): S_kappa=u,S_sigma=v}
   = sum_{S_kappa=u} W_kappa(v),  W_kappa = conv_i 1_[L_i(kappa),R_i(kappa)].
Tests:
 (E1) (A) <=> k(a,b) log-concave in a ;  (B) <=> P log-supermodular  [reformulations]
 (T1) kappa <= kappa' coordinatewise  =>  (W_kappa,W_kappa') TP2 ?
 (T2) S_kappa < S_kappa' (incomparable allowed) => TP2 ?
 (T3) does summing over each fibre {S_kappa=u} restore TP2 ?  (= (B), known true)
 (L2) does convolution preserve TP2 for PF2 rows?  (random test)
"""
import sys, itertools, random
from collections import defaultdict
sys.path.insert(0,'/home/clio/projects/proofs/code-q254')
sys.path.insert(0,'/home/clio/projects/proofs/code-2026-09-19')
sys.path.insert(0,'/home/clio/projects/proofs/code-q255')
import cyl as C
from winding import weights_counted
from falsify import cyl_shapes

def ext(x,i,n,m):
    q,r=divmod(i,m); return x[r]+q*n

def pieces(mu,lam,n,m):
    """yield (kappa, W_kappa as dict v->count)"""
    Ar=[range(mu[i], ext(mu,i+1,n,m)) for i in range(m)]
    for kap in itertools.product(*Ar):
        ivs=[]; ok=True
        for i in range(m):
            lo=max(kap[i], ext(lam,i-1,n,m)+1)
            hi=min(ext(kap,i+1,n,m)-1, lam[i])
            if lo>hi: ok=False; break
            ivs.append((lo,hi))
        if not ok: continue
        W={0:1}
        for (lo,hi) in ivs:
            W2=defaultdict(int)
            for s,c in W.items():
                for t in range(lo,hi+1): W2[s+t]+=c
            W=dict(W2)
        yield kap, W

def tp2_fail(f,g,vs):
    """rows f (lower) and g (upper): need f(v)g(v')>=f(v')g(v) for v<v'"""
    bad=[]
    for i,v in enumerate(vs):
        for v2 in vs[i+1:]:
            if f.get(v,0)*g.get(v2,0) < f.get(v2,0)*g.get(v,0): bad.append((v,v2))
    return bad

t1=t1f=t2=t2f=0
for n in range(2,6):
  for m in range(2,n+1):
    for mu in cyl_shapes(n,m):
      for lam in itertools.product(*[range(mu[i],mu[i]+n+1) for i in range(m)]):
        if not C.is_shape(lam,n,m) or not C.contains(lam,mu): continue
        d=C.size(mu,lam)
        if d<2 or d>7: continue
        Ws=list(pieces(mu,lam,n,m))
        if len(Ws)<2: continue
        vs=sorted({v for _,W in Ws for v in W})
        for (k1,W1),(k2,W2) in itertools.permutations(Ws,2):
            le_coord=all(k1[i]<=k2[i] for i in range(m))
            if le_coord:
                t1+=1
                if tp2_fail(W1,W2,vs): t1f+=1
            if sum(k1)<sum(k2) and not le_coord:
                t2+=1
                if tp2_fail(W1,W2,vs): t2f+=1
print("(T1) kappa<=kappa' coordinatewise : %d pairs, %d fail TP2"%(t1,t1f))
print("(T2) S_kappa<S_kappa', incomparable: %d pairs, %d fail TP2"%(t2,t2f))

# (T3) the fibre sums
t3=t3f=0
for n in range(2,6):
  for m in range(2,n+1):
    for mu in cyl_shapes(n,m):
      for lam in itertools.product(*[range(mu[i],mu[i]+n+1) for i in range(m)]):
        if not C.is_shape(lam,n,m) or not C.contains(lam,mu): continue
        d=C.size(mu,lam)
        if d<2 or d>7: continue
        F=defaultdict(lambda: defaultdict(int))
        for kap,W in pieces(mu,lam,n,m):
            for v,c in W.items(): F[sum(kap)][v]+=c
        us=sorted(F); vs=sorted({v for u in F for v in F[u]})
        for i,u in enumerate(us):
            for u2 in us[i+1:]:
                t3+=1
                if tp2_fail(F[u],F[u2],vs): t3f+=1
print("(T3) fibre sums F_u over S_kappa=u : %d pairs, %d fail TP2   [this is condition (B)]"%(t3,t3f))

# (L2) convolution preserves TP2 for PF2 rows?
rng=random.Random(5); tot=fail=0
def pf2(seq):
    ks=[i for i,x in enumerate(seq) if x>0]
    if not ks or ks!=list(range(ks[0],ks[-1]+1)): return False
    return all(seq[ks[i]]**2>=seq[ks[i-1]]*seq[ks[i+1]] for i in range(1,len(ks)-1))
def conv(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
def tp2(f,g):
    return all(f[i]*g[j]>=f[j]*g[i] for i in range(len(f)) for j in range(i+1,len(f)))
for t in range(300000):
    L=5
    f=[rng.randint(0,4) for _ in range(L)]; fp=[rng.randint(0,4) for _ in range(L)]
    g=[rng.randint(0,4) for _ in range(L)]; gp=[rng.randint(0,4) for _ in range(L)]
    if not(pf2(f) and pf2(fp) and pf2(g) and pf2(gp)): continue
    if not(tp2(f,fp) and tp2(g,gp)): continue
    tot+=1
    if not tp2(conv(f,g),conv(fp,gp)): fail+=1
print("(L2) convolution of TP2 pairs with PF2 rows: %d tested, %d fail TP2"%(tot,fail))
