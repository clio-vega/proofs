import sys, itertools, random
from collections import defaultdict
sys.path.insert(0,'/home/clio/projects/proofs/code-q254')
sys.path.insert(0,'/home/clio/projects/proofs/code-2026-09-19')
sys.path.insert(0,'/home/clio/projects/proofs/code-q255')
import cyl as C
from winding import weights_counted
from falsify import cyl_shapes
from cond import hessians, condA, condB
from route import pieces, tp2_fail

# ---- wider (T1)/(T2)
t1=t1f=t2=t2f=0; wit=[]
for n in range(2,8):
  for m in range(2,n+1):
    for mu in cyl_shapes(n,m):
      for lam in itertools.product(*[range(mu[i],mu[i]+n+1) for i in range(m)]):
        if not C.is_shape(lam,n,m) or not C.contains(lam,mu): continue
        d=C.size(mu,lam)
        if d<2 or d>9: continue
        Ws=list(pieces(mu,lam,n,m))
        if len(Ws)<2 or len(Ws)>60: continue
        vs=sorted({v for _,W in Ws for v in W})
        for (k1,W1),(k2,W2) in itertools.permutations(Ws,2):
            le=all(k1[i]<=k2[i] for i in range(m))
            if le:
                t1+=1
                if tp2_fail(W1,W2,vs): t1f+=1
            elif sum(k1)<sum(k2):
                t2+=1
                b=tp2_fail(W1,W2,vs)
                if b:
                    t2f+=1
                    if len(wit)<3: wit.append((n,m,mu,lam,k1,k2,dict(W1),dict(W2),b[0]))
print("(T1) WIDE kappa<=kappa' coordinatewise : %d pairs, %d fail TP2"%(t1,t1f))
print("(T2) WIDE incomparable, S_kappa<S_kappa': %d pairs, %d fail TP2"%(t2,t2f))
for w in wit: print("      witness n=%d m=%d mu=%s lam=%s kappa=%s kappa'=%s W=%s W'=%s at v-pair %s"%w)

# ---- (L2) with proper PF2 generators
def randpf2(rng,L):
    a=rng.randint(0,L-1); b=rng.randint(a,L-1)
    seq=[0]*L; v=rng.randint(1,6); r=rng.uniform(0.4,2.5)
    for k in range(a,b+1):
        seq[k]=max(1,int(round(v))); v*=r; r*=rng.uniform(0.35,1.0)
    for _ in range(80):
        ks=[i for i,x in enumerate(seq) if x>0]; ok=True
        for i in range(1,len(ks)-1):
            if seq[ks[i]]**2<seq[ks[i-1]]*seq[ks[i+1]]:
                seq[ks[i+1]]=max(1,seq[ks[i]]**2//max(1,seq[ks[i-1]])); ok=False
        if ok: break
    return seq
def pf2(s):
    ks=[i for i,x in enumerate(s) if x>0]
    if not ks or ks!=list(range(ks[0],ks[-1]+1)): return False
    return all(s[ks[i]]**2>=s[ks[i-1]]*s[ks[i+1]] for i in range(1,len(ks)-1))
def conv(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
def tp2(f,g): return all(f[i]*g[j]>=f[j]*g[i] for i in range(len(f)) for j in range(i+1,len(f)))
rng=random.Random(99); tot=fail=0
for t in range(600000):
    L=6
    f=randpf2(rng,L); fp=randpf2(rng,L); g=randpf2(rng,L); gp=randpf2(rng,L)
    if not(pf2(f) and pf2(fp) and pf2(g) and pf2(gp)): continue
    if not(tp2(f,fp) and tp2(g,gp)): continue
    tot+=1
    if not tp2(conv(f,g),conv(fp,gp)): fail+=1
print("(L2) convolution of TP2 pairs with PF2 rows: %d tested, %d fail"%(tot,fail))
# NEGATIVE control for (L2): drop PF2
rng=random.Random(100); tot2=fail2=0
for t in range(400000):
    L=5
    f=[rng.randint(0,4) for _ in range(L)]; fp=[rng.randint(0,4) for _ in range(L)]
    g=[rng.randint(0,4) for _ in range(L)]; gp=[rng.randint(0,4) for _ in range(L)]
    if not(tp2(f,fp) and tp2(g,gp)): continue
    if pf2(f) and pf2(fp) and pf2(g) and pf2(gp): continue
    tot2+=1
    if not tp2(conv(f,g),conv(fp,gp)): fail2+=1
print("(L2) NEGATIVE control (TP2 but not all PF2): %d tested, %d fail  (must be >0)"%(tot2,fail2))

# ---- (E1) the two reformulations
print()
bad_a=bad_b=tot_e=0
for n in range(2,7):
  for m in range(1,n+1):
    for mu in cyl_shapes(n,m):
      for lam in itertools.product(*[range(mu[i],mu[i]+n+1) for i in range(m)]):
        if not C.is_shape(lam,n,m) or not C.contains(lam,mu): continue
        d=C.size(mu,lam)
        if d<2 or d>9: continue
        K=weights_counted(mu,lam,n,m,3); K={a:c for a,c in K.items() if c>0}
        if not K: continue
        k=lambda a,b: K.get((a,b,d-a-b),0)
        # Hessian-level (A),(B)
        hA=any(condA(M,3) for _,M in hessians(K,3,d)); hB=any(condB(M,3) for _,M in hessians(K,3,d))
        # reformulated (A): log-concavity in a ; (B): log-submodularity of k
        rA=any(k(a,b)**2 < k(a+1,b)*k(a-1,b) for a in range(1,d) for b in range(0,d+1))
        rB=any(k(a+1,b)*k(a,b+1) < k(a,b)*k(a+1,b+1) for a in range(0,d+1) for b in range(0,d+1))
        tot_e+=1
        if hA!=rA: bad_a+=1
        if hB!=rB: bad_b+=1
print("(E1) %d instances: (A)<=>log-concave in a mismatches=%d ; (B)<=>log-submodular in (a,b) mismatches=%d"%(tot_e,bad_a,bad_b))
