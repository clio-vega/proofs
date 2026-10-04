"""Round 2: sharpen C1, C2, C3 and size the hedge."""
import itertools, random, sys
from lor import *
random.seed(41004)

print("=== C1b. WIDER: delete one interior support point, chain must refuse ===")
tot=ref=0; reasons={}
sizes=[]
for h in range(2,8):
    for l in range(0,h+1):
        c=Ptilde(l,h); sizes.append(len(c))
        for p in sorted(c):
            if len(c)<2: continue
            c2=dict(c); del c2[p]; tot+=1
            ok,w=is_N_lorentzian(c2,return_witness=True)
            if not ok:
                ref+=1; r=w if isinstance(w,str) else 'Hessian'
                reasons[r]=reasons.get(r,0)+1
print(f"  deletions enumerated {tot}, refusals {ref}; support sizes {min(sizes)}..{max(sizes)}")
print(f"  refusal reasons: {reasons}")

print("=== C2b. crossing family: a witness with INTERVAL support where LOG-CONCAVITY fails ===")
def k_cross(m,l,h,D,cross):
    out={}
    for x in itertools.product(*[range(0,h[i]+1) for i in range(m)]):
        for e in itertools.product(*[range(0,h[i]+1) for i in range(m)]):
            if sum(x)+sum(e)!=D: continue
            if any(not(l[i]<=x[i]+e[i]<=h[i]) for i in range(m)): continue
            clo,chi=cross
            if any(not(clo<=e[i]+x[(i+1)%m]<=chi) for i in range(m)): continue
            out[sum(x)]=out.get(sum(x),0)+1
    return [out.get(a,0) for a in range(0,D+1)] if out else []
best=None; nlc_interval=0; tot=0
for cross in [(0,1),(1,2),(0,2),(0,3),(1,3)]:
    for m in (2,3):
        for h in itertools.product(range(1,5),repeat=m):
            for l in itertools.product(*[range(0,hi+1) for hi in h]):
                for D in range(sum(l),sum(h)+1):
                    ks=k_cross(m,list(l),list(h),D,cross)
                    if not ks or all(v==0 for v in ks): continue
                    tot+=1
                    nz=[i for i,v in enumerate(ks) if v]
                    interval = all(ks[i]!=0 for i in range(nz[0],nz[-1]+1))
                    lcfail = any(ks[a]**2 < ks[a-1]*ks[a+1] for a in range(nz[0]+1,nz[-1]))
                    if interval and lcfail:
                        nlc_interval+=1
                        if best is None or sum(ks)>sum(best[-1]): best=(cross,m,l,h,D,ks)
print(f"  instances enumerated {tot}; NOT log-concave WITH interval support: {nlc_interval}")
print(f"  strongest such witness: cross={best[0]} m={best[1]} l={best[2]} h={best[3]} D={best[4]}")
print(f"     k = {best[5]}   (interval support, log-concavity FAILS)")

print("=== C3b. isolate the M-convexity hypothesis: all-1 coefficients, support NOT ===")
print("===      M-convex, yet EVERY Hessian condition PASSES.  Then only the support ===")
print("===      hypothesis stands between normalizedcoefficients and a false certificate ===")
def mono(d,n=3):
    def rec(k,rem):
        if k==1: yield (rem,);return
        for v in range(rem+1):
            for t in rec(k-1,rem-v): yield (v,)+t
    return list(rec(n,d))
found=[]
for d in (3,4):
    ms=mono(d)
    for r in range(2,min(6,len(ms))+1):
        if len(found)>=4: break
        for S in itertools.combinations(ms,r):
            c={a:1 for a in S}
            if is_Mconvex(list(S)): continue
            if is_N_lorentzian(c,check_support=False):
                found.append((d,S)); 
                if len(found)>=4: break
for d,S in found:
    c={a:1 for a in S}
    print(f"  degree {d}, support {S}: M-convex=False, all Hessian conditions PASS, "
          f"certified Lorentzian={is_N_lorentzian(c)}")
if not found: print("  none found in range -- the two hypotheses are not separable here")

print("=== HEDGE b. size the M-concavity measurement; report OBJECT SIZES not just counts ===")
def is_M_concave(c):
    supp=[a for a,v in c.items() if v!=0]; S=set(supp)
    if not is_Mconvex(supp): return False
    n=len(supp[0])
    for a in supp:
        for b in supp:
            for i in range(n):
                if a[i]>b[i]:
                    ok=False
                    for j in range(n):
                        if a[j]<b[j]:
                            a2=list(a);a2[i]-=1;a2[j]+=1;a2=tuple(a2)
                            b2=list(b);b2[i]+=1;b2[j]-=1;b2=tuple(b2)
                            if a2 in S and b2 in S and c[a]*c[b]<=c[a2]*c[b2]: ok=True;break
                    if not ok: return False
    return True
tot=yes=0; sizes=[]; degs=[]; firstno=None
for m in (1,2,3,4):
    hmax = 5 if m<=2 else (4 if m==3 else 3)
    for h in itertools.product(range(1,hmax+1),repeat=m):
        for l in itertools.product(*[range(0,hi+1) for hi in h]):
            f={(0,0,0):1}
            for i in range(m): f=poly_mult(f,Ptilde(l[i],h[i]))
            tot+=1; sizes.append(len(f)); degs.append(sum(h))
            if is_M_concave(f): yes+=1
            elif firstno is None: firstno=(m,l,h)
print(f"  products enumerated {tot}; nu M-concave in {yes}; first failure {firstno}")
print(f"  support sizes {min(sizes)}..{max(sizes)}; degrees H {min(degs)}..{max(degs)}; m up to 4")
