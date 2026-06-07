import sys, os
sys.path.insert(0,"/home/clio/projects/scratch/2026-06-07-prove")
from residual import (Gtab, hook_f, v2, vpi, core_quotient_4, beta_set, beta_to_part)
from collections import defaultdict

# (A) Quotient stability under r -> r+4 (convention robustness)
def core_quotient_4_r(lam, r):
    lam=[x for x in lam if x>0]; ell=len(lam)
    assert r%4==0 and r>=ell
    betas=beta_set(lam,r)
    runners={j:[] for j in range(4)}
    for b in betas: runners[b%4].append(b//4)
    quot=tuple(beta_to_part(runners[j]) for j in range(4))
    cb=[]
    for j in range(4):
        for k in range(len(runners[j])): cb.append(4*k+j)
    return beta_to_part(cb), quot

def parts(n):
    def gen(n,mx):
        if n==0: yield (); return
        for k in range(min(n,mx),0,-1):
            for rest in gen(n-k,k): yield (k,)+rest
    return list(gen(n,n))

stable=True
for n in range(0,15):
    for lam in parts(n):
        ell=len([x for x in lam if x>0])
        r0=ell
        while r0%4: r0+=1
        if r0==0: r0=4
        c0,q0=core_quotient_4_r(lam,r0)
        c1,q1=core_quotient_4_r(lam,r0+4)
        c2,q2=core_quotient_4_r(lam,r0+8)
        if not(q0==q1==q2 and c0==c1==c2):
            stable=False; print("UNSTABLE",lam,q0,q1,q2)
print("(A) ordered-quotient stability under r->r+4 (r=0 mod4):", "STABLE" if stable else "UNSTABLE")

# rebuild rows
rows=[]
for m in range(1,11):
    for lam,(re,im) in Gtab(m).items():
        vp=vpi(re,im); fl=hook_f(lam); v2f=v2(fl)
        core,quot=core_quotient_4(lam)
        nci=sum(core)
        gc=(1,0) if nci==0 else Gtab(nci//2).get(tuple(core),(0,0))
        vpc=vpi(*gc); fc=hook_f(core); v2fc=v2(fc)
        rows.append(dict(m=m,lam=lam,vp=vp,v2f=v2f,core=core,quot=quot,vpc=vpc,v2fc=v2fc,re=re,im=im))

def residual(r):
    if r['vpc']==float('inf') or r['vp']==float('inf'): return None
    return r['vp']-r['v2f']-r['vpc']+r['v2fc']

# (B) Smallest counterexample by n: same ordered quotient, different residual
groups=defaultdict(list)
for r in rows:
    res=residual(r)
    if res is None: continue
    groups[r['quot']].append(r)
cex=[]
for q,rs in groups.items():
    byres=defaultdict(list)
    for r in rs: byres[residual(r)].append(r)
    if len(byres)>1:
        # pick two rows with diff residual, smallest n
        items=sorted(rs,key=lambda r:sum(r['lam']))
        seen={}
        for r in items:
            res=residual(r)
            for r2res,r2 in seen.items():
                if r2res!=res:
                    cex.append((sum(r['lam']),q,r2,r)); break
            else:
                seen[res]=r; continue
            break
cex.sort(key=lambda x:x[0])
print("\n(B) SMALLEST counterexamples (same ordered 4-quotient, different residual):")
for n,q,a,b in cex[:5]:
    print(f"  n={n} quotient={q}")
    for r in (a,b):
        print(f"     lam={r['lam']} G=({r['re']},{r['im']}) vpi={r['vp']} v2f={r['v2f']} core={r['core']} vpc={r['vpc']} v2fc={r['v2fc']} -> residual={residual(r)}")

# (C) Does residual depend on (core, ordered quotient) jointly?
g2=defaultdict(set)
for r in rows:
    res=residual(r)
    if res is None: continue
    g2[(r['core'],r['quot'])].add(res)
nc2=[(k,v) for k,v in g2.items() if len(v)>1]
print(f"\n(C) group by (4-core, ordered 4-quotient): #non-constant groups = {len(nc2)} / {len(g2)}")
for k,v in nc2[:8]: print("    ",k," residuals",sorted(v))
