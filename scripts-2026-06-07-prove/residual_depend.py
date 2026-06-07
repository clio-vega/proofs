import sys
sys.path.insert(0,"/home/clio/projects/scratch/2026-06-07-prove")
from residual import Gtab, hook_f, v2, vpi, core_quotient_4
from collections import defaultdict

rows=[]
for m in range(1,11):
    for lam,(re,im) in Gtab(m).items():
        vp=vpi(re,im); fl=hook_f(lam); v2f=v2(fl)
        core,quot=core_quotient_4(lam); nci=sum(core)
        gc=(1,0) if nci==0 else Gtab(nci//2).get(tuple(core),(0,0))
        vpc=vpi(*gc); fc=hook_f(core); v2fc=v2(fc)
        if vpc==float('inf') or vp==float('inf'): res=None
        else: res=vp-v2f-vpc+v2fc
        rows.append(dict(lam=lam,vp=vp,v2f=v2f,core=core,quot=quot,vpc=vpc,v2fc=v2fc,res=res,
                         qsz=sum(sum(q) for q in quot), csz=nci))

# H1: residual function of (v_pi(G_core), v2(f_core), quotient)?  i.e. core only via its valuations
g=defaultdict(set)
for r in rows:
    if r['res'] is None: continue
    g[(r['vpc'],r['v2fc'],r['quot'])].add(r['res'])
nc=[(k,sorted(v)) for k,v in g.items() if len(v)>1]
print("H1 [residual = f(vpc, v2fc, quotient)]: non-constant groups =",len(nc),"/",len(g))
for k,v in nc[:6]: print("   ",k,v)

# H2: residual function of (core-size, quotient)?
g=defaultdict(set)
for r in rows:
    if r['res'] is None: continue
    g[(r['csz'],r['quot'])].add(r['res'])
nc=[(k,sorted(v)) for k,v in g.items() if len(v)>1]
print("\nH2 [residual = f(core-size, quotient)]: non-constant groups =",len(nc),"/",len(g))
for k,v in nc[:6]: print("   ",k,v)

# H3: For fixed quotient, tabulate core -> residual (a few simple quotients)
print("\nH3  core -> residual, at fixed quotient:")
for target in [((),(),(1,),()), ((),(),(),(1,)), ((1,),(),(),())]:
    m={}
    for r in rows:
        if r['res'] is None or r['quot']!=target: continue
        m.setdefault(r['core'],set()).add(r['res'])
    print(f"  quotient {target}:")
    for core in sorted(m,key=lambda c:(sum(c),c)):
        print(f"      core={core}: residual(s)={sorted(m[core])}")

# H4: Is residual ADDITIVE: res(lam) =? r_core + r_quot for some functions?
# Test: does res - (res at empty quotient, same core) depend only on quotient?
# Build base_core[core] = residual when quotient is empty (lam=core): that's res=0 by defn (lam=core => vp=vpc,v2f=v2fc => res=0). 
# So compare res(lam) to "res of some canonical rep with same quotient".
# Instead: fix core, see if res - f(core) is quotient-only. Try f(core)=0:
# already know res not quotient-only. Try res as function of (core, qsz):
g=defaultdict(set)
for r in rows:
    if r['res'] is None: continue
    g[(r['core'],r['qsz'])].add(r['res'])
nc=[(k,sorted(v)) for k,v in g.items() if len(v)>1]
print("\nH5 [residual = f(core, quotient-size)]: non-constant groups =",len(nc),"/",len(g))
for k,v in nc[:8]: print("   ",k,v)
