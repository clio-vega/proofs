"""
JOB B probe — pin down the d = v_pi(G)-mu predictor.

Headline from jobB_tower: d is single-valued on 34/39 distinct 4-cores.
Dig into the 5 exceptions; test refined predictors; stress the worked cases.
"""
import sys, math
from collections import defaultdict, Counter
from job1_tie_census import partitions, M_vector, v2
from jobB_tower import deflation_trace, surviving_index_set, is_box
from core_quotient import core_quotient

MMAX=int(sys.argv[1]) if len(sys.argv)>1 else 12

# collect ties
ties=[]
for m in range(1,MMAX+1):
    for lam in partitions(2*m):
        info=deflation_trace(lam,m)
        if len(info['J'])<2: continue
        core,quot=core_quotient(lam,4)
        ties.append((m,lam,info,tuple(core),tuple(tuple(q) for q in quot)))

print(f"ties m<={MMAX}: {len(ties)}")

# d by 4-core, list exceptions
d_by_core=defaultdict(list)
for (m,lam,info,core,quot) in ties:
    if info['d'] is not None:
        d_by_core[core].append((info['d'],lam,m,len(info['J'])))
print("\n=== 4-cores with NON-unique d (the exceptions) ===")
exc=0
for core in sorted(d_by_core):
    ds=set(x[0] for x in d_by_core[core])
    if len(ds)>1:
        exc+=1
        print(f"  core={core}: d in {sorted(ds)}   (n={len(d_by_core[core])})")
        # show the split by |J*| and quotient
        for (d,lam,m,sz) in sorted(d_by_core[core]):
            pass
print(f"  ({exc} exceptional cores)")

# refined predictor: d by (4-core, |J*|)?
def det(keyfun,label):
    by=defaultdict(set)
    for (m,lam,info,core,quot) in ties:
        if info['d'] is None: continue
        by[keyfun(m,lam,info,core,quot)].add(info['d'])
    single=sum(1 for v in by.values() if len(v)==1)
    print(f"  {label:45s}: {single}/{len(by)} single-valued")

print("\n=== refined d predictors ===")
det(lambda m,lam,info,core,quot:(core,), "by 4-core")
det(lambda m,lam,info,core,quot:(core,len(info['J'])), "by (4-core, |J*|)")
det(lambda m,lam,info,core,quot:(core,info['J'][0]%2), "by (4-core, parity)")
det(lambda m,lam,info,core,quot:(core,len(info['J']),info['J'][0]%2), "by (4-core,|J*|,parity)")
det(lambda m,lam,info,core,quot:(len(info['J']),), "by |J*| alone")
det(lambda m,lam,info,core,quot:(quot,), "by 4-quotient")

# Is d1 (within-J*) better predicted?
print("\n=== d1 (within-J* depth) predictors ===")
def det1(keyfun,label):
    by=defaultdict(set)
    for (m,lam,info,core,quot) in ties:
        if info['d1'] is None: continue
        by[keyfun(m,lam,info,core,quot)].add(info['d1'])
    single=sum(1 for v in by.values() if len(v)==1)
    print(f"  {label:45s}: {single}/{len(by)} single-valued")
det1(lambda m,lam,info,core,quot:(core,), "d1 by 4-core")
det1(lambda m,lam,info,core,quot:(len(info['J']),), "d1 by |J*|")

# relationship d vs |J*|: is d >= something like C(|J*| structure)?
print("\n=== d vs |J*| cross-tab ===")
ct=defaultdict(Counter)
for (m,lam,info,core,quot) in ties:
    if info['d'] is None: continue
    ct[len(info['J'])][info['d']]+=1
for sz in sorted(ct):
    print(f"  |J*|={sz}: min d={min(ct[sz])}, max d={max(ct[sz])}, dist={dict(sorted(ct[sz].items()))}")

# ---- stress cases ----
print("\n=== STRESS CASES ===")
def show(lam,m):
    info=deflation_trace(lam,m)
    core,quot=core_quotient(lam,4)
    surv=surviving_index_set(lam,m,info)
    print(f"\nlam={lam} m={m}")
    print(f"  M={info['Ms']}")
    print(f"  mu={info['mu']}  J*={info['J']} (|J*|={len(info['J'])})  vpi(G)={info['vpG']}  d={info['d']}  d1={info['d1']}")
    print(f"  4-core={core}  4-quot={quot}")
    print(f"  surviving-order index set={surv}  is_box={None if surv is None else is_box(surv)}")
    print(f"  deflation trace (level, new indices, v_pi(cumulative)):")
    for (v,newj,vp) in info['trace']:
        print(f"     level {v}: +{newj}  -> v_pi={vp}")

show((2,2),2)
show((3,3,1,1),4)
# find the largest-d tie
big=max((t for t in ties if t[2]['d'] is not None), key=lambda t:t[2]['d'])
show(big[1],big[0])
