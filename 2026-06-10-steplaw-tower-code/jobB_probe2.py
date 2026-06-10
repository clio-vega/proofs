"""
JOB B probe 2 — is the '4-core determines d' headline real or a low-sample artifact?
And what (if anything) predicts d?
"""
import sys
from collections import defaultdict, Counter
from job1_tie_census import partitions
from jobB_tower import deflation_trace
from core_quotient import core_quotient

MMAX=int(sys.argv[1]) if len(sys.argv)>1 else 12
ties=[]
for m in range(1,MMAX+1):
    for lam in partitions(2*m):
        info=deflation_trace(lam,m)
        if len(info['J'])<2: continue
        core,quot=core_quotient(lam,4)
        ties.append((m,lam,info,tuple(core),tuple(tuple(q) for q in quot)))

# sample size of single-valued vs multi-valued cores
d_by_core=defaultdict(list)
for (m,lam,info,core,quot) in ties:
    if info['d'] is not None:
        d_by_core[core].append(info['d'])
single_n=[len(v) for c,v in d_by_core.items() if len(set(v))==1]
multi_n=[len(v) for c,v in d_by_core.items() if len(set(v))>1]
print(f"single-valued cores: {len(single_n)}, their sample sizes: {sorted(single_n,reverse=True)}")
print(f"multi-valued cores:  {len(multi_n)}, their sample sizes: {sorted(multi_n,reverse=True)}")
print(f"ties on single-valued cores: {sum(single_n)} / {sum(single_n)+sum(multi_n)} total")
print(f"  => {sum(multi_n)} ties live on the 5 multi-valued cores")

# within fixed core (2,2): does 4-quotient or |J*| predict d?
print("\n=== within 4-core=(2,2): predictors of d ===")
sub=[(m,lam,info,quot) for (m,lam,info,core,quot) in ties if core==(2,2) and info['d'] is not None]
print(f"  {len(sub)} ties")
def det(keyfun,label):
    by=defaultdict(set)
    for (m,lam,info,quot) in sub:
        by[keyfun(m,lam,info,quot)].add(info['d'])
    s=sum(1 for v in by.values() if len(v)==1)
    print(f"    {label:35s}: {s}/{len(by)} single-valued")
det(lambda m,lam,info,quot:(quot,), "by 4-quotient")
det(lambda m,lam,info,quot:(len(info['J']),), "by |J*|")
det(lambda m,lam,info,quot:(m,), "by m")
det(lambda m,lam,info,quot:(len(info['J']),quot), "by (|J*|,4-quot)")

# global: best simple predictor of d
print("\n=== global predictors of d (all ties) ===")
def detg(keyfun,label):
    by=defaultdict(set)
    for (m,lam,info,core,quot) in ties:
        if info['d'] is None: continue
        by[keyfun(m,lam,info,core,quot)].add(info['d'])
    s=sum(1 for v in by.values() if len(v)==1)
    print(f"  {label:45s}: {s}/{len(by)} single-valued")
detg(lambda m,lam,info,core,quot:(len(info['J']),quot), "by (|J*|, 4-quot)")
detg(lambda m,lam,info,core,quot:(core,quot), "by (4-core,4-quot)=lam [vacuous]")
detg(lambda m,lam,info,core,quot:(info['d1'],), "by d1 (within-J* depth)")
detg(lambda m,lam,info,core,quot:(info['d1'],len(info['J'])), "by (d1,|J*|)")

# correlation d vs d1
print("\n=== d vs d1 ===")
eq=sum(1 for (m,lam,info,core,quot) in ties if info['d'] is not None and info['d']==info['d1'])
gt=sum(1 for (m,lam,info,core,quot) in ties if info['d'] is not None and info['d']>info['d1'])
lt=sum(1 for (m,lam,info,core,quot) in ties if info['d'] is not None and info['d']<info['d1'])
print(f"  d==d1 (within-J* is survivor): {eq}")
print(f"  d <d1 (next levels LOWER the survivor): {lt}")
print(f"  d >d1 (within-J* over-cancels, recovered higher): {gt}")
