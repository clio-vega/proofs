"""
JOB A probe 2 — non-vacuous determination tests for v2(M_j).

(core,quotient) TOGETHER = lam (bijection), so that determination is vacuous.
Test SINGLE invariants and additive decompositions instead.
"""
import sys, math
from collections import defaultdict, Counter
from job1_tie_census import partitions, M_vector, v2
from core_quotient import core_quotient

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10

def det(keyfun, label, restrict=None):
    by=defaultdict(set)
    for m in range(1,MMAX+1):
        for lam in partitions(2*m):
            Ms=M_vector(lam,m)
            for j in range(m+1):
                if Ms[j]==0: continue
                if restrict and not restrict(m,lam,j): continue
                by[keyfun(m,lam,j)].add(v2(Ms[j]))
    single=sum(1 for v in by.values() if len(v)==1)
    print(f"  {label:50s}: {single}/{len(by)} single-valued")

def cq(lam,d):
    core,quot=core_quotient(lam,d)
    return tuple(core), tuple(tuple(q) for q in quot)

print(f"=== single-invariant determination of v2(M_j) (m<={MMAX}) ===")
det(lambda m,lam,j:(m,j,cq(lam,2)[1]), "by (m,j,2-quot)")
det(lambda m,lam,j:(m,j,cq(lam,2)[0]), "by (m,j,2-core)")
det(lambda m,lam,j:(m,j,cq(lam,4)[0]), "by (m,j,4-core)")
det(lambda m,lam,j:(m,j,cq(lam,4)[1]), "by (m,j,4-quot)")
det(lambda m,lam,j:(m,j,cq(lam,2)[1],cq(lam,4)[0]), "by (m,j,2-quot,4-core)")

# additive decomposition of v2(M_0)=v2(f) over 2-quotient (Macdonald):
# classical: f^lam 2-adic valuation relates to quotient f-numbers.
print(f"\n=== v2(f^lam) additive over 2-quotient?  (m<={MMAX}) ===")
# f^lam = (n!/prod hooks). For d-quotient: prod over quotient parts of f-numbers,
# times multinomial.  Test v2(f) = v2(multinomial) + sum v2(f^{quot_r}) for d=2.
def fnum(lam):
    n=sum(lam)
    # hook length formula
    lam=list(lam); k=len(lam)
    hooks=1
    for i in range(k):
        for jj in range(lam[i]):
            arm=lam[i]-jj-1
            leg=sum(1 for r in range(i+1,k) if lam[r]>jj)
            hooks*=(arm+leg+1)
    return math.factorial(n)//hooks
ok=0;tot=0
for m in range(1,MMAX+1):
    for lam in partitions(2*m):
        core,quot=core_quotient(lam,2)
        if sum(core)!=0:  # 2-core nonempty -> f formula has core correction; skip pure test
            continue
        n=sum(lam)
        sizes=[sum(q) for q in quot]
        from math import comb
        # multinomial n!/(prod (2*size_r)!)? careful: each quotient part contributes 2*size boxes
        multinom=math.factorial(n)
        for s in sizes:
            multinom//=math.factorial(2*s)
        pred=v2(multinom) if multinom else None
        for q in quot:
            if sum(q)>0:
                pred+=v2(fnum(q))
        if pred==v2(fnum(lam)):
            ok+=1
        tot+=1
print(f"  v2(f)=v2(multinom)+sum v2(f^quot) (2-core empty): {ok}/{tot}")
