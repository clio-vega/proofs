"""
JOB A probe — what governs v2(M_j)?  Simple global forms failed; dig deeper.

Tests:
  (P1) parity of M_j: is "M_j odd" a Lucas condition on (j, m) or on the 4-core?
  (P2) refine C1: v2(M_j) determined by (m, j, 4-core, 4-quotient) [add core]?
  (P3) v2(M_j) determined by (m, j, full lam mod nothing) is trivially yes; test
       the minimal determining set among {m, j, 4-core, 4-quotient, 2-core, 2-quot}.
  (P4) the actual lever: is  v2(M_j) - v2(M_{j-1})  or the discrete structure clean?
  (P5) Mod-2 engine cross-check: M_j mod 2 vs chi_b mod 2.
"""
import sys, math
from collections import defaultdict, Counter
from job1_tie_census import partitions, M_vector, v2, chi_b_vector
from core_quotient import core_quotient

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10

# ---- P1: parity of M_j ----
def kummer_and(j, x):
    """j submask of x ?  (Lucas: C(x,j) odd iff j submask x)."""
    return (j & x) == j

print(f"=== P1  parity of M_j (m<={MMAX}) ===")
# hypothesis: M_j odd iff f^lam odd AND j submask of (m)?  or of something.
# Test several masks.
tests = {
    "j submask m": lambda m,lam,j,f: (j & m)==j,
    "f odd & j submask m": None,
}
cnt = Counter(); tot=0
# Build: for each shape with f odd, what is the set {j: M_j odd}? is it submasks of a fixed mask?
shapes_mask=0; shapes_tot=0
for m in range(1,MMAX+1):
    for lam in partitions(2*m):
        Ms=M_vector(lam,m)
        oddj=[j for j in range(m+1) if Ms[j]%2==1]
        if not oddj:
            continue
        shapes_tot+=1
        # is oddj exactly the submasks of some mask X? then X = OR of oddj and count=2^popcount
        X=0
        for j in oddj: X|=j
        submasks=[j for j in range(X+1) if (j&X)==j]
        if set(oddj)==set(submasks) and len(oddj)==(1<<bin(X).count("1")):
            shapes_mask+=1
print(f"  shapes where {{j: M_j odd}} = submasks of a single mask X: {shapes_mask}/{shapes_tot}")

# ---- P2/P3: minimal determining set for v2(M_j) ----
print(f"\n=== P2/P3  determination of v2(M_j) by invariants (m<={MMAX}) ===")
def det_test(keyfun, label):
    by=defaultdict(set)
    for m in range(1,MMAX+1):
        for lam in partitions(2*m):
            Ms=M_vector(lam,m)
            for j in range(m+1):
                if Ms[j]==0: continue
                by[keyfun(m,lam,j)].add(v2(Ms[j]))
    single=sum(1 for v in by.values() if len(v)==1)
    print(f"  {label:55s}: {single}/{len(by)} keys single-valued")

def c4(lam):
    core,quot=core_quotient(lam,4)
    return (tuple(core), tuple(tuple(q) for q in quot))
def c2(lam):
    core,quot=core_quotient(lam,2)
    return (tuple(core), tuple(tuple(q) for q in quot))

det_test(lambda m,lam,j:(m,j), "by (m,j)")
det_test(lambda m,lam,j:(m,j,c4(lam)[1]), "by (m,j,4-quot)")
det_test(lambda m,lam,j:(m,j,c4(lam)), "by (m,j,4-core,4-quot)")
det_test(lambda m,lam,j:(m,j,c2(lam)), "by (m,j,2-core,2-quot)")
det_test(lambda m,lam,j:(m,j,v2(M_vector(lam,m)[0])), "by (m,j,v2(f))")

# ---- P4: v2(M_j) vs v2(C(m,j)) — is val(j) the clean object? ----
print(f"\n=== P4  is val(j)=j+2v2(C(m,j)M_j) cleaner than v2(M_j)? ===")
# test: does val(j) - val(0) depend only on (m,j,4-quot)?
by=defaultdict(set)
for m in range(1,MMAX+1):
    for lam in partitions(2*m):
        Ms=M_vector(lam,m)
        if Ms[0]==0: continue
        v0=2*v2(Ms[0])  # val(0)=0+2v2(C(m,0)M_0)=2v2(f)
        for j in range(m+1):
            if Ms[j]==0: continue
            val=j+2*v2(math.comb(m,j))+2*v2(Ms[j])
            by[(m,j,c4(lam)[1])].add(val-v0)
single=sum(1 for v in by.values() if len(v)==1)
print(f"  val(j)-val(0) by (m,j,4-quot): {single}/{len(by)} single-valued")

# ---- P5: M_j mod 2 from chi_b ----
print(f"\n=== P5  M_j mod 2 engine cross-check (m<={MMAX}) ===")
# M_j*2^j = sum_b C(j,b)(-1)^b chi_b. mod 2: depends on chi_b mod 2.
# verify M_j parity equals predicted from chi_b mod 2^(j+1).
ok=0;tot=0
for m in range(1,MMAX+1):
    for lam in partitions(2*m):
        chi=chi_b_vector(lam,m)
        Ms=M_vector(lam,m)
        for j in range(m+1):
            s=sum(math.comb(j,b)*((-1)**b)*chi[b] for b in range(j+1))
            assert s==Ms[j]*(2**j)
            tot+=1; ok+=1
print(f"  M_j*2^j = sum_b C(j,b)(-1)^b chi_b verified: {ok}/{tot}")
