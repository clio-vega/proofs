"""Q430 PART 9: the holdings gap is CLOSED.  HKKO 2301.13117 is on disk at
/home/clio/data/fetched/arxiv-2301.13117.pdf (+ .txt).  Verify sign_hkko against
eq (3.4) VERBATIM, by an independent re-implementation transcribed off the source.

SOURCE QUOTE, HKKO 2301.13117v3, Theorem 3.3, eq (3.4), p.11 (line 1042-1060 of the
text layer):
   "Let h and w be positive integers. For a partition lambda with length <= 2h and
    lambda_1 - lambda_{2h} <= w, we define
      c^{+-}_{2h,w}(lambda) =
         1,   if lambda_{2i-1} = lambda_{2i} for all 1 <= i <= h,
        +-1,  if lambda_1 - lambda_{2h} = w and lambda_{2i} = lambda_{2i+1}
                                                        for all 1 <= i <= h-1,
         0,   otherwise."
DICTIONARY: HKKO's h  ==  my k ;  HKKO's 2h  ==  my h = 2k.  So HKKO's
"1 <= i <= h" is my "1 <= i <= k", and the determinant size h is my k.
"""
import sys, itertools
sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
import sympy as sp
from korff import paths
t=sp.Symbol('t')

def sign_hkko(la,k,w):                       # the implementation under audit
    h=2*k
    if all(la[2*i]==la[2*i+1] for i in range(k)): return 1
    if la[0]-la[h-1]==w and all(la[2*i+1]==la[2*i+2] for i in range(k-1)): return -1
    return 0

def sign_from_source(la, k, w):
    """Transcribed off eq (3.4) using 1-INDEXED lambda, HKKO's h = k.
    Deliberately written in HKKO's own indices to avoid inheriting an off-by-one."""
    H = k                                     # HKKO's h
    lam = {i+1: la[i] for i in range(len(la))}        # 1-indexed
    def L(i): return lam.get(i, 0)
    # case 1: lambda_{2i-1} = lambda_{2i} for all 1 <= i <= H
    if all(L(2*i-1) == L(2*i) for i in range(1, H+1)): return 1
    # case 2: lambda_1 - lambda_{2H} = w  and  lambda_{2i} = lambda_{2i+1}, 1 <= i <= H-1
    if L(1) - L(2*H) == w and all(L(2*i) == L(2*i+1) for i in range(1, H)): return -1
    return 0

print("PART 9 -- auditing sign_hkko against HKKO eq (3.4), now that the source is held.")
print()
bad=0; tot=0; nz=0
# exhaustive over every endpoint reachable in every fibre we used, plus a brute-force sweep
cases=[]
for k in (1,2,3,4):
    for ell in (0,1,2):
        w=2*ell+2; h=2*k
        for n in range(1, min(w,5)+1):
            for P in paths(h,w,tuple([1]*n)):
                cases.append((tuple(P[-1]),k,w))
# brute force: all weakly decreasing la in Par(2h,w) for small h,w
for k in (1,2,3):
    h=2*k
    for w in (2,4,6):
        for la in itertools.product(range(0, w+2), repeat=h):
            if all(la[i]>=la[i+1] for i in range(h-1)) and la[0]-la[h-1]<=w:
                cases.append((la,k,w))
cases=list(dict.fromkeys(cases))
for (la,k,w) in cases:
    a=sign_hkko(list(la),k,w); b=sign_from_source(list(la),k,w)
    tot+=1
    if a!=0: nz+=1
    if a!=b:
        bad+=1
        if bad<=10: print(f"   MISMATCH la={la} k={k} w={w}: impl={a} source={b}")
print(f" endpoints + brute-force partitions compared: {tot}")
print(f" of these, impl returns a NONZERO sign on: {nz}   (non-vacuous: {nz>0})")
print(f" MISMATCHES: {bad}")
print(f" VERDICT: sign_hkko agrees with eq (3.4) on all {tot} cases -> {bad==0}")
print()
print(" NEGATIVE CONTROL -- a deliberately corrupted transcription must DISAGREE:")
def sign_corrupt(la,k,w):
    H=k; lam={i+1:la[i] for i in range(len(la))}
    def L(i): return lam.get(i,0)
    if all(L(2*i-1)==L(2*i) for i in range(1,H+1)): return 1
    if L(1)-L(2*H)==w and all(L(2*i)==L(2*i+1) for i in range(1,H+1)): return -1   # H-1 -> H
    return 0
d=sum(1 for (la,k,w) in cases if sign_hkko(list(la),k,w)!=sign_corrupt(list(la),k,w))
print(f"   'for all 1<=i<=h' instead of 'h-1' in case 2: disagreements = {d}  (control fires: {d>0})")
def sign_corrupt2(la,k,w):
    H=k; lam={i+1:la[i] for i in range(len(la))}
    def L(i): return lam.get(i,0)
    if all(L(2*i-1)==L(2*i) for i in range(1,H+1)): return 1
    if L(1)-L(2*H)==w-1 and all(L(2*i)==L(2*i+1) for i in range(1,H)): return -1   # w -> w-1
    return 0
d2=sum(1 for (la,k,w) in cases if sign_hkko(list(la),k,w)!=sign_corrupt2(list(la),k,w))
print(f"   'lambda_1-lambda_{{2h}} = w-1' instead of 'w': disagreements = {d2}  (control fires: {d2>0})")
print()
print("SENTINEL: q430f.py part 9 complete.")
