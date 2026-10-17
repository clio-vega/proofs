"""C: the q=0 structure claim, and D: coverage of the stratum max(I) <= b.

C. CLAIM (derived by hand from Motegi's symmetrisation formula Thm 3.3, L491-504):
   at q=0 every w-block collapses --
       j<=r:  prod_{k>i_j}(u_j - q w_k) -> u_j^{N-i_j},   prod_{k<i_j}(u_j - w_k)  SURVIVES
       j>r :  prod_{k=1}^N (u_j - q w_k) -> u_j^N
   and  E^C|_{q=0} = 1.  Hence Atilde_{M,I}|_{q=0} depends on w_k for k not in M
   ONLY IF k < max(I).  (Pre-registered BEFORE measuring.)
   PLANT C1: the same claim with "k <= max(I)" must FAIL somewhere (else the test is blind).
"""
import sympy as sp, itertools
import q484_engine as E
from q484_engine import q, wvars

print("=== C0: E^C|_{q=0} = 1 ? ===")
for N in (2,3,4):
    w=wvars(N)
    for M in E.subsets(N):
        n=len(M); u=[w[m-1] for m in M]
        for r in range(0,n+1):
            v = sp.simplify(E.EC(N,n,r,u,w).subs(q,0))
            if v != 1: print(f"  !! N={N} M={sorted(M)} r={r}: E^C|_q=0 = {v}")
    print(f"  N={N}: checked")

print()
print("=== C: at q=0, does Atilde_{M,I} involve w_k (k not in M) only for k < max(I)? ===")
for N in (2,3,4):
    w=wvars(N)
    nclaim=nplant=0; claim_ok=plant_ok=True
    for M in E.subsets(N):
        for I in E.subsets(N):
            if not E.geq(M,I): continue
            A0 = sp.expand(sp.cancel(E.Atilde(N,w,M,I,X='C')).subs(q,0))
            if A0 == 0: continue
            mx = max(I) if I else 0
            for k in range(1,N+1):
                if k in M: continue
                present = A0.has(w[k-1])
                pred_claim = (k < mx)
                pred_plant = (k <= mx)
                if present != pred_claim:
                    claim_ok=False
                    print(f"  CLAIM FAILS: N={N} M={sorted(M)} I={sorted(I)} k={k} "
                          f"present={present} predicted={pred_claim}")
                else: nclaim+=1
                if present != pred_plant: nplant+=1; plant_ok=False
    print(f"  N={N}: claim agreements {nclaim}, claim holds = {claim_ok};  "
          f"PLANT C1 (k<=max I) disagreements = {nplant} (must be >0 to discriminate)")

print()
print("=== D: coverage of the stratum max(I) <= b, over the nonvanH domain ===")
print("    domain: a<b, M>=I, {a,b} cap M = empty.  (combinatorial only)")
for N in range(2,8):
    tot=cov=0
    for a,b in itertools.combinations(range(1,N+1),2):
        for M in E.subsets(N):
            if a in M or b in M: continue
            for I in E.subsets(N):
                if not E.geq(M,I): continue
                tot+=1
                if (max(I) if I else 0) <= b: cov+=1
    print(f"  N={N}: {cov}/{tot} = {100.0*cov/tot:.1f}%")
