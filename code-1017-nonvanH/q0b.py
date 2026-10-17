"""The DIRECTION I need: k >= max(I)  =>  w_k absent from Atilde_{M,I}|_{q=0}.
(The iff I pre-registered is FALSE -- extra vanishing occurs. Separate the two directions.)
PLANT: the reverse implication 'w_k absent => k >= max(I)' must FAIL (it does, by above),
       which is what shows the two directions are genuinely different.
Also D: coverage of the stratum max(I) <= b."""
import sympy as sp, itertools
import q484_engine as E
from q484_engine import q, wvars

print("=== C': k >= max(I)  =>  w_k absent  (q=0, k not in M, M>=I) ===")
for N in (2,3,4):
    w=wvars(N); need=0; fail=0; rev_fail=0
    for M in E.subsets(N):
        for I in E.subsets(N):
            if not E.geq(M,I): continue
            A0 = sp.expand(sp.cancel(E.Atilde(N,w,M,I,X='C')).subs(q,0))
            if A0 == 0: continue
            mx = max(I) if I else 0
            for k in range(1,N+1):
                if k in M: continue
                present = A0.has(w[k-1])
                if k >= mx:
                    need += 1
                    if present:
                        fail += 1
                        print(f"  !! N={N} M={sorted(M)} I={sorted(I)} k={k}: present but k>=max(I)")
                else:
                    if not present: rev_fail += 1
    print(f"  N={N}: implication holds in {need-fail}/{need};  "
          f"reverse implication fails {rev_fail} times (so the two directions differ)")

print()
print("=== D: coverage of stratum max(I) <= b over the nonvanH domain (a<b, M>=I, {a,b} cap M = 0) ===")
for N in range(2,9):
    tot=cov=0
    for a,b in itertools.combinations(range(1,N+1),2):
        for M in E.subsets(N):
            if a in M or b in M: continue
            for I in E.subsets(N):
                if not E.geq(M,I): continue
                tot+=1
                if (max(I) if I else 0) <= b: cov+=1
    print(f"  N={N}: {cov}/{tot} = {100.0*cov/tot:.1f}%")
