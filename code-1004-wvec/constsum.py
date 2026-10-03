"""Theorem R's load-bearing step: does the symmetric exchange axiom FORCE a constant
coordinate sum?  I could not close this from the local axiom alone (the all-(+-1) case
cycles), so test it exhaustively before relying on it."""
import mconv, itertools
tot=0; viol=[]
for m,rng in [(2,range(4)),(2,range(5)),(3,range(2)),(4,range(2)),(3,range(3))]:
    pts=[p for p in itertools.product(rng,repeat=m)]
    N=len(pts)
    if N>18:
        print(f" m={m} grid {list(rng)}^{m}: {N} pts -- SKIPPED"); continue
    n_mc=n_mc_multisum=0
    for mask in range(1,1<<N):
        S=[pts[k] for k in range(N) if mask>>k & 1]
        if mconv.exch_M(S)[0]:
            n_mc+=1
            if len({sum(p) for p in S})>1:
                n_mc_multisum+=1; viol.append((m,S))
    tot+=n_mc
    print(f" m={m} grid {list(rng)}^{m}: M-convex subsets {n_mc}; of those with NON-constant "
          f"coordinate sum: {n_mc_multisum}")
print(f" total M-convex subsets found {tot}; violations {len(viol)}")
for v in viol[:3]: print("  VIOLATION",v)
