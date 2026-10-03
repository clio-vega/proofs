"""Can an M^natural-convex set have a GAP in its set of coordinate sums?
Exhaustive over all subsets of grids small enough to enumerate exactly."""
import mconv, itertools
GRIDS=[(2,range(4)),(2,range(5)),(3,range(2)),(4,range(2)),(3,range(3))]
for m,rng in GRIDS:
    pts=[p for p in itertools.product(rng,repeat=m)]
    N=len(pts)
    if N>18:
        print(f" m={m} grid {list(rng)}^{m}: {N} points, 2^{N} subsets -- SKIPPED (too large)"); continue
    gapped=mnat_gapped=mnat_any=multi=0
    wit=None
    for mask in range(1,1<<N):
        S=[pts[k] for k in range(N) if mask>>k & 1]
        sums=sorted({sum(p) for p in S})
        if len(sums)<2: continue
        multi+=1
        isgap = sums!=list(range(sums[0],sums[-1]+1))
        mn = mconv.exch_Mnat(S)[0]
        if mn: mnat_any+=1
        if isgap:
            gapped+=1
            if mn:
                mnat_gapped+=1
                if wit is None: wit=(S,sums)
    print(f" m={m} grid {list(rng)}^{m}: subsets with >=2 distinct sums {multi}; of these "
          f"M^nat-convex {mnat_any} (positive control, nonzero => test not vacuous); "
          f"gapped {gapped}; gapped AND M^nat-convex {mnat_gapped}")
    if wit: print("    COUNTEREXAMPLE:",wit)
