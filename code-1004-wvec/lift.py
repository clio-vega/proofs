"""Lemma G needs: W M^nat-convex => {sum w} is an integer interval.
Two-line proof if I may use the LIFT definition (Murota): W is M^nat-convex iff
  tilde W = {(-sum w, w) : w in W}  is M-convex in Z^{m+1}.
Verify that the lift definition and the local M^nat exchange axiom agree, exhaustively."""
import mconv, itertools

def lift(W):
    return [(-sum(w),)+tuple(w) for w in W]

for m,rng in [(2,range(4)),(3,range(2)),(4,range(2)),(2,range(3))]:
    pts=[p for p in itertools.product(rng,repeat=m)]
    N=len(pts); agree=dis=0; wit=None
    for mask in range(1,1<<N):
        S=[pts[k] for k in range(N) if mask>>k & 1]
        a=mconv.exch_Mnat(S)[0]; b=mconv.exch_M(lift(S))[0]
        if a==b: agree+=1
        else:
            dis+=1
            if wit is None: wit=(S,a,b)
    print(f" m={m} grid {list(rng)}^{m}: {agree} subsets agree, {dis} disagree", wit or "")
