"""Lemma P (parity): on ANY slice, sum_i w_i(nu) = G+m-sum_i|y_i| and sum|y_i| = sigma (mod 2),
so all width vectors of one slice have coordinate sum of ONE fixed parity.
=> distinct half-widths differ by >=1, i.e. coordinate sums differ by >=2: a GAP.
Test (i) the parity statement, (ii) whether M^nat-convex sets can have gapped sum-sets."""
import gen, tent, mconv, itertools, random
from collections import Counter

print("=== (i) Lemma P on real slices ===")
st=Counter()
for m in range(2,7):
  for n in range(m,m+6):
    for (mu,lam) in gen.pairs(n,m,7):
        d=sum(lam)-sum(mu); G=n-m
        for b in range(0,d+1):
            target=sum(mu)+b; sig=None; par=set(); nw=0
            for nu in gen.box(mu,n,m):
                if sum(nu)!=target: continue
                lr=gen.LR(nu,lam,n,m); w=[R-L+1 for (L,R) in lr]
                if any(x<=0 for x in w): continue
                nw+=1
                y=[nu[i]-tent.lam_at(lam,i-1,n,m)-1 for i in range(m)]
                sig=sum(y); par.add(sum(w)%2)
                # the identity itself
                if sum(w)==G+m-sum(abs(t) for t in y): st['id_ok']+=1
                else: st['id_BAD']+=1
            if nw==0: continue
            st['slices']+=1
            if len(par)==1: st['one_parity']+=1
            else: st['MULTI_PARITY']+=1
            if len(par)==1 and (par.pop()==(G+m+sig)%2): st['parity_predicted']+=1
            else: st['parity_MISPREDICTED']+=1
print(" ",dict(st))

print("\n=== (ii) can an M^nat-convex set have a GAP in its set of coordinate sums? ===")
print("   exhaustive search over all subsets of small grids")
found=[]
for m,rng in [(2,3),(2,4),(3,3)]:
    pts=[p for p in itertools.product(range(rng),repeat=m)]
    N=len(pts)
    cnt=0
    for mask in range(1,1<<N):
        S=[pts[k] for k in range(N) if mask>>k & 1]
        sums=sorted({sum(p) for p in S})
        if len(sums)<2 or sums==list(range(sums[0],sums[-1]+1)): continue  # need a gap
        cnt+=1
        if mconv.exch_Mnat(S)[0]:
            found.append((m,rng,S,sums))
    print(f"   m={m} grid [0,{rng-1}]^{m}: {cnt} gapped subsets tested, "
          f"{sum(1 for f in found if f[0]==m and f[1]==rng)} of them M^nat-convex")
if found:
    print("   COUNTEREXAMPLES to 'M^nat => sums form an interval':")
    for f in found[:5]: print("    ",f)
else:
    print("   NONE: over every gapped subset of these grids, M^nat-convexity FAILS.")

print("\n=== (iii) positive control: M^nat-convex sets DO exist in these grids (test not vacuous) ===")
for m,rng in [(2,3),(3,3)]:
    pts=[p for p in itertools.product(range(rng),repeat=m)]
    N=len(pts); good=0; tot=0
    for mask in range(1,1<<N):
        S=[pts[k] for k in range(N) if mask>>k & 1]
        if len({sum(p) for p in S})<2: continue
        tot+=1
        if mconv.exch_Mnat(S)[0]: good+=1
    print(f"   m={m}: {good}/{tot} subsets with >=2 distinct sums are M^nat-convex")
