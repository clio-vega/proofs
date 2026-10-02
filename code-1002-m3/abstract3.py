"""Candidate abstract theorem:
  concentric box splines whose half-width MULTIPLICITY PROFILE M is the positive part of a
  concave function  ==>  the sum is PF2.
At m=1 this is exactly the proved m=2 cylindric architecture (M PF2 -> tail LC).
Test at m=2,3 and push past the birth range."""
import gen, layers
from itertools import combinations_with_replacement, product
def spl(w): return gen.conv_intervals(list(w))
def half(w): return (len(spl(w))-1)/2
def addc(cos):
    H=max(len(c)-1 for c in cos); out=[0]*(H+1)
    for c in cos:
        p=H-(len(c)-1)
        if p%2: return None
        for j,v in enumerate(c): out[p//2+j]+=v
    return out
def Mprofile(hs):
    lo,hi=min(hs),max(hs)
    return [sum(1 for h in hs if abs(h-(lo+k))<1e-9) for k in range(int(hi-lo)+1)]
for m,wmax,Ns in [(1,11,range(2,9)),(2,6,range(2,7)),(3,5,range(2,6))]:
    WS=[w for w in product(range(1,wmax+1),repeat=m)]
    for N in Ns:
        bad=[];tot=0
        for ms in combinations_with_replacement(WS,N):
            hs=[half(w) for w in ms]
            M=Mprofile(hs)
            if not layers.is_concave_pospart(M): continue
            g=addc([spl(w) for w in ms])
            if g is None: continue
            tot+=1
            if not gen.is_pf2(g): bad.append((ms,M,g))
        print(f"  m={m} N={N}: {tot} multisets with M concave; sum NOT PF2: {len(bad)}  e.g. {bad[:2]}", flush=True)
