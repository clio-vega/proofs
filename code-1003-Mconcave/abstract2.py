"""EXTEND the abstract search past its birth range (size<=3 gave 0 failures -- the tell)."""
import gen, layers, random
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
def gapfree(hs):
    s=sorted(set(hs)); return all(abs(s[i+1]-s[i]-1)<1e-9 for i in range(len(s)-1))

print("exhaustive, m=1, widths<=9, multiset size up to 6 (ALL gap-free)")
WS=[(w,) for w in range(1,10)]
for N in range(2,7):
    bad=[]; tot=0
    for ms in combinations_with_replacement(WS,N):
        hs=[half(w) for w in ms]
        if not gapfree(hs): continue
        g=addc([spl(w) for w in ms])
        if g is None: continue
        tot+=1
        if not gen.is_pf2(g): bad.append((ms,g))
    print(f"  N={N}: {tot} gap-free multisets, G not PF2: {len(bad)}  smallest e.g. {bad[:2]}")
print("\nexhaustive, m=2, widths<=5, size 4,5")
WS=[w for w in product(range(1,6),repeat=2)]
for N in (4,5):
    bad=[]; tot=0
    for ms in combinations_with_replacement(WS,N):
        hs=[half(w) for w in ms]
        if not gapfree(hs): continue
        g=addc([spl(w) for w in ms])
        if g is None: continue
        tot+=1
        if not gen.is_pf2(g): bad.append((ms,g))
    print(f"  N={N}: {tot} gap-free, G not PF2: {len(bad)}  e.g. {bad[:2]}")
