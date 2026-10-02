"""Is 'a multiset of concentric box splines with GAP-FREE half-widths' enough for
   (i) G = sum PF2, and (ii) beta = radial decrement log-concave?

If NOT, the counterexample locates exactly what more a general-m proof must use
beyond concentricity (Thm 1) + the tent lemma (Thm 2).
Exhaustive over small multisets; also the two-element certificate that gap-freeness kills.
"""
import gen, layers
from itertools import combinations_with_replacement, product

def spl(w):
    return gen.conv_intervals(list(w))

def half(w):
    c = spl(w)
    return (len(c)-1)/2

def add_concentric(cos):
    """add coefficient lists all centred at a common point"""
    H = max((len(c)-1) for c in cos)
    out=[0]*(H+1)
    for c in cos:
        pad=(H-(len(c)-1))//2
        if (H-(len(c)-1))%2: return None   # parity clash: cannot be concentric on one lattice
        for j,v in enumerate(c): out[pad+j]+=v
    return out

def gapfree(hs):
    s=sorted(set(hs))
    return all(abs(s[i+1]-s[i]-1)<1e-9 for i in range(len(s)-1))

print("baseline: the 2-element certificate WITHOUT gap-freeness")
cos=[spl((1,)), spl((5,))]
g=add_concentric(cos); print("  w=(1) and w=(5): halves",[half((1,)),half((5,))],"G =",g,"PF2 =",gen.is_pf2(g))

print("\nexhaustive search over multisets of concentric splines with GAP-FREE halves")
for m in (1,2,3):
    WS=[w for w in product(range(1,6),repeat=m)]
    bad_pf2=[]; bad_beta=[]; tested=0
    for N in (2,3):
        for ms in combinations_with_replacement(WS, N):
            hs=[half(w) for w in ms]
            if not gapfree(hs): continue
            cos=[spl(w) for w in ms]
            g=add_concentric(cos)
            if g is None: continue
            tested+=1
            if not gen.is_pf2(g): bad_pf2.append((ms,g))
            bet=layers.beta_of(layers.radial(g))
            if min(bet)<0 or not layers.is_lc(bet): bad_beta.append((ms,g,bet))
    print(f"  m={m}: {tested} gap-free multisets of size<=4")
    print(f"     G not PF2: {len(bad_pf2)}   e.g. {bad_pf2[:2]}")
    print(f"     beta not LC: {len(bad_beta)}  e.g. {bad_beta[:2]}")
