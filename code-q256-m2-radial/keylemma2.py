"""Sharper KEY LEMMA: beta PF2 (log-concave, interval support, finite) =>
gamma(y)=sum_{r>=y} beta(r) is log-concave.

Proof: with q = beta(y)/beta(y-1) < 1, log-concavity makes ratios non-increasing,
so beta(y+k) <= beta(y) q^k, hence gamma(y) <= beta(y)/(1-q) = beta(y-1)beta(y)/Delta.

REFUSAL CONTROL: drop log-concavity (keep interval support) and the instrument MUST fire."""
import random, itertools
def gamma_of(beta):
    g={}; acc=0
    for k in range(len(beta)-1,-1,-1):
        acc+=beta[k]; g[k]=acc
    return g
def tail_bad(beta):
    g=gamma_of(beta); S=sum(beta); n=len(beta)
    get=lambda y: (S if y<0 else g.get(y,0))
    return [(y,get(y),get(y-1),get(y+1)) for y in range(-1,n+2)
            if get(y)**2 < get(y-1)*get(y+1)]
def interval_support(b):
    nz=[i for i,v in enumerate(b) if v>0]
    return (not nz) or nz[-1]-nz[0]+1==len(nz)
def is_lc(b):
    return all(b[i]**2>=b[i-1]*b[i+1] for i in range(1,len(b)-1))

print("=== CONTROL: can the instrument refuse at all? ===")
for probe in [[1,0,5],[10,1,10],[1,0,1],[2,1,3]]:
    print(f"   beta={probe} intervalsupp={interval_support(probe)} LC={is_lc(probe)} -> bad={tail_bad(probe)[:2]}")

random.seed(7)
print()
print("=== (A) CONFIRM: beta PF2 (log-concave + interval support) ===")
tr=0; f=0
for _ in range(600000):
    L=random.randint(2,8)
    b=[random.randint(0,12) for _ in range(L)]
    if sum(b)==0 or not interval_support(b) or not is_lc(b): continue
    tr+=1
    bad=tail_bad(b)
    if bad:
        f+=1
        if f<=3: print("   FAIL",b,bad)
print(f"   trials={tr} failures={f}")

print()
print("=== (B) REFUSE: interval support but NOT log-concave  (must find failures) ===")
tr=0; f=0; shown=0
for _ in range(600000):
    L=random.randint(3,8)
    b=[random.randint(0,12) for _ in range(L)]
    if sum(b)==0 or not interval_support(b) or is_lc(b): continue
    tr+=1
    bad=tail_bad(b)
    if bad:
        f+=1
        if shown<3: print("   refused:",b,bad[:1]); shown+=1
print(f"   trials={tr} failures={f}  ({100*f/max(tr,1):.1f}% refused)")
assert f>0, "INSTRUMENT BLIND"

print()
print("=== (C) exhaustive: every PF2 b in [0,7]^L, L<=6 ===")
tr=0; f=0
for L in range(1,7):
    for b in itertools.product(range(8), repeat=L):
        b=list(b)
        if sum(b)==0 or not interval_support(b) or not is_lc(b): continue
        tr+=1
        if tail_bad(b): f+=1; print("   FAIL",b) if f<4 else None
print(f"   exhaustive PF2 trials={tr} failures={f}")

print()
print("=== (D) long/extreme PF2 sequences (geometric, near-equality cases) ===")
tr=0; f=0
cases=[[2**k for k in range(12)],[2**(11-k) for k in range(12)],
       [1]*30,[3**(9-k) for k in range(10)],
       [1,2,4,8,16,32,64,32,16,8,4,2,1]]
for b in cases:
    if not (interval_support(b) and is_lc(b)): print("   skip (not PF2)",b[:5]); continue
    tr+=1; bad=tail_bad(b)
    print(f"   len={len(b)} PF2 ok, failures={len(bad)} {bad[:1]}")
    if bad: f+=1
print(f"   extreme trials={tr} failures={f}")
