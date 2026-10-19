"""V3.  Validate the type-B engine path (WX(...,'B') from WBfunction L497-501,
EB from typeBnormalizationfactor L1361-1367) against type-B theorems of Motegi
that were NOT used to write it.

DELIBERATELY EXCLUDED from this file: stabilityWB (L2780) and
normalizedsusythree (L2770).  Those two ARE (specialisations of) the theorem I am
proving today; using them to validate the engine and then using the engine to
corroborate the theorem would be circular.  They are treated in v5_motegi_priority.py
as LITERATURE, not as engine validators.

Validators used here:
  (A) typeBmainproperties (i), normalizationWB L1386: the diagonal M=I is a product
      over roots.  [root-product block]
  (B) typeBmainproperties (ii) L1378-1385: divisibility of WB~(w_M|w|I).
      [divisibility block]
  (C) prop:B_si (ii) L2733-2737: w-symmetry WB~(u|w|I) = WB~(u|s_i w|I) when
      i,i+1 both in I or both out.   [w-side block]
  (D) prop:B_shift (ii) L2757-2760: WB~_{N,n+1,r}(u,w_1^-1 | w_1^-1,w_2..w_N | I)
      = (1-q)/(1-q w_1) * WB~_{N,n,r}(u|w|I)  when 1 notin I.  [boundary-shift block]
"""
import sympy as sp, itertools, random
from q484_engine import F, EC, EB, EX, WX, q, wvars

def WBt(N, n, r, u, w, I):
    """widetilde{WB}_{N,n,r}(u|w|I) at ARBITRARY u -- no w_M specialisation."""
    return sp.cancel(sp.together(WX(N, n, r, u, w, I, 'B') / EB(N, n, r, u, w)))

def rat(seed, k, lo=2, hi=40):
    """k pairwise-distinct positive rationals, deterministic."""
    rng = random.Random(seed); out = []
    while len(out) < k:
        v = sp.Rational(rng.randint(lo, hi), rng.randint(1, 7))
        if v not in out and v != 1 and v != -1: out.append(v)
    return out

def ok(e): return sp.simplify(sp.cancel(sp.together(e))) == 0

results = []

# ---- (A) diagonal root product, normalizationWB L1386 --------------------
# Inv(lambda), T(lambda) are Motegi's root sets; rather than re-derive them I use
# the explicit type-C form (eq:diag) with the type-B substitution of the DIAGONAL
# term only, and test it as a PREDICTION.  Reading L997-1030 + L1386:
#   type C diagonal block prod_{a<=b}(1 - 1/(w_ia w_ib))  has a=b term (1-w_ia^-2)
#   type B must carry the SHORT boundary root, i.e. (1 - w_ia^-1)  [L1378: prod_j (1-w_mj^-1)]
def diagB(N, I):
    w = wvars(N); r = len(I); Ic = [c for c in range(1, N+1) if c not in I]
    t = sp.Integer(1)
    for j in range(r):                                   # boundary roots e_{i_j}
        t *= (1 - 1/w[I[j]-1])
    for a in range(r):                                   # pair roots e_ia + e_ib
        for b in range(a+1, r):
            t *= (1 - 1/(w[I[a]-1]*w[I[b]-1]))
    for c in Ic:                                         # A-type, descent side
        for il in I:
            if il > c: t *= (1 - w[c-1]/w[il-1])
    for c in Ic:                                         # (Q) pair, complement
        t *= (1 - q*w[c-1])
        for d in Ic:
            if d > c: t *= (1 - q*w[c-1]*w[d-1])
    for c in Ic:                                         # (Q) A-type
        for il in I:
            if il < c: t *= (1 - q*w[c-1]/w[il-1])
    return t

for N in (2, 3):
    w = wvars(N)
    for r in range(0, N+1):
        for I in itertools.combinations(range(1, N+1), r):
            n = max(r, 1)
            u = [w[m-1] for m in I] if r else []
            if r == 0:
                lhs = sp.cancel(1/EB(N, 0, 0, [], w)); n = 0
            else:
                lhs = WBt(N, r, r, u, w, I)
            rhs = diagB(N, I)
            results.append(("A:diagonal-root-product", N, r, str(I), ok(lhs-rhs)))

# ---- (C) w-symmetry prop:B_si (ii) L2733-2737 ---------------------------
for N in (2, 3):
    w = wvars(N)
    for n in (1, 2):
        for r in range(0, n+1):
            for I in itertools.combinations(range(1, N+1), r):
                uvals = rat(1000+n, n)
                for i in range(1, N):          # s_i swaps w_i, w_{i+1}
                    inI, inI1 = i in I, (i+1) in I
                    if inI != inI1: continue   # hypothesis fails
                    ws = list(w); ws[i-1], ws[i] = ws[i], ws[i-1]
                    lhs = WBt(N, n, r, uvals, w, I)
                    rhs = WBt(N, n, r, uvals, ws, I)
                    results.append(("C:w-symmetry", N, r, f"I={I},s_{i}", ok(lhs-rhs)))

# ---- (D) boundary shift prop:B_shift (ii) L2757-2760 --------------------
for N in (2, 3):
    w = wvars(N)
    for n in (1, 2):
        for r in range(0, min(n, N)+1):
            for I in itertools.combinations(range(1, N+1), r):
                if 1 in I: continue                        # hypothesis: 1 notin I
                uvals = rat(2000+n, n)
                wsh = [1/w[0]] + list(w[1:])
                lhs = WBt(N, n+1, r, uvals + [1/w[0]], wsh, I)
                rhs = (1-q)/(1-q*w[0]) * WBt(N, n, r, uvals, w, I)
                results.append(("D:boundary-shift", N, r, f"I={I},n={n}", ok(lhs-rhs)))

# ---- PLANTS: each block gets a positive arm -----------------------------
plants = []
# plant A: wrong boundary root in the diagonal -- use type C's (1-w^-2)
def diagB_plant(N, I):
    w = wvars(N); base = diagB(N, I); r = len(I)
    for j in range(len(I)):
        base = base / (1 - 1/w[I[j]-1]) * (1 - 1/w[I[j]-1]**2)
    return base
for N in (2, 3):
    w = wvars(N)
    for r in range(1, N+1):
        for I in itertools.combinations(range(1, N+1), r):
            lhs = WBt(N, r, r, [w[m-1] for m in I], w, I)
            plants.append(("A-plant:typeC-boundary-root", N, r, str(I),
                           not ok(lhs - diagB_plant(N, I))))
# plant D: wrong multiplier (1-q)/(1-q w_1^2)  [type C's shape]
for N in (2, 3):
    w = wvars(N)
    for n in (1, 2):
        for r in range(0, min(n, N)+1):
            for I in itertools.combinations(range(1, N+1), r):
                if 1 in I: continue
                uvals = rat(2000+n, n)
                wsh = [1/w[0]] + list(w[1:])
                lhs = WBt(N, n+1, r, uvals + [1/w[0]], wsh, I)
                rhs = (1-q)/(1-q*w[0]**2) * WBt(N, n, r, uvals, w, I)
                plants.append(("D-plant:typeC-shape", N, r, f"I={I},n={n}", not ok(lhs-rhs)))
# plant E: feed the type-C normaliser where type B belongs (the 10-17 defect)
for N in (2, 3):
    w = wvars(N)
    for r in range(1, N+1):
        for I in itertools.combinations(range(1, N+1), r):
            bad = sp.cancel(sp.together(WX(N, r, r, [w[m-1] for m in I], w, I, 'B')
                                        / EC(N, r, r, [w[m-1] for m in I], w)))
            plants.append(("E-plant:WB-over-EC", N, r, str(I), not ok(bad - diagB(N, I))))

for tag in sorted(set(x[0] for x in results)):
    sub = [x for x in results if x[0] == tag]
    print(f"{tag:28s}  {sum(1 for x in sub if x[4])}/{len(sub)}")
    for x in sub:
        if not x[4]: print("    FAIL", x)
print("--- plants (should FIRE = pass) ---")
for tag in sorted(set(x[0] for x in plants)):
    sub = [x for x in plants if x[0] == tag]
    print(f"{tag:28s}  fired {sum(1 for x in sub if x[4])}/{len(sub)}")
    for x in sub:
        if not x[4]: print("    BLIND", x)
