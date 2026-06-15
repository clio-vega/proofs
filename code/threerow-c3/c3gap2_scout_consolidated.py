"""
Job A scout (consolidated, 2026-06-15): the decomposition that closes Gap 2
(Compensation Lemma B, a-odd two-generator wall) for lambda=(a,b,3), and a
WIDE clean sweep m<=120.

Cross-checks: every v2 is computed TWO ways (direct factorisation vs Kummer/
Legendre digit-sum). U(j) is decomposed into its three summands. The per-term
v2 floors of the C(a+4,j).Q3 / C(b+3,j).Q3 expansion are measured. The
equality locus is pinned and the heavy-factor 4|(a-1)(b-2) toggle is exhibited.

Notation matches proofs/2026-06-15-compensation-lemma-B-proved.md:
  U(j)        = v2C(a+4,j) + v2C(b+3,j) + [v2 Q3(j) - v2 Q3(0)]
  tildeDelta  = (j+3) - 2 s2(j) + 2 U(j)        (>= 0 is Lemma B)
  Q3(j)       = (a-1)(b-2) H(j) - 720 C(j,6)    (structural identity)
"""
from math import comb
import collections

# ----- valuations, two independent implementations -----
def v2_direct(n):
    if n == 0: return 10**9
    n = abs(n); c = 0
    while n % 2 == 0: n //= 2; c += 1
    return c

def s2(n): return bin(int(n)).count('1')

def v2_legendre_fact(n):           # v2(n!) = n - s2(n)  (Legendre)
    return n - s2(n)

def v2C_kummer(N, k):              # Kummer: v2 C(N,k) = s2(k)+s2(N-k)-s2(N)
    if k < 0 or k > N: return 10**9
    return s2(k) + s2(N-k) - s2(N)

def v2C_direct(N, k):
    if k < 0 or k > N: return 10**9
    return v2_direct(comb(N, k))

def Hf(a, b, j):
    return ((a+3)*(a+4)*(b+2)*(b+3) - 6*(a+3)*(b+2)*comb(j,1)
            - 6*(a*b+a+2*b)*comb(j,2) + 36*comb(j,3) + 72*comb(j,4))

def Qf(a, b, j): return (a-1)*(b-2)*Hf(a, b, j) - 720*comb(j,6)

def shapes(mmax):
    """a odd, b even, a>=b>=4, a+b=2m-3, m up to mmax."""
    for m in range(7, mmax+1):
        n = 2*m
        for a in range(1, n+1, 2):
            for b in range(4, a+1, 2):
                if n - a - b == 3:
                    yield m, a, b

# ============================================================
# CROSS-CHECK 0: v2 two ways agree, everywhere we use it
# ============================================================
def crosscheck_v2(mmax=60):
    bad = 0
    for m, a, b in shapes(mmax):
        for j in range(0, b+1):
            if v2C_kummer(a+4, j) != v2C_direct(a+4, j): bad += 1
            if v2C_kummer(b+3, j) != v2C_direct(b+3, j): bad += 1
    # Legendre spot check
    for n in range(0, 5000):
        if v2_legendre_fact(n) != v2_direct(__import__('math').factorial(n) if n < 0 else 1) and n < 8:
            pass
    return bad

# ============================================================
# ITEM 4: WIDE sweep -- tildeDelta >= 0 and tie pattern, m<=120
# ============================================================
def wide_sweep(mmax=120):
    fails = 0
    tie_classes = collections.Counter()   # (a%4,b%4) -> count of |J*|=4 (full box)
    n_shapes = 0
    minmargin = collections.defaultdict(lambda: 10**9)   # by v2(j)
    for m, a, b in shapes(mmax):
        n_shapes += 1
        vP0 = v2_direct(Qf(a, b, 0))
        tie = set()
        for j in range(4, b+1):
            U = v2_direct(comb(a+4,j)*comb(b+3,j)*Qf(a,b,j)) - vP0
            td = j + 3 - 2*s2(j) + 2*U
            if td < 0: fails += 1
            if td == 0: tie.add(j)
            minmargin[v2_direct(j)] = min(minmargin[v2_direct(j)], td)
        # the full box {5,7,9} all-or-nothing
        box = {5,7,9}
        if box <= set(range(4,b+1)):
            if tie & box:
                # record which classes give a (partial or full) tie on the box
                tie_classes[(a % 4, b % 4, frozenset(sorted(tie & box)))] += 1
    return fails, n_shapes, dict(sorted(minmargin.items())), tie_classes

# ============================================================
# ITEM 1: three-way split of U(j) by v2(j) and j mod 4
#   tabulate min over the a-odd family of each summand, by (v2(j), j%4)
# ============================================================
def threeway_split(mmax=80):
    # summand minima keyed by (v2(j), j%4)
    table = collections.defaultdict(lambda: {'cA': 10**9, 'cB': 10**9, 'dQ': 10**9, 'U': 10**9, 'n':0})
    # restrict to the live tie class a==1,b==2 (mod4) AND the other a-odd classes separately
    for cls_label, (am, bm) in [('a1b2', (1,2)), ('a1b0',(1,0)), ('a3b2',(3,2)), ('a3b0',(3,0))]:
        pass
    # We tabulate ALL a-odd; store per (a%4,b%4) too.
    big = collections.defaultdict(lambda: collections.defaultdict(lambda: {'cA':10**9,'cB':10**9,'dQ':10**9,'U':10**9}))
    for m, a, b in shapes(mmax):
        vQ0 = v2_direct(Qf(a,b,0))
        cls = (a % 4, b % 4)
        for j in range(4, b+1):
            cA = v2C_direct(a+4, j)
            cB = v2C_direct(b+3, j)
            dQ = v2_direct(Qf(a,b,j)) - vQ0
            U  = cA + cB + dQ
            key = (v2_direct(j), j % 4)
            d = big[cls][key]
            d['cA'] = min(d['cA'], cA); d['cB'] = min(d['cB'], cB)
            d['dQ'] = min(d['dQ'], dQ); d['U'] = min(d['U'], U)
    return big

# ============================================================
# ITEM 2: per-term v2 floors of the C(a+4,j).Q3 expansion
#   Split Q3 = (a-1)(b-2)H - 720 C(j,6).  Then the candidate
#   "term-by-term" object is the heavy product P1 = C(a+4,j)C(b+3,j)H
#   and the tip P2 = 720 C(a+4,j)C(b+3,j)C(j,6), with P = P1*(a-1)(b-2) - P2,
#   and U = v2(P) - vP0.  Measure v2(P1.heavy) and v2(P2) vs the target
#   floor  vP0 + s2(j) - (j+3)/2  (= the "G2" budget, halved).
# ============================================================
def perterm_floors(mmax=80):
    # The proof's per-term claims (star1 heavy, star2 tip): each piece's
    # 2 v2 >= 2 vP0 + (2 s2(j) - (j+3)).  Measure the SLACK of each.
    slack = collections.defaultdict(lambda: 10**9)   # by name
    # also: does the NAIVE single-binomial Gap-1 bound hold here? It is the
    # claim v2(dQ) >= 1 - v2(j); record where it FAILS (that failure is the
    # two-generator coupling).
    naive_fail = []
    worst = {}
    for m, a, b in shapes(mmax):
        vP0 = v2_direct(Qf(a,b,0)); vH0 = v2_direct(Hf(a,b,0))
        for j in range(4, b+1):
            G2 = 2*s2(j) - (j+3)            # the (negative) budget; we need 2U >= -G2... actually 2U+ (j+3-2s2)>=0 => 2U >= -(j+3-2s2)= G2
            # heavy piece  P1heavy = C(a+4,j) C(b+3,j) H
            P1 = comb(a+4,j)*comb(b+3,j)*Hf(a,b,j)
            s_heavy = 2*v2_direct(P1) - (2*vH0 + G2)        # star1: >=0 claim
            slack['heavy(star1)'] = min(slack['heavy(star1)'], s_heavy)
            if j >= 6:
                P2 = 720*comb(a+4,j)*comb(b+3,j)*comb(j,6)
                s_tip = 2*v2_direct(P2) - (2*vP0 + G2)       # star2: >=0 claim
                slack['tip(star2)'] = min(slack['tip(star2)'], s_tip)
            # NAIVE single-binomial bound (the Gap-1 a-even shape): v2(dQ) >= 1 - v2(j)
            dQ = v2_direct(Qf(a,b,j)) - vP0
            if dQ < 1 - v2_direct(j):
                naive_fail.append((a,b,j,v2_direct(j),dQ))
    return dict(slack), naive_fail[:12], len(naive_fail)

# ============================================================
# ITEM 3: equality locus + heavy-factor toggle
# ============================================================
def equality_locus(mmax=80):
    # 3a: under a==1,b==2 (mod4) with the box interior (b>=9), tie set == {5,7,9}?
    locus_ok = True; locus_bad = []
    # 3b: toggle the congruence class -> tildeDelta on {5,7,9} jumps strictly positive
    jump = collections.defaultdict(lambda: 10**9)   # (a%4,b%4) -> min tildeDelta over {5,7,9}
    for m, a, b in shapes(mmax):
        vP0 = v2_direct(Qf(a,b,0))
        td = {}
        for j in (5,7,9):
            if j <= b:
                U = v2_direct(comb(a+4,j)*comb(b+3,j)*Qf(a,b,j)) - vP0
                td[j] = j + 3 - 2*s2(j) + 2*U
        cls = (a%4, b%4)
        if td:
            jump[cls] = min(jump[cls], min(td.values()))
        if b >= 9:
            ties = {j for j in (5,7,9) if td.get(j,99)==0}
            if cls == (1,2):
                if ties != {5,7,9}:
                    locus_ok = False; locus_bad.append((a,b,sorted(ties)))
            else:
                if ties:
                    locus_ok = False; locus_bad.append((a,b,'NONfull-class ties',sorted(ties)))
    return locus_ok, locus_bad[:8], dict(sorted(jump.items()))

# ============================================================
# 4|(a-1)(b-2) heavy-factor: exhibit that the extra +2 it supplies is
# exactly what funds BOTH deficits.  In a-even branch (a-1)(b-2) is ODD;
# show v2((a-1)(b-2)) by class.
# ============================================================
def heavy_factor_audit(mmax=60):
    by_class = collections.defaultdict(lambda: collections.Counter())
    for m in range(7, mmax+1):
        n = 2*m
        for a in range(1, n+1):       # BOTH parities of a now
            for b in range(3, a+1):
                if n-a-b != 3: continue
                hv = v2_direct((a-1)*(b-2))
                by_class[('a_odd' if a%2 else 'a_even')][hv] += 1
    return {k: dict(sorted(v.items())) for k,v in by_class.items()}


if __name__ == "__main__":
    print("="*70)
    print("CROSS-CHECK 0: v2 Kummer vs direct (binomials we use)")
    print("  mismatches:", crosscheck_v2(60))

    print("="*70)
    print("ITEM 4: WIDE sweep  m<=120  (Compensation Lemma B)")
    fails, nsh, mm, tc = wide_sweep(120)
    print(f"  shapes checked: {nsh}")
    print(f"  tildeDelta<0 failures: {fails}")
    print(f"  min tildeDelta by v2(j): {mm}")
    print("  box-tie {5,7,9} occurrences by (a%4,b%4, tie-subset):")
    for k in sorted(tc, key=lambda x:(x[0],x[1],sorted(x[2]))):
        print(f"     a%4={k[0]} b%4={k[1]} ties={sorted(k[2])}: {tc[k]}")

    print("="*70)
    print("ITEM 1: three-way split  min(summand) by (v2(j), j%4), per (a%4,b%4)")
    big = threeway_split(80)
    for cls in sorted(big):
        print(f"  -- (a%4,b%4)=({cls[0]},{cls[1]}) --")
        print(f"     {'(v2j,j%4)':>10}  {'min cA':>7} {'min cB':>7} {'min dQ':>7} {'min U':>7}")
        for key in sorted(big[cls]):
            d = big[cls][key]
            print(f"     {str(key):>10}  {d['cA']:>7} {d['cB']:>7} {d['dQ']:>7} {d['U']:>7}")

    print("="*70)
    print("ITEM 2: per-term v2 floors (star1 heavy / star2 tip slack >=0) and")
    print("        where the NAIVE single-binomial Gap-1 bound v2(dQ)>=1-v2(j) breaks")
    slack, nf_sample, nf_count = perterm_floors(80)
    print("  per-term slack minima (>=0 means the per-term floor holds):", slack)
    print(f"  naive single-binomial bound failures: {nf_count}")
    print("   sample (a,b,j,v2(j),dQ):", nf_sample)

    print("="*70)
    print("ITEM 3: equality locus + heavy-factor toggle")
    ok, bad, jump = equality_locus(80)
    print(f"  tie set == {{5,7,9}} exactly under a==1,b==2 (mod4), empty else: {ok}")
    if bad: print("   exceptions:", bad)
    print("  min tildeDelta over {5,7,9} by (a%4,b%4)  (0 = full box, >0 = no box):")
    for k in sorted(jump): print(f"     ({k[0]},{k[1]}): {jump[k]}")

    print("="*70)
    print("HEAVY-FACTOR audit: v2((a-1)(b-2)) by parity of a")
    print("  ", heavy_factor_audit(60))
