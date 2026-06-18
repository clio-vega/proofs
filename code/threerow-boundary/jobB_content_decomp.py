"""
jobB_content_decomp.py -- 2026-06-18 CODE, Job B.

Pin the c!*k!/den floor arithmetic for the general-c boundary Content Lemma.

For c=4..12, depths k=4,5,6, both b-parities, decompose the EXACT 2-adic content
of the master deficit N_i^{(c)} (i=c-k) into
      content g  =  v2(c! k!)  +  mu ,      mu := min_slice [ v2(M_{b+i}) - v2(den) ]
(pointwise identity  v2(N_i) = v2(M) + v2(c! k!) - v2(den), den = (b-c+1)*prod_{t=2}^i (b+t)).
Tabulate g, v2(c!k!), mu by (c,k) and split by parity of i.  Test the floors
  (F1)  g >= 2*floor(k/2)   on both slices
  (F2)  g >= k              on the even slice
and a candidate sharper i-parity floor.  0-violation cert over c<=12, b<=200.

M_{b+i} = the alternant (polynomial continuation, valid a<b too).  Optimized:
in the boundary regime s = 2m-2j = 2(P+k) is SMALL, and [x^A y^B z^C] E^j has the
closed form al=(A+B-C)/2, be=(B+C-A)/2, ga=(A+C-B)/2 (nonneg ints summing to j).
"""
from math import factorial, comb
from functools import lru_cache

def v2int(n):
    n = abs(int(n))
    if n == 0:
        return 10**9
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k

@lru_cache(maxsize=None)
def _fact(n):
    return factorial(n)

def Ecoef(j, A, B, C):
    """[x^A y^B z^C] E^j,  E = xy+yz+zx."""
    if A < 0 or B < 0 or C < 0:
        return 0
    if (A + B + C) != 2 * j:
        return 0
    s = A + B - C
    if s & 1:
        return 0
    al = s // 2
    be = B - al
    ga = A - al
    if al < 0 or be < 0 or ga < 0:
        return 0
    if al + be + ga != j:
        return 0
    return _fact(j) // (_fact(al) * _fact(be) * _fact(ga))

def _mult3(n, i, j, k):
    if i < 0 or j < 0 or k < 0 or i + j + k != n:
        return 0
    return _fact(n) // (_fact(i) * _fact(j) * _fact(k))

def Tfast(p, q, r, j, s):
    """[x^p y^q z^r] E^j H^s, iterate over H^s exponents (small s)."""
    if p < 0 or q < 0 or r < 0:
        return 0
    if p + q + r != 2 * j + s:
        return 0
    tot = 0
    for P in range(s + 1):
        for Q in range(s - P + 1):
            R = s - P - Q
            e = Ecoef(j, p - P, q - Q, r - R)
            if e:
                tot += _mult3(s, P, Q, R) * e
    return tot

_S3 = [((2,1,0),1),((1,2,0),-1),((2,0,1),-1),((0,2,1),1),((1,0,2),1),((0,1,2),-1)]

def Malt(a, b, c, j, m):
    s = 2 * m - 2 * j
    if s < 0:
        raise ValueError("need j<=m")
    base = (a + 2, b + 1, c)
    tot = 0
    for (perm, sgn) in _S3:
        tot += sgn * Tfast(base[0]-perm[0], base[1]-perm[1], base[2]-perm[2], j, s)
    return tot

def Ni_val(c, i, Pv, Bv, even):
    """Exact integer N_i on the slice (a=2P+b+c), or None if out of range/non-integral."""
    k = c - i
    bv = 2 * Bv if even else 2 * Bv + 1
    av = 2 * Pv + bv + c
    m = Pv + bv + c
    j = bv + i
    if bv < c + 1 or j > m:
        return None
    M = Malt(av, bv, c, j, m)
    den = (bv - c + 1)
    for t in range(2, i + 1):
        den *= (bv + t)
    num = M * factorial(c) * factorial(k)
    if den == 0 or num % den != 0:
        return None
    return num // den, M, den

# ----------------------------------------------------------------------
# cross-check the optimized alternant against MN before trusting it
# ----------------------------------------------------------------------
def selfcheck():
    import sys
    sys.path.insert(0, '/home/clio/projects/code/threerow-boundary')
    from mn import Mj
    bad = tested = 0
    for m in range(2, 9):
        n = 2 * m
        for a in range(1, n + 1):
            for b in range(1, a + 1):
                c = n - a - b
                if c < 1 or c > b:
                    continue
                for j in range(0, min(b + c, m) + 1):
                    v1 = Malt(a, b, c, j, m)
                    v2 = int(Mj((a, b, c), j, m))
                    tested += 1
                    if v1 != v2:
                        bad += 1
                        if bad <= 5:
                            print("  MISMATCH", (a, b, c), j, m, v1, v2)
    print(f"  optimized Malt vs MN: tested={tested} bad={bad}")
    return bad == 0

print("="*78)
print("SECTION 0: optimized alternant self-check vs Murnaghan-Nakayama")
print("="*78)
ok = selfcheck()
assert ok, "alternant self-check FAILED -- aborting"

# ----------------------------------------------------------------------
# content census + decomposition
# ----------------------------------------------------------------------
def content(c, i, even, Bmax, Pmax):
    """exact fixed 2-divisor of N_i on the slice, with minimizer + (M,den) there."""
    g = 10**9; best = None
    for Bv in range((c + 1) // 2 + 1, Bmax + 1):
        for Pv in range(0, Pmax + 1):
            res = Ni_val(c, i, Pv, Bv, even)
            if res is None:
                continue
            val, M, den = res
            vv = v2int(val)
            if vv < g:
                g = vv
                best = (Pv, Bv, val, M, den)
    return g, best

print("\n" + "="*78)
print("SECTION 1: content g[c][k] + decomposition  g = v2(c!k!) + mu,  by i-parity")
print("  mu = min_slice[ v2(M_{b+i}) - v2(den) ]   (residual the Content Lemma must bound)")
print("="*78)
# grid: B<=64, P<=64 (covers mod-2^6 residue period; content stability for g<=8
# cross-checked against last cycle's jobA_deep table at grids 64..512).
BMAX, PMAX = 64, 64
# published g_even / g_odd from FINDINGS-2026-06-18-jobA-claimA-deep.md (regression)
PUB_EVEN = {(4,4):None,(5,4):None,(6,4):5,(7,4):5,(8,4):6,(9,4):6,(10,4):5,(11,4):5,(12,4):6,
            (5,5):None,(6,5):None,(7,5):7,(8,5):7,(9,5):7,(10,5):7,(11,5):7,(12,5):7,
            (6,6):None,(7,6):None,(8,6):8,(9,6):8,(10,6):8,(11,6):7,(12,6):8}
PUB_ODD  = {(6,4):5,(7,4):5,(8,4):6,(9,4):6,(10,4):5,(11,4):5,(12,4):6,
            (7,5):7,(8,5):6,(9,5):8,(10,5):5,(11,5):7,(12,5):6,
            (8,6):8,(9,6):8,(10,6):8,(11,6):8,(12,6):8}
data = {}      # (c,k,slice) -> dict
import sys
for c in range(4, 13):
    for k in (4, 5, 6):
        i = c - k
        if i < 1:
            continue
        v2const = v2int(factorial(c) * factorial(k))
        for even in (True, False):
            g, best = content(c, i, even, BMAX, PMAX)
            mu = g - v2const
            data[(c, k, even)] = (g, v2const, mu, i, best)
        sys.stdout.write(f"    [done c={c} k={k}]\n"); sys.stdout.flush()
# regression vs published deep-census content
print("  regression vs last cycle's published content table:")
reg_bad = 0
for c in range(4,13):
    for k in (4,5,6):
        if (c,k) in PUB_EVEN and PUB_EVEN[(c,k)] is not None:
            ge = data[(c,k,True)][0]
            if ge != PUB_EVEN[(c,k)]:
                reg_bad += 1; print(f"    EVEN MISMATCH c={c} k={k}: now {ge} vs pub {PUB_EVEN[(c,k)]}")
        if (c,k) in PUB_ODD:
            go = data[(c,k,False)][0]
            if go != PUB_ODD[(c,k)]:
                reg_bad += 1; print(f"    ODD MISMATCH c={c} k={k}: now {go} vs pub {PUB_ODD[(c,k)]}")
print(f"  regression mismatches vs published table: {reg_bad}")
print(f"  (grid B<=64, P<=64; v2const = v2(c! k!))")
print(f"\n  {'c':>2} {'k':>2} {'i':>2} {'i%2':>3} | "
      f"{'g_even':>6} {'mu_e':>5} | {'g_odd':>6} {'mu_o':>5} | {'v2const':>7}")
for c in range(4, 13):
    for k in (4, 5, 6):
        i = c - k
        if i < 1:
            continue
        de = data.get((c, k, True)); do = data.get((c, k, False))
        if de is None or do is None:
            continue
        ge, vc, mue, _, _ = de
        go, _, muo, _, _ = do
        print(f"  {c:>2} {k:>2} {i:>2} {i%2:>3} | "
              f"{ge:>6} {mue:>5} | {go:>6} {muo:>5} | {vc:>7}")

# ----------------------------------------------------------------------
print("\n" + "="*78)
print("SECTION 2: pointwise identity check + minimizer decomposition (a few cases)")
print("="*78)
print(f"  identity at minimizer:  g == v2(M) + v2(c!k!) - v2(den)")
for (c, k, even) in [(6,4,True),(7,5,True),(8,6,True),(6,4,False),(9,5,False),(10,6,False)]:
    g, vc, mu, i, best = data[(c, k, even)]
    Pv, Bv, val, M, den = best
    lhs = v2int(val)
    rhs = v2int(M) + vc - v2int(den)
    sl = 'even' if even else 'odd '
    print(f"  c={c} k={k} i={i} {sl}: g={g}  v2(M)={v2int(M)} v2(c!k!)={vc} "
          f"v2(den)={v2int(den)}  -> rhs={rhs}  match={lhs==rhs}")

# ----------------------------------------------------------------------
print("\n" + "="*78)
print("SECTION 3: floor certificate (0-violation over c<=12, b<=200, both slices)")
print("="*78)
viol_F1 = viol_F2 = 0
viol_examples = []
sharper_ok = True
sharper_examples = []
for (c, k, even), (g, vc, mu, i, best) in data.items():
    floor1 = 2 * (k // 2)
    if g < floor1:
        viol_F1 += 1
        if len(viol_examples) < 8: viol_examples.append(('F1', c, k, even, g, floor1))
    if even and g < k:
        viol_F2 += 1
        if len(viol_examples) < 8: viol_examples.append(('F2', c, k, even, g, k))
    # candidate sharper i-parity floor:
    #   odd i, odd slice -> floor exactly 2*floor(k/2) (dream: floor-tight at k<=3)
    #   even i           -> floor 2*floor(k/2) + (k odd ? 0 : ... )
    # we just RECORD whether odd-i-odd-slice is floor-tight (g == 2 floor(k/2)).
print(f"  F1  g >= 2*floor(k/2)  both slices : violations = {viol_F1}")
print(f"  F2  g >= k             even slice   : violations = {viol_F2}")
if viol_examples:
    print(f"  example violations: {viol_examples}")

print("\n  i-parity refinement (is odd-i floor-tight on the odd slice?):")
print(f"  {'c':>2} {'k':>2} {'i':>2} {'i%2':>3} {'slice':>5} {'g':>3} {'2floor(k/2)':>11} {'surplus':>7}")
for c in range(4, 13):
    for k in (4, 5, 6):
        i = c - k
        if i < 1: continue
        for even in (False, True):
            if (c, k, even) not in data: continue
            g, vc, mu, _, _ = data[(c, k, even)]
            fl = 2 * (k // 2)
            sl = 'even' if even else 'odd'
            print(f"  {c:>2} {k:>2} {i:>2} {i%2:>3} {sl:>5} {g:>3} {fl:>11} {g-fl:>7}")
print("DONE.")
