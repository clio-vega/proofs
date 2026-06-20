#!/usr/bin/env python3
"""Job B (2026-06-20 CODE): the beta'(c) constant-floor sequence + beta vs beta'.

PROVE asks: does the heavy quotient H_c carry a CONSTANT 2-adic floor 2^{beta'(c)}
for all c, and what is the sequence?  We have beta'(4)=4, beta'(5)=3 (non-monotonic).

Setup (verified empirically for c=6,7,8 before use):
  closed form   M_j = C(2(m-j),b-j)(a-b+1) Q_c / [ c! (a+c+1-j) prod_{i=1..c}(b+i-j) ]
  heavy/tip     Q_c = (a-(c-2))(b-(c-1)) H_c  -  (-1)^c (2c)! C(j,2c)
  =>            H_c(a,b,j) = [ Q_c + (-1)^c (2c)! C(j,2c) ] / [(a-(c-2))(b-(c-1))]

Q_c is read off MN; H_c is reconstructed SYMBOLICALLY by interpolation (deg_a=deg_b=c-1,
deg_j=2c-2) and the reconstruction is cross-checked against MN.  beta'(c) = min 2-adic
valuation of H_c over the valid-parity sublattice (a==b mod 2 for c even, a!=b mod 2 for
c odd; j free), found by a residue-grid scan mod 2^k with growing k until stable.

Also: beta(c) = (c-1) + v2((c-1)!)  for c=2..8, to settle the c=5 6-vs-7 question.
"""
import sys, time
sys.path.insert(0, '/home/clio/projects/code/threerow-c3')
from mn import Mj
from math import comb, factorial
from fractions import Fraction
import sympy as sp

a_s, b_s, j_s = sp.symbols('a b j')

def v2(n):
    n = abs(int(n))
    if n == 0:
        return float('inf')
    return (n & -n).bit_length() - 1

def s2(n):
    return bin(int(n)).count('1')

# ---------------------------------------------------------------------------
# Black-box integer H_c(A,B,J) from MN (any valid interior shape, any J).
# ---------------------------------------------------------------------------
def Qc_int(c, A, B, J):
    m = (A + B + c) // 2
    M = Mj((A, B, c), J, m)
    N = 2 * (m - J)
    df = factorial(c) * (A + c + 1 - J)
    for i in range(1, c + 1):
        df *= (B + i - J)
    num = M * df
    den = comb(N, B - J) * (A - B + 1)
    q = Fraction(num, den)
    assert q.denominator == 1, f"Q_{c} non-integer ({A},{B}) J={J}"
    return int(q)

def Hc_int(c, A, B, J):
    # Convention (verified): Q_c = heavy*H + (-1)^c (2c)! C(j,2c),  heavy=(a-(c-2))(b-(c-1)).
    #   c=4: Q4 = heavy*H + 8! C(j,8)   (sign +)
    #   c=5: Q5 = heavy*H - 10!C(j,10)  (sign -)
    # => H = (Q_c - (-1)^c (2c)! C(j,2c)) / heavy.
    q = Qc_int(c, A, B, J)
    heavy = (A - (c - 2)) * (B - (c - 1))
    Hnum = q - (-1)**c * factorial(2 * c) * comb(J, 2 * c)
    assert Hnum % heavy == 0, f"heavy does not divide ({A},{B}) J={J}"
    return Hnum // heavy

# ---------------------------------------------------------------------------
# Symbolic reconstruction of H_c via interpolation (generalises H5.py).
# deg_a = deg_b = c-1,  deg_j = 2c-2.
# ---------------------------------------------------------------------------
def interp_ab(c, jfix):
    """Interpolate H_c(a,b,jfix) as a 2D polynomial of degree (c-1,c-1)."""
    npts = c                       # c points -> degree c-1
    # b-values all > 2c (>= jfix since jfix <= 2c-2) and > c, distinct -> valid shapes
    Bvals = [2 * c + 2 + 2 * t for t in range(npts)]
    polysA = []
    for B0 in Bvals:
        ptsA = []
        A0 = B0 + 3
        seen = set()
        while len(ptsA) < npts:
            A = A0
            if (A + B0 + c) % 2 != 0:                 # need m integer
                A += 1
            if A > B0 and A not in seen:
                ptsA.append((A, Hc_int(c, A, B0, jfix)))
                seen.add(A)
            A0 = A + 1
        pA = sp.interpolate(ptsA, a_s)
        polysA.append((B0, sp.expand(pA)))
    # interpolate each a-coefficient across b
    maxdeg = max(sp.Poly(p, a_s).degree() for _, p in polysA)
    result = 0
    for deg in range(maxdeg + 1):
        ptsB = []
        for B0, p in polysA:
            P = sp.Poly(p, a_s)
            coeff = P.coeff_monomial(a_s**deg) if deg > 0 else P.coeff_monomial(1)
            ptsB.append((B0, coeff))
        pB = sp.interpolate(ptsB, b_s)
        result += sp.expand(pB) * a_s**deg
    return sp.expand(result)

def binom_poly(var, k):
    r = sp.Integer(1)
    for t in range(k):
        r *= (var - t)
    return r / factorial(k)

def build_Hc(c):
    """Return (H_c symbolic, fast integer evaluator, validation report)."""
    degj = 2 * c - 2
    # H_c(a,b,i) for i = 0..degj  (all i < 2c so tip vanishes there anyway)
    Hi = [interp_ab(c, i) for i in range(degj + 1)]
    hcoeff = []
    for k in range(degj + 1):
        hk = sum((-1)**(k - i) * comb(k, i) * Hi[i] for i in range(k + 1))
        hcoeff.append(sp.expand(hk))
    Hsym = sp.expand(sum(hcoeff[k] * binom_poly(j_s, k) for k in range(degj + 1)))
    # fast integer evaluator
    terms = [(tuple(int(e) for e in mono), int(co))
             for mono, co in sp.Poly(sp.expand(Hsym), a_s, b_s, j_s).terms()]
    def H_eval(av, bv, jv):
        return sum(co * av**ea * bv**eb * jv**ej for (ea, eb, ej), co in terms)
    # H(0) leading-product sanity:  expect prod_{t=2}^{c+1}(a+c+t-... ) -- just record factor
    H0fac = sp.factor(Hsym.subs(j_s, 0))
    return Hsym, H_eval, H0fac, terms

# ---------------------------------------------------------------------------
# Validate reconstruction vs MN, then floor scan.
# ---------------------------------------------------------------------------
def validate(c, H_eval, ntest=300):
    import random
    rng = random.Random(1234 + c)
    bad = 0; tested = 0
    for _ in range(ntest):
        B = rng.randrange(c + 1, c + 35)
        A = rng.randrange(B + 1, B + 45)
        if (A + B + c) % 2 != 0:
            A += 1
        if A <= B:
            continue
        J = rng.randrange(0, B + 1)
        if B - J < 0:
            continue
        mn = Hc_int(c, A, B, J)
        sy = H_eval(A, B, J)
        tested += 1
        if mn != sy:
            bad += 1
            if bad <= 4:
                print(f"   RECON MISMATCH c={c} ({A},{B}) J={J}: MN={mn} sym={sy}")
    return tested, bad

def content_floor(c, Hsym):
    """RIGOROUS beta'(c): the 2-adic content of H_c over the valid-parity sublattice.

    An integer-valued polynomial P(x) satisfies 2^t | P(n) for ALL integers n iff,
    expanded in the integer-valued basis prod C(x_i, m_i), every coefficient is
    divisible by 2^t.  The content (min t failing) = min 2-adic valuation of those
    basis coefficients = (Delta-table at 0).  We split into parity sheets:
       c even: (a,b) = (2u+r, 2v+r),  r in {0,1}   [a == b mod 2]
       c odd : (a,b) = (2u+r, 2v+1-r), r in {0,1}  [a != b mod 2]
    j is unconstrained (its own variable).  Returns (floor, sheet) -- EXACT, no scan."""
    u, vv = sp.symbols('u v')
    da = db = c - 1
    dj = 2 * c - 2
    sheets = [(0, 0), (1, 1)] if c % 2 == 0 else [(0, 1), (1, 0)]
    best = None
    best_sheet = None
    for (ra, rb) in sheets:
        P = sp.expand(Hsym.subs({a_s: 2 * u + ra, b_s: 2 * vv + rb}))
        # grid of integer values of P on [0..da]x[0..db]x[0..dj]
        G = [[[int(P.subs({u: i, vv: k, j_s: l})) for l in range(dj + 1)]
              for k in range(db + 1)] for i in range(da + 1)]
        # full forward-difference table: coeff[i][k][l] = (Delta_u^i Delta_v^k Delta_j^l P)(0)
        # difference along axis u
        Au = [[[ _fd([G[ii][k][l] for ii in range(da + 1)], i)
                 for l in range(dj + 1)] for k in range(db + 1)] for i in range(da + 1)]
        Av = [[[ _fd([Au[i][kk][l] for kk in range(db + 1)], k)
                 for l in range(dj + 1)] for k in range(db + 1)] for i in range(da + 1)]
        Aj = [[[ _fd([Av[i][k][ll] for ll in range(dj + 1)], l)
                 for l in range(dj + 1)] for k in range(db + 1)] for i in range(da + 1)]
        mn = min(v2(Aj[i][k][l]) for i in range(da + 1)
                 for k in range(db + 1) for l in range(dj + 1))
        if best is None or mn < best:
            best = mn
            best_sheet = (ra, rb)
    return best, best_sheet

def _fd(seq, order):
    """order-th forward difference of seq evaluated at index 0."""
    s = list(seq)
    for _ in range(order):
        s = [s[i + 1] - s[i] for i in range(len(s) - 1)]
    return s[0]

def small_witness(c, H_eval, floor, R=20):
    """Find a small integer (a,b,j) of valid parity attaining v2(H_c) = floor."""
    same = (c % 2 == 0)
    for av in range(-R, R + 1):
        for bv in range(-R, R + 1):
            if ((av - bv) % 2 == 0) != same:
                continue
            for jv in range(-2, R + 1):
                val = H_eval(av, bv, jv)
                if val != 0 and v2(val) == floor:
                    return (av, bv, jv)
    return None

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 72)
    print("Job B: beta'(c) constant-floor sequence")
    print("=" * 72)
    results = {}
    for c in (4, 5, 6, 7, 8):
        t0 = time.time()
        print(f"\n--- c = {c} ---")
        Hsym, H_eval, H0fac, terms = build_Hc(c)
        tested, bad = validate(c, H_eval)
        print(f"  symbolic reconstruction vs MN: tested={tested}  mismatches={bad}")
        if bad:
            print("  *** reconstruction FAILED, skipping floor scan ***")
            continue
        floor, sheet = content_floor(c, Hsym)
        w = small_witness(c, H_eval, floor)
        results[c] = (floor, w)
        print(f"  H_c(0) factor = {H0fac}")
        print(f"  beta'({c}) = {floor}  (RIGOROUS: binomial-basis content, parity sheet {sheet})")
        wval = H_eval(*w) if w else None
        print(f"  sharp witness (a,b,j)={w}  v2(H_c)={v2(wval) if w else 'n/a'}")
        print(f"  [c={c}: {time.time()-t0:.1f}s]")

    print("\n" + "=" * 72)
    print("beta'(c) SEQUENCE")
    print("=" * 72)
    seq = {4: 4, 5: 3}
    for c in sorted(results):
        seq[c] = results[c][0]
    print("  c      :", "  ".join(f"{c:2d}" for c in sorted(seq)))
    print("  beta'(c):", "  ".join(f"{seq[c]:2d}" for c in sorted(seq)))
    print("  witnesses:")
    for c in sorted(results):
        print(f"    beta'({c})={results[c][0]} at (a,b,j)={results[c][1]}")
    # known: beta'(4)=4 (from c4 census 16|H), beta'(5)=3 (8|H5 sharp)
    print(f"\n  sanity: my scan gives beta'(4)={results.get(4,('?',))[0]} (expect 4), "
          f"beta'(5)={results.get(5,('?',))[0]} (expect 3)")

    print("\n" + "=" * 72)
    print("beta(c) = (c-1) + v2((c-1)!)   [the RIGID/Number-Lemma exponent]")
    print("=" * 72)
    print("   c | (c-1)+v2((c-1)!) | 2(c-1)-s2(c-1) | match?")
    for c in range(2, 9):
        b1 = (c - 1) + v2(factorial(c - 1))
        b2 = 2 * (c - 1) - s2(c - 1)
        print(f"   {c} |        {b1:2d}        |      {b2:2d}        | {b1 == b2}")
    print()
    print("  RECONCILIATION for c=5:")
    b5 = (5 - 1) + v2(factorial(4))
    print(f"    beta(5) = 4 + v2(4!) = 4 + v2(24) = 4 + 3 = {b5}")
    print(f"    2(5-1) - s2(4) = 8 - 1 = {2*4 - s2(4)}")
    print(f"    => beta(5) = 7  (NOT 6).  The c=5 writeup's '6' is an error;")
    print(f"       Rick's prediction 7 and the formula agree.  beta' (the FLOOR) is 3;")
    print(f"       beta (the Number-Lemma RIGID exponent) is a DIFFERENT quantity = 7.")
