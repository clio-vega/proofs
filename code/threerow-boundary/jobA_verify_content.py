"""
jobA_verify_content.py -- sanity gate for the fixed-divisor content table.

Fast & rigorous, two checks:
 (1) GRID STABILITY: content = min v2(N_i) over (P,B) in [0,G)^2 (within parity).
     grid-min is an UPPER bound on the true content (min over a subset >= global
     min), so if it does NOT decrease from G=64 to 512 the value is the true
     content.  Uses the VALIDATED slice polynomial (fast int eval).
 (2) DIRECT CROSS-CHECK: on a small window (b<=120, P<=30) compute content
     straight from Malt and confirm it equals the poly-based content over the
     same window -> the poly reproduces N_i's 2-adics, not just its values.
Plus Malt-vs-MN on real partitions (transitive trust to Murnaghan-Nakayama).
"""
import sympy as sp
from math import factorial
from fast_alt import Malt
from mn import Mj

P, B = sp.symbols('P B')


def v2int(n):
    n = abs(int(n))
    if n == 0:
        return 10**9
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def Ni_val(c, i, Pv, Bv, even):
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
    return num // den


def fit_slice_poly(c, i, even, pdeg, bdeg):
    monos = [(p, q) for p in range(pdeg + 1) for q in range(bdeg + 1)]
    B0 = ((c + 1) // 2) + 2
    pts, rhs = [], []
    for Pv in range(0, pdeg + 1):
        for Bv in range(B0, B0 + bdeg + 1):
            v = Ni_val(c, i, Pv, Bv, even)
            if v is None:
                return None
            pts.append((Pv, Bv))
            rhs.append(sp.Integer(v))
    A = sp.Matrix([[sp.Integer(Pv)**p * sp.Integer(Bv)**q for (p, q) in monos]
                   for (Pv, Bv) in pts])
    sol = A.LUsolve(sp.Matrix(rhs))
    return sp.expand(sum(sol[t] * P**p * B**q for t, (p, q) in enumerate(monos)))


def int_coeffs(expr):
    poly = sp.Poly(sp.expand(expr), P, B)
    D = 1
    for cf in poly.coeffs():
        D = sp.ilcm(D, sp.Rational(cf).q)
    terms = [(int(m[0]), int(m[1]), int(sp.Rational(cf) * D))
             for m, cf in poly.terms()]
    return terms, v2int(D)


def poly_content(terms, vD, even, grid, c):
    qmax = max(t[1] for t in terms)
    g = 10**9
    B0 = ((c + 1) // 2) + 1
    for Bv in range(B0, grid):
        Bpow = [Bv**q for q in range(qmax + 1)]
        for Pv in range(0, grid):
            val = 0
            for (pp, qq, cc) in terms:
                val += cc * (Pv ** pp) * Bpow[qq]
            vv = v2int(val) - vD
            if vv < g:
                g = vv
    return g


def direct_content(c, i, even, Bmax, Pmax):
    g = 10**9
    B0 = ((c + 1) // 2) + 1
    for Bv in range(B0, Bmax + 1):
        for Pv in range(0, Pmax + 1):
            v = Ni_val(c, i, Pv, Bv, even)
            if v is None:
                continue
            vv = v2int(v)
            if vv < g:
                g = vv
    return g


def malt_vs_mn(c, i):
    bad = cnt = 0
    B0 = ((c + 1) // 2) + 1
    for Bv in range(B0, B0 + 25):
        for Pv in range(0, 6):
            bv = 2 * Bv
            av = 2 * Pv + bv + c
            m = Pv + bv + c
            j = bv + i
            if av < bv or j > m:
                continue
            if Malt(av, bv, c, j, m) != int(Mj((av, bv, c), j, m)):
                bad += 1
            cnt += 1
    return bad, cnt


CASES = [(6, 2, 4), (7, 2, 5), (7, 4, 3), (8, 2, 6), (8, 5, 3),
         (9, 4, 5), (9, 6, 3), (10, 6, 4), (11, 7, 4), (12, 9, 3),
         (5, 2, 3), (6, 3, 3), (10, 8, 2)]


def degrees_for(k):
    return (k // 2 + 3, k + 4)


if __name__ == "__main__":
    print("Malt-vs-MN (must be 0 bad):")
    for (c, i, _) in CASES[:5]:
        print("  c=%d i=%d : %d/%d bad" % ((c, i) + malt_vs_mn(c, i)))

    print("\nGrid stability + direct cross-check (poly content @ grids; direct small):")
    print("  c  i  k slc | g64 g128 g256 g512 | poly@win direct@win  match?")
    allok = True
    for (c, i, k) in CASES:
        pdeg, bdeg = degrees_for(k)
        for even in (True, False):
            expr = fit_slice_poly(c, i, even, pdeg, bdeg)
            terms, vD = int_coeffs(expr)
            gs = [poly_content(terms, vD, even, G, c) for G in (64, 128, 256, 512)]
            stable = (len(set(gs)) == 1)
            # window cross-check
            win = 60
            pw = poly_content(terms, vD, even, win, c)
            dw = direct_content(c, i, even, Bmax=win - 1, Pmax=win - 1)
            match = (pw == dw)
            if not (stable and match):
                allok = False
            print("  %-2d %-2d %-2d %-3s | %-3d %-3d %-3d %-3d | %-8d %-10d  %s%s"
                  % (c, i, k, "ev" if even else "od", gs[0], gs[1], gs[2], gs[3],
                     pw, dw, "OK" if match else "MISMATCH",
                     "" if stable else " <-DRIFT"))
    print("\nALL stable & direct-matched:", allok)
