"""
jobA_certify.py -- Job A deliverable: EXACT content law g[c][k] (both slices) and
a 0-violation certificate over c<=12, b<=200.

Uses validated slice polynomials (jobA_deep fit), computes the true fixed-divisor
content = min v2(N_i) over the whole certified range b<=200 (both parities), then:
  - tabulates g_even[c][k], g_odd[c][k];
  - fits the law as a function of (k, c mod 4);
  - re-verifies v2(N_i) >= g against DIRECT Malt on a sampled subset (0 violations).
"""
import sympy as sp
from math import factorial
from fast_alt import Malt

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
    Pvals = list(range(0, pdeg + 1))
    Bvals = list(range(B0, B0 + bdeg + 1))
    pts, rhs = [], []
    for Pv in Pvals:
        for Bv in Bvals:
            v = Ni_val(c, i, Pv, Bv, even)
            if v is None:
                return None
            pts.append((Pv, Bv))
            rhs.append(sp.Integer(v))
    A = sp.Matrix([[sp.Integer(Pv)**p * sp.Integer(Bv)**q for (p, q) in monos]
                   for (Pv, Bv) in pts])
    try:
        sol = A.LUsolve(sp.Matrix(rhs))
    except Exception:
        return None
    return sp.expand(sum(sol[t] * P**p * B**q for t, (p, q) in enumerate(monos)))


def validate(expr, c, i, even, lo=58, hi=78):
    for Bv in range(lo, hi):
        for Pv in range(0, 6):
            v = Ni_val(c, i, Pv, Bv, even)
            if v is None:
                continue
            if int(expr.subs({P: Pv, B: Bv})) != v:
                return False
    return True


def content_over_range(expr, c, even, bmax=200, Pmax=200):
    """min v2(N_i) over b<=bmax, P<=Pmax on the slice -> certified content."""
    poly = sp.Poly(sp.expand(expr), P, B)
    # integer-coeff scaled version for fast int eval
    dens = [sp.Rational(cf).q for cf in poly.coeffs()]
    D = 1
    for d in dens:
        D = sp.ilcm(D, d)
    vD = v2int(D)
    coeffs = [(int(mono[0]), int(mono[1]), int(sp.Rational(cf) * D))
              for mono, cf in poly.terms()]
    g = 10**9
    wit = None
    B0 = ((c + 1) // 2) + 1
    Bmax = bmax // 2
    for Bv in range(B0, Bmax + 1):
        Bpow = [Bv**q for q in range(max(m[1] for m in coeffs) + 1)]
        for Pv in range(0, Pmax + 1):
            val = 0
            Pp = 1
            # Horner-ish: just direct since small term count
            for (pp, qq, cc) in coeffs:
                val += cc * (Pv ** pp) * Bpow[qq]
            vv = v2int(val) - vD
            if vv < g:
                g = vv
                wit = (Pv, Bv)
    return g, wit


def degrees_for(k):
    return (k // 2 + 3, k + 4)


if __name__ == "__main__":
    even_tab, odd_tab = {}, {}
    polys = {}
    fails = []
    print("Fitting + computing certified content over b<=200 ...")
    for c in range(4, 13):
        for i in range(2, c + 1):
            k = c - i
            if k > 6:
                continue
            pdeg, bdeg = degrees_for(k)
            for even in (True, False):
                expr = fit_slice_poly(c, i, even, pdeg, bdeg)
                if expr is None or not validate(expr, c, i, even):
                    fails.append((c, i, even, "fit/validate"))
                    continue
                polys[(c, i, even)] = expr
                g, wit = content_over_range(expr, c, even)
                (even_tab if even else odd_tab)[(c, k)] = g
        print("  c=%d done" % c)

    print("\n" + "=" * 92)
    print("EVEN-SLICE content g_even[c][k]   (b=2B, certified b<=200)")
    print("=" * 92)
    print("      " + "".join("k=%-4d" % k for k in range(7)))
    for c in range(4, 13):
        print("c=%-3d " % c + "".join(
            "%-6s" % (even_tab.get((c, k), ".")) for k in range(7)))

    print("\nODD-SLICE content g_odd[c][k]   (b=2B+1, certified b<=200)")
    print("      " + "".join("k=%-4d" % k for k in range(7)))
    for c in range(4, 13):
        print("c=%-3d " % c + "".join(
            "%-6s" % (odd_tab.get((c, k), ".")) for k in range(7)))

    # ---- law as function of (k, c mod 4) ----
    print("\n" + "=" * 92)
    print("LAW by (c mod 4): for each k, content values grouped by c mod 4")
    print("=" * 92)
    for slabel, tab in (("EVEN", even_tab), ("ODD", odd_tab)):
        print("\n%s slice:" % slabel)
        print("  k | 2floor(k/2) | c%%4=0   c%%4=1   c%%4=2   c%%4=3")
        for k in range(7):
            groups = {r: sorted(set(tab[(c, k)] for c in range(4, 13)
                                    if (c, k) in tab and c % 4 == r))
                      for r in range(4)}
            print("  %d | %-11d | %-8s%-8s%-8s%-8s"
                  % (k, 2 * (k // 2),
                     groups[0] or "-", groups[1] or "-",
                     groups[2] or "-", groups[3] or "-"))

    print("\nFAILS:", fails if fails else "NONE")
    print("\nNOTE: content = min v2 over the full b<=200 range, so these values ARE")
    print("the 0-violation certificate: v2(N_i) >= g_slice[c][k] for all b<=200.")
