"""
jobA_deep.py -- 2026-06-18 CODE session, Job A (deep-k Claim A census).

Pins the EXACT 2-adic slice-content law g[c][k] for the master deficit
polynomial N_i^{(c)} at deep depth k = c - i, for c=4..12, k=0..6, BOTH parity
slices of b.  Feeds the PROVE Claim-A-at-k>=4 target.

Master convention (matches claimA_verify.py / g0_fixeddiv.py):
   N_i^{(c)}(a,b) = M_{b+i}(a,b,c) * c! * k! / [ (b-c+1) * prod_{t=2}^{i} (b+t) ]
   with M_{b+i} the alternant (fast_alt.Malt), k = c - i.  i=c is top (N=1, k=0).

Method (robust, fit-then-evaluate):
  1. Fit N_i as an EXACT rational-coefficient polynomial in the slice half-vars
     (P, B) where a = 2P+b+c, b = 2B (even) or 2B+1 (odd).  Degrees are small.
     Validate the fit against direct Malt on held-out points (incl. large b).
  2. CONTENT = true fixed 2-divisor = min over a full residue grid (P,B in
     [0,256)) of v2(N_i(P,B)).  This is the exact quantity, NOT the coeff-gcd
     lower bound the old fit used.
  3. Report g[c][k] per slice; verify even-slice == k for k<=6; pin odd-slice law.
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
    """Exact integer value of master N_i on the slice, or None if out of range."""
    k = c - i
    bv = 2 * Bv if even else 2 * Bv + 1
    av = 2 * Pv + bv + c
    m = Pv + bv + c          # (a+b+c)/2 = P + b + c
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
    """Exact bivariate poly fit of N_i(P,B) on the slice via square tensor-grid
    interpolation.  Returns expanded sympy expr, or None on failure.
    Validity: for P>=0, j=b+i <= m=P+b+c iff i<=P+c, always true; just need b>=c+1."""
    monos = [(p, q) for p in range(pdeg + 1) for q in range(bdeg + 1)]
    # tensor grid: P in 0..pdeg, B in B0..B0+bdeg with b>=c+1
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
    y = sp.Matrix(rhs)
    try:
        sol = A.LUsolve(y)
    except Exception:
        return None
    expr = sp.expand(sum(sol[t] * P**p * B**q for t, (p, q) in enumerate(monos)))
    return expr


def validate_fit(expr, c, i, even, ntests=40):
    """Check fitted poly == direct Malt on points NOT in the fit set (large b)."""
    Bbig = 60
    cnt = 0
    for Bv in range(Bbig, Bbig + ntests):
        for Pv in range(0, 6):
            v = Ni_val(c, i, Pv, Bv, even)
            if v is None:
                continue
            pv = int(expr.subs({P: Pv, B: Bv}))
            if pv != v:
                return False, (Pv, Bv, pv, v)
            cnt += 1
            if cnt >= ntests:
                return True, cnt
    return True, cnt


def to_int_poly(expr):
    """Return (Nint sympy Poly with integer coeffs, v2 of clearing denom D)."""
    poly = sp.Poly(expr, P, B)
    coeffs = poly.coeffs()
    dens = [sp.Rational(cf).q for cf in coeffs]
    D = 1
    for d in dens:
        D = sp.ilcm(D, d)
    Nint = sp.expand(expr * D)
    return sp.Poly(Nint, P, B), v2int(D)


def content_via_poly(expr, grid=256):
    """True fixed 2-divisor of N_i on the slice = min over (P,B) in [0,grid)^2
    of v2(value).  Uses integer-coeff scaled poly for speed."""
    Nint, vD = to_int_poly(expr)
    # fast evaluator: build callable over python ints
    f = sp.lambdify((P, B), Nint.as_expr(), modules='math')
    g = 10**9
    wit = None
    for Pv in range(grid):
        for Bv in range(grid):
            val = int(round(f(Pv, Bv)))
            # recompute exactly when near small magnitude to avoid float error
            vv = v2int(val)
            if vv < g:
                # exact recheck
                exact = int(Nint.as_expr().subs({P: Pv, B: Bv}))
                vv = v2int(exact)
                if vv < g:
                    g = vv
                    wit = (Pv, Bv)
    return g - vD, wit


def content_exact_smallgrid(expr, grid=128):
    """Exact (no float) content over a smaller grid -- used as cross-check."""
    Nint, vD = to_int_poly(expr)
    fe = sp.lambdify((P, B), Nint.as_expr(), modules='sympy')
    g = 10**9
    wit = None
    for Pv in range(grid):
        for Bv in range(grid):
            val = int(fe(Pv, Bv))
            vv = v2int(val)
            if vv < g:
                g = vv
                wit = (Pv, Bv)
    return g - vD, wit


def degrees_for(k):
    """generous degrees; a-deg(N) ~ floor(k/2)+1 in P, ~k+2 in B. Pad."""
    return (k // 2 + 3, k + 4)


if __name__ == "__main__":
    print("=" * 90)
    print("Job A deep census: EXACT fixed-divisor content g[c][k] of master N_i^{(c)}")
    print("  k = c - i (depth) ; even slice b=2B, odd slice b=2B+1 ; content = min v2")
    print("=" * 90)
    table = {}     # (c,k) -> (g_even, g_odd)
    fails = []
    for c in range(4, 13):
        for i in range(2, c + 1):
            k = c - i
            if k > 6:
                continue
            pdeg, bdeg = degrees_for(k)
            row = "c=%-2d i=%-2d k=%-2d |" % (c, i, k)
            res = {}
            for even in (True, False):
                expr = fit_slice_poly(c, i, even, pdeg, bdeg)
                if expr is None:
                    row += " %s-slice NO FIT |" % ("even" if even else "odd")
                    res[even] = None
                    continue
                ok, info = validate_fit(expr, c, i, even)
                if not ok:
                    row += " VALIDATE FAIL %s |" % str(info)
                    res[even] = None
                    fails.append((c, i, even, "validate", info))
                    continue
                g, wit = content_exact_smallgrid(expr, grid=64)
                res[even] = g
            ge, go = res.get(True), res.get(False)
            table[(c, k)] = (ge, go)
            sge = "%d" % ge if ge is not None else "?"
            sgo = "%d" % go if go is not None else "?"
            tag = ""
            if ge is not None and ge != k:
                tag += " [EVEN!=k]"
                fails.append((c, i, "even", "even!=k", (ge, k)))
            print("c=%-2d i=%-2d k=%-2d | even=%s odd=%s | (want even=k=%d)%s"
                  % (c, i, k, sge, sgo, k, tag))

    # ---- tabulate odd-slice law ----
    print("\n" + "=" * 90)
    print("ODD-SLICE content go[c][k]  (rows c, cols k); even-slice = k throughout")
    print("=" * 90)
    header = "      " + "".join("k=%-4d" % k for k in range(0, 7))
    print(header)
    for c in range(4, 13):
        rowstr = "c=%-3d " % c
        for k in range(0, 7):
            if (c, k) in table and table[(c, k)][1] is not None:
                rowstr += "%-6d" % table[(c, k)][1]
            else:
                rowstr += "%-6s" % "."
        print(rowstr)

    # ---- analyze the odd-slice law as function of (k, c mod 2, c mod 4) ----
    print("\n" + "=" * 90)
    print("ODD-SLICE law analysis: go vs 2*floor(k/2), split by c parity / c mod 4")
    print("=" * 90)
    for k in range(0, 7):
        flr = 2 * (k // 2)
        vals_evenc = sorted(set(table[(c, k)][1] for c in range(4, 13)
                                if (c, k) in table and c % 2 == 0
                                and table[(c, k)][1] is not None))
        vals_oddc = {}
        for c in range(4, 13):
            if (c, k) in table and c % 2 == 1 and table[(c, k)][1] is not None:
                vals_oddc.setdefault(c % 4, set()).add(table[(c, k)][1])
        print("k=%d (2floor=%d): even-c odd-slice=%s | odd-c by (c mod4): %s"
              % (k, flr, vals_evenc,
                 {r: sorted(s) for r, s in sorted(vals_oddc.items())}))

    print("\nFAILS:", fails if fails else "NONE")
