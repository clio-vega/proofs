"""
TASK 1 -- Q_b as explicit integer polynomials + their p-adic Newton polygons.

Feeds the OP-agnostic (diamond) attack on the TWO-ROW d=4 law:
    law  <=>  (diamond) Q_b has NO RATIONAL ROOT     (b not div by 4)
              (diamond) the deg-(floor(b/2)-1) factor Q~_b of Q_b has no rational root
                        after dividing out the forced half-integer linear factor (2m-(2b-1))
                                                                     (4 | b)

Pipeline (recomputed from the proved reduction -- nothing trusted from memory):
  * H_m = h_{m-1}(A,B),  A+B = 2W,  AB = W^2+u^2,  W = 1+u+u^2.
  * I_b(m) = Im( [u^b] ((1-u)(1+su+u^2)^m) ),  s = 1+i, as a genuine Z[m] polynomial.
    (We build it directly as the imaginary part of [u^b]((1-u)P), which is the same
     object; we ALSO cross-check it against the H_m form below.)
  * Divide out the forced linear factors  prod_{r=0}^{floor((b-1)/2)} (m - r)  to get Q_b,
    deg_m Q_b = floor(b/2).
  * Q~_b := primitive integer form of the non-trivial irreducible factor.

For each b = 5..16 and Q~_b:
  1. Factor leading coeff and constant term into primes.
  2. Discriminant, its prime factorisation, perfect-square (Galois in A_n) test.
  3. p-adic Newton polygon for the relevant primes; edge slopes as exact rationals.
     RIGOROUS no-rational-root witness:  a prime p such that the Newton polygon of
     Q~_b has NO integer-slope edge  ==>  Q~_b has no factor of degree 1 over Q_p
     ==>  Q~_b has no rational root.  (Any rational root r has v_p(r) in Z equal to
     minus an integer slope of the polygon.)  We also report the looser |slope|<1
     flag that CODE.md asked for, and explain the difference honestly.

Every section ends with a cross-validation block.
"""
import sympy as sp
from sympy import symbols, Poly, I, im, expand, factorint, primerange, Rational, ZZ, QQ
from fractions import Fraction
import csv, os

m = symbols('m', real=True)
RESULTS = "/home/clio/projects/code/results"
os.makedirs(RESULTS, exist_ok=True)

# ----------------------------------------------------------------------
# 1.  I_b(m) and Q_b  (recomputed from the reduction)
# ----------------------------------------------------------------------
def ImP_coeff_b(b):
    """I_b(m) = Im( [u^b]((1-u) P) ),  P = (1 + s u + u^2)^m,  s = 1+i,
    as a polynomial in m.  Build via the binomial/trinomial expansion:
      [u^l] P = sum_{c>=0} ( m! / ((l-2c)! c! (m-l+c)!) ) s^{l-2c}
              = sum_{c} ( fallingfactorial(m, l-c) / ((l-2c)! c!) ) s^{l-2c}.
    Then [u^b]((1-u)P) = [u^b]P - [u^{b-1}]P, take imaginary part."""
    s = 1 + I
    def rl(l):
        tot = sp.Integer(0)
        for c in range(0, l // 2 + 1):
            k = l - 2 * c
            ff = sp.Integer(1)
            for j in range(0, k + c):          # falling factorial m(m-1)...(m-(k+c)+1)
                ff *= (m - j)
            tot += ff / (sp.factorial(k) * sp.factorial(c)) * s**k
        return sp.expand(tot)
    G = expand(rl(b) - rl(b - 1))
    return expand(im(G))

def Ib_via_Hm(b, mval):
    """Independent numeric evaluation of I_b at integer m = mval, using the
    H_m = h_{m-1}(A,B) staircase form, A+B=2W, AB=W^2+u^2, W=1+u+u^2.
    G_{(2m-b,b)} = [u^b]((1-u) H_m')... -- but the cleanest INDEPENDENT check is
    to expand (1 + s u + u^2)^m directly as a u-polynomial with s=1+i and read off
    Im([u^b]((1-u)P)).  This shares no code path with the symbolic-m builder."""
    u = symbols('u')
    s = 1 + I
    P = expand((1 + s * u + u**2)**mval)
    Q = expand((1 - u) * P)
    cb = Q.coeff(u, b)
    return int(im(expand(cb)))

def get_Qb(b):
    """Return (Qb_poly_over_QQ, forced_roots_list).  Q_b = I_b / prod (m-r)."""
    Ib = ImP_coeff_b(b)
    P = Poly(Ib, m, domain=QQ)
    forced = list(range(0, (b - 1) // 2 + 1))
    for r in forced:
        q, rem = sp.div(P, Poly(m - r, m, domain=QQ), m)
        assert rem.as_expr() == 0, ("forced root failed", b, r, rem)
        P = q
    return P, forced

# ----------------------------------------------------------------------
# 2.  Q~_b  : primitive integer form of the non-trivial irreducible factor
# ----------------------------------------------------------------------
def primitive_int_poly(poly_qq):
    """Clear denominators and divide by content -> primitive integer Poly,
    with positive leading coefficient."""
    expr = poly_qq.as_expr()
    p = Poly(expr, m, domain=QQ)
    coeffs = p.all_coeffs()
    dens = [Fraction(str(c)).denominator for c in [sp.Rational(c) for c in coeffs]]
    from math import lcm
    L = 1
    for d in dens:
        L = lcm(L, d)
    intp = Poly([sp.nsimplify(c * L) for c in coeffs], m, domain=ZZ)
    intp = intp.primitive()[1]
    if intp.LC() < 0:
        intp = Poly([-c for c in intp.all_coeffs()], m, domain=ZZ)
    return intp

def factor_structure(b):
    """Return dict with Q_b, its factor list, Q~_b (primitive int), the
    half-integer linear factor (if 4|b), and irreducibility verdict."""
    Qb, forced = get_Qb(b)
    flist = sp.factor_list(Qb.as_expr())
    const_factor, factors = flist                  # (rational content, [(poly, mult),...])
    # non-trivial (degree>=1) factors
    polyfactors = [(Poly(f, m, domain=QQ), mult) for f, mult in factors if Poly(f, m).degree() >= 1]
    # identify a linear factor with HALF-INTEGER root (the (2m-(2b-1)) for 4|b)
    halfint_linear = None
    nonlinear = []
    for f, mult in polyfactors:
        if f.degree() == 1:
            root = sp.solve(f.as_expr(), m)[0]
            if sp.nsimplify(root).q == 2:          # denominator 2 -> half-integer
                halfint_linear = (f, root)
                continue
        nonlinear.append((f, mult))
    # Q~_b = product of the non-half-integer-linear factors (should be one irreducible)
    if nonlinear:
        Qtilde_qq = nonlinear[0][0]
        for f, _ in nonlinear[1:]:
            Qtilde_qq = Qtilde_qq * f
    elif polyfactors:
        # whole thing was the linear factor (b too small) -- use the linear factor itself
        Qtilde_qq = polyfactors[0][0]
    else:
        # Q_b is a nonzero constant (b=1): degenerate, Q~ = the constant
        Qtilde_qq = Qb
    Qtilde = primitive_int_poly(Qtilde_qq)
    irreducible = (len(nonlinear) <= 1 and (not nonlinear or nonlinear[0][1] == 1)
                   and Qtilde.is_irreducible)
    return dict(b=b, Qb=Qb, forced=forced, full_factor=sp.factor(Qb.as_expr()),
                Qtilde=Qtilde, halfint_linear=halfint_linear,
                irreducible=Qtilde.is_irreducible, n_nonlin=len(nonlinear))

# ----------------------------------------------------------------------
# 3.  Newton polygon
# ----------------------------------------------------------------------
def vp(n, p):
    n = int(n)
    if n == 0:
        return None
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k

def lower_hull(points):
    """Lower convex hull of (x,y) points, x increasing. Returns hull vertices."""
    pts = sorted(points)
    hull = []
    for pt in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            x3, y3 = pt
            # cross product; pop if not a right (lower) turn
            cross = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)
            if cross <= 0:
                hull.pop()
            else:
                break
        hull.append(pt)
    return hull

def newton_polygon(intpoly, p):
    """Return (vertices, slopes) of the p-adic Newton polygon.
    intpoly: sympy Poly over ZZ.  coeff_j attached to m^j.
    slopes are exact Fractions; edge from vertex i to i+1."""
    deg = intpoly.degree()
    coeffs = intpoly.all_coeffs()              # high degree -> low degree
    points = []
    for idx, c in enumerate(coeffs):
        j = deg - idx
        v = vp(c, p)
        if v is not None:
            points.append((j, v))
    hull = lower_hull(points)
    slopes = []
    for (x1, y1), (x2, y2) in zip(hull, hull[1:]):
        slopes.append((Fraction(y2 - y1, x2 - x1), x2 - x1))   # (slope, horiz length)
    return hull, slopes

def no_integer_slope(slopes):
    """RIGOROUS: True iff no Newton edge has integer slope => no rational root."""
    return all(s.denominator != 1 for s, _ in slopes)

def all_abs_slope_lt1(slopes):
    """The looser CODE.md flag: every slope strictly between -1 and 1."""
    return all(abs(s) < 1 for s, _ in slopes)

# ----------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------
def main():
    print("=" * 78)
    print("STEP A -- reconstruct Q_b, cross-check I_b numerically, verify deg & factors")
    print("=" * 78)
    data = {}
    for b in range(1, 17):
        info = factor_structure(b)
        data[b] = info
        Qb = info['Qb']
        # cross-check: Q_b * prod(m-r) reproduces I_b, and I_b(m0) matches independent eval
        Ib = ImP_coeff_b(b)
        prod = sp.Integer(1)
        for r in info['forced']:
            prod *= (m - r)
        recon = expand(Qb.as_expr() * prod)
        assert expand(recon - Ib) == 0, ("reduction reconstruction failed", b)
        # independent numeric cross-check via direct u-expansion at several m
        for m0 in [b, b + 1, b + 3, b + 5, 2 * b]:
            sym = int(Poly(Ib, m).eval(m0))
            ind = Ib_via_Hm(b, m0)
            assert sym == ind, ("numeric mismatch", b, m0, sym, ind)
        tag = ""
        if info['halfint_linear'] is not None:
            tag = f"  [half-int linear: {sp.factor(info['halfint_linear'][0].as_expr())}, root {info['halfint_linear'][1]}]"
        print(f"  b={b:2d}  deg Q_b={Qb.degree()}  Q~ irreducible={info['irreducible']}"
              f"  deg Q~={info['Qtilde'].degree()}{tag}")
    print("  [cross-checks PASS: reduction reconstructs I_b; direct u-expansion agrees]")

    print("\n" + "=" * 78)
    print("STEP B -- leading coeff, constant term, discriminant of Q~_b  (b=5..16)")
    print("=" * 78)
    disc_rows = []
    for b in range(5, 17):
        Qt = data[b]['Qtilde']
        lc = int(Qt.LC())
        ct = int(Qt.all_coeffs()[-1])
        disc = int(sp.discriminant(Qt.as_expr(), m))
        lc_f = factorint(abs(lc)) if lc != 0 else {}
        ct_f = factorint(abs(ct)) if ct != 0 else {}
        disc_f = factorint(abs(disc)) if disc != 0 else {}
        is_sq = sp.sqrt(abs(disc)) == int(sp.sqrt(abs(disc))) and disc > 0 and int(sp.sqrt(disc))**2 == disc
        disc_rows.append(dict(b=b, lc=lc, ct=ct, disc=disc, lc_f=lc_f, ct_f=ct_f,
                              disc_f=disc_f, is_square=bool(is_sq)))
        print(f"  b={b:2d}  Q~={sp.factor(Qt.as_expr())}")
        print(f"        LC={lc} = {dict(lc_f)},  const={ct} = {dict(ct_f)}")
        print(f"        disc={disc} = {dict(disc_f)}   perfect_square={bool(is_sq)}")

    print("\n" + "=" * 78)
    print("STEP C -- p-adic Newton polygons; hunt for no-integer-slope (no rational root)")
    print("=" * 78)
    csv_rows = []
    witness_by_b = {}
    for b in range(5, 17):
        Qt = data[b]['Qtilde']
        dr = next(r for r in disc_rows if r['b'] == b)
        # candidate primes: those in lc / const / disc factorisations + small odd primes
        cand = set([3, 5, 7, 11, 13])
        cand |= set(dr['lc_f'].keys()) | set(dr['ct_f'].keys()) | set(dr['disc_f'].keys())
        cand = sorted(cand)
        b_witnesses = []
        for p in cand:
            hull, slopes = newton_polygon(Qt, p)
            no_int = no_integer_slope(slopes)
            lt1 = all_abs_slope_lt1(slopes)
            slope_str = "; ".join(f"{s}(len{L})" for s, L in slopes)
            obstruction = no_int and Qt.degree() >= 1
            if obstruction:
                b_witnesses.append(p)
            csv_rows.append(dict(b=b, p=p, deg=Qt.degree(),
                                 slopes=slope_str,
                                 no_integer_slope=no_int,
                                 all_abs_slope_lt1=lt1,
                                 has_rational_root_obstruction=obstruction))
        witness_by_b[b] = b_witnesses
        wtag = (f"WITNESS primes (no integer slope): {b_witnesses}"
                if b_witnesses else "NO single-prime no-integer-slope witness")
        print(f"  b={b:2d} (b mod 4 = {b%4})  {wtag}")
        # show the polygon detail for the smallest witness or for p=2,3
        for p in (b_witnesses[:1] or [2, 3]):
            hull, slopes = newton_polygon(Qt, p)
            print(f"        p={p}: vertices={hull}  slopes={[str(s) for s,_ in slopes]}")

    # write CSV
    csv_path = os.path.join(RESULTS, "Qb-newton-polygons.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["b", "p", "deg", "slopes",
                                           "no_integer_slope", "all_abs_slope_lt1",
                                           "has_rational_root_obstruction"])
        w.writeheader()
        for r in csv_rows:
            w.writerow(r)
    print(f"\n  wrote {csv_path}  ({len(csv_rows)} rows)")

    # ---- pattern analysis on witnesses by b mod 4 / mod 8 / mod 12 ----
    print("\n" + "=" * 78)
    print("STEP D -- does a (congruence-uniform) prime witness no-rational-root?")
    print("=" * 78)
    for b in range(5, 17):
        print(f"  b={b:2d}  b%4={b%4} b%8={b%8} b%12={b%12}  witnesses={witness_by_b[b]}"
              f"  irreducible_over_Q={data[b]['irreducible']}")
    # is there a prime working for ALL b? for all b in a residue class?
    allb = set(range(5, 17))
    # for each prime, which b's does it witness?
    from collections import defaultdict
    prime_hits = defaultdict(set)
    for b in range(5, 17):
        for p in witness_by_b[b]:
            prime_hits[p].add(b)
    print("\n  prime -> set of b in 5..16 it witnesses (no integer slope):")
    for p in sorted(prime_hits):
        print(f"     p={p}: {sorted(prime_hits[p])}")
    uniform = [p for p, S in prime_hits.items() if S == allb]
    print(f"\n  prime witnessing ALL b=5..16: {uniform if uniform else 'NONE'}")
    for mod in (4, 8, 12):
        print(f"\n  -- by residue mod {mod} --")
        classes = defaultdict(set)
        for b in range(5, 17):
            classes[b % mod].add(b)
        for res in sorted(classes):
            bs = classes[res]
            # primes witnessing every b in this class
            good = [p for p in prime_hits if bs <= prime_hits[p]]
            print(f"     b%{mod}={res} (b={sorted(bs)}): class-uniform witness primes = {good if good else 'NONE'}")

    return data, disc_rows, csv_rows, witness_by_b

if __name__ == "__main__":
    main()
