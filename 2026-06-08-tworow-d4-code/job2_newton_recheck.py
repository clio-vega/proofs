"""
job2_newton_recheck.py  --  CORRECTED Newton-polygon obstruction for Job 2.

The right single-prime certificate for "no rational root" of an integer polynomial f:
  a rational root r is a p-adic root whose valuation v_p(r) equals (minus) an edge slope
  of the p-adic Newton polygon of f, hence an INTEGER.  Therefore

      every edge slope of the p-adic Newton polygon is NON-INTEGRAL
        ==>  f has no root in Q_p  ==>  f has no rational root.

(The earlier "all |slope| < 1" flag is wrong: slope 0 is an integer < 1 and corresponds to
unit roots, which can be rational.  We use the correct non-integral-slope test here.)

We re-examine, for b == 2,3 (mod 4):
  (M) the monic q_b itself -- does some prime p | const(q_b) give all-non-integral slopes?
      (this would OVERTURN last cycle's "monic always has a slope-0 edge" verdict);
  (R) the reciprocal m^deg q_b(1/m) (non-monic, leading = |const term|);
  (S) the scaled q_b(p*m).
For each we report any prime giving an all-non-integral-slope obstruction.
"""

import sympy as sp
from sympy import Poly, primerange
from fractions import Fraction
from dfour_tworow import qb_monic, m

PRIMES = list(primerange(2, 200))


def newton_slopes(poly, p):
    coeffs = poly.all_coeffs()[::-1]
    pts = []
    for i, a in enumerate(coeffs):
        a = int(a)
        if a == 0:
            continue
        v, aa = 0, a
        while aa % p == 0:
            aa //= p
            v += 1
        pts.append((i, v))
    pts.sort()
    hull = []
    for pt in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (x2 - x1) * (pt[1] - y1) - (y2 - y1) * (pt[0] - x1) <= 0:
                hull.pop()
            else:
                break
        hull.append(pt)
    edges = []
    for (x1, y1), (x2, y2) in zip(hull, hull[1:]):
        edges.append((Fraction(y2 - y1, x2 - x1), x2 - x1))  # (slope, horizontal length)
    return edges


def all_nonintegral(edges):
    return bool(edges) and all(s.denominator != 1 for s, _ in edges)


def reciprocal(poly):
    return Poly(poly.all_coeffs()[::-1], m)


def main():
    bs = [b for b in range(6, 43) if b % 4 in (2, 3)]
    print("=" * 78)
    print("CORRECTED Newton obstruction (all-slopes-NON-INTEGRAL) for b == 2,3 (mod 4)")
    print("=" * 78)
    print(f"{'b':>3} {'%4':>3} {'deg':>4}  monic q_b: primes w/ all-noninteg slopes |"
          f" reciprocal | scaled")
    any_monic = False
    for b in bs:
        q = qb_monic(b)
        rec = reciprocal(q)
        mon_p, rec_p, sc_p = [], [], []
        for p in PRIMES:
            if all_nonintegral(newton_slopes(q, p)):
                mon_p.append(p)
            if all_nonintegral(newton_slopes(rec, p)):
                rec_p.append(p)
            qs = Poly(q.as_expr().subs(m, p * m), m)
            if all_nonintegral(newton_slopes(qs, p)):
                sc_p.append(p)
        if mon_p:
            any_monic = True
        print(f"{b:3d} {b%4:3d} {q.degree():4d}  {str(mon_p):>34s} | {str(rec_p):>10s} | {sc_p}")
    print()
    print("INTERPRETATION:")
    print("  * A non-empty 'monic' entry = a single odd prime that ALONE proves (law for b)")
    print("    via a non-integral Newton polygon -- a genuine local certificate.")
    print("  *", "monic certificates DO exist (overturns last cycle's slope-0 verdict)."
          if any_monic else
          "no monic certificate; monic Newton polygon stays inconclusive (slope-0 edge).")
    print("  * reciprocal / scaled break monicity and may certify where the monic cannot.")


if __name__ == "__main__":
    main()
