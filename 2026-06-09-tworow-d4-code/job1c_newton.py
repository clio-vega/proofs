"""JOB 1 (Newton-polygon-of-COEFFICIENTS) — the direct recurrence-Newton lever.

Hajir / Filaseta-type irreducibility comes from the p-adic Newton polygon of the
COEFFICIENTS of Q_b (not the discriminant): if for some prime p the lower hull of
the points (i, v_p(a_i)) is a single segment from (0, v_p(a_0)) to (deg, 0) whose
only lattice points are the endpoints (gcd(v_p(a_0), deg) = 1, a "pure" / Dumas
slope), then Q_b is irreducible over Q_p, hence over Q.  (Eisenstein is the case
v_p(a_0)=1.)  More generally a Newton polygon all of whose segment lengths share
no common proper divisor restricts the factor degrees.

We scan small primes and report, per b, the best Newton-polygon structure found.
a_i are the coefficients of Q_b, a_deg = 1 (monic) so v_p(a_deg)=0 always.
"""
import sys
sys.set_int_max_str_digits(500000)
from math import gcd
from sympy import primerange
from _qbcore import load_cache

def vp(n, p):
    if n == 0:
        return None  # +inf
    k = 0
    n = abs(int(n))
    while n % p == 0:
        n //= p; k += 1
    return k

def lower_hull(points):
    """Lower convex hull of (x,y) points, x increasing, ignoring y=None (inf)."""
    pts = sorted((x, y) for x, y in points if y is not None)
    hull = []
    for p in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            x3, y3 = p
            # cross product to keep lower hull (turn must be counterclockwise-ish)
            if (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1) <= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    return hull

def newton_segments(coeffs, p):
    """coeffs high->low. Return list of (run_length, rise) segments of lower hull
    over points (i, v_p(a_i)) with i = degree of that term (0..deg)."""
    deg = len(coeffs) - 1
    pts = []
    for idx, c in enumerate(coeffs):
        i = deg - idx
        pts.append((i, vp(c, p)))
    hull = lower_hull(pts)
    segs = []
    for k in range(len(hull)-1):
        (x1, y1), (x2, y2) = hull[k], hull[k+1]
        segs.append((x2-x1, y1-y2))  # horizontal run, vertical drop
    return hull, segs

def is_pure_irreducible(coeffs, p):
    """Single segment from (0,v) to (deg,0) with gcd(v,deg)=1 -> irreducible."""
    deg = len(coeffs)-1
    if vp(coeffs[0], p) != 0:   # leading must be a p-unit (monic -> always true)
        return False
    v0 = vp(coeffs[-1], p)      # constant term
    if v0 is None:
        return False
    hull, segs = newton_segments(coeffs, p)
    if len(segs) != 1:
        return False
    run, rise = segs[0]
    return run == deg and rise == v0 and gcd(rise, deg) == 1

def main():
    bmax = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    cache = load_cache()
    primes = list(primerange(2, 500))
    print(f"{'b':>4} {'deg':>4} {'pure-irred prime':>16} {'min#segments(over p)':>22}")
    n_pure = 0
    rows = []
    for b in sorted(cache):
        if b > bmax:
            continue
        coeffs = [int(c) for c in cache[b].all_coeffs()]
        deg = len(coeffs)-1
        pure_p = None
        min_segs = 999
        for p in primes:
            if is_pure_irreducible(coeffs, p):
                pure_p = p
                break
            _, segs = newton_segments(coeffs, p)
            # count nontrivial structure: segments whose runs share no common divisor
            runs = [r for r, _ in segs]
            if runs:
                min_segs = min(min_segs, len(segs))
        if pure_p:
            n_pure += 1
        print(f"{b:>4} {deg:>4} {str(pure_p):>16} {min_segs:>22}")
        rows.append((b, deg, pure_p, min_segs))
    print(f"\nb with a PURE (Dumas/Eisenstein) irreducibility prime p<500: "
          f"{n_pure}/{len(rows)}")
    pures = [(b, pp) for b, d, pp, ms in rows if pp]
    print("examples:", pures[:15])

if __name__ == "__main__":
    main()
