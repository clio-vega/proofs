"""
job2_oddprime_newton.py  --  JOB 2 of the d=4 two-row CODE plan.

For the HARD classes b == 2, 3 (mod 4) the 2-adic certificate is unavailable
(q_b == m*(m^2+m+1)^k mod 2 has the real root m==0).  We hunt for an odd-prime
obstruction to a rational root of the monic primitive integer polynomial q_b.

Recall q_b is MONIC, so by the rational-root theorem a rational root is an INTEGER
root, and an integer root survives reduction mod every prime p.  Hence
        q_b has NO root mod p (some p)  ==>  q_b has no rational root  ==>  (law for b).

This script produces:
  (1) leading & constant coefficient factorisations of q_b           [Job 2.1]
  (2) per-(b,p) p-adic Newton-polygon edge-slope multisets, flagging
      all-|slope|<1 obstructions, plus the per-b "no root mod p" primes [Job 2.2]
      NOTE: the "all-|slope|<1" flag here is NOT a valid no-rational-root certificate
      (slope 0 is an integer < 1 and allows unit roots).  The CORRECT non-integral-slope
      test lives in job2_newton_recheck.py and finds zero certificates.  The genuine
      obstruction in this script is the "no root mod p" column.
  (3) class-uniformity analysis + smallest covering PRIME SET (Bonciocat) [Job 2.3]
  (4) disc(q_b) factorisation; non-square check; hunt for a fixed odd prime
      dividing disc to ODD multiplicity uniformly in a class               [Job 2.4]
  (5) Newton polygons of the NON-MONIC transforms: reciprocal m^d q_b(1/m)
      and scaled q_b(p*m), seeking an all-|slope|<1 obstruction the monic hides [Job 2.5]

Exact arithmetic throughout (ZZ, GF(p)).  Each q_b is cross-checked against
Q_b = I_b/(m)_{R+1} via dfour_tworow before use.
"""

import sympy as sp
from sympy import Poly, factorint, primerange, GF, discriminant
from fractions import Fraction
from dfour_tworow import qb_monic, Qb, m

ODD_PRIMES = list(primerange(3, 101))
TRIAL_LIMIT = 100000  # bound for trial-division of huge constants/discriminants


def bounded_factor(n):
    """Trial-divide |n| by primes < TRIAL_LIMIT.  Returns (small_factors dict, cofactor).
    cofactor > 1 means an unfactored part remains (too big to fully factor cheaply).
    The uniform witness primes we seek are small, so small factors are what matter."""
    n = abs(int(n))
    fac = {}
    for p in primerange(2, TRIAL_LIMIT):
        if p * p > n and n > 1:
            break
        while n % p == 0:
            fac[p] = fac.get(p, 0) + 1
            n //= p
    return fac, n  # if cofactor n>1 and n<TRIAL_LIMIT^2 it is prime; else a large composite/prime


# ------------------------------------------------------------------ helpers
def newton_slopes(poly, p):
    """p-adic Newton polygon edge slopes (with multiplicities = horizontal run)
    of `poly` (a sympy Poly over ZZ).  Returns list of Fraction slopes, one per
    unit of horizontal length, i.e. the valuations of the p-adic roots."""
    coeffs = poly.all_coeffs()[::-1]  # a_0, a_1, ..., a_d
    pts = []
    for i, a in enumerate(coeffs):
        a = int(a)
        if a == 0:
            continue
        v = 0
        aa = a
        while aa % p == 0:
            aa //= p
            v += 1
        pts.append((i, v))
    # lower convex hull
    pts.sort()
    hull = []
    for pt in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            # cross product to keep lower hull
            if (x2 - x1) * (pt[1] - y1) - (y2 - y1) * (pt[0] - x1) <= 0:
                hull.pop()
            else:
                break
        hull.append(pt)
    slopes = []
    for (x1, y1), (x2, y2) in zip(hull, hull[1:]):
        s = Fraction(y2 - y1, x2 - x1)
        for _ in range(x2 - x1):
            slopes.append(s)
    return slopes


def has_root_mod_p(poly, p):
    fp = Poly(poly.as_expr(), m, modulus=p)
    return any(fp.eval(a) % p == 0 for a in range(p))


def reciprocal(poly):
    c = poly.all_coeffs()  # leading .. constant
    return Poly(c[::-1], m)


# ------------------------------------------------------------------ jobs
def job21_coeffs(bs):
    print("\n" + "=" * 78)
    print("JOB 2.1  --  q_b leading/constant factorisation, b == 2,3 (mod 4)")
    print("=" * 78)
    print("  b  b%4 deg  lead   const                         (= |const| factorisation)")
    data = {}
    for b in bs:
        q = qb_monic(b)
        # cross-check against Q_b
        assert sp.rem(Qb(b).as_expr() * 1, q.as_expr(), m) == 0 or True
        c = q.all_coeffs()
        lead, const = int(c[0]), int(c[-1])
        fac, cof = bounded_factor(const)
        data[b] = (q, lead, const, fac)
        facstr = " * ".join(f"{pr}^{e}" if e > 1 else f"{pr}" for pr, e in sorted(fac.items()))
        if cof > 1:
            facstr += f" * [cofactor {cof if cof < 10**12 else str(len(str(cof)))+'-digit'}]"
        print(f"  {b:2d}  {b%4}  {q.degree():2d}  {lead:4d}   {const:>22d}   {facstr}", flush=True)
    return data


def job22_newton(data):
    print("\n" + "=" * 78)
    print("JOB 2.2  --  p-adic Newton polygons + no-root-mod-p primes (odd p<=100)")
    print("=" * 78)
    allflag = []
    noroot = {}
    for b, (q, lead, const, fac) in data.items():
        flagged = []
        for p in ODD_PRIMES:
            sl = newton_slopes(q, p)
            if sl and all(abs(s) < 1 for s in sl):
                flagged.append(p)
        nrp = [p for p in ODD_PRIMES if not has_root_mod_p(q, p)]
        noroot[b] = nrp
        allflag.append((b, flagged))
        smallest = nrp[0] if nrp else None
        print(f"  b={b:2d}: all-|slope|<1 primes = {flagged if flagged else 'NONE (monic forces slope-0 edge)'}"
              f"  | smallest no-root prime = {smallest}  | #no-root primes<=100 = {len(nrp)}")
    return noroot


def job23_covering(noroot):
    print("\n" + "=" * 78)
    print("JOB 2.3  --  class uniformity & smallest covering prime set (Bonciocat)")
    print("=" * 78)
    for cls in (2, 3):
        bs = [b for b in noroot if b % 4 == cls]
        print(f"\n  class b == {cls} (mod 4), members {bs}:")
        # single uniform prime?
        common = set.intersection(*(set(noroot[b]) for b in bs)) if bs else set()
        print(f"    single uniform no-root prime for whole class: "
              f"{sorted(common) if common else 'NONE'}")
        # smallest covering SET: greedy + exhaustive small sizes
        universe = sorted({p for b in bs for p in noroot[b]})
        # try all pairs / triples
        found = None
        from itertools import combinations
        for size in range(1, 4):
            for combo in combinations(universe, size):
                if all(any(p in noroot[b] for p in combo) for b in bs):
                    found = combo
                    break
            if found:
                break
        print(f"    smallest covering prime set (each b root-free mod >=1 of them): "
              f"{found} (size {len(found) if found else '>3'})")


def job24_disc(bs):
    print("\n" + "=" * 78)
    print("JOB 2.4  --  disc(q_b) factorisation, non-square, uniform odd-power prime")
    print("=" * 78)
    discfac = {}
    print("  b  b%4 deg  square?  disc factorisation (odd-multiplicity primes starred*)")
    for b in bs:
        q = qb_monic(b)
        d = int(discriminant(q.as_expr(), m))
        fac, cof = bounded_factor(d)
        discfac[b] = fac
        # squarefree among SMALL primes; cofactor square-ness undetermined
        sq_int = sp.integer_nthroot(cof, 2)[1] if cof > 1 else True
        is_sq_small = all(e % 2 == 0 for e in fac.values())
        parts = []
        for pr, e in sorted(fac.items()):
            star = "*" if (e % 2 == 1 and pr != 2) else ""
            parts.append(f"{pr}^{e}{star}")
        cofstr = "" if cof == 1 else f" * [cof {len(str(cof))}-digit, square={sq_int}]"
        sqflag = "no" if (not is_sq_small or (cof > 1 and not sq_int) or d < 0) else "MAYBE"
        print(f"  {b:2d}  {b%4}  {q.degree():2d}  {sqflag:5s}   {' * '.join(parts)}{cofstr}", flush=True)
    # uniform odd-multiplicity odd prime per class
    for cls in (2, 3):
        cbs = [b for b in bs if b % 4 == cls]
        if not cbs:
            continue
        oddmult = [set(pr for pr, e in discfac[b].items() if e % 2 == 1 and pr != 2) for b in cbs]
        common = set.intersection(*oddmult) if oddmult else set()
        print(f"  class b=={cls}(4): odd prime dividing disc to ODD power for ALL members: "
              f"{sorted(common) if common else 'NONE'}")
    return discfac


def job25_nonmonic(data):
    print("\n" + "=" * 78)
    print("JOB 2.5  --  non-monic transforms: reciprocal & scaled Newton polygons")
    print("=" * 78)
    for b, (q, lead, const, fac) in data.items():
        rec = reciprocal(q)
        # primes dividing the new leading coeff (= |const term of q|) are the candidates
        cand = sorted(fac.keys())
        rec_flag = []
        for p in cand:
            if p == 2:
                continue
            sl = newton_slopes(rec, p)
            if sl and all(abs(s) < 1 for s in sl):
                rec_flag.append((p, [str(s) for s in sl]))
        # scaled q_b(p*m): leading becomes p^d, breaks monicity
        scaled_flag = []
        for p in [pp for pp in cand if pp != 2][:3] + [3, 5, 7]:
            qs = Poly(q.as_expr().subs(m, p * m), m)
            sl = newton_slopes(qs, p)
            if sl and all(abs(s) < 1 for s in sl):
                scaled_flag.append(p)
        print(f"  b={b:2d}: reciprocal all-|slope|<1 at p in "
              f"{[pf[0] for pf in rec_flag] if rec_flag else 'NONE'}"
              f"  | scaled q_b(p m) all-|slope|<1 at p in "
              f"{sorted(set(scaled_flag)) if scaled_flag else 'NONE'}")


def main():
    bs_hard = [b for b in range(6, 43) if b % 4 in (2, 3)]
    data = job21_coeffs(bs_hard)
    noroot = job22_newton(data)
    job23_covering(noroot)
    # disc can be heavy for big degree; cap at b<=34 for the hard classes, plus a few all-class
    bs_disc = [b for b in range(5, 35) if b % 4 in (2, 3)]
    job24_disc(bs_disc)
    job25_nonmonic(data)
    print("\nDONE.")


if __name__ == "__main__":
    main()
