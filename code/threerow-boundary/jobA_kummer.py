"""
jobA_kummer.py -- Job A step 4: Kummer cross-check of the even-slice content.

Lemma 0 writes M_{b+i} as a signed S3 sum of products of two multinomials
(see fast_alt.T).  We expand the master N_i as that explicit alternating sum and,
on the even slice at the content-achieving minimal point, count the 2-adic carries
(Kummer) of each surviving term.  Goal: see whether the even-slice content matches
a clean carry-count on the surviving terms (Kummer route) or whether it only
appears after cancellation (factorization route).

We report, per (c,i): the measured fixed-divisor content vs sum-of-factor-content,
and the v2 of each S3 term at a generic minimal point -- showing how much content
is "free" (every term already divisible) vs "cancellation-born".
"""
from math import factorial
from fast_alt import Malt, _mult3, _S3


def v2int(n):
    n = abs(int(n))
    if n == 0:
        return 10**9
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def term_v2s(a, b, c, j, m):
    """v2 of each S3 alternating term of M_j (before summation)."""
    s = 2 * m - 2 * j
    base = (a + 2, b + 1, c)
    out = []
    for (perm, sgn) in _S3:
        p, q, r = base[0] - perm[0], base[1] - perm[1], base[2] - perm[2]
        # T(p,q,r,j,s) is itself a sum; we report v2 of the whole T-term
        tot = 0
        if p >= 0 and q >= 0 and r >= 0:
            for al in range(j + 1):
                for be in range(j - al + 1):
                    ga = j - al - be
                    Pp, Qq, Rr = p - al - ga, q - al - be, r - be - ga
                    if Pp < 0 or Qq < 0 or Rr < 0 or Pp + Qq + Rr != s:
                        continue
                    tot += _mult3(j, al, be, ga) * _mult3(s, Pp, Qq, Rr)
        out.append((sgn, tot, v2int(tot) if tot else None))
    return out


def Ni_at(c, i, Pv, Bv, even):
    k = c - i
    bv = 2 * Bv if even else 2 * Bv + 1
    av = 2 * Pv + bv + c
    m = Pv + bv + c
    j = bv + i
    if bv < c + 1 or j > m:
        return None
    M = Malt(av, bv, c, j, m)
    den = bv - c + 1
    for t in range(2, i + 1):
        den *= (bv + t)
    num = M * factorial(c) * factorial(k)
    if num % den:
        return None
    return num // den, av, bv, m, j, M


CASES = [(6, 2, 4), (8, 2, 6), (7, 2, 5)]

if __name__ == "__main__":
    for (c, i, k) in CASES:
        print("=" * 78)
        print("c=%d i=%d k=%d  Kummer/term analysis on EVEN slice" % (c, i, k))
        # find a content-achieving minimal even-slice point and a few neighbours
        best = None
        for Bv in range(c, c + 40):
            for Pv in range(0, 20):
                r = Ni_at(c, i, Pv, Bv, True)
                if r is None:
                    continue
                Nval = r[0]
                vv = v2int(Nval)
                if best is None or vv < best[0]:
                    best = (vv, Pv, Bv, r)
        vv, Pv, Bv, (Nval, av, bv, m, j, M) = best
        print("  minimal point P=%d B=%d -> a=%d b=%d m=%d j=%d" % (Pv, Bv, av, bv, m, j))
        print("  v2(N_i)=%d  v2(M_{b+i})=%d  v2(c!k!)=%d  v2(den)=%d"
              % (vv, v2int(M), v2int(factorial(c) * factorial(k)),
                 v2int(Nval and (factorial(c) * factorial(k) * M)) - v2int(Nval)))
        print("  S3 terms of M_{b+i}=Malt (sgn, value, v2):")
        for (sgn, tval, tv2) in term_v2s(av, bv, c, j, m):
            print("     sgn=%+d  v2=%s  val=%d" % (sgn, tv2, tval))
        minterm = min(v2int(t) for (_, t, _) in term_v2s(av, bv, c, j, m) if t)
        print("  min term v2 = %d ; v2(M)=%d  => cancellation lifts v2 by %d"
              % (minterm, v2int(M), v2int(M) - minterm))
