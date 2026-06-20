#!/usr/bin/env python3
"""Job C (2026-06-20 CODE): the LGV / Gaussian-binomial q-determinant of M_j.

Builds on 2026-06-19 route-1 negative: the dimension q-hook collapses to 0 at the
boundary; the correct q-lift is the Gaussian-binomial determinant of the CLOSED FORM
   M_j = Sum_{sigma in S_3} sgn(sigma) D_j(L0+sigma1, L1+sigma2, L2+sigma3),
   D_j(p,q,r) = [s^p t^q u^r] e1^{2(m-j)} e2^j
             = Sum_{i+k+l=j} multinom(j;i,k,l) * multinom(N; p-(i+k), q-(i+l), r-(k+l))
   (N = 2(m-j)),   L = (a-1, b-2, c-3).   (Verified vs MN in threerow-c1/dj.py.)

Each multinomial is a lattice-path count in a box (LGV), so M_j is a signed sum of
products of ordinary binomials.  Item 1 of the plan: q-lift every binomial to its
Gaussian analogue [n choose k]_q, check the q=1 value recovers M_j (sanity), and test
whether the result is q-POSITIVE coefficient-wise (LGV => no cancellation) -- or report
exactly where positivity fails.

We test two q-liftings:
  (A) NAIVE: replace each multinomial by its q-multinomial, no extra q-powers.
  (B) CHARGE: the same with the standard q-multinomial nesting (already a product of
      Gaussian binomials carrying the inversion statistic).  [A and B coincide here
      because the q-multinomial IS the charge-weighted box-path count -- so a single
      honest construction; we test its positivity.]
"""
import sys
sys.path.insert(0, '/home/clio/projects/code/threerow-c3')
from mn import Mj
import sympy as sp

q = sp.symbols('q')

def qint(n):
    if n <= 0:
        return sp.Integer(0) if n < 0 else sp.Integer(1) if n == 0 else None
    return sp.expand(sum(q**i for i in range(n)))

def qfact(n):
    r = sp.Integer(1)
    for t in range(1, n + 1):
        r = sp.expand(r * qint(t))
    return r

def qbinom(n, k):
    if k < 0 or k > n or n < 0:
        return sp.Integer(0)
    # Gaussian binomial is a polynomial; cancel performs the exact division.
    return sp.expand(sp.cancel(qfact(n) / (qfact(k) * qfact(n - k))))

def qmultinom(N, parts):
    """[N; parts]_q = product of nested Gaussian binomials (q-multinomial)."""
    if any(p < 0 for p in parts) or sum(parts) != N:
        return sp.Integer(0)
    r = sp.Integer(1)
    rem = N
    for p in parts[:-1]:
        r = sp.expand(r * qbinom(rem, p))
        rem -= p
    return r

def Dj_q(p, r, s, j, m):
    """q-lift of D_j(p,q,r); 's' here is the third exponent (avoid name clash with q)."""
    if p < 0 or r < 0 or s < 0 or p + r + s != 2 * m:
        return sp.Integer(0)
    N = 2 * (m - j)
    tot = sp.Integer(0)
    for i in range(j + 1):
        for k in range(j - i + 1):
            l = j - i - k
            coef = qmultinom(j, [i, k, l])
            tri = qmultinom(N, [p - (i + k), r - (i + l), s - (k + l)])
            if tri != 0 and coef != 0:
                tot = sp.expand(tot + coef * tri)
    return tot

S3 = [(1, (1, 2, 3)), (-1, (1, 3, 2)), (-1, (2, 1, 3)),
      (1, (2, 3, 1)), (1, (3, 1, 2)), (-1, (3, 2, 1))]

def Mj_q(lam, j, m):
    a, b, c = lam
    L = [a - 1, b - 2, c - 3]
    tot = sp.Integer(0)
    for sg, sigma in S3:
        idx = [L[0] + sigma[0], L[1] + sigma[1], L[2] + sigma[2]]
        tot = sp.expand(tot + sg * Dj_q(idx[0], idx[1], idx[2], j, m))
    return sp.expand(tot)

# ---------------------------------------------------------------------------
def v2(n):
    n = abs(int(n))
    return 10**9 if n == 0 else (n & -n).bit_length() - 1

def cyclotomic_tower_v2(poly):
    """Sum_{r>=1} mult_{Phi_{2^r}}(poly), i.e. 2-content visible in cyclotomic factors.
    = Sum over factors Phi_{2^r}^{e} of e (each Phi_{2^r}(1)=2). Returns (tower, leak)
    where leak = v2(poly(1)) - tower (2-content NOT in 2-power-cyclotomic factors)."""
    P = sp.Poly(sp.expand(poly), q)
    if P.is_zero:
        return (None, None)
    val1 = int(P.eval(1))
    if val1 == 0:
        return (None, None)
    tower = 0
    fac = sp.factor_list(P.as_expr())
    # iterate factors; detect Phi_{2^r}: cyclotomic poly of order 2^r is q^{2^{r-1}}+1
    const, factors = fac
    for base, e in factors:
        bp = sp.Poly(sp.expand(base), q)
        # check if base == q^{2^{r-1}} + 1 for some r
        d = bp.degree()
        # Phi_{2^r} has degree 2^{r-1}
        if d >= 1 and (d & (d - 1)) == 0:  # d is a power of 2
            cand = q**d + 1
            if sp.expand(base - cand) == 0:
                tower += e
    leak = v2(val1) - tower
    return tower, leak

if __name__ == "__main__":
    print("=" * 72)
    print("Job C: Gaussian-binomial (LGV) q-determinant of M_j")
    print("=" * 72)

    # ---- item 1a: q=1 sanity vs MN, INTERIOR and BOUNDARY indices ----
    print("\n[1a] q=1 recovery of M_j vs MN (interior j<=c AND boundary j=b+i):")
    bad = 0; tested = 0
    shapes = [(5, 4, 3), (6, 5, 4), (7, 5, 4), (8, 6, 4), (9, 7, 5), (10, 8, 6)]
    for lam in shapes:
        a, b, c = lam
        m = (a + b + c) // 2
        if (a + b + c) % 2 != 0:
            continue
        for j in range(0, b + 3):          # include boundary j=b+1,b+2
            mq = Mj_q(lam, j, m)
            v_at1 = int(sp.Poly(mq, q).eval(1)) if mq != 0 else 0
            mn = Mj(lam, j, m)
            tested += 1
            if v_at1 != mn:
                bad += 1
                if bad <= 6:
                    print(f"    q=1 MISMATCH {lam} j={j}: qdet(1)={v_at1} MN={mn}")
    print(f"    q=1 vs MN: tested={tested}  mismatches={bad}")

    # ---- item 1b: q-positivity test ----
    print("\n[1b] q-positivity (LGV: signed S_3 sum should telescope to >=0 coeffs):")
    pos_fail = 0; pos_ok = 0; examples = []
    for lam in shapes:
        a, b, c = lam
        m = (a + b + c) // 2
        if (a + b + c) % 2 != 0:
            continue
        for j in range(0, b + 1):
            mq = Mj_q(lam, j, m)
            if mq == 0:
                continue
            coeffs = sp.Poly(mq, q).all_coeffs()
            neg = [co for co in coeffs if co < 0]
            if neg:
                pos_fail += 1
                if len(examples) < 6:
                    examples.append((lam, j, min(coeffs)))
            else:
                pos_ok += 1
    print(f"    q-positive cases: {pos_ok}   NON-positive (has a negative coeff): {pos_fail}")
    if examples:
        print(f"    first positivity failures (lam,j,min-coeff): {examples}")

    # ---- item 2: cyclotomic tower vs v2 at the boundary deficit ----
    print("\n[2] cyclotomic 2^r-tower vs v2(M_{b+i})  (leak = non-cyclotomic 2-content):")
    print("    lam        j=b+i  M_j        v2(M)  tower  leak")
    for lam in [(6, 5, 4), (8, 6, 4), (9, 7, 5), (10, 8, 6)]:
        a, b, c = lam
        if (a + b + c) % 2 != 0:
            continue
        m = (a + b + c) // 2
        for i in (1, 2):
            j = b + i
            mq = Mj_q(lam, j, m)
            if mq == 0:
                print(f"    {lam}  j={j}: q-det == 0  (M_j={Mj(lam,j,m)})")
                continue
            mn = Mj(lam, j, m)
            tower, leak = cyclotomic_tower_v2(mq)
            print(f"    {lam}  j={j}   {mn:<10} {v2(mn):<5}  {tower}     {leak}")

    print("\n[done]")
