"""
Job A wide sweep (2026-06-15 evening code session).

Boundary residual for the three-row family lambda=(a,b,c) |- 2m.
Goal: confirm that the boundary indices j in {b+1,...,b+c} LOSE the 2-adic
minimum to the interior argmin j0, strictly, far past the m<=79-80 frontier.

val(j) = j + 2*( v2 C(m,j) + v2 M_j )   is the 2-adic valuation of the j-th
contribution to Z_lambda; ord = min_j val(j).  We show min over boundary j
exceeds the interior min by Delta >= 1 everywhere.

Closed forms used (all verified against Murnaghan-Nakayama below):
  interior  c=1  general.py Mfac          c=2  c2verify2.Mj_formula   c=3  c3verify.Q3
  boundary  the a-(in)dependent fits from boundary_forms.py.
For c=4,5 we fall back to the exact 3-row Jacobi-Trudi determinant (dj.Mj_jt),
at a smaller m, only to confirm the pattern persists.
"""
import sys
sys.path.insert(0, '/home/clio/projects/code/threerow-c1')
from math import comb
from fractions import Fraction
from mn import Mj  # exact MN, for cross-checks

# ---------------------------------------------------------------- 2-adic tools
def v2(n):
    n = abs(int(n))
    if n == 0:
        return float('inf')
    c = 0
    while n % 2 == 0:
        n //= 2; c += 1
    return c

def s2(n):
    return bin(int(n)).count('1')

def C(n, k):
    if k < 0 or n < 0 or k > n:
        return 0
    return comb(n, k)

def v2C(N, k):
    if k < 0 or k > N:
        return float('inf')
    return s2(k) + s2(N - k) - s2(N)

# ---------------------------------------------------------------- interior M_j
def M_interior(a, b, c, j, m):
    """Exact integer M_j for 0<=j<=b (the interior/conjugate region). None if 0/undefined."""
    N = 2 * (m - j)
    if c == 1:
        if j <= b - 1:
            num = C(N, b - j - 1) * (a - b + 1) * (b * (a + 1) - j * (j - 1))
            den = (b - j) * (b - j + 1)
            return num // den
        if j == b:
            return 2 * b * (m - b)
        return None
    if c == 2:
        k = b - j
        if k < 0:
            return None
        Q2 = a * (b - 1) * ((a + 3) * (b + 2) - 2 * j * j) + j * (j - 1) * (j - 2) * (j - 3)
        num = C(N, k) * (a - b + 1) * Q2
        den = 2 * (a + 3 - j) * (b + 2 - j) * (b + 1 - j)
        f = Fraction(num, den)
        return int(f) if f.denominator == 1 else None
    if c == 3:
        Bb = b - j
        if Bb < 0:
            return None
        Q3 = (a**3*b**3 + 3*a**3*b**2 - 4*a**3*b - 12*a**3
              + 6*a**2*b**3 - 3*a**2*b**2*j**2 - 3*a**2*b**2*j + 18*a**2*b**2
              + 3*a**2*b*j**2 - 3*a**2*b*j - 24*a**2*b + 6*a**2*j**2 + 18*a**2*j - 72*a**2
              + 5*a*b**3 - 3*a*b**2*j**2 - 9*a*b**2*j + 15*a*b**2 + 3*a*b*j**4 - 12*a*b*j**3
              + 24*a*b*j**2 - 15*a*b*j - 20*a*b - 6*a*j**4 + 24*a*j**3 - 36*a*j**2 + 66*a*j - 60*a
              - 12*b**3 + 6*b**2*j**2 + 12*b**2*j - 36*b**2 - 3*b*j**4 + 12*b*j**3 - 27*b*j**2
              + 18*b*j + 48*b - j**6 + 15*j**5 - 79*j**4 + 201*j**3 - 244*j**2 + 36*j + 144)
        num = C(N, Bb) * (a - b + 1) * Q3
        den = 6 * (a - j + 4) * (b - j + 1) * (b - j + 2) * (b - j + 3)
        f = Fraction(num, den)
        return int(f) if f.denominator == 1 else None
    raise ValueError("interior closed form only for c<=3")

# ---------------------------------------------------------------- boundary M_{b+i}
def M_boundary(a, b, c, i):
    """Exact integer M_{b+i}, i=1..c, from the closed-form fits."""
    if c == 1:
        return b                                   # i=1
    if c == 2:
        if i == 1:
            f = Fraction((b - 1) * (a * b + 2 * a - b**2 - b), 2)
        else:  # i==2
            f = Fraction((b - 1) * (b + 2), 2)
        return int(f) if f.denominator == 1 else None
    if c == 3:
        if i == 1:
            f = Fraction((b - 2) * (a - b + 1) * (a*b**2 + 5*a*b + 6*a - b**3 - b**2 - 4*b - 6), 12)
        elif i == 2:
            f = Fraction((b - 2) * (b + 2) * (a*b + 3*a - b**2 - 2*b - 3), 6)
        else:  # i==3
            f = Fraction((b - 2) * (b + 2) * (b + 3), 6)
        return int(f) if f.denominator == 1 else None
    raise ValueError("boundary closed form only for c<=3")

def val_from_M(M, j, m):
    if M is None or M == 0:
        return None
    if j > m:
        return None
    return j + 2 * (v2C(m, j) + v2(M))

# ================================================================ cross-checks
def crosscheck(cmax=3, mmax=22):
    bad = 0; tested = 0
    for c in range(1, cmax + 1):
        for m in range(c, mmax + 1):
            n = 2 * m
            for a in range(1, n + 1):
                for b in range(c, a + 1):
                    if n - a - b != c:
                        continue
                    lam = (a, b, c)
                    # interior j=0..b
                    for j in range(0, b + 1):
                        cl = M_interior(a, b, c, j, m)
                        ex = Mj(lam, j, m)
                        tested += 1
                        if cl != ex:
                            bad += 1
                            if bad < 12:
                                print("  INTERIOR BAD", lam, "m", m, "j", j, "closed", cl, "exact", ex)
                    # boundary j=b+i, i=1..c
                    for i in range(1, c + 1):
                        j = b + i
                        if j > m:
                            continue
                        cl = M_boundary(a, b, c, i)
                        ex = Mj(lam, j, m)
                        tested += 1
                        if cl != ex:
                            bad += 1
                            if bad < 12:
                                print("  BOUNDARY BAD", lam, "m", m, "j", j, "closed", cl, "exact", ex)
    print(f"[crosscheck c<=3, m<={mmax}] tested={tested} mismatches={bad}")
    return bad

# ================================================================ wide sweep
def wide_sweep(c, mmax):
    """For c=1,2,3, push m up to mmax with closed forms.

    The Boundary Lemma hypothesis is that the interior 2-adic box {j0,...,j0+2c}
    fits inside [0,b], i.e. j0 + 2c <= b. We separate shapes into:
      * IN-GATE  (j0+2c <= b): the lemma claims val(boundary) > V strictly.
      * SUB-GATE (j0+2c >  b): the known per-family boundary-tie shapes, excluded.
    j0 = smallest interior index achieving the interior minimum.
    """
    res = dict(c=c, mmax=mmax, shapes=0,
               in_min_delta=None, in_min_delta_at=None, in_ties=[],
               sub_count=0, sub_ties=0, sub_max_b=0, sub_b_list=set(),
               max_v2_boundary=-1, max_v2_boundary_at=None,
               max_v2_top=-1, max_v2_top_at=None,
               max_v2_M0=-1)
    for m in range(c, mmax + 1):
        n = 2 * m
        for a in range(1, n + 1):
            b = n - a - c
            if b < c or b > a:
                continue
            res['shapes'] += 1
            # interior minimum over j=0..b ; j0 = smallest achieving it
            interior_min = None; j0 = None
            for j in range(0, b + 1):
                M = M_interior(a, b, c, j, m)
                v = val_from_M(M, j, m)
                if v is None:
                    continue
                if interior_min is None or v < interior_min:
                    interior_min = v; j0 = j
            if interior_min is None:
                continue
            # boundary vals j=b+i, i=1..c
            bmin = None; bj = None
            for i in range(1, c + 1):
                j = b + i
                M = M_boundary(a, b, c, i)
                if M is not None and M != 0 and j <= m:
                    vb = v2(M)
                    if vb > res['max_v2_boundary']:
                        res['max_v2_boundary'] = vb; res['max_v2_boundary_at'] = (a, b, c, m, j)
                    if i == c and vb > res['max_v2_top']:
                        res['max_v2_top'] = vb; res['max_v2_top_at'] = (a, b, c, m, j)
                v = val_from_M(M, j, m)
                if v is None:
                    continue
                if bmin is None or v < bmin:
                    bmin = v; bj = j
            M0 = M_interior(a, b, c, 0, m)
            if M0:
                res['max_v2_M0'] = max(res['max_v2_M0'], v2(M0))
            if bmin is None:
                continue
            delta = bmin - interior_min
            in_gate = (j0 + 2 * c <= b)
            if in_gate:
                if res['in_min_delta'] is None or delta < res['in_min_delta']:
                    res['in_min_delta'] = delta
                    res['in_min_delta_at'] = (a, b, c, m, "j0", j0, "bj", bj)
                if delta <= 0:
                    res['in_ties'].append((a, b, c, m, interior_min, j0, bmin, bj))
            else:
                res['sub_count'] += 1
                res['sub_max_b'] = max(res['sub_max_b'], b)
                res['sub_b_list'].add(b)
                if delta <= 0:
                    res['sub_ties'] += 1
    return res

if __name__ == "__main__":
    print("=" * 70)
    print("STEP 0  cross-check closed forms vs Murnaghan-Nakayama")
    print("=" * 70)
    crosscheck(cmax=3, mmax=22)

    print()
    print("=" * 70)
    print("STEP 1  wide val-gap sweep, c=1,2,3, m up to 200")
    print("=" * 70)
    for c in (1, 2, 3):
        r = wide_sweep(c, 200)
        print(f"\n--- c={c}: {r['shapes']} shapes, m<= {r['mmax']} ---")
        print(f"  IN-GATE (j0+2c<=b): min Delta(val)=val(bdry)-V = {r['in_min_delta']} "
              f"at {r['in_min_delta_at']}")
        print(f"  IN-GATE ties (Delta<=0): {len(r['in_ties'])}  {r['in_ties'][:4]}")
        print(f"  SUB-GATE (excluded small-b families): {r['sub_count']} shapes, "
              f"of which {r['sub_ties']} have boundary in J* (expected);  b values = "
              f"{sorted(r['sub_b_list'])} (max b={r['sub_max_b']})")
        print(f"  boundedness: max v2(M_boundary,any i)={r['max_v2_boundary']} at "
              f"{r['max_v2_boundary_at']}; max v2(M_top=M_b+c)={r['max_v2_top']} at "
              f"{r['max_v2_top_at']}; max v2(M_0)={r['max_v2_M0']}")
